function [Ycorrect,max_err,frobenius] = algorithm(I,V,Y)
%{
    This function ingests a list (of length M) of current and voltage meausements
    (each measurement length N) generated from the true admittance matrix Y and
    outputs an estimate Ytilde of the admittance matrix computed from only the
    measurements themselves.
    
    Note that there are several threshhold parameters specified to defaults
    internally and in subfunctions.
%}

%{
    Version v2 differs from the original and has a modified ranking
function that enforces additional structure properties on Y, namely
    that the real parts of off-diagonal entries are non-positive and the
    imaginary parts of off-diagonal entries are non-negative. There are
    also additional outputs from the 'errors' function to report the
    percent error per entry in the recovered Ytilde.
%}

% Problem Parameters
%get the network size (number of rows and columns in Ytilde)
[~,N] = size(I);

Ytilde = complex(zeros(N,N)); %this matrix will be reduced in dimension at each iteration
Ycorrect = complex(zeros(N,N));  %this matrix will stay the same size and be filled in
fixedrows = []; %list will keep track of which rows have been fixed in Ycorrect



%initialize live plotting
figure(1)
drawnow

fprintf('\nIteration ')

for k = 1:N
    fprintf('%i, ',k)
    
    subplot(1,3,1)
    img1 = imagesc(abs(Y));
    set(img1,'Alphadata',~isnan(abs(Y)))
    title('Y')
    drawnow
    
    % Step 1 (Run Compressed Sensing On Each Row of Ytilde) %%
    
    [M,n] = size(I); %get number of measurements in current system
    
    %check if the number of measurements is larger than the system
    %size. If so, then pick a random subset of N-1 rows.
    
    if M > n
        idxs = randperm(M,n);
        
        %redefine I and V
        I = I(idxs,:);
        V = V(idxs,:);
    end
    
    
    [Ytilde] = CS(Ytilde,I,V);
    
    subplot(1,3,3)
    img2 = imagesc(abs(Ytilde));
    set(img2,'Alphadata',~isnan(abs(Ytilde)))
    title('Ytilde')
    drawnow
    
    
    % Step 2 (Evaluate Consistency of Ytilde and compute row rankings) %%
    [scores] = compute_scores(Ytilde,10); %c=10
    
    % Step 3 (Select Highest Rank) %%
    
    % select the **lowest??** ranked row
    [~,bestrow] = min(scores);
    
    % save "correct" row of Ytilde to final Ytilde (which is always N x N)
    [Ycorrect, fixedrows] = update_Ycorrect(Ycorrect, Ytilde, fixedrows, bestrow);
    
    subplot(1,3,2)
    img3 = imagesc(abs(Ycorrect));
    set(img3,'Alphadata',~isnan(abs(Ycorrect)))
    title('Ycorr')
    drawnow
    
    % Step 4 (System Reduction) %%
    
    % (i) update Ytilde (needs input of bestrow index)
    
    [Ytilde] = update_Ytilde(bestrow,Ytilde);
    
    % (ii) update current measurement matrix I (not the same as the original!)
    [I] = update_Itilde(bestrow,I,Ytilde,V);
    
    % (iii) reduce Ytilde (removes bestrow row AND column)
    [Ytilde] = reduce_matrix(Ytilde,bestrow,bestrow);
    
    % (iv) reduce I (updated-valued current measurement matrix) by removing bestrow column
    [I] = reduce_matrix(I,bestrow);
    
    % (v) reduce V (original-valued voltage measurement matrix) by removing bestrow column
    [V] = reduce_matrix(V,bestrow);
    
    % Step 5 (Loop Back!) %%
end

% Check the correctness of Ytilde against the ground truth Y
[frobenius,max_err] = errors(Y, Ycorrect);
% print()
fprintf('\nThe frobenius norm error is %d\n', frobenius)
fprintf('\nThe max entry deviation is %d percent\n', 100*max_err)

end %end CSalgorithm() function

%% Sub-functions of algorithm()

function [Ynew] = CS(Ytilde,I,V)
%{
    This function runs compressed sensing on all the rows in the (potentially reduced) Ytilde matrix and gets an updated Ytilde
%}

%------------------
Ynew = Ytilde;
N = length(Ynew);
idx = eye(N,N);
% compute a row of Y tilde using l-1 minimization
cvx_begin quiet
cvx_solver sedumi
variable x(N,N) complex symmetric
minimize(norm(x,1))
subject to

norm(I-V*x,2) <= 1e-4;
real(x(~idx)) <= 0;
imag(x(~idx)) >= 0;
real(diag(x)) >= 0;
% norm(I(:,i) - V*x,2) <= 1e-10;
%                         I(:,i) == V*x
%                         real(x([1:i-1 i+1:end])) <= 0;
%                         imag(x([1:i-1 i+1:end])) >= 0;
%                       real(x(i)) >= 0;
%I(:,i) == V*x %what is the tolerance on this equality
%constraint?
cvx_end
% store the row
Ynew= x;  %cvx seems to output the conjugate of the correct numbers for some solutions??

%------------------

end %end CS() function

function [scores] = compute_scores(Ytilde,c)
%{
    This function ingests the Ytilde computed by CS and evaluates its symmetric consistency.
    The output is a list of ranks of length len(Ytilde), with each element corresponding to a row of Ytilde.
    c is the constant that forces the score to be nonzero (default specified)
%}

%take care of default values
if nargin<2
    c = 10; %this is the parameter for term2
end

%------------------
N = length(Ytilde);  % number of rows in Ytilde
scores = zeros(N,1);  % pre-allocate scores (one for each row)

% weightings of score terms
lm1 = 1;
lm2 = 1;

% compute ranks of all rows of Ytilde
for i = 1:N  % rows
    for j = 1:N  % cols
        
        %symmetry
        term1 = abs(Ytilde(i,j) - Ytilde(j,i)); %could save this to not calculate twice??
        
        %sparsity (c is set to default value)
        term2 = 1.0/(abs(Ytilde(i,j))*abs(Ytilde(j,i)) + c);
        
        % add values from each column (j index) to score for the row
        scores(i) = scores(i) + lm1*term1 + lm2*term2;
    end
end
%------------------

end %end function

function [mag] = RealNeg_ImagPos(y)
%{
    Calculates the magnitude of the real part of the input ONLY IF the real
    part is greater than zero and sums it with the magnitude of
    the imaginary part of the input ONLY IF the imaginary part is less than
    zero. This rule means that if re(y) <=0 and im(y)>=0, then the function
    will return 0, which is what we want.
%}

%------------------
if (real(y) > 0 || imag(y) < 0)
    mag = abs(real(y)) + abs(imag(y));
else
    mag=0;
end
%------------------
end %end function

function [Ycorrect,fixedrows] = update_Ycorrect(Ycorrect,Ytilde,fixedrows,bestrow)
%{
    This function plugs in the values from the smaller Ytilde into the final estimate Ycorrect
%}

%------------------
N = length(Ycorrect);

proxyrows = [];

for row = 1:N
    if ~ismember(row,fixedrows)
        proxyrows = [proxyrows row];
    end
end

actualrowidx = proxyrows(bestrow);
fixedrows = [fixedrows actualrowidx];

%specify the row
counter = 1;
for j = 1:N
    if ismember(j,proxyrows)
        Ycorrect(actualrowidx,j) = Ytilde(bestrow,counter);
        counter = counter + 1;
    end
end

%specify the column now
Ycorrect(:,actualrowidx) = Ycorrect(actualrowidx,:);
%------------------
end %end function

function A = reduce_matrix(A,col,row)
%{
    This function takes in a matrix and removes the specified column (required argument) and removes the specified row (optional argument).
    A is a numpy array
    col is a non-negative integer
    row is a non-negative integer if it's specified, otherwise it's -1
%}

%take care of default values
if nargin<3
    row = -1;
end

%delete row if it's called for
if row ~=-1
    A(row,:) = [];
end

%always delete column
A(:,col) = [];

end %end function


function [Ytilde] = update_Ytilde(idx,Ytilde)
%{
    This function takes in the index of the "correct" row of Ytilde and fixes the corresponding entries
    Ytilde is a numpy array
    idx is an integer 0,...,N
%}
Ytemp = Ytilde;
Ytemp(idx,:) = Ytilde(idx,:); %enter row
Ytemp(:,idx) = Ytilde(idx,:).'; %enter col
Yiltde = Ytemp;

end %end function


function [Itilde] = update_Itilde(idx,I,Ytilde,V)
%{
    This function updates I to Itilde using the formula in step 4(ii) of the algorithm
    All input matrices are numpy arrays
    idx is the index of the to-be-removed row of Ytilde
%}

[~,N] = size(I);
Itilde = I;

for i =1:N %loop through columns of I
    Itilde(:,i) = I(:,i) - Ytilde(idx,i)*V(:,idx);
end

end %end function

function [frobenius,max_err] = errors(Ytruth, Ytilde)
%{
    This function the difference between the estimate Ytilde and the ground truth.
    Takes Frobenius norm of the difference between the two matrices. Also
    returns the "percent error" of the maximum deviation from the ground
    truth.
%}

% frobenius norm normalized by the number of matrix entries to get an
% "average" error per entry
frobenius = norm(Ytruth-Ytilde, 'fro')/(length(Ytruth)^2);

% max deviation entry
N = length(Ytruth); %system size
max_err = 0;
for i = 1:N
    for j = 1:N
        if Ytruth(i,j) ~=0 && abs(Ytruth(i,j) - Ytilde(i,j)) > max_err
            max_err = abs(Ytruth(i,j) - Ytilde(i,j))/Ytruth(i,j);
        end
    end
end


end %end function