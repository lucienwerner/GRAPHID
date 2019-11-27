function [x,errors,time] = algorithm_v3(I,V,Y,N_sub,ew,lm1,lm2,gamma)
    plotting = 0;
    %{
    This function ingests a list (of length M) of current and voltage meausements 
    (each measurement length N) generated from the true admittance matrix Y and
    outputs an estimate Ytilde of the admittance matrix computed from only the
    measurements themselves.
    
    Note that there are several threshhold parameters specified to defaults 
    internally and in subfunctions.
    %}
    
    %{
    Version v3 differs from the original and has a modified ranking
    function that enforces additional structure properties on Y, namely
    that the real parts of off-diagonal entries are non-positive and the
    imaginary parts of off-diagonal entries are non-negative. There are
    also additional outputs from the 'errors' function to report the
    percent error per entry in the recovered Ytilde. 

    Also, made compatible with CVX_opt algorithm, can be run from
    pipeline_v1.m
    %}

    %Start Timing
    tic;

    % Problem Parameters 
    
    %get the network size (number of rows and columns in Ytilde)
    [~,N] = size(I);
    
    %set number of iterations of the algorithm
    %N_sub = 1; %size of subset of rows to pick in each iteration
    n_iter = floor(N/N_sub);
    if mod(N,N_sub)~=0
        n_iter = n_iter + 1;
    end
    
    Ytilde = complex(zeros(N,N)); %this matrix will be reduced in dimension at each iteration
    x = complex(zeros(N,N));  %this matrix will stay the same size and be filled in
    row_inds = 1:1:N; %starting row indices
   
    %initialize live plotting
    if plotting
        figure(1)
        subplot(1,3,1)
        img1 = imagesc(abs(Y));
        set(img1,'Alphadata',~isnan(abs(Y)))
        title('Y')
        drawnow
    end
   
    fprintf('\nSystem Order: ')
    
    for k = 1:n_iter
        [M,n] = size(I); %get size of current system
        
        fprintf('%i, ',n)
        if mod(k,10)==0
            fprintf('\n              ')
        end

        
        % Step 1 (Run Compressed Sensing To Get Each Row Of Ytilde) %

        %check if the number of measurements is larger than the system
        %size. If so, then pick a random subset of N-1 rows.
        
        if M > n
            idxs = randperm(M,n);
            
            %redefine I and V
            I = I(idxs,:);
            V = V(idxs,:);
        end    
        
        [Ytilde] = CSConstrained(I,V,ew,gamma); %entrywise constraints, need to set lm3 in scoring function to 0!
        
        if sum(isnan(Ytilde(:))) > 0
            time = Inf;
            errors(1:3) = Inf;
            return
        end

        % Step 2 (Evaluate Consistency of Ytilde and compute row rankings) %%
        
        [scores] = compute_scores(Ytilde,10,lm1,lm2); %c=10
        
        % Step 3 (Select Best Scoring) %
        
        % select the **lowest** scoring rows
        if n >= N_sub
            pickrows = N_sub;
        else
            pickrows = rem(N,N_sub);
        end
        
        [~,inds] = sort(scores,'ascend');
        bestrows = inds(1:pickrows);
       
        
        % save "correct" row(s) of Ytilde to final x (which is always N x N)
        for i=1:pickrows
            bestrow = bestrows(1);
            [x,row_inds] = update_x(x,Ytilde,bestrow,row_inds);
            
            if plotting
                subplot(1,3,2)
                img3 = imagesc(abs(x));
                set(img3,'Alphadata',~isnan(abs(x)))
                title('Ycorr')
                drawnow
            end
        
        % Step 4 (System Reduction) %
            
            % (i) update Ytilde (needs input of bestrows indices)
            if plotting
                subplot(1,3,3)
                img2 = imagesc(abs(Ytilde));
                set(img2,'Alphadata',~isnan(abs(Ytilde)))
                title('Ytilde')
                drawnow
            end
            [Ytilde] = update_Ytilde(bestrow,Ytilde);

            % (ii) update current measurement matrix I (not the same as the original!)
            [I] = update_Itilde(bestrow,I,Ytilde,V);

            % (iii) reduce Ytilde (removes bestrow row AND column)
            [Ytilde] = reduce_matrix(Ytilde,bestrow,bestrow);

            % (iv) reduce I (updated-valued current measurement matrix) by removing bestrow column
            [I] = reduce_matrix(I,bestrow);

            % (v) reduce V (original-valued voltage measurement matrix) by removing bestrow column
            [V] = reduce_matrix(V,bestrow);    
            
            % (vi) update remaining bestrows indices
            [bestrows] = update_bestrows(bestrows,bestrow);
        end
        
        % Step 5 (Loop Back!) %
        %plot empty Ytilde on the last iteration (just for visualization)
%         if k == n_iter
%             subplot(1,3,3)
%             img2 = imagesc(abs(Ytilde));
%             set(img2,'Alphadata',~isnan(abs(Ytilde)))
%             title('Ytilde')
%             drawnow
%             fprintf('Done ')
%         end
    end
    fprintf('\n')
    
    %End timing
    time = toc;
    
    errors(1) = norm(Y - x,'fro')/(N^2); %normalized frobenius "per entry"
    errors(2) = max(max(abs(Y - x))); %largest per-entry error
    errors(3) = topology_error(Y,x,1e-2); %number of topology errors, thresh = 1e-5

end %end algorithm

    
%% Sub-functions of algorithm()
    
function [Ytilde] = CS(I,V)
    %{
    This function runs compressed sensing on all the rows in the (potentially reduced) Ytilde matrix and gets an updated Ytilde
    %}

    [~,N] = size(V);
    Ytilde = zeros(N,N);
    
    for i = 1:N
        % compute a row of Ytilde using l1 minimization
        cvx_begin quiet 
            cvx_solver gurobi
            variable x(N) complex
            minimize(norm(x,1))
            subject to
                I(:,i) == V*x; %what is the tolerance on this equality?
        cvx_end

        % store the row
        Ytilde(i,:) = x.'; 
    end
    %------------------
    
end %end CS() function

function [Ytilde] = CSConstrained(I,V,ew,gamma)
    %{
    This function runs compressed sensing on all the rows in the (potentially reduced) Ytilde matrix and gets an updated Ytilde
    
    'ew' argument is a boolean to turn on/off the entrywise constraints in
    the cvx problem
    %}

    [~,N] = size(V);
    Ytilde = zeros(N,N);
    %gamma = 1e-4;
    for i = 1:N
        
        onehot = zeros(N,1);
        onehot(i) = 1;
        onehot = logical(onehot); %onehot encoding of each column
        
        % compute a row of Ytilde using l1 minimization
        cvx_begin quiet 
            cvx_solver sedumi
            variable x(N) complex
%             minimize(0.5*power(2,norm(I(:,i) - V*x,2)) + 0.2*norm(x,1))
            minimize(norm(x,1))
            subject to
                I(:,i) == V*x; %linear system
                %norm(I(:,i) - V*x,2) <= gamma; %LS
                if ew==1
                    real(x(~onehot)) <= 0;
                    imag(x(~onehot)) >= 0;
                    real(onehot) >= 0;
                end
        cvx_end

        % store the row
        Ytilde(i,:) = x.'; 
    end
    %------------------
    
end %end CSConstrained() function
    
function [scores] = compute_scores(Ytilde,c,lm1,lm2)
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
    %lm1 = 1;
    %lm2 = 1;
    lm3 = 0;
    
    % compute ranks of all rows of Ytilde
    for i = 1:N  % rows
        for j = 1:N  % cols
            
            %symmetry
            term1 = abs(Ytilde(i,j) - Ytilde(j,i)); %could save this to not calculate twice??
            
            %sparsity (c is set to default value)
            term2 = 1.0/(abs(Ytilde(i,j))*abs(Ytilde(j,i)) + c); 

            %real/imaginary signs
            if (i ~= j)
                term3 = RealNeg_ImagPos(Ytilde(i,j));
            else
                term3 = 0;
            end    
             
            % add values from each column (j index) to score for the row
            scores(i) = scores(i) + lm1*term1 + lm2*term2 + lm3*term3;
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

function [x,row_inds] = update_x(x,Ytilde,bestrow,row_inds)

%bestrow is in Ytilde indexing
%row_inds is in original indexing
%fixedrows is in original indexing

n = length(Ytilde);
 
actualrowind = row_inds(bestrow);

for j = 1:n
    x_ind = row_inds(j);
    x(actualrowind,x_ind) = Ytilde(bestrow,j);
    x(x_ind,actualrowind) = x(actualrowind,x_ind);
end

row_inds = setdiff(row_inds,row_inds(bestrow)); 

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

function [bestrows] = update_bestrows(bestrows,bestrow)

bestrows = bestrows(2:end); %remove bestrow (i.e., first element) from list

%update remaining indices
for j = 1:length(bestrows)
    if bestrows(j) > bestrow
        bestrows(j) = bestrows(j)-1;
    end
end
end %end function

function [err3] = topology_error(Y,x,zero_thresh)

if nargin < 3
    zero_thresh = 1e-5; 
end

Y = logical(abs(Y) > zero_thresh);
x = logical(abs(x) > zero_thresh);
err3 = sum(sum(xor(x,Y)));
end %end function 

