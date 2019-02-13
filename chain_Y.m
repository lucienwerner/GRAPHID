function [Y] = chain_Y(n,plotting)

if nargin < 2
    plotting = 0;
end

N = linspace(1,n,n);
E = linspace(1,n-1,n-1);
Z_line = zeros(length(E),4);
Z_self = zeros(length(N),1); %only complex

% pick both real and imag values uniformly in [-100,100] for all line and self
% impedances, but set some equal to ~zero (1e-7)
magnitude=100;

% pick a subset from N to have non-zero Z_self
pzero = 0.1; %proportion of zeros
zs = randperm(n,ceil(n*pzero)); %bus indices (from N) that will be zero

% pick a subset from E to have non-zero Re(Z_line) (should be most)
pzero = 0.02; %proportion of zeros
zlr = randperm(n,ceil(n*pzero)); %line indices (from E) that will be zero real
zli = randperm(n,ceil(n*pzero)); %line indices (from E) that will be zero imag

for i=1:n
    Z_self(i) = (2*rand-1) * magnitude;
    if i > 1
        Z_line(i-1,1) = i-1;
        Z_line(i-1,2) = i;
        Z_line(i-1,3) = magnitude * rand;
        Z_line(i-1,4) = -magnitude * rand;
    end
end

Z_self(zs) = 1e-7;
Z_line(zlr,3) = 1e-7;
Z_line(zli,4) = -1e-7;

Y = zeros(n,n);

for e = 1:length(E)
    i = Z_line(e,1);
    j = Z_line(e,2);
    Y(i,j) = -1 * (Z_line(e,3) + 1i*Z_line(e,4));
    Y(j,i) = Y(i,j);
end

for i = 1:n
    Y(i,i) = 1i*Z_self(i) - sum(Y(i,:));
end
            
% Plotting
if plotting
    figure
    imagesc(abs(Y))
    title('Y')
end

end %end function
