function [Idata,Vdata,Y] = simul_data(Y,M,xi,plotting)

%{
Function generates simulated (random) complex-valued data for a given  Y matrix
%}
if nargin < 3
    xi = 0; %noise set to zero by default
end

if nargin < 4
    plotting =0; %plotting off by default
end

Vlims = [0.95,1.05]; 
n = length(Y);

% Generate random voltage data within PU limits. Centered around 1.0pu
Vdata = (Vlims(2)-Vlims(1)) * rand(M,n) + Vlims(1);

% Generate mean-zero uniform random additive noise Z. Size scaled by xi
Zdata = xi * (rand(M,n) - 0.5); %uniform
Zdata = xi * randn(M,n);

%Compute currents
Idata = Vdata*Y + Zdata;

% Plotting
if plotting
    figure
    subplot(1,3,1)
    imagesc(abs(Y))
    title('Y')
    subplot(1,3,2)
    imagesc(abs(Idata))
    title('I')
    xlabel('Bus')
    ylabel('Datapoint')
    subplot(1,3,3)
    imagesc(abs(Vdata))
    title('V')
    xlabel('Bus')
    ylabel('Datapoint')
end

end %end function