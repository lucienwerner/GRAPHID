function [Idata,Vdata,Y] = simul_data(Y,M,plotting)

%{
Function generates simulated (random) complex-valued data for a given  Y matrix
%}
if nargin < 3
    plotting =0; %plotting off by default
end

Vlims = [0.95,1.05]; 
n = length(Y);

Vdata = (Vlims(2)-Vlims(1)) * rand(M,n) + Vlims(1);
Idata = Vdata*Y;

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