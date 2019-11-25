montecarlo = 3;

xis = linspace(1e-8,1e-2,3); %error scaling

casename = 'case24_ieee_rts';
mpc = ext2int(loadcase(casename)); %load case
Y = makeYbus(mpc); %Pre-compute the admittance matrix Y
N = length(Y);
errors = zeros(montecarlo,N,length(xis));

for m = 1:montecarlo
    for M=1:N
        for xi = 1:length(xis)
        
            % Get data (M samples)
            [Idata,Vdata,~] = simul_data(Y,M,xis(xi)); 

            % Run algorithm
            %[~,urrs,~] = CVX_opt(Idata,Vdata,Y); %OPT
            
            s = 2;
            Nsub = ceil(N/s); %use subset size of half
            [~,urrs,~] = algorithm_v3(Idata,Vdata,Y,Nsub,0,0,1); %OUR ALGORITHM
            
            errors(m,M,xi) = urrs(1);
            fprintf('mc = %i, M= %i\n, xi = %i\n',m,M,xi)
        end
    end
end
  

%% do averaging

min_samples = zeros(montecarlo,length(xis));

for mc = 1: montecarlo
    for xi = 1:length(xis)
        n=1;
        while n < 24 && errors(mc,n,xi)> .01
            n = n+1;
         end
        min_samples(mc,xi) = n;
    end
end
        

%average over MC
noisydata = mean(min_samples,1);


%% Plotting

figure (2)

plot(xis,noisydata)
%% Topology error plotting

mcerr = mean(errors,1);
figure(1)
for i = 1:N
    err = errors(1,i,:);
    hold on
    plot(xis,err(:),'Color',[(i/N),(i/N),(i/N)])
end
hold off

%%
i=14
err = errors(1,i,:);
plot(xis,err(:))
