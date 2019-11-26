montecarlo = 50; %# of montecarlo iterations

xis = [0,1e-9,1e-7,1e-6,1e-5,1e-3,1e-2]; %error variance

casename = 'case24_ieee_rts';
mpc = ext2int(loadcase(casename)); %load case
Y = makeYbus(mpc); %Pre-compute the admittance matrix Y
N = length(Y);
num_errors = 3;
errors = zeros(montecarlo,N,length(xis),num_errors);

for m = 1:montecarlo
    for M=1:N
        for xi = 1:length(xis)
        
            % Get data (M samples)
            [Idata,Vdata,~] = simul_data(Y,M,xis(xi)); 

            % Run algorithm
            %[~,urrs,~] = CVX_opt(Idata,Vdata,Y); %OPT
            gamma = sqrt(N) * sqrt(xi);
            s = 2;
            Nsub = ceil(N/s); %use subset size of half
            [~,errors(m,M,xi,:),~] = algorithm_v3(Idata,Vdata,Y,Nsub,0,1,0,gamma); %OUR ALGORITHM
            fprintf('mc = %i, M= %i\n, xi = %d\n',m,M,xis(xi))
        end
    end
end
  

%% do averaging

min_samples = zeros(50,length(xis));

for mc = 1: montecarlo
    for xi = 1:length(xis)
        n=1;
        while n < N && myerrors(mc,n,xi,1)> .7
            n = n+1;
         end
        min_samples(mc,xi) = n;
    end
end
        

%average over MC
noisydata2 = mean(min_samples,1);
figure (2)
semilogx(xis,noisydata2)
xlabel('Noise Variance')
ylabel('Number of samples ($m$)', 'Interpreter','latex')

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

load('errors.mat') %called 'myerrors'
fro_raw = myerrors(1:50,:,:,1);
max_raw = myerrors(1:50,:,:,2);
top_raw = myerrors(1:50,:,:,3);

%do averaging over MC iterations
fro_errors = mean(fro_raw,1);
max_errors = mean(max_raw,1);
top_errors = mean(top_raw,1);


figure(1)
colors = linspace(0,.8,length(xis));
for i = 1:length(xis)
    hold on
    plot(fro_errors(1,:,i),'-+','Color',[colors(i),colors(i),colors(i)],'LineWidth',3,'MarkerSize',8)
end
hold off
xlim([1,24])
ylim([0.35,0.9])
xlabel('Samples ($m$)','Interpreter','latex')
ylabel('$\frac{1}{n^2}||$\boldmath$X -  $\boldmath$Y||$','Interpreter','latex')
%ylabel('Num topology errors','Interpreter','latex')
legend({'$\sigma_N^2 = 0$','$\sigma_N^2 = 10^{-9}$','$\sigma_N^2 = 10^{-7}$','$\sigma_N^2 = 10^{-6}$','$\sigma_N^2 = 10^{-5}$','$\sigma_N^2 = 10^{-3}$','$\sigma_N^2 = 10^{-2}$'},'Location','southwest','NumColumns',1,'Interpreter','latex')
set(gcf,'color','w')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
grid on
box on
set(gcf,'position',[700 900 800 500])
hfig= tightfig(gcf);


    
