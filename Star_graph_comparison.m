<<<<<<< Updated upstream
montecarlo = 1000;

N=25; %number of nodes

errors1 = zeros(montecarlo,N,3);
errors2 = zeros(montecarlo,N,3);
errors3 = zeros(montecarlo,N,3);
errors4 = zeros(montecarlo,N,3);

=======
montecarlo = 1;

N=24; %number of nodes
Nsub=12; %size of subset
errors1 = zeros(montecarlo,N,3);
errors2 = zeros(montecarlo,N,3);
errors3 = zeros(montecarlo,N,3);
>>>>>>> Stashed changes

for m = 1:montecarlo

    %Get Y matrix for the graph
    Y = star_Y(N);
    for M=1:N
        % Get data (M samples)
<<<<<<< Updated upstream
        [Idata,Vdata,~] = simul_data(Y,M); 
        % Run algorithms
        [~,errors1(m,M,:),~] = CVX_opt_LIN(Idata,Vdata,Y); %only linear system
        fprintf('mc = %i, M= %i\n',m,M)
    end
end
  

%do averaging
err1_1 = mean(errors1(:,:,1), 1);
err2_1 = mean(errors2(:,:,1), 1);
err3_1 = mean(errors3(:,:,1), 1);
err4_1 = mean(errors4(:,:,1), 1);

err1_2 = mean(errors1(:,:,2), 1);
err2_2 = mean(errors2(:,:,2), 1);
err3_2 = mean(errors3(:,:,2), 1);
err4_2 = mean(errors4(:,:,2), 1);

err1_3 = mean(errors1(:,:,3), 1);
err2_3 = mean(errors2(:,:,3), 1);
err3_3 = mean(errors3(:,:,3), 1);
err4_3 = mean(errors4(:,:,3), 1);

%% Plotting

figure (2)

subplot(1,3,1)
plot(linspace(1,N,N),err1_1,linspace(1,N,N),err2_1,linspace(1,N,N),err3_1,linspace(1,N,N),err4_1)

subplot(1,3,2)
plot(linspace(1,N,N),err1_2,linspace(1,N,N),err2_2,linspace(1,N,N),err3_2,linspace(1,N,N),err4_2)

subplot(1,3,3)
plot(linspace(1,N,N),err1_3,linspace(1,N,N),err2_3,linspace(1,N,N),err3_3,linspace(1,N,N),err4_3)


=======
        % Noiseless
        [Idata,Vdata,~] = simul_data(Y,M); 
        % Run algorithms
%         [~,errors1(m,M,:),~] = CVX_opt_LIN(Idata,Vdata,Y); %only linear system
        [~,errors2(m,M,:),~] = CVX_opt_SYM(Idata,Vdata,Y); %linear system + symmetry
%         [~,errors3(m,M,:),~] = algorithm_v3(Idata,Vdata,Y,Nsub,0,1,0); %OUR ALGORITHM
        fprintf('mc = %i, M= %i\n',m,M)
    end
end


%% Plotting

%do averaging
err1 = mean(errors1(:,:,1), 1);
err2 = mean(errors2(:,:,1), 1);
err3 = mean(errors3(:,:,1), 1);
%% Plotting

figure



plot(linspace(1,N,N),err2)
% plot(linspace(1,N,N),err1, linspace(1,N,N),err2, linspace(1,N,N),err3)
% legend('Linear','Symmetry','Iterative')
>>>>>>> Stashed changes




