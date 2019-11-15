function [x,errors,time] = CVX_opt_LIN(Idata,Vdata,Y)


tic; %start timing
N = length(Y);
idx = eye(N); %mask for entry-wise constraints

cvx_begin quiet
% cvx_solver gurobi %mosek
variable x(N,N) complex 
minimize(norm(x,1))
    subject to
        Idata == Vdata*x;
cvx_end
time = toc;

errors(1) = norm(Y - x,'fro')/(N^2); %normalized frobenius "per entry"
errors(2) = max(max(abs(Y - x))); %largest per-entry error
errors(3) = topology_error(Y,x); %number of topology errors, thresh = 1e-5
fprintf('cvx_optval = %d\n',cvx_optval)

end %end function


function [err3] = topology_error(Y,x)
zero_thresh = 1e-5; 
Y = logical(abs(Y) > zero_thresh);
x = logical(abs(x) > zero_thresh);
err3 = sum(sum(xor(x,Y)));
end %end function 