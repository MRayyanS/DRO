function [xout,sigmaout,muout,x, sigma, mu, min_F, sup_F, primal_sub_optimality, duality_gap, FW_gap, Ropt, Keps] = min_variance_DRO_ellipsoidal_support(XI, Q, rho, epsilon, x_regulariser, sigma0, mu0)
[n,N] = size(XI);

Keps = compute_Keps(epsilon,rho);

min_F = zeros(2*Keps + 2, 1);
FW_gap = zeros(2*Keps + 2, 1);
sup_F = zeros(2*Keps + 2, 1);
duality_gap = zeros(2*Keps + 2, 1);
primal_sub_optimality = zeros(2*Keps + 2, 1);



%% dummy placeholder values for xout, sigmaout
sigmaout = zeros(n,n);
muout = zeros(n,1);
xout = zeros(n,1);


%% initialization for minimization over x for the first iteration, it updates itself for later iterations

%x0 = (1/n).*ones(n,1);
x0 = rand(n,1);
x0 = (1/norm(x0,1)).*x0;

%% initialization for FW algorithm

x = x0;

sigma = sigma0;
mu = mu0;

%% Frank-Wolfe algorithm for first regime of diminishing step-size

for k = 0:Keps
    
    stepsize = 2/(k + 2);

    var = sigma - mu*mu';
    
  % solving minimization over x at kth iteration  
    
    [x, min_F(k+1) ] = solve_for_x(var, x_regulariser, n);

  % computing worst case cost for the current decision x
    
    sup_F(k+1) = compute_worst_case_cost(x, mu, sigma, Q, XI, rho, N, n, var) ;
  
  % computing Duality gap
    
    duality_gap(k+1) = sup_F(k+1) - min_F(k+1);
   
  % solving the FW linear minimization problem

     v = mu;
     
    [mu_FW, sigma_FW, FW_gap(k+1)] = FWOracle_min_variance_ellipsoidal_suuport(Q, XI, rho, x, v, N, n, var);

  % FW update
  
    sigma = sigma + stepsize.*(sigma_FW - sigma);
    mu = mu + stepsize.*(mu_FW - mu) ;

end

%% Frank-Wolfe procedure for second regime of constant step-size
for k = Keps+1:2*Keps+1
%   k = Keps+1:2*Keps+1

    stepsize = 2/(Keps + 2);

    var = sigma - mu*mu';

  % solving minimization over x at kth iteration  
   
    [x, min_F(k+1) ] = solve_for_x(var, x_regulariser, n);

  % computing worst case cost for the current decision x

    sup_F(k+1) = compute_worst_case_cost(x, mu, sigma, Q, XI, rho, N, n, var) ;
    
  % computing Duality gap

    duality_gap(k+1) = sup_F(k+1) - min_F(k+1);
  
  % Solving the FW linear minimization problem

     v = mu;
     
    [mu_FW, sigma_FW, FW_gap(k+1)] = FWOracle_min_variance_ellipsoidal_suuport(Q, XI, rho, x, v, N, n, var);

  % if the FW_gap is small enough, then define the outputs

    if FW_gap(k+1) <= epsilon
       sigmaout = sigma;
       muout = mu;
       xout = x;
    end
    
  % if the FW_gap is not small enough, then update the iterates

    sigma = sigma + stepsize.*(sigma_FW - sigma);
    mu = mu + stepsize.*(mu_FW - mu) ;

end

%% running the FW iterations to finally compute the optimal value and a saddle point by running "extra_iter" many iterations

extra_iter = Keps;
for k = 0:extra_iter
    stepsize = 2/(Keps+ k + 2);

    var = sigma - mu*mu';

  % solving minimization over x at kth iteration  
    
    [x, Ropt] = solve_for_x(var, x_regulariser, n);
   
%% Solving the FW linear minimization problem

     v = mu;
     
    [mu_FW, sigma_FW] = FWOracle_min_variance_ellipsoidal_suuport(Q, XI, rho, x, v, N, n, var);

    sigma = sigma + stepsize.*(sigma_FW - sigma);
    mu = mu + stepsize.*(mu_FW - mu) ;

end



%% computing primal sub-optimality using Ropt computed from extra iterations

for k = 1:2*Keps + 2
    primal_sub_optimality(k) = Ropt - min_F(k);
end

end





function Keps = compute_Keps(epsilon,rho)
%Keps = ceil(10*rho/epsilon) - 2;
Keps = 75;
end