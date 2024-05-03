function worst_case_cost = compute_worst_case_cost(x, mu0, sigma0, Q, XI, rho, N, n, var)

max_iter = 20;

sigma = sigma0;
mu = mu0;


for k = 0:max_iter
    stepsize = 2/(k + 2);
  
% Solving the FW linear minimization problem

     v = mu;
     
    [mu_FW, sigma_FW] = FWOracle_min_variance_ellipsoidal_suuport(Q, XI, rho, x, v, N, n, var);


    % FW update
    sigma = sigma + stepsize.*(sigma_FW - sigma);
    mu = mu + stepsize.*(mu_FW - mu) ;

end

var = sigma - mu*mu' ;
worst_case_cost = x'*var*x ;

end