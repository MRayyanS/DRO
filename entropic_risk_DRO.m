function [x, Q_points, p_prob, q_prob, min_F, sup_F, primal_sub_optimality, duality_gap, FW_gap, Ropt, Keps, conv_iter, e_flag_x, e_flag_FW] = entropic_risk_DRO(Xi,rho,c,theta,regularizer,epsilon)
[n,T] = size(Xi) ;

Keps      = 350 ;
conv_iter = 2*Keps+1 ; % a dummy conv_iter place holder
% Keps = ceil(10*rho/epsilon) - 2;

max_iter   = Keps ;
extra_iter = Keps ;
total_iter = max_iter + extra_iter ;

% quantities to plot
min_F                 = zeros(2*Keps + 1, 1);
FW_gap                = zeros(2*Keps + 1, 1);
sup_F                 = zeros(2*Keps + 1, 1);
duality_gap           = zeros(2*Keps + 1, 1);
primal_sub_optimality = zeros(2*Keps + 1, 1);

e_flag_FW             = zeros(2*Keps + 1, 1) ;
e_flag_x              = zeros(2*Keps + 1, 1) ;

%% initialization for minimization over x for the first iteration, it updates itself for later iterations
%x0 = (1/n).*ones(n,1);
x0 = rand(n,1);
x0 = (1/norm(x0,1)).*x0;

x = x0;

%% initialization for FW procedure
p_prob = 1 ;
q_prob = zeros(total_iter,1) ;

Q_points = zeros(n,T,total_iter) ;


%% Frank-Wolfe procedure for first regime of diminishing step-size

for k = 1:Keps+1
    stepsize = 2/(k + 1) ; % since the iterate starts from 1 and not 0
    
%% The two oracles

  % solving minimization over x at kth iteration
    [x, min_F(k), EP_expterm, e_flag_x(k)] = solve_for_x(theta, Xi, Q_points, p_prob, q_prob, x, regularizer, k) ;
    
  % Solving the FW linear minimization problem
    [Xi_perturbed, EQ_expterm, e_flag_FW(k)] = FW_oracle_entropic_risk(Xi, c, rho, x, theta) ;

%% Computing different gaps
  % computing the FW gap
    FW_gap(k) = compute_FW_gap(EP_expterm, EQ_expterm, theta) ;

  % computing worst case cost for the current decision x
  % the FW optimal distribution is the worst-case distribution for the
  % non-linear problem as well

    sup_F(k) = (1./theta)'*log(EQ_expterm) ;

  % computing duality-gap
    duality_gap(k) = sup_F(k) - min_F(k) ;
 
%% Update of the Fw step
  % storing all the perturbed points
    Q_points(:,:,k) = Xi_perturbed ;

  % update of Probabilities
    p_prob    = (1 - stepsize)*p_prob ;
    q_prob    = (1 - stepsize)*q_prob ; % scales all the previous probabilities by a factor of (1 - stepsize_k)
    q_prob(k) = stepsize ;

end


%% Frank-Wolfe procedure for second regime of constant step-size

for k = Keps+2:2*Keps+1
    stepsize = 2/(Keps + 1) ; % constant stepsize
  
%% The two oracles
  % solving minimization over x at kth iteration
    [x, min_F(k), EP_expterm, e_flag_x(k)] = solve_for_x(theta, Xi, Q_points, p_prob, q_prob, x, regularizer, k) ;
    
  % Solving the FW linear minimization problem
    [Xi_perturbed, EQ_expterm, e_flag_FW(k)] = FW_oracle_entropic_risk(Xi, c, rho, x, theta) ;

%% Computing different gaps
  % computing the FW gap
    FW_gap(k) = compute_FW_gap(EP_expterm, EQ_expterm, theta) ;

  % computing worst case cost for the current decision x
    sup_F(k) = (1./theta)'*log(EQ_expterm) ;

  % computing duality-gap
    duality_gap(k) = sup_F(k) - min_F(k) ;

  %% Chech for convergence, if not update of the FW step
  % if the FW_gap is small enough, then break
    if FW_gap(k) <= epsilon
        conv_iter = k ;
    end
    
  % if the FW_gap is not small enough, then update the iterates
  % storing all the perturbed points
    Q_points(:,:,k) = Xi_perturbed ;

  % update of Probabilities
    p_prob    = (1 - stepsize)*p_prob ;
    q_prob    = (1 - stepsize)*q_prob ; % scales all the previous probabilities by a factor of (1 - stepsize_k)
    q_prob(k) = stepsize ;

end


%% running the extra FW iterations to finally compute the optimal value and a saddle point by running "extra_iter" many iterations

for k = 1:extra_iter
    stepsize = 2/(k + 1) ; % since the iterate starts from 1 and not 0
  
%% The two oracles
  % solving minimization over x at kth iteration
    [x, Ropt] = solve_for_x(theta, Xi, Q_points, p_prob, q_prob, x, regularizer, k) ;
    
  % Solving the FW linear minimization problem
    [Xi_perturbed] = FW_oracle_entropic_risk(Xi, c, rho, x, theta) ;

  %% Update of the Fw step
  % storing all the perturbed points
    Q_points(:,:,2*Keps+1+k) = Xi_perturbed ;

  % update of Probabilities
    p_prob    = (1 - stepsize)*p_prob ;
    q_prob    = (1 - stepsize)*q_prob ; % scales all the previous probabilities by a factor of (1 - stepsize_k)
    q_prob(k) = stepsize ;

end


%% computing primal sub-optimality using Ropt computed from extra iterations

for k = 1:2*Keps + 1
    primal_sub_optimality(k) = Ropt - min_F(k);
end

% the entropic_risk_DRO() function ends
end



%% function for computing FW_gap

function FW_gap = compute_FW_gap(EP_expterm, EQ_expterm, theta)
n = size(EP_expterm, 1) ;
FW_gap = 0 ;
for j = 1:n
    FW_gap_j = ( EQ_expterm(j) - EP_expterm(j) )/( theta(j)*EP_expterm(j) ) ;
    FW_gap   = FW_gap + FW_gap_j ;
end
end