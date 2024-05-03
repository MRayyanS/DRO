function [Xi_perturbed, EQ_expterm, e_flag_FW] = FW_oracle_entropic_risk(Xi, c, rho, x, theta)

n = size(x,1);

Xi_perturbed = zeros(size(Xi)) ;
EQ_expterm   = zeros(n,1) ;
e_flag_FW    = 0 ; 

for j = 1:n
    % individually solving the problem for each index j = 1,2,...,n    
    thetaj = theta(j) ;
    xj     = x(j) ;
    Xij    = Xi(j,:) ;

    [Xi_perturbed(j,:), EQ_expterm(j), e_flag_FW_j]  = FW_oracle_scalar_problem(thetaj, xj, Xij, c, rho) ;
    e_flag_FW = e_flag_FW + e_flag_FW_j ;
end

end

function [xi_perturbed, EQ_expterm, e_flag_FW] = FW_oracle_scalar_problem(theta, x, xi, c, rho)

T = size(xi, 2) ;

r_factor = (c - theta*x)/c ; % occurs multiple times
etocrho = exp(c*rho) ;
thetax = theta*x ;

if thetax == 0
    xi_perturbed = xi ;
    EQ_expterm   = 1 ;
    e_flag_FW    = 0 ;

else
    exp_seq    = exp(-(thetax/r_factor).*xi) ;
    sorted_seq = sort(exp_seq, 'ascend') ;
    tailsum    = norm(sorted_seq,1) ; 

 for t = 1:T
    tailsum   = tailsum - sorted_seq(t) ;
    check_sum = tailsum/(T*etocrho - t) ;
    if exp_seq(t) >= check_sum % check for der. w.r.t eta being >= 0
       % if yes, definining eta_opt
       eta_opt = (thetax / c)*check_sum^r_factor ;
       break ;
    end
 end

% defining optimal perturbed points
xi_perturbed = zeros(1,T) ;

% delete this later
exp_perturbation = zeros(1, T) ;

logterm = log( (c*eta_opt) / (thetax)) ;
for t = 1:T
    second_term = ( c*xi(t) + logterm )/(c - thetax) ;
    xi_perturbed(t) = min( xi(t) , second_term ) ;

    exp_perturbation(t) = exp( c*(xi(t) - xi_perturbed(t)) ) ;
end

% computing expectation w.r.t. Q - the perturbed points
dummy = exp( -thetax.*xi_perturbed )*ones(T,1) ;
EQ_expterm = (dummy/T) ;

e_flag_FW = mean(exp_perturbation,'all') - exp(c*rho) ;

end

end