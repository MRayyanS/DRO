function [x, risk_val, EP_expterm, e_flag] = solve_for_x(theta, Xi, Q_points, p_prob, q_prob, x0, regularizer, current_iter)

max_iter         = 100;
e_flag_threshold = 0.005 ;

stepbase = 0.025;

% initialization of x if not taken as input
% x0 = (1/n)*ones(n,1) ;

% initialization for accelerated algrithm
y = x0;
xprev = x0 ;
tprev = 1 ;

for s = 1:max_iter
    %% computing the gradients

    [grad, risk_val, EP_expterm] = compute_ER_gradient(y, theta, Xi, Q_points, p_prob, q_prob, regularizer, current_iter) ;

    %% FISTA
    % stepsize at iteration t
    stepsize = (stepbase*10)/(10+s) ;

    % gradient descent step
    descent_point = y - stepsize*grad ;

    % FISTA steps
    proj_simplex_array = @(y) max(bsxfun(@minus,y,max(bsxfun(@rdivide,cumsum(sort(y,1,'descend'),1)-1,(1:size(y,1))'),[],1)),0) ;
    x = proj_simplex_array(descent_point) ;

    % compute exit_flag
    e_flag = norm(y - x, 2) ;

    t = ( 1 + sqrt(1 + 4*tprev^2) )/2 ;
    y = x + ((tprev - 1)/t)*(x - xprev)  ;

    % update xprev and tprev for next interates
    xprev = x ;
    tprev = t ;

    % break from the for loop if e_flag is small enough
    if e_flag <= e_flag_threshold
        break ;
    end
end

% end of function solve_for_x(theta, Xi, Q_points, p_prob, q_prob, x0, regularizer, current_iter)
end

function [grad, risk_val, EP_expterm] = compute_ER_gradient( y, theta, Xi, Q_points, p_prob, q_prob, regularizer, current_iter )

[n,T] = size(Xi) ;
grad = zeros(size(y)) ;

risk_val = 0 ;
EP_expterm = zeros(n,1) ;

for j = 1:n

    thetaj = theta(j) ;
    yj     = y(j) ;
    thetayj = thetaj*yj ;

    gradj_num = 0 ;
    gradj_den = 0 ;

    for t = 1:T
        % compute the num and den terms for each t
        exp_term_t = exp( - thetayj*Xi(j,t) ) ;
        num_term_t = p_prob*Xi(j,t)*exp_term_t ;
        den_term_t = p_prob*exp_term_t ;

        % computing sum over l
        for l = 1:(current_iter-1)
            exp_term_l = exp( - thetayj*Q_points(j,t,l) ) ;

            num_term_t = num_term_t + q_prob(l)*Q_points(j,t,l)*exp_term_l ;
            den_term_t = den_term_t + q_prob(l)*exp_term_l ; 
        end

        % update teh grad terms for each t
        gradj_num = gradj_num + num_term_t ;
        gradj_den = gradj_den + den_term_t ;
    end

    % compute the jth componenet of the gradient with regularizer
    gradj = (regularizer*yj) - (gradj_num/gradj_den) ;
    grad(j) = gradj ;

    % the expectation w.r.t. P of the exponent term is the denominator of the gradient
    EP_expterm(j) = (1/T)*gradj_den ;
    
    risk_val      = risk_val + ( log( EP_expterm(j) ) )/(thetaj) ;

end

end
