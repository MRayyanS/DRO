function [xout, min_F ] = solve_for_x(var, x_regulariser, n)

e = ones(n,1);


%% optimization solver to compute minimizing x
cvx_begin

cvx_quiet true

variable x(n,1)

 J = x'*var*x + 0.5*x_regulariser*(x'*x);  


minimize (J)

subject to

% x in simplex
x >= 0;
e'*x == 1;

cvx_end

%% solver ends
% output the minimizing x
xout = x ;

%% optimization solver to compute minimum over x without regularizer
cvx_begin

cvx_quiet true

variable x(n,1)

 J = x'*var*x ;  


minimize (J)

subject to

% x in simplex
x >= 0;
e'*x == 1;

cvx_end

%% solver ends

%output the min over x value (without regularizer) for plotting purpose
min_F = J ;


end