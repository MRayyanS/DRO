function [mu_FW, sigma_FW, FW_gap] = FWOracle_min_variance_ellipsoidal_suuport(Q, XI, rho, x, v, N, n, var)
 
%% Defining all relevant matrices for the SDP 

% Definig Mx
Mx = [-x*x' (x'*v).*x ; (x'*v).*x' -(x'*v)^2];

% Defining Mlambda
Mlambda = [Q zeros(n,1); zeros(1,n) -1];

% Defining Mgamma
Mgamma = [zeros(n,n) zeros(n,1); zeros(1,n) -1];

%Definig matrices Meta(i)
Meta = zeros(n+1,n+1,N);
for i = 1:N
    Meta(:,:,i) = [eye(n) -XI(:,i); -XI(:,i)' norm(XI(:,i),2)^2] ;
end

%% Deining other required quantites

% definig the e vector
e = (1/N).*ones(N,1);





%% SDP begins
cvx_begin sdp

cvx_quiet true

variable eta
variable gam(N,1)
variable lambda_mul(N,1)

J = eta*rho^2 - e'*gam ;

minimize (J)

subject to

eta >= 0 ;
lambda_mul >= 0 ;


for i = 1:N
    Mx + lambda_mul(i).*Mlambda + gam(i).*Mgamma + eta.*Meta(:,:,i) >= 0 ;
end

cvx_end




% computing q(eta*, xi_i) and FW points

FW_points = zeros(n,N);

A0 = eta.*eye(n) - x*x' ;
A1 = Q;

for i = 1:N

    b0 = (x'*v).*x - eta.*XI(:,i);

    FW_points(:,i) = -(A0 + lambda_mul(i).*A1)\b0;

end

% Computing mean and variances of FW_points

mu_FW = (1/N).*FW_points*ones(N,1) ;
sigma_FW = (1/N).*(FW_points*FW_points') ;

% Computing the FW gap

FW_gap = J - x'*var*x ;   % note: no x_regulariser*(x'*x) here








end
