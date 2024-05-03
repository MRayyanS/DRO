function [Q,XI] = generate_data_ellipsoidal_support(n,N)


% generate the data
rng('shuffle');
XI = randn(n,N) ;  

%  The the variance matrix randomly 
rng('shuffle');
U = rand(n,n+10);
U = orth(U);
spectrum = rand(n,1) ;

Var = U'*diag(spectrum)*U;
sqrtVar = sqrtm(Var);

% comoute a random mean and add 
mu = randn(n,1) + 0.5.*ones(n,1);

% transform the data for correct variance and mean
XI = sqrtVar*XI + mu*ones(1,N);

%% Characterising the ellipsoid and normalising the data

% Generate the Q matrix
rng('shuffle');
U = rand(n,n+10);
U = orth(U);
spectrum = rand(n,1) + 0.25.*ones(n,1) ;

Q = U'*diag(spectrum)*U;

normalise_factor = 0.125 ;

Q = normalise_factor.*Q ;


% normalise the data to fit in the ellipsoid
normalisation_vec = diag(XI'*Q*XI) ;

randon_norm_vec = ones(N,1) ; 

for i = 1:N
    XI(:,i) = (randon_norm_vec(i)/sqrt(normalisation_vec(i))).*XI(:,i) ; % normalised data points
end






end