%% generating data and parameters for the NDRO problem and calling the solver

%% data related to the empirical distribution

n = 250; % number of assets

T = 2*n; % data points

Xi = generate_data(n,T);


%% entropic risk parameters

theta = rand(n, 1) ;

regularizer = 0 ;

%% accuracy parameters

epsilon = 0.0125; %accuracy in duality gap and bound on worst case performance

%% wasserstein distance parameters

rho = 15; % radius of ambiguity set

c = 1; % wasserstein distance transportation cost exponent


%% running the FW algorithm to solve the DRO problem

[x, Q_points, p_prob, q_prob, min_F, sup_F, primal_sub_optimality, duality_gap, FW_gap, Ropt, Keps, conv_iter, e_flag_x, e_flag_FW] = entropic_risk_DRO(Xi,rho,c,theta,regularizer,epsilon) ;


% run plotting-n_pdf to plot the results and generate pdf figures
plot_n_pdf
