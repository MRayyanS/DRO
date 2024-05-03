%% generating data and parameters for the NDRO problem and calling the solver

function [xout, sigmaout, muout, x,sigma,mu,Ropt, data,Q,rho, x_regulariser, Keps, epsilon, min_F, sup_F, primal_sub_optimality, duality_gap, FW_gap] = portfolio_optimization_DRO_ellipsoidal_support(Q,data)

%% wasserstein distance parameters

rho = 1.5 ; % radius of ambiguity set

% m = 2; % order of wasserstein distance

% p = 2; % norm used for wasserstein distance

x_regulariser = 0.3 ; % regulariser if no strong convexity

%% accuracy parameters

epsilon = 0.25; %accuracy in duality gap and bound on worst case performance

% delta = 0 ; % accuracy of FW-oracle


%% data related to the empirical distribution

 n = 25; % number of assets

 N = 2*n; % data points

 %[Q,data] = generate_data_ellipsoidal_support(n,N);


sigma0 = (1/N).*data*data';
mu0 = (1/N).*data*ones(N,1);

%% equidistribution startegy as a  base strategy

% x_equi = (1/n).*ones(n,1);

% alpha_min = x_equi'*mu0 - rho*norm(x_equi, 2) ;
% alpha_min = 0.5;

%% running the FW algorithm to solve the DRO problem

[xout,sigmaout,muout,x, sigma, mu, min_F, sup_F, primal_sub_optimality, duality_gap, FW_gap, Ropt, Keps] = min_variance_DRO_ellipsoidal_support(data, Q, rho, epsilon, x_regulariser, sigma0, mu0) ;

%% plotting the results

k = 1:2*Keps+2 ;

figure(1);
plot(k,min_F,'-o','MarkerSize',4);
hold on
grid on;
plot(k,sup_F,'-*','MarkerSize',4);
axis padded;
legend('$\ \min\limits_{x \in \Delta^{n}} \ F(x,P_k)$', '$\sup\limits_{P \in {W}_2 (\widehat{P},\rho) } F \big( x_k,P \big)$','FontSize',15,'Interpreter','latex','Location','southeast');
xlabel('Iterations, $k$', 'FontSize',15,'Interpreter','latex');
hold off
% Get current figure handle
fig = gcf;
% Synchronise units for screen and paper
fig.Units      = 'centimeters';
fig.PaperUnits = 'centimeters';
% Get width and heigth of figure (in centimeters)
fig_width  = fig.Position(3);
fig_height = fig.Position(4);
% Activate auto-positioning of figure on paper
fig.PaperPositionMode = 'auto';
% Set paper size to figure size
fig.PaperSize = [fig_width, fig_height];


figure(2);
semilogy(k,primal_sub_optimality,'-o','MarkerSize',4);
hold on;
grid on;
semilogy(k,duality_gap,'-*','MarkerSize',4);
axis padded;
legend('Primal sub optimality', 'Duality gap','FontSize',12,'Interpreter','latex','Location','northeast');
xlabel('Iterations, $k$', 'FontSize',15,'Interpreter','latex');
hold off;
% Get current figure handle
fig = gcf;
% Synchronise units for screen and paper
fig.Units      = 'centimeters';
fig.PaperUnits = 'centimeters';
% Get width and heigth of figure (in centimeters)
fig_width  = fig.Position(3);
fig_height = fig.Position(4);
% Activate auto-positioning of figure on paper
fig.PaperPositionMode = 'auto';
% Set paper size to figure size
fig.PaperSize = [fig_width, fig_height];












end






