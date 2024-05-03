%% plotting the results

k = 1:2*Keps+1 ; % iteration index


%% plotting the primal and dual functions

figure(1)
fig_distance = tiledlayout(1,1,'Padding','tight');
fig_distance.Units = 'inches';
fig_distance.OuterPosition = [0.25 0.25 5 5];
nexttile;

plot(k,min_F,'-o','MarkerSize',3);
hold on
plot(k,sup_F,'-x','MarkerSize',3);
grid on;
axis padded;
xlabel('Iterations, $k$', 'FontSize',10,'Interpreter','latex');
legend('$ \min\limits_{x \in \Delta^{n}} \mathcal{E} (x,{P}_k)$', '$\sup\limits_{ {P} \in \mathcal{P} } \, \mathcal{E} ( x_k,P )$','FontSize',18,'Interpreter','latex','Location','southeast');
im = gcf;
exportgraphics(im,'pd_func_rhoalpha0.pdf','ContentType','vector') ;
hold off;



%% plotting primal sub-optimality and duality gap

figure(2)
fig_distance = tiledlayout(1,1,'Padding','tight');
fig_distance.Units = 'inches';
fig_distance.OuterPosition = [0.25 0.25 5 5];
nexttile;

semilogy(k,primal_sub_optimality,'-o','MarkerSize',3);
hold on
semilogy(k,duality_gap,'-x','MarkerSize',3);
grid on;
axis padded;
xlabel('Iterations, $k$', 'FontSize',10,'Interpreter','latex');
legend('Primal sub optimality', 'Duality gap','FontSize',15,'Interpreter','latex','Location','northeast');
im = gcf;
exportgraphics(im,'sub_opt_rhoalpha0.pdf','ContentType','vector') ;
hold off;