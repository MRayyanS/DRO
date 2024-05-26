# NDRO - minimum variance portfolio selection

This repository contains the MATLAB code for the DRO problem of portfolio selection to minimize the variance.

$$
\min_{x \in \Delta^n} \sup_{\mathbb{P} \in \mathcal{P}} \quad V ({x},{\mathbb{P}}) \coloneqq x^\top \big( \Sigma_{\mathbb{P}} - \mu_{\mathbb{P}} \mu_{\mathbb{P}}^\top \big) x ,
$$

where $` \Sigma_{\mathbb{P}} \coloneqq \mathbb{E}_{\mathbb{P}} [\xi \xi^\top] `$ and $` \mu_{\mathbb{P}} \coloneqq \mathbb{E}_{\mathbb{P}} [\xi] `$. For more information about the problem, see [Section 7, 1], particularly page no. 31.

1. The "portfolio_optimization_DRO_ellipsoidal_support.m" is the main file

2. The "generate_data_ellipsoidal_support.m" file generates synthetic data for the problem as discussed in [1, Section 7.3]

3. The "min_variance_DRO_ellipsoidal_support.m" file takes the prolem data and solves the DRO problem

4. The "solve_for_x.m" file implements $` \min_{x \in  \Delta^n} V({x},{\mathbb{P}}) `$ using cvx

5. The "FWOracle_min_variance_ellipsoidal_suuport.m" is the FW oracle for the minimum variance potfolio selection problem [Lemma 7.6, 1]

6. The "compute_worst_case_cost.m" file evaluates $` \sup_{\mathbb{P} \in \mathcal{P}} V({x},{\mathbb{P}}) `$ for plotting purpose

7. The "plot_n_pdf.m" file generates the convergence plots for the algorithm and creates individual pdfs


%%%==========================================================

References

[1] M. R. Sheriff, P. Mohajerin Esfahani, "Nonlinear Distributionally Robust Optimization", arXiv: 2306.03202
