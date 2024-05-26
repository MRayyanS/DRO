# NDRO - minimum variance portfolio selection

This repository contains the MATLAB code for the DRO problem of portfolio selection to minimize the Entropic risk.

$$
\min_{x \in \Delta^n} \sup_{\mathbb{P} \in \mathcal{P}} \quad V ({x},{\mathbb{P}}) \coloneqq x^\top \big( \Sigma_{\mathb{P}} - \mu_{\mathb{P}} {\mathb{P}}^\top \big) x .
$$

For more information about the problem, see [1], particularly Section 6 therein.

1. The "portfolio_optimization_DRO.m" is the main file

2. The "generate_data.m" file generates synthetic data for the problem as discussed in [1, Section 6.3]

3. The "entropic_risk_DRO.m" file takes the prolem data and solves the DRO problem

4. The "solve_for_x.m" file implements $` \min_{x \in  \Delta^n} \mathcal{E}({x},{\mathbb{P}}) `$ using FISTA algorith, [2]

5. The "FW_oracle_entropic_risk.m" is the FW oracle for the Entropic risk DRO problem [Lemma 6.4, 1]

6. The "compute_worst_case_cost.m" file evaluates $` \sup_{\mathbb{P} \in \mathcal{P}} \mathcal{E}({x},{\mathbb{P}}) `$ for plotting purpose

7. The "plot_n_pdf.m" file generates the convergence plots for the algorithm and creates individual pdfs


%%%==========================================================

References

[1] M. R. Sheriff, P. Mohajerin Esfahani, "Nonlinear Distributionally Robust Optimization", arXiv: 2306.03202

[2] A. Beck, and M. Teboulle. "A fast iterative shrinkage-thresholding algorithm for linear inverse problems." SIAM journal on imaging sciences 2.1 (2009): 183-202.
