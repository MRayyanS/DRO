This repository contains the MATLAB code for the DRO problem of portfolio selection to minimize the Entropic risk.

$$
\min_{x \in \Delta^n} \sup_{\mathbb{P} \in \mathcal{P}} \quad \mathcal{E}({x},{\mathbb{P}}) \coloneqq \Sigma_{j = 1}^{n} {\frac{1}{\theta_j} \log \left( \mathbb{E}_{\mathbb{P}_j} [e^{-\theta_j x_j \xi_j}] \right) }  .
$$

For more information about the problem, see [1], particularly Section 6 therein.

1. The "portfolio_optimization_DRO.m" is the main file

2. The "generate_data.m" file generates synthetic data for the problem as discussed in [1, Section 6.3]

3. The "entropic_risk_DRO.m" file takes the prolem data and solves the DRO problem

4. The "solve_for_x.m" file implements $` \min_{x \in  \Delta^n} \mathcal{E}({x},{\mathbb{P}}) `$ using FISTA algorith, [2]

5. The "FW_oracle_entropic_risk.m" is the FW oracle for the Entropic risk DRO problem [Lemma 6.4, 1]

6. The "compute_worst_case_cost.m" file evaluates $` \sup_{\mathbb{P} \in \mathcal{P}} \quad \mathcal{E}({x},{\mathbb{P}}) `$ for plotting purpose

7. The "plot_n_pdf.m" file generates the convergence plots for the algorithm and creates individual pdfs


%%%==========================================================

References

[1] NDRO paper

[2] FISTA paper
