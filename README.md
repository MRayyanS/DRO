This repository contains the MATLAB code for the DRO problem of portfolio selection to minimize the Entropic risk.
$$
\min_{x \in \Delta^n} \sup_{\mathbb{P} \in \mathcal{P}} \quad \mathcal{E}({x},{\mathbb{P}}) \coloneqq \Sigma_{j = 1}^{n} {\frac{1}{\theta_j} \log \left( \mathbb{E}_{\mathbb{P}_j} [e^{-\theta_j x_j \xi_j}] \right) }  .
$$
For more information about the problem, see [1], particularly Section 6 therein.
