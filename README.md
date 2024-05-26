This repository contains the MATLAB code for the DRO problem of portfolio selection to minimize the Entropic risk. See [1], particularly Section 6 therein for a more information about the problem.

$$
\min_{x \in \Delta^n} \sup_{\mathbb{P} \in \mathcal{P}} \quad \mathcal{E}({x},{\mathbb{P}}) \coloneqq \Sigma\limits_{j = 1}^{n} {\frac{1}{\theta_j} \log \left( \mathbb{E}_{\mathbb{P}_j} [e^{-\theta_j x_j \xi_j}] \right) }  .
$$

