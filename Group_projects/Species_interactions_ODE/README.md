# Estimating species interactions from dynamic data - Florian Hartig and Frank Pennekamp

This repository contains data and information for the Bayesian thinking project on estimating species interactions from dynamic data

## Project goals

The project aims to infer species interactions from dynamic data using ordinary differential equations fitted to time-series. Refer to the Species_interactions_ODE_intro.pdf for further details. More specific goals will be defined as soon as the group for part B has formed. 

## Data
There are currently three datasets from the following publications:

Daugaard, U., Petchey, O. L., & Pennekamp, F. (2019). Warming can destabilize predator–prey interactions by shifting the functional response from Type III to Type II. Journal of Animal Ecology, 88(10), 1575–1586. https://doi.org/10.1111/1365-2656.13053

Uszko, W., Diehl, S., Englund, G., & Amarasekare, P. (2017). Effects of warming on predator–prey interactions – a resource-based approach and a theoretical synthesis. Ecology Letters, 20(4), 513–523. https://doi.org/10.1111/ele.12755

DeLong, J. P., & Lyon, S. (2020). Temperature alters the shape of predator–prey cycles through effects on underlying mechanisms. PeerJ, 8, e9377. https://doi.org/10.7717/peerj.9377

The first two datasets are from classic functional response experiments, whereas the third provides time-series. All three have temperature treatments included.

## Previous work

### Predator-prey models implemented in Stan
https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html

https://www.magesblog.com/post/2021-02-08-fitting-multivariate-ode-models-with-brms/

## Relevant literature

Adams, M. P., Sisson, S. A., Helmstedt, K. J., Baker, C. M., Holden, M. H., Plein, M., Holloway, J., Mengersen, K. L., & McDonald‐Madden, E. (n.d.). Informing management decisions for ecological networks, using dynamic models calibrated to noisy time-series data. Ecology Letters, n/a(n/a). https://doi.org/10.1111/ele.13465

Boersch‐Supan, P. H., Ryan, S. J., & Johnson, L. R. (2017). deBInfer: Bayesian inference for dynamical models of biological systems in R. Methods in Ecology and Evolution, 8(4), 511–518. https://doi.org/10.1111/2041-210X.12679

Pascual, M. A., & Kareiva, P. (1996). Predicting the Outcome of Competition Using Experimental Data: Maximum Likelihood and Bayesian Approaches. Ecology, 77(2), 337–349. https://doi.org/10.2307/2265613

Rosenbaum, B., Raatz, M., Weithoff, G., Fussmann, G. F., & Gaedke, U. (2019). Estimating Parameters From Multiple Time Series of Population Dynamics Using Bayesian Inference. Frontiers in Ecology and Evolution, 6. https://doi.org/10.3389/fevo.2018.00234

Rosenbaum, B., & Rall, B. C. (2018). Fitting functional responses: Direct parameter estimation by simulating differential equations. Methods in Ecology and Evolution, 9(10), 2076–2090. https://doi.org/10.1111/2041-210X.13039



