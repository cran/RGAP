# RGAP: Output Gap Estimation in R

'RGAP' provides tools for modeling and estimating the bivariate unobserved component models involved in the European Commission's Cobb-Douglas production function methodology to estimate potential output and the output gap.

[KOF Working Paper](http://hdl.handle.net/20.500.11850/552089)

If you use 'RGAP' in your paper, please cite it properly, see `citation("RGAP")` in R, or above link to the paper.

## Details

The output gap indicates the percentage difference between the actual output of an economy and its potential. Since potential output is a latent process, the estimation of the output gap poses a challenge and numerous filtering techniques have been proposed. 'RGAP' facilitates the estimation of a Cobb-Douglas production function type output gap, as suggested by the European Commission ([Havik et al. 2014](https://ideas.repec.org/p/euf/ecopap/0535.html)). To that end, the non-accelerating wage rate of unemployment (NAWRU) and the trend of total factor productivity (TFP) can be estimated in two bivariate unobserved component (UC) models by means of Kalman filtering and smoothing. 'RGAP' features a flexible modeling framework for the appropriate state-space models and offers frequentist as well as Bayesian estimation techniques. Additional functionalities include direct access to the 'AMECO' database and automated model selection procedures.

## Main features

- Data fetching from 'AMECO' database
- Data pre processing
- Modeling bivariate state-space models for the NARWU and the TFP trend
- Estimation of defined models via the Kalman filter and smoother and
  - maximum likelihood estimation or
  - bayesian estimation via Gibbs procedure
- Output gap computation
- Prediction
- Alternative approaches: HP-filter and bivariate UC model by Kuttner (1994)
- Additional features: automated model selection

## Install the package
You can install the package from ‘Github’ using the **install_github** function from the **devtools** package.
``` 
library(devtools)
install_github('sinast3000/RGAP')
```

***

Kuttner, K. N. (1994), Estimating potential output as a latent variable, Journal of Business & Economic Statistics 12(3), 361–368.

Havik, K., Mc Morrow, K., Orlandi, F., Planas, C., Raciborski, R., Roeger, W., Rossi, A., Thum-Thysen, A. & Vandermeulen, V. (2014), The production function methodology for calculating potential growth rates & output gaps, Technical report, Directorate
General Economic and Financial Affairs (DG ECFIN), European Commission.


