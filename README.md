# structmcmc

structmcmc is a set of tools for performing structural inference for Bayesian Networks using MCMC in [`R`][R].
The widely-used MC<sup>3</sup> ([Madigan & Raftery, 1995][Madigan:1995p10499]) is implemented, as well as a number of variants of the algorithm.
The MC<sup>3</sup> algorithm is a Metropolis-Hastings sampler for which the target distribution is the posterior distribution of Bayesian networks.
Tools for exact solutions are also available, but for networks with more than, say, 6 nodes, these will be prohibitively slow.

The implementation allows the local conditional distributions to be multinomial or Gaussian, using standard priors.
Arbitrary structural priors for the Bayesian network can be specified.
The main difficulty in sampling Bayesian networks efficiently is ensuring the acyclicity constraint is not violated.
The package implements the cycle-checking methods introduced by [King & Sagert (2002)][King:2002gt], which is an alternative to the method introduced by [Giudici & Castelo (2003)][Giudici:2003cn].
To enable convergence to be assessed, a number of tools for creating diagnostic plots are included.

Interfaces to a number of other `R` packages for Bayesian networks are available, including [`deal`][cran:deal] (hill-climbing and heuristic search), [`bnlearn`][cran:bnlearn] (a number of constraint-based and score-based algorithms) and [`pcalg`][cran:pcalg] (PC-algorithm).
An interface to [`gRain`][cran:gRain] is also included to allow its probability propagation routines to be used easily.

# Discrete data
Each random variable has a Multinomial distribution, with the conjugate Dirichlet priors.

``` r
# Multinomial-Dirichlet model
x1 <- factor(c(1, 1, 0, 0, 1, 1, 0, 1, 0))
x2 <- factor(c(1, 0, 1, 1, 0, 1, 1, 0, 0))
x3 <- factor(c(1, 0, 0, 0, 0, 0, 0, 1, 1))
x <- data.frame(x1 = x1, x2 = x2, x3 = x3)

# MC^3 algorithm
set.seed(1234)
mcmc <- posterior(data = x, method = "mh-mcmc")
epmcmc <- ep(mcmc)
levelplot(epmcmc)

# Exact computation by exhaustive enumeration
exact <- posterior(data = x, method = "exact")
epexact <- ep(exact)
levelplot(epexact)

mcmc2 <- posterior(data = x, method = "mh-mcmc")
epmcmc2 <- ep(mcmc2)
levelplot(epmcmc2)

levelplot(ep.list(exact = epexact, mcmc = epmcmc))

cumep(list(mcmc, mcmc2))
```

# Continuous data
Each random variable is assumed to be Normally-distributed, with a G-prior...

...

# Installation

# 


[R]: http://www.r-project.org "The R Project for Statistical Computing"
[Madigan:1995p10499]: http://www.jstor.org/stable/1403615  "Madigan, D., & York, J. C. (1995). Bayesian Graphical Models for Discrete Data. International Statistical Review / Revue Internationale de Statistique, 63(2), 215-232."
[King:2002gt]: http://dx.doi.org/10.1006/jcss.2002.1883 "King, V., & Sagert, G. (2002). A Fully Dynamic Algorithm for Maintaining the Transitive Closure. Journal of Computer and System Sciences, 65(1), 150-167."
[Giudici:2003cn]: http://dx.doi.org/10.1023/A:1020202028934 "Giudici, P., & Castelo, R. (2003). Improving Markov Chain Monte Carlo Model Search for Data Mining. Machine Learning, 50, 127-158."
[cran:deal]: http://cran.r-project.org/web/packages/deal/ "deal: Learning Bayesian Networks with Mixed Variables"
[cran:bnlearn]: http://cran.r-project.org/web/packages/bnlearn/ "bnlearn: Bayesian network structure learning, parameter learning and inference"
[cran:pcalg]: http://cran.r-project.org/web/packages/pcalg/ "pcalg: Estimation of CPDAG/PAG and causal inference using the IDA algorithm"
[cran:gRain]: http://cran.r-project.org/web/packages/gRain "gRain: Graphical Independence Networks"