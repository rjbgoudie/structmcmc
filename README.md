structmcmc
==========

structmcmc is a set of tools for performing structural inference for Bayesian Networks using MCMC in [`R`][R], a free software environment for statistical computing and graphics.

The widely-used MC<sup>3</sup> algorithm ([Madigan & Raftery, 1995][Madigan:1995p10499]) is implemented, as well as a number of variants of the algorithm. The MC<sup>3</sup> algorithm is a Metropolis-Hastings sampler for which the target distribution is the posterior distribution of Bayesian networks. Tools for exact solutions are also available, but for networks with more than, say, 6 nodes, these will be prohibitively slow.

The implementation allows the local conditional distributions to be multinomial or Gaussian, using standard priors. Arbitrary structural priors for the Bayesian network can be specified. The main difficulty in sampling Bayesian networks efficiently is ensuring the acyclicity constraint is not violated. The package implements the cycle-checking methods introduced by [King & Sagert (2002)][King:2002gt], which is an alternative to the method introduced by [Giudici & Castelo (2003)][Giudici:2003cn]. To enable convergence to be assessed, a number of tools for creating diagnostic plots are included.

Interfaces to a number of other `R` packages for Bayesian networks are available, including [`deal`][cran:deal] (hill-climbing and heuristic search), [`bnlearn`][cran:bnlearn] (a number of constraint-based and score-based algorithms) and [`pcalg`][cran:pcalg] (PC-algorithm). An interface to [`gRain`][cran:gRain] is also included to allow its probability propagation routines to be used easily.

Basic operation, discrete data
------------------------------
[View script as file](https://gist.github.com/970279)

Each random variable has a Multinomial distribution, with the conjugate Dirichlet priors.

Data must be supplied as a [`data.frame`][rdoc:data.frame] with `p` columns (corresponding to `p` random variables) and `n` columns (corresponding to the `n` samples). Each column must be a [`factor`][rdoc:factor] variable.

``` r
x1 <- factor(c("a", "a", "g", "c", "c", "a", "g", "a", "a"))
x2 <- factor(c(2, 2, 4, 3, 1, 4, 4, 4, 1))
x3 <- factor(c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
x <- data.frame(x1 = x1, x2 = x2, x3 = x3)
```

Draw samples from the posterior using MC<sup>3</sup>.

``` r
set.seed(1234)
initial <- bn(c(), c(), c())
mcmc <- posterior(data = x, method = "mc3", nSamples = 10000, nBurnin = 1000, initial = initial)
```

Compute and plot estimated edge probabilities given by the MCMC run

``` r
epmcmc <- ep(mcmc)
levelplot(epmcmc)
```

Since this is a problem with `p = 3`, we can compute the posterior edge probabilies by exhaustive enumeration. This is only feasible for `p <= 6` or so.

``` r
exact <- posterior(x, "exact")
epexact <- ep(exact)
levelplot(epexact)
```

Comparing multiple MCMC runs

``` r
mcmc2 <- posterior(data = x, method = "mc3", nSamples = 10000, nBurnin = 1000, initial = initial)
epmcmc2 <- ep(mcmc2)
levelplot(epmcmc2)
```

Compare the final edge probabilities between runs

``` r
splom(bnpostmcmc.list(mcmc, mcmc2))
levelplot(ep.list(exact = epexact, mcmc = epmcmc))
```

Plot how the cumulative edge probabilities change as samples are drawn.

``` r
xyplot(cumep(bnpostmcmc.list(mcmc, mcmc2)))
```

Plot how the moving averaging edge probabilities change as samples are drawn.

``` r
xyplot(mwep(bnpostmcmc.list(mcmc, mcmc2)))
```

Plot how the cumulative total variance distance changes as samples are drawn

``` r
exactgp <- gp(exact)
xyplot(cumtvd(exactgp, bnpostmcmc.list(mcmc, mcmc2)))
```


Basic operation, continuous data
--------------------------------
[View script as file](https://gist.github.com/974390)

Each random variable is assumed to be Normally-distributed, with a G-prior.

Data must be supplied as a [`matrix`][rdoc:matrix] with `p` columns (corresponding to `p` random variables) and `n` columns (corresponding to the `n` samples). Each column must be a [`numeric`][rdoc:numeric] variable.

``` r
x1 <- rnorm(20)
x2 <- rnorm(20)
x3 <- rnorm(20)
x <- matrix(c(x1, x2, x3), ncol = 3)
```

Draw samples from the posterior using MC<sup>3</sup>.

``` r
set.seed(1234)
initial <- bn(c(), c(), c())
mcmc <- posterior(data = x, method = "mc3",
                  logScoreFUN = logScoreZellnerFUN(),
                  nSamples = 10000, nBurnin = 1000, initial = initial)
```

Compute and plot estimated edge probabilities given by the MCMC run

``` r
epmcmc <- ep(mcmc)
levelplot(epmcmc)
```

Since this is a problem with `p = 3`, we can compute the posterior edge probabilies by exhaustive enumeration. This is only feasible for `p <= 6` or so.

``` r
exact <- posterior(x, "exact", logScoreFUN = logScoreZellnerFUN())
epexact <- ep(exact)
levelplot(epexact)
```

Comparing multiple MCMC runs

``` r
mcmc2 <- posterior(data = x, method = "mc3",
                   logScoreFUN = logScoreZellnerFUN(),
                   nSamples = 10000, nBurnin = 1000, initial = initial)
epmcmc2 <- ep(mcmc2)
levelplot(epmcmc2)
```

Compare the final edge probabilities between runs

``` r
splom(bnpostmcmc.list(mcmc, mcmc2))
levelplot(ep.list(exact = epexact, mcmc = epmcmc))
```

Plot how the cumulative edge probabilities change as samples are drawn.

``` r
xyplot(cumep(bnpostmcmc.list(mcmc, mcmc2)))
```

Plot how the moving averaging edge probabilities change as samples are drawn.

``` r
xyplot(mwep(bnpostmcmc.list(mcmc, mcmc2)))
```

Plot how the cumulative total variance distance changes as samples are drawn

``` r
exactgp <- gp(exact)
xyplot(cumtvd(exactgp, bnpostmcmc.list(mcmc, mcmc2)))
```

Installation
------------
Download the current version, and `unzip` the file. Then install in `R` using the following, where `rjbgoudie-structmcmc-XXXXX` is the name of the `unzip`ped directory/folder, and `path/to/rjbgoudie-structmcmc-XXXXX` is the path to this folder.

``` r
install.packages("path/to/rjbgoudie-structmcmc-XXXXX", repos = NULL, type = "source")
```

The package depends on [`parental`][pkg:parental], which provides a very lightweight directed graph object and basic manipulation tools for R.

The package also depends on [`lattice`][cran:lattice], and `grid`, both of which are included with R.

[`gRain`][cran:gRain], [`nnet`][cran:nnet], [`reshape`][cran:reshape], [`zoo`][cran:zoo] are also recommended, and can be installed from CRAN, e.g. using

``` r
install.packages("reshape")
```

Contact
-------


[R]: http://www.r-project.org "The R Project for Statistical Computing"
[Madigan:1995p10499]: http://www.jstor.org/stable/1403615  "Madigan, D., & York, J. C. (1995). Bayesian Graphical Models for Discrete Data. International Statistical Review / Revue Internationale de Statistique, 63(2), 215-232."
[King:2002gt]: http://dx.doi.org/10.1006/jcss.2002.1883 "King, V., & Sagert, G. (2002). A Fully Dynamic Algorithm for Maintaining the Transitive Closure. Journal of Computer and System Sciences, 65(1), 150-167."
[Giudici:2003cn]: http://dx.doi.org/10.1023/A:1020202028934 "Giudici, P., & Castelo, R. (2003). Improving Markov Chain Monte Carlo Model Search for Data Mining. Machine Learning, 50, 127-158."
[cran:deal]: http://cran.r-project.org/web/packages/deal/ "deal: Learning Bayesian Networks with Mixed Variables"
[cran:bnlearn]: http://cran.r-project.org/web/packages/bnlearn/ "bnlearn: Bayesian network structure learning, parameter learning and inference"
[cran:pcalg]: http://cran.r-project.org/web/packages/pcalg/ "pcalg: Estimation of CPDAG/PAG and causal inference using the IDA algorithm"
[cran:gRain]: http://cran.r-project.org/web/packages/gRain "gRain: Graphical Independence Networks"
[cran:lattice]: http://cran.r-project.org/web/packages/lattice "lattice: Lattice Graphics"
[cran:nnet]: http://cran.r-project.org/web/packages/nnet "nnet: Feed-forward Neural Networks and Multinomial Log-Linear Models"
[cran:reshape]: http://cran.r-project.org/web/packages/reshape "reshape: Flexibly reshape data"
[cran:zoo]: http://cran.r-project.org/web/packages/zoo "zoo: Z's ordered observations"
[rdoc:factor]: http://stat.ethz.ch/R-manual/R-devel/library/base/html/factor.html "R Documentation: Factors"
[rdoc:data.frame]: http://stat.ethz.ch/R-manual/R-devel/library/base/html/data.frame.html "R Documentation: Data Frames"
[rdoc:matrix]: http://stat.ethz.ch/R-manual/R-devel/library/base/html/matrix.html "R Documentation: Matrices"
[rdoc:numeric]: http://stat.ethz.ch/R-manual/R-devel/library/base/html/numeric.html "R Documentation: Numeric"
[pkg:parental]: https://github.com/rjbgoudie/parental "parental: a very lightweight directed graph object and basic manipulation tools for R"
