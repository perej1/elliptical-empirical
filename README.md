# elliptical-empirical

Empirical example about estimation of multivariate elliptical extreme quantiles in a financial context. Data consists of three time series of market price indices from 2.7.2001 to 29.6.2007:

- Standard and Poors S&P 500,
- Financial Times Stock Exchange FTSE 100 and
- Nikkei 225.

We compute returns $Y_t = \log(X_{t+1} / X_t)$ and fit $\mathrm{EGARCH}(1, 1)$ model to time series of returns for each index. Then we use innovations for computing extreme quantile regions for several probabilities $p\in(1/2000, 1/5000, 1/10000)$.

## Requirements

1. Clone or unzip the repository.
    ```
    git clone git@github.com:perej1/elliptical-sim.git
    ```

2. Install required packages by running the following R command in the project's root folder (R package `renv` has to be installed).
    ```
    renv::restore()
    ```

## Running the code

Emprical example involves three steps:

1. `modify.R` - Fitting of $\mathrm{EGARCH}(1, 1)$ models and computation of innovations.

2. `test-assumptions.R` - Test assumptions of independence, ellipticity and regular variation.

3. `analyze.R` - Estimation of extreme quantile regions.

Lastly, the script `main.R` is responsible for argument parsing. In order to perform all the three parts of the empirical example, just run the following.
```
Rscript main.R
```

If you wish to skip some parts of the empirical example, you can run something like below.
```
Rscript main.R --modify TRUE --test FALSE --analyze TRUE
```
That is, above snippet of code skips statistical tests.