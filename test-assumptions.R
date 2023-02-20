library(readr)
library(dplyr)
library(purrr)
library(ellipticalsymmetry)
library(robustbase)

#' Compute PE Test Statistic
#'
#' @param r Double vector. Sample from a distribution.
#' @param eta Double. Free parameter, see (Hüsler and Li 2006) for details.
#' @param m Integer. Discretization level in integration.
#' @param k Integer. Number of tail observations.
#'
#' @return Double. Approximate test statistic.
compute_test_stat_regvar <- function(r, eta = 0.5, m = 10000, k = 80) {
  n <- length(r)
  t <- seq(0, 1, length.out = m + 1)[-1]
  r <- sort(r, decreasing = FALSE)
  d <- log(r[n - floor(k * t)]) - log(r[n - k])
  gamma <- mean((log(r[(n - k):n]) - log(r[n - k]))[-1])
  k * mean((d / gamma + log(t))^2 * t^eta)
}

#' Simulate Observation from Asymptotic Distribution of the PE Test Statistic
#'
#' @param eta Double. Free parameter, see (Hüsler and Li 2006) for details.
#' @param m Integer. Discretization level in integration.
#'
#' @return Double. Approximate value for an observation from asymptotic
#'  distribution 
sim_asymp_test_stat_regvar <- function(eta = 0.5, m = 10000) {
  t <- seq(0, 1, length.out = m + 1)[-1]
  bm <- cumsum(rnorm(m, 0, sqrt(1 / m)))
  bm_bridge <- bm - t * bm[m]
  integral_1 <- mean((1 / t) * bm_bridge)
  mean(((1 / t) * bm_bridge + log(t) * integral_1)^2 * t^eta)
}

#' Computer P-Value for the PE test
#'
#' @param r Double vector. Sample from a distribution.
#' @param eta Double. Free parameter, see (Hüsler and Li 2006) for details.
#' @param m Integer. Discretization level in integration.
#' @param k Integer. Number of tail observations.
#' @param mm Integer. Sample size from asymptotic distribution.
#'
#' @return Double. Approximate p-value for the PE-test.
test_regvar <- function(r, eta = 0.5, m = 10000, k = 80, mm  = 10000) {
  test_stat <- compute_test_stat_regvar(r, eta, m, k)
  fn <- ecdf(replicate(mm, sim_asymp_test_stat_regvar(eta, m)))
  1 - fn(test_stat)
}

stock <- read_csv("data/stock-price-mod.csv", col_names = TRUE)

# Ljung-box test for returns
box_test_return <- stock %>%
  select(ends_with("return")) %>%
  map(~ Box.test(., lag = 20, type = "Ljung-Box"))

box_test_return$us_return$p.value
box_test_return$uk_return$p.value
box_test_return$jpn_return$p.value

# Ljung-box test for innovations
box_test_innovation <- stock %>%
  select(contains("innovation")) %>%
  map(~ Box.test(., lag = 20, type = "Ljung-Box"))

box_test_innovation$us_return$p.value
box_test_innovation$uk_return$p.value
box_test_innovation$jpn_return$p.value

# Test for ellipticity (Huffer and Park 2007)
ellipticity_test <- stock %>%
  select(contains("innovation")) %>%
  as.matrix %>%
  HufferPark(c = 3, R = 1000, sector = "orthants", nJobs = -1)

ellipticity_test$p.value

# Visual inspection of the ellipticity assumption
plot(stock$us_return_innovation, stock$uk_return_innovation)
plot(stock$us_return_innovation, stock$jpn_return_innovation)
plot(stock$uk_return_innovation, stock$jpn_return_innovation)

# Test regular variation of the Mahalanobis distance
sigma <- stock %>%
  select(contains("innovation")) %>%
  covMcd(alpha = 0.5)

reg_test <- stock %>%
  select(contains("innovation")) %>%
  as.matrix %>%
  stats::mahalanobis(center = sigma$center, cov = sigma$cov) %>%
  test_regvar(eta = 0.5, m = 10000, k = 80, mm = 10000)
