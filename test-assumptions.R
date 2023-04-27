suppressPackageStartupMessages(library(dplyr))


#' Compute PE test statistic
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


#' Simulate observation from asymptotic distribution of the PE test statistic
#'
#' @param eta Double. Free parameter, see (Hüsler and Li 2006) for details.
#' @param m Integer. Discretization level in integration.
#'
#' @return Double. Approximate value for an observation from asymptotic
#'  distribution.
sim_asymp_test_stat_regvar <- function(eta = 0.5, m = 10000) {
  t <- seq(0, 1, length.out = m + 1)[-1]
  bm <- cumsum(rnorm(m, 0, sqrt(1 / m)))
  bm_bridge <- bm - t * bm[m]
  integral_1 <- mean((1 / t) * bm_bridge)
  mean(((1 / t) * bm_bridge + log(t) * integral_1)^2 * t^eta)
}


#' Computer p-value for the PE test
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


# Read data
stock <- readr::read_csv("data/stock-price-mod.csv", col_names = TRUE,
                  col_types = readr::cols())

# Ljung-Box test for returns
box_test_return <- stock %>%
  select(ends_with("return")) %>%
  purrr::map(~ Box.test(., lag = 20, type = "Ljung-Box")) %>%
  purrr::map(~ broom::glance(.)$p.value) %>%
  as_tibble() %>%
  rename(us = us_return, uk = uk_return, jpn = jpn_return)

# Ljung-Box test for innovations
box_test_innovation <- stock %>%
  select(ends_with("innovation")) %>%
  purrr::map(~ Box.test(., lag = 20, type = "Ljung-Box")) %>%
  purrr::map(~ broom::glance(.)$p.value) %>%
  as_tibble() %>%
  rename(us = us_innovation, uk = uk_innovation, jpn = jpn_innovation)

# Combine results for Ljung-Box test
box_test <- bind_rows(box_test_return, box_test_innovation) %>%
  mutate(type = c("return", "innovation"))

# Test for ellipticity (Huffer and Park 2007)
ellipticity_test <- stock %>%
  select(contains("innovation")) %>%
  as.matrix() %>%
  ellipticalsymmetry::HufferPark(c = 3, R = 1000, sector = "orthants",
                                 nJobs = -1) %>%
  broom::glance() %>%
  pull(p.value)

# Test regular variation of the estimated generating variate
sigma <- stock %>%
  select(ends_with("innovation")) %>%
  robustbase::covMcd(alpha = 0.5)

reg_test <- stock %>%
  select(ends_with("innovation")) %>%
  as.matrix %>%
  stats::mahalanobis(center = sigma$center, cov = sigma$cov) %>%
  sqrt() %>%
  test_regvar(eta = 0.5, m = 10000, k = 80, mm = 10000)

# Combine test results for ellipticity and regular variation
test_res <- tibble(ellipticity_test = ellipticity_test,
                   regular_variation_test = reg_test)

cli::cli_h3("P-values of the the Ljung-Box test")
box_test

cli::cli_h3("P-values for the tests of ellipticity and regular variation")
test_res

# Save p-values
readr::write_csv(box_test, "data/box_test.csv")
readr::write_csv(test_res, "data/test_res.csv")
