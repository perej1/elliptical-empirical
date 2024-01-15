# Test assumptions of independence, ellipticity and regular variation
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)


#' Compute value of the Hill estimator
#'
#' @param data Double vector of data.
#' @param k Integer, number of tail observations.
#'
#' @return Double, estimate for the extreme value index
hill <- function(data, k) {
  n <- length(data)
  data <- sort(data, decreasing = FALSE)
  mean((log(data[(n - k):n]) - log(data[n - k]))[-1])
}


#' Compute PE test statistic
#'
#' @param r Double vector representing a sample from a distribution.
#' @param eta Double, free parameter, see (Hüsler and Li 2006) for details.
#' @param m Integer, discretization level in integration.
#' @param k Integer, number of tail observations.
#'
#' @return Double representing approximate test statistic.
compute_test_stat_regvar <- function(r, eta = 0.5, m = 10000, k = 160) {
  n <- length(r)
  t <- seq(0, 1, length.out = m + 1)[-1]
  r <- sort(r, decreasing = FALSE)
  d <- log(r[n - floor(k * t)]) - log(r[n - k])
  gamma <- mean((log(r[(n - k):n]) - log(r[n - k]))[-1])
  k * mean((d / gamma + log(t))^2 * t^eta)
}


#' Simulate observation from asymptotic distribution of the PE test statistic
#'
#' @param eta Double, free parameter, see (Hüsler and Li 2006) for details.
#' @param m Integer, discretization level in integration.
#'
#' @return Double, approximate value for an observation from asymptotic
#'  distribution.
sim_asymp_test_stat_regvar <- function(eta = 0.5, m = 10000) {
  t <- seq(0, 1, length.out = m + 1)[-1]
  bm <- cumsum(rnorm(m, 0, sqrt(1 / m)))
  bm_bridge <- bm - t * bm[m]
  integral_1 <- mean((1 / t) * bm_bridge)
  mean(((1 / t) * bm_bridge + log(t) * integral_1)^2 * t^eta)
}


#' Computer approximate p-value for the PE test
#'
#' @param r Double vector representing a sample from a distribution.
#' @param eta Double, free parameter, see (Hüsler and Li 2006) for details.
#' @param m Integer, discretization level in integration.
#' @param k Integer, number of tail observations.
#' @param mm Integer, sample size from asymptotic distribution.
#'
#' @return Double, approximate p-value for the PE-test.
test_regvar <- function(r, eta = 0.5, m = 10000, k = 160, mm  = 10000) {
  test_stat <- compute_test_stat_regvar(r, eta, m, k)
  fn <- ecdf(replicate(mm, sim_asymp_test_stat_regvar(eta, m)))
  1 - fn(test_stat)
}


# Read data
stock <- readr::read_csv("data/stock-price-mod.csv", col_names = TRUE,
                         col_types = readr::cols())

innovations <- stock %>%
  select(ends_with("innovation")) %>%
  rename_with(~ gsub("_innovation", "", .))

# Ljung-Box test for returns
box_test_return <- stock %>%
  select(ends_with("return")) %>%
  purrr::map(~ Box.test(., lag = 20, type = "Ljung-Box")) %>%
  purrr::map(~ broom::glance(.)$p.value) %>%
  as_tibble() %>%
  rename_with(~ gsub("_return", "", .))

# Ljung-Box test for innovations
box_test_innovation <- innovations %>%
  purrr::map(~ Box.test(., lag = 20, type = "Ljung-Box")) %>%
  purrr::map(~ broom::glance(.)$p.value) %>%
  as_tibble()

# Combine results for Ljung-Box test
box_test <- bind_rows(box_test_return, box_test_innovation) %>%
  mutate(type = c("return", "innovation"))

# Histograms of marginals of innovations
for (country in c("us", "uk", "jpn")) {
  ggplot(innovations, aes(x = !!ggplot2::sym(country))) +
    geom_histogram(bins = 30)
  ggsave(stringr::str_c("figures/", country, "_hist", ".jpg"), width = 7,
         height = 7, dpi = 1000)
}

# QQ-plots and scatter plots
innovations_comb <- innovations %>%
  combn(2, simplify = FALSE)
for (data in innovations_comb) {
  countries <- names(data)
  names(data) <- c("x", "y")
  qq_data <- as_tibble(qqplot(data[[1]], data[[2]], plot.it = FALSE))
  ggplot(qq_data, aes(x = x, y = y)) +
    geom_point() +
    xlab(countries[1]) + ylab(countries[2])
  ggsave(stringr::str_c("figures/", stringr::str_c(countries, collapse = "_"),
                        "_qq", ".jpg"), width = 7, height = 7, dpi = 1000)
  
  ggplot(data, aes(x = x, y = y)) +
    geom_point() +
    xlab(countries[1]) + ylab(countries[2])
  ggsave(stringr::str_c("figures/", stringr::str_c(countries, collapse = "_"),
                        "_scatter", ".jpg"), width = 7, height = 7, dpi = 1000)
}

# Test for ellipticity (Huffer and Park 2007)
ellipticity_test <- innovations %>%
  as.matrix() %>%
  ellipticalsymmetry::HufferPark(c = 3, R = 1000, sector = "orthants",
                                 nJobs = -1) %>%
  broom::glance() %>%
  pull(p.value)

# Compute extreme value indices for the tails
gamma_left <- purrr::map(-innovations, ~ hill(., 80))
gamma_right <- purrr::map(innovations, ~ hill(., 80))

gamma <- bind_rows(gamma_left, gamma_right) %>%
  mutate(type = c("left", "right"))


# Test regular variation of the estimated generating variate
sigma <- robustbase::covMcd(innovations, alpha = 0.5)

reg_test <- innovations %>%
  as.matrix() %>%
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

cli::cli_h3("Extreme value indices")
gamma

# Save diagnostics
readr::write_csv(box_test, "data/box_test.csv")
readr::write_csv(test_res, "data/test_res.csv")
readr::write_csv(gamma, "data/gamma.csv")
