library(readr)
library(dplyr)
library(purrr)
library(ellipticalsymmetry)

compute_test_stat_radial <- function(r, eta = 0.5, m = 10000, k = 80) {
  n <- length(r)
  t <- seq(0, 1, length.out = m + 1)[-1]
  r <- sort(r, decreasing = FALSE)
  d <- log(r[n - floor(k * t)]) - log(r[n - k])
  gamma <- mean((log(r[(n - k):n]) - log(r[n - k]))[-1])
  k * mean((d / gamma + log(t))^2 * t^eta)
}

sim_asymp_test_stat_radial <- function(eta = 0.5, m = 10000) {
  t <- seq(0, 1, length.out = m + 1)[-1]
  bm <- cumsum(rnorm(m, 0, sqrt(1 / m)))
  bm_bridge <- bm - t * bm
  integral_1 <- mean((1 / t) * bm_bridge)
  mean(((1 / t) * bm_bridge + log(t) * integral_1)^2 * t^eta)
}

test_radial <- function(r, eta = 0.5, m = 10000, k = 80, mm  = 10000) {
  test_stat <- compute_test_stat_radial(r, eta, m, k)
  fn <- ecdf(replicate(mm, sim_asymp_test_stat_radial(eta, m)))
  1 - fn(test_stat)
}

stock <- read_csv("data/stock-price-mod.csv", col_names = TRUE)

# Ljung-box test for innovations
box_test <- stock %>%
  select(contains("innovation")) %>%
  map(~ Box.test(., lag = 20, type = "Ljung-Box"))

# Test for ellipticity (Huffer and Park 2007)
ellipticity_test <- stock %>%
  select(contains("innovation")) %>%
  as.matrix %>%
  HufferPark(c = 3, R = 1000, sector = "orthants", nJobs = -1)

# Visual inspection of the ellipticity assumption
plot(stock$us_return_innovation, stock$uk_return_innovation)
plot(stock$us_return_innovation, stock$jpn_return_innovation)
plot(stock$uk_return_innovation, stock$jpn_return_innovation)
