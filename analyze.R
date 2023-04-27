library(readr)
library(dplyr)
library(purrr)
library(ellipticalsymmetry)
library(robustbase)
library(purrr)
library(tidyr)
library(reshape2)
library(stringr)
library(ggplot2)


#' Generate m equally spaced points from a circle
#'
#' @param m Integer, number of points.
#'
#' @return Double matrix of points, one row represents one point.
get_ball_mesh <- function(m) {
  w <- seq(0, 2 * pi, length.out = m)
  cbind(cos(w), sin(w))
}


#' Compute square root of a positive definite matrix
#'
#' More precisely, function computes matrix lambda
#' s.t. sigma = lambda %*% lambda.
#'
#' @param sigma Double matrix, positive definite matrix.
#'
#' @return Double matrix, square root of the matrix.
sqrtmat <- function(sigma) {
  eigenval <- eigen(sigma)$values
  if (any(eigenval <= 0) || any(sigma != t(sigma))) {
    rlang::abort("`sigma` must be a symmetric positive definite matrix.")
  }
  eigenvec <- eigen(sigma)$vectors
  eigenvec %*% diag(eigenval^0.5) %*% t(eigenvec)
}


#' Estimate elliptical extreme quantile region
#'
#' Consistent for heavy-tailed elliptical distributions under some other
#' technical assumptions.
#'
#' @param data Double matrix of observations, each row represents one
#'   observation.
#' @param mu_est Double vector, estimate of the location.
#' @param sigma_est Double matrix, estimate of the scatter.
#' @param p Double, probability in quantile region.
#' @param k Integer, threshold for the sample from the tail.
#' @param m Integer, number of points to return.
#'
#' @return Double matrix, m points from the boundary of the quantile region.
elliptical_extreme_qregion <- function(data, mu_est, sigma_est, p, k, m) {
  n <- nrow(data)
  w <- get_ball_mesh(m)
  
  # Center data
  data <- sweep(data, 2, mu_est, "-")
  
  # Approximate generating variate
  radius <- sqrt(stats::mahalanobis(data, FALSE, sigma_est, inverted = FALSE))
  radius_sort <- sort(radius, decreasing = FALSE)
  
  # Estimate extreme value index
  gamma_est <- mean((log(radius_sort[(n - k):n]) - log(radius_sort[n - k]))[-1])
  
  # Estimate extreme quantile of generating variate
  r_hat <- radius_sort[n - k] * (k / (n * p))^gamma_est
  
  # Estimate extreme quantile region
  lambda <- sqrtmat(r_hat^2 * sigma_est)
  sweep(w %*% t(lambda), 2, mu_est, "+")
}


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


# Read data
stock <- read_csv("data/stock-price-mod.csv", col_names = TRUE)
parameters <- read_csv("data/model-parameters.csv", col_names = TRUE)

# Ljung-box test for returns
box_test_return <- stock %>%
  select(ends_with("return")) %>%
  map(~ Box.test(., lag = 20, type = "Ljung-Box"))

box_test_return$us_return$p.value
box_test_return$uk_return$p.value
box_test_return$jpn_return$p.value

# Ljung-box test for innovations
box_test_innovation <- stock %>%
  select(ends_with("innovation")) %>%
  map(~ Box.test(., lag = 20, type = "Ljung-Box"))

box_test_innovation$us_innovation$p.value
box_test_innovation$uk_innovation$p.value
box_test_innovation$jpn_innovation$p.value

# Test for ellipticity (Huffer and Park 2007)
ellipticity_test <- stock %>%
  select(contains("innovation")) %>%
  as.matrix %>%
  HufferPark(c = 3, R = 1000, sector = "orthants", nJobs = -1)

ellipticity_test$p.value

# Visual inspection of the ellipticity assumption
plot(stock$us_innovation, stock$uk_innovation)
plot(stock$us_innovation, stock$jpn_innovation)
plot(stock$uk_innovation, stock$jpn_innovation)

# Test regular variation of the Mahalanobis distance
sigma <- stock %>%
  select(ends_with("innovation")) %>%
  covMcd(alpha = 0.5)

reg_test <- stock %>%
  select(ends_with("innovation")) %>%
  as.matrix %>%
  stats::mahalanobis(center = sigma$center, cov = sigma$cov) %>%
  sqrt() %>%
  test_regvar(eta = 0.5, m = 10000, k = 80, mm = 10000)

# Compute extreme quantile estimates for innovations
innovation_comb <- stock %>%
  select(ends_with("innovation")) %>%
  combn(2, simplify = FALSE)

prediction_comb <- stock %>%
  select(ends_with("prediction")) %>%
  combn(2, simplify = FALSE)

m <- 1000
k <- 160
p <- c(low = 1 / 2000, medium = 1 / 5000, high = 1 / 10000)
labels <- c(us = "S&P 500", uk = "FTSE 100", jpn = "Nikkei 225")

for (i in 1:length(innovation_comb)) {
  countries <- c(str_split(names(innovation_comb[[i]]), "_",
                           simplify = TRUE))[1:2]
  sigma <- parameters %>%
    filter(country %in% countries) %>%
    pull(sigma_pred) %>%
    diag()
  
  offset <- parameters %>%
    filter(country %in% countries) %>%
    pull(offset)
  
  mcd_est <- covMcd(innovation_comb[[i]], alpha = 0.5)
  estimates <- map(p, ~ elliptical_extreme_qregion(innovation_comb[[i]],
                                                           mcd_est$center,
                                                           mcd_est$cov,
                                                           ., k, m)) %>%
    map(~ sweep(. %*% t(sigma), 2, offset, "+")) %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(group = rep(c("low", "medium", "high"), each = m)) %>%
    rename(x = V1, y = V2)
  
  colnames(prediction_comb[[i]]) <- c("x", "y")
  g <- ggplot(prediction_comb[[i]], aes(x = x, y = y)) +
    geom_point() +
    geom_path(data = estimates,
              aes(x = x, y = y, group = group, linetype = group),
              show.legend = FALSE) +
    coord_fixed() +
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 15),
          axis.text = element_text(size = 15)) +
    xlab(labels[countries[1]]) + ylab(labels[countries[2]]) +
    scale_linetype_manual(name = NULL,
                          values = c("low" = "solid", "medium" = "dashed",
                                     "high" = "dotted"))
  
  # Save figure
  filename <- str_c("figures/", str_c(countries, collapse = "_"), ".jpg")
  ggsave(filename, plot = g, width = 7, height = 7,
         dpi = 1000)
}
