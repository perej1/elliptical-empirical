# Estimation of extreme quantile regions
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)


#' Convert a spherical coordinate to a Cartesian coordinate
#'
#' @param radius Double, radius of the point.
#' @param theta Double vector, d - 1 angular coordinates.
#'
#' @return Double vector, d coordinates giving the location of the point in
#'   space.
spherical_to_cartesian <- function(radius, theta) {
  d <- length(theta) + 1
  cartesian <- rep(NA, d)
  cartesian[1] <- cos(theta[1])
  cartesian[d] <- prod(sin(theta))
  
  if (d > 2) {
    for (i in 2:(d - 1)) {
      cartesian[i] <- prod(sin(theta[1:(i - 1)])) * cos(theta[i])
    }
  }
  radius * cartesian
}


#' Generate m points from a (d - 1)-sphere
#'
#' @param d Integer, dimension of the ball, must be equal to or greater than 2.
#' @param m_angle Integer, number of points to return.
#'
#' @return List of two double matrices, the first gives the ball in spherical
#'   and the second in Cartesian coordinates. One row represents one point.
get_ball_mesh <- function(d, m_angle) {
  m <- ceiling(m_angle^(1 / (d - 1)))
  angle_list <- vector("list", d - 1)
  angle_list[[d - 1]] <- seq(0, 2 * pi, length.out = m)
  if (d > 2) {
    for (i in 1:(d - 2)) {
      angle_list[[i]] <- seq(0, pi, length.out = m)
    }
  }
  spherical <- as.matrix(expand.grid(angle_list))
  list(spherical = spherical,
       cartesian = t(apply(spherical, 1, spherical_to_cartesian, r = 1)),
       m_effective = nrow(spherical))
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
#' Consistent for heavy-tailed elliptical distributions under some technical
#' assumptions.
#'
#' @param data Double matrix of observations, each row represents one
#'   observation.
#' @param mu_est Double vector, estimate of the location.
#' @param sigma_est Double matrix, estimate of the scatter.
#' @param p Double, probability in quantile region.
#' @param k Integer, threshold for the sample from the tail.
#' @param m_angle Integer, number of points to return.
#'
#' @return Double matrix, m_angle points from the boundary of the quantile
#'   region.
elliptical_extreme_qregion <- function(data, mu_est, sigma_est, p, k, m_angle) {
  n <- nrow(data)
  d <- ncol(data)
  w <- get_ball_mesh(d, m_angle)$cartesian
  
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
  data <- sweep(w %*% t(lambda), 2, mu_est, "+")
  list(data = data, r_hat = r_hat)
}


# Read data
stock <- readr::read_csv("data/stock-price-mod.csv", col_names = TRUE,
                         col_types = readr::cols())
parameters <- readr::read_csv("data/model-parameters.csv", col_names = TRUE,
                              col_types = readr::cols())

# Get bivariate combinations for innovations and predictions
innovation_comb <- stock %>%
  select(ends_with("innovation")) %>%
  combn(2, simplify = FALSE)

prediction_comb <- stock %>%
  select(ends_with("prediction")) %>%
  combn(2, simplify = FALSE)

# Compute predicted bivariate extreme quantile regions
m <- 1000
k <- 160
p <- c(low = 1 / 2000, medium = 1 / 5000, high = 1 / 10000)
labels <- c(us = "S&P 500", uk = "FTSE 100", jpn = "Nikkei 225")

for (i in seq_along(innovation_comb)) {
  countries <- c(stringr::str_split(names(innovation_comb[[i]]), "_",
                           simplify = TRUE))[1:2]
  sigma <- parameters %>%
    filter(country %in% countries) %>%
    pull(sigma_pred) %>%
    diag()

  offset <- parameters %>%
    filter(country %in% countries) %>%
    pull(offset)

  mcd_est <- robustbase::covMcd(innovation_comb[[i]], alpha = 0.5)
  estimates <- purrr::map(p, ~ elliptical_extreme_qregion(innovation_comb[[i]],
                                                          mcd_est$center,
                                                          mcd_est$cov,
                                                          ., k, m)$data) %>%
    purrr::map(~ sweep(. %*% t(sigma), 2, offset, "+")) %>%
    do.call(rbind, .) %>%
    as_tibble(.name_repair = "universal") %>%
    mutate(group = rep(c("low", "medium", "high"), each = m)) %>%
    rename(x = `...1`, y = `...2`)

  colnames(prediction_comb[[i]]) <- c("x", "y")
  g <- ggplot(prediction_comb[[i]], aes(x = x, y = y)) +
    geom_point() +
    geom_path(data = estimates,
              aes(x = x, y = y, group = group, linetype = group),
              show.legend = FALSE) +
    coord_fixed() + xlim(-0.07, 0.07) + ylim(-0.07, 0.07) +
    theme(axis.title = element_text(size = 15, face = "bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 15),
          axis.text = element_text(size = 15)) +
    xlab(labels[countries[1]]) + ylab(labels[countries[2]]) +
    scale_linetype_manual(name = NULL,
                          values = c("low" = "solid", "medium" = "solid",
                                     "high" = "solid"))

  # Save figure
  filename <- stringr::str_c("figures/",
                             stringr::str_c(countries, collapse = "_"),
                             "_estimate", ".jpg")
  ggsave(filename, plot = g, width = 7, height = 7, dpi = 1000)
}

# Outlier detection with three dimensional extreme quantile regions
innovations <- stock %>%
  select(ends_with("innovation"))

mcd_est <- robustbase::covMcd(innovations, alpha = 0.5)

r_hats <- purrr::map_dbl(p,
                         ~ elliptical_extreme_qregion(innovations,
                                                      mu_est = mcd_est$center,
                                                      sigma_est = mcd_est$cov,
                                                      .,
                                                      k = k,
                                                      m = 10000)$r_hat)

outliers <- purrr::map(r_hats,
                       ~ mahalanobis(innovations,
                                     mcd_est$center,
                                     mcd_est$cov) >= .^2) %>%
  purrr::map(~ stock$date[.])

# Print results of outlier detection
cli::cli_alert_info("Estimates of the location and scatter:")
mcd_est$center
mcd_est$cov

cli::cli_h3("Results of outlier detection for different levels of p when k = {k}")
cli::cli_alert_info("p = {p['low']}:
                    \t - r_hat = {r_hats['low']}
                    \t - {length(outliers$low)} outlier{?s} found, {outliers$low}")

cli::cli_alert_info("p = {p['medium']}:
                    \t - r_hat = {r_hats['medium']}
                    \t - {length(outliers$medium)} outlier{?s} found, {outliers$medium}")

cli::cli_alert_info("p = {p['high']}:
                    \t - r_hat = {r_hats['high']}
                    \t - {length(outliers$high)} outlier{?s} found, {outliers$high}")
