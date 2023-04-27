suppressPackageStartupMessages(library(dplyr))
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
                                                          ., k, m)) %>%
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
    coord_fixed() +
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
                          values = c("low" = "solid", "medium" = "dashed",
                                     "high" = "dotted"))

  # Save figure
  filename <- stringr::str_c("figures/",
                             stringr::str_c(countries, collapse = "_"),
                             ".jpg")
  ggsave(filename, plot = g, width = 7, height = 7, dpi = 1000)
}
