library(readr)
library(dplyr)
library(rugarch)

#' Compute Innovations corresponding to the Garch Model
#'
#' @param spec Specification of the garch model as uGARCHspec.
#' @param x Data as a double vector.
#'
#' @return Double vector of innovations
get_innovations <- function(spec, x) {
  fit <- ugarchfit(spec, x)
  as.double((x - coef(fit)["mu"]) / sigma(fit))
}

#' Three Dimensional Polar Coordinate from Cartesian Coordinate
#'
#' @param x Cartesian coordinate as a double vector
#'
#' @return Named vector of polar coordinates
car3pol <- function(x) {
  x <- unname(x)
  r <- norm(x, type = "2")
  inc <- asin(x[3] / r)
  azi <- sign(x[2]) * acos(x[1] / norm(x[1:2], type = "2"))
  c(radius = r, inclination = inc, azimuth = azi)
}

# Specification for the EGARCH(1, 1) model
garch_spec <- ugarchspec(variance.model = list(model = "eGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0)),
                         distribution.model = "std")

# Read data, calculate daily returns and innovations
stock <- read_csv("data/stock-price-raw.csv", col_names = TRUE) %>%
  rename(date = ...1,
         us = `US-SP500`,
         uk = `UK-FTSE100`,
         jpn = `JPN-Nikke225`) %>%
  mutate(across(c(us, uk, jpn),
                ~ log(. / lag(.)),
                .names = "{col}_return")) %>%
  slice(-1) %>%
  mutate(across(c(us_return, uk_return, jpn_return),
                ~ get_innovations(garch_spec, .),
                .names = "{col}_innovation"))

# Calculate polar coordinates for innovations
polar <- stock %>%
  select(ends_with("innovation")) %>%
  pmap(., lift_vd(car3pol)) %>%
  transpose() %>%
  map(flatten_dbl) %>%
  as_tibble()

# Write data
bind_cols(stock, polar) %>%
  write_csv("data/stock-price-mod.csv")
