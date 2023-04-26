library(readr)
library(dplyr)
library(rugarch)
library(rlang)
library(stringr)
library(glue)
library(purrr)


#' Compute Innovations corresponding to the Garch model
#'
#' @param fit Fitted GARCH model.
#' @param x Data as a double vector.
#'
#' @return Double vector of innovations
get_innovations <- function(fit, x) {
  as.double((x - coef(fit)["mu"]) / sigma(fit))
}


#' Compute one step predictions for the GARCH model
#'
#' @param fit Fitted GARCH model.
#' @param x Data as a double vector.
#'
#' @return Double vector of predictions.
get_predictions <- function(fit, x) {
  sigma_pred <- as.double(sigma(ugarchforecast(fit, n.ahead = 1)))
  as.double(x * sigma_pred + coef(fit)["mu"])
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


# Calculate daily returns
stock <- read_csv("data/stock-price-raw.csv", col_names = TRUE) %>%
  rename(date = ...1,
         us = `US-SP500`,
         uk = `UK-FTSE100`,
         jpn = `JPN-Nikke225`) %>%
  mutate(across(c(us, uk, jpn),
                ~ log(. / lag(.)),
                .names = "{col}_return")) %>%
  slice(-1)

# Specification for the EGARCH(1, 1) model
garch_spec <- ugarchspec(variance.model = list(model = "eGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0)),
                         distribution.model = "std")

# Fit EGARCH for each country
fits <- stock %>%
  select(ends_with("return")) %>%
  purrr::map(~ ugarchfit(garch_spec, .x))

# Compute innovations and predictions
for (country in c("us", "uk", "jpn")) {
  cols <- paste0(country, c("_return", "_innovation", "_prediction"))
  stock <- stock %>%
    mutate("{cols[2]}" := get_innovations(fits[[cols[1]]], get(cols[1])),
           "{cols[3]}" := get_predictions(fits[[cols[1]]], get(cols[2])))
}

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
