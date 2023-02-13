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

# Specification for the EGARCH(1, 1) model
garch_spec <- ugarchspec(variance.model = list(model = "eGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0)),
                         distribution.model = "std")

# Read data, calculate daily returns and innovations, write data
read_csv("data/stock-price-raw.csv", col_names = TRUE) %>%
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
                .names = "{col}_innovation")) %>%
  write_csv("data/stock-price-mod.csv")
