# Compute returns, fit EGARCH(1, 1) models and compute innovations
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rugarch))


#' Compute innovations corresponding to the Garch model
#'
#' @param fit Fitted GARCH model of class uGARCHfit.
#' @param x Double vector of raw returns.
#'
#' @return Double vector of innovations.
get_innovations <- function(fit, x) {
  as.double((x - coef(fit)["mu"]) / sigma(fit))
}


#' Compute predicted returns
#'
#' @param fit Fitted GARCH model of class uGARCHfit.
#' @param x Double vector of innovations.
#'
#' @return Double vector of predicted returns.
get_predictions <- function(fit, x) {
  sigma_pred <- as.double(sigma(ugarchforecast(fit, n.ahead = 1)))
  as.double(x * sigma_pred + coef(fit)["mu"])
}


# Calculate daily returns
stock <- readr::read_csv("data/stock-price-raw.csv", col_names = TRUE,
                         col_types = readr::cols()) %>%
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

# Write data
stock %>%
  readr::write_csv("data/stock-price-mod.csv")

# Write EGARCH model parameters
offset <- purrr::map_dbl(fits, ~ coef(.)["mu"])
sigma_pred <- purrr::map_dbl(fits,
                             ~ as.double(
                               sigma(ugarchforecast(., n.ahead = 1))
                               )
                             )
country <- purrr::map_chr(names(fits),
                          ~ stringr::str_split(., "_", simplify = TRUE)[1])
tibble(offset, sigma_pred, country) %>%
  readr::write_csv("data/model-parameters.csv")
