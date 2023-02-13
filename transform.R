library(readr)
library(dplyr)
library(rugarch)

# Read data and calculate daily returns
stock <- read_csv("data/stock_price.csv", col_names = TRUE) %>%
  rename(date = ...1,
         us = `US-SP500`,
         uk = `UK-FTSE100`,
         jpn = `JPN-Nikke225`) %>%
  mutate(across(-date, ~ log(. / lag(.)))) %>%
  slice(-1)


garch_spec <- ugarchspec(variance.model = list(model = "eGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0)),
                         distribution.model = "std")

fit_us <- ugarchfit(spec = garch_spec, data = stock$us)
res_us <- (stock$us - coef(fit_us)["mu"]) / sigma(fit_us)

View(stock)
