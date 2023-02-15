library(readr)
library(dplyr)
library(purrr)
library(ellipticalsymmetry)

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
