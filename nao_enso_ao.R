########################   Data   ########################  
# NAO
# https://crudata.uea.ac.uk/cru/data//nao/index.htm
# 
# The NAO is traditionally defined as the normalized pressure difference between a station on the Azores
# and one on Iceland. [...] Here we give data for SW Iceland (Reykjavik), Gibraltar and Ponta Delgada (Azores). 
# cite: Jones, P.D., Jónsson, T. and Wheeler, D., 1997: Extension to the North Atlantic Oscillation using early instrumental pressure observations from Gibraltar and South-West Iceland. Int. J. Climatol. 17, 1433-1450. doi: 10.1002/(SICI)1097-0088(19971115)17:13<1433::AID-JOC203>3.0.CO;2-P
#
#
# ENSO
# https://bmcnoldy.rsmas.miami.edu/tropics/oni/
#
# ONI ≥ 0.5°C indicate El Niño, values ≤ -0.5°C indicate La Niña, and values in between indicate a Neutral phase.
#
#
# AO
# https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/ao.shtml
# https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/history/method.shtml


########################   Links   ########################
# https://psl.noaa.gov/data/correlation/
# https://ore.exeter.ac.uk/repository/handle/10871/34601 (Skillful long-range prediction of European and North American winters)
# Previous studies have shown that the El Niño–Southern Oscillation can drive interannual variations
# in the NAO [Brönnimann et al., 2007] and hence Atlantic and European winter climate via the
# stratosphere [Bell et al., 2009].Figures 2b and 2c confirm that this teleconnection to the tropical
# Pacific is active in our experiments, with forecasts initialized in El Niño/La Niña conditions in
# November tending to be followed by negative/positive NAO conditions in winter. 


library(torch)
library(tidyverse)
library(fable)
library(tsibble)
library(feasts)


########################   NAO   ########################

# start 3 years later, in 1824
# last valid value is 2022-9
nao <-
  read_table(
    "data/nao.dat",
    col_names = FALSE,
    na = "-99.99",
    skip = 3
  ) %>%
  select(-X1,-X14) %>%
  as.matrix() %>% 
  t() %>%
  as.vector() %>%
  .[1:(length(.) - 3)] %>%
  tibble(x = seq.Date(
    from = as.Date("1824-01-01"),
    to = as.Date("2022-09-01"),
    by = "months"
  ),
  nao = .) %>%
  mutate(x = yearmonth(x)) %>%
  fill(nao) %>%
  as_tsibble(index = x) 

fft <- torch_fft_fft(nao$nao)
below_nyquist <- fft[1:(length(fft)/2)]

### tbd correct!!!! ###

years <- length(below_nyquist)/12
x <- seq(0, years, length.out = length(below_nyquist))/years

df <- data.frame(f = x, y = as.numeric(below_nyquist$abs()))
df %>% ggplot(aes(f, y)) + geom_line() +
  xlab("frequency (years)")

nao %>%
  model(STL(nao ~ season(period = 12 * 7))) %>%
  components() %>%
  autoplot()
nao %>%
  filter(x >= yearmonth("1950-01")) %>%
  model(STL(nao ~ season(period = 12 * 7))) %>%
  components() %>%
  autoplot()


########################   tbd adjust   ########################
nao_valid %>% features(y, feat_stl)
feat_stl(nao_valid$y, .period = 12, s.window = 7) %>% round(2)

# ACF
nao_valid %>% features(y, feat_acf)
nao_valid %>% ACF(y) %>% autoplot()

# other features
# rate at which autocorrelations decrease as the lag between pairs of values increases
# > 0.5: long-term positive autocorrelations
# < 0.5: mean-reverting
# 0.5: random walk
nao_valid %>% features(y, coef_hurst) 
nao_valid %>% features(y, feat_spectral) #[0, 1]
nao_valid %>% features(y, feat_acf)
nao_valid %>% features(y, feat_acf)
nao_valid %>% features(y, feat_acf)


# fit

fit <- nao_train %>% model(
  # Error ={A,M}, Trend ={N,A,Ad} and Seasonal ={N,A,M}.
  ets = ETS(y ~ season(method = "A", gamma = 0.5))
)
fit

fit <- nao_train %>% model(
  ets = ETS(y ~ season(method = "A", gamma = 0.1)), # 0: seasonal pattern will not change
  ets2 = ETS(y ~ season(method = "A", gamma = 0.5)),
  ets3 = ETS(y ~ season(method = "A", gamma = 0.9)), # 1: seasonality will have no memory of past periods
  arima = ARIMA(y),
  snaive = SNAIVE(y)
) 

fc <- fit %>%
  forecast(h = "2 years") 

fc %>% 
  autoplot(filter(nao_valid, x < yearmonth("1992-01")), level = NULL)

accuracy(fc, filter(nao_valid, x < yearmonth("1992-01")))

fit %>% select(ets3) %>% report()
fit %>% report()





########################   ENSO   ########################

enso <- read_table2("data/ONI_NINO34_1854-2020.txt", skip = 9) %>%
  mutate(x = yearmonth(as.Date(paste0(YEAR, "-", `MON/MMM`, "-01")))) %>%
  select(x, enso = NINO34_MEAN) %>%
  filter(x >= yearmonth("1854-01"), x <= yearmonth("2020-07")) %>%
  as_tsibble(index = x) 




