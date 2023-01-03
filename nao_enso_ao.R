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

fft <- torch_fft_fft(as.numeric(scale(nao$nao)))

num_samples <- nrow(nao)
nyquist_cutoff <- floor(num_samples/2) 
bins_below_nyquist <- 0:nyquist_cutoff

sampling_rate <- 12 # per year
frequencies_per_bin <- sampling_rate / num_samples
frequencies <- frequencies_per_bin * bins_below_nyquist

df <- data.frame(f = frequencies, y = as.numeric(fft[1:(nyquist_cutoff + 1)]$abs()))
df %>% ggplot(aes(f, y)) + geom_line() +
  xlab("frequency (per year)")

torch_topk(fft[1:(nyquist_cutoff/2)]$abs(), 9)
# dominant frequencies: 1/year, 2/year, 0.5/year, 0.08/year ...
# [[1]]
# torch_tensor
# 356.4380
# 318.9236
# 281.5578
# 230.8849
# 225.3260
# 212.9847
# 207.1933
# 205.6112
# 205.2927
# [ CPUFloatType{9} ]
# 
# [[2]]
# torch_tensor
# 200
# 399
# 398
# 269
# 356
# 424
# 89
# 400
# 16
# [ CPULongType{9} ]

nao %>%
  model(STL(nao ~ season(period = 144) + season(period = 24) + season(period = 12) + season(period = 6))) %>%
  components() %>%
  autoplot()

nao %>%
  model(STL(nao ~ season(period = 12 * 7))) %>%
  components() %>%
  autoplot()


nao %>% features(nao, feat_stl) %>% glimpse()
feat_stl(nao$nao, .period = 12) %>% round(2)
feat_stl(nao$nao, .period = 6) %>% round(2)
feat_stl(nao$nao, .period = 24) %>% round(2)

nao %>% ACF(nao) %>% autoplot()





########################   ENSO   ########################

enso <- read_table("data/ONI_NINO34_1854-2022.txt", skip = 9) %>%
  mutate(x = yearmonth(as.Date(paste0(YEAR, "-", `MON/MMM`, "-01")))) %>%
  select(x, enso = NINO34_MEAN) %>%
  filter(x >= yearmonth("1854-01"), x <= yearmonth("2022-11")) %>%
  as_tsibble(index = x) 

fft <- torch_fft_fft(as.numeric(scale(enso$enso)))

num_samples <- nrow(enso)
nyquist_cutoff <- floor(num_samples / 2)
bins_below_nyquist <- 0:nyquist_cutoff

sampling_rate <- 12 # per year
frequencies_per_bin <- sampling_rate / num_samples
frequencies <- frequencies_per_bin * bins_below_nyquist

df <- data.frame(f = frequencies, y = as.numeric(fft[1:(nyquist_cutoff + 1)]$abs()))
df %>% ggplot(aes(f, y)) + geom_line() +
  xlab("frequency (per year)")

torch_topk(fft[1:(nyquist_cutoff/2)]$abs(), 9)

########################   AO    ########################


