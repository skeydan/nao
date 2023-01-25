########################   Data   ########################  
# NAO
# https://crudata.uea.ac.uk/cru/data//nao/index.htm
# 
# ENSO
# https://bmcnoldy.rsmas.miami.edu/tropics/oni/
#
# AO
# https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/ao.shtml


########################   Links   ########################
# https://psl.noaa.gov/data/correlation/
# https://ore.exeter.ac.uk/repository/handle/10871/34601 
#

library(torch)
library(tidyverse)
library(tsibble)
library(feasts)


########################   NAO   ########################

download.file(
  "https://crudata.uea.ac.uk/cru/data//nao/nao.dat",
  destfile = "data/nao.dat"
)

# last valid value is 2022-9
# for all 3, use 1950-01 to 2022-09

use_months <- seq.Date(
  # from = as.Date("1824-01-01"),
  from = as.Date("1950-01-01"),
  to = as.Date("2022-09-01"),
  by = "months"
)
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
  .[1:length(use_months)] %>%
  tibble(x = use_months,
  nao = .) %>%
  mutate(x = yearmonth(x)) %>%
  fill(nao) %>%
  as_tsibble(index = x) 

fft <- torch_fft_fft(as.numeric(scale(nao$nao)))

num_samples <- nrow(nao)
nyquist_cutoff <- ceiling(num_samples / 2)
bins_below_nyquist <- 0:nyquist_cutoff

sampling_rate <- 12 # per year
frequencies_per_bin <- sampling_rate / num_samples
frequencies <- frequencies_per_bin * bins_below_nyquist

df <- data.frame(f = frequencies, y = as.numeric(fft[1:(nyquist_cutoff + 1)]$abs()))
df %>% ggplot(aes(f, y)) +
  geom_line() +
  xlab("frequency (per year)") +
  ylab("magnitude") +
  ggtitle("Spectrum of NAO data.")

strongest <- torch_topk(fft[1:(nyquist_cutoff/2)]$abs(), 6)
strongest
# [[1]]
# torch_tensor
# 102.7191
# 80.5129
# 76.1179
# 75.9949
# 72.9086
# 60.8281
# [ CPUFloatType{6} ]
# 
# [[2]]
# torch_tensor
# 147
# 99
# 146
# 59
# 33
# 78
# [ CPULongType{6} ]

important_freqs <- frequencies[as.numeric(strongest[[2]])]
important_freqs
# [1] 2.0068729 1.3470790 1.9931271 0.7972509 0.4398625 1.0584192
num_observations_in_season <- 12/important_freqs  
num_observations_in_season
# [1]  5.979452  8.908163  6.020690 15.051724 27.281250 11.337662
nao %>%
  model(STL(nao ~ season(period = 6) + season(period = 9) + season(period = 15) + season(period = 27) + season(period = 12))) %>%
  components() %>%
  autoplot()

nao %>%
  autoplot()

nao %>%
  filter(x >= yearmonth("2000-01")) %>%
  autoplot()

nao %>% ACF(nao) %>% autoplot()





########################   ENSO   ########################
download.file(
  "https://bmcnoldy.rsmas.miami.edu/tropics/oni/ONI_NINO34_1854-2022.txt",
  destfile = "ONI_NINO34_1854-2022.txt"
)

enso <- read_table("data/ONI_NINO34_1854-2022.txt", skip = 9) %>%
  mutate(x = yearmonth(as.Date(paste0(YEAR, "-", `MON/MMM`, "-01")))) %>%
  select(x, enso = NINO34_MEAN) %>%
  filter(x >= yearmonth("1950-01"), x <= yearmonth("2022-09")) %>%
  as_tsibble(index = x) 

fft <- torch_fft_fft(as.numeric(scale(enso$enso)))

num_samples <- nrow(enso)
nyquist_cutoff <- ceiling(num_samples / 2)
bins_below_nyquist <- 0:nyquist_cutoff

sampling_rate <- 12 # per year
frequencies_per_bin <- sampling_rate / num_samples
frequencies <- frequencies_per_bin * bins_below_nyquist

df <- data.frame(f = frequencies, y = as.numeric(fft[1:(nyquist_cutoff + 1)]$abs()))
df %>% ggplot(aes(f, y)) +
  geom_line() +
  xlab("frequency (per year)") +
  ylab("magnitude") +
  ggtitle("Spectrum of Niño 3.4 data")

strongest <- torch_topk(fft[1:(nyquist_cutoff/2)]$abs(), 3)
strongest
# [[1]]
# torch_tensor
# 233.9855
# 172.2784
# 142.3784
# [ CPUFloatType{9} ]
# 
# [[2]]
# torch_tensor
# 74
# 21
# 7
# [ CPULongType{9} ]
important_freqs <- frequencies[as.numeric(strongest[[2]])]
important_freqs
# [1] 1.00343643 0.27491409 0.08247423 
num_observations_in_season <- 12/important_freqs  
num_observations_in_season
# [1] 11.95890  43.65000 145.50000  

enso %>%
  model(STL(enso ~ season(period = 12) + season(period = 44) + season(period = 145))) %>%
  components() %>%
  autoplot()

enso %>%
  autoplot()

enso %>%
  filter(x >= yearmonth("2000-01")) %>%
  autoplot()

enso %>% ACF(enso) %>% autoplot()



########################   AO    ########################

download.file(
  "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/monthly.ao.index.b50.current.ascii.table",
  destfile = "data/ao.dat"
)

ao <-
  read_table(
    "data/ao.dat",
    col_names = FALSE,
    skip = 1
  ) %>%
  select(-X1) %>%
  as.matrix() %>% 
  t() %>%
  as.vector() %>%
  .[1:length(use_months)] %>%
  tibble(x = use_months,
         ao = .) %>%
  mutate(x = yearmonth(x)) %>%
  fill(ao) %>%
  as_tsibble(index = x) 

fft <- torch_fft_fft(as.numeric(scale(ao$ao)))

num_samples <- nrow(ao)
nyquist_cutoff <- ceiling(num_samples / 2)
bins_below_nyquist <- 0:nyquist_cutoff

sampling_rate <- 12 # per year
frequencies_per_bin <- sampling_rate / num_samples
frequencies <- frequencies_per_bin * bins_below_nyquist

df <- data.frame(f = frequencies, y = as.numeric(fft[1:(nyquist_cutoff + 1)]$abs()))
df %>% ggplot(aes(f, y)) +
  geom_line() +
  xlab("frequency (per year)") +
  ylab("magnitude") +
  ggtitle("Spectrum ofAO data")

strongest <- torch_topk(fft[1:(nyquist_cutoff/2)]$abs(), 5)
strongest
# [[1]]
# torch_tensor
# 77.3358
# 73.1349
# 72.3889
# 72.3759
# 72.2528
# [ CPUFloatType{5} ]
# 
# [[2]]
# torch_tensor
# 2
# 27
# 130
# 94
# 6
# [ CPULongType{5} ]

important_freqs <- frequencies[as.numeric(strongest[[2]])]
important_freqs
# [1] 0.01374570 0.35738832 1.77319588 1.27835052 0.06872852
num_observations_in_season <- 12/important_freqs  
num_observations_in_season
# [1] 873.000000  33.576923   6.767442   9.387097 174.600000 
ao %>%
  model(STL(ao ~ season(period = 33) + season(period = 7) + season(period = 9) + season(period = 174))) %>%
  components() %>%
  autoplot()

ao %>%
  autoplot()

ao %>%
  filter(x >= yearmonth("2000-01")) %>%
  autoplot()

ao %>% ACF(ao) %>% autoplot()


########################   wavelet analysis    ########################

library(torchwavelets)

### nao ###

nao_idx <- nao$nao %>% as.numeric() %>% torch_tensor()
dt <- 1/12 # same interval as for ENSO
wtf <- wavelet_transform(length(nao_idx), dt = dt)
power_spectrum <- wtf$power(nao_idx)

times <- lubridate::year(nao$x) + lubridate::month(nao$x)/12 # also same
scales <- as.numeric(wtf$scales) # will be same because both series have same length

df <- as_tibble(as.matrix(power_spectrum$t()), .name_repair = "universal") %>%
  mutate(time = times) %>%
  pivot_longer(!time, names_to = "scale", values_to = "power") %>%
  mutate(scale = scales[scale %>% 
                          str_remove("[\\.]{3}") %>%
                          as.numeric()])

coi <- wtf$coi(times[1], times[length(nao_idx)])
coi_df <- data.frame(x = as.numeric(coi[[1]]), y = as.numeric(coi[[2]]))

labeled_scales <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64) # same since scales are same 
labeled_frequencies <- round(as.numeric(wtf$fourier_period(labeled_scales)), 1)

ggplot(df) +
  scale_y_continuous(
    trans = scales::compose_trans(scales::log2_trans(), scales::reverse_trans()),
    breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64),
    limits = c(max(scales), min(scales)),
    expand = c(0,0),
    sec.axis = dup_axis(
      labels = scales::label_number(labeled_frequencies),
      name = "Fourier period (years)")
  ) + 
  ylab("scale (years)") +
  scale_x_continuous(breaks = seq(1950, 2020, by = 5), expand = c(0,0)) +
  xlab("year") +
  geom_contour_filled(aes(time, scale, z = power), show.legend = FALSE) +
  scale_fill_viridis_d(option = "turbo") +
  geom_ribbon(data = coi_df, aes(x = x, ymin = y, ymax = max(scales)), fill = "black", alpha = 0.6) +
  theme(legend.position = "none")


### enso ###


enso_idx <- enso$enso %>% as.numeric() %>% torch_tensor()
dt <- 1/12
wtf <- wavelet_transform(length(enso_idx), dt = dt)
power_spectrum <- wtf$power(enso_idx)

times <- lubridate::year(enso$x) + lubridate::month(enso$x)/12
scales <- as.numeric(wtf$scales)

df <- as_tibble(as.matrix(power_spectrum$t()), .name_repair = "universal") %>%
  mutate(time = times) %>%
  pivot_longer(!time, names_to = "scale", values_to = "power") %>%
  mutate(scale = scales[scale %>% 
                          str_remove("[\\.]{3}") %>%
                          as.numeric()])
coi <- wtf$coi(times[1], times[length(enso_idx)])
coi_df <- data.frame(x = as.numeric(coi[[1]]), y = as.numeric(coi[[2]]))

labeled_scales <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)
labeled_frequencies <- round(as.numeric(wtf$fourier_period(labeled_scales)), 1)

ggplot(df) +
  scale_y_continuous(
    trans = scales::compose_trans(scales::log2_trans(), scales::reverse_trans()),
    breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64),
    limits = c(max(scales), min(scales)),
    expand = c(0,0),
    sec.axis = dup_axis(
      labels = scales::label_number(labeled_frequencies),
      name = "Fourier period (years)")
  ) + 
  ylab("scale (years)") +
  scale_x_continuous(breaks = seq(1950, 2020, by = 5), expand = c(0,0)) +
  xlab("year") +
  geom_contour_filled(aes(time, scale, z = power), show.legend = FALSE) +
  scale_fill_viridis_d(option = "turbo") +
  geom_ribbon(data = coi_df, aes(x = x, ymin = y, ymax = max(scales)), fill = "black", alpha = 0.6) +
  theme(legend.position = "none")



### ao ###

ao_idx <- ao$ao %>% as.numeric() %>% torch_tensor()
dt <- 1/12 # same interval as for ENSO and NAO
wtf <- wavelet_transform(length(ao_idx), dt = dt)
power_spectrum <- wtf$power(ao_idx)

times <- lubridate::year(ao$x) + lubridate::month(ao$x)/12 # also same
scales <- as.numeric(wtf$scales) # will be same because all series have same length

df <- as_tibble(as.matrix(power_spectrum$t()), .name_repair = "universal") %>%
  mutate(time = times) %>%
  pivot_longer(!time, names_to = "scale", values_to = "power") %>%
  mutate(scale = scales[scale %>% 
                          str_remove("[\\.]{3}") %>%
                          as.numeric()])

coi <- wtf$coi(times[1], times[length(ao_idx)])
coi_df <- data.frame(x = as.numeric(coi[[1]]), y = as.numeric(coi[[2]]))

labeled_scales <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64) # same since scales are same 
labeled_frequencies <- round(as.numeric(wtf$fourier_period(labeled_scales)), 1)

ggplot(df) +
  scale_y_continuous(
    trans = scales::compose_trans(scales::log2_trans(), scales::reverse_trans()),
    breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64),
    limits = c(max(scales), min(scales)),
    expand = c(0,0),
    sec.axis = dup_axis(
      labels = scales::label_number(labeled_frequencies),
      name = "Fourier period (years)")
  ) + 
  ylab("scale (years)") +
  scale_x_continuous(breaks = seq(1950, 2020, by = 5), expand = c(0,0)) +
  xlab("year") +
  geom_contour_filled(aes(time, scale, z = power), show.legend = FALSE) +
  scale_fill_viridis_d(option = "turbo") +
  geom_ribbon(data = coi_df, aes(x = x, ymin = y, ymax = max(scales)), fill = "black", alpha = 0.6) +
  theme(legend.position = "none")


