library(lubridate)
library(forecast)
library(exuber)
library(dplyr)

radf_mc_cv2 <- function(tb, minw0, nrep, seed) {
  set.seed(seed)
  a <- 1 / tb
  distr <- numeric(nrep)
  
  for (i in 1:nrep) {
    y <- cumsum(rnorm(tb + 1) + a)
    distr[i] <- radf(y, minw0, lag = 0)$gsadf
  }
  
  quantile(distr, c(0.90, 0.95, 0.99))
}

locate_episodes <- function(bsadf, cv, dates, min_duration = 1) {
  dates <- as_date(dates)
  signal <- bsadf > cv
  
  start <- signal & !lag(signal, default = FALSE)
  end <- signal & !lead(signal, default = FALSE)
  starts <- dates[start]
  ends <- dates[end]
  
  episodes <- data.frame(start = starts,
                         end = ends,
                         duration = as.integer(ends - starts) + 1)
  
  episodes[episodes$duration >= min_duration, ]
}

data <- read.csv("close_price.csv")
data$datetime <- as_datetime(data$datetime) 
data[, 2:4] <- log(data[, 2:4])

auto.arima(data$btc, ic = "bic", max.q = 0, max.d = 1) 
auto.arima(data$eth, ic = "bic", max.q = 0, max.d = 1)
auto.arima(data$ltc, ic = "bic", max.q = 0, max.d = 1)

N <- nrow(data)
r0 <- 0.01 + 1.8 / sqrt(N)
minw0 <- floor(r0 * N)
tb <- minw0 + 180 - 1

bsadf_btc <- radf(data$btc, minw0, lag = 0)$bsadf
gsadf_btc <- max(bsadf_btc)

bsadf_eth <- radf(data$eth, minw0, lag = 0)$bsadf
gsadf_eth <- max(bsadf_eth)

bsadf_ltc <- radf(data$ltc, minw0, lag = 0)$bsadf
gsadf_ltc <- max(bsadf_ltc)

options(exuber.parallel = TRUE)
options(exuber.ncores = parallel::detectCores() - 1)

gsadf_mc_cv <- radf_mc_cv(N, minw0, nrep = 2000, seed = 123)$gsadf_cv
gsadf_wb_cv_btc <- radf_wb_cv(data$btc, minw0, nboot = 2000, seed = 123)$gsadf_cv
gsadf_wb_cv_eth <- radf_wb_cv(data$eth, minw0, nboot = 2000, seed = 123)$gsadf_cv
gsadf_wb_cv_ltc <- radf_wb_cv(data$ltc, minw0, nboot = 2000, seed = 123)$gsadf_cv

bsadf_mc_cv <- radf_mc_cv2(tb = tb, minw0, nrep = 2000, seed = 123)
bsadf_cb_cv_btc <- radf_wb_cv2(data$btc, minw0, nboot = 5000, tb = tb, seed = 123)$gsadf_cv
bsadf_cb_cv_eth <- radf_wb_cv2(data$eth, minw0, nboot = 5000, tb = tb, seed = 123)$gsadf_cv
bsadf_cb_cv_ltc <- radf_wb_cv2(data$ltc, minw0, nboot = 5000, tb = tb, seed = 123)$gsadf_cv

data.frame(datetime = data$datetime,
           bsadf_btc = c(rep(NA, minw0 - 1), bsadf_btc),
           bsadf_eth = c(rep(NA, minw0 - 1), bsadf_eth),
           bsadf_ltc = c(rep(NA, minw0 - 1), bsadf_ltc)) %>%
  write.csv(file = "results/bsadf.csv", row.names = FALSE)

saveRDS(list(gsadf = list(gsadf_mc_cv = as.matrix(gsadf_mc_cv), 
                          gsadf_wb_cv_btc = gsadf_wb_cv_btc, 
                          gsadf_wb_cv_eth = gsadf_wb_cv_eth,
                          gsadf_wb_cv_ltc = gsadf_wb_cv_ltc),
             bsadf = list(bsadf_mc_cv = as.matrix(bsadf_mc_cv), 
                          bsadf_cb_cv_btc = bsadf_cb_cv_btc, 
                          bsadf_cb_cv_eth = bsadf_cb_cv_eth,
                          bsadf_cb_cv_ltc = bsadf_cb_cv_ltc)), 
        file = "results/radf_cv.rds")

cv <- readRDS("results/radf_cv.rds")
bsadf <- read.csv("results/bsadf.csv")
bsadf$datetime <- as_datetime(bsadf$datetime) 
bsadf <- bsadf %>% na.omit()

episodes_btc_mc <- locate_episodes(bsadf$bsadf_btc, cv = cv$bsadf$bsadf_mc_cv[2], bsadf$datetime)
episodes_btc_cb <- locate_episodes(bsadf$bsadf_btc, cv = cv$bsadf$bsadf_cb_cv_btc[2], bsadf$datetime)

episodes_eth_mc <- locate_episodes(bsadf$bsadf_eth, cv = cv$bsadf$bsadf_mc_cv[2], bsadf$datetime)
episodes_eth_cb <- locate_episodes(bsadf$bsadf_eth, cv = cv$bsadf$bsadf_cb_cv_eth[2], bsadf$datetime)

episodes_ltc_mc <- locate_episodes(bsadf$bsadf_ltc, cv = cv$bsadf$bsadf_mc_cv[2], bsadf$datetime)
episodes_ltc_cb <- locate_episodes(bsadf$bsadf_ltc, cv = cv$bsadf$bsadf_cb_cv_ltc[2], bsadf$datetime)

saveRDS(list(btc = list(mc = episodes_btc_mc, cb = episodes_btc_cb),
             eth = list(mc = episodes_eth_mc, cb = episodes_eth_cb),
             ltc = list(mc = episodes_ltc_mc, cb = episodes_ltc_cb)),
        file = "results/episodes.rds")
