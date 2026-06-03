library(lubridate)
library(forecast)
library(dplyr)

estimate_variance <- function(y, k, b1, b2, sd0) {
  N <- length(y)
  dy <- c(NA, diff(y))
  i <- seq(N)
  
  h1 <- b1 * N^(-1/3) / log(N)
  h2 <- b2 * N^(-1/4)
  trunc_thresh <- sd0 * log(N)
  
  e_hat <- rep(NA, N)
  
  for (t in (k + 2):N) {
    idx <- which(abs((i - t) / (N * h1)) <= 1 & i >= k + 2)
    
    X <- cbind(y[idx - 1], 
               if (k > 0) sapply(seq(k), function(j) dy[idx - j]))
    
    fit <- lm(dy[idx] ~ X - 1)
    e_hat[t] <- residuals(fit)[which(idx == t)]
  }
  
  e_trunc <- ifelse(abs(e_hat) < trunc_thresh, e_hat, NA)
  
  var_hat <- rep(NA, N)
  
  for (t in 1:N) {
    idx <- which(abs((i - t) / (N * h2)) <= 1)
    
    var_hat[t] <- sum(e_trunc[idx]^2, na.rm = TRUE) / sum(!is.na(e_trunc[idx]))
  }
  
  return(list(var_hat = var_hat, e_trunc = e_trunc))
}

compute_CV1 <- function(y, k, b1) {
  N <- length(y)
  dy <- c(NA, diff(y))
  i <- seq(N)
  
  h1 <- b1 * N^(-1/3) / log(N)
  e_hat_loo <- rep(NA, N)
  
  for (t in (k + 2):N) {
    idx <- which(abs((i - t) / (N * h1)) <= 1 & i >= k + 2 & i != t)
    
    X <- cbind(y[idx - 1], 
               if (k > 0) sapply(seq(k), function(j) dy[idx - j]))
    fit <- lm(dy[idx] ~ X - 1)
    x_t <- matrix(c(y[t - 1], 
                    if (k > 0) sapply(seq(k), function(j) dy[t - j])), nrow = 1)
    
    e_hat_loo[t] <- dy[t] - x_t %*% coef(fit)
  }
  
  CV1 <- mean(e_hat_loo^2, na.rm = TRUE)
  return(CV1)
}

compute_CV2_loo <- function(y, k, b1, b2, sd0) {
  N <- length(y)
  dy <- c(NA, diff(y))
  i <- seq(N)
  
  h1 <- b1 * N^(-1/3) / log(N)
  h2 <- b2 * N^(-1/4)
  trunc_thresh <- sd0 * log(N)
  
  e_hat <- rep(NA, N)
  
  for (t in (k + 2):N) {
    idx <- which(abs((i - t) / (N * h1)) <= 1 & i >= k + 2)
    
    X <- cbind(y[idx - 1], 
               if (k > 0) sapply(seq(k), function(j) dy[idx - j]))
    
    fit <- lm(dy[idx] ~ X - 1)
    e_hat[t] <- residuals(fit)[which(idx == t)]
  }
  
  e_trunc <- ifelse(abs(e_hat) < trunc_thresh, e_hat, NA)
  
  var_hat_loo <- rep(NA, N)
  
  for (t in 1:N) {
    idx <- which(abs((i - t) / (N * h2)) <= 1 & i != t)
    
    var_hat_loo[t] <- sum(e_trunc[idx]^2, na.rm = TRUE) / sum(!is.na(e_trunc[idx]))
  }
  
  CV2 <- mean((e_trunc^2 - var_hat_loo)^2, na.rm = TRUE)
  return(CV2)
}

compute_CV2_lpo <- function(y, k, b1, b2, sd0, p) {
  N <- length(y)
  dy <- c(NA, diff(y))
  i <- seq(N)
  
  h1 <- b1 * N^(-1/3) / log(N)
  h2 <- b2 * N^(-1/4)
  trunc_thresh <- sd0 * log(N)
  
  e_hat <- rep(NA, N)
  
  for (t in (k + 2):N) {
    idx <- which(abs((i - t) / (N * h1)) <= 1 & i >= k + 2)
    
    X <- cbind(y[idx - 1], 
               if (k > 0) sapply(seq(k), function(j) dy[idx - j]))
    
    fit <- lm(dy[idx] ~ X - 1)
    e_hat[t] <- residuals(fit)[which(idx == t)]
  }
  
  e_trunc <- ifelse(abs(e_hat) < trunc_thresh, e_hat, NA)
  
  var_hat_lpo <- rep(NA, N)
  
  for (t in 1:N) {
    idx <- which(abs((i - t) / (N * h2)) <= 1 & abs(i - t) > p)
    
    var_hat_lpo[t] <- sum(e_trunc[idx]^2, na.rm = TRUE) / sum(!is.na(e_trunc[idx]))
  }
  
  CV2 <- mean((e_trunc^2 - var_hat_lpo)^2, na.rm = TRUE)
  return(CV2)
}

data <- read.csv("close_price.csv")
data[, 2:4] <- log(data[, 2:4])

auto.arima(data$btc, ic = "bic", max.q = 0, max.d = 1) 
auto.arima(data$eth, ic = "bic", max.q = 0, max.d = 1)
auto.arima(data$ltc, ic = "bic", max.q = 0, max.d = 1)

sd0 <- data %>%
  filter(datetime >= as_date("2016-04-01") & 
           datetime <= as_date("2016-12-31")) %>%
  summarise(btc = sd(diff(btc)),
            eth = sd(diff(eth)),
            ltc = sd(diff(ltc)))

N <- nrow(data)
b1_grid <- seq(30, 360, by = 10) / (2 * N^(-1/3) / log(N) * N)
b2_grid <- seq(30, 360, by = 10) / (2 * N^(-1/4) * N)

b1_scores_btc <- sapply(b1_grid, function(x) 
  compute_CV1(data$btc, k = 0, b1 = x))
b1_opt_btc <- b1_grid[which.min(b1_scores_btc)] 
b2_scores_btc <- sapply(b2_grid, function(x) 
  compute_CV2_lpo(data$btc, k = 0, b1 = b1_opt_btc, b2 = x, sd0 = sd0$btc, p = 7))
b2_opt_btc <- b2_grid[which.min(b2_scores_btc)]

btc_var <- estimate_variance(data$btc, k = 0, 
                             b1 = b1_opt_btc, b2 = b2_opt_btc, sd0 = sd0$btc)

b1_scores_eth <- sapply(b1_grid, function(x) 
  compute_CV1(data$eth, k = 0, b1 = x))
b1_opt_eth <- b1_grid[which.min(b1_scores_eth)] 
b2_scores_eth <- sapply(b2_grid, function(x) 
  compute_CV2_lpo(data$eth, k = 0, b1 = b1_opt_eth, b2 = x, sd0 = sd0$eth, p = 7))
b2_opt_eth <- b2_grid[which.min(b2_scores_eth)]

eth_var <- estimate_variance(data$eth, k = 0, 
                             b1 = b1_opt_eth, b2 = b2_opt_eth, sd0 = sd0$eth)

b1_scores_ltc <- sapply(b1_grid, function(x) 
  compute_CV1(data$ltc, k = 0, b1 = x))
b1_opt_ltc <- b1_grid[which.min(b1_scores_ltc)] 
b2_scores_ltc <- sapply(b2_grid, function(x) 
  compute_CV2_lpo(data$ltc, k = 0, b1 = b1_opt_ltc, b2 = x, sd0 = sd0$ltc, p = 10))
b2_opt_ltc <- b2_grid[which.min(b2_scores_ltc)]

ltc_var <- estimate_variance(data$ltc, k = 0, 
                             b1 = b1_opt_ltc, b2 = b2_opt_ltc, sd0 = sd0$ltc)

pacf(btc_var$e_trunc^2, na.action = na.pass, lag.max = 30)
pacf(eth_var$e_trunc^2, na.action = na.pass, lag.max = 30)
pacf(ltc_var$e_trunc^2, na.action = na.pass, lag.max = 30)

rw_btc <- rollapply(c(NA, diff(data$btc)^2), 
                    width = round(2 * b2_opt_btc * N^(-1/4) * N), 
                    FUN = function(x) mean(x, na.rm = TRUE),
                    fill = NA, align = "center", partial = TRUE)

rw_eth <- rollapply(c(NA, diff(data$eth)^2), 
                    width = round(2 * b2_opt_eth * N^(-1/4) * N), 
                    FUN = function(x) mean(x, na.rm = TRUE),
                    fill = NA, align = "center", partial = TRUE)

rw_ltc <- rollapply(c(NA, diff(data$ltc)^2), 
                    width = round(2 * b2_opt_ltc * N^(-1/4) * N), 
                    FUN = function(x) mean(x, na.rm = TRUE),
                    fill = NA, align = "center", partial = TRUE)

data.frame(datetime = data$datetime,
           hlz_btc = sqrt(btc_var$var_hat),
           rw_btc = sqrt(rw_btc),
           hlz_eth = sqrt(eth_var$var_hat),
           rw_eth = sqrt(rw_eth),
           hlz_ltc = sqrt(ltc_var$var_hat),
           rw_ltc = sqrt(rw_ltc)) %>%
  write.csv(file = "results/volatility.csv", row.names = FALSE)