# Cryptocurrency market speculative bubbles

Code, data and results for a bachelor thesis on cryptocurrency market bubbles using the Phillips-Shi-Yu (2015) approach and the Harvey-Leybourne-Zu (2023) nonparametric variance function estimator.

## Files

**`close_price.csv`** — daily closing prices of Bitcoin, Ethereum and Litecoin from 2016-01-01 to 2025-12-31, sourced from [CoinGecko](https://www.coingecko.com).

**`hlz_volatility.R`** — implements the nonparametric variance function estimator of Harvey et al. (2023) with a uniform kernel, including bandwidth selection via leave-one-out and leave-p-out cross-validation.

**`radf.R`** — uses the `exuber` package for recursive ADF-based explosiveness tests and implements Monte Carlo critical value estimation for the date-stamping procedure accounting for the multiple testing problem, and the `locate_episodes` function for date-stamping explosive periods.

## Results

Results are stored in the `results/` directory: `volatility.csv`, `bsadf.csv`, `radf_cv.rds`, `episodes.rds`.
