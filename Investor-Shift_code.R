library(quantmod)
library(moments)
library(rugarch)
library(xts)
library(ggplot2)
library(naniar)  # can delete
library(dplyr)  # fix British spelling issue
library(tidyr)              
library(zoo)
library(reshape2)
library(scales)  # for date formatting
library(lubridate)
library(lmtest)
library(car)
library(sandwich)
library(broom)
# sp500 data download
sp500 <- getSymbols("^GSPC", src = "yahoo", auto.assign = FALSE, from = "2010-01-01", to = "2025-04-17")

# gold price download
gold <- getSymbols("GC=F", src = "yahoo", auto.assign = FALSE, from = "2010-01-01", to = "2025-04-17")

# EDA: exploratory data analysis+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sp500_adj <- Ad(sp500)
gold_adj  <- Ad(gold)
# merge 
data_merge <- merge(sp500_adj, gold_adj, join = "inner")
colnames(data_merge) <- c("SP500", "Gold")
# summary of NA values
sapply(data_merge, function(x) sum(is.na(x)))
# Drop rows with any NA values
data_merge <- na.omit(data_merge)
sapply(data_merge, function(x) sum(is.na(x)))
data_df <- data.frame(date = index(data_merge), coredata(data_merge))
# basic statsitics
summary_stats <- data_df %>%
  summarise(across(c(SP500, Gold), list(
    mean = ~mean(.x, na.rm = TRUE),
    sd   = ~sd(.x, na.rm = TRUE),
    min  = ~min(.x, na.rm = TRUE),
    max  = ~max(.x, na.rm = TRUE),
    skew = ~moments::skewness(.x, na.rm = TRUE),
    kurt = ~moments::kurtosis(.x, na.rm = TRUE)
  )))

print(summary_stats)

# Transform to long format
data_long <- pivot_longer(data_df, cols = c(SP500, Gold), names_to = "Asset", values_to = "Price")

# Plot
ggplot(data_long, aes(x = date, y = Price, color = Asset)) +
  geom_line() +
  labs(
    title = "Daily Adjusted Prices: S&P 500 vs Gold",
    x = "Date", y = "Adjusted Close"
  ) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# returns ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# daily log returns
data_df <- data_df %>%
  mutate(
    SP500_ret = c(NA, diff(log(SP500))),
    Gold_ret  = c(NA, diff(log(Gold)))
  )

# plot histograms
data_ret_long <- pivot_longer(data_df, cols = c(SP500_ret, Gold_ret), names_to = "Asset", values_to = "LogReturn")

ggplot(data_ret_long, aes(x = LogReturn, fill = Asset)) +
  geom_histogram(bins = 100, alpha = 0.6, position = "identity") +
  facet_wrap(~Asset, scales = "free") +
  theme_minimal()

# Calculate summary statistics for returns
summary_stats <- data_df %>%
  summarise(
    mean_SP500  = mean(SP500_ret, na.rm = TRUE),
    sd_SP500    = sd(SP500_ret, na.rm = TRUE),
    min_SP500   = min(SP500_ret, na.rm = TRUE),
    max_SP500   = max(SP500_ret, na.rm = TRUE),
    skew_SP500  = skewness(SP500_ret, na.rm = TRUE),
    kurt_SP500  = kurtosis(SP500_ret, na.rm = TRUE),
    
    mean_Gold   = mean(Gold_ret, na.rm = TRUE),
    sd_Gold     = sd(Gold_ret, na.rm = TRUE),
    min_Gold    = min(Gold_ret, na.rm = TRUE),
    max_Gold    = max(Gold_ret, na.rm = TRUE),
    skew_Gold   = skewness(Gold_ret, na.rm = TRUE),
    kurt_Gold   = kurtosis(Gold_ret, na.rm = TRUE)
  )

# Print the results
print(summary_stats)

# drop na for returns
returns_clean <- na.omit(data_df[, c("SP500_ret", "Gold_ret")])

# entire correlation
cor(returns_clean)

# scatter plot
ggplot(returns_clean, aes(x = SP500_ret, y = Gold_ret)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Scatter Plot of Daily Returns: SP500 vs Gold", x = "S&P 500 Return", y = "Gold Return") +
  theme_minimal()



# Step 1: Extract year and month
monthly_counts <- data_df %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(trading_days = n()) %>%
  ungroup()

# Step 2: Calculate average
average_days <- mean(monthly_counts$trading_days)

# Print result
print(average_days)


# the realized correlation calculation++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# realized variance & covariance (21-day rolling sum of squares and cross-products)
data_df <- data_df %>%
  mutate(
    RV_sp500 = rollapply(SP500_ret^2, width = 21, FUN = sum, fill = NA, align = "right"),
    RV_gold  = rollapply(Gold_ret^2, width = 21, FUN = sum, fill = NA, align = "right"),
    RCOV     = rollapply(SP500_ret * Gold_ret, width = 21, FUN = sum, fill = NA, align = "right")
  ) %>%
  # Step 2: Realized correlation
  mutate(
    rho = RCOV / sqrt(RV_sp500 * RV_gold)
  )

# view last few rows of results
tail(data_df[, c("date", "RV_sp500", "RV_gold", "RCOV", "rho")], 10)

# plot to observe results
# reshape 
var_plot_df <- data_df %>%
  select(date, RV_sp500, RV_gold) %>%
  pivot_longer(cols = c(RV_sp500, RV_gold), names_to = "Asset", values_to = "Variance")

# Clean up asset names
var_plot_df$Asset <- recode(var_plot_df$Asset,
                            RV_sp500 = "S&P 500",
                            RV_gold  = "Gold")
# plot realized variance (volatility)

ggplot(var_plot_df, aes(x = date, y = Variance, color = Asset)) +
  geom_line(size = 1) +
  labs(
    title = "21-Day Rolling Realized Variance",
    x = "Date",
    y = "Variance",
    color = "Asset"
  ) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# plot realized correlation
ggplot(data_df, aes(x = date, y = rho)) +
  geom_line(color = "steelblue", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "21-Day Rolling Realized Correlation: S&P 500 vs Gold",
    x = "Date",
    y = "Realized Correlation (ρ)"
  ) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



# period for shock+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set up search window
dates_to_test <- seq(as.Date("2020-03-20"), as.Date("2020-06-20"), by = "5 days") #start date difference

# Create empty data frame to collect results
pval_grid <- expand.grid(start = dates_to_test, end = dates_to_test) %>%
  filter(end > start + 20)  

# Run regression for each window
# Run regression for each window and store both estimate and p-value
scan_results <- pval_grid %>%
  rowwise() %>%
  mutate(
    model_stats = {
      df <- data_df %>%
        mutate(
          dummy = ifelse(date >= start & date <= end, 1, 0),
          rho_lag = lag(rho)
        )
      reg_data <- na.omit(df[, c("rho", "rho_lag", "dummy")])
      if (nrow(reg_data) < 100) return(tibble(estimate = NA_real_, p_value = NA_real_))
      model <- lm(rho ~ rho_lag + dummy, data = reg_data)
      tidy(model) %>%
        filter(term == "dummy") %>%
        select(estimate, p.value)
    }
  ) %>%
  unnest(model_stats) %>%
  ungroup()


# Pivot for heatmap
heatmap_data <- scan_results %>%
  pivot_wider(names_from = end, values_from = p.value)

# Convert to matrix and format rownames/colnames
heatmap_mat <- as.matrix(heatmap_data[,-1])
rownames(heatmap_mat) <- as.character(scan_results$start[seq(1, nrow(heatmap_data), length.out = nrow(heatmap_data))])

# Plot heatmap
heatmap_long <- melt(heatmap_mat, varnames = c("start", "end"), value.name = "p_value")

ggplot(heatmap_long, aes(x = end, y = start, fill = p_value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = 1, na.value = "gray", limits = c(0, 0.3)) +
  labs(title = "Heatmap of Dummy Regression p-values",
       subtitle = "Each cell: p-value of dummy (start to end window)",
       x = "End Date", y = "Start Date", fill = "p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Show significant dummy windows
scan_results %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  arrange(p.value) %>%
  print(n = Inf)

# regression for Covid, a beginning regression tryout+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
covid_window_df <- data_df %>%
  filter(date >= as.Date("2019-06-01") & date <= as.Date("2021-12-31")) %>%
  mutate(
    covid_dummy = ifelse(date >= as.Date("2020-04-30") & date <= as.Date("2020-06-29"), 1, 0),
    rho_lag = lag(rho)
  )

model_covid <- lm(rho ~ rho_lag + covid_dummy, data = na.omit(covid_window_df))
summary(model_covid)


# Heteroskedasticity test (Breusch-Pagan)
bptest(model_covid)

# Autocorrelation test (Breusch-Godfrey, 1st and 2nd lag)
bgtest(model_covid, order = 1)
bgtest(model_covid, order = 2)


# regression for Covid, a beginning regression tryout+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
covid_window_df <- data_df %>%
  filter(date >= as.Date("2010-01-01") & date <= as.Date("2025-04-17")) %>%
  mutate(
    covid_dummy = ifelse(date >= as.Date("2020-04-30") & date <= as.Date("2020-06-29"), 1, 0),
    rho_lag = lag(rho)
  )

model_covid1 <- lm(rho ~ rho_lag + covid_dummy, data = na.omit(covid_window_df))
summary(model_covid1)
bptest(model_covid1)

# Autocorrelation test (Breusch-Godfrey, 1st and 2nd lag)
bgtest(model_covid1, order = 1)
bgtest(model_covid1, order = 2)

# solve the autocorrelation issue
covid_window_df <- data_df %>%
  mutate(
    covid_dummy = ifelse(date >= as.Date("2020-04-30") & date <= as.Date("2020-06-29"), 1, 0),
    rho_lag = lag(rho),
    RV_sp500_lag = lag(RV_sp500),
    RV_gold_lag = lag(RV_gold)
  )

model_covid_corrected <- lm(rho ~ rho_lag + RV_sp500_lag + RV_gold_lag + covid_dummy, data = na.omit(covid_window_df))
summary(model_covid_corrected)

bgtest(model_covid_corrected, order = 1)
bgtest(model_covid_corrected, order = 2) # the issue still there

coeftest(model_covid_corrected, vcov. = NeweyWest) # this prove the significance of dummy variable coefficient


# period for shock tariff+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set up search window
dates_to_test <- seq(as.Date("2025-01-01"), as.Date("2025-04-17"), by = "5 days") #start date difference

# Create empty data frame to collect results
pval_grid <- expand.grid(start = dates_to_test, end = dates_to_test) %>%
  filter(end > start + 20)  

# Run regression for each window
# Run regression for each window and store both estimate and p-value
scan_results <- pval_grid %>%
  rowwise() %>%
  mutate(
    model_stats = {
      df <- data_df %>%
        mutate(
          dummy = ifelse(date >= start & date <= end, 1, 0),
          rho_lag = lag(rho)
        )
      reg_data <- na.omit(df[, c("rho", "rho_lag", "dummy")])
      if (nrow(reg_data) < 100) return(tibble(estimate = NA_real_, p_value = NA_real_))
      model <- lm(rho ~ rho_lag + dummy, data = reg_data)
      tidy(model) %>%
        filter(term == "dummy") %>%
        select(estimate, p.value)
    }
  ) %>%
  unnest(model_stats) %>%
  ungroup()


# Pivot for heatmap
heatmap_data <- scan_results %>%
  pivot_wider(names_from = end, values_from = p.value)

# Convert to matrix and format rownames/colnames
heatmap_mat <- as.matrix(heatmap_data[,-1])
rownames(heatmap_mat) <- as.character(scan_results$start[seq(1, nrow(heatmap_data), length.out = nrow(heatmap_data))])

# Plot heatmap

heatmap_long <- melt(heatmap_mat, varnames = c("start", "end"), value.name = "p_value")

ggplot(heatmap_long, aes(x = end, y = start, fill = p_value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = 1, na.value = "gray", limits = c(0, 0.3)) +
  labs(title = "Heatmap of Dummy Regression p-values",
       subtitle = "Each cell: p-value of dummy (start to end window)",
       x = "End Date", y = "Start Date", fill = "p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Show significant dummy windows
scan_results %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  arrange(p.value) %>%
  print(n = Inf)

# regression for tariff+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
tariff_window_df <- data_df %>%
  filter(date >= as.Date("2024-06-01") & date <= as.Date("2025-04-17")) %>%
  mutate(
    tariff_dummy = ifelse(date >= as.Date("2025-02-20") & date <= as.Date("2025-03-28"), 1, 0),
    rho_lag = lag(rho)
  )

model_tariff <- lm(rho ~ rho_lag + tariff_dummy, data = na.omit(tariff_window_df))
summary(model_tariff)

bptest(model_tariff)

# Autocorrelation test (Breusch-Godfrey, 1st and 2nd lag)
bgtest(model_tariff, order = 1)
bgtest(model_tariff, order = 2)


# regression for tariff+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
tariff_window_df <- data_df %>%
  filter(date >= as.Date("2010-01-01") & date <= as.Date("2025-04-17")) %>%
  mutate(
    tariff_dummy = ifelse(date >= as.Date("2025-02-20") & date <= as.Date("2025-03-27"), 1, 0),
    rho_lag = lag(rho)
  )

model_tariff1 <- lm(rho ~ rho_lag + tariff_dummy, data = na.omit(tariff_window_df))
summary(model_tariff1)

bptest(model_tariff1)

# Autocorrelation test (Breusch-Godfrey, 1st and 2nd lag)
bgtest(model_tariff1, order = 1)
bgtest(model_tariff1, order = 2)

#solve for autocorrelation
# Add lagged variables
tariff_window_df <- data_df %>%
  mutate(
    tariff_dummy = ifelse(date >= as.Date("2025-02-20") & date <= as.Date("2025-03-27"), 1, 0),
    rho_lag = lag(rho),
    RV_sp500_lag = lag(RV_sp500),
    RV_gold_lag = lag(RV_gold)
  )

# New model with lagged regressors
model_tariff_corrected <- lm(rho ~ rho_lag + RV_sp500_lag + RV_gold_lag + tariff_dummy, 
                             data = na.omit(tariff_window_df))

summary(model_tariff_corrected)

# Test again for autocorrelation
bptest(model_tariff_corrected)
bgtest(model_tariff_corrected, order = 1)
bgtest(model_tariff_corrected, order = 2)

coeftest(model_tariff_corrected, vcov. = NeweyWest)


# period for non shock period+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set up search window
dates_to_test <- seq(as.Date("2019-01-01"), as.Date("2009-05-01"), by = "5 days") #start date difference

# Create empty data frame to collect results
pval_grid <- expand.grid(start = dates_to_test, end = dates_to_test) %>%
  filter(end > start + 20)  

# Run regression for each window
# Run regression for each window and store both estimate and p-value
scan_results <- pval_grid %>%
  rowwise() %>%
  mutate(
    model_stats = {
      df <- data_df %>%
        mutate(
          dummy = ifelse(date >= start & date <= end, 1, 0),
          rho_lag = lag(rho)
        )
      reg_data <- na.omit(df[, c("rho", "rho_lag", "dummy")])
      if (nrow(reg_data) < 100) return(tibble(estimate = NA_real_, p_value = NA_real_))
      model <- lm(rho ~ rho_lag + dummy, data = reg_data)
      tidy(model) %>%
        filter(term == "dummy") %>%
        select(estimate, p.value)
    }
  ) %>%
  unnest(model_stats) %>%
  ungroup()


# Pivot for heatmap
heatmap_data <- scan_results %>%
  pivot_wider(names_from = end, values_from = p.value)

# Convert to matrix and format rownames/colnames
heatmap_mat <- as.matrix(heatmap_data[,-1])
rownames(heatmap_mat) <- as.character(scan_results$start[seq(1, nrow(heatmap_data), length.out = nrow(heatmap_data))])

# Plot heatmap

heatmap_long <- melt(heatmap_mat, varnames = c("start", "end"), value.name = "p_value")

ggplot(heatmap_long, aes(x = end, y = start, fill = p_value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = 1, na.value = "gray", limits = c(0, 0.3)) +
  labs(title = "Heatmap of Dummy Regression p-values",
       subtitle = "Each cell: p-value of dummy (start to end window)",
       x = "End Date", y = "Start Date", fill = "p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Show significant dummy windows
scan_results %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  arrange(p.value) %>%
  print(n = Inf)


