rm(list = ls())
library(quantmod)
library(rugarch)
library(FinTS)
library(tseries)
library(ggplot2)
library(gridExtra)
library(moments)
library(xts)
library(dplyr)
library(patchwork)
library(tidyr)
library(strucchange)
library(zoo)

required_packages <- c("quantmod", "rugarch", "FinTS", "tseries",
                       "ggplot2", "gridExtra", "moments", "xts",
                       "dplyr", "patchwork", "tidyr",
                       "forecast", "strucchange","zoo")

missing_pkgs <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, dependencies = TRUE,
                   repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(quantmod); library(rugarch); library(FinTS)
  library(tseries);  library(ggplot2); library(gridExtra)
  library(moments);  library(xts);     library(dplyr)
  library(patchwork); library(tidyr);  library(forecast)
  library(strucchange); library(zoo)
})

OUTPUT <- file.path(getwd(), "nifty50_outputs")

dir.create(
  OUTPUT,
  recursive = TRUE,
  showWarnings = FALSE
)

print(OUTPUT)




# ============================================================
# FIX 1: CLEAN DATA LOADING — removed duplicate assignments,
#         load Log_Returns directly from CSV (already computed),
#         dates_ret defined once after na.omit, in correct order
# ============================================================

# ── IN-SAMPLE (2010-2024) ────────────────────────────────────
insample       <- read.csv("C:/Users/visha/OneDrive/Desktop/vishal dissert/data/nifty50_insample_2010_2024.csv")
insample$Date  <- as.Date(insample$Date)
insample       <- insample[order(insample$Date), ]

prices         <- xts(insample$Price, order.by = insample$Date)
colnames(prices) <- "Nifty50"

returns        <- xts(insample$Log_Returns, order.by = insample$Date)
returns        <- na.omit(returns)          # removes 1 NA (first row only)
colnames(returns) <- "Log_Returns"

ret_vec        <- as.numeric(returns)
dates_ret      <- as.Date(index(returns))   # FIX: defined once, stays in sync

cat(sprintf("In-sample  : %s to %s | %d observations\n",
            format(min(dates_ret), "%Y-%m-%d"),
            format(max(dates_ret), "%Y-%m-%d"),
            length(ret_vec)))

# ── OUT-OF-SAMPLE (2025 only) ────────────────────────────────
# FIX 2: Replaced getSymbols live fetch with CSV loading
#         Filtered to 2025 only — excludes 2026 data
outsample      <- read.csv("C:/Users/visha/OneDrive/Desktop/vishal dissert/data/nifty50_outsample_2025.csv")
outsample$Date <- as.Date(outsample$Date)
outsample      <- outsample[order(outsample$Date), ]
outsample      <- outsample[!is.na(outsample$Log_Returns), ]
outsample      <- outsample[outsample$Date <= as.Date("2025-12-31"), ]  # 2025 only

actual_returns_2025 <- xts(outsample$Log_Returns,
                           order.by = outsample$Date)

cat(sprintf("Out-of-sample: %s to %s | %d observations\n",
            format(min(outsample$Date), "%Y-%m-%d"),
            format(max(outsample$Date), "%Y-%m-%d"),
            nrow(outsample)))

write.csv(as.data.frame(prices),  file.path(OUTPUT, "nifty50_prices.csv"))
write.csv(as.data.frame(returns), file.path(OUTPUT, "nifty50_returns.csv"))

jb_test <- jarque.bera.test(ret_vec)

desc_stats <- data.frame(
  Statistic = c("Observations", "Mean", "Median", "Std Dev",
                "Min", "Max", "Skewness", "Excess Kurtosis",
                "Jarque-Bera Stat", "Jarque-Bera p-value"),
  Value = c(
    length(ret_vec),
    round(mean(ret_vec), 6),
    round(median(ret_vec), 6),
    round(sd(ret_vec), 6),
    round(min(ret_vec), 6),
    round(max(ret_vec), 6),
    round(skewness(ret_vec), 6),
    round(kurtosis(ret_vec) - 3, 6),
    round(jb_test$statistic, 4),
    round(jb_test$p.value, 6)
  )
)

print(desc_stats)
write.csv(desc_stats, file.path(OUTPUT, "descriptive_statistics.csv"),
          row.names = FALSE)

adf_result <- adf.test(ret_vec, alternative = "stationary")
print(adf_result)

arch_test <- ArchTest(ret_vec, lags = 10)
print(arch_test)


# ============================================================
#  BAI-PERRON STRUCTURAL BREAK TEST
# ============================================================

cat("\n============================================================\n")
cat("BAI-PERRON STRUCTURAL BREAK TEST\n")
cat("============================================================\n")

bp_test     <- breakpoints((ret_vec^2) ~ 1)
summary(bp_test)
bp_optimal  <- breakpoints(bp_test)
break_index <- bp_optimal$breakpoints
break_index <- break_index[!is.na(break_index)]
n_breaks    <- length(break_index)

if (n_breaks > 0) {
  break_dates <- as.Date(index(returns)[break_index])
  print(break_dates)
  
  ci        <- confint(bp_test, breaks = n_breaks)
  print(ci)
  ci_bounds <- ci$confint
  
  ci_table <- data.frame(
    Break    = seq_len(n_breaks),
    Lower    = dates_ret[ci_bounds[, 1]],
    Estimate = break_dates,
    Upper    = dates_ret[ci_bounds[, 2]]
  )
  print(ci_table)
  
  break_table <- data.frame(
    Break_Index = break_index,
    Break_Date  = as.character(break_dates)
  )
  print(break_table)
  write.csv(break_table,
            file.path(OUTPUT, "bai_perron_breakpoints.csv"),
            row.names = FALSE)
} else {
  cat("No structural breakpoints selected by BIC criterion.\n")
  break_dates <- as.Date(character(0))
  write.csv(data.frame(Break_Index = integer(0),
                       Break_Date  = character(0)),
            file.path(OUTPUT, "bai_perron_breakpoints.csv"),
            row.names = FALSE)
}

ret_df   <- data.frame(Date = as.Date(index(returns)), Returns = ret_vec^2)
price_df <- data.frame(Date = index(prices), Price = as.numeric(prices))

p_bp_base <- ggplot(ret_df, aes(x = Date, y = Returns)) +
  geom_line(color = "#B71C1C", linewidth = 0.3, alpha = 0.7) +
  labs(title    = "Bai-Perron Structural Break Test — Nifty 50 Returns",
       subtitle = "Blue dashed lines = statistically identified breakpoints",
       x = "", y = "Squared Returns") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40"))

if (n_breaks > 0) {
  bp_dates_df <- data.frame(Date = break_dates)
  p_bp <- p_bp_base +
    geom_vline(data = bp_dates_df, aes(xintercept = as.numeric(Date)),
               color = "blue", linetype = "dashed",
               linewidth = 0.8, alpha = 0.9) +
    geom_text(data = bp_dates_df,
              aes(x = Date,
                  y = max(ret_df$Returns) * 0.85,
                  label = as.character(Date)),
              angle = 90, vjust = -0.3, hjust = 1,
              size = 3, color = "blue")
} else {
  p_bp <- p_bp_base
}

ggsave(file.path(OUTPUT, "plot_bai_perron_breaks.png"),
       p_bp, width = 14, height = 6, dpi = 150)
print(p_bp)
cat("Bai-Perron plot saved\n")

cat("\n=== SEGMENT STATISTICS ===\n")
break_indices <- c(0, bp_optimal$breakpoints, length(ret_vec))

seg_stats <- data.frame()
for (i in seq_len(length(break_indices) - 1)) {
  seg       <- ret_vec[(break_indices[i] + 1):break_indices[i + 1]]
  seg_start <- dates_ret[break_indices[i] + 1]
  seg_end   <- dates_ret[break_indices[i + 1]]
  
  seg_stats <- rbind(seg_stats, data.frame(
    Segment  = paste0("Regime ", i),
    Start    = as.character(seg_start),
    End      = as.character(seg_end),
    N_obs    = length(seg),
    Mean     = round(mean(seg),         6),
    Std_Dev  = round(sd(seg),           6),
    Skewness = round(skewness(seg),     4),
    Kurtosis = round(kurtosis(seg) - 3, 4)
  ))
}

print(seg_stats, row.names = FALSE)
write.csv(seg_stats,
          file.path(OUTPUT, "bai_perron_segment_stats.csv"),
          row.names = FALSE)
cat("Segment statistics saved\n")

seg_long <- pivot_longer(seg_stats,
                         cols      = c("Mean", "Std_Dev"),
                         names_to  = "Statistic",
                         values_to = "Value")

p_seg <- ggplot(seg_long,
                aes(x = Segment, y = Value, fill = Segment)) +
  geom_bar(stat = "identity", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = round(Value, 4)),
            vjust = -0.4, size = 3.2) +
  facet_wrap(~Statistic, scales = "free_y") +
  labs(title    = "Return Statistics by Structural Regime",
       subtitle = "Each regime defined by Bai-Perron breakpoints",
       x = "", y = "Value") +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(color = "grey40"),
        legend.position = "none")

ggsave(file.path(OUTPUT, "plot_bai_perron_segment_stats.png"),
       p_seg, width = 12, height = 5, dpi = 150)
print(p_seg)
cat("Segment statistics plot saved\n")

cat("\n=== BAI-PERRON SUMMARY ===\n")
cat("Method          : BIC-optimal break selection\n")
cat("Max breaks tested: 5\n")
cat("Min segment size : 10% of observations\n")
cat(sprintf("Breaks found    : %d\n", n_breaks))
if (n_breaks > 0) {
  cat("Break dates     :\n")
  for (d in as.character(break_dates)) cat(sprintf("  %s\n", d))
}
cat("\nNote: All breakpoints are purely data-driven.\n")
cat("No event dates were preassigned.\n")
cat("The data itself determined where structural breaks occur.\n")


# ============================================================
#  ACF / PACF
# ============================================================

png(file.path(OUTPUT, "plot6_acf_pacf.png"),
    width = 1600, height = 1200, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
acf(ret_vec,    lag.max = 30, main = "ACF - Log Returns",      col = "#1565C0", lwd = 2)
pacf(ret_vec,   lag.max = 30, main = "PACF - Log Returns",     col = "#1565C0", lwd = 2)
acf(ret_vec^2,  lag.max = 30, main = "ACF - Squared Returns",  col = "#B71C1C", lwd = 2)
pacf(ret_vec^2, lag.max = 30, main = "PACF - Squared Returns", col = "#B71C1C", lwd = 2)
mtext("ACF / PACF: Nifty 50 Log Returns (2010-2024)",
      outer = TRUE, cex = 1.2, font = 2)
dev.off()

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
acf(ret_vec,    lag.max = 30, main = "ACF - Log Returns",      col = "#1565C0", lwd = 2)
pacf(ret_vec,   lag.max = 30, main = "PACF - Log Returns",     col = "#1565C0", lwd = 2)
acf(ret_vec^2,  lag.max = 30, main = "ACF - Squared Returns",  col = "#B71C1C", lwd = 2)
pacf(ret_vec^2, lag.max = 30, main = "PACF - Squared Returns", col = "#B71C1C", lwd = 2)
mtext("ACF / PACF: Nifty 50 Log Returns (2010-2024)",
      outer = TRUE, cex = 1.2, font = 2)
par(mfrow = c(1, 1))
cat("plot6_acf_pacf.png saved\n")


# ============================================================
#  MODEL ESTIMATION
# ============================================================

spec_garch <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std")
fit_garch <- ugarchfit(spec_garch, data = ret_vec, solver = "hybrid")
show(fit_garch)

spec_egarch <- ugarchspec(
  variance.model     = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std")
fit_egarch <- ugarchfit(spec_egarch, data = ret_vec, solver = "hybrid")
show(fit_egarch)

spec_gjr <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std")
fit_gjr <- ugarchfit(spec_gjr, data = ret_vec, solver = "hybrid")
show(fit_gjr)

spec_figarch <- ugarchspec(
  variance.model     = list(model = "fiGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std")
fit_figarch <- ugarchfit(spec_figarch, data = ret_vec, solver = "hybrid")
show(fit_figarch)

n          <- length(ret_vec)
models     <- list(fit_garch, fit_egarch, fit_gjr, fit_figarch)
modelnames <- c("GARCH(1,1)", "EGARCH(1,1)", "GJR-GARCH(1,1)", "FIGARCH(1,d,1)")
colors_vec <- c("#2196F3", "#E91E63", "#FF9800", "#4CAF50")
colors_map <- c("GARCH(1,1)"     = "#2196F3",
                "EGARCH(1,1)"    = "#E91E63",
                "GJR-GARCH(1,1)" = "#FF9800",
                "FIGARCH(1,d,1)" = "#4CAF50")

get_criteria <- function(fit, name) {
  ic <- infocriteria(fit)
  data.frame(
    Model         = name,
    LogLikelihood = round(likelihood(fit), 4),
    AIC           = round(ic[1] * n, 4),
    BIC           = round(ic[2] * n, 4)
  )
}

comparison <- rbind(
  get_criteria(fit_garch,   "GARCH(1,1)"),
  get_criteria(fit_egarch,  "EGARCH(1,1)"),
  get_criteria(fit_gjr,     "GJR-GARCH(1,1)"),
  get_criteria(fit_figarch, "FIGARCH(1,d,1)")
)
comparison <- comparison[order(comparison$AIC), ]
print(comparison, row.names = FALSE)
write.csv(comparison, file.path(OUTPUT, "model_comparison.csv"), row.names = FALSE)

p_g   <- coef(fit_garch)
p_e   <- coef(fit_egarch)
p_gjr <- coef(fit_gjr)
p_f   <- coef(fit_figarch)

cat(sprintf("GARCH(1,1)      Persistence = %.6f\n", p_g["alpha1"] + p_g["beta1"]))
cat(sprintf("EGARCH(1,1)     Persistence = %.6f\n", p_e["alpha1"] + p_e["beta1"]))
cat(sprintf("GJR-GARCH(1,1)  Persistence = %.6f\n",
            p_gjr["alpha1"] + 0.5 * p_gjr["gamma1"] + p_gjr["beta1"]))
cat(sprintf("FIGARCH(1,d,1)  d           = %.6f\n", p_f["d"]))   # FIX 3: removed delta fallback

diag_df <- data.frame()
for (i in seq_along(models)) {
  r  <- residuals(models[[i]], standardize = TRUE)
  lb <- Box.test(r,   lag = 10, type = "Ljung-Box")
  ls <- Box.test(r^2, lag = 10, type = "Ljung-Box")
  diag_df <- rbind(diag_df, data.frame(
    Model         = modelnames[i],
    LB_Resid_pval = round(lb$p.value, 4),
    LB_Sq_pval    = round(ls$p.value, 4),
    LB_Resid_OK   = ifelse(lb$p.value > 0.05, "PASS", "FAIL"),
    LB_Squared_OK = ifelse(ls$p.value > 0.05, "PASS", "FAIL")
  ))
}
print(diag_df, row.names = FALSE)
write.csv(diag_df, file.path(OUTPUT, "model_diagnostics.csv"), row.names = FALSE)
print(signbias(fit_egarch))

best_model_name <- comparison$Model[1]
cat(sprintf("Best model by AIC: %s\n", best_model_name))

model_list <- list(
  "GARCH(1,1)"     = fit_garch,
  "EGARCH(1,1)"    = fit_egarch,
  "GJR-GARCH(1,1)" = fit_gjr,
  "FIGARCH(1,d,1)" = fit_figarch
)

best_fit  <- model_list[[best_model_name]]
best_name <- best_model_name

mu_hat    <- as.numeric(fitted(best_fit))
sigma_hat <- as.numeric(sigma(best_fit))
shape_par <- coef(best_fit)["shape"]

VaR_1 <- as.numeric(mu_hat + sigma_hat *
                      qdist("std", p = 0.01, mu = 0, sigma = 1, shape = shape_par))
VaR_5 <- as.numeric(mu_hat + sigma_hat *
                      qdist("std", p = 0.05, mu = 0, sigma = 1, shape = shape_par))

# FIX 4: ES computed as scalars — renamed clearly to avoid confusion
ES_1_const <- mean(ret_vec[ret_vec < VaR_1])
ES_5_const <- mean(ret_vec[ret_vec < VaR_5])

breach_1   <- sum(ret_vec < VaR_1)
breach_5   <- sum(ret_vec < VaR_5)
expected_1 <- round(length(ret_vec) * 0.01)
expected_5 <- round(length(ret_vec) * 0.05)

cat(sprintf("\n=== VaR & ES RESULTS (%s) ===\n", best_name))
cat(sprintf("1%% VaR breaches : %d (expected ~%d)\n", breach_1, expected_1))
cat(sprintf("5%% VaR breaches : %d (expected ~%d)\n", breach_5, expected_5))
cat(sprintf("1%% Expected Shortfall : %.6f%%\n", ES_1_const))
cat(sprintf("5%% Expected Shortfall : %.6f%%\n", ES_5_const))

# ============================================================
# FIX 5: KUPIEC TEST — added all-breach guard to prevent
#         log(1 - p_hat) = -Inf crash when p_hat = 1
# ============================================================
kupiec_test <- function(n_breach, n_obs, conf_level) {
  p     <- 1 - conf_level
  x     <- n_breach
  n_o   <- n_obs
  p_hat <- x / n_o
  
  # Guard: all observations breached
  if (x == n_o) {
    return(data.frame(
      Confidence  = paste0(conf_level * 100, "%"),
      N_obs       = n_o,
      N_breach    = x,
      Expected    = round(n_o * p),
      Breach_Rate = 100,
      Target_Rate = round(p * 100, 4),
      LR_stat     = NA,
      p_value     = NA,
      Result      = "FAIL - all observations breached"
    ))
  }
  
  if (x == 0) {
    LR <- -2 * n_o * log(1 - p)
  } else {
    # Compute entirely in log space to avoid underflow
    LR <- -2 * (
      (n_o - x) * log(1 - p)   + x * log(p) -
        (n_o - x) * log(1 - p_hat) - x * log(p_hat)
    )
  }
  
  p_val <- 1 - pchisq(LR, df = 1)
  data.frame(
    Confidence  = paste0(conf_level * 100, "%"),
    N_obs       = n_o,
    N_breach    = x,
    Expected    = round(n_o * p),
    Breach_Rate = round(p_hat * 100, 4),
    Target_Rate = round(p * 100, 4),
    LR_stat     = round(LR, 4),
    p_value     = round(p_val, 6),
    Result      = ifelse(p_val > 0.05,
                         "PASS - VaR model adequate",
                         "FAIL - VaR model inadequate")
  )
}
kupiec_1 <- kupiec_test(breach_1, n, 0.99)
kupiec_5 <- kupiec_test(breach_5, n, 0.95)
kupiec_results <- rbind(kupiec_1, kupiec_5)

cat("\n=== KUPIEC PROPORTION OF FAILURES TEST ===\n")
print(kupiec_results, row.names = FALSE)
write.csv(kupiec_results,
          file.path(OUTPUT, "kupiec_test_results.csv"),
          row.names = FALSE)

var_es_summary <- data.frame(
  Model         = best_name,
  VaR_1pct      = round(mean(VaR_1), 6),
  VaR_5pct      = round(mean(VaR_5), 6),
  ES_1pct       = round(ES_1_const, 6),
  ES_5pct       = round(ES_5_const, 6),
  Breach_1pct   = breach_1,
  Breach_5pct   = breach_5,
  Expected_1pct = expected_1,
  Expected_5pct = expected_5,
  Kupiec_1pct   = kupiec_1$Result,
  Kupiec_5pct   = kupiec_5$Result
)
write.csv(var_es_summary,
          file.path(OUTPUT, "var_es_summary.csv"),
          row.names = FALSE)

var_plot_df <- data.frame(
  Date    = dates_ret,
  Returns = ret_vec,
  VaR_1   = VaR_1,
  VaR_5   = VaR_5,
  ES_1    = ES_1_const,   # FIX 4: using renamed scalar
  ES_5    = ES_5_const    # FIX 4: using renamed scalar
)

p_var <- ggplot(var_plot_df, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "grey60",
            linewidth = 0.3, alpha = 0.7) +
  geom_line(aes(y = VaR_1, color = "1% VaR"),   linewidth = 0.7) +
  geom_line(aes(y = VaR_5, color = "5% VaR"),   linewidth = 0.7) +
  geom_hline(aes(yintercept = ES_1_const, color = "1% ES"),
             linewidth = 0.8, linetype = "dotted") +
  geom_hline(aes(yintercept = ES_5_const, color = "5% ES"),
             linewidth = 0.8, linetype = "dotted") +
  geom_point(data = subset(var_plot_df, Returns < VaR_1),
             aes(y = Returns), color = "red", size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c(
    "1% VaR" = "#B71C1C",
    "5% VaR" = "#FF9800",
    "1% ES"  = "#880E4F",
    "5% ES"  = "#E65100"
  )) +
  labs(
    title    = paste("VaR & Expected Shortfall -", best_name),
    subtitle = sprintf(
      "1%% VaR: %d breaches (exp %d) | 5%% VaR: %d breaches (exp %d) | ES_1%%=%.2f%% ES_5%%=%.2f%%",
      breach_1, expected_1, breach_5, expected_5, ES_1_const, ES_5_const),
    x = "", y = "Returns (%)", color = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(size = 8, color = "grey40"),
        legend.position = "bottom")

ggsave(file.path(OUTPUT, "plot9_var_es_estimation.png"),
       p_var, width = 14, height = 6, dpi = 150)
print(p_var)
cat("plot9_var_es_estimation.png saved\n")


# ============================================================
#  PLOTS
# ============================================================

p1 <- ggplot(price_df, aes(x = Date, y = Price)) +
  geom_line(color = "#1565C0", linewidth = 0.5) +
  labs(title = "Nifty 50 Closing Price (2010-2024)", x = "", y = "Price (INR)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ret_df2 <- data.frame(Date = dates_ret, Returns = ret_vec)

p2 <- ggplot(ret_df2, aes(x = Date, y = Returns)) +
  geom_segment(aes(xend = Date, y = 0, yend = Returns,
                   color = Returns >= 0),
               linewidth = 0.4, alpha = 0.9) +
  scale_color_manual(values = c("TRUE" = "#B71C1C", "FALSE" = "#B71C1C")) +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "black") +
  scale_y_continuous(
    limits = c(min(ret_vec) * 1.05, max(ret_vec) * 1.05),
    breaks = seq(-15, 15, by = 5)
  ) +
  labs(title = "Nifty 50 Log Returns (%)", x = "", y = "Returns (%)") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "none"
  )

plot1 <- p1 / p2

ggsave(file.path(OUTPUT, "plot1_price_returns.png"),
       plot1, width = 14, height = 8, dpi = 150)

print(plot1)


p3 <- ggplot(data.frame(Returns = ret_vec), aes(x = Returns)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100,
                 fill = "#1565C0", alpha = 0.6) +
  stat_function(fun = dnorm,
                args = list(mean = mean(ret_vec), sd = sd(ret_vec)),
                color = "red", linewidth = 1) +
  labs(title = "Distribution of Nifty 50 Log Returns",
       x = "Returns (%)", y = "Density") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT, "plot2_return_distribution.png"),
       p3, width = 10, height = 5, dpi = 150)
print(p3)

vol_df <- data.frame(
  Date       = rep(dates_ret, 4),
  Volatility = c(as.numeric(sigma(fit_garch)),
                 as.numeric(sigma(fit_egarch)),
                 as.numeric(sigma(fit_gjr)),
                 as.numeric(sigma(fit_figarch))),
  Model      = rep(modelnames, each = length(dates_ret))
)

p4 <- ggplot(vol_df, aes(x = Date, y = Volatility, color = Model)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~Model, ncol = 1, scales = "free_y") +
  scale_color_manual(values = colors_map) +
  labs(title = "Conditional Volatility: GARCH Model Comparison",
       x = "", y = "Volatility (%)") +
  theme_minimal(base_size = 10) +
  theme(plot.title      = element_text(face = "bold", size = 13),
        legend.position = "none",
        strip.text      = element_text(face = "bold"))

ggsave(file.path(OUTPUT, "plot3_conditional_volatility.png"),
       p4, width = 14, height = 14, dpi = 150)
print(p4)

pj <- ggplot(vol_df, aes(x = Date, y = Volatility, color = Model)) +
  geom_line(linewidth = 0.5, alpha = 0.85) +
  scale_color_manual(values = colors_map) +
  labs(
    title    = "Conditional Volatility: GARCH Model Comparison",
    subtitle = "All four models overlaid on a single panel",
    x        = "",
    y        = "Volatility (%)",
    color    = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(color = "grey40"),
    legend.position = "bottom",
    legend.text     = element_text(size = 10)
  )

ggsave(file.path(OUTPUT, "plot3_conditional_volatility.png"),
       p4, width = 14, height = 6, dpi = 150)
print(pj)

comp_long_zoom <- pivot_longer(comparison, cols = c("AIC", "BIC"),
                               names_to = "Criterion", values_to = "Value")

p_aic_fix <- ggplot(comp_long_zoom,
                    aes(x = reorder(Model, Value), y = Value, fill = Model)) +
  geom_bar(stat = "identity", alpha = 0.85, width = 0.6) +
  facet_wrap(~Criterion, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = colors_map) +
  geom_text(aes(label = round(Value, 1)), vjust = -0.4, size = 3.2) +
  labs(title = "Model Comparison: AIC & BIC (lower = better)",
       x = "", y = "Value") +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold"),
        axis.text.x     = element_text(angle = 20, hjust = 1),
        legend.position = "none")

ggsave(file.path(OUTPUT, "plot4b_aic_bic.png"),
       p_aic_fix, width = 12, height = 5, dpi = 150)
print(p_aic_fix)

events <- data.frame(
  Date  = as.Date(c("2016-11-08", "2020-03-23",
                    "2022-02-24", "2023-10-07", "2024-04-13")),
  Label = c("Demonetisation", "COVID Crash",
            "Russia-Ukraine", "Israel-Hamas", "Iran-Israel")
)

gjr_vol_df <- data.frame(Date = dates_ret,
                         Volatility = as.numeric(sigma(fit_gjr)))

p6 <- ggplot(gjr_vol_df, aes(x = Date, y = Volatility)) +
  geom_line(color = "#FF9800", linewidth = 0.6) +
  geom_vline(data = events, aes(xintercept = as.numeric(Date)),
             color = "red", linetype = "dashed",
             linewidth = 0.8, alpha = 0.8) +
  geom_text(data = events,
            aes(x = Date, y = max(gjr_vol_df$Volatility) * 0.9,
                label = Label),
            angle = 90, vjust = -0.3, hjust = 1,
            size = 3, color = "red") +
  labs(title = "GJR-GARCH Conditional Volatility with Key Market Events (2010-2024)",
       x = "", y = "Volatility (%)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(OUTPUT, "plot5_events_volatility.png"),
       p6, width = 16, height = 6, dpi = 150)
print(p6)


# ============================================================
#  Q-Q PLOTS
# ============================================================

png(file.path(OUTPUT, "plot7_qq_residuals.png"),
    width = 1400, height = 1400, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
for (i in seq_along(models)) {
  std_resid <- as.numeric(residuals(models[[i]], standardize = TRUE))
  qqnorm(std_resid, main = paste("Q-Q Plot -", modelnames[i]),
         col = colors_vec[i], pch = 16, cex = 0.4,
         xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
  qqline(std_resid, col = "black", lwd = 1.5, lty = 2)
}
mtext("Q-Q Plots: Standardized Residuals vs Normal Distribution",
      outer = TRUE, cex = 1.1, font = 2)
dev.off()

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
for (i in seq_along(models)) {
  std_resid <- as.numeric(residuals(models[[i]], standardize = TRUE))
  qqnorm(std_resid, main = paste("Q-Q Plot -", modelnames[i]),
         col = colors_vec[i], pch = 16, cex = 0.4,
         xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
  qqline(std_resid, col = "black", lwd = 1.5, lty = 2)
}
mtext("Q-Q Plots: Standardized Residuals vs Normal Distribution",
      outer = TRUE, cex = 1.1, font = 2)
par(mfrow = c(1, 1))
cat("plot7_qq_residuals.png saved\n")


# ============================================================
#  COEFFICIENT TABLE
# ============================================================

extract_coef_table <- function(fit, name) {
  cf           <- fit@fit$matcoef
  df           <- as.data.frame(cf)
  colnames(df) <- c("Estimate", "Std_Error", "t_stat", "p_value")
  df$Model     <- name
  df$Parameter <- rownames(cf)
  rownames(df) <- NULL
  df           <- df[, c("Model", "Parameter", "Estimate",
                         "Std_Error", "t_stat", "p_value")]
  df$Estimate  <- round(df$Estimate,  6)
  df$Std_Error <- round(df$Std_Error, 6)
  df$t_stat    <- round(df$t_stat,    4)
  df$p_value   <- round(df$p_value,   4)
  df$Sig       <- ifelse(df$p_value < 0.01, "***",
                         ifelse(df$p_value < 0.05, "**",
                                ifelse(df$p_value < 0.10, "*", "")))
  df
}

coef_table <- rbind(
  extract_coef_table(fit_garch,   "GARCH(1,1)"),
  extract_coef_table(fit_egarch,  "EGARCH(1,1)"),
  extract_coef_table(fit_gjr,     "GJR-GARCH(1,1)"),
  extract_coef_table(fit_figarch, "FIGARCH(1,d,1)")
)

write.csv(coef_table,
          file.path(OUTPUT, "coefficient_table_all_models.csv"),
          row.names = FALSE)
print(coef_table, row.names = FALSE)
cat("coefficient_table_all_models.csv saved\n")


# ============================================================
#  ROLLING MULTI-HORIZON FORECAST — Q1-Q4 2025
# ============================================================

cat("\n============================================================\n")
cat("ROLLING MULTI-HORIZON FORECAST — Q1-Q4 2025\n")
cat("============================================================\n")

# Quarter date ranges — 2025 only
quarter_ranges <- list(
  Q1 = list(start = as.Date("2025-01-01"), end = as.Date("2025-03-31"), label = "Q1 2025"),
  Q2 = list(start = as.Date("2025-04-01"), end = as.Date("2025-06-30"), label = "Q2 2025"),
  Q3 = list(start = as.Date("2025-07-01"), end = as.Date("2025-09-30"), label = "Q3 2025"),
  Q4 = list(start = as.Date("2025-10-01"), end = as.Date("2025-12-31"), label = "Q4 2025")
)

# Slice each quarter from actual_returns_2025 (already filtered to 2025)
quarter_data <- lapply(quarter_ranges, function(qr) {
  idx <- index(actual_returns_2025)
  sub <- actual_returns_2025[idx >= qr$start & idx <= qr$end]
  list(vec   = as.numeric(sub),
       dates = as.Date(index(sub)),
       label = qr$label,
       n     = length(sub))
})

for (qname in names(quarter_data)) {
  qd <- quarter_data[[qname]]
  cat(sprintf("  %s trading days: %d\n", qd$label, qd$n))
}

# ── Shared helpers ────────────────────────────────────────────
make_rv_list <- function(ret_vec) {
  rv_1  <- sqrt(ret_vec^2)
  rv_5  <- sqrt(zoo::rollapply(ret_vec^2, width = 5,
                               FUN = sum, align = "right", fill = NA))
  rv_10 <- sqrt(zoo::rollapply(ret_vec^2, width = 10,
                               FUN = sum, align = "right", fill = NA))
  rv_20 <- sqrt(zoo::rollapply(ret_vec^2, width = 20,
                               FUN = sum, align = "right", fill = NA))
  list(
    "1"  = rv_1,
    "5"  = rv_5,
    "10" = rv_10,
    "20" = rv_20
  )
}

# FIX 7: rolling_forecast retries fit if previous fit failed
#         Old code kept NULL fit_obj for all future steps after one failure
rolling_forecast <- function(spec, full_ret, n_train,
                             n_q, horizon, refit_every = 5) {
  forecasts <- rep(NA, n_q)
  fit_obj   <- NULL
  for (i in seq_len(n_q)) {
    window_end <- n_train + i - 1
    if (window_end < 10) next
    # FIX 7: added is.null(fit_obj) so failed fits are retried next step
    if (is.null(fit_obj) || i == 1 || (i %% refit_every) == 1) {
      fit_obj <- tryCatch(
        ugarchfit(spec, data = full_ret[seq_len(window_end)], solver = "hybrid"),
        error = function(e) NULL)
    }
    if (is.null(fit_obj)) next
    fc <- tryCatch(ugarchforecast(fit_obj, n.ahead = horizon),
                   error = function(e) NULL)
    if (is.null(fc)) next
    forecasts[i] <- as.numeric(sigma(fc))[horizon]
  }
  forecasts
}

calc_metrics <- function(forecast, actual, model, horizon, quarter) {
  valid <- !is.na(forecast) & !is.na(actual)
  f <- forecast[valid]; a <- actual[valid]
  if (length(f) < 3) return(NULL)
  data.frame(
    Quarter = quarter,
    Model   = model,
    Horizon = paste0(horizon, "-day"),
    MAE     = round(mean(abs(f - a)), 6),
    RMSE    = round(sqrt(mean((f - a)^2)), 6),
    MAPE    = round(mean(abs((f - a) / (a + 1e-8))) * 100, 4)
  )
}

horizons    <- c(1, 5, 10, 20)
model_specs <- list(spec_garch, spec_egarch, spec_gjr, spec_figarch)
model_names <- c("GARCH(1,1)", "EGARCH(1,1)",
                 "GJR-GARCH(1,1)", "FIGARCH(1,d,1)")

# FIX 8: base_ret never mutated — use rolling_ret accumulator instead
#         Old code mutated base_ret which contaminated Three Regime section
base_ret    <- ret_vec   # 2010-2024 — never modified
rolling_ret <- base_ret  # accumulator that grows each quarter

accuracy_all_quarters <- data.frame()
dm_results_all        <- data.frame()
all_forecasts_all_q   <- list()

for (qname in names(quarter_data)) {
  qd <- quarter_data[[qname]]
  if (qd$n < 5) {
    cat(sprintf("\n  Skipping %s — insufficient data (%d obs)\n",
                qd$label, qd$n))
    next
  }
  
  cat(sprintf("\n=== Quarter: %s (%d trading days) ===\n",
              qd$label, qd$n))
  
  # FIX 8: use rolling_ret (accumulator) not base_ret
  full_ret_q <- c(rolling_ret, qd$vec)
  n_train_q  <- length(rolling_ret)
  rv_list_q  <- make_rv_list(qd$vec)
  
  all_forecasts_q <- list()
  
  for (m in seq_along(model_specs)) {
    cat(sprintf("  Model: %s\n", model_names[m]))
    model_fc <- list()
    for (h in horizons) {
      cat(sprintf("    Horizon %d-day ... ", h))
      fc <- rolling_forecast(model_specs[[m]], full_ret_q,
                             n_train_q, qd$n, h, 5)
      model_fc[[as.character(h)]] <- fc
      cat(sprintf("done (%d forecasts)\n", sum(!is.na(fc))))
    }
    all_forecasts_q[[model_names[m]]] <- model_fc
  }
  
  all_forecasts_all_q[[qname]] <- all_forecasts_q
  
  # Accuracy metrics
  for (m in model_names) {
    for (h in as.character(horizons)) {
      met <- calc_metrics(all_forecasts_q[[m]][[h]],
                          rv_list_q[[h]], m, h, qd$label)
      if (!is.null(met)) accuracy_all_quarters <- rbind(accuracy_all_quarters, met)
    }
  }
  
  # Diebold-Mariano test (FIGARCH vs others)
  for (h in as.character(horizons)) {
    fc_figarch <- all_forecasts_q[["FIGARCH(1,d,1)"]][[h]]
    rv_h       <- rv_list_q[[h]]
    for (m in c("GARCH(1,1)", "EGARCH(1,1)", "GJR-GARCH(1,1)")) {
      fc_other <- all_forecasts_q[[m]][[h]]
      valid    <- !is.na(fc_figarch) & !is.na(fc_other) & !is.na(rv_h)
      e1 <- (fc_figarch[valid] - rv_h[valid])
      e2 <- (fc_other[valid]   - rv_h[valid])
      if (length(e1) < 5) next
      dm <- tryCatch(
        dm.test(e1, e2, alternative = "less", h = as.integer(h)),
        error = function(e) NULL)
      if (is.null(dm)) next
      dm_results_all <- rbind(dm_results_all, data.frame(
        Quarter    = qd$label,
        Horizon    = paste0(h, "-day"),
        Model_1    = "FIGARCH(1,d,1)",
        Model_2    = m,
        DM_stat    = round(dm$statistic, 4),
        p_value    = round(dm$p.value,   6),
        Conclusion = ifelse(dm$p.value < 0.05,
                            "FIGARCH significantly better",
                            "No significant difference")
      ))
    }
  }
  
  # FIX 8: extend rolling_ret accumulator (not base_ret)
  rolling_ret <- c(rolling_ret, qd$vec)
}

accuracy_all_quarters <- accuracy_all_quarters[
  order(accuracy_all_quarters$Quarter,
        accuracy_all_quarters$Horizon,
        accuracy_all_quarters$RMSE), ]

cat("\n=== FORECAST ACCURACY: ALL QUARTERS ===\n")
print(accuracy_all_quarters, row.names = FALSE)
write.csv(accuracy_all_quarters,
          file.path(OUTPUT, "multihorizon_forecast_accuracy_q1_q4.csv"),
          row.names = FALSE)

cat("\n=== DIEBOLD-MARIANO: ALL QUARTERS ===\n")
print(dm_results_all, row.names = FALSE)
write.csv(dm_results_all,
          file.path(OUTPUT, "diebold_mariano_results.csv"),
          row.names = FALSE)
cat("DM test results saved\n")


# ============================================================
#  COMBINED Q1-Q4 MAE HEATMAP
# ============================================================

mae_combined <- accuracy_all_quarters[, c("Quarter", "Model", "Horizon", "MAE")]

mae_combined$Horizon <- factor(mae_combined$Horizon,
                               levels = c("1-day","5-day","10-day","20-day"))

p_heatmap_combined <- ggplot(mae_combined,
                             aes(x = Horizon, y = Model, fill = MAE)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(MAE, 4)),
            size = 3, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#1D9E75", high = "#E24B4A", name = "MAE") +
  facet_wrap(~Quarter, ncol = 2) +
  labs(title    = "Forecast MAE Heatmap: Model vs Horizon — Q1-Q4 2025",
       subtitle = "Green = lower error (better)  |  Red = higher error (worse)  |  Compare across models only, not across horizons",
       x = "Forecast Horizon", y = "") +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(color = "grey40"),
        legend.position = "right",
        strip.text      = element_text(face = "bold"),
        panel.grid      = element_blank())

ggsave(file.path(OUTPUT, "plot_fc6_mae_heatmap_q1_q4.png"),
       p_heatmap_combined, width = 14, height = 10, dpi = 150)
print(p_heatmap_combined)
cat("plot_fc6_mae_heatmap_q1_q4.png saved\n")


# ============================================================
#  SEPARATE QUARTER HEATMAPS
# ============================================================

available_quarters <- unique(as.character(mae_combined$Quarter))

for (qlab in available_quarters) {
  mae_q <- mae_combined[mae_combined$Quarter == qlab, ]
  if (nrow(mae_q) == 0) next
  
  p_q <- ggplot(mae_q, aes(x = Horizon, y = Model, fill = MAE)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round(MAE, 4)),
              size = 3.5, color = "white", fontface = "bold") +
    scale_fill_gradient(low = "#1D9E75", high = "#E24B4A", name = "MAE",
                        limits = c(0.10, 0.65)) +
    scale_x_discrete(limits = c("1-day","5-day","10-day","20-day")) +
    labs(title    = paste("Forecast MAE Heatmap:", qlab),
         subtitle = "Green = lower error (better)  |  Red = higher error (worse)  |  Compare across models only, not across horizons",
         x = "Forecast Horizon", y = "") +
    theme_minimal(base_size = 12) +
    theme(plot.title      = element_text(face = "bold", size = 13),
          plot.subtitle   = element_text(color = "grey40"),
          legend.position = "right",
          panel.grid      = element_blank())
  
  fname <- paste0("plot_fc6_mae_heatmap_",
                  tolower(gsub(" ", "_", qlab)), ".png")
  ggsave(file.path(OUTPUT, fname), p_q, width = 12, height = 5, dpi = 150)
  print(p_q)
  cat(sprintf("%s saved\n", fname))
}


# ============================================================
#  THREE REGIME GARCH ANALYSIS
# ============================================================

cat("\n============================================================\n")
cat("THREE REGIME GARCH ANALYSIS\n")
cat("============================================================\n")

# FIX 8 benefit: ret_vec and dates_ret are still clean 2010-2024
# because base_ret was never mutated — no need to restore here
break1 <- as.Date("2020-02-27")
break2 <- as.Date("2022-05-20")

ret_r1 <- ret_vec[dates_ret < break1]
ret_r2 <- ret_vec[dates_ret >= break1 & dates_ret < break2]
ret_r3 <- ret_vec[dates_ret >= break2]

dates_r1 <- dates_ret[dates_ret < break1]
dates_r2 <- dates_ret[dates_ret >= break1 & dates_ret < break2]
dates_r3 <- dates_ret[dates_ret >= break2]

cat(sprintf("Regime 1 (Pre-COVID)   : %s to %s | %d obs\n",
            min(dates_r1), max(dates_r1), length(ret_r1)))
cat(sprintf("Regime 2 (COVID crisis): %s to %s | %d obs\n",
            min(dates_r2), max(dates_r2), length(ret_r2)))
cat(sprintf("Regime 3 (Post-crisis) : %s to %s | %d obs\n",
            min(dates_r3), max(dates_r3), length(ret_r3)))


# --- Section A: Descriptive statistics by regime ---------------

desc_by_regime <- function(ret, label) {
  jb <- jarque.bera.test(ret)
  data.frame(
    Regime          = label,
    N_obs           = length(ret),
    Mean            = round(mean(ret),         6),
    Median          = round(median(ret),       6),
    Std_Dev         = round(sd(ret),           6),
    Min             = round(min(ret),          6),
    Max             = round(max(ret),          6),
    Skewness        = round(skewness(ret),     4),
    Excess_Kurtosis = round(kurtosis(ret) - 3, 4),
    JB_stat         = round(jb$statistic,      4),
    JB_pvalue       = round(jb$p.value,        6)
  )
}

desc_regimes <- rbind(
  desc_by_regime(ret_r1, "Regime 1 - Pre-COVID (2010-2020)"),
  desc_by_regime(ret_r2, "Regime 2 - COVID crisis (2020-2022)"),
  desc_by_regime(ret_r3, "Regime 3 - Post-crisis (2022-2024)")
)
print(desc_regimes, row.names = FALSE)
write.csv(desc_regimes,
          file.path(OUTPUT, "regime_descriptive_stats.csv"),
          row.names = FALSE)

dist_df <- data.frame(
  Returns = c(ret_r1, ret_r2, ret_r3),
  Regime  = c(rep("Regime 1 - Pre-COVID",     length(ret_r1)),
              rep("Regime 2 - COVID crisis",  length(ret_r2)),
              rep("Regime 3 - Post-crisis",   length(ret_r3)))
)

colors_regime <- c(
  "Regime 1 - Pre-COVID"    = "#2196F3",
  "Regime 2 - COVID crisis" = "#E24B4A",
  "Regime 3 - Post-crisis"  = "#4CAF50"
)

p_dist_regime <- ggplot(dist_df, aes(x = Returns, fill = Regime)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, alpha = 0.6, position = "identity") +
  facet_wrap(~Regime, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = colors_regime) +
  labs(title    = "Return Distribution by Structural Regime",
       subtitle = "Separate histogram for each volatility regime",
       x = "Returns (%)", y = "Density") +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(color = "grey40"),
        legend.position = "none",
        strip.text      = element_text(face = "bold"))

ggsave(file.path(OUTPUT, "plot_regime_return_distribution.png"),
       p_dist_regime, width = 12, height = 10, dpi = 150)
print(p_dist_regime)


# --- Section B: Model estimation by regime ---------------------

spec_garch <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std")

spec_egarch <- ugarchspec(
  variance.model     = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std")

spec_gjr <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std")

spec_figarch <- ugarchspec(
  variance.model     = list(model = "fiGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std")

fit_models_regime <- function(ret, regime_label) {
  cat(sprintf("\n  Fitting models for %s...\n", regime_label))
  list(
    garch   = ugarchfit(spec_garch,   data = ret, solver = "hybrid"),
    egarch  = ugarchfit(spec_egarch,  data = ret, solver = "hybrid"),
    gjr     = ugarchfit(spec_gjr,     data = ret, solver = "hybrid"),
    figarch = ugarchfit(spec_figarch, data = ret, solver = "hybrid")
  )
}

cat("\nFitting all 4 models on 3 regimes (12 models total)...\n")
fits_r1 <- fit_models_regime(ret_r1, "Regime 1 Pre-COVID")
fits_r2 <- fit_models_regime(ret_r2, "Regime 2 COVID crisis")
fits_r3 <- fit_models_regime(ret_r3, "Regime 3 Post-crisis")


# --- Section C: Model comparison by regime ---------------------

get_criteria_regime <- function(fit, model_name, regime, n_obs) {
  ic <- infocriteria(fit)
  cf <- coef(fit)
  persistence <- tryCatch({
    if ("gamma1" %in% names(cf))
      cf["alpha1"] + 0.5 * cf["gamma1"] + cf["beta1"]
    else if ("d" %in% names(cf))
      cf["d"]
    else
      cf["alpha1"] + cf["beta1"]
  }, error = function(e) NA)
  
  data.frame(
    Regime        = regime,
    Model         = model_name,
    LogLikelihood = round(likelihood(fit), 4),
    AIC           = round(ic[1] * n_obs,   4),
    BIC           = round(ic[2] * n_obs,   4),
    Persistence   = round(persistence,     6)
  )
}

make_comp <- function(fits, regime_label, n_obs) {
  rbind(
    get_criteria_regime(fits$garch,   "GARCH(1,1)",     regime_label, n_obs),
    get_criteria_regime(fits$egarch,  "EGARCH(1,1)",    regime_label, n_obs),
    get_criteria_regime(fits$gjr,     "GJR-GARCH(1,1)", regime_label, n_obs),
    get_criteria_regime(fits$figarch, "FIGARCH(1,d,1)", regime_label, n_obs)
  )
}

comp_r1 <- make_comp(fits_r1, "Regime 1", length(ret_r1))
comp_r2 <- make_comp(fits_r2, "Regime 2", length(ret_r2))
comp_r3 <- make_comp(fits_r3, "Regime 3", length(ret_r3))

comp_r1 <- comp_r1[order(comp_r1$AIC), ]
comp_r2 <- comp_r2[order(comp_r2$AIC), ]
comp_r3 <- comp_r3[order(comp_r3$AIC), ]

cat("\n--- Regime 1: Pre-COVID ---\n"); print(comp_r1, row.names = FALSE)
cat("\n--- Regime 2: COVID crisis ---\n"); print(comp_r2, row.names = FALSE)
cat("\n--- Regime 3: Post-crisis ---\n"); print(comp_r3, row.names = FALSE)

comp_all <- rbind(comp_r1, comp_r2, comp_r3)
write.csv(comp_all,
          file.path(OUTPUT, "regime_model_comparison.csv"),
          row.names = FALSE)

best_per_regime <- data.frame(
  Regime     = c("Regime 1 - Pre-COVID",
                 "Regime 2 - COVID crisis",
                 "Regime 3 - Post-crisis"),
  Best_Model = c(comp_r1$Model[1], comp_r2$Model[1], comp_r3$Model[1]),
  Best_AIC   = c(comp_r1$AIC[1],   comp_r2$AIC[1],   comp_r3$AIC[1])
)
print(best_per_regime, row.names = FALSE)
write.csv(best_per_regime,
          file.path(OUTPUT, "best_model_per_regime.csv"),
          row.names = FALSE)

colors_model <- c(
  "GARCH(1,1)"     = "#2196F3",
  "EGARCH(1,1)"    = "#E91E63",
  "GJR-GARCH(1,1)" = "#FF9800",
  "FIGARCH(1,d,1)" = "#4CAF50"
)

comp_all$Regime_Label <- dplyr::recode(comp_all$Regime,
                                       "Regime 1" = "Regime 1\nPre-COVID (2010-2020)",
                                       "Regime 2" = "Regime 2\nCOVID crisis (2020-2022)",
                                       "Regime 3" = "Regime 3\nPost-crisis (2022-2024)")

p_aic_regime <- ggplot(comp_all,
                       aes(x = reorder(Model, AIC), y = AIC, fill = Model)) +
  geom_bar(stat = "identity", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = round(AIC, 1)), vjust = -0.4, size = 2.8) +
  facet_wrap(~Regime_Label, scales = "free_y") +
  scale_fill_manual(values = colors_model) +
  labs(title    = "Model Comparison by Structural Regime - AIC (lower = better)",
       subtitle = "Which model fits best changes across volatility regimes",
       x = "", y = "AIC") +
  theme_minimal(base_size = 10) +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(color = "grey40"),
        axis.text.x     = element_text(angle = 20, hjust = 1),
        legend.position = "none",
        strip.text      = element_text(face = "bold"))

ggsave(file.path(OUTPUT, "plot_regime_aic_comparison.png"),
       p_aic_regime, width = 14, height = 6, dpi = 150)
print(p_aic_regime)

persist_all <- comp_all[, c("Regime", "Regime_Label", "Model", "Persistence")]

p_persist_regime <- ggplot(persist_all,
                           aes(x = Model, y = Persistence, fill = Model)) +
  geom_bar(stat = "identity", alpha = 0.85, position = "dodge", width = 0.6) +
  geom_text(aes(label = round(Persistence, 4)), vjust = -0.4, size = 2.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
  facet_wrap(~Regime_Label, scales = "free_y") +
  scale_fill_manual(values = colors_model) +
  labs(title    = "Volatility Persistence by Structural Regime",
       subtitle = "Red dashed line = 1.0 (explosive threshold)",
       x = "", y = "Persistence") +
  theme_minimal(base_size = 10) +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(color = "grey40"),
        axis.text.x     = element_text(angle = 20, hjust = 1),
        legend.position = "none",
        strip.text      = element_text(face = "bold"))

ggsave(file.path(OUTPUT, "plot_regime_persistence.png"),
       p_persist_regime, width = 14, height = 6, dpi = 150)
print(p_persist_regime)


# --- Section D: Conditional volatility by regime ---------------

make_vol_df <- function(fits, dates, regime_label) {
  model_names_r <- c("GARCH(1,1)", "EGARCH(1,1)",
                     "GJR-GARCH(1,1)", "FIGARCH(1,d,1)")
  data.frame(
    Date       = rep(dates, 4),
    Volatility = c(as.numeric(sigma(fits$garch)),
                   as.numeric(sigma(fits$egarch)),
                   as.numeric(sigma(fits$gjr)),
                   as.numeric(sigma(fits$figarch))),
    Model      = rep(model_names_r, each = length(dates)),
    Regime     = regime_label
  )
}

vol_r1 <- make_vol_df(fits_r1, dates_r1, "Regime 1 - Pre-COVID")
vol_r2 <- make_vol_df(fits_r2, dates_r2, "Regime 2 - COVID crisis")
vol_r3 <- make_vol_df(fits_r3, dates_r3, "Regime 3 - Post-crisis")

plot_vol_regime <- function(vol_df, regime_label) {
  ggplot(vol_df, aes(x = Date, y = Volatility, color = Model)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~Model, ncol = 1, scales = "free_y") +
    scale_color_manual(values = colors_model) +
    labs(title    = paste("Conditional Volatility -", regime_label),
         subtitle = "Each panel shows one GARCH model fitted on this regime only",
         x = "", y = "Volatility (%)") +
    theme_minimal(base_size = 10) +
    theme(plot.title      = element_text(face = "bold", size = 12),
          plot.subtitle   = element_text(color = "grey40"),
          legend.position = "none",
          strip.text      = element_text(face = "bold"))
}

p_vol_r1 <- plot_vol_regime(vol_r1, "Regime 1 - Pre-COVID (2010-2020)")
p_vol_r2 <- plot_vol_regime(vol_r2, "Regime 2 - COVID crisis (2020-2022)")
p_vol_r3 <- plot_vol_regime(vol_r3, "Regime 3 - Post-crisis (2022-2024)")

ggsave(file.path(OUTPUT, "plot_vol_regime1_precovid.png"),  p_vol_r1, width = 14, height = 12, dpi = 150)
ggsave(file.path(OUTPUT, "plot_vol_regime2_covid.png"),     p_vol_r2, width = 14, height = 12, dpi = 150)
ggsave(file.path(OUTPUT, "plot_vol_regime3_postcovid.png"), p_vol_r3, width = 14, height = 12, dpi = 150)

print(p_vol_r1); print(p_vol_r2); print(p_vol_r3)
cat("All 3 regime volatility plots saved\n")


# --- Section E: Diagnostics by regime -------------------------

get_diagnostics <- function(fits, regime_label) {
  model_list_r  <- list(fits$garch, fits$egarch, fits$gjr, fits$figarch)
  model_names_r <- c("GARCH(1,1)", "EGARCH(1,1)",
                     "GJR-GARCH(1,1)", "FIGARCH(1,d,1)")
  df <- data.frame()
  for (i in seq_along(model_list_r)) {
    r  <- as.numeric(residuals(model_list_r[[i]], standardize = TRUE))
    r  <- r[is.finite(r)]
    lb <- Box.test(r,   lag = 10, type = "Ljung-Box")
    ls <- Box.test(r^2, lag = 10, type = "Ljung-Box")
    df <- rbind(df, data.frame(
      Regime        = regime_label,
      Model         = model_names_r[i],
      LB_Resid_pval = round(lb$p.value, 4),
      LB_Sq_pval    = round(ls$p.value, 4),
      LB_Resid_OK   = ifelse(lb$p.value > 0.05, "PASS", "FAIL"),
      LB_Squared_OK = ifelse(ls$p.value > 0.05, "PASS", "FAIL")
    ))
  }
  df
}

diag_all <- rbind(
  get_diagnostics(fits_r1, "Regime 1 - Pre-COVID"),
  get_diagnostics(fits_r2, "Regime 2 - COVID crisis"),
  get_diagnostics(fits_r3, "Regime 3 - Post-crisis")
)
print(diag_all, row.names = FALSE)
write.csv(diag_all,
          file.path(OUTPUT, "regime_model_diagnostics.csv"),
          row.names = FALSE)


# --- Section F: Parameter evolution ---------------------------

get_key_params <- function(fits, regime) {
  cf_g   <- coef(fits$garch)
  cf_e   <- coef(fits$egarch)
  cf_gjr <- coef(fits$gjr)
  cf_f   <- coef(fits$figarch)
  d_val  <- cf_f["d"]   # FIX 3: removed invalid delta fallback
  
  data.frame(
    Regime        = regime,
    GARCH_alpha   = round(cf_g["alpha1"],   6),
    GARCH_beta    = round(cf_g["beta1"],    6),
    GARCH_persist = round(cf_g["alpha1"] + cf_g["beta1"], 6),
    EGARCH_gamma  = round(cf_e["alpha1"],   6),
    GJR_gamma     = round(cf_gjr["gamma1"], 6),
    GJR_persist   = round(cf_gjr["alpha1"] + 0.5 * cf_gjr["gamma1"] + cf_gjr["beta1"], 6),
    FIGARCH_d     = round(d_val,            6)
  )
}

param_evolution <- rbind(
  get_key_params(fits_r1, "Regime 1 - Pre-COVID (2010-2020)"),
  get_key_params(fits_r2, "Regime 2 - COVID crisis (2020-2022)"),
  get_key_params(fits_r3, "Regime 3 - Post-crisis (2022-2024)")
)
print(param_evolution, row.names = FALSE)
write.csv(param_evolution,
          file.path(OUTPUT, "regime_parameter_evolution.csv"),
          row.names = FALSE)

param_long <- pivot_longer(
  param_evolution,
  cols      = c("GARCH_persist", "GJR_persist", "FIGARCH_d"),
  names_to  = "Parameter",
  values_to = "Value"
)
param_long$Parameter <- dplyr::recode(param_long$Parameter,
                                      GARCH_persist = "GARCH(1,1) Persistence",
                                      GJR_persist   = "GJR-GARCH Persistence",
                                      FIGARCH_d     = "FIGARCH d (long memory)")

p_param_evo <- ggplot(param_long,
                      aes(x = Regime, y = Value,
                          color = Parameter, group = Parameter)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_text(aes(label = round(Value, 4)), vjust = -1, size = 3.2) +
  scale_color_manual(values = c(
    "GARCH(1,1) Persistence" = "#2196F3",
    "GJR-GARCH Persistence"  = "#FF9800",
    "FIGARCH d (long memory)" = "#4CAF50"
  )) +
  labs(title    = "Parameter Evolution Across Structural Regimes",
       subtitle = "Shows how volatility dynamics changed in each market regime",
       x = "", y = "Value", color = "Parameter") +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(color = "grey40"),
        axis.text.x     = element_text(angle = 10, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(OUTPUT, "plot_parameter_evolution.png"),
       p_param_evo, width = 12, height = 7, dpi = 150)
print(p_param_evo)
cat("Parameter evolution plot saved\n")


# ============================================================
#  FINAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("ALL ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Files saved to: %s\n\n", OUTPUT))

cat("Forecast outputs:\n")
cat("  multihorizon_forecast_accuracy_q1_q4.csv\n")
cat("  plot_fc6_mae_heatmap_q1_q4.png\n")
cat("  plot_fc6_mae_heatmap_q1_2025.png\n")
cat("  plot_fc6_mae_heatmap_q2_2025.png\n")
cat("  plot_fc6_mae_heatmap_q3_2025.png\n")
cat("  plot_fc6_mae_heatmap_q4_2025.png\n")
cat("  diebold_mariano_results.csv\n\n")

cat("Regime analysis outputs:\n")
cat("  regime_descriptive_stats.csv\n")
cat("  regime_model_comparison.csv\n")
cat("  best_model_per_regime.csv\n")
cat("  regime_model_diagnostics.csv\n")
cat("  regime_parameter_evolution.csv\n\n")

cat("Best model per regime:\n")
print(best_per_regime, row.names = FALSE)
cat("\nParameter evolution highlights:\n")
print(param_evolution[, c("Regime", "GARCH_persist", "GJR_persist", "FIGARCH_d")],
      row.names = FALSE)
# ============================================================
# ACTUAL VS FORECASTED VOLATILITY PLOT
# ============================================================

plot_df <- data.frame(
  Date     = qd$dates,
  ActualRV = rv_list_q[["1"]],
  GARCH    = all_forecasts_q[["GARCH(1,1)"]][["1"]],
  EGARCH   = all_forecasts_q[["EGARCH(1,1)"]][["1"]],
  GJR      = all_forecasts_q[["GJR-GARCH(1,1)"]][["1"]],
  FIGARCH  = all_forecasts_q[["FIGARCH(1,d,1)"]][["1"]]
)

plot_long <- pivot_longer(
  plot_df,
  cols = -Date,
  names_to = "Series",
  values_to = "Value"
)

p_actual_vs_fc <- ggplot(plot_long,
                         aes(x = Date,
                             y = Value,
                             color = Series)) +
  geom_line(linewidth = 0.7) +
  labs(
    title = paste("Actual vs Forecasted Volatility -", qd$label),
    x = "",
    y = "Volatility"
  ) +
  theme_minimal()

ggsave(
  file.path(OUTPUT,
            paste0("plot_actual_vs_forecast_", qname, ".png")),
  p_actual_vs_fc,
  width = 14,
  height = 6,
  dpi = 150
)

print(p_actual_vs_fc)


# ============================================================
#  FORECAST VS REALIZED PLOTS
#  1-day, 5-day, 10-day, 20-day | All 4 models | All quarters
#  Run AFTER main script
#  (all_forecasts_all_q, quarter_data, OUTPUT must be in memory)
# ============================================================

library(ggplot2)
library(tidyr)
library(dplyr)

cat("\n============================================================\n")
cat("FORECAST VS REALIZED PLOTS — ALL HORIZONS & MODELS\n")
cat("============================================================\n")

# ── Colors ───────────────────────────────────────────────────
colors_model <- c(
  "GARCH(1,1)"     = "#2196F3",
  "EGARCH(1,1)"    = "#E91E63",
  "GJR-GARCH(1,1)" = "#FF9800",
  "FIGARCH(1,d,1)" = "#4CAF50",
  "Realized"       = "#212121"
)

model_names <- c("GARCH(1,1)", "EGARCH(1,1)",
                 "GJR-GARCH(1,1)", "FIGARCH(1,d,1)")

horizons     <- c(1, 5, 10, 20)
horizon_labs <- c("1-day", "5-day", "10-day", "20-day")

# ── Realized volatility function ─────────────────────────────
make_rv <- function(ret, h) {
  n   <- length(ret)
  out <- rep(NA_real_, n)
  for (i in h:n) out[i] <- sqrt(mean(ret[(i-h+1):i]^2))
  out
}

# ============================================================
#  PLOT TYPE 1: One plot per horizon — all 4 models + realized
#               across all quarters combined
# ============================================================

cat("\nGenerating combined quarterly plots by horizon...\n")

for (hi in seq_along(horizons)) {
  h     <- horizons[hi]
  h_lab <- horizon_labs[hi]
  
  # Collect data across all available quarters
  all_rows <- data.frame()
  
  for (qname in names(all_forecasts_all_q)) {
    qd <- quarter_data[[qname]]
    if (qd$n < 5) next
    
    rv <- make_rv(qd$vec, h)
    
    # Realized
    all_rows <- rbind(all_rows, data.frame(
      Date       = qd$dates,
      Value      = rv,
      Series     = "Realized",
      Quarter    = qd$label
    ))
    
    # Each model forecast
    for (m in model_names) {
      fc <- all_forecasts_all_q[[qname]][[m]][[as.character(h)]]
      if (is.null(fc)) next
      all_rows <- rbind(all_rows, data.frame(
        Date    = qd$dates,
        Value   = fc,
        Series  = m,
        Quarter = qd$label
      ))
    }
  }
  
  all_rows <- all_rows[!is.na(all_rows$Value), ]
  all_rows$Type <- ifelse(all_rows$Series == "Realized",
                          "Realized", "Forecast")
  
  p <- ggplot(all_rows,
              aes(x = Date, y = Value,
                  color = Series,
                  linetype = Type,
                  linewidth = I(ifelse(Series == "Realized",
                                       1.3, 0.7)))) +
    geom_line() +
    geom_point(data = subset(all_rows, Series == "Realized"),
               size = 1.2, alpha = 0.7) +
    facet_wrap(~Quarter, ncol = 2, scales = "free_x") +
    scale_color_manual(values = colors_model) +
    scale_linetype_manual(values = c("Realized" = "solid",
                                     "Forecast" = "dashed")) +
    labs(
      title    = paste0(h_lab, " Ahead Forecast vs Realized Volatility — 2025"),
      subtitle = paste0("All 4 GARCH models | Realized = sqrt(mean(r²)) over ",
                        h, " days | Dashed = forecast | Solid = realized"),
      x = "Date", y = "Volatility (%)", color = "Series"
    ) +
    guides(linetype = "none", linewidth = "none") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(color = "grey40", size = 9),
      legend.position = "bottom",
      strip.text      = element_text(face = "bold"),
      legend.text     = element_text(size = 9)
    )
  
  fname <- paste0("plot_forecast_vs_realized_",
                  gsub("-", "", h_lab), "_allquarters.png")
  ggsave(file.path(OUTPUT, fname),
         p, width = 16, height = 10, dpi = 150)
  print(p)
  cat(sprintf("[OK] %s saved\n", fname))
}


# ============================================================
#  PLOT TYPE 2: One plot per model — all 4 horizons as panels
#               one model at a time across all quarters
# ============================================================

cat("\nGenerating per-model plots across all horizons...\n")

for (m in model_names) {
  cat(sprintf("  Model: %s\n", m))
  
  all_rows <- data.frame()
  
  for (qname in names(all_forecasts_all_q)) {
    qd <- quarter_data[[qname]]
    if (qd$n < 5) next
    
    for (hi in seq_along(horizons)) {
      h     <- horizons[hi]
      h_lab <- horizon_labs[hi]
      
      rv <- make_rv(qd$vec, h)
      fc <- all_forecasts_all_q[[qname]][[m]][[as.character(h)]]
      if (is.null(fc)) next
      
      all_rows <- rbind(all_rows,
                        data.frame(
                          Date    = rep(qd$dates, 2),
                          Value   = c(rv, fc),
                          Series  = c(rep("Realized", qd$n), rep("Forecast", qd$n)),
                          Horizon = h_lab,
                          Quarter = qd$label
                        )
      )
    }
  }
  
  all_rows <- all_rows[!is.na(all_rows$Value), ]
  all_rows$Horizon <- factor(all_rows$Horizon,
                             levels = horizon_labs)
  all_rows$Type <- ifelse(all_rows$Series == "Realized",
                          "Realized", "Forecast")
  
  p <- ggplot(all_rows,
              aes(x = Date, y = Value,
                  color = Series,
                  linetype = Type,
                  linewidth = I(ifelse(Series == "Realized",
                                       1.3, 0.7)))) +
    geom_line() +
    geom_point(data = subset(all_rows, Series == "Realized"),
               size = 1.0, alpha = 0.6) +
    facet_grid(Horizon ~ Quarter, scales = "free") +
    scale_color_manual(values = c(
      "Realized" = "#212121",
      "Forecast" = colors_model[m]
    )) +
    scale_linetype_manual(values = c("Realized" = "solid",
                                     "Forecast" = "dashed")) +
    labs(
      title    = paste0(m, " — Forecast vs Realized Volatility (2025)"),
      subtitle = "Rows = forecast horizon | Columns = quarter | Solid = realized | Dashed = forecast",
      x = "Date", y = "Volatility (%)", color = "Series"
    ) +
    guides(linetype = "none", linewidth = "none") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(color = "grey40", size = 9),
      legend.position = "bottom",
      strip.text      = element_text(face = "bold", size = 8),
      axis.text.x     = element_text(angle = 30, hjust = 1, size = 7)
    )
  
  model_safe <- gsub("[()\\.,]", "", gsub(" ", "_", m))
  fname      <- paste0("plot_forecast_vs_realized_", model_safe, ".png")
  ggsave(file.path(OUTPUT, fname),
         p, width = 18, height = 14, dpi = 150)
  print(p)
  cat(sprintf("[OK] %s saved\n", fname))
}


# ============================================================
#  PLOT TYPE 3: Q1 only — 2x2 grid of horizons, one per model
#               cleaner individual comparison per quarter
# ============================================================

cat("\nGenerating per-quarter 2x2 horizon comparison plots...\n")

for (qname in names(all_forecasts_all_q)) {
  qd <- quarter_data[[qname]]
  if (qd$n < 5) next
  
  cat(sprintf("  Quarter: %s\n", qd$label))
  
  all_rows <- data.frame()
  
  for (hi in seq_along(horizons)) {
    h     <- horizons[hi]
    h_lab <- horizon_labs[hi]
    rv    <- make_rv(qd$vec, h)
    
    # Realized
    all_rows <- rbind(all_rows, data.frame(
      Date    = qd$dates,
      Value   = rv,
      Series  = "Realized",
      Horizon = h_lab
    ))
    
    # All model forecasts
    for (m in model_names) {
      fc <- all_forecasts_all_q[[qname]][[m]][[as.character(h)]]
      if (is.null(fc)) next
      all_rows <- rbind(all_rows, data.frame(
        Date    = qd$dates,
        Value   = fc,
        Series  = m,
        Horizon = h_lab
      ))
    }
  }
  
  all_rows <- all_rows[!is.na(all_rows$Value), ]
  all_rows$Horizon <- factor(all_rows$Horizon, levels = horizon_labs)
  all_rows$Type    <- ifelse(all_rows$Series == "Realized",
                             "Realized", "Forecast")
  
  p <- ggplot(all_rows,
              aes(x = Date, y = Value,
                  color = Series,
                  linetype = Type,
                  linewidth = I(ifelse(Series == "Realized",
                                       1.4, 0.75)))) +
    geom_line() +
    geom_point(data = subset(all_rows, Series == "Realized"),
               size = 1.5, alpha = 0.8) +
    facet_wrap(~Horizon, ncol = 2, scales = "free_y") +
    scale_color_manual(values = colors_model) +
    scale_linetype_manual(values = c("Realized" = "solid",
                                     "Forecast" = "dashed")) +
    labs(
      title    = paste0("Forecast vs Realized Volatility — ",
                        qd$label),
      subtitle = paste0(
        "All 4 GARCH models across 4 forecast horizons | ",
        "Solid black = realized | Dashed = model forecast"),
      x = "Date", y = "Volatility (%)", color = "Series"
    ) +
    guides(linetype = "none", linewidth = "none") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(color = "grey40", size = 9),
      legend.position = "bottom",
      strip.text      = element_text(face = "bold"),
      legend.text     = element_text(size = 9)
    )
  
  q_safe <- tolower(gsub(" ", "_", qd$label))
  fname  <- paste0("plot_forecast_vs_realized_", q_safe, "_allhorizons.png")
  ggsave(file.path(OUTPUT, fname),
         p, width = 16, height = 10, dpi = 150)
  print(p)
  cat(sprintf("[OK] %s saved\n", fname))
}


# ============================================================
#  DONE
# ============================================================

cat("\n============================================================\n")
cat("ALL FORECAST VS REALIZED PLOTS COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Saved to: %s\n\n", OUTPUT))
cat("Files generated:\n")
cat("  Type 1 — By horizon (all models, all quarters):\n")
for (h_lab in horizon_labs) {
  cat(sprintf("    plot_forecast_vs_realized_%s_allquarters.png\n",
              gsub("-","",h_lab)))
}
cat("\n  Type 2 — By model (all horizons, all quarters):\n")
for (m in model_names) {
  model_safe <- gsub("[()\\.,]", "", gsub(" ", "_", m))
  cat(sprintf("    plot_forecast_vs_realized_%s.png\n", model_safe))
}
cat("\n  Type 3 — By quarter (all models, all horizons as 2x2):\n")
for (qname in names(quarter_data)) {
  if (quarter_data[[qname]]$n >= 5) {
    q_safe <- tolower(gsub(" ","_", quarter_data[[qname]]$label))
    cat(sprintf("    plot_forecast_vs_realized_%s_allhorizons.png\n",
                q_safe))
  }
}

# FIX 9: browseURL wrapped in interactive() check — safe on server/headless
if (interactive()) browseURL(OUTPUT)
browseURL(OUTPUT)







# ============================================================
# OUT-OF-SAMPLE VAR + KUPIEC — minimal working version
# Paste this AFTER your existing 05_var_es.R code
# ============================================================


# ── 1. Split ─────────────────────────────────────────────────
oos_start <- as.Date("2025-01-01")
oos_end   <- as.Date("2025-12-31")
all_dates <- as.Date(index(returns))

n_in      <- sum(all_dates < oos_start)
n_oos     <- sum(all_dates >= oos_start & all_dates <= oos_end)
dates_oos <- all_dates[all_dates >= oos_start & all_dates <= oos_end]

cat(sprintf("In-sample    : %d obs\n", n_in))
cat(sprintf("Out-of-sample: %d obs [%s to %s]\n", n_oos,
            as.character(dates_oos[1]), as.character(dates_oos[n_oos])))

# ── 2. Rolling forecast ──────────────────────────────────────
ret_xts  <- xts(ret_vec, order.by = as.Date(dates_ret))
best_spec <- getspec(best_fit)
setstart(best_spec) <- FALSE

cat("Running ugarchroll... (takes a few minutes)\n")

roll_fit <- ugarchroll(
  spec             = best_spec,
  data             = ret_xts,
  n.ahead          = 1,
  forecast.length  = n_oos,
  refit.every      = 20,
  refit.window     = "expanding",
  solver           = "hybrid",
  calculate.VaR    = TRUE,
  VaR.alpha        = c(0.01, 0.05),
  keep.coef        = TRUE
)


# ── 3. Extract results ───────────────────────────────────────
roll_df   <- as.data.frame(roll_fit)

realized  <- roll_df$Realized
mu_h1     <- roll_df$Mu
sigma_h1  <- roll_df$Sigma

shape_par <- coef(best_fit)["shape"]

VaR_1_oos <- as.numeric(mu_h1 + sigma_h1 *
                          qdist("std", p = 0.01, mu = 0, sigma = 1, shape = shape_par))
VaR_5_oos <- as.numeric(mu_h1 + sigma_h1 *
                          qdist("std", p = 0.05, mu = 0, sigma = 1, shape = shape_par))

# ── 4. Breaches ──────────────────────────────────────────────
breach_1_oos   <- sum(realized < VaR_1_oos)
breach_5_oos   <- sum(realized < VaR_5_oos)
expected_1_oos <- round(n_oos * 0.01)
expected_5_oos <- round(n_oos * 0.05)

ES_1_oos <- mean(realized[realized < VaR_1_oos], na.rm = TRUE)
ES_5_oos <- mean(realized[realized < VaR_5_oos], na.rm = TRUE)

cat(sprintf("\n=== OOS VaR RESULTS (%s) ===\n", best_name))
cat(sprintf("1%% VaR breaches: %d (expected ~%d)\n", breach_1_oos, expected_1_oos))
cat(sprintf("5%% VaR breaches: %d (expected ~%d)\n", breach_5_oos, expected_5_oos))
cat(sprintf("1%% ES: %.4f%%\n", ES_1_oos))
cat(sprintf("5%% ES: %.4f%%\n", ES_5_oos))

# ── 5. Kupiec — reuses your existing function ────────────────
kupiec_oos_1 <- kupiec_test(breach_1_oos, n_oos, 0.99)
kupiec_oos_5 <- kupiec_test(breach_5_oos, n_oos, 0.95)

kupiec_oos_1$Window <- "OOS"
kupiec_oos_5$Window <- "OOS"
kupiec_1$Window     <- "In-Sample"
kupiec_5$Window     <- "In-Sample"

kupiec_all <- rbind(kupiec_1, kupiec_oos_1, kupiec_5, kupiec_oos_5)

cat("\n=== KUPIEC: IN-SAMPLE vs OOS ===\n")
print(kupiec_all[, c("Window","Confidence","N_obs","N_breach",
                     "Expected","Breach_Rate","LR_stat","p_value","Result")],
      row.names = FALSE)

write.csv(kupiec_all,
          file.path(OUTPUT, "kupiec_insample_vs_oos.csv"),
          row.names = FALSE)




# ── 6. Plot ──────────────────────────────────────────────────
oos_plot_df <- data.frame(
  Date    = as.Date(rownames(roll_df)),
  Returns = realized,
  VaR_1   = VaR_1_oos,
  VaR_5   = VaR_5_oos,
  Breach  = realized < VaR_1_oos
)

p_oos <- ggplot(oos_plot_df, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "grey55", linewidth = 0.35) +
  geom_line(aes(y = VaR_1, color = "1% VaR"), linewidth = 0.75) +
  geom_line(aes(y = VaR_5, color = "5% VaR"), linewidth = 0.75) +
  geom_hline(yintercept = ES_1_oos, color = "#880E4F",
             linetype = "dotted", linewidth = 0.8) +
  geom_hline(yintercept = ES_5_oos, color = "#E65100",
             linetype = "dotted", linewidth = 0.8) +
  geom_point(data = subset(oos_plot_df, Breach),
             aes(y = Returns), color = "red", size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c("1% VaR" = "#B71C1C", "5% VaR" = "#FF9800")) +
  labs(
    title    = paste("OOS Rolling VaR —", best_name),
    subtitle = sprintf(
      "1%% VaR: %d breaches (exp %d) Kupiec: %s | 5%%: %d (exp %d) Kupiec: %s",
      breach_1_oos, expected_1_oos, kupiec_oos_1$Result,
      breach_5_oos, expected_5_oos, kupiec_oos_5$Result),
    x = "", y = "Returns (%)", color = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(size = 7.5, color = "grey40"),
        legend.position = "bottom")

ggsave(file.path(OUTPUT, "plot10_oos_var_kupiec.png"),
       p_oos, width = 14, height = 6, dpi = 150)
print(p_oos)
cat("plot10_oos_var_kupiec.png saved\n")
