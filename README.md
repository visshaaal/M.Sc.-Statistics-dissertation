# Volatility Modelling of the Nifty 50 Index Using GARCH-Family Models

M.Sc. Statistics dissertation comparing GARCH(1,1), EGARCH, GJR-GARCH, and FIGARCH models for volatility forecasting and risk estimation on the Nifty 50 index, with a regulatory (Basel III/IV) risk-management framing.

## Overview

This project analyzes ~3,700 daily Nifty 50 observations (2010–2024) to compare how different GARCH-family specifications capture volatility clustering, leverage effects, and long memory in returns. The mean equation is modeled as ARMA(3,2), selected via `auto.arima()`, with volatility dynamics estimated using four competing GARCH specifications.

**Key components:**
- In-sample model fitting and diagnostic testing
- Out-of-sample multi-horizon volatility forecasting
- Value-at-Risk (VaR) and Expected Shortfall (ES) estimation
- Backtesting via Kupiec (unconditional coverage) and Christoffersen (conditional coverage) tests
- Structural break analysis (Bai-Perron)
- Regime-based model comparison (pre-COVID, COVID, post-COVID)
- Diebold-Mariano tests for forecast accuracy comparison

## Methodology

| Component | Approach |
|---|---|
| Mean equation | ARMA(3,2), selected via `auto.arima()` |
| Volatility models | GARCH(1,1), EGARCH, GJR-GARCH, FIGARCH |
| Risk metrics | VaR, Expected Shortfall (multiple horizons) |
| Backtesting | Kupiec test, Christoffersen test |
| Regulatory framing | Basel III/IV capital adequacy context |
| Structural breaks | Bai-Perron multiple breakpoint test |

## Repository Structure

```
├── chapter1.tex – chapter5.tex   # Dissertation chapters
├── main.tex                       # Master LaTeX document
├── abstract.tex, frontpage.tex    # Front matter
├── dissert code                   # R implementation (GARCH modeling, backtesting)
├── *.csv                          # Model outputs (diagnostics, comparisons, backtest results)
├── plot_*.png                     # Forecast, regime, and diagnostic visualizations
├── refrences.bib                  # Bibliography
```

## Tools & Packages

- **R** — `rugarch`, `forecast` (for `auto.arima()`)
- **LaTeX** (Overleaf) — dissertation write-up

## Author

Vishal Verma — M.Sc. Statistics, Central University of Haryana

---
*This repository accompanies an academic dissertation submitted in partial fulfillment of the M.Sc. Statistics program.*
