# Investor Shift from Risky Assets to Safe Haven in USA During Major Shocks

## 📌 Overview
This project investigates whether investors in the U.S. shift from risky assets (S&P 500 index) to safe haven assets (gold futures) during major shocks. Using **R** for econometric modeling, the study analyzes two events:
- The **COVID-19 pandemic (2020)**
- The **2025 U.S. tariff disturbances**

The project applies realized correlation, rolling regression, and dummy variable testing to quantify investor behavior during turbulent periods.

---

## 📊 Key Findings
- Investors shifted to gold during major shocks, with **negative and statistically significant correlation** observed.  
- The **2025 tariff disturbance** had a **larger impact** than COVID-19, reflecting stronger market uncertainty.  
- The safe-haven effect was **short-lived and lagged**, rather than persistent.  
- No significant safe-haven dynamics were found during peaceful, non-shock periods.  
- Broader windows diluted the effect, showing that the reaction is sharp but temporary.  

---

## 🧪 Methods
1. **Data Collection**  
   - Downloaded daily S&P 500 (^GSPC) and gold futures (GC=F) data from Yahoo Finance (2010–2025).  
   - Used `quantmod::getSymbols()`.  

2. **Exploratory Data Analysis (EDA)**  
   - Computed descriptive statistics (mean, std, skewness, kurtosis).  
   - Plotted price trends, histograms of log returns, and scatter plots.  

3. **Realized Correlation**  
   - Calculated 21-day rolling realized variance, covariance, and correlation.  
   - Visualized rolling correlation and volatility.  

4. **Event Study with Regression**  
   - Model:  
     ```math
     ρ_t = α + βρ_{t-1} + δ·ShockDummy_t + u_t
     ```  
   - Dummy windows tested via **p-value heatmaps**.  
   - Separate models for COVID-19 and tariff shocks.  
   - Controlled for autocorrelation/heteroskedasticity using **BG, BP tests, and Newey–West robust errors**.  

---

## 📂 Repository Structure
- **Investor-Shift_code.R** – Complete R code for data collection, EDA, rolling correlation, regressions, and plots.  
- **Investor-Shift_report.pdf** – Academic report with methodology, figures, regression tables, and conclusions.  

---

## ⚙️ Requirements
This project uses **R ≥ 4.0** with the following libraries:

```r
install.packages(c(
  "quantmod", "moments", "rugarch", "xts", "ggplot2",
  "dplyr", "tidyr", "zoo", "reshape2", "scales",
  "lubridate", "lmtest", "car", "sandwich", "broom"
))
```

Optional: `naniar` (for NA checks, not essential).

---

## 🚀 Usage
1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/Investor-Shift-from-Risky-Assets-to-Safe-Haven.git
   ```
2. Open **Investor-Shift_code.R** in R or RStudio.  
3. Run sections sequentially:
   - Data download & cleaning  
   - EDA & return analysis  
   - Realized correlation plots  
   - Event regression models (COVID & Tariff)  
4. Review results and figures in the console/plots.  
5. See **Investor-Shift_report.pdf** for the full write-up.  

---

