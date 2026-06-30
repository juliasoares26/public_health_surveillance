set.seed(42)
suppressWarnings(suppressMessages(library(rpart)))

data_url <- "https://raw.githubusercontent.com/weinbergerlab/PCV-mortality-Brazil/master/national.csv"
data_file <- "national.csv"

if (!file.exists(data_file)) {
  cat("Downloading national.csv from Weinberger Lab GitHub repository...\n")
  download.file(data_url, destfile = data_file, mode = "wb", quiet = TRUE)
}
stopifnot(file.exists(data_file))

# 01. Data Preparation
# Aggregate full ICD-10 chapters (J00-J99 respiratory,
# I00-I99 circulatory), restrict to the broadest, non-nested
# age band (3-<60mo, i.e. under-5), and construct covariates

df <- read.csv(data_file, stringsAsFactors = FALSE)
df$date <- as.Date(df$date, format = "%m/%d/%Y")

j_cols <- c("J00_J06", "J09_J18", "J20_J22", "J30_J39", "J40_J47",
            "J60_J70", "J80_J84", "J85_J86", "J90_J94", "J95_J99")
i_cols <- c("I00_I02", "I05_I09", "I10_I15", "I20_I25", "I26_I28",
            "I30_I52", "I60_69", "I70_I79", "I80_I89", "I95_I99")

stopifnot(all(j_cols %in% names(df)))
stopifnot(all(i_cols %in% names(df)))

df$deaths_respiratory <- rowSums(df[, j_cols])
df$deaths_circulatory <- rowSums(df[, i_cols])

# Additional real covariates already present in the source data
keep <- c("age_group", "date", "pop",
          "deaths_respiratory", "deaths_circulatory",
          "J12_J18",  # all-cause pneumonia (PCV-relevant subset)
          "NoPCV")    # control cause-of-death bucket from source study

mort_full <- df[, keep]
mort_full$year <- as.numeric(format(mort_full$date, "%Y"))
mort_full$quarter <- as.numeric(format(mort_full$date, "%m"))
mort_full$quarter <- ceiling(mort_full$quarter / 3)
mort_full$time_index <- as.numeric(factor(mort_full$date))

# Rates per 100,000 population (standard epidemiological scaling)
mort_full$resp_rate <- mort_full$deaths_respiratory / mort_full$pop * 1e5
mort_full$circ_rate <- mort_full$deaths_circulatory / mort_full$pop * 1e5

write.csv(mort_full, "brazil_child_mortality_resp_circ.csv", row.names = FALSE)

cat("\n Data preparation \n")
cat("Rows:", nrow(mort_full), "\n")
cat("Age groups:", paste(unique(mort_full$age_group), collapse = ", "), "\n")
cat("Date range:", as.character(min(mort_full$date)), "to",
    as.character(max(mort_full$date)), "\n\n")
cat("Summary, deaths_respiratory:\n"); print(summary(mort_full$deaths_respiratory))
cat("\nSummary, deaths_circulatory:\n"); print(summary(mort_full$deaths_circulatory))
cat("\nSummary, resp_rate (per 100k):\n"); print(summary(mort_full$resp_rate))
cat("\nSummary, circ_rate (per 100k):\n"); print(summary(mort_full$circ_rate))

# 02. Modeling
# Restrict to the broadest, non-nested age band (3-<60mo) to
# avoid the artificial variance inflation that arises from the source data's nested age bands (3-<12mo subset of 3-<24mo subset of 3-<60mo

# Models:
# 1 Poisson GLM (classical baseline; assumes mean = variance)
# 2 Quasi-Poisson GLM (classical correction for overdispersion)
# 3 Bagged regression-tree ensemble (rpart), an AI ensemble built from R base + rpart only, with bootstrap-based empirical prediction intervals

mort <- subset(mort_full, age_group == "3-<60mo")
mort$date <- as.Date(mort$date)
mort <- mort[order(mort$date), ]
mort$log_pop <- log(mort$pop)
n <- nrow(mort)
cat("\n Modeling \n")
cat("n (under-5, quarterly observations, Brazil 2004-2014):", n, "\n\n")

# Poisson GLM 
glm_pois <- glm(deaths_respiratory ~ time_index + factor(quarter) +
                  offset(log_pop),
                family = poisson(link = "log"), data = mort)

# 2 Quasi-Poisson GLM
glm_qp <- glm(deaths_respiratory ~ time_index + factor(quarter) +
                offset(log_pop),
              family = quasipoisson(link = "log"), data = mort)

pearson_resid <- residuals(glm_pois, type = "pearson")
dispersion <- sum(pearson_resid^2) / glm_pois$df.residual
cat("Pearson dispersion statistic (Poisson model):", round(dispersion, 3), "\n\n")

cat("Poisson GLM summary \n")
print(summary(glm_pois))

# 3 Bagged regression-tree ensemble (manual random forest) 
bag_predict <- function(train, newdata, n_trees = 300, seed_offset = 0) {
  preds <- matrix(NA_real_, nrow = nrow(newdata), ncol = n_trees)
  for (b in seq_len(n_trees)) {
    set.seed(2000 * seed_offset + b)
    idx <- sample(seq_len(nrow(train)), size = nrow(train), replace = TRUE)
    boot_data <- train[idx, ]
    tree <- rpart(deaths_respiratory ~ time_index + quarter + log_pop + circ_rate,
                  data = boot_data,
                  control = rpart.control(cp = 0.01, minsplit = 6, maxdepth = 4))
    preds[, b] <- predict(tree, newdata = newdata)
  }
  preds
}

n_trees <- 500
bag_preds_full <- bag_predict(mort, mort, n_trees = n_trees, seed_offset = 0)
mort$bag_fit <- rowMeans(bag_preds_full)
mort$bag_lwr <- apply(bag_preds_full, 1, quantile, probs = 0.025)
mort$bag_upr <- apply(bag_preds_full, 1, quantile, probs = 0.975)

pois_pred <- predict(glm_pois, type = "response", se.fit = TRUE)
mort$pois_fit <- pois_pred$fit
mort$pois_lwr <- exp(log(pois_pred$fit) - 1.96 * pois_pred$se.fit / pois_pred$fit)
mort$pois_upr <- exp(log(pois_pred$fit) + 1.96 * pois_pred$se.fit / pois_pred$fit)

qp_pred <- predict(glm_qp, type = "response", se.fit = TRUE)
mort$qp_fit <- qp_pred$fit
mort$qp_lwr <- exp(log(qp_pred$fit) - 1.96 * qp_pred$se.fit / qp_pred$fit)
mort$qp_upr <- exp(log(qp_pred$fit) + 1.96 * qp_pred$se.fit / qp_pred$fit)

# Leave-one-out cross-validation: accuracy + interval calibration 
loo_pois <- loo_qp <- loo_bag <- numeric(n)
loo_pois_w <- loo_qp_w <- loo_bag_w <- numeric(n)
loo_pois_cov <- loo_qp_cov <- loo_bag_cov <- logical(n)

for (i in seq_len(n)) {
  train <- mort[-i, ]
  test  <- mort[i, , drop = FALSE]
  
  gp <- glm(deaths_respiratory ~ time_index + factor(quarter) +
              offset(log_pop), family = poisson(link = "log"), data = train)
  pp <- predict(gp, newdata = test, type = "response", se.fit = TRUE)
  loo_pois[i] <- pp$fit
  pl <- exp(log(pp$fit) - 1.96 * pp$se.fit / pp$fit)
  pu <- exp(log(pp$fit) + 1.96 * pp$se.fit / pp$fit)
  loo_pois_w[i] <- pu - pl
  loo_pois_cov[i] <- test$deaths_respiratory >= pl & test$deaths_respiratory <= pu
  
  gq <- glm(deaths_respiratory ~ time_index + factor(quarter) +
              offset(log_pop), family = quasipoisson(link = "log"), data = train)
  qp <- predict(gq, newdata = test, type = "response", se.fit = TRUE)
  loo_qp[i] <- qp$fit
  ql <- exp(log(qp$fit) - 1.96 * qp$se.fit / qp$fit)
  qu <- exp(log(qp$fit) + 1.96 * qp$se.fit / qp$fit)
  loo_qp_w[i] <- qu - ql
  loo_qp_cov[i] <- test$deaths_respiratory >= ql & test$deaths_respiratory <= qu
  
  bp <- bag_predict(train, test, n_trees = 200, seed_offset = i)
  loo_bag[i] <- mean(bp)
  bl <- quantile(bp, 0.025); bu <- quantile(bp, 0.975)
  loo_bag_w[i] <- bu - bl
  loo_bag_cov[i] <- test$deaths_respiratory >= bl & test$deaths_respiratory <= bu
}

results <- data.frame(
  Model = c("Poisson GLM", "Quasi-Poisson GLM", "Bagged trees (AI ensemble)"),
  RMSE = c(sqrt(mean((mort$deaths_respiratory - loo_pois)^2)),
           sqrt(mean((mort$deaths_respiratory - loo_qp)^2)),
           sqrt(mean((mort$deaths_respiratory - loo_bag)^2))),
  MAE = c(mean(abs(mort$deaths_respiratory - loo_pois)),
          mean(abs(mort$deaths_respiratory - loo_qp)),
          mean(abs(mort$deaths_respiratory - loo_bag))),
  Coverage_95 = c(mean(loo_pois_cov), mean(loo_qp_cov), mean(loo_bag_cov)),
  Mean_Interval_Width = c(mean(loo_pois_w), mean(loo_qp_w), mean(loo_bag_w))
)

cat("\n Leave-one-out cross-validation results (under-5 respiratory deaths) \n")
print(results, digits = 4)

write.csv(results, "loo_results_table.csv", row.names = FALSE)
saveRDS(list(mort = mort, results = results,
             loo_pois = loo_pois, loo_qp = loo_qp, loo_bag = loo_bag,
             loo_pois_cov = loo_pois_cov, loo_qp_cov = loo_qp_cov, loo_bag_cov = loo_bag_cov,
             glm_pois = glm_pois, glm_qp = glm_qp, dispersion = dispersion),
        "model_results.rds")
write.csv(mort, "mort_with_predictions.csv", row.names = FALSE)
cat("\nSaved loo_results_table.csv, model_results.rds, mort_with_predictions.csv\n")

# 03. XAI Interpretability
# Variable importance (permutation-free, RSS-reduction) and
# partial dependence profiles for the bagged-tree ensemble.
# Uses only rpart and base R — no external XAI library required.

cat("\n XAI: Ensemble Interpretability \n")

# 3a. Variable importance
# For each bootstrap tree, rpart stores $variable.importance: the total
# reduction in RSS attributed to each predictor across all splits.
# We average these scores over all B trees and normalise to 100%.

vi_predictors <- c("time_index", "quarter", "log_pop", "circ_rate")
vi_mat <- matrix(0, nrow = n_trees, ncol = length(vi_predictors),
                 dimnames = list(NULL, vi_predictors))

set.seed(42)
for (b in seq_len(n_trees)) {
  set.seed(2000 * 0 + b)
  idx <- sample(seq_len(n), size = n, replace = TRUE)
  tree <- rpart(deaths_respiratory ~ time_index + quarter + log_pop + circ_rate,
                data = mort[idx, ],
                control = rpart.control(cp = 0.01, minsplit = 6, maxdepth = 4))
  vi <- tree$variable.importance
  for (nm in names(vi)) {
    if (nm %in% vi_predictors) vi_mat[b, nm] <- vi[nm]
  }
}

mean_vi   <- colMeans(vi_mat)                   # mean RSS reduction per tree
norm_vi   <- mean_vi / sum(mean_vi) * 100       # normalised to 100 %

vi_df <- data.frame(
  Predictor           = c("Time index (t)", "Log population (log N_t)",
                          "Circulatory rate", "Calendar quarter (Q_t)"),
  Variable            = vi_predictors,
  Mean_RSS_reduction  = round(mean_vi, 1),
  Normalized_pct      = round(norm_vi, 1),
  stringsAsFactors    = FALSE
)
vi_df <- vi_df[order(-vi_df$Normalized_pct), ]

cat("Variable importance (permutation-free, normalised to 100%):\n")
print(vi_df[, c("Predictor", "Mean_RSS_reduction", "Normalized_pct")],
      row.names = FALSE, digits = 4)

write.csv(vi_df, "xai_variable_importance.csv", row.names = FALSE)

# 3b. Partial dependence profiles 
# PDP for time_index: set time_index to each observed value,
# marginalise over all other predictors (observed distribution),
# predict with the full 500-tree ensemble, take grand mean.

cat("\nComputing partial dependence profiles (this may take ~30 s) ...\n")

pdp_time <- numeric(n)
for (i in seq_len(n)) {
  tmp             <- mort
  tmp$time_index  <- mort$time_index[i]
  pdp_time[i]     <- mean(rowMeans(
    sapply(seq_len(n_trees), function(b) {
      set.seed(2000 * 0 + b)
      idx  <- sample(seq_len(n), n, replace = TRUE)
      tree <- rpart(deaths_respiratory ~ time_index + quarter +
                      log_pop + circ_rate,
                    data = mort[idx, ],
                    control = rpart.control(cp = 0.01, minsplit = 6,
                                            maxdepth = 4))
      predict(tree, newdata = tmp)
    })
  ))
}

# PDP for quarter: group mean ensemble prediction by quarter
# (quarter takes only 4 integer values, so a simple group mean
#  of the already-computed bag_preds_full suffices here)
pdp_quarter <- tapply(rowMeans(bag_preds_full), mort$quarter, mean)

cat("PDP time_index — mean predicted deaths:\n")
cat("  First 5 quarters (2004):", round(mean(pdp_time[1:5])), "\n")
cat("  Last  5 quarters (2014):", round(mean(pdp_time[40:44])), "\n")
cat("PDP calendar quarter — mean predicted deaths:\n")
for (q in 1:4) cat("  Q", q, "=", round(pdp_quarter[q]), "\n")

# Save PDP data
pdp_time_df    <- data.frame(time_index = mort$time_index,
                             date       = mort$date,
                             pdp_pred   = round(pdp_time, 1))
pdp_quarter_df <- data.frame(quarter  = 1:4,
                             pdp_pred = round(as.numeric(pdp_quarter), 1))
write.csv(pdp_time_df,    "xai_pdp_time.csv",    row.names = FALSE)
write.csv(pdp_quarter_df, "xai_pdp_quarter.csv", row.names = FALSE)

# 3c. XAI figures 

# Figure 5: Variable importance bar chart
png("fig5_variable_importance.png", width = 1400, height = 900, res = 150)
par(mar = c(6, 5, 3, 1))
bp_vi <- barplot(vi_df$Normalized_pct,
                 names.arg = rep("", nrow(vi_df)),
                 col = c("#1b4f72", "#2980b9", "#e67e22", "#c0392b"),
                 ylim = c(0, 50),
                 ylab = "Normalised variable importance (%)",
                 main = "Variable Importance: Bagged-Tree Ensemble\n(mean RSS reduction, averaged over 500 trees)",
                 cex.main = 0.95)
text(x = bp_vi,
     y = par("usr")[3] - 2,
     labels = vi_df$Predictor,
     xpd = NA, srt = 15, adj = c(1, 0.5), cex = 0.82)
text(bp_vi, vi_df$Normalized_pct + 1.5,
     labels = paste0(round(vi_df$Normalized_pct, 1), "%"), cex = 0.85)
dev.off()

# Figure 6: Partial dependence — time trend
png("fig6_pdp_time.png", width = 1400, height = 900, res = 150)
par(mar = c(4.5, 4.5, 3, 1))
plot(mort$date, pdp_time, type = "l", lwd = 2.5, col = "#1b4f72",
     xlab = "Quarter", ylab = "Partial dependence\n(mean predicted deaths, under-5)",
     main = "Partial Dependence on Time Trend\n(ensemble mean, all other predictors marginalised)",
     cex.main = 0.95, ylim = c(400, 900))
grid(col = "grey88")
dev.off()

# Figure 7: Partial dependence — calendar quarter
png("fig7_pdp_quarter.png", width = 1100, height = 900, res = 150)
par(mar = c(5, 5, 3, 1))
barplot(pdp_quarter, names.arg = paste0("Q", 1:4),
        col = "#2980b9", ylim = c(0, 850),
        ylab = "Partial dependence\n(mean predicted deaths, under-5)",
        main = "Partial Dependence on Calendar Quarter\n(bagged-tree ensemble)",
        cex.main = 0.95)
text(x = c(0.7, 1.9, 3.1, 4.3),
     y = pdp_quarter + 25,
     labels = round(pdp_quarter), cex = 0.88)
dev.off()

cat("XAI outputs saved: xai_variable_importance.csv, xai_pdp_time.csv,",
    "xai_pdp_quarter.csv\n")
cat("XAI figures: fig5_variable_importance.png, fig6_pdp_time.png,",
    "fig7_pdp_quarter.png\n")

# 04. Figures 

cat("\n Generating figures \n")

# Figure 1: Time series of under-5 respiratory deaths, Brazil 
png("fig1_time_series.png", width = 1400, height = 900, res = 150)
par(mar = c(4.5, 4.5, 2.5, 1))
plot(mort$date, mort$deaths_respiratory, type = "b", pch = 16, cex = 0.8,
     col = "#1b4f72", lwd = 1.5,
     xlab = "Quarter", ylab = "Respiratory deaths (under-5, Brazil)",
     main = "Quarterly Under-5 Respiratory Mortality, Brazil 2004-2014",
     cex.main = 1.0)
grid(col = "grey85")
dev.off()

# Figure 2: Model fit comparison (Poisson vs Bagged) with 95% intervals 
png("fig2_model_comparison.png", width = 1600, height = 1000, res = 150)
par(mar = c(4.5, 4.5, 2.5, 1))
ylim_range <- range(c(mort$deaths_respiratory, mort$pois_lwr, mort$pois_upr,
                      mort$bag_lwr, mort$bag_upr), na.rm = TRUE)
plot(mort$date, mort$deaths_respiratory, type = "p", pch = 16, cex = 0.9,
     col = "black", ylim = ylim_range,
     xlab = "Quarter", ylab = "Respiratory deaths (under-5)",
     main = "Observed Deaths vs. Model Fits with 95% Intervals",
     cex.main = 1.0)

lines(mort$date, mort$pois_fit, col = "#c0392b", lwd = 2)
polygon(c(mort$date, rev(mort$date)), c(mort$pois_lwr, rev(mort$pois_upr)),
        col = adjustcolor("#c0392b", alpha.f = 0.15), border = NA)

lines(mort$date, mort$bag_fit, col = "#1b4f72", lwd = 2, lty = 2)
polygon(c(mort$date, rev(mort$date)), c(mort$bag_lwr, rev(mort$bag_upr)),
        col = adjustcolor("#1b4f72", alpha.f = 0.15), border = NA)

points(mort$date, mort$deaths_respiratory, pch = 16, cex = 0.9, col = "black")

legend("topright", legend = c("Observed", "Poisson GLM fit (95% CI)",
                              "Bagged-tree ensemble fit (95% PI)"),
       col = c("black", "#c0392b", "#1b4f72"), lty = c(NA, 1, 2),
       pch = c(16, NA, NA), lwd = c(NA, 2, 2), bty = "n", cex = 0.85)
grid(col = "grey90")
dev.off()

# Figure 3: Coverage / calibration bar chart 
png("fig3_coverage.png", width = 1400, height = 1000, res = 150)
par(mar = c(9, 4.5, 2.5, 1))
short_names <- c("Poisson\nGLM", "Quasi-Poisson\nGLM", "Bagged trees\n(AI ensemble)")
bp <- barplot(results$Coverage_95 * 100, names.arg = rep("", 3),
              col = c("#c0392b", "#e67e22", "#1b4f72"),
              ylim = c(0, 100), ylab = "Empirical 95% interval coverage (%)",
              main = "Calibration of Uncertainty Intervals (Leave-One-Out CV)",
              cex.main = 1.0)
text(x = bp, y = par("usr")[3] - 5, labels = short_names, xpd = NA,
     srt = 0, adj = c(0.5, 1), cex = 0.85)
abline(h = 95, lty = 2, col = "grey30", lwd = 1.5)
text(bp, results$Coverage_95 * 100 + 4, labels = paste0(round(results$Coverage_95 * 100, 1), "%"),
     cex = 0.85)
text(max(bp) * 0.55, 97, "Nominal target (95%)", cex = 0.75, col = "grey30")
dev.off()

# Figure 4: Respiratory vs circulatory rates over time 
png("fig4_resp_circ_rates.png", width = 1400, height = 900, res = 150)
par(mar = c(4.5, 4.5, 2.5, 4.5))
plot(mort$date, mort$resp_rate, type = "l", col = "#1b4f72", lwd = 2,
     xlab = "Quarter", ylab = "Respiratory death rate (per 100,000, under-5)",
     main = "Respiratory vs. Circulatory Mortality Rates, Under-5, Brazil",
     cex.main = 1.0)
par(new = TRUE)
plot(mort$date, mort$circ_rate, type = "l", col = "#c0392b", lwd = 2, lty = 2,
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4)
mtext("Circulatory death rate (per 100,000, under-5)", side = 4, line = 3)
legend("topright", legend = c("Respiratory (J00-J99)", "Circulatory (I00-I99)"),
       col = c("#1b4f72", "#c0392b"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
dev.off()

cat("Figures written: fig1_time_series.png, fig2_model_comparison.png,",
    "fig3_coverage.png, fig4_resp_circ_rates.png\n")

cat("\n Pipeline complete \n")