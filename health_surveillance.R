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

keep <- c("age_group", "date", "pop",
          "deaths_respiratory", "deaths_circulatory",
          "J12_J18",
          "NoPCV")

mort_full <- df[, keep]
mort_full$year <- as.numeric(format(mort_full$date, "%Y"))
mort_full$quarter <- as.numeric(format(mort_full$date, "%m"))
mort_full$quarter <- ceiling(mort_full$quarter / 3)
mort_full$time_index <- as.numeric(factor(mort_full$date))

mort_full$resp_rate <- mort_full$deaths_respiratory / mort_full$pop * 1e5
mort_full$circ_rate <- mort_full$deaths_circulatory / mort_full$pop * 1e5

write.csv(mort_full, "brazil_child_mortality_resp_circ.csv", row.names = FALSE)

cat("\n Data preparation \n")
cat("Rows:", nrow(mort_full), "\n")
cat("Age groups:", paste(unique(mort_full$age_group), collapse = ", "), "\n")
cat("Date range:", as.character(min(mort_full$date)), "to",
    as.character(max(mort_full$date)), "\n\n")

# 02. Modeling
mort <- subset(mort_full, age_group == "3-<60mo")
mort$date <- as.Date(mort$date)
mort <- mort[order(mort$date), ]
mort$log_pop <- log(mort$pop)
n <- nrow(mort)
cat("\n Modeling \n")
cat("n (under-5, quarterly observations, Brazil 2004-2014):", n, "\n\n")

glm_pois <- glm(deaths_respiratory ~ time_index + factor(quarter) +
                  offset(log_pop),
                family = poisson(link = "log"), data = mort)

glm_qp <- glm(deaths_respiratory ~ time_index + factor(quarter) +
                offset(log_pop),
              family = quasipoisson(link = "log"), data = mort)

pearson_resid <- residuals(glm_pois, type = "pearson")
dispersion <- sum(pearson_resid^2) / glm_pois$df.residual
cat("Pearson dispersion statistic (Poisson model):", round(dispersion, 4), "\n\n")

# raw variance-to-mean ratio
vmr <- var(mort$deaths_respiratory) / mean(mort$deaths_respiratory)
cat("Raw variance-to-mean ratio (deaths_respiratory):", round(vmr, 3), "\n\n")

# ==================================================================
# NEW BLOCK 1: correlation between time_index and log_pop
# (the actual predictor variables, not the VI matrix)
# ==================================================================
cat("\n=== Correlation diagnostic: time_index vs log_pop (predictor variables) ===\n")
spearman_tp <- cor(mort$time_index, mort$log_pop, method = "spearman")
pearson_tp  <- cor(mort$time_index, mort$log_pop, method = "pearson")
rank_inversions <- sum(rank(mort$time_index) != rank(-mort$log_pop))

cat("Spearman correlation (time_index vs log_pop):", spearman_tp, "\n")
cat("Pearson correlation (time_index vs log_pop):", pearson_tp, "\n")
cat("Number of rank inversions:", rank_inversions, "\n\n")

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

# Leave-one-out cross-validation
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
print(results, digits = 6)

cat("\nFold coverage counts:\n")
cat("Poisson covered:", sum(loo_pois_cov), "/", n, "\n")
cat("Quasi-Poisson covered:", sum(loo_qp_cov), "/", n, "\n")
cat("Bagged trees covered:", sum(loo_bag_cov), "/", n, "\n")

# ==================================================================
# NEW BLOCK 2: exact (Clopper-Pearson) 95% CIs for LOOCV coverage,
# computed directly from the fold-level coverage indicator vectors
# ==================================================================
cat("\n=== Exact (Clopper-Pearson) 95% CIs for LOOCV coverage ===\n")

ci_pois <- binom.test(sum(loo_pois_cov), n)$conf.int
ci_qp   <- binom.test(sum(loo_qp_cov), n)$conf.int
ci_bag  <- binom.test(sum(loo_bag_cov), n)$conf.int

coverage_ci <- data.frame(
  Model = c("Poisson GLM", "Quasi-Poisson GLM", "Bagged trees (AI ensemble)"),
  Successes = c(sum(loo_pois_cov), sum(loo_qp_cov), sum(loo_bag_cov)),
  N = n,
  Coverage_pct = round(100 * c(mean(loo_pois_cov), mean(loo_qp_cov), mean(loo_bag_cov)), 1),
  CI_lower_pct = round(100 * c(ci_pois[1], ci_qp[1], ci_bag[1]), 1),
  CI_upper_pct = round(100 * c(ci_pois[2], ci_qp[2], ci_bag[2]), 1)
)
print(coverage_ci, row.names = FALSE)
write.csv(coverage_ci, "loocv_coverage_ci.csv", row.names = FALSE)

write.csv(results, "loo_results_table.csv", row.names = FALSE)
saveRDS(list(mort = mort, results = results,
             loo_pois = loo_pois, loo_qp = loo_qp, loo_bag = loo_bag,
             loo_pois_cov = loo_pois_cov, loo_qp_cov = loo_qp_cov, loo_bag_cov = loo_bag_cov,
             glm_pois = glm_pois, glm_qp = glm_qp, dispersion = dispersion),
        "model_results.rds")
write.csv(mort, "mort_with_predictions.csv", row.names = FALSE)
cat("\nSaved loo_results_table.csv, model_results.rds, mort_with_predictions.csv, loocv_coverage_ci.csv\n")

# ==================================================================
# NEW BLOCK 3: rolling-origin / expanding-window robustness check
# Minimum training window = 20 quarters, strictly past-only training.
# ==================================================================
cat("\n=== Rolling-origin / expanding-window robustness check ===\n")

min_train <- 20
test_idx <- (min_train + 1):n   # origins available to test
n_test <- length(test_idx)
cat("Minimum training window:", min_train, "quarters\n")
cat("Number of evaluable origins (n_test):", n_test, "\n\n")

ro_pois <- ro_qp <- ro_bag <- numeric(n_test)
ro_pois_w <- ro_qp_w <- ro_bag_w <- numeric(n_test)
ro_pois_cov <- ro_qp_cov <- ro_bag_cov <- logical(n_test)
ro_actual <- numeric(n_test)

for (k in seq_along(test_idx)) {
  t <- test_idx[k]
  train <- mort[1:(t - 1), ]      # strictly past observations only
  test  <- mort[t, , drop = FALSE]
  ro_actual[k] <- test$deaths_respiratory

  gp <- glm(deaths_respiratory ~ time_index + factor(quarter) +
              offset(log_pop), family = poisson(link = "log"), data = train)
  pp <- predict(gp, newdata = test, type = "response", se.fit = TRUE)
  ro_pois[k] <- pp$fit
  pl <- exp(log(pp$fit) - 1.96 * pp$se.fit / pp$fit)
  pu <- exp(log(pp$fit) + 1.96 * pp$se.fit / pp$fit)
  ro_pois_w[k] <- pu - pl
  ro_pois_cov[k] <- test$deaths_respiratory >= pl & test$deaths_respiratory <= pu

  gq <- glm(deaths_respiratory ~ time_index + factor(quarter) +
              offset(log_pop), family = quasipoisson(link = "log"), data = train)
  qp <- predict(gq, newdata = test, type = "response", se.fit = TRUE)
  ro_qp[k] <- qp$fit
  ql <- exp(log(qp$fit) - 1.96 * qp$se.fit / qp$fit)
  qu <- exp(log(qp$fit) + 1.96 * qp$se.fit / qp$fit)
  ro_qp_w[k] <- qu - ql
  ro_qp_cov[k] <- test$deaths_respiratory >= ql & test$deaths_respiratory <= qu

  bp <- bag_predict(train, test, n_trees = 200, seed_offset = 10000 + k)
  ro_bag[k] <- mean(bp)
  bl <- quantile(bp, 0.025); bu <- quantile(bp, 0.975)
  ro_bag_w[k] <- bu - bl
  ro_bag_cov[k] <- test$deaths_respiratory >= bl & test$deaths_respiratory <= bu
}

rolling_results <- data.frame(
  Model = c("Poisson GLM", "Quasi-Poisson GLM", "Bagged trees (AI ensemble)"),
  RMSE = c(sqrt(mean((ro_actual - ro_pois)^2)),
           sqrt(mean((ro_actual - ro_qp)^2)),
           sqrt(mean((ro_actual - ro_bag)^2))),
  MAE = c(mean(abs(ro_actual - ro_pois)),
          mean(abs(ro_actual - ro_qp)),
          mean(abs(ro_actual - ro_bag))),
  Coverage_pct = round(100 * c(mean(ro_pois_cov), mean(ro_qp_cov), mean(ro_bag_cov)), 1),
  Coverage_frac = c(sprintf("%d/%d", sum(ro_pois_cov), n_test),
                     sprintf("%d/%d", sum(ro_qp_cov), n_test),
                     sprintf("%d/%d", sum(ro_bag_cov), n_test)),
  Mean_Interval_Width = c(mean(ro_pois_w), mean(ro_qp_w), mean(ro_bag_w))
)

cat("Rolling-origin / expanding-window results:\n")
print(rolling_results, digits = 6, row.names = FALSE)
write.csv(rolling_results, "rolling_origin_results.csv", row.names = FALSE)

# 03. XAI Interpretability
cat("\n XAI: Ensemble Interpretability \n")

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

mean_vi   <- colMeans(vi_mat)
norm_vi   <- mean_vi / sum(mean_vi) * 100

predictor_labels <- c(
  time_index = "Time index (t)",
  quarter    = "Calendar quarter (Q_t)",
  log_pop    = "Log population (log N_t)",
  circ_rate  = "Circulatory rate"
)

vi_df <- data.frame(
  Predictor           = unname(predictor_labels[vi_predictors]),
  Variable            = vi_predictors,
  Mean_RSS_reduction  = round(mean_vi, 3),
  Normalized_pct      = round(norm_vi, 3),
  stringsAsFactors    = FALSE
)
vi_df <- vi_df[order(-vi_df$Normalized_pct), ]

cat("Variable importance (permutation-free, normalised to 100%):\n")
print(vi_df[, c("Predictor", "Mean_RSS_reduction", "Normalized_pct")],
      row.names = FALSE, digits = 8)

# Diagnostic: how many trees actually split on each predictor,
# and raw (non-rounded) means
cat("\nDiagnostic - raw (unrounded) mean VI per predictor:\n")
print(mean_vi, digits = 10)
cat("\nDiagnostic - number of trees (of", n_trees, ") with nonzero importance per predictor:\n")
print(colSums(vi_mat > 0))
cat("\nDiagnostic - are time_index and log_pop columns identical across all trees?\n")
print(identical(vi_mat[, "time_index"], vi_mat[, "log_pop"]))
cat("Number of trees where they differ:", sum(vi_mat[, "time_index"] != vi_mat[, "log_pop"]), "\n")
cat("Correlation between the two columns:", cor(vi_mat[, "time_index"], vi_mat[, "log_pop"]), "\n")

write.csv(vi_df, "xai_variable_importance.csv", row.names = FALSE)
write.csv(as.data.frame(vi_mat), "xai_vi_matrix_raw.csv", row.names = FALSE)

cat("\n Pipeline complete \n")
