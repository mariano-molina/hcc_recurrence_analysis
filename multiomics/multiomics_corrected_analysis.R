# ============================================================
# Liver tumour and adjacent tissue recurrence analysis
# Batch-corrected multiomics (mRNA + miRNA)
#
# Author: MAM
# Date: April 2026
#
# Description:
# This script reproduces the batch-corrected multiomics analyses
# for liver tumour, healthy/adjacent tissue, and all-sample settings.
# It integrates mRNA and miRNA data across KI and TCGA, applies
# ComBat-based batch correction, and performs:
#
# - preprocessing and metadata matching
# - mRNA + miRNA multi-omics integration
# - t-SNE visualization after batch correction
# - repeated within-cohort RF/PLS-DA cross-validation
# - cross-cohort transfer analyses
# - DE-informed recurrence score analyses
# - clinical association analyses
# - multivariable logistic regression and forest plots
#
# Input:
# - data/mRNA_merge_csv.csv
# - data/merge_mirna.csv
# - data/mRNA_miRNA_patient_sample_info.csv
# - data/SurroundingCategories_260408.xlsx
#
# Output:
# - Figures (PDF) in results/figures/
# - Tables (.csv) in results/
# - Processed objects (.rds) in results/objects/
#
# Notes:
# - One predefined tumour outlier sample is excluded:
#   X27.A033.303.Tumor
# - mRNA preprocessing follows the submission-ready mRNA logic
# - miRNA preprocessing follows the submission-ready miRNA logic
# - Batch variable is cohort/source (KI vs TCGA)
# - ComBat protects recurrence status during correction
# ============================================================

# ============================================================
# 0. Setup
# ============================================================

suppressPackageStartupMessages({
  library(sva)
  library(Rtsne)
  library(FactoMineR)
  library(ggplot2)
  library(patchwork)
  library(mixOmics)
  library(randomForest)
  library(pROC)
  library(scales)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(limma)
  library(readxl)
})

set.seed(999)

# Project directories
input_dir   <- "data"
output_dir  <- "results"
figure_dir  <- file.path(output_dir, "figures")
object_dir  <- file.path(output_dir, "objects")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(object_dir, showWarnings = FALSE, recursive = TRUE)

# Input files
mrna_file    <- file.path(input_dir, "mRNA_merge_csv.csv")
mirna_file   <- file.path(input_dir, "merge_mirna.csv")
metadata_file <- file.path(input_dir, "mRNA_miRNA_patient_sample_info.csv")
viral_file   <- file.path(input_dir, "SurroundingCategories_260408.xlsx")

# Predefined QC exclusion
bad_sample <- "X27.A033.303.Tumor"

# ============================================================
# 1. Helper functions
# ============================================================

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", sprintf("%.3f", p))
}

get_y_top <- function(x, expand = 0.06) {
  rng <- range(x, na.rm = TRUE)
  rng[2] + expand * diff(rng)
}

add_n_labels <- function(df, var) {
  tab <- table(df[[var]])
  labs <- paste0(names(tab), "\n(n=", as.integer(tab), ")")
  names(labs) <- names(tab)
  labs
}

normalize_ids <- function(x) {
  x <- trimws(x)
  x <- gsub("\\.", "_", x)
  x <- gsub("-", "_", x)
  x <- gsub("__", "_", x)
  x
}

remove_leading_X <- function(x) {
  gsub("^X", "", x)
}

fix_ki_mrna_ids <- function(x) {
  x <- gsub("N_RNA_", "", x)
  x <- gsub("T_RNA_", "", x)
  x <- gsub("([0-9]+)[NT]_([A-Za-z]+)$", "\\1_\\2", x)
  x <- gsub("([0-9]+)(Normal)$", "\\1_Normal", x)
  x <- gsub("([0-9]+)(Tumor)$", "\\1_Tumor", x)
  x
}

fix_tcga_ids <- function(x) {
  gsub("(TCGA_[A-Z0-9]+_[A-Z0-9]+)T(_[0-9][0-9]$)", "\\1\\2", x)
}

make_clean_id <- function(x) {
  x <- trimws(x)
  x <- gsub("_", ".", x)
  x <- gsub("-", ".", x)
  x
}

plot_roc_curve <- function(roc_obj, auc_val, ci_vals = NULL, title_text, color_value) {
  if (!is.null(ci_vals)) {
    label_text <- sprintf(
      "AUC = %.3f\n95%% CI: %.3f-%.3f",
      auc_val,
      ci_vals[1],
      ci_vals[3]
    )
  } else {
    label_text <- sprintf("AUC = %.3f", auc_val)
  }
  
  ggroc(roc_obj, legacy.axes = TRUE, color = color_value, linewidth = 1.2) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_minimal(base_size = 14) +
    labs(
      title = title_text,
      x = "Specificity",
      y = "Sensitivity"
    ) +
    annotate(
      "label",
      x = 0.78,
      y = 0.18,
      label = label_text,
      size = 4.5,
      color = "black",
      fill = alpha("white", 0.85),
      label.r = unit(0.1, "lines"),
      hjust = 1
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

plot_auc_distribution <- function(auc_values, title_text, fill_color) {
  auc_df <- data.frame(AUC = auc_values)
  
  ggplot(auc_df, aes(x = "", y = AUC)) +
    geom_boxplot(width = 0.25, fill = fill_color, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.08, size = 2, alpha = 0.6) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal(base_size = 14) +
    labs(
      title = title_text,
      x = "",
      y = "AUC"
    ) +
    annotate(
      "label",
      x = 1,
      y = 0.95,
      label = paste0(
        "Mean AUC = ", round(mean(auc_values), 3),
        "\nSD = ", round(sd(auc_values), 3),
        "\nRepeats = ", length(auc_values)
      ),
      size = 4,
      fill = alpha("white", 0.7)
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

plot_mdg <- function(importance_df, title_text, fill_color, top_n = 10, feature_col = "Feature") {
  plot_df <- importance_df[order(-importance_df$MeanDecreaseGini), , drop = FALSE]
  plot_df <- plot_df[1:min(top_n, nrow(plot_df)), , drop = FALSE]
  plot_df[[feature_col]] <- factor(plot_df[[feature_col]], levels = rev(plot_df[[feature_col]]))
  
  ggplot(plot_df, aes_string(x = feature_col, y = "MeanDecreaseGini")) +
    geom_col(fill = fill_color) +
    coord_flip() +
    theme_minimal(base_size = 14) +
    labs(
      title = title_text,
      x = feature_col,
      y = "Mean Decrease Gini"
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

run_plsda_rf_cv_repeated <- function(X, y, n_repeats = 20, k = 5,
                                     ncomp = 2, n_top = 100, ntree = 500,
                                     seed = 999) {
  set.seed(seed)
  
  X <- as.matrix(X)
  y <- factor(y, levels = c("Remission", "Recidivism"))
  
  stopifnot(nrow(X) == length(y))
  stopifnot(!is.null(rownames(X)))
  
  auc_values <- c()
  rf_importance_list <- list()
  last_all_probs <- NULL
  last_all_true  <- NULL
  
  for (r in seq_len(n_repeats)) {
    fold_assignments <- rep(NA_integer_, length(y))
    names(fold_assignments) <- rownames(X)
    
    for (cls in levels(y)) {
      cls_ids <- sample(rownames(X)[y == cls])
      fold_ids <- rep(seq_len(k), length.out = length(cls_ids))
      fold_assignments[cls_ids] <- fold_ids
    }
    
    all_probs <- c()
    all_true <- c()
    
    for (i in seq_len(k)) {
      test_ids  <- names(fold_assignments)[fold_assignments == i]
      train_ids <- setdiff(rownames(X), test_ids)
      
      X_train_full <- X[train_ids, , drop = FALSE]
      X_test_full  <- X[test_ids, , drop = FALSE]
      y_train <- factor(y[train_ids], levels = c("Remission", "Recidivism"))
      y_test  <- factor(y[test_ids], levels = c("Remission", "Recidivism"))
      
      if (length(unique(y_train)) < 2) next
      
      nzv <- apply(X_train_full, 2, sd, na.rm = TRUE) > 0
      X_train_plsda <- X_train_full[, nzv, drop = FALSE]
      X_test_plsda  <- X_test_full[, nzv, drop = FALSE]
      
      if (ncol(X_train_plsda) < 2) next
      
      colnames(X_train_plsda) <- make.unique(colnames(X_train_plsda))
      colnames(X_test_plsda)  <- colnames(X_train_plsda)
      
      plsda_fit <- plsda(X_train_plsda, y_train, ncomp = ncomp)
      load1 <- plsda_fit$loadings$X[, 1]
      load2 <- plsda_fit$loadings$X[, 2]
      loading_score <- abs(load1) + abs(load2)
      
      ranked_features <- names(sort(loading_score, decreasing = TRUE))
      selected_features <- ranked_features[1:min(n_top, length(ranked_features))]
      selected_features <- intersect(selected_features, colnames(X_train_plsda))
      selected_features <- intersect(selected_features, colnames(X_test_plsda))
      
      if (length(selected_features) < 2) next
      
      X_train <- X_train_plsda[, selected_features, drop = FALSE]
      X_test  <- X_test_plsda[, selected_features, drop = FALSE]
      
      rf_fit <- randomForest(
        x = X_train,
        y = y_train,
        ntree = ntree,
        importance = TRUE,
        strata = y_train,
        sampsize = rep(min(table(y_train)), 2)
      )
      
      rf_importance_list[[length(rf_importance_list) + 1]] <- data.frame(
        Feature = rownames(rf_fit$importance),
        MeanDecreaseGini = rf_fit$importance[, "MeanDecreaseGini"],
        Repeat = r,
        Fold = i,
        stringsAsFactors = FALSE
      )
      
      y_prob <- predict(rf_fit, newdata = X_test, type = "prob")[, "Recidivism"]
      all_probs <- c(all_probs, y_prob)
      all_true  <- c(all_true, as.character(y_test))
    }
    
    if (length(unique(all_true)) < 2) next
    
    roc_obj_repeat <- roc(
      response = all_true,
      predictor = all_probs,
      levels = c("Remission", "Recidivism"),
      direction = "<"
    )
    
    auc_values <- c(auc_values, as.numeric(auc(roc_obj_repeat)))
    last_all_probs <- all_probs
    last_all_true  <- all_true
  }
  
  if (length(auc_values) == 0) stop("No repeated CV AUC values were generated.")
  if (length(rf_importance_list) == 0) stop("No RF importance values were collected across repeats/folds.")
  
  importance_df <- do.call(rbind, rf_importance_list)
  importance_summary <- aggregate(MeanDecreaseGini ~ Feature, data = importance_df, FUN = mean)
  importance_summary <- importance_summary[order(-importance_summary$MeanDecreaseGini), ]
  
  roc_last <- roc(
    response = factor(last_all_true, levels = c("Remission", "Recidivism")),
    predictor = last_all_probs,
    levels = c("Remission", "Recidivism"),
    direction = "<"
  )
  
  list(
    auc_values = auc_values,
    mean_auc = mean(auc_values),
    sd_auc = sd(auc_values),
    roc_last = roc_last,
    importance = importance_summary,
    ranked_features = importance_summary$Feature
  )
}

run_plsda_rf_transfer <- function(X_train_full, y_train, X_test_full, y_test,
                                  ncomp = 2, n_plsda = 100, n_rf = 100,
                                  ntree_rank = 1000, ntree_final = 1000,
                                  seed = 999) {
  set.seed(seed)
  
  X_train_full <- as.matrix(X_train_full)
  X_test_full  <- as.matrix(X_test_full)
  y_train <- factor(y_train, levels = c("Remission", "Recidivism"))
  y_test  <- factor(y_test, levels = c("Remission", "Recidivism"))
  
  nzv <- apply(X_train_full, 2, sd, na.rm = TRUE) > 0
  X_train_plsda <- X_train_full[, nzv, drop = FALSE]
  X_test_plsda  <- X_test_full[, nzv, drop = FALSE]
  
  if (ncol(X_train_plsda) < 2) stop("Too few non-zero variance features available.")
  
  colnames(X_train_plsda) <- make.unique(colnames(X_train_plsda))
  colnames(X_test_plsda)  <- colnames(X_train_plsda)
  
  plsda_fit <- plsda(X_train_plsda, y_train, ncomp = ncomp)
  
  load1 <- plsda_fit$loadings$X[, 1]
  load2 <- plsda_fit$loadings$X[, 2]
  
  top_c1 <- names(sort(abs(load1), decreasing = TRUE)[1:min(n_plsda, length(load1))])
  top_c2 <- names(sort(abs(load2), decreasing = TRUE)[1:min(n_plsda, length(load2))])
  top_plsda_features <- union(top_c1, top_c2)
  
  X_train_rank <- X_train_plsda[, top_plsda_features, drop = FALSE]
  X_test_rank  <- X_test_plsda[, top_plsda_features, drop = FALSE]
  
  rf_rank <- randomForest(
    x = X_train_rank,
    y = y_train,
    importance = TRUE,
    ntree = ntree_rank
  )
  
  importance_df <- data.frame(
    Feature = rownames(rf_rank$importance),
    MeanDecreaseGini = rf_rank$importance[, "MeanDecreaseGini"],
    stringsAsFactors = FALSE
  )
  importance_df <- importance_df[order(-importance_df$MeanDecreaseGini), ]
  ranked_features <- importance_df$Feature
  
  selected_features <- ranked_features[1:min(n_rf, length(ranked_features))]
  selected_features <- intersect(selected_features, colnames(X_train_rank))
  selected_features <- intersect(selected_features, colnames(X_test_rank))
  
  if (length(selected_features) < 2) stop("Too few selected features remained.")
  
  X_train <- X_train_rank[, selected_features, drop = FALSE]
  X_test  <- X_test_rank[, selected_features, drop = FALSE]
  
  rf_final <- randomForest(
    x = X_train,
    y = y_train,
    ntree = ntree_final,
    importance = TRUE
  )
  
  y_pred <- predict(rf_final, X_test)
  y_prob <- predict(rf_final, X_test, type = "prob")[, "Recidivism"]
  
  roc_obj <- roc(
    response = y_test,
    predictor = y_prob,
    levels = c("Remission", "Recidivism"),
    direction = "<"
  )
  
  list(
    roc = roc_obj,
    auc = as.numeric(auc(roc_obj)),
    ci = ci.auc(roc_obj),
    accuracy = mean(y_pred == y_test),
    confusion_matrix = table(True = y_test, Predicted = y_pred),
    importance = importance_df,
    ranked_features = ranked_features,
    selected_features = selected_features,
    plsda = plsda_fit
  )
}

run_de_signature_transfer <- function(expr_train, y_train, expr_test, y_test,
                                      train_name, test_name, fdr_cut = 0.25) {
  y_train <- factor(y_train, levels = c("Remission", "Recidivism"))
  y_test  <- factor(y_test, levels = c("Remission", "Recidivism"))
  
  design_train <- model.matrix(~ y_train)
  fit_train <- eBayes(lmFit(t(expr_train), design_train))
  res_train <- topTable(fit_train, coef = 2, number = Inf)
  res_train$Feature <- rownames(res_train)
  
  design_test <- model.matrix(~ y_test)
  fit_test <- eBayes(lmFit(t(expr_test), design_test))
  res_test <- topTable(fit_test, coef = 2, number = Inf)
  res_test$Feature <- rownames(res_test)
  
  merged <- merge(
    res_train[, c("Feature", "logFC", "adj.P.Val")],
    res_test[, c("Feature", "logFC", "adj.P.Val")],
    by = "Feature",
    suffixes = c(paste0("_", train_name), paste0("_", test_name))
  )
  
  fc_train <- paste0("logFC_", train_name)
  fc_test  <- paste0("logFC_", test_name)
  fdr_train <- paste0("adj.P.Val_", train_name)
  fdr_test  <- paste0("adj.P.Val_", test_name)
  
  merged$SameDir <- sign(merged[[fc_train]]) == sign(merged[[fc_test]])
  
  signature_features <- merged[
    merged$SameDir &
      merged[[fdr_train]] < fdr_cut &
      merged[[fdr_test]] < fdr_cut,
    "Feature"
  ]
  
  up_features <- merged$Feature[merged[[fc_train]] > 0 & merged$SameDir]
  down_features <- merged$Feature[merged[[fc_train]] < 0 & merged$SameDir]
  
  up_features <- intersect(up_features, colnames(expr_train))
  up_features <- intersect(up_features, colnames(expr_test))
  down_features <- intersect(down_features, colnames(expr_train))
  down_features <- intersect(down_features, colnames(expr_test))
  
  score_fun <- function(mat, up, down) {
    up_score <- if (length(up) > 0) rowMeans(mat[, up, drop = FALSE], na.rm = TRUE) else rep(0, nrow(mat))
    down_score <- if (length(down) > 0) rowMeans(mat[, down, drop = FALSE], na.rm = TRUE) else rep(0, nrow(mat))
    up_score - down_score
  }
  
  score_train <- score_fun(expr_train, up_features, down_features)
  score_test  <- score_fun(expr_test, up_features, down_features)
  
  roc_train <- roc(response = y_train, predictor = score_train, levels = c("Remission", "Recidivism"))
  roc_test  <- roc(response = y_test, predictor = score_test, levels = c("Remission", "Recidivism"))
  
  list(
    merged = merged,
    signature_features = signature_features,
    up_features = up_features,
    down_features = down_features,
    score_train = score_train,
    score_test = score_test,
    roc_train = roc_train,
    roc_test = roc_test,
    auc_train = as.numeric(auc(roc_train)),
    auc_test = as.numeric(auc(roc_test))
  )
}

run_rf_cv_signature <- function(X, y, ntree = 500, k = 5, seed = 999) {
  set.seed(seed)
  
  X <- as.matrix(X)
  y <- factor(y, levels = c("Remission", "Recidivism"))
  
  folds <- rep(NA_integer_, length(y))
  names(folds) <- rownames(X)
  
  for (cls in levels(y)) {
    ids <- rownames(X)[y == cls]
    ids <- sample(ids)
    folds[ids] <- rep(seq_len(k), length.out = length(ids))
  }
  
  all_probs <- c()
  all_true <- c()
  
  for (f in seq_len(k)) {
    test_ids  <- names(folds)[folds == f]
    train_ids <- setdiff(rownames(X), test_ids)
    
    X_train <- X[train_ids, , drop = FALSE]
    X_test  <- X[test_ids, , drop = FALSE]
    y_train <- y[train_ids]
    y_test  <- y[test_ids]
    
    if (length(unique(y_train)) < 2) next
    
    rf_fit <- randomForest(
      x = X_train,
      y = y_train,
      ntree = ntree,
      importance = TRUE
    )
    
    y_prob <- predict(rf_fit, X_test, type = "prob")[, "Recidivism"]
    all_probs <- c(all_probs, y_prob)
    all_true <- c(all_true, as.character(y_test))
  }
  
  roc_obj <- roc(
    response = all_true,
    predictor = all_probs,
    levels = c("Remission", "Recidivism"),
    direction = "<"
  )
  
  list(
    roc = roc_obj,
    auc = as.numeric(auc(roc_obj))
  )
}

make_forest_df <- function(model, model_type = c("full", "viral")) {
  model_type <- match.arg(model_type)
  
  coef_table <- summary(model)$coefficients
  ci <- suppressMessages(confint(model))
  
  forest_df <- data.frame(
    Variable = rownames(coef_table),
    Estimate = coef_table[, "Estimate"],
    SE = coef_table[, "Std. Error"],
    P_value = coef_table[, "Pr(>|z|)"],
    CI_low = ci[, 1],
    CI_high = ci[, 2],
    stringsAsFactors = FALSE
  )
  
  forest_df <- forest_df %>%
    mutate(
      OR = exp(Estimate),
      OR_low = exp(CI_low),
      OR_high = exp(CI_high)
    ) %>%
    filter(Variable != "(Intercept)")
  
  if (model_type == "full") {
    forest_df$Variable <- recode(
      forest_df$Variable,
      "Score_z" = "Adjacent score",
      "CirrosisYes" = "Cirrhosis",
      "MilanIn" = "Milan criteria",
      "PS1" = "PS 1",
      "PS2" = "PS 2",
      "PS3" = "PS 3"
    )
    
    forest_df$Variable <- factor(
      forest_df$Variable,
      levels = rev(c("Adjacent score", "Cirrhosis", "Milan criteria", "PS 1", "PS 2", "PS 3"))
    )
  }
  
  if (model_type == "viral") {
    forest_df$Variable <- recode(
      forest_df$Variable,
      "Score_z" = "Adjacent score",
      "CirrosisYes" = "Cirrhosis",
      "MilanIn" = "Milan criteria",
      "PS1" = "PS 1",
      "PS2" = "PS 2",
      "PS3" = "PS 3",
      "ViralStatusViral" = "Viral etiology"
    )
    
    forest_df$Variable <- factor(
      forest_df$Variable,
      levels = rev(c("Adjacent score", "Cirrhosis", "Milan criteria", "PS 1", "PS 2", "PS 3", "Viral etiology"))
    )
  }
  
  forest_df %>%
    mutate(Significant = P_value < 0.05)
}

plot_forest <- function(forest_df, title_text) {
  ggplot(forest_df, aes(x = OR, y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_errorbarh(
      aes(xmin = OR_low, xmax = OR_high, color = Significant),
      height = 0.16,
      linewidth = 0.8
    ) +
    geom_point(aes(color = Significant), size = 3) +
    scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "black"), guide = "none") +
    scale_x_log10(breaks = c(0.3, 1, 3, 10, 30)) +
    coord_cartesian(xlim = c(0.3, 30)) +
    labs(
      title = title_text,
      x = "Odds ratio (log scale)",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(face = "bold", color = "black", size = 11)
    )
}

# ============================================================
# 2. Load and preprocess mRNA data
# ============================================================

message("Loading mRNA matrix and metadata...")
mrna_raw <- read.csv(mrna_file, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)

mrna_clean <- mrna_raw[!(is.na(mrna_raw[[1]]) | mrna_raw[[1]] == ""), ]
mrna_clean <- mrna_clean[rowSums(mrna_clean[, -1] != 0) > 0, ]

mrna_sample_ids <- colnames(mrna_clean)[-1]
sample_id_column <- "mRNA"

mrna_sample_ids_clean <- trimws(mrna_sample_ids)
metadata_ids_clean <- trimws(as.character(metadata[[sample_id_column]]))

mrna_sample_ids_clean <- gsub("-", ".", mrna_sample_ids_clean)
metadata_ids_clean <- gsub("-", ".", metadata_ids_clean)

mrna_sample_ids_match <- ifelse(
  grepl("^TCGA", mrna_sample_ids_clean),
  mrna_sample_ids_clean,
  paste0("X", mrna_sample_ids_clean)
)

recurrence_lookup <- setNames(metadata$Recidiv, metadata_ids_clean)
tissue_lookup     <- setNames(metadata$Tumour, metadata_ids_clean)

mrna_recurrence_labels <- sapply(mrna_sample_ids_match, function(id) {
  value <- recurrence_lookup[id]
  if (is.na(value)) return(NA_character_)
  if (value == 1) return("Recidivism")
  if (value == 0) return("Remission")
  NA_character_
})

mrna_tissue_labels <- sapply(mrna_sample_ids_match, function(id) {
  value <- tissue_lookup[id]
  if (is.na(value)) return(NA_character_)
  if (value == 1) return("Tumour")
  if (value == 0) return("Healthy")
  NA_character_
})

mrna_expr <- apply(mrna_clean[, -1], 2, as.numeric)
mrna_expr <- t(mrna_expr)
rownames(mrna_expr) <- mrna_sample_ids_clean
colnames(mrna_expr) <- make.unique(mrna_clean[[1]])

valid_mrna <- !is.na(mrna_recurrence_labels) & !is.na(mrna_tissue_labels)

mrna_expr <- mrna_expr[valid_mrna, , drop = FALSE]
mrna_recurrence_labels <- factor(mrna_recurrence_labels[valid_mrna], levels = c("Remission", "Recidivism"))
mrna_tissue_labels     <- factor(mrna_tissue_labels[valid_mrna], levels = c("Healthy", "Tumour"))

mrna_expr <- mrna_expr[, apply(mrna_expr, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
mode(mrna_expr) <- "numeric"
mrna_expr <- log2(mrna_expr + 1)

names(mrna_recurrence_labels) <- rownames(mrna_expr)
names(mrna_tissue_labels)     <- rownames(mrna_expr)

bad_sample_mrna <- gsub("-", ".", bad_sample)
if (bad_sample_mrna %in% rownames(mrna_expr)) {
  mrna_expr <- mrna_expr[rownames(mrna_expr) != bad_sample_mrna, , drop = FALSE]
  mrna_recurrence_labels <- mrna_recurrence_labels[names(mrna_recurrence_labels) != bad_sample_mrna]
  mrna_tissue_labels     <- mrna_tissue_labels[names(mrna_tissue_labels) != bad_sample_mrna]
}

# Harmonize mRNA rownames to integration format
mrna_ids <- rownames(mrna_expr)
mrna_ids <- normalize_ids(mrna_ids)
mrna_ids <- remove_leading_X(mrna_ids)
mrna_ids <- fix_ki_mrna_ids(mrna_ids)
mrna_ids <- fix_tcga_ids(mrna_ids)
rownames(mrna_expr) <- mrna_ids

names(mrna_recurrence_labels) <- rownames(mrna_expr)
names(mrna_tissue_labels)     <- rownames(mrna_expr)

# ============================================================
# 3. Load and preprocess miRNA data
# ============================================================

message("Loading miRNA matrix...")

mirna_raw <- read.csv(mirna_file, stringsAsFactors = FALSE)

mirna_clean <- mirna_raw[!(is.na(mirna_raw[[1]]) | mirna_raw[[1]] == ""), ]
mirna_clean <- mirna_clean[rowSums(mirna_clean[, -1] != 0) > 0, ]

mirna_names <- mirna_clean[[1]]
mirna_expr <- apply(mirna_clean[, -1], 2, as.numeric)
rownames(mirna_expr) <- mirna_names
mirna_expr <- t(mirna_expr)

mirna_ids <- rownames(mirna_expr)
mirna_ids <- normalize_ids(mirna_ids)
mirna_ids <- remove_leading_X(mirna_ids)
mirna_ids <- fix_ki_mrna_ids(mirna_ids)
mirna_ids <- fix_tcga_ids(mirna_ids)
rownames(mirna_expr) <- mirna_ids

mirna_expr <- mirna_expr[, apply(mirna_expr, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
colnames(mirna_expr) <- make.unique(colnames(mirna_expr))
mode(mirna_expr) <- "numeric"
mirna_expr <- log2(mirna_expr + 1)

# ============================================================
# 4. Merge mRNA and miRNA into batch-corrected multi-omics
# ============================================================

common_samples <- intersect(rownames(mrna_expr), rownames(mirna_expr))

cat("Matched multi-omics samples:", length(common_samples), "\n")
cat("TCGA matched:", sum(grepl("^TCGA", common_samples)), "\n")
cat("KI matched:", sum(!grepl("^TCGA", common_samples)), "\n")

message("Matched multi-omics samples: ", length(common_samples))
if (length(common_samples) == 0) {
  stop("No overlapping samples remained between mRNA and miRNA matrices.")
}

mrna_multi <- mrna_expr[common_samples, , drop = FALSE]
mirna_multi <- mirna_expr[common_samples, , drop = FALSE]

multi_omics_matrix <- cbind(mrna_multi, mirna_multi)

recidiv_multi <- mrna_recurrence_labels[common_samples]
tumour_multi  <- mrna_tissue_labels[common_samples]

source_multi <- factor(
  ifelse(grepl("^TCGA", common_samples), "TCGA", "KI"),
  levels = c("KI", "TCGA")
)
names(source_multi) <- common_samples

stopifnot(nrow(multi_omics_matrix) == length(recidiv_multi))
stopifnot(nrow(multi_omics_matrix) == length(tumour_multi))

# Batch correction
message("Running ComBat batch correction...")
combat_input <- t(multi_omics_matrix)
batch <- droplevels(factor(source_multi))
recidiv_multi <- droplevels(factor(recidiv_multi, levels = c("Remission", "Recidivism")))

cat("Batch distribution:\n")
print(table(batch, useNA = "ifany"))

cat("\nOutcome distribution:\n")
print(table(recidiv_multi, useNA = "ifany"))

if (nlevels(batch) < 2) {
  stop("Batch correction requires at least two batches. Only one cohort is present after sample matching.")
}

if (nlevels(recidiv_multi) >= 2) {
  mod <- model.matrix(~ recidiv_multi)
  expr_combat <- ComBat(
    dat = combat_input,
    batch = batch,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
} else {
  warning("Only one recurrence class remained after matching. Running ComBat without biological covariate protection.")
  expr_combat <- ComBat(
    dat = combat_input,
    batch = batch,
    par.prior = TRUE,
    prior.plots = FALSE
  )
}

multi_omics_bc <- t(expr_combat)
multi_omics_corrected <- multi_omics_bc

saveRDS(multi_omics_matrix, file.path(object_dir, "multi_omics_matrix_raw.rds"))
saveRDS(multi_omics_bc, file.path(object_dir, "multi_omics_matrix_batch_corrected.rds"))
saveRDS(recidiv_multi, file.path(object_dir, "multi_omics_recurrence_labels.rds"))
saveRDS(tumour_multi, file.path(object_dir, "multi_omics_tissue_labels.rds"))
saveRDS(source_multi, file.path(object_dir, "multi_omics_source_labels.rds"))

sample_summary <- data.frame(
  SampleID = rownames(multi_omics_bc),
  Cohort = source_multi,
  Tissue = tumour_multi,
  Outcome = recidiv_multi,
  stringsAsFactors = FALSE
)
write.csv(sample_summary, file.path(output_dir, "multiomics_sample_summary.csv"), row.names = FALSE)

# ============================================================
# 5. t-SNE on batch-corrected multi-omics
# ============================================================

message("Running PCA + t-SNE...")
set.seed(42)

pca_res <- PCA(multi_omics_bc, ncp = min(100, nrow(multi_omics_bc) - 1), graph = FALSE)
pca_features <- pca_res$ind$coord

tsne_result <- Rtsne(
  pca_features,
  perplexity = min(30, floor((nrow(pca_features) - 1) / 3)),
  verbose = TRUE,
  max_iter = 1500
)

tsne_df <- data.frame(
  TSNE1 = tsne_result$Y[, 1],
  TSNE2 = tsne_result$Y[, 2],
  Outcome = recidiv_multi,
  Status = tumour_multi,
  Source = source_multi,
  stringsAsFactors = FALSE
)

p_tsne_outcome <- ggplot(tsne_df, aes(TSNE1, TSNE2, fill = Outcome)) +
  geom_point(size = 2.5, shape = 21, color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("Recidivism" = "#E41A1C", "Remission" = "#377EB8")) +
  labs(title = "SF-7A", x = "t-SNE1", y = "t-SNE2") +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"))

p_tsne_tissue <- ggplot(tsne_df, aes(TSNE1, TSNE2, fill = Status)) +
  geom_point(size = 2.5, shape = 21, color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("Tumour" = "orange", "Healthy" = "#4DAF4A")) +
  labs(title = "Fig. 5A", x = "t-SNE1", y = "t-SNE2") +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"))

p_tsne_source <- ggplot(tsne_df, aes(TSNE1, TSNE2, fill = Source)) +
  geom_point(size = 2.5, shape = 21, color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("TCGA" = "#93C47D", "KI" = "#912F6C")) +
  labs(title = "SF-7B", x = "t-SNE1", y = "t-SNE2") +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"))

ggsave(file.path(figure_dir, "multiomics_tsne_outcome.pdf"), p_tsne_outcome, width = 5, height = 4)
ggsave(file.path(figure_dir, "multiomics_tsne_tissue.pdf"), p_tsne_tissue, width = 5, height = 4)
ggsave(file.path(figure_dir, "multiomics_tsne_source.pdf"), p_tsne_source, width = 5, height = 4)

# ============================================================
# 6. Define tumour, healthy, and all-sample matrices
# ============================================================

tumour_ids  <- names(tumour_multi)[tumour_multi == "Tumour"]
healthy_ids <- names(tumour_multi)[tumour_multi == "Healthy"]

expr_tumour <- multi_omics_bc[tumour_ids, , drop = FALSE]
expr_healthy <- multi_omics_bc[healthy_ids, , drop = FALSE]

y_tumour <- factor(recidiv_multi[tumour_ids], levels = c("Remission", "Recidivism"))
y_healthy <- factor(recidiv_multi[healthy_ids], levels = c("Remission", "Recidivism"))

src_tumour <- source_multi[tumour_ids]
src_healthy <- source_multi[healthy_ids]

# Tumour split
X_tcga_tumour <- expr_tumour[src_tumour == "TCGA", , drop = FALSE]
y_tcga_tumour <- y_tumour[src_tumour == "TCGA"]

X_ki_tumour <- expr_tumour[src_tumour == "KI", , drop = FALSE]
y_ki_tumour <- y_tumour[src_tumour == "KI"]

# Healthy split
X_tcga_healthy <- expr_healthy[src_healthy == "TCGA", , drop = FALSE]
y_tcga_healthy <- y_healthy[src_healthy == "TCGA"]

X_ki_healthy <- expr_healthy[src_healthy == "KI", , drop = FALSE]
y_ki_healthy <- y_healthy[src_healthy == "KI"]

# All split
X_tcga_all <- multi_omics_bc[source_multi == "TCGA", , drop = FALSE]
y_tcga_all <- factor(recidiv_multi[source_multi == "TCGA"], levels = c("Remission", "Recidivism"))

X_ki_all <- multi_omics_bc[source_multi == "KI", , drop = FALSE]
y_ki_all <- factor(recidiv_multi[source_multi == "KI"], levels = c("Remission", "Recidivism"))


# ============================================================
# 7. Within-cohort CV (tumour multi-omics, single CV + 95% CI)
# ============================================================

message("Running CV: tumour multi-omics...")

rf_tcga_tumour <- run_rf_cv_signature(X_tcga_tumour, y_tcga_tumour)
rf_ki_tumour   <- run_rf_cv_signature(X_ki_tumour, y_ki_tumour)

# Compute CI
ci_tcga <- ci.auc(rf_tcga_tumour$roc)
ci_ki   <- ci.auc(rf_ki_tumour$roc)

# Save summary
write.csv(
  data.frame(
    Model = c("TCGA tumour CV", "KI tumour CV"),
    AUC = c(rf_tcga_tumour$auc, rf_ki_tumour$auc),
    CI_lower = c(ci_tcga[1], ci_ki[1]),
    CI_upper = c(ci_tcga[3], ci_ki[3])
  ),
  file.path(output_dir, "tumour_multiomics_cv_summary.csv"),
  row.names = FALSE
)

# Plot ROC
ggsave(
  file.path(figure_dir, "tumour_multiomics_roc_tcga_cv.pdf"),
  plot_roc_curve(rf_tcga_tumour$roc, rf_tcga_tumour$auc, ci_tcga,
                 "SF-7D", "springgreen4"),
  width = 5, height = 4
)

ggsave(
  file.path(figure_dir, "tumour_multiomics_roc_ki_cv.pdf"),
  plot_roc_curve(rf_ki_tumour$roc, rf_ki_tumour$auc, ci_ki,
                 "SF-7C", "violetred4"),
  width = 5, height = 4
)

# ============================================================
# 7B. Within-cohort CV (healthy multi-omics, single CV + 95% CI)
# ============================================================

message("Running CV: healthy multi-omics...")

rf_tcga_healthy <- run_rf_cv_signature(X_tcga_healthy, y_tcga_healthy)
rf_ki_healthy   <- run_rf_cv_signature(X_ki_healthy, y_ki_healthy)

# Compute CI
ci_tcga_healthy <- ci.auc(rf_tcga_healthy$roc)
ci_ki_healthy   <- ci.auc(rf_ki_healthy$roc)

# Save summary
write.csv(
  data.frame(
    Model = c("TCGA healthy CV", "KI healthy CV"),
    AUC = c(rf_tcga_healthy$auc, rf_ki_healthy$auc),
    CI_lower = c(ci_tcga_healthy[1], ci_ki_healthy[1]),
    CI_upper = c(ci_tcga_healthy[3], ci_ki_healthy[3])
  ),
  file.path(output_dir, "healthy_multiomics_cv_summary.csv"),
  row.names = FALSE
)

# Plot ROC
ggsave(
  file.path(figure_dir, "healthy_multiomics_roc_tcga_cv.pdf"),
  plot_roc_curve(
    rf_tcga_healthy$roc,
    rf_tcga_healthy$auc,
    ci_tcga_healthy,
    "SF-7H",
    "springgreen4"
  ),
  width = 5,
  height = 4
)

ggsave(
  file.path(figure_dir, "healthy_multiomics_roc_ki_cv.pdf"),
  plot_roc_curve(
    rf_ki_healthy$roc,
    rf_ki_healthy$auc,
    ci_ki_healthy,
    "SF-7G",
    "violetred4"
  ),
  width = 5,
  height = 4
)

# ============================================================
# 7C. Within-cohort CV (all-sample multi-omics, single CV + 95% CI)
# ============================================================

message("Running CV: all-sample multi-omics...")

rf_tcga_all <- run_rf_cv_signature(X_tcga_all, y_tcga_all)
rf_ki_all   <- run_rf_cv_signature(X_ki_all, y_ki_all)

# Compute CI
ci_tcga_all <- ci.auc(rf_tcga_all$roc)
ci_ki_all   <- ci.auc(rf_ki_all$roc)

# Save summary
write.csv(
  data.frame(
    Model = c("TCGA all-sample CV", "KI all-sample CV"),
    AUC = c(rf_tcga_all$auc, rf_ki_all$auc),
    CI_lower = c(ci_tcga_all[1], ci_ki_all[1]),
    CI_upper = c(ci_tcga_all[3], ci_ki_all[3])
  ),
  file.path(output_dir, "all_multiomics_cv_summary.csv"),
  row.names = FALSE
)

# Plot ROC
ggsave(
  file.path(figure_dir, "all_multiomics_roc_tcga_cv.pdf"),
  plot_roc_curve(
    rf_tcga_all$roc,
    rf_tcga_all$auc,
    ci_tcga_all,
    "SF-7J",
    "springgreen4"
  ),
  width = 5,
  height = 4
)

ggsave(
  file.path(figure_dir, "all_multiomics_roc_ki_cv.pdf"),
  plot_roc_curve(
    rf_ki_all$roc,
    rf_ki_all$auc,
    ci_ki_all,
    "SF-7I",
    "violetred4"
  ),
  width = 5,
  height = 4
)

# ============================================================
# 8. Cross-cohort transfer (tumour multi-omics)
# ============================================================

message("Running transfer: tumour multi-omics...")
res_tcga_to_ki_tumour <- run_plsda_rf_transfer(X_tcga_tumour, y_tcga_tumour, X_ki_tumour, y_ki_tumour)
res_ki_to_tcga_tumour <- run_plsda_rf_transfer(X_ki_tumour, y_ki_tumour, X_tcga_tumour, y_tcga_tumour)

write.csv(
  data.frame(
    Model = c("TCGA -> KI tumour", "KI -> TCGA tumour"),
    AUC = c(res_tcga_to_ki_tumour$auc, res_ki_to_tcga_tumour$auc),
    CI_lower = c(res_tcga_to_ki_tumour$ci[1], res_ki_to_tcga_tumour$ci[1]),
    CI_upper = c(res_tcga_to_ki_tumour$ci[3], res_ki_to_tcga_tumour$ci[3]),
    Accuracy = c(res_tcga_to_ki_tumour$accuracy, res_ki_to_tcga_tumour$accuracy)
  ),
  file.path(output_dir, "tumour_multiomics_transfer_summary.csv"),
  row.names = FALSE
)

write.csv(res_tcga_to_ki_tumour$importance, file.path(output_dir, "tumour_multiomics_importance_tcga_to_ki.csv"), row.names = FALSE)
write.csv(res_ki_to_tcga_tumour$importance, file.path(output_dir, "tumour_multiomics_importance_ki_to_tcga.csv"), row.names = FALSE)

ggsave(
  file.path(figure_dir, "tumour_multiomics_roc_tcga_to_ki.pdf"),
  plot_roc_curve(res_tcga_to_ki_tumour$roc, res_tcga_to_ki_tumour$auc, res_tcga_to_ki_tumour$ci, "TCGA -> KI tumour", "springgreen4"),
  width = 5, height = 4
)
ggsave(
  file.path(figure_dir, "tumour_multiomics_roc_ki_to_tcga.pdf"),
  plot_roc_curve(res_ki_to_tcga_tumour$roc, res_ki_to_tcga_tumour$auc, res_ki_to_tcga_tumour$ci, "KI -> TCGA tumour", "violetred4"),
  width = 5, height = 4
)

ggsave(
  file.path(figure_dir, "tumour_multiomics_mdg_tcga_to_ki.pdf"),
  plot_mdg(res_tcga_to_ki_tumour$importance, "TCGA -> KI tumour", "springgreen4", top_n = 10),
  width = 5, height = 4
)
ggsave(
  file.path(figure_dir, "tumour_multiomics_mdg_ki_to_tcga.pdf"),
  plot_mdg(res_ki_to_tcga_tumour$importance, "KI -> TCGA tumour", "violetred4", top_n = 10),
  width = 5, height = 4
)

table(y_tcga_tumour)
table(y_ki_tumour)
levels(y_tcga_tumour)
levels(y_ki_tumour)

# ============================================================
# 9. DE-informed transfer analyses
# ============================================================

message("Running DE-informed transfer analyses...")

# Tumour
de_tumour_tcga_to_ki <- run_de_signature_transfer(
  expr_train = X_tcga_tumour, y_train = y_tcga_tumour,
  expr_test = X_ki_tumour, y_test = y_ki_tumour,
  train_name = "TCGA", test_name = "KI"
)

de_tumour_ki_to_tcga <- run_de_signature_transfer(
  expr_train = X_ki_tumour, y_train = y_ki_tumour,
  expr_test = X_tcga_tumour, y_test = y_tcga_tumour,
  train_name = "KI", test_name = "TCGA"
)

# Healthy
de_healthy_tcga_to_ki <- run_de_signature_transfer(
  expr_train = X_tcga_healthy, y_train = y_tcga_healthy,
  expr_test = X_ki_healthy, y_test = y_ki_healthy,
  train_name = "TCGA", test_name = "KI"
)

de_healthy_ki_to_tcga <- run_de_signature_transfer(
  expr_train = X_ki_healthy, y_train = y_ki_healthy,
  expr_test = X_tcga_healthy, y_test = y_tcga_healthy,
  train_name = "KI", test_name = "TCGA"
)

# All
de_all_tcga_to_ki <- run_de_signature_transfer(
  expr_train = X_tcga_all, y_train = y_tcga_all,
  expr_test = X_ki_all, y_test = y_ki_all,
  train_name = "TCGA", test_name = "KI"
)

de_all_ki_to_tcga <- run_de_signature_transfer(
  expr_train = X_ki_all, y_train = y_ki_all,
  expr_test = X_tcga_all, y_test = y_tcga_all,
  train_name = "KI", test_name = "TCGA"
)

# 95% CI for external ROC curves
ci_de_tumour_tcga_to_ki <- ci.auc(de_tumour_tcga_to_ki$roc_test)
ci_de_tumour_ki_to_tcga <- ci.auc(de_tumour_ki_to_tcga$roc_test)

ci_de_healthy_tcga_to_ki <- ci.auc(de_healthy_tcga_to_ki$roc_test)
ci_de_healthy_ki_to_tcga <- ci.auc(de_healthy_ki_to_tcga$roc_test)

ci_de_all_tcga_to_ki <- ci.auc(de_all_tcga_to_ki$roc_test)
ci_de_all_ki_to_tcga <- ci.auc(de_all_ki_to_tcga$roc_test)

de_summary <- data.frame(
  Analysis = c(
    "Tumour TCGA -> KI", "Tumour KI -> TCGA",
    "Healthy TCGA -> KI", "Healthy KI -> TCGA",
    "All TCGA -> KI", "All KI -> TCGA"
  ),
  Internal_AUC = c(
    de_tumour_tcga_to_ki$auc_train, de_tumour_ki_to_tcga$auc_train,
    de_healthy_tcga_to_ki$auc_train, de_healthy_ki_to_tcga$auc_train,
    de_all_tcga_to_ki$auc_train, de_all_ki_to_tcga$auc_train
  ),
  External_AUC = c(
    de_tumour_tcga_to_ki$auc_test, de_tumour_ki_to_tcga$auc_test,
    de_healthy_tcga_to_ki$auc_test, de_healthy_ki_to_tcga$auc_test,
    de_all_tcga_to_ki$auc_test, de_all_ki_to_tcga$auc_test
  ),
  CI_lower = c(
    ci_de_tumour_tcga_to_ki[1], ci_de_tumour_ki_to_tcga[1],
    ci_de_healthy_tcga_to_ki[1], ci_de_healthy_ki_to_tcga[1],
    ci_de_all_tcga_to_ki[1], ci_de_all_ki_to_tcga[1]
  ),
  CI_upper = c(
    ci_de_tumour_tcga_to_ki[3], ci_de_tumour_ki_to_tcga[3],
    ci_de_healthy_tcga_to_ki[3], ci_de_healthy_ki_to_tcga[3],
    ci_de_all_tcga_to_ki[3], ci_de_all_ki_to_tcga[3]
  ),
  Signature_Size = c(
    length(de_tumour_tcga_to_ki$signature_features), length(de_tumour_ki_to_tcga$signature_features),
    length(de_healthy_tcga_to_ki$signature_features), length(de_healthy_ki_to_tcga$signature_features),
    length(de_all_tcga_to_ki$signature_features), length(de_all_ki_to_tcga$signature_features)
  ),
  stringsAsFactors = FALSE
)
write.csv(de_summary, file.path(output_dir, "de_informed_transfer_summary.csv"), row.names = FALSE)

write.csv(de_tumour_tcga_to_ki$merged, file.path(output_dir, "de_tumour_tcga_to_ki_merged.csv"), row.names = FALSE)
write.csv(de_tumour_ki_to_tcga$merged, file.path(output_dir, "de_tumour_ki_to_tcga_merged.csv"), row.names = FALSE)
write.csv(de_healthy_tcga_to_ki$merged, file.path(output_dir, "de_healthy_tcga_to_ki_merged.csv"), row.names = FALSE)
write.csv(de_healthy_ki_to_tcga$merged, file.path(output_dir, "de_healthy_ki_to_tcga_merged.csv"), row.names = FALSE)
write.csv(de_all_tcga_to_ki$merged, file.path(output_dir, "de_all_tcga_to_ki_merged.csv"), row.names = FALSE)
write.csv(de_all_ki_to_tcga$merged, file.path(output_dir, "de_all_ki_to_tcga_merged.csv"), row.names = FALSE)

# Example ROC outputs for DE-informed models
ggsave(
  file.path(figure_dir, "de_tumour_tcga_to_ki_external_roc.pdf"),
  plot_roc_curve(
    de_tumour_tcga_to_ki$roc_test,
    de_tumour_tcga_to_ki$auc_test,
    ci_de_tumour_tcga_to_ki,
    "SF-7F",
    "springgreen4"
  ),
  width = 5, height = 4
)

ggsave(
  file.path(figure_dir, "de_healthy_tcga_to_ki_external_roc.pdf"),
  plot_roc_curve(
    de_healthy_tcga_to_ki$roc_test,
    de_healthy_tcga_to_ki$auc_test,
    ci_de_healthy_tcga_to_ki,
    "Fig. 5C",
    "springgreen4"
  ),
  width = 5, height = 4
)

ggsave(
  file.path(figure_dir, "de_all_tcga_to_ki_external_roc.pdf"),
  plot_roc_curve(
    de_all_tcga_to_ki$roc_test,
    de_all_tcga_to_ki$auc_test,
    ci_de_all_tcga_to_ki,
    "SF-7K",
    "springgreen4"
  ),
  width = 5, height = 4
)

ggsave(
  file.path(figure_dir, "de_tumour_ki_to_tcga_external_roc.pdf"),
  plot_roc_curve(
    de_tumour_ki_to_tcga$roc_test,
    de_tumour_ki_to_tcga$auc_test,
    ci_de_tumour_ki_to_tcga,
    "SF-7E",
    "violetred4"
  ),
  width = 5, height = 4
)

ggsave(
  file.path(figure_dir, "de_healthy_ki_to_tcga_external_roc.pdf"),
  plot_roc_curve(
    de_healthy_ki_to_tcga$roc_test,
    de_healthy_ki_to_tcga$auc_test,
    ci_de_healthy_ki_to_tcga,
    "Fig. 5B",
    "violetred4"
  ),
  width = 5, height = 4
)

ggsave(
  file.path(figure_dir, "de_all_ki_to_tcga_external_roc.pdf"),
  plot_roc_curve(
    de_all_ki_to_tcga$roc_test,
    de_all_ki_to_tcga$auc_test,
    ci_de_all_ki_to_tcga,
    "Fig. 5D",
    "violetred4"
  ),
  width = 5, height = 4
)

# ============================================================
# 10. RF CV on DE-informed signatures
# ============================================================

message("Running RF CV on DE-informed signatures...")

# Tumour
rf_sig_tcga_tumour <- run_rf_cv_signature(
  X_tcga_tumour[, intersect(de_tumour_tcga_to_ki$signature_features, colnames(X_tcga_tumour)), drop = FALSE],
  y_tcga_tumour
)
rf_sig_ki_tumour <- run_rf_cv_signature(
  X_ki_tumour[, intersect(de_tumour_ki_to_tcga$signature_features, colnames(X_ki_tumour)), drop = FALSE],
  y_ki_tumour
)

# Healthy
rf_sig_tcga_healthy <- run_rf_cv_signature(
  X_tcga_healthy[, intersect(de_healthy_tcga_to_ki$signature_features, colnames(X_tcga_healthy)), drop = FALSE],
  y_tcga_healthy
)
rf_sig_ki_healthy <- run_rf_cv_signature(
  X_ki_healthy[, intersect(de_healthy_ki_to_tcga$signature_features, colnames(X_ki_healthy)), drop = FALSE],
  y_ki_healthy
)

# All
rf_sig_tcga_all <- run_rf_cv_signature(
  X_tcga_all[, intersect(de_all_tcga_to_ki$signature_features, colnames(X_tcga_all)), drop = FALSE],
  y_tcga_all
)
rf_sig_ki_all <- run_rf_cv_signature(
  X_ki_all[, intersect(de_all_ki_to_tcga$signature_features, colnames(X_ki_all)), drop = FALSE],
  y_ki_all
)

rf_sig_summary <- data.frame(
  Analysis = c(
    "TCGA tumour DE-signature RF CV",
    "KI tumour DE-signature RF CV",
    "TCGA healthy DE-signature RF CV",
    "KI healthy DE-signature RF CV",
    "TCGA all DE-signature RF CV",
    "KI all DE-signature RF CV"
  ),
  AUC = c(
    rf_sig_tcga_tumour$auc,
    rf_sig_ki_tumour$auc,
    rf_sig_tcga_healthy$auc,
    rf_sig_ki_healthy$auc,
    rf_sig_tcga_all$auc,
    rf_sig_ki_all$auc
  ),
  stringsAsFactors = FALSE
)
write.csv(rf_sig_summary, file.path(output_dir, "de_signature_rf_cv_summary.csv"), row.names = FALSE)

# ============================================================
# 11. Adjacent recurrence score distribution
# ============================================================

up_adj <- intersect(de_healthy_tcga_to_ki$up_features, colnames(expr_healthy))
down_adj <- intersect(de_healthy_tcga_to_ki$down_features, colnames(expr_healthy))

plot_df_adj <- data.frame(
  Sample = rownames(expr_healthy),
  Cohort = as.character(src_healthy),
  Outcome = factor(as.character(y_healthy), levels = c("Remission", "Recidivism")),
  stringsAsFactors = FALSE
)

up_score <- if (length(up_adj) > 0) rowMeans(expr_healthy[, up_adj, drop = FALSE], na.rm = TRUE) else rep(0, nrow(expr_healthy))
down_score <- if (length(down_adj) > 0) rowMeans(expr_healthy[, down_adj, drop = FALSE], na.rm = TRUE) else rep(0, nrow(expr_healthy))
plot_df_adj$RecurrenceScore <- up_score - down_score

pvals_adj <- plot_df_adj %>%
  group_by(Cohort) %>%
  summarise(
    p_value = tryCatch(wilcox.test(RecurrenceScore ~ Outcome)$p.value, error = function(e) NA_real_),
    ymax = max(RecurrenceScore, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(is.na(p_value), "p = NA", paste0("p = ", signif(p_value, 2))),
    x = 1.5,
    y = ymax + 0.05 * diff(range(plot_df_adj$RecurrenceScore, na.rm = TRUE))
  )

write.csv(plot_df_adj, file.path(output_dir, "adjacent_recurrence_score_distribution.csv"), row.names = FALSE)

p_adj_score <- ggplot(plot_df_adj, aes(x = Outcome, y = RecurrenceScore, fill = Outcome)) +
  geom_violin(trim = FALSE, alpha = 0.85, color = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.75, color = "black") +
  facet_wrap(~ Cohort, scales = "free_y") +
  scale_fill_manual(values = c("Remission" = "#377EB8", "Recidivism" = "#E41A1C")) +
  geom_text(data = pvals_adj, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 4) +
  labs(
    title = "Fig. 5F",
    x = NULL,
    y = "Recurrence score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black")
  )

ggsave(file.path(figure_dir, "adjacent_recurrence_score_distribution.pdf"), p_adj_score, width = 7, height = 5)

# ============================================================
# 12. Clinical associations of adjacent recurrence score
# ============================================================

meta <- metadata
meta$Sample <- normalize_ids(meta$mRNA)
meta$Sample <- remove_leading_X(meta$Sample)
meta$Sample <- fix_ki_mrna_ids(meta$Sample)
meta$Sample <- fix_tcga_ids(meta$Sample)

clinical_df <- plot_df_adj %>%
  left_join(meta, by = "Sample")

clinical_df$Cirrosis <- factor(clinical_df$Cirrosis, levels = c(0, 1), labels = c("No", "Yes"))
clinical_df$Milan    <- factor(clinical_df$Milan, levels = c(0, 1), labels = c("Out", "In"))
clinical_df$AFP20    <- factor(clinical_df$AFP20, levels = c(0, 1), labels = c("<20", ">=20"))
clinical_df$Sex      <- factor(clinical_df$Sex, levels = c(0, 1), labels = c("Female", "Male"))
clinical_df$TMN      <- droplevels(factor(clinical_df$TMN))
clinical_df$PS       <- factor(clinical_df$PS)
clinical_df$DAS      <- factor(clinical_df$DAS)
clinical_df$S123     <- factor(clinical_df$S123)

theme_clin <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(color = "black", size = 8.5),
    axis.text.y = element_text(color = "black", size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  )

# AFP
df_afp <- clinical_df %>%
  filter(
    !is.na(AFP),
    is.finite(AFP),
    AFP > 0,
    !is.na(RecurrenceScore),
    is.finite(RecurrenceScore)
  )

if (nrow(df_afp) < 3) {
  warning("Not enough finite observations for AFP correlation.")
  
  p_afp <- ggplot() +
    annotate("text", x = 1, y = 1, label = "Not enough AFP data") +
    xlim(0, 2) + ylim(0, 2) +
    labs(title = "AFP", x = NULL, y = NULL) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  
} else {
  afp_test <- cor.test(df_afp$AFP, df_afp$RecurrenceScore, method = "spearman")
  afp_label <- paste0(
    "rho = ", round(unname(afp_test$estimate), 2), "\n",
    format_p(afp_test$p.value)
  )
  
  p_afp <- ggplot(df_afp, aes(x = AFP, y = RecurrenceScore, color = Outcome)) +
    geom_point(size = 2.2, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    annotate(
      "text",
      x = min(df_afp$AFP, na.rm = TRUE) * 1.3,
      y = get_y_top(df_afp$RecurrenceScore, expand = 0.01),
      label = afp_label,
      hjust = 0,
      size = 3.3
    ) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
    scale_color_manual(values = c("Remission" = "#377EB8", "Recidivism" = "#E41A1C")) +
    labs(
      title = "AFP",
      x = "AFP (log scale)",
      y = "Adjacent recurrence score"
    ) +
    theme_clin
}

# Cirrhosis
df_cirr <- clinical_df %>% filter(!is.na(Cirrosis))
cirr_labels <- add_n_labels(df_cirr, "Cirrosis")
cirr_test <- wilcox.test(RecurrenceScore ~ Cirrosis, data = df_cirr)

p_cirr <- ggplot(df_cirr, aes(x = Cirrosis, y = RecurrenceScore, fill = Cirrosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.62, linewidth = 0.7) +
  geom_jitter(width = 0.10, alpha = 0.55, size = 1.7, color = "black") +
  annotate(
    "text",
    x = 1.5,
    y = get_y_top(df_cirr$RecurrenceScore, 0.015),
    label = format_p(cirr_test$p.value),
    size = 3.3
  ) +
  scale_x_discrete(labels = cirr_labels) +
  scale_fill_manual(values = c("No" = "#377EB8", "Yes" = "#E41A1C")) +
  labs(title = "Cirrhosis", x = NULL, y = "Adjacent recurrence score") +
  theme_clin +
  theme(legend.position = "none")

# Milan
df_milan <- clinical_df %>% filter(!is.na(Milan))
milan_labels <- add_n_labels(df_milan, "Milan")
milan_test <- wilcox.test(RecurrenceScore ~ Milan, data = df_milan)

p_milan <- ggplot(df_milan, aes(x = Milan, y = RecurrenceScore, fill = Milan)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.62, linewidth = 0.7) +
  geom_jitter(width = 0.10, alpha = 0.55, size = 1.7, color = "black") +
  annotate(
    "text",
    x = 1.5,
    y = get_y_top(df_milan$RecurrenceScore, 0.015),
    label = format_p(milan_test$p.value),
    size = 3.3
  ) +
  scale_x_discrete(labels = milan_labels) +
  scale_fill_manual(values = c("Out" = "#377EB8", "In" = "#E41A1C")) +
  labs(title = "Milan criteria", x = NULL, y = "Adjacent recurrence score") +
  theme_clin +
  theme(legend.position = "none")

# AFP20
df_afp20 <- clinical_df %>% filter(!is.na(AFP20))
afp20_labels <- add_n_labels(df_afp20, "AFP20")
afp20_test <- wilcox.test(RecurrenceScore ~ AFP20, data = df_afp20)

p_afp20 <- ggplot(df_afp20, aes(x = AFP20, y = RecurrenceScore, fill = AFP20)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.62, linewidth = 0.7) +
  geom_jitter(width = 0.10, alpha = 0.55, size = 1.7, color = "black") +
  annotate(
    "text",
    x = 1.5,
    y = get_y_top(df_afp20$RecurrenceScore, 0.015),
    label = format_p(afp20_test$p.value),
    size = 3.3
  ) +
  scale_x_discrete(labels = afp20_labels) +
  scale_fill_manual(values = c("<20" = "#377EB8", ">=20" = "#E41A1C")) +
  labs(title = "AFP20", x = NULL, y = "Adjacent recurrence score") +
  theme_clin +
  theme(legend.position = "none")

# TMN
df_tmn <- clinical_df %>%
  filter(!is.na(TMN)) %>%
  filter(!TMN %in% c("1000000", "100000")) %>%
  droplevels()

tmn_test <- kruskal.test(RecurrenceScore ~ TMN, data = df_tmn)

p_tmn <- ggplot(df_tmn, aes(x = TMN, y = RecurrenceScore, fill = TMN)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.62, linewidth = 0.7) +
  geom_jitter(width = 0.10, alpha = 0.55, size = 1.7, color = "black") +
  annotate(
    "text",
    x = ceiling(length(unique(df_tmn$TMN)) / 2),
    y = get_y_top(df_tmn$RecurrenceScore, 0.015),
    label = format_p(tmn_test$p.value),
    size = 3.3
  ) +
  labs(title = "TMN", x = NULL, y = "Adjacent recurrence score") +
  theme_clin +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 9)
  )

# PS
df_ps <- clinical_df %>% filter(!is.na(PS))
ps_labels <- add_n_labels(df_ps, "PS")
ps_test <- kruskal.test(RecurrenceScore ~ PS, data = df_ps)

p_ps <- ggplot(df_ps, aes(x = PS, y = RecurrenceScore, fill = PS)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.62, linewidth = 0.7) +
  geom_jitter(width = 0.10, alpha = 0.55, size = 1.7, color = "black") +
  annotate(
    "text",
    x = ceiling(length(unique(df_ps$PS)) / 2),
    y = get_y_top(df_ps$RecurrenceScore, 0.015),
    label = format_p(ps_test$p.value),
    size = 3.3
  ) +
  scale_x_discrete(labels = ps_labels) +
  labs(title = "PS", x = NULL, y = "Adjacent recurrence score") +
  theme_clin +
  theme(legend.position = "none")

p_clinical <- (p_afp | p_cirr | p_milan) / (p_afp20 | p_tmn | p_ps)

ggsave(
  file.path(figure_dir, "adjacent_score_clinical_associations.pdf"),
  p_clinical,
  width = 12,
  height = 8
)

# ============================================================
# 13. Viral status integration and multivariable models
# ============================================================

viral_info <- read_excel(viral_file)

viral_info$Sample <- normalize_ids(viral_info$mRNA)
viral_info$Sample <- remove_leading_X(viral_info$Sample)
viral_info$Sample <- fix_ki_mrna_ids(viral_info$Sample)
viral_info$Sample <- fix_tcga_ids(viral_info$Sample)

viral_info2 <- viral_info %>%
  dplyr::select(Sample, `Healthy/Viral/Non-viral`) %>%
  dplyr::rename(ViralStatus = `Healthy/Viral/Non-viral`)

clinical_df2 <- clinical_df %>%
  left_join(viral_info2, by = "Sample")

model_df_full <- clinical_df2 %>%
  dplyr::select(Outcome, RecurrenceScore, Cirrosis, Milan, PS) %>%
  na.omit()

model_df_full$Outcome <- factor(model_df_full$Outcome, levels = c("Remission", "Recidivism"))
model_df_full$Score_z <- as.numeric(scale(model_df_full$RecurrenceScore))

clinical_df_viral <- clinical_df2 %>%
  filter(ViralStatus %in% c("Non-viral", "Viral")) %>%
  mutate(ViralStatus = factor(ViralStatus, levels = c("Non-viral", "Viral")))

model_df_viral <- clinical_df_viral %>%
  dplyr::select(Outcome, RecurrenceScore, Cirrosis, Milan, PS, ViralStatus) %>%
  na.omit()

model_df_viral$Outcome <- factor(model_df_viral$Outcome, levels = c("Remission", "Recidivism"))
model_df_viral$Score_z <- as.numeric(scale(model_df_viral$RecurrenceScore))

model_full <- glm(
  Outcome ~ Score_z + Cirrosis + Milan + PS,
  data = model_df_full,
  family = binomial
)

model_viral <- glm(
  Outcome ~ Score_z + Cirrosis + Milan + PS + ViralStatus,
  data = model_df_viral,
  family = binomial
)

forest_df_full <- make_forest_df(model_full, "full")
forest_df_viral <- make_forest_df(model_viral, "viral")

write.csv(
  forest_df_full,
  file.path(output_dir, "adjacent_score_full_model_forest_df.csv"),
  row.names = FALSE
)

write.csv(
  forest_df_viral,
  file.path(output_dir, "adjacent_score_viral_model_forest_df.csv"),
  row.names = FALSE
)

p_forest_full <- plot_forest(
  forest_df_full,
  "Figure 5F"
)

p_forest_viral <- plot_forest(
  forest_df_viral,
  "Adjusted odds ratios for recurrence\nRestricted to annotated viral status"
)

ggsave(
  file.path(figure_dir, "adjacent_score_forest_full.pdf"),
  p_forest_full,
  width = 6,
  height = 4.5
)

ggsave(
  file.path(figure_dir, "adjacent_score_forest_viral.pdf"),
  p_forest_viral,
  width = 6,
  height = 5
)

# Interaction model
model_interaction <- glm(
  Outcome ~ Score_z * ViralStatus + Cirrosis + Milan + PS,
  data = model_df_viral,
  family = binomial
)

interaction_or <- as.data.frame(
  exp(cbind(OR = coef(model_interaction), confint(model_interaction)))
)
interaction_or$Variable <- rownames(interaction_or)

write.csv(
  interaction_or,
  file.path(output_dir, "adjacent_score_interaction_model_or.csv"),
  row.names = FALSE
)

# Viral status boxplot
df_plot_viral <- clinical_df2 %>%
  filter(ViralStatus %in% c("Non-viral", "Viral")) %>%
  mutate(ViralStatus = factor(ViralStatus, levels = c("Non-viral", "Viral")))

viral_labels <- add_n_labels(df_plot_viral, "ViralStatus")
viral_test <- wilcox.test(RecurrenceScore ~ ViralStatus, data = df_plot_viral)

p_viral <- ggplot(df_plot_viral, aes(x = ViralStatus, y = RecurrenceScore, fill = ViralStatus)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.62, linewidth = 0.7) +
  geom_jitter(width = 0.10, alpha = 0.55, size = 1.7, color = "black") +
  annotate(
    "text",
    x = 1.5,
    y = get_y_top(df_plot_viral$RecurrenceScore, expand = 0.015),
    label = format_p(viral_test$p.value),
    size = 3.3
  ) +
  scale_x_discrete(labels = viral_labels) +
  scale_fill_manual(values = c("Non-viral" = "#377EB8", "Viral" = "#E41A1C")) +
  labs(
    title = "SF-6D",
    x = NULL,
    y = "Adjacent recurrence score"
  ) +
  theme_clin +
  theme(legend.position = "none")

ggsave(
  file.path(figure_dir, "adjacent_score_by_viral_status.pdf"),
  p_viral,
  width = 5,
  height = 4
)

write.csv(
  df_plot_viral %>% count(ViralStatus),
  file.path(output_dir, "viral_status_counts.csv"),
  row.names = FALSE
)