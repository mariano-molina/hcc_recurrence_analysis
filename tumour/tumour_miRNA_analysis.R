# ============================================================
# Liver tumour recurrence analysis - miRNA
#
# Author: MAM
# Date: April 2026
#
# Description:
# This script reproduces the main tumour-based miRNA analyses presented
# in the manuscript, including preprocessing, cohort structure
# visualization (t-SNE), within-cohort cross-validation,
# cross-cohort prediction, and feature ranking.
#
# Input:
# - merge_mirna.csv
# - mRNA_miRNA_patient_sample_info.csv
#
# Output:
# - Figures (PDF) in results/figures/
# - Processed objects (.rds) in results/objects/
# - Summary tables (.csv) in results/
#
# Notes:
# - Input files are not included in this repository due to data sharing
#   restrictions.
# - For miRNA data, TCGA sample IDs are kept as-is, whereas non-TCGA
#   metadata sample IDs are converted from dot-separated to
#   underscore-separated format to match the expression matrix.
# ============================================================

# ============================================================
# 0. Setup
# ============================================================

suppressPackageStartupMessages({
  library(Rtsne)
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
  library(UpSetR)
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
expression_file <- file.path(input_dir, "merge_mirna.csv")
metadata_file   <- file.path(input_dir, "mRNA_miRNA_patient_sample_info.csv")

# Optional predefined QC exclusion
# Set to NA_character_ if not needed
bad_sample <- NA_character_

# ============================================================
# 1. Load and preprocess data
# ============================================================

message("Loading miRNA matrix and metadata...")

# Keep default column-name conversion for miRNA expression data,
# as in the original script
expr_raw <- read.csv(expression_file, stringsAsFactors = FALSE)

# Keep metadata names unchanged
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)

# Remove miRNAs with missing/empty identifiers
expr_clean <- expr_raw[!(is.na(expr_raw[[1]]) | expr_raw[[1]] == ""), ]

# Remove miRNAs with zero expression across all samples
expr_clean <- expr_clean[rowSums(expr_clean[, -1] != 0) > 0, ]

# Extract sample IDs from expression matrix
sample_ids <- colnames(expr_clean)[-1]

# Standardize sample IDs before matching
sample_ids_clean <- trimws(sample_ids)

# Use the same metadata ID logic as in the original miRNA script:
# the first metadata column contains the matching sample identifiers
metadata_ids_clean <- trimws(as.character(metadata[[1]]))

# Original miRNA-specific harmonization:
# keep TCGA IDs as-is, convert non-TCGA metadata IDs from dots to underscores
metadata_ids_match <- ifelse(
  grepl("^TCGA", metadata_ids_clean),
  metadata_ids_clean,
  gsub("\\.", "_", metadata_ids_clean)
)

sample_ids_match <- sample_ids_clean

# Build lookup tables from metadata
recurrence_lookup <- setNames(metadata$Recidiv, metadata_ids_match)
tissue_lookup     <- setNames(metadata$Tumour, metadata_ids_match)

# Match sample-level labels to the expression matrix
recurrence_labels <- sapply(sample_ids_match, function(id) {
  value <- recurrence_lookup[id]
  if (is.na(value)) return(NA_character_)
  if (value == 1) return("Recidivism")
  if (value == 0) return("Remission")
  return(NA_character_)
})

tissue_labels <- sapply(sample_ids_match, function(id) {
  value <- tissue_lookup[id]
  if (is.na(value)) return(NA_character_)
  if (value == 1) return("Tumour")
  if (value == 0) return("Healthy")
  return(NA_character_)
})

message("Number of expression samples: ", length(sample_ids_clean))
message("Matched recurrence labels: ", sum(!is.na(recurrence_labels)))
message("Matched tissue labels: ", sum(!is.na(tissue_labels)))

# Build samples x miRNAs numeric matrix
expr_matrix <- apply(expr_clean[, -1], 2, as.numeric)
expr_matrix <- t(expr_matrix)
rownames(expr_matrix) <- sample_ids_clean
colnames(expr_matrix) <- make.unique(expr_clean[[1]])

# Remove samples with missing labels
valid_samples <- !is.na(recurrence_labels) & !is.na(tissue_labels)

if (sum(valid_samples) == 0) {
  stop("No samples remained after matching miRNA data to metadata. Check sample ID formatting and metadata sample ID column.")
}

expr_matrix <- expr_matrix[valid_samples, , drop = FALSE]
recurrence_labels <- factor(recurrence_labels[valid_samples], levels = c("Remission", "Recidivism"))
tissue_labels     <- factor(tissue_labels[valid_samples], levels = c("Healthy", "Tumour"))

# Remove zero-variance miRNAs across retained samples
expr_matrix <- expr_matrix[, apply(expr_matrix, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]

# Log-transform expression values
mode(expr_matrix) <- "numeric"
log_expr_matrix <- log2(expr_matrix + 1)

# Assign names to label vectors
names(recurrence_labels) <- rownames(log_expr_matrix)
names(tissue_labels)     <- rownames(log_expr_matrix)

# Remove predefined outlier sample from final analyses, if specified
if (!is.na(bad_sample) && bad_sample %in% rownames(log_expr_matrix)) {
  log_expr_matrix <- log_expr_matrix[rownames(log_expr_matrix) != bad_sample, , drop = FALSE]
  recurrence_labels <- recurrence_labels[names(recurrence_labels) != bad_sample]
  tissue_labels     <- tissue_labels[names(tissue_labels) != bad_sample]
}

# Restrict tumour-only dataset
tumour_samples <- names(tissue_labels)[tissue_labels == "Tumour"]
expr_tumour    <- log_expr_matrix[tumour_samples, , drop = FALSE]
labels_tumour  <- recurrence_labels[tumour_samples]

# Save processed objects
saveRDS(log_expr_matrix, file.path(object_dir, "mirna_log_expr_matrix.rds"))
saveRDS(recurrence_labels, file.path(object_dir, "mirna_recurrence_labels.rds"))
saveRDS(tissue_labels, file.path(object_dir, "mirna_tissue_labels.rds"))
saveRDS(expr_tumour, file.path(object_dir, "mirna_expr_tumour.rds"))
saveRDS(labels_tumour, file.path(object_dir, "mirna_labels_tumour.rds"))

# ============================================================
# 2. Exploratory t-SNE for cohort structure and QC
# ============================================================

message("Running t-SNE...")
set.seed(42)

tsne_input <- as.matrix(log_expr_matrix)
n_samples <- nrow(tsne_input)
perplexity_value <- min(30, floor((n_samples - 1) / 3) - 1)

if (perplexity_value < 2) {
  stop("Too few samples for t-SNE after filtering.")
}

message("Using t-SNE perplexity: ", perplexity_value)

tsne_result <- Rtsne(
  tsne_input,
  perplexity = perplexity_value,
  verbose = TRUE,
  max_iter = 1000
)

tsne_df <- data.frame(
  TSNE1 = tsne_result$Y[, 1],
  TSNE2 = tsne_result$Y[, 2],
  Outcome = recurrence_labels,
  Status = tissue_labels,
  Source = ifelse(grepl("^TCGA", rownames(tsne_input)), "TCGA", "KI"),
  SampleID = rownames(tsne_input),
  stringsAsFactors = FALSE
)

plot_tsne <- function(data, fill_var, title_text, palette) {
  ggplot(data, aes(x = TSNE1, y = TSNE2, fill = .data[[fill_var]])) +
    geom_point(size = 2.5, shape = 21, color = "black", stroke = 0.5) +
    scale_fill_manual(values = palette) +
    labs(title = title_text, x = "t-SNE1", y = "t-SNE2") +
    theme_minimal(base_size = 13) +
    theme(
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

p_tsne_outcome <- plot_tsne(
  tsne_df, "Outcome", "Fig. 2A-4",
  c("Recidivism" = "#E41A1C", "Remission" = "#377EB8")
)

p_tsne_tissue <- plot_tsne(
  tsne_df, "Status", "Fig. 2A-3",
  c("Tumour" = "orange", "Healthy" = "#4DAF4A")
)

p_tsne_source <- plot_tsne(
  tsne_df, "Source", "SF-3A",
  c("TCGA" = "#93C47D", "KI" = "#912F6C")
)

ggsave(file.path(figure_dir, "mirna_tsne_outcome.pdf"), p_tsne_outcome, width = 5, height = 4)
ggsave(file.path(figure_dir, "mirna_tsne_tissue.pdf"), p_tsne_tissue, width = 5, height = 4)
ggsave(file.path(figure_dir, "mirna_tsne_source.pdf"), p_tsne_source, width = 5, height = 4)

# ============================================================
# 3. Utility functions for modelling
# ============================================================

run_plsda_rf_cv <- function(X, y, k = 5, ncomp = 2, n_top = 100, ntree = 500, seed = 999) {
  set.seed(seed)
  
  X <- as.matrix(X)
  y <- factor(y, levels = c("Remission", "Recidivism"))
  
  stopifnot(nrow(X) == length(y))
  stopifnot(!is.null(rownames(X)))
  
  fold_assignments <- rep(NA_integer_, length(y))
  names(fold_assignments) <- rownames(X)
  
  for (cls in levels(y)) {
    cls_samples <- sample(rownames(X)[y == cls])
    fold_ids <- rep(seq_len(k), length.out = length(cls_samples))
    fold_assignments[cls_samples] <- fold_ids
  }
  
  all_probs <- c()
  all_true <- c()
  rf_importance_list <- list()
  
  for (i in seq_len(k)) {
    test_samples <- names(fold_assignments)[fold_assignments == i]
    train_samples <- setdiff(rownames(X), test_samples)
    
    X_train_full <- X[train_samples, , drop = FALSE]
    X_test_full  <- X[test_samples, , drop = FALSE]
    y_train <- factor(y[train_samples], levels = c("Remission", "Recidivism"))
    y_test  <- factor(y[test_samples], levels = c("Remission", "Recidivism"))
    
    if (length(unique(y_train)) < 2) next
    
    non_zero_var <- apply(X_train_full, 2, sd, na.rm = TRUE) > 0
    X_train_plsda <- X_train_full[, non_zero_var, drop = FALSE]
    X_test_plsda  <- X_test_full[, non_zero_var, drop = FALSE]
    
    if (ncol(X_train_plsda) < 2) next
    
    colnames(X_train_plsda) <- make.unique(colnames(X_train_plsda))
    colnames(X_test_plsda)  <- colnames(X_train_plsda)
    
    plsda_fit <- plsda(X_train_plsda, y_train, ncomp = ncomp)
    
    loadings_comp1 <- plsda_fit$loadings$X[, 1]
    loadings_comp2 <- plsda_fit$loadings$X[, 2]
    loading_score <- abs(loadings_comp1) + abs(loadings_comp2)
    
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
      importance = TRUE
    )
    
    rf_importance_list[[length(rf_importance_list) + 1]] <- data.frame(
      Feature = rownames(rf_fit$importance),
      MeanDecreaseGini = rf_fit$importance[, "MeanDecreaseGini"],
      Fold = i,
      stringsAsFactors = FALSE
    )
    
    y_prob <- predict(rf_fit, newdata = X_test, type = "prob")[, "Recidivism"]
    all_probs <- c(all_probs, y_prob)
    all_true  <- c(all_true, as.character(y_test))
  }
  
  if (length(all_probs) == 0 || length(all_true) == 0) {
    stop("No predictions were generated during cross-validation.")
  }
  
  if (length(rf_importance_list) == 0) {
    stop("No variable importance values were collected during cross-validation.")
  }
  
  all_true <- factor(all_true, levels = c("Remission", "Recidivism"))
  
  roc_obj <- roc(
    response = all_true,
    predictor = all_probs,
    levels = c("Remission", "Recidivism"),
    direction = "<"
  )
  
  auc_val <- as.numeric(auc(roc_obj))
  ci_vals <- ci.auc(roc_obj)
  
  importance_df <- do.call(rbind, rf_importance_list)
  importance_summary <- aggregate(
    MeanDecreaseGini ~ Feature,
    data = importance_df,
    FUN = mean
  )
  importance_summary <- importance_summary[order(-importance_summary$MeanDecreaseGini), ]
  
  list(
    roc = roc_obj,
    auc = auc_val,
    ci = ci_vals,
    importance = importance_summary,
    ranked_features = importance_summary$Feature
  )
}

run_plsda_rf_transfer <- function(X_train_full, y_train, X_test_full, y_test,
                                  ncomp = 2, n_plsda = 100, n_rf = 100,
                                  ntree_rank = 1000, ntree_final = 500, seed = 999) {
  set.seed(seed)
  
  X_train_full <- as.matrix(X_train_full)
  X_test_full  <- as.matrix(X_test_full)
  y_train <- factor(y_train, levels = c("Remission", "Recidivism"))
  y_test  <- factor(y_test, levels = c("Remission", "Recidivism"))
  
  stopifnot(nrow(X_train_full) == length(y_train))
  stopifnot(nrow(X_test_full) == length(y_test))
  
  non_zero_var <- apply(X_train_full, 2, sd, na.rm = TRUE) > 0
  X_train_plsda <- X_train_full[, non_zero_var, drop = FALSE]
  X_test_plsda  <- X_test_full[, non_zero_var, drop = FALSE]
  
  if (ncol(X_train_plsda) < 2) {
    stop("Too few non-zero variance features available for transfer analysis.")
  }
  
  colnames(X_train_plsda) <- make.unique(colnames(X_train_plsda))
  colnames(X_test_plsda)  <- colnames(X_train_plsda)
  
  plsda_fit <- plsda(X_train_plsda, y_train, ncomp = ncomp)
  
  loadings_comp1 <- plsda_fit$loadings$X[, 1]
  loadings_comp2 <- plsda_fit$loadings$X[, 2]
  
  top_features_c1 <- names(sort(abs(loadings_comp1), decreasing = TRUE)[1:min(n_plsda, length(loadings_comp1))])
  top_features_c2 <- names(sort(abs(loadings_comp2), decreasing = TRUE)[1:min(n_plsda, length(loadings_comp2))])
  top_plsda_features <- union(top_features_c1, top_features_c2)
  
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
  
  if (length(selected_features) < 2) {
    stop("Too few selected features remained for transfer analysis.")
  }
  
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
  
  auc_val <- as.numeric(auc(roc_obj))
  ci_vals <- ci.auc(roc_obj)
  conf_mat <- table(True = y_test, Predicted = y_pred)
  acc <- mean(y_pred == y_test)
  
  list(
    roc = roc_obj,
    auc = auc_val,
    ci = ci_vals,
    accuracy = acc,
    confusion_matrix = conf_mat,
    importance = importance_df,
    ranked_features = ranked_features,
    selected_features = selected_features
  )
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
      fill = scales::alpha("white", 0.85),
      label.r = grid::unit(0.1, "lines"),
      hjust = 1
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

plot_mdg <- function(importance_df, title_text, fill_color, top_n = 10) {
  plot_df <- importance_df[order(-importance_df$MeanDecreaseGini), , drop = FALSE]
  plot_df <- plot_df[1:min(top_n, nrow(plot_df)), , drop = FALSE]
  plot_df$Feature <- factor(plot_df$Feature, levels = rev(plot_df$Feature))
  
  ggplot(plot_df, aes(x = Feature, y = MeanDecreaseGini)) +
    geom_col(fill = fill_color) +
    coord_flip() +
    theme_minimal(base_size = 14) +
    labs(
      title = title_text,
      x = "miRNA",
      y = "Mean Decrease Gini"
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

plot_plsda_scores <- function(plsda_fit, y, title_text) {
  scores_df <- data.frame(
    Comp1 = plsda_fit$variates$X[, 1],
    Comp2 = plsda_fit$variates$X[, 2],
    Outcome = factor(y, levels = c("Remission", "Recidivism"))
  )
  
  ggplot(scores_df, aes(x = Comp1, y = Comp2, fill = Outcome)) +
    geom_point(shape = 21, size = 3, color = "black", stroke = 0.6) +
    scale_fill_manual(values = c("Remission" = "#377EB8", "Recidivism" = "#E41A1C")) +
    theme_minimal(base_size = 14) +
    labs(
      title = title_text,
      x = "PLS-DA Component 1",
      y = "PLS-DA Component 2"
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    )
}

# ============================================================
# 4. Tumour-only cohort splits
# ============================================================

X_tcga <- expr_tumour[grepl("^TCGA", rownames(expr_tumour)), , drop = FALSE]
y_tcga <- factor(labels_tumour[grepl("^TCGA", rownames(expr_tumour))], levels = c("Remission", "Recidivism"))

X_ki <- expr_tumour[!grepl("^TCGA", rownames(expr_tumour)), , drop = FALSE]
y_ki <- factor(labels_tumour[!grepl("^TCGA", rownames(expr_tumour))], levels = c("Remission", "Recidivism"))

stopifnot(nrow(X_tcga) == length(y_tcga))
stopifnot(nrow(X_ki) == length(y_ki))

# ============================================================
# 4B. PLS-DA visualization (training cohorts)
# ============================================================

# TCGA
non_zero_tcga <- apply(X_tcga, 2, sd, na.rm = TRUE) > 0
X_tcga_pls <- X_tcga[, non_zero_tcga, drop = FALSE]
plsda_tcga <- plsda(X_tcga_pls, y_tcga, ncomp = 2)

p_plsda_tcga <- plot_plsda_scores(plsda_tcga, y_tcga, "SF-3F")
ggsave(file.path(figure_dir, "mirna_plsda_tcga_tumour.pdf"), p_plsda_tcga, width = 5, height = 4)

# KI
non_zero_ki <- apply(X_ki, 2, sd, na.rm = TRUE) > 0
X_ki_pls <- X_ki[, non_zero_ki, drop = FALSE]
plsda_ki <- plsda(X_ki_pls, y_ki, ncomp = 2)

p_plsda_ki <- plot_plsda_scores(plsda_ki, y_ki, "miRNA PLS-DA: KI tumour")
ggsave(file.path(figure_dir, "mirna_plsda_ki_tumour.pdf"), p_plsda_ki, width = 5, height = 4)

# ============================================================
# 5. Within-cohort cross-validation
# ============================================================

message("Running TCGA 5-fold CV...")
res_tcga_cv <- run_plsda_rf_cv(X_tcga, y_tcga, k = 5, ncomp = 2, n_top = 100, ntree = 500)
ranked_mirna_tcga_cv <- res_tcga_cv$ranked_features
p_roc_tcga_cv <- plot_roc_curve(
  res_tcga_cv$roc, res_tcga_cv$auc, res_tcga_cv$ci,
  "SF-3C", "springgreen4"
)
ggsave(file.path(figure_dir, "mirna_roc_tcga_cv.pdf"), p_roc_tcga_cv, width = 5, height = 4)

message("Running KI 5-fold CV...")
res_ki_cv <- run_plsda_rf_cv(X_ki, y_ki, k = 5, ncomp = 2, n_top = 100, ntree = 500)
ranked_mirna_ki_cv <- res_ki_cv$ranked_features
p_roc_ki_cv <- plot_roc_curve(
  res_ki_cv$roc, res_ki_cv$auc, res_ki_cv$ci,
  "SF-3B", "violetred4"
)
ggsave(file.path(figure_dir, "mirna_roc_ki_cv.pdf"), p_roc_ki_cv, width = 5, height = 4)

# ============================================================
# 6. Cross-cohort transfer analyses
# ============================================================

message("Running TCGA to KI transfer analysis...")
res_tcga_ki <- run_plsda_rf_transfer(X_tcga, y_tcga, X_ki, y_ki)
ranked_mirna_tcga_ki <- res_tcga_ki$ranked_features
p_roc_tcga_ki <- plot_roc_curve(
  res_tcga_ki$roc, res_tcga_ki$auc, res_tcga_ki$ci,
  "Fig. 2B-4", "springgreen4"
)
ggsave(file.path(figure_dir, "mirna_roc_tcga_to_ki.pdf"), p_roc_tcga_ki, width = 5, height = 4)

p_mdg_tcga_ki <- plot_mdg(
  res_tcga_ki$importance,
  "SF-3G",
  "springgreen4",
  top_n = 10
)
ggsave(file.path(figure_dir, "mirna_mdg_tcga_to_ki.pdf"), p_mdg_tcga_ki, width = 5, height = 4)

message("Running KI to TCGA transfer analysis...")
res_ki_tcga <- run_plsda_rf_transfer(X_ki, y_ki, X_tcga, y_tcga)
ranked_mirna_ki_tcga <- res_ki_tcga$ranked_features
p_roc_ki_tcga <- plot_roc_curve(
  res_ki_tcga$roc, res_ki_tcga$auc, res_ki_tcga$ci,
  "Fig. 2B-3", "violetred4"
)
ggsave(file.path(figure_dir, "mirna_roc_ki_to_tcga.pdf"), p_roc_ki_tcga, width = 5, height = 4)

p_mdg_ki_tcga <- plot_mdg(
  res_ki_tcga$importance,
  "SF-3E",
  "violetred4",
  top_n = 10
)
ggsave(file.path(figure_dir, "mirna_mdg_ki_to_tcga.pdf"), p_mdg_ki_tcga, width = 5, height = 4)

# Save summary metrics
auc_summary <- data.frame(
  Model = c("TCGA CV", "KI CV", "TCGA -> KI", "KI -> TCGA"),
  AUC = c(res_tcga_cv$auc, res_ki_cv$auc, res_tcga_ki$auc, res_ki_tcga$auc),
  CI_lower = c(
    res_tcga_cv$ci[1],
    res_ki_cv$ci[1],
    res_tcga_ki$ci[1],
    res_ki_tcga$ci[1]
  ),
  CI_upper = c(
    res_tcga_cv$ci[3],
    res_ki_cv$ci[3],
    res_tcga_ki$ci[3],
    res_ki_tcga$ci[3]
  ),
  stringsAsFactors = FALSE
)
write.csv(auc_summary, file.path(output_dir, "mirna_model_auc_summary.csv"), row.names = FALSE)

# ============================================================
# 7. Overlap dot plot
# ============================================================

top_n <- 20

mirna_lists <- list(
  "TCGA -> KI" = ranked_mirna_tcga_ki[1:min(top_n, length(ranked_mirna_tcga_ki))],
  "KI -> TCGA" = ranked_mirna_ki_tcga[1:min(top_n, length(ranked_mirna_ki_tcga))],
  "TCGA CV"    = ranked_mirna_tcga_cv[1:min(top_n, length(ranked_mirna_tcga_cv))],
  "KI CV"      = ranked_mirna_ki_cv[1:min(top_n, length(ranked_mirna_ki_cv))]
)

all_mirnas <- unique(unlist(mirna_lists))

dot_df <- expand.grid(
  miRNA = all_mirnas,
  Model = names(mirna_lists),
  stringsAsFactors = FALSE
)

dot_df$Rank <- mapply(function(g, m) {
  idx <- match(g, mirna_lists[[m]])
  ifelse(is.na(idx), NA, idx)
}, dot_df$miRNA, dot_df$Model)

dot_df$Present <- !is.na(dot_df$Rank)

mirna_counts <- dot_df %>%
  group_by(miRNA) %>%
  summarise(N_models = sum(Present), .groups = "drop")

dot_df <- dot_df %>%
  left_join(mirna_counts, by = "miRNA") %>%
  filter(N_models >= 2)

dot_df$PointSize <- ifelse(dot_df$Present, top_n - dot_df$Rank + 1, NA)

mirna_order <- mirna_counts %>%
  filter(N_models >= 2) %>%
  arrange(N_models, miRNA) %>%
  pull(miRNA)

dot_df$miRNA <- factor(dot_df$miRNA, levels = mirna_order)
dot_df$Model <- factor(dot_df$Model, levels = c("TCGA -> KI", "KI -> TCGA", "TCGA CV", "KI CV"))

p_dot <- ggplot(dot_df %>% filter(Present), aes(x = Model, y = miRNA, size = PointSize)) +
  geom_point(aes(fill = Model), shape = 21, color = "black", stroke = 0.3) +
  scale_size(range = c(2.5, 8), guide = guide_legend(title = "Higher rank")) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste("Shared top", top_n, "miRNAs across models"),
    x = "Model",
    y = "miRNA"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(figure_dir, "mirna_overlap_dotplot.pdf"), p_dot, width = 7, height = 6)

# ============================================================
# 8. Candidate-miRNA violin plot
# ============================================================

mirnas_of_interest <- c("hsa-miR-2355-5p")

missing_mirnas <- setdiff(mirnas_of_interest, colnames(expr_tumour))
if (length(missing_mirnas) > 0) {
  warning("The following candidate miRNAs were not found and will be skipped: ",
          paste(missing_mirnas, collapse = ", "))
}

mirnas_of_interest <- intersect(mirnas_of_interest, colnames(expr_tumour))

if (length(mirnas_of_interest) > 0) {
  plot_metadata <- data.frame(
    SampleID = rownames(expr_tumour),
    Cohort = ifelse(grepl("^TCGA", rownames(expr_tumour)), "TCGA", "KI"),
    Outcome = factor(labels_tumour, levels = c("Remission", "Recidivism")),
    stringsAsFactors = FALSE
  )
  
  expr_sub <- expr_tumour[, mirnas_of_interest, drop = FALSE]
  expr_sub <- data.frame(SampleID = rownames(expr_sub), expr_sub, stringsAsFactors = FALSE, check.names = FALSE)
  
  plot_long <- expr_sub %>%
    pivot_longer(cols = all_of(mirnas_of_interest), names_to = "miRNA", values_to = "Expression") %>%
    left_join(plot_metadata, by = "SampleID")
  
  p_mirna_violin <- ggplot(plot_long, aes(x = Outcome, y = Expression, fill = Outcome)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.4, color = "black") +
    facet_grid(miRNA ~ Cohort, scales = "free_y") +
    stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, size = 5) +
    scale_fill_manual(values = c("Remission" = "#377EB8", "Recidivism" = "#E41A1C")) +
    theme_minimal(base_size = 14) +
    labs(title = "Fig. 2C-4", x = "", y = "log2 expression") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
  ggsave(file.path(figure_dir, "mirna_candidate_violin_plot.pdf"), p_mirna_violin, width = 7, height = 6)
}
