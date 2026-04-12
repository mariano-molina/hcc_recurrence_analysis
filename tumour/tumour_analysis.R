# ============================================================
# Liver tumour recurrence analysis
#
# Author: MAM
# Date: April 2026
#
# Description:
# This script reproduces the main tumour-based analyses presented in the manuscript,
# including preprocessing, cohort structure visualization (t-SNE),
# within-cohort cross-validation, cross-cohort prediction, and feature ranking.
#
# Input:
# - mRNA_merge_csv.csv
# - mRNA_miRNA_patient_sample_info.csv
#
# Output:
# - Figures (PDF) in results/figures/
# - Processed objects (.rds) in results/objects/
# - Summary tables (.csv) in results/
#
# Notes:
# - One predefined outlier sample (X27.A033.303.Tumor) is excluded
#   based on QC and PLSDA inspection.

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
expression_file <- file.path(input_dir, "mRNA_merge_csv.csv")
metadata_file   <- file.path(input_dir, "mRNA_miRNA_patient_sample_info.csv")

# Predefined QC exclusion used in final manuscript analyses (t-SNE and data inspection)
bad_sample <- "27.A033.303.Tumor"

# ============================================================
# 1. Load and preprocess data
# ============================================================

message("Loading expression matrix and metadata...")
expr_raw <- read.csv(expression_file, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)

# Remove genes with missing/empty identifiers
expr_clean <- expr_raw[!(is.na(expr_raw[[1]]) | expr_raw[[1]] == ""), ]

# Remove genes with zero expression across all samples
expr_clean <- expr_clean[rowSums(expr_clean[, -1] != 0) > 0, ]

# Extract sample IDs from expression matrix
sample_ids <- colnames(expr_clean)[-1]

# Use the metadata column containing mRNA sample identifiers
sample_id_column <- "mRNA"

# Standardize sample IDs before matching
sample_ids_clean <- trimws(sample_ids)
metadata_ids_clean <- trimws(as.character(metadata[[sample_id_column]]))

# Harmonize delimiters across datasets
sample_ids_clean <- gsub("-", ".", sample_ids_clean)
metadata_ids_clean <- gsub("-", ".", metadata_ids_clean)

# Add leading X to non-TCGA sample IDs to match metadata format
sample_ids_match <- ifelse(
  grepl("^TCGA", sample_ids_clean),
  sample_ids_clean,
  paste0("X", sample_ids_clean)
)

# Build lookup tables from metadata
recurrence_lookup <- setNames(metadata$Recidiv, metadata_ids_clean)
tissue_lookup     <- setNames(metadata$Tumour, metadata_ids_clean)

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

# Build samples x genes numeric matrix
expr_matrix <- apply(expr_clean[, -1], 2, as.numeric)
expr_matrix <- t(expr_matrix)
rownames(expr_matrix) <- sample_ids_clean
colnames(expr_matrix) <- make.unique(expr_clean[[1]])

# Remove samples with missing labels
valid_samples <- !is.na(recurrence_labels) & !is.na(tissue_labels)

if (sum(valid_samples) == 0) {
  stop("No samples remained after matching expression data to metadata. Check sample ID formatting and metadata sample ID column.")
}

expr_matrix <- expr_matrix[valid_samples, , drop = FALSE]
recurrence_labels <- factor(recurrence_labels[valid_samples], levels = c("Remission", "Recidivism"))
tissue_labels     <- factor(tissue_labels[valid_samples], levels = c("Healthy", "Tumour"))

# Remove zero-variance genes across retained samples
expr_matrix <- expr_matrix[, apply(expr_matrix, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]

# Log-transform expression values
mode(expr_matrix) <- "numeric"
log_expr_matrix <- log2(expr_matrix + 1)

# Assign names to label vectors
names(recurrence_labels) <- rownames(log_expr_matrix)
names(tissue_labels)     <- rownames(log_expr_matrix)

# Remove predefined outlier sample from final analyses
bad_sample_clean <- gsub("-", ".", bad_sample)

if (bad_sample_clean %in% rownames(log_expr_matrix)) {
  log_expr_matrix <- log_expr_matrix[rownames(log_expr_matrix) != bad_sample_clean, , drop = FALSE]
  recurrence_labels <- recurrence_labels[names(recurrence_labels) != bad_sample_clean]
  tissue_labels     <- tissue_labels[names(tissue_labels) != bad_sample_clean]
}

# Restrict tumour-only dataset
tumour_samples <- names(tissue_labels)[tissue_labels == "Tumour"]
expr_tumour    <- log_expr_matrix[tumour_samples, , drop = FALSE]
labels_tumour  <- recurrence_labels[tumour_samples]

# Save processed objects
saveRDS(log_expr_matrix, file.path(object_dir, "log_expr_matrix.rds"))
saveRDS(recurrence_labels, file.path(object_dir, "recurrence_labels.rds"))
saveRDS(tissue_labels, file.path(object_dir, "tissue_labels.rds"))
saveRDS(expr_tumour, file.path(object_dir, "expr_tumour.rds"))
saveRDS(labels_tumour, file.path(object_dir, "labels_tumour.rds"))

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
  tsne_df, "Outcome", "Fig. 2A-2",
  c("Recidivism" = "#E41A1C", "Remission" = "#377EB8")
)

p_tsne_tissue <- plot_tsne(
  tsne_df, "Status", "Fig. 2A-1",
  c("Tumour" = "orange", "Healthy" = "#4DAF4A")
)

p_tsne_source <- plot_tsne(
  tsne_df, "Source", "SF-2A",
  c("TCGA" = "#93C47D", "KI" = "#912F6C")
)

ggsave(file.path(figure_dir, "tumour_tsne_outcome.pdf"), p_tsne_outcome, width = 5, height = 4)
ggsave(file.path(figure_dir, "tumour_tsne_tissue.pdf"), p_tsne_tissue, width = 5, height = 4)
ggsave(file.path(figure_dir, "tumour_tsne_source.pdf"), p_tsne_source, width = 5, height = 4)

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
    
    ranked_genes <- names(sort(loading_score, decreasing = TRUE))
    selected_genes <- ranked_genes[1:min(n_top, length(ranked_genes))]
    selected_genes <- intersect(selected_genes, colnames(X_train_plsda))
    selected_genes <- intersect(selected_genes, colnames(X_test_plsda))
    
    if (length(selected_genes) < 2) next
    
    X_train <- X_train_plsda[, selected_genes, drop = FALSE]
    X_test  <- X_test_plsda[, selected_genes, drop = FALSE]
    
    rf_fit <- randomForest(
      x = X_train,
      y = y_train,
      ntree = ntree,
      importance = TRUE
    )
    
    rf_importance_list[[length(rf_importance_list) + 1]] <- data.frame(
      Gene = rownames(rf_fit$importance),
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
    MeanDecreaseGini ~ Gene,
    data = importance_df,
    FUN = mean
  )
  importance_summary <- importance_summary[order(-importance_summary$MeanDecreaseGini), ]
  
  list(
    roc = roc_obj,
    auc = auc_val,
    ci = ci_vals,
    importance = importance_summary,
    ranked_genes = importance_summary$Gene
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
    stop("Too few non-zero variance genes available for transfer analysis.")
  }
  
  colnames(X_train_plsda) <- make.unique(colnames(X_train_plsda))
  colnames(X_test_plsda)  <- colnames(X_train_plsda)
  
  plsda_fit <- plsda(X_train_plsda, y_train, ncomp = ncomp)
  
  loadings_comp1 <- plsda_fit$loadings$X[, 1]
  loadings_comp2 <- plsda_fit$loadings$X[, 2]
  
  top_genes_c1 <- names(sort(abs(loadings_comp1), decreasing = TRUE)[1:min(n_plsda, length(loadings_comp1))])
  top_genes_c2 <- names(sort(abs(loadings_comp2), decreasing = TRUE)[1:min(n_plsda, length(loadings_comp2))])
  top_plsda_genes <- union(top_genes_c1, top_genes_c2)
  
  X_train_rank <- X_train_plsda[, top_plsda_genes, drop = FALSE]
  X_test_rank  <- X_test_plsda[, top_plsda_genes, drop = FALSE]
  
  rf_rank <- randomForest(
    x = X_train_rank,
    y = y_train,
    importance = TRUE,
    ntree = ntree_rank
  )
  
  importance_df <- data.frame(
    Gene = rownames(rf_rank$importance),
    MeanDecreaseGini = rf_rank$importance[, "MeanDecreaseGini"],
    stringsAsFactors = FALSE
  )
  importance_df <- importance_df[order(-importance_df$MeanDecreaseGini), ]
  ranked_genes <- importance_df$Gene
  
  selected_genes <- ranked_genes[1:min(n_rf, length(ranked_genes))]
  selected_genes <- intersect(selected_genes, colnames(X_train_rank))
  selected_genes <- intersect(selected_genes, colnames(X_test_rank))
  
  if (length(selected_genes) < 2) {
    stop("Too few selected genes remained for transfer analysis.")
  }
  
  X_train <- X_train_rank[, selected_genes, drop = FALSE]
  X_test  <- X_test_rank[, selected_genes, drop = FALSE]
  
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
    ranked_genes = ranked_genes,
    selected_genes = selected_genes
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

# MDG plots

plot_mdg <- function(importance_df, title_text, fill_color, top_n = 10) {
  plot_df <- importance_df[order(-importance_df$MeanDecreaseGini), , drop = FALSE]
  plot_df <- plot_df[1:min(top_n, nrow(plot_df)), , drop = FALSE]
  plot_df$Gene <- factor(plot_df$Gene, levels = rev(plot_df$Gene))
  
  ggplot(plot_df, aes(x = Gene, y = MeanDecreaseGini)) +
    geom_col(fill = fill_color) +
    coord_flip() +
    theme_minimal(base_size = 14) +
    labs(
      title = title_text,
      x = "Gene",
      y = "Mean Decrease Gini"
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
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

# TCGA
non_zero_tcga <- apply(X_tcga, 2, sd, na.rm = TRUE) > 0
X_tcga_pls <- X_tcga[, non_zero_tcga, drop = FALSE]
plsda_tcga <- plsda(X_tcga_pls, y_tcga, ncomp = 2)

p_plsda_tcga <- plot_plsda_scores(plsda_tcga, y_tcga, "SF-2F")
ggsave(file.path(figure_dir, "plsda_tcga_tumour.pdf"), p_plsda_tcga, width = 5, height = 4)

# KI
non_zero_ki <- apply(X_ki, 2, sd, na.rm = TRUE) > 0
X_ki_pls <- X_ki[, non_zero_ki, drop = FALSE]
plsda_ki <- plsda(X_ki_pls, y_ki, ncomp = 2)

p_plsda_ki <- plot_plsda_scores(plsda_ki, y_ki, "SF-2D")
ggsave(file.path(figure_dir, "plsda_ki_tumour.pdf"), p_plsda_ki, width = 5, height = 4)

# ============================================================
# 5. Within-cohort cross-validation
# ============================================================

message("Running TCGA 5-fold CV...")
res_tcga_cv <- run_plsda_rf_cv(X_tcga, y_tcga, k = 5, ncomp = 2, n_top = 100, ntree = 500)
ranked_genes_tcga_cv <- res_tcga_cv$ranked_genes
p_roc_tcga_cv <- plot_roc_curve(res_tcga_cv$roc, res_tcga_cv$auc, res_tcga_cv$ci, "SF-2C",
  "springgreen4")
ggsave(file.path(figure_dir, "roc_tcga_cv.pdf"), p_roc_tcga_cv, width = 5, height = 4)

message("Running KI 5-fold CV...")
res_ki_cv <- run_plsda_rf_cv(X_ki, y_ki, k = 5, ncomp = 2, n_top = 100, ntree = 500)
ranked_genes_ki_cv <- res_ki_cv$ranked_genes
p_roc_ki_cv <- plot_roc_curve(res_ki_cv$roc, res_ki_cv$auc, res_ki_cv$ci, "SF-2B", "violetred4")
ggsave(file.path(figure_dir, "roc_ki_cv.pdf"), p_roc_ki_cv, width = 5, height = 4)

# ============================================================
# 6. Cross-cohort transfer analyses
# ============================================================

message("Running TCGA to KI transfer analysis...")
res_tcga_ki <- run_plsda_rf_transfer(X_tcga, y_tcga, X_ki, y_ki)
ranked_genes_tcga_ki <- res_tcga_ki$ranked_genes
p_roc_tcga_ki <- plot_roc_curve(res_tcga_ki$roc, res_tcga_ki$auc, res_tcga_ki$ci, "Fig. 2B-2", "springgreen4")
ggsave(file.path(figure_dir, "roc_tcga_to_ki.pdf"), p_roc_tcga_ki, width = 5, height = 4)

p_mdg_tcga_ki <- plot_mdg(
  res_tcga_ki$importance,
  "SF-2G",
  "springgreen4",
  top_n = 10
)

ggsave(file.path(figure_dir, "mdg_tcga_to_ki.pdf"), p_mdg_tcga_ki, width = 5, height = 4)

message("Running KI to TCGA transfer analysis...")
res_ki_tcga <- run_plsda_rf_transfer(X_ki, y_ki, X_tcga, y_tcga)
ranked_genes_ki_tcga <- res_ki_tcga$ranked_genes
p_roc_ki_tcga <- plot_roc_curve(res_ki_tcga$roc, res_ki_tcga$auc, res_ki_tcga$ci, "Fig. 2B-1", "violetred4")
ggsave(file.path(figure_dir, "roc_ki_to_tcga.pdf"), p_roc_ki_tcga, width = 5, height = 4)

p_mdg_ki_tcga <- plot_mdg(
  res_ki_tcga$importance,
  "SF-2E",
  "violetred4",
  top_n = 10
)

ggsave(file.path(figure_dir, "mdg_ki_to_tcga.pdf"), p_mdg_ki_tcga, width = 5, height = 4)

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
write.csv(auc_summary, file.path(output_dir, "tumour_model_auc_summary.csv"), row.names = FALSE)


# ============================================================
# 7. Candidate-gene boxplots
# ============================================================

genes_of_interest <- c("EFNA3", "NKAIN2")

plot_metadata <- data.frame(
  SampleID = rownames(expr_tumour),
  Cohort = ifelse(grepl("^TCGA", rownames(expr_tumour)), "TCGA", "KI"),
  Outcome = factor(labels_tumour, levels = c("Remission", "Recidivism")),
  stringsAsFactors = FALSE
)

expr_sub <- expr_tumour[, genes_of_interest, drop = FALSE]
expr_sub <- data.frame(SampleID = rownames(expr_sub), expr_sub, stringsAsFactors = FALSE)

plot_long <- expr_sub %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression") %>%
  left_join(plot_metadata, by = "SampleID")

p_genes <- ggplot(plot_long, aes(x = Outcome, y = Expression, fill = Outcome)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.4, color = "black") +
  facet_grid(Gene ~ Cohort, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, size = 5) +
  scale_fill_manual(values = c("Remission" = "#377EB8", "Recidivism" = "#E41A1C")) +
  theme_minimal(base_size = 14) +
  labs(title = "Expression of overlapping genes across cohorts", x = "", y = "log2 expression") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(figure_dir, "candidate_gene_boxplots.pdf"), p_genes, width = 7, height = 6)

# ============================================================
# 9. Candidate-gene ROC comparison
# ============================================================

y_all <- factor(labels_tumour, levels = c("Remission", "Recidivism"))

roc_efna3 <- roc(response = y_all, predictor = expr_tumour[, "EFNA3"], levels = c("Remission", "Recidivism"), direction = "<")
roc_nkain2 <- roc(response = y_all, predictor = expr_tumour[, "NKAIN2"], levels = c("Remission", "Recidivism"), direction = "<")

two_gene_df <- data.frame(
  Outcome = y_all,
  EFNA3 = expr_tumour[, "EFNA3"],
  NKAIN2 = expr_tumour[, "NKAIN2"]
)

glm_two_gene <- glm(Outcome ~ EFNA3 + NKAIN2, data = two_gene_df, family = binomial)
two_gene_prob <- predict(glm_two_gene, type = "response")
roc_two_gene <- roc(response = y_all, predictor = two_gene_prob, levels = c("Remission", "Recidivism"), direction = "<")

auc_df <- data.frame(
  Model = c("TCGA CV", "KI CV", "TCGA -> KI", "KI -> TCGA", "EFNA3", "NKAIN2", "EFNA3 + NKAIN2"),
  AUC = c(
    res_tcga_cv$auc,
    res_ki_cv$auc,
    res_tcga_ki$auc,
    res_ki_tcga$auc,
    as.numeric(auc(roc_efna3)),
    as.numeric(auc(roc_nkain2)),
    as.numeric(auc(roc_two_gene))
  ),
  stringsAsFactors = FALSE
)

auc_df$Model <- factor(auc_df$Model, levels = auc_df$Model[order(auc_df$AUC, decreasing = TRUE)])

p_auc <- ggplot(auc_df, aes(x = Model, y = AUC)) +
  geom_col(fill = "steelblue4", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.35, size = 5) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 14) +
  labs(title = "SF-2H", x = "", y = "AUC") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(file.path(figure_dir, "auc_comparison_across_models.pdf"), p_auc, width = 7, height = 5)
