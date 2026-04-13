# ============================================================
# Liver healthy-tissue recurrence analysis - mRNA
#
# Author: MAM
# Date: April 2026
#
# Description:
# This script reproduces the healthy tissue–based recurrence analyses
# presented in the manuscript, including preprocessing, cohort
# structure visualization, PLS-DA, repeated within-cohort cross-
# validation, cross-cohort transfer learning, pathway enrichment,
# and pathway-level program scoring.
#
# Input:
# - data/mRNA_merge_csv.csv
# - data/mRNA_miRNA_patient_sample_info.csv
#
# Output:
# - Figures (PDF) in results/figures/
# - Processed objects (.rds) in results/objects/
# - Summary tables (.csv) in results/
#
# Notes:
# - One predefined outlier sample (27.A033.303.Tumor) is excluded
#   based on QC and downstream inspection.
# - This script follows the same preprocessing and sample-matching
#   logic as the tumour and miRNA submission-ready scripts.
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
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(GSVA)
  library(GSEABase)
  library(AnnotationDbi)
  library(GO.db)
  library(rstatix)
  library(stringr)
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

# Predefined QC exclusion used in final manuscript analyses
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

# Use metadata column containing mRNA sample identifiers
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

# Restrict healthy-only dataset
healthy_samples <- names(tissue_labels)[tissue_labels == "Healthy"]
expr_healthy    <- log_expr_matrix[healthy_samples, , drop = FALSE]
labels_healthy  <- recurrence_labels[healthy_samples]

# Save processed objects
saveRDS(log_expr_matrix, file.path(object_dir, "healthy_log_expr_matrix.rds"))
saveRDS(recurrence_labels, file.path(object_dir, "healthy_recurrence_labels.rds"))
saveRDS(tissue_labels, file.path(object_dir, "healthy_tissue_labels.rds"))
saveRDS(expr_healthy, file.path(object_dir, "expr_healthy.rds"))
saveRDS(labels_healthy, file.path(object_dir, "labels_healthy.rds"))

# ============================================================
# 1B. Cohort summary table
# ============================================================

sample_summary <- data.frame(
  SampleID = rownames(log_expr_matrix),
  Cohort = ifelse(grepl("^TCGA", rownames(log_expr_matrix)), "TCGA", "KI"),
  Tissue = tissue_labels,
  Outcome = recurrence_labels,
  stringsAsFactors = FALSE
)

sample_summary$PatientID <- ifelse(
  sample_summary$Cohort == "TCGA",
  sapply(strsplit(gsub("\\.", "-", sample_summary$SampleID), "-"), function(x) paste(x[1:3], collapse = "-")),
  NA_character_
)

write.csv(sample_summary, file.path(output_dir, "healthy_sample_summary.csv"), row.names = FALSE)

sample_counts <- sample_summary %>%
  count(Cohort, Tissue, Outcome)
write.csv(sample_counts, file.path(output_dir, "healthy_sample_counts.csv"), row.names = FALSE)

paired_tcga_patients <- sample_summary %>%
  filter(Cohort == "TCGA", !is.na(PatientID)) %>%
  group_by(PatientID) %>%
  summarise(TissueSet = paste(sort(unique(as.character(Tissue))), collapse = "_"), .groups = "drop") %>%
  filter(TissueSet == "Healthy_Tumour")

write.csv(paired_tcga_patients, file.path(output_dir, "tcga_paired_patients.csv"), row.names = FALSE)

message("Initial TCGA samples in retained dataset: ", sum(sample_summary$Cohort == "TCGA"))
message("Unique TCGA patients: ", length(unique(na.omit(sample_summary$PatientID))))
message("TCGA patients with paired Healthy + Tumour samples: ", nrow(paired_tcga_patients))

# ============================================================
# 2. Utility functions for modelling and plotting
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

plot_auc_distribution <- function(auc_values, title_text, fill_color) {
  auc_df <- data.frame(AUC = auc_values)
  mean_auc <- mean(auc_values)
  sd_auc <- sd(auc_values)
  
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
        "Mean AUC = ", round(mean_auc, 3),
        "\nSD = ", round(sd_auc, 3),
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
        importance = TRUE,
        strata = y_train,
        sampsize = rep(min(table(y_train)), 2)
      )
      
      rf_importance_list[[length(rf_importance_list) + 1]] <- data.frame(
        Gene = rownames(rf_fit$importance),
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
  
  if (length(auc_values) == 0) {
    stop("No repeated CV AUC values were generated.")
  }
  if (length(rf_importance_list) == 0) {
    stop("No RF importance values were collected across repeats/folds.")
  }
  
  importance_df <- do.call(rbind, rf_importance_list)
  importance_summary <- aggregate(MeanDecreaseGini ~ Gene, data = importance_df, FUN = mean)
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
    ranked_genes = importance_summary$Gene
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
  
  stopifnot(nrow(X_train_full) == length(y_train))
  stopifnot(nrow(X_test_full) == length(y_test))
  
  nzv <- apply(X_train_full, 2, sd, na.rm = TRUE) > 0
  X_train_plsda <- X_train_full[, nzv, drop = FALSE]
  X_test_plsda  <- X_test_full[, nzv, drop = FALSE]
  
  if (ncol(X_train_plsda) < 2) {
    stop("Too few non-zero variance genes available for transfer analysis.")
  }
  
  colnames(X_train_plsda) <- make.unique(colnames(X_train_plsda))
  colnames(X_test_plsda)  <- colnames(X_train_plsda)
  
  plsda_fit <- plsda(X_train_plsda, y_train, ncomp = ncomp)
  
  load1 <- plsda_fit$loadings$X[, 1]
  load2 <- plsda_fit$loadings$X[, 2]
  
  top_c1 <- names(sort(abs(load1), decreasing = TRUE)[1:min(n_plsda, length(load1))])
  top_c2 <- names(sort(abs(load2), decreasing = TRUE)[1:min(n_plsda, length(load2))])
  top_plsda_genes <- union(top_c1, top_c2)
  
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
  
  list(
    roc = roc_obj,
    auc = as.numeric(auc(roc_obj)),
    ci = ci.auc(roc_obj),
    accuracy = mean(y_pred == y_test),
    confusion_matrix = table(True = y_test, Predicted = y_pred),
    importance = importance_df,
    ranked_genes = ranked_genes,
    selected_genes = selected_genes,
    plsda = plsda_fit
  )
}

get_top_loadings <- function(plsda_obj, model_name, top_n = 15) {
  load_mat <- as.matrix(plsda_obj$loadings$X[, 1:2, drop = FALSE])
  
  data.frame(
    Gene = rownames(load_mat),
    Comp1 = as.numeric(load_mat[, 1]),
    Comp2 = as.numeric(load_mat[, 2]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      LoadingScore = abs(Comp1) + abs(Comp2),
      Direction = ifelse(Comp1 >= 0, "Positive", "Negative"),
      Model = model_name
    ) %>%
    arrange(desc(LoadingScore)) %>%
    slice_head(n = top_n)
}

run_plsda_gsea <- function(plsda_fit, y_train, prefix) {
  y_train <- factor(y_train[rownames(plsda_fit$variates$X)], levels = c("Remission", "Recidivism"))
  scores <- plsda_fit$variates$X[, 1]
  mean_scores <- tapply(scores, y_train, mean, na.rm = TRUE)
  loadings_c1 <- plsda_fit$loadings$X[, 1]
  
  if (mean_scores["Recidivism"] > mean_scores["Remission"]) {
    gene_rank <- loadings_c1
  } else {
    gene_rank <- -loadings_c1
  }
  
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  gene_rank <- gene_rank[!is.na(gene_rank)]
  gene_rank <- gene_rank[!duplicated(names(gene_rank))]
  
  gsea_res <- gseGO(
    geneList = gene_rank,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
  gsea_df <- as.data.frame(gsea_res)
  recid_pathways <- gsea_df %>% filter(NES > 0, p.adjust < 0.05) %>% arrange(desc(NES))
  rem_pathways   <- gsea_df %>% filter(NES < 0, p.adjust < 0.05) %>% arrange(NES)
  
  write.csv(gsea_df, file.path(output_dir, paste0(prefix, "_gsea_go_bp.csv")), row.names = FALSE)
  write.csv(recid_pathways, file.path(output_dir, paste0(prefix, "_recidivism_pathways.csv")), row.names = FALSE)
  write.csv(rem_pathways, file.path(output_dir, paste0(prefix, "_remission_pathways.csv")), row.names = FALSE)
  
  list(gsea_df = gsea_df, recid_pathways = recid_pathways, rem_pathways = rem_pathways)
}

simplify_pathway_labels <- function(x) {
  case_when(
    x == "ribosome biogenesis" ~ "Ribosome biogenesis",
    x == "rRNA processing" ~ "rRNA processing",
    x == "rRNA metabolic process" ~ "rRNA metabolism",
    x == "ribonucleoprotein complex biogenesis" ~ "RNP complex biogenesis",
    x == "ribosomal small subunit biogenesis" ~ "Small ribosomal subunit biogenesis",
    x == "ribosomal large subunit biogenesis" ~ "Large ribosomal subunit biogenesis",
    x == "ncRNA processing" ~ "ncRNA processing",
    x == "tRNA metabolic process" ~ "tRNA metabolism",
    x == "MHC class II protein complex assembly" ~ "MHC II complex assembly",
    x == "peptide antigen assembly with MHC class II protein complex" ~ "MHC II antigen assembly",
    x == "immunoglobulin production" ~ "Immunoglobulin production",
    x == "immunoglobulin mediated immune response" ~ "Ig-mediated immune response",
    x == "B cell mediated immunity" ~ "B cell-mediated immunity",
    x == "production of molecular mediator of immune response" ~ "Immune mediator production",
    x == "lymphocyte mediated immunity" ~ "Lymphocyte-mediated immunity",
    x == "leukocyte mediated immunity" ~ "Leukocyte-mediated immunity",
    x == "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains" ~ "Adaptive immune response",
    TRUE ~ x
  )
}

extract_core <- function(gsea_df, pathways_keep) {
  gsea_df %>%
    filter(Description %in% pathways_keep) %>%
    select(Description, core_enrichment) %>%
    mutate(Gene = strsplit(core_enrichment, "/")) %>%
    tidyr::unnest(Gene)
}

# ============================================================
# 3. Healthy-only cohort splits
# ============================================================

X_tcga <- expr_healthy[grepl("^TCGA", rownames(expr_healthy)), , drop = FALSE]
y_tcga <- factor(labels_healthy[grepl("^TCGA", rownames(expr_healthy))], levels = c("Remission", "Recidivism"))

X_ki <- expr_healthy[!grepl("^TCGA", rownames(expr_healthy)), , drop = FALSE]
y_ki <- factor(labels_healthy[!grepl("^TCGA", rownames(expr_healthy))], levels = c("Remission", "Recidivism"))

stopifnot(nrow(X_tcga) == length(y_tcga))
stopifnot(nrow(X_ki) == length(y_ki))

# ============================================================
# 4. PLS-DA visualization (training cohorts)
# ============================================================

non_zero_tcga <- apply(X_tcga, 2, sd, na.rm = TRUE) > 0
X_tcga_pls <- X_tcga[, non_zero_tcga, drop = FALSE]
plsda_tcga <- plsda(X_tcga_pls, y_tcga, ncomp = 2)
p_plsda_tcga <- plot_plsda_scores(plsda_tcga, y_tcga, "Fig. 3B")
ggsave(file.path(figure_dir, "plsda_tcga_healthy.pdf"), p_plsda_tcga, width = 5, height = 4)

non_zero_ki <- apply(X_ki, 2, sd, na.rm = TRUE) > 0
X_ki_pls <- X_ki[, non_zero_ki, drop = FALSE]
plsda_ki <- plsda(X_ki_pls, y_ki, ncomp = 2)
p_plsda_ki <- plot_plsda_scores(plsda_ki, y_ki, "Fig. 3A")
ggsave(file.path(figure_dir, "plsda_ki_healthy.pdf"), p_plsda_ki, width = 5, height = 4)

# ============================================================
# 5. Within-cohort repeated cross-validation
# ============================================================

message("Running TCGA repeated 5-fold CV...")
res_tcga_cv <- run_plsda_rf_cv_repeated(X_tcga, y_tcga, n_repeats = 20, k = 5, ncomp = 2, n_top = 100, ntree = 500)
ranked_genes_tcga_cv <- res_tcga_cv$ranked_genes
p_auc_tcga_cv <- plot_auc_distribution(res_tcga_cv$auc_values, "SF-4D", "springgreen4")
ggsave(file.path(figure_dir, "auc_distribution_tcga_healthy_cv.pdf"), p_auc_tcga_cv, width = 5, height = 4)
p_roc_tcga_cv <- plot_roc_curve(res_tcga_cv$roc_last, res_tcga_cv$mean_auc, NULL, "SF-4C", "springgreen4")
ggsave(file.path(figure_dir, "roc_tcga_healthy_cv_last_repeat.pdf"), p_roc_tcga_cv, width = 5, height = 4)

message("Running KI repeated 5-fold CV...")
res_ki_cv <- run_plsda_rf_cv_repeated(X_ki, y_ki, n_repeats = 20, k = 5, ncomp = 2, n_top = 100, ntree = 500)
ranked_genes_ki_cv <- res_ki_cv$ranked_genes
p_auc_ki_cv <- plot_auc_distribution(res_ki_cv$auc_values, "SF-4B", "violetred4")
ggsave(file.path(figure_dir, "auc_distribution_ki_healthy_cv.pdf"), p_auc_ki_cv, width = 5, height = 4)
p_roc_ki_cv <- plot_roc_curve(res_ki_cv$roc_last, res_ki_cv$mean_auc, NULL, "SF-4A", "violetred4")
ggsave(file.path(figure_dir, "roc_ki_healthy_cv_last_repeat.pdf"), p_roc_ki_cv, width = 5, height = 4)

cv_summary <- data.frame(
  Model = c("TCGA repeated CV", "KI repeated CV"),
  Mean_AUC = c(res_tcga_cv$mean_auc, res_ki_cv$mean_auc),
  SD_AUC = c(res_tcga_cv$sd_auc, res_ki_cv$sd_auc),
  Repeats = c(length(res_tcga_cv$auc_values), length(res_ki_cv$auc_values)),
  stringsAsFactors = FALSE
)
write.csv(cv_summary, file.path(output_dir, "healthy_repeated_cv_summary.csv"), row.names = FALSE)

# ============================================================
# 6. Cross-cohort transfer analyses
# ============================================================

message("Running TCGA to KI transfer analysis...")
res_tcga_ki <- run_plsda_rf_transfer(X_tcga, y_tcga, X_ki, y_ki)
ranked_genes_tcga_ki <- res_tcga_ki$ranked_genes
p_roc_tcga_ki <- plot_roc_curve(res_tcga_ki$roc, res_tcga_ki$auc, res_tcga_ki$ci, "SF-4G", "springgreen4")
ggsave(file.path(figure_dir, "roc_tcga_to_ki_healthy.pdf"), p_roc_tcga_ki, width = 5, height = 4)
p_mdg_tcga_ki <- plot_mdg(res_tcga_ki$importance, "SF-4H", "springgreen4", top_n = 10)
ggsave(file.path(figure_dir, "mdg_tcga_to_ki_healthy.pdf"), p_mdg_tcga_ki, width = 5, height = 4)

message("Running KI to TCGA transfer analysis...")
res_ki_tcga <- run_plsda_rf_transfer(X_ki, y_ki, X_tcga, y_tcga)
ranked_genes_ki_tcga <- res_ki_tcga$ranked_genes
p_roc_ki_tcga <- plot_roc_curve(res_ki_tcga$roc, res_ki_tcga$auc, res_ki_tcga$ci, "SF-4E", "violetred4")
ggsave(file.path(figure_dir, "roc_ki_to_tcga_healthy.pdf"), p_roc_ki_tcga, width = 5, height = 4)
p_mdg_ki_tcga <- plot_mdg(res_ki_tcga$importance, "SF-4F", "violetred4", top_n = 10)
ggsave(file.path(figure_dir, "mdg_ki_to_tcga_healthy.pdf"), p_mdg_ki_tcga, width = 5, height = 4)

transfer_summary <- data.frame(
  Model = c("TCGA -> KI", "KI -> TCGA"),
  AUC = c(res_tcga_ki$auc, res_ki_tcga$auc),
  CI_lower = c(res_tcga_ki$ci[1], res_ki_tcga$ci[1]),
  CI_upper = c(res_tcga_ki$ci[3], res_ki_tcga$ci[3]),
  Accuracy = c(res_tcga_ki$accuracy, res_ki_tcga$accuracy),
  stringsAsFactors = FALSE
)
write.csv(transfer_summary, file.path(output_dir, "healthy_transfer_summary.csv"), row.names = FALSE)

# ============================================================
# 7. Top loading genes from cross-cohort PLS-DA models
# ============================================================

df_tcga_ki <- get_top_loadings(res_tcga_ki$plsda, "TCGA -> KI", top_n = 15)
df_ki_tcga <- get_top_loadings(res_ki_tcga$plsda, "KI -> TCGA", top_n = 15)

write.csv(df_tcga_ki, file.path(output_dir, "top_loadings_tcga_to_ki_healthy.csv"), row.names = FALSE)
write.csv(df_ki_tcga, file.path(output_dir, "top_loadings_ki_to_tcga_healthy.csv"), row.names = FALSE)

plot_top_loadings <- function(df, title_text) {
  df <- df %>% mutate(Gene = reorder(Gene, LoadingScore))
  
  ggplot(df, aes(x = LoadingScore, y = Gene, fill = Direction)) +
    geom_col(color = "black", width = 0.7) +
    scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) +
    theme_minimal(base_size = 14) +
    labs(
      title = title_text,
      x = "Combined loading score",
      y = NULL
    ) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}

p_load_tcga_ki <- plot_top_loadings(df_tcga_ki, "Fig. 4C-1")
p_load_ki_tcga <- plot_top_loadings(df_ki_tcga, "Fig. 4C-2")

ggsave(
  file.path(figure_dir, "top_loading_genes_cross_cohort_healthy.pdf"),
  p_load_tcga_ki + p_load_ki_tcga,
  width = 10,
  height = 6
)

# ============================================================
# 8. Pathway enrichment from cross-cohort PLS-DA models
# ============================================================

gsea_ki_tcga <- run_plsda_gsea(res_ki_tcga$plsda, y_ki, prefix = "ki_to_tcga_healthy")
gsea_tcga_ki <- run_plsda_gsea(res_tcga_ki$plsda, y_tcga, prefix = "tcga_to_ki_healthy")

shared_recid_pathways <- intersect(
  gsea_ki_tcga$recid_pathways$Description,
  gsea_tcga_ki$recid_pathways$Description
)

shared_rem_pathways <- intersect(
  gsea_ki_tcga$rem_pathways$Description,
  gsea_tcga_ki$rem_pathways$Description
)

write.csv(
  data.frame(Description = shared_recid_pathways),
  file.path(output_dir, "shared_recidivism_pathways.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Description = shared_rem_pathways),
  file.path(output_dir, "shared_remission_pathways.csv"),
  row.names = FALSE
)

shared_all <- union(shared_recid_pathways, shared_rem_pathways)

ki_df <- bind_rows(
  gsea_ki_tcga$recid_pathways %>% select(Description, NES),
  gsea_ki_tcga$rem_pathways %>% select(Description, NES)
) %>%
  filter(Description %in% shared_all) %>%
  mutate(Model = "KI -> TCGA")

tcga_df <- bind_rows(
  gsea_tcga_ki$recid_pathways %>% select(Description, NES),
  gsea_tcga_ki$rem_pathways %>% select(Description, NES)
) %>%
  filter(Description %in% shared_all) %>%
  mutate(Model = "TCGA -> KI")

plot_df_heat <- bind_rows(ki_df, tcga_df) %>%
  mutate(Description = simplify_pathway_labels(Description))

path_order <- plot_df_heat %>%
  group_by(Description) %>%
  summarise(meanNES = mean(NES), .groups = "drop") %>%
  arrange(meanNES) %>%
  pull(Description)

plot_df_heat$Description <- factor(plot_df_heat$Description, levels = path_order)

p_shared_heatmap <- ggplot(plot_df_heat, aes(x = Model, y = Description, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0
  ) +
  labs(
    title = "Fig. 3E",
    x = NULL,
    y = NULL,
    fill = "NES"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(figure_dir, "healthy_shared_pathway_heatmap.pdf"),
  p_shared_heatmap,
  width = 7,
  height = 8
)

# ============================================================
# 9. Recurrent leading-edge genes across shared pathways
# ============================================================

lead_recid_ki   <- extract_core(gsea_ki_tcga$gsea_df, shared_recid_pathways)
lead_rem_ki     <- extract_core(gsea_ki_tcga$gsea_df, shared_rem_pathways)
lead_recid_tcga <- extract_core(gsea_tcga_ki$gsea_df, shared_recid_pathways)
lead_rem_tcga   <- extract_core(gsea_tcga_ki$gsea_df, shared_rem_pathways)

recid_counts <- bind_rows(
  lead_recid_ki %>% mutate(Model = "KI_to_TCGA"),
  lead_recid_tcga %>% mutate(Model = "TCGA_to_KI")
) %>%
  group_by(Gene) %>%
  summarise(
    n_occurrences = n(),
    n_pathways = n_distinct(Description),
    n_models = n_distinct(Model),
    .groups = "drop"
  ) %>%
  filter(n_models == 2) %>%
  arrange(desc(n_occurrences), desc(n_pathways))

rem_counts <- bind_rows(
  lead_rem_ki %>% mutate(Model = "KI_to_TCGA"),
  lead_rem_tcga %>% mutate(Model = "TCGA_to_KI")
) %>%
  group_by(Gene) %>%
  summarise(
    n_occurrences = n(),
    n_pathways = n_distinct(Description),
    n_models = n_distinct(Model),
    .groups = "drop"
  ) %>%
  filter(n_models == 2) %>%
  arrange(desc(n_occurrences), desc(n_pathways))

write.csv(recid_counts, file.path(output_dir, "healthy_recurrent_leading_edge_recidivism.csv"), row.names = FALSE)
write.csv(rem_counts, file.path(output_dir, "healthy_recurrent_leading_edge_remission.csv"), row.names = FALSE)

top_recid_drivers <- recid_counts %>% slice_head(n = 10) %>% mutate(Gene = factor(Gene, levels = rev(Gene)))
p_recid_drivers <- ggplot(top_recid_drivers, aes(x = n_occurrences, y = Gene)) +
  geom_segment(aes(x = 0, xend = n_occurrences, y = Gene, yend = Gene), color = "grey70", linewidth = 1) +
  geom_point(size = 4, color = "#E41A1C") +
  theme_minimal(base_size = 12) +
  labs(title = "Top recurrent leading-edge genes: Recidivism, Fig. 3D", x = "Recurrent leading-edge genes", y = NULL)

top_rem_drivers <- rem_counts %>% slice_head(n = 10) %>% mutate(Gene = factor(Gene, levels = rev(Gene)))
p_rem_drivers <- ggplot(top_rem_drivers, aes(x = n_occurrences, y = Gene)) +
  geom_segment(aes(x = 0, xend = n_occurrences, y = Gene, yend = Gene), color = "grey70", linewidth = 1) +
  geom_point(size = 4, color = "#377EB8") +
  theme_minimal(base_size = 12) +
  labs(title = "Top recurrent leading-edge genes: Remission, SF-4I", x = "Recurrent leading-edge genes", y = NULL)

ggsave(file.path(figure_dir, "healthy_recurrent_leading_edge_genes.pdf"), p_recid_drivers + p_rem_drivers, width = 10, height = 5)

# ============================================================
# 10. Sample-level pathway scores with ssGSEA
# ============================================================

expr_mat <- t(expr_healthy)  # genes x samples
sample_info <- data.frame(
  Sample = colnames(expr_mat),
  Outcome = labels_healthy[colnames(expr_mat)],
  Cohort = ifelse(grepl("^TCGA", colnames(expr_mat)), "TCGA", "KI"),
  stringsAsFactors = FALSE
)

keep <- !is.na(sample_info$Outcome)
expr_mat <- expr_mat[, keep, drop = FALSE]
sample_info <- sample_info[keep, , drop = FALSE]
sample_info$Outcome <- factor(sample_info$Outcome, levels = c("Remission", "Recidivism"))

pathway_terms <- c(
  "ribosome biogenesis",
  "rRNA processing",
  "MHC class II protein complex assembly",
  "immunoglobulin production"
)

go_map <- AnnotationDbi::select(
  GO.db,
  keys = pathway_terms,
  columns = c("GOID", "TERM"),
  keytype = "TERM"
)
go_map <- unique(go_map[, c("GOID", "TERM")])
go_map <- go_map[!is.na(go_map$GOID), ]

gene_sets_list <- lapply(seq_len(nrow(go_map)), function(i) {
  genes_df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = go_map$GOID[i],
    columns = c("SYMBOL"),
    keytype = "GOALL"
  )
  unique(na.omit(genes_df$SYMBOL))
})
names(gene_sets_list) <- go_map$TERM

gene_sets_list <- lapply(gene_sets_list, function(gs) intersect(gs, rownames(expr_mat)))
gene_sets_list <- gene_sets_list[sapply(gene_sets_list, length) >= 10]

if (length(gene_sets_list) == 0) {
  stop("No pathway gene sets remained after filtering. Check gene names and GO mappings.")
}

ssgsea_scores <- gsva(
  expr = as.matrix(expr_mat),
  gset.idx.list = gene_sets_list,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = FALSE,
  verbose = FALSE
)

score_df <- as.data.frame(t(ssgsea_scores))
score_df$Sample <- rownames(score_df)
score_df <- left_join(score_df, sample_info, by = "Sample")

plot_df <- score_df %>%
  pivot_longer(cols = all_of(names(gene_sets_list)), names_to = "Pathway", values_to = "Score") %>%
  mutate(
    Pathway = recode(
      Pathway,
      "ribosome biogenesis" = "Ribosome biogenesis",
      "rRNA processing" = "rRNA processing",
      "MHC class II protein complex assembly" = "MHC II assembly",
      "immunoglobulin production" = "Immunoglobulin production"
    ),
    Pathway = factor(
      Pathway,
      levels = c("Ribosome biogenesis", "rRNA processing", "MHC II assembly", "Immunoglobulin production")
    )
  )

write.csv(plot_df, file.path(output_dir, "healthy_ssgsea_pathway_scores_long.csv"), row.names = FALSE)

p_all <- ggplot(plot_df, aes(x = Outcome, y = Score, fill = Outcome)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.7) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Remission" = "#377EB8", "Recidivism" = "#E41A1C")) +
  theme_minimal(base_size = 12) +
  labs(title = "SF-4J", x = NULL, y = "ssGSEA score") +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave(file.path(figure_dir, "healthy_ssgsea_pathway_scores.pdf"), p_all, width = 8, height = 6)

stats_df <- plot_df %>%
  group_by(Pathway) %>%
  wilcox_test(Score ~ Outcome) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
write.csv(stats_df, file.path(output_dir, "healthy_ssgsea_wilcox_stats.csv"), row.names = FALSE)

program_df <- plot_df %>%
  mutate(
    Program = case_when(
      Pathway %in% c("Ribosome biogenesis", "rRNA processing") ~ "Ribosome",
      Pathway %in% c("MHC II assembly", "Immunoglobulin production") ~ "Immune",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Program)) %>%
  group_by(Sample, Outcome, Cohort, Program) %>%
  summarise(Score = mean(Score), .groups = "drop") %>%
  pivot_wider(names_from = Program, values_from = Score)

write.csv(program_df, file.path(output_dir, "healthy_program_scores_mRNA.csv"), row.names = FALSE)

program_tcga <- program_df %>% filter(Cohort == "TCGA")
program_ki   <- program_df %>% filter(Cohort == "KI")

cor_tcga <- cor.test(program_tcga$Ribosome, program_tcga$Immune, method = "spearman")
cor_ki   <- cor.test(program_ki$Ribosome, program_ki$Immune, method = "spearman")

cor_summary <- data.frame(
  Cohort = c("TCGA", "KI"),
  Spearman_rho = c(unname(cor_tcga$estimate), unname(cor_ki$estimate)),
  P_value = c(cor_tcga$p.value, cor_ki$p.value),
  stringsAsFactors = FALSE
)
write.csv(cor_summary, file.path(output_dir, "healthy_program_correlations.csv"), row.names = FALSE)

make_program_corr_plot <- function(df, cor_obj, title_text, label_left = FALSE) {
  x_rng <- range(df$Ribosome, na.rm = TRUE)
  y_rng <- range(df$Immune, na.rm = TRUE)
  
  ggplot(df, aes(x = Ribosome, y = Immune, color = Outcome)) +
    geom_point(size = 3.2, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.1) +
    scale_color_manual(values = c("Remission" = "#377EB8", "Recidivism" = "#E41A1C")) +
    annotate(
      "text",
      x = if (label_left) x_rng[1] + 0.02 * diff(x_rng) else x_rng[2] - 0.02 * diff(x_rng),
      y = y_rng[2] - 0.02 * diff(y_rng),
      hjust = if (label_left) 0 else 1,
      vjust = 1,
      size = 4.5,
      label = paste0(
        "Spearman rho = ", round(unname(cor_obj$estimate), 2),
        "\np = ", signif(cor_obj$p.value, 2)
      )
    ) +
    labs(title = title_text, x = "Ribosome program score", y = "Immune program score") +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

p_tcga_corr <- make_program_corr_plot(program_tcga, cor_tcga, "TCGA, Fig. 3F", label_left = FALSE)
p_ki_corr   <- make_program_corr_plot(program_ki, cor_ki, "KI, SF-4K", label_left = TRUE)

ggsave(file.path(figure_dir, "healthy_program_correlation_tcga.pdf"), p_tcga_corr, width = 6, height = 5)
ggsave(file.path(figure_dir, "healthy_program_correlation_ki.pdf"), p_ki_corr, width = 6, height = 5)
