# ============================================================
#  Liver tumour and adjacent tissue recurrence analysis
#  Uncorrected multiomics (mRNA + miRNA)
#
# Author: MAM
# Date: April 2026
#
# Description:
# This script reproduces Figure 4 in the manuscript.
#
# Input:
# - data/mRNA_merge_csv.csv
# - data/merge_mirna.csv
# - data/mRNA_miRNA_patient_sample_info.csv
#
# Output:
# - results/figures
# - results/objects/
# - results/tables/
#
# Notes:
# - Uses harmonized matched mRNA + miRNA samples
# - Uses raw log-transformed expression

# ============================================================

# ============================================================
# 0. Setup
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(limma)
  library(mixOmics)
  library(GSVA)
  library(GSEABase)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(clusterProfiler)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(multiMiR)
})

set.seed(999)

input_dir  <- "data"
output_dir  <- "results"
figure_dir  <- file.path(output_dir, "figures")
table_dir  <- file.path(output_dir, "tables")
object_dir <- file.path(output_dir, "objects")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(object_dir, showWarnings = FALSE, recursive = TRUE)

mrna_file     <- file.path(input_dir, "mRNA_merge_csv.csv")
mirna_file    <- file.path(input_dir, "merge_mirna.csv")
metadata_file <- file.path(input_dir, "mRNA_miRNA_patient_sample_info.csv")

bad_sample <- "X27.A033.303.Tumor"

# ============================================================
# 1. Helper functions
# ============================================================

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

get_patient_id <- function(x) {
  x <- gsub("_Tumor", "", x)
  x <- gsub("_Normal", "", x)
  x <- gsub("_01$", "", x)
  x <- gsub("_11$", "", x)
  x
}

format_r_p <- function(test_obj) {
  paste0(
    "R = ", round(unname(test_obj$estimate), 2),
    ", p = ", signif(test_obj$p.value, 2)
  )
}

safe_cor_test <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  }
  
  suppressWarnings(cor.test(x, y, method = method))
}

theme_fig4 <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold")
    )
}

# ============================================================
# 2. Load and preprocess mRNA data
# ============================================================

message("Loading mRNA data...")
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

message("Loading miRNA data...")
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
# 4. Match mRNA and miRNA samples
# ============================================================

common_samples <- intersect(rownames(mrna_expr), rownames(mirna_expr))

message("Matched multi-omics samples: ", length(common_samples))
if (length(common_samples) == 0) {
  stop("No overlapping samples remained between mRNA and miRNA matrices.")
}

log_mrna <- mrna_expr[common_samples, , drop = FALSE]
mirna_subset <- mirna_expr[common_samples, , drop = FALSE]

recidiv_multi <- mrna_recurrence_labels[common_samples]
tumour_multi  <- mrna_tissue_labels[common_samples]
source_multi  <- factor(ifelse(grepl("^TCGA", common_samples), "TCGA", "KI"),
                        levels = c("KI", "TCGA"))
names(source_multi) <- common_samples

write.csv(
  data.frame(
    Sample = common_samples,
    Cohort = source_multi,
    Tissue = tumour_multi,
    Outcome = recidiv_multi
  ),
  file.path(table_dir, "sample_summary.csv"),
  row.names = FALSE
)

saveRDS(log_mrna, file.path(object_dir, "log_mrna.rds"))
saveRDS(mirna_subset, file.path(object_dir, "log_mirna.rds"))

# ============================================================
# 5. Pathway programs and tumour-adjacent correlations
# ============================================================

message("Running GSVA...")

expr_mat <- t(log_mrna)  # genes x samples

sample_info <- data.frame(
  Sample = colnames(expr_mat),
  TumourStatus = tumour_multi[colnames(expr_mat)],
  Outcome = recidiv_multi[colnames(expr_mat)],
  Cohort = source_multi[colnames(expr_mat)],
  stringsAsFactors = FALSE
)

keep <- !is.na(sample_info$TumourStatus) & !is.na(sample_info$Outcome)
expr_mat <- expr_mat[, keep, drop = FALSE]
sample_info <- sample_info[keep, , drop = FALSE]

sample_info$TumourStatus <- factor(sample_info$TumourStatus, levels = c("Healthy", "Tumour"))
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
  go_id <- go_map$GOID[i]
  genes_df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = go_id,
    columns = c("SYMBOL"),
    keytype = "GOALL"
  )
  unique(na.omit(genes_df$SYMBOL))
})
names(gene_sets_list) <- go_map$TERM
gene_sets_list <- lapply(gene_sets_list, function(gs) intersect(gs, rownames(expr_mat)))
gene_sets_list <- gene_sets_list[sapply(gene_sets_list, length) >= 10]

if (length(gene_sets_list) == 0) {
  stop("No pathway gene sets remained after filtering.")
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

ribosome_cols <- intersect(c("ribosome biogenesis", "rRNA processing"), colnames(score_df))
immune_cols   <- intersect(c("MHC class II protein complex assembly", "immunoglobulin production"), colnames(score_df))

if (length(ribosome_cols) == 0) stop("No ribosome pathways found in score_df.")
if (length(immune_cols) == 0) stop("No immune pathways found in score_df.")

score_df$RibosomeProgram <- rowMeans(score_df[, ribosome_cols, drop = FALSE], na.rm = TRUE)
score_df$ImmuneProgram   <- rowMeans(score_df[, immune_cols, drop = FALSE], na.rm = TRUE)

tumor_df <- score_df %>%
  filter(TumourStatus == "Tumour") %>%
  mutate(Patient = get_patient_id(Sample))

adj_df <- score_df %>%
  filter(TumourStatus == "Healthy") %>%
  mutate(Patient = get_patient_id(Sample))

common_patients <- intersect(tumor_df$Patient, adj_df$Patient)
message("Matched tumour-adjacent pairs: ", length(common_patients))

tumor_match <- tumor_df[match(common_patients, tumor_df$Patient), ]
adj_match   <- adj_df[match(common_patients, adj_df$Patient), ]

stopifnot(all(tumor_match$Patient == adj_match$Patient))

corr_df <- data.frame(
  Patient = common_patients,
  Outcome = factor(tumor_match$Outcome, levels = c("Remission", "Recidivism")),
  Cohort = tumor_match$Cohort,
  tumor_ribo = tumor_match$RibosomeProgram,
  adj_ribo   = adj_match$RibosomeProgram,
  tumor_imm  = tumor_match$ImmuneProgram,
  adj_imm    = adj_match$ImmuneProgram,
  stringsAsFactors = FALSE
)

write.csv(corr_df, file.path(table_dir, "correlation_input.csv"), row.names = FALSE)

cor_ribo_rem <- safe_cor_test(
  corr_df$tumor_ribo[corr_df$Outcome == "Remission"],
  corr_df$adj_ribo[corr_df$Outcome == "Remission"]
)
cor_ribo_rec <- safe_cor_test(
  corr_df$tumor_ribo[corr_df$Outcome == "Recidivism"],
  corr_df$adj_ribo[corr_df$Outcome == "Recidivism"]
)

cor_imm_rem <- safe_cor_test(
  corr_df$tumor_imm[corr_df$Outcome == "Remission"],
  corr_df$adj_imm[corr_df$Outcome == "Remission"]
)
cor_imm_rec <- safe_cor_test(
  corr_df$tumor_imm[corr_df$Outcome == "Recidivism"],
  corr_df$adj_imm[corr_df$Outcome == "Recidivism"]
)

lab_ribo_rem <- format_r_p(cor_ribo_rem)
lab_ribo_rec <- format_r_p(cor_ribo_rec)
lab_imm_rem  <- format_r_p(cor_imm_rem)
lab_imm_rec  <- format_r_p(cor_imm_rec)

x_ribo_lab <- max(corr_df$tumor_ribo, na.rm = TRUE) - 0.01 * diff(range(corr_df$tumor_ribo, na.rm = TRUE))
y_ribo_lab1 <- min(corr_df$adj_ribo, na.rm = TRUE) + 0.12 * diff(range(corr_df$adj_ribo, na.rm = TRUE))
y_ribo_lab2 <- min(corr_df$adj_ribo, na.rm = TRUE) + 0.06 * diff(range(corr_df$adj_ribo, na.rm = TRUE))

x_imm_lab <- max(corr_df$tumor_imm, na.rm = TRUE) - 0.01 * diff(range(corr_df$tumor_imm, na.rm = TRUE))
y_imm_lab1 <- min(corr_df$adj_imm, na.rm = TRUE) + 0.12 * diff(range(corr_df$adj_imm, na.rm = TRUE))
y_imm_lab2 <- min(corr_df$adj_imm, na.rm = TRUE) + 0.06 * diff(range(corr_df$adj_imm, na.rm = TRUE))

p_ribo <- ggplot(corr_df, aes(x = tumor_ribo, y = adj_ribo)) +
  geom_point(aes(color = Outcome), size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
  annotate("text", x = x_ribo_lab, y = y_ribo_lab1, label = lab_ribo_rem,
           hjust = 1, color = "#377EB8", size = 4) +
  annotate("text", x = x_ribo_lab, y = y_ribo_lab2, label = lab_ribo_rec,
           hjust = 1, color = "#E41A1C", size = 4) +
  scale_color_manual(values = c("Remission" = "#377EB8", "Recidivism" = "#E41A1C")) +
  labs(
    title = "",
    x = "Tumor ribosome program score",
    y = "Adjacent tissue ribosome program score"
  ) +
  theme_fig4(14) +
  theme(legend.position = "none")

p_immune <- ggplot(corr_df, aes(x = tumor_imm, y = adj_imm)) +
  geom_point(aes(color = Outcome), size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
  annotate("text", x = x_imm_lab, y = y_imm_lab1, label = lab_imm_rem,
           hjust = 1, color = "#377EB8", size = 4) +
  annotate("text", x = x_imm_lab, y = y_imm_lab2, label = lab_imm_rec,
           hjust = 1, color = "#E41A1C", size = 4) +
  scale_color_manual(values = c("Remission" = "#377EB8", "Recidivism" = "#E41A1C")) +
  labs(
    title = "",
    x = "Tumor immune program score",
    y = "Adjacent tissue immune program score"
  ) +
  theme_fig4(14) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "ribosome_correlation.pdf"), p_ribo, width = 5, height = 4)
ggsave(file.path(output_dir, "immune_correlation.pdf"), p_immune, width = 5, height = 4)

# ============================================================
# 6. Integrated multi-omic network
# ============================================================

message("Building network...")

tumor_samples <- names(tumour_multi)[tumour_multi == "Tumour"]
adj_samples   <- names(tumour_multi)[tumour_multi == "Healthy"]

tumor_patients <- get_patient_id(tumor_samples)
adj_patients   <- get_patient_id(adj_samples)

common_patients_net <- intersect(tumor_patients, adj_patients)

tumor_idx <- match(common_patients_net, tumor_patients)
adj_idx   <- match(common_patients_net, adj_patients)

tumor_ids_matched <- tumor_samples[tumor_idx]
adj_ids_matched   <- adj_samples[adj_idx]

tumor_mrna <- log_mrna[tumor_ids_matched, , drop = FALSE]
adj_mrna   <- log_mrna[adj_ids_matched, , drop = FALSE]
adj_mirna  <- mirna_subset[adj_ids_matched, , drop = FALSE]

rownames(tumor_mrna) <- common_patients_net
rownames(adj_mrna)   <- common_patients_net
rownames(adj_mirna)  <- common_patients_net

network_ribosome_genes <- c("NPM1", "RIOK2", "EXOSC2", "NOL6", "EIF4A3")
network_immune_genes   <- c("CCR6", "CD28", "IL2", "TBX21", "HLA-DRA")
key_mirnas <- c("hsa-miR-5699-3p", "hsa-miR-5589-3p", "hsa-miR-301a-5p", "hsa-miR-1-3p")

network_ribosome_genes <- intersect(network_ribosome_genes, colnames(adj_mrna))
network_immune_genes   <- intersect(network_immune_genes, colnames(adj_mrna))
key_mirnas             <- intersect(key_mirnas, colnames(adj_mirna))

if (length(network_ribosome_genes) == 0) stop("No ribosome genes available for network.")
if (length(network_immune_genes) == 0) stop("No immune genes available for network.")
if (length(key_mirnas) == 0) stop("No selected miRNAs available for network.")

adj_ribo_score   <- rowMeans(adj_mrna[, network_ribosome_genes, drop = FALSE], na.rm = TRUE)
adj_immune_score <- rowMeans(adj_mrna[, network_immune_genes, drop = FALSE], na.rm = TRUE)
tumor_ribo_score <- rowMeans(tumor_mrna[, network_ribosome_genes, drop = FALSE], na.rm = TRUE)
tumor_immune_score <- rowMeans(tumor_mrna[, network_immune_genes, drop = FALSE], na.rm = TRUE)

net_df <- data.frame(
  Patient = common_patients_net,
  Adj_Ribosome   = adj_ribo_score,
  Adj_Immune     = adj_immune_score,
  Tumor_Ribosome = tumor_ribo_score,
  Tumor_Immune   = tumor_immune_score,
  stringsAsFactors = FALSE
)

for (g in network_ribosome_genes) {
  net_df[[paste0("Adj_", g)]] <- adj_mrna[, g]
}
for (g in network_immune_genes) {
  net_df[[paste0("Adj_", g)]] <- adj_mrna[, g]
}
for (m in key_mirnas) {
  net_df[[m]] <- adj_mirna[, m]
}

num_df <- net_df[, setdiff(colnames(net_df), "Patient"), drop = FALSE]
vars <- colnames(num_df)

cor_mat <- matrix(NA_real_, nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
p_mat   <- matrix(NA_real_, nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))

for (i in seq_along(vars)) {
  for (j in seq_along(vars)) {
    test <- safe_cor_test(num_df[[vars[i]]], num_df[[vars[j]]], method = "spearman")
    cor_mat[i, j] <- unname(test$estimate)
    p_mat[i, j]   <- test$p.value
  }
}

edge_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
colnames(edge_df) <- c("from", "to", "rho")

p_df <- as.data.frame(as.table(p_mat), stringsAsFactors = FALSE)
colnames(p_df) <- c("from", "to", "p_value")

edge_df <- edge_df %>%
  left_join(p_df, by = c("from", "to")) %>%
  filter(from != to) %>%
  mutate(pair = paste(pmin(from, to), pmax(from, to), sep = "__")) %>%
  distinct(pair, .keep_all = TRUE)

main_edges <- edge_df %>%
  filter(
    abs(rho) >= 0.40,
    !is.na(p_value),
    p_value < 0.05,
    (
      (from %in% c("Adj_Ribosome", "Adj_Immune") & grepl("^Adj_", to)) |
        (to %in% c("Adj_Ribosome", "Adj_Immune") & grepl("^Adj_", from)) |
        (from %in% c("Tumor_Ribosome", "Tumor_Immune") & to %in% c("Adj_Ribosome", "Adj_Immune")) |
        (to %in% c("Tumor_Ribosome", "Tumor_Immune") & from %in% c("Adj_Ribosome", "Adj_Immune"))
    )
  )

mirna_edges <- edge_df %>%
  filter(
    abs(rho) >= 0.30,
    !is.na(p_value),
    p_value < 0.05,
    (
      (grepl("^hsa-", from) & grepl("^Adj_", to)) |
        (grepl("^hsa-", to) & grepl("^Adj_", from)) |
        (grepl("^hsa-", from) & to %in% c("Adj_Ribosome", "Adj_Immune")) |
        (grepl("^hsa-", to) & from %in% c("Adj_Ribosome", "Adj_Immune"))
    )
  )

forced_edges <- data.frame(
  from = c("Tumor_Ribosome", "Tumor_Immune"),
  to   = c("Adj_Ribosome",   "Adj_Immune"),
  rho  = c(
    suppressWarnings(cor(num_df$Tumor_Ribosome, num_df$Adj_Ribosome, method = "spearman", use = "pairwise.complete.obs")),
    suppressWarnings(cor(num_df$Tumor_Immune, num_df$Adj_Immune, method = "spearman", use = "pairwise.complete.obs"))
  ),
  p_value = c(
    safe_cor_test(num_df$Tumor_Ribosome, num_df$Adj_Ribosome)$p.value,
    safe_cor_test(num_df$Tumor_Immune, num_df$Adj_Immune)$p.value
  ),
  stringsAsFactors = FALSE
)

edge_df_keep <- bind_rows(main_edges, mirna_edges, forced_edges) %>%
  mutate(pair = paste(pmin(from, to), pmax(from, to), sep = "__")) %>%
  distinct(pair, .keep_all = TRUE)

write.csv(edge_df_keep, file.path(table_dir, "network_edges.csv"), row.names = FALSE)

node_names <- unique(c(edge_df_keep$from, edge_df_keep$to))
node_df <- data.frame(name = node_names, stringsAsFactors = FALSE)

node_df$layer <- dplyr::case_when(
  node_df$name %in% c("Adj_Ribosome", "Adj_Immune") ~ "Adjacent program",
  node_df$name %in% c("Tumor_Ribosome", "Tumor_Immune") ~ "Tumor program",
  grepl("^Adj_", node_df$name) ~ "Adjacent gene",
  grepl("^hsa-", node_df$name) ~ "miRNA",
  TRUE ~ "Other"
)

node_df$size <- dplyr::case_when(
  node_df$layer %in% c("Adjacent program", "Tumor program") ~ 13,
  node_df$layer == "miRNA" ~ 8,
  TRUE ~ 6
)

node_df$label <- node_df$name
node_df$label[node_df$label == "Adj_Ribosome"]   <- "Adjacent ribosome"
node_df$label[node_df$label == "Adj_Immune"]     <- "Adjacent immune"
node_df$label[node_df$label == "Tumor_Ribosome"] <- "Tumor ribosome"
node_df$label[node_df$label == "Tumor_Immune"]   <- "Tumor immune"
node_df$label <- gsub("^Adj_", "", node_df$label)

write.csv(node_df, file.path(table_dir, "network_nodes.csv"), row.names = FALSE)

g <- graph_from_data_frame(
  d = edge_df_keep[, c("from", "to", "rho", "p_value")],
  vertices = node_df,
  directed = FALSE
)

set.seed(123)
p_network <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = abs(rho), color = rho), alpha = 0.80) +
  scale_edge_width(range = c(0.8, 2.5)) +
  scale_edge_color_gradient2(
    low = "#377EB8",
    mid = "grey80",
    high = "#E41A1C",
    midpoint = 0
  ) +
  geom_node_point(aes(size = size, color = layer), alpha = 0.95) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3.2) +
  scale_size_identity() +
  scale_color_manual(values = c(
    "Adjacent program" = "#6BAED6",
    "Tumor program"    = "#FB6A4A",
    "Adjacent gene"    = "#9ECAE1",
    "miRNA"            = "#74C476",
    "Other"            = "grey60"
  )) +
  theme_void(base_size = 13) +
  labs(
    title = "C",
    edge_color = "Spearman rho",
    edge_width = "|rho|",
    color = "Node type"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "network.pdf"), p_network, width = 7, height = 5.5)

# ============================================================
# 7. Pathway convergence across omics
# ============================================================

message("Running pathway convergence...")

adj_samples <- names(tumour_multi)[tumour_multi == "Healthy"]

expr_adj_mrna <- log_mrna[adj_samples, , drop = FALSE]
expr_adj_mirna <- mirna_subset[adj_samples, , drop = FALSE]

sample_info_adj <- data.frame(
  Sample = adj_samples,
  Outcome = recidiv_multi[adj_samples],
  Cohort = ifelse(grepl("^TCGA", adj_samples), "TCGA", "KI"),
  stringsAsFactors = FALSE
)

keep <- !is.na(sample_info_adj$Outcome)
sample_info_adj <- sample_info_adj[keep, , drop = FALSE]
sample_info_adj$Outcome <- factor(sample_info_adj$Outcome, levels = c("Remission", "Recidivism"))

expr_adj_mrna <- expr_adj_mrna[sample_info_adj$Sample, , drop = FALSE]
expr_adj_mirna <- expr_adj_mirna[sample_info_adj$Sample, , drop = FALSE]

run_limma_by_cohort <- function(expr_mat, sample_info, cohort_name) {
  idx <- sample_info$Cohort == cohort_name
  sub_info <- sample_info[idx, , drop = FALSE]
  sub_expr <- expr_mat[sub_info$Sample, , drop = FALSE]
  
  if (nrow(sub_expr) < 3 || length(unique(sub_info$Outcome)) < 2) {
    stop(paste("Not enough samples/classes for cohort:", cohort_name))
  }
  
  design <- model.matrix(~ sub_info$Outcome)
  fit <- limma::lmFit(t(sub_expr), design)
  fit <- limma::eBayes(fit)
  
  res <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res$Feature <- rownames(res)
  res$Cohort <- cohort_name
  rownames(res) <- NULL
  res
}

mrna_res_ki <- run_limma_by_cohort(expr_adj_mrna, sample_info_adj, "KI") %>%
  dplyr::arrange(P.Value) %>%
  dplyr::mutate(rank_ki = dplyr::row_number()) %>%
  dplyr::select(Feature, logFC_ki = logFC, P_ki = P.Value, rank_ki)

mrna_res_tcga <- run_limma_by_cohort(expr_adj_mrna, sample_info_adj, "TCGA") %>%
  dplyr::arrange(P.Value) %>%
  dplyr::mutate(rank_tcga = dplyr::row_number()) %>%
  dplyr::select(Feature, logFC_tcga = logFC, P_tcga = P.Value, rank_tcga)

mrna_rank_df <- dplyr::full_join(mrna_res_ki, mrna_res_tcga, by = "Feature") %>%
  dplyr::mutate(
    score_ki = ifelse(!is.na(logFC_ki) & !is.na(P_ki), sign(logFC_ki) * -log10(P_ki + 1e-300), NA_real_),
    score_tcga = ifelse(!is.na(logFC_tcga) & !is.na(P_tcga), sign(logFC_tcga) * -log10(P_tcga + 1e-300), NA_real_)
  ) %>%
  dplyr::mutate(combined_score = rowMeans(cbind(score_ki, score_tcga), na.rm = TRUE)) %>%
  dplyr::filter(is.finite(combined_score)) %>%
  dplyr::filter(grepl("^[A-Z0-9]+$", Feature))

mrna_rank_mapped <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(mrna_rank_df$Feature),
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "SYMBOL"
) %>%
  dplyr::filter(!is.na(ENTREZID), !is.na(SYMBOL)) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::inner_join(mrna_rank_df, by = c("SYMBOL" = "Feature")) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_max(order_by = abs(combined_score), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

geneList_mrna <- mrna_rank_mapped$combined_score
names(geneList_mrna) <- mrna_rank_mapped$ENTREZID
geneList_mrna <- sort(geneList_mrna, decreasing = TRUE)

gsea_mrna <- clusterProfiler::gseGO(
  geneList = geneList_mrna,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  verbose = FALSE
)

mrna_enrich_df <- as.data.frame(gsea_mrna)

key_mirnas_panel <- c(
  "hsa-miR-5699-3p", "hsa-miR-6723-5p", "hsa-miR-376b-5p", "hsa-miR-4705",
  "hsa-miR-370-5p", "hsa-miR-203b-3p", "hsa-miR-6516-5p", "hsa-miR-4767",
  "hsa-miR-5589-3p", "hsa-miR-4685-3p", "hsa-miR-487a-3p", "hsa-miR-642a-5p",
  "hsa-miR-301a-5p", "hsa-miR-5590-3p", "hsa-miR-570-3p", "hsa-miR-1-3p"
)

key_mirnas_panel <- intersect(key_mirnas_panel, colnames(expr_adj_mirna))
if (length(key_mirnas_panel) == 0) {
  stop("No miRNAs available for pathway convergence panel.")
}

safe_get_multimir_one <- function(mirna_name, tbl = "validated") {
  tryCatch({
    res <- multiMiR::get_multimir(
      mirna = mirna_name,
      table = tbl,
      summary = FALSE,
      legacy.out = FALSE
    )
    df <- res@data
    if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) return(NULL)
    df$queried_mirna <- mirna_name
    df
  }, error = function(e) {
    message("multiMiR failed for ", mirna_name, ": ", e$message)
    NULL
  })
}

mirna_target_list <- lapply(key_mirnas_panel, safe_get_multimir_one, tbl = "validated")
names(mirna_target_list) <- key_mirnas_panel
mirna_target_list <- mirna_target_list[!sapply(mirna_target_list, is.null)]

if (length(mirna_target_list) == 0) {
  stop("No validated targets found for queried miRNAs.")
}

mirna_targets_df <- dplyr::bind_rows(mirna_target_list)

if ("target_symbol" %in% colnames(mirna_targets_df)) {
  mirna_target_genes <- mirna_targets_df %>%
    dplyr::filter(!is.na(target_symbol), target_symbol != "") %>%
    dplyr::pull(target_symbol) %>%
    unique()
} else if ("target_entrez" %in% colnames(mirna_targets_df)) {
  entrez_ids <- mirna_targets_df %>%
    dplyr::filter(!is.na(target_entrez)) %>%
    dplyr::pull(target_entrez) %>%
    as.character() %>%
    unique()
  
  mapped_targets <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = entrez_ids,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "ENTREZID"
  ) %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
  
  mirna_target_genes <- unique(mapped_targets$SYMBOL)
} else {
  stop("Could not find target_symbol or target_entrez in multiMiR output.")
}

mirna_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = mirna_target_genes,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "SYMBOL"
) %>%
  dplyr::filter(!is.na(ENTREZID), !is.na(SYMBOL)) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

ego_mirna <- clusterProfiler::enrichGO(
  gene = mirna_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

mirna_enrich_df <- as.data.frame(ego_mirna)

pathway_patterns <- c(
  "ribosom", "rRNA", "RNA", "translat", "splice",
  "immune", "lymph", "antigen", "interferon", "MHC"
)

keep_pathways_gsea <- function(df, source_name) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  df %>%
    dplyr::filter(grepl(paste(pathway_patterns, collapse = "|"), Description, ignore.case = TRUE)) %>%
    dplyr::mutate(
      Source = source_name,
      neglog10FDR = -log10(p.adjust),
      Count = setSize
    ) %>%
    dplyr::select(Description, Count, p.adjust, neglog10FDR, Source)
}

keep_pathways_ora <- function(df, source_name) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  df %>%
    dplyr::filter(grepl(paste(pathway_patterns, collapse = "|"), Description, ignore.case = TRUE)) %>%
    dplyr::mutate(
      Source = source_name,
      neglog10FDR = -log10(p.adjust)
    ) %>%
    dplyr::select(Description, Count, p.adjust, neglog10FDR, Source)
}

plot_mrna <- keep_pathways_gsea(mrna_enrich_df, "mRNA")
plot_mirna <- keep_pathways_ora(mirna_enrich_df, "miRNA targets")

combined_plot_df <- dplyr::bind_rows(plot_mrna, plot_mirna)

combined_plot_df$PathwayGroup <- dplyr::case_when(
  grepl("ribosom|translat", combined_plot_df$Description, ignore.case = TRUE) ~ "Ribosome / translation",
  grepl("rRNA|RNA|splice", combined_plot_df$Description, ignore.case = TRUE) ~ "RNA processing",
  grepl("immune|lymph|antigen|interferon|MHC", combined_plot_df$Description, ignore.case = TRUE) ~ "Immune signaling",
  TRUE ~ combined_plot_df$Description
)

combined_plot_df <- combined_plot_df %>%
  dplyr::group_by(Source, PathwayGroup) %>%
  dplyr::arrange(p.adjust, dplyr::desc(Count)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

combined_plot_df$PathwayGroup <- factor(
  combined_plot_df$PathwayGroup,
  levels = c("Ribosome / translation", "RNA processing", "Immune signaling")
)

combined_plot_df$Source <- factor(
  combined_plot_df$Source,
  levels = c("mRNA", "miRNA targets")
)

write.csv(
  combined_plot_df,
  file.path(table_dir, "pathway_convergence.csv"),
  row.names = FALSE
)

p_convergence <- ggplot(
  combined_plot_df,
  aes(x = Source, y = PathwayGroup, size = Count, color = neglog10FDR)
) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "grey80", high = "#d7301f") +
  scale_size(range = c(4.5, 11)) +
  scale_x_discrete(expand = expansion(add = 0.25)) +
  labs(
    title = "D",
    x = NULL,
    y = NULL,
    color = "-log10(FDR)",
    size = "Gene count"
  ) +
  theme_fig4(14) +
  theme(
    axis.text.y = element_text(face = "bold")
  )

ggsave(
  file.path(output_dir, "convergence.pdf"),
  p_convergence,
  width = 5,
  height = 4.5
)

# ============================================================
# 8. Paired tumour-adjacent program plot
# ============================================================

message("Running paired tumour-adjacent program plot...")

paired_prog_df <- rbind(
  data.frame(
    Patient = corr_df$Patient,
    Outcome = corr_df$Outcome,
    Cohort = corr_df$Cohort,
    Program = "Ribosome",
    Tumour = corr_df$tumor_ribo,
    Adjacent = corr_df$adj_ribo,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Patient = corr_df$Patient,
    Outcome = corr_df$Outcome,
    Cohort = corr_df$Cohort,
    Program = "Immune",
    Tumour = corr_df$tumor_imm,
    Adjacent = corr_df$adj_imm,
    stringsAsFactors = FALSE
  )
)

paired_prog_long <- paired_prog_df %>%
  tidyr::pivot_longer(
    cols = c("Tumour", "Adjacent"),
    names_to = "Tissue",
    values_to = "Score"
  )

paired_prog_long$Tissue <- factor(
  paired_prog_long$Tissue,
  levels = c("Tumour", "Adjacent")
)

paired_prog_long$Program <- factor(
  paired_prog_long$Program,
  levels = c("Ribosome", "Immune")
)

write.csv(
  paired_prog_long,
  file.path(output_dir, "paired_program_plot_data.csv"),
  row.names = FALSE
)

# ============================================================
# Paired statistical tests (Tumour vs Adjacent)
# ============================================================

paired_tests <- paired_prog_df %>%
  dplyr::group_by(Program) %>%
  dplyr::summarise(
    p_value = wilcox.test(Tumour, Adjacent, paired = TRUE)$p.value,
    y = max(c(Tumour, Adjacent), na.rm = TRUE) + 0.04 * diff(range(c(Tumour, Adjacent), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = dplyr::case_when(
      is.na(p_value) ~ "p = NA",
      p_value < 0.001 ~ "p < 0.001",
      TRUE ~ paste0("p = ", sprintf("%.3f", p_value))
    ),
    x = 1.5
  )

print(paired_tests)

write.csv(
  paired_tests,
  file.path(output_dir, "paired_program_plot_stats.csv"),
  row.names = FALSE
)

p_paired_program <- ggplot(
  paired_prog_long,
  aes(x = Tissue, y = Score, group = Patient, color = Outcome)
) +
  geom_line(alpha = 0.45, linewidth = 0.7) +
  geom_point(size = 2.2, alpha = 0.9) +
  facet_wrap(~ Program, scales = "free_y") +
  geom_text(
    data = paired_tests,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_color_manual(values = c("Remission" = "#377EB8", "Recidivism" = "#E41A1C")) +
  labs(
    title = "Paired tumour-adjacent program scores",
    x = NULL,
    y = "Program score"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  )

ggsave(
  file.path(output_dir, "paired_program_plot.pdf"),
  p_paired_program,
  width = 6.5,
  height = 4.5
)
