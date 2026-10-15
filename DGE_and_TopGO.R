#SLime Gland vs Other Tissues DGE

#Import Packages
library("tximport")
library("readr")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("tidyverse")
library("dplyr")
library("tidyverse")
library("ggplot2")

#Specify data location
dir <- "/projects/health_sciences/bms/biochemistry/dearden_lab/Taylor/stringtie/salmon_quant"
samples <- read.csv(file.path(dir, "RNASeq_Sample_Info.csv"), header=TRUE)

exclude_samples <- c("A2", "C2", "O1", "V2") 
samples <- samples[!(samples$Batch %in% exclude_samples), ]
samples <- samples[!is.na(samples$Condition) & samples$Condition != "", ]

# Set rownames to Batch
rownames(samples) <- samples$Batch

# Create file paths to quant.sf files using Batch as folder names
files <- file.path(dir, samples$Batch, "quant.sf")
names(files) <- samples$Batch
samples$Batch

# Read the file 
tx2gene <- read.delim(
  "/projects/health_sciences/bms/biochemistry/dearden_lab/Taylor/stringtie/gene_trans_map.tsv",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Swap columns to be 1. transcript and 2. gene
tx2gene <- tx2gene[, c(2,1)]
colnames(tx2gene) <- c("transcript", "gene")


# Tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
samples$Condition <- factor(
  samples$Condition,
  levels = c(
    "Antennae",
    "Brain",
    "Carcass",
    "Ovary",
    "Slime glands",
    "Tail segment",
    "Ventral nerve cord"
  ),
  labels = c(
    "Antennae",
    "Brain",
    "Carcass",
    "Ovary",
    "Slime glands",
    "Tail segment",
    "Ventral nerve cord"
  )
)

txi$counts <- txi$counts[, rownames(samples)]
txi$abundance <- txi$abundance[, rownames(samples)]
txi$length <- txi$length[, rownames(samples)]

# Convert Condition column to factor for groupings
samples$Condition <- gsub(" ", "_", samples$Condition)

samples$Individual <- factor(samples$Individual)
samples$Condition <- factor(samples$Condition)

# Create a DDS using the txi object and sample info
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ Individual + Condition)
nrow(ddsTxi)

# Prefiltering 
# Smallest group size = 6 biological replicates
smallestGroupSize <- 6
keep <- rowSums(counts(ddsTxi) >= 10) >= smallestGroupSize
ddsTxi <- ddsTxi[keep,]

nrow(ddsTxi)
#cut from 39,000 transcripts down to ~30,000

#Run DESeq
cat("1. Running DESeq2 on unified model (all tissues distinct)...\n")
dds_unified <- DESeq(ddsTxi)

#Run pairwise contrasts
# Define all baseline tissues
other_tissues <- c("Antennae", "Brain", "Carcass", "Ovary", "Tail_segment", "Ventral_nerve_cord")

padj_cutoff <- 0.01
lfc_cutoff  <- 2

sig_gene_lists <- list()
all_res_list   <- list()

for (bg in other_tissues) {
  # Extract contrast: Slime vs X tissue
  res_pair <- results(dds_unified, 
                      contrast = c("Condition", "Slime_glands", bg),
                      alpha = padj_cutoff)
  
  # Store the full results table for this contrast
  all_res_list[[bg]] <- as.data.frame(res_pair)
  
  # Identify significant genes for this specific contrast
  sig_genes <- rownames(subset(res_pair, !is.na(padj) & padj < padj_cutoff & log2FoldChange >= lfc_cutoff))
  sig_gene_lists[[bg]] <- sig_genes
}
#Identify strinct intersection (slime DEGs that are exclusively DE in slime after comparison to everything else)
# gene must be upregulated against every single tissue to pass
candidate_gene_ids <- Reduce(intersect, sig_gene_lists)

cat("\n======================================================\n")
cat("Total genes tested:                ", nrow(dds_unified), "\n")
cat("Gold-Standard Slime-Specific Genes:", length(candidate_gene_ids), "\n")
cat("======================================================\n\n")


#create composit results DF
cat("4. Building composite res_df and fetching TPMs...\n")

res_df <- data.frame(row.names = rownames(dds_unified))
res_df$baseMean <- all_res_list[[1]]$baseMean

# Extract all LFCs and padjs into matrices
lfc_matrix  <- sapply(all_res_list, function(x) x$log2FoldChange)
padj_matrix <- sapply(all_res_list, function(x) x$padj)

# Average Fold Change across the 6 contrasts
res_df$log2FoldChange <- rowMeans(lfc_matrix, na.rm = TRUE)

# Maximum (Worst) padj (ensures conservative significance reporting
res_df$padj <- apply(padj_matrix, 1, max, na.rm = TRUE)
res_df$padj[is.infinite(res_df$padj)] <- NA

#  TPM Abundances directly from tximport
tpm_matrix <- txi$abundance[rownames(res_df), ]

# Calculate Mean TPM for Slime vs the average of all other samples
slime_cols <- rownames(samples[samples$Condition == "Slime_glands", ])
other_cols <- rownames(samples[samples$Condition != "Slime_glands", ])

res_df$mean_tpm_slime <- rowMeans(tpm_matrix[, slime_cols], na.rm = TRUE)
res_df$mean_tpm_other <- rowMeans(tpm_matrix[, other_cols], na.rm = TRUE)

#export candidates
slime_candidates <- res_df[candidate_gene_ids, ]

# Sort by highest Slime TPM (or Fold Change)
slime_candidates <- slime_candidates[order(slime_candidates$mean_tpm_slime, decreasing = TRUE), ]

# View the top hits
cat("Top 10 Slime Candidates:\n")
print(head(slime_candidates[, c("log2FoldChange", "padj", "mean_tpm_slime", "mean_tpm_other")], 10))

# Export for GO and Trinotate
write.csv(slime_candidates, "Slime_Gland_Enriched_Genes_Pairwise.csv", row.names = TRUE)



### VOLCANO PLOT ###
library(ggplot2)
library(ggrepel)

#prepare from results DF
plot_df <- as.data.frame(res_df)
plot_df$gene_id <- rownames(plot_df)

# Remove any NA rows
plot_df <- plot_df[!is.na(plot_df$padj) & !is.na(plot_df$log2FoldChange), ]

# Prevent -Inf on y-axis for p-values of 0
min_p <- min(plot_df$padj[plot_df$padj > 0], na.rm = TRUE)
plot_df$padj[plot_df$padj == 0] <- min_p * 1e-2
plot_df$y_val <- -log10(plot_df$padj)

# Assign significance categories
plot_df$status <- "Not Significant"

# Define Downregulated (Average LFC <= -2, Max padj < 0.01)
plot_df$status[plot_df$padj < 0.01 & plot_df$log2FoldChange <= -2] <- "Downregulated in Slime"

# Only label Upregulated if passed the rigorous intersection!
plot_df$status[plot_df$gene_id %in% rownames(slime_candidates)] <- "Strict Slime Marker"

# Set factor levels for the legend
plot_df$status <- factor(plot_df$status, levels = c("Strict Slime Marker", "Downregulated in Slime", "Not Significant"))


# Attach trinotate gene names
plot_df$label <- plot_df$gene_id

if (exists("trinotate")) {
  g_col <- grep("gene_id", colnames(trinotate), ignore.case = TRUE, value = TRUE)[1]
  n_col <- grep("Preferred_name", colnames(trinotate), ignore.case = TRUE, value = TRUE)[1]
  
  if (!is.na(g_col) && !is.na(n_col)) {
    match_idx <- match(plot_df$gene_id, trinotate[[g_col]])
    names_found <- trinotate[[n_col]][match_idx]
    has_valid_name <- !is.na(names_found) & names_found != "" & names_found != "."
    plot_df$label[has_valid_name] <- names_found[has_valid_name]
  }
}

#label the top 10 upregulated slime genes with trinotate name if available
up_genes <- plot_df[plot_df$status == "Strict Slime Marker", ]

# Order by statistical significance and fold-change
up_genes <- up_genes[order(up_genes$y_val, up_genes$log2FoldChange, decreasing = TRUE), ]
labels_to_plot <- head(up_genes, 10)


#plot volcano plot
p_volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = y_val)) +
  
  geom_point(aes(color = status), alpha = 0.45, size = 1.4) +
  
  scale_color_manual(values = c(
    "Strict Slime Marker"    = "firebrick2",  
    "Downregulated in Slime" = "steelblue3",
    "Not Significant"        = "grey75" #may be DEG for one comparison but has to be all to be red
  )) +
  
  # Threshold guide lines
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  
  geom_text_repel(data = labels_to_plot, 
                  aes(label = label),
                  size = 3.3, 
                  fontface = "bold", 
                  max.overlaps = 50,
                  box.padding = 0.4,
                  point.padding = 0.3) +
  
  theme_bw(base_size = 11) +
  labs(
    title    = "Volcano Plot: True Slime Gland Markers",
    subtitle = "Upregulated vs. ALL other tissues (Intersected Pairwise Contrasts)",
    x        = expression(Average~Log[2]~Fold~Change),
    y        = expression(Maximum~-Log[10]~Adjusted~italic(P)-value),
    color    = ""
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13)
  )

print(p_volcano)
ggsave("Slime_Gland_Volcano_Plot_Pairwise.png", plot = p_volcano, width = 8.5, height = 7, dpi = 300)



##### TOP GO ANALYSIS #######

gene2GO <- readMappings("/projects/health_sciences/bms/biochemistry/dearden_lab/Taylor/trinotate/gene2go.tsv", sep = "\t", IDsep = ",")

#build topGO universe
cat("1. Building universe and selection vectors...\n")

# All genes evaluated in DESeq2 (from res_df)
all_tested_genes <- rownames(res_df)

# Universe = Genes tested in DESeq2 that ALSO have GO annotations
all_universe_genes <- intersect(names(gene2GO), all_tested_genes)

# Target genes = Slime gland candidates from the pairwise intersection
target_candidates <- candidate_gene_ids

# Named binary factor: 1 for Target (Slime Candidate), 0 for Background
geneList <- factor(as.integer(all_universe_genes %in% target_candidates))
names(geneList) <- all_universe_genes

cat(paste(" -> Total Universe Size:", length(geneList), "annotated genes\n"))
cat(paste(" -> Annotated Slime Candidates:", sum(geneList == 1), "genes\n"))

# Function for topGO to identify interesting genes
selection_function <- function(x) { return(x == 1) }


# run topGO (weight01 + Fisher)
cat("2. Running topGO with weight01 algorithm...\n")

# Build topGOdata object
chosen_ontology <- "CC"   # Change this to "MF" or "CC" for the other plots

GOdata <- new(
  "topGOdata",
  ontology = chosen_ontology,
  allGenes = geneList,
  geneSel  = selection_function,
  annot    = annFUN.gene2GO,
  gene2GO  = gene2GO,
  nodeSize = 3 # Excludes extremely obscure terms supported by <3 genes
)

# Run the topology-aware weight01 test
results_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

# Extract results table (numChar = 1000 prevents label truncation)
top_results <- GenTable(
  GOdata,
  weight01_p = results_weight01,
  orderBy    = "weight01_p",
  topNodes   = 40,
  numChar    = 1000
)

# Format p-values, calculate Fold Enrichment and -log10(p-value)
top_results <- top_results %>%
  mutate(
    weight01_p      = as.numeric(gsub("<", "", weight01_p)),
    neg_log10_p     = round(-log10(weight01_p), 2),
    Fold_Enrichment = round(Significant / Expected, 2)
  ) %>%
  filter(weight01_p < 0.05) %>%
  arrange(desc(Fold_Enrichment))


#generate figures
cat("3. Generating figure...\n")

plot_data <- top_results %>% head(43)

p_topgo <- ggplot2::ggplot(plot_data, aes(x = Fold_Enrichment, y = reorder(Term, Fold_Enrichment))) +
  geom_point(aes(size = Significant, color = neg_log10_p)) +
  scale_color_viridis_c(option = "magma", name = "-log10(p-value)") +
  scale_size_continuous(name = "Gene Count", range = c(3, 8)) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_bw(base_size = 11) +
  labs(
    title    = paste("Top Enriched", chosen_ontology, "Terms in Slime Gland"),
    subtitle = "Strict Pairwise Intersection | topGO (weight01 graph algorithm)",
    x        = "Fold Enrichment (Observed / Expected)",
    y        = ""
  ) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 9, color = "grey30"),
    axis.text.y     = element_text(size = 9.5, face = "bold", color = "black"),
    legend.position = "right"
  )

print(p_topgo)
ggsave(paste0("Slime_Gland_topGO_", chosen_ontology, "_Pairwise.png"), plot = p_topgo, width = 9.5, height = 7, dpi = 300)

#identify genes driving terms
cat("4. Extracting driving genes for each term...\n")

# Pull all significant genes directly from GOdata
sig_genes_vec <- sigGenes(GOdata)

# Attach the driving gene IDs to every row in top_results
top_results$Driving_Genes <- sapply(top_results$GO.ID, function(id) {
  in_term <- genesInTerm(GOdata, id)[[1]]
  paste(intersect(in_term, sig_genes_vec), collapse = ", ")
})

# Reorder columns cleanly
top_results_clean <- top_results %>%
  dplyr::select(GO.ID, Term, Significant, Expected, Fold_Enrichment, weight01_p, neg_log10_p, Driving_Genes)

# View in spreadsheet
View(top_results_clean)

# Save to CSV
write.csv(top_results_clean, paste0("topGO_", chosen_ontology, "_Enriched_Terms_Summary_Pairwise.csv"), row.names = FALSE)



