
library(genoPlotR)
library(tidyverse)
library(Cairo)

#load data
p_raw <- read.table("/weka/health_sciences/bms/biochemistry/dearden_lab/galta815/blast/additional_slime_protein_tblastn.txt", header=FALSE)
colnames(p_raw) <- c("q", "subject", "sstart", "send", "e", "qcov", "len", "bit")

e_raw <- read.table("/weka/health_sciences/bms/biochemistry/dearden_lab/galta815/blast/epiperipatus_broadwayi_tblastn.txt", header=FALSE)
colnames(e_raw) <- c("q", "subject", "sstart", "send", "e", "qcov", "len", "bit")

#specify windows
p_chr <- "Scaffold_18__2_contigs__length_300523913"
e_chr <- "JAQFVV010012593.1"

p_win1 <- c(204100000, 204650000) # Sub-cluster A
p_win2 <- c(204650001, 205300000) # Sub-cluster B (Captures the 205.20Mb gene)
e_win  <- c(800000, 1150000)      # Ebro Ancestral Anchor

make_final_seg <- function(raw_df, target_chr, window) {
  sub <- raw_df[raw_df$subject == target_chr, ]
  sub$actual_start <- pmin(sub$sstart, sub$send)
  sub$actual_end   <- pmax(sub$sstart, sub$send)
  sub <- sub[sub$actual_start >= window[1] & sub$actual_end <= window[2], ]
  
  if(nrow(sub) == 0) return(NULL)
  
  # Merge fragments within 500bp to define distinct gene loci
  sub <- sub[order(sub$actual_start), ]
  loci <- list()
  curr_start <- sub$actual_start[1]; curr_end <- sub$actual_end[1]
  for(i in 2:nrow(sub)) {
    if(sub$actual_start[i] < (curr_end + 500)) { 
      curr_end <- max(curr_end, sub$actual_end[i])
    } else {
      loci[[length(loci)+1]] <- c(curr_start, curr_end)
      curr_start <- sub$actual_start[i]; curr_end <- sub$actual_end[i]
    }
  }
  loci[[length(loci)+1]] <- c(curr_start, curr_end)
  
  # Gold genes for a high-impact "Precious Locus" look
  df <- do.call(rbind, lapply(1:length(loci), function(i) {
    data.frame(name="", start=loci[[i]][1], end=loci[[i]][2],
               strand=1, col="darkgoldenrod", fill="gold", gene_type="arrows")
  }))
  return(dna_seg(df))
}

seg_pote1 <- make_final_seg(p_raw, p_chr, p_win1)
seg_ebro  <- make_final_seg(e_raw, e_chr, e_win)
seg_pote2 <- make_final_seg(p_raw, p_chr, p_win2)

#link up
make_fan_comp <- function(s1, s2) {
  links <- expand.grid(i1 = 1:nrow(s1), i2 = 1:nrow(s2))
  comparison(data.frame(
    start1 = s1$start[links$i1], end1 = s1$end[links$i1],
    start2 = s2$start[links$i2], end2 = s2$end[links$i2],
    col = rgb(0.8, 0.6, 0, alpha = 0.2) # Translucent gold ribbons
  ))
}

compA <- make_fan_comp(seg_pote1, seg_ebro)
compB <- make_fan_comp(seg_ebro, seg_pote2)

#plot
#CairoPDF("Figure_PanelC_Gold_Synteny_Final.pdf", width = 14, height = 12)

plot_gene_map(
  dna_segs = list("P. otepoti (12 Gb) - Sub-cluster A" = seg_pote1, 
                  "E. broadwayi (5.6 Gb) - Ancestral Anchor" = seg_ebro, 
                  "P. otepoti (12 Gb) - Sub-cluster B" = seg_pote2),
  comparisons = list(compA, compB),
  xlims = list(p_win1, e_win, p_win2),
  main = "Tandem Expansion and Intergenic Bloat of Slime Protein Loci",
  gene_type = "arrows",
  arrow_head_len = 4000, 
  dna_seg_line = "grey40",
  scale = TRUE,
  dna_seg_scale = TRUE,
  main_pos = "left"
)


#dev.off()

cat("\nSUCCESS: Gold synteny figure generated.")
cat("\nTotal genes in P. otepoti window:", nrow(seg_pote1) + nrow(seg_pote2), "\n")