BiocManager::install(c("Gviz", "rtracklayer", "GenomicFeatures"))

library(Gviz)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(txdbmaker)
library(dplyr)

#Use custom chromosome names
options(ucscChromosomeNames = FALSE)

#Load in BRAKER gff3
txdb <- txdbmaker::makeTxDbFromGFF("/weka/users/guhjo98p/velvetworm/source_invest/braker_busco_rerun_2026-05-21_existingbam_v3/braker_clean_etp/braker.gff3", format="gff3")

#Load in BLAST hits
blast_data <- read.table("/weka/health_sciences/bms/biochemistry/dearden_lab/galta815/blast/slime_protein_genome_blast_hits.txt", header=FALSE)
colnames(blast_data) <- c("query", "subject", "sstart", "send", "evalue", "qcovs", "length", "bitscore")

#BLAST to GRanges object
blast_gr <- GRanges(
  seqnames = blast_data$subject,
  ranges = IRanges(start = pmin(blast_data$sstart, blast_data$send), 
                   end = pmax(blast_data$sstart, blast_data$send)),
  strand = ifelse(blast_data$sstart < blast_data$send, "+", "-"),
  query_id = blast_data$query
)

# Read in BigWig RNA-seq files
bw_dir <- "/weka/health_sciences/bms/biochemistry/dearden_lab/galta815/rna-seq/star/trimgalore_mapped_reads_post_removal" 
bw_files <- list.files(bw_dir, pattern = "\\.bw$", full.names = TRUE)

tissue_map <- c("C" = "Carcass", "S" = "Slime glands", "V" = "VNC", 
                "A" = "Antennae", "B" = "Brain", "O" = "Ovaries", "T" = "Tail segment")

tissue_colors <- c("Carcass" = "brown", "Slime glands" = "blue", "VNC" = "magenta",
                   "Antennae" = "red", "Brain" = "orange", "Ovaries" = "green", "Tail segment" = "purple")

#Make RNA-seq metadata
metadata <- data.frame(file_path = bw_files, file_name = basename(bw_files)) %>%
  mutate(code = substr(file_name, 1, 1),
         tissue = tissue_map[code])

#Create list of files grouped by tissue
tissue_list <- split(metadata$file_path, metadata$tissue)

#Plotting parameters (change these for Scaff2 vs Scaff 18)
target_chr <- "Scaffold_18__2_contigs__length_300523913" 
target_start <- 102039000
target_end   <- 102050000
max_coverage <- 400000 

#Create tracks
#RNA-seq Tracks
merged_tracks <- lapply(names(tissue_list), function(t_name) {
  files <- tissue_list[[t_name]]
  
  individual_replicate_tracks <- lapply(files, function(f) {
    DataTrack(range = f, 
              type = "histogram", 
              name = t_name,
              fill.histogram = tissue_colors[t_name], 
              col.histogram = tissue_colors[t_name],
              ylim = c(0, max_coverage), 
              
              #Style
              background.title = "grey90",  
              fontcolor.title = "black",    
              col.axis = "black",           
              fontcolor.axis = "black",     
              cex.axis = 0.4,               
              cex.title = 0.5,              
              showAxis = TRUE           
    )
  })
  
  OverlayTrack(trackList = individual_replicate_tracks, 
               name = t_name, 
               ylim = c(0, max_coverage))
})

#BLAST track
blast_track <- AnnotationTrack(blast_gr, 
                               name = "Slime BLAST", 
                               fill = "gold", 
                               stacking = "squish",
                               #Style
                               background.title = "grey90", 
                               fontcolor.title = "black",   
                               fontcolor.group = "black",
                               cex.group = 0.5)

#BRAKER track
braker_track <- GeneRegionTrack(txdb, 
                                chromosome = target_chr, 
                                name = "BRAKER3",
                                fill = "darkblue", 
                                col = "black", 
                                transcriptAnnotation = "transcript",
                                #Style
                                background.title = "grey90", 
                                fontcolor.title = "black",   
                                fontcolor.group = "black",
                                cex.group = 0.5)

# Plotting
plotTracks(c(list(GenomeAxisTrack(col="black", fontcolor="black", cex=0.7)), 
             merged_tracks, 
             list(braker_track, blast_track)),
           chromosome = target_chr, 
           from = target_start, 
           to = target_end,
           sizes = c(1, rep(2, 7), 1, 2), 
           
           #Style
           title.width = 1.2,        
           margin = 40,              
           innerMargin = 4,          
           main = paste0("Slime Gene BLAST Hits - Chromosome 18"),
           cex.main = 1.2,
           fontcolor.main = "black")
