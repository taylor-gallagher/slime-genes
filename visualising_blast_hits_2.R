#Use custom chromosome names
options(ucscChromosomeNames = FALSE)
library(dplyr)
library(txdbmaker)
library(Gviz)

#Load in BRAKER gff3
txdb <- txdbmaker::makeTxDbFromGFF("/weka/users/guhjo98p/velvetworm/source_invest/braker_busco_rerun_2026-05-21_existingbam_v3/braker_clean_etp/braker.gff3", format="gff3")

#Load in BLAST hits
blast_data <- read.table("/projects/health_sciences/bms/biochemistry/dearden_lab/Taylor/blast/additional_slime_protein_tblastn.txt", header=FALSE)
colnames(blast_data) <- c("query", "subject", "sstart", "send", "evalue", "qcovs", "length", "bitscore")

#BLAST to GRanges object
blast_gr <- GRanges(
  seqnames = blast_data$subject,
  ranges = IRanges(start = pmin(blast_data$sstart, blast_data$send), 
                   end = pmax(blast_data$sstart, blast_data$send)),
  strand = ifelse(blast_data$sstart < blast_data$send, "+", "-"),
  id = blast_data$query,
  group = blast_data$query
)

# Read in BigWig RNA-seq files
bw_dir <- "/projects/health_sciences/bms/biochemistry/dearden_lab/Taylor/Trimgalore_Remapped_VW/bigwig" 
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

#Plotting parameters
target_chr <- "Scaffold_2__1_contigs__length_717158983" 
target_start <- 275420000    
target_end   <- 275495000

#Create tracks
#RNA-seq Tracks
rna_tracks <- lapply(names(tissue_list), function(t_name) {
  files <- tissue_list[[t_name]]
  track_label <- paste0(t_name, "\nBPM")
  
  OverlayTrack(trackList = lapply(files, function(f) {
    DataTrack(range = f, type = "histogram", name = track_label, 
              fill = tissue_colors[t_name], col = tissue_colors[t_name],
              lwd = 0, alpha = 0.4, transformation = function(x) log2(x + 1),
              ylim = c(0, 17), baseline = 0, showAxis = TRUE, col.axis = "black",      
              cex.axis = 0.4, fontcolor.title = "black", col.title = "black",       
              cex.title = 0.5, background.title = "#E0E0E0", rotation.title = 0)
  }), name = track_label) 
})

# StringTie Track
stringtie_txdb <- txdbmaker::makeTxDbFromGFF(
  "/projects/health_sciences/bms/biochemistry/dearden_lab/Taylor/stringtie/P_otepoti_master_curated.gtf", 
  format = "gtf"
)

# THE FIX: Added showId = FALSE to kill the text rendering engine completely
trinity_tr <- GeneRegionTrack(stringtie_txdb,
                              chromosome = target_chr,
                              start = target_start, # Passed start/end here saves huge amounts of memory
                              end = target_end,
                              name = "StringTie",
                              fill = "#2E7D32", 
                              col = "#1B5E20",
                              shape = "arrow",
                              collapseTranscripts = "meta",
                              showId = FALSE,      # <-- KILL SWITCH 1
                              background.title = "#E0E0E0",
                              rotation.title = 0)

#BLAST track
blast_tr <- AnnotationTrack(blast_gr, 
                            name = "Slime\nBLAST", 
                            fill = "darkgoldenrod", 
                            col = "darkgoldenrod",
                            stacking = "dense",
                            showId = FALSE,        # <-- KILL SWITCH 2
                            background.title = "#E0E0E0",
                            rotation.title = 0)

#BRAKER track
# THE FIX: Removed transcriptAnnotation = "transcript" and added showId = FALSE
braker_tr <- GeneRegionTrack(txdb, 
                             chromosome = target_chr, 
                             start = target_start,
                             end = target_end,
                             name = "BRAKER3",
                             fill = "darkblue", 
                             col = "darkblue",
                             collapseTranscripts = "meta",
                             stacking = "dense",
                             shape = "arrow", 
                             showId = FALSE,       # <-- KILL SWITCH 3
                             background.title = "#E0E0E0",
                             rotation.title = 0)

# Plotting
track_list <- c(
  list(GenomeAxisTrack(col="black", fontcolor="black", cex=0.7)), 
  rna_tracks, 
  list(braker_tr, trinity_tr, blast_tr)
)

track_sizes <- c(1, rep(1.5, 7), 1.5, 4, 1)

CairoPDF("Chromosome_2.pdf", width = 14, height = 16)
p <- plotTracks(track_list,
                chromosome = target_chr, 
                from = target_start, 
                to = target_end,
                sizes = track_sizes, 
                title.width = 2.5,       
                main = "Chromosome 2 - Candidate Slime Genes",
                cex.main = 1.2,
                fontcolor.title = "black", 
                col.title = "black",
                innerMargin = 1,          
                margin = 20,
                extend.left = 0.01,
                collapseTranscripts = "meta",
                stackHeight = 0.5,
                showId = FALSE)           
dev.off()