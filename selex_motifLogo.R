library(tidyverse)
library(Biostrings)
library(msa)
library(ggseqlogo)
library(ips)
library(argparse)

get_motif_logo <- function(input, output, motif_file, align_mode, n) {
  
  if (align_mode %in% c("MAFFT", "Clustal")) {
    print(paste("Processing alignment using the", align_mode, "algorithm..."))
  } else {
    stop("Wrong alignment mode, please use either MAFFT or Clustal!")
  }
  
  input <- read_csv(input) %>%
    #filter(cycle == max(.$cycle)) %>%
    #Only selects the highest hits at the last cycle
    arrange(desc(value)) %>%
    dplyr::slice(1:n) %>%
    mutate(relative_value = round((value / sum(.$value) * 1000), 0))
  
  #Calculates relative fold change by summing up all of the fold change values, and then dividing each fold change to this value. Finally, multiple this by a 1000.
  #This will make weighting each motif much more efficient since they will be multipled by this value. 
  
  rows <- nrow(input)
  
  weighted_kmer <- c()
  
  for (i in seq(rows)) {
    kmer <- rep(input$kmer[i], input$relative_value[i])
    weighted_kmer <- c(weighted_kmer, kmer)
  }
  
  weighted_kmer_DNA <- DNAStringSet(weighted_kmer)
  
  if (align_mode == "Clustal") {
    
    print(paste("Running", align_mode, "alignment..."))
    
    dna_string_set <- msaClustalOmega(weighted_kmer_DNA)
    
  } else if (align_mode == "MAFFT") {
    
    print(paste("Running", align_mode, "alignment..."))
    
    weighted_kmer_DNAbin <- as.DNAbin(weighted_kmer_DNA)
    
    aligned_kmer <- as.character(mafft(weighted_kmer_DNAbin, method = "auto", maxiterate = 1000, op = 2,
                                       ep = 0.2, options = c("--adjustdirectionaccurately"), thread = -1, exec = "/usr/local/bin/mafft"))
    
    # Convert the character matrix to a character vector
    aligned_kmer_vector <- apply(aligned_kmer, 1, paste, collapse = "")
    
    # Convert the character vector to a DNAStringSet object
    dna_string_set <- DNAStringSet(aligned_kmer_vector)
    
  }
  
  consensus <- consensusString(dna_string_set, ambiguityMap = "N", threshold = 0.5)
  
  print(paste("Consensus motif for", protein, "with k-mer size", k, "is:", consensus))
  
  write_lines(consensus, motif_file)
  
  PSSM <- consensusMatrix(dna_string_set, baseOnly = T, as.prob = T)
  
  pwm <- PSSM[DNA_BASES, ]
  
  motif_plot <- ggseqlogo(pwm, method='custom')
  
  ggsave(output, motif_plot)
  
  #seqLogo::seqLogo(pwm, ic.scale = F)
  
  #pwm_sum <- colSums(pwm)
  
  #pwm_fixed <- t(apply(pwm, 1, function(x) x/pwm_sum))
  
  return(motif_plot)
}

# Create a new parser object
parser <- ArgumentParser()

# Add arguments
parser$add_argument('--input_file', required=TRUE, help='path to .csv file with information')
parser$add_argument('--output_folder', required = TRUE, help='Path to output folder')
parser$add_argument('--kmer_min', required = TRUE, help = 'Minimum kmer size')
parser$add_argument('--kmer_max', required = TRUE, help = 'Maximum kmer size')
parser$add_argument('--align_mode', default = "MAFFT", choices = c('MAFFT', 'Clustal'), help = 'Algorithm to align motifs. Either "MAFFT", or "Clustal"')

# Parse command line arguments
args <- parser$parse_args()


kmer_list <- seq(as.numeric(args$kmer_min), as.numeric(args$kmer_max))
#List of kmers to work with

output_folder <- args$output_folder
#Output folder for processed files

protein_list <- read.csv(args$input_file) %>%
  filter(protein != 'None') %>%
  .$protein %>%
  unique()

for (protein in protein_list) {
  for (k in kmer_list) {
    print(paste("Plotting motif logo for:", protein, "with k-mer size", k))
    get_motif_logo(input = paste0(output_folder,protein,"/",protein,"_kmer_",k,"_selex_py.csv"),
                   output = paste0(output_folder,protein,"/",protein,"_kmer_",k,"_selex_py_top10_motifLogo.pdf"),
                   motif_file = paste0(output_folder,protein,"/",protein,"_kmer_",k,"_selex_py_top10_motif.txt"),
                   align_mode = args$align_mode,
                   n = 10)
    
  }
}


