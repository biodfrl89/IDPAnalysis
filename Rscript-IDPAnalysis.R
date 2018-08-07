#!/usr/bin/env Rscript
library(optparse)

option_list <-  list(
  make_option( opt_str = c("-f", "--filename"), type = "character", default = NULL, 
               help = "The fasta file with the sequence(s) to be analized.", metavar = "filename"),
  make_option( opt_str = c("-o", "--output"), type = "character", default = "output.txt", 
               help = "The filename of the output [default = %default]", metavar = "filename")
) 

opt_parser <-  OptionParser(option_list = option_list)
opt <-  parse_args(opt_parser)

if (is.null(opt$file)) {
  stop("At least one argument must be supplied (input file).\n\n Look the help section with the argument -h or --help \n\n", call.=FALSE)
}

s <- 0

library(seqinr)
fasta <- read.fasta(file = opt$filename, seqtype = "AA", as.string = FALSE, set.attributes = FALSE)
fasta <- read.fasta(file = "LEA4_fasta.fasta", seqtype = "AA", as.string = FALSE, set.attributes = FALSE)

AA <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
        "M", "N", "P", "Q", "R", "S", "T", "V", "Y", "W")

df <- data.frame(matrix(vector(), 
                        0, 42, 
                        dimnames=list(c(), c("GeneID", paste0("cont_",AA), "Long_prot", paste0("prop_",AA)))
                        ),
                 stringsAsFactors=F)
i = 1 
while(i <= length(fasta)) {
  prot_sample <- fasta[[i]]
  prot_sample_name <- names(fasta[i])
  length_prot <- length(fasta[[i]])
  for (j in AA) {
    assign( paste0("cont_",j), length( which(fasta[[i]] == j)))
    assign( paste0("prop_",j), length(which(fasta[[i]] == j)) * 100 / length(fasta[[i]]))
  }
  df[i,] <- c(prot_sample_name,
              cont_A, cont_C, cont_D, cont_E, cont_F, cont_G, cont_H, cont_I, cont_K, cont_L,
              cont_M, cont_N, cont_P, cont_Q, cont_R, cont_S, cont_T, cont_V, cont_Y, cont_W,
              length_prot,
              round(prop_A, 3), round(prop_C,3), round(prop_D,3), round(prop_E,3), round(prop_F,3),
              round(prop_G,3), round(prop_H,3), round(prop_I,3), round(prop_K,3), round(prop_L,3),
              round(prop_M,3), round(prop_N,3), round(prop_P,3), round(prop_Q,3), round(prop_R,3),
              round(prop_S,3), round(prop_T,3), round(prop_V,3), round(prop_Y,3), round(prop_W,3)
  )
  i = i + 1
}
rm(list=eval(ls()[grep(pattern = "prop_", ls())]))
rm(list=eval(ls()[grep(pattern = "cont_", ls())]))
rm(list=c("prot_sample", "prot_sample_name", "i","j", "length_prot"))

if (s == 0) stop("Test Done", call. = TRUE)
