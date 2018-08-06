#!/usr/bin/env Rscript
library(optparse)

option_list <-  list(
  make_option( opt_str = c("-f", "--file"), type = "character", default = NULL, 
               help = "The fasta file with the sequence(s) to be analized.", metavar = "filename"),
  make_option( opt_str = c("-o", "--out"), type = "character", default = "out.txt", 
               help = "The filename of the output [default = %default]", metavar = "filename"),
  make_option( opt_str = c("-db", "--database"), type = "character", default = "place", 
               help = "A string defining the database that is going to be used to make the analisys [default = %default]. It can be modified to 'atcis'", 
               metavar = "String")
) 

opt_parser <-  OptionParser(option_list = option_list)
opt <-  parse_args(opt_parser)

if (is.null(opt$file)) {
  stop("At least one argument must be supplied (input file).\n\n Look the help section with the argument -h or --help \n\n", call.=FALSE)
}
