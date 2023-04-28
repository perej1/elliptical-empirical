# Argument parsing
library(optparse)

cli::cli_h1("Emprical example starting")

option_list <- list(
  make_option("--modify", type = "logical", default = TRUE,
              help = "Perform initial data modification or not"),
  make_option("--test", type = "logical", default = TRUE,
              help = "Perform statistical tests or not"),
  make_option("--analyze", type = "logical", default = TRUE,
              help = "Perform analysis or not")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (opt$modify) {
  cli::cli_h2("Data modification starting")
  system("Rscript modify.R")
  cli::cli_h2("Data modification done")
}

if (opt$test) {
  cli::cli_h2("Statistical testing starting")
  system("Rscript test-assumptions.R")
  cli::cli_h2("Statistical tests done")
}

if (opt$analyze) {
  cli::cli_h2("Analysis starting")
  system("Rscript analyze.R")
  cli::cli_h2("Analysis done")
}

cli::cli_h1("Emprical example done")
