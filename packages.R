## library() calls go here
library(conflicted)
library(dotenv)
library(targets)
library(tarchetypes)
library(dplyr)
library(purrr)
#library(mirai)
#daemons(3)
library(quarto)
library(ggplot2)
theme_set(cowplot::theme_cowplot())

suppressMessages(library(SummarizedExperiment))

tar_option_set(memory = "transient")
options(tflow.report_dir = "docs")
paint::mask_print()
