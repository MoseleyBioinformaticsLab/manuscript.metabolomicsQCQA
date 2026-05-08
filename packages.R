## library() calls go here
library(conflicted)
library(dotenv)
library(targets)
library(tarchetypes)
library(dplyr)
library(purrr)
#library(mirai)
#daemons(3)

tar_option_set(memory = "transient")
options(tflow.report_dir = "docs")
paint::mask_print()
library(quarto)
