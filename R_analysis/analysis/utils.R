if (!require("pacman")) install.packages("pacman")
pkgs = c(
    'multidplyr',
    'dplyr',
    'arrow',
    'haven',
    'readr',
    'fixest',
    'glue',
    'ggplot2',
    'foreach',
    'doParallel',
    'stringr',
    'ggpubr',
    'purrr',
    'memoise'
)
pacman::p_load(pkgs, character.only = TRUE)
here::i_am('R_analysis/analysis/utils.R')
CODE = here::here()



