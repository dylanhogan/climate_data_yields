suppressMessages(here::i_am('R_analysis/init.R'))

# Load required packages
if (!require("pacman")) install.packages("pacman")
pkgs = c(
    'dplyr',
    'data.table',
    'tidyr',
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
    'memoise',
    'marginaleffects',
    'xtable',
    'ggh4x',
    'sf',
    'doRNG',
    'cowplot'
)
pacman::p_load(pkgs, character.only = TRUE)

# Set data path
DATA = read_file(here::here('data_path.txt'))

# Set cache size
cm = cachem::cache_mem(max_size = 5000 * 1024^2)

local_load = function(rel_path) {
    source(glue("{here::here()}/R_analysis/{rel_path}"))
}
    
init = TRUE