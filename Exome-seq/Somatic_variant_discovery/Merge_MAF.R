#!/bin/bash


args <- commandArgs(trailingOnly = TRUE)
# 1 args - final vcf path


library(maftools)
library(tidyverse)
maf_path_funco <- list.files(path = args[1], pattern = "*anno.maf", full.names = T)
maf_path_vep <- list.files(path = args[1], pattern = "*.vep.maf", full.names = T)

maf_list_vep <- lapply(X = maf_path_vep, FUN = function(path){
  read.maf(maf = path) %>% return()
})
maf_list_funco <- lapply(X = maf_path_funco, FUN = function(path){
  read.maf(maf = path) %>% return()
})

merged_maf_vep <- merge_mafs(maf_list_vep)
merged_maf_funco <- merge_mafs(maf_list_funco)

write.mafSummary(maf = merged_maf_vep, basename = paste0(args[1], '/merged_vef'))
write.mafSummary(maf = merged_maf_funco, basename = paste0(args[1], '/merged_funco'))