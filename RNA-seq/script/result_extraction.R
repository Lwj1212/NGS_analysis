library(tidyverse)
library(data.table)

setwd("~/RNA_SEQ/work/kallisto/kallisto_result/PUBLISH/EBI/0513/")
kallistoResult <- function(root_path, transcript_id){
  setwd(root_path) # kallisto result folder
  folder_dir <- list.files(full.names = T)
  
  temp_result <- lapply(X = folder_dir, FUN = function(dir){
    name <- dir %>% str_extract(string = ., pattern = "(?<=./)[a-zA-Z0-9]+")
    abundance <- fread(file = paste0(dir, "/abundance.tsv"))
    
    temp_id_DF <- tibble()
    
    for(id in transcript_id){
      temp_id_DF <- temp_id_DF %>% 
        bind_rows(., abundance %>% filter(str_detect(string = target_id, pattern = id)))
    }
    
    base_DF <- temp_id_DF %>% select(target_id, est_counts, tpm) %>% 
      mutate(tpm_log = log(x = (tpm+1), base = 2))
    
    count_DF <- base_DF %>% select(target_id, value = est_counts) %>% 
      bind_rows(., base_DF %>% select(target_id, value = tpm)) %>% 
      bind_rows(., base_DF %>% select(target_id, value = tpm_log))
    
    names(count_DF) <- c("target_id", name)
    return(count_DF)
    
  })
  
  return(temp_result %>% bind_cols())
}
cufflinkResult <- function(root_path, transcript_id){
  setwd(root_path) # kallisto result folder
  folder_dir <- list.files(full.names = T)
  
  temp_result <- lapply(X = folder_dir, FUN = function(dir){
    name <- dir %>% str_extract(string = ., pattern = "(?<=./)[a-zA-Z0-9]+")
    cuff <- fread(file = paste0(dir, "/isoforms.fpkm_tracking"))
    
    temp_id_DF <- tibble()
    cuff_fpkm_sum <- sum(cuff$FPKM)
    
    for(id in transcript_id){
      temp_id_DF <- temp_id_DF %>% 
        bind_rows(., cuff %>% filter(str_detect(string = tracking_id, pattern = id)))
    }
    
    base_DF <- temp_id_DF %>% select(tracking_id, gene_id, gene_short_name, FPKM) %>% 
      mutate(tpm = (FPKM / cuff_fpkm_sum) * 1e6) %>% 
      mutate(tpm_log = log(x = (tpm+1), base = 2))
    
    count_DF <- base_DF %>% select(tracking_id, value = FPKM) %>% 
      bind_rows(., base_DF %>% select(tracking_id, value = tpm)) %>% 
      bind_rows(., base_DF %>% select(tracking_id, value = tpm_log))
    
    colnames(count_DF) <- c("target_id", name)
    
    return(count_DF)
    
  })
  return(temp_result %>% bind_cols())
}



kallistoResult(root_path = "~/RNA_SEQ/work/kallisto/kallisto_result/PUBLISH/EBI/0513/",
               transcript_id = c("ENST00000296474.8", "RON_del5611", "RON_del56","ENST00000344206")) %>% 
 write_delim(delim = "\t", file = "../temp.txt")

cufflinkResult(root_path = "~/gcloud/GSE_cuff/cufflinks_result/",
               transcript_id = c("ENST00000296474.8", "ENST00000344206.8", "ENST00000256078", "ENST00000311936","ENST00000523873", "ENST00000523950", "ENST00000518824",
                                 "ENST00000307677", "ENST00000376055")) %>% write_delim(delim = "\t", file = "temp.txt")
