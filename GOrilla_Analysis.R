library(tidyverse)
gene_info_2 <- read_csv("Data/gene_biological_info_gff_and_RBH_PVAP_MLUC_PALE_RAEG.csv")
View(gene_info_2)


getwd()
gene_list <- gene_info_2[,c(1,38)] %>% 
  filter( EHELV_gene_ID %in% final_count_list$EHELV_gene_ID)
background_gene_list <- gene_list  %>% 
  filter(!is.na(gene_list$merged_gene_name) )
background_gene_list <- background_gene_list[,c(2)] 
  row.names(background_gene_list) = NULL

write.csv( background_gene_list, file = "background_gene_list_cluster_4w.csv", row.names = FALSE)
View(background_gene_list)
#28 comes from net$colors
for (i in 0:31){
module_list <- final_count_list %>% 
  filter(classification_WGCNA == i )
  module_list <- module_list[,c(1)] %>% 
    as.data.frame()
  colnames(module_list)[1] <- "EHELV_gene_ID"
   
  assign(paste("module_list",i, sep = "_"), module_list)

target_gene_list <-  gene_info_2[,c(1,38)] %>% 
  filter( EHELV_gene_ID %in% module_list$EHELV_gene_ID)
target_gene_list <- target_gene_list  %>% 
  filter(!is.na(target_gene_list$merged_gene_name) )
target_gene_list <- target_gene_list[,c(2)]

assign(paste("target_gene_list",i, sep = "_"),target_gene_list)
write.csv( target_gene_list, file = paste(paste("target_gene_list_cluster_4w",i, sep = "_"), "csv", sep = "."), row.names = FALSE)
}

