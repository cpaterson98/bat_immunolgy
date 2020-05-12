library(biomaRt)


    
GO_gene_info <-  read_csv("Data/gene_biological_info_gff_and_RBH_PVAP_MLUC_PALE_RAEG.csv") %>% 
  as.data.frame() 

GO_info <- GO_gene_info %>% 
  filter(EHELV_gene_ID %in% final_count_list$EHELV_gene_ID) 
GO_info <- GO_info[,c(1,38)]
 View(final_count_list)
  View(GO_info)
immune_genes <- final_count_list[c(1,14)]

immune_genes <- inner_join(immune_genes,GO_info,c("EHELV_gene_ID")) 

  
immune_genes_1 <- immune_genes    
  
immune_genes <- immune_genes %>% 
  filter(grepl("TLR", merged_gene_name))

immune_genes_to_add <- immune_genes_1 %>% 
  filter(grepl("DDX58", merged_gene_name))

immune_genes <- rbind(immune_genes,immune_genes_to_add)


immune_genes_to_add <- immune_genes_1 %>% 
  filter(grepl("NLR", merged_gene_name))

immune_genes <- rbind(immune_genes,immune_genes_to_add)

immune_genes_to_add <- immune_genes_1 %>% 
  filter(grepl("NOD", merged_gene_name))

immune_genes <- rbind(immune_genes,immune_genes_to_add)

immune_genes_to_add <- immune_genes_1 %>% 
  filter(grepl("IL", merged_gene_name))

immune_genes <- rbind(immune_genes,immune_genes_to_add)

immune_genes_to_add <- immune_genes_1 %>% 
  filter(grepl("IFI", merged_gene_name))

immune_genes <- rbind(immune_genes,immune_genes_to_add)


immune_genes_to_add <- immune_genes_1 %>% 
  filter(grepl("RSAD2", merged_gene_name))

immune_genes <- rbind(immune_genes,immune_genes_to_add)



#Take Care to remove accidental genes included such as COIL CARMIL NRKBIL etc.
View(immune_genes)

immune_genes <- immune_genes[-c(12,14,15,19:21,23,24,33,36,37,38,41,45:47),]

View(immune_genes)





