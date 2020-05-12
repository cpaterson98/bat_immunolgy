

m_o_i_1w <- read_csv("module_of_interest_cluster_1w.csv")
m_o_i_4w <- read_csv("module_of_interest_cluster_4w.csv")
m_o_i_5w_out <- read_csv("module_of_interest_cluster_5w_out.csv")


View(m_o_i_1w)
View(m_o_i_4w)
View(m_o_i_5w_out)


m_o_i_1w <- m_o_i_1w[,-c(1)]
m_o_i_4w <- m_o_i_4w[,-c(1)]
m_o_i_5w_out <- m_o_i_5w_out[,-c(1)]


GO_list <- rbind(m_o_i_1w,m_o_i_4w)
GO_list <- rbind(GO_list,m_o_i_5w_out)
View(GO_list)
GO_list <- GO_list %>% 
  distinct(EHELV_gene_ID, .keep_all = TRUE)

GO_all <- rbind(m_o_i_1w,m_o_i_4w,m_o_i_5w_out)
View(GO_all)

  GO_list <- GO_list %>% 
    mutate("1_w" = EHELV_gene_ID %in% m_o_i_1w$EHELV_gene_ID) %>% 
    mutate("4_w" = EHELV_gene_ID %in% m_o_i_4w$EHELV_gene_ID) %>% 
  mutate("5_w" = EHELV_gene_ID %in% m_o_i_5w_out$EHELV_gene_ID)

  
  unannotated_genes <- GO_list %>% 
    filter(is.na(merged_gene_name))
View(unannotated_genes)  
  
annotated_genes <- GO_list %>% 
  filter(!is.na(merged_gene_name))
View(annotated_genes)  

  
  GO_list_totals <- GO_list %>% 
    mutate("sum" = rowSums(GO_list[,as.numeric(c(3,4,5))]))
View(GO_list_totals)  


multiple_modules <- GO_list_totals %>% 
  filter( sum == 3)
View(multiple_modules)
multiple_modules <- multiple_modules[,c(2)]

unannotated_multiple_modules <- multiple_modules[,c(1:5)] %>% 
  filter(is.na(merged_gene_name))

unannotated_multiple_modules <- unannotated_multiple_modules[,c(1,3,4,5)]
View(unannotated_multiple_modules)

annotated_multiple <- multiple_modules[,c(1:5)] %>% 
  filter(!is.na(merged_gene_name))

annotated_multiple <- annotated_multiple[,c(1,2)] %>% 
  column_to_rownames("EHELV_gene_ID")
View(annotated_multiple)


gene_counts_immune <- 

