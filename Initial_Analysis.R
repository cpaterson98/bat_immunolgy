# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
workingDir = "/Users/charlottepaterson/Project/Sys_bio/Project/Analysis/Tester"; 
setwd(workingDir); 
#Firstly, I downloaded the relevant libraries and data sets. 

library(tidyverse)
library(ggfortify)
library(DESeq2)
library(ggrepel)
library(knitr)
library(WGCNA)



## Read in Data ##
sample_info <- read_csv("Data/colData_all_samples_included_27.02.2020.csv")
counts <- read_csv("Data/counttable_all_samples_included_27.02.2020.csv")
gene_info <- read_csv("Data/gene_biological_info_gff_and_RBH_PVAP.csv")

#I then proceded to change the formats of the tables to make them more suitable for the analyses and filtered out lowly expressed genes.


### Manipulate format of the datasets 
counts <- counts %>% 
  column_to_rownames( "X1") %>%  # this has genes as rownames
  rename_all(str_remove, "_ReadsPerGene.out")



sample_info <- sample_info %>% 
  column_to_rownames("X1") %>% 
  t() %>% 
  as.data.frame() %>% 
  rename_all(str_remove, "_ReadsPerGene.out") %>% 
  t()
sample_info <- sample_info %>% 
  as.data.frame()

row.names.remove <- c("EHELV_G00013549","EHELV_G00013550", "EHELV_G00013551", "EHELV_G00013552", "EHELV_G00015016","EHELV_G00027280", "EHELV_G00031256", "EHELV_G00032092", "EHELV_G00032093", "EHELV_G00034890", "EHELV_G00035000")

counts <-  counts[!(row.names(counts) %in% row.names.remove), ]





### Filter the counts for lowly represented genes 
keep <-  rowSums(counts) >600
counts <- counts[keep,]
librarySizes <- colSums(counts)


librarySizes <- librarySizes %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") 

colnames(librarySizes)[2] <- "lib_size"
librarySizes <- librarySizes %>% 
  mutate( "Individual" = sample_info$individual )



#This barplot shows us the gene expression levels for each sample, which I used as a quality check on the reads of each sample.

ggplot(data=librarySizes, aes(x=Sample , y=lib_size, fill = Individual)) +
  geom_bar(stat="identity") + labs( title = "Library Size of Each Sample",
                                    x = "Samples",
                                    y = "Library Size") + theme(axis.text.x=element_blank())



### Outlier

sample_outliers <- librarySizes$lib_size > 1.0e07

counts_outliers <- counts %>% 
  rownames_to_column("Gene") %>% 
  mutate("average" = rowSums(counts)/130) %>% 
  select(Gene,ebola_A0710_17_5900117, ebola_A0711_09_5900117,average)

counts_outliers <- counts_outliers %>% 
  mutate("ratio_710" = ebola_A0710_17_5900117/average ) %>% 
  mutate("ratio_711" = ebola_A0711_09_5900117/average ) %>% 
  select(Gene, ratio_710, ratio_711)



counts_outlier_genes <- counts_outliers %>% 
  filter(ratio_710 > 10 | ratio_711 >10)


#need gene_info_2 for this 


outlier_gene_list <- gene_info_2[,c(1,38)] %>% 
  filter( EHELV_gene_ID %in% counts_outlier_genes$Gene) %>% 
  filter(!is.na(merged_gene_name) )

outlier_gene_list <- outlier_gene_list[,c(2)]
row.names(outlier_gene_list) = NULL
write.csv( outlier_gene_list, file = "outlier_gene_list.csv", row.names = FALSE)










#Ebola Group Library Size

librarySizes_ebola <- librarySizes %>% 
  as.data.frame() %>% 
  mutate("group" = sample_info$treatment) %>% 
  filter(group == "ebola")
  names(librarySizes_ebola)[2] <- "library_size"

  ggplot(data=librarySizes_ebola, aes(x=Sample , y=library_size, fill = Individual)) +
    geom_bar(stat="identity") + labs( title = "Library Size of Each Sample in the Ebola Group",
                                      x = "Samples",
                                      y = "Library Size") + theme(axis.text.x=element_blank())
  
#As we can see, although there is some variation, I didnt see any reason that any of the samples should be ignored at this stage.

#I then proceded to make a master table, displaying together the sample information and the counts of each gene. I would use this to filter my data later on.

  librarySizes_sample <- librarySizes %>% 
    as.data.frame() %>% 
    mutate("group" = sample_info$treatment) %>%
    mutate("day" = sample_info$day) %>% 
    filter(group == "ebola") %>% 
    filter(day == "A0711" | day == "A0807")
  
  librarySizes_sample <- librarySizes_sample[-c(1,14),]
  names(librarySizes_sample)[2] <- "library_size"

  attach(librarySizes_sample)
  librarySizes <- librarySizes_sample[order(Individual),] 
  detach(librarySizes_sample)

  
  ggplot(data=librarySizes_sample, aes(x=Sample , y=library_size, fill = Individual)) +
    geom_bar(stat="identity") + labs( title = "Library Size of Each Sample in the Clustering Group",
                                      x = "Samples",
                                      y = "Library Size") + theme(axis.text.x=element_blank())
  
  
  
### Create Master Table 
# counts_t has samples as the rows

counts_t <- counts %>% 
  t()


big_table <- merge(sample_info, counts_t)


#To start the Analysis on the genes, I created lists filtering for certain characteristics, including all non-control bats post vaccination, all bats in the ebola group and samples post ebola vaccination.



#Group 1
#All bats in the sample
all_bats_sample_info <- sample_info %>% 
 rownames_to_column("Long_Sample")

all_bat_counts <- counts_t%>% 
  as.data.frame() %>% 
  rownames_to_column("Long_Sample") %>% 
  filter(Long_Sample %in% all_bats_sample_info$Long_Sample) %>% 
  column_to_rownames("Long_Sample") %>% 
  t()


#Group 2
#List of the Ebola Group 

ebola_group_bats <- sample_info %>% 
  as.data.frame()
ebola_group_bats <- ebola_group_bats[ebola_group_bats$treatment == "ebola",]
ebola_group_bats <- ebola_group_bats %>% 
  rownames_to_column("Long_Sample")




ebola_group_counts <- counts_t%>% 
  as.data.frame() %>% 
  rownames_to_column("Long_Sample") %>% 
  filter(Long_Sample %in% ebola_group_bats$Long_Sample) %>% 
  column_to_rownames("Long_Sample") %>% 
  t()





#Sample Bat Group

Before_and_After_infected_ebola_bats <- sample_info %>% 
  as.data.frame()
Before_and_After_infected_ebola_bats <- Before_and_After_infected_ebola_bats[Before_and_After_infected_ebola_bats$treatment == "ebola",]
Before_and_After_infected_ebola_bats <- Before_and_After_infected_ebola_bats[Before_and_After_infected_ebola_bats$day == "A0711"| Before_and_After_infected_ebola_bats$day == "A0814"  ,]
Before_and_After_infected_ebola_bats <- Before_and_After_infected_ebola_bats[-c(7,13),]
Before_and_After_infected_ebola_bats <- Before_and_After_infected_ebola_bats %>% 
  rownames_to_column("Long_Sample")

Before_and_After_infected_ebola_counts <- counts_t %>% 
  as.data.frame() %>% 
  rownames_to_column("Long_Sample") %>% 
  filter(Long_Sample %in% Before_and_After_infected_ebola_bats$Long_Sample) %>% 
  column_to_rownames("Long_Sample") %>% 
  t()  

View(Before_and_After_infected_ebola_bats)

### Manipulate and Normalise Data to get rid of mean dependant variance ###
#enter in place of counts the count data for
#whichever set of data you wish to run PCA on


### All Bats Boxplot

counter <- counts

logcounts <- log2(counter+1)

logcounts_long <- logcounts %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  mutate("Individual" = sample_info$individual )

long_format <- gather(logcounts_long, Gene, Log_Count, 2:13267, factor_key=TRUE)

ggplot(long_format, aes(x=Sample, y=Log_Count, fill = Individual)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=1,
               outlier.size=2) + ggtitle( 'Boxplot showing variation in count reads in each Sample') + theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_blank())


### Clusterin group Boxplot


counter <- Before_and_After_infected_ebola_counts
logcounts <- log2(counter +1)




logcounts_long <- logcounts %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  mutate(Individual = Before_and_After_infected_ebola_bats$individual ) %>% 
  mutate(day = Before_and_After_infected_ebola_bats$day )

long_format <- gather(logcounts_long, Gene, Log_Count, 2:13267, factor_key=TRUE)


ggplot(long_format, aes(x=Sample, y=Log_Count, fill = Individual, colour = day)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=1,
               outlier.size=2) + ggtitle( 'Boxplot showing variation in count reads in each sample of our Clustering Group') + theme(axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5)) 

# Ebola group Boxplot

counter <- ebola_group_counts
logcounts <- log2(counter +1)




logcounts_long <- logcounts %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  mutate("Individual" = ebola_group_bats$individual )

long_format <- gather(logcounts_long, Gene, Log_Count, 2:13267, factor_key=TRUE)


ggplot(long_format, aes(x=Sample, y=Log_Count, fill = Individual)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=1,
               outlier.size=2) + ggtitle( 'Boxplot showing variation in count reads in each Ebola Group Sample') + theme(axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))


rlogcounts <- vst(as.matrix(counts))
rlog_ebola_group_counts <- vst(as.matrix(ebola_group_counts))
rlog_cluster_group_counts <- vst(as.matrix(Before_and_After_infected_ebola_counts))

View(Before_and_After_infected_ebola_counts)
#I then ran PCA on the groups of all samples, ebola post vaccination and all the ebola group which showed some time dependent variation, with the greatest variance shown between pre and post first vaccination.


### Run PCA 
pcDat <- prcomp(t(rlogcounts))
pcDat2 <- prcomp(t(rlog_ebola_group_counts))
pcDat3 <- prcomp(t(rlog_cluster_group_counts))



# plot PCA
sample_info$time <- factor(sample_info$time, levels = c("before", "1w_ps", "2w_ps", "3w_ps", "4w_ps", "5w_ps_1w_ps", "6w_ps_2w_ps", "7w_ps_3w_ps", "8w_ps_4w_ps"))

plot <- autoplot(pcDat,
                 data = sample_info,
                 colour="time", 
                 shape = "treatment",
                 size=5)
plot + ggtitle( 'PCA Analysis for all Bats') + theme(plot.title = element_text(hjust = 0.5))

View(infected_ebola_bats)



ebola_group_bats$time <- factor(ebola_group_bats$time, levels = c("before", "1w_ps", "2w_ps", "3w_ps", "4w_ps", "5w_ps_1w_ps", "6w_ps_2w_ps", "7w_ps_3w_ps", "8w_ps_4w_ps"))

plot2 <- autoplot(pcDat2,
                  data = ebola_group_bats,
                  colour="time", 
                  size=5, label = FALSE, label.size =3)
plot2 + ggtitle( 'PCA Analysis for Ebola Group Bats') + theme(plot.title = element_text(hjust = 0.5)) + geom_text(label="1",x= -0.22,y=0.10,size = 5, fontface = "plain") + geom_text(label="2",  x= -0.16,y=-0.2,size = 5, fontface = "plain")




View(Before_and_After_infected_ebola_bats)
Before_and_After_infected_ebola_bats$time <- factor(Before_and_After_infected_ebola_bats$time, levels = c("before", "5w_ps_1w_ps"))

plot3 <- autoplot(pcDat3,
                  data = Before_and_After_infected_ebola_bats,
                  shape ="time", colour = "individual", 
                  size=5, label = FALSE, label.size =3)
plot3 + ggtitle( 'PCA Analysis for our Cluster Group Bats') + theme(plot.title = element_text(hjust = 0.5))


#I saved the loadings to try and extract highly expressed genes that could identify the genes responsible for the variance in PC1 and PC2.


pc_loadings <- pcDat2$rotation
pc_loadings_cluster <- pcDat3$rotation







