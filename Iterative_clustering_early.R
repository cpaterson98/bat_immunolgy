library(mclust)
library(tidyverse)
library(WGCNA)

sample_info <- read_csv("Data/colData_all_samples_included_27.02.2020.csv")
counts <- read_csv("Data/counttable_all_samples_included_27.02.2020.csv")
gene_info <- read_csv("Data/gene_biological_info_gff_and_RBH_PVAP_MLUC_PALE_RAEG.csv")

View(counts_b)
counts_m <- counts %>% 
  as.data.frame() %>% 
  column_to_rownames("X1") %>% 
  as.matrix() 
total_counts <- as.data.frame(rowSums(counts_m))

counts <- counts %>% 
  as.data.frame() %>% 
  column_to_rownames("X1")

counts_2 <- cbind(counts, total_counts)


row.names.remove <- c("EHELV_G000003486","EHELV_G00013550", "EHELV_G00013551", "EHELV_G00013552", "EHELV_G00015016","EHELV_G00027280", "EHELV_G00031256", "EHELV_G00032092", "EHELV_G00032093", "EHELV_G00034890", "EHELV_G00035000")

counts <-  counts[!(row.names(counts) %in% row.names.remove), ]

keep <-  rowSums(counts) >600
counts <- counts[keep,]
librarySizes <- colSums(counts)

sample_normalised_counts <- counts/librarySizes
sample_normalised_counts <- sample_normalised_counts*1000000

norm_counts <- log2(sample_normalised_counts +1)
counts <- norm_counts %>% 
  as.data.frame() %>% 
  rename_all(str_remove, "_ReadsPerGene.out")

sample_info <- sample_info %>% 
  column_to_rownames("X1") %>% 
  t() %>% 
  as.data.frame() %>% 
  rename_all(str_remove, "_ReadsPerGene.out") %>% 
  t()
sample_info <- sample_info %>% 
  as.data.frame()

counts_t <- t(counts)


Before_and_After_infected_ebola_bats <- sample_info %>% 
  as.data.frame()
Before_and_After_infected_ebola_bats <- Before_and_After_infected_ebola_bats[Before_and_After_infected_ebola_bats$treatment == "ebola",]
Before_and_After_infected_ebola_bats <- Before_and_After_infected_ebola_bats[Before_and_After_infected_ebola_bats$day == "A0711"| Before_and_After_infected_ebola_bats$day == "A0807"  ,]

View(Before_and_After_infected_ebola_bats)
Before_and_After_infected_ebola_bats <- Before_and_After_infected_ebola_bats[-c(1,14),]

Before_and_After_infected_ebola_bats <- Before_and_After_infected_ebola_bats %>% 
  rownames_to_column("Long_Sample")


Before_and_After_infected_ebola_counts <- counts_t%>% 
  as.data.frame() %>% 
  rownames_to_column("Long_Sample") %>% 
  filter(Long_Sample %in% Before_and_After_infected_ebola_bats$Long_Sample) %>% 
  column_to_rownames("Long_Sample") %>% 
  t()

counts <- Before_and_After_infected_ebola_counts


counts_m <- counts %>% 
  as.matrix() 
total_counts <- as.data.frame(rowSums(counts_m))


counts_2 <- cbind(counts, total_counts)















## Mcluster 1

target_group_counts <- Before_and_After_infected_ebola_counts
target_group_bats <- Before_and_After_infected_ebola_bats
target_group_counts <- target_group_counts %>% 
  as.data.frame() %>% 
  rename_all(str_remove, "ebola_")

initial_cluster <- Mclust(target_group_counts)
summary(initial_cluster)
class <- initial_cluster$classification
X <-  target_group_counts
clPairs(X, class, cex.labels = 0.5)

#BIC <- mclustBIC(X)
#plot(BIC)
#summary(BIC)
#mod1 <- Mclust(X, x = BIC)
#summary(mod1, parameters = TRUE)
#plot(mod1, what = "classification")





## WGCNA Cluster 1

datExpr0 <- as.data.frame(t(target_group_counts)); 

sample_data <- target_group_bats


#Then we run the process:

### Checking for Excessive Missing Values and identifying of outliers
gsg = goodSamplesGenes(datExpr0, verbose = 3); 
gsg$allOK
if (!gsg$allOK) {
  
  # Optionally, print the gene and sample names that were removed: if (sum(!gsg$goodGenes)>0)
  
  printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", "))); if (sum(!gsg$goodSamples)>0)
    
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", "))); # Remove the offending genes and samples from the data:
  
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = fastcluster::hclust(dist(datExpr0), method = "average")

sizeGrWindow(12,9) 

par(cex = 0.6); 
par(mar = c(0,4,2,0)) 
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 300, col = "red");


#Determine the cut height if there are any outliers

cut_height <- 400

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut_height, minSize = 2)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]

nGenes = ncol(datExpr0) 
nSamples = nrow(datExpr0)
traitData <- sample_data %>% 
  column_to_rownames("Long_Sample")
traitRows<-  rownames(datExpr0)
datTraits <- traitData[traitRows,6]



datTraits <- traitData[traitRows,6]


collectGarbage();

# Re-cluster samples 
sampleTree2 <-  fastcluster::hclust(dist(datExpr0), method = "average") 
# Convert traits to a color representation: white means low, red means high, grey means missing entry 
traitColors <-  numbers2colors(as.numeric(datTraits), signed = FALSE); 
# Plot the sample dendrogram and the colors underneath.








options(stringsAsFactors = FALSE); # Allow multi-threading within WGCNA. This helps speed up certain calculations.


# Choose a set of soft-thresholding powers 
powers = c(c(1:10), seq(from = 12, to=40, by=2)) 
# Call the network topology analysis function 
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5) 
# Plot the results:

sizeGrWindow(24, 24) 
par(mfrow = c(1,2)); 
cex1 = 0.9; 
# Scale-free topology fit index as a function of the soft-thresholding power 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")); text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                                               labels=powers,cex=cex1,col="red"); 
# this line corresponds to using an R^2 cut-off of h 
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power 
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



cor <- WGCNA::cor
net = blockwiseModules(datExpr0, power = 20,
                       TOMType = "unsigned", minModuleSize = 80, maxModuleSize = 1000, 
                       reassignThreshold = 0, mergeCutHeight = 0.1, 
                       numericLabels = TRUE, pamRespectsDendro = FALSE, 
                       saveTOMs = TRUE, saveTOMFileBase = "BatGenesTOM", 
                       verbose = 3)

table(net$colors)

net$colors
# Convert labels to colors for plotting 
mergedColors <-  labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath 
sizeGrWindow(12, 9) 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "DendronTree for Round 1 of WGCNA")

cutree(net$dendrograms[[1]] , h = 0.9)
net$dendrograms
moduleLabels <-  net$colors 
moduleColors <-  labels2colors(net$colors) 
MEs <-  net$MEs; 
geneTree <-  net$dendrograms[[1]]; 
save(MEs, moduleLabels, moduleColors, geneTree, 
  file = "cluster_4w_autonetwork.RData")



## Get genes of high module membership 

geneModuleMembership <-  as.data.frame(cor(datExpr0, MEs, use = "p"))

geneModuleMembership[, "max"] <- apply(geneModuleMembership[,], 1, max)

filtered_genes <- target_group_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  mutate(classification_mclust = initial_cluster$classification) %>% 
  mutate(classification_WGCNA = net$colors) 


filtered_genes_WCCNA <- filtered_genes %>% 
  filter(geneModuleMembership$max >= 0.6)
dim(filtered_genes_WCCNA)


filtered_genes_mclust <-  filtered_genes %>% 
  filter(initial_cluster$uncertainty <= 0.3)
dim(filtered_genes_mclust)

target_group_bats_rownames <- target_group_bats %>% 
  column_to_rownames("Long_Sample")


row.names(target_group_bats_rownames)



target_group_counts_2 <- inner_join(filtered_genes_mclust, filtered_genes_WCCNA, c("gene", "A0711_10_5900074", "A0711_11_5896894", "A0711_12_5896826",
                                                                                   "A0711_16_5897207", "A0711_17_5903551", "A0711_18_5895929",
                                                                                    "A0807_09_5900074", "A0807_10_5896826", "A0807_11_5896894",
                                                                                    "A0807_16_5895929", "A0807_17_5897207", "A0807_18_5903551", "classification_WGCNA", "classification_mclust"))
#"A0711_09_5900117"
View(target_group_counts_2)
dim(target_group_counts_2)


target_group_counts_2 <- target_group_counts_2[ -c(14,15)]
target_group_counts_2 <- target_group_counts_2 %>% 
  column_to_rownames("gene")






### ROUND 2 ####









### Mclust 2

secondary_cluster <- Mclust(target_group_counts_2)
summary(secondary_cluster)
class <- secondary_cluster$classification
table(class)
X <-  target_group_counts_2
clPairs(X, class, cex.labels = 0.5)
#BIC <- mclustBIC(X)
#plot(BIC)
#summary(BIC)
#mod1 <- Mclust(X, x = BIC)
#summary(mod1, parameters = TRUE)
#plot(mod1, what = "classification")




## WGCNA 2
datExpr0 <- as.data.frame(t(target_group_counts_2)); 
sample_data <- target_group_bats

#Then we run the process:

### Checking for Excessive Missing Values and identifying of outliers
gsg = goodSamplesGenes(datExpr0, verbose = 3); 
gsg$allOK
if (!gsg$allOK) {
  
  # Optionally, print the gene and sample names that were removed: if (sum(!gsg$goodGenes)>0)
  
  printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", "))); if (sum(!gsg$goodSamples)>0)
    
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", "))); # Remove the offending genes and samples from the data:
  
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = fastcluster::hclust(dist(datExpr0), method = "average")

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches 
# The user should change the dimensions if the window is too large or too small. 
#sizeGrWindow(12,9) 
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9); 
#par(cex = 0.6); 
#par(mar = c(0,4,2,0)) 
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
#    cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
#abline(h = 0.3, col = "red");


#Determine the cut height if there are any outliers

#cut_height <- 0.3

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut_height, minSize = 2)
#table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]

nGenes = ncol(datExpr0) 
nSamples = nrow(datExpr0)

traitData <- sample_data %>% 
  column_to_rownames("Long_Sample")
traitRows<-  rownames(datExpr0)
datTraits <- traitData[,6]

collectGarbage();

# Re-cluster samples 
sampleTree2 <-  fastcluster::hclust(dist(datExpr0), method = "average") 
# Convert traits to a color representation: white means low, red means high, grey means missing entry 
traitColors <-  numbers2colors(as.numeric(datTraits), signed = FALSE); 
# Plot the sample dendrogram and the colors underneath.

#plotDendroAndColors(sampleTree2, traitColors, 
# groupLabels =  names(datTraits), 
#main = " Inital Dendrogram and Variable heatmap")

#save(datExpr0, datTraits, file = "LateInfectedEbolaBatGeneExpression.RData")


# The following setting is important, do not omit.

options(stringsAsFactors = FALSE); # Allow multi-threading within WGCNA. This helps speed up certain calculations.


# Load the data saved in the first part 
#lnames = load(file = "LateInfectedEbolaBatGeneExpression.RData"); 

# Choose a set of soft-thresholding powers 
powers = c(c(1:10), seq(from = 12, to=40, by=2)) 
# Call the network topology analysis function 
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5) 
# Plot the results:

sizeGrWindow(24, 24) 
par(mfrow = c(1,2)); 
cex1 = 0.9; 
# Scale-free topology fit index as a function of the soft-thresholding power 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")); text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red"); 
 # this line corresponds to using an R^2 cut-off of h 
 abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power 
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

                                              

                                               
                                               
  cor <- WGCNA::cor
 net = blockwiseModules(datExpr0, power = 20,
  TOMType = "unsigned", minModuleSize = 80, maxModuleSize = 1000, 
 reassignThreshold = 0, mergeCutHeight = 0.1, 
 numericLabels = TRUE, pamRespectsDendro = FALSE, 
 saveTOMs = TRUE, saveTOMFileBase = "BatGenesTOM", 
verbose = 3)
                                               
 table(net$colors)
 
 
 # Convert labels to colors for plotting 
 mergedColors <-  labels2colors(net$colors)
 # Plot the dendrogram and the module colors underneath 
 sizeGrWindow(12, 9) 
 plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],  "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "DendronTree for Round 2 of WGCNA")
 
 
 moduleLabels <-  net$colors 
 moduleColors <-  labels2colors(net$colors) 
 MEs <-  net$MEs; 
 geneTree <-  net$dendrograms[[1]]; 
 save(MEs, moduleLabels, moduleColors, geneTree, 
  file = "cluster_4w_autonetwork_2.RData")
 
 
 
 geneModuleMembership2 <-  as.data.frame(cor(datExpr0, MEs, use = "p"))
 
 geneModuleMembership2[, "max"] <- apply(geneModuleMembership2[, ], 1, max)
 
 
 filtered_genes_2 <- target_group_counts_2 %>% 
   as.data.frame() %>% 
   rownames_to_column("gene") %>% 
   mutate(classification_WGCNA = net$colors) %>% 
   mutate(classification_mclust = secondary_cluster$classification)
 
 
 filtered_genes_WCCNA_2 <- filtered_genes_2 %>% 
   filter(geneModuleMembership2$max >= 0.6)
 
 filtered_genes_mclust_2 <-  filtered_genes_2 %>% 
   filter(secondary_cluster$uncertainty <= 0.3)
 
 target_group_counts_3 <- inner_join(filtered_genes_mclust_2, filtered_genes_WCCNA_2, c("gene", "A0711_10_5900074", "A0711_11_5896894", "A0711_12_5896826",
                                                                                        "A0711_16_5897207", "A0711_17_5903551", "A0711_18_5895929",
                                                                                        "A0807_09_5900074", "A0807_10_5896826", "A0807_11_5896894",
                                                                                        "A0807_16_5895929", "A0807_17_5897207", "A0807_18_5903551", "classification_WGCNA", "classification_mclust"))
 

 dim(target_group_counts_3)
 View(target_group_counts_3)
 
 
 
 
 View(target_group_counts_3)
 
 
 
 
 
 
 ##### FINAL SORTING #######
 
 View(datExpr0)
 datExpr0 <- target_group_counts_3[,-c(14,15)] %>% 
   column_to_rownames("gene") %>% 
   t() %>% 
   as.data.frame()
   
 sample_data <- target_group_bats
 
 #Then we run the process:
 
 ### Checking for Excessive Missing Values and identifying of outliers
 gsg = goodSamplesGenes(datExpr0, verbose = 3); 
 gsg$allOK
 if (!gsg$allOK) {
   
   # Optionally, print the gene and sample names that were removed: if (sum(!gsg$goodGenes)>0)
   
   printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", "))); if (sum(!gsg$goodSamples)>0)
     
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", "))); # Remove the offending genes and samples from the data:
   
   datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
 }

 
 nGenes = ncol(datExpr0) 
 nSamples = nrow(datExpr0)
 
 traitData <- sample_data %>% 
   column_to_rownames("Long_Sample")
 traitRows<-  rownames(datExpr0)
 datTraits <- traitData[,6]
 View(datTraits)
 View(traitData)
 collectGarbage();
 
 # Re-cluster samples 
 sampleTree2 <-  fastcluster::hclust(dist(datExpr0), method = "average") 
 # Convert traits to a color representation: white means low, red means high, grey means missing entry 
 traitColors <-  numbers2colors(as.numeric(datTraits), signed = FALSE); 
 
 
 cor <- WGCNA::cor
 net = blockwiseModules(datExpr0, power = 20,
                        TOMType = "unsigned", minModuleSize = 80, maxModuleSize = 1000, 
                        reassignThreshold = 0, mergeCutHeight = 0.1, 
                        numericLabels = TRUE, pamRespectsDendro = FALSE, 
                        saveTOMs = TRUE, saveTOMFileBase = "BatGenesTOM", 
                        verbose = 3)
 
 table(net$colors)
 
 
 # Convert labels to colors for plotting 
 mergedColors <-  labels2colors(net$colors)
 # Plot the dendrogram and the module colors underneath 
 sizeGrWindow(12, 9) 
 plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],  "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "DendronTree for Round 3 of WGCNA")
 table(net$colors)
 
 moduleLabels <-  net$colors 
 moduleColors <-  labels2colors(net$colors) 
 MEs <-  net$MEs; 
 geneTree <-  net$dendrograms[[1]]; 
 save(MEs, moduleLabels, moduleColors, geneTree, 
      file = "cluster_4w_autonetwork_final.RData")
 
 
 
 geneModuleMembership2 <-  as.data.frame(cor(datExpr0, MEs, use = "p"))
 
 geneModuleMembership2[, "max"] <- apply(geneModuleMembership2[, ], 1, max)
 
 View(target_group_counts_3)
 filtered_genes_3 <- target_group_counts_3[,-c(15)] %>% 
   as.data.frame() %>% 
   column_to_rownames("gene")
   

 
 
 final_count_list <- target_group_counts_3[,-c(14,15)] %>% 
   as.data.frame() %>% 
   column_to_rownames("gene") %>% 
   rownames_to_column("EHELV_gene_ID") %>% 
   mutate("classification_WGCNA" = net$colors)
 
 

 

 
 View(final_count_list)


 

 
                                               
               
                                               
                                               