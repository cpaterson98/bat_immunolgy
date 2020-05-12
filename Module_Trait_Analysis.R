# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
datTraits <- as.data.frame(datTraits)
View(datTraits)
names(datTraits)[names(datTraits) == "datTraits"] <- "time_sample_taken"
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#Since we have a moderately large number of modules and traits, a suitable graphical representation will help in reading the table. We color code each association by the correlation value:
  sizeGrWindow(25,25)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3, 12, 3, 16))
# Display the correlation values within a heatmap plot 
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
yColorWidth = 0.08,
yColorOffset = 0.01,
colorLabels = TRUE,
keepLegendSpace = FALSE,
cex.lab.x = 1, cex.lab.y = 0.7,xLabelsAngle = 0,
colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.8,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
View(moduleTraitCor)
