#!/usr/bin/Rscript --no-save --no-restore

# Prepare a set of plots describing the SOM results...

library(kohonen2)
library(RColorBrewer)

# Read in the SOM dataset
load("/data/mouseENCODE.new/simulations/Workspace.0.RData")

data.som = lines.deep.gt1.som

#Extract data:
data.som.info<-matrix(unlist(strsplit(rownames(data.som$data),":")),ncol=4,byrow=TRUE)
data.som.info.classes <- sort(unique(data.som.info[,1]))
data.som.bins <- apply(data.som$data, 1, sum)
som.x <- data.som$grid$xdim
som.y <- data.som$grid$ydim


############################################################################
# "Off the shelf" plots...

## General statistics on the bins:
# Plot the total number of shared TFs in each feature for each class
#for(i in 1:length(data.som.info.classes)) {
#    svg(paste(data.som.info.classes[i], ".barplot.svg", sep = ""))
#    barplot(table(data.som.bins[which(data.som.info[,1]==data.som.info.classes#[i])])[1:max(data.som.bins)], xlab="Number of Shared TFs", ylab="Number of Bin#s", col="cadetblue3", border="blue", main=data.som.info.classes[i])
#    dev.off()
#}
#rm(i)
## And plot for all classes combined
#svg("All-Classes.barplot.svg")
#barplot(table(data.som.bins)[1:max(data.som.bins)], xlab="Number of Shared TFs#", ylab="Number of Bins", col="cadetblue3", border="blue")
#dev.off()
#
##Plot standard SOM plots for information purposes
#svg("changes.svg")
#plot(data.som, type="changes")
#dev.off()
svg("counts.svg")
plot(data.som, type="counts")
dev.off()
#svg("codes.svg")
#plot(data.som, type="codes", codeRendering="segments")
#dev.off()
svg("dist.neighbours.svg")
plot(data.som, type="dist.neighbours")
dev.off()
##Quantization Error
#mean(data.som$distances)


#############################################################################
# Specialized plots...

#Build the classification vectors
data.som.info.classf = factor(data.som$unit.classif, levels=seq(1:(som.x*som.y)))
data.som.classes <- matrix(nrow=(som.x*som.y), ncol=length(data.som.info.classes))
for(i in 1:(som.x*som.y)) {
  data.som.info.class.factor <- factor(data.som.info[which(data.som.info.classf==i),1], levels=data.som.info.classes)
  temp<-table(data.som.info.class.factor)
  data.som.classes[i,] = as.vector(temp)
}
rm(temp)
rm(i)
colnames(data.som.classes)<-data.som.info.classes

# Load pattern counts from significant patterns (will be used to filter total data)
pattern.data = read.table("/home/adadiehl/mouseENCODE.new/conservation/SOM_Selection-Pattern_Data.ALL.txt", header = TRUE, row.names = "Neuron")
data.som.patterns = as.data.frame(matrix(nrow=(som.x*som.y), ncol=length(pattern.data)))
for(i in 1:(som.x*som.y)) {
      row = c(rep(NA, length(pattern.data)))

      if (any(row.names(pattern.data) == i)) {

          data.som.patterns[i,] = as.vector(pattern.data[which(row.names(pattern.data) == i),])

      } else {
          data.som.patterns[i,] = as.vector(row)
      }
}
colnames(data.som.patterns)<-colnames(pattern.data)

# Prepare filtered classes matrix (only significant patterns)
data.som.classes.sig = matrix(nrow=(som.x*som.y), ncol=ncol(data.som.classes))
data.som.classes.sig[(which(data.som.patterns$Pattern != "NA")),] = data.som.classes[ (which(data.som.patterns$Pattern != "NA")), ]
colnames(data.som.classes.sig)<-data.som.info.classes

# Plot a density of location on the map for each class. ie - The density of each cell type for broad and deep
for(i in 1:length(data.som.info.classes)) {
  data.som.ratio <- data.som.classes.sig[,i] / rowSums(data.som.classes.sig)
  data.som.ratio[is.na(data.som.ratio)] <- -1
  svg(paste("density.class.", colnames(data.som.classes.sig)[i], "-ALL.svg", sep = ""))
  plot(data.som, type="property", property=as.vector(data.som.ratio), main=paste("Ratio of",colnames(data.som.classes.sig)[i], "to All", sep=" "), zlim=c(0, 1), palette.name=colorRampPalette(brewer.pal(8,"Blues")), contin=TRUE)
  dev.off()
}
rm(i)


# Plot Selected Fraction in divided hexagons
svg("Selected-Fraction.svg")
plot.kohprop2(data.som, property1=as.vector(data.som.patterns[,7]), property2=as.vector(data.som.patterns[,8]), main="Selected Fraction", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

# Plot Coverage by PhastCons Elements
svg("PhastCons-Coverage.svg")
plot.kohprop2(data.som, property1=as.vector(data.som.patterns[,9]), property2=as.vector(data.som.patterns[,10]), main="PhastCons Coverage", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

# Plot Average, Max and Min PhastCons Scores
svg("PhastCons-Average.svg")
plot.kohprop2(data.som, property1=as.vector(data.som.patterns[,11]), property2=as.vector(data.som.patterns[,12]), main="PhastCons Average Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

svg("PhastCons-Max.svg")
plot.kohprop2(data.som, property1=as.vector(data.som.patterns[,13]), property2=as.vector(data.som.patterns[,14]), main="PhastCons Maximum Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

svg("PhastCons-Min.svg")
plot.kohprop2(data.som, property1=as.vector(data.som.patterns[,15]), property2=as.vector(data.som.patterns[,16]), main="PhastCons Minimum Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

# Plot Intersection with PhastCons Elements
svg("PhastCons-Int.svg")
plot.kohprop2(data.som, property1=as.vector(data.som.patterns[,25]), property2=as.vector(data.som.patterns[,26]), main="Average Num PhastCons Elements/Peak", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

# Plot Average, Max and Min PhyloP Scores
svg("PhyloP-Average.svg")
plot.kohprop2(data.som, property1=as.vector(data.som.patterns[,17]), property2=as.vector(data.som.patterns[,18]), main="PhyloP Average Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

svg("PhyloP-Max.svg")
plot.kohprop2(data.som, property1=as.vector(data.som.patterns[,19]), property2=as.vector(data.som.patterns[,20]), main="PhyloP Maximum Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

svg("PhyloP-Min.svg")
plot.kohprop2(data.som, property1=as.vector(data.som.patterns[,21]), property2=as.vector(data.som.patterns[,22]), main="PhyloP Minimum Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()


# Plot filtered SOM Pattern usage in RGB color-blend plot
svg("Pattern-Usage.svg")
plot.kohpropblend(data.som, property1 = as.vector(data.som.patterns[,3]), property2 = as.vector(data.som.patterns[,4]), main = "Relative Pattern Usage in Human and Mouse", color1="blue", color2="red", contin = TRUE, contin2 = TRUE, keepMargins = TRUE, heatkeywidth = 0.5, heatkeytype = "2d", heatkeyname = "Human", heatkeyname2 = "Mouse", alpha = TRUE, alphatrans = TRUE)
dev.off()


# Plot the pie chart showing ratio of each class in each neuron
svg("Pattern-Usage-Pie.svg")
plot.kohclasses(data.som, main = "Pattern Usage in Four Mammalian Cell Types", keepMargins = TRUE, bgcol=NULL, codeRendering = NULL, property = data.som.classes.sig, scale = FALSE, fg = "white", fgPie = "transparent", bg = "white", palette.name = colorRampPalette(c(rgb(135,206,250,153,maxColorValue=255), rgb(0,0,255,153,maxColorValue=255), rgb(240,128,128,153,maxColorValue=255), rgb(255,0,0,153,maxColorValue=255) )))
dev.off()


# Plot the density location of each code seperately. ie - The density of each TF
par(ask=T) 
for(i in 1:length(colnames(data.som$codes))) {
    code.data = matrix(nrow=(som.x*som.y), ncol=1)
    code.data[(which(data.som.patterns$Pattern != "NA")),] = data.som$codes[(which(data.som.patterns$Pattern != "NA")), colnames(data.som$codes)[i]]    
    svg(paste(colnames(data.som$codes)[i], ".svg", sep = ""))
    plot(data.som, type="property", property=code.data, main=paste(colnames(data.som$codes)[i], " Enrichment", sep=" "), zlim=c(0, 1), palette.name=colorRampPalette(brewer.pal(8,"Blues")), contin=TRUE)
    dev.off()
}
par(ask=F) 
rm(i)


# Plot a combined plot of ratio of two classes only
if(dim(data.som.classes.sig)[2] == 2) { # Two classes only!
    svg(paste("ratio.", colnames(data.som.classes.sig)[1], '.', colnames(data.som.classes.sig)[2], '.svg', sep = ""))
    plot(data.som, type="property", property=as.vector(data.som.ratio), main=paste("Ratio of",colnames(data.som.classes.sig)[1], "to", colnames(data.som.classes.sig)[2], sep=" "), zlim=c(0, 1), palette.name=colorRampPalette(brewer.pal(11,"RdYlBu")), contin=TRUE, keepMargins=TRUE)
    dev.off()
}


# Plot TF's per neuron
data.som.ntf = as.data.frame(matrix(nrow=(som.x*som.y)), ncol = 1)
for(i in 1:(som.x*som.y)) {
      if (!is.na(data.som.patterns[i,1])) {
      	   data.som.ntf[i,] = length(colnames(data.som$codes)[which(data.som$codes[i,] > .9)])
      }
}
svg("TF-count.svg")
plot(data.som, type="property", property=as.vector(data.som.ntf[,1]), main="TFs per Neuron", contin = TRUE, keepMargins = TRUE, heatkey = TRUE)
dev.off()


## TSS Distance Plots
tss_dist_avg = log2(data.som.patterns[,23] + data.som.patterns[,24] / 2)
tss_dist_avg[ which(tss_dist_avg == -Inf) ] = 0
svg("TssDist_Median.svg")
# Overall Average
plot(data.som, type="property", property=as.vector(tss_dist_avg), main="Median Log2 Distance to Nearest TSS", palette.name=colorRampPalette(brewer.pal(8,"Blues")), contin=TRUE)
dev.off()
# Species-wise averages in divided plots
v1 = log2(data.som.patterns[,23])
v1[ which(v1 == -Inf) ] = 0
v2 = log2(data.som.patterns[,24])
v2[ which(v2 == -Inf) ] = 0
svg("TssDist_HgMm_Median.svg")
plot.kohprop2(data.som, property1=as.vector(v1), property2=as.vector(v2), main="Median Log2 Distance to Nearest TSS", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()


###
# Histone Modification Plots
hist_data = read.table("/home/adadiehl/mouseENCODE.new/histone_mod/som-plot-data.tsv", header = TRUE)

svg("H3K9ac.svg")
plot.kohprop2(data.som, property1=as.vector(hist_data$H3K9ac_hg19), property2=as.vector(hist_data$H3K9ac_mm9), main="H3K9 Acetylation Median Max Log-2 Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

svg("H3K79me2.svg")
plot.kohprop2(data.som, property1=as.vector(hist_data$H3K79me2_hg19), property2=as.vector(hist_data$H3K79me2_mm9), main="H3K79 Dimethylation Median Max Log-2 Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

svg("H3K4me3.svg")
plot.kohprop2(data.som, property1=as.vector(hist_data$H3K4me3_hg19), property2=as.vector(hist_data$H3K4me3_mm9), main="H3K4 Trimethylation Median Max Log-2 Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

svg("H3K4me1.svg")
plot.kohprop2(data.som, property1=as.vector(hist_data$H3K4me1_hg19), property2=as.vector(hist_data$H3K4me1_mm9), main="H3K4 Monomethylation Median Max Log-2 Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

svg("H3K36me3.svg")
plot.kohprop2(data.som, property1=as.vector(hist_data$H3K36me3_hg19), property2=as.vector(hist_data$H3K36me3_mm9), main="H3K36 Trimethylation Median Max Log-2 Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

svg("H3K27me3.svg")
plot.kohprop2(data.som, property1=as.vector(hist_data$H3K27me3_hg19), property2=as.vector(hist_data$H3K27me3_mm9), main="H3K27 Trimethylation Median Max Log-2 Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()

svg("H3K27ac.svg")
plot.kohprop2(data.som, property1=as.vector(hist_data$H3K27ac_hg19), property2=as.vector(hist_data$H3K27ac_mm9), main="H3K27 Acetylation Median Max Log-2 Score", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "Human", heatkeyname2 = "Mouse")
dev.off()


###
# Gene Expression Plots

# TPM data
expr_data = read.table("/home/adadiehl/mouseENCODE.new/gene_expr/neurons.gene_expr.kallisto.tsv", header = TRUE)

svg("gene-expr_GM12878.fpkm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_GM12878+1)), main="GM12878 Median Log2 TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_CH12.fpkm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_CH12+1)), main="CH12 Median Log2 TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_K562.fpkm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_K562+1)), main="K562 Median Log2 TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_MEL.fpkm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_MEL+1)), main="MEL Median Log2 TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_CH12-MEL.fpkm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_CH12+1)), property2=as.vector(log2(expr_data$FPKM_MEL+1)), main="Median Log2 TPM, CH12 vs MEL", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "CH12", heatkeyname2 = "MEL")
dev.off()

svg("gene-expr_GM12878-K562.fpkm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_GM12878+1)), property2=as.vector(log2(expr_data$FPKM_K562+1)), main="Median Log2 TPM, GM12878 vs K562", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "GM12878", heatkeyname2 = "K562")
dev.off()

svg("gene-expr_K562-MEL.fpkm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_K562+1)), property2=as.vector(log2(expr_data$FPKM_MEL+1)), main="Median Log2 TPM, K562 vs MEL", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "K562", heatkeyname2 = "MEL")
dev.off()

svg("gene-expr_GM12878-CH12.fpkm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_GM12878+1)), property2=as.vector(log2(expr_data$FPKM_CH12+1)), main="Median Log2 TPM, GM12878 vs CH12", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "GM12878", heatkeyname2 = "CH12")
dev.off()

svg("gene-expr_fc_K562-GM12878.fpkm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_GM12878_K562), main="Median Log2 TPM Fold Change, K562/GM12878", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_GM12878_K562), heatkeyname = "K562", heatkeyname2 = "GM12878")
dev.off()

svg("gene-expr_fc_MEL-CH12.fpkm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_CH12_MEL), main="Median Log2 TPM Fold Change, MEL/CH12", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_CH12_MEL), heatkeyname = "CH12", heatkeyname2 = "MEL")
dev.off()

svg("gene-expr_fc_CH12-GM12878.fpkm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_CH12_GM12878), main="Median Log2 TPM Fold Change, CH12/GM12878", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_CH12_GM12878), heatkeyname = "GM12878", heatkeyname2 = "CH12")
dev.off()

svg("gene-expr_fc_MEL-K562.fpkm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_MEL_K562), main="Median Log2 TPM Fold Change, MEL/K562", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_MEL_K562), heatkeyname = "K562", heatkeyname2 = "MEL")
dev.off()


# Batch normalized data
expr_data = read.table("/home/adadiehl/mouseENCODE.new/gene_expr/neurons.gene_expr.bnorm.tsv", header = TRUE)

svg("gene-expr_GM12878.bnorm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_GM12878+1)), main="GM12878 Median Log2 Batch Normalized TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_CH12.bnorm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_CH12+1)), main="CH12 Median Log2 Batch Normalized TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_K562.bnorm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_K562+1)), main="K562 Median Log2 Batch Normalized TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_MEL.bnorm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_MEL+1)), main="MEL Median Log2 Batch Normalized TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_CH12-MEL.bnorm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_CH12+1)), property2=as.vector(log2(expr_data$FPKM_MEL+1)), main="Median Log2 Batch Norm TPM, CH12 vs MEL", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "CH12", heatkeyname2 = "MEL")
dev.off()

svg("gene-expr_GM12878-K562.bnorm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_GM12878+1)), property2=as.vector(log2(expr_data$FPKM_K562+1)), main="Median Log2 Batch Norm TPM, GM12878 vs K562", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "GM12878", heatkeyname2 = "K562")
dev.off()

svg("gene-expr_K562-MEL.bnorm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_K562+1)), property2=as.vector(log2(expr_data$FPKM_MEL+1)), main="Median Log2 Batch Norm TPM, K562 vs MEL", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "K562", heatkeyname2 = "MEL")
dev.off()

svg("gene-expr_GM12878-CH12.bnorm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_GM12878+1)), property2=as.vector(log2(expr_data$FPKM_CH12+1)), main="Median Log2 Batch Norm TPM, GM12878 vs CH12", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "GM12878", heatkeyname2 = "CH12")
dev.off()

svg("gene-expr_fc_K562-GM12878.bnorm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_GM12878_K562), main="Median Log2 Bnorm TPM Fold Change, K562/GM12878", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_GM12878_K562), heatkeyname = "K562", heatkeyname2 = "GM12878")
dev.off()

svg("gene-expr_fc_MEL-CH12.bnorm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_CH12_MEL), main="Median Log2 Bnorm TPM Fold Change, MEL/CH12", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_CH12_MEL), heatkeyname = "CH12", heatkeyname2 = "MEL")
dev.off()

svg("gene-expr_fc_CH12-GM12878.bnorm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_CH12_GM12878), main="Median Log2 Bnorm TPM Fold Change, CH12/GM12878", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_CH12_GM12878), heatkeyname = "GM12878", heatkeyname2 = "CH12")
dev.off()

svg("gene-expr_fc_MEL-K562.bnorm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_MEL_K562), main="Median Log2 Bnorm TPM Fold Change, MEL/K562", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_MEL_K562), heatkeyname = "K562", heatkeyname2 = "MEL")
dev.off()

# Quantile normalized data
expr_data = read.table("/home/adadiehl/mouseENCODE.new/gene_expr/neurons.gene_expr.kallisto.qnorm.tsv", header = TRUE)

svg("gene-expr_GM12878.qnorm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_GM12878+1)), main="GM12878 Median Log2 Quantile Normalized TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_CH12.qnorm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_CH12+1)), main="CH12 Median Log2 Quantile Normalized TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_K562.qnorm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_K562+1)), main="K562 Median Log2 Quantile Normalized TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_MEL.qnorm.svg")
plot(data.som, type="property", property=as.vector(log2(expr_data$FPKM_MEL+1)), main="MEL Median Log2 Quantile Normalized TPM", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off()

svg("gene-expr_CH12-MEL.qnorm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_CH12+1)), property2=as.vector(log2(expr_data$FPKM_MEL+1)), main="Median Log2 Quantile Norm TPM, CH12 vs MEL", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "CH12", heatkeyname2 = "MEL")
dev.off()

svg("gene-expr_GM12878-K562.qnorm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_GM12878+1)), property2=as.vector(log2(expr_data$FPKM_K562+1)), main="Median Log2 Quantile Norm TPM, GM12878 vs K562", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "GM12878", heatkeyname2 = "K562")
dev.off()

svg("gene-expr_K562-MEL.qnorm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_K562+1)), property2=as.vector(log2(expr_data$FPKM_MEL+1)), main="Median Log2 Quantile Norm TPM, K562 vs MEL", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "K562", heatkeyname2 = "MEL")
dev.off()

svg("gene-expr_GM12878-CH12.qnorm.svg")
plot.kohprop2(data.som, property1=as.vector(log2(expr_data$FPKM_GM12878+1)), property2=as.vector(log2(expr_data$FPKM_CH12+1)), main="Median Log2 Quantile Norm TPM, GM12878 vs CH12", palette.name1=colorRampPalette(brewer.pal(8,"Blues")), palette.name2=colorRampPalette(brewer.pal(8,"Reds")), heatkey = TRUE, contin1=TRUE, contin2=TRUE, keepMargins = TRUE, heatkeywidth=0.5, heatkeyname1 = "GM12878", heatkeyname2 = "CH12")
dev.off()

svg("gene-expr_fc_K562-GM12878.qnorm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_GM12878_K562), main="Median Log2 Qnorm TPM Fold Change, K562/GM12878", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_GM12878_K562), heatkeyname = "K562", heatkeyname2 = "GM12878")
dev.off()

svg("gene-expr_fc_MEL-CH12.qnorm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_CH12_MEL), main="Median Log2 Qnorm TPM Fold Change, MEL/CH12", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_CH12_MEL), heatkeyname = "CH12", heatkeyname2 = "MEL")
dev.off()

svg("gene-expr_fc_CH12-GM12878.qnorm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_CH12_GM12878), main="Median Log2 Qnrom TPM Fold Change, CH12/GM12878", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_CH12_GM12878), heatkeyname = "GM12878", heatkeyname2 = "CH12")
dev.off()

svg("gene-expr_fc_MEL-K562.qnorm.svg")
plot(data.som, type="property", property=as.vector(expr_data$FC_MEL_K562), main="Median Log2 Qnorm TPM Fold Change, MEL/K562", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(rev(brewer.pal(8,"RdBu"))), zlim = get_zero_centered_zlim(expr_data$FC_MEL_K562), heatkeyname = "K562", heatkeyname2 = "MEL")
dev.off()


# GWAS Association plots
gwas_data = read.table("/home/adadiehl/mouseENCODE.new/gwas/neurons-gwas.txt", header = TRUE)

svg("gwas_count.svg")
plot(data.som, type="property", property=as.vector(log2(gwas_data$gwas_count+1)), main="Log2 GWAS Count", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off();

svg("gwas_immune-count.svg")
plot(data.som, type="property", property=as.vector(log2(gwas_data$gwas_immune_count+1)), main="Log2 GWAS Immune-Related Count", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off();

svg("gwas_norm-count.svg")
plot(data.som, type="property", property=as.vector(log2(gwas_data$gwas_norm_count+1)), main="Log2 GWAS Count / N Human Peaks", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off();

svg("gwas_norm-count-immune.svg")
plot(data.som, type="property", property=as.vector(log2(gwas_data$gwas_norm_immune_count+1)), main="Log2 GWAS Immune Count / N Human Peaks", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")))
dev.off();

# DNase Hypersensitivity plots
dnase_data = read.table("/home/adadiehl/mouseENCODE.new/dnase1/dnase_stats.neurons.txt", header = TRUE)

svg("dnase_GM12878.svg")
plot(data.som, type="property", property=as.vector(dnase_data$GM12878), main="Human GM12878 Avg DNase Sensitivity Coverage", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

svg("dnase_K562.svg")
plot(data.som, type="property", property=as.vector(dnase_data$K562), main="Human K562 Avg DNase Sensitivity Coverage", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

svg("dnase_CH12.svg")
plot(data.som, type="property", property=as.vector(dnase_data$CH12), main="Mouse CH12 Avg DNase Sensitivity Coverage", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

svg("dnase_MEL.svg")
plot(data.som, type="property", property=as.vector(dnase_data$MEL), main="Mouse MEL Avg DNase Sensitivity Coverage", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

svg("dnase_mm9-GM12878-orth.svg")
plot(data.som, type="property", property=as.vector(dnase_data$mm9_GM12878_orth), main="Mouse Orthologs in GM12878 Avg DNase Sens. Covg", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

svg("dnase_mm9-K562-orth.svg")
plot(data.som, type="property", property=as.vector(dnase_data$mm9_K562_orth), main="Mouse Orthologs in K562 Avg DNase Sens. Covg", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

svg("dnase_hg19-CH12-orth.svg")
plot(data.som, type="property", property=as.vector(dnase_data$hg19_CH12_orth), main="Human Orthologs in CH12 Avg DNase Sens. Covg", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

svg("dnase_hg19-MEL-orth.svg")
plot(data.som, type="property", property=as.vector(dnase_data$hg19_MEL_orth), main="Human Orthologs in MEL Avg DNase Sens. Covg", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

svg("dnase_GM12878_all.svg")
plot(data.som, type="property", property=as.vector(dnase_data$GM12878_all), main="GM12878 Composite Avg DNase Sens. Covg", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

svg("dnase_K562_all.svg")
plot(data.som, type="property", property=as.vector(dnase_data$K562_all), main="K562 Composite Avg DNase Sens. Covg", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

svg("dnase_CH12_all.svg")
plot(data.som, type="property", property=as.vector(dnase_data$CH12_all), main="CH12 Composite Avg DNase Sens. Covg", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

svg("dnase_MEL_all.svg")
plot(data.som, type="property", property=as.vector(dnase_data$MEL_all), main="MEL Composite Avg DNase Sens. Covg", contin = TRUE, keepMargins = TRUE, heatkey = TRUE, palette.name=colorRampPalette(brewer.pal(8,"Blues")), zlim=c(0,1))
dev.off();

# Plot of SOM neuron classifications
svg("som_neuron_classes.svg")
plot.somclasses(data.som, property, main="Pattern Classifications", contin=TRUE, heatkey=TRUE, palette.name=colorRampPalette(c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAEDAB", "#F7EE55", "#F6C141", "#F1932D", "#EB601C", "#DC050C")))
dev.off()
