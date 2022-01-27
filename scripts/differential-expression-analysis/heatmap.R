#Load libraries
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)


"
texpr:  extract transcript-level expression measurements from ballgown objects.
gexpr:  extract gene-level expression measurements from ballgown objects.
rowVars: Variance estimates for each row (column) in a matrix.
stattest: statistical tests for differential expression in ballgown Description.
hclust: Hierarchical Clustering.
dist: Distance Matrix Computation.
"

#If X11 not available, open a pdf device for output of all plots
pdf(file="Heatmap-0.05.pdf")

#Import expression and differential expression results from the HISAT2/StringTie/Ballgown pipeline
load('bg.rda')

# load gene names for lookup later
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

# Pull the gene expression data frame from the ballgown object
gene_expression = as.data.frame(gexpr(bg))

#Assign colors to each.  You can specify color by RGB, Hex code, or name
#To get a list of color names:
data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3")

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Calculate the differential expression results including significance
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

#get only those that are significant according to Ballgown
sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])

# Set the columns for finding FPKM and create shorter names for figures
data_columns=c(1:6)
short_names=c("Ad-Rbm10_rep1", "Ad-Rbm10_rep2", "Ad-Rbm10_rep3", "Ad-Gfb_rep1", "Ad-Gfb_rep2", "Ad-Gfb_rep3")

#### Plot #11 - Create a heatmap to vizualize expression differences between the eight samples

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

main_title="sig DE Transcripts"
par(cex.main=0.8)

sig_genes=results_genes[sig,"id"]
sig_gene_names=results_genes[sig,"gene_name"]

data=log2(as.matrix(gene_expression[sig_genes,data_columns])+1)

heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both",
          margins=c(6,7), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none",
          trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names,col=rev(heat.colors(75)))

dev.off()