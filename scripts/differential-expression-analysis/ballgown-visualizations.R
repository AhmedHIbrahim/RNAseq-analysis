# Load libraries needed for this analysis
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

#########
#NOTE::
# ONLY for ccne1_gene, change target_gene_index wisely
#> geneNames<-ballgown::geneNames(bg)
#> ccne2_trans_index <- names(geneNames[geneNames=="CCNE2"])
#>
#> transcriptNames[ccne2_trans_index]
#########


# Define a path for the output PDF to be written
outfile="./Final_Part2_ballgown_outpu_ccne2_gene.pdf"

target_gene_index=208874

# Load phenotype data
pheno_data = read.csv("Ad-Rbm10_vs_Ad-Gfb.csv")

# Display the phenotype data
pheno_data

# Load the ballgown object from file
load('bg.rda')

# The load command, loads an R object from a file into memory in our R session. 
# You can use ls() to view the names of variables that have been loaded
ls()

# Print a summary of the ballgown object
bg

# Open a PDF file where we will save some plots. We will save all figures and then view the PDF at the end
pdf(file=outfile)

# Extract FPKM values from the 'bg' object
fpkm = texpr(bg,meas="FPKM")

# View the last several rows of the FPKM table
tail(fpkm)

# Transform the FPKM values by adding 1 and convert to a log2 scale
fpkm = log2(fpkm+1)

# View the last several rows of the transformed FPKM table
tail(fpkm)

# Create boxplots to display summary statistics for the FPKM values for each sample
boxplot(fpkm,las=2,ylab='log2(FPKM+1)')

# Display the transcript ID for a single row of data
ballgown::transcriptNames(bg)[target_gene_index]

# Display the gene name for a single row of data 
ballgown::geneNames(bg)[target_gene_index]

# Create a BoxPlot comparing the expression of a single gene for all replicates of both conditions
plot(fpkm[target_gene_index,] ~ as.factor(pheno_data$type), border=c(2,3), main=paste(ballgown::geneNames(bg)[target_gene_index],' : ', ballgown::transcriptNames(bg)[target_gene_index]),pch=19, xlab="Type", ylab='log2(FPKM+1)')

# Add the FPKM values for each sample onto the plot 
points(fpkm[target_gene_index,] ~ jitter(as.numeric(as.factor(pheno_data$type))), col=as.numeric(as.factor(pheno_data$type))+1, pch=16)

# Create a plot of transcript structures observed in each replicate and color transcripts by expression level
plotTranscripts(ballgown::geneIDs(bg)[target_gene_index], bg, main=c('Gene in sample Ad-Rbm10_Rep1'), sample=c('Ad-Rbm10_rep1'), labelTranscripts = TRUE)
plotTranscripts(ballgown::geneIDs(bg)[target_gene_index], bg, main=c('Gene in sample Ad-Rbm10_Rep2'), sample=c('Ad-Rbm10_rep2'), labelTranscripts = TRUE)
plotTranscripts(ballgown::geneIDs(bg)[target_gene_index], bg, main=c('Gene in sample Ad-Rbm10_Rep3'), sample=c('Ad-Rbm10_rep3'), labelTranscripts = TRUE)
plotTranscripts(ballgown::geneIDs(bg)[target_gene_index], bg, main=c('Gene in sample Ad-Gfb_Rep1'), sample=c('Ad-Gfb_rep1'), labelTranscripts = TRUE)
plotTranscripts(ballgown::geneIDs(bg)[target_gene_index], bg, main=c('Gene in sample Ad-Gfb_Rep2'), sample=c('Ad-Gfb_rep2'), labelTranscripts = TRUE)
plotTranscripts(ballgown::geneIDs(bg)[target_gene_index], bg, main=c('Gene in sample Ad-Gfb_Rep3'), sample=c('Ad-Gfb_rep3'), labelTranscripts = TRUE)

#plotMeans('TST',bg,groupvar="type",legend=FALSE)

# Close the PDF device where we have been saving our plots
dev.off()
