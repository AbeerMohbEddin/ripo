# Set working directory
setwd("/mnt/c/BGAproject/data")

# Load the required libraries
library(edgeR)
library(biomaRt)
library(enrichR)
library(ggrepel)
library(msigdbr)
library(tidyverse)
library(DOSE)
library(clusterProfiler)
library(AnnotationDbi)

# ---- Normalization Step ----

# Load and normalize data for NP
countfile_NP <- read.csv("CGE_NP.csv", row.names = 1)
counts_NP <- as.matrix(countfile_NP)
y_NP <- DGEList(counts = counts_NP)
dge_NP <- calcNormFactors(y_NP)
normcounts_log2cpm_NP <- cpm(dge_NP, log = TRUE, prior.count = 0.25)
write.csv(normcounts_log2cpm_NP, file = "normcounts_log2cpm_NP.csv")

# Load and normalize data for PE
countfile_PE <- read.csv("CGE_PE.csv", row.names = 1)
counts_PE <- as.matrix(countfile_PE)
y_PE <- DGEList(counts = counts_PE)
dge_PE <- calcNormFactors(y_PE)
normcounts_log2cpm_PE <- cpm(dge_PE, log = TRUE, prior.count = 0.25)
write.csv(normcounts_log2cpm_PE, file = "normcounts_log2cpm_PE.csv")

# ---- Pathway Analysis Step ----

# Load the NP data file
df <- read.csv("NP_MEVHIGHANDLOWNEW.csv")

# Perform pathway analysis
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
NP <- getBM(filters = 'ensembl_gene_id', attributes = c('ensembl_gene_id',"entrezgene_id",'hgnc_symbol', 'description','gene_biotype'), values = df['ENSEMBL.ID'], mart=mart, uniqueRows = TRUE)
write.csv(NP, "GENELIST.CSV")

myGene = NP$hgnc_symbol
setEnrichrSite("Enrichr")
websiteLive <- TRUE
mydbs = c("KEGG_2021_Human", "GO_Cellular_Component_2023","GO_Biological_Process_2023")

if(websiteLive) {
  enriched <- enrichr(myGene, mydbs)
}

Np_Pathtable = if (websiteLive) enriched[["KEGG_2021_Human"]]
Np_Pathtable = as.data.frame(Np_Pathtable)
write.csv(Np_Pathtable, "NP_Pathtable.csv")

# Repeat the process for PE data file
df <- read.csv("PE_MEVHIGHANDLOWNEW.csv")

# Perform pathway analysis
PE <- getBM(filters = 'ensembl_gene_id', attributes = c('ensembl_gene_id',"entrezgene_id",'hgnc_symbol', 'description','gene_biotype'), values = df['ENSEMBL.ID'], mart=mart, uniqueRows = TRUE)
write.csv(PE, "GENELIST_PE.CSV")

myGene = PE$hgnc_symbol
setEnrichrSite("Enrichr")
websiteLive <- TRUE
mydbs = c("KEGG_2021_Human", "GO_Cellular_Component_2023","GO_Biological_Process_2023")

if(websiteLive) {
  enriched <- enrichr(myGene, mydbs)
}

PE_Pathtable = if (websiteLive) enriched[["KEGG_2021_Human"]]
PE_Pathtable = as.data.frame(PE_Pathtable)
write.csv(PE_Pathtable, "PE_Pathtable.csv")

# ---- Gene Ontology Step ----

# Assign the gene ids to a variable
genes_to_testNP <- rownames(sigsNP)
genes_to_testPE <- rownames(sigsPE)

# Perform the Gene Ontology
GO_resultsNP <- enrichGO(gene = genes_to_testNP, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
GO_resultsPE <- enrichGO(gene = genes_to_testPE, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

# Visualize the result in barplot
fit_NP <- barplot(GO_resultsNP, showCategory = 10)
fit_PE <- barplot(GO_resultsPE, showCategory = 10)

# Save the plots as png files
png("BP_NP.png", res = 250, width = 1400, height = 1800)
print(fit_NP)
dev.off()

png("BP_PE.png", res = 250, width = 1400, height = 1800)
print(fit_PE)
dev.off()

# ---- Visualization Step ----

# Read Highly and Lowly Regulated data NP
HLRData <- read.csv("ConsistentlyExpGenesPE.csv", row.names = 1)

# Convert the p-value into a -log10(p-value)
ggplot(data=HLRData, aes(x=mean, y=-log10(RP.values.up.))) +
  geom_point() +
  geom_smooth()

# Change the confidence interval fill color
ggplot(data=HLRData, aes(x=mean, y=-log10(RP.values.up.))) +
  geom_point(shape=18, color="blue") +
  geom_smooth(method=lm, linetype="dashed", color="darkred", fill="blue") +
  geom_text(label=rownames(HLRData))

# Load NormalisedData for heatmap 
NormalisedData <- read.csv("", row.names = 1)
# Create the heatmap
pheatmap(NormalisedData, color = colorRampPalette(c("blue", "white", "red"))(100), 
         show_rownames = TRUE, show_colnames = TRUE)

# Visualization for PE
Data2 <- read.csv("ConsistelyExpGenesPE.csv", row.names = 1)
# Scatter plot between RP.values.up and mean
ggplot(data=Data2, aes(x=std.dev, y=RP.values.up.)) + geom_point()

# Loess method
ggplot(data=Data2, aes(x=mean, y=-log10(RP.values.up.))) +
  geom_point() +
  geom_smooth()

# Change the confidence interval fill color
ggplot(data=Data2, aes(x=mean, y=-log10(RP.values.up.))) +
  geom_point(shape=18, color="blue") +
  geom_smooth(method=lm, linetype="dashed", color="darkred", fill="blue") +
  geom_text(label=rownames(Data2))

# Create the heatmap
data3 <- Data2[ -c(3) ] # removing the third column
pheatmap(data3, color = colorRampPalette(c("blue", "white", "red"))(100), 
         show_rownames = TRUE, show_colnames = TRUE)
