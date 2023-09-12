# script to analyse RNAseq data for ovarian cancer cell lines
# Carolin Sauer

# load libraries
require(tidyverse)
require(DESeq2)
require(fgsea)

# Note that this script was written for and run with DESeq2 v1.26.0. The code might have to be adjusted with altered parameters to be compatible with newer versions of DESeq2 and to yield the same results.

source("cms_pal.R")

# load and prepare data
## Load R object (from tximport - salmon data)
load("../data/cell_line_RNAseq_data/salmon_gene_txi.rda")

## Sample info / metadata
sample_order <- as.factor(colnames(txi$counts))

meta <- read_csv("../data/cell_line_RNAseq_data/CellLines_RNAextract_all.csv") %>% 
  filter(RNAJBLABnumber != "JBLAB-25681") %>% #JBLAB-25681 has not produced any reads and was not analysed, so we can remove this sample from the meta data
  arrange(factor(RNAJBLABnumber, levels = sample_order)) %>% #rownames of sample table must align with colnames of txi$counts
  mutate(SampleID = RNAJBLABnumber) %>% 
  column_to_rownames(var = "RNAJBLABnumber") %>%  #rownames of sample table must align with colnames of txi$counts
  dplyr::select(cell_line, SampleID, passage, Media, OxygenConditions, Type, comment) %>% 
  mutate(CellType = case_when(Type == "HGSOC" ~ "HGSOC",
                              str_detect(Type, "normal") ~ "Normal",
                              str_detect(Type, "lymphoblastoid") ~ "Normal",
                              Type == "unknown" ~ "unknown")) %>% 
  mutate(CellType = if_else(is.na(CellType), "OtherOV", CellType))

## Centrosome/MN data
IDjoin <- read.csv("../data/cell_line_RNAseq_data/cellCharac_metaData_forRNAseqAnalysis.csv")[,c("Cell.Type","SampleID", "subtype", "subtype2")]

IF_data <- read.csv("../data/cell_line_imaging_data/CellLine_summary_data_stain1-4_V2_CMS20220523.csv") %>% 
    mutate(CA_group = if_else(X..CA.cells > median(X..CA.cells), "high", "low"), #Median split CA and MN
         MN_group = if_else(X..MN.cells > median(X..MN.cells), "high", "low")) %>% 
  left_join(IDjoin, by = "Cell.Type")  

CL_info <- meta[,c(1:6)] %>% 
  left_join(IF_data, by = "SampleID")

CL_info_s <- CL_info %>% 
  filter(SampleID != "JBLAB-25876") %>% #CIOV2 sequenced twice, only include one (p29)
  filter(!is.na(Cell.Type)) # only include samples for which both RNA and imaging data available (n=60)


## Prep deseq2 data object
CL_info_s$CA_group <- factor(CL_info_s$CA_group, levels = c("low", "high"))
CL_info_s$subtype <- factor(CL_info_s$subtype, levels = c("HGSOC", "other", "normal", "unknown"))

## Generate txi object containing selected samples only (i.e. those with complete dataset)
txi[[1]]<-txi[[1]][,CL_info_s$SampleID]
txi[[2]]<-txi[[2]][,CL_info_s$SampleID]
txi[[3]]<-txi[[3]][,CL_info_s$SampleID]

## Build DESeq2 DataSet
dds <- DESeqDataSetFromTximport(txi, CL_info_s, design = ~CA_group+subtype) #additive model

## Filter Genes
cnts <- counts(dds)
## Only keep genes that have at least 25 reads in at least 2 or more samples
keep <- apply(cnts, 1, function(x){ sum(x > 25) }) > 1
dds <- dds[keep,]

# Run DESeq2
ddsObj <- DESeq(dds, fitType = "local")

#Checking log2 distribution of counts
cnts <- counts(ddsObj, normalized=TRUE)
logCnts <- log2(cnts + 1)
#head(logCnts)
logCnts %>% 
  as.data.frame() %>%
  rownames_to_column("Gene") %>% 
  pivot_longer(names_to="Sample", values_to="Counts", starts_with("JB")) %>% 
  ggplot(aes(x=Counts)) +
  geom_density(aes(colour=Sample)) +
  guides(colour=FALSE)

###Checking (can be skipped)

modelMatrix <- model.matrix(as.formula(~ CA_group + OxygenConditions), data = CL_info_s)
#modelMatrix
countdata <- as.matrix(txi$counts)
logcounts <- log2(countdata + 1)

limma::plotMA(logcounts)
abline(h=0, col="red")

normalizedCounts <- counts(ddsObj, normalized=TRUE) 
logNormalizedCounts <- log2(normalizedCounts + 1)

limma::plotMA(logNormalizedCounts)
abline(h=0, col="red")

plotDispEsts(ddsObj)


#Generate results table
resultsNames(ddsObj)

res_CA_HighVsLow <-  DESeq2::results(ddsObj, alpha = 0.05, name =  "CA_group_high_vs_low")
res_normal_vs_HGSOC <-  DESeq2::results(ddsObj, alpha = 0.05, name =  "subtype_normal_vs_HGSOC")

sum(res_CA_HighVsLow$padj < 0.05 ,na.rm = TRUE)
sum(res_CA_HighVsLow$padj < 0.05 & res_CA_HighVsLow$log2FoldChange >0,na.rm = TRUE)
#193 genes were differentially expressed in additive model,and 103 with log2FC>0

sum(res_normal_vs_HGSOC$padj < 0.05 ,na.rm = TRUE)
sum(res_normal_vs_HGSOC$padj < 0.05 & res_normal_vs_HGSOC$log2FoldChange >0,na.rm = TRUE)
#2078 genes were differentially expressed in additive model,and 530 with log2FC>0

sig_DGE_genes <- as.data.frame(res_CA_HighVsLow) %>% 
  dplyr::filter(!is.na(padj)) %>% 
  dplyr::filter(padj < 0.05) %>% 
  arrange(padj) %>% 
  rownames_to_column(var = "gene")


### TEST if additive model (taking into account the histological subtype) is better than simple model ###
# likelihood ratio test

# create the simpler model
design.reduced <- as.formula(~ CA_group )

ddsObjC <- DESeq(ddsObj, test="LRT", reduced=design.reduced)
res_ddsObjC <- results(ddsObjC)


sum(res_ddsObjC$padj < 0.05, na.rm=TRUE)
# For only 2153 genes the more complex model fits the data better

sum(res_ddsObjC$padj < 0.05, na.rm=TRUE)/length(res_ddsObjC$padj)*100
#This equals a fraction of 11 % of genes!


# Add annotations to DESeq2 results
#Annotation File
anns <- read.table("../data/cell_line_RNAseq_data/Annotation_EnsDb.Hsapiens.v99_CMS20210420.tsv", sep = "\t", header = T)

###Add annotations
res_CA_HighVsLow_annot <- as.data.frame(res_CA_HighVsLow) %>% 
  rownames_to_column(var = "gene_id") %>% 
  left_join(anns, by = "gene_id")

res_normal_vs_HGSOC_annot <- as.data.frame(res_normal_vs_HGSOC) %>% 
  rownames_to_column(var = "gene_id") %>% 
  left_join(anns, by = "gene_id")

################
# VISUALISATION
################

#contrast CA high vs low
ddsShrink <- lfcShrink(ddsObj, coef = "CA_group_high_vs_low")

shrink_res <- as.data.frame(ddsShrink) %>%
  rownames_to_column("gene_id") %>% 
  left_join(anns, "gene_id") 

cutoff <- sort(shrink_res$pvalue)[10]
shrink_res <- shrink_res %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoff, symbol, ""))

filtTab <- shrink_res %>% 
  # dplyr::filter(!is.na(padj)) %>%
  mutate(`-log10(pvalue)` = -log10(pvalue))

topgenes <- filtTab %>% 
  filter(padj < 0.05,
         abs(log2FoldChange)>=log2(1.5)) %>% 
  arrange(desc(abs(log2FoldChange)))


filtTab %>%
  mutate(symbol = as.character(symbol)) %>% 
  filter(!is.na(padj)) %>% 
  mutate(TopGeneLabel=ifelse(gene_id %in% topgenes[1:25,]$gene_id, symbol, "")) %>% 
  ggplot(aes(x = log2FoldChange, y=`-log10(pvalue)`)) + 
  geom_point(aes(colour=padj < 0.1), size=0.25) +
  scale_colour_manual(values = c(CMSgrey, CMSlightblue), name = "FDR < 0.1") +
  geom_vline(xintercept = 0, colour = CMSdarkgrey, alpha = 0.7, linetype = "longdash") +
  # geom_vline(xintercept = c(-log2(1.5),log2(1.5)), colour = CMSdarkgrey, alpha = 0.8, linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label=TopGeneLabel), size = 1.5, position=position_jitter(),segment.size = 0.25) +
  xlim(-3,3) +
  theme_bw(7) +
  theme(legend.position = c(0.85, 0.85),
        legend.key.size = unit(0.05, "in"))
# ggsave("../figures/CAOV_paper_figure7a.pdf", height = 2.5, width = 3.5)
# write.csv(filtTab, "../source_data_figures/source_data_figure7a.csv", row.names = F)

# Gene set testing
gseaDat <- dplyr::filter(shrink_res, !is.na(entrezid), !is.na(log2FoldChange), !is.na(pvalue))

rankData <- -log10(gseaDat$pvalue) * sign(gseaDat$log2FoldChange)
names(rankData) <- gseaDat$entrezid
head(rankData)

#Load pathways
load("../data/cell_line_RNAseq_data/human_H_v5p2.rdata")

set.seed(47)

#conduct GSEA analysis
pathwaysH <- Hs.H
fgseaRes <- fgsea(pathwaysH, 
                  rankData, 
                  minSize=15, 
                  maxSize = 500, 
                  nperm=1000)

#top 10 results
fgseaRes %>% 
  arrange(desc(abs(NES))) %>%
  filter(padj < 0.05) 

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>% 
  mutate(Pathway = str_replace_all(pathway, "_", " ")) %>%
  select(-leadingEdge)

fgseaResTidy %>% 
  filter(pval < 0.05, padj < 0.1) %>%
  # top_n(20, abs(NES)) %>% 
  ggplot(aes(x = reorder(Pathway, NES), y = NES)) +
  geom_col(aes(alpha=padj), fill = CMSlightblue)  +
  geom_hline(yintercept = 0, colour = CMSgrey, size = 0.5, linetype = "dashed") +
  coord_flip() +
  scale_alpha_continuous(range = c(1,0.4), name = "FDR") +
  labs(y="Normalized Enrichment Score",
       title="% CA cells - high vs. low \n(median split)") + 
  theme_minimal(6) +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "in")) +
  guides(alpha =guide_legend(nrow=2))
# ggsave("../figures/CAOV_paper_figure7c.pdf", width = 3.5, height = 2.5)
# write.csv(fgseaResTidy, "../source_data_figures/source_data_figure7c.csv", row.names = F)

# Enrichment score plots for all significant pathways
plotting_data <- fgseaRes %>% 
  arrange(desc(abs(NES))) %>%
  filter(padj < 0.1)

#Parameters
ticksSize = 0.2
gseaParam = 1
#prep rank data
rnk <- rank(-rankData)
ord <- order(rnk)
statsAdj <- rankData[ord]
statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
statsAdj <- statsAdj/max(abs(statsAdj))

for (i in 1:nrow(plotting_data)) {
  #Define
  pw <- plotting_data[i]$pathway
  title <- str_replace_all(pw, "_", " ")
  pval <- plotting_data[i]$pval %>% round(digits = 3)
  fdr <- plotting_data[i]$padj %>% round(digits = 3)
  nes <- plotting_data[i]$NES %>% round(digits = 1)
  
  pathway <- pathwaysH[[pw]]
  
  #Prep
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  
  #Plot
  ggplot(toPlot, aes(x = x, y = y)) + 
    geom_point(color = CMSlightblue, size = 0.1) + 
    geom_hline(yintercept = 0, colour = CMSdarkgrey) + 
    geom_hline(yintercept = c(max(tops),min(bottoms)), colour = CMSorange, linetype = "dashed") + 
    geom_line(color = CMSlightblue) + 
    theme_bw(6) + 
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), 
                 size = ticksSize, colour = CMSdarkgrey) +
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank()) + 
    labs(x = "rank", y = "enrichment score") +
    ggtitle(paste0(title, "\nNES=", nes, ", p.value=", pval, ", FDR=", fdr)) 
  
  # ggsave(paste0("../figures/enrichm_plots_rnaseq/", pw, "_addMod_CA.pdf"), height = 2, width = 2.75)
}

