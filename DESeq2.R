library(DESeq2)
setwd("~/Project/AlfredCHENG/HDAC8/")

# 1.Prep data -----
meta_info <- read.csv("4.DESeq2/metaInfo.csv")
meta_info$Group <- factor(meta_info$Group, levels = c("H8f","H8KO"))
meta_info$TimePoint <- factor(meta_info$TimePoint, levels = c("Day_0", "Day_7", "Day_14","Day_21"))
str(meta_info)

dat <- read.table("4.DESeq2/countFile/all_count_idTrans.txt",header = T)
deseqDat <- dat[,-c(1:6)]
rownames(deseqDat) <- deseqDat$gene_name
deseqDat <- deseqDat[,-17]

# 2.Make deseq and run -----
## 2.1 Day_0 H8f&H8KO
dds <- DESeqDataSetFromMatrix(countData = select(deseqDat, starts_with("D0")), 
                              colData = meta_info[which(meta_info$TimePoint == "Day_0"),], design = ~Group)
dds <- DESeq(dds)
res <- results(dds, name = "Group_H8KO_vs_H8f")
res <- res[order(res$pvalue), ] 
summary(res)
table(res$padj < 0.05)
diff_gene_deseq <- subset(res, padj < 0.1 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq)
write.csv(diff_gene_deseq, file = "4.DESeq2/DEGs/Day_0-H8f_H8KO-0.1_1.csv")

## 2.2 Day_7 H8f&H8KO
dds <- DESeqDataSetFromMatrix(countData = select(deseqDat, starts_with("Day_7")), 
                              colData = meta_info[which(meta_info$TimePoint == "Day_7"),], design = ~Group)
dds <- DESeq(dds)
res <- results(dds, name = "Group_H8KO_vs_H8f")
res <- res[order(res$pvalue), ] 
summary(res)
table(res$padj < 0.1)
diff_gene_deseq <- subset(res, padj < 0.1 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq)
write.csv(diff_gene_deseq, file = "4.DESeq2/DEGs/Day_7-H8f_H8KO-0.1_1.csv")

## 2.3 Day_14 H8f&H8KO
dds <- DESeqDataSetFromMatrix(countData = select(deseqDat, starts_with("Day_14")), 
                              colData = meta_info[which(meta_info$TimePoint == "Day_14"),], design = ~Group)
dds <- DESeq(dds)
res <- results(dds, name = "Group_H8KO_vs_H8f")
res <- res[order(res$pvalue), ] 
summary(res)
table(res$padj < 0.1)
diff_gene_deseq <- subset(res, padj < 0.1 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq)
write.csv(diff_gene_deseq, file = "4.DESeq2/DEGs/Day_14-H8f_H8KO-0.1_1.csv")

## 2.4 Day_21 H8f&H8KO
dds <- DESeqDataSetFromMatrix(countData = select(deseqDat, starts_with("Day_21")), 
                              colData = meta_info[which(meta_info$TimePoint == "Day_21"),], design = ~Group)
dds <- DESeq(dds)
res <- results(dds, name = "Group_H8KO_vs_H8f")
res <- res[order(res$pvalue), ] 
summary(res)
table(res$padj < 0.1)
diff_gene_deseq <- subset(res, padj < 0.1 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq)
write.csv(diff_gene_deseq, file = "4.DESeq2/DEGs/Day_21-H8f_H8KO-0.1_1.csv")


# 3.Extract results -----
resultsNames(dds)
## "Group_H8KO_vs_H8f"

setwd("~/Project/AlfredCHENG/HDAC8/4.DESeq2/countFile/")
dat <- read.table("all_count.txt",header = T)
dim(dat)

# 4. remove Day14 bad sample, Day_14 H8f&H8KO
dds <- DESeqDataSetFromMatrix(countData = select(deseqDat, one_of("Day_14_H8f.1A","Day_14_H8KO.1A","Day_14_H8KO.2A")), 
                              colData = meta_info[which(meta_info$Sample %in% c("Day_14_H8f.1A","Day_14_H8KO.1A","Day_14_H8KO.2A")),], design = ~Group)
dds <- DESeq(dds)
res <- results(dds, name = "Group_H8KO_vs_H8f")
res <- res[order(res$pvalue), ] 
summary(res)
table(res$padj < 0.1)
diff_gene_deseq <- subset(res, padj < 0.1 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq)
write.csv(res, file = "4.DESeq2/DEGs/withD14-H8f-1A_Day_14-H8f_H8KO.csv")

