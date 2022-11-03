library(DESeq2)
setwd("~/Project/AlfredCHENG/HDAC8/")

# 1.Prep data -----
meta_info <- read.csv("4.DESeq2/metaInfo2.csv")
meta_info$Group <- factor(meta_info$Group, levels = c("H8f","H8KO"))
meta_info$TimePoint <- factor(meta_info$TimePoint, levels = c("Day_0", "Day_7", "Day_14","Day_21"))
str(meta_info)

dat <- read.table("4.DESeq2/countFile/all_count_idTrans.txt",header = T)
deseqDat <- dat[,-c(1:6)]
rownames(deseqDat) <- deseqDat$gene_name
deseqDat <- deseqDat[,-17]

# 2. Perform DESeq -----
dds <- DESeqDataSetFromMatrix(countData = deseqDat,
                              colData = meta_info, design = ~Group)
dds <- DESeq(dds)
dds

# remoce D14-H8f1A
dds <- DESeqDataSetFromMatrix(countData = select(deseqDat,-Day_14_H8f.1A),
                              colData = meta_info[-which(meta_info$Sample == "Day_14_H8f.1A"),], design = ~Group)
dds <- DESeq(dds)
dds

# 3. Extract DESeq normalized data -----
rld <- rlogTransformation(dds)
vsd <- vst(dds)
exprSet_vsd <- assay(vsd)
write.table(exprSet_vsd, file="4.DESeq2/countFile/all_DESeq2.vsd.normalization.rm.D14h8f.1A.txt", sep="\t",quote=F)

head(exprSet_vsd)

pcaData <- plotPCA(vsd, intgroup=c('Group'), returnData=TRUE)
pcaData$simpName <- c("Day0-1","Day0-2","Day0-1","Day0-2","Day_14-1","Day_14-2","Day_14-1","Day_14-2","Day_21-1","Day_21-2","Day_21-1","Day_21-2","Day_7-1","Day_7-2","Day_7-1","Day_7-2")
pcaData$TimePoint <- c(rep("Day0",4),rep("Day14",4),rep("Day21",4),rep("Day7",4))

# 4. Plot -----
# A basic scatterplot with color depending on Species
ggplot(pcaData, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point(size=6) +
  ggtitle("Removed D14-H8f-1A") +
  theme_ipsum() + 
  geom_text(
    label=pcaData$simpName, 
    size = 2.5,
    nudge_x = 0.25, nudge_y = 0.75
  )

# Time-Point PCA
dds <- DESeqDataSetFromMatrix(countData = select(deseqDat, starts_with("Day_7")), 
                              colData = meta_info[which(meta_info$TimePoint == "Day_7"),], design = ~Group)
dds <- DESeq(dds)
rld <- rlogTransformation(dds)
vsd <- vst(dds)
exprSet_vsd <- assay(vsd)
pcaData <- plotPCA(vsd, intgroup=c('Group'), returnData=TRUE)
pcaData$simpName <- c("Day7-1","Day7-2","Day7-1","Day7-2")

Day7 <- ggplot(pcaData, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point(size=6) +
  ggtitle("Day_7") +
  theme_ipsum() + 
  stat_ellipse(level = 0.9) + 
  geom_text(
    label=pcaData$simpName, 
    size = 2.5,
    nudge_x = 0.25, nudge_y = 0.75
  ) 
ggsave("4.DESeq2/PCA.plot/Day7.png",Day7, width = 5, height = 6, units = "in")

# remove 14-2 and Day_21 PCA
dds <- DESeqDataSetFromMatrix(countData = select(deseqDat, one_of("D0_H8f.1A","D0_H8f.2A","D0_H8KO.1A","D0_H8KO.2A","Day_14_H8f.1A","Day_14_H8KO.1A",
                                                                  "Day_7_H8f.1A","Day_7_H8f.2A","Day_7_H8KO.1A","Day_7_H8KO.2A")), 
                              colData = meta_info[which(meta_info$Sample %in% c("D0_H8f.1A","D0_H8f.2A","D0_H8KO.1A","D0_H8KO.2A","Day_14_H8f.1A","Day_14_H8KO.1A",
                                                                                "Day_7_H8f.1A","Day_7_H8f.2A","Day_7_H8KO.1A","Day_7_H8KO.2A")),], design = ~Group)
dds <- DESeq(dds)
rld <- rlogTransformation(dds)
vsd <- vst(dds)
exprSet_vsd <- assay(vsd)
pcaData <- plotPCA(vsd, intgroup=c('Group'), returnData=TRUE)
pcaData$simpName <- c("Day0-1","Day0-2","Day0-1","Day0-2","Day_14-1","Day_14-1","Day_7-1","Day_7-2","Day_7-1","Day_7-2")

ggplot(pcaData, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point(size=6) +
  ggtitle("Remove 14-2 and Day_21 samples") +
  theme_ipsum() + 
  geom_text(
    label=pcaData$simpName, 
    size = 2.5,
    nudge_x = 0.25, nudge_y = 0.75
  ) 
ggsave("4.DESeq2/PCA.plot/rm14-2-Day21.png", width = 5, height = 6, units = "in")
