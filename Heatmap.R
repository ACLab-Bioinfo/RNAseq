setwd("~/Project/AlfredCHENG/HDAC8/")

dds <- DESeqDataSetFromMatrix(countData = deseqDat,
                              colData = meta_info, design = ~Group)
dds <- DESeq(dds)
dds

# 3. Extract DESeq normalized data -----
rld <- rlogTransformation(dds)
vsd <- vst(dds)
exprSet_vsd <- assay(vsd)
write.table(exprSet_vsd, file="4.DESeq2/countFile/all_DESeq2.vsd.normalization.txt", sep="\t",quote=F)

head(exprSet_vsd)

expr.corr <- cor(exprSet_vsd)

library(pheatmap)

pheatmap(expr.corr,cluster_cols = T,
              show_rownames = T,
              show_colnames = T,
              scale = "none", 
              fontsize_row = 8,
              fontsize_col = 8,
              border = FALSE,
              color = colorRampPalette(colors = c("blue","white","red"))(100))
