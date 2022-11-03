library(clusterProfiler)
library(ggnewscale)
library(org.Mm.eg.db)

degs <- read.csv("4.DESeq2/DEGs/Day_7-H8f_H8KO.csv",header = T)
diff_gene_deseq <- subset(degs, pvalue < 0.05 & abs(log2FoldChange) > 0.5)
dim(diff_gene_deseq)

gene.GO <- enrichGO(gene = diff_gene_deseq$X,
                    OrgDb = org.Mm.eg.db,
                    keyType = "SYMBOL",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.5,
                    readable = T)

write.table(gene.GO@result, file = "src/resFile/geneGO_cell.txt",sep = "\t", row.names = T,col.names = T)

gene.KEGG <- enrichKEGG(gene = D,
                        organism = "hsa",
                        keyType = "ncbi-geneid",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.01,
                        qvalueCutoff = 0.05 )
write.table(gene.KEGG@result, file = "src/resFile/geneKEGG.txt",sep = "\t", row.names = T,col.names = T)

# plotting -----
library(enrichplot)
p1 <- barplot(gene.GO, showCategory = 15)
p2 <- dotplot(gene.GO, showCategory = 15)
plot_grid(p1,p2,ncol = 2)

p3 <- barplot(gene.KEGG, showCategory = 20)
p4 <- dotplot(gene.KEGG, showCategory = 20)
plot_grid(p3,p4,ncol = 2)

cnetplot(gene.KEGG, showCategory = 5)
cnetplot(gene.GO, showCategory = 10)

pw <- pairwise_termsim(gene.GO)
emapplot(pw, showCategory = 20) 

pw <- pairwise_termsim(gene.GO)
emapplot(pw, showCategory = 100) 

goplot(gene.GO)
