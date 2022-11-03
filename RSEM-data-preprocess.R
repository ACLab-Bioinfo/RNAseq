setwd("~/Project/AlfredCHENG/HDAC8/")

library(stringr)

nameList <- list.files("3.Count/rsem_out/",pattern = "*.genes.results")
sampName <- str_split_fixed(nameList,"\\.",n = 2)[,1]

rsem_list <- list()
for (i in 1:length(nameList)){
  dfnew <- read.table(paste0("3.Count/rsem_out/",nameList[i]), header = T)[,c(1,6)]
  colnames(dfnew) <- c("gene_id",paste0(sampName[i],"_TPM"))
  rsem_list[[i]] <- dfnew 
}

rsem_all <- Reduce(function(x, y) merge(x, y, by = "gene_id", allow.cartesian=TRUE),rsem_list)

# -------- ID Trans --------
library(rtracklayer)
gff <- readGFF("/home/data/ssy033/References/genome/mm10-GRCm38.p6/gencode.vM25.annotation.gtf")
mapid <- gff[gff$type == "gene", c("gene_id", "gene_name")]
df <- merge(rsem_all, mapid, by="gene_id")

# remove duplicate gene name(different ENSEMBL version but same name)
table(duplicated(df$gene_name))
df <- df[!duplicated(df$gene_name),]

dat <- df
write.table(dat, file = "4.DESeq2/countFile/all_RSEM_TPM_idTrans.txt")
# --------------------------


