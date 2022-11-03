library(dplyr)
library(plyr)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(reshape2)
library(stringr)
library(ggsci)
library(ggpubr)
library(ggsignif)

plotDat <- read.table("4.DESeq2/countFile/all_count_tpm.txt", header = T)
colnames(plotDat)[7:22] <- c("D0_H8f-1A","D0_H8f-2A","D0_H8KO-1A","D0_H8KO-2A","Day_14_H8f-1A","Day_14_H8f-2A","Day_14_H8KO-1A","Day_14_H8KO-2A","Day_21_H8f-1A","Day_21_H8f-2A","Day_21_H8KO-1A","Day_21_H8KO-2A","Day_7_H8f-1A","Day_7_H8f-2A","Day_7_H8KO-1A","Day_7_H8KO-2A")

# -------- ID Trans --------
library(rtracklayer)
gff <- readGFF("/home/data/ssy033/References/genome/mm10-GRCm38.p6/gencode.vM25.annotation.gtf")
mapid <- gff[gff$type == "gene", c("gene_id", "gene_name")]
plot.df <- merge(plotDat, mapid, by.x="Geneid", by.y="gene_id")

# remove duplicate gene name(different ENSEMBL version but same name)
table(duplicated(plot.df$gene_name))
plot.df <- plot.df[!duplicated(plot.df$gene_name),]
# --------------------------
plot.df <- plot.df[,-c(1:6)]
write.table(plot.df, file = "4.DESeq2/countFile/all_tpm_idTrans_for_plot.txt")
# plot.df <- read.table("4.DESeq2/countFile/all_tpm_idTrans_for_plot.txt",header = T)

# memory-associated molecules -----
# Tcf7, Eomes, Foxo1, Zeb1, Lef1, Bach2, Bcl2, Il7r, Cd27, Cxcr3,Hdac8
dat1 <- filter(plot.df, gene_name %in% c("Tcf7", "Eomes", "Foxo1", "Zeb1", "Lef1", "Bach2", "Bcl2", "Il7r", "Cd27", "Cxcr3","Hdac8"))
dat1 <- melt(dat1, id.vars = c("gene_name"), 
             measure.vars = c("D0_H8f.1A","D0_H8f.2A","D0_H8KO.1A","D0_H8KO.2A","Day_14_H8f.1A","Day_14_H8f.2A","Day_14_H8KO.1A","Day_14_H8KO.2A","Day_21_H8f.1A","Day_21_H8f.2A","Day_21_H8KO.1A","Day_21_H8KO.2A","Day_7_H8f.1A","Day_7_H8f.2A","Day_7_H8KO.1A","Day_7_H8KO.2A"), 
             na.rm = T, variable.name = "Sample", value.name = "TPM")

dat1[which(dat1$Sample == str_extract(dat1$Sample,"^D0.+")),"TimePoint"] <- "Day_0"
dat1[which(dat1$Sample == str_extract(dat1$Sample,"^Day_14.+")),"TimePoint"] <- "Day_14"
dat1[which(dat1$Sample == str_extract(dat1$Sample,"^Day_21.+")),"TimePoint"] <- "Day_21"
dat1[which(dat1$Sample == str_extract(dat1$Sample,"^Day_7.+")),"TimePoint"] <- "Day_7"
dat1$TimePoint <- factor(dat1$TimePoint, levels = c("Day_0","Day_7","Day_14","Day_21"))

dat1[which(ifelse(str_extract(dat1$Sample,"H8f") == "H8f",TRUE,FALSE)),"Cell"] <- "H8f"
dat1[which(ifelse(str_extract(dat1$Sample,"H8KO") == "H8KO",TRUE,FALSE)),"Cell"] <- "H8KO"

dat1 <- ddply(dat1,.(gene_name,TimePoint,Cell),function(sub){data.frame(meanTPM=mean(sub$TPM))})

ggplot(dat1, aes(fill=Cell, y=meanTPM, x=TimePoint)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Memory-associated Molecules") +
  facet_wrap(~gene_name, scale="free") +
  theme_ipsum() +
  theme(legend.position="right") +
  stat_compare_means(aes(group = Cell), label = "p.format") + 
  scale_fill_manual(labels = c("H8f","H8KO"),values = c("#F1C076","#75A6F0")) + 
  xlab("")

# TE-associated molecules -----
# Tbx21, Zeb2, Gzma, Gzmb, Klrg1, Cx3cr1
dat2 <- filter(plot.df, gene_name %in% c("Tbx21", "Zeb2", "Gzma", "Gzmb", "Klrg1", "Cx3cr1"))
dat2 <- melt(dat2, id.vars = c("gene_name"), 
             measure.vars = c("D0_H8f-1A","D0_H8f-2A","D0_H8KO-1A","D0_H8KO-2A","Day_14_H8f-1A","Day_14_H8f-2A","Day_14_H8KO-1A","Day_14_H8KO-2A","Day_21_H8f-1A","Day_21_H8f-2A","Day_21_H8KO-1A","Day_21_H8KO-2A","Day_7_H8f-1A","Day_7_H8f-2A","Day_7_H8KO-1A","Day_7_H8KO-2A"), 
             na.rm = T, variable.name = "Sample", value.name = "TPM")

dat2[which(dat2$Sample == str_extract(dat2$Sample,"^D0.+")),"TimePoint"] <- "Day_0"
dat2[which(dat2$Sample == str_extract(dat2$Sample,"^Day_14.+")),"TimePoint"] <- "Day_14"
dat2[which(dat2$Sample == str_extract(dat2$Sample,"^Day_21.+")),"TimePoint"] <- "Day_21"
dat2[which(dat2$Sample == str_extract(dat2$Sample,"^Day_7.+")),"TimePoint"] <- "Day_7"
dat2$TimePoint <- factor(dat2$TimePoint, levels = c("Day_0","Day_7","Day_14","Day_21"))

dat2[which(ifelse(str_extract(dat2$Sample,"H8f") == "H8f",TRUE,FALSE)),"Cell"] <- "H8f"
dat2[which(ifelse(str_extract(dat2$Sample,"H8KO") == "H8KO",TRUE,FALSE)),"Cell"] <- "H8KO"

dat2 <- ddply(dat2,.(gene_name,TimePoint,Cell),function(sub){data.frame(meanTPM=mean(sub$TPM))})

ggplot(dat2, aes(fill=Cell, y=meanTPM, x=TimePoint)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("TE-associated Molecules") +
  facet_wrap(~gene_name, scale="free") +
  theme_ipsum() +
  theme(legend.position="right") +
  scale_fill_manual(labels = c("H8f","H8KO"),values = c("#E1A1BA","#A2BAE2")) + 
  xlab("")

# TFs and surface makers -----
# Tbx21, Gata3, Rorγ, Stat3, Foxp3, Eomes, Id3, Bcl6, Tcf7, Il7r, Klrg1
dat3 <- filter(plot.df, gene_name %in% c("Tbx21","Gata3","Rorc","Stat3", "Foxp3", "Eomes", "Id3", "Bcl6", "Tcf7", "Il7r", "Klrg1","Rorc"))
dat3 <- melt(dat3, id.vars = c("gene_name"), 
             measure.vars = c("D0_H8f-1A","D0_H8f-2A","D0_H8KO-1A","D0_H8KO-2A","Day_14_H8f-1A","Day_14_H8f-2A","Day_14_H8KO-1A","Day_14_H8KO-2A","Day_21_H8f-1A","Day_21_H8f-2A","Day_21_H8KO-1A","Day_21_H8KO-2A","Day_7_H8f-1A","Day_7_H8f-2A","Day_7_H8KO-1A","Day_7_H8KO-2A"), 
             na.rm = T, variable.name = "Sample", value.name = "TPM")

dat3[which(dat3$Sample == str_extract(dat3$Sample,"^D0.+")),"TimePoint"] <- "Day_0"
dat3[which(dat3$Sample == str_extract(dat3$Sample,"^Day_14.+")),"TimePoint"] <- "Day_14"
dat3[which(dat3$Sample == str_extract(dat3$Sample,"^Day_21.+")),"TimePoint"] <- "Day_21"
dat3[which(dat3$Sample == str_extract(dat3$Sample,"^Day_7.+")),"TimePoint"] <- "Day_7"
dat3$TimePoint <- factor(dat3$TimePoint, levels = c("Day_0","Day_7","Day_14","Day_21"))

dat3[which(ifelse(str_extract(dat3$Sample,"H8f") == "H8f",TRUE,FALSE)),"Cell"] <- "H8f"
dat3[which(ifelse(str_extract(dat3$Sample,"H8KO") == "H8KO",TRUE,FALSE)),"Cell"] <- "H8KO"

dat3 <- ddply(dat3,.(gene_name,TimePoint,Cell),function(sub){data.frame(meanTPM=mean(sub$TPM))})

ggplot(dat3, aes(fill=Cell, y=meanTPM, x=TimePoint)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("TFs and Surface Makers") +
  facet_wrap(~gene_name, scale="free") +
  theme_ipsum() +
  theme(legend.position="right") +
  scale_fill_npg(alpha = 0.7) + 
  xlab("")

# Cytokines -----
# Il2, Il7, Il15, Il21, Il23, Il17a, Il10, Ifnγ
dat4 <- filter(plot.df, gene_name %in% c("Il2", "Il7", "Il15", "Il21", "Il23a", "Il17a", "Il10", "Ifng"))
dat4 <- melt(dat4, id.vars = c("gene_name"), 
             measure.vars = c("D0_H8f.1A","D0_H8f.2A","D0_H8KO.1A","D0_H8KO.2A","Day_14_H8f.1A","Day_14_H8f.2A","Day_14_H8KO.1A","Day_14_H8KO.2A","Day_21_H8f.1A","Day_21_H8f.2A","Day_21_H8KO.1A","Day_21_H8KO.2A","Day_7_H8f.1A","Day_7_H8f.2A","Day_7_H8KO.1A","Day_7_H8KO.2A"), 
             na.rm = T, variable.name = "Sample", value.name = "TPM")

dat4[which(dat4$Sample == str_extract(dat4$Sample,"^D0.+")),"TimePoint"] <- "Day_0"
dat4[which(dat4$Sample == str_extract(dat4$Sample,"^Day_14.+")),"TimePoint"] <- "Day_14"
dat4[which(dat4$Sample == str_extract(dat4$Sample,"^Day_21.+")),"TimePoint"] <- "Day_21"
dat4[which(dat4$Sample == str_extract(dat4$Sample,"^Day_7.+")),"TimePoint"] <- "Day_7"
dat4$TimePoint <- factor(dat4$TimePoint, levels = c("Day_0","Day_7","Day_14","Day_21"))

dat4[which(ifelse(str_extract(dat4$Sample,"H8f") == "H8f",TRUE,FALSE)),"Cell"] <- "H8f"
dat4[which(ifelse(str_extract(dat4$Sample,"H8KO") == "H8KO",TRUE,FALSE)),"Cell"] <- "H8KO"

dat4 <- ddply(dat4,.(gene_name,TimePoint,Cell),function(sub){data.frame(meanTPM=mean(sub$TPM))})

ggplot(dat4, aes(fill=Cell, y=meanTPM, x=TimePoint)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Cytokines") +
  facet_wrap(~gene_name, scale="free") +
  stat_compare_means(aes(group = Cell), label = "p.signif") + 
  theme_ipsum() +
  theme(legend.position="right") +
  scale_fill_jama(alpha = 0.7) + 
  xlab("")

# Jing's validate
# "Il2","Rorc","Il17a","Il15","Bcl6","Gata3","Id3","Ifng","Il7r","Klrg1","Tbx21","Il23a","Il7","Hdac8"
dat5 <- filter(plot.df, gene_name %in% c("Il2","Rorc","Il17a","Il15","Bcl6","Gata3","Id3","Ifng","Il7r","Klrg1","Tbx21","Il23a","Il7","Hdac8"))
dat5 <- melt(dat5, id.vars = c("gene_name"), 
             measure.vars = c("D0_H8f.1A","D0_H8f.2A","D0_H8KO.1A","D0_H8KO.2A","Day_14_H8f.1A","Day_14_H8f.2A","Day_14_H8KO.1A","Day_14_H8KO.2A","Day_21_H8f.1A","Day_21_H8f.2A","Day_21_H8KO.1A","Day_21_H8KO.2A","Day_7_H8f.1A","Day_7_H8f.2A","Day_7_H8KO.1A","Day_7_H8KO.2A"), 
             na.rm = T, variable.name = "Sample", value.name = "TPM")

dat5[which(dat5$Sample == str_extract(dat5$Sample,"^D0.+")),"TimePoint"] <- "Day_0"
dat5[which(dat5$Sample == str_extract(dat5$Sample,"^Day_14.+")),"TimePoint"] <- "Day_14"
dat5[which(dat5$Sample == str_extract(dat5$Sample,"^Day_21.+")),"TimePoint"] <- "Day_21"
dat5[which(dat5$Sample == str_extract(dat5$Sample,"^Day_7.+")),"TimePoint"] <- "Day_7"
dat5$TimePoint <- factor(dat5$TimePoint, levels = c("Day_0","Day_7","Day_14","Day_21"))

dat5[which(ifelse(str_extract(dat5$Sample,"H8f") == "H8f",TRUE,FALSE)),"Cell"] <- "H8f"
dat5[which(ifelse(str_extract(dat5$Sample,"H8KO") == "H8KO",TRUE,FALSE)),"Cell"] <- "H8KO"

# dat5$Cell <- factor(dat5$Cell)
dat5$gene_name <- factor(dat5$gene_name, levels = c("Il2","Rorc","Il17a","Il15","Bcl6","Gata3","Id3","Ifng","Il7r","Klrg1","Tbx21","Il23a","Il7","Hdac8"))

dat5.mean <- ddply(dat5,.(gene_name,TimePoint,Cell),function(sub){data.frame(meanTPM=mean(sub$TPM))})

ggplot(dat5, aes(fill=Cell, y=TPM, x=TimePoint)) + 
  geom_bar(position="dodge", stat="summary", fun = "mean") + 
  stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black",
               width = 0.15,position = position_dodge(1)) + 
  stat_compare_means(aes(label = ..p.format..), method = "t.test", method.args = list(alternative = "two.sided"), size = 3) + 
  ggtitle("Validated Genes") +
  facet_wrap(~gene_name, scale="free") +
  theme_ipsum() +
  theme(legend.position="right") +
  scale_fill_jama(alpha = 0.7) + 
  xlab("") 
ggsave("3.Count/validate_gene_check.png", width = 14, height = 12)
# stat_compare_means(aes(label = ..p.signif..), method = "t.test", method.args = list(alternative = "two.sided")) +
