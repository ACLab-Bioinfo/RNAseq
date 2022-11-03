setwd("~/Project/AlfredCHENG/HDAC8/")

library(stringr)

# HDAC8
# "ENSMUST00000087916.10","ENSMUST00000154872.7","ENSMUST00000137785.7","ENSMUST00000148481.7","ENSMUST00000154024.1"

nameList <- list.files("3.Count/rsem_out/",pattern = "*.isoforms.results")
sampName <- str_split_fixed(nameList,"\\.",n = 2)[,1]

rsem_list <- list()
for (i in 1:length(nameList)){
  dfnew <- read.table(paste0("3.Count/rsem_out/",nameList[i]), header = T)[,c(1,6)]
  colnames(dfnew) <- c("transcript_id",paste0(sampName[i],"_TPM"))
  rsem_list[[i]] <- dfnew 
}

rsem_all <- Reduce(function(x, y) merge(x, y, by = "transcript_id", allow.cartesian=TRUE),rsem_list)

Hdac8.dat <- filter(rsem_all, transcript_id %in% c("ENSMUST00000087916.10","ENSMUST00000154872.7","ENSMUST00000137785.7","ENSMUST00000148481.7","ENSMUST00000154024.1"))
Hdac8.dat.plot <- melt(Hdac8.dat, id.vars = c("transcript_id"), 
             measure.vars = paste0(sampName,"_TPM"), 
             na.rm = T, variable.name = "Sample", value.name = "TPM")
Hdac8.dat.plot[which(ifelse(str_extract(Hdac8.dat.plot$Sample,"H8f") == "H8f",TRUE,FALSE)),"Cell"] <- "H8f"
Hdac8.dat.plot[which(ifelse(str_extract(Hdac8.dat.plot$Sample,"H8KO") == "H8KO",TRUE,FALSE)),"Cell"] <- "H8KO"
Hdac8.dat.plot[which(Hdac8.dat.plot$Sample == str_extract(Hdac8.dat.plot$Sample,"^D0.+")),"TimePoint"] <- "Day_0"
Hdac8.dat.plot[which(Hdac8.dat.plot$Sample == str_extract(Hdac8.dat.plot$Sample,"^Day_14.+")),"TimePoint"] <- "Day_14"
Hdac8.dat.plot[which(Hdac8.dat.plot$Sample == str_extract(Hdac8.dat.plot$Sample,"^Day_21.+")),"TimePoint"] <- "Day_21"
Hdac8.dat.plot[which(Hdac8.dat.plot$Sample == str_extract(Hdac8.dat.plot$Sample,"^Day_7.+")),"TimePoint"] <- "Day_7"
Hdac8.dat.plot$TimePoint <- factor(Hdac8.dat.plot$TimePoint, levels = c("Day_0","Day_7","Day_14","Day_21"))

ggplot(Hdac8.dat.plot[which(Hdac8.dat.plot$TimePoint =="Day_21"),], aes(fill=Sample, y=TPM, x=transcript_id)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Hdac8 Transcript--Day_21") +
  theme_ipsum() +
  theme(legend.position="right") +
  scale_fill_simpsons() +
  xlab("")
ggsave("3.Count/Hdac8.trans.check/Hdac8-trans-Day21.png", width = 12, height = 5)


