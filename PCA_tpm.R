library(ggpubr)
library(ggthemes)
library(gmodels)
library(export)

#使用prcomp函数进行PCA分析
# 数据转置，转换成行为样本，列为基因的矩阵
tpmDat <- read.csv("3.Count/tpmDat.csv", row.names = "Geneid")[,-1]
tpm.t <- as.data.frame(t(tpmDat))
data.pca <- prcomp(tpm.t)
# 查看PCA的分析结果
summary(data.pca)
# 绘制主成分的碎石图
screeplot(data.pca, npcs = 10, type = "lines")

#查看主成分的结果
pca.scores <- as.data.frame(data.pca$x)
head(pca.scores)

#用gmodels里面自带的fast.prcomp对PC进行计算。
pca.info <- fast.prcomp(tpmDat)
#查看PCA的结果。
head(pca.info)
summary(pca.info) #能看到每个PC具体的解释比例
head(pca.info$rotation) #特征向量，回归系数
head(pca.info$sdev) #特征值的开方
head(pca.info$x) #样本得分score

#将PCA的结果和数据的表型构建为一个数据矩阵。
##这个里面的Type=c(rep("control",50),rep("case",374))和列表里面样本的名字是对应的，根据自己数据的情况灵活选择。
pca.data <- data.frame(sample = rownames(pca.info$rotation), 
                       Type=factor(c(rep("D7_H8f",3),rep("D7_H8KO",3),rep("D14_H8f",3),rep("D14_H8KO",3),rep("D21_H8f",3),rep("D21_H8KO",3)), 
                                   levels = c("D7_H8f","D7_H8KO","D14_H8f","D14_H8KO","D21_H8f","D21_H8KO"),
                                   ordered = TRUE),
                       pca.info$rotation)

cairo_pdf(file = "3.Count/pca_tpm.pdf", width = 5, height = 5)
ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=TRUE, ellipse.type="confidence",main="PCA plot", legend = "right",
          xlab = "PC1 - 97.38%", ylab = "PC2 - 2.10%")
dev.off()

