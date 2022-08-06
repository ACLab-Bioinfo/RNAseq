# This is code blocks for RNA-seq analysis
# workflow: QC(fastp)--Mapping(STAR/HISAT2)--Count(HTSeq)

# Reference download
# 下载的小鼠基因组
cd ~/reference
mkdir -p  genome/mm10  && cd genome/mm10
nohup wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz  &
tar zvfx chromFa.tar.gz
cat *.fa > mm10.fa
rm chr*.fa

# 下载hg19：
cd ~/reference
mkdir -p genome/hg19  && cd genome/hg19
nohup wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz &
tar zvfx chromFa.tar.gz
cat *.fa > hg19.fa
rm chr*.fa

# 下载hg38
cd ~/reference
mkdir -p genome/hg38  && cd genome/hg38
nohup wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz  &

# bowtie软件建立索引文件
cd ~/reference
mkdir -p index/bowtie && cd index/bowtie
nohup time ~/biosoft/bowtie/bowtie2-2.2.9/bowtie2-build  ~/reference/genome/hg19/hg19.fa  ~/reference/index/bowtie/hg19 1>hg19.bowtie_index.log 2>&1 &
nohup time ~/biosoft/bowtie/bowtie2-2.2.9/bowtie2-build  ~/reference/genome/hg38/hg38.fa  ~/reference/index/bowtie/hg38 1>hg38.bowtie_index.log 2>&1 &
nohup time ~/biosoft/bowtie/bowtie2-2.2.9/bowtie2-build  ~/reference/genome/mm10/mm10.fa  ~/reference/index/bowtie/mm10 1>mm10.bowtie_index.log 2>&1 &

# bwa软件建立索引文件
cd ~/reference
mkdir -p index/bwa && cd index/bwa
nohup time ~/biosoft/bwa/bwa-0.7.15/bwa index   -a bwtsw   -p ~/reference/index/bwa/hg19  ~/reference/genome/hg19/hg19.fa 1>hg19.bwa_index.log 2>&1   &
nohup time ~/biosoft/bwa/bwa-0.7.15/bwa index   -a bwtsw   -p ~/reference/index/bwa/hg38  ~/reference/genome/hg38/hg38.fa 1>hg38.bwa_index.log 2>&1   &
nohup time ~/biosoft/bwa/bwa-0.7.15/bwa index   -a bwtsw   -p ~/reference/index/bwa/mm10  ~/reference/genome/mm10/mm10.fa 1>mm10.bwa_index.log 2>&1   &

# hisat软件建立索引文件
cd ~/reference
mkdir -p index/hisat && cd index/hisat
nohup wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg19.tar.gz  &
nohup wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg38.tar.gz  &
nohup wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38.tar.gz &
tar zxvf hg19.tar.gz
tar zxvf grcm38.tar.gz
tar zxvf hg38.tar.gz

#######################



## -------------------------- QC --------------------------
# fastqc 
fastqc -t 24 -o ./FastQC_result.fix/ -q ./fix.fastq/test_*.gz &

# fastp
# 参数:
# -z 输出压缩格式 -w 线程数
# -q 设置低质量的标准，设置Q20
# -f, –trim_front1 裁剪read1前多少个碱基，默认0
# -F, –trim_front2 裁剪read2前多少个碱基，默认0，如果没有指定，将保持与read1相同设置
# -t, –trim_tail1 裁剪read1末尾多少个碱基，默认0
# -T, –trim_tail2 裁剪read2末尾多少个碱基，默认0，如果没有指定，将保持与read1相同设置

for i in 1 2 3; do
    {
    fastp -i rawdata-trim/c${i}_1_s.fastq -o cleandata-trim/c${i}_1_s.fastq.gz \
        -I rawdata-trim/c${i}_2_s.fastq -O cleandata-trim/c${i}_2_s.fastq.gz \
        --thread=3 \
        > log/fastp.log/c${i}.fastp.log 2>&1 &
    }
done

for i in 1 2 3; do
    {
    fastp -i rawdata-trim/t${i}_1_s.fastq -o cleandata-trim/t${i}_1_s.fastq.gz \
        -I rawdata-trim/t${i}_2_s.fastq -O cleandata-trim/t${i}_2_s.fastq.gz \
        --thread=3 \
        > log/fastp.log/t${i}.fastp.log 2>&1 &
    }
done

# cutadapt
cutadapt -j 6 --times 1  -e 0.1  -O 3  --quality-cutoff 25  -m 55 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-o fix.fastq/test_R1_cutadapt.temp.fq.gz \
-p fix.fastq/test_R2_cutadapt.temp.fq.gz \
raw.fastq/test_R1.fq.gz \
raw.fastq/test_R2.fq.gz > fix.fastq/test_cutadapt.temp.log 2>&1 &

# trim_galore
# –length 去除长度小于参数值的reads
# -e 允许的最大误差
# –stringency 设置与接头重叠的序列
# –paired 双端文件
trim_galore -q 6 --phred33 \
--length 35 -e 0.1 --stringency 4 --paired -o clean/ \
fqFile/CD4_Rep1_1.fastq.gz fqFile/CD4_Rep1_2.fastq.gz \
> log/CD4_Rep1.trim.log 2>&1 &
# 批处理
ls fqFile/*_1.fastq.gz > fqFile/1 
ls fqFile/*_2.fastq.gz > fqFile/2
paste raw/srr.list fqFile/1 fqFile/2 > fqFile/config.raw

cat fqFile/config.raw | while read id; \
do echo $id
arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
nohup time trim_galore -q 6 --phred33 --length 35 -e 0.1 --stringency 4 \
--paired -o clean/ $fq1 $fq2 > log/trim.log/${sample}_trim.log 2>&1 & \
done

# MultiQC v1.9  
# -i 输出的前缀
multiqc tmp.res/qcReport/cleanqc -o tmp.res/qcReport -i clean

## --------------------------------------------------------



############################################################
# hisat2 build index and mapping
############################################################
# download some resource
## SNP
# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/

## GTF

# hisat2-2.1.0
# -x ：参考基因组索引的basename，即前缀名
# -1 ：双端测序的read1 list
# -2 ：双端测序的read2 list
# -S ：SAM写入的文件名，默认写入到标准输出中

# make exon 
hisat2_extract_exons.py hg38_refseq.gtf > hg38_refseq.exon &

# make splice site
hisat2_extract_splice_sites.py hg38_refseq.gtf > hg38_refseq.ss &

# make snp and haplotype
hisat2_extract_snps_haplotypes_UCSC.py ref_hg38.fa snp151Common.txt snp151Common &

# build index
hisat2-build -p 6 --snp snp151Common.snp --haplotype snp151Common.haplotype --exon hg38_refseq.exon  --ss hg38_refseq.ss ref_hg38.fa ref_hg38.fa.snp_gtf > hisat2_build.log 2>&1 & 

# mapping
hisat2 -p 6 \
-x /Users/meng/ngs_course/reference/hisat2_index/ref_hg38.fa.snp_gtf \
-1 ./fix.fastq/test_R1_cutadapt.fq.gz \
-2 ./fix.fastq/test_R2_cutadapt.fq.gz \
-S ./bam/test_hisat2.sam > ./bam/test_hisat2.log 2>&1 &

# 批处理
cd clean
ls *_1_clean.fq.gz > 1  
ls *_2_clean.fq.gz > 2
 
paste ../raw/fq.list 1 2 > config.clean

cat clean/config.clean | while read id; \
do
echo $id
arr=(${id})
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
hisat2 -p 6 \
-x /public/workspace/lincs/lab7/ref/hisat2/genome_tran \
-1 clean/${fq1} \
-2 clean/${fq2} \
-S HISAT2_aligned/sam/${sample}_hisat2.sam > log/hisat2.log/${sample}_hisat2.log 2>&1 & \
done


############################################################
# STAR build index and mapping
############################################################
#  STAR-2.7.1a
# –runMode mapping 如果需要建index，参数为genomeGenerate
# –genomeDir index所在路径
# –readFilesCommand 如果read是gz压缩的，要写zcat
# –outFileNamePrefix 输出的文件夹
# –outSAMtype 输出文件类型
# –outFilterIntronMotifs RemoveNoncanonical 删掉非典型的

time STAR \
--runMode alignReads \
--genomeDir  /public/workspace/lincs/lab7/ref/star_271a/ \
--runThreadN 6 \
--readFilesIn clean/c1_1_clean.fq.gz clean/c1_2_clean.fq.gz \
--readFilesCommand zcat \
--outFileNamePrefix bam/c1_STAR \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--outSAMattributes All \
--outFilterIntronMotifs RemoveNoncanonical > bam/c1_STAR.log 2>&1 &

# 批处理
cd clean
ls *_1_clean.fq.gz > 1  
ls *_2_clean.fq.gz > 2
 
paste ../raw/fq.list 1 2 > config.clean

cat clean/config.clean | while read id; \
do
echo $id
arr=(${id})
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
time STAR \
--runMode alignReads \
--genomeDir  /public/workspace/lincs/lab7/ref/star_271a/ \
--runThreadN 4 \
--readFilesIn ./clean/${fq1} ./clean/${fq2} \
--readFilesCommand zcat \
--outFileNamePrefix ./bam \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--outSAMattributes All \
--outFilterIntronMotifs RemoveNoncanonical > bam/${sample}_STAR.log 2>&1 \
done


############################################################
# HTSeq
############################################################
mamba install htseq
mamba install subread

# HTSeq
# -f 输入的格式 -r 文件是按坐标排序还是名字排序(pos,name)
# –stranded 分链测序还是不分链（reads无论比对到正链还是反链都计数)
# –minaqual 默认10代表1%的错误，一般20/30
# –type count exon –idattr 输出后以什么为行名
# –nonunique 如果不是unique mapping （一条reads比对到多个），none都不选
# –secondary-alignments 次优的不考虑
# –supplementary-alignments 补充的比对不考虑
htseq-count -f bam -r pos \
--max-reads-in-buffer 1000000 \
--stranded no \
--minaqual 10 \
--type exon --idattr gene_id \
--mode union --nonunique none --secondary-alignments ignore --supplementary-alignments ignore \
--counts_output ./count_result/test_count.tsv \
--nprocesses 1 \
./bam/test_hisat2.sort.bam  ./reference/gtf/hg38_refseq_from_ucsc.rm_XM_XR.fix_name.gtf > ./count_result/test_count.HTSeq.log  2>&1 & 


# HTSeq批处理流程
# 用*替代
nohup time htseq-count -f bam -r pos \
--max-reads-in-buffer 1000000 \
--stranded no \
--minaqual 10 \
--type exon --idattr gene_id \
--mode intersection-nonempty --nonunique none --secondary-alignments ignore \
--supplementary-alignments ignore \
STAR_aligned/coordsrt/*_STARAligned_coordsrt.out.bam \
/public/workspace/lincs/lab7/ref/Homo_sapiens.GRCh38.98.gtf \
> STAR_HTSeq/*_STARAligned_HTSeq_counts.txt &

############################################################
# feature count
############################################################
# subread-featureCounts
# -a 输入GTF/GFF基因组注释文件(annotation file)
# -F 指定区间注释文件的格式，默认是GTF
# -o 输出文件:可输出raw counts的txt文本及raw counts的summary文本
# -p 针对paired-end数据
# -g 指定统计的meta-feature名称，默认是gene_id
# -t 指定统计的feature名称， 默认是exon
# -T 指定线程数

# 对所有bam一起处理，自动merge成counts文件
featureCounts -T 2 -p \
-t exon -g gene_id \
-a /public/workspace/lincs/lab7/ref/Homo_sapiens.GRCh38.98.gtf \
-o HISAT2_featureCounts/all.id.txt \
HISAT2_aligned/hisat2.bam/*_hisat2.bam > ./count_result/test_count.featureCounts.log  2>&1 & 

# 使用multiqc对featureCounts统计结果进行可视化
multiqc all.id.txt.summary


############################################################
# Salmon
############################################################

# index
# salmon v1.4.0
# -i arg salmon index(索引)
# -t [–transcripts] arg 转录本fasta 文件
salmon index -t /public/workspace/lincs/lab7/ref/Homo_sapiens.GRCh38.cdna.all.fa.gz \
-i salmon.res/GRCh38_salmon_index

# quant
# 参数:
# -l #arg 字符串类型描述文库类型
# -i #arg salmon index(索引)
# -r #arg 文件中包括不匹配的序列（如：单端测序序列）
# -1 #arg 文件中包含#1匹配
# -2 #arg 文件中包含#2匹配

time salmon quant -i salmon.res/GRCh38_salmon_index \
-l A -p 6 \
-1 clean/t1_1_clean.fq.gz \
-2 clean/t1_2_clean.fq.gz \
-o salmon.res/t1_salmon > log/salmon.log/t1_salmon.log 2>&1 &

# 生成文件（以样本t1为例） salmon.res/t1_salmon
# aux_info 辅助文件夹，内含多个文件
# aux_info/fld.gz：在辅助文件夹中，该文件记录的是观察到的片段长度分布的近似值
# cmd_info.json 记录salmon程序运行的命令和参数
# lib_format_counts.json 记录有关文库格式和reads比对的情况
# libParams
# logs
# quant.sf

# quant.sf文件有5列，分别是Name，Length，EffectiveLength，TPM和NumReads。分别表示的含义如下所述：
# Name — target transcript 名称， 由输入的 transcript database (FASTA file)所提供。 Length — target transcript 长度，即有多少个核苷酸.
# EffectiveLength — target transcript 计算的有效长度。
# TPM — 估计转录本的表达量。
# NumReads — 估计比对到每个转录本的reads数（期望值）。

# Salmon批处理流程
cat clean/config.clean | while read id; \
do
echo $id
arr=(${id})
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
time salmon quant -i salmon.res/GRCh38_salmon_index -l A -p 6 \
-1 clean/${fq1} -2 clean/${fq2} \
-o salmon.res/${sample}_salmon > log/salmon.log/${sample}_salmon.log 2>&1 & \
done
