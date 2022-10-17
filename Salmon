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
