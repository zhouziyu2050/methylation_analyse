## docker操作
```
## 镜像拉取
docker pull docker.anyhub.us.kg/continuumio/anaconda3

## 创建容器
docker run -itd \
  -p 22224:22 \
  --name methylationv0.4 \
  -v ~/methylation:/methylation \
  --cpus=62 \
  --memory=120g \
  --memory-swap=120g \
  --oom-kill-disable \
  methylation:v0.4


## 进入容器
docker exec -it methylationv0.3 /bin/bash

## 退出容器
exit

## 进入甲基化分析数据目录
cd /methylation/F24A080000424_MUSekgzH_20240805100100/
```

## conda换源

```
conda config --set show_channel_urls yes
vim ~/.condarc
```

```
channels:
  - https://mirrors.pku.edu.cn/anaconda/pkgs/main
  - https://mirrors.pku.edu.cn/anaconda/pkgs/r
  - conda-forge
  - bioconda
  - defaults
show_channel_urls: true
channel_priority: disabled
```

## conda操作

```
# 创建环境
conda create -n methylathion python==3.10.14
```

## 依赖安装
```
conda install shuaizhou::soapnuke==2.1.9  # 安装后调用名称为SOAPnuke
# conda install bioconda::bsmap
conda install bioconda::bismark  # 主要包含bowtie2、samtools、perl、libgcc-ng等，其中libgcc-ng来自conda-forge源
conda install bioconda::samtools==1.20  # 安装指定版本的samtools
# DMR注释
conda install bioconda::bedtools
# 用于gtf2bed
conda install bioconda::bedops
# 用于DMR注释
# conda install bioconda::metilene
# 用于DMRs导出bed文件、读取gtf文件
conda install bioconda::bioconductor-rtracklayer

conda install bioconda::homer

# 安装R及依赖
conda install r-base==4.3.3
# 用于计算DMR
conda install bioconda::bioconductor-dmrcaller
# 用于基因富集
conda install bioconda::bioconductor-clusterprofiler
conda install conda-forge::libxml2 # conda安装clusterprofiler时可能遗漏该包，手动安装
# 数据处理工具集
# conda install conda-forge::r-tidyverse
# 用于将ensembl gene id转gene_ids
conda install bioconda::bioconductor-biomart
# 安装注释库
conda install bioconda::bioconductor-txdb.mmusculus.ucsc.mm39.knowngene
# 安装ggbio
conda install bioconda::bioconductor-ggbio 或 BiocManager::install("ggbio")
# 安装devtools（装transPlotR）
conda install conda-forge::r-devtools

# 以下在R交互式命令行中安装
install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db") # 小鼠基因注释数据库
install.packages("circlize")  # 绘制Circos图
devtools::install_github("junjunlab/transPlotR") # 绘制基因转录本位置
install.packages("dplyr") # 数据处理

# 在jupyter notebook中使用R语言
install.packages('IRkernel')  # 安装IRkernel
IRkernel::installspec(user = FALSE)  # 在jupyter notebook中安装

# python依赖安装
pip install jupyter pandas matplotlib seaborn

# C语言脚本编译依赖
apt-get install gcc zlib1g-dev
apt-get install g++ libjsoncpp-dev
```

## 自定义脚本编译
```
gcc -o refseq2chr refseq2chr.c -lz
gcc -o methylation_coverage_analyse methylation_coverage_analyse.c -lz
gcc -o methylation_distribution_analysis methylation_distribution_analysis.c -lm -lz
gcc -o methylation_depth_analysis methylation_depth_analysis.c -lm -lz
```

## bashrc修改
```
PS1="\u@\h:\W\$ "
conda activate methylation
cd /methylation
export PATH="/methylation/utils:$PATH"

# Check if the SSH service is running
if ! pgrep -x "sshd" > /dev/null
then
    service ssh start
fi
```

## locals优化

```
locale  # 查看当前设置的locale变量
locale -a  # 查看当前已安装的locale，发现设置的en_US.UTF-8未安装
apt-get install locales  # 安装locales
dpkg-reconfigure locales  # 选择安装en_US.UTF-8
```

## Homer使用
在 http://homer.ucsd.edu/homer/download.html 下载```configureHomer.pl```文件，用于homer及其数据集的下载、安装
```
# 查看所有内置的物种
perl configureHomer.pl --list
# 下载mm39注释数据集（需切换到下载目录）
perl configureHomer.pl  -install mm39
```


## bedGraphToBigWig操作步骤
```
# 依赖安装
conda install bioconda::ucsc-bedgraphtobigwig
# 构建基因组索引
gzip -dk /methylation/genome/mm39/GCF_000001635.27_GRCm39_genomic.fa.gz
samtools faidx /methylation/genome/mm39/GCF_000001635.27_GRCm39_genomic.fa
# 转换 (已集成脚本utils/bedgraph2bigwig.sh)
gzip -dk 13A/output/bismark_methylation/13A_1_bismark_bt2_pe.deduplicated.bedGraph.gz
bedGraphToBigWig  13A/output/bismark_methylation/13A_1_bismark_bt2_pe.deduplicated.bedGraph  /methylation/genome/mm39/GCF_000001635.27_GRCm39_genomic.fa.fai  13A/output/bismark_methylation/13A_1_bismark_bt2_pe.deduplicated.bigwig
```
