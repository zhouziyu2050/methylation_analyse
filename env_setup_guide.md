## docker操作
```
## 镜像拉取
docker pull docker.anyhub.us.kg/continuumio/anaconda3

## 创建容器
docker run -itd \
  -p 22224:22 \
  --name methylation \
  -v ~/methylation:/methylation \
  --cpus=62 \
  --memory=120g \
  --memory-swap=120g \
  --oom-kill-disable \
  methylation:v1.0


## 进入容器
docker exec -it methylation /bin/bash

## 进入甲基化分析数据目录
cd /methylation/F24A080000424_MUSekgzH_20240805100100/
```



## conda操作

**conda换源**

需注意：必须先换源再安装虚拟环境

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

**conda虚拟环境安装** 

conda的相关依赖及其注释被保存在[environment.yml](environment.yml)

```
conda env create -f environment.yml -n methylation
```

**其他依赖安装**

```
# pip换源
pip config set global.index-url https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple

# python依赖安装
pip install jupyter pandas matplotlib seaborn

# C语言脚本编译依赖
apt-get install gcc zlib1g-dev
apt-get install g++ libjsoncpp-dev

# 以下在R交互式命令行中安装
devtools::install_github("junjunlab/transPlotR") # 绘制基因转录本位置
IRkernel::installspec(user = FALSE)  # 将IRkernel注册到jupyter notebook
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
