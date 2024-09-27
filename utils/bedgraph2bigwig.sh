#!/bin/bash

# 样本列表
samples=("13A" "32A" "50A" "26A" "28A" "38A" "74A" "25A" "75A" "88A" "97A")

# 基础目录
base_dir="/methylation/F24A080000424_MUSekgzH_20240805100100/"
# 参考文件
genome_fai="/methylation/genome/mm39/GCF_000001635.27_GRCm39_genomic.fa.fai"

# 遍历样本列表
for sample in "${samples[@]}"; do
    # 动态构建路径
    sample_base_dir="${base_dir}${sample}/output/bismark_methylation"
    bedgraph_gz="${sample_base_dir}/${sample}_1_bismark_bt2_pe.deduplicated.bedGraph.gz"
    bedgraph="${sample_base_dir}/${sample}_1_bismark_bt2_pe.deduplicated.bedGraph"
    bigwig="${sample_base_dir}/${sample}_1_bismark_bt2_pe.deduplicated.bigwig"

    # 解压 bedGraph 文件
    echo "Processing sample ${sample}..."
    echo "Unzipping ${bedgraph_gz}..."
    gzip -dk "$bedgraph_gz"

    # 转换 bedGraph 到 bigWig
    echo "Converting ${bedgraph} to bigWig..."
    bedGraphToBigWig "$bedgraph" "$genome_fai" "$bigwig"

    # 删除中间文件
    echo "Cleaning up..."
    rm "$bedgraph"

    echo "Sample ${sample} processed."
    echo "---------------------------------"
done
