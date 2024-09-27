#!/bin/bash

# 样本名称列表
sample_names=("13A" "32A" "50A" "26A" "28A" "38A" "74A" "25A" "75A" "88A" "97A")

# 基础目录
base_dir="/methylation/F24A080000424_MUSekgzH_20240805100100/"

# 循环遍历每个样本
for sample in "${sample_names[@]}"; do
    echo "Processing sample $sample..."

    # 排序 BAM 文件
    samtools sort -@8 -m 4G ${base_dir}${sample}/output/bismark_deduplicate/${sample}_1_bismark_bt2_pe.deduplicated.bam -o ${base_dir}${sample}/output/bismark_deduplicate/${sample}_1_bismark_bt2_pe.deduplicated.sort.bam
    echo "$sample sorted."

    # 索引排序后的 BAM 文件
    samtools index -@8 ${base_dir}${sample}/output/bismark_deduplicate/${sample}_1_bismark_bt2_pe.deduplicated.sort.bam
    echo "$sample indexed."

    echo "$sample processing completed."
done

echo "All samples processed."
