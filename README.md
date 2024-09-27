# methylation_analyse
基于bismark的全流程甲基化分析，包括数据清洗、序列比对、甲基化提取及DMR分析、GO/KEGG等


# 关键分析步骤
## bismark分析及质控
分析操作程序（耗时长）：```methylation_analyse.py```
质控报告生成程序：``qc_report.ipynb``

## DMR分析及绘图
DMR分析程序（耗时长）：```DMR_analyse.R```
DMR绘图程序:```DMR_plot.ipynb```

## GO & KEGG分析
分析程序：```GO_and_KEGG_analyse.ipynb```

# 运行环境部署过程
查看：```Environmental setup.md```

# report输出结构
```
report 报告文件夹
│
├── 样本文件夹
│   ├── output 输出数据的文件夹
│   │   ├── *_1_bismark_bt2_PE_report.txt                      # 甲基化比对报告
│   │   ├── *_1_bismark_bt2_pe.deduplicated.sort.bam           # 去重后的bam文件
│   │   ├── *_1_bismark_bt2_pe.deduplicated.sort.bam.bai       # 去重后bam文件的索引
│   │   ├── *_1_bismark_bt2_pe.deduplication_report.txt        # 去重报告
│   │   ├── *_1_bismark_bt2_pe.deduplicated_splitting_report.txt # 甲基化分析报告
│   │   ├── *_1_bismark_bt2_pe.deduplicated.bedGraph.gz        # 甲基化结果的bedGraph格式文件
│   │   └── *_1_bismark_bt2_pe.deduplicated.bigwig             # 甲基化结果的bigWig格式文件
│   │
│   ├── Coverage of Corresponding Depth in Sample.jpg          # 甲基化深度分布图，横坐标是覆盖深度（按1%统计），纵坐标是占比
│   ├── Coverage Rate Group By Chromosome.tsv                  # 每条染色体上的甲基化区域覆盖度
│   ├── Cumulative Coverage of Corresponding Depth in Sample.jpg # 甲基化深度累积分布图，横坐标是覆盖深度（按1%统计），纵坐标是占比
│   ├── M-bias R1.jpg                                          # M-bias分析图（读数R1）
│   ├── M-bias R2.jpg                                          # M-bias分析图（读数R2）
│   ├── Methylation Level Distribution.jpg                     # 甲基化水平分布图，横坐标是甲基化水平（按10%统计），纵坐标是占比
│   ├── Methylation Level Groupp By Chromosome.tsv             # 每条染色体上的甲基化水平
│   └── The proportion of every type mC.jpg                    # 不同context的占比图
│
├── 组间对照文件夹
│   ├── circos_plot.png                                        # DMR区域及相关基因在染色体上的分布环形图
│   ├── DMR_gene_all.tsv                                       # DMR分析的结果及完整的基因注释
│   ├── DMR_genes.tsv                                          # 仅保留基因的关键字段并去重
│   ├── DMR_summary.tsv                                        # 不同染色体上的DMR数量统计
│   └── GO/KEGG富集-(all/gain/loss).(png/tsv)                  # GO/KEGG富集分析结果
│
├── Coverage Rate Group By Chromosome.tsv                      # 每条染色体上的甲基化区域覆盖度
├── Methylation Level Groupp By Chromosome.tsv                 # 每条染色体上的甲基化水平
├── Proportions of Three Types of Methylated Cytosine.tsv      # 每个样本不同context的甲基化比例
├── QC quality control table for each sample (deduplicated).tsv # 去重后的质控数据
└── QC quality control table for each sample.tsv               # 质控数据
```

# 染色体名称对应关系

| Accession   | Chromosome |
|-------------|------------|
| NC_000067.7 | chr1       |
| NC_000068.8 | chr2       |
| NC_000069.7 | chr3       |
| NC_000070.7 | chr4       |
| NC_000071.7 | chr5       |
| NC_000072.7 | chr6       |
| NC_000073.7 | chr7       |
| NC_000074.7 | chr8       |
| NC_000075.7 | chr9       |
| NC_000076.7 | chr10      |
| NC_000077.7 | chr11      |
| NC_000078.7 | chr12      |
| NC_000079.7 | chr13      |
| NC_000080.7 | chr14      |
| NC_000081.7 | chr15      |
| NC_000082.7 | chr16      |
| NC_000083.7 | chr17      |
| NC_000084.7 | chr18      |
| NC_000085.7 | chr19      |
| NC_000086.8 | chrX       |
| NC_000087.8 | chrY       |
| NC_005089.1 | chrMT      |
