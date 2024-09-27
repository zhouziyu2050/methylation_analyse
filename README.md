# methylation_analyse
基于bismark的全流程甲基化分析，包括数据清洗、序列比对、甲基化提取及DMR分析、GO/KEGG等

# 运行环境部署过程

部署过程说明：[env_setup_guide.md](env_setup_guide.md)

# 关键分析步骤

## 1. bismark分析及质控

分析操作程序（耗时长）：[methylation_analyse.py](methylation_analyse.py)

质控报告生成程序（Python）：[qc_report.ipynb](qc_report.ipynb)

## 2. DMR分析及绘图

DMR分析程序（耗时长）：[DMR_analyse.R](DMR_analyse.R)

DMR绘图程序（R语言）：[DMR_plot.ipynb](DMR_plot.ipynb)

## 3. GO & KEGG分析

分析程序（R语言）：[GO_and_KEGG_analyse.ipynb](GO_and_KEGG_analyse.ipynb)

# report输出结构
```
report 报告文件夹
│
├── 样本文件夹
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
参考：[Chromosome comparison table.ipynb](Chromosome comparison table.ipynb)
