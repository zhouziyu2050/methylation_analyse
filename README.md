# methylation_analyse
基于bismark的全流程甲基化分析，包括数据清洗、序列比对、甲基化提取及DMR分析、GO/KEGG等

# 运行环境部署过程

部署过程说明：[env_setup_guide.md](env_setup_guide.md)

# 关键步骤

## 1. 甲基化分析及质控

**甲基化分析程序**（耗时长）：[methylation_analyse.py](methylation_analyse.py)

关键参数描述：

| 参数                           | 默认值                           | 描述                                                   |
|--------------------------------|-----------------------------------|-------------------------------------------------------|
| `-h, --help`                   |                                   | 显示帮助信息并退出                                     |
| `--config <file>`              |                                   | 添加配置文件（示例文件：[config.json](config.json)）    |
| `--genome_folder <folder>`     |                                   | 参考基因组文件所在文件夹的路径（必传）                  |
| `--utils_folder <folder>`      | {当前文件夹}                      | utils文件夹的路径，默认值为当前文件夹                   |
| `--skip_filter`                |                                   | 添加该参数以跳过数据清洗步骤                            |
| `--parallel_num <num>`         | 30                                | 最大使用线程数                                          |
| `--parallel_alignment <num>`   | 6                                 | 对齐比对的线程数，线程过多容易内存溢出                  |
| `--sample_name <name>`         |                                   | 样本名（必传）                                          |
| `--group_name <name>`          |                                   | 样本所属分组（必传）                                    |
| `--input_1 <path>`             | `{sample_name}/{sample_name}_1.fq.gz` | 测序文件1的路径                                     |
| `--input_2 <path>`             | `{sample_name}/{sample_name}_2.fq.gz` | 测序文件2的路径                                     |
| `--output_dir <folder>`        | `{input_1所在文件夹}/output`    | 输出的中间文件存放路径                                     |
| `--log_dir <folder>`           | `{input_1所在文件夹}/log`       | 日志文件夹                                                 |
| `--report_dir <folder>`        | `{input_1所在文件夹}/report`    | 报告输出路径                                               |

注：

1 使用config配置文件可以同时传入多个样本，否则只能每次传入一个样本。

2 如果使用了配置文件，则其它参数均不生效。

3 为了方便使用，配置文件中可以使用```//```和```/* */```注释，程序解析时会自动忽略注释内容。

参考基因组文件下载地址：[mm39小鼠基因组](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/) , [其他基因组](https://www.ncbi.nlm.nih.gov/datasets/genome/)

主要分析步骤：

| 步骤 | 程序来源 | 调用程序  | 步骤描述   |
|------|----------|------------|------------|
| 1        | bismark  | [bismark_genome_preparation](https://felixkrueger.github.io/Bismark/options/genome_preparation/) | 创建参考基因组的索引文件（检测到文件已存在则自动跳过）|
| 2        | SOAPnuke | [soapnuke_filter](https://github.com/BGI-flexlab/SOAPnuke) | 进行数据过滤|
| 3        | bismark  | [bismark_alignment](https://felixkrueger.github.io/Bismark/options/alignment/) | 执行序列比对|
| 4        | bismark  | [bismark_deduplicate](https://felixkrueger.github.io/Bismark/options/deduplication/) | 去除重复片段|
| 5        | bismark  | [bismark_methylation_extractor](https://felixkrueger.github.io/Bismark/options/methylation_extraction/) | 提取甲基化信息|
| 6        | C语言脚本 | [methylation_depth_analysis](utils/methylation_depth_analysis) | 输出甲基化测序深度信息|
| 7        | C语言脚本 | [methylation_coverage_analyse](utils/methylation_coverage_analyse) | 输出基于染色体和context的甲基化覆盖度信息|
| 8        | C语言脚本 | [methylation_distribution_analysis](utils/methylation_distribution_analysis) | 输出基于染色体和context的甲基化分布信息|


**质控报告生成程序**（Python）：[qc_report.ipynb](qc_report.ipynb)


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
参考对照表：[Chromosome_comparison_table.md](Chromosome_comparison_table.md)
