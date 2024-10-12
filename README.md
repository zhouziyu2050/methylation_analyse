# methylation_analyse
基于bismark的全流程甲基化分析，包括数据清洗、序列比对、甲基化提取及DMR分析、GO/KEGG等

# 运行环境部署

```
# 拉取docker镜像
docker push zhouziyu2050/methylation:latest

# 创建容器（映射ssh端口、映射基因文件夹、设置资源上限、设置禁止kill）
docker run -itd \
  -p 2222:22 \
  --name methylation \
  -v ~/methylation:/methylation \
  --cpus=62 \
  --memory=120g \
  --memory-swap=120g \
  --oom-kill-disable \
  zhouziyu2050/methylation:latest

# 进入容器
docker exec -it methylation /bin/bash
```

环境搭建完整步骤查看：[env_setup_guide.md](env_setup_guide.md)

# 关键步骤

## 1. 甲基化分析及质控

**1.1 甲基化分析程序**（耗时长）：[methylation_analyse.py](methylation_analyse.py)

参数描述：

| 参数                           | 默认值                           | 描述                                                   |
|--------------------------------|-----------------------------------|-------------------------------------------------------|
| `-h`, `--help`                 |                                   | 显示帮助信息                                           |
| **全局参数**                   |                                   |                                                        |
| `--config <file>`              | `NULL`                            | 配置文件路径（从配置文件读取所有参数，json格式，示例文件：[config.json](config.json)）    |
| `--genome_folder <folder>`     | `NULL`                            | 参考基因组文件所在文件夹的路径（必传）                  |
| `--utils_folder <folder>`      | `{当前文件夹}`                    | utils文件夹的路径，默认值为当前文件夹                   |
| `--skip_filter`                | `false`                           | 添加该参数以跳过数据清洗步骤                            |
| `--parallel_num <num>`         | `30`                              | 最大使用线程数                                          |
| `--parallel_alignment <num>`   | `6`                               | 对齐比对的线程数，线程过多容易内存溢出                  |
| `--samples_file <file>`        | `NULL`                            | 样本配置文件路径（从配置文件读取样本参数，支持csv/tsv/excel格式，示例文件：[config_samples.tsv](config_samples.tsv)） |
| **样本参数**                   |                                   | 可从命令行中输入单个样本的参数                          |
| `--sample_name <name>`         | `NULL`                            | 样本名（必传）                                          |
| `--group_name <name>`          | `NULL`                            | 样本所属分组（必传）                                    |
| `--input_1 <path>`             | `{sample_name}/{sample_name}_1.fq.gz` | 测序文件1的路径                                     |
| `--input_2 <path>`             | `{sample_name}/{sample_name}_2.fq.gz` | 测序文件2的路径                                     |
| `--output_dir <folder>`        | `{input_1所在文件夹}/output`    | 输出的中间文件存放路径                                     |
| `--log_dir <folder>`           | `{input_1所在文件夹}/log`       | 日志文件夹                                                 |
| `--report_dir <folder>`        | `{input_1所在文件夹}/report`    | 报告输出路径                                               |

注：

- 使用config中的samples参数或samples_file配置文件可以传入多个样本的参数，通过命令行只能传入单个样本的参数。

- 若设置了配置文件`config`，其他所有参数都仅从配置文件读取。推荐使用`config`文件配置参数，后续步骤可以复用。

- 若设置了样本配置文件`samples_file`，所有样本参数都仅从该配置文件读取。

- 为了方便阅读，配置文件中可以使用```//```和```/* */```注释，程序解析时会自动忽略注释内容。

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


**1.2 质控报告生成程序**（Python）：[qc_report.ipynb](qc_report.ipynb)

质控数据图表：

| 名称                     | 文件名                                               | 描述                                         |
|--------------------------|------------------------------------------------------|-----------------------------------------------|
| 测序质量基本统计表        | Clean_stat.xls                                       | 数据清洗后的基本统计数据                      |
| QC质控表                  | QC quality control table for each sample.tsv         | 与参考基因组比对结果的质控统计信息            |
| 测序深度分布图            | Coverage of Corresponding Depth in Sample.jpg        | 每个样本中不同C碱基在不同测序深度的分布情况    |
| 测序深度累积分布图        | Cumulative Coverage of Corresponding Depth in Sample.jpg | 每个样本中不同C碱基在不同测序深度的累积分布情况 |
| C位点覆盖度统计表         | Coverage Rate Group By Chromosome.tsv                | 样品在全基因组各染色体上的C位点覆盖度          |
| QC质控表(去重后)          | QC quality control table for each sample (deduplicated).tsv | 去重后的QC质控分析                     |
| 染色体甲基化水平统计表    | Methylation Level Groupp By Chromosome.tsv           | 样品全基因组及各染色体的平均甲基化水平         |
| M-bias图                  | M-bias Read {1/2}.jpg                                |                                               |
| 不同context比例统计表     | Proportions of Three Types of Methylated Cytosine.tsv | mCG、mCHG和mCHH三种类型甲基化胞嘧啶的比例统计 |
| 不同context比例饼图       | The proportion of every type mC.jpg                  | mCG、mCHG和mCHH三种类型甲基化胞嘧啶的比例饼图   |
| 甲基化水平分布图          | Methylation Level Distribution.jpg                   | 不同甲基化比例的分布情况（按每10%分组）         |
| 甲基化水平累积分布图      | Cumulative Methylation Level Distribution.jpg        | 不同甲基化比例的累积分布情况（按每1%分组）       |

## 2. DMR分析及绘图

**2.1 DMR分析程序**（耗时长）：[DMR_analyse.R](DMR_analyse.R)

参数描述：

| 参数                  | 默认值                    | 描述                           |
|-----------------------|--------------------------|--------------------------------|
| `-h`, `--help`        |                          | 显示帮助信息                    |
| `-c`, `--config`      | `NULL`                   | 配置文件路径                    |
| `-o`, `--output_dir`  | `{当前文件夹}/output`     | 中间文件的输出文件夹            |
| `-r`, `--report_dir`  | `{当前文件夹}/report`     | 报告的输出文件夹                |
| `-a`, `--group_a`     | `NULL`                   | DMR的组A名称 (必传)             |
| `-b`, `--group_b`     | `NULL`                   | DMR的组B名称 (必传)             |
| `-f`, `--samples_file`| `NULL`                   | 以tsv/csv/excel文件传入样本参数 |
| `-g`, `--gtf_file`    | `NULL`                   | gtf注释文件路径，支持gtf/gtf.gz格式（必传） |

注：
- 同时使用配置文件和命令行参数时，命令行参数会覆盖配置文件的参数。
- `config`及`samples_file`参数，推荐复用第1步中的文件。
- gtf注释文件下载地址：[https://www.gencodegenes.org/](https://www.gencodegenes.org/)

**2.2 DMR绘图程序**：[DMR_plot.R](DMR_plot.R)

参数描述

| 参数名                 | 默认值                   | 描述                                                       |
|------------------------|---------------------------|------------------------------------------------------------|
| **全局参数**           |                           |                                                           |
| `-h`, `--help`         |                           | 显示帮助信息                                               |
| `-c`, `--config`       | `NULL`                    | 配置文件路径                                               |
| `-o`, `--output_dir`   | `{当前文件夹}/output`     | 中间文件的输出文件夹                                       |
| `-r`, `--report_dir`   | `{当前文件夹}/report`     | 报告的输出文件夹                                           |
| `-a`, `--group_a`      | `NULL`                    | DMR的组A名称                                               |
| `-b`, `--group_b`      | `NULL`                    | DMR的组B名称                                               |
| `-f`, `--samples_file` | `NULL`                    | 以tsv/csv/excel文件传入样本参数                            |
| **DMR和甲基化位置分布图**  |                       |                                                           |
| `-p`, `--plot_type`    | `NULL`                    | DMR和甲基化位置分布图的绘制形式，可选值为`line/bar/point`，不传则不绘制此图 |
| `-n`, `--seqname`      | `NULL`                    | 绘制DMR和甲基化位置分布图的染色体名称                     |
| `-s`, `--start`        | `NULL`                    | 绘制DMR和甲基化位置分布图的起始位置                       |
| `-e`, `--end`          | `NULL`                    | 绘制DMR和甲基化位置分布图的结束位置                       |
| `-g`, `--gtf_file`     | `NULL`                    | 绘制DMR和甲基化位置分布图所使用的gtf注释文件路径，支持gtf/gtf.gz格式 |
| **DMR环形分布图**      |                           |                                                          |
| `-y`, `--cytoband_path`| `NULL`                    | cytoband文件路径，不传则不绘制DMR环形分布图               |
| `-t`, `--text_num`     | `88`                      | 在环内显示标签的数量（微调该参数使标签数量刚好铺满整个环） |

注：
- 同时使用配置文件和命令行参数时，命令行参数会覆盖配置文件的参数。
- `config`及`samples_file`参数，推荐复用第1步中的文件。
- gtf注释文件下载地址：[https://www.gencodegenes.org/](https://www.gencodegenes.org/)
- cytoband文件下载地址：[https://hgdownload.cse.ucsc.edu/goldenPath/mm39/database/cytoBandIdeo.txt](https://hgdownload.cse.ucsc.edu/goldenPath/mm39/database/cytoBandIdeo.txt)

## 3. GO & KEGG分析

GO & KEGG分析程序：[GO_and_KEGG_analyse.R](GO_and_KEGG_analyse.R)

参数描述：

| 参数                        | 默认值               | 描述                                          |
|-----------------------------|---------------------|------------------------------------------------|
| `-h`, `--help`              |                     | 显示帮助信息                                   |
| `-g`, `--genes`             | `NULL`              | DMR输出的基因文件路径，可以使用相对路径或绝对路径（必传） |
| `-r`, `--report_dir`        | `{genes所在文件夹}`  | 输出报告的文件夹路径，可选             |
| `-p`, `--pathways_selected` | `NULL`              | 指定通路，多个通路以英文逗号连接，如 'GO:0007015,GO:0007264'，可选 |

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
│   ├── DMR_gene_all.tsv                                       # DMR分析的结果及完整的基因注释
│   ├── DMR_genes.tsv                                          # 仅保留基因的关键字段并去重
│   ├── DMR_summary.tsv                                        # 不同染色体上的DMR数量统计
│   ├── Circos plot of DMR.png                                 # DMR环形分布图
│   ├── Position of DMR and methylation.png                    # DMR和甲基化位置分布图
│   └── GO/KEGG富集-{all/gain/loss}.{png/tsv}                  # GO/KEGG富集分析结果
│
├── Coverage Rate Group By Chromosome.tsv                      # 每条染色体上的甲基化区域覆盖度
├── Methylation Level Groupp By Chromosome.tsv                 # 每条染色体上的甲基化水平
├── Proportions of Three Types of Methylated Cytosine.tsv      # 每个样本不同context的甲基化比例
├── QC quality control table for each sample (deduplicated).tsv # 去重后的质控数据
└── QC quality control table for each sample.tsv               # 质控数据
```

# 染色体名称对应关系
参考对照表：[Chromosome_comparison_table.md](Chromosome_comparison_table.md)
