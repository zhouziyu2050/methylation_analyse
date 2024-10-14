import argparse
from io import StringIO
import json
import os
import re
import sys
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.ticker import FuncFormatter

##################################################################
# 定义相关的工具方法
##################################################################


# 读取并解析带注释的JSON文件
def jsonload(file_path):
    with open(file_path, "r") as file:
        json_string = file.read()

    # 处理注释和多余的逗号
    json_string = re.sub(r"//.*", "", json_string)  # 删除单行注释
    json_string = re.sub(r"/\*.*?\*/", "", json_string, flags=re.DOTALL)  # 删除多行注释
    json_string = re.sub(r",(\s*[\]}])", r"\1", json_string)  # 删除末尾的逗号
    # 返回解析的JSON数据
    return json.loads(json_string)


class DotDict(dict):
    def __getattr__(self, item, default_value=None):
        if item in self:
            return self[item]
        else:
            return default_value

    def __setattr__(self, key, value):
        self[key] = value


# 默认值
DEFAULTS = {
    # 样本的默认参数
    "sample_name": None,
    "group_name": None,
    "input_1": "{sample_name}/{sample_name}_1.fq.gz",
    "input_2": "{sample_name}/{sample_name}_2.fq.gz",
    "output_dir": "{sample_dir}/output",
    "log_dir": "{sample_dir}/log",
    "report_dir": "{sample_dir}/report",
}


# 解析样本参数
def parse_sample_config(data):
    if "sample_name" not in data:
        raise NameError(f"sample_name不可缺失")

    sample = DotDict()
    sample["sample_name"] = data.get("sample_name", DEFAULTS["sample_name"])
    sample["group_name"] = data.get("group_name", DEFAULTS["group_name"])

    # 处理样本输入文件路径
    for key in ["input_1", "input_2"]:
        sample[key] = data.get(key, DEFAULTS[key]).format(sample_name=data["sample_name"])

    # 检查输入文件是否存在
    if not os.path.exists(sample.input_1):
        raise FileNotFoundError(f"输入文件不存在: {sample.input_1}")
    if not os.path.exists(sample.input_2):
        raise FileNotFoundError(f"输入文件不存在: {sample.input_2}")

    # 获取文件前缀，默认值为输入文件1的文件名，该参数暂时不支持手动设置
    sample["prefix"] = os.path.basename(sample["input_1"]).split(".")[0]

    # 处理输出、日志、报告目录（默认值为input_1所在文件夹）
    for key in ["output_dir", "log_dir", "report_dir"]:
        # 如果参数不存在则使用默认参数
        sample[key] = sample.get(
            key, DEFAULTS[key].format(sample_dir=os.path.dirname(sample["input_1"]))
        ).rstrip("/")
        # 如果文件夹不存在则创建
        if not os.path.exists(sample[key]):
            os.makedirs(sample[key])

    return sample


# 自动识别格式并读取samples文件为pandas格式
def read_samples_file(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError("samples文件不存在")
    if file_path.endswith(".csv"):
        # 读取 CSV 文件
        samples = pd.read_csv(file_path)
    elif file_path.endswith(".xls") or file_path.endswith(".xlsx"):
        # 读取 Excel 文件
        samples = pd.read_excel(file_path)
    elif file_path.endswith(".tsv"):
        # 读取 TSV 文件 (tab分隔符)
        samples = pd.read_csv(file_path, sep="\t")
    else:
        # 默认为tsv格式
        samples = pd.read_csv(file_path, sep="\t")
    return samples


##################################################################
# 设置命令行参数
##################################################################


parser = argparse.ArgumentParser(description="甲基化分析参数描述")
parser.add_argument("-c", "--config", type=str, help="添加配置文件（配置文件中的参数可以被命令行参数覆盖）")
parser.add_argument(
    "-r",
    "--report_dir",
    type=str,
    help="全局报告输出文件夹（不宜放在样本文件夹中的报告），默认值为：{当前文件夹}/report",
)
parser.add_argument(
    "f",
    "--samples_file",
    type=str,
    help="样本配置文件路径（从配置文件所有的样本参数，支持csv/tsv/excel格式）",
)
args = DotDict(vars(parser.parse_args()))


##################################################################
# 参数解析
##################################################################

print("参数解析...")

# 交互式命令时使用以下代码加载config文件
# args = DotDict({"config":"config.json"})

# 初始化一个空的配置参数字典
config = DotDict()
# 如果提供了配置文件，则加载配置文件的参数
if args.config:
    # 读取配置文件
    config = DotDict(jsonload(args.config))

# 用命令行参数覆盖配置文件的参数
for k, v in args.items():
    if v != None:
        config[k] = v

# 解析报告文件夹
config.report_dir = config.get("report_dir", "./report").rstrip("/")
# 如果文件夹不存在则创建
if not os.path.exists(config.report_dir):
    os.makedirs(config.report_dir)

# print(config)

# 解析样本参数
if config.samples_file:
    # 传入表格文件，则使用pandas读取表格
    samples_pd = read_samples_file(config.samples_file)
    samples = [parse_sample_config(sample) for i, sample in samples_pd.iterrows()]
elif isinstance(config.samples, list) and len(config.samples) > 0:
    # 遍历json解析样本参数
    samples = [parse_sample_config(sample) for sample in config.samples]
else:
    raise Exception("sample_file参数和config文件中的samples参数都不存在，无法读取样本信息")

print("参数解析完成")

##################################################################
# 报告生成
##################################################################

#  1. 数据基本处理与质控
# 将下机数据进行过滤，包括去污染，去测序接头和低质量碱基比例过高的reads，得到clean data。


# 表2 QC质控表
print("生成QC质控表")


def calc_qc(sample):
    data = {}
    path = f"{sample.output_dir}/bismark_alignment/{sample.prefix}_bismark_bt2_PE_report.txt"
    with open(path, "r") as file:
        for line in file:
            if line.count(":") == 1:  # 匹配包含1个冒号的行
                parts = line.split(":")
                try:
                    data[parts[0].strip()] = int(parts[1].strip())
                except:
                    pass
    item = {}
    item["Sample Id"] = sample.sample_name
    item["Clean Reads"] = data["Sequence pairs analysed in total"]
    item["Uniquely Mapped Reads"] = data["Number of paired-end alignments with a unique best hit"]
    item["Uniquely Mapping Rate"] = round(
        data["Number of paired-end alignments with a unique best hit"]
        / data["Sequence pairs analysed in total"]
        * 100,
        2,
    )
    item["Bisulfite Conversion Rate"] = round(
        (data["Total unmethylated C's in CHG context"] + data["Total unmethylated C's in CHH context"])
        / (
            data["Total unmethylated C's in CHG context"]
            + data["Total unmethylated C's in CHH context"]
            + data["Total methylated C's in CHG context"]
            + data["Total methylated C's in CHH context"]
            + data["Total methylated C's in Unknown context"]
        )
        * 100,
        2,
    )
    return item


report_qc = []
for sample in samples:
    item = calc_qc(sample)
    item = report_qc.append(item)
report_qc = pd.DataFrame(report_qc)
report_qc.to_csv(f"{config.report_dir}QC quality control table for each sample.tsv", sep="\t", index=False)


# 图3 测序深度分布图
print("绘制测序深度分布图")


def plot_methylation_depth_distribution(sample):
    input_file = f"{sample.output_dir}/{sample.sample_name}_methylation_depth_report.txt"
    output_file = f"{sample.report_dir}/Coverage of Corresponding Depth in Sample.jpg"
    df = pd.read_csv(input_file, sep="\t")
    df["C"] = df.sum(axis=1)
    df = df.set_index("Depth")
    df = df[["C", "CG", "CHG", "CHH"]]
    df = df / df.sum(axis=0)
    df.plot()
    plt.title(f"Coverage of Corresponding Depth in Sample:{sample.sample_name}")
    plt.ylabel("Fraction of Covered (%)")
    plt.legend(title="Context")
    plt.savefig(output_file, dpi=1000)
    plt.close()


for sample in samples:
    df_coverage = plot_methylation_depth_distribution(sample)


# 图4 C碱基测序深度的累积分布图
print("绘制C碱基测序深度的累积分布图")


def plot_cumulative_methylation_depth_distribution(sample):
    input_file = f"{sample.output_dir}/{sample.sample_name}_methylation_depth_report.txt"
    output_file = f"{sample.report_dir}/Cumulative Coverage of Corresponding Depth in Sample.jpg"
    df = pd.read_csv(input_file, sep="\t")
    df["C"] = df.sum(axis=1)
    df = df.set_index("Depth")
    df = df[["C", "CG", "CHG", "CHH"]]
    df = df / df.sum(axis=0)
    df = df.cumsum()  # 计算累积分布
    df.plot()
    plt.title(f"Cumulative Coverage of Corresponding Depth in Sample:{sample.sample_name}")
    plt.ylabel("Fraction of Covered (%)")
    plt.legend(title="Context")
    plt.savefig(output_file, dpi=1000)
    plt.close()


for sample in samples:
    df_coverage = plot_cumulative_methylation_depth_distribution(sample)


# 表3、表4 样品在全基因组各染色体上的C位点覆盖度统计表
print("生成C位点覆盖度统计表")


def calc_coverage_rate_by_chromosome(sample):
    input_file = f"{sample.output_dir}/{sample.sample_name}_methylation_coverage_report.txt"
    output_file = f"{sample.report_dir}/Coverage Rate Group By Chromosome.tsv"
    df = pd.read_csv(input_file, sep="\t")
    # 过滤只保留Chromosome列中以NC开头的行
    df = df[df["Chromosome"].str.startswith("NC")]
    # 计算每个 context 的覆盖率
    df["CoverageRate"] = round(df["covered"] / df["Count"] * 100, 2)

    # 计算每个染色体的整体覆盖率
    overall_coverage = (
        df.groupby("Chromosome")
        .apply(lambda x: round(x["covered"].sum() / x["Count"].sum() * 100, 2), include_groups=False)
        .reset_index()
    )
    overall_coverage.columns = ["Chromosome", "C (%)"]

    # 创建以染色体为index、context为表头的透视表，并添加总体覆盖率列
    pivot_table = df.pivot(index="Chromosome", columns="Context", values="CoverageRate")

    # 将整体覆盖率加入透视表
    pivot_table = overall_coverage.merge(pivot_table, on="Chromosome")

    pivot_table.columns = ["Chr", "C (%)", "CG (%)", "CHG (%)", "CHH (%)"]

    # 保存结果到文件
    pivot_table.to_csv(output_file, sep="\t", index=False)
    return pivot_table


df_coverage_list = []
for sample in samples:
    df_coverage = calc_coverage_rate_by_chromosome(sample)
    df_coverage.insert(0, "sample_name", sample.sample_name)
    df_coverage_list.append(df_coverage)
df_coverage_list = pd.concat(df_coverage_list, axis=0).reset_index(drop=True)
df_coverage_list.to_csv(f"{config.report_dir}/Coverage Rate Group By Chromosome.tsv", sep="\t", index=False)


# 表7 各样品QC质控表(去重之后的数据)
print("生成QC质控表(去重后)")


def calc_qc_deduplicated(sample):
    PE_report = {}
    path = f"{sample.output_dir}/bismark_alignment/{sample.prefix}_bismark_bt2_PE_report.txt"
    with open(path, "r") as file:
        for line in file:
            if line.count(":") == 1:  # 匹配包含1个冒号的行
                parts = line.split(":")
                try:
                    PE_report[parts[0].strip()] = int(parts[1].strip())
                except:
                    pass

    # 查找deduplication_percentage
    percentage_pattern = re.compile(r"Total number duplicated alignments removed:.*?(\d+\.\d+)%")
    path = f"{sample.output_dir}/bismark_deduplicate/{sample.prefix}_bismark_bt2_pe.deduplication_report.txt"
    with open(path, "r") as file:
        for line in file:
            match = percentage_pattern.search(line)
            if match:
                deduplication_percentage = match.group(1)
                break

    # 查找去重后的甲基化、非甲基化数据
    splitting_report = {}
    path = f"{sample.output_dir}/bismark_methylation/{sample.prefix}_bismark_bt2_pe.deduplicated_splitting_report.txt"
    with open(path, "r") as file:
        for line in file:
            if line.count(":") == 1:  # 匹配包含1个冒号的行
                parts = line.split(":")
                try:
                    splitting_report[parts[0].strip()] = int(parts[1].strip())
                except:
                    pass

    # 计算平均深度
    path = f"{sample.output_dir}/{sample.sample_name}_methylation_depth_report.txt"
    depth_report = pd.read_csv(path, sep="\t")
    depth_report["C"] = depth_report[["CG", "CHG", "CHH"]].sum(axis=1)
    depth_report["CxDepth"] = depth_report["C"] * depth_report["Depth"]
    avg_depth = (depth_report["CxDepth"].sum() / depth_report["C"].sum()).round(2)

    # 计算Coverage
    path = f"{sample.output_dir}/{sample.sample_name}_methylation_coverage_report.txt"
    coverage_report = pd.read_csv(path, sep="\t")
    coverage = (coverage_report["covered"].sum() / coverage_report["Count"].sum() * 100).round(2)

    item = {}
    item["Sample ID"] = sample.sample_name
    item["Clean Data Size (bp)"] = PE_report["Sequence pairs analysed in total"]
    item["Mapping Rate (%)"] = round(
        PE_report["Number of paired-end alignments with a unique best hit"]
        / PE_report["Sequence pairs analysed in total"]
        * 100,
        2,
    )
    item["Bisulfite Conversion Rate"] = round(
        (
            splitting_report["Total C to T conversions in CpG context"]
            + splitting_report["Total C to T conversions in CHH context"]
        )
        / (
            splitting_report["Total C to T conversions in CHG context"]
            + splitting_report["Total C to T conversions in CHH context"]
            + splitting_report["Total methylated C's in CHG context"]
            + splitting_report["Total methylated C's in CHH context"]
        )
        * 100,
        2,
    )
    item["Duplication Rate (%)"] = deduplication_percentage
    item["Average Depth (X)"] = avg_depth
    item["Coverage (%)"] = coverage
    return item


report_qc_deduplicated = []
for sample in samples:
    item = calc_qc_deduplicated(sample)
    report_qc_deduplicated.append(item)
report_qc_deduplicated = pd.DataFrame(report_qc_deduplicated)
report_qc_deduplicated.to_csv(
    f"{config.report_dir}/QC quality control table for each sample (deduplicated).tsv", sep="\t", index=False
)


# 2. 全基因组甲基化水平分析
# 用于分析的DNA样品为多细胞样品，因此C碱基的甲基化水平是一个0% ～ 100%范围内的数值，等于该C碱基上覆盖到的支持mC的序列数除以有效覆盖的序列总数，通常CG甲基化存在于基因和重复序列中，在基因表达调控过程中起到非常重要的作用。非CG类型的序列（CHG和CHH）在基因中十分少见，主要存在于基因间区和富含重复序列的区域，在沉默转座子过程中起关键作用。


# 表8、表9 样品全基因组及各染色体的平均甲基化水平
print("生成染色体甲基化水平统计表")


def calc_methylation_level_by_chromosome(sample_name):
    input_file = f"{sample.output_dir}/{sample.sample_name}_methylation_coverage_report.txt"
    output_file = f"{sample.report_dir}/Methylation Level Groupp By Chromosome.tsv"

    # 读取表格文件
    df = pd.read_csv(input_file, sep="\t")

    # 过滤只保留Chromosome列中以NC开头的行
    df = df[df["Chromosome"].str.startswith("NC")]

    # 计算每个 context 的甲基化水平
    df["MethylationLevel"] = round(df["totalReadsM"] / df["totalReadsN"] * 100, 2)

    # 创建以染色体为index、context为表头的透视表
    pivot_table = df.pivot(index="Chromosome", columns="Context", values="MethylationLevel")

    # 计算每个染色体的总体甲基化水平
    overall_methylation = (
        df.groupby("Chromosome")
        .apply(
            lambda x: round(x["totalReadsM"].sum() / x["totalReadsN"].sum() * 100, 2), include_groups=False
        )
        .reset_index()
    )
    overall_methylation.columns = ["Chromosome", "C (%)"]

    # 将总体甲基化水平加入透视表
    pivot_table = overall_methylation.merge(pivot_table, on="Chromosome")

    pivot_table = pivot_table.fillna(0)

    pivot_table.columns = ["Chr", "C (%)", "CG (%)", "CHG (%)", "CHH (%)"]

    # 保存结果到文件
    pivot_table.to_csv(output_file, sep="\t", index=False)

    return pivot_table


df_methylation_level_list = []
for sample in samples:
    df_methylation_level = calc_methylation_level_by_chromosome(sample)
    df_methylation_level.insert(0, "sample_name", sample.sample_name)
    df_methylation_level_list.append(df_methylation_level)
df_methylation_level_list = pd.concat(df_methylation_level_list, axis=0).reset_index(drop=True)
df_methylation_level_list.to_csv(
    f"{config.report_dir}/Methylation Level Groupp By Chromosome.tsv", sep="\t", index=False
)

# 绘制M-bias图
print("绘制M-bias图")


def plot_mbias(mbias_file, sample_name, report_dir):
    # 读取整个文件内容
    with open(mbias_file, "r") as file:
        content = file.read()

    # 定义正则表达式查找不同上下文表格
    # 正则表达式将匹配上下文标题、分隔符、表头和表数据
    pattern = r"(\w+) context \((\w+)\)\n=+\n([\s\S]+?)(?=\n[A-Z]|$)"
    matches = re.findall(pattern, content)

    df = []
    # 解析匹配到的表格
    for context, read, table_data in matches:
        # 将表格转换为 DataFrame
        item = pd.read_csv(StringIO(table_data), sep="\t")
        # 添加上下文和read标签列到 DataFrame
        item.insert(0, "context", context)
        item.insert(1, "read", read)
        df.append(item)

    df = pd.concat(df)
    for read_name in df["read"].unique():
        # 创建一个图表
        fig, ax1 = plt.subplots(figsize=(8, 4))

        # 创建组合字段
        df["read_context"] = df["read"] + " - " + df["context"]

        # 绘制甲基化率（左侧 Y 轴）
        sns.lineplot(
            data=df[df["read"] == read_name],
            x="position",
            y="% methylation",
            hue="context",
            linewidth=3,
            ax=ax1,
        )

        ax1.grid(True, linestyle="--", linewidth=0.7, color="#dddddd")
        ax1.set_xlabel("Position in Read [bp]")
        # 设置 x 和 y 轴从 0 开始
        ax1.set_xlim(left=0, right=df["position"].max())
        ax1.set_ylim(bottom=0, top=100)
        # ax1.set_ylabel("% Methylation", color='tab:blue')
        # ax1.tick_params(axis='y', labelcolor='tab:blue')
        ax1.set_ylabel("% Methylation")

        # 创建第二个 Y 轴（右侧）
        ax2 = ax1.twinx()  # 共享 x 轴

        # 绘制甲基化数量（右侧 Y 轴）
        sns.lineplot(
            data=df[df["read"] == read_name],
            x="position",
            y="count methylated",
            hue="context",
            linewidth=1,
            linestyle="--",
            ax=ax2,
        )

        # 设置第二个 Y 轴
        # ax2.set_ylabel("# Methylation Calls", color='tab:orange')
        # ax2.tick_params(axis='y', labelcolor='tab:orange')
        ax2.set_ylabel("# Methylation Calls")
        ax2.set_ylim(bottom=0)

        # 自定义格式化函数
        def format_func(value, tick_number):
            if value >= 10_000_000:
                return f"{value / 1_000_000:.0f}M"
            elif value >= 1_000_000:
                return f"{value / 1_000_000:.1f}M"
            elif value >= 10_000:
                return f"{value / 1_000:.0f}K"
            elif value >= 1_000:
                return f"{value / 1_000:.1f}K"
            else:
                return f"{value:.0f}"

        ax2.yaxis.set_major_formatter(FuncFormatter(format_func))

        read_fulname = {"R1": "Read 1", "R2": "Read 2"}[read_name]
        plt.title(f"{sample_name} {read_fulname}")

        # 显示图例
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        labels = [label + " methylation" for label in labels]
        labels2 = [label + " total calls" for label in labels2]
        ax1.legend(lines + lines2, labels + labels2, loc="upper right")
        ax2.get_legend().remove()

        plt.savefig(f"{report_dir}/M-bias {read_name}.jpg", dpi=1000)
        # plt.show()
        plt.close()


for sample in samples:
    mbias_file = (
        f"{sample.output_dir}/bismark_methylation/{sample.prefix}_bismark_bt2_pe.deduplicated.M-bias.txt"
    )
    plot_mbias(mbias_file, sample.sample_name, sample.report_dir)

# 3. 甲基化 C碱基中 CG, CHG 与CHH的分布比例
# mCG，mCHG和mCHH三种碱基类型的构成比例在不同物种中，甚至在同一物种不同样品中都存在很大差异。因此，不同时间、空间、生理条件下的样品会表现出不同的甲基化图谱，各类型mC( mCG、mCHG和mCHH )的数目，及其在全部mC的位点中所占的比例，在一定程度上反映了特定物种的全基因组甲基化图谱的特征。mCG、mCHG和mCHH分别表示表示甲基化CG、甲基化CHG和甲基化CHH。三种碱基类型占比总和为100%，甲基化C鉴定方法依据Lister的文章描述进行。


# 表12 mCG、mCHG和mCHH三种类型甲基化胞嘧啶的比例
print("生成不同context比例统计表")


def calc_context_proportion(sample):
    data = {}
    path = f"{sample.output_dir}/bismark_methylation/{sample.prefix}_bismark_bt2_pe.deduplicated_splitting_report.txt"
    with open(path, "r") as file:
        for line in file:
            if line.count(":") == 1:  # 匹配包含1个冒号的行
                parts = line.split(":")
                try:
                    data[parts[0].strip()] = int(parts[1].strip())
                except:
                    pass
    item = {}
    item["Sample Id"] = sample.sample_name
    item["mCG"] = data["Total methylated C's in CpG context"]
    item["mCHG"] = data["Total methylated C's in CHG context"]
    item["mCHH"] = data["Total methylated C's in CHH context"]
    total_mC = item["mCG"] + item["mCHG"] + item["mCHH"]
    item["mCG proportion (%)"] = round(item["mCG"] / total_mC * 100, 2)
    item["mCHG proportion (%)"] = round(item["mCHG"] / total_mC * 100, 2)
    item["mCHH proportion (%)"] = round(item["mCHH"] / total_mC * 100, 2)

    return item


df_context_proportion = []
for sample in samples:
    item = calc_context_proportion(sample)
    df_context_proportion.append(item)
df_context_proportion = pd.DataFrame(df_context_proportion).set_index("Sample Id")
df_context_proportion.to_csv(
    f"{config.report_dir}/Proportions of Three Types of Methylated Cytosine.tsv", sep="\t"
)


# 图5 不同序列类型甲基化C碱基的分布比例
print("绘制不同context比例统计图")


def plot_context_proportion(sample, proportion):
    proportion.plot(kind="pie", autopct="%1.1f%%", figsize=(5, 5))
    plt.title(f"The proportion of every type mC (Sample: {sample.sample_name})")
    plt.legend()
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(f"{sample.report_dir}/The proportion of every type mC.jpg", dpi=1000)
    plt.close()


for sample in samples:
    proportion = df_context_proportion.loc[sample.sample_name]
    plot_context_proportion(sample, proportion)


# 4. 甲基化 CG、CHG和CHH的甲基化水平分布
# 不同类型的C碱基(mCG、mCHG和mCHH )，其甲基化水平在不同物种间，甚至同一物种不同细胞类型不同条件下其甲基化水平都存在差异。此图统计每种类型( CG、CHG和CHH )甲基化C的甲基化水平分布，反映了该物种DNA甲基化特征


# 绘制甲基化水平分布图（按每10%分组）
print("绘制甲基化水平分布图")


def plot_methylation_level_distribution(sample):
    input_file = f"{sample.output_dir}/{sample.sample_name}_methylation_distribution_report.txt"
    output_file = f"{sample.report_dir}/Methylation Level Distribution.jpg"

    # 读取TSV格式的输出文件
    df = pd.read_csv(input_file, sep="\t")
    df = df[df["methylation_level"] > 0]

    # 确保数据按 context 和 methylation_level 排序
    df.sort_values(by=["context", "methylation_level"], inplace=True)

    # 将数据分组到指定的区间
    df = bin_methylation_levels(df, range_size=10)

    # 从长格式转为宽格式
    df = df.pivot(index="methylation_level_range", columns="context", values="count")
    df.columns.name = None
    df.index.name = None

    # 计算总甲基化数量
    df["C"] = df[["CG", "CHG", "CHH"]].sum(axis=1)
    df = df[["C", "CG", "CHG", "CHH"]]
    df = df / df.sum() * 100  # 计算不同甲基化水平的占比

    # 创建绘图
    df.plot(marker="o", figsize=(6, 4))
    plt.xlabel("Methylation Level")
    plt.ylabel("Percentage (%)")
    plt.title(f"Methylation Level Distribution (Sample: {sample.sample_name})")
    plt.legend(title="Context")
    plt.grid(True, linewidth=1, color="#f6f6f6")
    plt.tight_layout()
    plt.savefig(output_file, dpi=1000)
    plt.close()


def bin_methylation_levels(df, range_size=10):
    """将数据分组到区间"""
    # 创建 bin 的区间
    bin_edges = np.arange(0, 101, range_size)
    bin_labels = [f"{i+range_size}" for i in bin_edges[:-1]]
    df["methylation_level_range"] = pd.cut(
        df["methylation_level"], bins=bin_edges, labels=bin_labels, right=False
    )
    # 计算每个 methylation_level_range 中的 count 总和
    range_df = df.groupby(["context", "methylation_level_range"], observed=False)["count"].sum().reset_index()
    return range_df


for sample in samples:
    plot_methylation_level_distribution(sample)


# 绘制甲基化水平累积分布图（按每1%分组）
print("绘制甲基化水平累积分布图")


def plot_cumulative_methylation_level_distribution(sample):
    input_file = f"{sample.output_dir}/{sample.sample_name}_methylation_distribution_report.txt"
    output_file = f"{sample.report_dir}/Cumulative Methylation Level Distribution.jpg"
    # 读取TSV格式的输出文件
    df = pd.read_csv(input_file, sep="\t")
    df = df[df["methylation_level"] > 0]

    # 确保数据按 context 和 methylation_level 排序
    df.sort_values(by=["context", "methylation_level"], inplace=True)

    # 从长格式转为宽格式
    df = df.pivot(index="methylation_level", columns="context", values="count")
    df.columns.name = None
    df.index.name = None

    # 计算总甲基化数量
    df["C"] = df[["CG", "CHG", "CHH"]].sum(axis=1)
    df = df[["C", "CG", "CHG", "CHH"]]
    df = df / df.sum() * 100  # 计算不同甲基化水平的占比
    df = df.cumsum()  # 计算累积分布

    # 创建绘图
    df.plot(figsize=(6, 4))
    plt.xlabel("Methylation Level")
    plt.ylabel("Percentage (%)")
    plt.title(f"Cumulative Methylation Level Distribution (Sample: {sample.sample_name})")
    plt.legend(title="Context")
    plt.grid(True, linewidth=1, color="#f6f6f6")
    plt.tight_layout()
    plt.savefig(output_file, dpi=1000)
    plt.close()


for sample in samples:
    plot_cumulative_methylation_level_distribution(sample)

print("质控相关报告已全部生成")
