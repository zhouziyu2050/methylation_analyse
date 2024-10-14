import argparse
import json
import re
import subprocess
import datetime
import os
import signal
import sys
import time


# 用于将字典参数构造成命令字符串
def dict2cmd(prefix, params):
    cmd = prefix
    for param, value in params.items():
        if value == None or value == "":
            cmd += f" {param}"
        else:
            cmd += f" {param} {value}"
    # print("执行命令：",cmd)
    return cmd


# 创建中间文件输出目录
def mkdirs(sample, config):
    directories = [
        "bismark_alignment",
        "bismark_alignment/temp",
        "bismark_methylation",
        "bismark_deduplicate",
    ]
    if not config.skip_filter:
        directories.append("soapnuke")
    directories = [f"{sample.output_dir}/{x}" for x in directories]
    directories.append(sample.report_dir)
    directories.append(sample.log_dir)

    params = {
        "-p": "",  # 创建多级目录
        "-m": 777,  # 目录权限
        " ".join(directories): "",  # 目录列表
    }
    cmd = dict2cmd("mkdir", params)
    return cmd


# 1.创建参考基因组的索引文件（每个参考基因组只需要执行一次）
def bismark_genome_preparation(config):
    # 文档地址：https://felixkrueger.github.io/Bismark/options/genome_preparation/
    # 基因组下载地址：https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/

    # 定义参数字典
    params = {
        "--bowtie2": "",  # 使用 Bowtie2 作为比对工具
        "--parallel": config.parallel_num // 2,  # 每个线程实际使用2个核心，所以这里减半
        config.genome_folder: "",  # 参考基因组文件夹
    }
    # 构造命令字符串
    cmd = dict2cmd("bismark_genome_preparation", params)
    return cmd


# 2.使用SOAPnuke做数据过滤（已过滤则不需要这一步）
def soapnuke_filter(sample, config):
    # 文档地址：https://github.com/BGI-flexlab/SOAPnuke/blob/master/Readme.md
    input_file_1 = sample.input_1
    input_file_2 = sample.input_2
    output_file_1 = f"{os.path.basename(sample.input_1)}"
    output_file_2 = f"{os.path.basename(sample.input_2)}"
    output_dir = f"{sample.output_dir}/soapnuke/"
    # 定义参数字典
    params = {
        "-1": input_file_1,  # 输入的第一个（正向）读段文件
        "-2": input_file_2,  # 输入的第二个（反向）读段文件
        "-C": output_file_1,  # 输出的第一个清理后的（正向）读段文件
        "-D": output_file_2,  # 输出的第二个清理后的（反向）读段文件
        "-o": output_dir,  # 输出结果的文件夹名称
        "-l": 5,  # 过滤器的最小长度阈值
        "-q": 0.5,  # 过滤器的最小质量阈值（0到1之间）
        "-n": 0.1,  # 允许的最大错误率
        "-T": config.parallel_num,  # 使用的线程数
    }
    # 构造命令字符串
    cmd = dict2cmd("SOAPnuke filter", params)
    return cmd


# 3.序列比对（20小时）
def bismark_alignment(sample, config):
    # 文档地址：https://felixkrueger.github.io/Bismark/options/alignment/
    input_file_1 = f"{sample.output_dir}/{'soapnuke/' if not sample.skip_filter else '/'}{os.path.basename(sample.input_1)}"
    input_file_2 = f"{sample.output_dir}/{'soapnuke/' if not sample.skip_filter else '/'}{os.path.basename(sample.input_2)}"
    temp_dir = f"{sample.output_dir}/bismark_alignment/temp/"
    output_dir = f"{sample.output_dir}/bismark_alignment/"
    # 定义参数字典
    params = {
        "--genome": config.genome_folder,  # 指定参考基因组文件夹
        "-N": 0,  # 允许最多 N（0 或 1）个错配，默认值 0
        "-1": input_file_1,  # 输入的第一个（正向）读段文件
        "-2": input_file_2,  # 输入的第二个（反向）读段文件
        # "-un": "",                           # 保存未比对的读段
        "--bowtie2": "",  # 使用 Bowtie2 作为比对工具
        "--bam": "",  # 输出文件为 BAM 格式
        "--parallel": config.parallel_alignment,  # 线程数，需注意程序会额外占用几个线程，注意容易内存溢出
        "--temp_dir": temp_dir,  # 临时文件目录
        # "--non_directional":"",                # 测序库以非链特异性的方式构建
        "-o": output_dir,  # 指定输出文件夹
    }
    # 构造命令字符串
    cmd = dict2cmd("bismark", params)
    return cmd


# 4.去除重复片段
def bismark_deduplicate(sample, config):
    # 文档地址：https://felixkrueger.github.io/Bismark/options/deduplication/
    input_filename = f"{sample.output_dir}/bismark_alignment/{sample.prefix}_bismark_bt2_pe.bam"
    output_dir = f"{sample.output_dir}/bismark_deduplicate/"
    # output_filename="output/bismark_deduplicate/clean_1_bismark_deduplicate_bt2_pe.bam"
    # 定义参数字典
    params = {
        "-p": "",  # paired-end，即指定为双端
        "--bam": "",  # 输出bam格式
        "--output_dir": output_dir,  # 输出的目录，用于存储去重后的结果及report
        input_filename: "",  # 输入的 BAM 文件所在目录
    }
    # 构造命令字符串
    cmd = dict2cmd("deduplicate_bismark", params)
    return cmd


# 5.提取甲基化信息，并将测序数据的覆盖度转换为细胞碱基甲基化数据
def bismark_methylation_extractor(sample, config):
    # 文档地址：https://felixkrueger.github.io/Bismark/options/methylation_extraction/
    input_file = f"{sample.output_dir}/bismark_deduplicate/{sample.prefix}_bismark_bt2_pe.deduplicated.bam"
    output_dir = f"{sample.output_dir}/bismark_methylation/"
    # 定义参数字典
    params = {
        "--bedGraph": "",  # 生成 bedGraph 文件
        "--CX": "",  # 计算 CpG、CHG 和 CHH 位点的甲基化水平
        "--gzip": "",  # 对输出进行 gzip 压缩
        "--multicore": config.parallel_num
        // 3,  # 实际为3倍线程（1个用于甲基化提取器本身，1个用于Samtools流，1个用于GZIP流）
        "--buffer_size": "30%",  # 设置缓冲区大小为 32GB、总内存的30%等
        "-o": output_dir,  # 指定输出目录
        # 同时输出cytosine_report报告的参数
        "--cytosine_report": "",  # 输出cytosine_report报告
        "--genome_folder": config.genome_folder,  # 参考基因的文件夹，这里必须是绝对路径
        "--split_by_chromosome": "",  # 将输出按染色体分割，每个染色体生成一个文件
        input_file: "",  # 输入文件的路径
    }
    # 构造命令字符串
    cmd = dict2cmd("bismark_methylation_extractor", params)
    return cmd


# 使用自定义脚本1输出基于染色体的甲基化测序深度信息
def methylation_depth_analysis(sample, config):
    input_file = f'"{sample.output_dir}/bismark_methylation/{sample.prefix}_bismark_bt2_pe.deduplicated.CX_report.txt*.gz"'
    output_file = f"{sample.output_dir}/{sample.sample_name}_methylation_depth_report.txt"

    # 定义参数字典
    params = {
        input_file: "",  # 输入文件的路径
        output_file: "",  # 输出文件的路径
    }
    # 使用bash避免权限问题
    cmd = dict2cmd(f"bash {config.utils_folder}/utils/methylation_depth_analysis", params)
    return cmd


# 使用自定义脚本2输出基于染色体和context的甲基化覆盖的统计信息
def methylation_coverage_analyse(sample):
    input_file = f'"{sample.output_dir}/bismark_methylation/{sample.prefix}_bismark_bt2_pe.deduplicated.CX_report.txt*.gz"'
    output_file = f"{sample.output_dir}/{sample.sample_name}_methylation_coverage_report.txt"

    # 定义参数字典
    params = {
        input_file: "",  # 输入文件的路径
        output_file: "",  # 输出文件的路径
    }
    # 使用bash避免权限问题
    cmd = dict2cmd(f"bash {config.utils_folder}/utils/methylation_coverage_analyse", params)
    return cmd


# 使用自定义脚本3输出基于染色体和context的甲基化分布信息（按百分比）
def methylation_distribution_analysis(sample):
    input_file = f'"{sample.output_dir}/bismark_methylation/{sample.prefix}_bismark_bt2_pe.deduplicated.CX_report.txt*.gz"'
    output_file = f"{sample.output_dir}/{sample.sample_name}_methylation_distribution_report.txt"
    # 定义参数字典
    params = {
        input_file: "",  # 输入文件的路径
        output_file: "",  # 输出文件的路径
    }
    # 使用bash避免权限问题
    cmd = dict2cmd(f"bash {config.utils_folder}/utils/methylation_distribution_analysis", params)
    return cmd


# 执行命令，并将结果重定向到log文件
def execute_shell_command(command, log_dir="./log/"):
    # 检查并创建日志目录
    if not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)

    # 从命令的第一个词提取程序名称
    command_args = command.split()
    if command.split()[0].lower() in ["python", "rscript", "bash"]:
        program_name = command_args[0] + "_" + os.path.basename(command_args[1])
    else:
        program_name = os.path.basename(command_args[0])

    # 设置东八区时区（当前进程及子进程有效）
    os.environ["TZ"] = "Asia/Shanghai"
    time.tzset()

    # 获取当前时间作为日志文件名的一部分
    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file_name = f"{log_dir}/{current_time}_{program_name}.log"

    # 定义一个用于处理进程结束时的信号
    def kill_child_processes(signum, frame):
        print("Received signal to terminate, killing all child processes...")
        # 发送信号给整个进程组
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        sys.exit(1)

    # 注册信号处理函数，处理进程终止的信号
    signal.signal(signal.SIGINT, kill_child_processes)  # Ctrl+C
    signal.signal(signal.SIGTERM, kill_child_processes)  # kill 命令

    # 打开日志文件用于写入
    with open(log_file_name, "w") as log_file:
        # 获取当前时间，格式为 HH:MM:SS
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # 先将完整的命令写入日志
        log_file.write(f"[{current_time}] Executing command: {command}\n")
        log_file.flush()
        # 创建一个子进程
        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1,  # 设置行缓冲
            preexec_fn=os.setsid,  # 将子进程放入新的进程组
            executable="/bin/bash",
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
        )

        # 实时读取子进程的输出并写入日志
        for line in process.stdout:
            # 获取当前时间，格式为 HH:MM:SS
            timestamp = datetime.datetime.now().strftime("%H:%M:%S")
            # 写入日志文件，每行以时间开头
            log_file.write(f"[{timestamp}] {line}")
            log_file.flush()  # 立即写入文件

            # 同时也可以选择打印输出到控制台
            print(f"[{timestamp}] {line}", end="")

        # 等待进程结束
        process.wait()

        # 检查退出状态码，非 0 表示执行失败
        if process.returncode != 0:
            raise RuntimeError(f"Command '{command}' failed with return code {process.returncode}")

    # 返回子进程的退出代码
    return process.returncode


# 连接多层文件夹，并自动处理文件夹分隔符
def path_join(parent_path, *child_paths):
    # 去除 parent_path 末尾的斜杠
    parent_path = parent_path.rstrip("/")

    # 去除每个 child_path 两侧的斜杠并连接
    child_paths = [child_path.strip("/") for child_path in child_paths]

    return os.path.join(parent_path, *child_paths)


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


# 默认值
DEFAULTS = {
    # 公共默认参数
    "genome_folder": None,
    "utils_folder": ".",
    "skip_filter": False,
    "parallel_num": 30,
    "parallel_alignment": 4,
    # 样本的默认参数
    "sample_name": None,
    "group_name": None,
    "input_1": "{sample_name}/{sample_name}_1.fq.gz",
    "input_2": "{sample_name}/{sample_name}_2.fq.gz",
    "output_dir": "{sample_dir}/output",
    "log_dir": "{sample_dir}/log",
    "report_dir": "{sample_dir}/report",
}


class DotDict(dict):
    def __getattr__(self, item, default_value=None):
        if item in self:
            return self[item]
        else:
            return default_value

    def __setattr__(self, key, value):
        self[key] = value


# 解析公共参数
def parse_public_config(data):
    if "genome_folder" not in data:
        raise "genome_folder不可缺失"

    config = DotDict()
    config.genome_folder = data.get("genome_folder").rstrip("/")
    config.utils_folder = data.get("utils_folder", DEFAULTS["utils_folder"]).rstrip("/")
    config.skip_filter = data.get("skip_filter", DEFAULTS["skip_filter"])
    config.parallel_num = data.get("parallel_num", DEFAULTS["parallel_num"])
    config.parallel_alignment = data.get("parallel_alignment", DEFAULTS["parallel_alignment"])

    # 将相对路径转为绝对路径
    if not os.path.isabs(config.genome_folder):
        config.genome_folder = os.path.abspath(config.genome_folder)

    if not os.path.exists(config.genome_folder):
        raise FileNotFoundError(f"参考基因组文件夹不存在: {config.genome_folder}")

    if not os.path.exists(config.utils_folder):
        raise FileNotFoundError(f"utils文件夹不存在: {config.utils_folder}")

    return config


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
        sample[key] = sample.get(
            key, DEFAULTS[key].format(sample_dir=os.path.dirname(sample["input_1"]))
        ).rstrip("/")

    return sample


# 自动识别格式并读取samples文件为pandas格式
def read_samples_file(file_path):
    import pandas as pd

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="甲基化分析参数描述")
    # 添加config文件路径参数
    parser.add_argument("--config", type=str, help="添加配置文件（如果设置了该参数，其他参数都不生效）")
    # 添加公共参数
    parser.add_argument("--genome_folder", type=str, help="参考基因组文件所在文件夹（必传）")
    parser.add_argument("--skip_filter", action="store_true", help="添加该参数以跳过数据清洗步骤")
    parser.add_argument("--parallel_num", type=int, default=30, help="最大使用线程数，默认值为30")
    parser.add_argument(
        "--parallel_alignment",
        type=int,
        default=6,
        help="比对使用的线程数，容易内存溢出，默认值为6",
    )
    # 添加样本参数
    parser.add_argument("--sample_name", type=str, help="样本名（必传）")
    parser.add_argument("--group_name", type=str, help="样本所属分组（必传）")
    parser.add_argument(
        "--input_1", type=str, help="测序文件1的路径，默认值为: {sample_name}/{sample_name}_1.fq.gz"
    )
    parser.add_argument(
        "--input_2", type=str, help="测序文件2的路径，默认值为: {sample_name}/{sample_name}_2.fq.gz"
    )
    parser.add_argument(
        "--output_dir", type=str, help="输出的中间文件存放路径，默认值为:{input_1所在文件夹}/output"
    )
    parser.add_argument("--log_dir", type=str, help="日志文件夹，默认值为:{input_1所在文件夹}/log")
    parser.add_argument("--report_dir", type=str, help="报告输出路径，默认值为:{input_1所在文件夹}/report")

    args = parser.parse_args()

    # 如果提供了配置文件，则加载参数
    if args.config:
        # 读取配置文件
        data = jsonload(args.config)
        # 解析公共参数
        config = parse_public_config(data)
        # 解析样本参数
        # if isinstance(data["samples"], str):
        if "samples_file" in data:
            # 传入表格文件，则使用pandas读取表格
            samples_pd = read_samples_file(data["samples_file"])
            samples = [parse_sample_config(sample) for i, sample in samples_pd.iterrows()]
        else:
            # 遍历json解析样本参数
            samples = [parse_sample_config(sample) for sample in data["samples"]]
    else:  # 否则从args读取参数
        # 解析公共参数
        config = parse_public_config(args)
        # 解析样本参数
        samples = [parse_sample_config(args)]

    # 创建参考基因组的索引文件
    if not os.path.exists(config.genome_folder + "/Bisulfite_Genome/"):
        cmd = bismark_genome_preparation(config)
        print("创建参考基因组的索引文件: ", cmd)
        execute_shell_command(cmd, samples[0].log_dir)
    else:
        print("检测到参考基因组的索引文件已存在，跳过索引构建")

    for sample in samples:
        print(f"开始处理样本{sample.sample_name}...")
        # 创建输出目录
        cmd = mkdirs(sample, config)
        print("-----------------------")
        print("创建输出目录: ", cmd)
        execute_shell_command(cmd, sample.log_dir)

        # 使用SOAPnuke做数据过滤
        if not config.skip_filter:
            cmd = soapnuke_filter(sample, config)
            print("-----------------------")
            print("使用SOAPnuke做数据过滤: ", cmd)
            execute_shell_command(cmd, sample.log_dir)

        # 序列比对（12~16小时）
        cmd = bismark_alignment(sample, config)
        print("-----------------------")
        print("序列比对: ", cmd)
        execute_shell_command(cmd, sample.log_dir)

        # 去除重复片段（5小时）
        cmd = bismark_deduplicate(sample, config)
        print("-----------------------")
        print("去除重复片段: ", cmd)
        execute_shell_command(cmd, sample.log_dir)

        # 提取甲基化信息，并将测序数据的覆盖度转换为细胞碱基甲基化数据（20小时）
        cmd = bismark_methylation_extractor(sample, config)
        print("-----------------------")
        print("提取甲基化信息: ", cmd)
        execute_shell_command(cmd, sample.log_dir)

        # 使用自定义脚本1输出基于染色体的甲基化测序深度信息（10分钟）
        cmd = methylation_depth_analysis(sample)
        print("-----------------------")
        print("输出甲基化测序深度信息: ", cmd)
        execute_shell_command(cmd, sample.log_dir)

        # 使用自定义脚本2输出基于染色体和context的甲基化覆盖的统计信息（10分钟）
        cmd = methylation_coverage_analyse(sample)
        print("-----------------------")
        print("输出甲基化覆盖度信息: ", cmd)
        execute_shell_command(cmd, sample.log_dir)

        # 使用自定义脚本3输出基于染色体和context的甲基化分布信息（按百分比）（10分钟）
        cmd = methylation_distribution_analysis(sample)
        print("-----------------------")
        print("输出甲基化分布信息: ", cmd)
        execute_shell_command(cmd, sample.log_dir)

        print(f"样本{sample.sample_name}处理完成")
        print("=======================")
