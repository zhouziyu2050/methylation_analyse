import subprocess
import datetime
import os
import signal
import sys
import time

# 执行命令，并将结果重定向到log文件
def execute_shell_command(command,log_dir="./log/"):
    # 检查并创建日志目录
    if not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)

    # 从命令的第一个词提取程序名称
    program_name = command.split()[0]
    # 如果 program_name 是路径，则提取文件名
    if os.path.isfile(program_name):
        program_name = os.path.basename(program_name)

    # 设置东八区时区（当前进程及子进程有效）
    os.environ['TZ'] = 'Asia/Shanghai'
    time.tzset()  # 更新系统时区设置

    # 获取当前时间作为日志文件名的一部分
    current_time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    log_file_name = f"{log_dir}{current_time}_{program_name}.log"
    
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
    with open(log_file_name, 'w') as log_file:
        # 获取当前时间，格式为 HH:MM:SS
        current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        # 先将完整的命令写入日志
        log_file.write(f"[{current_time}] Executing command: {command}\n")
        log_file.flush()  # 立即写入文件
        # 创建一个子进程，设置进程组 ID (os.setsid 将子进程放入一个新的进程组)
        process = subprocess.Popen(
            command, 
            shell=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            preexec_fn=os.setsid  # 将子进程放入新的进程组
        )
        
        # 实时读取子进程的输出并写入日志
        for line in process.stdout:
            # 获取当前时间，格式为 HH:MM:SS
            timestamp = datetime.datetime.now().strftime('%H:%M:%S')
            # 写入日志文件，每行以时间开头
            log_file.write(f"[{timestamp}] {line}")
            log_file.flush()  # 立即写入文件
            
            # 同时也可以选择打印输出到控制台
            # print(f"[{timestamp}] {line}", end="")
        
        # 等待进程结束
        process.wait()

        # 检查退出状态码，非 0 表示执行失败
        if process.returncode != 0:
            raise RuntimeError(f"Command '{command}' failed with return code {process.returncode}")

    # 返回子进程的退出代码
    return process.returncode

# 用于将字典参数构造成命令字符串
def dict2cmd(prefix,params):
    cmd = prefix
    for param, value in params.items():
        if value==None or value=="":
            cmd += f" {param}"
        else:
            cmd += f" {param} {value}"
    # print("执行命令：",cmd)
    return cmd

# 创建输出目录
def mkdirs():
    directories = [
        "output/bismark_alignment",
        "output/bismark_alignment/temp",
        "output/bismark_methylation",
        "output/bismark_deduplicate",
        # "report"
    ]
    if genome_filter:
        directories.append("output/soapnuke")
    directories=[f"{root_dir}{sample_name}/{x}" for x in directories]
    directories.append(f"{root_dir}report/{sample_name}")

    params = {
        "-p":"",                       # 创建多级目录
        "-m":777,                      # 目录权限
        " ".join(directories):"",      # 目录列表
    }
    cmd = dict2cmd("mkdir",params)
    return cmd

# 1.创建参考基因组的索引文件（每个参考基因组只需要执行一次）
def bismark_genome_preparation():
    # 文档地址：https://felixkrueger.github.io/Bismark/options/genome_preparation/
    
    # 定义参数字典
    params = {
        "--bowtie2": "",          # 使用 Bowtie2 作为比对工具
        "--parallel": parallel_num//2,   # 每个线程实际使用2个核心，所以这里减半
        genome_folder: ""         # 参考基因组文件夹
    }
    # 构造命令字符串
    cmd = dict2cmd("bismark_genome_preparation",params)
    return cmd

# 2.使用SOAPnuke做数据过滤（已过滤则不需要这一步）
def soapnuke_filter():
    # 文档地址：https://github.com/BGI-flexlab/SOAPnuke/blob/master/Readme.md
    input_file_1=f"{root_dir}{sample_name}/{sample_name}_1.fq.gz"
    input_file_2=f"{root_dir}{sample_name}/{sample_name}_2.fq.gz"
    output_file_1=f"{root_dir}{sample_name}/output/soapnuke/{sample_name}_1.fq.gz"
    output_file_2=f"{root_dir}{sample_name}/output/soapnuke/{sample_name}_2.fq.gz"
    output_dir=f"{root_dir}{sample_name}/output/soapnuke/"
    # 定义参数字典
    params = {
        "-1": input_file_1,         # 输入的第一个（正向）读段文件
        "-2": input_file_2,         # 输入的第二个（反向）读段文件
        "-C": output_file_1,        # 输出的第一个清理后的（正向）读段文件
        "-D": output_file_2,        # 输出的第二个清理后的（反向）读段文件
        "-o": output_dir,           # 输出结果的文件夹名称
        "-l": 5,                    # 过滤器的最小长度阈值
        "-q": 0.5,                  # 过滤器的最小质量阈值（0到1之间）
        "-n": 0.1,                  # 允许的最大错误率
        "-T": parallel_num          # 使用的线程数
    }
    # 构造命令字符串
    cmd = dict2cmd("SOAPnuke filter",params)
    return cmd

# 3.序列比对（20小时）
def bismark_alignment():
    # 文档地址：https://felixkrueger.github.io/Bismark/options/alignment/
    input_file_1=f"{root_dir}{sample_name}/{'output/soapnuke/' if genome_filter else ''}{sample_name}_1.fq.gz"
    input_file_2=f"{root_dir}{sample_name}/{'output/soapnuke/' if genome_filter else ''}{sample_name}_2.fq.gz"
    temp_dir=f"{root_dir}{sample_name}/output/bismark_alignment/temp/"
    output_dir=f"{root_dir}{sample_name}/output/bismark_alignment/"
    # 定义参数字典
    params = {
        "--genome": genome_folder,             # 指定参考基因组文件夹
        "-N": 0,                               # 允许最多 N（0 或 1）个错配，默认值 0
        "-1": input_file_1,                    # 输入的第一个（正向）读段文件
        "-2": input_file_2,                    # 输入的第二个（反向）读段文件
        # "-un": "",                           # 保存未比对的读段
        "--bowtie2": "",                       # 使用 Bowtie2 作为比对工具
        "--bam": "",                           # 输出文件为 BAM 格式
        "--parallel": parallel_alignment,      # 线程数，需注意程序会额外占用几个线程，注意容易内存溢出
        "--temp_dir": temp_dir,                # 临时文件目录
        # "--non_directional":"",                # 测序库以非链特异性的方式构建
        "-o": output_dir,                      # 指定输出文件夹
    }
    # 构造命令字符串
    cmd = dict2cmd("bismark",params)
    return cmd

# 4.去除重复片段
def bismark_deduplicate():
    # 文档地址：https://felixkrueger.github.io/Bismark/options/deduplication/
    input_filename=f"{root_dir}{sample_name}/output/bismark_alignment/{sample_name}_1_bismark_bt2_pe.bam"
    output_dir=f"{root_dir}{sample_name}/output/bismark_deduplicate/"
    # output_filename="output/bismark_deduplicate/clean_1_bismark_deduplicate_bt2_pe.bam"
    # 定义参数字典
    params = {
        "-p": "",                  # paired-end，即指定为双端
        "--bam":"",                # 输出bam格式
        "--output_dir":output_dir, # 输出的目录，用于存储去重后的结果及report
        input_filename:"",         # 输入的 BAM 文件所在目录
    }
    # 构造命令字符串
    cmd = dict2cmd("deduplicate_bismark",params)
    return cmd

# 5.提取甲基化信息，并将测序数据的覆盖度转换为细胞碱基甲基化数据
def bismark_methylation_extractor():
    # 文档地址：https://felixkrueger.github.io/Bismark/options/methylation_extraction/
    input_file = f"{root_dir}{sample_name}/output/bismark_deduplicate/{sample_name}_1_bismark_bt2_pe.deduplicated.bam"
    output_dir = f"{root_dir}{sample_name}/output/bismark_methylation/"
    # 定义参数字典
    params = {
        "--bedGraph": "",                     # 生成 bedGraph 文件
        "--CX": "",                           # 计算 CpG、CHG 和 CHH 位点的甲基化水平
        "--gzip": "",                         # 对输出进行 gzip 压缩
        "--multicore":parallel_num//3,        # 实际为3倍线程（1个用于甲基化提取器本身，1个用于Samtools流，1个用于GZIP流）
        "--buffer_size": "30%",               # 设置缓冲区大小为 32GB、总内存的30%等
        "-o": output_dir,                     # 指定输出目录
        
        # 同时输出cytosine_report报告的参数
        "--cytosine_report": "",              # 输出cytosine_report报告
        "--genome_folder": genome_folder,     # 参考基因的文件夹，这里必须是绝对路径
        '--split_by_chromosome': '',          # 将输出按染色体分割，每个染色体生成一个文件
        
        input_file:""                         # 输入文件的路径
    }
    # 构造命令字符串
    cmd = dict2cmd("bismark_methylation_extractor",params)
    return cmd

# # 单独执行coverage2cytosine
# def coverage2cytosine():
#     input_file = f"{root_dir}{sample_name}/output/bismark_methylation/{sample_name}_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
#     output_dir = f"{root_dir}{sample_name}/output/bismark_coverage2cytosine/{sample_name}_1_bismark_bt2_pe.deduplicated.CX_report.txt"
#     # 定义参数字典
#     params = {
#         "--CX": "",                           # 计算 CpG、CHG 和 CHH 位点的甲基化水平
#         "--gzip": "",                         # 对输出进行 gzip 压缩
#         "-o": output_dir,                     # 指定输出目录
#         "--genome_folder": genome_folder,     # 参考基因的文件夹，这里必须是绝对路径
#         '--split_by_chromosome': '',          # 将输出按染色体分割，每个染色体生成一个文件
        
#         input_file:""                         # 输入文件的路径
#     }
#     # 构造命令字符串
#     cmd = dict2cmd("coverage2cytosine",params)
#     return cmd

# 使用自定义脚本1输出基于染色体的甲基化测序深度信息
def methylation_depth_analysis():
    input_file = f"\"{root_dir}{sample_name}/output/bismark_methylation/{sample_name}_1_bismark_bt2_pe.deduplicated.CX_report.txt*.gz\""
    output_file = f"{root_dir}{sample_name}/output/{sample_name}_methylation_depth_report.txt"
    # 定义参数字典
    params = {
        input_file:"",                         # 输入文件的路径
        output_file:"",                        # 输出文件的路径
    }
    # 构造命令字符串
    cmd = dict2cmd("/methylation/utils/methylation_depth_analysis",params)
    return cmd

# 使用自定义脚本2输出基于染色体和context的甲基化覆盖的统计信息
def methylation_coverage_analyse():
    input_file = f"\"{root_dir}{sample_name}/output/bismark_methylation/{sample_name}_1_bismark_bt2_pe.deduplicated.CX_report.txt*.gz\""
    output_file = f"{root_dir}{sample_name}/output/{sample_name}_methylation_coverage_report.txt"
    # 定义参数字典
    params = {
        input_file:"",                         # 输入文件的路径
        output_file:"",                        # 输出文件的路径
    }
    # 构造命令字符串
    cmd = dict2cmd("/methylation/utils/methylation_coverage_analyse",params)
    return cmd

# 使用自定义脚本3输出基于染色体和context的甲基化分布信息（按百分比）
def methylation_distribution_analysis():
    input_file = f"\"{root_dir}{sample_name}/output/bismark_methylation/{sample_name}_1_bismark_bt2_pe.deduplicated.CX_report.txt*.gz\""
    output_file = f"{root_dir}{sample_name}/output/{sample_name}_methylation_distribution_report.txt"
    # 定义参数字典
    params = {
        input_file:"",                         # 输入文件的路径
        output_file:"",                        # 输出文件的路径
    }
    # 构造命令字符串
    cmd = dict2cmd("/methylation/utils/methylation_distribution_analysis",params)
    return cmd


if __name__=="__main__":
    # 全局参数设置
    parallel_num=30                                 # 最大线程数
    parallel_alignment=6                            # 对齐比对的线程数，需注意程序会额外占用几个线程，注意容易内存溢出
    genome_folder="/methylation/genome/mm39/"       # 文件夹中需包含xxx.fa或xxx.fa.gz
    genome_filter=False                             # 是否需要清洗数据
    root_dir="/methylation/F24A080000424_MUSekgzH_20240805100100/"  # 基因文件的根目录
    sample_names=["13A","32A","50A","26A","28A","38A","74A","25A","75A","88A","97A"] # 样本名称
    log_dir = "./log/"                              # 日志输出文件夹

    # 1.创建参考基因组的索引文件
    if not os.path.exists(genome_folder+"Bisulfite_Genome/"):
        cmd = bismark_genome_preparation()
        print("创建参考基因组的索引文件: ", cmd)
        execute_shell_command(cmd, log_dir)
    else:
        print("检测到参考基因组的索引文件已存在，跳过索引构建")

    for sample_name in sample_names:
        print(f"开始处理: {sample_name}...")

        # 将日志放到样本目录中
        log_dir = f"{root_dir}{sample_name}/log/"

        # 创建输出目录
        cmd = mkdirs()
        print("-----------------------")
        print("创建输出目录: ", cmd)
        execute_shell_command(cmd, log_dir)

        # 使用SOAPnuke做数据过滤
        if genome_filter:
            cmd = soapnuke_filter()
            print("-----------------------")
            print("使用SOAPnuke做数据过滤: ", cmd)
            execute_shell_command(cmd, log_dir)

        # 序列比对（12~16小时）
        cmd = bismark_alignment()
        print("-----------------------")
        print("序列比对: ", cmd)
        execute_shell_command(cmd, log_dir)

        # 去除重复片段（5小时）
        cmd = bismark_deduplicate()
        print("-----------------------")
        print("去除重复片段: ", cmd)
        execute_shell_command(cmd, log_dir)

        # 提取甲基化信息，并将测序数据的覆盖度转换为细胞碱基甲基化数据（20小时）
        cmd = bismark_methylation_extractor()
        print("-----------------------")
        print("提取甲基化信息: ", cmd)
        execute_shell_command(cmd, log_dir)

        # 使用自定义脚本1输出基于染色体的甲基化测序深度信息（10分钟）
        cmd = methylation_depth_analysis()
        print("-----------------------")
        print("输出甲基化测序深度信息: ", cmd)
        execute_shell_command(cmd, log_dir)

        # 使用自定义脚本2输出基于染色体和context的甲基化覆盖的统计信息（10分钟）
        cmd = methylation_coverage_analyse()
        print("-----------------------")
        print("输出甲基化覆盖度信息: ", cmd)
        execute_shell_command(cmd, log_dir)

        # 使用自定义脚本3输出基于染色体和context的甲基化分布信息（按百分比）（10分钟）
        cmd = methylation_distribution_analysis()
        print("-----------------------")
        print("输出甲基化分布信息: ", cmd)
        execute_shell_command(cmd, log_dir)
        

        print(f"{sample_name}处理完成")
        print("=======================")
