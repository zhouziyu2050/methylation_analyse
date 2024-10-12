##################################################################
# 加载读取参数的包及必要的function
##################################################################

# 加载读取数据的包
library(optparse)
library(jsonlite)
library(readr) # 用于读取CSV/TSV文件
library(readxl) # 用于读取Excel文件

# 读取JSON文件并移除注释
read_json <- function(config_file) {
  json_string <- readLines(config_file, warn = FALSE)
  json_string <- paste(json_string, collapse = "\n")

  # 去除单行注释 // 之后的内容
  json_string <- gsub("//.*", "", json_string, perl = TRUE)
  # 移除多行注释 (/* ... */)，非贪婪模式匹配
  json_string <- gsub("/\\*[\\s\\S]*?\\*/", "", json_string, perl = TRUE)
  # 删除 `}` 或 `]` 前面的逗号
  json_string <- gsub(",(\\s*[\\]}])", "\\1", json_string, perl = TRUE)

  # 将处理后的字符串解析为JSON
  return(fromJSON(json_string))
}

# 读取表格数据
read_data <- function(file_path) {
  # 获取文件扩展名
  file_ext <- tools::file_ext(file_path)
  # 根据文件扩展名选择读取方法
  if (file_ext == "csv") {
    data <- read_csv(file_path, show_col_types = FALSE) # 使用 readr 包读取 CSV
  } else if (file_ext == "tsv") {
    data <- read_tsv(file_path, show_col_types = FALSE) # 使用 readr 包读取 TSV
  } else if (file_ext == "xlsx" || file_ext == "xls") {
    data <- read_excel(file_path) # 使用 readxl 包读取 Excel 文件
  } else {
    data <- read_tsv(file_path) # 默认使用tsv格式
  }
  return(as.data.frame(data))
}

##################################################################
# 设置和读取参数（如果使用R交互式命令窗口则不运行这部分）
##################################################################

# 定义命令行参数
option_list <- list(
  make_option(c("-c", "--config"),
    type = "character", default = NULL,
    help = "配置文件路径", metavar = "character"
  ),
  make_option(c("-o", "--output_dir"),
    type = "character", default = NULL,
    help = "中间文件的输出文件夹，默认值为：{当前文件夹}/output", metavar = "character"
  ),
  make_option(c("-r", "--report_dir"),
    type = "character", default = NULL,
    help = "报告的输出文件夹，默认值为：{当前文件夹}/report", metavar = "character"
  ),
  make_option(c("-a", "--group_a"),
    type = "character", default = NULL,
    help = "DMR的组A名称", metavar = "character"
  ),
  make_option(c("-b", "--group_b"),
    type = "character", default = NULL,
    help = "DMR的组B名称", metavar = "character"
  ),
  make_option(c("-f", "--samples_file"),
    type = "character", default = NULL,
    help = "以tsv/csv/excel文件传入样本参数", metavar = "character"
  ),
  make_option(c("-g", "--gtf_file"),
    type = "character", default = NULL,
    help = "gtf注释文件路径（下载地址：https://www.gencodegenes.org/）",
    metavar = "character"
  )
)
# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 创建一个空的config向量
config <- list()

# 如果有config参数则从配置文件读取参数
if (!is.null(opt$config)) {
  # 读取配置文件
  config <- read_json(opt$config)
}

# 如果命令行传入了某参数，则更新config中的对应参数
for (k in names(opt)) {
  if (!is.null(opt[[k]]) && opt[[k]] != "") {
    config[[k]] <- opt[[k]]
  }
}

##################################################################
# 校验和调整参数
##################################################################

# 若使用R交互式命令行，则运行以下代码加载config文件
# config <- read_json("config.json")

# 确保必填参数不为空
if (is.null(config$group_a) || is.null(config$group_b)) {
  print_help(opt_parser)
  stop("缺少必填参数：group_a 和/或 group_b")
}
if (is.null(config$samples_file) && is.null(config$samples)) {
  print_help(opt_parser)
  stop("samples_file参数和配置文件中的samples参数不能同时为空")
}

# 如果提供了samples_file参数，则从文件读取samples数据，否则从config$samples读取
if (!is.null(config$samples_file)) {
  samples <- read_data(config$samples_file)
} else {
  samples <- as.data.frame(config$samples)
}

# 设置默认输出文件夹
if (is.null(config$output_dir)) {
  config$output_dir <- "./output"
}
if (is.null(config$report_dir)) {
  config$report_dir <- "./report"
}

# 检查group_a, group_b是否存在
if (!(config$group_a %in% samples$group_name)) {
  stop(paste0("组A (", config$group_a, ") 在 samples 数据中不存在"))
}
if (!(config$group_b %in% samples$group_name)) {
  stop(paste0("组B (", config$group_b, ") 在 samples 数据中不存在"))
}

# 按传入的两个分组过滤数据，并重置行号
samples <- samples[samples$group_name %in% c(config$group_a, config$group_b), ]
row.names(samples) <- NULL

# 检查sample_name和group_names是否存在空值
if (any(is.na(samples$sample_name)) || any(nchar(samples$sample_name) == 0)) {
  stop("所有样本的sample_name不能为空")
}
if (any(is.na(samples$group_name)) || any(nchar(samples$group_name) == 0)) {
  stop("所有样本的group_name不能为空")
}

# 设置样本输出文件的默认路径和prefix（用于寻找CX_report文件）
samples$output_dir <- ifelse(
  is.na(samples$output_dir) | samples$output_dir == "",
  paste0(samples$sample_name, "/output"),
  samples$output_dir
)
samples$prefix <- ifelse(!is.na(samples$input_1) & samples$input_1 != "",
  sub("\\..*$", "", basename(samples$input_1)), # 提取文件名并去掉后缀
  paste0(samples$sample_name, "_1") # 使用 sample_name 和 _1
)

# 读取sample_names和group_names
sample_names <- samples$sample_name
group_names <- samples$group_name


# 打印以确认参数
cat("中间文件输出文件夹:", config$output_dir, "\n")
cat("报告输出文件夹:", config$report_dir, "\n")
cat("组A名称:", config$group_a, "\n")
cat("组B名称:", config$group_b, "\n")
cat("样本名称:", paste(samples$sample_name, collapse = ", "), "\n")
cat("样本分组:", paste(samples$group_name, collapse = ", "), "\n")

##################################################################
# 补充参数及文件夹检查
##################################################################

# 需要处理的染色体名称
seqnames <- c(
  "NC_000067.7", "NC_000068.8", "NC_000069.7", "NC_000070.7", "NC_000071.7",
  "NC_000072.7", "NC_000073.7", "NC_000074.7", "NC_000075.7", "NC_000076.7",
  "NC_000077.7", "NC_000078.7", "NC_000079.7", "NC_000080.7", "NC_000081.7",
  "NC_000082.7", "NC_000083.7", "NC_000084.7", "NC_000085.7", "NC_000086.8",
  "NC_000087.8", "NC_005089.1"
)

# accession和chromosome的映射关系
accession2chromosome <- c(
  "NC_000067.7" = "chr1", "NC_000068.8" = "chr2", "NC_000069.7" = "chr3",
  "NC_000070.7" = "chr4", "NC_000071.7" = "chr5", "NC_000072.7" = "chr6",
  "NC_000073.7" = "chr7", "NC_000074.7" = "chr8", "NC_000075.7" = "chr9",
  "NC_000076.7" = "chr10", "NC_000077.7" = "chr11", "NC_000078.7" = "chr12",
  "NC_000079.7" = "chr13", "NC_000080.7" = "chr14", "NC_000081.7" = "chr15",
  "NC_000082.7" = "chr16", "NC_000083.7" = "chr17", "NC_000084.7" = "chr18",
  "NC_000085.7" = "chr19", "NC_000086.8" = "chrX", "NC_000087.8" = "chrY",
  "NC_005089.1" = "chrM"
)


# 给输出文件夹追加分组信息
config$output_dir <- paste0(
  sub("/$", "", config$output_dir), "/",
  config$group_a, "_and_", config$group_b
)
config$report_dir <- paste0(
  sub("/$", "", config$report_dir), "/",
  config$group_a, "_and_", config$group_b
)

# 检查并创建对应的输出文件夹
if (!dir.exists(config$output_dir)) {
  dir.create(config$output_dir, recursive = TRUE)
}
if (!dir.exists(config$report_dir)) {
  dir.create(config$report_dir, recursive = TRUE)
}


##################################################################
# 加载库
##################################################################

# library(DMRcaller)
# library(betareg)
# library(tibble)
# library(data.table)

# 静默加载
suppressPackageStartupMessages(library(DMRcaller))
suppressPackageStartupMessages(library(betareg))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(data.table))


##################################################################
# DMR区域计算
##################################################################


DMRsReplicatesBinsList <- list()
for (seqname in seqnames) {
  # seqname <- "NC_000085.7"
  cat("开始处理染色体：", seqname, "\n")
  # Initialize an empty list to store the data
  methylationDataList <- list()
  # Loop through each sample name and read the corresponding file
  for (i in seq_len(nrow(samples))) {
    file_path <- paste0(
      sub("/$", "", samples[i, "output_dir"]), "/bismark_methylation/",
      samples[i, "prefix"], "_bismark_bt2_pe.deduplicated.CX_report.txt.chr",
      seqname, ".CX_report.txt.gz"
    )
    methylationDataList[[i]] <- readBismark(file_path)
  }

  # Initialize methylationData with the first element of methylationDataList
  methylationData <- methylationDataList[[1]]

  # Loop through the rest of the list and merge each element into methylationData
  for (i in 2:length(methylationDataList)) {
    methylationData <- joinReplicates(methylationData, methylationDataList[[i]])
  }

  cat("DMR分析...\n")
  # 使用 computeDMRsReplicates 函数计算差异甲基化区域 (DMRs)
  tryCatch(
    {
      # CG
      DMRsReplicatesBins <- computeDMRsReplicates(
        methylationData = methylationData, # 输入的甲基化数据
        condition = group_names, # 样本的条件（分组信息）
        regions = NULL, # 目标区域，如果为 NULL 则使用默认的全基因组区域
        context = c("CG"), # 分析的上下文类型（例如 CG）
        method = "bins", # 用于检测 DMR 的方法，这里使用 "bins" 方法
        binSize = 200, # 分箱大小，单位为碱基对，这里设定为 200bp
        test = "betareg", # 统计测试方法，这里选择了贝塔回归 (betareg)
        pValueThreshold = 0.05, # 设定的 p 值阈值，用于判断显著性，默认值为 0.05
        minCytosinesCount = 5, # 分析时要求每个区域中至少有 5 个胞嘧啶
        minReadsPerCytosine = 5, # 每个胞嘧啶至少需要 5 个读数（reads）
        minProportionDifference = 0.2, # 最小甲基化比例差异，只有大于此差异的区域才会被视为 DMR。CG为0.2，CHG为0.1，CHH类型为0.05
        minGap = 0, # 区域之间的最小间隔，设置为 0 表示没有间隔
        cores = 1, # 使用的 CPU 核心数
        # minSize = 5,                       # 最小区域大小，单位为碱基对，这里设定为 50bp
        # pseudocountM = 1,                # 添加到甲基化读数 (M count) 的伪计数，默认值为 0
        # pseudocountN = 2,                # 添加到非甲基化读数 (N count) 的伪计数，默认值为 0
      )
      if (length(DMRsReplicatesBins) > 0) {
        DMRsReplicatesBinsList[[paste0(seqname, "CG")]] <- DMRsReplicatesBins
      }
    },
    error = function(e) {
      message(seqname, "Error encountered: ", e$message)
    },
    warning = function(w) {
      message(seqname, "Warning: ", w$message)
    }
  )
  tryCatch(
    {
      # CHG
      DMRsReplicatesBins <- computeDMRsReplicates(
        methylationData = methylationData, # 输入的甲基化数据
        condition = group_names, # 样本的条件（分组信息）
        regions = NULL, # 目标区域，如果为 NULL 则使用默认的全基因组区域
        context = c("CHG"), # 分析的上下文类型（例如 CG）
        method = "bins", # 用于检测 DMR 的方法，这里使用 "bins" 方法
        binSize = 200, # 分箱大小，单位为碱基对，这里设定为 200bp
        test = "betareg", # 统计测试方法，这里选择了贝塔回归 (betareg)
        pValueThreshold = 0.05, # 设定的 p 值阈值，用于判断显著性，默认值为 0.05
        minCytosinesCount = 5, # 分析时要求每个区域中至少有 5 个胞嘧啶
        minReadsPerCytosine = 5, # 每个胞嘧啶至少需要 5 个读数（reads）
        minProportionDifference = 0.1, # 最小甲基化比例差异，只有大于此差异的区域才会被视为 DMR。CG为0.2，CHG为0.1，CHH类型为0.05
        minGap = 0, # 区域之间的最小间隔，设置为 0 表示没有间隔
        cores = 1, # 使用的 CPU 核心数
        # minSize = 5,                       # 最小区域大小，单位为碱基对，这里设定为 50bp
        # pseudocountM = 1,                # 添加到甲基化读数 (M count) 的伪计数，默认值为 0
        # pseudocountN = 2,                # 添加到非甲基化读数 (N count) 的伪计数，默认值为 0
      )
      if (length(DMRsReplicatesBins) > 0) {
        DMRsReplicatesBinsList[[paste0(seqname, "CHG")]] <- DMRsReplicatesBins
      }
    },
    error = function(e) {
      message(seqname, "Error encountered: ", e$message)
    },
    warning = function(w) {
      message(seqname, "Warning: ", w$message)
    }
  )

  tryCatch(
    {
      # CHH
      DMRsReplicatesBins <- computeDMRsReplicates(
        methylationData = methylationData, # 输入的甲基化数据
        condition = group_names, # 样本的条件（分组信息）
        regions = NULL, # 目标区域，如果为 NULL 则使用默认的全基因组区域
        context = c("CHH"), # 分析的上下文类型（例如 CG）
        method = "bins", # 用于检测 DMR 的方法，这里使用 "bins" 方法
        binSize = 200, # 分箱大小，单位为碱基对，这里设定为 200bp
        test = "betareg", # 统计测试方法，这里选择了贝塔回归 (betareg)
        pValueThreshold = 0.05, # 设定的 p 值阈值，用于判断显著性，默认值为 0.05
        minCytosinesCount = 5, # 分析时要求每个区域中至少有 5 个胞嘧啶
        minReadsPerCytosine = 5, # 每个胞嘧啶至少需要 5 个读数（reads）
        minProportionDifference = 0.05, # 最小甲基化比例差异，CG为0.2，CHG为0.1，CHH为0.05
        minGap = 0, # 区域之间的最小间隔，设置为 0 表示没有间隔
        cores = 1, # 使用的 CPU 核心数
        # minSize = 5,                       # 最小区域大小，单位为碱基对，这里设定为 50bp
        # pseudocountM = 1,                # 添加到甲基化读数 (M count) 的伪计数，默认值为 0
        # pseudocountN = 2,                # 添加到非甲基化读数 (N count) 的伪计数，默认值为 0
      )
      if (length(DMRsReplicatesBins) > 0) {
        DMRsReplicatesBinsList[[paste0(seqname, "CHH")]] <- DMRsReplicatesBins
      }
    },
    error = function(e) {
      message(seqname, "Error encountered: ", e$message)
    },
    warning = function(w) {
      message(seqname, "Warning: ", w$message)
    }
  )
}
cat("DMR结果合并导出...\n")
# 清除DMRsReplicatesBinss列表元素的名称（不清除名称无法合并）
names(DMRsReplicatesBinsList) <- NULL
# 使用 c() 函数将多个 GRanges 对象拼接成一个，并转为frame
DMRsReplicatesBinsCombined <- as.data.frame(do.call(c, DMRsReplicatesBinsList))
# 将accession转为chromosome
DMRsReplicatesBinsCombined$seqnames <- accession2chromosome[DMRsReplicatesBinsCombined$seqnames]
# 将 DMRs 导出为文本文件
# 为了方便寻找重叠区域未输出表头，其表头为[seqnames start end width strand sumReadsM1 sumReadsN1 proportion1 sumReadsM2 sumReadsN2 proportion2 cytosinesCount context direction pValue regionType]
write.table(DMRsReplicatesBinsCombined, file = paste0(config$output_dir, "/DMRsReplicatesBins.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 统计DMR信息（每条染色体、不同context所包含DMR的数量）
cat("输出DMR统计信息...\n")
# 生成频数表
DMRsummary <- table(DMRsReplicatesBinsCombined$seqnames, DMRsReplicatesBinsCombined$context)
# 将频数表转换为数据框
DMRsummary <- as.data.frame.matrix(DMRsummary)
# 将Seqnames转换为普通列（rownames_to_column函数依赖于tibble包）
DMRsummary <- rownames_to_column(DMRsummary, var = "Seqnames")
# 输出到文件
write.table(DMRsummary,
  file = paste0(config$report_dir, "DMR_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)


##################################################################
# DMR注释
##################################################################

cat("DMR注释...\n")

gtf_bed_file <- paste0(config$gtf_file, ".bed")
# 判断并生成对应的bed文件
if (!file.exists(gtf_bed_file)) {
  if (grepl("\\.gz$", config$gtf_file)) {
    # gz格式，则先解压再使用gtf2bed转为bed格式
    command <- paste0("zcat ", config$gtf_file, " | gtf2bed > ", gtf_bed_file)
  } else {
    # 不是gz格式，则直接使用gtf2bed转为bed格式
    command <- paste0("gtf2bed < ", config$gtf_file, " > ", gtf_bed_file)
  }
  system(command)
}

# 寻找重叠区域并注释
command <- paste0(
  "bedtools intersect -a ", config$output_dir, "/DMRsReplicatesBins.txt",
  " -b ", gtf_bed_file,
  " -wa -wb > ", config$output_dir, "/DMR_gene_association.bed"
)
system(command)


##################################################################
# 统计DMR基因信息
##################################################################

cat("整理DMR注释...\n")

# 注释解析函数
parse_key_value_pairs <- function(text) {
  # 将字符串按分号分割为键值对
  pairs <- unlist(strsplit(text, ";"))

  # 去除空格
  pairs <- trimws(pairs)

  # 初始化一个空的列表用于存储键值对
  key_value_list <- list()

  # 遍历每个键值对
  for (pair in pairs) {
    # 将键和值按空格分开
    parts <- unlist(strsplit(pair, " ", 2))

    if (length(parts) == 2) {
      # 提取键和值
      key <- trimws(parts[1])
      value <- trimws(parts[2])

      # 去除引号
      value <- gsub('"', "", value)

      # 将键值对添加到列表
      key_value_list[[key]] <- value
    }
  }

  return(key_value_list)
}

# 将多行的 extra 列解析为数据框
parse_extra_column <- function(original_df) {
  # 解析 extra 列
  parsed_list <- lapply(original_df$extra, parse_key_value_pairs)

  # 获取所有可能的列名
  all_colnames <- unique(unlist(lapply(parsed_list, names)))

  # 将列表转换为数据框，并填充缺失的列
  df_extra <- rbindlist(lapply(parsed_list, function(x) {
    # 创建一个具有所有可能列名的空行
    row <- setNames(as.list(rep(NA, length(all_colnames))), all_colnames)
    # 填充现有列的值
    row[names(x)] <- unlist(x)
    return(row)
  }), fill = TRUE)

  # 删除 extra 列
  original_df <- original_df[, !names(original_df) %in% "extra"]

  # 合并到原数据框
  combined_df <- cbind(original_df, df_extra)
  return(combined_df)
}

# 获取去重结果
filtered_genes <- function(df) {
  # 使用 strsplit 按点分隔 gene_id 列
  gene_split <- strsplit(df$gene_id, "\\.")

  # 将分隔结果转换为两列
  gene_split_df <- do.call(rbind, lapply(gene_split, function(x) {
    if (length(x) == 2) {
      return(c(x[1], x[2])) # gene_id 和 version
    } else {
      return(c(x[1], NA)) # 如果没有 version，则填充 NA
    }
  }))

  # 设置列名
  colnames(gene_split_df) <- c("gene_id", "version")
  # 删除原 gene_id 列
  df <- df[, !names(df) %in% "gene_id"]
  # 将分隔的两列合并到原数据框中
  df <- cbind(df, gene_split_df)

  # 保留指定列
  selected_columns <- c(
    "seqnames", "context", "regionType",
    "gene_type", "gene_name", "gene_id", "version"
  )
  dmr_genes <- df[, selected_columns, drop = FALSE]

  # 所有列联合去重
  dmr_genes <- unique(dmr_genes)
  # 按 gene_id 列去重
  # dmr_genes <- dmr_genes[!duplicated(dmr_genes$gene_id), ]
  # 按 gene_id和regionType 列去重
  # dmr_genes <- dmr_genes[!duplicated(dmr_genes[, c("gene_id", "regionType")]), ]
  # 排序、重置行号
  dmr_genes <- dmr_genes[order(
    dmr_genes$seqnames, dmr_genes$context, dmr_genes$regionType,
    dmr_genes$gene_type, dmr_genes$gene_name, dmr_genes$gene_id
  ), ]
  rownames(dmr_genes) <- NULL

  return(dmr_genes)
}


# 读取数据
file_path <- paste0(config$output_dir, "/DMR_gene_association.bed")
data <- read.table(file_path,
  sep = "\t", header = FALSE,
  stringsAsFactors = FALSE
)

# 设置列名
colnames(data) <- c(
  "seqnames", "start", "end", "width", "strand", "sumReadsM1", "sumReadsN1",
  "proportion1", "sumReadsM2", "sumReadsN2", "proportion2", "cytosinesCount",
  "context", "direction", "pValue", "regionType",
  "gene_seqnames", "gene_start", "gene_end", "gene_id", "unkown1",
  "strand", "source", "feature_type", "unkown2",
  "extra"
)

# 删除无效列和重复列
data <- data[, !(names(data) %in% c(
  "gene_seqnames", "gene_id", "unkown1", "unkown2"
))]

# 解析注释内容
dmr_gene_all <- parse_extra_column(data)

# 导出到文件
write.table(dmr_gene_all,
  file = paste0(config$report_dir, "DMR_gene_all.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = ""
)

# 分隔gene_id并去重
dmr_genes <- filtered_genes(dmr_gene_all)
# 导出到文件
write.table(dmr_genes,
  file = paste0(config$report_dir, "DMR_genes.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("分析结束！\n")
