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

  # 绘制DMR甲基化位置分布图相关参数
  make_option(c("-p", "--plot_type"),
    type = "character", default = NULL,
    help = "DMR和甲基化位置分布图的绘制形式，可选值为line/bar/point，不传则不绘制此图", metavar = "character"
  ),
  make_option(c("-n", "--seqname"),
    type = "character", default = NULL,
    help = "绘制DMR和甲基化位置分布图的染色体名称", metavar = "character"
  ),
  make_option(c("-s", "--start"),
    type = "integer", default = NULL,
    help = "绘制DMR和甲基化位置分布图的起始位置", metavar = "integer"
  ),
  make_option(c("-e", "--end"),
    type = "integer", default = NULL,
    help = "绘制DMR和甲基化位置分布图的结束位置", metavar = "integer"
  ),
  make_option(c("-g", "--gtf_file"),
    type = "character", default = NULL,
    help = "绘制DMR和甲基化位置分布图所使用的gtf注释文件路径（下载地址：https://www.gencodegenes.org/）", metavar = "character"
  ),

  # 绘制DMR环形分布图相关参数
  make_option(c("-y", "--cytoband_path"),
    type = "character", default = NULL,
    help = "cytoband文件路径，不传则不绘制DMR环形分布图（下载地址：https://hgdownload.cse.ucsc.edu/goldenPath/mm39/database/cytoBandIdeo.txt）", metavar = "character"
  ),
  
  make_option(c("-t", "--text_num"),
    type = "integer", default = NULL,
    help = "在环内显示标签的数量（微调该参数使标签数量刚好铺满整个环）", metavar = "integer"
  ),
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
if (!is.null(config$plot_type)) {
  if (!config$plot_type %in% c("line", "bar", "point")) {
    stop("plot_type 参数无效，应为 'line'、'bar'、'point' 之一")
  }
  if (is.null(config$start) || is.null(config$end) || is.null(config$seqname) || is.null(config$gtf_file)) {
    stop("绘制DMR甲基化位置分布图时，'start'、'end'、'seqname'、'gtf_file'参数都不能为空")
  }
}


# 如果提供了samples_file参数，则从文件读取samples数据，否则从config$samples读取
if (!is.null(config$samples_file)) {
  samples <- read_data(config$samples_file)
} else {
  samples <- as.data.frame(config$samples)
}

# 设置默认值
if (is.null(config$output_dir)) {
  config$output_dir <- "./output"
}
if (is.null(config$report_dir)) {
  config$report_dir <- "./report"
}
if (is.null(config$text_num)) {
  config$text_num <- 88
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

# seqnames的类型（accession或chromosome，如果设置为accession会在生成dmr报告时自动转换为chromosome）
seqnames_type <- "accession"

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
# 绘制甲基化位置分布图
##################################################################

if (!is.null(config$plot_type)) {
  # 加载库
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(dplyr))
  # 读取 GTF注释 文件
  suppressPackageStartupMessages(library(rtracklayer))
  # 绘制基因位置
  suppressPackageStartupMessages(library(transPlotR))
  # 合并子图
  suppressPackageStartupMessages(library(cowplot))

  # 遍历读取甲基化分布数据
  df_list <- list()
  region <- GRanges(config$seqname, IRanges(config$start, config$end))
  # methylationDataFiltered <- subsetByOverlaps(methylationData, region)
  for (i in seq_len(nrow(samples))) {
    # 读取数据
    file_path <- paste0(
      sub("/$", "", samples[i, "output_dir"]), "/bismark_methylation/",
      samples[i, "prefix"], "_bismark_bt2_pe.deduplicated.CX_report.txt.chr",
      config$seqname, ".CX_report.txt.gz"
    )
    methylationData <- readBismark(file_path)
    # 过滤指定区域
    methylationDataFiltered <- subsetByOverlaps(methylationData, region)
    # 计算甲基化率
    methylation_proportion <- methylationDataFiltered$readsM / methylationDataFiltered$readsN

    df_list[[i]] <- data.frame(
      position = start(methylationDataFiltered),
      methylation_proportion = methylation_proportion,
      group_name = group_names[i],
      sample = sample_names[i],
      # label = paste(group_names[i], ":", sample_names[i]),
      stringsAsFactors = FALSE
    )
  }
  # 合并所有样本的数据框
  df <- do.call(rbind, df_list)
  # 删除0或者NA
  df <- subset(df, methylation_proportion != 0 & !is.na(methylation_proportion))
  # 将 sample 列转换为因子，并指定级别顺序，用于分面时排序
  df$sample <- factor(df$sample, levels = sample_names)

  # 创建颜色映射
  color_palette <- c("#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#D55E00")
  names(color_palette) <- unique(group_names)

  # 读取DMR结果文件
  print(paste0(config$output_dir, "/DMRsReplicatesBins.txt"))
  DMRsReplicatesBinsCombined <- read.table(
    paste0(config$output_dir, "/DMRsReplicatesBins.txt"),
    sep = "\t", header = FALSE, stringsAsFactors = FALSE
  )
  # 设置列名
  colnames(DMRsReplicatesBinsCombined) <- c(
    "seqnames", "start", "end", "width", "strand", "sumReadsM1", "sumReadsN1",
    "proportion1", "sumReadsM2", "sumReadsN2", "proportion2", "cytosinesCount",
    "context", "direction", "pValue", "regionType"
  )
  # 筛选指定范围和 seqnames 的行
  start_range_val <- config$start
  end_range_val <- config$end
  specified_seqnames <- accession2chromosome[config$seqname]
  filtered_DMRs <- DMRsReplicatesBinsCombined %>%
    filter(start >= start_range_val & end <= end_range_val & seqnames == specified_seqnames)


  # 读取 GTF注释 文件
  gtf_data <- import(config$gtf_file, format = "gtf")

  # 绘散点图
  if (config$plot_type == "point") {
    options(repr.plot.width = 12, repr.plot.height = 6)
    plot1 <- ggplot(df, aes(x = position, y = methylation_proportion)) +
      geom_point(aes(
        color = group_name,
        shape = sample
      ), size = 2) +
      scale_color_manual(values = color_palette, labels = names(color_palette)) +
      scale_shape_manual(values = seq_along(sample_names), labels = sample_names) +
      labs(
        x = paste0("Position in ", accession2chromosome[config$seqname], " (bp)"),
        y = "Methylation proportion",
        color = "Group",
        shape = "Sample",
      ) +
      # ylim(0, 1) +
      theme_minimal()
  }

  # 绘制折线图
  if (config$plot_type == "line") {
    options(repr.plot.width = 12, repr.plot.height = 6)
    plot1 <- ggplot(data = df, mapping = aes(
      x = position,
      y = methylation_proportion,
      color = group_name
    )) +
      geom_line() +
      labs(
        x = paste0("Position in ", accession2chromosome[config$seqname], " (bp)"),
        y = "Methylation Proportion",
        fill = "Group"
      ) +
      scale_color_manual(values = color_palette) + # 设置 group 的颜色
      theme_minimal() +
      theme(
        strip.text.y.right = element_text(angle = 270), # 调整标题文本方向
        strip.placement = "outside" # 将标题放在外部
      ) +
      facet_grid(sample ~ .) # 按 sample 列分面，竖直排列并切换到右侧
  }

  # 绘制柱状图
  if (config$plot_type == "bar") {
    options(repr.plot.width = 12, repr.plot.height = 6)
    plot1 <- ggplot(data = df, mapping = aes(
      x = position,
      y = methylation_proportion,
      fill = group_name
    )) +
      geom_bar(stat = "identity", position = "dodge") + # 使用柱状图
      labs(
        x = paste0("Position in ", accession2chromosome[config$seqname], " (bp)"),
        y = "Methylation Proportion",
        fill = "Group"
      ) +
      scale_fill_manual(values = color_palette) + # 设置 group 的颜色
      theme_minimal() +
      theme(
        strip.text.y.right = element_text(angle = 270), # 调整文本方向
        strip.placement = "outside" # 将标签放在外部
      ) +
      facet_grid(sample ~ .) # 按 sample 列分面，竖直排列并切换到右侧
  }

  # 绘制DMR区域高亮
  if (nrow(filtered_DMRs) > 0) { # 检查 filtered_DMRs 是否有数据
    # 创建一个数据框用于绘制矩形
    rect_data <- data.frame(
      xmin = filtered_DMRs$start,
      xmax = filtered_DMRs$end,
      ymin = -Inf,
      ymax = Inf
    )

    # 添加矩形层
    plot1 <- plot1 +
      geom_rect(
        data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "lightblue", alpha = 0.5, inherit.aes = FALSE
      )
  }

  # 绘制基因位置
  # 教程地址：https://blog.csdn.net/qazplm12_3/article/details/125903915
  options(repr.plot.width = 12, repr.plot.height = 2)
  plot2 <- trancriptVis(
    gtfFile = as.data.frame(gtf_data), # 将 GTF 数据文件转换为数据框，并用于绘制转录本图
    Chr = accession2chromosome[config$seqname], # 将染色体的 accession 转换为实际染色体名称
    posStart = config$start, # 可视化区域的起始位置
    posEnd = config$end, # 可视化区域的结束位置
    textLabel = "gene_name", # 标签文字类型，可选 gene_name, transcript_name, transcript_id, gene_id
    collapse = TRUE, # 当多个转录本重叠时，是否将其折叠显示为一行
    textLabelSize = 4, # 标签文字大小
    relTextDist = -0.22, # 标签文字距离箭头的相对距离
    # textLabelColor = "red",            # 设置标签文字颜色（默认颜色为黑色）
    # exonColorBy = "transcript_name",   # 根据转录本名称为外显子着色（此行注释掉）
    # arrowType = "closed",              # 设置箭头的类型为闭合箭头（默认为开放箭头）
    arrowCol = "#e69f00", # 设置箭头的颜色为橙色
    exonFill = "#56b4e9", # 设置外显子区域的填充颜色为蓝色
    arrowBreak = 0.5 # 设置箭头的断开位置
  ) +
    labs(x = paste0("Position in ", accession2chromosome[config$seqname], " (bp)")) +
    ylim(0.7, 1.2) +
    scale_x_continuous(
      limits = c(config$start, config$end) # 设置 x 轴的显示范围为 start 和 end 之间
    ) # + theme_void()           # （可选）隐藏所有坐标轴和背景网格线

  # 合并子图
  options(repr.plot.width = 12, repr.plot.height = 8)
  combined_plot <- plot_grid(
    plot1 +
      theme(
        plot.background = element_rect(fill = "white", colour = "white"), # 白色背景
        axis.title.x = element_blank(), # 移除 x 轴标题
        axis.ticks.x = element_blank(), # 移除 x 轴刻度线
        axis.text.x = element_blank(), # 移除 x 轴刻度标签
        panel.border = element_blank(), # 移除面板边框
        plot.margin = margin(5, 5, 0, 5) # 设置图形周围的边距
      ),
    plot2 +
      theme(
        plot.margin = margin(0, 0, 10, 0) # 设置图形周围的边距
      ),
    ncol = 1, # 设置为 1 列，垂直排列
    align = "v", # 垂直对齐
    axis = "lrtb", # 设置轴参数以确保对齐
    rel_heights = c(6, 2) # 设置相对高度
  )

  # 保存绘图到文件
  ggsave(
    paste0(config$report_dir, "/Position of DMR and methylation.png"),
    plot = combined_plot,
    width = 12, height = 8
  )
}

# 绘制DMR环形分布图
if (!is.null(config$cytoband_path)) {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(circlize))
  # 参考教程：https://www.jianshu.com/p/3b73b029ee8e

  # 读取DMR结果
  colname_lists <- c("chr", "start", "end", "width", "strand", "sumReadsM1", "sumReadsN1", "proportion1", "sumReadsM2", "sumReadsN2", "proportion2", "cytosinesCount", "context", "direction", "pValue", "regionType")
  DMRs <- read.table(paste0(config$output_dir, "/DMRsReplicatesBins.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(DMRs) <- colname_lists

  # 读取DMR注释结果
  DMR_gene <- read.table(paste0(config$report_dir, "/DMR_gene_all.tsv"), sep = "\t", header = 1, stringsAsFactors = FALSE)
  DMR_gene <- DMR_gene[DMR_gene$feature_type == "gene", ]
  DMR_gene <- DMR_gene %>%
    group_by(seqnames, gene_name) %>%
    summarise(
      start = first(gene_start), # 保留该组的第一个 gene_start
      end = first(gene_end), # 保留该组的第一个 gene_end
      proportion_ratio = sum(proportion1) / sum(proportion2),
      .groups = "drop" # 去掉分组信息
    ) %>%
    arrange(desc(proportion_ratio)) # 按 value 降序排列
  DMR_gene <- DMR_gene[c("seqnames", "start", "end", "gene_name")]

  # 开始绘图
  par(bg = "white") # 设置背景为白色
  # 设置画布大小
  # options(repr.plot.width = 8, repr.plot.height = 8)
  # 设置输出文件名和尺寸
  png(paste0(config$report_dir, "/circos plot of DMR.png"), width = 2000, height = 2000, pointsize = 50)

  circos.par(track.margin = c(0.00, 0.01)) # 调整上边距和下边距（单位是相对半径）

  # 初始化坐标轴，并绘制染色体分布图（红线表示着丝粒）
  # circos.initializeWithIdeogram(
  #   cytoband = cytoband_path,
  #   plotType = c("ideogram", "labels", "axis"),
  #   chromosome.index = accession2chromosome,
  #   # track.height = 0.2,
  #   ideogram.height = 0.05
  # )

  # 初始化坐标轴
  circos.initializeWithIdeogram(
    cytoband = config$cytoband_path,
    plotType = c("labels"),
    chromosome.index = accession2chromosome,
    labels.cex = 0.6 * par("cex"), # 染色体标签大小
  )

  # 第一圈：绘制染色体分布图（红线表示着丝粒）
  circos.genomicIdeogram(
    cytoband = config$cytoband_path,
    track.height = 0.03,
    track.margin = c(0.01, 0.04) # 修改染色体图与染色体标签的距离
  )

  # 在染色体外层绘制坐标轴刻度标签
  circos.track(
    track.index = get.current.track.index(),
    panel.fun = function(x, y) {
      circos.genomicAxis(
        h = "top",
        direction = "outside",
        labels.cex = 0.4 * par("cex"), # 坐标轴标签大小
      )
    },
    track.height = 0.03
  )

  # 第二圈：绘制DMR分布密度图（只需要"chr", "start", "end"列的数据）
  circos.genomicDensity(
    DMRs,
    col = c("#67a9cf80"),
    track.height = 0.1
  )
  # 第三圈：绘制DMR基因名称
  circos.genomicLabels(
    head(DMR_gene, config$text_num), # 需要调整head数量控制标签刚好铺满一圈
    labels.column = 4,
    side = "inside",
    # col = "#ef8a62",
    # labels_height = 0.3, # 调整文字的高度
    connection_height = 0.08, # 调整连接线的长度
  )
  circos.clear()
  dev.off() # 关闭绘图设备，保存绘图文件
}
