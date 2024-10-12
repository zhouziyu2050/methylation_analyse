############################################
# 参数读取
############################################

# 加载 optparse 包
library(optparse)

# 定义命令行参数的规格
option_list <- list(
  make_option(c("--species", "-s"), metavar = "<string>", type = "character", default = "mouse", help = "物种，可选值为human/mouse，默认值为mouse"),
  make_option(c("--genes", "-g"), metavar = "<file>", type = "character", default = NULL, help = "DMR输出的基因文件路径，可以使用相对路径或绝对路径"),
  make_option(c("--report_dir", "-r"), metavar = "<folder>", type = "character", default = NULL, help = "输出报告的文件夹路径，可选，默认值为genes所在文件夹"),
  make_option(c("--pathways_selected", "-p"), metavar = "<string>", type = "character", default = NULL, help = "指定通路，如'GO:0007015,GO:0007264'，多个通路以英文逗号连接")
  # GO:0007015,GO:0007264,GO:1902903,GO:0032970
)

# 解析命令行参数
parser <- OptionParser(option_list = option_list)
options <- parse_args(parser, commandArgs(TRUE))

# 检查是否提供了 --genes 参数
if (is.null(options$genes)) {
  stop("genes参数不可为空。\n")
}

# 提取参数值
DMR_File <- options$genes
report_dir <- options$report_dir
species <- options$species
cat("当前所选物种：", species, "\n")
# 检查文件是否存在
if (file.exists(DMR_File)) {
  # 输出文件存在信息
  cat("输入 DMR gene 文件：", DMR_File, "\n")
} else {
  # 文件不存在，输出错误信息
  stop("文件不存在：", DMR_File, "\n")
}

# 如果提供了 --report_dir 参数，则使用自定义的报告目录
if (!is.null(options$report_dir)) {
  report_dir <- options$report_dir
} else {
  # 否则从文件路径中提取目录部分
  report_dir <- dirname(DMR_File)
}
cat("输出报告文件夹路径：", report_dir, "\n")

# 如果提供了 --pathways_selected 参数，则将其格式化为向量
if (!is.null(options$pathways_selected)) {
  pathways_selected <- options$pathways_selected
  pathways_selected <- gsub(" ", "", pathways_selected)
  pathways_selected <- gsub("，", ",", pathways_selected)
  pathways_selected <- unlist(strsplit(pathways_selected, ","))
  cat("指定GO通路：", pathways_selected, "\n")
}

############################################
# 加载包
############################################

# 静默加载依赖包，避免显示过多提示
cat("依赖包加载中...", "\n")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(biomaRt))
if (species == "mouse") {
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  dataset <- "mmusculus_gene_ensembl"
  organism <- "mmu"
  OrgDb <- org.Mm.eg.db
} else if (species == "human") {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  dataset <- "hsapiens_gene_ensembl"
  organism <- "hsa"
  OrgDb <- org.Hs.eg.db
}
cat("依赖包加载完成", "\n")


############################################
# GO/KEGG分析
# 参考文档地址：https://lishensuo.github.io/posts/bioinfo/056clusterprofiler%E5%8C%85%E5%AF%8C%E9%9B%86%E5%88%86%E6%9E%90%E4%B8%8E%E5%8F%AF%E8%A7%86%E5%8C%96/#2gsea%e6%89%93%e5%88%86
############################################


# 读取DMR导出的表
d <- read.csv(DMR_File, sep = "\t")

region_types <- c("gain", "loss", "all")
# 循环处理不同的 regionType
for (region_type in region_types) {
  cat("开始绘制", region_type, "类型...", "\n")
  if (region_type == "all") {
    # 去重后的gene_id列
    gene_ids <- unique(d$gene_id)
  } else {
    # 过滤regionType并去重
    gene_ids <- unique(d$gene_id[d$regionType == region_type])
  }

  # go富集分析
  cat("GO富集分析...", "\n")
  ego <- enrichGO(
    gene          = gene_ids, # 输入基因列表
    keyType       = "ENSEMBL", # 指定基因ID类型为 Ensembl 基因 ID
    OrgDb         = OrgDb, # 使用小鼠基因数据库
    ont           = "ALL", # 指定 GO 类别：CC（细胞组分）、BP（生物过程）、MF（分子功能）
    pAdjustMethod = "BH", # 多重假设检验校正方法，使用 Benjamini-Hochberg 方法
    pvalueCutoff  = 0.01, # p 值阈值
    qvalueCutoff  = 0.05, # q 值阈值
    readable      = TRUE # 是否将结果转换为可读的基因符号
  )
  # 导出结果为 TSV 文件
  write.table(as.data.frame(ego),
    file = paste0(report_dir, "/GO富集-", region_type, ".tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  # 绘图
  barplot(
    ego,
    split = "ONTOLOGY",
    showCategory = 12,
    label_format = 50,
  ) + facet_grid(ONTOLOGY ~ ., scale = "free")
  ggsave(
    paste0(report_dir, "/GO富集-", region_type, ".png"),
    width = 8, height = 6
  )

  # 指定GO通路
  if (exists("pathways_selected") && !is.null(pathways_selected)) {
    cat("GO富集分析(指定通路)...", "\n")
    # 检查pathways_selected是否存在于ego结果中
    found_pathways <- intersect(pathways_selected, rownames(ego@result))
    not_found_pathways <- setdiff(pathways_selected, found_pathways) # 找到未匹配的通路
    if (length(found_pathways) > 0) {
      if (length(not_found_pathways) > 0) {
        cat("\033[31m警告：部分指定的通路未找到：", paste(not_found_pathways, collapse = ", "), "\033[0m\n")
      }
      # 绘制匹配到的通路
      barplot(ego, showCategory = ego@result$Description[
        which(rownames(ego@result) %in% found_pathways)
      ])
      ggsave(
        paste0(report_dir, "/GO富集(指定通路)-", region_type, ".png"),
        width = 8, height = 6
      )
    } else {
      cat("\033[31m警告：所有指定的通路均未找到\033[0m\n")
    }
  }


  # 使用 clusterProfiler 包中的 simplify 函数对富集分析结果进行去冗余处理
  cat("GO富集分析(去冗余)...", "\n")
  ego_sim <- clusterProfiler::simplify(
    ego, # 输入的富集分析结果对象
    cutoff = 0.7, # 去冗余的阈值。相似度大于这个值的 GO term 将被合并
    measure = "Wang", # 相似度计算方法，这里指定为 "Wang"。Wang 方法基于信息内容来计算 GO term 的相似度
    by = "p.adjust", # 按哪个字段进行去冗余操作。表示将相似度高的 GO term 合并时，保留调整后的 p 值最低的 GO term
    select_fun = min # 选择保留 GO term 的标准，这里指定为 min。表示选择 p 值最小的 GO term
  )
  # 导出结果为 TSV 文件
  write.table(
    as.data.frame(ego_sim),
    file = paste0(report_dir, "/GO富集(去冗余)-", region_type, ".tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  # 绘图
  barplot(
    ego_sim,
    split = "ONTOLOGY",
    showCategory = 12,
    label_format = 50,
  ) + facet_grid(ONTOLOGY ~ ., scale = "free")
  ggsave(
    paste0(report_dir, "/GO富集(去冗余)-", region_type, ".png"),
    width = 8, height = 6
  )

  # 从 Ensembl Gene ID 转换到 Entrez Gene ID
  # mmusculus_gene_ensembl 小鼠数据集,hsapiens_gene_ensembl 人类数据集
  ensembl <- useMart("ensembl", dataset = dataset)
  ensembl_gene <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = ensembl
  )
  entrezgene_ids <- ensembl_gene$entrezgene_id

  # KEGG
  cat("KEGG分析(该步骤需要联网获取数据)...", "\n")
  ekg <- enrichKEGG(
    gene = entrezgene_ids, # 输入的差异表达基因。
    keyType = "kegg", # 可选值：kegg, ncbi-geneid, ncbi-proteinid, uniprot
    organism = organism, # 物种标识符，hsa:人类; mmu:小鼠
    pvalueCutoff = 0.05 # p 值的阈值，用于筛选富集的 KEGG 路径。
  )
  # 将 ekg 结果设置为可读格式
  ekg <- setReadable(ekg, OrgDb = OrgDb, keyType = "ENTREZID")

  # 导出结果为 TSV 文件
  write.table(
    as.data.frame(ekg),
    file = paste0(report_dir, "/KEGG富集-", region_type, ".tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  # 绘制气泡图。颜色映射P值，大小映射交集基因数(差异基因与通路基因集)，横轴表示比例(count/geneset)
  dotplot(ego, showCategory = 20, label_format = 50)
  ggsave(
    paste0(report_dir, "/KEGG富集-", region_type, ".png"),
    width = 8, height = 6
  )
}
