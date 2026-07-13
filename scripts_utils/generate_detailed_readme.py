#!/usr/bin/env python3

from __future__ import annotations

import csv
import subprocess
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
README_PATH = PROJECT_ROOT / "README.md"


SCRIPT_DESCRIPTIONS = {
    "01_download_and_prepare_scRNA.R": (
        "下载并整理GSE215403单细胞RNA-seq原始数据，"
        "核对GEO补充文件，建立按样本组织的10x表达矩阵目录。"
    ),
    "02_read_and_QC_scRNA.R": (
        "读取各样本10x表达矩阵，创建Seurat对象，计算nFeature、nCount和线粒体比例，"
        "完成基础质控检查并保存过滤前对象。"
    ),
    "03_QC_reproduction_candidate.R": (
        "按照论文复现目标重新评估质控阈值，对过滤前后细胞数和QC指标进行比较，"
        "保存候选QC对象和统计结果。"
    ),
    "04_standard_Seurat_PCA_UMAP_resolution_scan.R": (
        "执行标准Seurat流程，包括归一化、高变基因选择、ScaleData、PCA、Harmony、"
        "邻居图、UMAP以及多resolution聚类扫描。"
    ),
    "05_manual_annotation_and_target_gene_summary.R": (
        "基于经典marker、cluster差异基因和样本构成进行人工细胞类型注释，"
        "并总结SASH1、COL1A1、EMP1和MYH11表达。"
    ),
    "06_malignant_candidate_diagnostic.R": (
        "针对上皮和肿瘤相关cluster进行进一步诊断，结合marker、状态评分、"
        "样本构成和轨迹信息筛选恶性候选群体。"
    ),
    "07_CopyKAT_malignant_call.R": (
        "使用CopyKAT推断细胞拷贝数状态，将aneuploid预测作为恶性细胞判断的辅助证据。"
    ),
    "07_CopyKAT_malignant_call_batch.sh": (
        "脚本07的批处理启动脚本，用于长时间运行CopyKAT并记录日志。"
    ),
    "08_finalize_malignant_call.R": (
        "整合人工注释、cluster信息、上皮状态和CopyKAT结果，生成最终恶性细胞标签。"
    ),
    "08_finalize_malignant_call_batch.sh": (
        "脚本08的批处理启动脚本。"
    ),
    "09_malignant_cell_composition_check.R": (
        "检查恶性细胞在不同样本、cluster和细胞类型中的数量及比例，"
        "评估后续轨迹分析的样本可用性。"
    ),
    "09_malignant_cell_composition_check_batch.sh": (
        "脚本09的批处理启动脚本。"
    ),
    "10_manual_review_epithelial_core.R": (
        "对上皮细胞和核心肿瘤相关cluster进行人工复核，比较不同状态评分、"
        "样本组成和核心基因表达。"
    ),
    "11a_malignant_focus_cluster_audit.R": (
        "审查核心、扩展和候选恶性focus clusters，综合CopyKAT支持率、marker、"
        "样本组成和轨迹位置形成最终cluster层面判断。"
    ),
    "11b_core_extended_malignant_overlay_rebuild.R": (
        "重建核心恶性focus和扩展恶性focus覆盖标签，整理Monocle3轨迹坐标、"
        "主图边和顶点信息，供脚本16使用。"
    ),
    "12_spatial_download_QC_gene_maps.R": (
        "下载并读取GSE252265空间转录组表达矩阵和坐标文件，完成spot质控、"
        "barcode匹配检查及核心基因空间表达图。"
    ),
    "13_spatial_SASH1_COL1A1_neighborhood_analysis.R": (
        "定义SASH1-high和COL1A1-high spots，计算重叠比例、最近邻距离及空间置换检验。"
    ),
    "14_spatial_domain_annotation_and_core_gene_maps.R": (
        "根据空间marker module score近似注释Fibrotic Stroma、Inflammatory Zone、"
        "不同肿瘤状态和CSC-like Niche，并生成Fig.6-like图。"
    ),
    "15_scRNA_core_gene_expression_Figure4_like.R": (
        "在主要细胞类型层面生成SASH1、COL1A1、EMP1和MYH11的UMAP、DotPlot、"
        "FeaturePlot和VlnPlot，形成Fig.4-like组合图。"
    ),
    "15b_scRNA_core_gene_expression_Figure4_cluster_like.R": (
        "在cluster与cell type组合层面生成更细粒度的Fig.4-like核心基因表达图。"
    ),
    "16_scRNA_malignant_trajectory_Figure5_like.R": (
        "将恶性focus标签、轨迹坐标和核心基因表达结合，生成Fig.5-like轨迹、"
        "gene-high状态和拟时序或拟时序代理趋势图。"
    ),
    "17_scTenifold_virtual_KO_OE.R": (
        "选取cluster4、6和11中的7353个恶性focus细胞，使用scTenifoldKnk和"
        "scTenifoldNet对SASH1、COL1A1、EMP1和MYH11执行虚拟KO与虚拟OE分析。"
    ),
    "17_scTenifold_virtual_KO_OE_batch.sh": (
        "脚本17的批处理启动脚本，适合后台运行耗时较长的scTenifold网络分析。"
    ),
}


STAGE_SUMMARIES = {
    "01": "单细胞原始数据下载、文件核对和10x目录整理",
    "02": "单细胞表达矩阵读取和基础QC",
    "03": "论文复现候选QC过滤",
    "04": "Seurat标准降维、Harmony整合和聚类",
    "05": "人工细胞类型注释和目标基因表达总结",
    "06": "上皮及恶性候选细胞诊断",
    "07": "CopyKAT拷贝数推断和恶性支持分析",
    "08": "最终恶性细胞标签整合",
    "09": "恶性细胞样本和cluster组成检查",
    "10": "上皮核心群体人工复核",
    "11": "恶性focus cluster审查和轨迹输入重建",
    "12": "空间转录组下载、QC和核心基因空间图",
    "13": "SASH1-high与COL1A1-high空间关系分析",
    "14": "空间结构域近似注释和Fig.6-like图",
    "15": "Fig.4-like单细胞核心基因表达",
    "16": "Fig.5-like恶性细胞轨迹",
    "17": "scTenifold虚拟KO和虚拟OE",
}


def git_tracked_files() -> list[str]:
    result = subprocess.run(
        ["git", "ls-files"],
        cwd=PROJECT_ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    return sorted(line.strip() for line in result.stdout.splitlines() if line.strip())


def stage_from_name(filename: str) -> str | None:
    for stage in [
        "01", "02", "03", "04", "05", "06", "07", "08", "09",
        "10", "11", "12", "13", "14", "15", "16", "17",
    ]:
        if filename.startswith(stage):
            return stage
    return None


def describe_result_file(path_str: str) -> str:
    path = Path(path_str)
    name = path.name
    suffix = path.suffix.lower()
    stage = stage_from_name(name)

    if "sessionInfo" in name:
        return "记录该分析阶段的R版本、操作系统和已加载R包版本，用于环境追踪和可重复性检查。"

    if "output_file_check" in name:
        return "检查该阶段预期输出文件或目录是否真实存在。"

    if "run_status_summary" in name:
        return "汇总运行输入、细胞数量以及成功、失败和跳过的分析数量。"

    if "runtime_config" in name:
        return "记录脚本实际使用的输入对象、细胞筛选来源、目标基因、网络参数和软件版本。"

    if "summary" in name.lower():
        base = "汇总该分析阶段的主要统计结果"
    elif "metadata" in name.lower():
        base = "保存细胞或spot层面的metadata"
    elif "marker" in name.lower():
        base = "保存marker基因或marker表达统计"
    elif "expression" in name.lower():
        base = "保存基因表达结果"
    elif "trajectory" in name.lower():
        base = "保存轨迹分析相关结果"
    elif "copykat" in name.lower():
        base = "保存CopyKAT预测或支持率结果"
    elif "permutation" in name.lower():
        base = "保存空间置换检验结果"
    elif "distance" in name.lower():
        base = "保存空间距离统计结果"
    elif "threshold" in name.lower():
        base = "保存分析中使用的阈值"
    elif "check" in name.lower():
        base = "保存文件、基因或分析条件检查结果"
    elif "all_results" in name.lower():
        base = "保存全部差异调控或统计结果"
    elif "significant_results" in name.lower():
        base = "保存达到显著性阈值的结果"
    else:
        base = "保存该分析阶段生成的结果"

    if suffix == ".csv":
        file_type = "，文件格式为CSV，可用于后续统计或绘图。"
    elif suffix == ".txt":
        file_type = "，文件格式为纯文本。"
    elif suffix == ".pdf":
        file_type = "，文件格式为PDF矢量图。"
    elif suffix == ".png":
        file_type = "，文件格式为PNG预览图。"
    elif suffix == ".rds":
        file_type = "，文件格式为RDS对象，供后续R脚本直接读取。"
    else:
        file_type = "。"

    if stage:
        return f"由脚本{stage}生成，{base}{file_type}"

    return f"{base}{file_type}"


def describe_figure(path_str: str) -> str:
    name = Path(path_str).name
    stage = stage_from_name(name)

    keywords = [
        ("volcano", "差异调控火山图"),
        ("top_genes_barplot", "Top差异调控基因条形图"),
        ("UMAP", "UMAP降维可视化"),
        ("DotPlot", "marker或目标基因DotPlot"),
        ("VlnPlot", "基因表达小提琴图"),
        ("FeaturePlot", "基因表达FeaturePlot"),
        ("heatmap", "热图"),
        ("trajectory", "轨迹或拟时序相关图"),
        ("spatial", "空间坐标或空间表达图"),
        ("colocalization", "空间共定位统计图"),
        ("permutation", "空间置换检验图"),
        ("CopyKAT", "CopyKAT预测或支持率图"),
        ("QC", "质量控制图"),
        ("PCA", "PCA相关图"),
        ("elbow", "PCA肘部图"),
        ("composition", "细胞组成图"),
        ("expression", "基因表达图"),
        ("marker", "marker基因图"),
    ]

    plot_type = "分析结果图"
    lower_name = name.lower()

    for keyword, description in keywords:
        if keyword.lower() in lower_name:
            plot_type = description
            break

    if stage:
        return f"由脚本{stage}生成的{plot_type}。"

    return f"项目生成的{plot_type}。"


def describe_object(path_str: str) -> str:
    path = Path(path_str)
    name = path.name
    stage = stage_from_name(name)

    if ".parts" in path_str or ".part_" in name:
        return "大型RDS对象的分片文件，需要使用scripts_utils中的重建脚本合并后读取。"

    if "OE_networks" in path_str:
        return "脚本17保存的WT网络、虚拟过表达网络和manifold alignment结果。"

    if "17_scTenifold_analysis_input" in name:
        return "脚本17实际用于虚拟KO/OE分析的7353个cluster4、6和11细胞Seurat对象。"

    if "spatial" in name.lower():
        return "空间转录组分析生成的Seurat中间对象，保存表达矩阵、空间坐标和新增metadata。"

    if "trajectory" in name.lower():
        return "轨迹分析相关Seurat或Monocle中间对象。"

    if stage:
        return f"脚本{stage}生成的Seurat或分析中间对象，供后续脚本继续读取。"

    return "分析流程生成的RDS中间对象。"


def describe_data_file(path_str: str) -> str:
    path = Path(path_str)
    name = path.name

    if name == "barcodes.tsv.gz":
        return "10x格式barcode列表，对应表达矩阵中的细胞列。"
    if name == "features.tsv.gz":
        return "10x格式feature列表，记录基因ID、基因符号和feature类型。"
    if name == "matrix.mtx.gz":
        return "10x格式稀疏表达矩阵，行为基因，列为细胞。"
    if name.endswith("_barcodes.tsv.gz"):
        return "从GEO补充文件整理得到的样本barcode列表。"
    if name.endswith("_features.tsv.gz"):
        return "从GEO补充文件整理得到的样本feature注释表。"
    if name.endswith("_matrix.mtx.gz"):
        return "从GEO补充文件整理得到的样本稀疏表达矩阵。"
    if name.endswith("_RAW.tar"):
        return "GEO提供的补充原始文件压缩包，由Git LFS管理。"
    if "filtered_feature_bc_matrix.h5" in name:
        return "过滤后的feature-barcode表达矩阵，HDF5格式。"
    if "tissue_positions" in name:
        return "空间spot组织坐标表，用于将barcode映射到空间位置。"
    if "aggregation" in name:
        return "GSE252265聚合信息文件。"
    if "barcodes" in name:
        return "空间转录组barcode列表。"
    if "features" in name:
        return "空间转录组feature注释表。"
    if "matrix" in name:
        return "空间转录组Matrix Market表达矩阵。"
    return "项目输入或整理后的数据文件。"


def describe_file(path_str: str) -> str:
    if path_str == ".gitattributes":
        return "定义Git LFS追踪规则，用于管理RDS、压缩包和其他大型文件。"
    if path_str == ".gitignore":
        return "定义不进入Git版本控制的临时文件、本地缓存和可重建文件。"
    if path_str.endswith(".Rproj"):
        return "RStudio项目文件，用于从项目根目录打开和运行分析。"
    if path_str == "README.md":
        return "本项目的总体说明、运行顺序、分析结果、文件清单和限制。"

    if path_str.startswith("config/"):
        return "样本metadata或项目配置文件。"

    if path_str.startswith("data/"):
        return describe_data_file(path_str)

    if path_str.startswith("scripts/"):
        return SCRIPT_DESCRIPTIONS.get(
            Path(path_str).name,
            "主分析或批处理脚本。"
        )

    if path_str.startswith("scripts_utils/"):
        if "reconstruct_split_rds_objects" in path_str:
            return "重建脚本11、15、15b和16产生的大型分片RDS对象。"
        if "reconstruct_split_files" in path_str:
            return "根据manifest重建通用大型分片文件。"
        if "generate_detailed_readme" in path_str:
            return "自动扫描Git追踪文件并生成详细README。"
        return "项目辅助工具脚本。"

    if path_str.startswith("split_file_manifest/"):
        return "记录大型文件分片路径、大小或重建信息。"

    if path_str.startswith("results/figures/"):
        return describe_figure(path_str)

    if path_str.startswith("results/tables/"):
        return describe_result_file(path_str)

    if path_str.startswith("results/objects/"):
        return describe_object(path_str)

    return "项目文件。"


def read_sc_tenifold_summary() -> list[dict[str, str]]:
    summary_path = (
        PROJECT_ROOT
        / "results"
        / "tables"
        / "17_scTenifold_KO_OE_summary.csv"
    )

    if not summary_path.exists():
        return []

    with summary_path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def markdown_table(rows: list[tuple[str, str]]) -> list[str]:
    output = [
        "| 文件 | 作用 |",
        "|---|---|",
    ]

    for file_name, description in rows:
        safe_description = description.replace("|", "\\|")
        output.append(f"| `{file_name}` | {safe_description} |")

    return output


def build_stage_output_mapping(files: list[str]) -> list[str]:
    lines: list[str] = [
        "## 01–17脚本与输出文件对应关系",
        "",
        "以下清单根据当前Git已追踪文件自动生成。",
        "",
    ]

    result_files = [
        path for path in files
        if path.startswith("results/")
    ]

    for stage, summary in STAGE_SUMMARIES.items():
        matched = [
            path for path in result_files
            if Path(path).name.startswith(stage)
        ]

        lines.extend([
            f"### 阶段{stage}：{summary}",
            "",
        ])

        matching_scripts = [
            path for path in files
            if path.startswith("scripts/")
            and Path(path).name.startswith(stage)
        ]

        if matching_scripts:
            lines.append("对应代码：")
            lines.append("")
            for script in matching_scripts:
                lines.append(
                    f"- `{script}`："
                    f"{SCRIPT_DESCRIPTIONS.get(Path(script).name, '分析脚本。')}"
                )
            lines.append("")

        if matched:
            lines.append("当前已追踪输出：")
            lines.append("")
            for output_file in matched:
                lines.append(
                    f"- `{output_file}`：{describe_file(output_file)}"
                )
        else:
            lines.append("当前Git中没有以该编号命名的独立结果文件。")

        lines.append("")

    return lines


def build_inventory(files: list[str]) -> list[str]:
    grouped: dict[str, list[str]] = defaultdict(list)

    for path_str in files:
        path = Path(path_str)

        if len(path.parts) == 1:
            group = "根目录"
        elif path_str.startswith("data/processed/"):
            group = "data/processed/"
        elif path_str.startswith("data/raw/GSE252265/"):
            group = "data/raw/GSE252265/"
        elif path_str.startswith("data/raw/"):
            group = "data/raw/"
        elif path_str.startswith("results/figures/17_scTenifold_KO/"):
            group = "results/figures/17_scTenifold_KO/"
        elif path_str.startswith("results/figures/17_scTenifold_OE/"):
            group = "results/figures/17_scTenifold_OE/"
        elif path_str.startswith("results/figures/"):
            group = "results/figures/"
        elif path_str.startswith("results/tables/17_scTenifold_KO/"):
            group = "results/tables/17_scTenifold_KO/"
        elif path_str.startswith("results/tables/17_scTenifold_OE/"):
            group = "results/tables/17_scTenifold_OE/"
        elif path_str.startswith("results/tables/"):
            group = "results/tables/"
        elif path_str.startswith("results/objects/17_scTenifold_OE_networks/"):
            group = "results/objects/17_scTenifold_OE_networks/"
        elif path_str.startswith("results/objects/"):
            group = "results/objects/"
        else:
            group = path.parts[0] + "/"

        grouped[group].append(path_str)

    order = [
        "根目录",
        "config/",
        "data/raw/",
        "data/raw/GSE252265/",
        "data/processed/",
        "scripts/",
        "scripts_utils/",
        "split_file_manifest/",
        "results/figures/",
        "results/figures/17_scTenifold_KO/",
        "results/figures/17_scTenifold_OE/",
        "results/tables/",
        "results/tables/17_scTenifold_KO/",
        "results/tables/17_scTenifold_OE/",
        "results/objects/",
        "results/objects/17_scTenifold_OE_networks/",
    ]

    lines = [
        "## Git已追踪文件详细清单",
        "",
        "本节由`scripts_utils/generate_detailed_readme.py`自动生成。",
        "只列出当前Git实际追踪的文件，因此不会把`.DS_Store`、临时日志或未追踪文件写入README。",
        "",
    ]

    used_groups = set()

    for group in order:
        if group not in grouped:
            continue

        used_groups.add(group)
        rows = [
            (path_str, describe_file(path_str))
            for path_str in sorted(grouped[group])
        ]

        lines.extend([
            f"### {group}",
            "",
            *markdown_table(rows),
            "",
        ])

    for group in sorted(set(grouped) - used_groups):
        rows = [
            (path_str, describe_file(path_str))
            for path_str in sorted(grouped[group])
        ]

        lines.extend([
            f"### {group}",
            "",
            *markdown_table(rows),
            "",
        ])

    return lines


def main() -> None:
    files = git_tracked_files()
    sc_summary = read_sc_tenifold_summary()

    lines: list[str] = [
        "# HNSCC/OSCC SASH1多组学论文复现",
        "",
        "本仓库用于复现HNSCC/OSCC中SASH1相关研究的单细胞RNA-seq、"
        "空间转录组和基因调控网络虚拟扰动分析。",
        "",
        "当前01–17流程已经覆盖单细胞数据处理、人工注释、CopyKAT辅助恶性判断、"
        "恶性focus群体定义、Fig.4-like核心基因表达、Fig.5-like轨迹、"
        "GSE252265空间分析、Fig.6-like空间结构域以及scTenifold虚拟KO/OE。",
        "",
        "需要明确的是：本仓库已完成论文的单细胞、空间转录组和虚拟扰动相关模块，"
        "但尚未完整覆盖bulk RNA-seq、机器学习特征筛选、TCGA预后模型和外部验证队列。",
        "",
        "## 当前完成状态",
        "",
        "- GSE215403单细胞RNA-seq原始数据下载和10x矩阵整理",
        "- 单细胞QC、Harmony整合、PCA、UMAP和多resolution聚类",
        "- 基于marker和cluster差异基因的人工细胞类型注释",
        "- 上皮及肿瘤相关细胞群体复核",
        "- CopyKAT辅助恶性细胞判断",
        "- 核心恶性focus clusters：6和11",
        "- 扩展恶性focus cluster：4",
        "- 候选恶性相关肿瘤clusters：2和3",
        "- Fig.4-like单细胞核心基因表达图",
        "- Fig.5-like恶性细胞轨迹和核心基因表达图",
        "- GSE252265空间转录组下载、坐标匹配和spot QC",
        "- Fig.6-like空间结构域近似注释和核心基因空间图",
        "- SASH1-high与COL1A1-high共定位、最近邻距离和置换检验",
        "- SASH1、COL1A1、EMP1和MYH11虚拟KO与虚拟OE",
        "- 8个scTenifold分析成功，0失败，0跳过",
        "",
        "## 数据集",
        "",
        "### GSE215403",
        "",
        "用于单细胞RNA-seq分析。原始补充文件和整理后的10x矩阵分别位于：",
        "",
        "```text",
        "data/raw/scRNA_GSE215403/",
        "data/processed/scRNA_GSE215403/",
        "```",
        "",
        "### GSE252265",
        "",
        "用于空间转录组分析。表达矩阵、barcode、feature和空间坐标位于：",
        "",
        "```text",
        "data/raw/GSE252265/",
        "```",
        "",
        "## 仓库目录结构",
        "",
        "```text",
        "config/",
        "data/",
        "  raw/",
        "  processed/",
        "results/",
        "  figures/",
        "  tables/",
        "  objects/",
        "scripts/",
        "scripts_utils/",
        "split_file_manifest/",
        "README.md",
        "HNSCC_SASH1_reproduction.Rproj",
        "```",
        "",
        "- `config/`：样本metadata和项目配置。",
        "- `data/raw/`：GEO原始补充文件和空间原始输入。",
        "- `data/processed/`：按样本整理后的10x表达矩阵。",
        "- `scripts/`：01–17主分析脚本和对应batch脚本。",
        "- `scripts_utils/`：大型RDS重建和README生成工具。",
        "- `results/figures/`：PDF和PNG结果图。",
        "- `results/tables/`：CSV统计结果、运行参数、输出检查和sessionInfo。",
        "- `results/objects/`：Seurat对象、轨迹对象、网络对象和大型对象分片。",
        "- `split_file_manifest/`：大型分片文件的重建记录。",
        "",
        "## 01–17主流程",
        "",
        "| 脚本 | 作用 |",
        "|---|---|",
    ]

    for script_name, description in SCRIPT_DESCRIPTIONS.items():
        lines.append(f"| `scripts/{script_name}` | {description} |")

    lines.extend([
        "",
        "## 推荐运行顺序",
        "",
        "在仓库根目录依次运行：",
        "",
        "```bash",
        "Rscript scripts/01_download_and_prepare_scRNA.R",
        "Rscript scripts/02_read_and_QC_scRNA.R",
        "Rscript scripts/03_QC_reproduction_candidate.R",
        "Rscript scripts/04_standard_Seurat_PCA_UMAP_resolution_scan.R",
        "Rscript scripts/05_manual_annotation_and_target_gene_summary.R",
        "Rscript scripts/06_malignant_candidate_diagnostic.R",
        "Rscript scripts/07_CopyKAT_malignant_call.R",
        "Rscript scripts/08_finalize_malignant_call.R",
        "Rscript scripts/09_malignant_cell_composition_check.R",
        "Rscript scripts/10_manual_review_epithelial_core.R",
        "Rscript scripts/11a_malignant_focus_cluster_audit.R",
        "Rscript scripts/11b_core_extended_malignant_overlay_rebuild.R",
        "Rscript scripts/12_spatial_download_QC_gene_maps.R",
        "Rscript scripts/13_spatial_SASH1_COL1A1_neighborhood_analysis.R",
        "Rscript scripts/14_spatial_domain_annotation_and_core_gene_maps.R",
        "Rscript scripts/15_scRNA_core_gene_expression_Figure4_like.R",
        "Rscript scripts/15b_scRNA_core_gene_expression_Figure4_cluster_like.R",
        "Rscript scripts/16_scRNA_malignant_trajectory_Figure5_like.R",
        "Rscript scripts/17_scTenifold_virtual_KO_OE.R",
        "```",
        "",
        "耗时较长的步骤可使用：",
        "",
        "```bash",
        "bash scripts/07_CopyKAT_malignant_call_batch.sh",
        "bash scripts/08_finalize_malignant_call_batch.sh",
        "bash scripts/09_malignant_cell_composition_check_batch.sh",
        "bash scripts/17_scTenifold_virtual_KO_OE_batch.sh",
        "```",
        "",
        "## 恶性细胞focus定义",
        "",
        "```text",
        "核心恶性focus clusters：6和11",
        "扩展恶性focus cluster：4",
        "候选恶性相关肿瘤clusters：2和3",
        "```",
        "",
        "该定义综合考虑人工注释、上皮/肿瘤marker、样本组成、轨迹位置和CopyKAT支持率。",
        "CopyKAT仅作为辅助证据，不作为单独的恶性判定标准。",
        "",
        "脚本17使用cluster4、6和11，共7353个细胞作为虚拟扰动输入。",
        "",
        "## 空间转录组分析说明",
        "",
        "当前空间分析使用GSE252265公开表达矩阵和坐标文件。",
        "由于公开文件缺少标准Seurat空间对象所需的完整H&E图像和明确的按样本图像拆分信息，"
        "当前分析采用基于坐标的空间表达和marker-score结构域近似注释。",
        "",
        "当前限制：",
        "",
        "```text",
        "空间样本ID统一记为All_spots",
        "没有H&E病理图像叠加",
        "空间结构域属于marker-score近似注释",
        "不等同于原文作者基于病理图像完成的人工区域标注",
        "```",
        "",
        "当前空间分析显示SASH1-high和COL1A1-high spots重叠有限；"
        "在现有置换检验框架下，没有观察到明显空间共定位富集。",
        "",
        "## scTenifold虚拟KO/OE分析",
        "",
        "脚本17使用：",
        "",
        "```text",
        "scTenifoldKnk 1.1",
        "scTenifoldNet 1.4",
        "RSpectra 0.16.2",
        "```",
        "",
        "目标基因：",
        "",
        "```text",
        "SASH1",
        "COL1A1",
        "EMP1",
        "MYH11",
        "```",
        "",
        "每个基因分别进行虚拟KO和虚拟OE，输出全部差异调控结果、显著结果、"
        "火山图、Top基因条形图及OE网络对象。",
        "",
        "| 基因 | 扰动 | 状态 | 表达细胞数 | 表达比例 | 显著基因数 |",
        "|---|---|---|---:|---:|---:|",
    ])

    if sc_summary:
        for row in sc_summary:
            expression_pct = row.get("target_expression_pct", "")
            try:
                expression_pct_text = f"{float(expression_pct):.2f}%"
            except (TypeError, ValueError):
                expression_pct_text = expression_pct

            lines.append(
                f"| {row.get('target_gene', '')} "
                f"| {row.get('perturbation', '')} "
                f"| {row.get('status', '')} "
                f"| {row.get('target_expressing_cells', '')} "
                f"| {expression_pct_text} "
                f"| {row.get('significant_genes', '')} |"
            )
    else:
        lines.append("| 未读取到结果 | - | - | - | - | - |")

    lines.extend([
        "",
        "当前结果中，仅MYH11虚拟OE检测到2个FDR显著差异调控基因。",
        "其余分析在当前细胞选择、网络参数和FDR阈值下未达到显著性。",
        "",
        "虚拟KO/OE属于基因调控网络计算模拟，不等同于CRISPR敲除、转染过表达、"
        "动物模型或其他湿实验验证。无FDR显著结果也不等于相应基因没有生物学作用。",
        "",
        "## 大型RDS对象和Git LFS",
        "",
        "仓库中的大型输入文件、Seurat对象、网络对象和RDS分片由Git LFS管理。",
        "",
        "完整RDS对象包括：",
        "",
        "```text",
        "results/objects/12_spatial_tissue_spots_Seurat.rds",
        "results/objects/13_spatial_SASH1_COL1A1_high_annotated_Seurat.rds",
        "results/objects/14_spatial_domain_annotated_paper_style_Seurat.rds",
        "results/objects/17_scTenifold_analysis_input_Seurat.rds",
        "results/objects/17_scTenifold_OE_networks/",
        "```",
        "",
        "采用分片方式保存的对象包括：",
        "",
        "```text",
        "results/objects/11_global_trajectory_Seurat.rds.parts/",
        "results/objects/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds.parts/",
        "results/objects/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds.parts/",
        "results/objects/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds.parts/",
        "```",
        "",
        "重建大型RDS：",
        "",
        "```bash",
        "bash scripts_utils/reconstruct_split_rds_objects.sh",
        "```",
        "",
        "通用分片文件重建：",
        "",
        "```bash",
        "bash scripts_utils/reconstruct_split_files.sh",
        "```",
        "",
        "克隆后下载LFS对象：",
        "",
        "```bash",
        "git lfs install",
        "git clone https://github.com/ssshuoshuo/HNSCC_SASH1_multiomics_reproduction.git",
        "cd HNSCC_SASH1_multiomics_reproduction",
        "git lfs pull",
        "```",
        "",
        "## 可重复性说明",
        "",
        "各阶段尽量保存以下信息：",
        "",
        "- 分析参数和实际输入对象",
        "- PDF或PNG结果图",
        "- CSV统计表",
        "- Seurat或网络RDS对象",
        "- `sessionInfo.txt`",
        "- 输出文件存在性检查表",
        "",
        "部分脚本包含本地绝对路径。换电脑或移动项目后，需要修改脚本顶部的`project_dir`。",
        "",
        "脚本17保留了缺失依赖自动安装逻辑。Apple Silicon macOS首次运行前还需要安装：",
        "",
        "```text",
        "Apple Command Line Tools",
        "GNU Fortran 14.2",
        "```",
        "",
        "已经安装依赖后，可以把脚本17中的：",
        "",
        "```r",
        "install_missing_packages <- TRUE",
        "```",
        "",
        "改为：",
        "",
        "```r",
        "install_missing_packages <- FALSE",
        "```",
        "",
        "以避免每次运行都重复检查依赖。",
        "",
        "## 项目范围和限制",
        "",
        "- 当前仓库重点复现单细胞RNA-seq、空间转录组和虚拟扰动模块。",
        "- 尚未纳入bulk RNA-seq差异分析、机器学习筛选、TCGA预后模型和外部队列验证。",
        "- 恶性细胞注释综合依赖marker、人工判断、cluster结构和CopyKAT。",
        "- 空间结构域为marker-score近似注释。",
        "- 公开空间数据缺少完整H&E图像和明确样本拆分。",
        "- 拟时序或拟时序代理结果不能替代真实时间序列。",
        "- 虚拟KO/OE不能替代真实湿实验。",
        "",
    ])

    lines.extend(build_stage_output_mapping(files))
    lines.extend(build_inventory(files))

    lines.extend([
        "## README自动更新",
        "",
        "新增、删除或重命名Git追踪文件后，可重新运行：",
        "",
        "```bash",
        "python3 scripts_utils/generate_detailed_readme.py",
        "```",
        "",
        "该脚本会重新扫描`git ls-files`并更新本README中的脚本说明、"
        "脚本与输出对应关系以及完整文件清单。",
        "",
    ])

    README_PATH.write_text(
        "\n".join(lines).rstrip() + "\n",
        encoding="utf-8",
    )

    print(f"README已生成：{README_PATH}")
    print(f"Git追踪文件数量：{len(files)}")
    print(f"README行数：{len(lines)}")


if __name__ == "__main__":
    main()
