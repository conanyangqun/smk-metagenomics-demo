# =====================================================================
# 宏基因组物种分类流程 (Demo)
#   输入: NGS 单端 / NGS 双端 / NGS 双端交错单文件 / 三代(长读长) 数据 (fastq 或 fastq.gz)
#   流程: Kraken2 物种分类 -> Bracken 丰度重估 -> 汇总丰度表
#
# 使用模式: 程序目录(Snakefile 所在目录)仅存放代码；
#           在你自己的分析目录下准备 数据/配置文件 并执行:
#
#   cd /path/to/analysis_dir
#   snakemake -s /path/to/smk-metagenomics-demo/Snakefile \
#             --configfile config.yaml --use-conda --cores 8
#
#   - config.yaml / samples.tsv 模板见程序目录 config/ 下
#   - 配置与样本清单中的相对路径均相对于【执行 snakemake 时所在的目录】
#   - results/ 与 logs/ 生成在执行目录（可用 results_dir/logs_dir 修改）
# =====================================================================

import csv
import os
import re

try:
    from snakemake.exceptions import WorkflowError
except ImportError:  # 兼容旧版本
    class WorkflowError(Exception):
        pass


if not config:
    raise WorkflowError(
        "未指定配置文件。请在你的分析目录下执行:\n"
        "  snakemake -s <本Snakefile路径> --configfile <你的config.yaml> ...\n"
        "配置模板见程序目录 config/config.yaml，样本清单模板 config/samples.tsv"
    )

SAMPLES_FILE = config.get("samples_file", "samples.tsv")
try:
    with open(SAMPLES_FILE, encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []   # 直接解析表头行，文件无数据行时也能准确报告
        SAMPLES_TBL = list(reader)
except FileNotFoundError:
    raise WorkflowError(
        f"找不到样本清单文件: {SAMPLES_FILE}\n"
        "(相对路径相对于执行 snakemake 命令时所在的目录，建议使用绝对路径)"
    )
except IsADirectoryError:
    raise WorkflowError(f"样本清单路径是一个目录而不是文件: {SAMPLES_FILE}")
except PermissionError:
    raise WorkflowError(f"样本清单文件没有读取权限: {SAMPLES_FILE}")
except UnicodeDecodeError:
    raise WorkflowError(
        f"样本清单文件无法按 UTF-8 解码，请确认文件编码（应为 UTF-8 文本）: {SAMPLES_FILE}"
    )

# 校验样本清单: 表头列名 / 样本名唯一且非空 / layout 合法 / paired 必须提供 fastq2
VALID_LAYOUTS = ("paired", "single", "long", "interleaved")
REQUIRED_COLS = ("sample", "fastq1", "fastq2", "layout")
missing_cols = [c for c in REQUIRED_COLS if c not in header]
if missing_cols:
    if not header:
        hint = "文件为空，缺少表头行(表头列名: " + " / ".join(REQUIRED_COLS) + ")"
    else:
        hint = "请确认第 1 行为表头(列名: " + " / ".join(REQUIRED_COLS) + ")，且使用 Tab 分隔"
    raise WorkflowError(
        f"样本清单 {SAMPLES_FILE} 表头缺少列: {missing_cols}\n"
        f"当前表头: {list(header)}\n" + hint
    )
if not SAMPLES_TBL:
    raise WorkflowError(
        f"样本清单 {SAMPLES_FILE} 没有任何样本数据行(表头下至少需要 1 行样本记录)"
    )
problems, seen = [], set()
problems_truncated = False
for i, row in enumerate(SAMPLES_TBL, start=2):   # 第1行为表头
    # DictReader 将超出表头的多余字段放入键 None，需检查是否多敲了 Tab
    extra = row.get(None)
    if extra:
        problems.append(f"第{i}行多出 {len(extra)} 列(检查是否多了一个Tab分隔符): {extra}")
        continue
    # DictReader 对字段不足的行会以 None 填充(键仍存在)，需显式检查 None
    missing = [c for c in REQUIRED_COLS if row.get(c) is None]
    if missing:
        problems.append(f"第{i}行缺少列(检查是否少了一个Tab分隔符): {missing}")
        continue
    s = row["sample"].strip()
    layout = row["layout"].strip()
    fq1 = row["fastq1"].strip()
    fq2 = (row.get("fastq2") or "").strip()
    if not s:
        problems.append(f"第{i}行 sample 为空")
    elif s in seen:
        problems.append(f"第{i}行样本名重复: {s}")
    elif not re.fullmatch(r"[A-Za-z0-9_-]+", s):
        problems.append(f"第{i}行样本名 {s!r} 含非法字符(仅允许字母/数字/_/-)")
    if s:
        seen.add(s)
    if layout not in VALID_LAYOUTS:
        problems.append(f"第{i}行 layout 必须为 {VALID_LAYOUTS} 之一，当前: {layout!r}")
    else:
        if not fq1:
            problems.append(f"第{i}行 fastq1 为空")
        if layout == "paired" and not fq2:
            problems.append(f"第{i}行 layout=paired 但 fastq2 为空")
        if layout in ("single", "long", "interleaved") and fq2:
            problems.append(f"第{i}行 layout={layout} 但 fastq2 已填写(仅 layout=paired 填写 fastq2，其余留空)")
        # kraken2 的 --gzip-compressed 作用于全部输入，混合压缩会误读 gz 文件
        read_list = [fq1, fq2] if layout == "paired" else [fq1]
        gz = [f.endswith(".gz") for f in read_list]
        if any(gz) and not all(gz):
            problems.append(f"第{i}行 fastq 压缩状态不一致（.gz 与非 .gz 混用），请先统一压缩格式")
    if problems and len(problems) >= 20:
        problems_truncated = True
        break
if problems:
    raise WorkflowError(
        f"样本清单 {SAMPLES_FILE} 存在以下问题:\n  - " + "\n  - ".join(problems)
        + ("\n  (问题较多，仅显示前 20 条，修正后重新运行以检查其余问题)" if problems_truncated else "")
    )

SAMPLES = [r["sample"].strip() for r in SAMPLES_TBL]
LAYOUT = {r["sample"].strip(): r["layout"].strip() for r in SAMPLES_TBL}
FQ1 = {r["sample"].strip(): r["fastq1"].strip() for r in SAMPLES_TBL}
FQ2 = {r["sample"].strip(): (r.get("fastq2") or "").strip() for r in SAMPLES_TBL}

def _cfg(key):
    """读取必需配置项，缺失时给出友好报错"""
    if key not in config:
        raise WorkflowError(f"config.yaml 缺少必需参数: {key}")
    return config[key]


LEVELS = config.get("levels", config.get("level"))
if isinstance(LEVELS, str):
    LEVELS = [LEVELS]
if not isinstance(LEVELS, list) or not LEVELS:
    raise WorkflowError(
        "config.yaml 中 levels 必须为非空列表，如 levels: [\"S\", \"G\"]\n"
        "(也兼容旧的单值写法 level: \"S\")"
    )
invalid_levels = [lv for lv in LEVELS
                  if not isinstance(lv, str) or len(lv) != 1 or lv not in "PCOFGS"]
if invalid_levels:
    raise WorkflowError(
        f"config.yaml 中 levels 内的层级必须为 P/C/O/F/G/S 之一，非法值: {invalid_levels}\n"
        "(Bracken 不支持 D/K 水平的丰度重估)"
    )
LEVELS = [lv.strip().upper() for lv in LEVELS]
RES = config.get("results_dir", "results")
LOGS = config.get("logs_dir", "logs")

# --------------------- 数据库路径预校验 ------------------------------
# 在 DAG 构建(dry-run)阶段即检查数据库文件，避免到运行时才报错
K2_DB = _cfg("kraken2_db")
BR_DB = _cfg("bracken_db")
if not isinstance(K2_DB, str) or not K2_DB.strip():
    raise WorkflowError("config.yaml 中 kraken2_db 不能为空，请填写数据库目录路径")
if not isinstance(BR_DB, str) or not BR_DB.strip():
    raise WorkflowError("config.yaml 中 bracken_db 不能为空，请填写数据库目录路径")
READ_LEN = str(_cfg("read_len"))

k2_missing = [f for f in ("hash.k2d", "opts.k2d", "taxo.k2d")
              if not os.path.exists(os.path.join(K2_DB, f))]
if k2_missing:
    raise WorkflowError(
        f"Kraken2 数据库目录不完整或不存在: {K2_DB}\n"
        f"缺少文件: {k2_missing}\n"
        "(相对路径相对于执行 snakemake 命令时所在的目录，建议使用绝对路径)"
    )
if not os.path.exists(os.path.join(BR_DB, f"database{READ_LEN}mers.kmer_distrib")):
    raise WorkflowError(
        f"Bracken 数据库目录 {BR_DB} 中找不到 database{READ_LEN}mers.kmer_distrib\n"
        "请确认 bracken_db 指向含 kmer 分布文件的目录，且 config.yaml 的 read_len 与库内分布文件匹配"
        "(官方库一般含 50/75/100/150/200/250/300)"
    )


THRESHOLD = config.get("threshold", 10)
try:
    THRESHOLD = int(THRESHOLD)
except (TypeError, ValueError):
    raise WorkflowError(
        f"config.yaml 中 threshold 必须为非负整数，当前: {config.get('threshold')!r}"
    )
if THRESHOLD < 0:
    raise WorkflowError(f"config.yaml 中 threshold 必须为非负整数，当前: {THRESHOLD}")


K2_THREADS = config.get("threads", 8)
try:
    K2_THREADS = int(K2_THREADS)
except (TypeError, ValueError):
    raise WorkflowError(
        f"config.yaml 中 threads 必须为正整数，当前: {config.get('threads')!r}"
    )
if K2_THREADS < 1:
    raise WorkflowError(f"config.yaml 中 threads 必须为正整数，当前: {K2_THREADS}")

# 预建 bracken 环境目录: 锚定到本 Snakefile 所在目录(程序目录)的绝对路径，
# 这样在任意分析目录下执行 snakemake 都能正确定位（相对路径会按执行目录解析）
BRACKEN_ENV_DIR = os.path.join(workflow.basedir, "envs", "bracken-2.6.1-osx64")
if not os.path.isdir(BRACKEN_ENV_DIR) or not os.path.isdir(
    os.path.join(BRACKEN_ENV_DIR, "conda-meta")
):
    raise WorkflowError(
        f"找不到预建的 bracken osx-64 环境: {BRACKEN_ENV_DIR}\n"
        "(bracken 3.x 无 macOS 构建，需预建 osx-64 环境。预建方法:\n"
        "  cd <程序目录> && CONDA_SUBDIR=osx-64 conda env create -y \\\n"
        "      -f envs/bracken.yaml -p envs/bracken-2.6.1-osx64 )"
    )


def sample_reads(wc):
    """paired 返回双端读段，其余(interleaved/single/long)返回单个文件"""
    if LAYOUT[wc.sample] == "paired":
        return [FQ1[wc.sample], FQ2[wc.sample]]
    return [FQ1[wc.sample]]


def _gz_flag(wc):
    """全部输入均为 .gz 时启用 --gzip-compressed（gz 混用已在样本清单校验阶段拦截）"""
    return "--gzip-compressed" if all(f.endswith(".gz") for f in sample_reads(wc)) else ""


rule all:
    input:
        expand(f"{RES}/kraken2/{{sample}}.k2report", sample=SAMPLES),
        expand(f"{RES}/bracken/{{sample}}.{{level}}.bracken",
               sample=SAMPLES, level=LEVELS),
        expand(f"{RES}/bracken/{{sample}}.{{level}}.bracken.report",
               sample=SAMPLES, level=LEVELS),
        expand(f"{RES}/abundance_table.{{level}}.tsv", level=LEVELS),


# ------------------------- Kraken2 物种分类 -------------------------
rule kraken2_classify:
    input:
        reads=sample_reads,
    output:
        report=f"{RES}/kraken2/{{sample}}.k2report",
    log:
        f"{LOGS}/kraken2/{{sample}}.log",
    params:
        db=K2_DB,
        paired=lambda wc: "--paired" if LAYOUT[wc.sample] == "paired" else "",
        interleaved=lambda wc: "--interleaved-input" if LAYOUT[wc.sample] == "interleaved" else "",
        gz=_gz_flag,
    threads: K2_THREADS
    conda: "envs/kraken2.yaml"
    shell:
        r"""
        kraken2 --db "{params.db}" --threads {threads} \
            {params.paired} {params.interleaved} {params.gz} \
            --output /dev/null --report "{output.report}" \
            {input.reads:q} 2> "{log}"
        """


# ----------------------- Bracken 丰度重估 ---------------------------
rule bracken_abundance:
    input:
        report=f"{RES}/kraken2/{{sample}}.k2report",
    output:
        abundance=f"{RES}/bracken/{{sample}}.{{level}}.bracken",
        report=f"{RES}/bracken/{{sample}}.{{level}}.bracken.report",
    log:
        f"{LOGS}/bracken/{{sample}}.{{level}}.log",
    params:
        db=BR_DB,
        read_len=READ_LEN,
        threshold=THRESHOLD,
    threads: 1
    # bracken 3.x 无 macOS 构建(osx-64 最高 2.6.1)，此处指向预建的 osx-64 环境；
    # 预建方法: CONDA_SUBDIR=osx-64 conda env create -y \
    #             -f envs/bracken.yaml -p envs/bracken-2.6.1-osx64
    conda: BRACKEN_ENV_DIR
    shell:
        r"""
        bracken -d "{params.db}" -i "{input.report}" \
            -o "{output.abundance}" -w "{output.report}" \
            -l {wildcards.level} -r {params.read_len} \
            -t {params.threshold} > "{log}" 2>&1
        """


# --------------------- 多样本丰度表汇总 ------------------------------
rule abundance_table:
    input:
        expand(f"{RES}/bracken/{{sample}}.{{level}}.bracken",
               sample=SAMPLES, level=LEVELS),
    output:
        f"{RES}/abundance_table.{{level}}.tsv",
    log:
        f"{LOGS}/abundance_table.{{level}}.log",
    run:
        level = str(wildcards.level).upper()
        # 读取各样本 bracken 结果，按 (name, taxonomy_id) 合并为丰度矩阵
        order = []   # taxon 首次出现顺序
        seen_keys = set()   # 配合 order 做去重，避免对 list 做线性查找
        data = {s: {} for s in SAMPLES}

        for s in SAMPLES:
            with open(f"{RES}/bracken/{s}.{level}.bracken") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    key = (row["name"], row["taxonomy_id"])
                    data[s][key] = [
                        int(float(row["new_est_reads"])),
                        float(row["fraction_total_reads"]),
                    ]
                    if key not in seen_keys:
                        seen_keys.add(key)
                        order.append(key)

        with open(str(output), "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            header = ["name", "taxonomy_id"]
            for s in SAMPLES:
                header += [f"{s}_reads", f"{s}_abundance"]
            writer.writerow(header)

            # 按平均丰度降序输出
            def mean_abund(key):
                return sum(data[s].get(key, [0, 0.0])[1] for s in SAMPLES) / len(SAMPLES)

            for name, taxid in sorted(order, key=mean_abund, reverse=True):
                line = [name, taxid]
                for s in SAMPLES:
                    reads, abund = data[s].get((name, taxid), [0, 0.0])
                    line += [reads, f"{abund:.6f}"]
                writer.writerow(line)

        # run 规则的日志文件需主动写入才会生成
        with open(str(log), "w") as lf:
            lf.write(f"样本数: {len(SAMPLES)}\n")
            lf.write(f"分类单元数: {len(order)}\n")
            lf.write(f"输出: {output}\n")
