# smk-metagenomics-demo

基于 [Snakemake](https://snakemake.readthedocs.io/) 的宏基因组物种分类演示流程（Demo）。

输入 NGS 单端、NGS 双端（R1/R2 分文件或 R1/R2 交错存于单文件）或三代（长读长）测序数据，使用 **Kraken2** 进行物种分类学鉴定，使用 **Bracken** 重估物种丰度，最终输出各样本的物种及相对丰度汇总表。

> 本流程仅做物种分类学鉴定（taxonomic profiling），面向演示用途，不包含质控、组装、功能注释等步骤。

**使用模式**：本程序目录只存放流程代码；数据和结果不放在程序目录下。使用时在自己的分析目录中准备数据与配置文件，从分析目录调用程序目录中的 `Snakefile` 执行分析，结果与日志生成在分析目录中。

---

## 1. 目录结构

```text
smk-metagenomics-demo/            # 程序目录（仅代码，运行时不写入任何数据）
├── Snakefile                     # 主流程文件
├── README.md                     # 本文档
├── config/
│   ├── config.yaml               # 配置模板（复制到分析目录后修改使用）
│   └── samples.tsv               # 样本清单模板（复制到分析目录后修改使用）
└── envs/
    ├── kraken2.yaml              # kraken2 conda 环境
    └── bracken.yaml              # bracken conda 环境

~/my_analysis/                    # 示例分析目录（自建，数据与结果都在这里）
├── config.yaml                   # 从程序目录复制并修改的配置
├── samples.tsv                   # 从程序目录复制并修改的样本清单
├── data/*.fq.gz                  # 测序数据
├── results/                      # 运行后自动生成：结果输出
└── logs/                         # 运行后自动生成：日志
```

## 2. 环境要求

- Python 3.11+（Snakemake 8/9 的要求；如只有旧版 Python 可使用 Snakemake 7）
- Snakemake ≥ 7（推荐 ≥ 8）
- Conda / Mamba（用于 `--use-conda` 自动创建 kraken2、bracken 运行环境）

安装 Snakemake（推荐使用 conda）：

```bash
conda create -n smk-meta -c conda-forge -c bioconda snakemake
conda activate smk-meta
```

## 3. 准备 Kraken2 / Bracken 数据库

流程需要一份同时包含 Kraken2 数据库文件和 Bracken kmer 分布文件的数据库目录，例如官方预构建的 standard 数据库（建议放在独立数据库目录，不要放在程序目录或分析目录中）：

```bash
mkdir -p ~/databases/k2_standard
cd ~/databases/k2_standard
# 8GB 精简版标准库（演示推荐，完整库请见索引页；文件名带日期后缀，可到索引页获取最新版本）
curl -LO https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz
tar -xzf k2_standard_08gb_20250402.tar.gz && rm k2_standard_08gb_20250402.tar.gz
```

解压后目录内应包含：

| 文件 | 用途 |
|---|---|
| `hash.k2d`, `opts.k2d`, `taxo.k2d` | Kraken2 数据库 |
| `database{N}mers.kmer_distrib` | Bracken 丰度重估所需（`N` 与读长相关） |

自建数据库的方法请参考 [Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual) 与 [Bracken manual](https://github.com/jenniferlu713/Bracken)。

## 4. 准备分析目录

### 4.1 建立分析目录并复制模板

```bash
mkdir -p ~/my_analysis && cd ~/my_analysis
cp <程序目录>/config/config.yaml  ./config.yaml
cp <程序目录>/config/samples.tsv  ./samples.tsv
mkdir data   # 将 fastq 放入（路径也可指向任意其他位置）
```

以下所有命令均在**分析目录**中执行。

### 4.2 样本清单 `samples.tsv`

制表符分隔，包含 4 列：`sample`、`fastq1`、`fastq2`、`layout`。

| 列 | 说明 |
|---|---|
| `sample` | 样本名（唯一，仅含字母/数字/`_`/`-`） |
| `fastq1` | R1、交错单文件，或单端/长读长 fastq 路径（未压缩或 `.gz` 均可） |
| `fastq2` | R2 fastq 路径，**仅 `layout=paired` 时填写**，其余留空 |
| `layout` | 数据类型：`paired`（NGS 双端，R1/R2 分文件）/ `interleaved`（NGS 双端，R1/R2 交错存于同一文件，流程自动以 `--interleaved-input` 分类）/ `single`（NGS 单端）/ `long`（三代长读长） |

示例（对应程序目录自带模板，请替换为真实路径）：

```text
sample	fastq1	fastq2	layout
demo_pe	data/demo_pe_1.fq.gz	data/demo_pe_2.fq.gz	paired
demo_il	data/demo_il.fq.gz		interleaved
demo_se	data/demo_se.fq.gz		single
demo_ont	data/demo_ont.fq.gz		long
```

注意：

- 列之间必须使用**制表符（Tab）**分隔；`layout=paired` 以外的样本 `fastq2` 列留空但仍需保留分隔符。
- fastq 路径相对于分析目录（执行 snakemake 的目录），**建议使用绝对路径**避免歧义。
- fastq 支持**未压缩**与 **`.gz` 压缩**两种形式，流程对每个样本自动识别（`.gz` 文件自动加 `--gzip-compressed`）；同一 `paired` 样本的 R1/R2 压缩状态需一致，不同样本之间可以混用两种格式。

### 4.3 参数配置 `config.yaml`

| 参数 | 说明 | 默认 |
|---|---|---|
| `samples_file` | 样本清单路径（相对分析目录） | `samples.tsv` |
| `results_dir` | 结果输出目录（相对分析目录，可用绝对路径） | `results` |
| `logs_dir` | 日志目录（相对分析目录，可用绝对路径） | `logs` |
| `kraken2_db` | Kraken2 数据库目录 | 需修改 |
| `bracken_db` | Bracken 数据库目录（通常与 kraken2_db 相同） | 需修改 |
| `level` | 分类层级：`S`=种, `G`=属, `F`=科, `O`=目, `C`=纲, `P`=门（Bracken 不支持 `D`/`K`） | `S` |
| `read_len` | Bracken 读长，必须与库内 `database{N}mers.kmer_distrib` 文件匹配 | `150` |
| `threshold` | Bracken 重估的最小 reads 数阈值 | `10` |
| `threads` | 每个样本 Kraken2 分类线程数 | `8` |

> 关于三代数据：Kraken2 可直接分类任意读长；但 Bracken 的重估依赖数据库中与 `read_len` 匹配的 kmer 分布。若数据库没有长读长对应的分布文件，可降低 `read_len`（如 100/150），长读长场景下丰度估计精度有限，请谨慎解读。

## 5. 运行

在**分析目录**中执行（`PIPELINE_DIR` 替换为程序目录的绝对路径，如 `~/src/smk-metagenomics-demo`）：

```bash
# 预演（dry-run），检查将要执行的任务
snakemake -s $PIPELINE_DIR/Snakefile --configfile config.yaml -n --use-conda

# 正式运行（首次会自动创建 conda 环境）
snakemake -s $PIPELINE_DIR/Snakefile --configfile config.yaml --use-conda --cores 8

# 指定其他分类层级（例如属水平）运行
snakemake -s $PIPELINE_DIR/Snakefile --configfile config.yaml --config level=G --use-conda --cores 8

# 生成流程 DAG 图（可选，需要 graphviz 的 dot 命令，如 conda install graphviz）
snakemake -s $PIPELINE_DIR/Snakefile --configfile config.yaml --dag | dot -Tpdf > dag.pdf
```

要点：

- `-s` 指定程序目录中的 `Snakefile`，`--configfile` 指定**你分析目录里的**配置文件（必须提供）。
- 相对路径（fastq、`samples_file`、`results_dir`、`logs_dir`）均相对于执行命令时所在的目录。
- 结果与日志生成在分析目录的 `results/`、`logs/` 下，程序目录保持只读；同一程序目录可同时服务多个不同的分析目录。

常用选项：

- `--cores N`：总并行核数
- `--use-conda`：自动为每个规则创建独立 conda 环境（推荐始终开启）
- `--rerun-incomplete`：中断续跑
- 删除分析目录中的 `results/` 后重跑可强制重算

## 6. 流程与输出

流程步骤：

```text
fastq (paired/single/long)
        │
        ▼
  kraken2_classify      Kraken2 物种分类
        │
        ▼
  bracken_abundance     Bracken 物种丰度重估
        │
        ▼
  abundance_table       多样本丰度表汇总
```

输出文件（分析目录 `results/` 下）：

| 文件 | 说明 |
|---|---|
| `kraken2/{sample}.k2report` | Kraken2 分类报告（各层级分类单元及 reads 比例） |
| `bracken/{sample}.{level}.bracken` | Bracken 重估后的物种丰度（reads 数与相对丰度） |
| `bracken/{sample}.{level}.bracken.report` | Bracken 重估后的分类报告 |
| `abundance_table.{level}.tsv` | **最终汇总表**：行 = 物种（按平均丰度降序），列 = 各样本 reads 数与相对丰度 |

汇总表格式示例（`level=S`）：

```text
name	taxonomy_id	demo_pe_reads	demo_pe_abundance	demo_se_reads	demo_se_abundance
Escherichia coli	562	8123	0.812300	4321	0.432100
Staphylococcus aureus	1280	1000	0.100000	2100	0.210000
```

日志位于分析目录 `logs/`（每个任务一个文件），任务失败时优先查看对应日志排查。

## 7. 常见问题

**Q: 报错"未指定配置文件"**
流程要求通过 `--configfile` 显式指定你自己分析目录中的配置文件，模板见程序目录 `config/config.yaml`。

**Q: 报错"样本清单 samples.tsv 存在以下问题"**
按提示检查对应行：`layout` 只能取 `paired` / `interleaved` / `single` / `long`；`layout=paired` 时 `fastq2` 必须填写，其余 layout 的 `fastq2` 留空；样本名不能重复且只能含字母/数字/`_`/`-`；各列需用 Tab 分隔（某行列数过多或过少也会报错）。

**Q: 启动时提示"Kraken2 数据库目录不完整"或"找不到 database{N}mers.kmer_distrib"**
数据库路径写错或 `read_len` 与库内 kmer 分布文件不匹配。检查 `kraken2_db`/`bracken_db` 目录（应含 `hash.k2d`/`opts.k2d`/`taxo.k2d` 及 `database{N}mers.kmer_distrib`），并把 `read_len` 改为库中存在的长度（官方库一般含 50/75/100/150/200/250/300）。该检查在 dry-run 阶段即会执行。

**Q: 报错 `layout` 不合法**
`samples.tsv` 的 `layout` 列只能取 `paired` / `single` / `long`。

**Q: 报错 MissingInputException**
`samples.tsv` 中 fastq 路径不存在。相对路径是相对于执行 snakemake 命令时所在的分析目录，建议改用绝对路径。

**Q: 如何只跑单个样本？**
临时编辑 `samples.tsv` 只保留该样本行，跑完后再恢复完整清单。注意 `--until`/`--omit-from` 等选项是按规则名截断流程，不能用来限制样本。

**Q: 多个分析能共用这份流程吗？**
可以。程序目录只读不写，每个分析目录有自己的 `config.yaml`、`samples.tsv`、数据与 `results/`。

## 8. 局限性说明（Demo）

- 未包含数据质控（fastp/fastqc）、宿主去除、去冗余等步骤，正式分析建议在分类前完成质控与宿主过滤。
- 未包含功能注释、多样性分析等下游分析。
- Bracken 对三代长读长的丰度重估精度有限（见 4.3 节说明）。
