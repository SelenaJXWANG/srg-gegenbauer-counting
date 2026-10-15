# 通过 Gegenbauer 正性与强正则图计数证明不存在性

[English](README.md) · [繁體中文](README_zh-TW.md) · [正文](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_2026-10-14.pdf) · [405 补充材料](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_Supplement_405_2026-10-14.pdf) · [arXiv:2610.20259](https://arxiv.org/abs/2610.20259)

**Zhen-Qi Liao、Jiaxin Wang、Wei-Hsuan Yu、Jianfeng Cai**

文章日期：**2026-10-14** · **[arXiv:2610.20259v1 [math.CO]](https://arxiv.org/abs/2610.20259v1)**

仓库正文 PDF 与 arXiv v1 的官方 PDF 逐字节一致。正文 TeX、405 补充 TeX 及全部 64 份固定证书／核验输入也与 arXiv 源文件包一致。详见[版本核对记录](paper/README.md)与 [BibTeX 引用](paper/citation.bib)。

本项目展示结合 Gegenbauer 正性与强正则图计数信息的统一方法，用于同时处理多组参数，并探讨半正定规划与强正则图计数之间的联系。数值 SDP 解为候选组合提供指引；本文给出的有理系数经过单独重建与精确核验，再代入对偶形式得到明确的矛盾式。读者可以直接阅读数学推导，无需重新运行求解器。精确核验程序提供对有限算术计算的独立检查。

这项工作是 Liao–Yu 的自然延伸，共处理以下十组参数。“Liao–Yu 已排除”指此前的不存在性结果；本篇的证书与推导不必与此前相同。

| 参数 $(v,k,\lambda,\mu)$ | 与 Liao–Yu 的关系 | 本篇证明使用的额外结构信息 |
|---|---|---|
| 351, 140, 73, 44 | Liao–Yu 已排除 | 不需要 |
| 375, 102, 45, 21 | 本篇另行处理的参数 | 使用（第 4 节） |
| 405, 132, 63, 33 | 本篇另行处理的参数 | 不需要 |
| 441, 88, 35, 13 | 本篇另行处理的参数 | 使用（第 4 节） |
| 476, 133, 60, 28 | 本篇另行处理的参数 | 不需要 |
| 540, 147, 66, 30 | 本篇另行处理的参数 | 不需要 |
| 550, 162, 75, 36 | Liao–Yu 已排除 | 不需要 |
| 575, 112, 45, 16 | 本篇另行处理的参数 | 使用（第 4 节） |
| 703, 182, 81, 35 | Liao–Yu 已排除 | 不需要 |
| 1344, 221, 88, 26 | Liao–Yu 已排除 | 不需要 |

[arXiv TeX 原始文件包](paper/arxiv-sources.zip)

## 阅读证明

- [正文 PDF](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_2026-10-14.pdf) 与 [TeX 源代码](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_2026-10-14.tex)。
- [405 补充 PDF](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_Supplement_405_2026-10-14.pdf) 与 [TeX 源代码](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_Supplement_405_2026-10-14.tex)。
- [记号与归一化](docs/notation.md)、[程序与论文对照](docs/code_article_alignment.md)、[精确核验说明](docs/verification.md)。

正文对每组参数都给出明确矛盾。补充材料收录篇幅较长的 405 函数与系数表，以及五点模型在所有次数上的比较。405 的两份记录分开存放：`ordinary_counting_certificate/` 是本篇使用的普通计数证明；`earlier_structural_certificate/` 是包含额外零计数来源的早期替代证书。两者不能混用，也不算作两组排除结果。

本资料为 **2026-10-14 证书快照**，已公开于固定提交 [85cd39b](https://github.com/SelenaJXWANG/srg-gegenbauer-counting/tree/85cd39b74c9ef8a61ab746deb1102af3eeae7e10)。正文与补充文件已加入公开来源链接，证书数值表保持不变；预印本编号为 [arXiv:2610.20259](https://arxiv.org/abs/2610.20259)，DOI 为 [10.48550/arXiv.2610.20259](https://doi.org/10.48550/arXiv.2610.20259)。固定输入清单为 [extras/certificate_inputs.sha256](extras/certificate_inputs.sha256)，SHA-256 为 `58ebb057e5d840a1571d503bf2c3a4fc0d1854fe6c7336615eb1fc7bc0398b9b`。涵盖范围与证书入口详见[复现记录](docs/reproducibility_snapshot.md)。

## 精确核验

只需 **Python 3.9 及以上版本**的标准库，无需半正定规划求解器或网络连接：

```sh
python3 -S code/srg.py verify
```

也可以直接运行：

```sh
python3 -S code/verify.py
```

五个阶段依次重建十组证书，核对普通计数的 405 证明与所有次数的比较，读取实际 TeX 中的短证书与完整证书表，再比对文章的参数表、结构来源与严格有理数下界。程序在临时目录中运行，论文与证书原文件保持不变；结果写入 `results/latest_verification/`。请勿使用会关闭断言的 Python `-O` 选项。

## 可选数值探索

数值探索用于寻找新证书，与已有证明的精确核验分开。所需软件包列于 [code/discovery/requirements.txt](code/discovery/requirements.txt)。

```sh
python3 code/srg.py discover --check-inputs
python3 code/srg.py discover --cases 540_147_66_30 --output results/runs/example_540_147_66_30
```

六组五点或六点模型从参数重建；四组七点模型则使用保留且已核验的计数模型，必要时加入已证明的结构信息。新求解不读取已有证书的方向或权重。安装、限制以及不求解时的准备方式见[数值探索说明](docs/numerical_discovery.md)。

## 文件夹结构

```text
code/                        公开入口、核验与探索程序
results/certificates/        固定证书输入，按完整四元参数分类
results/catalogue/           共用小图目录与标准标号映射
results/release_checks/      本次版本实际完成的精确核验
results/runs/                新数值运行结果（生成后不纳入版本控制）
paper/                       正文与补充材料的 PDF、TeX
docs/                        数学约定与使用说明
extras/                      来源对照、沿革与文件哈希
```

结构参考 [Liao–Yu 的复现项目](https://github.com/Wideep886/srg_sdp)（`de824a2`），另加入 `paper/`。公开文件名直接说明用途。为避免重命名影响已经核验的数学算法，原核验引擎与数学数据逐字节保留，由临时路径适配器处理历史相对路径；[来源沿革](extras/source_history.md)和[机器可读对照](extras/file_mapping.json)列出来源及其对应关系。本版本不包含旧文章或其他数值运行的替代结果。

## 致谢与使用

数学想法、研究设计、约束选择及结果解释由作者提出；GPT（OpenAI）在作者指导下协助程序实现。作者对数学、计算与呈现负责；技术引用详见文章。

代码采用 [MIT](LICENSE)；论文与证书采用 [CC BY 4.0](LICENSE-CC-BY-4.0.txt)。各目录的适用范围见 [许可证说明](LICENSES.md)。

[共同作者上传说明](docs/collaborator_upload.md)。

[复现快照与十组参数的证书入口](docs/reproducibility_snapshot.md)。

[arXiv 投稿说明](docs/arxiv_submission.txt) · [投稿表单文字](docs/arxiv_metadata.txt)。
