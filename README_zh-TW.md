# 透過 Gegenbauer 正性與強正則圖計數證明不存在性

[English](README.md) · [简体中文](README_zh-CN.md) · [主文](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_2026-10-14.pdf) · [405 補充資料](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_Supplement_405_2026-10-14.pdf) · [arXiv:2610.20259](https://arxiv.org/abs/2610.20259)

**Zhen-Qi Liao、Jiaxin Wang、Wei-Hsuan Yu、Jianfeng Cai**

文章日期：**2026-10-14** · **[arXiv:2610.20259v1 [math.CO]](https://arxiv.org/abs/2610.20259v1)**

倉庫主文 PDF 與 arXiv v1 的官方 PDF 逐位元一致。主文 TeX、405 補充 TeX 及全部 64 份固定證書／核驗輸入也與 arXiv 原始檔包一致。詳見[版本核對記錄](paper/README.md)與 [BibTeX 引用](paper/citation.bib)。

本專案展示結合 Gegenbauer 正性與強正則圖計數資訊的共同方法，用來同時處理多組參數，並探討半正定規劃與強正則圖計數的關聯。數值 SDP 解為候選組合提供指引；本文給出的有理數係數經過另行重建與精確核驗，再代入對偶形式，得到明確的矛盾式。讀者可以直接閱讀數學推導，不需要重新執行求解器。精確核驗程式提供有限算術的獨立核對。

這項工作是 Liao–Yu 的自然延伸，共處理以下十組參數。「Liao–Yu 已排除」指先前的不存在性結果；本篇的證書與推導不必與先前相同。

| 參數 $(v,k,\lambda,\mu)$ | 與 Liao–Yu 的關係 | 本篇證明所用的額外結構資訊 |
|---|---|---|
| 351, 140, 73, 44 | Liao–Yu 已排除 | 不需要 |
| 375, 102, 45, 21 | 本篇另處理的參數 | 使用（第 4 節） |
| 405, 132, 63, 33 | 本篇另處理的參數 | 不需要 |
| 441, 88, 35, 13 | 本篇另處理的參數 | 使用（第 4 節） |
| 476, 133, 60, 28 | 本篇另處理的參數 | 不需要 |
| 540, 147, 66, 30 | 本篇另處理的參數 | 不需要 |
| 550, 162, 75, 36 | Liao–Yu 已排除 | 不需要 |
| 575, 112, 45, 16 | 本篇另處理的參數 | 使用（第 4 節） |
| 703, 182, 81, 35 | Liao–Yu 已排除 | 不需要 |
| 1344, 221, 88, 26 | Liao–Yu 已排除 | 不需要 |

[arXiv TeX 原始檔包](paper/arxiv-sources.zip)

## 閱讀證明

- [主文 PDF](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_2026-10-14.pdf) 與 [TeX 原始碼](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_2026-10-14.tex)。
- [405 補充 PDF](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_Supplement_405_2026-10-14.pdf) 與 [TeX 原始碼](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_Supplement_405_2026-10-14.tex)。
- [記號與正規化](docs/notation.md)、[程式與論文對照](docs/code_article_alignment.md)、[精確核驗說明](docs/verification.md)。

主文對每組參數都給出明確矛盾。補充資料收錄篇幅較長的 405 函數與係數表，以及五點模型的所有次數比較。405 的兩份記錄分開存放：`ordinary_counting_certificate/` 是本篇使用的普通計數證明；`earlier_structural_certificate/` 是含額外零計數來源的早期替代證書。兩者不能混用，也不另算兩組排除結果。

本資料為 **2026-10-14 證書快照**，已公開於固定提交 [85cd39b](https://github.com/SelenaJXWANG/srg-gegenbauer-counting/tree/85cd39b74c9ef8a61ab746deb1102af3eeae7e10)。主文與補充文件已加入公開來源連結，證書數值表維持不變；預印本編號為 [arXiv:2610.20259](https://arxiv.org/abs/2610.20259)，DOI 為 [10.48550/arXiv.2610.20259](https://doi.org/10.48550/arXiv.2610.20259)。固定輸入清單是 [extras/certificate_inputs.sha256](extras/certificate_inputs.sha256)，SHA-256 為 `58ebb057e5d840a1571d503bf2c3a4fc0d1854fe6c7336615eb1fc7bc0398b9b`。涵蓋範圍及證書入口詳見[重現記錄](docs/reproducibility_snapshot.md)。

## 精確核驗

只需 **Python 3.9 以上版本**的標準函式庫，不需半正定規劃求解器或網路：

```sh
python3 -S code/srg.py verify
```

也可以直接執行：

```sh
python3 -S code/verify.py
```

五個階段依序重建十組證書，核對普通計數的 405 證明與所有次數比較，讀取實際 TeX 中的短證書與完整證書表，再比對文章的參數表、結構來源與嚴格有理數下界。程式在暫存目錄執行，論文與證書原檔不變；結果寫入 `results/latest_verification/`。請勿使用會關閉斷言的 Python `-O` 選項。

## 可選數值探索

數值探索用於尋找新證書，與既有證明的精確核驗分開。所需套件列於 [code/discovery/requirements.txt](code/discovery/requirements.txt)。

```sh
python3 code/srg.py discover --check-inputs
python3 code/srg.py discover --cases 540_147_66_30 --output results/runs/example_540_147_66_30
```

六組五點或六點模型從參數重建；四組七點模型則使用保留且已核驗的計數模型，必要時加上已證明的結構資訊。新求解不讀取既有證書的方向或權重。安裝、限制與不求解的準備方式見[數值探索說明](docs/numerical_discovery.md)。

## 資料夾配置

```text
code/                        公開入口、核驗與探索程式
results/certificates/        固定證書輸入，以完整四元參數分類
results/catalogue/           共用小圖目錄與標準標號映射
results/release_checks/      本次版本實際完成的精確核驗
results/runs/                新數值執行結果（產生後不納入版本控制）
paper/                       主文與補充資料的 PDF、TeX
docs/                        數學慣例與使用說明
extras/                      來源對照、沿革與檔案雜湊
```

配置參考 [Liao–Yu 的重現專案](https://github.com/Wideep886/srg_sdp)（`de824a2`），另加入 `paper/`。公開檔名直接描述用途。為避免重新命名影響已核驗的數學演算法，原核驗引擎與數學資料逐位元保留，由暫存路徑轉接器處理歷史相對路徑；[來源沿革](extras/source_history.md)及[機器可讀對照](extras/file_mapping.json)列出來源與對應。本版本沒有混入舊文章或其他數值執行的替代結果。

## 致謝與使用

數學想法、研究設計、限制式選擇及結果詮釋由作者發展；GPT（OpenAI）在作者指導下協助程式實作。作者對數學、計算與呈現負責；技術引用詳見文章。

程式碼採用 [MIT](LICENSE)；論文與證書採用 [CC BY 4.0](LICENSE-CC-BY-4.0.txt)。各目錄的適用範圍見 [授權說明](LICENSES.md)。

[共同作者上傳說明](docs/collaborator_upload.md)。

[重現快照與十組參數的證書入口](docs/reproducibility_snapshot.md)。

[arXiv 投稿說明](docs/arxiv_submission.txt) · [投稿欄位文字](docs/arxiv_metadata.txt)。
