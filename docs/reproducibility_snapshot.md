# Reproducibility snapshot / 重現資料快照

This is the **2026-10-14 local certificate snapshot**, a local pre-release. No public repository revision or DOI has yet been assigned. The fixed certificate material is identified by [certificate_inputs.sha256](../extras/certificate_inputs.sha256). Its SHA-256 digest is:

```text
58ebb057e5d840a1571d503bf2c3a4fc0d1854fe6c7336615eb1fc7bc0398b9b
```

The manifest lists **64 fixed files**: mathematical JSON and graph-atlas data, exact verification engines, and the exact command/path adapter. It excludes paper files, documentation, generated reports, and optional discovery code. The article can print this digest without a circular dependency on its own PDF or TeX. The complete distribution, including papers, has the separate release checksum list `extras/SHA256SUMS`.

## Entry points for all ten parameter tuples

The common directory prefix is `results/certificates/`. The machine-readable [certificate index](../extras/certificate_index.json) gives full paths and the corresponding presentation checkers.

| Parameter directory | Entry used by the article or its reconstruction |
|---|---|
| `351_140_73_44/` | `raw_certificate.json`; labelled-count conversion checked against the printed tables |
| `375_102_45_21/` | `exact_certificate.json` and accompanying structural-source JSON |
| `405_132_63_33/` | `ordinary_counting_certificate/exact_certificate.json` and `counting_model.json` |
| `441_88_35_13/` | `exact_certificate.json` and accompanying structural-source JSON |
| `476_133_60_28/` | `raw_certificate.json`; labelled-count conversion checked against the printed tables |
| `540_147_66_30/` | Actual short-certificate tables in the main TeX; independently reconstructed by the short-table checker |
| `550_162_75_36/` | Actual short-certificate tables in the main TeX; independently reconstructed by the short-table checker |
| `575_112_45_16/` | `exact_certificate.json` and accompanying structural-source JSON |
| `703_182_81_35/` | Actual short-certificate tables in the main TeX; independently reconstructed by the short-table checker |
| `1344_221_88_26/` | Actual main-TeX tables and the construction in `code/verification/verify_1344_certificate.py`; no separate certificate JSON is required |

The short-table checker is `code/verification/verify_displayed_tables.py`. It reads the actual rational functions and weights from the main TeX. This matters because the earlier retained JSON can describe a different valid certificate, notably for 540. For 405, the earlier structural certificate and the article's ordinary-counting certificate also remain separate.

## Exact command and expected result

Use Python 3.9 or later, standard library only:

```sh
python3 -S code/verify.py
```

The expected result is five successful stages followed by:

```text
PASS: five exact verification stages; fixed inputs unchanged.
```

The stages check the ten-case sources, the ordinary-counting 405 proof and all-degree comparison, actual short TeX certificates, actual expanded TeX certificates, and the article's parameter/structural tables and strict bounds. The exact identities and positivity arguments in the paper are the mathematical proofs; this command independently checks the listed finite constructions and arithmetic. No SDP solve is required.

Optional fresh search uses `python3 code/srg.py discover` and the dependencies in `code/discovery/requirements.txt`. It proposes new certificates, which may differ from the supplied functions or weights. A new numerical run is not the verification command for the certificate printed in this article.

## 繁體中文

本資料稱為 **2026-10-14 本機證書快照**，屬尚未發布的本機版本，沒有虛構公開專案版本或 DOI。上述 SHA-256 識別固定證書輸入清單，涵蓋 64 份數學資料、精確核驗引擎及入口／路徑轉接器；不涵蓋論文、文件、產生的報告或可選數值搜尋。因而文章可印出此雜湊，不會形成 PDF 自己雜湊自己的循環。整份發行資料另有 `extras/SHA256SUMS`。

十組參數的入口列於上表與機器可讀索引。540、550、703、1344 的短證書直接從實際主文 TeX 讀取，不能只把較早 JSON 當成本文係數；1344 的獨立構造就在專用核驗程式，並非遺漏 JSON。405 普通計數證書與早期結構證書也不能混用。

執行 `python3 -S code/verify.py` 只需 Python 3.9 以上標準函式庫，預期五階段全部 PASS，並確認固定輸入未變。這是對文章指定的有限構造與精確算術進行獨立核對，不需要重新求解 SDP。文章中的恆等式與正性論證才是數學證明。

可選新搜尋使用獨立數值套件，可能找到不同方向與權重；它不等同於核驗本文已印出的證書。
