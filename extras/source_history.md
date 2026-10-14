# Source history and file mapping / 來源沿革與檔案對照

The repository layout follows the public Liao–Yu reproduction project at commit `de824a2`: `code/`, parameter-indexed `results/certificates/`, recorded release checks, ignored new runs, and linked English/Traditional Chinese READMEs. This snapshot supplies English, Traditional Chinese, and Simplified Chinese READMEs. This article's PDF and TeX files are additionally collected under `paper/`. The reference project is https://github.com/Wideep886/srg_sdp . Its license has not been copied or assumed.

The source of this preparation is the previously verified article-support package. `file_mapping.json` records each original relative path, its descriptive public path, a temporary compatibility path, and its SHA-256 digest. It contains no private machine path. Mathematical JSON/atlas inputs and verification/search engine bytes are preserved. Duplicate graph catalogues, duplicate degree-two/degree-six counting models, and the shared reconstructed 1344 verifier are stored once publicly and staged to their original required locations.

The descriptive public names do not alter the mathematical engines. `code/repository_layout.py` checks every pinned file and copies it into a temporary tree before execution. This explicit adapter preserves the existing internal imports and relative paths, rather than silently editing audited mathematics during renaming. Three internal manifests are generated for the staged retained inputs. Historical prose-only provenance, old standalone theory TeX, redundant verification summaries and obsolete reports are omitted. Their omission changes neither certificate arithmetic nor algorithms. The current 405 supplement contains the all-degree comparison argument.

The article and supplement are copied with title-and-full-date filenames. The map is refreshed with their final hashes before release. `SHA256SUMS` lists all distributed files except itself. Generated local run directories are outside this release manifest.

## Important distinction: two 405 certificates

`results/certificates/405_132_63_33/ordinary_counting_certificate/` contains the article's later seven-point proof, using ordinary counting equations and Gegenbauer positivity. `earlier_structural_certificate/` retains the older alternative that also uses zero-count sources. The two records have different functions and coefficients and must not be mixed.

## 繁體中文

資料夾配置參考 Liao–Yu 專案 `de824a2`，提供英文、繁體中文與簡體中文三份 README，並新增主文與補充資料的 `paper/`。沒有直接沿用或假設該專案的授權。

來源為先前已核驗的文章重現資料。`file_mapping.json` 記錄原相對路徑、公開直觀檔名、暫存相容路徑與 SHA-256，沒有私人電腦絕對路徑。數學 JSON、圖目錄與核驗／探索演算法逐位元保留；相同資料公開只保存一份，在暫存樹中重建必要的重複位置。

路徑轉接器先核對雜湊，再建立暫存樹，讓檔名整理不必改動已核驗演算法的內部匯入與相對路徑。暫存清單只列實際保留的輸入；過期純文字沿革、舊獨立理論 TeX、重複核驗摘要與舊報告不收錄。這些刪除不影響算術。所有次數比較的論證已在目前 405 補充資料中。

405 普通計數證書與早期含結構零計數的證書分開分類，方向與權重不能交叉使用。文章檔名包含完整標題與日期；正式核驗前會更新對應雜湊。
