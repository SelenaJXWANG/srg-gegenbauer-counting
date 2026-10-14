# Exact verification / 精確核驗

The article presents rational contradiction identities derived from the semidefinite dual. Its mathematical proofs are the primary presentation; the programs provide an independent finite-arithmetic check. No numerical optimizer is needed to verify the supplied results.

```sh
python3 -S code/verify.py
```

The five checks are:

1. Reconstruct the retained ten certificates, including graph types, count equations, kernel rows, structural zero counts, completion equations, integer moments, and residual bounds.
2. Reconstruct the article's ordinary-counting seven-point 405 proof, and its separate all-degree five-point feasible relaxation witness.
3. Read the actual main TeX file and check the four short contradiction certificates in labelled counts.
4. Read both TeX files and check every displayed six-point and 405 function/coefficient entry and rational contradiction bound.
5. Compare the ten parameter rows, spectra, dimensions, certificate orders, structural classification, three source inventories, six quartics, and three strict margins with the exact inputs.

The 405 five-point witness is a point of a real relaxation, not a graph construction. The structural cases do not claim that every larger counting/kernel model must require the same additional information. The structural lemmas and the all-degree tail argument are mathematical arguments in the documents; their finite ingredients are checked here. This is not a proof-assistant formalization of every sentence.

Fixed inputs are checked against `extras/file_mapping.json`, staged into an isolated temporary tree, and checked again after the run. The adapter preserves the audited engine bytes and their historical internal paths while giving public files descriptive names. No source is rewritten. `results/latest_verification/` receives generated reports. For full logs, add `--keep-logs`; for another location, add `--reports PATH`.

The release report in `results/release_checks/summary.json` records the actual complete run for this release. SHA-256 values identify inputs; hashes alone do not establish the mathematical claims.

## 繁體中文

文章透過半正定對偶寫出有理數矛盾式。正文中的數學推導是主要呈現；程式另提供有限算術的獨立核對，不需要數值最佳化求解器。

五個階段分別核驗十組證書來源、普通計數的 405 七點證書與所有次數的五點比較、主文四組短證書、兩份 TeX 中所有六點與 405 完整資料，以及參數表、結構來源表、六個四次式與三個嚴格餘量。

405 五點可行見證不是圖的構造。三個結構案例也不宣稱每個更高階模型都必須用到相同資訊。結構引理及尾項論證由文章中的數學推導證明，程式核對其有限成分；這並非全文的形式化證明。

程式依來源對照核對雜湊，在暫存樹中執行保留的核驗引擎，完成後確認固定資料未變。公開檔名皆描述用途；歷史相對路徑只在暫存轉接層使用。新報告寫入 `results/latest_verification/`，本次正式核驗記錄位於 `results/release_checks/`。雜湊用來識別輸入，不能單獨證明數學結論。
