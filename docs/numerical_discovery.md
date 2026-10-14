# Optional numerical discovery / 可選數值探索

Exact verification uses only Python's standard library. Numerical search is optional and has separate dependencies:

```sh
python3 -m pip install -r code/discovery/requirements.txt
python3 code/srg.py discover --check-inputs
python3 code/srg.py discover --cases 540_147_66_30 --output results/runs/example_540_147_66_30
```

Use a separate Python environment for the numerical dependencies. The recorded search environment used Python 3.12. The requirements preserve the recorded package versions; preparing this release did not install them or rerun optimization.

The six models for 351, 476, 540, 550, 703, and 1344 reconstruct five/six-point count equations and kernels from parameters. The models for 375, 405, 441, and 575 use the retained verified seven-point source models. The three structural cases also use the proved structural information. A fresh solve does not read saved certificate directions or multipliers. Every candidate must pass exact rational reconstruction and a strict contradiction test.

Use full tuples with `--cases`. Omitting that option selects all ten. The output directory must be new; prior results are never overwritten. The run contains the staged engines and model inputs, `run_configuration.json`, and new results under `fresh/`. Staged internal filenames follow the retained research engines; the public repository paths are listed in the source map.

To inspect preparation without optimization or installed numerical dependencies:

```sh
python3 code/srg.py discover --cases 405_132_63_33 --prepare-only --output results/runs/prepared_405_132_63_33
```

A numerical status such as `optimal_inaccurate` describes a search candidate, not a proof. Accepted output needs an exact rational identity or positive margin. New optimization can take substantial time and return different directions or weights from those printed in the paper.

## 繁體中文

精確核驗只需標準函式庫；可選數值探索另使用 `code/discovery/requirements.txt` 中的套件，建議安裝於獨立環境。原記錄的探索環境為 Python 3.12。本次整理未重新安裝套件或執行數值最佳化。

六個五點或六點模型由參數重建；375、405、441、575 四個七點模型使用保留且已核驗的來源模型，三個結構案例另加入已證明的資訊。新求解不讀取保存的方向與權重，候選必須通過精確有理數重建與嚴格矛盾檢查。

`--cases` 使用完整四元參數；省略則選擇十組。`--output` 必須是新資料夾。程式保存暫存引擎、模型輸入、執行設定及 `fresh/` 中的新結果；既有證書不變。暫存引擎保留歷史內部名稱，與公開直觀檔名的關係可查來源對照。

`--check-inputs` 只查語法與資料位置，`--prepare-only` 可在不安裝數值套件、不求解的情況下準備輸入。數值求解器狀態不能當作證明，最終須有精確有理數恆等式或嚴格正餘量。新求解可能耗時，也可能得到與文章不同的方向或權重。
