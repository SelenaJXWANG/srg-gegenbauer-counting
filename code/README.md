# Code entry points / 程式入口

- `srg.py verify`: run all five exact checks.
- `verify.py`: equivalent direct exact entry point.
- `srg.py discover`: optional fresh numerical search using full parameter tuples.
- `repository_layout.py`: hash-checked staging adapter for preserved engines.
- `verification/`: exact mathematical engines and article-data parsers.
- `discovery/`: optional numerical engines and requirements.

Run the public entry points above. Individual retained engine files expect the isolated compatibility tree created by the adapter, rather than the public layout. Their descriptive filenames and original source hashes are documented in `extras/file_mapping.json`. This allows file organization to change without changing audited mathematical algorithms.

請使用上述公開入口。`verification/` 與 `discovery/` 內部保留引擎由轉接器在相容暫存樹執行，不直接假設目前公開檔案配置。來源對照記錄每個直觀檔名與原始雜湊，避免整理名稱時改變數學演算法。
