# Release validation / 發行核驗

The final release must pass `python3 -S code/verify.py` from a clean extracted archive, with all five stages successful and fixed inputs unchanged. The result is recorded under `results/release_checks/`. Mathematical input and engine bytes are compared against their source using `extras/file_mapping.json`. Paper files use the final title-and-date copies.

Other release checks cover command help, syntax, optional discovery input resolution, internal Markdown links, public path readability, absence of stale article copies and private machine paths, and the release SHA-256 manifest. Optional numerical optimization is not rerun as part of packaging.

本版由乾淨壓縮檔解壓後執行五階段精確核驗，確認固定輸入不變；另查入口說明、語法、可選探索輸入、文件連結、直觀檔名、舊稿及私人路徑，並建立 SHA-256 清單。整理資料夾不重跑數值最佳化。
