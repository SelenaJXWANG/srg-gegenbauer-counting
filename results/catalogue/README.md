# Shared graph catalogue / 共用小圖目錄

`graph_types.txt` has rows `id h b`: zero-based catalogue identifier, vertex count, and edge mask. In the integer mask, edge `ij` with `i < j` occupies bit `j(j-1)/2+i`. Catalogue identifiers in JSON are not edge masks. `canonical_graph_map_N_vertices.bin` maps labelled N-vertex masks to canonical catalogue identifiers; the retained verifiers decode the binary arrays explicitly. `five_and_six_point_theory.json` retains the source theory metadata and count-coordinate conventions.

圖目錄每列為 `id h b`：圖型編號、頂點數及邊遮罩。JSON 圖型編號不可直接當邊遮罩。二進位映射將標號小圖轉為目錄標準圖型；實際解碼由保留核驗程式處理。記號及正規化換算另見 `docs/notation.md`。
