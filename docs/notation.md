# Notation and normalization / 記號與正規化

## English

The article's notation is authoritative. A graph $F$ has fixed labels, and $N_F$ counts injective induced embeddings of $F$ into the hypothetical strongly regular graph. It counts labelled embeddings, rather than unlabelled vertex subsets.

| Symbol | Meaning |
|---|---|
| $d$ | Dimension of the positive restricted eigenspace used for the spherical embedding |
| $\alpha,\beta$ | Inner products for distinct nonadjacent and adjacent vertices, respectively |
| $H,m$ | A labelled reference graph and its number of vertices |
| $D=d-m$ | Residual dimension when the reference Gram matrix is positive definite |
| $f:2^{[m]}\to\mathbb Q$ | Rational function of the adjacency set to the reference vertices |
| $S(H,\ell,f)$ | Sum over all labelled reference embeddings and ordered exterior pairs, including the diagonal |
| $A N=b$ | The counting system; $A$ does not denote the graph adjacency matrix |
| $\xi$ | Selected original free counts after exact elimination |

The profile order is binary lexicographic. For three references it is `000,001,010,011,100,101,110,111`. Canonical edge strings use `01,02,...,0(m-1),12,...`. The seven-point atlas instead stores edge $ij$, $i<j$, at integer bit $j(j-1)/2+i$. These are different encodings; the verifiers explicitly decode them.

The older five/six-point certificates store

$$X_F=\frac{N_F}{v|\operatorname{Aut}F|},\qquad S_{\rm old}(H,\ell,f)=\frac{S(H,\ell,f)}{v|\operatorname{Aut}H|}.$$

Some retained JSON fields and diagnostic strings write `X_F=N_F/v`. In those historical strings only, `N_F` denotes the old **unlabelled vertex-subset count**. With the current article's labelled `N_F`, the correct conversion is the formula above. These original metadata strings are preserved for provenance; they do not redefine the article's notation.

The paper and printed-table checkers include the corresponding automorphism factors. The seven-point data store

$$y_F=\frac{2^{\binom h2}N_F}{(v)_h},\qquad h=|F|,\quad (v)_h=v(v-1)\cdots(v-h+1).$$

In `results/catalogue/graph_types.txt`, each line is `id h b`: a zero-based catalogue identifier, vertex count, and integer edge mask. JSON fields such as `anchor` and `graph_id` contain **catalogue identifiers**, not edge masks. Look up the matching line before decoding `b`. The supplement prints `(h,b)` directly. For the 405 count substitution, `counting_model.json` lists core identifiers in coordinate order, starting with the singleton. For each structural case, `exact_counting_elimination.json` lists `free_core_ids` in the order of coordinates $x_1,\ldots,x_{110}$; coordinate zero is the constant one.

For the structural cases, omitted impossible profiles are filled with zero. The centered-vector kernels in the stored data differ from the unit-vector kernels by the positive factors $1$, $3-\rho$, and $(3-\rho)^2(D-1)/D$ in degrees $0,1,2$, with $\rho=(k+3)/v$. The article explains the additional positive row scaling. These conventions must be converted together; coefficients from differently normalized certificates cannot be mixed.

## 繁體中文

記號以論文為準。$N_F$ 計算有固定標號的圖 $F$ 至假設強正則圖的單射誘導嵌入數，不是無標號頂點子集數。$d$ 是球面表示維度；$\alpha$ 對應非相鄰頂點，$\beta$ 對應相鄰頂點。$A N=b$ 中的 $A$ 是計數方程矩陣。

$f$ 是定義於參考點鄰接集合的有理數函數。座標採二進位字典序；三個參考點時為 `000,001,010,011,100,101,110,111`。$S$ 對所有標號參考嵌入及有序外部頂點對求和，包含對角項，且不除以參考圖自同構數。

舊五點、六點資料的 $X_F$ 與平均核和依上述公式換算；七點資料則使用 $y_F$。七點 atlas 的整數 bit 順序與正文邊字串不同，程式明確處理兩者的轉換。結構證書省略的不可能 profile 以零補齊；中心化向量的核再依上述正比例因子換成正文單位向量的核。使用證書時，方向、權重與正規化必須成套換算。


保留資料內某些歷史 `normalization` 欄位或核驗訊息寫成 `X_F=N_F/v`；那裡的舊 `N_F` 指**無標號頂點子集數**。換成目前論文的標號嵌入數，正確公式是上面的 $X_F=N_F/(v|\operatorname{Aut}F|)$。保留原始欄位是為了追溯來源，不代表論文改用舊記號。

七點 atlas 每列是 `id h b`，分別為從零起算的圖型編號、頂點數及整數邊遮罩。JSON 的 `anchor`、`graph_id` 是**圖型編號**，不能直接當作邊遮罩解碼；應先查該列的 `b`。補充文件直接列出 $(h,b)$。405 的 `cores` 依座標順序排列，首項是單點圖；結構案例的 `free_core_ids` 依 $x_1,\ldots,x_{110}$ 排列，座標零是常數一。

See [the code–article mapping](code_article_alignment.md) for the exact implementation correspondence. 詳細公式與實作對應見該文件。
