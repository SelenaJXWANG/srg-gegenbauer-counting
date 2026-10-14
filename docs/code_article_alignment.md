# Code–article alignment / 程式與論文對照

## English

The source review checked conventions and formulas as well as executing the certificates. It found no change needed to the retained mathematical inputs. The corresponding equations are identified by their subject below, so ordinary changes in page or equation numbers do not change this mapping.

| Article calculation | Implementation and conversion |
|---|---|
| Labelled induced counts and ordinary extensions | `code/discovery/build_counting_kernel_model.py` constructs one coefficient per labelled extension and the totals $\sum N_F/|\operatorname{Aut}F|=\binom vh$. `verify_certificate.py` reconstructs them independently. |
| Older five/six-point count coordinates | `code/verification/verify_raw_counting_certificate.py` uses $X_F=N_F/(v|\operatorname{Aut}F|)$ and $S_{\rm old}=S/(v|\operatorname{Aut}H|)$. Its child/anchor automorphism ratios implement this conversion. |
| Printed short certificates | `code/verification/verify_displayed_tables.py` converts both dependent and free coordinates back to labelled counts, reads all 19 printed functions and 20 coefficients, and reconstructs the identities. The printed five-function 540 proof is the alternative recovered certificate, not the four-function certificate retained in the older JSON. |
| Printed six-point certificates | `code/verification/verify_expanded_certificates.py` checks all 33 functions and 48 integer multipliers. It verifies $\omega_i/L_v=w_i^{\rm old}/|\operatorname{Aut}H_i|$ and the analogous count factor; the article then uses $\omega_i/(vL_v)$. |
| Unit-vector kernels | The residual dimension is $D=d-m$; degree two is $(Dz^2-hh')/(D-1)$ and its diagonal is $h^2$. The raw, pure-405 and fresh six-case verifiers use this convention. |
| Seven-point count scaling | The stored coordinates are $y_F=2^{\binom h2}N_F/(v)_h$. A profile-kernel row uses a single-vertex scale $2^{\binom{m+1}2}/(v)_{m+1}$ and a distinct-pair/single ratio $(v-m-1)/2^{m+1}$. No anchor automorphism divisor occurs. |
| Structural kernels | `code/verification/verify_structural_kernel_completion.py` uses centered vectors $u=\sqrt{3-\rho}\,p$ and kernels $1,z,z^2-hh'/D$. Multiply the article's normalized energies by $\chi_\ell=1,3-\rho,(3-\rho)^2(D-1)/D$, respectively. Omitted impossible profiles are zero-filled, and proved zero graph counts are eliminated. |
| Structural count elimination | `verify_structural_counting_equations.py` rebuilds the ordinary deletion equations, 58 zero-count witnesses and 76 completion equations with 273 witnesses. Identity free rows and an exact rank argument establish all 110 free counts plus the constant. |
| Integer common-neighbour quartics | `verify_integer_moments_and_margins.py` expands falling-factorial moments and verifies the adjacent-pair roots of every active quartic. Their six printed instances are also read by `verify_article_metadata.py`. |
| Residual payment | The seven-point sources have sign $R_i\le0$. The exact margin is $c_0+\sum_{c_j<0}c_jU_j>0$. Every leftover coefficient is paid using a proved count bound; no floating residual is discarded. |

### The two 405 records

The main article and supplement use `results/certificates/405_132_63_33/ordinary_counting_certificate/exact_certificate.json`: 67 count sources and 107 kernel sources (54 degree zero, 50 degree one, three degree two). Its only equalities are ordinary counting equations. Its printed coefficient is the stored `weight/scale`, already reduced; do not divide by the row scale again.

The earlier `results/certificates/405_132_63_33/earlier_structural_certificate/raw_certificate.json` has additional zero-count sources. It is retained for the older replay, but its functions and weights must not be substituted into the article's later proof. Both exclude the same parameters.

The independent counting-model audit reconstructs all 1,252 graph types, 4,135 ordinary extension equations and 1,078 triangular eliminations. The all-degree five-point comparison separately checks a rational real-count witness and a rigorous tail bound. It is not a graph construction or a claim about the intermediate six-point level.

### Encodings and source identifiers

See [notation](notation.md): profile vectors are binary lexicographic; string edge codes and seven-point integer masks have different bit orders. Seven-point JSON `anchor` and `graph_id` fields contain catalogue identifiers, not edge masks. Resolve the `(id,h,b)` atlas row first. Free coordinate order comes from `cores` or `free_core_ids`, not sorted edge masks.

Some original JSON strings say `X_F=N_F/v` using the older meaning “unlabelled subset count” for `N_F`. The current article's labelled count requires the automorphism factor above. Those source strings are preserved, and the conversion is explicit here.

### What is checked automatically

`code/verify.py` performs five stages: the retained ten-case replay; the ordinary-counting 405 and all-degree comparison; actual short certificate tables; actual expanded certificate tables; and article metadata alignment. The last stage checks all ten parameter rows, restricted eigenvalues, dimensions, displayed orders, additional-structure classification, three structural source inventories, six quartics and the three strict margin claims. It complements the source reconstruction; it does not automatically prove every sentence or replace the structural theorem's proof.

The three structural inventories, before residual count-bound terms, are respectively $(13,97,1)$, $(21,88,2)$ and $(60,48,3)$ for counts, kernels and quartics. The strict margins exceed $7/10^9$, $9/10^7$ and $3/10^6$. Positive stored source weights are divided by the positive stored reduced-row normalization exactly once.

### Optional fresh discovery

The six five/six-point models rebuild ordinary counts and kernels from their parameters. Their numerical dual eigenvectors only propose directions; rational reconstruction accepts a certificate only after obtaining nonnegative exact coefficients and the identity $-1$.

The four seven-point models use retained verified count/elimination and structural data, then perform fresh optimization. The 405 model loads neither forbidden counts nor integer-moment sources. Every selected numerical row is replaced by its exact rational reconstruction, and a strictly positive margin after residual payment is mandatory. New runs can produce different functions, source counts and margins. The curated folder includes the printed certificates and discovery implementation; archived alternative numerical-run outputs are not required for exact replay and are not bundled here.

## 繁體中文

本次檢查不只重播證書，也直接比對數學公式、座標與實作。未發現需要修改既存數學輸入的差異。

- **標號計數與自同構因子：**新五點、六點模型直接使用論文的標號誘導嵌入數。舊資料使用 $X_F=N_F/(v|\operatorname{Aut}F|)$，核和另除以 $v|\operatorname{Aut}H|$。印出資料核驗程式同時換算被消去的計數與自由座標，因此不會漏掉自同構因子。舊 JSON 內 `X_F=N_F/v` 的 `N_F` 指當時的無標號子集數。
- **核與尺度：**所有殘餘維度均為 $D=d-m$。405 使用單位球面向量的正規化核；三個結構案例使用中心化向量，所以次數零、一、二分別相差正因子 $1$、$3-\rho$、$(3-\rho)^2(D-1)/D$。七點的 $y_F$ 尺度與有序點對係數依上表換算。
- **編碼：**方向依二進位字典序排列；舊邊字串與七點整數遮罩的 bit 順序不同。七點 JSON 的 `anchor`、`graph_id` 是圖型編號，必須先查 atlas 的 `(id,h,b)`，不能把編號直接當邊遮罩。自由座標順序取自 `cores` 或 `free_core_ids`。
- **405 的兩份證書：**論文採用普通計數的後期證書，含 67 個計數來源及 107 個核來源，其中次數分布為 54、50、3。補充文件已將 `weight/scale` 印成完整係數，不可再除一次尺度。早期 405 含額外零計數來源，保留供舊流程重播，不能與後期方向或係數混用，也不另算一組結果。
- **結構資料：**程式核對 58 個零計數見證、76 個補全方程與 273 個見證，並證明消去後留下 110 個自由計數與常數。選用的六個四次式及來源數量均和正文相符。
- **矛盾符號：**七點來源採 $R_i\le0$；正組合的殘餘係數以計數上界逐項處理，得到嚴格正的有理數常數。沒有忽略任何浮點殘差。

頂層精確核驗有五個階段。除了獨立重建來源與實際印出的完整證書，新增對照程式也讀取十列參數、特徵值、維度、階數及結構分類，並核對三個結構來源表、六個四次式與嚴格下界。這是指定資料的自動核對，不能取代結構定理的數學證明，也不宣稱自動驗證文章每句敘述。

可選數值發現流程的六個五點或六點模型從參數重建；四個七點模型使用保留且已核驗的來源與消去資料，再重新求解。數值解只提供候選方向，必須通過精確有理數重建及嚴格矛盾後才接受。重新執行可能產生不同方向與權重；本資料夾收錄論文證書及發現程式，沒有混入其他數值執行的替代結果。
