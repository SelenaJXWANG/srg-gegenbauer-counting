# Manuscript and exact certificates - proposed v1.0.0

Authors: Zhen-Qi Liao, Jiaxin Wang, Wei-Hsuan Yu, Jianfeng Cai.

This release accompanies *Nonexistence of Strongly Regular Graphs via Gegenbauer
Positivity and SRG Counting*, manuscript dated 14 October 2026, and its 405
supplement. It provides exact certificates for the ten parameter tuples listed
in the manuscript, including new counting certificates for four cases previously
excluded by Liao and Yu.

- Article and supplement in PDF and TeX.
- Fixed rational certificates and shared small-graph catalogues.
- Exact verifier for all ten cases and the actual displayed manuscript tables.
- Optional numerical discovery programs, kept separate from exact verification.
- Coauthor citation metadata, checksums and provenance.

Verification:

```sh
python3 -S code/verify.py
```

Expected result: `PASS: five exact verification stages; fixed inputs unchanged.`
No SDP solver or internet access is needed for verification.

The fixed-input manifest SHA-256 is
`58ebb057e5d840a1571d503bf2c3a4fc0d1854fe6c7336615eb1fc7bc0398b9b`.
Code uses MIT; the paper and certificates use CC BY 4.0. The public repository is `SelenaJXWANG/srg-gegenbauer-counting`. These remain draft release notes until a v1.0.0 tag/release is published; check the Actions tab for hosted verification results.
