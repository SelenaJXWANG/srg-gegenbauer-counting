# Nonexistence of Strongly Regular Graphs via Gegenbauer Positivity and SRG Counting

[繁體中文](README_zh-TW.md) · [简体中文](README_zh-CN.md) · [Main article](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_2026-10-14.pdf) · [405 supplement](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_Supplement_405_2026-10-14.pdf)

**Zhen-Qi Liao, Jiaxin Wang, Wei-Hsuan Yu, Jianfeng Cai**

Article date: **2026-10-14**

This repository accompanies a common method for excluding strongly regular graph parameters and studying the relationship between semidefinite programming and strongly regular graph counting. It combines Gegenbauer positivity with exact counting information. Numerical SDP solutions guide candidate combinations; the supplied rational coefficients are reconstructed and checked separately before being used in explicit dual contradiction identities. The printed arguments can be read as mathematics without running a solver. The optional exact verifier independently checks their finite arithmetic.

The work naturally extends Liao–Yu. It gives certificates for the following ten parameter tuples. “Previously excluded” refers to Liao–Yu's result; the present certificates and derivations need not coincide with theirs.

| Parameters $(v,k,\lambda,\mu)$ | Relation to Liao–Yu | Additional structural information in this proof |
|---|---|---|
| 351, 140, 73, 44 | Previously excluded by Liao–Yu | Not needed |
| 375, 102, 45, 21 | Additional case treated here | Used (Section 4) |
| 405, 132, 63, 33 | Additional case treated here | Not needed |
| 441, 88, 35, 13 | Additional case treated here | Used (Section 4) |
| 476, 133, 60, 28 | Additional case treated here | Not needed |
| 540, 147, 66, 30 | Additional case treated here | Not needed |
| 550, 162, 75, 36 | Previously excluded by Liao–Yu | Not needed |
| 575, 112, 45, 16 | Additional case treated here | Used (Section 4) |
| 703, 182, 81, 35 | Previously excluded by Liao–Yu | Not needed |
| 1344, 221, 88, 26 | Previously excluded by Liao–Yu | Not needed |

[arXiv TeX source ZIP](paper/arxiv-sources.zip)

## Read the proofs

- [Main article PDF](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_2026-10-14.pdf) and [editable TeX](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_2026-10-14.tex).
- [405 supplementary PDF](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_Supplement_405_2026-10-14.pdf) and [editable TeX](paper/Nonexistence_of_Strongly_Regular_Graphs_via_Gegenbauer_Positivity_and_SRG_Counting_Supplement_405_2026-10-14.tex).
- [Notation and normalization](docs/notation.md), [code–article alignment](docs/code_article_alignment.md), and [verification guide](docs/verification.md).

The main article states an explicit contradiction for every parameter tuple. The supplement records the long 405 coefficient/function tables and the all-degree five-point comparison. The two 405 records are separated: `ordinary_counting_certificate/` is the proof used by this article; `earlier_structural_certificate/` is a retained alternative with extra zero-count sources. They are not interchangeable and do not count as two parameter exclusions.

The **2026-10-14 certificate snapshot** is publicly available at fixed commit [85cd39b](https://github.com/SelenaJXWANG/srg-gegenbauer-counting/tree/85cd39b74c9ef8a61ab746deb1102af3eeae7e10). The article and supplement now link to this snapshot; their displayed certificate coefficients are unchanged. No DOI or arXiv identifier has been assigned. The fixed-input manifest is [extras/certificate_inputs.sha256](extras/certificate_inputs.sha256), with SHA-256 `58ebb057e5d840a1571d503bf2c3a4fc0d1854fe6c7336615eb1fc7bc0398b9b`. Detailed scope and the certificate entry points are in the [reproducibility record](docs/reproducibility_snapshot.md).

## Exact verification

Python **3.9 or later**, standard library only; no SDP solver or network access is required:

```sh
python3 -S code/srg.py verify
```

Equivalent direct command:

```sh
python3 -S code/verify.py
```

The five stages reconstruct the ten certificates, check the ordinary-counting 405 proof and its all-degree comparison, read the actual short and expanded TeX certificate tables, and compare the article's parameter/structural tables and strict rational bounds. The programs run in a temporary directory and leave the paper and certificate files unchanged. Output is written to `results/latest_verification/`. Do not use Python's `-O` option, which disables assertions.

## Optional numerical discovery

Discovery searches for new certificates; it is separate from verifying the supplied proofs. Its dependencies are in [code/discovery/requirements.txt](code/discovery/requirements.txt).

```sh
python3 code/srg.py discover --check-inputs
python3 code/srg.py discover --cases 540_147_66_30 --output results/runs/example_540_147_66_30
```

Six five/six-point models are rebuilt from parameters. Four seven-point models use retained, verified counting models and, when applicable, proved structural sources. Fresh optimization does not load saved certificate directions or multipliers. See [the discovery guide](docs/numerical_discovery.md) for setup, limitations, and preparation without optimization.

## Repository layout

```text
code/                        Public commands and verification/search engines
results/certificates/        Fixed certificate inputs, grouped by full parameter tuple
results/catalogue/           Shared small-graph types and canonical-label maps
results/release_checks/      Exact verification recorded for this release
results/runs/                New numerical runs (generated and ignored)
paper/                       Article and supplement, PDF and TeX
docs/                        Mathematical conventions and usage guides
extras/                      Source mapping, provenance, and release hashes
```

The organization follows the [Liao–Yu reproduction repository](https://github.com/Wideep886/srg_sdp) at commit `de824a2`, with `paper/` added here. Public filenames describe their purpose. The original exact verification engines and mathematical inputs are byte-preserved behind a temporary path adapter; [source history](extras/source_history.md) and [the machine-readable map](extras/file_mapping.json) identify each input and its origin. There are no archived alternative numerical runs or old article copies in this release.

## Acknowledgement and reuse

The authors developed the mathematical ideas, research design, constraints, and interpretation. GPT (OpenAI) assisted with coding implementations under the authors' direction. The authors take responsibility for the mathematics, computations, and presentation. Detailed technical references are in the article.

Software is licensed under [MIT](LICENSE); the paper and certificates are licensed under [CC BY 4.0](LICENSE-CC-BY-4.0.txt). See [license scope](LICENSES.md) for the directory-level distinction.

[Collaborator upload guide](docs/collaborator_upload.md).

[Reproducibility snapshot and all ten certificate entry points](docs/reproducibility_snapshot.md).

[arXiv submission notes](docs/arxiv_submission.txt) · [Submission-form metadata](docs/arxiv_metadata.txt).

## GitHub collaboration and release

The author list and order are recorded in [CITATION.cff](CITATION.cff).
The [collaboration guide](CONTRIBUTING.md) describes branches, pull requests and
versioned corrections. A [publication guide in Chinese](docs/github_publication_zh.md)
covers repository creation, coauthor access, exact verification, releases and DOI
archiving. [Draft release notes](docs/release_notes_v1.0.0.md) are also provided.

GitHub Actions is configured to run the exact verifier on pushes to `main` and
pull requests. The [initial publication run](https://github.com/SelenaJXWANG/srg-gegenbauer-counting/actions/runs/38571633757) passed all five exact stages on Python 3.9 and 3.12. Each later revision has its own workflow result. Repository
destination is [SelenaJXWANG/srg-gegenbauer-counting](https://github.com/SelenaJXWANG/srg-gegenbauer-counting), a public personal repository.
Coauthor access remains to be configured using the authors' GitHub usernames.
The bundled local exact checks passed; consult the Actions tab for hosted results.
