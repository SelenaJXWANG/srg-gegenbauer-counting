# Collaboration and changes

The author list and order follow the accompanying manuscript. GitHub commit
authors record implementation contributions; they do not replace the manuscript
author list in `CITATION.cff`.

Create a branch for each change and open a pull request against `main`.
Describe any change to mathematical constraints, normalization, coefficients or
verification scope. Run `python3 -S code/verify.py` and include the result.

Paper PDFs and TeX sources should describe the same manuscript revision. Changes
to displayed coefficients require the corresponding table checks to pass.
Changes to fixed mathematical inputs require an explicit new snapshot and an
updated provenance/manifest record, not merely a replacement checksum.

Keep released tags stable. Publish corrections as a new release with a note
explaining which statements, data and verification results changed. Record
numerical experiments separately from the certificates used in the manuscript.

The optional discovery code is distinct from solver-free verification. A solver
status, a negative floating-point margin or a timeout is not a proof.

Repository membership and release permissions are configured on GitHub by the
repository owner. The intended owner is SelenaJXWANG and repository name is
`srg-gegenbauer-counting`. Preserve the manuscript author order. Code uses MIT;
paper and certificates use CC BY 4.0, as detailed in `LICENSES.md`.
