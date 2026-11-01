# Privacy scan

Status: **PASS**. Scanned 152 text files; parsed 64 Python source files. Secret scan: gitleaks 8.30.1 (pass). No patient data were read by this static scanner.

The checks cover secrets/tokens, private keys, absolute workstation paths, database connection strings, email addresses, literal patient-ID filters, forbidden aggregate columns, full-precision aggregate values, file types and notebook outputs. Matched secret values are never written to this report. The machine-readable companion contains rule locations and file hashes.

Notebooks selected: 0. No outputs or execution counts were retained. The documented nbstripout workflow is required if notebooks are added; this release contains none.

Identifiers such as subject_id, stay_id and patientunitstayid occur in extraction SQL and analytical grouping code as schema names. They are not ID observations. Static matches were classified by file type; aggregate CSV columns were separately checked. Dictionary item identifiers are concept keys, not patient identifiers. Actual patient or pair records and literal patient-ID selections are prohibited.

The private denominator audit verified the minimum contributing population before rounding and export. Cells below ten patients and cells without sufficient denominator evidence were withheld. No exact patient/pair counts, marginal count totals, individual-contribution percentages, or bootstrap draws are present. Full-precision ratios were replaced with manuscript-display precision. Cross-table disclosure review included code assertions and captions, not only CSV columns. The same safeguards apply even while this repository is private.

This scan is supplemented by an independent pre-upload review. It does not certify unrestricted rights to the source databases or replace the signed data-use agreements. Original row-level audit files remain outside the repository.
