# Software and validation scope

| Component | Verified version | Use |
|---|---|---|
| Python | 3.12 (exact patch in reproduction receipt) | Original analysis and aggregate replot |
| NumPy | 2.3.5 | Original analysis and aggregate replot |
| pandas | 3.0.1 | Original analysis |
| DuckDB | 1.4.3 | File-based extraction and private intermediate queries |
| PyArrow | 26.0.0 | Original private intermediate files |
| SciPy | 1.18.1 | Original analysis |
| Matplotlib | 3.11.2 | Original figures |
| Matplotlib | 3.10.8 | Isolated aggregate reproduction environment |
| R | 4.6.0 (2026-04-24 UCRT) | Models and pooling |
| data.table | 1.18.4 | R data handling |
| lme4 | 2.0.1 | Mixed effects |
| metafor | 5.2.1 | Random effects pooling |
| sandwich | 3.1.1 | Robust standard errors |
| nbstripout | 0.9.1 | Notebook output check; no notebooks selected |
| gitleaks | 8.30.1 | Secret scan |

Original package versions were read from installed package metadata; the original models were not rerun for repository preparation. The aggregate reproduction environment is pinned in `requirements.txt`, and its fresh-clone result is recorded separately. PostgreSQL server versions and full raw-data runtime are not newly certified by this release.

The five study data-quality domains remain source-verification concerns: variable missingness, cross-database rates, ascertainment coverage, conserved cohort counts, and independent source-record review. Existing study source checks were reused as provenance; the repository release adds aggregate source/rounding comparisons, suppressions, syntax validation, a secret/path/identifier scan and an isolated reproduction run. A passing repository check is not a new source-record audit.
