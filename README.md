# Same cuff, different meaning: vasopressor dose and arterial hypotension behind an acceptable cuff reading in four ICU databases

Research code and disclosure-controlled aggregate results for a retrospective paired-measurement study in MIMIC-IV, eICU, AmsterdamUMCdb and SICdb.

**Software release:** v1.0.0. This repository contains the research code and approved disclosure-controlled aggregates accompanying the manuscript named above. Manuscript submission, acceptance and publication status are separate from this software release.

The v5 supplementary display uses database-specific disclosure grids. See [v5 display preparation](docs/brief_report_v5.md) for the current figure export code and retained S6 layer. The public release preserves the primary analytical definitions and the aggregate disclosure rules.

The clinical question is how often a cuff reading at or above the treatment threshold accompanies arterial reference hypotension in patients receiving vasopressors with both measurements available. The study does not establish the effect of a cuff-first monitoring strategy or transportability to patients monitored only by cuff.

## Reproduce the two Brief Report tables and Figure 1

This path needs no database credentials. It reads only the rounded, disclosure-controlled CSV files in `outputs/aggregate/`. Python 3.12 was used for the verified environment.

```bash
python -m venv .venv
# Windows PowerShell
.venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python analysis/reproduce_aggregate.py --out reproduced
```

On Linux or macOS, use `.venv/bin/python` in place of `.venv/Scripts/python`. The outputs are `Table1.csv`, `Table2.csv` (also Markdown and HTML), `Figure1.png`, `Figure1.svg`, `Figure1.pdf`, a figure source trace, and a receipt containing input hashes and software versions.

Table 1 contains overall missed, deep-missed and false-alarm rates. Table 2 contains the five dose strata, including crude and depth-standardized missed rates and cuff-conditioned probabilities. Figure 1 is the four-panel residual-hypotension figure used by the Brief Report, **not** the full-length manuscript's cohort flow diagram. It includes observed and common-relation expected probabilities in the top panels, and cuff-interval probabilities in the bottom panels.

Counts are deliberately withheld. The script reproduces rate components and redraws the figure from display-precision summaries; it does not reproduce all manuscript table columns or refit models. Current display corrections are applied, including `Not pooled` for uninterpretable pooled false-alarm estimates. See [aggregate release policy](docs/aggregate_release_policy.md) and [input/output map](docs/script_output_map.md).

## Obtain source data

Each investigator must obtain their own permitted access; this repository grants none.

| Database | Access and documentation |
|---|---|
| MIMIC-IV | [PhysioNet MIMIC-IV](https://physionet.org/content/mimiciv/3.1/): credentialing, training and a data-use agreement |
| eICU | [PhysioNet eICU](https://physionet.org/content/eicu-crd/2.0/): credentialing, training and a data-use agreement |
| AmsterdamUMCdb | [Official access application](https://amsterdammedicaldatascience.nl/#request-access): personal approval and an End User License; original-table data are required by these extraction scripts |
| SICdb | [PhysioNet SICdb](https://physionet.org/content/sicdb/1.0.8/): credentialing, research-project approval and the contributor-review data license |

The [license review](docs/data_use_review.md) records the official terms checked for this repository. MIT applies to the study's original code, not to any source database or third-party software. No raw records, paired observations, patient identifiers, bootstrap draws, source-audit rows, vendor libraries, or credentials are included.

## AmsterdamUMCdb

This study uses **AmsterdamUMCdb v1.0.2 (legacy CSV release)**, with the original `admissions`, `drugitems`, `listitems`, `numericitems` and `processitems` tables. The local admissions-file fingerprints match the archived official download review dated 2026-09-29; the data release number is distinct from the archive portal's record-version number. The original-table schema must not be substituted with a different release without reviewing the extraction definitions.

| Code | AmsterdamUMCdb work |
|---|---|
| `analysis/legacy/05_stage_amsterdam_csv.py`, `analysis/legacy/amsterdam_stage_inner.py` | Stage the original CSV and nested numeric archive in the private workspace |
| `analysis/legacy/05_scan_amsterdam.py`, `analysis/legacy/02_prepare_maps.py` | Scan numeric items and resolve pressure and drug concepts from source dictionaries |
| `analysis/legacy/08_core.py`, `08_build_stays_drugs.py`, `09_build_bp.py`, `10_pair_and_counts.py` in the same directory | Construct the source-specific cohort, treatment episodes, blood-pressure series and pairs |
| `analysis/legacy/phase1/` and `analysis/legacy/phase2/` | Run the shared four-database rate, agreement, pooling and supplementary-result pipeline; see the script map for individual entry points and database-specific exceptions |
| `analysis/legacy/phase4/01_source_extract.py`, `analysis/legacy/phase4/02_analyse.py` | Extract documented arterial catheter sites and estimate the supplementary site comparisons |

Database citation: Thoral PJ, Peppink JM, Driessen RH, et al. Sharing ICU Patient Data Responsibly Under the Society of Critical Care Medicine/European Society of Intensive Care Medicine Joint Data Science Collaboration: The Amsterdam University Medical Centers Database (AmsterdamUMCdb) Example. *Crit Care Med.* 2021;49:e563–e577. [doi:10.1097/CCM.0000000000004916](https://doi.org/10.1097/CCM.0000000000004916).

Code is shared in accordance with the AmsterdamUMCdb End User License. The investigator supplied the applicable requirement to deposit all code behind publications in a public repository no later than publication. This main repository is the release location for the complete study code across all four databases. The authors authorized public release of this repository on 6 November 2026. Source records remain subject to the database agreement, and aggregate results follow the disclosure policy above.

## Run the full analysis with approved data access

The legacy scripts are **templates**. Prepare a workspace outside this checkout before running them:

```powershell
$env:CUFF_WORKSPACE = '/your/private/workspace'
$env:AMSTERDAM_ROOT = '/your/approved/amsterdam/source'
$env:SICDB_ROOT = '/your/approved/sicdb/source'
.venv/Scripts/python analysis/prepare_workspace.py
```

Use forward-slash paths valid on your operating system. The preparer never reads patient data, never overwrites changed files, and refuses a destination inside the repository. It replaces path tokens, installs dictionary mappings and parameter files, and creates private intermediate directories. Do not add the generated workspace to Git.

Install the additional packages in `requirements-analysis.txt`; R and PostgreSQL are also needed. PostgreSQL schemas are `mimiciv_icu`, `mimiciv_hosp` (and available derived tables), and eICU `public`. Use standard PostgreSQL environment/service configuration or a private password file. Do not place credentials in code. Execute `.venv/Scripts/python analysis/run_private.py --workspace /your/private/workspace --database mimiciv --script 03_export_mimic.sql` (or the corresponding eICU script). The launcher validates the external workspace before enabling SQL extraction. Use it for Python/R stages as well; specify the script relative to the prepared `sql/` directory and use `--rscript` for an R executable not on PATH. Unprepared Python/R templates reject execution, and SQL templates reject calls without the launcher's prepared-workspace flag. Database names may be replaced by a local service configuration. File-based analyses require the original AmsterdamUMCdb CSV tables and SICdb version 1.0.8 tables.

Run stages in the following dependency order. Extraction utilities and bounded source checks are listed individually in [the script map](docs/script_output_map.md).

1. **Dictionary and extraction:** `sql/00_inventory.sql`, dictionary SQL, and `sql/02_prepare_maps.py`; review concept mappings. Run the PostgreSQL exports and file-database staging utilities.
2. **Cohort, dose and pairing:** `sql/08_build_stays_drugs.py`, `09_build_bp.py`, `09_decode_sicdb_vectorized.py`, `10_pair_and_counts.py`, followed by the source-unit corrections and aggregate checks. The SICdb decoder supports stored vector observations. Do not infer absence of a concept from one item identifier.
3. **Primary analyses:** `sql/phase1/01_prepare.py`, `02_rates.py`, `03_models.R`; recovery sensitivity and figure scripts follow their aggregate inputs.
4. **Explanatory analyses:** `sql/phase2/01_analyses.py`, `02_pool.R`, the standardization module and metadata/validation utilities; baseline extraction/building is a separate branch. The existing positional-waveform sensitivity uses source flags generated by `sql/phase2/audit_mimic_positional_all.py` before `08_positional_sensitivity.py`; these stages require approved source access and the original private paired observations.
5. **Population restrictions:** `sql/phase3/population_prepare.py`, source SQL, `population_analyse.py`.
6. **Catheter site and bias direction:** `sql/phase4/01_source_extract.py`, `02_analyse.py`, `04_figure.py`.
7. **Common measurement relation:** `sql/phase5/01_analyse.py`, `04_figure.py`.
8. **Time-window and selection decomposition:** `sql/phase6/01_analyse.py`, `03_figure.py`.
9. **Additional cuff threshold and likelihood ratios:** `sql/icm_rev4/01_analyse.py`, `audit_bootstrap.py`. Manuscript-specific authoring files and frozen audit receipts are not distributed. The old figure script requires its historical aggregate-export layout; use `analysis/reproduce_aggregate.py` for this release's Figure 1.

10. **SICdb interval diagnostics:** `sql/rev8d/sicdb_bootstrap_diagnostics.py`; all outputs, including contribution and deletion diagnostics, stay under the external workspace's `restricted/rev8d_bootstrap/`. Existing full-precision primary aggregate files are required for its independent reproduction comparisons.

11. **Full-length manuscript displays:** `sql/phase3/figures_tables.py` and `figures_supplement.py` retain the full-length manuscript figure/table implementations. `sql/revision8a/prepare_aggregate_tables.py` organizes the later four-database tables and supplementary depth/cuff-threshold/likelihood-ratio summaries. These templates require the original private aggregates and produce outputs in the external workspace. The credential-free Brief Report reproduction remains the separate first workflow above.

These raw-data stages require private intermediate files and are not claimed to have been rerun end to end during repository preparation. Source-level checks from the study remain distinct from this repository's syntax, privacy and aggregate-reproduction checks. Inspect the dependency notes before running historical validation scripts that expect private audit receipts.

The script map includes the general independent numerical-recomputation, sampling and source-query algorithms used to verify the study. Run them against their corresponding private stage outputs under the same approved access. Case-specific generated SQL and row-level audit evidence remain in the private workspace; they are not distributed as code or aggregate data.

Place stage-specific Python or R arguments after `--` in the launcher command; they are passed through unchanged. PostgreSQL extraction uses the launcher's own database and script options.

## Methods, reproducibility and limits

The primary pairing is a cuff-anchored arterial median within ±5 minutes during vasopressor treatment. The root random seed is **20261027**, with **2,000 patient-cluster bootstrap replicates**. Main intervals are percentile intervals; pooled proportions use the logit scale, REML and modified Hartung–Knapp inference. Independent database streams used by joint standardization are documented in its parameter file. Original frozen definitions are retained in `config/parameters/`.

eICU dose is missing for **72.6%** of paired observations. Its infusion intervals are inferred with a bounded continuation rule. Dose-stratified results use original known-dose observations; recovery analyses remain sensitivity analyses. Conditioning on arterial hypotension and conditioning on an acceptable cuff reading answer different questions. Selected paired observations do not constitute a cuff-only monitoring cohort.

Cuff-conditioned probabilities, low-pressure-depth standardization, population restrictions, catheter-site/bias analyses, the common measurement relation, common-anchor timing analyses and the later likelihood-ratio/cuff-threshold analyses have distinct post hoc timing. See [analysis chronology](docs/analysis_chronology.md). No post hoc analysis is relabeled as prespecified.

The verified aggregate environment uses Python 3.12, NumPy 2.3.5 and Matplotlib 3.10.8. Original analysis versions and validation scope are in [software versions](docs/software_versions.md). The [privacy scan](docs/privacy_scan.md) and [code manifest](docs/code_manifest.csv) describe the release checks and source adaptation.

## Citation and archive

`CITATION.cff` and `.zenodo.json` provide the software title, manuscript-aligned author list, affiliations, version and license. The versioned GitHub release is the archival source. See [the archive checklist](docs/zenodo_checklist.md) for the publication and verification steps. A reserved DOI is not described as a published archive until the record and its files are publicly accessible.

## Data-use checks before transfer

Read `DATA_USE_RULES.md` and install the local pre-commit and pre-push hooks with `python tools/install_data_use_hooks.py`. The hooks and `tools/package_with_gate.py` block transfer when the current payload fails `tools/small_cell_gate.py`. Run `python tools/test_small_cell_gate.py` for synthetic boundary tests; XLSX and Parquet checks require openpyxl and pyarrow and fail closed if a reader is unavailable.

The withdrawn unsuppressed S6 density export is not distributed. The replacement exporter retains the original grid but releases only bins with at least 10 pairs and at least 10 distinct patients. Full-data bias and limits of agreement are unchanged. Single-patient contribution and deletion outputs are not returned by diagnostic stdout or included in delivery payloads. The private launcher redirects complete stdout and stderr to an external workspace's `restricted/agent_run_logs/` and returns only a status and log location.

JSON schema-clock parameter strings are source field definitions, not patient timestamps. A passing data-format check does not replace review of Word, PDF or image contents; these require the documented figure-source and document checks. Do not bypass hooks or rename fields to avoid a failure.
