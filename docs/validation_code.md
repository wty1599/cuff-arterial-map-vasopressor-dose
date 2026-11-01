# Research validation code

The repository includes independent numerical reconstruction, resampling,
sampling and source-query algorithms as well as the production analysis code.
All are private-workspace templates. Inclusion of code does not mean that these
checks were rerun during repository preparation.

Prepare the external workspace with `analysis/prepare_workspace.py`, install
the analysis requirements, and use `analysis/run_private.py`. Data access,
original intermediate inputs and the corresponding study-stage aggregates are
required. The PostgreSQL source-query wrappers use `CUFF_PSQL` (default `psql`)
and standard PostgreSQL authentication. The queries and source definitions are
unchanged; the original workstation-specific WSL launcher was replaced by the
configured local client.

Pass a script's original command-line arguments after `--`, for example:

```bash
python analysis/run_private.py --workspace /your/private/workspace --script audit_file_databases.py -- amsterdam
```

Insert `--check` before `--` to validate the path and arguments without executing
the script. Argument forwarding applies to Python/R scripts; PostgreSQL options
use the launcher's named arguments.

| Prepared path, relative to `workspace/sql/` | Research validation implemented |
|---|---|
| `audit_source_recheck.py`, `audit_file_databases.py`, `audit_missing_sources.py`, `audit_missing_file_sources.py`, `audit_recheck_tolerances.py` | Independent source reconstruction of blood pressure, drug conversion, treatment intervals and missing-observation checks for the initial cohort |
| `phase1/audit_prepare.py`, `audit_pg_source.py`, `audit_file_source.py`, `audit_compare_clean.py`, `audit_quality.py`, `audit_quality_checks.py`, `recovery_audit.py` | Sample selection, source reconstruction, quality-range checks, explicitly synthetic branch checks and recovered-dose validation |
| `phase2/audit_prepare.py`, `audit_pg_source.py`, `audit_file_source.py`, `audit_compare_clean.py`, `audit_quality.py`, `audit_artifact_evidence.py` | Deep-missed source audit and descriptive waveform/device evidence checks |
| `phase2/stats_review_checks.py`, `stats_review_joint_C.py` | Independent pressure-band counts, standardization and declared patient-bootstrap stream replay |
| `phase2/table1_audit_prepare.py`, `table1_audit.py` | Independent baseline and peak-dose sampling/reconstruction |
| `phase3/population_audit.py` | Source-field and population-restriction count verification |
| `phase4/audit_independent.py` | Independent catheter-site classification, source intervals, pressure-range bias and saved-bootstrap interval checks |
| `phase5/audit_prepare.py`, `audit_points.R`, `audit_compare.py` | Independent R comparison of the common measurement relation and a patient resample, plus saved-replicate checks |
| `phase6/audit_source.py`, `audit_compute.py`, `audit_report.py` | Common-anchor source values, independent SQL rate/interval reconstruction and timing/value-change counts; the last file retains only its numerical validation block |
| `icm_rev4/audit_prepare.py`, `audit_bootstrap.py`, `audit_source_local.py`, `audit_finalize.py` | Additional cuff-threshold and likelihood-ratio recount, bootstrap replay and source-dose comparison; the finalizer retains its general verification blocks |
| `revision8a/independent_numeric_audit.py`, `sicdb_false_alarm_source_audit.py` | Independent frozen-pair recount and bounded SICdb source reconstruction for false alarms |

The path names after the first file in a cell share that cell's directory.
`docs/code_manifest.csv` maps each implementation to its original source and
records the adaptations. The original seeds, sample sizes, error tolerances,
selection rules, conversion formulas and estimands are retained. The SICdb
audit's historical observed frame count is read from the original private S2
aggregate instead of being embedded in distributed code.

PostgreSQL `.sql` templates have the prepared-workspace guard. DuckDB queries
that depend on tables registered by their Python wrapper are retained as
non-executable `.sql.txt` companions. They are source code for inspection, not
standalone commands; the owning Python implementation constructs the private
inputs and executes the corresponding query. The launcher intentionally refuses
to execute a companion text file. No selected patient IDs are embedded in these
companions. The sample tables are private inputs.

The historical ICM finalizer expects
`phase3/icm_revision4/analysis/_restricted/validation_expected_checks.json` in
the private workspace. It contains the original `quality.A1_expected_counts`
and `quality.A2_expected_rounded` fields from the author's private validation
parameters. Those observed check values are not supplied in the repository.
Other validation scripts likewise require their original stage-specific
intermediates and receipts. The public aggregate command has no such dependency.

Patient-level source observations, selected-case files, bootstrap arrays,
comparison rows and all unreviewed outputs must stay in the external workspace.
Manuscript/Word/reference/figure-layout checks and author narrative/package
builders are not research-analysis implementations and are not distributed.
The complete local source inventory records these exclusions individually.
