# Script-to-output map

| Entry point | Output | Input requirement |
|---|---|---|
| `analysis/reproduce_aggregate.py` | Brief Report Table 1 and Table 2 rate components; Figure 1 PNG/SVG/PDF; source trace | Only shipped rounded aggregates |
| `analysis/run_private.py` | Validated execution of one external prepared Python/R/SQL script | Absolute private workspace outside checkout |
| `analysis/prepare_workspace.py` | Portable legacy code and dictionary/parameter setup outside this checkout | User-configured private paths; no records read |
| `analysis/privacy_scan.py` | Privacy checklist and machine-readable scan receipt | This repository; optional gitleaks executable |
| `sql/00_inventory.sql` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `sql/01_eicu_dictionary.sql` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `sql/01_mimic_dictionary.sql` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/02_prepare_maps.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `sql/03_export_eicu.sql` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `sql/03_export_mimic.sql` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/03_generate_pg_exports.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/05_scan_amsterdam.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/05_scan_sicdb.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/05_stage_amsterdam_csv.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/06_probe_and_prepare.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `sql/07_eicu_assessment.sql` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `sql/07_metadata_followup.sql` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `sql/07_reextract_eicu_drugs.sql` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/08_build_stays_drugs.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/08_core.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/09_build_bp.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/09_decode_sicdb_vectorized.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/10_pair_and_counts.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/11_signal_and_audit_sample.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/12_verify_synthetic.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `sql/13_eicu_source_metadata.sql` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `sql/13_mimic_source_metadata.sql` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/13_source_metadata.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/14_fix_mimic_height.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/14_mimic_height_probe.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/15_prepare_missing_check.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/16_complete_aggregates.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/17_finalize_mappings.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/18_final_metadata.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/19_rounding_and_qc.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/20_probe_sicdb_covariates.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/21_correct_sicdb_units_af.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/22_final_aggregate_checks.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/23_quality_summary_fix.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `sql/26_mimic_careunit_probe.sql` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/27_mimic_careunit_bound.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/amsterdam_stage_inner.py` | Cohort, exposure and paired measurement inputs for all tables/figures | Credentialed source/private intermediate files |
| `analysis/legacy/phase1/common.py` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase1/01_prepare.py` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase1/02_rates.py` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase1/03_models.R` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase1/04_recovery_sensitivity.py` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase1/05_figures.py` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase1/06_boundary_labels.py` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase1/07_validate.py` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase1/08_validate_ba.R` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase1/eicu_recovery.py` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `sql/phase1/covariates_eicu.sql` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `sql/phase1/covariates_mimic.sql` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `sql/phase1/recovery_export_medication.sql` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase1/recovery_probe.py` | Primary/false-alarm rates, repeated-measurement agreement and dose model; BR Tables 1–2 and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/common.py` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/01_analyses.py` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/02_pool.R` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/03_validate.py` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/04_figures.py` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/05_complete_metadata.py` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/06_fill_same_minute_prevalence.py` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/08_positional_sensitivity.py` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/standardized.py` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/table1_build.py` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/table1_finalize.py` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `sql/phase2/table1_export_eicu.sql` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `sql/phase2/table1_export_mimic.sql` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `sql/phase2/table1_probe_eicu.sql` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `sql/phase2/table1_probe_mimic.sql` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase2/table1_probe.py` | Deep missed rates, cuff-conditioned probabilities, standardization and baseline characteristics; BR Tables 1–2 and Figure 1 lower panels | Credentialed source/private intermediate files |
| `analysis/legacy/phase3/population_prepare.py` | Population-restriction sensitivity tables in Additional file 1 | Credentialed source/private intermediate files |
| `analysis/legacy/phase3/population_analyse.py` | Population-restriction sensitivity tables in Additional file 1 | Credentialed source/private intermediate files |
| `sql/phase3/population_source_eicu.sql` | Population-restriction sensitivity tables in Additional file 1 | Credentialed source/private intermediate files |
| `sql/phase3/population_source_mimic.sql` | Population-restriction sensitivity tables in Additional file 1 | Credentialed source/private intermediate files |
| `analysis/legacy/phase4/01_source_extract.py` | Catheter-site and pressure-range-bias supplementary tables/figure | Credentialed source/private intermediate files |
| `analysis/legacy/phase4/02_analyse.py` | Catheter-site and pressure-range-bias supplementary tables/figure | Credentialed source/private intermediate files |
| `analysis/legacy/phase4/04_figure.py` | Catheter-site and pressure-range-bias supplementary tables/figure | Credentialed source/private intermediate files |
| `sql/phase4/source_eicu.sql` | Catheter-site and pressure-range-bias supplementary tables/figure | Credentialed source/private intermediate files |
| `sql/phase4/source_mimic.sql` | Catheter-site and pressure-range-bias supplementary tables/figure | Credentialed source/private intermediate files |
| `analysis/legacy/phase5/01_analyse.py` | Common measurement relation; BR Figure 1 upper panels and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase5/04_figure.py` | Common measurement relation; BR Figure 1 upper panels and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase5/audit_points.R` | Common measurement relation; BR Figure 1 upper panels and supplement | Credentialed source/private intermediate files |
| `analysis/legacy/phase6/01_analyse.py` | Common-anchor time-window/selection supplementary tables/figure | Credentialed source/private intermediate files |
| `analysis/legacy/phase6/03_figure.py` | Common-anchor time-window/selection supplementary tables/figure | Credentialed source/private intermediate files |
| `analysis/legacy/icm_rev4/01_analyse.py` | Additional cuff threshold and likelihood-ratio supplementary tables; historical Figure 1 code | Credentialed source/private intermediate files |
| `analysis/legacy/icm_rev4/02_figure.py` | Additional cuff threshold and likelihood-ratio supplementary tables; historical Figure 1 code | Credentialed source/private intermediate files |
| `analysis/legacy/icm_rev4/audit_bootstrap.py` | Additional cuff threshold and likelihood-ratio supplementary tables; historical Figure 1 code | Credentialed source/private intermediate files |

`analysis/legacy/rev8d/sicdb_bootstrap_diagnostics.py` implements the later SICdb percentile/BCa/logit-Wald interval and influence diagnostics. It requires approved local paired observations plus original primary aggregates. Every result is directed to the private workspace, including summaries; none is included in this aggregate release.

Legacy scripts are executed after workspace preparation. Shared `common.py` and standardization/recovery modules are dependencies, not standalone stages. Extraction SQL creates only private outputs in that workspace. Historical reports, manuscript/DOCX generators, environment inventories, literature files, frozen source-review records and package-build scripts were intentionally excluded. Some historical validation/finalization utilities expect those private receipts; they are not dependencies of aggregate reproduction. No patient-data file is an input to the public reproduction path.

The aggregate filenames retain their supplement lineage. `S2_*` supplies overall and dose-stratified error rates; `S3_cuff_conditioned_probabilities.csv` supplies the cuff perspective and lower figure panels; `S5_standardized_missed_hypotension.csv` supplies depth-standardized rates; `phase5_S5_counterfactual.csv` supplies upper figure panels. Source CSV names and original row locations are retained in the released aggregate rows. This metadata identifies summary rows, not patients.

## Additional source and presentation implementations

These historical scripts are prepared into the external workspace and require
authorized source access or the original private aggregate inputs. They are not
part of the public rounded-aggregate reproduction command. Their outputs can
contain restricted counts and must stay in that workspace until independently
reviewed for disclosure.

| Prepared script (relative to workspace/sql) | Implementation and dependency | Outputs |
|---|---|---|
| `phase2/audit_mimic_positional_all.py` | MIMIC waveform-label source query for every existing primary pair; uses `source_query_helpers.py`; PostgreSQL client via `CUFF_PSQL` (default `psql`), database via `CUFF_MIMIC_DATABASE` (default `mimiciv`) | Private all-pair positional-waveform flag file consumed by `phase2/08_positional_sensitivity.py`; private coverage and source receipts |
| `phase2/audit_mimic_positional_all.sql` | Static source query companion; the Python entrypoint materializes the input key file and renders the equivalent SQL | Private waveform observations only |
| `phase3/figures_tables.py` | Original phase 1/2 frozen aggregates, phase 2 flow/baseline and parameters, verified private `phase3/sources/everdac_dose_review.csv` | Full-manuscript Figures 1–3, historical compact letter figure, Tables 1–2, baseline detail and numerical trace |
| `phase3/figures_supplement.py` | Phase 1 repeated-measures agreement/density and fitted model curve outputs | Historical supplement Bland–Altman and dose-association figures and source traces |
| `revision8a/prepare_aggregate_tables.py` | Existing `phase3/revision4/data/*.csv` and `phase3/icm_revision4/analysis/A1_cuff_ge70.csv`, `A2_likelihood_ratios.csv` | Four-database Table 2 source organization, false-alarm display, hypotension-depth distribution, Supplementary Tables 17–18; no model refit |

The static positional SQL is generated by its Python entrypoint; it is supplied
for inspectability, not a second estimator. Historical manuscript figure captions
and precision labels retain their original versions. Current Brief Report
display rules are implemented by `analysis/reproduce_aggregate.py` and do not
change the frozen estimates. The old density plot's pair-count display threshold
is not a patient disclosure rule: no density output is shipped in this repository.
All new plots from private inputs require the same release review as other
private analytical outputs.

## Independent research validation

General source-query, sampling, independent recount, quality-range and bootstrap
verification implementations are also included. See [research validation
code](validation_code.md) for the phase-specific entry points, private input
requirements and non-executable DuckDB query companions. Their historical
sampling and statistical rules remain intact. Source-review records and
patient-specific outputs are not included; none of these checks was rerun as
part of the repository-completeness update.
# Brief Report v5 display additions

| Script | Inputs | Outputs and scope |
|---|---|---|
| `analysis/v5/regenerate_s6.py` | Approved local pair Parquet files and existing full-data reference lines | Dual-threshold retained density bins, whole-database range coverage and Supplementary Figure 1; no estimate refit |
| `analysis/v5/export_figure1.py` | Existing repository display-precision summaries | Four-panel Figure 1 at 170 mm and 600 dpi; no model fitting |

