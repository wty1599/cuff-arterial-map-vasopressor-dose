# Brief Report v5 display preparation

The primary estimands, dose definitions, patient-cluster bootstrap and frozen
result estimates are unchanged. Manuscript-specific Word files, internal number
traces and author declarations are not distributed in this repository.

## Bland–Altman display layer

The v5 supplementary density display uses 2×2 mmHg bins for MIMIC-IV and eICU.
For AmsterdamUMCdb and SICdb, 5×5 and 10×10 mmHg candidates were compared using
a rule recorded before regeneration. The smallest width covering at least 80%
of the full-data central 95% difference interval was retained. Coverage is the
union length of visible retained difference-bin intervals, divided by the full
central-interval length. Both supplementary databases use 10×10 mmHg bins.

Every published bin contains at least 10 paired measurements **and** 10 distinct
patients. Other bins are entirely omitted, including coordinates and counts.
Bias and repeated-measurement limits of agreement reuse the complete-data
estimates. Suppression changes the graphical layer only.

Approved layers and coverage summaries are in `outputs/aggregate/v5_s6/`.
To rebuild with locally approved access and existing paired files:

```bash
python analysis/v5/regenerate_s6.py --pairs-dir /approved/private/pairs \
  --reference-csv outputs/aggregate/v5_s6/S6_reference_lines.csv \
  --output /approved/private/display-output
```

Install the analysis dependencies, including pandas and PyArrow, before running
this private-data route. It exports retained aggregate bins only. Original pair
records and unsuppressed bins are not included in the repository.

## Brief Report Figure 1

```bash
python analysis/v5/export_figure1.py --out reproduced_v5
```

This route redraws the four-panel figure from the existing rounded aggregate
summaries at a full-page width of 170 mm and 600 dpi. The submission Figure 1
was exported from the original producer using the original full-precision
aggregates; no analysis or model was rerun. The two routes have the same
clinical series, dose groups and cuff intervals, within source display rounding.

All transfer and packaging still require the repository data gate. The review
repository remains private; this update does not create a release or DOI.
