# R figure redraws

These scripts redraw existing disclosure-controlled aggregate estimates. They do not refit models or recompute confidence intervals.

Requirements: R 4.6, grid, svglite, ragg, jsonlite and systemfonts; Arial. Figure 1 requires the project aggregate sources and Python for the existing disclosure gate. Figure 2 takes the existing `01_dose_contrasts.csv`; only the six specified pair-weighted, main-target highest-minus-lowest contrasts are drawn.

Run the project small-cell gate on the actual aggregate input before plotting and on the generated output afterward. Figure 1 also invokes this gate internally. Use the guarded project packager before sharing output. No patient-level inputs are required by these scripts.

```
Rscript draw_figure1.R PROJECT_ROOT OUTPUT_ROOT PYTHON_FOR_GATE
Rscript draw_figure2.R AGGREGATE_CONTRAST_CSV FIGURES_OUTPUT_DIRECTORY
```

Both figures use manually placed R grid primitives, Arial body text at 7.5 pt, 0.8 pt axes and 3.5 pt outward ticks. Outputs are editable SVG, embedded-font PDF, and 600 dpi TIFF/PNG. Cairo PDF page dimensions are rounded down to whole PDF points (less than 0.4 mm difference from the nominal 170 mm width).
