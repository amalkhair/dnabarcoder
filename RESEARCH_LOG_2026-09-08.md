# Research Log — 2026-09-08

## Goal for today
Create a clean, interpretable set of fungal biomarker results with reproducible commands.

## Outputs targeted
- Master biomarker comparison table
- Taxonomic coverage summary
- Coverage figure
- Coverage vs resolving power figure
- TEF1 variation result
- Figure interpretation notes
- Updated GitHub branch


## Step 2 — Inspected summary tables

### Characterization summary
File:
`results/final_summary/marker_characterization_summary.csv`

Contains:
- sequence counts
- length statistics
- indistinguishable-group statistics

### Threshold summary
File:
`results/final_summary/marker_threshold_summary.csv`

Contains:
- number of predicted cutoffs
- mean / median / min / max cutoff
- mean / median confidence
- total sequences in thresholded genera

Important correction:
The current threshold summary contains **confidence**, not F-measure. We must verify the source of the resolving-power figure before interpreting it as F-measure.

