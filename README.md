# CHARLS questionnaire treatment-state transitions

This repository contains the final analysis code, aggregate source data, and figures for:

> *Transitions Between Reported Questionnaire Treatment States Among Middle-Aged and Older Chinese Adults With Self-Reported Digestive Disease: A National Panel Study, 2011–2018*

The analysis describes adjacent-wave transitions between treatment-response profiles in CHARLS and estimates standardized associations by joint hukou–residence group. A `none-listed` response means that none of the questionnaire-listed treatment options was selected; it is not interpreted as untreated disease, nonadherence, or unmet need.

## Contents

- `analysis/01_prepare_descriptive_outputs.py`: constructs wave and person-interval samples and produces descriptive transition outputs.
- `analysis/02_run_primary_and_sensitivity_models.R`: produces the primary GEE estimates, standardized probabilities and risk differences, unadjusted estimates, and sensitivity analyses.
- `analysis/03_make_figures.R`: creates the two submission-ready TIFF figures.
- `source_data/`: aggregate machine-readable data supporting the manuscript tables and figures.
- `figures/`: final figures at 300 dpi.

No individual-level CHARLS data are included.

## Data access

CHARLS data are available to registered users from the [official CHARLS data portal](https://charls.pku.edu.cn/). Place the following files in `data/`, or set `CHARLS_DATA_DIR` to the directory containing them:

```text
harmonized.dta
w1_health.dta
w2_health.dta
w3_health.dta
w4_health.dta
```

The harmonized file is the RAND/Harmonized CHARLS longitudinal release used in the analysis. Raw data remain subject to the CHARLS data-use terms.

## Reproduce the analysis

The verified environment used Python 3.12.14 and R 4.6.1. Python package versions are pinned in `requirements.txt`; required R package versions are listed in `R-packages.txt`.

From PowerShell, optionally set the data directory and run:

```powershell
$env:CHARLS_DATA_DIR = "D:\path\to\CHARLS"
.\run_all.ps1
```

The scripts write aggregate outputs to `source_data/` and figures to `figures/`. The final complete-case primary models included 2,782 person-intervals (1,355 events) for none-listed to listed and 5,105 person-intervals (1,372 events) for listed to none-listed.

## Citation

Use the stable concept DOI [10.5281/zenodo.21966963](https://doi.org/10.5281/zenodo.21966963), which resolves to the latest archived release. Citation metadata are also provided in `CITATION.cff`.

## License

Code is released under the MIT License. The license does not apply to CHARLS source data, which are not redistributed here.
