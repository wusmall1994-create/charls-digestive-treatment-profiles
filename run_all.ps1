$ErrorActionPreference = "Stop"

python analysis/01_prepare_descriptive_outputs.py
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Rscript analysis/02_run_primary_and_sensitivity_models.R
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Rscript analysis/03_make_figures.R
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "Analysis completed. See source_data/ and figures/."
