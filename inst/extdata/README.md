# Vignette Fixtures

This directory contains compact fixtures used by the package vignettes.

The `vignette_fit_*.rds` and `vignette_sensitivity_fits.rds` files are
xz-compressed `dpmirt_fit` objects containing evenly thinned retained
posterior draws. They preserve the fields needed for vignette summaries,
plots, diagnostics, draw extraction, and WAIC comparison, but omit
session-bound compiled NIMBLE objects and raw sample storage.

The compact fit fixtures use 250 retained posterior draws selected from the
original 8,000 post-burn-in draws. They are intended for documentation and
examples, not for independent convergence assessment.
