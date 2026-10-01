# SPRING 0.99.0

* Initial Bioconductor submission.
* Replaced `mixedCCA` with `latentcor` for latent correlation estimation.
* Removed the dependency on `pulsar` (archived on CRAN). StARS stability
  selection is now implemented internally and gives results identical to
  `pulsar::pulsar()` / `pulsar::refit()` for the same `seed`. The structure of
  the `SPRING()` output is unchanged.
* Removed the unused `SpiecEasi` import (which also depends on `pulsar`).
* `SPRING()` no longer resets the global RNG state when `seed` is given.
* `SPRING()` now stops early with an informative error if `lambdaseq` is not a
  decreasing sequence of positive numbers.
* Added `QMP` dataset with taxonomy-level column labels.
