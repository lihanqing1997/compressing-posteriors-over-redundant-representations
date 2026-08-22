# Compressing Posteriors over Redundant Representations

This repository contains the code and exact-reference audit artifacts for
certified regional compression of finite discrete posteriors.

The method selects a minimum-declared-cost union of disjoint regions under an
omitted-mass tolerance. It reports a feasible retained union, an omitted-mass
bound, lower and upper bounds on the constrained optimum, and a certified cost
gap. The included study evaluates the method in sparse Gaussian Bayesian model
averaging with aligned, perturbed, and crossed region trees.

## Contents

`supplement/` is the verified reproducibility package accompanying the
manuscript. It contains:

- the active R implementation and exact Gaussian BMA helpers;
- the frozen exact-study protocol and deterministic seed rules;
- six canonical result tables with 4,800 registered method-setting runs;
- derived summaries and the seven-check audit;
- release, protocol, and unit-test scripts.

Historical experiments, pilot outputs, caches, logs, manuscript sources, and
superseded density-ratio reporting code are intentionally excluded.

## Requirements

- R 4.4.2 or a compatible recent R release
- the CRAN package `yaml`
- `certutil` on Windows, or `sha256sum`/`shasum` on Unix-like systems

## Verify the released artifact

Run the following commands from `supplement/`:

```text
Rscript sim/verify_release.R
Rscript sim/validate_tmlr_revision_protocol.R
Rscript sim/run_region_restricted_smc_unit_tests.R
Rscript sim/run_constrained_region_search_unit_tests.R
Rscript sim/run_factorized_region_mass_unit_tests.R
```

The release verifier checks the 23-file allowlist, SHA-256 hashes, 4,800
registered rows, protocol grid, seed formulas, budget accounting, summaries,
and reported audit results.

## Recompute derived outputs

From `supplement/`, write recomputed files outside the verified artifact:

```text
Rscript sim/audit_constrained_exact_production.R ../tmlr_recomputed/constrained_exact_gate1_audit.csv
Rscript sim/summarize_tmlr_revision_results.R ../tmlr_recomputed/tmlr_constrained_exact_summary.csv
Rscript sim/make_tmlr_revision_figures.R ../tmlr_recomputed/figures
```

Full production is computationally expensive and is not required for artifact
verification. See `supplement/REPRODUCTION.txt` for the registered single-run
example and additional details.

## Scope

The package supports the reported exact optimization, cost-sandwich,
omitted-mass, stochastic interval-coverage, and grouping-sensitivity results.
It does not claim automatic region discovery, polynomial worst-case runtime,
universal speedups, or superiority to full-posterior samplers.

## License

The software is released under the MIT License. See `LICENSE`.
