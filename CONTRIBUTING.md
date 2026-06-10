# Contributing

## Development Workflow

- Install package dependencies with `rv sync`.
- Run `Rscript tests/testthat.R` for a local test pass.
- Run a CRAN-style package check from outside the repository tree:
  `R CMD build --no-manual .`, then `R CMD check --as-cran --no-manual`.
- Regenerate documentation with a temporary script that calls
  `roxygen2::roxygenise()` after updating roxygen comments.
- For rOpenSci readiness, run `srr::srr_stats_pre_submit(".")` and
  `pkgcheck::pkgcheck(".")` in an environment where those packages are
  installed.

## Pull Requests

- Keep behavior changes covered by `testthat`.
- Prefer targeted patches over large formatting-only rewrites.
