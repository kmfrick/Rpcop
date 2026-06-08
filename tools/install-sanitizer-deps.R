lib <- "/tmp/Rlibs"
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

required <- c("Rcpp", "testthat", "princurve")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  install.packages(
    missing,
    repos = "https://cloud.r-project.org",
    Ncpus = 2
  )
}

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing R packages after installation: ", paste(missing, collapse = ", "), call. = FALSE)
}
