# Principal Curves of Oriented Points

This repository provides an implementation of principal curves on oriented points that can be compiled on modern systems, using Rcpp.

The original implementation, of which this repository is an adaptation, is provided [here](https://www-eio.upc.es/~delicado/PCOP/index.html) by the authors.

Principal curves on oriented points are introduced in Delicado and Huerta ([2003](https://link.springer.com/article/10.1007/s001800300145)).

## Tests

Use these commands from repository root:

```bash
# Run package tests with your default compiler flags while ignoring ~/.R/Makevars
R_MAKEVARS_USER=/dev/null R_LIBS=/tmp/Rlibs-san \
  /opt/R-devel-san/bin/Rscript -e "testthat::test_local('pkg')"
```

```bash
# Run full sanitizer check with ASan/UBSan instrumentation from tools/sanitize.Makevars
ASAN_OPTIONS='halt_on_error=1' \
UBSAN_OPTIONS='print_stacktrace=1:halt_on_error=1' \
R_MAKEVARS_USER=$PWD/tools/sanitize.Makevars \
R_LIBS=/tmp/Rlibs-san \
  /opt/R-devel-san/bin/R CMD check pkg
```
