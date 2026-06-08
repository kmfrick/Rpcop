# Principal Curves of Oriented Points

This repository provides an implementation of principal curves on oriented points that can be compiled on modern systems, using Rcpp.

The original implementation, of which this repository is an adaptation, is provided [here](https://www-eio.upc.es/~delicado/PCOP/index.html) by the authors.

Principal curves on oriented points are introduced in Delicado and Huerta ([2003](https://link.springer.com/article/10.1007/s001800300145)).

## Tests

Use these commands from repository root:

```bash
# Build the source tarball, then run a CRAN-style check outside the repo tree
pkgdir=$PWD
mkdir -p /private/tmp/rpcop-check
cd /private/tmp/rpcop-check
R CMD build --no-manual --no-build-vignettes "$pkgdir"
R_MAKEVARS_USER=/dev/null R CMD check --as-cran --no-manual --output=/private/tmp/rpcop-check Rpcop_1.2.1.tar.gz
```

```bash
# Run Linux ASAN/UBSAN diagnostics
vagrant up --provision
```
