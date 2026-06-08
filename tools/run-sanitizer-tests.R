.libPaths(c("/tmp/Rlibs", .libPaths()))

library(Rpcop)

source("tests/testthat/helper-stress-data.R")

sanitizer_n_for_dim <- function(p) {
  50L * as.integer(p)
}

sanitizer_backend_cases <- function() {
  list(
    list(
      id = "gaussian_p2",
      p = 2L,
      x = stress_gen_gaussian(sanitizer_n_for_dim(2L), 2L, seed = 1002L),
      Cd = 0.30,
      Ch = 1.00
    ),
    list(
      id = "anisotropic_p3",
      p = 3L,
      x = stress_gen_anisotropic(sanitizer_n_for_dim(3L), 3L, seed = 2003L),
      Cd = 0.25,
      Ch = 1.50
    ),
    list(
      id = "rank2_p5",
      p = 5L,
      x = stress_gen_rank2(sanitizer_n_for_dim(5L), 5L, seed = 3005L),
      Cd = 0.50,
      Ch = 0.50
    )
  )
}

sanitizer_wrapper_cases <- function() {
  list(
    list(
      id = "wrapper_parabola_p2",
      p = 2L,
      x = local({
        set.seed(5001L)
        n <- sanitizer_n_for_dim(2L)
        t <- runif(n, -2, 2)
        cbind(t, t^2 + rnorm(n, sd = 0.07))
      }),
      Cd = 0.30,
      Ch = 1.50
    ),
    list(
      id = "wrapper_scurve_p3",
      p = 3L,
      x = local({
        set.seed(5002L)
        n <- sanitizer_n_for_dim(3L)
        t <- runif(n, -pi, pi)
        cbind(sin(t), t / pi, sign(t) * (1 - cos(t)) + rnorm(n, sd = 0.05))
      }),
      Cd = 0.25,
      Ch = 1.00
    ),
    list(
      id = "wrapper_rank2_p5",
      p = 5L,
      x = stress_gen_rank2(sanitizer_n_for_dim(5L), 5L, seed = 5003L),
      Cd = 0.50,
      Ch = 0.50
    )
  )
}

check_matrix <- function(x, label) {
  if (!is.matrix(x)) {
    stop(label, " did not return a matrix", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(label, " returned an empty matrix", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop(label, " returned non-finite values", call. = FALSE)
  }
}

run_backend_case <- function(case) {
  label <- paste0("backend/", case$id)
  message("START ", label)
  out <- Rpcop:::pcop_backend(case$x, case$Cd, case$Ch)
  check_matrix(out, label)
  message("OK ", label)
}

run_wrapper_case <- function(case) {
  label <- paste0("wrapper/", case$id)
  message("START ", label)
  fit <- pcop(case$x, Ch = case$Ch, Cd = case$Cd, plot.true = FALSE)
  if (!is.list(fit) || !all(c("pcop.f1", "pcop.f2") %in% names(fit))) {
    stop(label, " returned a malformed fit", call. = FALSE)
  }
  if (!is.data.frame(fit$pcop.f1) || nrow(fit$pcop.f1) < 1L) {
    stop(label, " returned malformed pcop.f1", call. = FALSE)
  }
  if (!identical(colnames(fit$pcop.f1), stress_expected_cols(case$p))) {
    stop(label, " returned unexpected pcop.f1 columns", call. = FALSE)
  }
  if (!all(is.finite(as.matrix(fit$pcop.f1)))) {
    stop(label, " returned non-finite pcop.f1 values", call. = FALSE)
  }
  check_matrix(fit$pcop.f2$s, paste0(label, "/pcop.f2$s"))
  message("OK ", label)
}

for (case in sanitizer_backend_cases()) {
  run_backend_case(case)
}
for (case in sanitizer_wrapper_cases()) {
  run_wrapper_case(case)
}

message("DONE sanitizer stress cases")
