test_that("hard backend stress cases remain finite and non-empty", {
  cases <- stress_backend_cases()
  expect_identical(length(cases), 4L)

  for (case in cases) {
    msg <- sprintf("backend case: %s", case$id)
    out <- tryCatch(
      Rpcop:::pcop_backend(case$x, case$Cd, case$Ch),
      error = identity
    )
    expect_false(inherits(out, "error"), info = msg)
    if (inherits(out, "error")) {
      next
    }
    expect_true(is.matrix(out), info = msg)
    expect_true(nrow(out) > 0L, info = msg)
    expect_true(ncol(out) > 0L, info = msg)
    expect_true(all(is.finite(out)), info = msg)
  }
})

test_that("hard wrapper stress cases return well-formed outputs", {
  cases <- stress_wrapper_cases()
  expect_identical(length(cases), 4L)

  for (case in cases) {
    msg <- sprintf("wrapper case: %s", case$id)
    fit <- tryCatch(
      pcop(case$x, Ch = case$Ch, Cd = case$Cd, plot.true = FALSE),
      error = identity
    )
    expect_false(inherits(fit, "error"), info = msg)
    if (inherits(fit, "error")) {
      next
    }

    expect_true(is.list(fit), info = msg)
    expect_true(all(c("pcop.f1", "pcop.f2") %in% names(fit)), info = msg)

    expect_true(is.data.frame(fit$pcop.f1), info = msg)
    expect_true(nrow(fit$pcop.f1) > 0L, info = msg)
    expect_identical(colnames(fit$pcop.f1), stress_expected_cols(case$p), info = msg)
    expect_true(all(is.finite(as.matrix(fit$pcop.f1))), info = msg)

    expect_true(is.list(fit$pcop.f2), info = msg)
    expect_true(is.matrix(fit$pcop.f2$s), info = msg)
    expect_true(nrow(fit$pcop.f2$s) > 0L, info = msg)
    expect_true(all(is.finite(fit$pcop.f2$s)), info = msg)
  }
})
