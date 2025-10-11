#' Principal Curve of Oriented Points
#'
#' @description Computes a principal curve as defined in Delicado (2001)
#' \doi{10.1007/s001800300145}.
#' @param x A numeric matrix of \eqn{n} points in dimension \eqn{p}.
#' @param Ch The smoothing parameter \eqn{h} is \eqn{C_H} times the value given by
#'   the normal reference rule. Default value \eqn{1.5}. Constraints
#'   \eqn{0.5 \le C_H \le 1.5}.
#' @param Cd The distance between two consecutive principal oriented points in a
#'   PCOP is about \eqn{C_D} times the value of the smoothing parameter
#'   \eqn{h}. Default value \eqn{0.3}. Constraints \eqn{0.25 \le C_D \le 0.5}.
#' @param plot.true If \code{TRUE}, produce a plot of the resulting curve.
#' @param ... Additional parameters passed to \code{lines}.
#' @return
#' A list with two elements:
#' \describe{
#'   \item{pcop.f1}{Data frame storing the principal curve of oriented points in
#'     the original format, with columns \code{param}, \code{dens},
#'     \code{span}, \code{orth.var}, \code{pop1}, \code{pop2}, \ldots,
#'     \code{pr.dir1}, \code{pr.dir2}, \ldots}
#'   \item{pcop.f2}{List conforming to the format used in \pkg{princurve}; see that
#'     package for details.}
#' }
#' @examples
#' n <- 500
#' p <- 3
#' x <- matrix(rnorm(n * p), ncol = p) %*% diag(p:1)
#' pcop(x, plot.true = TRUE, lwd = 4, col = 2)
#'
#' x <- runif(100, -1, 1)
#' x <- cbind(x, x ^ 2 + rnorm(100, sd = 0.1))
#' pcop(x, plot.true = TRUE, lwd = 4, col = 2)
#' @importFrom princurve project_to_curve
#' @importFrom graphics lines
#' @export
pcop <- function(x, Ch = 1.5, Cd = 0.3, plot.true = FALSE, ...) {
  proyec <- as.data.frame(pcop_backend(x, Cd, Ch))
  pr.curve <- proyec[, -1]

  p <- ncol(x)
  names(pr.curve) <- c(
    "param", "dens", "span", "orth.var",
    paste0("pop", seq_len(p)),
    paste0("pr.dir", seq_len(p))
  )

  s <- as.matrix(pr.curve[, 5:(4 + p)])
  pc <- project_to_curve(x, s)

  if (plot.true) {
    adjust_range <- function(values, cte = 1) {
      rgx <- range(values)
      extra <- diff(rgx) * (cte - 1) / 2
      c(rgx[1] - extra, rgx[2] + extra)
    }
    plot(
      x[, 1:2],
      xlim = adjust_range(x[, 1], 1.1),
      ylim = adjust_range(x[, 2], 1.1),
      col = 8
    )
    lines(pc$s[pc$ord, 1:2], ...)
  }

  list(pcop.f1 = pr.curve, pcop.f2 = pc)
}
