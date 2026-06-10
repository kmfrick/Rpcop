extern "C" {
#include <stdlib.h>
}
#include "espai.h"
#include <Rcpp.h>
#include <memory>

//' @title pcop_backend
//' @name pcop_backend
//' @description Internal backend used by \code{pcop()} to compute the principal curve defined in Delicado and Huerta (2003) \doi{10.1007/s001800300145}.
//' @param x Numeric matrix of input points; see \code{pcop()}.
//' @param c_d Distance/scaling parameter passed from \code{pcop()}.
//' @param c_h Bandwidth scaling parameter passed from \code{pcop()}.
//' @return Numeric matrix consumed by \code{pcop()}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix pcop_backend(const Rcpp::NumericMatrix &x, float c_d,
                                 float c_h) {
  int profreq = 1;
  int nparts = 4;
  int i;

  // Default initialization values (may be overridden in setup).
  // Modified 2002-04-16.
  // PROFREQ = 1;
  // NPARTs = 4;
  // C_h = 0.75;
  // C_d = 0.4; // always < 0.5

  int Dim = x.ncol();
  std::unique_ptr<ll_p> ll_pt(new ll_p(Dim));
  for (i = 0; i < x.nrow(); i++) {
    float *d_punt = new float[Dim + 1];
    d_punt[0] = 1;
    d_punt++; // coordinate -1 holds the weight.
    for (int j = 0; j < x.ncol(); j++) {
      d_punt[j] = x(i, j);
    }
    ll_pt->add_ordX_principal(d_punt); // ###
  }
  if (ll_pt->n_punts() < NPTMIN * Dim) {
    Rcpp::stop("Warning: Not enough points in data matrix. At least %d points "
               "are needed for dimension %d.\n",
               NPTMIN * Dim, Dim);
  }
  float **Ma = new float *[Dim];
  if (Ma == NULL) {
    Rcpp::stop("Could not allocate Ma.\n");
  }
  int ma_rows = 0;
  auto ma_deleter = [&](float **matrix) {
    if (!matrix) {
      return;
    }
    for (int row = 0; row < ma_rows; row++) {
      delete[] matrix[row];
    }
    delete[] matrix;
  };
  std::unique_ptr<float *, decltype(ma_deleter)> ma_guard(Ma, ma_deleter);

  for (i = 0; i < Dim; i++) {
    Ma[i] = new float[Dim]();
    Ma[i][i] = 1;
    ma_rows++;
  }
  std::unique_ptr<float[]> mx(new float[Dim]());
  if (mx == NULL) {
    Rcpp::stop("Could not allocate mx.\n");
  }
  std::unique_ptr<espai> psp(new espai(ll_pt.release(), Dim, 0));
  // Modified 2002-04-16.
  psp->inicializar_nparts_ch_cd(profreq, nparts, c_h, c_d);
  psp->rebre_M_a(new M_a(Dim, 0, ma_guard.release(), mx.release()));
  float *returned_xomig = NULL;
  static_cast<void>(
      psp->obtenir_VTG(&returned_xomig)); // static_cast to ignore return value
  int nrow, ncol;
  int max_output_rows = (2 * x.nrow()) + 1;
  int output_cols = 2 * Dim + 5;
  std::unique_ptr<float[]> out(new float[output_cols * max_output_rows]);
  if (out == NULL) {
    Rcpp::stop("Could not allocate out.\n");
  }
  psp->obtenir_data(out.get(), &nrow, &ncol, max_output_rows);

  Rcpp::NumericMatrix proy(nrow, ncol, out.get());
  proy = transpose(proy);
  psp.reset();
  std::unique_ptr<float[]> returned_xomig_guard(returned_xomig);
  return proy;
}
