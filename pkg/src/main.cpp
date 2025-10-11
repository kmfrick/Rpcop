extern "C" {
#include <stdlib.h>
}
#include "espai.h"
#ifdef __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <Rcpp.h>
#ifdef __clang__
# pragma clang diagnostic pop
#endif

//' @title pcop_backend
//' @name pcop_backend
//' @description Internal backend used by \code{pcop()} to compute the principal curve defined in Delicado (2001) \doi{10.1007/s001800300145}.
//' @param x Numeric matrix of input points; see \code{pcop()}.
//' @param c_d Distance scaling parameter passed from \code{pcop()}.
//' @param c_h Bandwidth scaling parameter passed from \code{pcop()}.
//' @return Numeric matrix consumed by \code{pcop()}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix pcop_backend(const Rcpp::NumericMatrix &x, float c_d,
                                 float c_h) {
  int profreq = 1;
  int nparts = 4;
  espai *psp;
  float **Ma;
  int i;
  ll_p *ll_pt;
  float *mx;

  // inicializacin valores por defecto, despues en el fichero de setup puede
  // que cambien modificado 16/4/2002
  // PROFREQ =1;
  // NPARTs = 4;
  // C_h = 0.75;
  // C_d = 0.4;		// siempre menor que 0.5

  int Dim = x.ncol();
  ll_pt = new ll_p(Dim);
  for (i = 0; i < x.nrow(); i++) {
    float *d_punt = new float[Dim + 1];
    d_punt[0] = 1;
    d_punt++; // coordenada -1 pel pes.
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
  Ma = new float *[Dim];
  if (Ma == NULL) {
    Rcpp::stop("Could not allocate Ma.\n");
  }

  for (i = 0; i < Dim; i++) {
    Ma[i] = new float[Dim]();
    Ma[i][i] = 1;
  }
  mx = new float[Dim]();
  if (mx == NULL) {
    Rcpp::stop("Could not allocate mx.\n");
  }
  psp = new espai(ll_pt, Dim, 0);
  // modificacion 16/4/2002
  psp->inicializar_nparts_ch_cd(profreq, nparts, c_h, c_d);
  // fin
  psp->rebre_M_a(new M_a(Dim, 0, Ma, mx));
  static_cast<void>(psp->obtenir_VTG(&mx)); // static_cast to ignore return value
  int nrow, ncol;
  float *out = new float[(2 * Dim + 5) * x.nrow()];
  if (out == NULL) {
    Rcpp::stop("Could not allocate out.\n");
  }
  psp->obtenir_data(out, &nrow, &ncol);

  delete psp;
  Rcpp::NumericMatrix proy(nrow, ncol, out);
  proy = transpose(proy);
  return proy;
}
