extern "C" {
#include <stdlib.h>
}
#include "espai.h"
#include <Rcpp.h>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
Rcpp::NumericMatrix pcop (Rcpp::NumericMatrix x, double c_d, double c_h, int profreq = 1, int nparts = 4){
	espai *psp;
	double **Ma;
	int i;
	ll_p *ll_pt;
	double vtg;
	double *mx;

	// inicialización valores por defecto, despues en el fichero de setup puede que cambien
	// modificado 16/4/2002
	//PROFREQ =1;
	//NPARTs = 4;
	//C_h = 0.75;
	//C_d = 0.4;		// siempre menor que 0.5

	int Dim = x.ncol();
	ll_pt = new ll_p(Dim);
	for (i = 0; i < x.nrow(); i++) {
		if (&(x.row(i)[0]) == NULL) {
			Rcpp::stop("0-th element of row is null\n");
		}
		ll_pt->add_ordX_principal(&(x.row(i)[0])); //###
	}

	Ma = (double **) malloc(Dim*sizeof(double *));
	if (Ma == NULL) {
		Rcpp::stop("Could not allocate Ma.\n");
	}

	for (i=0;i<Dim;i++){
		Ma[i] = (double *)calloc(Dim,sizeof(double));
		Ma[i][i]=1;
	}
	mx = (double *) calloc(Dim,sizeof(double));
	if (mx == NULL) {
		Rcpp::stop("Could not allocate mx.\n");
	}
	psp = new espai(ll_pt,Dim,0);
	//modificacion 16/4/2002
	psp->inicializar_nparts_ch_cd(profreq,nparts,c_h,c_d);
	//fin
	psp->rebre_M_a(new M_a(Dim,0,Ma,mx));
	vtg = psp->obtenir_VTG(&mx);
	int nrow, ncol;
	double *out = (double*)malloc((2 * Dim + 5) * x.nrow() * sizeof(float));
	if (out == NULL) {
		Rcpp::stop("Could not allocate out.\n");
	}
	psp->obtenir_data(out, &nrow, &ncol);

	delete psp;
	return Rcpp::NumericMatrix(nrow, ncol, out); // FIXME: Probable memory leak

}
