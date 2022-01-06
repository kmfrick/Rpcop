#include "ll_pnt.h"
#include "ll_flt.h"
// #include "Ma.h"
#include "mb.h"
#include "ll_p.h"

extern "C" {
#include <memory.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
}



//#define PROF_REQ 1  // 0= calculamos curva sobre STV, 1= sobre las curvas de los subespacios,2= sobre las superficies de los subspacios.....



#define C_EPS 0.05   // Sin recursividad, el valor 0.01 fue muy bien
#define NPTMIN 50    // n puntos mínimo para calcular la curva
#define LD   0.5

#define PINZA_MAX PI/4
#define PINZA_MIN PI/((4*NPARTS)+1)  // 90/4 

class espai {
	private :

		int   Dim;
		int   profundidad;
		//###	 M_a   *Ma;
		ll_p  *ll_pt;

		double suma_d;  // suma distancies del mstree.
		double h_tail;  // distancia al pla.
		double delta;   // advance
		double diam;    // longitud de la diagonal del cluster. 
		double *eps_x;  // eps_x conte les distancies minimes entre el xo del cluster i el xmean del cluster perque aquest sigui valid. Es unic per tots els cluster de l'espai actual.
		double *xomig;  // l'xomig es necesitara per obtenir_STV() i per obtenir el xo i bo_opt inicial de una corba.
		// quan calculem la corba en sentit contrari no partirem d'aquest xomig si no del corresponent pop.   

		int bficorba;  // boolean que controla si hi ha aparegut un nou punt en el calcul del nou pop.

		typedef struct m_d_s{
			double *xmean;
			double span;
			double density;
		} m_d_s;

		struct opt{
			double  VTG;
			M_b    *Mb;
			M_b    *Mb_ant;
			espai  *espai;
			m_d_s  mds;
		}optims;

		struct x{
			double *act;
			double *ant;
		}xo;

		struct pop{
			double *alpha;
			double I;
			double *b_ast;
			double var_k;
			double span;
			double density;
			espai *espai;         
		};

		//	 ll_pnt *ll_pop; ###

		double Var_PC;
		double Var_res;
		double GTV;

		int    dist_al_pla(double *n_punt);
		void   calcular_htail_delta_xomig_epsx();
		espai  *obtenir_cluster(M_b *Mb,m_d_s *mds);
		double  *treure_coord(double *n_pnt);  // fem la projeccio sobre el pla per pasar al subspai de dimensio inferior
		int    fi_corba(double *n_pnt);
		int    no_creua_corba(double *pop);
		void   calcular_Mb(int ejegir,M_b *Mb,double porcion_pinza);
		double  calcular_corba_en_un_sentit();
		double  calcular_corba_en_sentit_contrari();
		double  finalitzacio();       /* retorna la VTG de la corba */
		double  *allargar(double *bopt);
		double  obtenir_STV();
		double  *obtenir_bo_inicial(double *alfa);
		double  Bmst();
		double  gammln(double xx);
		double  kernel(double d);
		int    major(double *v1,double *v2);
		void   eigens(double *A,double *RR,double *E,int N );  //Copyright 1973, 1991 by Stephen L. Moshier

		// vect ops
		double  distancia(double *pnt1,double *pnt2);
		double *mult_esc(double e,double *v);
		double  mult_v(double *v1,double *v2);
		double *sum_v(double *v1,double *v2);
		double *dif_v (double *v2,double *v1);
		double *norma_v(double *v);

		// modificado 16/4/2002 declaramos las variables de forma statica para que no varien para los diferentes espais
		static int  PROF_REQ;
		static int  NPARTS;
		static double C_H;  
		static double C_D;     // siempre menor que 0.5
		// fin

	public:

		espai(ll_p *ll_punts,int d,int p);
		~espai();
		double   obtenir_VTG(double **xm);
		void    rebre_M_a(M_a *n_Ma);      // li pasem el nou Ma al subespai.


		M_a   *Ma;      //###
		ll_pnt *ll_pop; //###



		// modificado 16/4/2002 inicializamos las variables
		// operacions inicialització i extracció d'informació només per el 1er espai.
		void  inicializar_nparts_ch_cd(int profreq,int nparts,double c_h,double c_d);		
		void  obtenir_data(double*, int*, int*);
		// fin

};
