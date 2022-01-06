#include "ma.h"


class M_b{
	private :

		int      Dim;
		double    *xo;
		double    **Mb;
		double    **MId;
		double    **MInv;

		// vect ops

		double **inv(double **M);
		double *Mxv(double **M,double *v);
		double **MxM(double **M1,double **M2);
		double *mult_esc(double e,double *v);
		double  mult_v(double *v1,double *v2);
		double *sum_v(double *v1,double *v2);
		double *dif_v(double *v2,double *v1);
		double *norma_v(double *v);

	public:

		// constructoras

		M_b(int Dim,double *b);   // reusamos el vector b pasado por parametro.
		M_b(int Dim,double **n_M,double *n_xo);
		~M_b();

		M_b  *girar(int eix,double angle); // s'aplicara un gir a Mb per l'eix o dimensio donat. L'angle d'aquests girs sera el que diferenci un objecte matriu resultant d'un altre
		M_b  *replicar();                 // fará una copia de la M_b; 
		void  calcular_la_inversa(); 	// calculamos la inversa  realizados los giros de Mb y antes de aplicarla sobre los puntos. Si el Mb resulta optimo se calculara mas de 1 vez, sino una sola vez
		M_a   *donar_M_a(M_a *Ma);     // Ma es la matriu que els espais inferiors  necesitaran(una per cada subspai)
		// per pasar els punts a les coordenades originals
		void  rebre_xo(double *punt);
		double *aplicar(double *punt);  // aplica Mb al punt y dona el punt per l'espai inferior
		double *desaplicar(double *punt); // extraer para las coordenadas originales. op inversa a aplicar
		double *donar_bopt();

};
