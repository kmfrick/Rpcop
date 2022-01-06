// contiene la dimension del espacio original, la matriz y el xo  necesaria para
// transformar los puntos al sistema de coordenadas original

class M_a{
	private:
		int 		 Dim;
		int      profundidad;
		double	 **Ma;
		double    *xa;

		// vect ops

		double    *Mxv(double **M1,double *v);
		double    **MxM(double **M1,double **M2);
		double    *sum_v(double *v1,double *v2);

	public:
		// constructores
		M_a(int d,int p,double **M,double *x);
		~M_a();

		double *aplicar_Ma_punt(double *punt);
		double *aplicar_Ma_vect(double *vect);
		M_a   *donar_M_a(double **Mbopt,double *xo); // proporciona el M_a  per posar els punts en les cordenades de l'espai inicial, al nou subspai
};
