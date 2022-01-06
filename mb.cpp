#include "mb.h"
#include <cstdlib>
#include <cmath>
extern "C" {
#include <memory.h>
}


M_b::M_b(int d,double *b){
	int i;
	int j;
	double *v_dif,*v_acum1,*v_acum2,*mesc;

	Dim =d;
	/* creamos e inicializamos Mb y MId */
	Mb    = (double**)malloc(Dim*sizeof(double *));
	MId   = (double**)malloc(Dim*sizeof(double *));
	MInv  = NULL;

	for (i=0;i<Dim;i++){
		Mb[i] = (double*)calloc(Dim,sizeof(double));
		MId[i] = (double*)calloc(Dim,sizeof(double));
	}
	for (i=0;i<Dim;i++){
		Mb[i][i] = 1;
		MId[i][i] = 1;
	}

	/* insertamos b_act en matriz identidad */

	i = 0;
	while(!b[i]) i++;
	if (i){
		for(j=i-1;j;j--)Mb[j+1]=Mb[j]; 
		for(j=Dim-2;j>i+1;j--)Mb[j+1]=Mb[j];
	}  

	Mb[0] = b;   // reusamos el vector b pasado por parametro.


	/* proceso de ortonormalizaci�n de Gram-Schmidt */

	v_acum1 = (double *)calloc(Dim,sizeof(double));
	Mb[0] = norma_v(Mb[0]);
	for (i=1;i<Dim;i++){
		delete v_acum1;
		v_acum1 = (double *)calloc(Dim,sizeof(double)); 
		for (j=0;j<i;j++){
			mesc = mult_esc(mult_v(Mb[i],Mb[j]),Mb[j]);
			v_acum2 = sum_v(v_acum1,mesc);
			delete v_acum1;
			v_acum1 = v_acum2;
		}
		v_dif = dif_v(Mb[i],v_acum1);
		Mb[i]= norma_v(v_dif);
	}

}

M_b::M_b(int d,double **n_M,double *n_xo){
	int i;
	Dim = d;
	Mb = n_M;
	xo = n_xo;

	/* creamos e inicializamos  MId */
	MId   = (double**)malloc(Dim*sizeof(double *));
	MInv  = NULL;

	for (i=0;i<Dim;i++){
		MId[i] = (double*)calloc(1,Dim*sizeof(double));
	}
	for (i=0;i<Dim;i++){
		MId[i][i] = 1;
	}

}

M_b::~M_b(){
	int i;

	if (MInv)
		for (i=0;i<Dim;i++){
			free(Mb[i]);
			free(MInv[i]);
			free(MId[i]);
		}
	else
		for (i=0;i<Dim;i++){
			free(Mb[i]);
			free(MId[i]);
		}

	free(Mb);
	free(MInv);
	free(MId);

	// free(xo) xo  es gestiona desde espai (calcular_corba_en_sentit(), calcular_corba_en_sentit_contrari()).
}


M_b *M_b::girar(int eix,double angle){
	// 0 < eix < Dim                    giro bo i  eix en sentit positiu

	double **Mbaux;

	// insertamos el giro en la matriz identidad
	MId[0][0] = cos(angle);
	MId[0][eix] = sin(angle);
	MId[eix][eix] = cos(angle);
	MId[eix][0] = -1*sin(angle);

	Mbaux = MxM(Mb,MId);

	// corregimos la matriz identidad
	MId[0][0]  = 1;
	MId[0][eix]= 0;
	MId[eix][eix] = 1;
	MId[eix][0]= 0;

	return new M_b(Dim,Mbaux,xo);
}

M_b *M_b::replicar(){
	double **c_Mb;
	int i,j;

	c_Mb = (double **) malloc(Dim*sizeof(double *));
	for (i=0;i<Dim;i++)
		c_Mb[i] = (double *)malloc(Dim*sizeof(double));


	for (i=0;i<Dim;i++) memmove(c_Mb[i],Mb[i],Dim*sizeof(double));
	//for (j=0;j<Dim;j++) c_Mb[i][j] = Mb[i][j];  


	return new M_b(Dim,c_Mb,xo);
}


void M_b::calcular_la_inversa(){
	int i,j;
	double **c_Mb;
	// calculamos la inversa  realizados los giros de Mb y antes de aplicarla sobre puntos. Si el Mb resulta optimo se calculara mas de 1 vez

	if (MInv){
		for (i=0;i<Dim;i++)  free(MInv[i]);
		free(MInv);
	}

	/* copiar la Mb  */ // ya que sera modificada por inv()

	c_Mb = (double **) malloc(Dim*sizeof(double *));
	for (i=0;i<Dim;i++)
		c_Mb[i] = (double *)malloc(Dim*sizeof(double));


	for (i=0;i<Dim;i++)  memmove(c_Mb[i],Mb[i],Dim*sizeof(double));
	//for (j=0;j<Dim;j++) c_Mb[i][j] = Mb[i][j];  

	/* calcular inversa */

	MInv = inv(c_Mb);

	/* borrar c_Mb  */
	for (i=0;i<Dim;i++)  free(c_Mb[i]);
	free(c_Mb);
}


double  *M_b::aplicar(double *punt){  /* aplica Mb al punt y dona el punt per l'espai inferior */
	double *p2;
	double *p3;
	p2 = dif_v(punt,xo);
	p3 = Mxv(MInv,p2);     /* se pasa la inversa */

	// la primera coordenada sera la distancia al pla, les seguents, les del punt
	// al pla inferior.

	free(p2);
	return p3;
}

double  *M_b::desaplicar(double *punt){  /* operaci�n inversa a aplicar */
	double *p2;
	double *p3;

	p2 = Mxv(Mb,punt);     
	p3 = sum_v(p2,xo);

	// la primera coordenada sera la distancia al pla, les seguents, les del punt
	// al pla inferior.

	free(p2);
	return p3;
}


void M_b::rebre_xo(double *punt){
	// se ejecutara una vez por avance_cluster() y todas las veces que el pop candidato cruce el hiperplano del pop anterior.
	//free xo;    //  no fa falta alliverar l'espai. O b� ha estat eliminat al crear el nou xo, o es un xo compartir amb la matriu del pop anterior.
	xo = punt;  

}

/* Ma es la matriu que l'espai inferior  necesacitar� per pasar*/
/* els seus punts i vectors finals a les coordenades originals    */

M_a *M_b::donar_M_a(M_a *Ma){
	// se ejecutara como m�ximo una vez si resulta ser el M_b optimo.

	return Ma->donar_M_a(Mb,xo);  // xo, en aquest punt haur� de ser de tamany Dim+prof.

}


double *M_b::donar_bopt(){
	return Mb[0];
}


// PRIVATE //////////////////////////////////////////////////////////////////

double **M_b::inv(double **M){
	int i,j,ii,aux;
	double **Inv,Mji;

	Inv = (double**)malloc(Dim*sizeof(double *));
	for (i=0;i<Dim;i++)
		Inv[i] = (double *)calloc(Dim,sizeof(double));
	for (i=0;i<Dim;i++)
		Inv[i][i] = 1;


	for (i=0;i<Dim;i++){
		for(j=(i+1)%Dim;j!=i;j=(j+1)%Dim){
			Mji = M[j][i];
			for(ii=0;ii<Dim;ii++){
				Inv[j][ii]= Inv[j][ii]*M[i][i]-Inv[i][ii]*Mji;
				M[j][ii]= M[j][ii]*M[i][i]-M[i][ii]*Mji;
			}
		}
	}


	/*  for (i=0;i<Dim;i++){
	    j=(i+1)%Dim;
	    aux = j;            // si no se vuelca la  j, el % no funciona correctamente borlan c++ !!��???
	    while (j!=i){
	    Mji = M[j][i];
	    for(ii=0;ii<Dim;ii++){
	    Inv[j][ii]= (Inv[j][ii]*M[i][i])-(Inv[i][ii]*Mji);   // no pot modificarse M avans de Inv
	    M[j][ii]= (M[j][ii]*M[i][i])-(M[i][ii]*Mji);
	    }
	    j=(j+1)%Dim;
	    }
	    }
	 */
	for (j=0;j<Dim;j++){
		for(ii=0;ii<Dim;ii++)
			Inv[j][ii]= Inv[j][ii]/M[j][j];
	}
	return Inv;
}

double *M_b::Mxv(double **M1,double *v){
	// vxM, trabajamos con vectores fila.
	int i,j;
	double sum;
	double *v3 = (double *) malloc(Dim*sizeof(double));

	for(i=0;i<Dim;i++){
		sum = 0;
		for(j=0;j<Dim;j++){
			sum += v[j]*M1[j][i];
		}
		v3[i] = sum;
	}
	return v3;
}


double **M_b::MxM(double **M1,double **M2){
	// vxM, trabajos con vectores fila.
	int i,ii,j;
	double sum;
	double **M3;

	M3 = (double**)malloc(Dim*sizeof(double *));

	for (i=0;i<Dim;i++)
		M3[i] = (double *)calloc(1,Dim*sizeof(double));

	for(i=0;i<Dim;i++){
		for(ii=0;ii<Dim;ii++){
			sum = 0;
			for(j=0;j<Dim;j++){
				sum += M1[i][j]*M2[j][ii];
			}
			M3[i][ii] = sum;
		}
	}
	return M3;
}

double *M_b::mult_esc(double e,double *v){
	int i;
	double *v3;

	v3 = (double *)malloc(Dim* sizeof(double));
	for(i=0;i<Dim;i++)  v3[i] = v[i]*e;
	return v3;
}

double M_b::mult_v(double *v1,double *v2){
	int i;
	double sum = 0.0;
	for(i=0;i<Dim;i++){
		sum += v1[i]* v2[i];
	}
	return sum;
}

double *M_b::sum_v (double *v1,double *v2){
	int i;
	double *v3;

	v3 = (double *)malloc(Dim* sizeof(double));
	for(i=0;i<Dim;i++)  v3[i] = v1[i]+v2[i];
	return v3;
}

double *M_b::dif_v (double *v2,double *v1){
	int i;
	double *v3;

	v3 = (double *)malloc(Dim* sizeof(double));
	for(i=0;i<Dim;i++)  v3[i] = v2[i]-v1[i];
	return v3;
}

double *M_b::norma_v(double *v){
	// devuelve el vector normalizado
	int i;
	double nrm =0.0;
	double *v3;

	v3 = (double *)malloc(Dim* sizeof(double));
	for(i=0;i<Dim;i++)  nrm += pow(v[i],2);
	nrm = sqrt(nrm);
	for(i=0;i<Dim;i++)  v3[i] = v[i]/nrm;
	return v3;
}





