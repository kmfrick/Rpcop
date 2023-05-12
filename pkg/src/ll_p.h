
#define FALSE 0
#define TRUE 1
#define ESQUERRA 0
#define DRETA 1
#define X 0
#define INF 9999

extern "C" {
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <string.h>
}
#include "ll_q.h"
#include "pila.h"

class ll_p {

private:
  int Dim;
  float dmax;
  int orcluster; // primero damos los candidatos a la derecha y luego a la
                 // izquierda de xo sobre la coordenada X
  int numcl;     // lleva la cuenta del n de rondas
  float sum_w;   // suma de los pesos de los puntos
  float suma_d;

  typedef struct node {
    float *coord;
    int marca;
    node *seg[2];
    void *noin[2]; // 1:enlaza nodos fuera del spamming tree. 2:enlaza
                   // satelites del nodo
  } node;

  typedef struct node_satelit {
    node *ptnode;
    node_satelit *seg;
  } node_satelit;

  node *xorig;   // punto medio de la curba
  node *xoant;   // punto del cluster cercano a un xo dado
  node *semilla; // node a partir del que obtenemos sus satelites para usarlos
                 // como candidats
  node_satelit *candidat; // posible punto del cluster enviado a validar

  pila p_n;

  int vn_punts;
  node *topright;
  node *topleft;
  float *min;
  float *max;
  float *x_mean; // per calcular el xmig

  // modificadoras

  void
  mstinsertar(node *pt); // marca com dins del msptree i treu de  la llista noin
  void add_satelit(int, node *, node *); // enlaza satelites

  // consultoras
  int mstinsertat(node *pt); // mira si esta insertat al min. spanning tree */

  // inicialitzacio
  // float *calcular_xomig_corba(); // cerca el punt origen de la corba (el ms
  // proper a xmean)
  void calcular_max_min_cluster();

  void obtener_quartiles(
      ll_q *ll_qt); // calculamos los quartiles sobre las distancias obtenidas
                    // del minium spaming tree de los puntos.
  float *obtener_satelites();

  // ops vect

  float *mult_esc(float e, float *v);
  float distancia(float *pnt1, float *pnt2);
  float *sum_v(float *v1, float *v2);

public:
  // constructuras
  ll_p(int d);
  ~ll_p();

  void add_ordX_principal(float *vect);

  // inicializacion
  float inicialitzacio_principal();
  void inicialitzacio_final();
  void tornar_a_xomig(); // tornem al punt origen per ferla en sentit contrari
  // als espais finals, tornarem el xmean ponderat per l'htail del cluster

  // modificadoras
  void
  trobar_primer_candidat_clt(float *xo); // cerca el primer candidat al cluster.
  float *canviar_orientacio_clt(); // comena la cerca en sentit contrari, torna
                                   // el 1er candidat

  // consultoras

  float *primer_candidat_clt(); // retorna el primer candidat al cluster.
  float *seguent_candidat_clt(int validacio); // cerca el seguent candidat i li
                                              // pasan la validacio del ultim
  int n_punts();                              // n punts del cluster
  void donar_max_min_xomig(float **mx, float **mn, float **xm, float *s_d);

  // consultoras  amb punter

  node* resetpt() { return topleft->seg[DRETA]; }
  void *noend(node *pt);
  float *llpt(node *pt);
  void advpt(node **pt);

};
