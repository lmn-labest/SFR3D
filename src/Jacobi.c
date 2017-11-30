#include<Jacobi.h>
/***********************************************************************
 * Data de criacao    : 14/05/2016                                     *
 * Data de modificaco : 28/11/2017                                     * 
 * ------------------------------------------------------------------  * 
 * CYCLIC_JACOBI : Calculo de autovalores e autovetores de matrizes    * 
 * simetricas pelo metodo de JACOBI.                                   *
 * ------------------------------------------------------------------  * 
 * Parametros de entrada:                                              *
 * ------------------------------------------------------------------  * 
 * a(n,n) - matriz simetrica                                           *
 * x(n,n) - nao definido                                               *
 * w(n)   - nao definido                                               *
 * ip(n)  - arranjo auxiliar (nao definido)                            *
 * n      - dimensao da matriz                                         *
 * tol    - tolerancia de convergencia                                 *
 * maxciclos - numero maximo de ciclos                                 *
 * cflag     -.true . ordena os autovalores e vetores em ordem         *
 *             crescente                                               *
 *            .false.ordena os autovalores e vetores em ordem          *
 *             decrescente                                             *
 * ------------------------------------------------------------------  * 
 * Parametros de saida:                                                *
 * ------------------------------------------------------------------  *
 * a - matriz espectral: a(i,i) = i-esimo autovalor                    *
 * x - matriz de autovetores: x(i,j),(i=1,n) = autovetor               *
 *                            associado ao j-esimo autovalor           *
 * w - w(i) = (a(i,i) = i-esimo autovalor)                             *
 * ------------------------------------------------------------------  * 
 * OBS:                                                                *
 * ------------------------------------------------------------------  * 
 * Versao: Cyclic-by-row Jacobi                                        *
 * fonte: Matrix Computation Golub e Van Loan pag. 480 Ed.4            * 
 ***********************************************************************/
void cyclic_jacobi(DOUBLE *RESTRICT a   , DOUBLE *RESTRICT x
                 , DOUBLE *RESTRICT w   , INT *RESTRICT ip
                 , INT const        n   , DOUBLE const tol
                 , short const maxciclos, bool const cflag) {
  short k;
  bool fail_conv = true;
  INT p,q;
  DOUBLE off,c,s,conv;
  DOUBLE i1,i2,i3;


/*... matriz indentidade*/
  indenty(x,n);
/*.....................................................................*/

  off = off_forbenius_norm_sym(a,n);
  if( off != 0.e0){
    conv = tol*forbenius_norm_sym(a,n) ;
/*... Ciclos de Jacobi:*/
    for(k=0;k<maxciclos;k++){
        for(p=0;p<n;p++){
/*...*/
          for(q=p+1;q<n;q++){
/* ... Calculo da rotacao para o coeficiente a(p,q): 
       2-by-2 Symmetric Schur Decomposition
*/
            i1 = MAT2D(q,q,a,n);
            i2 = MAT2D(p,p,a,n);
            i3 = MAT2D(p,q,a,n);
            symSchur2(i1,i2,i3,&c,&s);
/*.....................................................................*/
                 
/* ... Produto Xt A X:*/
            xax(a,c,s,n,p,q);

/* ... Produto (Xi)(Xi+1):*/                
            xx(x,c,s,n,p,q);
          }
/*.....................................................................*/
        }
/*.....................................................................*/
        off =  off_forbenius_norm_sym(a,n);
        if (off<conv) {
          fail_conv = false;
          break;
        }

      }
/*.....................................................................*/

/*...*/
      if(fail_conv){
        printf("CYCLIC_JACOBI: Sem convergecia\n %d: %e > %e!\n"
              ,k,off,conv);
        exit(EXIT_FAILURE);
      }
/*.....................................................................*/
  }

/*...*/
  ordena(a,x,w,ip,n,cflag);
/*.....................................................................*/
}
/***********************************************************************/

/***********************************************************************/
void ordena(DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
           ,DOUBLE *RESTRICT w,INT *RESTRICT ip
           ,INT const n       ,bool const fAsc ){

  INT i,j,k,iaux;
  DOUBLE aux;
  bool itroca;

/*...*/
  for( i = 0; i < n; i++){
    ip[i] = i;
    w[i]  = MAT2D(i,i,a,n);
  }
/*.....................................................................*/

/*...*/
  if(fAsc){
   do{
      itroca = false;
      for(i = 1; i < n; i++){
        j = i-1;
        if( w[j] > w[i]){
          aux   = w[i];
          w[i]  = w[j];
          w[j]  = aux;
          iaux  = ip[i];
          ip[i] = ip[j];
          ip[j] = iaux;
          for(k = 0; k < n; k++){
            aux            = MAT2D(k,i,x,n);
            MAT2D(k,i,x,n) = MAT2D(k,j,x,n);
            MAT2D(k,j,x,n) = aux;
          }
          itroca = true;
        }
      }
    }while(itroca);
  }
/*.....................................................................*/

/*...*/
  else{
    do{
      itroca = false;
      for(i = 1; i < n; i++){
        j = i-1;
        if( w[j] < w[i]){
          aux   = w[i];
          w[i]  = w[j];
          w[j]  = aux;
          iaux  = ip[i];
          ip[i] = ip[j];
          ip[j] = iaux;
          for(k = 0; k < n; k++){
            aux            = MAT2D(k,i,x,n);
            MAT2D(k,i,x,n) = MAT2D(k,j,x,n);
            MAT2D(k,j,x,n) = aux;
          }
          itroca = true;
        }
      }
    }while(itroca);
  }
/*.....................................................................*/
}
/***********************************************************************/

void indenty(DOUBLE *RESTRICT x, INT const n) {

  INT i,j;

  for (i = 0; i < n; i++) {
    MAT2D(i,i,x,n) = 1.e0;
    for (j = i+1; j < n; j++) {
      MAT2D(i,j,x,n) = 0.e0;
      MAT2D(j,i,x,n) = 0.e0;
    }
  }
}

void xax(DOUBLE *RESTRICT a, DOUBLE const c, DOUBLE const s
         , INT const n       , INT const p   , INT const q){

  INT i;
  DOUBLE a1,a2;

/*... Produto AX:*/
  for (i=0;i<n;i++){
    a1             = MAT2D(i,p,a,n);
    a2             = MAT2D(i,q,a,n);
    MAT2D(i,p,a,n) = a1*c - a2*s;
    MAT2D(i,q,a,n) = a1*s + a2*c;
  }
/*....................................................................*/              

/*... Produto Xt(AX):*/
  for (i=0;i<n;i++){
    a1             = MAT2D(p,i,a,n);
    a2             = MAT2D(q,i,a,n);
    MAT2D(p,i,a,n) = a1*c - a2*s;
    MAT2D(q,i,a,n) = a1*s + a2*c;
  }
/*....................................................................*/
}
/**********************************************************************/

void xx(DOUBLE *RESTRICT x, DOUBLE const c, DOUBLE const s
      , INT const n       , INT const p   , INT const q){

  INT i;
  DOUBLE a1,a2;

/*...*/
  for (i=0;i<n;i++){
    a1             = MAT2D(i,p,x,n);
    a2             = MAT2D(i,q,x,n);
    MAT2D(i,p,x,n) = a1*c - a2*s;
    MAT2D(i,q,x,n) = a1*s + a2*c;
  }
/*....................................................................*/
}
/**********************************************************************/

DOUBLE forbenius_norm_sym(DOUBLE *RESTRICT a,INT const n){


  INT i,j;
  DOUBLE x,norm,diag;
      
  norm = 0.e0;
  diag = 0.e0;
/*...*/
  for(i=0;i<n;i++){
    x     = MAT2D(i,i,a,n);
    diag += x*x;
    for(j=0;j<i;j++){
      x     = MAT2D(i,j,a,n);
      norm += x*x;
    }
  }
/*...................................................................*/

  return sqrt(diag+2.e0*norm);
}

DOUBLE off_forbenius_norm_sym(DOUBLE *RESTRICT a,INT const n){


  INT i,j;
  DOUBLE norm,x;
      
  norm = 0.e0;
/*...*/
  for(i=0;i<n;i++){
    for (j = 0; j < i; j++) {
      x     = MAT2D(i,j,a,n);
      norm += x*x;
    }
  }
/*...................................................................*/

  return sqrt(2.e0*norm);
}
/**********************************************************************/

 static void symSchur2(DOUBLE const aqq,DOUBLE const app,DOUBLE const apq
               ,DOUBLE *c       ,DOUBLE *s){

  DOUBLE r,d,t;
  t  = 0.e0;
  *c = 1.e0;
  *s = 0.e0;
  if(apq != 0.e0){
    r = 0.5e0*(aqq-app)/apq;
    d = sqrt(1.e0+r*r);
    if( r >= 0.e0)
      t = 1.e0/(r+d);
    else
      t = 1.e0/(r-d);
        
    *c = 1.e0/sqrt(1.e0+t*t);
    *s = t*(*c);
  } 
}
/**********************************************************************/