#include<CellLoop.h>
/********************************************************************* 
 * CELLDIF2D: Celula 2D para difusao pura                            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lx        -> coordenadas dos nos da celula central e seus viznhos * 
 * lnFace    -> numero de faces da celula central e seus vizinhos    * 
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lprop     -> propriedade fisicas das celulas                      * 
 * lViz      -> viznhos da celula central                            * 
 * xc        -> centroides das celulas                               * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * mEta      -> modulo do vetor eta                                  * 
 * normal    -> vetores normais as faces das celulas                 * 
 * area      -> area da celula central                               * 
 * xm        -> pontos medios das faces da celula cenral             * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * mkm       -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca      -> menor distancia do centroide central a faces desta   *
 *              celula                                               * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * u0        -> solucao conhecida                                    * 
 * faceR     -> restricoes por elmento                               * 
 * faceS     -> carga por elemento                                   * 
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellDif2D(double *restrict lx                        
              ,short *restrict lGeomType,double *restrict prop
              ,INT *restrict lViz       ,double *restrict xc                           
              ,double *restrict ksi     ,double *restrict mKsi
              ,double *restrict eta     ,double *restrict mEta
              ,double *restrict normal  ,double *restrict volume
              ,double *restrict xm      ,double *restrict xmcc
              ,double *restrict dcca    ,double *restrict mkm
              ,double *restrict lA      ,double *restrict lB
              ,short  *restrict lFaceR  ,double *restrict lFaceS
              ,double *restrict u0      ,const short nen     
              ,short const nFace        ,const short ndm
              ,INT const nel)
{ 

  double coefDifC,coefDif,coefDifV;
  double p,aP,sP,alpha,nk,dfd;
  double dPViz,v[2];
  short idCell = nFace;
  short nAresta;
  INT vizNel;

/*... propriedades da celula*/
  coefDifC = MAT2D(idCell,0,prop,MAXPROP);
/*...................................................................*/

  p  = 0.0e0;
  aP = 0.0e0;
  sP = 0.0e0;
  for(nAresta=0;nAresta<nFace;nAresta++){
    vizNel = lViz[nAresta];
/*... dominio*/
    if( vizNel  > -1 ){
/*...*/
      v[0]  = mkm[nAresta]*MAT2D(nAresta,0,eta,ndm) 
            + MAT2D(nAresta,0,xmcc,ndm);
      v[1]  = mkm[nAresta]*MAT2D(nAresta,1,eta,ndm) 
            + MAT2D(nAresta,1,xmcc,ndm);
      dPViz = sqrt(v[0]*v[0] + v[1]*v[1]);
      alpha = dPViz/mKsi[nAresta];
/*...................................................................*/

/*... media harmonica*/
      coefDifV =  MAT2D(nAresta,0,prop,MAXPROP); 
      coefDif  = alpha/coefDifC + (1.0e0-alpha)/coefDifV;
      coefDif  = 1.0e0/coefDif;
/*...................................................................*/

/*... produtos internos*/
      nk = MAT2D(nAresta,0,normal,ndm) * MAT2D(nAresta,0,ksi,ndm)
         + MAT2D(nAresta,1,normal,ndm) * MAT2D(nAresta,1,ksi,ndm);
/*...................................................................*/

/*... difusao direta*/
      dfd = (coefDif*mEta[nAresta])/(nk*mKsi[nAresta]);
/*...................................................................*/

/*...*/
      lA[nAresta] = dfd;
/*...................................................................*/
    }
/*... contorno*/
    else{
      lA[nAresta] = 0.0e0;
/*... fluxo prescrito*/
      if(lFaceR[nAresta] == 0){ 
        p +=  mEta[nAresta]*lFaceS[nAresta];
      }
/*... temperatura prescrita*/
      else{
        aP  = coefDifC*mEta[nAresta]/dcca[nAresta];
        sP += aP;
        p  += aP*lFaceS[nAresta]; 
      }
    }
/*...................................................................*/
  }

/*...*/
  lA[idCell] = sP;
/*residuo: R = F - KpUp*/ 
  p -= lA[idCell]*u0[idCell]; 
  for(nAresta=0;nAresta<nFace;nAresta++){
    lA[idCell] += lA[nAresta];
/*residuo (R = F-KvizUviz ) e valores prescritos por elemento*/
    p         += lA[nAresta]*u0[nAresta];  
  }
/*...................................................................*/

/*...*/
  lB[0] = p;
/*...................................................................*/
}
