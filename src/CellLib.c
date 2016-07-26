#include<CellLoop.h>
/********************************************************************* 
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 05/07/2016                                   * 
 *-------------------------------------------------------------------* 
 * CELLLIBSIMPLEVEl: chamada de bibliotecas de celulas para          *
 * problema de escoamento de fluidos (VEL)                           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loadsVel  -> definicoes de cargas de velocidades                  * 
 * loadsPres -> definicoes de cargas de pressao                      * 
 * advVel    -> tecnica da discretizacao do termo advecao            * 
 * typeSimple-> tipo do metodo simple                                *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            * 
 * lId       -> equa da celula                                       * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area das faces                                       * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> volume celula central                                * 
 * xm        -> pontos medios das faces da celula central            * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               * 
 * lDensity  -> massa especifica com variacao temporal               * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceVelR  -> restricoes por elemento de velocidades               * 
 * faceVelL  -> carga por elemento de velocidades                    * 
 * facePresR -> restricoes por elemento de pressao                   * 
 * facePresL -> carga por elemento de pressao                        * 
 * pres      -> campo de pressao conhecido                           * 
 * gradPes   -> gradiente reconstruido da pressao                    * 
 * vel       -> campo de velocidade conhecido                        * 
 * dField    -> matriz D do metodo simple                            * 
 * cc        -> centroides da celula centra e seus vizinhos          *
 * underU    -> fator underrelaxtion sinple                          * 
 * sPressure -> reconstrucao de segunda ordem para pressoes nas      *
 *              faces                                                *
 * nEn       -> numero de nos da celula central                      * 
 * nFace     -> numero de faces da celula central                    * 
 * ndm       -> numero de dimensoes                                  * 
 * lib       -> numero da biblioteca                                 * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellLibSimpleVel(Loads *loadsVel    ,Loads *loadsPres        
             ,Advection  advVel          ,short const typeSimple 
             ,short *restrict lGeomType  ,DOUBLE *restrict lprop
             ,INT   *restrict lViz       ,INT *restrict lId  
             ,DOUBLE *restrict ksi       ,DOUBLE *restrict mKsi
             ,DOUBLE *restrict eta       ,DOUBLE *restrict fArea
             ,DOUBLE *restrict normal    ,DOUBLE *restrict volume
             ,DOUBLE *restrict xm        ,DOUBLE *restrict xmcc
             ,DOUBLE *restrict dcca      ,DOUBLE *restrict lDensity
             ,DOUBLE *restrict vSkew     ,DOUBLE *restrict mvSkew
             ,DOUBLE *restrict lA        ,DOUBLE *restrict lB
             ,DOUBLE *restrict lRcell    ,Temporal const ddt
             ,short  *restrict lFaceVelR ,short  *restrict lFaceVelL       
             ,short  *restrict lFacePresR,short  *restrict lFacePresL             
             ,DOUBLE *restrict pres      ,DOUBLE *restrict gradPres 
             ,DOUBLE *restrict vel       ,DOUBLE *restrict gradVel
             ,DOUBLE *restrict dField    ,DOUBLE *restrict cc
             ,DOUBLE const underU        ,const bool sPressure
             ,short const nEn            ,short  const nFace     
             ,short const ndm            ,short const lib    
             ,INT const nel)
{

/*... */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2){
      cellSimpleVel2D(loadsVel,loadsPres   
                 ,advVel      ,typeSimple
                 ,lGeomType   ,lprop
                 ,lViz        ,lId
                 ,ksi         ,mKsi
                 ,eta         ,fArea
                 ,normal      ,volume
                 ,xm          ,xmcc
                 ,dcca        ,lDensity 
                 ,vSkew       ,mvSkew
                 ,lA          ,lB
                 ,lRcell      ,ddt 
                 ,lFaceVelR   ,lFaceVelL
                 ,lFacePresR  ,lFacePresL
                 ,pres        ,gradPres 
                 ,vel         ,gradVel
                 ,dField      ,cc
                 ,underU      ,sPressure                         
                 ,nEn         ,nFace 
                 ,ndm         ,nel);  
    }
/*..................................................................*/   

/*... 3D*/
    else if(ndm == 3){
      cellSimpleVel3D(loadsVel,loadsPres   
                 ,advVel      ,typeSimple
                 ,lGeomType   ,lprop
                 ,lViz        ,lId
                 ,ksi         ,mKsi
                 ,eta         ,fArea
                 ,normal      ,volume
                 ,xm          ,xmcc
                 ,dcca        ,lDensity 
                 ,vSkew       ,mvSkew
                 ,lA          ,lB
                 ,lRcell      ,ddt 
                 ,lFaceVelR   ,lFaceVelL
                 ,lFacePresR  ,lFacePresL
                 ,pres        ,gradPres 
                 ,vel         ,gradVel
                 ,dField      ,cc
                 ,underU      ,sPressure                         
                 ,nEn         ,nFace 
                 ,ndm         ,nel);  
    } 
/*..................................................................*/   
  }

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 01/07/2016                                   *
 * Data de modificaco : 09/07/2016                                   * 
 *-------------------------------------------------------------------* 
 * CELLLIBSIMPLEPRES : chamada de bibliotecas de celulas para        *
 * problema de escoamento de fluidos (PRES)                          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loadsVel  -> definicoes de cargas de velocidades                  * 
 * loadsPres -> definicoes de cargas de pressao                      * 
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            * 
 * lId       -> equa da celula                                       * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area das faces                                       * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> volume celula central                                * 
 * xm        -> pontos medios das faces da celula central            * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               * 
 * lDensity  -> massa especifica com variacao temporal               * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceVelR  -> restricoes por elemento de velocidades               * 
 * faceVelL  -> carga por elemento de velocidades                    * 
 * facePresR -> restricoes por elemento de pressao                   * 
 * facePresL -> carga por elemento de pressao                        * 
 * pres      -> campo de pressao conhecido                           * 
 * gradPes   -> gradiente reconstruido da pressao                    * 
 * vel       -> campo de velocidade conhecido                        * 
 * dField    -> matriz D do metodo simple                            * 
 * nEn       -> numero de nos da celula central                      * 
 * nFace     -> numero de faces da celula central                    * 
 * ndm       -> numero de dimensoes                                  * 
 * lib       -> numero da biblioteca                                 * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellLibSimplePres(Loads *loadsVel    ,Loads *loadsPres          
               ,short *restrict lGeomType,DOUBLE *restrict lprop
               ,INT   *restrict lViz     ,INT *restrict lId  
               ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
               ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea
               ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
               ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
               ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
               ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
               ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
               ,DOUBLE *restrict lRcell  ,Temporal const ddt
               ,short  *restrict lFaceVelR ,short  *restrict lFaceVelL       
               ,short  *restrict lFacePresR,short  *restrict lFacePresL             
               ,DOUBLE *restrict pres    ,DOUBLE *restrict gradPres 
               ,DOUBLE *restrict vel     ,DOUBLE *restrict dField 
               ,short const nEn          ,short  const nFace     
               ,short const ndm          ,short const lib    
               ,INT const nel)
{

/*... */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2){
      cellSimplePres2D(loadsVel,loadsPres
                 ,lGeomType,lprop
                 ,lViz     ,lId
                 ,ksi      ,mKsi
                 ,eta      ,fArea
                 ,normal   ,volume
                 ,xm       ,xmcc
                 ,dcca     ,lDensity 
                 ,vSkew    ,mvSkew
                 ,lA       ,lB
                 ,lRcell    
                 ,lFaceVelR   ,lFaceVelL
                 ,lFacePresR  ,lFacePresL
                 ,pres     ,gradPres 
                 ,vel      ,dField                                     
                 ,nEn      ,nFace  
                 ,ndm      ,nel);  
    }
/*..................................................................*/   

/*... 3D*/
    else if(ndm == 3){
      cellSimplePres3D(loadsVel,loadsPres
                 ,lGeomType,lprop
                 ,lViz     ,lId
                 ,ksi      ,mKsi
                 ,eta      ,fArea
                 ,normal   ,volume
                 ,xm       ,xmcc
                 ,dcca     ,lDensity 
                 ,vSkew    ,mvSkew
                 ,lA       ,lB
                 ,lRcell    
                 ,lFaceVelR   ,lFaceVelL
                 ,lFacePresR  ,lFacePresL
                 ,pres     ,gradPres 
                 ,vel      ,dField                                     
                 ,nEn      ,nFace  
                 ,ndm      ,nel);  
    } 
/*..................................................................*/   
  }

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * CELLLIBTRANS : chamada de bibliotecas de celulas para o problema  *
 * de transporte.                                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads     -> definicoes de cargas                                 * 
 * advT      -> tecnica da discretizacao do termo advecao            * 
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            * 
 * lId       -> equa da celula                                       * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area das faces                                       * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> volume celula central                                * 
 * xm        -> pontos medios das faces da celula central            * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               * 
 * lDensity  -> massa especifica com variacao temporal               * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceR     -> restricoes por elemento                              * 
 * faceL     -> carga por elemento                                   * 
 * u0        -> solucao conhecida                                    * 
 * gradU0    -> gradiente rescontruido da solucao conhecida          * 
 * vel       -> campo de velocidade conhecido                        * 
 * nEn       -> numero de nos da celula central                      * 
 * nFace     -> numero de faces da celula central                    * 
 * ndm       -> numero de dimensoes                                  * 
 * lib       -> numero da biblioteca                                 * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellLibTrans(Loads *loads           ,Advection advT 
               ,short *restrict lGeomType,DOUBLE *restrict lprop
               ,INT   *restrict lViz     ,INT *restrict lId  
               ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
               ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea
               ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
               ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
               ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
               ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
               ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
               ,DOUBLE *restrict lRcell  ,Temporal const ddt
               ,short  *restrict lFaceR  ,short  *restrict lFaceL       
               ,DOUBLE *restrict u0      ,DOUBLE *restrict gradU0
               ,DOUBLE *restrict vel                                 
               ,short const nEn          ,short  const nFace     
               ,short const ndm          ,short const lib    
               ,INT const nel)
{

/*... */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2){
      cellTrans2D(loads    ,advT
                 ,lGeomType,lprop
                 ,lViz     ,lId
                 ,ksi      ,mKsi
                 ,eta      ,fArea
                 ,normal   ,volume
                 ,xm       ,xmcc
                 ,dcca     ,lDensity 
                 ,vSkew    ,mvSkew
                 ,lA       ,lB
                 ,lRcell   ,ddt 
                 ,lFaceR   ,lFaceL
                 ,u0       ,gradU0      
                 ,vel                   
                 ,nEn      ,nFace 
                 ,ndm      ,nel);
    }
/*..................................................................*/   

/*... 3D*/
    else if(ndm == 3){
      cellTrans3D(loads    ,advT
                 ,lGeomType,lprop
                 ,lViz     ,lId
                 ,ksi      ,mKsi
                 ,eta      ,fArea
                 ,normal   ,volume
                 ,xm       ,xmcc
                 ,dcca     ,lDensity 
                 ,vSkew    ,mvSkew
                 ,lA       ,lB
                 ,lRcell   ,ddt 
                 ,lFaceR   ,lFaceL
                 ,u0       ,gradU0      
                 ,vel                   
                 ,nEn      ,nFace 
                 ,ndm      ,nel);
    } 
/*..................................................................*/   
  }

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * CELLLIBDIF : chamada de bibliotecas de celulas para o problema    *
 * de difusao.                                                       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads     -> definicoes de cargas                                 * 
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            * 
 * lId       -> equa da celula                                       * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area das faces                                       * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> volume celula central                                * 
 * xm        -> pontos medios das faces da celula central            * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               * 
 * lDensity  -> massa especifica com variacao temporal               * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * u0        -> solucao conhecida                                    * 
 * gradU0    -> gradiente rescontruido da solucao conhecida          * 
 * faceR     -> restricoes por elemento                              * 
 * faceLd1   -> carga por elemento                                   * 
 * nEn       -> numero de nos da celula central                      * 
 * nFace     -> numero de faces da celula central                    * 
 * ndm       -> numero de dimensoes                                  * 
 * lib       -> numero da biblioteca                                 * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellLibDif(Loads *loads
               ,short *restrict lGeomType,DOUBLE *restrict lprop
               ,INT   *restrict lViz     ,INT *restrict lId  
               ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
               ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea
               ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
               ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
               ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
               ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
               ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
               ,DOUBLE *restrict lRcell  ,Temporal const ddt
               ,short  *restrict lFaceR  ,short  *restrict lFaceL       
               ,DOUBLE *restrict u0      ,DOUBLE *restrict gradU0
               ,short const nEn          ,short  const nFace     
               ,short const ndm          ,short const lib    
               ,INT const nel)
{

/*... difusao pura */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2){
      cellDif2D(loads
               ,lGeomType,lprop
               ,lViz     ,lId
               ,ksi      ,mKsi
               ,eta      ,fArea
               ,normal   ,volume
               ,xm       ,xmcc
               ,dcca     ,lDensity 
               ,vSkew    ,mvSkew
               ,lA       ,lB
               ,lRcell   ,ddt 
               ,lFaceR   ,lFaceL
               ,u0       ,gradU0      
               ,nEn      ,nFace 
               ,ndm      ,nel);
    }
/*..................................................................*/   

/*... 3D*/
    else if(ndm == 3){
      cellDif3D(loads     
               ,lGeomType,lprop
               ,lViz     ,lId
               ,ksi      ,mKsi
               ,eta      ,fArea
               ,normal   ,volume
               ,xm       ,xmcc
               ,dcca     ,lDensity
               ,vSkew    ,mvSkew
               ,lA       ,lB
               ,lRcell   ,ddt
               ,lFaceR   ,lFaceL 
               ,u0       ,gradU0      
               ,nEn      ,nFace 
               ,ndm      ,nel);
    } 
/*..................................................................*/   
  }

}
/*********************************************************************/
 
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * CELLGEOM2D : calculo geometrico de propriedade de celula 2D       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lx        -> coordenadas dos nos da celula central e seus viznhos * 
 * lnFace    -> numero da faces da celula central e seus vizinhos    * 
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * xc        -> indefenido                                           * 
 * ksi       -> indefinido                                           * 
 * mksi      -> indefinido                                           * 
 * eta       -> indefinido                                           * 
 * meta      -> indefinido                                           * 
 * normal    -> indefinido                                           * 
 * area      -> indefinido                                           * 
 * xm        -> indefinido                                           * 
 * xmcc      -> indefinido                                           * 
 * dcca      -> indefinido                                           * 
 * vSkew     -> indefinido                                           * 
 * mvSkew    -> indefinido                                           * 
 * sn        -> numeracao dos nos por aresta                         * 
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * numel     -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * xc     -> centroide das celulas                                   * 
 * Ksi    -> vetores que unem centroide da celula central aos        *
 *            vizinhos destas                                        * 
 * mKsi   -> modulo do vetor ksi                                     * 
 * eta    -> vetores paralelos as faces das celulas                  * 
 * fArea  -> modulo do vetor eta                                     * 
 * normal -> vetores normais as faces das celulas                    * 
 * area   -> area da celula central                                  * 
 * xm     -> pontos medios das faces da celula cenral                * 
 * xmcc   -> vetores que unem o centroide aos pontos medios das      * 
 *            faces da celula central                                * 
 * vSkew  -> vetor entre o ponto medio a intersecao que une os       * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew -> distacia entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca   -> menor distacia do centroide central a faces desta celula* 
 *-------------------------------------------------------------------* 
 * OBS: xc | x1 x2 x3 | -> celula vizinha da aresta 1                * 
 *         | x1 x2 x3 | -> celula vizinha da aresta 2                * 
 *         | x1 x2 x3 | -> celula vizinha da aresta 3                * 
 *         |   ...    | ->   ...                                     * 
 *         | x1 x2 x3 | -> celula central                            * 
 *                                                                   * 
 *     lx(vizinho,numero do no do vizinho, dimensao)                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellGeom2D(DOUBLE *restrict lx       ,short *restrict lnFace
               ,short  *restrict lGeomType,DOUBLE *restrict xc
               ,DOUBLE *restrict ksi      ,DOUBLE *restrict mksi
               ,DOUBLE *restrict eta      ,DOUBLE *restrict meta
               ,DOUBLE *restrict normal   ,DOUBLE *restrict area
               ,DOUBLE *restrict xm       ,DOUBLE *restrict xmcc
               ,DOUBLE *restrict dcca     
               ,DOUBLE *restrict vSkew    ,DOUBLE *restrict mvSkew
               ,short  *restrict sn
               ,short const maxNo         ,short const maxViz
               ,short const ndm           ,INT const nel)
{
  short i,j,k,cCell = maxViz;
  short no1,no2;
  double x1,x2;

  for(i=0;i<lnFace[cCell];i++){
    mksi[i]  = 0.0e0;
    meta[i]  = 0.0e0;
    dcca[i]  = 0.0e0;
    MAT2D(i,0,ksi,ndm)= 0.0e0; 
    MAT2D(i,1,ksi,ndm)= 0.0e0; 
    MAT2D(i,0,eta,ndm)= 0.0e0; 
    MAT2D(i,1,eta,ndm)= 0.0e0; 
  }
  
/*... calculo do centro geometrico das celulas*/      
  for(i=0;i<=maxViz;i++){
    for(j=0;j<ndm;j++){
      MAT2D(i,j,xc,ndm) = 0.0e0;
      if(lnFace[i]){
        for(k=0;k<lnFace[i];k++){
          MAT2D(i,j,xc,ndm) += MAT3D(i,k,j,lx,maxNo,ndm);
        }
        MAT2D(i,j,xc,ndm)*=(1.0e0/((DOUBLE)lnFace[i]));
      }
    }
  }
/*...................................................................*/


/*... vetor que une o centroide da celula ao vizinho(ksi)*/
  for(i=0;i<lnFace[cCell];i++){
    if(lnFace[i]){
      for(j=0;j<ndm;j++){
        MAT2D(i,j,ksi,ndm) = 
        MAT2D(i,j,xc,ndm) - MAT2D(cCell,j,xc,ndm);
        mksi[i] += MAT2D(i,j,ksi,ndm)*MAT2D(i,j,ksi,ndm); 
      }
      mksi[i] = sqrt(mksi[i]);
    }   
  }
/*...................................................................*/


/*... vetor paralelo as arestas (eta)*/
  for(i=0;i<lnFace[cCell];i++){
    no1 = MAT2D(i,0,sn,2); 
    no2 = MAT2D(i,1,sn,2);
    for(j=0;j<ndm;j++){
      x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
      x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
      MAT2D(i,j,eta,ndm) = x2 - x1;
    }
    *area = areaCell(eta,lGeomType[cCell],ndm,nel);
  }
/*...................................................................*/

/*... modulo do vetor das arestas (eta)*/
  for(i=0;i<lnFace[cCell];i++){
    for(j=0;j<ndm;j++){
      meta[i] += MAT2D(i,j,eta,ndm)*MAT2D(i,j,eta,ndm); 
    }
    meta[i] = sqrt(meta[i]);
  }
/*...................................................................*/

/*... calculo dos vetores unitarios*/
  for(i=0;i<lnFace[cCell];i++){
    for(j=0;j<ndm;j++){
      MAT2D(i,j,eta,ndm) = MAT2D(i,j,eta,ndm) / meta[i];
      if(lnFace[i]) 
        MAT2D(i,j,ksi,ndm) = MAT2D(i,j,ksi,ndm) / mksi[i];
    } 
  }
/*...................................................................*/

/*... calculo do vetor normal unitario*/
  for(i=0;i<lnFace[cCell];i++){
    MAT2D(i,0,normal,ndm) =  MAT2D(i,1,eta,ndm);
    MAT2D(i,1,normal,ndm) = -MAT2D(i,0,eta,ndm);
  }
/*...................................................................*/  

/*... ponto medio da aresta(xm) e vetor que une este ponto ao 
      centroide(mxcc) */
  for(i=0;i<lnFace[cCell];i++){
    no1 = MAT2D(i,0,sn,2); 
    no2 = MAT2D(i,1,sn,2);
    for(j=0;j<ndm;j++){
      x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
      x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
/*...*/
      MAT2D(i,j,xm  ,ndm) = 0.5e0*(x1+x2);
/*...................................................................*/  

/*...*/
      MAT2D(i,j,xmcc,ndm) = MAT2D(i,j,xm,ndm) - MAT2D(cCell,j,xc,ndm);   
/*...................................................................*/  

/*...*/
      dcca[i] += MAT2D(i,j,xmcc,ndm)*MAT2D(i,j,normal,ndm); 
/*...................................................................*/  
    }
  }
/*...................................................................*/  

/*... calculo do centro da celula fantasma*/      
  for(i=0;i<=maxViz;i++){
    for(j=0;j<ndm;j++){
      if(!lnFace[i]){
        MAT2D(i,j,xc,ndm) = MAT2D(cCell,j,xc,ndm) 
                          + 2.0e0*dcca[i]*MAT2D(i,j,normal,ndm);
      }
    }
  }
/*...................................................................*/

/*... vetor que une o centroide da celula ao vizinho fantasma(ksi)*/
  for(i=0;i<lnFace[cCell];i++){
    if(!lnFace[i]){
      for(j=0;j<ndm;j++){
        MAT2D(i,j,ksi,ndm) = MAT2D(i,j,xmcc,ndm);
        mksi[i] += MAT2D(i,j,ksi,ndm)*MAT2D(i,j,ksi,ndm); 
      }
      mksi[i] = sqrt(mksi[i]);
      for(j=0;j<ndm;j++)
        MAT2D(i,j,ksi,ndm) /= mksi[i];
    }   
  }
/*...................................................................*/

/*...*/
  vectorKm2d(lx     , xc
            ,xm     
            ,vSkew  ,mvSkew
            ,sn     ,lnFace[cCell]
            ,maxViz ,maxNo       
            ,ndm    ,nel);
/*...................................................................*/  
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * CELLGEOM3D : calculo geometrico de propriedade de celula 3D       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lx        -> coordenadas dos nos da celula central e seus viznhos * 
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lnEn      -> numero da de nos da celula central e seus vizinhos   * 
 * lnFace    -> numero da faces da celula central e seus vizinhos    * 
 * xc        -> indefenido                                           * 
 * ksi       -> indefinido                                           * 
 * mksi      -> indefinido                                           * 
 * eta       -> indefinido                                           * 
 * fArea     -> indefinido                                           * 
 * normal    -> indefinido                                           * 
 * area      -> indefinido                                           * 
 * xm        -> indefinido                                           * 
 * xmcc      -> indefinido                                           * 
 * dcca      -> indefinido                                           * 
 * mvSkew    -> indefinido                                           * 
 * sn        -> numeracao dos nos por aresta                         * 
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * numel     -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * xc     -> centroide das celulas                                   * 
 * Ksi    -> vetores que unem centroide da celula central aos        *
 *            vizinhos destas                                        * 
 * mKsi   -> modulo do vetor ksi                                     * 
 * eta    -> vetores paralelos as faces das celulas                  * 
 * fArea  -> area da face                                            * 
 * normal -> vetores normais as faces das celulas                    * 
 * volume -> volume da celula central                                * 
 * xm     -> pontos medios das faces da celula cenral                * 
 * xmcc   -> vetores que unem o centroide aos pontos medios das      * 
 *            faces da celula central                                * 
 * vSkew  -> vetor entre o ponto medio a intersecao que une os       * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew -> distacia entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca   -> menor distacia do centroide central a faces desta celula* 
 *-------------------------------------------------------------------* 
 * OBS: xc | x1 x2 x3 | -> celula vizinha da aresta 1                * 
 *         | x1 x2 x3 | -> celula vizinha da aresta 2                * 
 *         | x1 x2 x3 | -> celula vizinha da aresta 3                * 
 *         |   ...    | ->   ...                                     * 
 *         | x1 x2 x3 | -> celula central                            * 
 *                                                                   * 
 *     lx(vizinho,numero do no do vizinho, dimensao)                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellGeom3D(DOUBLE *restrict lx       ,short  *restrict lGeomType
               ,short *restrict lnFace    ,short *restrict lnEn      
               ,DOUBLE *restrict xc        
               ,DOUBLE *restrict ksi      ,DOUBLE *restrict mksi
               ,DOUBLE *restrict eta      ,DOUBLE *restrict fArea
               ,DOUBLE *restrict normal   ,DOUBLE *restrict volume
               ,DOUBLE *restrict xm       ,DOUBLE *restrict xmcc
               ,DOUBLE *restrict dcca     
               ,DOUBLE *restrict vSkew    ,DOUBLE *restrict mvSkew
               ,short  *restrict sn                         
               ,short const maxNo         ,short const maxViz
               ,short const ndm           ,INT const nel)
{
  short i,j,k,cCell = maxViz;
  short no1,no2,no3,no4,nenFaces=0;
  DOUBLE x1,x2,x3,x4,v1[3],v2[3],n[3],dot;

  for(i=0;i<lnFace[cCell];i++){
    mksi[i]  = 0.0e0;
    fArea[i] = 0.0e0;
    dcca[i]  = 0.0e0;
    MAT2D(i,0,ksi,ndm)= 0.0e0; 
    MAT2D(i,1,ksi,ndm)= 0.0e0; 
    MAT2D(i,2,ksi,ndm)= 0.0e0; 
    MAT2D(i,0,eta,ndm)= 0.0e0; 
    MAT2D(i,1,eta,ndm)= 0.0e0; 
    MAT2D(i,2,eta,ndm)= 0.0e0; 
  }
 
  if( lGeomType[cCell] == HEXACELL )
    nenFaces = 4;  
  else if( lGeomType[cCell] == TETRCELL )  
    nenFaces = 3;  
  

/*... calculo do centro geometrico das celulas*/    
  for(i=0;i<=maxViz;i++){
    for(j=0;j<ndm;j++){
      MAT2D(i,j,xc,ndm) = 0.0e0;
      if(lnFace[i]){
        for(k=0;k<lnEn[i];k++){
          MAT2D(i,j,xc,ndm) += MAT3D(i,k,j,lx,maxNo,ndm);
        }
        MAT2D(i,j,xc,ndm)*=(1.0e0/((DOUBLE)lnEn[i]));
      }
    }
  }
/*...................................................................*/

/*... vetor que une o centroide da celula ao vizinho(ksi)*/
  for(i=0;i<lnFace[cCell];i++){
    if(lnFace[i]){
      for(j=0;j<ndm;j++){
        MAT2D(i,j,ksi,ndm) = 
        MAT2D(i,j,xc,ndm) - MAT2D(cCell,j,xc,ndm);
        mksi[i] += MAT2D(i,j,ksi,ndm)*MAT2D(i,j,ksi,ndm); 
      }
      mksi[i] = sqrt(mksi[i]);
      for(j=0;j<ndm;j++)
        MAT2D(i,j,ksi,ndm) /=mksi[i];
    }   
  }
/*...................................................................*/

/*... calculo do vetor normal e area da face*/
  if( lGeomType[cCell] == HEXACELL ){
    for(i=0;i<lnFace[cCell];i++){
      no1 = MAT2D(i,0,sn,nenFaces);
      no2 = MAT2D(i,1,sn,nenFaces);
      no3 = MAT2D(i,2,sn,nenFaces);
      no4 = MAT2D(i,3,sn,nenFaces);
      for(j=0;j<ndm;j++){
        x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
        x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
        x3 = MAT3D(cCell,no3,j,lx,maxNo,ndm);
        x4 = MAT3D(cCell,no4,j,lx,maxNo,ndm);
        v1[j] = MAT2D(0,j,eta,ndm) = x2 - x1;
        MAT2D(1,j,eta,ndm) = x3 - x2;
        MAT2D(2,j,eta,ndm) = x4 - x3;
        v2[j] = MAT2D(3,j,eta,ndm) = x1 - x4;
      }
/*... area da face*/
      fArea[i] = areaCell(eta,QUADCELL,ndm,nel);
/*... normal a face*/
      prodVet(v2,v1,n);
/*... normal do vetor normal*/
      dot = 0.e0;
      for(j=0;j<ndm;j++){
        dot += n[j]*n[j];                             
      }
      dot = sqrt(dot);
/*... vetor unitario normal*/
      for(j=0;j<ndm;j++)
        MAT2D(i,j,normal,ndm) =  n[j]/dot;
    }
  }
/*...................................................................*/

/*... calculo do vetor normal e area da face*/
  else if( lGeomType[cCell] == TETRCELL ){
    for(i=0;i<lnFace[cCell];i++){
      no1 = MAT2D(i,0,sn,nenFaces);
      no2 = MAT2D(i,1,sn,nenFaces);
      no3 = MAT2D(i,2,sn,nenFaces);
      for(j=0;j<ndm;j++){
        x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
        x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
        x3 = MAT3D(cCell,no3,j,lx,maxNo,ndm);
        v1[j] = MAT2D(0,j,eta,ndm) = x2 - x1;
        v2[j] = MAT2D(1,j,eta,ndm) = x3 - x1;
      }
/*... area da face*/
      fArea[i] = areaCell(eta,TRIACELL,ndm,nel);
/*... normal a face*/
      prodVet(v1,v2,n);
/*... normal do vetor normal*/
      dot = 0.e0;
      for(j=0;j<ndm;j++){
        dot += n[j]*n[j];                             
      }
      dot = sqrt(dot);
/*... vetor unitario normal*/
      for(j=0;j<ndm;j++)
        MAT2D(i,j,normal,ndm) =  n[j]/dot;
    }
  }
/*...................................................................*/


/*... ponto medio da aresta(xm) e vetor que une este ponto ao 
      centroide(mxcc) */
  if( lGeomType[cCell] == HEXACELL ){
    for(i=0;i<lnFace[cCell];i++){
      no1 = MAT2D(i,0,sn,nenFaces);
      no2 = MAT2D(i,1,sn,nenFaces);
      no3 = MAT2D(i,2,sn,nenFaces);
      no4 = MAT2D(i,3,sn,nenFaces);
      for(j=0;j<ndm;j++){
        x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
        x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
        x3 = MAT3D(cCell,no3,j,lx,maxNo,ndm);
        x4 = MAT3D(cCell,no4,j,lx,maxNo,ndm);
/*...*/
        MAT2D(i,j,xm  ,ndm) = 0.25e0*(x1+x2+x3+x4);
/*...................................................................*/  

/*...*/
        MAT2D(i,j,xmcc,ndm) = MAT2D(i,j,xm,ndm) - MAT2D(cCell,j,xc,ndm);   
/*...................................................................*/  

/*...*/
        dcca[i] += MAT2D(i,j,xmcc,ndm)*MAT2D(i,j,normal,ndm); 
/*...................................................................*/  
      }
    }
/*...................................................................*/  
  }
/*...................................................................*/  

  
  else if( lGeomType[cCell] == TETRCELL ){
    for(i=0;i<lnFace[cCell];i++){
      no1 = MAT2D(i,0,sn,nenFaces);
      no2 = MAT2D(i,1,sn,nenFaces);
      no3 = MAT2D(i,2,sn,nenFaces);
      for(j=0;j<ndm;j++){
        x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
        x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
        x3 = MAT3D(cCell,no3,j,lx,maxNo,ndm);
/*...*/
        MAT2D(i,j,xm  ,ndm) =  oneDivTree*(x1+x2+x3);
/*...................................................................*/  

/*...*/
        MAT2D(i,j,xmcc,ndm) = MAT2D(i,j,xm,ndm) - MAT2D(cCell,j,xc,ndm);   
/*...................................................................*/  

/*...*/
        dcca[i] += MAT2D(i,j,xmcc,ndm)*MAT2D(i,j,normal,ndm); 
/*...................................................................*/  
      }
    }
/*...................................................................*/  
  }
/*...................................................................*/  
    
/*... calculo do centro da celula fantasma*/      
/*  for(i=0;i<=maxViz;i++){
    for(j=0;j<ndm;j++){
      if(!lnFace[i]){
        MAT2D(i,j,xc,ndm) = MAT2D(cCell,j,xc,ndm) 
                          + 2.0e0*dcca[i]*MAT2D(i,j,normal,ndm);
      }
    }
  }*/
/*...................................................................*/

/*... vetor que une o centroide da celula ao vizinho fantasma(ksi)*/
  for(i=0;i<lnFace[cCell];i++){
    if(!lnFace[i]){
      mksi[i] = 0.e0;
      for(j=0;j<ndm;j++){
        MAT2D(i,j,ksi,ndm) = MAT2D(i,j,xmcc,ndm);
        mksi[i] += MAT2D(i,j,ksi,ndm)*MAT2D(i,j,ksi,ndm); 
      }
      mksi[i] = sqrt(mksi[i]);
      for(j=0;j<ndm;j++)
        MAT2D(i,j,ksi,ndm) /= mksi[i];
    }   
  }
/*...................................................................*/

/*... volume da celula*/
  *volume = volume3DGreenGauss(xm,normal,fArea,lnFace[cCell]);

/*...*/
   vectorKm3d(xc            ,xm
             ,ksi           ,normal
             ,vSkew         ,mvSkew     
             ,lnFace[cCell] ,ndm
             ,maxViz        ,nel);
/*...................................................................*/  
}
/*********************************************************************/ 
    
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * CELLRCGRAD : chamada de bibliotecas de celulas para a reconstrucao*
 * de gradiente.                                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads     -> definicoes de cargas                                 * 
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lSquare -> matriz para a reconstrucao least Square                * 
 * lSquareR-> fatoracao R (RCLSQUAREQR)                              * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area da face                                         * 
 * normal    -> vetores normais as faces das celulas                 * 
 * area      -> area da celula central                               * 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face                      * 
 * xm        -> pontos medios das faces das celulas                  * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               * 
 * u         -> solucao conhecida na celula                          * 
 * gradU     -> gradiente rescontruido da solucao conhecida          * 
 * nU        -> solucao conhecida no no                              * 
 * faceR     -> restricoes por elmento                               * 
 * faceL     -> carga por elemento                                   * 
 * nFace     -> numero de faces                                      * 
 * ndm       -> numero de dimensoes                                  * 
 * lib       -> tipo de reconstrucao de gradiente                    * 
 * ndf       -> grauss de liberdade                                  * 
 * isNod     -> numeracao dos nos por aresta                         * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * gradU     -> gradiente calculodo                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellLibRcGrad(Loads *loads
                 ,INT   *restrict lViz    ,DOUBLE *restrict lProp    
                 ,DOUBLE *restrict lLsquare,DOUBLE *restrict lLsquareR
                 ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
                 ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea
                 ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
                 ,DOUBLE *restrict vSkew   
                 ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
                 ,DOUBLE *restrict lDcca 
                 ,short  *restrict lFaceR  ,short *restrict lFaceL
                 ,DOUBLE *restrict u       ,DOUBLE *restrict gradU 
                 ,DOUBLE *restrict lnU     ,short const ty                
                 ,short const nFace        ,short const ndm      
                 ,short const lib          ,short const ndf
                 ,short *restrict  isNod   ,INT const nel){
  long aux;
    
  switch(lib){
/*... green-Gauss linear baseado na celula*/  
    case RCGRADGAUSSC:
      greenGaussCell(loads
                    ,lViz    ,mKsi
                    ,lProp   ,lDcca
                    ,eta     ,fArea
                    ,normal  ,volume
                    ,vSkew   
                    ,xm      ,xmcc
                    ,lFaceR  ,lFaceL
                    ,u       ,gradU 
                    ,nFace   ,ndm   
                    ,ndf     ,nel);
    break;
/*...................................................................*/ 
    
/*... green-Gauss linear baseado no no*/  
    case RCGRADGAUSSN:
       greenGaussNode(lViz    ,fArea
                     ,normal  ,volume
                     ,lnU     ,gradU 
                     ,isNod   
                     ,nFace   ,ndm   
                     ,ndf     ,ty);
    break;
/*...................................................................*/ 
   
/*... minimos quadrados*/  
    case RCLSQUARE:
      leastSquare(loads
                 ,lLsquare,lViz
                 ,xm      ,xmcc
                 ,lProp   ,lDcca
                 ,u       ,gradU
                 ,lFaceR  ,lFaceL
                 ,nFace   ,ndf
                 ,ndm     ,nel);
    break;
/*...................................................................*/ 

/*... minimos quadrados QR*/  
    case RCLSQUAREQR:
      leastSquareQR(loads
                   ,lLsquare,lLsquareR
                   ,lProp   ,lDcca
                   ,lViz    
                   ,xm      ,xmcc
                   ,u       ,gradU
                   ,lFaceR  ,lFaceL
                   ,nFace   ,ndf
                   ,ndm);
    break;
/*...................................................................*/ 

/*...*/
    default: 
      aux = (long) nel;
      printf("Erro: tipo de reconstrucao nao exitentes. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
               ,aux,__FILE__,__func__);
      exit(EXIT_FAILURE);
    break; 
/*...................................................................*/ 

  }  

}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * GRREENGAUSSCELL: reconstrucao de gradiente green-gauss linear por * 
 * celula                                                            *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads     -> definicoes de cargas                                 * 
 * lViz      -> viznhos da celula central                            * 
 * lProp     -> propriedades dos material                            *
 * lDcca     -> menor distancia do centroide a faces desta celula    * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area da face                                         * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> area da celula central                               * 
 * vSkew      -> vetor entre o ponto medio a intersecao que une os   * 
 *            centrois compartilhado nessa face                      * 
 * xm        -> pontos medios das faces das celulas                  * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * lFaceR    -> restricoes por elmento                               * 
 * lFaceS    -> carga por elemento                                   * 
 * u         -> solucao conhecida                                    * 
 * gradU     -> gradiente rescontruido da solucao conhecida          * 
 * nFace     -> numero vizinhos por celula maximo da malha           * 
 * ndf       -> grauss de liberdade                                  * 
 * ndm       -> numero de dimensoes                                  * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * gradU     -> gradiente calculado                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void greenGaussCell(Loads *loads
               ,INT *restrict lViz   ,DOUBLE *restrict mKsi
               ,DOUBLE *restrict lProp   ,DOUBLE *restrict lDcca 
               ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea
               ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
               ,DOUBLE *restrict vSkew   
               ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc 
               ,short  *restrict lFaceR  ,short *restrict lFaceL
               ,DOUBLE *restrict u       ,DOUBLE *restrict gradU 
               ,short const nFace        ,short const ndm   
               ,short const ndf          ,INT const nel)
{
  DOUBLE v,dPviz,coefDif,tmp;
  DOUBLE uf[MAX_NUM_FACE*MAX_NDF],uC[MAX_NDF];
  DOUBLE lModKsi,alpha,alphaMenosUm,invVol,xx[3];
  INT vizNel;
  short idCell = nFace,nCarg,type;
  short i,j,k;

	invVol = 1.e0/volume[idCell];
/*...*/

/*... */
  if(ndf == 1){
    uC[0] = u[idCell];
    for(i=0;i<nFace;i++){
      vizNel = lViz[i];
/*... dominio*/
      if(vizNel > -1){
        lModKsi    = mKsi[i];
/*...*/
        dPviz = 0.e0;
        for(j=0;j<ndm;j++){
          v          = MAT2D(i,j,vSkew,ndm) + MAT2D(i,j,xmcc,ndm);
          dPviz     += v*v;
        }
/*...................................................................*/

/*...*/
        dPviz        = sqrt(dPviz);
        alpha        = dPviz/lModKsi;
        alphaMenosUm = 1.e0 - alpha;
/*...................................................................*/

/*...*/
        uf[i] = alphaMenosUm*uC[0] + alpha*u[i];
/*...................................................................*/
      }
/*...................................................................*/

/*... contorno*/
      else{
/*... temperatura prescrita na face(extrapolacao linear)*/
        if(lFaceR[i]){
          nCarg= lFaceL[i]-1;
          type = loads[nCarg].type;
/*... valor prescrito*/
          if( type == DIRICHLETBC || type == INLET){
            uf[i] = loads[nCarg].par[0];
          }
/*...................................................................*/

/*... fluxo prescrito*/
          else if (type == NEUMANNBC ){
            coefDif = lProp[COEFDIF];
            if(coefDif != 0.e0  )
              uf[i] = uC[0] - (loads[nCarg].par[0]/coefDif)*lDcca[i];
          }
/*...................................................................*/

/*... condicao de robin*/
          else if(type == ROBINBC){
            uf[i] =uC[0];
            for(j=0;j<ndm;j++)
              uf[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
          }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
          else if (type == OUTLET){
            uf[i] = uC[0];
          }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
          else if(type == SINBC){
            if(ndm == 2){
              xx[0] = MAT2D(i,0,xm,ndm);
              xx[1] = MAT2D(i,1,xm,ndm);
              xx[2] = 0.e0;
            }
            else{              
              xx[0] = MAT2D(i,0,xm,ndm);
              xx[1] = MAT2D(i,1,xm,ndm);
              xx[2] = MAT2D(i,2,xm,ndm);
             }
            loadSenProd(&uf[i],loads[nCarg].par,xx);
          }
/*...................................................................*/
        } 
/*...................................................................*/

/*... fluxo nulo*/
        else{
          uf[i] =uC[0];
          for(j=0;j<ndm;j++)
            uf[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
    for(k=0;k<ndm;k++){
      tmp = 0.e0;
      for(i=0;i<nFace;i++){
        tmp += uf[i]*fArea[i]*MAT2D(i,k,normal,ndm);; 
      gradU[k] = tmp*invVol; 
      }
    }
/*...................................................................*/
  } 
/*...................................................................*/

/*... graus de liberdade maior que 1*/
  else{
/*...*/
    for(k=0;k<ndf;k++)
      uC[k] = MAT2D(idCell,k,u,ndf);
/*...................................................................*/
    for(i=0;i<nFace;i++){
      vizNel = lViz[i];
/*... dominio*/
      if(vizNel > -1){
        lModKsi    = mKsi[i];
/*...*/
        dPviz = 0.e0;
        for(j=0;j<ndm;j++){
          v          = MAT2D(i,j,vSkew,ndm) + MAT2D(i,j,xmcc,ndm);
          dPviz     += v*v;
        }
/*...................................................................*/

/*...*/
        dPviz        = sqrt(dPviz);
        alpha        = dPviz/lModKsi;
        alphaMenosUm = 1.e0 - alpha;
/*...................................................................*/

/*...*/
        for(k=0;k<ndf;k++)
          MAT2D(i,k,uf,ndf) = alphaMenosUm*uC[k] 
                            + alpha*MAT2D(i,k,u,ndf);
/*...................................................................*/
      }
/*...................................................................*/

/*... contorno*/
      else{
/*... temperatura prescrita na face(extrapolacao linear)*/
        if(lFaceR[i] > 0){
          nCarg= lFaceL[i]-1;
          type = loads[nCarg].type;
/*... valor prescrito*/
          if( type == DIRICHLETBC || type == INLET || type == MOVEWALL){
            for(k=0;k<ndf;k++)
              MAT2D(i,k,uf,ndf) = loads[nCarg].par[k];
            
          }
/*...................................................................*/

/*... fluxo prescrito*/
          else if (type == NEUMANNBC ){
            coefDif = lProp[COEFDIF];
            if(coefDif != 0.e0  )
              for(k=0;k<ndf;k++)
                MAT2D(i,k,uf,ndf) = uC[k] 
                             - (loads[nCarg].par[k]/coefDif)*lDcca[i];
          }
/*...................................................................*/

/*... condicao de robin*/
          else if(type == ROBINBC){
            ERRO_GERAL(__FILE__,__func__,__LINE__
            ,"Condicao de robin na implementada");
          }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
          else if (type == OUTLET){
            for(k=0;k<ndf;k++)
              MAT2D(i,k,uf,ndf) = uC[k];
          }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
          else if(type == SINBC){
            ERRO_GERAL(__FILE__,__func__,__LINE__
            ,"Condicao SINBC robin nao implementada");
          }
/*...................................................................*/
        } 
/*...................................................................*/

/*... parede static impermeavel*/
        else if(lFaceR[i]==STATICWALL){
          for(k=0;k<ndf;k++)
            MAT2D(i,k,uf,ndf) = 0.e0;
        }
/*...................................................................*/

/*... fluxo nulo*/
        else{
          for(k=0;k<ndf;k++){
            MAT2D(i,k,uf,ndf) = uC[k];
            for(j=0;j<ndm;j++)
              MAT2D(i,k,uf,ndf) += 
              MAT2D(k,j,gradU,ndm)*MAT2D(i,j,xmcc,ndm);
          }
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*...*/    
    for(k=0;k<ndf;k++){
      for(j=0;j<ndm;j++){
        tmp = 0.e0;
        for(i=0;i<nFace;i++)
          tmp += MAT2D(i,k,uf,ndf)*fArea[i]*MAT2D(i,j,normal,ndm); 
        MAT2D(k,j,gradU,ndm) = tmp*invVol; 
      }
    }
/*...................................................................*/
  }
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * GRREENGAUSSNODE: reconstrucao de gradiente green-gauss linear por * 
 * celula                                                            *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lViz      -> viznhos da celula central                            * 
 * fArea     -> area do face                                         * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> area da celula central                               * 
 * lFaceR    -> restricoes por elmento                               * 
 * lFaceS    -> carga por elemento                                   * 
 * u         -> solucao conhecida nodal                              * 
 * gradU     -> gradiente rescontruido da solucao conhecida por      *
 *              celula                                               * 
 * isNod     -> numeracao dos nos por aresta                         * 
 * nFace     -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * ndf       -> grauss de liberdade                                  * 
 * ty        -> tipo geometrico da celula                            * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * gradU     -> gradiente calculado                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void greenGaussNode(INT *restrict lViz   ,DOUBLE *restrict fArea
               ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
               ,DOUBLE *restrict u       ,DOUBLE *restrict gradU 
               ,short *restrict isNod       
               ,short const nFace        ,short const ndm   
               ,short const ndf          ,short const ty)
{
  DOUBLE lNormal[MAX_NDM],aux1;
  DOUBLE uf[MAX_NDF];
  DOUBLE invVol;
  INT no[4];
  short idCell = nFace;
  short i,j,k;

  invVol = 1.e0/volume[idCell];
  
  for(i=0;i<ndf*ndm;i++)
    gradU[i]  = 0.e0;
  
/*...*/
  for(i=0;i<nFace;i++){
    
    for(j=0;j<ndm;j++)
      lNormal[j] = MAT2D(i,j,normal,ndm);
/*... triangulos e quadrilateros*/    
    if(ty == TRIACELL || ty == QUADCELL){
      no[0] = MAT2D(i,0,isNod,2);
      no[1] = MAT2D(i,1,isNod,2);
      for(j=0;j<ndf;j++)
        uf[j] = 0.5e0*(MAT2D(no[0],j,u,ndf) + MAT2D(no[1],j,u,ndf));
    }
/*...................................................................*/

/*... tetraedro*/    
    else if(ty == TETRCELL){
      no[0] = MAT2D(i,0,isNod,3);
      no[1] = MAT2D(i,1,isNod,3);
      no[2] = MAT2D(i,2,isNod,3);
      for(j=0;j<ndf;j++)
        uf[j] = oneDivTree*(MAT2D(no[0],j,u,ndf) 
                           +MAT2D(no[1],j,u,ndf)
                           +MAT2D(no[2],j,u,ndf));
    }
/*...................................................................*/

/*... hexaedro*/    
    else if(ty == HEXACELL){
      no[0] = MAT2D(i,0,isNod,4);
      no[1] = MAT2D(i,1,isNod,4);
      no[2] = MAT2D(i,2,isNod,4);
      no[3] = MAT2D(i,3,isNod,4);
      for(j=0;j<ndf;j++)
        uf[j] = 0.25e0*(MAT2D(no[0],j,u,ndf) 
                       +MAT2D(no[1],j,u,ndf)
                       +MAT2D(no[2],j,u,ndf)
                       +MAT2D(no[3],j,u,ndf));
    }
/*...................................................................*/

/*...*/
    for(j=0;j<ndf;j++){
      aux1 = uf[j]*fArea[i]; 
      for(k=0;k<ndm;k++)
        MAT2D(j,k,gradU,ndm) += aux1*lNormal[k]; 
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  for(i=0;i<ndf;i++)
    for(j=0;j<ndm;j++)
      MAT2D(i,j,gradU,ndm) *= invVol; 
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 09/07/2016                                   * 
 *-------------------------------------------------------------------* 
 * LEASTSQUARE : calcula o gradiente por minimos quadrados           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * loads     -> definicoes de cargas                                 * 
 * lLsquare  -> matriz para a reconstrucao least Square              * 
 * xm        -> pontos medios das faces das celulas                  * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * lProp     -> propriedades dos material                            *
 * lDcca     -> menor distancia do centroide a faces desta celula    * 
 * lViz      -> viznhos da celula central                            * 
 * u         -> solucao conhecida                                    * 
 * gradU     -> gradiente rescontruido da solucao conhecida          * 
 * lFaceR    -> restricoes por elmento                               * 
 * lFaceS    -> carga por elemento                                   * 
 * nFace     -> numero vizinhos por celula maximo da malha           * 
 * ndf       -> grauss de liberdade                                  * 
 * ndm       -> numero de dimensoes                                  * 
 * nel       -> numero do elemento                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lSquare   -> matriz para a reconstrucao least Square              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void  leastSquare(Loads *loads
                 ,DOUBLE *restrict lLsquare,INT *restrict lViz
                 ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc 
                 ,DOUBLE *restrict lProp   ,DOUBLE *restrict lDcca 
                 ,DOUBLE *restrict u       ,DOUBLE *restrict gradU
                 ,short  *restrict lFaceR  ,short *restrict lFaceL
                 ,short const nFace        ,short const ndf
                 ,short const ndm          ,INT const nel){

  DOUBLE du[MAX_NUM_FACE*MAX_NDF],uC[MAX_NDF];
  DOUBLE uT[MAX_NDF],xx[3];
  DOUBLE tmp,coefDif;
  INT vizNel;
  short idCell = nFace,nCarg,type;
  short i,j,k,l;

  for(l=0;l<2;l++){
/*... um grau de liberdade*/  
    if(ndf == 1){
      uC[0] = u[idCell];  
      for(i=0;i<nFace;i++){
        vizNel = lViz[i];
/*... dominio*/
        if(vizNel > -1)
          du[i] = u[i] - uC[0];
/*... contorno*/
        else{
/*... temperatura prescrita na face(extrapolacao linear)*/
          if(lFaceR[i]){
            nCarg=lFaceL[i]-1;
            type = loads[nCarg].type;
/*... valor prescrito*/
            if( type == DIRICHLETBC || type == INLET){
              du[i] = loads[nCarg].par[0] - uC[0];
             }
/*...................................................................*/

/*... fluxo prescrito*/
            else if (type == NEUMANNBC){
              coefDif = lProp[COEFDIF];
              if(coefDif != 0.e0  )
                du[i] = (loads[nCarg].par[0]/coefDif)*lDcca[i];
            }
/*...................................................................*/

/*... condicao de robin*/
            else if( type == ROBINBC ){
              du[i] =0.e0;
              for(j=0;j<ndm;j++)
                du[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
            }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
            else if (type == OUTLET){
              du[i] = 0.e0;
            }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
            else if( type == SINBC ){
              if(ndm == 2){
                xx[0] = MAT2D(i,0,xm,ndm);
                xx[1] = MAT2D(i,1,xm,ndm);
                xx[2] = 0.e0;
              }
              else{              
                xx[0] = MAT2D(i,0,xm,ndm);
                xx[1] = MAT2D(i,1,xm,ndm);
                xx[2] = MAT2D(i,2,xm,ndm);
              }
              loadSenProd(uT,loads[nCarg].par,xx);
              du[i] = uT[0] - uC[0];
            }
/*...................................................................*/
          } 
/*...................................................................*/

/*... fluxo nulo*/
          else{
            du[i] =0.e0;
            for(j=0;j<ndm;j++)
              du[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/

/*... gradU = G du*/
      for(i=0;i<ndm;i++){
        tmp = 0.e0;
        for(j=0;j<nFace;j++)
          tmp += MAT2D(i,j,lLsquare,nFace)*du[j]; 
        gradU[i] = tmp;
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... graus de liberdade maior que 1 */
    else{
      for(k=0;k<ndf;k++)
        uC[k] = MAT2D(idCell,k,u,ndf);
/*... loop nas faces*/
      for(i=0;i<nFace;i++){
        vizNel = lViz[i];
/*... dominio*/
        if(vizNel > -1)
          for(k=0;k<ndf;k++)
            MAT2D(i,k,du,ndf) =  MAT2D(i,k,u,ndf) - uC[k];
/*... contorno*/
        else{
/*... potencial prescrito na face(extrapolacao linear)*/
          if(lFaceR[i] > 0){
            nCarg= lFaceL[i]-1;
            type = loads[nCarg].type;
/*... valor prescrito*/
            if( type == DIRICHLETBC || type == INLET 
             || type == MOVEWALL){
              for(k=0;k<ndf;k++)
                MAT2D(i,k,du,ndf) = loads[nCarg].par[k] - uC[k];
            }
/*...................................................................*/

/*... condicao de robin*/
            else if(type == ROBINBC){
              ERRO_GERAL(__FILE__,__func__,__LINE__
              ,"Condicao de robin na implementada");
            }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
            else if(type == SINBC){
              ERRO_GERAL(__FILE__,__func__,__LINE__
              ,"Condicao SINBC robin nao implementada");
            }
/*...................................................................*/
         
/*... derivada nula (condicao localmente parabolica saida)*/
            else if (type == OUTLET)
              for(k=0;k<ndf;k++)
                MAT2D(i,k,du,ndf) = 0.e0;                          
/*...................................................................*/
          } 
/*...................................................................*/

/*... parede static impermeavel*/
          else if(lFaceR[i]==STATICWALL){
            for(k=0;k<ndf;k++)
              MAT2D(i,k,du,ndf) = - uC[k];
           }
/*...................................................................*/

/*... fluxo nulo*/
          else{
            for(k=0;k<ndf;k++){
              MAT2D(i,k,du,ndf) = 0.e0;
              for(j=0;j<ndm;j++)
                MAT2D(i,k,du,ndf) += 
                MAT2D(k,j,gradU,ndm)*MAT2D(i,j,xmcc,ndm);
            }
          }
/*...................................................................*/
        } 
/*...................................................................*/
      }
/*...................................................................*/

/*... gradU = G du*/
      for(k=0;k<ndf;k++)
        for(i=0;i<ndm;i++){
          tmp = 0.e0;
          for(j=0;j<nFace;j++)
            tmp += MAT2D(i,j,lLsquare,nFace)*MAT2D(j,k,du,ndf); 
          MAT2D(k,i,gradU,ndm) = tmp;
        }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
} 
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 09/07/2016                                   * 
 *-------------------------------------------------------------------* 
 * LEASTSQUAREQR: calcula o gradiente por minimos quadrados           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * loads     -> definicoes de cargas                                 * 
 * lLsQt     -> matriz para a reconstrucao least Square              * 
 * lLsR      -> fatoracao R (RCLSQUAREQR)                            * 
 * lProp     -> propriedades dos material                            *
 * lDcca     -> menor distancia do centroide a faces desta celula    * 
 * xm        -> pontos medios das faces das celulas                  * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * lViz      -> viznhos da celula central                            * 
 * u         -> solucao conhecida                                    * 
 * gradU     -> gradiente rescontruido da solucao conhecida          * 
 * lFaceR    -> restricoes por elmento                               * 
 * lFaceS    -> carga por elemento                                   * 
 * nFace     -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * ndf       -> grauss de liberdade                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lSquare   -> matriz para a reconstrucao least Square              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void  leastSquareQR(Loads *loads
                   ,DOUBLE *restrict lLsQt   ,DOUBLE *restrict lLsR
                   ,DOUBLE *restrict lProp   ,DOUBLE *restrict lDcca 
                   ,INT *restrict lViz       
                   ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
                   ,DOUBLE *restrict u       ,DOUBLE *restrict gradU
                   ,short  *restrict lFaceR  ,short *restrict lFaceL
                   ,short const nFace        ,short const ndf
                   ,short const ndm){

  DOUBLE du[MAX_NUM_FACE*MAX_NDF],uC[MAX_NDF],coefDif;
  DOUBLE uT[MAX_NDF],xx[3];
  DOUBLE b[MAX_NUM_FACE*MAX_NDF],b0,b1,b2;
  INT vizNel;
  short idCell = nFace,type;
  short i,j,k,l,nCarg;

  for(l=0;l<1;l++){
/*... um grau de liberdade*/  
    if(ndf == 1){
      uC[0] = u[idCell];  
      for(i=0;i<nFace;i++){
        vizNel = lViz[i];
/*... dominio*/
        if(vizNel > -1)
          du[i] = u[i] - uC[0];
/*... contorno*/
        else{
/*... temperatura prescrita na face(extrapolacao linear)*/
          if(lFaceR[i]){
            nCarg=lFaceL[i]-1;
            type = loads[nCarg].type;
/*... valor prescrito*/
            if( type == DIRICHLETBC || type == INLET ){
              du[i] = loads[nCarg].par[0] - uC[0];
             }
/*...................................................................*/

/*... fluxo prescrito*/
            else if (type == NEUMANNBC){
              coefDif = lProp[COEFDIF];
              if(coefDif != 0.e0  )
                du[i] = (loads[nCarg].par[0]/coefDif)*lDcca[i];
            }
/*...................................................................*/

/*... condicao de robin*/
            else if( type == ROBINBC ){
              du[i] =0.e0;
              for(j=0;j<ndm;j++)
                du[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
            }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
            else if (type == OUTLET){
              du[i] = 0.e0;
            }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
            else if( type == SINBC ){
              if(ndm == 2){
                xx[0] = MAT2D(i,0,xm,ndm);
                xx[1] = MAT2D(i,1,xm,ndm);
                xx[2] = 0.e0;
              }
              else{              
                xx[0] = MAT2D(i,0,xm,ndm);
                xx[1] = MAT2D(i,1,xm,ndm);
                xx[2] = MAT2D(i,2,xm,ndm);
              }
              loadSenProd(uT,loads[nCarg].par,xx);
              du[i] = uT[0] - uC[0];
            }
/*...................................................................*/
          } 
/*...................................................................*/

/*... fluxo nulo*/
          else{
            du[i] =0.e0;
            for(j=0;j<ndm;j++)
              du[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/

/*... b = Qtdu*/
      for(i=0;i<ndm;i++){
        b[i]      = 0.e0;
        for(j=0;j<nFace;j++){
          b[i] += MAT2D(i,j,lLsQt,nFace)*du[j]; 
        } 
      }
/*...................................................................*/
    
/*... R grad = Qtdu*/
      if(ndm == 2){
/*... retrosubstituicao*/
        gradU[1] = b[1]/lLsR[2];
        gradU[0] = (b[0] - lLsR[1]*gradU[1])/lLsR[0];
      }
      else if(ndm == 3){
/*... retrosubstituicao*/
        gradU[2] = b[2]/lLsR[5];
        gradU[1] = (b[1] - lLsR[4]*gradU[2])/lLsR[3];
        gradU[0] = (b[0] 
                 - lLsR[1]*gradU[1] 
                 - lLsR[2]*gradU[2])/lLsR[0];
      }
/*...................................................................*/
    }
/*...................................................................*/
  
/*... graus de liberdade maior que 1 */
    else{
      for(k=0;k<ndf;k++)
        uC[k] = MAT2D(idCell,k,u,ndf);
/*... loop nas faces*/
      for(i=0;i<nFace;i++){
        vizNel = lViz[i];
/*... dominio*/
        if(vizNel > -1)
          for(k=0;k<ndf;k++)
            MAT2D(i,k,du,ndf) =  MAT2D(i,k,u,ndf) - uC[k];
/*... contorno*/
        else{
/*... potencial prescrito na face(extrapolacao linear)*/
          if(lFaceR[i] > 0){
            nCarg= lFaceL[i]-1;
            type = loads[nCarg].type;
/*... valor prescrito*/
            if( type == DIRICHLETBC || type == INLET 
             || type == MOVEWALL){
              for(k=0;k<ndf;k++)
                MAT2D(i,k,du,ndf) = loads[nCarg].par[k] - uC[k];
            }
/*...................................................................*/

/*... condicao de robin*/
            else if(type == ROBINBC){
              ERRO_GERAL(__FILE__,__func__,__LINE__
              ,"Condicao de robin na implementada");
            }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
            else if(type == SINBC){
              ERRO_GERAL(__FILE__,__func__,__LINE__
              ,"Condicao SINBC robin nao implementada");
            }
/*...................................................................*/
         
/*... derivada nula (condicao localmente parabolica saida)*/
            else if (type == OUTLET)
              for(k=0;k<ndf;k++)
                MAT2D(i,k,du,ndf) = 0.e0;                          
/*...................................................................*/
          } 
/*...................................................................*/

/*... parede static impermeavel*/
          else if(lFaceR[i]==STATICWALL){
            for(k=0;k<ndf;k++)
              MAT2D(i,k,du,ndf) = - uC[k];
           }
/*...................................................................*/

/*... fluxo nulo*/
          else{
            for(k=0;k<ndf;k++){
              MAT2D(i,k,du,ndf) = 0.e0;
              for(j=0;j<ndm;j++)
                MAT2D(i,k,du,ndf) += 
                MAT2D(k,j,gradU,ndm)*MAT2D(i,j,xmcc,ndm);
            }
          }
/*...................................................................*/
        } 
/*...................................................................*/
      }
/*... b = Qtdu*/
      for(k=0;k<ndf;k++)
        for(i=0;i<ndm;i++){
          MAT2D(k,i,b,ndf) = 0.e0;
          for(j=0;j<nFace;j++)
            MAT2D(k,i,b,ndf) += MAT2D(i,j,lLsQt,nFace)*MAT2D(j,k,du,ndf); 
        }
/*...................................................................*/
    
/*... R grad = Qtdu*/
      for(k=0;k<ndf;k++){
        if(ndm == 2){
          b0 =  MAT2D(k,0,b,ndf);
          b1 =  MAT2D(k,1,b,ndf);
/*... retrosubstituicao*/
          MAT2D(k,1,gradU,ndf) = b1/lLsR[2];
          MAT2D(k,0,gradU,ndf) = 
                           (b0  -lLsR[1]*MAT2D(k,1,gradU,ndf))/lLsR[0];
        }
        else if(ndm == 3){
          b0 =  MAT2D(k,0,b,ndf);
          b1 =  MAT2D(k,1,b,ndf);
          b2 =  MAT2D(k,2,b,ndf);
/*... retrosubstituicao*/
          MAT2D(k,2,gradU,ndf) = b2/lLsR[5];
          MAT2D(k,1,gradU,ndf)= (b1 
                   - lLsR[4]* MAT2D(k,2,gradU,ndf))/lLsR[3];
          MAT2D(k,0,gradU,ndf) = (b0 
                   - lLsR[1]* MAT2D(k,1,gradU,ndf) 
                   - lLsR[2]* MAT2D(k,2,gradU,ndf))/lLsR[0];
        }
/*...................................................................*/
      }
/*...................................................................*/
    } 
/*...................................................................*/
  }
/*...................................................................*/
} 


/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * LEASTSQUAREMATRIX : calcula a matriz de minimos quadrados para a  * 
 * rescontrucao de gradiente                                         *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * lKsi      -> vetores que unem centroide da celula central aos     *
 * mKsi      -> modulo do vetor ksi                                  * 
 * lSquare   -> nao definido                                         * 
 * lSquareR  -> nao definido                                         * 
 * type      -> tecnica de resolucao                                 * 
 * lnFace    -> numero vizinhos da celula                            * 
 * ndm       -> numero de dimensoes                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lSquare   -> matriz para a reconstrucao least Square              * 
 * lSquareR  -> fatoracao R (RCLSQUAREQR)                            * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 void leastSquareMatrix(DOUBLE *restrict lKsi    ,DOUBLE *restrict lmKsi
                  ,DOUBLE *restrict lLsquare,DOUBLE *restrict lLsquareR
                  ,short const type        
                  ,short const lnFace       ,short const ndm){

  DOUBLE dx[MAX_NUM_FACE*MAX_NDM],w[MAX_NUM_FACE*MAX_NDM];
  DOUBLE aTwA[2*MAX_NDM],inv[2*MAX_NDM];
  DOUBLE r[MAX_NDM*MAX_NDM];
  DOUBLE aux1,detA;
  short i,j,k,nf;
 
  if(ndm == 2){
    aTwA[0] = 0.e0;
    aTwA[1] = 0.e0;
    aTwA[2] = 0.e0;
  }
  else{
    aTwA[0] = 0.e0;
    aTwA[1] = 0.e0;
    aTwA[2] = 0.e0;
    aTwA[3] = 0.e0;
    aTwA[4] = 0.e0;
    aTwA[5] = 0.e0;
  }
  
  for(i=0;i<2*MAX_NDM;i++){
    inv[i] = 0.e0;
    r[i]   = 0.e0;
  }
   

/*...G=inv(aTwA)wAt */
  if(type == RCLSQUARE){
    for(nf=0;nf<lnFace;nf++){
      w[nf]  = 1.e0/lmKsi[nf];
      w[nf] *= w[nf];
/*... dimensao 2*/
      if(ndm == 2){
        MAT2D(nf,0,dx,ndm) = MAT2D(nf,0,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,1,dx,ndm) = MAT2D(nf,1,lKsi,ndm)*lmKsi[nf];
/*atTwa(0,0)*/
        aTwA[0] +=  w[nf]*MAT2D(nf,1,dx,ndm)*MAT2D(nf,1,dx,ndm);
/*atTwa(0,1)*/
        aTwA[1] +=  w[nf]*MAT2D(nf,0,dx,ndm)*MAT2D(nf,1,dx,ndm);
/*atTwa(1,1)*/
        aTwA[2] +=  w[nf]*MAT2D(nf,0,dx,ndm)*MAT2D(nf,0,dx,ndm);
      }
/*...................................................................*/

/*... dimensao 3*/
      else if(ndm == 3){
        MAT2D(nf,0,dx,ndm) = MAT2D(nf,0,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,1,dx,ndm) = MAT2D(nf,1,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,2,dx,ndm) = MAT2D(nf,2,lKsi,ndm)*lmKsi[nf];
/*atTwa(0,0)*/
        aTwA[0] +=  w[nf]*MAT2D(nf,0,dx,ndm)*MAT2D(nf,0,dx,ndm);
/*atTwa(0,1)*/
        aTwA[1] +=  w[nf]*MAT2D(nf,0,dx,ndm)*MAT2D(nf,1,dx,ndm);
/*atTwa(0,2)*/
        aTwA[2] +=  w[nf]*MAT2D(nf,0,dx,ndm)*MAT2D(nf,2,dx,ndm);
/*atTwa(1,1)*/
        aTwA[3] +=  w[nf]*MAT2D(nf,1,dx,ndm)*MAT2D(nf,1,dx,ndm);
/*atTwa(1,2)*/
        aTwA[4] +=  w[nf]*MAT2D(nf,1,dx,ndm)*MAT2D(nf,2,dx,ndm);
/*atTwa(2,2)*/
        aTwA[5] +=  w[nf]*MAT2D(nf,2,dx,ndm)*MAT2D(nf,2,dx,ndm);
      }
/*...................................................................*/
    }

/*... dimensao 2*/
    if(ndm == 2){
      detA   = 1.e0/(aTwA[2]*aTwA[0] - aTwA[1]*aTwA[1]);
      inv[0] =  aTwA[2]*detA;
      inv[1] = -aTwA[1]*detA;
      inv[2] =  aTwA[0]*detA;
    }
/*...................................................................*/
 
/*... dimensao 3*/
    else if(ndm == 3){
      detA = aTwA[0]*aTwA[3]*aTwA[5] + 2.e0*aTwA[1]*aTwA[4]*aTwA[2] 
             -( aTwA[2]*aTwA[2]*aTwA[3] 
              + aTwA[4]*aTwA[4]*aTwA[0]    
              + aTwA[1]*aTwA[1]*aTwA[5]); 
      detA = 1.e0/detA;
/*... inversao com cofatores*/
      inv[0] = (aTwA[3]*aTwA[5] - aTwA[4]*aTwA[4])*detA;
      inv[1] =-(aTwA[1]*aTwA[5] - aTwA[2]*aTwA[4])*detA;
      inv[2] = (aTwA[1]*aTwA[4] - aTwA[2]*aTwA[3])*detA;

      inv[3] = (aTwA[0]*aTwA[5] - aTwA[2]*aTwA[2])*detA;
      inv[4] =-(aTwA[0]*aTwA[4] - aTwA[1]*aTwA[2])*detA;

      inv[5] = (aTwA[0]*aTwA[3] - aTwA[1]*aTwA[1])*detA;
    }

    for(nf=0;nf<lnFace;nf++){
/*... dimensao 2*/
      if(ndm == 2){
        aux1 = w[nf];
        MAT2D(0,nf,lLsquare,lnFace) = aux1*(MAT2D(nf,0,dx,ndm)*inv[0] 
                                          + MAT2D(nf,1,dx,ndm)*inv[1]);
        MAT2D(1,nf,lLsquare,lnFace) = aux1*(MAT2D(nf,0,dx,ndm)*inv[1] 
                                          + MAT2D(nf,1,dx,ndm)*inv[2]);
      }
/*...................................................................*/

/*... dimensao 3*/

      else if(ndm == 3){
        aux1 = w[nf];
        MAT2D(0,nf,lLsquare,lnFace) = aux1*(MAT2D(nf,0,dx,ndm)*inv[0] 
                                        + MAT2D(nf,1,dx,ndm)*inv[1]
                                        + MAT2D(nf,2,dx,ndm)*inv[2]);
        MAT2D(1,nf,lLsquare,lnFace) = aux1*(MAT2D(nf,0,dx,ndm)*inv[1] 
                                        + MAT2D(nf,1,dx,ndm)*inv[3]
                                        + MAT2D(nf,2,dx,ndm)*inv[4]);
        MAT2D(2,nf,lLsquare,lnFace) = aux1*(MAT2D(nf,0,dx,ndm)*inv[2] 
                                        + MAT2D(nf,1,dx,ndm)*inv[4]
                                        + MAT2D(nf,2,dx,ndm)*inv[5]);
      }
/*...................................................................*/
    }
  }
/*...................................................................*/

/*... fatoracao QR-MGS*/
  else if(type == RCLSQUAREQR){
    for(nf=0;nf<lnFace;nf++){
      w[nf]  = 1.e0/lmKsi[nf];
      w[nf] *= w[nf];
/*... dimensao 2*/
      if(ndm == 2){
        MAT2D(nf,0,dx,ndm) = MAT2D(nf,0,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,1,dx,ndm) = MAT2D(nf,1,lKsi,ndm)*lmKsi[nf];
      }
/*...................................................................*/

/*... dimensao 3*/
      else if(ndm == 3){
        MAT2D(nf,0,dx,ndm) = w[nf]*MAT2D(nf,0,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,1,dx,ndm) = w[nf]*MAT2D(nf,1,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,2,dx,ndm) = w[nf]*MAT2D(nf,2,lKsi,ndm)*lmKsi[nf];
      }
/*...................................................................*/
    }

/*... MGS*/
    for(k=0;k<ndm;k++){
      aux1 = 0.e0;
      for(i=0;i<lnFace;i++)
        aux1 += MAT2D(i,k,dx,ndm)* MAT2D(i,k,dx,ndm);
      
      MAT2D(k,k,r,ndm) = sqrt(aux1);
      
      for(i=0;i<lnFace;i++)
        MAT2D(i,k,dx,ndm) /= MAT2D(k,k,r,ndm); 
      
      for(j=k+1;j<ndm;j++){
        aux1 = 0.e0;
        for(i=0;i<lnFace;i++)
          aux1 += MAT2D(i,k,dx,ndm)*MAT2D(i,j,dx,ndm);
        
        MAT2D(k,j,r,ndm) = aux1; 
        
        for(i=0;i<lnFace;i++)
          MAT2D(i,j,dx,ndm) -=  MAT2D(i,k,dx,ndm)*MAT2D(k,j,r,ndm); 
      }
    }
/*...................................................................*/

/*... dimensao 2*/
    if(ndm == 2){
        lLsquareR[0] = MAT2D(0,0,r,ndm);
        lLsquareR[1] = MAT2D(0,1,r,ndm);
        lLsquareR[2] = MAT2D(1,1,r,ndm);
    }
/*...................................................................*/

/*... R dimensao 3*/
    else if(ndm == 3){
        lLsquareR[0] = MAT2D(0,0,r,ndm);
        lLsquareR[1] = MAT2D(0,1,r,ndm);
        lLsquareR[2] = MAT2D(0,2,r,ndm);
        lLsquareR[3] = MAT2D(1,1,r,ndm);
        lLsquareR[4] = MAT2D(1,2,r,ndm);
        lLsquareR[5] = MAT2D(2,2,r,ndm);
    }
/*...................................................................*/

    for(nf=0;nf<lnFace;nf++){
/*... dimensao 2*/
      if(ndm == 2){
        MAT2D(0,nf,lLsquare,lnFace) = MAT2D(nf,0,dx,ndm);
        MAT2D(1,nf,lLsquare,lnFace) = MAT2D(nf,1,dx,ndm);
      }
/*...................................................................*/

/*... Qt - dimensao 3*/
      else if(ndm == 3){
        MAT2D(0,nf,lLsquare,lnFace) = w[nf]*MAT2D(nf,0,dx,ndm);
        MAT2D(1,nf,lLsquare,lnFace) = w[nf]*MAT2D(nf,1,dx,ndm);
        MAT2D(2,nf,lLsquare,lnFace) = w[nf]*MAT2D(nf,2,dx,ndm);
      }
/*...................................................................*/
    }
/*...................................................................*/
  } 
/*...................................................................*/
   
  }
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * SN: numera dos nos das faces das celulas                          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * s   -> indefinido                                                 * 
 * ty  -> numero do tipo geometrico da celula                        * 
 * nel -> numero da celula                                           * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * s   -> numeracao dos nos por aresta                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
short sn(short *s,short ty,INT nel)
{
    
  long aux;
  short i;
  short isNodTria[6] = {0,1
                       ,1,2
                       ,2,0};
  short isNodQuad[8] = {0,1
                       ,1,2
                       ,2,3
                       ,3,0};
  short isNodTetr[12] = {1,2,3
                        ,0,3,2
                        ,0,1,3
                        ,0,2,1};
  short isNodHex[24] = {0,3,2,1
                       ,4,5,6,7
                       ,0,1,5,4
                       ,1,2,6,5
                       ,2,3,7,6
                       ,3,0,4,7};

  switch(ty){

/* ... triangulo*/  
    case TRIACELL:
      for(i=0;i< 6;i++)
        s[i] = isNodTria[i];
      return 2;
    break;
/*...................................................................*/

/* ... quadrilateros*/  
    case QUADCELL:
      for(i=0;i< 8;i++)
        s[i] = isNodQuad[i];
      return 2;
    break;
/*...................................................................*/

/* ... tetraedros*/  
    case TETRCELL:
      for(i=0;i< 12;i++)
        s[i] = isNodTetr[i];
      return 3;
    break;
/*...................................................................*/

/* ... hexaedros*/  
    case HEXACELL:
      for(i=0;i<24;i++)
        s[i] = isNodHex[i];
      return 4;
    break;
/*...................................................................*/

/*...*/    
    default: 
      aux = (long) nel;
      printf("Erro: tipo de celula nao exitentes. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
               ,aux,__FILE__,__func__);
      exit(EXIT_FAILURE);
    break; 
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * VOLUMECELL : volume(area) da celula                               * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * eta -> vetores paralelos as faces da celula                       * 
 * ty  -> tipo geometrico da celula                                  * 
 * ndm -> numero de dimensao                                         * 
 * nel -> numero da celula                                           * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * retorna o valores da volume (area) da celula                      * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE areaCell(DOUBLE *eta,short ty,short ndm,INT nel)
{
  
  long aux;

  switch(ty){

/* ... triangulo*/  
    case TRIACELL: 
      return areaTriaCell(eta,ndm);
    break;
/*...................................................................*/

/* ... quadrilateros*/  
    case QUADCELL:
      return areaQuadCell(eta,ndm);
    break;
/*...................................................................*/

/*...*/    
    default: 
      aux = (long) nel;
			printf("Erro: tipo de celula nao exitentes. Nel = %ld.\n"
				"Arquivo fonte:  \"%s\".\n"
				"Nome da funcao: \"%s\".\n"
				, aux, __FILE__, __func__);
      exit(EXIT_FAILURE);
    break; 
/*...................................................................*/
  }
}
/*********************************************************************/ 


/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * AREATRIACELL: calculo da area da celula triangular                * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * eta -> vetores paralelos as arestas da celula                     * 
 * ndm -> numero de dimensao                                         * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * retorna o valores da area da celula                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE areaTriaCell(DOUBLE *restrict eta, short ndm)
{
  double v1[3],v2[3],v3[3],a,dot;

/*... aresta 1*/
  v1[0] = MAT2D(0,0,eta,ndm);
  v1[1] = MAT2D(0,1,eta,ndm);
  if( ndm == 2)
    v1[2] = 0.0e0;               
  else
    v1[2] = MAT2D(0,2,eta,ndm);
  

/*... aresta 2*/
  v2[0] = MAT2D(1,0,eta,ndm);
  v2[1] = MAT2D(1,1,eta,ndm);
  if( ndm == 2)
    v2[2] = 0.0e0;               
  else
    v2[2] = MAT2D(1,2,eta,ndm);

/*...*/
  prodVet(v1,v2,v3);
  dot = v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2];
  a = 0.5e0*sqrt(dot); 
  return a;
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * AREAQUADCELL: calculo da area da celula quadrangular              * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * eta -> vetores paralelos as arestas da celula                     * 
 * ndm -> numero de dimensao                                         * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * retorna o valores da area da celula                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE areaQuadCell(DOUBLE *restrict eta,short ndm)
{

  double v1[3],v2[3],c[3];
  double a,dot;

/*...*/
  v1[0] = MAT2D(0,0,eta,ndm);
  v1[1] = MAT2D(0,1,eta,ndm);
  if(ndm == 2) 
    v1[2] = 0.0e0;
  else
    v1[2] = MAT2D(0,2,eta,ndm);
    
  v2[0] = MAT2D(1,0,eta,ndm);
  v2[1] = MAT2D(1,1,eta,ndm);
  if(ndm == 2) 
    v2[2] = 0.0e0;
  else
    v2[2] = MAT2D(1,2,eta,ndm);
  
  prodVet(v1,v2,c);
  dot = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
  a = 0.5e0*sqrt(dot);
/*...................................................................*/

/*...*/
  v1[0] = MAT2D(2,0,eta,ndm);
  v1[1] = MAT2D(2,1,eta,ndm);
  if(ndm == 2) 
    v1[2] = 0.0e0;
  else
    v1[2] = MAT2D(2,2,eta,ndm);
  

  v2[0] = MAT2D(3,0,eta,ndm);
  v2[1] = MAT2D(3,1,eta,ndm);
  if(ndm == 2) 
    v2[2] = 0.0e0;
  else
    v2[2] = MAT2D(3,2,eta,ndm);
  
  prodVet(v1,v2,c);
  dot = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
  a += 0.5e0*sqrt(dot);
/*...................................................................*/


/*...*/
  v1[0] = MAT2D(0,0,eta,ndm);
  v1[1] = MAT2D(0,1,eta,ndm);
  if(ndm == 2) 
    v1[2] = 0.0e0;
  else
    v1[2] = MAT2D(0,2,eta,ndm);
  
  v2[0] = MAT2D(3,0,eta,ndm);
  v2[1] = MAT2D(3,1,eta,ndm);
  if(ndm == 2) 
    v2[2] = 0.0e0;
  else
    v2[2] = MAT2D(3,2,eta,ndm);
  
  prodVet(v1,v2,c);
  dot = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
  a += 0.5e0*sqrt(dot);
/*...................................................................*/


/*...*/
  v1[0] = MAT2D(1,0,eta,ndm);
  v1[1] = MAT2D(1,1,eta,ndm);
  if(ndm == 2) 
    v1[2] = 0.0e0;
  else
    v1[2] = MAT2D(1,2,eta,ndm);
  
  v2[0] = MAT2D(2,0,eta,ndm);
  v2[1] = MAT2D(2,1,eta,ndm);
  if(ndm == 2) 
    v2[2] = 0.0e0;
  else
    v2[2] = MAT2D(2,2,eta,ndm);
  
  prodVet(v1,v2,c);
  dot = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
  a += 0.5e0*sqrt(dot);
/*...................................................................*/
  return 0.5e0*a;
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * VECTORKM2D : calculo da distancia entre a intersecao entre a reta *
 * que une o centroide da celulas e a face compartilhada             *  
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x      -> coordenadas dos nos da celula central e seus viznhos    * 
 * xc     -> centroide das celulas                                   * 
 * xm     -> pontos medios das faces da celula cenral                * 
 * vSkew  -> indefinido                                              * 
 * mvSkew -> indefinido                                              * 
 * sn     -> numeracao dos nos por aresta                            * 
 * nFace  -> numero de face da celula                                * 
 * sn     -> numeracao dos nos por aresta                            * 
 * maxNo  -> numero de nos por celula maximo da malha                * 
 * maxViz -> numero vizinhos por celula maximo da malha              * 
 * ndm    -> numero de dimensoes                                     * 
 * nel    -> numero da celula                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * vSkew  -> vetor entre o ponto medio a intersecao que une os       * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew -> distacia entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void vectorKm2d(DOUBLE *restrict x     ,DOUBLE *restrict xc
               ,DOUBLE *restrict xm
               ,DOUBLE *restrict vSkew ,DOUBLE *restrict mvSkew
               ,short  *restrict sn    ,short const nFace
               ,short const maxViz     ,short const maxNo       
               ,short const ndm        ,INT const nel)
{
  short i,cCell=maxViz;
  double a1,b1,c1,a2,b2,c2;
  double a[2][2],f[2],xi[2],lvSkew[2];
  short no1,no2;
  double det = 0.0;
  long aux;
/*... reta 1: a1x + b1y + c1 = 0
      reta 2: a2x + b2y + c2 = 0

      | a1 b1 | |x| = |-c1| 
      | a2 b2 | |y| = |-c2| 
    
  a1 = i e b1 = cCell (celula central)
  a2 e b2 pontos da aresta da celula central 
.......................................................................*/
  for(i=0;i<nFace;i++){
/*...ya - yb*/
    a1 = MAT2D(i,1,xc,ndm)     - MAT2D(cCell,1,xc,ndm);  
/*...xb - xa*/
    b1 = MAT2D(cCell,0,xc,ndm) - MAT2D(i,0,xc,ndm);  
/*...xayb - xbya*/
    c1 = MAT2D(i    ,0,xc,ndm)*MAT2D(cCell,1,xc,ndm) 
       - MAT2D(cCell,0,xc,ndm)*MAT2D(i    ,1,xc,ndm);
/*.....................................................................*/

    no1 = MAT2D(i,0,sn,ndm);
    no2 = MAT2D(i,1,sn,ndm);

/*...ya - yb*/
    a2 = MAT3D(cCell,no1,1,x,maxNo,ndm) - MAT3D(cCell,no2,1,x,maxNo,ndm);  
/*...xb - xa*/
    b2 = MAT3D(cCell,no2,0,x,maxNo,ndm) - MAT3D(cCell,no1,0,x,maxNo,ndm);  
/*...xayb - xbya*/
    c2 = MAT3D(cCell,no1,0,x,maxNo,ndm)*MAT3D(cCell,no2,1,x,maxNo,ndm)  
       - MAT3D(cCell,no2,0,x,maxNo,ndm)*MAT3D(cCell,no1,1,x,maxNo,ndm);  
/*.....................................................................*/

/*... verifica se ha intersecao entre duas retas*/
    det = a1*b2 - a2*b1;
    if(fabs(det) < sqrt(a1*a1+a2*a2)*1.e-16){
      aux = (long) nel;
      printf("Erro: As retas nao sao concorrentes. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
               "Determinante  : \"%e\".\n"
               ,aux,__FILE__,__func__,fabs(det));
      exit(EXIT_FAILURE);
    }
/*... matriz inversa A*/
    a[0][0] = (1.0/det)*b2;
    a[0][1] =-(1.0/det)*b1;
    a[1][0] =-(1.0/det)*a2;
    a[1][1] = (1.0/det)*a1;
/*... f*/
    f[0] = -c1;
    f[1] = -c2;
/*... x = A-1f*/
    xi[0] = a[0][0]*f[0] + a[0][1]*f[1];
    xi[1] = a[1][0]*f[0] + a[1][1]*f[1];
/*...vertor vSkew*/
    lvSkew[0] = xi[0] - MAT2D(i,0,xm,ndm);
    lvSkew[1] = xi[1] - MAT2D(i,1,xm,ndm);
/*...................................................................*/
    mvSkew[i] = sqrt(lvSkew[0]*lvSkew[0] + lvSkew[1]*lvSkew[1]);
/*...*/
    MAT2D(i,0,vSkew,ndm) = lvSkew[0];
    MAT2D(i,1,vSkew,ndm) = lvSkew[1];
  }
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * VECTORKM3D : calculo da distancia entre a intersecao entre a reta *
 * que une o centroide da celulas e a face compartilhada             *  
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x      -> coordenadas dos nos da celula central e seus viznhos    * 
 * xc     -> centroide das celulas                                   * 
 * xm     -> pontos medios das faces da celula cenral                * 
 * normal -> vetores normais as faces das celulas                    * 
 * ksi    -> vetores que unem centroide da celula central aos        *
 *           vizinhos destas                                         * 
 * vSkew  -> indefinido                                              * 
 * mvSkew -> indefinido                                              * 
 * nFace  -> numero de faces da celula                               * 
 * ndm    -> numero de dimensoes                                     * 
 * maxViz -> numero vizinhos por celula maximo da malha              * 
 * nel    -> numero da celula                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void vectorKm3d(DOUBLE *restrict xc   ,DOUBLE *restrict xm
               ,DOUBLE *restrict ksi  ,DOUBLE *restrict normal
               ,DOUBLE *restrict vSkew,DOUBLE *restrict mvSkew
               ,short const nFace     ,short const ndm
               ,short const maxViz    ,INT nel)
{
  short i,cCell=maxViz;
  DOUBLE n[3],lKsi[3],xP[3],xR[3],xi[3],lvSkew[3],dot1,dot2,t;
  long aux;
/*... eq do plano :  nx*(x-xm) + ny*(y-ym) nz*(z-zm) = 0
      eq parametrica (x,y,z) = (xc,yc,zc) + t*(ksix,kisy,kisz)

  t = (nx(xm-xc) + ny(ym - yc) + nz(zm - zc))/(nx*ksix + ny*ksiy + nz*ksiz)

  t = n * dx / n*ksi  
  
.......................................................................*/
  xR[0] = MAT2D(cCell,0,xc,ndm);
  xR[1] = MAT2D(cCell,1,xc,ndm);
  xR[2] = MAT2D(cCell,2,xc,ndm);
  for(i=0;i<nFace;i++){
/*... vetor normal que define o plano e ponto*/
    n[0]  = MAT2D(i,0,normal,ndm);
    n[1]  = MAT2D(i,1,normal,ndm);
    n[2]  = MAT2D(i,2,normal,ndm);
    xP[0] = MAT2D(i,0,xm,ndm);
    xP[1] = MAT2D(i,1,xm,ndm);
    xP[2] = MAT2D(i,2,xm,ndm);
    
/*... vetor paralelo a reta*/
    lKsi[0] = MAT2D(i,0,ksi,ndm);
    lKsi[1] = MAT2D(i,1,ksi,ndm);
    lKsi[2] = MAT2D(i,2,ksi,ndm);
/*... */
    dot2 = n[0]*(xP[0]-xR[0]) +  n[1]*(xP[1]-xR[1]) + n[2]*(xP[2]-xR[2]);

/*... verifica se ha intersecao entre duas retas*/
    dot1 = n[0]*lKsi[0] + n[1]*lKsi[1] + n[2]*lKsi[2];
    if(fabs(dot1) < fabs(dot2)*1.e-16){
      aux = (long) nel;
      printf("Erro: A reta e o plano nao cruzam. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
               "Determinante  : \"%e\".\n"
               ,aux,__FILE__,__func__,fabs(dot1));
      exit(EXIT_FAILURE);
    }

/*... valor de t*/
    t = dot2/dot1;
/*...*/
    xi[0] = xR[0] + t*lKsi[0];
    xi[1] = xR[1] + t*lKsi[1];
    xi[2] = xR[2] + t*lKsi[2];

/*...vertor vSkew*/
    lvSkew[0] = xi[0] - MAT2D(i,0,xm,ndm);
    lvSkew[1] = xi[1] - MAT2D(i,1,xm,ndm);
    lvSkew[2] = xi[2] - MAT2D(i,2,xm,ndm);
/*...................................................................*/
    mvSkew[i] = 
     sqrt(lvSkew[0]*lvSkew[0]
        + lvSkew[1]*lvSkew[1] 
        + lvSkew[2]*lvSkew[2]);
/*...*/
    MAT2D(i,0,vSkew,ndm) = lvSkew[0];
    MAT2D(i,1,vSkew,ndm) = lvSkew[1];
    MAT2D(i,2,vSkew,ndm) = lvSkew[2];
  }
/*...................................................................*/
}
/*********************************************************************/
  
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * VOLUME3DGREENGAUSS : volume calculado pelo teoreima do divergente *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * xm     -> pontos medios das faces da celula cenral                * 
 * normal -> vetor normal as faces                                   * 
 * fArea  -> area das faces                                          * 
 * nFace  -> numero de face da celula                                * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE volume3DGreenGauss(DOUBLE *restrict xm,DOUBLE *restrict normal
                         ,DOUBLE *restrict fArea
                         ,short const nFace)
{
  
  DOUBLE dot = 0.e0;
  short i;
 
  for(i=0;i<nFace;i++){
    dot += fArea[i]
         * (MAT2D(i,0,xm,3)*MAT2D(i,0,normal,3) 
         +  MAT2D(i,1,xm,3)*MAT2D(i,1,normal,3) 
         +  MAT2D(i,2,xm,3)*MAT2D(i,2,normal,3)); 
  }   
  
  return dot*oneDivTree;

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 24/07/2016                                   * 
 *-------------------------------------------------------------------* 
 * pLoadSimple : condicao de contorno para velocidades               *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sP         -> termo da diagonal                                   * 
 * p          -> forca local                                         * 
 * tA         -> velocidade na face                                  * 
 * VelC       -> velocidade na centro do elemento no passo (n-1)     * 
 * n          -> normal                                              * 
 * gradVel    -> dradiente da velocidade                             * 
 * xmcc       -> vetores que unem o centroide aos pontos medios das  * 
 *              faces da celula central                              * 
 * viscosityC -> coeficiente de viscosidade dinamica                 * 
 * densityC   -> massa especifica                                    * 
 * fArea      -> area da face                                        * 
 * dcca       -> menor distancia do centroide central a face desta   *
 *               celula                                              * 
 * ld         -> definicao da carga                                  * 
 * ndm        -> numero de dimensoes                                 * 
 * fCalVel    -> true - atualizada sP e p pela equacao de velocidades* 
 *               false- nao atualizada sP e p                        * 
 * fCalPres   -> true - atualizada sP e p pela da pressao            * 
 *               false- nao atualizada sP e p                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal atualizado                           * 
 * p       -> forcas locais atualizada                               * 
 * tA      -> valor pescrito na face                                 * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void pLoadSimple(DOUBLE *restrict sP  ,DOUBLE *restrict p
          ,DOUBLE *restrict tA        ,DOUBLE *restrict velC
          ,DOUBLE *restrict n       
          ,DOUBLE *restrict gradVel   ,DOUBLE *restrict xmcc
          ,DOUBLE const viscosityC    ,DOUBLE const densityC
          ,DOUBLE const fArea         ,DOUBLE const dcca
          ,Loads ld                   ,short  const ndm 
          ,bool const fCalVel         ,bool const fCalPres){

  DOUBLE aP,wfn,m,tmp[3],gradVelFace[9];

/*... parade impermeavel movel*/
  if( ld.type == MOVEWALL){
    tA[0]   = ld.par[0];
    tA[1]   = ld.par[1];
    if( ndm == 3 )  tA[2]   = ld.par[2];
/*...*/
    if(fCalVel){
      aP   = viscosityC*fArea/dcca;
      *sP += aP;
      if( ndm == 2) {
/*... x*/
        p[0]  += aP*(tA[0]*(1.0-n[0]*n[0]) 
              + (velC[1]-tA[1])*n[1]*n[0]
              + velC[0]*n[0]*n[0]);
/*... y*/
        p[1]  += aP*(tA[1]*(1.0-n[1]*n[1]) 
              + (velC[0]-tA[0])*n[0]*n[1]
              + velC[1]*n[1]*n[1]);
      }
/*...................................................................*/

/*...*/
      else if( ndm == 3 ){
/*... x*/
        p[0]  += aP*(tA[0]*(1.0-n[0]*n[0]) 
              + (velC[1]-tA[1])*n[1]*n[0]
              + (velC[2]-tA[2])*n[2]*n[0]
              + velC[0]*n[0]*n[0]);
/*... y*/
        p[1]  += aP*(tA[1]*(1.0-n[1]*n[1]) 
              + (velC[0]-tA[0])*n[0]*n[1]
              + (velC[2]-tA[2])*n[2]*n[1]
              + velC[1]*n[1]*n[1]);
/*... z*/
        p[2]  += aP*(tA[2]*(1.0-n[2]*n[2]) 
              + (velC[0]-tA[0])*n[0]*n[2]
              + (velC[1]-tA[1])*n[1]*n[2]
              + velC[2]*n[2]*n[2]);
      }
/*...................................................................*/
    } 
/*...................................................................*/

  }
/*...................................................................*/

/*... entrada de massa*/
  else if( ld.type ==  INLET){
    if( ndm == 2) {
      tA[0] = ld.par[0];
      tA[1] = ld.par[1];
      wfn   = tA[0]*n[0] + tA[1]*n[1];
    }
    else{
      tA[0] = ld.par[0];
      tA[1] = ld.par[1];
      tA[2] = ld.par[2];
      wfn   = tA[0]*n[0] + tA[1]*n[1] + tA[2]*n[2] ;
    }
/*...*/
    m  = densityC*wfn*fArea;
    if(fCalVel){
      p[0] -= m*tA[0];
      p[1] -= m*tA[1];
      if(ndm == 3) p[2] -= m*tA[2];
    }
    if(fCalPres){
      p[0] -= m;
    }
  } 
/*...................................................................*/

/*... saida (derivada nula)*/
  else if( ld.type ==  OUTLET){
    if ( ndm == 2 ) 
      wfn = velC[0]*n[0] + velC[1]*n[1];
    else 
      wfn = velC[0]*n[0] + velC[1]*n[1] + velC[2]*n[2];

    m    = densityC*wfn*fArea;
/*... vb = vc + Grad(Vc)*r*/
    if(fCalVel){
/*...*/
      if( ndm == 2){
        gradVelFace[0] = 0.0e0;
        gradVelFace[1] = 0.0e0;
        gradVelFace[2] = 0.0e0;
        gradVelFace[3] = 0.0e0;
        
        gradFaceNull(gradVelFace,gradVel,xmcc,ndm);
        
        tmp[0] = MAT2D(0,0,gradVelFace,2)*xmcc[0] 
               + MAT2D(0,1,gradVelFace,2)*xmcc[1];      
         
        tmp[1] = MAT2D(1,0,gradVelFace,2)*xmcc[0] 
               + MAT2D(1,1,gradVelFace,2)*xmcc[1];      

        p[0] -= m*tmp[0];
        p[1] -= m*tmp[1];
      }
/*...................................................................*/

/*...*/
      else{
        gradVelFace[0] = 0.0e0;
        gradVelFace[1] = 0.0e0;
        gradVelFace[2] = 0.0e0;
        gradVelFace[3] = 0.0e0;
        gradVelFace[4] = 0.0e0;
        gradVelFace[5] = 0.0e0;
        gradVelFace[6] = 0.0e0;
        gradVelFace[7] = 0.0e0;
        gradVelFace[8] = 0.0e0;
        
        gradFaceNull(gradVelFace,gradVel,xmcc,ndm);
        
        tmp[0] = MAT2D(0,0,gradVelFace,3)*xmcc[0] 
               + MAT2D(0,1,gradVelFace,3)*xmcc[1]       
               + MAT2D(0,2,gradVelFace,3)*xmcc[2];      
        
        tmp[1] = MAT2D(1,0,gradVelFace,3)*xmcc[0] 
               + MAT2D(1,1,gradVelFace,3)*xmcc[1]       
               + MAT2D(1,2,gradVelFace,3)*xmcc[2];      
        
        tmp[2] = MAT2D(2,0,gradVelFace,3)*xmcc[0] 
               + MAT2D(2,1,gradVelFace,3)*xmcc[1]       
               + MAT2D(2,2,gradVelFace,3)*xmcc[2]; 
     
        p[0] -= m*tmp[0];
        p[1] -= m*tmp[1];
        p[2] -= m*tmp[2];
      }
/*...................................................................*/
      *sP += m;
    } 
/*...*/
    if(fCalPres){
      m     = densityC*wfn*fArea;
      p[0] -=  m;
    }
/*...................................................................*/
  }
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 01/07/2016                                   *
 * Data de modificaco : 09/07/2016                                   * 
 *-------------------------------------------------------------------* 
 * pLoadSimplePres : condicao de contorno para pressao               *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sP         -> termo da diagonal                                   * 
 * p          -> forca local                                         * 
 * tA         -> nao definido                                        * 
 * dField     -> coeficiente de D na face                            * 
 * densityC   -> massa especifica                                    * 
 * wfn        -> velocidade normal a face                            * 
 * xm         -> coordenada do ponto medio da face                   * 
 * fArea      -> area da face                                        * 
 * dcca       -> menor distancia do centroide central a face desta   *
 *               celula                                              * 
 * ld         -> definicao da carga                                  * 
 * fCal       -> true - atualizada sP e p                            * 
 *               false- atualizada sP e p                            * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal atualizado                           * 
 * p       -> forcas locais atualizada                               * 
 * tA      -> valor pescrito na face                                 * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void pLoadSimplePres(DOUBLE *restrict sP  ,DOUBLE *restrict p
          ,DOUBLE *restrict tA
          ,DOUBLE const dField    ,DOUBLE const densityC
          ,DOUBLE const wfn                                              
          ,DOUBLE const fArea     ,DOUBLE const dcca
          ,Loads ld               ,bool const fCal){

/*... pressao prescrita*/
  if( ld.type == DIRICHLETBC){
    tA[0]   = ld.par[0];
    if(fCal){
      *sP += densityC*fArea*dField/dcca;
    }
  }
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * pLoad : cargas                                                    *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal                                      * 
 * p       -> forca local                                            * 
 * tA      -> nao definido                                           * 
 * coefDifC-> coeficiente de difusao                                 * 
 * densityC-> densidade                                              * 
 * wfn     -> velocidade normal a face                               * 
 * xm      -> coordenada do ponto medio da face                      * 
 * fArea   -> area da face                                           * 
 * dcca    -> menor distancia do centroide central a face desta      *
 *            celula                                                 * 
 * ld      -> definicao da carga                                     * 
 * fCal    -> true - atualizada sP e p                               * 
 *            false- atualizada sP e p                               * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal atualizado                           * 
 * p       -> forcas locais atualizada                               * 
 * tA      -> valor pescrito na face                                 * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void pLoad(DOUBLE *restrict sP  ,DOUBLE *restrict p
          ,DOUBLE *restrict tA
          ,DOUBLE const coefDifC,DOUBLE const densityC
          ,DOUBLE const wfn     ,DOUBLE *restrict xm                   
          ,DOUBLE const fArea   ,DOUBLE const dcca
          ,Loads ld             ,bool const fCal){

  DOUBLE aP,h;


/*... potencial prescrito*/
  if( ld.type == DIRICHLETBC){
    tA[0]   = ld.par[0];
/*...*/
    if(fCal){
      aP   = coefDifC*fArea/dcca;
      *sP += aP;
      *p  += aP*tA[0];
    } 
/*...................................................................*/
  }
/*...................................................................*/

/*... lei de resfriamento de newton*/
  else if( ld.type == ROBINBC){
    h     = ld.par[0];
    tA[0] = ld.par[1];
/*...*/
    if(fCal){
      aP  = ((coefDifC*h)/(coefDifC+h*dcca))*fArea;
      *sP += aP;
      *p  += aP*tA[0];
    }    
/*...................................................................*/
  }
/*...................................................................*/

/*... fluxo prestrito diferente de zero*/
   else if( ld.type == NEUMANNBC){
     tA[0]   = ld.par[0];
/*...*/
     if(fCal)
       *p += fArea*tA[0];
/*...................................................................*/
   }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
   else if( ld.type == SINBC){
     loadSenProd(tA,ld.par,xm); 
/*...*/
     if(fCal){ 
       aP   = coefDifC*fArea/dcca;
       *sP += aP;
       *p  += aP*tA[0];
     }  
/*...................................................................*/
   }
/*...................................................................*/

/*... potencial prescrito (entra)*/
   else if( ld.type == INLET){
     tA[0]   = ld.par[0];
/*...*/
     if(fCal)
       *p -= wfn*densityC*fArea*tA[0];
   }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
   else if( ld.type == OUTLET){
/*...*/
     if(fCal)
       *sP += wfn*densityC*fArea;
/*...................................................................*/
   }
/*...................................................................*/


}
/*********************************************************************/

     
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * LIMITFACEBASE : funcao limitadora de fluxo baseado na face        *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * r       -> parametro utilizado                                    * 
 * iCod    -> tipo de funcao limitadora                              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE limitFaceBase(DOUBLE const r,short const iCod)
{
  #define NFUNCLIMTEB 5
  DOUBLE a,b,c;
  short i;
  char word[][WORD_SIZE]=
       {"VanLeer","VanAlbada"
       ,"MidMod ","Osher"   ,"SuperBee"};

  switch(iCod){

/*... Van Leer - TVD*/
    case VANLEERFACE: 
      a   = fabs(r);
      return (r + a ) / (1.0e0 + a);
    break;
/*...................................................................*/

/*... Van Albada  - TVD*/
    case VANALBADAFACE:
      if( r > 0.e0){ 
        a = r*r;
        return (r + a ) / (1.0e0 + a);
      }
      else 
        return 0.e0;
    break;
/*...................................................................*/

/*... Mid-Mod - TVD*/
    case MIDMODFACE:
      return max(0.e0,min(r,1.e0));
    break;
/*...................................................................*/

/*... OSHER - TVD*/
    case OSHERFACE:
      return  max(0.e0,min(r,2.e0));
    break;
/*...................................................................*/

/*... SUPERBEE - TVD*/
    case SUPERBEEFACE:
      a   = min(2.e0*r,1.e0);
      b   = min(r,2.e0);
      c   = max(a,b);
      return max(0.e0,c);
    break;
/*...................................................................*/

/*...*/
    default:
      printf("Erro: tipo de funcao limiradora fluxo invalida.\n"
             "Arquivo fonte:  \"%s\".\n"
             "Nome da funcao: \"%s\".\n"
             "Linha         : \"%d\".\n"
             ,__FILE__,__func__,__LINE__);
      printf("Funcoes disponiveis:\n");
      for(i=0;i<NFUNCLIMTEB;i++)
        printf("%s\n",word[i]);
      exit(EXIT_FAILURE);
/*...................................................................*/
  }
  return 0.e0;
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * FACEBASETVD : metodo TVD baseado na face                          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE faceBaseTvd(short const nAresta    ,short const idCell
                 ,DOUBLE *restrict u0
                 ,DOUBLE *restrict gradUv,DOUBLE *restrict gradUp
                 ,DOUBLE *restrict lKsi  ,DOUBLE const lModKsi 
                 ,DOUBLE const cv
                 ,short const iCod       ,short const ndm)
{                    
  
  DOUBLE du,r,cvc=0.0e0,gf;
  DOUBLE eps=1.e-14;
 
  if(ndm == 2)
    if( cv < 0.0e0) {
      du    = u0[idCell] - u0[nAresta] + eps;
      gf    = -(gradUv[0]*lKsi[0] + gradUv[1]*lKsi[1])*lModKsi;
      r     = 2.0e0*gf/du - 1.0e0;
      cvc   = 0.5e0*limitFaceBase(r ,iCod)*du;
    }
    else{
      du    = u0[nAresta] - u0[idCell] + eps;
      gf    = (gradUp[0]*lKsi[0] + gradUp[1]*lKsi[1])*lModKsi;
      r     = 2.0e0*gf/du - 1.0e0;
      cvc   = 0.5e0*limitFaceBase(r,iCod)*du;
    }
  else if(ndm ==3){
    if( cv < 0.0e0) {
      du    = u0[idCell] - u0[nAresta] + eps;
      gf    = -( gradUv[0]*lKsi[0] 
               + gradUv[1]*lKsi[1]
               + gradUv[2]*lKsi[2])*lModKsi;
      r     = 2.0e0*gf/du - 1.0e0;
      cvc   = 0.5e0*limitFaceBase(r ,iCod)*du;
    }
    else{
      du    = u0[nAresta] - u0[idCell] + eps;
      gf    = (gradUp[0]*lKsi[0] 
             + gradUp[1]*lKsi[1]
             + gradUp[2]*lKsi[2])*lModKsi;
      r     = 2.0e0*gf/du - 1.0e0;
      cvc   = 0.5e0*limitFaceBase(r,iCod)*du;
    }
  }
 
  return cvc;
} 
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 04/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * FACEBASETVDV1 : metodo TVD baseado na face                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * uC     -> valor central                                           * 
 * uV     -> valor do vizinho                                        * 
 * gradUc -> gradiente central                                       * 
 * gradUv -> gradiente do vizinho                                    * 
 * lKsi   -> vetor unitario entre os centroides                      * 
 * lModKsi-> distancia entre os centroides                           * 
 * cv     ->                                                         * 
 * iCod   -> tecnica TVD                                             * 
 * ndm    -> numero de dimensoes                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE faceBaseTvdV1(DOUBLE const uC     ,DOUBLE const uV
                 ,DOUBLE *restrict gradUc,DOUBLE *restrict gradUv
                 ,DOUBLE *restrict lKsi  ,DOUBLE const lModKsi 
                 ,DOUBLE const cv
                 ,short const iCod       ,short const ndm)
{                    
  
  DOUBLE du,r,cvc=0.0e0,gf;
  DOUBLE eps=1.e-14;
 
  if(ndm == 2)
    if( cv < 0.0e0) {
      du    = uC - uV + eps;
      gf    = -(gradUv[0]*lKsi[0] + gradUv[1]*lKsi[1])*lModKsi;
      r     = 2.0e0*gf/du - 1.0e0;
      cvc   = 0.5e0*limitFaceBase(r ,iCod)*du;
    }
    else{
      du    = uV - uC + eps;
      gf    = (gradUc[0]*lKsi[0] + gradUc[1]*lKsi[1])*lModKsi;
      r     = 2.0e0*gf/du - 1.0e0;
      cvc   = 0.5e0*limitFaceBase(r,iCod)*du;
    }
  else if(ndm ==3){
    if( cv < 0.0e0) {
      du    = uC - uV + eps;
      gf    = -( gradUv[0]*lKsi[0] 
               + gradUv[1]*lKsi[1]
               + gradUv[2]*lKsi[2])*lModKsi;
      r     = 2.0e0*gf/du - 1.0e0;
      cvc   = 0.5e0*limitFaceBase(r ,iCod)*du;
    }
    else{
      du    = uV - uC + eps;
      gf    = (gradUc[0]*lKsi[0] 
             + gradUc[1]*lKsi[1]
             + gradUc[2]*lKsi[2])*lModKsi;
      r     = 2.0e0*gf/du - 1.0e0;
      cvc   = 0.5e0*limitFaceBase(r,iCod)*du;
    }
  }
 
  return cvc;
} 
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 16/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * UPWINDLINEARV1: upwind linear com resconstrucao linear            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * uC     -> valor central                                           * 
 * uV     -> valor do vizinho                                        * 
 * gradUc -> gradiente central                                       * 
 * gradUv -> gradiente do vizinho                                    * 
 * r      -> distancia ate o ponto central da face                   * 
 * wfn    -> velocidade normal                                       * 
 * iCod   -> tecnica TVD                                             * 
 * ndm    -> numero de dimensoes                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 * Tface = TfaceUp + grad(T)*r                                       *   
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE upwindLinearV1(DOUBLE const uC     ,DOUBLE const uV
                 ,DOUBLE *restrict gradUc,DOUBLE *restrict gradUv
                 ,DOUBLE *restrict r     ,DOUBLE const wfn
                 ,short const iCod       ,short const ndm)
{                    
  
  DOUBLE cvc=0.0e0;
 
  if(ndm == 2){
    if( wfn < 0.0e0) 
      cvc = gradUv[0]*r[0] +gradUv[1]*r[1];
    else 
      cvc = gradUc[0]*r[0] +gradUc[1]*r[1];
  }
  else if(ndm ==3){
    if( wfn < 0.0e0) 
      cvc = gradUv[0]*r[0] +gradUv[1]*r[1] + gradUv[2]*r[2] ;
    else 
      cvc = gradUc[0]*r[0] +gradUc[1]*r[1] + gradUc[2]*r[2] ;
  }
 
  return cvc;
} 
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 16/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * DEFERREDCD : Diferencao central com correcao atrasada             * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * uC     -> valor central                                           * 
 * uV     -> valor do vizinho                                        * 
 * wfn    -> velocidade normal                                       * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------*
 * valor na face = (Upwind)implicito + (Central - Upwind)explicito   * 
 *********************************************************************/
DOUBLE deferredCd(DOUBLE const uC,DOUBLE const uV,DOUBLE const wfn)
{                    
  
  DOUBLE cvc=0.0e0;
 
  if(wfn < 0.0e0) {
    cvc   = 0.5e0*(uC+uV) - uV;
  }
  else{
    cvc   = 0.5e0*(uC+uV) - uC;
  }
 
  return cvc;
} 
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * SETFACEBASE :                                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * word ->                                                           * 
 * iCod -> nao definido                                              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * iCod -> codigo da tecnica do termo advectivo                      * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void setFaceBase(char *word,short *iCod)
{
  #define NFUNCLIMTFACE 7
  short i;
  char fBase[][WORD_SIZE]=
       {"FoUp"   ,"VanLeer" ,"VanAlbada"
       ,"MidMod ","Osher"   ,"SuperBee"
       ,"CD"     ,"SoUp" };
/*...*/  
  if(!strcmp(word,"FoUp")){
    *iCod = FOUP; 
    if(!mpiVar.myId ) printf("iCod  : FOUP\n");
  }
/*...................................................................*/
  
/*...*/
  else if(!strcmp(word,"CD")){
    *iCod = CD; 
    if(!mpiVar.myId ) printf("iCod  : CD\n");
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,"VanLeer")){
    *iCod = VANLEERFACE; 
    if(!mpiVar.myId ) printf("iCod  : VanLeer\n");
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,"VanAlbada")){
    *iCod =  VANALBADAFACE; 
    if(!mpiVar.myId ) printf("iCod  : VanAlbade\n");
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,"MidMod")){
    *iCod =  MIDMODFACE; 
    if(!mpiVar.myId ) printf("iCod  : MidMod\n");
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,"Osher")){
    *iCod =  OSHERFACE; 
    if(!mpiVar.myId ) printf("iCod  : Osher\n");
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,"SuperBee")){
    *iCod =  SUPERBEEFACE; 
    if(!mpiVar.myId ) printf("iCod  : SuperBee\n");
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,"SuperBee")){
    *iCod =  SUPERBEEFACE; 
    if(!mpiVar.myId ) printf("iCod  : SuperBee\n");
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,"SoUp")){
    *iCod =  SOUP; 
    if(!mpiVar.myId ) printf("iCod  : SoUp\n");
  }
/*...................................................................*/

/*...*/
  else{
    printf("Erro: tipo de funcao limitadora fluxo invalida.\n"
           "Arquivo fonte:  \"%s\".\n"
           "Nome da funcao: \"%s\".\n"
           "Linha         : \"%d\".\n"
           ,__FILE__,__func__,__LINE__);
    printf("Funcoes disponiveis:\n");
    for(i=0;i<NFUNCLIMTFACE;i++)
      printf("%s\n",fBase[i]);
    exit(EXIT_FAILURE);
  }
/*...................................................................*/
}
/*********************************************************************/ 

       
/********************************************************************* 
 * Data de criacao    : 18/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * SIZECAR: tamanho caracteristico da celula                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * volume -> area ou voluma                                          * 
 * ndm    -> dimensao                                                * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE sizeCar(DOUBLE const volume,short const ndm)
{

  if(ndm == 2)      return sqrt(volume);
  else if(ndm == 3) return pow(volume,oneDivTree);

  return 0.0;
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 22/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * GRADFACENULL :  garante o gradiente nulo da face usando uma extra * 
 * polacao no centroide                                              * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void gradFaceNull(DOUBLE *restrict gradVelFace
                 ,DOUBLE *restrict gradVelCell
                 ,DOUBLE *restrict xmcc       ,short const ndm)
{

  DOUBLE prod[3],tensor[3][3],eFace[3],mod;

/*...................................................................*/
  if(ndm == 2){
/*... vetor unitario na direcao do centroide ate o ponto medio da face
      eFace */
    mod      = sqrt(xmcc[0]*xmcc[0] + xmcc[1]*xmcc[1]);
    eFace[0] = xmcc[0]/mod;
    eFace[1] = xmcc[1]/mod;
/*...................................................................*/

/*... (Grad(VelP),eFace)*/
    prod[0] = MAT2D(0,0,gradVelCell,2)*eFace[0]
            + MAT2D(0,1,gradVelCell,2)*eFace[1];

    prod[1] = MAT2D(1,0,gradVelCell,2)*eFace[0]
            + MAT2D(1,1,gradVelCell,2)*eFace[1];
/*...................................................................*/

/*... produto tensorial entre dois vetores*/
    tensor[0][0] = prod[0]*eFace[0];
    tensor[0][1] = prod[0]*eFace[1];
    
    tensor[1][0] = prod[1]*eFace[0];
    tensor[1][1] = prod[1]*eFace[1];

/*...................................................................*/

/*... Grad(VelFace) = Grad(VelP) - (Grad(VelP),eFace)eFace*/
     MAT2D(0,0,gradVelFace,2) =  MAT2D(0,0,gradVelCell,2)
                              -  tensor[0][0];        

     MAT2D(0,1,gradVelFace,2) =  MAT2D(0,1,gradVelCell,2)
                              -  tensor[0][1];        

     MAT2D(1,0,gradVelFace,2) =  MAT2D(1,0,gradVelCell,2)
                              -  tensor[1][0];        

     MAT2D(1,1,gradVelFace,2) =  MAT2D(1,1,gradVelCell,2)
                              -  tensor[1][1];        
/*...................................................................*/

  }
/*...................................................................*/

/*...*/  
  else if(ndm == 3){
/*... vetor unitario na direcao do centroide ate o ponto medio da face
      eFace */
    mod      = sqrt(xmcc[0]*xmcc[0] 
                  + xmcc[1]*xmcc[1] 
                  + xmcc[2]*xmcc[2]);
    eFace[0] = xmcc[0]/mod;
    eFace[1] = xmcc[1]/mod;
    eFace[2] = xmcc[2]/mod;
/*...................................................................*/

/*... (Grad(VelP),eFace)*/
    prod[0] = MAT2D(0,0,gradVelCell,3)*eFace[0]
            + MAT2D(0,1,gradVelCell,3)*eFace[1] 
            + MAT2D(0,2,gradVelCell,3)*eFace[2];

    prod[1] = MAT2D(1,0,gradVelCell,3)*eFace[0]
            + MAT2D(1,1,gradVelCell,3)*eFace[1] 
            + MAT2D(1,2,gradVelCell,3)*eFace[2];
    
    prod[2] = MAT2D(2,0,gradVelCell,3)*eFace[0]
            + MAT2D(2,1,gradVelCell,3)*eFace[1] 
            + MAT2D(2,2,gradVelCell,3)*eFace[2];
/*...................................................................*/

/*... produto tensorial entre dois vetores*/
    tensor[0][0] = prod[0]*eFace[0];
    tensor[0][1] = prod[0]*eFace[1];
    tensor[0][2] = prod[0]*eFace[2];
    
    tensor[1][0] = prod[1]*eFace[0];
    tensor[1][1] = prod[1]*eFace[1];
    tensor[1][2] = prod[1]*eFace[2];
    
    tensor[2][0] = prod[2]*eFace[0];
    tensor[2][1] = prod[2]*eFace[1];
    tensor[2][2] = prod[2]*eFace[2];

/*...................................................................*/

/*... Grad(VelFace) = Grad(VelP) - (Grad(VelP),eFace)eFace*/
     MAT2D(0,0,gradVelFace,3) =  MAT2D(0,0,gradVelCell,3)
                              -  tensor[0][0];        
     MAT2D(0,1,gradVelFace,3) =  MAT2D(0,1,gradVelCell,3)
                              -  tensor[0][1];        
     MAT2D(0,2,gradVelFace,3) =  MAT2D(0,2,gradVelCell,3)
                              -  tensor[0][2];        

     MAT2D(1,0,gradVelFace,3) =  MAT2D(1,0,gradVelCell,3)
                              -  tensor[1][0];        
     MAT2D(1,1,gradVelFace,3) =  MAT2D(1,1,gradVelCell,3)
                              -  tensor[1][1];        
     MAT2D(1,2,gradVelFace,3) =  MAT2D(1,2,gradVelCell,3)
                              -  tensor[1][2];        
     
     MAT2D(2,0,gradVelFace,3) =  MAT2D(2,0,gradVelCell,3)
                              -  tensor[2][0];        
     MAT2D(2,1,gradVelFace,3) =  MAT2D(2,1,gradVelCell,3)
                              -  tensor[2][1];        
     MAT2D(2,2,gradVelFace,3) =  MAT2D(2,2,gradVelCell,3)
                              -  tensor[2][2];        
/*...................................................................*/

  }


}
/*********************************************************************/ 
