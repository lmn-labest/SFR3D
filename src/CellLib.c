#include<CellLoop.h>
/********************************************************************* 
 * CELLLIBDIF : chamada de bibliotecas de celulas para o problema    *
 * de difusao.                                                       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
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
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * u0        -> solucao conhecida                                    * 
 * gradU0    -> gradiente rescontruido da solucao conhecida          * 
 * faceR     -> restricoes por elmento                               * 
 * faceS     -> carga por elemento                                   * 
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
void cellLibDif(short *restrict lGeomType,DOUBLE *restrict lprop
               ,INT   *restrict lViz     ,INT *restrict lId  
               ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
               ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea
               ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
               ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
               ,DOUBLE *restrict dcca    
               ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
               ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
               ,DOUBLE *restrict lRcell 
               ,short  *restrict lFaceR  ,DOUBLE *restrict lFaceS
               ,DOUBLE *restrict u0      ,DOUBLE *restrict gradU0
               ,short const nEn          ,short  const nFace     
               ,short const ndm          ,short const lib    
               ,INT const nel)
{

/*... difusao pura */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2)
      cellDif2D(lGeomType,lprop
               ,lViz     ,lId
               ,ksi      ,mKsi
               ,eta      ,fArea
               ,normal   ,volume
               ,xm       ,xmcc
               ,dcca     
               ,vSkew    ,mvSkew
               ,lA       ,lB
               ,lRcell 
               ,lFaceR   ,lFaceS
               ,u0       ,gradU0      
               ,nEn      ,nFace 
               ,ndm      ,nel);
/*..................................................................*/   

/*... 3D*/
    else if(ndm == 3){
      cellDif3D(lGeomType,lprop
               ,lViz     ,lId
               ,ksi      ,mKsi
               ,eta      ,fArea
               ,normal   ,volume
               ,xm       ,xmcc
               ,dcca     
               ,vSkew    ,mvSkew
               ,lA       ,lB
               ,lRcell 
               ,lFaceR   ,lFaceS
               ,u0       ,gradU0      
               ,nEn      ,nFace 
               ,ndm      ,nel);
    } 
/*..................................................................*/   
  }

}
/*********************************************************************/
 
/********************************************************************* 
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
 * meta   -> modulo do vetor eta                                     * 
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
               ,short maxNo               ,short maxViz
               ,short ndm                 ,INT nel)
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

/*... calculo do vetor normal e da area da face*/
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

/*... calculo do vetor normal e da area da face*/
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
 * CELLRCGRAD : chamada de bibliotecas de celulas para a reconstrucao*
 * de gradiente.                                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lViz      -> viznhos da celula central                            *
 * lSquare -> matriz para a reconstrucao least Square                * 
 * lSquareR-> fatoracao R (RCLSQUAREQR)                              * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * meta      -> modulo do vetor eta                                  * 
 * normal    -> vetores normais as faces das celulas                 * 
 * area      -> area da celula central                               * 
 * vSkew      -> vetor entre o ponto medio a intersecao que une os   * 
 *            centrois compartilhado nessa face                      * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * u         -> solucao conhecida na celula                          * 
 * gradU     -> gradiente rescontruido da solucao conhecida          * 
 * nU        -> solucao conhecida no no                              * 
 * faceR     -> restricoes por elmento                               * 
 * faceS     -> carga por elemento                                   * 
 * nFace     -> carga por elemento                                   * 
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
void cellLibRcGrad(INT   *restrict lViz    
                 ,DOUBLE *restrict lLsquare,DOUBLE *restrict lLsquareR
                 ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
                 ,DOUBLE *restrict eta     ,DOUBLE *restrict mEta
                 ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
                 ,DOUBLE *restrict vSkew   ,DOUBLE *restrict xmcc 
                 ,short  *restrict lFaceR  ,DOUBLE *restrict lFaceS
                 ,DOUBLE *restrict u       ,DOUBLE *restrict gradU 
                 ,DOUBLE *restrict lnU     ,short const ty                
                 ,short const nFace        ,short const ndm      
                 ,short const lib          ,short const ndf
                 ,short *restrict  isNod   ,INT const nel){
  long aux;
    
  switch(lib){
/*... green-Gauss linear baseado na celula*/  
    case RCGRADGAUSSC:
      greenGaussCell(lViz    ,mKsi
                    ,eta     ,mEta
                    ,normal  ,volume
                    ,vSkew   ,xmcc
                    ,lFaceR  ,lFaceS
                    ,u       ,gradU 
                    ,nFace   ,ndm   
                    ,ndf);
    break;
/*...................................................................*/ 
    
/*... green-Gauss linear baseado no no*/  
    case RCGRADGAUSSN:
       greenGaussNode(lViz    ,mEta
                     ,normal  ,volume
                     ,lnU     ,gradU 
                     ,isNod   
                     ,nFace   ,ndm   
                     ,ndf     ,ty);
    break;
/*...................................................................*/ 
   
/*... minimos quadrados*/  
    case RCLSQUARE:
      leastSquare(lLsquare,lViz
                 ,u       ,gradU
                 ,lFaceR  ,lFaceS
                 ,nFace   ,ndf
                 ,ndm);
    break;
/*...................................................................*/ 

/*... minimos quadrados QR*/  
    case RCLSQUAREQR:
      leastSquareQR(lLsquare,lLsquareR
                  ,lViz
                  ,u       ,gradU
                  ,lFaceR  ,lFaceS
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
 * GRREENGAUSSCELL: reconstrucao de gradiente green-gauss linear por * 
 * celula                                                            *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lViz      -> viznhos da celula central                            * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * meta      -> modulo do vetor eta                                  * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> area da celula central                               * 
 * vSkew      -> vetor entre o ponto medio a intersecao que une os   * 
 *            centrois compartilhado nessa face                      * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * lFaceR    -> restricoes por elmento                               * 
 * lFaceS    -> carga por elemento                                   * 
 * u         -> solucao conhecida                                    * 
 * gradU     -> gradiente rescontruido da solucao conhecida          * 
 * nFace     -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * ndf       -> grauss de liberdade                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * gradU     -> gradiente calculado                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void greenGaussCell(INT *restrict lViz   ,DOUBLE *restrict mKsi
               ,DOUBLE *restrict eta     ,DOUBLE *restrict mEta
               ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
               ,DOUBLE *restrict vSkew   ,DOUBLE *restrict xmcc 
               ,short  *restrict lFaceR  ,DOUBLE *restrict lFaceS
               ,DOUBLE *restrict u       ,DOUBLE *restrict gradU 
               ,short const nFace        ,short const ndm   
               ,short const ndf)
{
  DOUBLE lNormal[3],v,dPviz,aux1;
  DOUBLE uf[MAX_NDF],uC[MAX_NDF];
  DOUBLE lModKsi,alpha,alphaMenosUm,invVol;
  INT vizNel;
  short idCell = nFace;
  short i,j,k;

  invVol = 1.e0/volume[idCell];
/*...*/

  for(i=0;i<ndf;i++)
    uC[i] = MAT2D(idCell,i,u,ndf);
   
  for(i=0;i<ndf*ndm;i++)
    gradU[i]  = 0.e0;
      
/*... */
  for(i=0;i<nFace;i++){
    vizNel = lViz[i];
    for(j=0;j<ndm;j++)
      lNormal[j] = MAT2D(i,j,normal,ndm);
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
      dPviz = sqrt(dPviz);
      alpha        = dPviz/lModKsi;
      alphaMenosUm = 1.0 - alpha;
/*...................................................................*/

/*...*/
      for(j=0;j<ndf;j++)
        uf[j] = alphaMenosUm*uC[j] + alpha*MAT2D(i,j,u,ndf);
/*...................................................................*/
    }
/*...................................................................*/

/*... contorno*/
    else{
/*... temperatura prescrita na face*/
      if(lFaceR[i])
        for(j=0;j<ndf;j++)
          uf[j] = MAT2D(i,j,lFaceS,ndf); 

/*... temperatura prescrita na celula*/
      else if(lFaceR[nFace]==VPES)
        for(j=0;j<ndf;j++)
          uf[j] = MAT2D(idCell,j,lFaceS,ndf); 
      
/*... fluxo prescrito*/
      else 
        for(j=0;j<ndf;j++)
          uf[j] = MAT2D(idCell,j,u,ndf); 
    }
/*...................................................................*/
    
/*...*/
    for(j=0;j<ndf;j++){
      aux1 = uf[j]*mEta[i]; 
      for(k=0;k<ndm;k++)
        MAT2D(j,k,gradU,ndm) += aux1*lNormal[k]; 
    }
/*...................................................................*/
  }
/*...................................................................*/
 
  for(i=0;i<ndf;i++)
    for(j=0;j<ndm;j++)
      MAT2D(i,j,gradU,ndm) *= invVol; 
}
/*********************************************************************/

/********************************************************************* 
 * GRREENGAUSSNODE: reconstrucao de gradiente green-gauss linear por * 
 * celula                                                            *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lViz      -> viznhos da celula central                            * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * meta      -> modulo do vetor eta                                  * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> area da celula central                               * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * mvSkew       -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
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
void greenGaussNode(INT *restrict lViz   ,DOUBLE *restrict mEta
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
      aux1 = uf[j]*mEta[i]; 
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
 * LEASTSQUARE : calcula o gradiente por minimos quadrados           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * lLsquare  -> matriz para a reconstrucao least Square              * 
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
void  leastSquare(DOUBLE *restrict lLsquare,INT *restrict lViz 
                 ,DOUBLE *restrict u       ,DOUBLE *restrict gradU
                 ,short  *restrict lFaceR  ,DOUBLE *restrict lFaceS
                 ,short const nFace        ,short const ndf
                 ,short const ndm){

  DOUBLE du[MAX_NUM_FACE*MAX_NDF],uC[MAX_NDF];
  INT vizNel;
  short idCell = nFace;
  short i,j;

  for(i=0;i<ndf*ndm;i++){
    gradU[i]  = 0.e0;
  }

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
        if(lFaceR[i])
          du[i] = lFaceS[i] - uC[0]; 

/*... temperatura prescrita na celula*/
        else if(lFaceR[nFace]==VPES)
          du[i] = lFaceS[idCell] - uC[0]; 
      
/*... fluxo prescrito*/
        else {
          du[i] = lFaceS[i];
        }
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... gradU = G du*/
    for(i=0;i<ndm;i++)
      for(j=0;j<nFace;j++){
        gradU[i] += MAT2D(i,j,lLsquare,nFace)*du[j]; 
      } 
  }
/*...................................................................*/

/*... */  
//else
//  for(k=0;k<ndf;k++)
//    for(i=0;i<ndm;i++)
//      for(j=0;j<nFace;i++)
//        MAT2D(k,i,grad,ndf) 
//        += MAT2D(i,j,lLsquare,nFace)*MAT2D(k,j,du,nFace);  
/*...................................................................*/

} 
/********************************************************************* 
 * LEASTSQUAREQR: calcula o gradiente por minimos quadrados           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * lLsquareQt-> matriz para a reconstrucao least Square              * 
 * lLsquareR -> fatoracao R (RCLSQUAREQR)                            * 
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
void  leastSquareQR(DOUBLE *restrict lLsquareQt
                   ,DOUBLE *restrict lLsquareR
                   ,INT *restrict lViz 
                   ,DOUBLE *restrict u       ,DOUBLE *restrict gradU
                   ,short  *restrict lFaceR  ,DOUBLE *restrict lFaceS
                   ,short const nFace        ,short const ndf
                   ,short const ndm){

  DOUBLE du[MAX_NUM_FACE*MAX_NDF],uC[MAX_NDF];
  DOUBLE  b[MAX_NUM_FACE*MAX_NDF];
  INT vizNel;
  short idCell = nFace;
  short i,j;

  for(i=0;i<ndf*ndm;i++){
    gradU[i]  = 0.e0;
    b[i]      = 0.e0;
  }

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
        if(lFaceR[i])
          du[i] = lFaceS[i] - uC[0]; 

/*... temperatura prescrita na celula*/
        else if(lFaceR[nFace]==VPES)
          du[i] = lFaceS[idCell] - uC[0]; 
      
/*... fluxo prescrito*/
        else {
          du[i] = lFaceS[i];
        }
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... b = Qtdu*/
    for(i=0;i<ndm;i++)
      for(j=0;j<nFace;j++){
        b[i] += MAT2D(i,j,lLsquareQt,nFace)*du[j]; 
      } 
/*...................................................................*/
    
/*... R grad = Qtdu*/
    if(ndm == 2){
/*... retrosubstituicao*/
      gradU[1] = b[1]/lLsquareR[2];
      gradU[0] = (b[0] - lLsquareR[1]*gradU[1])/lLsquareR[0];
    }
    else if(ndm == 3){
/*... retrosubstituicao*/
      gradU[2] = b[2]/lLsquareR[5];
      gradU[1] = (b[1] - lLsquareR[4]*gradU[2])/lLsquareR[3];
      gradU[0] = (b[0] 
               - lLsquareR[1]*gradU[1] 
               - lLsquareR[2]*gradU[2])/lLsquareR[0];
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... */  
//else
//  for(k=0;k<ndf;k++)
//    for(i=0;i<ndm;i++)
//      for(j=0;j<nFace;i++)
//        MAT2D(k,i,grad,ndf) 
//        += MAT2D(i,j,lLsquare,nFace)*MAT2D(k,j,du,nFace);  
/*...................................................................*/

} 


/********************************************************************* 
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
/*... dimensao 2*/
      if(ndm == 2){
        MAT2D(nf,0,dx,ndm) = MAT2D(nf,0,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,1,dx,ndm) = MAT2D(nf,1,lKsi,ndm)*lmKsi[nf];
      }
/*...................................................................*/

/*... dimensao 3*/
      else if(ndm == 3){
        MAT2D(nf,0,dx,ndm) = MAT2D(nf,0,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,1,dx,ndm) = MAT2D(nf,1,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,2,dx,ndm) = MAT2D(nf,2,lKsi,ndm)*lmKsi[nf];
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
/*...................................................................*/
    }

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
        MAT2D(0,nf,lLsquare,lnFace) = inv[0]*MAT2D(nf,0,dx,ndm);
        MAT2D(1,nf,lLsquare,lnFace) = inv[2]*MAT2D(nf,1,dx,ndm);
      }
/*...................................................................*/

/*... Qt - dimensao 3*/
      else if(ndm == 3){
        MAT2D(0,nf,lLsquare,lnFace) = MAT2D(nf,0,dx,ndm);
        MAT2D(1,nf,lLsquare,lnFace) = MAT2D(nf,1,dx,ndm);
        MAT2D(2,nf,lLsquare,lnFace) = MAT2D(nf,2,dx,ndm);
      }
/*...................................................................*/
    }
/*...................................................................*/
  } 
/*...................................................................*/
   
  }
/*********************************************************************/

/********************************************************************* 
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
inline short sn(short *s,short ty,INT nel)
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
inline double areaCell(double *eta,short ty,short ndm,INT nel)
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
               ,aux,__FILE__,__func__);
      exit(EXIT_FAILURE);
    break; 
/*...................................................................*/
  }
}
/*********************************************************************/ 


/********************************************************************* 
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
inline double areaTriaCell(double *restrict eta, short ndm)
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
inline double areaQuadCell(double *restrict eta,short ndm)
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
 * nFace  -> numero de face da celula                                * 
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
