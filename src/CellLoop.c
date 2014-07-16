#include<CellLoop.h>
/********************************************************************* 
 * PGEOMFORM: calculo da propriedades geometricas das celulas        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x       -> cordenadas dos pontos                                  * 
 * el      -> conetividade dos celulas                               * 
 * nelcon  -> vizinhos dos elementos                                 * 
 * nen     -> numero de nos por celulas                              * 
 * nFace   -> numero de faces por celulas                            * 
 * geomType-> tipo geometrico das celulas                            * 
 * gCc     -> indefinido                                             * 
 * gKsi    -> indefinido                                             * 
 * gmKsi   -> indefinido                                             * 
 * gEta    -> indefinido                                             * 
 * gmEta   -> indefinido                                             * 
 * gNormal -> indefinido                                             * 
 * gVolume -> indefinido                                             * 
 * gXm     -> indefinido                                             * 
 * gXmcc   -> indefinido                                             * 
 * gmkm    -> indefinido                                             * 
 * gDcca   -> indefinido                                             * 
 * maxNo   -> numero de nos por celula maximo da malha               * 
 * maxViz  -> numero vizinhos por celula maximo da malha             * 
 * ndm     -> numero de dimensoes                                    * 
 * numel   -> numero de toral de celulas                             * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * gCc     -> centroide das celulas                                  * 
 * gKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        * 
 * gmKsi   -> modulo do vetor ksi                                    * 
 * gEta    -> vetores paralelos as faces das celulas                 * 
 * gmEta   -> modulo do vetor eta                                    * 
 * gNormal -> vetores normais as faces das celulas                   * 
 * gVolume -> volumes das celulas                                    * 
 * gXm     -> pontos medios das faces das celulas                    * 
 * gXmcc   -> vetores que unem o centroide aos pontos medios das     * 
 *            faces                                                  * 
 * gmkm    -> distacia entre o ponto medio a intersecao que une os   * 
 *            centrois compartilhado nessa face                      * 
 * gDcca   -> menor distancia do centroide a faces desta celula      * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void pGeomForm(double *restrict x      ,INT    *restrict el
              ,INT    *restrict nelcon ,short  *restrict nen    
              ,short  *restrict nFace  ,short  *restrict geomType
              ,double *restrict gCc
              ,double *restrict gKsi   ,double *restrict gmKsi
              ,double *restrict gEta   ,double *restrict gmEta
              ,double *restrict gNormal,double *restrict gVolume
              ,double *restrict gXm    ,double *restrict gXmcc   
              ,double *restrict gmKm   ,double *restrict gDcca                 
              ,short maxNo             ,short maxViz
              ,short ndm               ,INT numel)
{
  INT nel,no,vizNel;
  short i,j,k;
/*... variavel local */
  double lx[MAX_NUM_PONT];
  double lCc[(MAX_NUM_FACE+1)*3];
  double lKsi[MAX_NUM_FACE*3],lmKsi[MAX_NUM_FACE];
  double lEta[MAX_NUM_FACE*3],lmEta[MAX_NUM_FACE];
  double lNormal[MAX_NUM_FACE*3],lVolume;
  double lXm[MAX_NUM_FACE*3],lXmcc[MAX_NUM_FACE*3];
  double lDcca[MAX_NUM_FACE];
  double lmKm[MAX_NUM_FACE];  
  short  lnFace[MAX_NUM_FACE+1],lGeomType[MAX_NUM_FACE+1],ty;
  short  s[MAX_SN];
/*...................................................................*/
  
/*... loop nas celulas*/
  
  for(nel=0;nel<numel;nel++){

/*... zerando vetores*/
    for(i=0;i<MAX_NUM_PONT;i++){
      lx[i] = 0.0;
    }
    
    for(i=0;i<MAX_SN;i++){
      s[i] = 0.0;
    }

    for(i=0;i<MAX_NUM_FACE+1;i++){
      lnFace[i]    = 0;
      lGeomType[i] = 0;     
    }

/*... loop na celula central*/    
    lnFace[maxViz]    = nFace[nel];
    lGeomType[maxViz] = geomType[nel];
    for(j=0;j<nen[nel];j++){
/*...*/
      no = MAT2D(nel,j,el,maxNo) - 1;
      for(k=0;k<ndm;k++){
        MAT3D(maxViz,j,k,lx,maxNo,ndm) = MAT2D(no,k,x,ndm);
      }
    }

/*... loop na celulas vizinhas*/    
      for(i=0;i<nFace[nel];i++){
        vizNel = MAT2D(nel,i,nelcon,maxViz) - 1;
        if( vizNel != -2) {
          lnFace[i]    = nFace[vizNel];
          lGeomType[i] = geomType[vizNel];
          for(j=0;j<nen[vizNel];j++){
            no = MAT2D(vizNel,j,el,maxNo) - 1;
            for(k=0;k<ndm;k++){
              MAT3D(i,j,k,lx,maxNo,ndm) 
              =  MAT2D(no,k,x,ndm);
            }
          }
        }
      }
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
    ty = geomType[nel];
    if(ty == 2 || ty == 3){
      sn(s,ty,nel); 
      cellGeom2D(lx       ,lnFace
                ,lGeomType,lCc
                ,lKsi     ,lmKsi 
                ,lEta     ,lmEta 
                ,lNormal  ,&lVolume
                ,lXm      ,lXmcc
                ,lDcca    ,lmKm 
                ,s
                ,maxNo    ,maxViz
                ,ndm      ,nel);
    }
/*...................................................................*/

/*... atribuido valores aos celula global*/
    gVolume[nel] = lVolume;
    for(i=0;i<ndm;i++){
      MAT2D(nel,i,gCc ,ndm) = MAT2D(maxViz,i,lCc  ,ndm);
    }
    for(i=0;i<nFace[nel];i++) {
      MAT2D(nel,i,gmKsi,maxViz)        = lmKsi[i];
      MAT2D(nel,i,gmEta,maxViz)        = lmEta[i];
      MAT2D(nel,i,gDcca,maxViz)        = lDcca[i];
      MAT2D(nel,i,gmKm ,maxViz)        = lmKm[i];
      for(j=0;j<ndm;j++){
        MAT3D(nel,i,j,gKsi   ,maxViz,ndm) = MAT2D(i,j,lKsi   ,ndm);
        MAT3D(nel,i,j,gEta,maxViz,ndm)    = MAT2D(i,j,lEta   ,ndm);
        MAT3D(nel,i,j,gNormal,maxViz,ndm) = MAT2D(i,j,lNormal,ndm);
        MAT3D(nel,i,j,gXm    ,maxViz,ndm) = MAT2D(i,j,lXm    ,ndm);
        MAT3D(nel,i,j,gXmcc  ,maxViz,ndm) = MAT2D(i,j,lXmcc  ,ndm);
      }
    }
/*...................................................................*/
  }
/*...................................................................*/
}
