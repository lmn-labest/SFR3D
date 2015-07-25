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
 * gfArea  -> indefinido                                             * 
 * gNormal -> indefinido                                             * 
 * gVolume -> indefinido                                             * 
 * gXm     -> indefinido                                             * 
 * gXmcc   -> indefinido                                             * 
 * vSkew   -> indefinido                                             * 
 * mvSkew  -> indefinifo                                             * 
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
 * gfArea  -> area da face                                           * 
 * gNormal -> vetores normais as faces das celulas                   * 
 * gVolume -> volumes das celulas                                    * 
 * gXm     -> pontos medios das faces das celulas                    * 
 * gXmcc   -> vetores que unem o centroide aos pontos medios das     * 
 *            faces 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
 * gDcca   -> menor distancia do centroide a faces desta celula      * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void pGeomForm(DOUBLE *restrict x      ,INT    *restrict el
              ,INT    *restrict nelcon ,short  *restrict nen    
              ,short  *restrict nFace  ,short  *restrict geomType
              ,DOUBLE *restrict gCc    ,DOUBLE *restrict gKsi   
              ,DOUBLE *restrict gmKsi  ,DOUBLE *restrict gEta   
              ,DOUBLE *restrict gfArea ,DOUBLE *restrict gNormal
              ,DOUBLE *restrict gVolume,DOUBLE *restrict gXm   
              ,DOUBLE *restrict gXmcc  
              ,DOUBLE *restrict gvSkew ,DOUBLE *restrict gmvSkew
              ,DOUBLE *restrict gDcca                 
              ,short maxNo             ,short maxViz
              ,short ndm               ,INT numel)
{
  INT nel,no,vizNel;
  short i,j,k;
/*... variavel local */
  DOUBLE lx[MAX_NUM_PONT];
  DOUBLE lCc[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM],lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM],lVolume;
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM],lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE],lvSkew[MAX_NUM_FACE*MAX_NDM];  
  short  lnFace[MAX_NUM_FACE+1],lGeomType[MAX_NUM_FACE+1],ty;
  short  lnEn[MAX_NUM_FACE+1];
  short  isnod[MAX_SN];
/*...................................................................*/
  
/*... loop nas celulas*/
  
  for(nel=0;nel<numel;nel++){

/*... zerando vetores*/
    for(i=0;i<MAX_NUM_PONT;i++){
      lx[i] = 0.0;
    }
    
    for(i=0;i<MAX_SN;i++){
      isnod[i] = 0.0;
    }

    for(i=0;i<MAX_NUM_FACE+1;i++){
      lnFace[i]    = 0;
      lnEn[i]      = 0;
      lGeomType[i] = 0;     
    }

/*... loop na celula central*/    
    lnFace[maxViz]    = nFace[nel];
    lnEn[maxViz]      = nen[nel];
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
          lnEn[i]      = nen[vizNel];
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
/*... triangulos e quadrilateros*/
    if(ty == TRIACELL || ty == QUADCELL){
      sn(isnod,ty,nel); 
      cellGeom2D(lx       ,lnFace
                ,lGeomType,lCc
                ,lKsi     ,lmKsi 
                ,lEta     ,lfArea
                ,lNormal  ,&lVolume
                ,lXm      ,lXmcc
                ,lDcca    
                ,lvSkew   ,lmvSkew 
                ,isnod
                ,maxNo    ,maxViz
                ,ndm      ,nel);
    }
/*...................................................................*/

/*... tetraedros e hexaedros*/
    else if(ty == TETRCELL || ty == HEXACELL){
      sn(isnod,ty,nel); 
      cellGeom3D(lx       ,lGeomType
                ,lnFace   ,lnEn
                ,lCc
                ,lKsi     ,lmKsi 
                ,lEta     ,lfArea
                ,lNormal  ,&lVolume
                ,lXm      ,lXmcc
                ,lDcca    
                ,lvSkew   ,lmvSkew 
                ,isnod
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
      MAT2D(nel,i,gmKsi   ,maxViz)        = lmKsi[i];
      MAT2D(nel,i,gfArea  ,maxViz)        = lfArea[i];
      MAT2D(nel,i,gDcca   ,maxViz)        = lDcca[i];
      MAT2D(nel,i,gmvSkew ,maxViz)        = lmvSkew[i];
      for(j=0;j<ndm;j++){
        MAT3D(nel,i,j,gKsi   ,maxViz,ndm) = MAT2D(i,j,lKsi   ,ndm);
        MAT3D(nel,i,j,gEta,maxViz,ndm)    = MAT2D(i,j,lEta   ,ndm);
        MAT3D(nel,i,j,gNormal,maxViz,ndm) = MAT2D(i,j,lNormal,ndm);
        MAT3D(nel,i,j,gXm    ,maxViz,ndm) = MAT2D(i,j,lXm    ,ndm);
        MAT3D(nel,i,j,gXmcc  ,maxViz,ndm) = MAT2D(i,j,lXmcc  ,ndm);
        MAT3D(nel,i,j,gvSkew ,maxViz,ndm) = MAT2D(i,j,lvSkew ,ndm);
      }
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * SYSTFOMDIF : calculo do sistema de equacoes para problemas        * 
 * difusao (Ax=b)                                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * el      -> conetividade dos celulas                               * 
 * nelcon  -> vizinhos dos elementos                                 * 
 * nen     -> numero de nos por celulas                              * 
 * nFace   -> numero de faces por celulas                            * 
 * calType -> tipo de calculo das celulas                            * 
 * geomType-> tipo geometrico das celulas                            * 
 * gKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        * 
 * gmKsi   -> modulo do vetor ksi                                    * 
 * gEta    -> vetores paralelos as faces das celulas                 * 
 * gfArea   -> modulo do vetor eta                                   * 
 * gNormal -> vetores normais as faces das celulas                   * 
 * gVolume -> volumes das celulas                                    * 
 * gXm     -> pontos medios das faces das celulas                    * 
 * gXmcc   -> vetores que unem o centroide aos pontos medios das     * 
 *            faces                                                  * 
 * gvSkew  -> vetor entre o ponto medio a intersecao que une os      * 
 *            centrois compartilhado nessa face                      * 
 * gmvSkew -> distacia entre o ponto medio a intersecao que une os   * 
 *            centrois compartilhado nessa face                      * 
 * gDcca   -> menor distancia do centroide a faces desta celula      * 
 * ia      -> ponteiro para as linhas da matriz esparsa              * 
 * ja      -> ponteiro para as colunas da matriz esparsa             * 
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 * al      -> matriz de coeficientes esparsa                         * 
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- triangular inferior               )       *
 * b       -> vetor de forcas                                        * 
 * id      -> numera das equacoes                                    * 
 * u0      -> solucao conhecida                                      * 
 * rCell   -> nao definido                                           * 
 * faceR   -> restricoes por elmento                                 * 
 * faceS   -> carga por elemento                                     * 
 * maxNo   -> numero de nos por celula maximo da malha               * 
 * maxViz  -> numero vizinhos por celula maximo da malha             * 
 * ndm     -> numero de dimensoes                                    * 
 * numel   -> numero de toral de celulas                             * 
 * ndf     -> graus de liberdade                                     * 
 * storage -> tecnica de armazenamento da matriz esparsa             * 
 * forces  -> mantagem no vetor de forcas                            * 
 * matrix  -> mantagem da matriz de coeficientes                     * 
 * calRcell-> calculo do residuo de celula                           * 
 * unsym   -> matiz nao simetrica                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * au,a,al   -> coeficiente da linha i     (matriz = true)           *
 * b         -> vetor de forca da linha i  (forces = true)           *
 * rCell   -> residuo por celula                                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void systFormDif(INT    *restrict el     ,INT    *restrict nelcon 
               ,short  *restrict nen     ,short  *restrict nFace
               ,short  *restrict geomType,DOUBLE *restrict prop 
               ,short  *restrict calType ,short  *restrict mat     
               ,DOUBLE *restrict gKsi    ,DOUBLE *restrict gmKsi 
               ,DOUBLE *restrict gEta    ,DOUBLE *restrict gfArea 
               ,DOUBLE *restrict gNormal ,DOUBLE *restrict gVolume
               ,DOUBLE *restrict gXm     ,DOUBLE *restrict gXmcc 
               ,DOUBLE *restrict gvSkew  ,DOUBLE *restrict gmvSkew 
               ,DOUBLE *restrict gDcca 
               ,INT    *restrict ia      ,INT    *restrict ja   
               ,DOUBLE *restrict ad      ,DOUBLE *restrict al
               ,DOUBLE *restrict b       ,INT    *restrict id
               ,short  *restrict faceR   ,DOUBLE *restrict faceS  
               ,DOUBLE *restrict u0      ,DOUBLE *restrict gradU0 
               ,DOUBLE *restrict rCell                                   
               ,short const maxNo        ,short const maxViz
               ,short const ndm          ,INT const numel
               ,short const ndf          ,short const storage
               ,bool  const forces       ,bool const matrix 
               ,bool const calRcell      ,bool  const  unsym) 
{
  INT nel,vizNel;
  short i,j;
/*... variavel local */
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM],lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM],lVolume[MAX_NUM_FACE+1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM],lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE],lvSkew[MAX_NUM_FACE*MAX_NDM];
  DOUBLE dum;  
  short  lGeomType[MAX_NUM_FACE+1];
  short  lib;
  short  lFaceR[MAX_NUM_FACE+1];
  DOUBLE lFaceS[(MAX_NUM_FACE+1)*MAX_NDF];
  DOUBLE lA[(MAX_NUM_FACE+1)*MAX_NDF],lB[MAX_NDF];
  DOUBLE lProp[(MAX_NUM_FACE+1)*MAXPROP];
  DOUBLE lu0[(MAX_NUM_FACE+1)*MAX_NDF];
  DOUBLE lGradU0[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lRcell[MAX_NDF];
  INT    lId[(MAX_NUM_FACE+1)*MAX_NDF],lViz[MAX_NUM_FACE];
  short  aux1,aux2,lMat;

/*... loop nas celulas*/
  aux2    = maxViz+1;
  for(nel=0;nel<numel;nel++){
/*...*/
    if(calRcell)
      for(j=0;j<ndf;j++)
        MAT2D(nel,j,rCell ,ndf)   = 0.e0;;
/*...*/
    aux1    = nFace[nel];
/*... elementos com equacoes*/
    if(MAT2D(nel,aux1,faceR ,aux2) != 1){

/*... zerando vetores*/
      for(j=0;j<(MAX_NUM_FACE+1)*MAX_NDF;j++){
        lId[j] = -1;
        lu0[j] = 0.e0;    
      }
      
      for(j=0;j<(MAX_NUM_FACE+1)*MAXPROP;j++)
        lProp[j] = 0.0e0;    
      
      for(j=0;j<MAX_NUM_FACE+1;j++){
        lGeomType[j] = 0;     
      }
      

/*... loop na celula central*/    
      lMat            = mat[nel]-1;
      lib             = calType[lMat];
      lVolume[aux1]   = gVolume[nel]; 
      lGeomType[aux1] = geomType[nel];
      lFaceR[aux1]    = MAT2D(nel,aux1,faceR ,aux2);
      
      for(j=0;j<ndf;j++){
        MAT2D(aux1,j,lu0   ,ndf) = MAT2D(nel,j,u0   ,ndf);
        MAT2D(aux1,j,lFaceS,ndf) = MAT3D(nel,aux1,j,faceS ,aux2,ndf);
        MAT2D(aux1,j,lId   ,ndf) = MAT2D(nel,j,id   ,ndf) - 1;
      }
      
      for(j=0;j<MAXPROP;j++)
        MAT2D(aux1,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
      
      for(j=0;j<ndm;j++)
        MAT2D(aux1,j,lGradU0,ndm) = MAT2D(nel,j,gradU0,ndm);
       

      for(i=0;i<aux1;i++){
        lmKsi[i]   = MAT2D(nel,i,gmKsi   ,maxViz);
        lfArea[i]   = MAT2D(nel,i,gfArea   ,maxViz);
        lDcca[i]   = MAT2D(nel,i,gDcca   ,maxViz);
        lmvSkew[i] = MAT2D(nel,i,gmvSkew ,maxViz);
        aux2       = (maxViz+1);
        lFaceR[i]  = MAT2D(nel,i,faceR ,aux2);
        for(j=0;j<ndf;j++){
          MAT2D(i,j,lFaceS,ndf) = MAT3D(nel,i,j,faceS ,aux2,ndf);
        }
        for(j=0;j<ndm;j++){
          MAT2D(i,j,lKsi   ,ndm) = MAT3D(nel,i,j,gKsi   ,maxViz,ndm);
          MAT2D(i,j,lEta   ,ndm) = MAT3D(nel,i,j,gEta   ,maxViz,ndm);
          MAT2D(i,j,lNormal,ndm) = MAT3D(nel,i,j,gNormal,maxViz,ndm);
          MAT2D(i,j,lXm    ,ndm) = MAT3D(nel,i,j,gXm    ,maxViz,ndm); 
          MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,gXmcc  ,maxViz,ndm);
          MAT2D(i,j,lvSkew ,ndm) = MAT3D(nel,i,j,gvSkew ,maxViz,ndm);
        }
      }


/*... loop na celulas vizinhas*/    
      for(i=0;i<aux1;i++){
        vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
        lViz[i] = vizNel;
        if( vizNel != -2) {
          lVolume[i]   = gVolume[vizNel]; 
          lGeomType[i] = geomType[vizNel];
          lMat = mat[vizNel]-1;
          for(j=0;j<ndf;j++){
            MAT2D(i,j,lu0 ,ndf)   = MAT2D(vizNel,j,u0   ,ndf);
            MAT2D(i,j,lId ,ndf)   = MAT2D(vizNel,j,id   ,ndf) - 1;
          }
          for(j=0;j<ndm;j++)
            MAT2D(i,j,lGradU0,ndm)   = MAT2D(vizNel,j,gradU0,ndm);
          for(j=0;j<MAXPROP;j++)
            MAT2D(i,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
          }
        }
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
      cellLibDif(lGeomType ,lProp 
                ,lViz      ,lId           
                ,lKsi      ,lmKsi
                ,lEta      ,lfArea 
                ,lNormal   ,lVolume
                ,lXm       ,lXmcc
                ,lDcca     
                ,lvSkew    ,lmvSkew
                ,lA        ,lB
                ,lRcell
                ,lFaceR    ,lFaceS
                ,lu0       ,lGradU0     
                ,nen[nel]  ,nFace[nel] 
                ,ndm       ,lib   
                ,nel);
/*...................................................................*/

/*... residuo da celula*/
      if(calRcell)
        for(j=0;j<ndf;j++)
          MAT2D(nel,j,rCell ,ndf)   = lRcell[j];
/*...................................................................*/
      
/*...*/
      assbly(ia          ,ja
            ,&dum        ,ad 
            ,al          ,b  
            ,lId 
            ,lA          ,lB 
            ,nFace[nel]  ,ndf 
            ,storage     ,forces
            ,matrix      ,unsym); 
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * CELLPLOAD : Carregamento prescrito por centro da celula           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * faceR     -> restricoes por elmento                               * 
 * faceS     -> carga por elemento                                   * 
 * volume    -> volume das celulas                                   * 
 * id        -> numera das equacoes                                  * 
 * u         -> solucao                                              * 
 * f         -> vetor de forcas                                      * 
 * numel     -> numero de elementos                                  * 
 * ndf       -> graus de liberade                                    * 
 * maxViz    -> numero maximo de vizinhos                            * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * u0        -> atualizado com a restricao                           * 
 * f         -> atualizado com as forcas de volume                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellPload(short  *restrict faceR ,DOUBLE *restrict faceS
              ,DOUBLE *restrict volume,INT *restrict id 
              ,DOUBLE *restrict u     ,DOUBLE *restrict f
              ,INT const numel        ,short const ndf
              ,short const maxViz)
{
  INT nel,lNeq;
  short carg,j;
  short col = maxViz + 1;

/*...*/  
  for(nel = 0; nel < numel;nel++){
    carg = MAT2D(nel,maxViz,faceR,col);
/*... variavel prescrita no dominio*/
    if( carg == VPES){
      for(j = 0; j< ndf;j++)
        MAT2D(nel,j,u,ndf) = MAT3D(nel,maxViz,j,faceS,col,ndf);
    }
/*...................................................................*/

/*... carga */
    else if( carg ==  CARGCONST){
      for(j = 0; j< ndf;j++){
        lNeq = MAT2D(nel,j,id,ndf) - 1;
        if( lNeq > -1)
          MAT2D(lNeq,j,f,ndf) 
          = volume[nel]*MAT3D(nel,maxViz,j,faceS,col,ndf);
      }
    }
/*...................................................................*/
  }
/*...................................................................*/
}  
/*********************************************************************/ 

/********************************************************************* 
 * UPDATECELLVALUE:atualizacao dos valores das variaveis das celulas *
 * com os valores das respectivas equacoes                           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * u       -> variavel nas celulas                                   * 
 * x       -> solucao do sistema                                     * 
 * id      -> numera das equacoes                                    * 
 * numel   -> numero de elementos                                    * 
 * ndf     -> graus de liberdade                                     * 
 * fAdd    -> true add false sobreescreve                            * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * u      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 void updateCellValue(DOUBLE *restrict u,DOUBLE *restrict x
                 ,INT *restrict id 
                 ,INT const numel   ,short const ndf
                 ,bool const fAdd)
{
  INT nel,lNeq;
  short jNdf;
  if(fAdd)    
    for(nel=0;nel<numel;nel++){
      for(jNdf = 0;jNdf<ndf;jNdf++){ 
        lNeq = MAT2D(nel,jNdf,id,ndf) - 1;
        if( lNeq > -1)
          MAT2D(nel,jNdf,u,ndf) += MAT2D(lNeq,jNdf,x,ndf);
      }
    }
  else
    for(nel=0;nel<numel;nel++){
      for(jNdf = 0;jNdf<ndf;jNdf++){ 
        lNeq = MAT2D(nel,jNdf,id,ndf) - 1;
        if( lNeq > -1)
          MAT2D(nel,jNdf,u,ndf) = MAT2D(lNeq,jNdf,x,ndf);
      }
    }
}
/*********************************************************************/

 
/********************************************************************* 
 * INTERCELLNODE: interpolacao dos valores das celulas para o no da  *
 * malha                                                             *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m       -> vetor de memoria principal                             * 
 * noU     -> nao definido                                           * 
 * elU     -> valores nas celulas                                    * 
 * el      -> conectividades das celulas                             * 
 * geomType-> tipo geometrico das celulas                            * 
 * cc      -> centroide das celulas                                  * 
 * x       -> coordenadas                                            * 
 * nen     -> numero de nos por celulas                              * 
 * faceR   -> restricoes por elmento                                 * 
 * faceS   -> carga por elemento                                     * 
 * numel   -> numero de elementos                                    * 
 * nnode   -> numero de nos                                          * 
 * maxNo   -> numero de nos por celula maximo da malha               * 
 * ndf     -> graus de liberdade                                     * 
 * ndm     -> numero de dimensao                                     * 
 * type    -> tipo de interpolacao                                   * 
 *            1 - media simples                                      * 
 *            2 - media ponderada                                    * 
 * fBc     -> forca condicao de controno conhecida                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * u      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 void interCellNode(Memoria *m             
                   ,DOUBLE *restrict noU   ,DOUBLE *restrict elU
                   ,INT *restrict el       ,short  *restrict geomType 
                   ,DOUBLE *restrict cc    ,DOUBLE *restrict x 
                   ,short *restrict nen    ,short *restrict nFace
                   ,short  *restrict faceR ,DOUBLE *restrict faceS  
                   ,INT const numel        ,INT const nnode
                   ,short const maxNo      ,short const maxViz     
                   ,short const ndf        ,short const ndm
                   ,bool const fBc         ,short const type)

{
  int *md=NULL;
  DOUBLE *mdf=NULL;
  bool *flag=NULL;
  DOUBLE dist,dx;
  short i,j,k,n,nodeFace,aux=maxViz+1;
  INT nel,no1,no[4];
  short  isNod[MAX_SN],ty;
  
  switch(type){
/*... media simple*/
    case 1:
/*...*/
      HccaAlloc(int,m,md,nnode,"md",false);
      zero(md,nnode,"int");
      zero(noU,ndf*nnode,DOUBLEC);
/*...................................................................*/

/*...*/
      for(nel = 0; nel < numel; nel++){
        for(j = 0; j < nen[nel];j++){
          no1 = MAT2D(nel,j,el,maxNo) - 1;
          for(k = 0; k   < ndf;k++)
            MAT2D(no1,k,noU,ndf) += MAT2D(nel,k,elU,ndf);
          md[no1]++;
        }
      }
/*...................................................................*/

/*...*/
      for(no1 = 0; no1 < nnode; no1++)
        for(k = 0; k < ndf; k++)
          MAT2D(no1,k,noU,ndf) /= md[no1];
/*...................................................................*/
          
/*...*/
      HccaDealloc(m,md,"md",false);
/*...................................................................*/
    break;
/*...................................................................*/

/*... media ponderada*/
    case 2:
/*...*/
      HccaAlloc(DOUBLE,m,mdf,nnode,"mdf",false);
      zero(mdf,nnode,"double");
      zero(noU,ndf*nnode,DOUBLEC);
/*...................................................................*/

/*...*/
      for(nel = 0; nel < numel; nel++){
        for(j = 0; j < nen[nel];j++){
          no1 = MAT2D(nel,j,el,maxNo) - 1;
          dist = 0.e0;
          for(k = 0; k   < ndm;k++){
            dx = MAT2D(no1,k,x,ndm) - MAT2D(nel,k,cc,ndm);
            dist += dx*dx;
          }
          dist = 1.e0/sqrt(dist);
          for(k = 0; k   < ndf;k++)
            MAT2D(no1,k,noU,ndf) += MAT2D(nel,k,elU,ndf)*dist;
          mdf[no1]+=dist;
        }
      }
/*...................................................................*/

/*...*/
      for(no1 = 0; no1 < nnode; no1++)
        for(k = 0; k < ndf; k++)
          MAT2D(no1,k,noU,ndf) /= mdf[no1];
/*...................................................................*/
          
/*...*/
      HccaDealloc(m,mdf,"mdf",false);
/*...................................................................*/
    break;
/*...................................................................*/

/*... media simple*/
    default:
      ERRO_OP(__FILE__,__func__,type);
    break;
  }
/*...................................................................*/

/*...*/
  if(fBc){
    HccaAlloc(int ,m,md  ,nnode,"md",false);
    HccaAlloc(bool,m,flag,nnode,"flag",false);
    zero(md  ,nnode,"int");
    zero(flag,nnode,"char");
    for(nel = 0; nel < numel; nel++)
      for(i = 0; i < nFace[nel]; i++)
        if(MAT2D(nel,i,faceR,aux)){
          ty = geomType[nel];
          nodeFace =  sn(isNod,ty,nel); 
          for(n=0;n<nodeFace ;n++){
            no[n] = MAT2D(i,n,isNod,nodeFace );
            no[n] = MAT2D(nel,no[n],el,maxNo) - 1;
            if(flag[no[n]] == false){
              flag[no[n]] = true;
              for(k = 0; k   < ndf;k++)
                MAT2D(no[n],k,noU,ndf) = 0.e0;
            }
            for(k = 0; k   < ndf;k++) 
              MAT2D(no[n],k,noU,ndf) += MAT3D(nel,i,k,faceS,aux,ndf);
            md[no[n]]++;
          }
        }
/*...................................................................*/

/*...*/
    for(no1 = 0; no1 < nnode; no1++)
      if(flag[no1])
        for(k = 0; k < ndf; k++)
          MAT2D(no1,k,noU,ndf) /= md[no1];
/*...................................................................*/
    HccaDealloc(m,flag,"flag",false);
    HccaDealloc(m,md  ,"md"  ,false);
  }
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * RCGRADU: calculo do gradiente de um campo escalar ou vetorial     * 
 * conhecido.                                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m       -> vetor de memoria principal                             * 
 * el      -> conetividade dos celulas                               * 
 * cc      -> centroide das celulas                                  * 
 * x       -> cordenadas dos pontos                                  * 
 * nelcon  -> vizinhos dos elementos                                 * 
 * nen     -> numero de nos por celulas                              * 
 * nFace   -> numero de faces por celulas                            * 
 * geomType-> tipo geometrico das celulas                            * 
 * lSquare -> matriz para a reconstrucao least Square                * 
 * lSquareR-> fatoracao R (RCLSQUAREQR)                              * 
 * gKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        * 
 * gmKsi   -> modulo do vetor ksi                                    * 
 * gEta    -> vetores paralelos as faces das celulas                 * 
 * gfArea  -> areas das faces                                        * 
 * gNormal -> vetores normais as faces das celulas                   * 
 * gVolume -> volumes das celulas                                    * 
 * gvSkew  -> vetor entre o ponto medio a intersecao que une os      * 
 *            centrois compartilhado nessa face                      * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * faceR   -> restricoes por elmento                                 * 
 * faceS   -> carga por elemento                                     * 
 * u       -> solucao conhecida por celula (atualizado)              * 
 * gradU   -> gradiente da solucao         (desatualizado)           * 
 * nU      -> solucao conhecida por no     (desatualizado)           * 
 * lib     -> tipo de reconstrucao de gradiente                      * 
 * maxNo   -> numero de nos por celula maximo da malha               * 
 * maxViz  -> numero vizinhos por celula maximo da malha             * 
 * ndm     -> numero de dimensoes                                    * 
 * numel   -> numero de toral de celulas                             * 
 * ndf     -> graus de liberdade                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*
 *         |gradU1|   | dU1/dx1 dU1/dx2 dU1/dx3 |                    *
 * gradU = | ...  | = |         ...             |                    *
 *         |gradUn|   | dUn/dx1 dUn/dx2 dUn/dx3 |                    *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void rcGradU(Memoria *m
            ,INT    *restrict el      ,INT    *restrict nelcon 
            ,DOUBLE *restrict cc      ,DOUBLE *restrict x     
            ,short  *restrict nen     ,short  *restrict nFace
            ,short  *restrict geomType
            ,DOUBLE *restrict lSquare ,DOUBLE *restrict lSquareR            
            ,DOUBLE *restrict gKsi    ,DOUBLE *restrict gmKsi 
            ,DOUBLE *restrict gEta    ,DOUBLE *restrict gfArea 
            ,DOUBLE *restrict gNormal ,DOUBLE *restrict gVolume
            ,DOUBLE *restrict gvSkew  ,DOUBLE *restrict gXmcc  
            ,short  *restrict faceR   ,DOUBLE *restrict faceS  
            ,DOUBLE *restrict u       ,DOUBLE *restrict gradU              
            ,DOUBLE *restrict nU      ,short const lib 
            ,short const maxNo        ,short const maxViz
            ,short const ndf          ,short const ndm
            ,INT const numel          ,INT const nNode)
{
  INT nel,no,vizNel;
  short i,j;
/*... variavel local */
  DOUBLE lKsi[MAX_NUM_FACE*3],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*3],lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*3],lVolume[MAX_NUM_FACE+1];
  DOUBLE lvSkew[MAX_NUM_FACE*3];
  DOUBLE lXmcc[MAX_NUM_FACE*3];
  short  lFaceR[MAX_NUM_FACE+1];
  DOUBLE lFaceS[(MAX_NUM_FACE+1)*MAX_NDF];
  DOUBLE lu[(MAX_NUM_FACE+1)*MAX_NDF];
  DOUBLE lnU[(MAX_NUM_NODE)*MAX_NDF];
  DOUBLE lGradU[MAX_NDM*MAX_NDF];
  DOUBLE lLsquare[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lLsquareR[2*MAX_NDM];
  INT    lViz[MAX_NUM_FACE];
  short  aux1,aux2;
  short  isNod[MAX_SN],ty;

/*... reconstrucao de gradiente Green-Gauss nodal*/
  if(lib ==  RCGRADGAUSSN){
    interCellNode(m
                 ,nU        ,u
                 ,el        ,geomType
                 ,cc        ,x
                 ,nen       ,nFace
                 ,faceR     ,faceS   
                 ,numel     ,nNode    
                 ,maxNo     ,maxViz
                 ,ndf       ,ndm   
                 ,true      ,2);
  }
/*.....................................................................*/


/*... */
  aux2    = maxViz+1;
  for(nel=0;nel<numel;nel++){
    aux1    = nFace[nel];

/*... loop na celula central*/    
    lVolume[aux1]   = gVolume[nel]; 
    lFaceR[aux1]    = MAT2D(nel,aux1,faceR ,aux2);
      
    for(i=0;i<ndf;i++){
      MAT2D(aux1,i,lu    ,ndf) = MAT2D(nel,i,u    ,ndf);
      MAT2D(aux1,i,lFaceS,ndf) = MAT3D(nel,aux1,i,faceS ,aux2,ndf);
    }

/*... valor da funcao nodal nodias*/    
    if(lib ==  RCGRADGAUSSN){
      for(i=0;i<nen[nel];i++){
        no = MAT2D(nel,i,el,maxNo)-1;
        for(j=0;j<ndf;j++)
          MAT2D(i,j,lnU   ,ndf) = MAT2D(no,j,nU    ,ndf);
      }
    }

/*... leastSquare*/
    if(lib ==  RCLSQUARE)
      for(i=0;i<ndm;i++)
        for(j=0;j<aux1;j++)
          MAT2D(i,j,lLsquare,aux1) 
          = MAT3D(nel,i,j,lSquare,ndm,maxViz); 
/*...................................................................*/

/*... leastSquare-QR*/
    if(lib ==  RCLSQUAREQR){
      for(i=0;i<ndm;i++)
        for(j=0;j<aux1;j++)
          MAT2D(i,j,lLsquare,aux1) 
          = MAT3D(nel,i,j,lSquare,ndm,maxViz);
/*... R*/
       if(ndm == 2){ 
         lLsquareR[0] = MAT2D(nel,0,lSquareR,3); 
         lLsquareR[1] = MAT2D(nel,1,lSquareR,3); 
         lLsquareR[2] = MAT2D(nel,2,lSquareR,3); 
       }
       else if(ndm == 3){ 
         lLsquareR[0] = MAT2D(nel,0,lSquareR,6); 
         lLsquareR[1] = MAT2D(nel,1,lSquareR,6); 
         lLsquareR[2] = MAT2D(nel,2,lSquareR,6); 
         lLsquareR[3] = MAT2D(nel,3,lSquareR,6); 
         lLsquareR[4] = MAT2D(nel,4,lSquareR,6); 
         lLsquareR[5] = MAT2D(nel,5,lSquareR,6); 
       }
    }
/*...................................................................*/



/*...*/
    for(i=0;i<ndf;i++)
      for(j=0;j<ndm;j++)
        MAT2D(i,j,lGradU   ,ndf) =  MAT3D(nel,i,j,gradU,ndf,ndm);
/*...................................................................*/
      

/*...*/      
    for(i=0;i<aux1;i++){
      lmKsi[i]  = MAT2D(nel,i,gmKsi  ,maxViz);
      lfArea[i] = MAT2D(nel,i,gfArea ,maxViz);
      lFaceR[i] = MAT2D(nel,i,faceR  ,aux2);
      for(j=0;j<ndf;j++){
        MAT2D(i,j,lFaceS,ndf) = MAT3D(nel,i,j,faceS ,aux2,ndf);
      }
      for(j=0;j<ndm;j++){
        MAT2D(i,j,lKsi   ,ndm) = MAT3D(nel,i,j,gKsi   ,maxViz,ndm);
        MAT2D(i,j,lEta   ,ndm) = MAT3D(nel,i,j,gEta   ,maxViz,ndm);
        MAT2D(i,j,lNormal,ndm) = MAT3D(nel,i,j,gNormal,maxViz,ndm);
        MAT2D(i,j,lvSkew ,ndm) = MAT3D(nel,i,j,gvSkew ,maxViz,ndm);
        MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,gXmcc  ,maxViz,ndm);
      }
/*... loop na celulas vizinhas*/    
      vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
      lViz[i] = vizNel;
      if( vizNel != -2) {
        lVolume[i] = gVolume[vizNel]; 
        for(j=0;j<ndf;j++)
          MAT2D(i,j,lu ,ndf)   = MAT2D(vizNel,j,u   ,ndf);
      }

    }
/*...................................................................*/
    

/*... chamando a biblioteca de celulas*/
    ty = geomType[nel];
    sn(isNod,ty,nel); 
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
    cellLibRcGrad(lViz      
                 ,lLsquare  ,lLsquareR
                 ,lKsi      ,lmKsi
                 ,lEta      ,lfArea 
                 ,lNormal   ,lVolume
                 ,lvSkew    ,lXmcc 
                 ,lFaceR    ,lFaceS
                 ,lu        ,lGradU
                 ,lnU       ,ty     
                 ,nFace[nel],ndm 
                 ,lib       ,ndf
                 ,isNod     ,nel);  
/*...................................................................*/
  
/*...*/
    for(i=0;i<ndf;i++)
      for(j=0;j<ndm;j++)
        MAT3D(nel,i,j,gradU,ndf,ndm)  = MAT2D(i,j,lGradU   ,ndf);
/*...................................................................*/
  }
}
/*********************************************************************/

/********************************************************************* 
 * RCLEASTSQUARE : calcula a matriz dos minimos quadrados            *
 * para reconstrucao de gradiente                                    *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * gKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        * 
 * gmKsi   -> modulo do vetor ksi                                    * 
 * lSquare   -> nao definido                                         * 
 * lSquareR  -> nao definido                                         * 
 * numel     -> numero de elementos                                  * 
 * maxViz    -> numero maximo de vizinhos                            * 
 * nFace     -> numero vizinhos por celula maximo da malha           * 
 * type      -> tecnica de resolucao                                 * 
 * ndm       -> numero de dimensoes                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lSquare   -> matriz para a reconstrucao least Square              *
 * lSquareR  -> fatoracao R (RCLSQUAREQR)                            * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void rcLeastSquare(DOUBLE *restrict gKsi    ,DOUBLE *restrict gmKsi
                  ,DOUBLE *restrict lSquare ,DOUBLE *restrict lSquareR
                  ,short *restrict nFace       
                  ,INT const numel          ,short const maxViz
                  ,short const type         ,short const ndm){

  INT nEl;
  DOUBLE lLsquare[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lLsquareR[2*MAX_NDM];
  DOUBLE lKsi[MAX_NUM_FACE*3],lmKsi[MAX_NUM_FACE];
  short lnFace,i,j;

  
  for(nEl=0;nEl<numel;nEl++){
    lnFace  = nFace[nEl];


    for(i=0;i<lnFace;i++){
      lmKsi[i] = MAT2D(nEl,i,gmKsi ,maxViz);
      for(j=0;j<ndm;j++){
        MAT2D(i,j,lKsi   ,ndm) = MAT3D(nEl,i,j,gKsi   ,maxViz,ndm);
      }
    }
    
    leastSquareMatrix(lKsi    ,lmKsi
                     ,lLsquare,lLsquareR
                     ,type
                     ,lnFace  ,ndm);

    for(i=0;i<ndm;i++)
      for(j=0;j<lnFace;j++){ 
        MAT3D(nEl,i,j,lSquare,ndm,maxViz) 
         = MAT2D(i,j,lLsquare,lnFace);
      }

    if(type == RCLSQUAREQR){
      if(ndm == 2){
        MAT2D(nEl,0,lSquareR,3) = lLsquareR[0];
        MAT2D(nEl,1,lSquareR,3) = lLsquareR[1];
        MAT2D(nEl,2,lSquareR,3) = lLsquareR[2];
      }
      else if(ndm == 3){
        MAT2D(nEl,0,lSquareR,6) = lLsquareR[0];
        MAT2D(nEl,1,lSquareR,6) = lLsquareR[1];
        MAT2D(nEl,2,lSquareR,6) = lLsquareR[2];
        MAT2D(nEl,3,lSquareR,6) = lLsquareR[3];
        MAT2D(nEl,4,lSquareR,6) = lLsquareR[4];
        MAT2D(nEl,5,lSquareR,6) = lLsquareR[5];
      }
    } 
  }
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * CONVTEMPFORKELVIN : conversao entre C e Kelvin                    * 
 * conhecido.                                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * u       -> temperatura em C/K                                     * 
 * n       -> tamanho do arranjo                                     * 
 * fKevin  -> true/false                                             * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*
 * u       -> temperatura em K/C                                     *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
  void convTempForKelvin(DOUBLE *restrict u,INT const n
                        ,bool const fKelvin){
    int i;
    if(fKelvin)
      for(i=0;i<n;i++)
        u[i] += KELVINCONV; 
    else
      for(i=0;i<n;i++)
        u[i] += -KELVINCONV; 
    
  }
/*********************************************************************/

/********************************************************************* 
 * MESHQUALITY: calulo das propriedades da malha                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * gVolume -> volumes das celulas                                    * 
 * gKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        * 
 * gNormal -> vetores normais as faces das celulas                   * 
 * gmvSkew -> distacia entre o ponto medio a intersecao que une os   * 
 *            centrois compartilhado nessa face da celula central    * 
 * maxViz  -> numero vizinhos por celula maximo da malha             * 
 * ndm     -> numero de dimensoes                                    * 
 * numel   -> numero de toral de celulas                             * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * volume total                                                      *
 * volar maximo da nao-ortogonalidade em graus                       *
 * volar medio da nao-ortogonalidade em graus                        *
 * skewness medio                                                    *
 * skewness maximo                                                   *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void meshQuality(MeshQuality *mq
                ,short  *restrict nFace   ,DOUBLE *restrict volume
                ,DOUBLE *restrict gKsi    ,DOUBLE *restrict gNormal
                ,DOUBLE *restrict gmvSkew
                ,short const maxViz      ,short const ndm
                ,INT const numel)        
{
  INT nEl,k=0;
  DOUBLE volumeTotal = 0.e0,nk,nkMin=1.e0,nkMed=0.e0;
  DOUBLE skewMax=0.0e0, skewMed=0.0e0,teta;
  short nf,j;

  for(nEl=0;nEl<numel;nEl++){
/*... Volume total da malha*/
    volumeTotal += volume[nEl];
/*... nao-ortoganilidade*/
    for(nf=0;nf<nFace[nEl];nf++){
      nk = 0.e0;
/*... k * normal*/
      for(j=0;j<ndm;j++)
        nk += MAT3D(nEl,nf,j,gKsi   ,maxViz,ndm)
             *MAT3D(nEl,nf,j,gNormal,maxViz,ndm);
/*...................................................................*/ 
      
      nkMed += nk;
      nkMin = min(nkMin,nk); 
      k++;    
    }
/*...................................................................*/ 

/*... sKew*/  
    skewMed += gmvSkew[nEl];
    skewMax  = max(skewMax,gmvSkew[nEl]);
  }

  skewMed    /= numel;
  nkMed      /= k;

  teta           = acos(nkMed);
  teta           = radToDeg(teta);
  mq->nonOrthMed = teta; 

  teta           = acos(nkMin);
  teta           = radToDeg(teta);
  mq->nonOrthMax = teta; 
  
  mq->volume  = volumeTotal; 
  mq->skewMed = skewMed; 
  mq->skewMax = skewMax; 

}
/*********************************************************************/
