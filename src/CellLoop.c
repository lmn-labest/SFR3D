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
void pGeomForm(DOUBLE *restrict x      ,INT    *restrict el
              ,INT    *restrict nelcon ,short  *restrict nen    
              ,short  *restrict nFace  ,short  *restrict geomType
              ,DOUBLE *restrict gCc    ,DOUBLE *restrict gKsi   
              ,DOUBLE *restrict gmKsi  ,DOUBLE *restrict gEta   
              ,DOUBLE *restrict gmEta  ,DOUBLE *restrict gNormal
              ,DOUBLE *restrict gVolume,DOUBLE *restrict gXm   
              ,DOUBLE *restrict gXmcc  ,DOUBLE *restrict gmKm 
              ,DOUBLE *restrict gDcca                 
              ,short maxNo             ,short maxViz
              ,short ndm               ,INT numel)
{
  INT nel,no,vizNel;
  short i,j,k;
/*... variavel local */
  DOUBLE lx[MAX_NUM_PONT];
  DOUBLE lCc[(MAX_NUM_FACE+1)*3];
  DOUBLE lKsi[MAX_NUM_FACE*3],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*3],lmEta[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*3],lVolume;
  DOUBLE lXm[MAX_NUM_FACE*3],lXmcc[MAX_NUM_FACE*3];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmKm[MAX_NUM_FACE];  
  short  lnFace[MAX_NUM_FACE+1],lGeomType[MAX_NUM_FACE+1],ty;
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
      sn(isnod,ty,nel); 
      cellGeom2D(lx       ,lnFace
                ,lGeomType,lCc
                ,lKsi     ,lmKsi 
                ,lEta     ,lmEta 
                ,lNormal  ,&lVolume
                ,lXm      ,lXmcc
                ,lDcca    ,lmKm 
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
/*********************************************************************/ 

/********************************************************************* 
 * SYSTFOM : calculo do sistema de equacoes (Ax=b)                   * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x       -> cordenadas dos pontos                                  * 
 * el      -> conetividade dos celulas                               * 
 * nelcon  -> vizinhos dos elementos                                 * 
 * nen     -> numero de nos por celulas                              * 
 * nFace   -> numero de faces por celulas                            * 
 * calType -> tipo de calculo das celulas                            * 
 * geomType-> tipo geometrico das celulas                            * 
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
 * unsym   -> matiz nao simetrica                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * au,a,al   -> coeficiente da linha i     (matriz = true)           *
 * b         -> vetor de forca da linha i  (forces = true)           *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void systForm(DOUBLE *restrict x      ,INT    *restrict el
             ,INT    *restrict nelcon ,short  *restrict nen    
             ,short  *restrict nFace  ,short  *restrict geomType
             ,DOUBLE *restrict prop   ,short  *restrict calType
             ,short  *restrict mat    ,DOUBLE *restrict gCc 
             ,DOUBLE *restrict gKsi   ,DOUBLE *restrict gmKsi 
             ,DOUBLE *restrict gEta   ,DOUBLE *restrict gmEta 
             ,DOUBLE *restrict gNormal,DOUBLE *restrict gVolume
             ,DOUBLE *restrict gXm    ,DOUBLE *restrict gXmcc 
             ,DOUBLE *restrict gmKm   ,DOUBLE *restrict gDcca 
             ,INT    *restrict ia     ,INT    *restrict ja   
             ,DOUBLE *restrict ad     ,DOUBLE *restrict al
             ,DOUBLE *restrict b      ,INT    *restrict id
             ,short  *restrict faceR  ,DOUBLE *restrict faceS  
             ,DOUBLE *restrict u0                                           
             ,short const maxNo       ,short const maxViz
             ,short const ndm         ,INT const numel
             ,short const ndf         ,short const storage
             ,bool  const forces      ,bool const matrix 
             ,bool  const  unsym) 
{
  INT nel,no,vizNel;
  short i,j,k;
/*... variavel local */
  DOUBLE lx[MAX_NUM_PONT];
  DOUBLE lCc[(MAX_NUM_FACE+1)*3];
  DOUBLE lKsi[MAX_NUM_FACE*3],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*3],lmEta[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*3],lVolume[MAX_NUM_FACE+1];
  DOUBLE lXm[MAX_NUM_FACE*3],lXmcc[MAX_NUM_FACE*3];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmKm[MAX_NUM_FACE];
  DOUBLE dum;  
  short  lGeomType[MAX_NUM_FACE+1];
  short  lib;
  short  lFaceR[MAX_NUM_FACE+1];
  DOUBLE lFaceS[(MAX_NUM_FACE+1)*MAX_NDF];
  DOUBLE lA[(MAX_NUM_FACE+1)*MAX_NDF],lB[MAX_NDF];
  DOUBLE lProp[(MAX_NUM_FACE+1)*MAXPROP];
  DOUBLE lu0[(MAX_NUM_FACE+1)*MAX_NDF];
  INT    lId[(MAX_NUM_FACE+1)*MAX_NDF],lViz[MAX_NUM_FACE];
  short  aux1,aux2,lMat;
/*... loop nas celulas*/
  
  for(nel=0;nel<numel;nel++){
    aux2    = maxViz+1;
    aux1    = nFace[nel];
/*... elementos com equacoes*/
    if(MAT2D(nel,aux1,faceR ,aux2) != 1){

/*... zerando vetores*/
      for(i=0;i<MAX_NUM_PONT;i++){
        lx[i] = 0.0e0;
      }
      
      for(i=0;i<(MAX_NUM_FACE+1)*MAX_NDF;i++)
        lId[i] = -1;    
      
      for(i=0;i<(MAX_NUM_FACE+1)*MAXPROP;i++)
        lProp[i] = 0.0e0;    
      
      for(i=0;i<MAX_NUM_FACE+1;i++){
        lGeomType[i] = 0;     
      }

/*... loop na celula central*/    
      lMat    = mat[nel]-1;
      lib     = calType[lMat];
      lVolume[aux1]   = gVolume[nel]; 
      lGeomType[aux1] = geomType[nel];
      lFaceR[aux1] = MAT2D(nel,aux1,faceR ,aux2);
      for(i=0;i<ndm;i++)
        MAT2D(aux1,i,lCc ,ndm) = MAT2D(nel,i,gCc  ,ndm);
      
      for(i=0;i<ndf;i++){
        MAT2D(aux1,i,lu0   ,ndf) = MAT2D(nel,i,u0   ,ndf);
        MAT2D(aux1,i,lFaceS,ndf) = MAT3D(nel,aux1,i,faceS ,aux2,ndf);
        MAT2D(aux1,i,lId   ,ndf) = MAT2D(nel,i,id   ,ndf) - 1;
      }
      
      for(i=0;i<MAXPROP;i++)
        MAT2D(aux1,i,lProp,MAXPROP) = MAT2D(lMat,i,prop,MAXPROP);
       

      for(i=0;i<nFace[nel];i++){
        lmKsi[i] = MAT2D(nel,i,gmKsi ,maxViz);
        lmEta[i] = MAT2D(nel,i,gmEta ,maxViz);
        lDcca[i] = MAT2D(nel,i,gDcca ,maxViz);
        lmKm[i]  = MAT2D(nel,i,gmKm  ,maxViz);
        aux2     = (maxViz+1);
        lFaceR[i] = MAT2D(nel,i,faceR ,aux2);
        for(j=0;j<ndf;j++){
          MAT2D(i,j,lFaceS,ndf) = MAT3D(nel,i,j,faceS ,aux2,ndf);
        }
        for(j=0;j<ndm;j++){
          MAT2D(i,j,lKsi   ,ndm) = MAT3D(nel,i,j,gKsi   ,maxViz,ndm);
          MAT2D(i,j,lEta   ,ndm) = MAT3D(nel,i,j,gEta   ,maxViz,ndm);
          MAT2D(i,j,lNormal,ndm) = MAT3D(nel,i,j,gNormal,maxViz,ndm);
          MAT2D(i,j,lXm    ,ndm) = MAT3D(nel,i,j,gXm    ,maxViz,ndm); 
          MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,gXmcc  ,maxViz,ndm);
        }
      }

      for(j=0;j<nen[nel];j++){
/*...*/
        no = MAT2D(nel,j,el,maxNo) - 1;
        for(k=0;k<ndm;k++){
          MAT3D(maxViz,j,k,lx,maxNo,ndm) = MAT2D(no,k,x,ndm);
        }
      }

/*... loop na celulas vizinhas*/    
      for(i=0;i<nFace[nel];i++){
        vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
        lViz[i] = vizNel;
        if( vizNel != -2) {
          lVolume[i]   = gVolume[vizNel]; 
          lGeomType[i] = geomType[vizNel];
          aux1    = nFace[nel] +1;
          aux2    = maxViz     +1;
          lMat = mat[vizNel]-1;
          for(j=0;j<ndf;j++){
            MAT2D(i,j,lu0 ,ndf)   = MAT2D(vizNel,j,u0   ,ndf);
            MAT2D(i,j,lId ,ndf)   = MAT2D(vizNel,j,id   ,ndf) - 1;
          }
          for(j=0;j<MAXPROP;j++)
            MAT2D(i,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
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
      cellLib(lx       
             ,lGeomType ,lProp
             ,lViz      ,lCc           
             ,lKsi      ,lmKsi
             ,lEta      ,lmEta 
             ,lNormal   ,lVolume
             ,lXm       ,lXmcc
             ,lDcca     ,lmKm
             ,lA        ,lB
             ,lFaceR    ,lFaceS
             ,lu0       ,nen[nel] 
             ,nFace[nel],ndm
             ,lib       ,nel);
/*...................................................................*/
/*   printf("%d %lf %lf %lf %lf %lf %lf\n"
          ,nel+1,lA[0],lA[1],lA[2],lA[3],lA[4],lB[0]);
     printf("%d %d %d %d %d %d\n"
          ,nel+1,lId[0],lId[1],lId[2],lId[3],lId[4]);*/
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
 * u0        -> solucao                                              * 
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
              ,DOUBLE *restrict volume
              ,DOUBLE *restrict u     ,DOUBLE *restrict f
              ,INT const numel        ,short const ndf
              ,short const maxViz)
{
  INT nel;
  short carg,j;
  short col = maxViz + 1;

/*...*/  
  for(nel = 0; nel < numel;nel++){
    carg = MAT2D(nel,maxViz,faceR,col);
/*... carga constante*/
    if( carg == 1){
      for(j = 0; j< ndf;j++)
        MAT2D(nel,j,u,ndf) = MAT3D(nel,maxViz,j,faceS,col,ndf);
    }
/*...................................................................*/

/*...*/
    else if( carg == 2){
      for(j = 0; j< ndf;j++)
        MAT2D(nel,j,f,ndf) 
        = volume[nel]*MAT3D(nel,maxViz,j,faceS,col,ndf);
    }
/*...................................................................*/
  }
/*...................................................................*/
}  
/*********************************************************************/ 

/********************************************************************* 
 * UPDATECELLU : atualizacao dos valores das variaveis das celulas   *
 * com os valores das respectivas equacoes                           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * u       -> variavel nas celulas                                   * 
 * x       -> carga por elemento                                     * 
 * id      -> numera das equacoes                                    * 
 * numel   -> numero de elementos                                    * 
 * ndf     -> graus de liberdade                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * u      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 void updateCellU(DOUBLE *restrict u,DOUBLE *restrict x
                 ,INT *restrict id 
                 ,INT const numel   ,short const ndf)
{
  INT nel,lNeq;
  short jNdf;
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
 * m       -> variavel nas celulas                                   * 
 * noU     -> nao definido                                           * 
 * elU     -> valores nas celulas                                    * 
 * el      -> conectividades das celulas                             * 
 * nen     -> numero de nos por celulas                              * 
 * numel   -> numero de elementos                                    * 
 * nnode   -> numero de nos                                          * 
 * maxNo   -> numero de nos por celula maximo da malha               * 
 * ndf     -> graus de liberdade                                     * 
 * type    -> tipo de interpolacao                                   * 
 *            1 - media simples                                      * 
 *            2 - media ponderada                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * u      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 void interCellNode(Memoria *m
                   ,DOUBLE *restrict noU,DOUBLE *restrict elU
                   ,INT *restrict el                 
                   ,short *restrict nen
                   ,INT const numel      ,INT const nnode
                   ,short const maxNo    ,short const ndf
                   ,short const type)

{
  int *md=NULL;
  short j,k;
  INT nel,no;
  
  switch(type){
/*... media simple*/
    case 1:
/*...*/
      HccaAlloc(int,m,md,nnode,"md",false);
      zero(md,nnode,"int");
      zero(noU,ndf*nnode,DOUBLEC);
/*...................................................................*/

/*...*/
      for(nel = 0; nel < numel; nel++)
        for(j = 0; j   < nen[nel];j++){
          no = MAT2D(nel,j,el,maxNo) - 1;
          if( no > -1){
            for(k = 0; k   < ndf;k++)
              MAT2D(no,k,noU,ndf) += MAT2D(nel,k,elU,ndf);
            md[no]++;
          }
        }
/*...................................................................*/

/*...*/
      for(no = 0; no < nnode; no++)
        for(k = 0; k < ndf; k++)
          MAT2D(no,k,noU,ndf) = MAT2D(no,k,noU,ndf)/md[no];
/*...................................................................*/
          
/*...*/
      HccaDealloc(m,md,"md",false);
/*...................................................................*/
    break;
/*...................................................................*/

/*... media simple*/
    default:
      ERRO_OP(__FILE__,__func__,type);
    break;
  }
/*...................................................................*/
}
/*********************************************************************/
