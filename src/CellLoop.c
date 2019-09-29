#include<CellLoop.h>
/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 19/07/2018                                   * 
 *-------------------------------------------------------------------* 
 * PGEOMFORM: calculo da propriedades geometricas das celulas        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x       -> cordenadas dos pontos                                  * 
 * el      -> conetividade dos celulas                               * 
 * nelcon  -> vizinhos dos elementos                                 * 
 * nFace   -> numero de faces por celulas                            * 
 * geomType-> tipo geometrico das celulas                            * 
 * nen     -> numero de nos por celulas                              *
 * cellFace-> faces que componhe a celula                            *
 * geom    -> geometria por celula                                   *
 * face    -> geometria por face                                     *
 * iCel    -> interface de comunicao dos elmentos
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
 *            faces                                                  *
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
 * gDcca   -> menor distancia do centroide a faces desta celula      * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void pGeomForm(DOUBLE *RESTRICT x       ,INT    *RESTRICT el
              ,INT    *RESTRICT nelcon  ,short  *RESTRICT nFace
              ,short  *RESTRICT geomType,short *RESTRICT nen
              ,INT *RESTRICT cellFace
              ,Geom *RESTRICT geom      ,Face *RESTRICT face
              ,Interface *iCel
              ,short maxNo              ,short maxViz
              ,short ndm                ,INT numel)
{
  short i, j, k, aux1, aux2;
  INT nel,no,vizNel,idFace,idCell;
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

    aux1 = nFace[nel];

/*... zerando vetores*/
    for(i=0;i<MAX_NUM_PONT;i++){
      lx[i] = 0.0;
    }
    
    for(i=0;i<MAX_SN;i++){
      isnod[i] = 0;
    }

    for(i=0;i<MAX_NUM_FACE+1;i++){
      lnFace[i]    = 0;
      lnEn[i]      = 0;
      lGeomType[i] = 0;     
    }

/*... loop na celula central*/    
    lnFace[aux1]    = nFace[nel];
    lnEn[aux1]      = nen[nel];
    lGeomType[aux1] = geomType[nel];
    for(j=0;j<nen[nel];j++){
/*...*/
      no = MAT2D(nel,j,el,maxNo) - 1;
      for(k=0;k<ndm;k++){
        MAT3D(aux1,j,k,lx,maxNo,ndm) = MAT2D(no,k,x,ndm);
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
    else if(ty == TETRCELL || ty == HEXACELL || ty == PIRACELL){
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
    geom->volume[nel] = lVolume;
    for (i = 0; i<ndm; i++)  
      MAT2D(nel, i, geom->cc, ndm) = MAT2D(aux1, i, lCc, ndm);
/*...................................................................*/

/*... atribuindo valores por face*/
    for(i=0;i<nFace[nel];i++)
    {
      idFace = MAT2D(nel,i , cellFace, maxViz) - 1;
      idCell = MAT2D(idFace, 0, face->owner, 2) - 1;
      if (idCell == nel)
      {
        face->mksi[idFace]   = lmKsi[i];
        face->area[idFace]   = lfArea[i];
        face->mvSkew[idFace] = lmvSkew[i];
        for (j = 0; j < ndm; j++)
        {
          MAT2D(idFace,j,face->ksi,ndm)=MAT2D(i,j,lKsi,ndm);
          MAT2D(idFace,j,face->eta,ndm)=MAT2D(i,j,lEta,ndm);
          MAT2D(idFace,j,face->normal,ndm)=MAT2D(i,j,lNormal,ndm);
          MAT2D(idFace,j,face->xm,ndm)=MAT2D(i,j,lXm,ndm);
          MAT2D(idFace,j,face->vSkew,ndm)=MAT2D(i,j,lvSkew,ndm);
        }
      }  
    }  
/*...................................................................*/

/*... atribuido valores aos celula global*/
    for (i = 0; i<nFace[nel]; i++)
    {
      MAT2D(nel, i, geom->dcca, maxViz) = lDcca[i];
      for (j = 0; j<ndm; j++)
        MAT3D(nel, i, j, geom->xmcc, maxViz, ndm)
          = MAT2D(i, j, lXmcc, ndm);
    }
/*...................................................................*/

  }
/*...................................................................*/

/*...*/
  if(mpiVar.nPrcs > 1)
    comunicateCel(iCel,geom->cc,ndm,1);
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 26/09/2019                                   *
 *-------------------------------------------------------------------*
 * SYSTFOMDIF : calculo do sistema de equacoes para problemas        * 
 * difusao (Ax=b)                                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads   -> definicoes de cargas                                   * 
 * diff    -> tecnica da discretizacao do termo difusivo             *
 * dModel   -> configuracoes do modelo difusivo                      *
 * el      -> conetividade dos celulas                               * 
 * nelcon  -> vizinhos dos elementos                                 * 
 * nen     -> numero de nos por celulas                              * 
 * nFace   -> numero de faces por celulas                            * 
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * gVolume -> volumes das celulas                                    *
 * gDcca   -> menor distancia do centroide a faces desta celula      *
 * gXmCc   -> vetores que unem o centroide aos pontos medios das     *
 * gModKsi -> modulo do vetor ksi                                    *
 * fModKsi -> o mudolo do vetores que unem centroide da celula       * 
 *            central aos vizinhos destas                            *
 * fKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        * 
 * fEta    -> vetores paralelos as faces das celulas                 *
 * fArea   -> modulo do vetor eta                                    *
 * fNormal -> vetores normais as faces das celulas                   *
 * fXm     -> pontos medios das faces das celulas                    *
 * fModvSkew -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face                      * 
 * fvSkew  -> vetor entre o ponto medio a intersecao que une os      *
 *            centrois compartilhado nessa face                      *
 * calType -> tipo de calculo das celulas                            * 
 * geomType-> tipo geometrico das celulas                            *  
 * mat     -> material por celula                                    * 
 * density -> massa especifica com variacao temporal                 *
 * cDiffD  -> coeficiente de difusao com variacao temporal           *
 * ia      -> ponteiro para as linhas da matriz esparsa              * 
 * ja      -> ponteiro para as colunas da matriz esparsa             * 
 * a       -> matriz de coeficientes esparsa                         * 
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- fora da diagonal principal        )       *
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 * b       -> vetor de forcas                                        * 
 * id      -> numera das equacoes                                    * 
 * faceR   -> restricoes por elemento                                * 
 * faceLd1 -> carga por elemento                                     * 
 * u0      -> solucao conhecida                                      * 
 * gradU0  -> gradiente da solucao conhecido                         * 
 * rCell   -> nao definido                                           * 
 * ddt     -> discretizacao temporal                                 *
 * nEq     -> numero de equacoes                                     *
 * neqNov  -> numero de equacoes nao sobrepostas                     *
 * nAd     -> numero de termos nao nulos                             *
 * nAdR    -> numero de termos nao nulos na parte retangular         *
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
 * rCell     -> residuo por celula                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void systFormDif(Loads *loads             ,Diffusion *diff
               ,DiffModel *dModel
               ,INT    *RESTRICT el       ,INT    *RESTRICT nelcon 
               ,short  *RESTRICT nen      ,short  *RESTRICT nFace
               ,INT *RESTRICT cellFace    ,INT *RESTRICT fOwner
               ,DOUBLE *RESTRICT gVolume  ,DOUBLE *RESTRICT gDcca
               ,DOUBLE *RESTRICT gXmCc    ,DOUBLE *RESTRICT fModksi 
               ,DOUBLE *RESTRICT fKsi     ,DOUBLE *RESTRICT fEta
               ,DOUBLE *RESTRICT fArea    ,DOUBLE *RESTRICT fNormal 
               ,DOUBLE *RESTRICT fXm      ,DOUBLE *RESTRICT fModvSkew
               ,DOUBLE *RESTRICT fvSkew   ,short  *RESTRICT geomType 
               ,short  *RESTRICT calType  ,short  *RESTRICT mat     
               ,DOUBLE *RESTRICT density  ,DOUBLE *RESTRICT cDiffD
               ,INT    *RESTRICT ia       ,INT    *RESTRICT ja   
               ,DOUBLE *RESTRICT a        ,DOUBLE *RESTRICT ad
               ,DOUBLE *RESTRICT b        ,INT    *RESTRICT id
               ,short  *RESTRICT faceR    ,short  *RESTRICT faceLd1        
               ,DOUBLE *RESTRICT u0       ,DOUBLE *RESTRICT gradU0 
               ,DOUBLE *RESTRICT rCell    ,Temporal *ddt 
               ,INT nEq                   ,INT nEqNov     
               ,INT nAd                   ,INT nAdR                     
               ,short maxNo               ,short maxViz
               ,short ndm                 ,INT numel
               ,short ndf                 ,short storage
               ,bool forces               ,bool matrix 
               ,bool calRcell             ,bool unsym) 
{
  short i,j;
  INT nel, vizNel;
/*... variavel local */
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM],lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM],lVolume[MAX_NUM_FACE+1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM],lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE],lvSkew[MAX_NUM_FACE*MAX_NDM];
  short  lGeomType[MAX_NUM_FACE+1];
  short  lib;
  short  lFaceR[MAX_NUM_FACE+1];
  short  lFaceL[MAX_NUM_FACE+1];
  DOUBLE lDensity,lCoefDiffD[MAX_NUM_FACE+1];
  DOUBLE lA[(MAX_NUM_FACE+1)*MAX_NDF],lB[MAX_NDF];
  DOUBLE lu0[(MAX_NUM_FACE+1)*MAX_NDF];
  DOUBLE lGradU0[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lRcell[MAX_NDF];
  INT    lId[(MAX_NUM_FACE+1)*MAX_NDF],lViz[MAX_NUM_FACE];
  INT idFace, cellOwner,ch;
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
    if(MAT2D(nel,aux1,faceR ,aux2) != PCCELL){

/*... loop na celula central*/    
      lMat            = mat[nel]-1;
      lib             = calType[lMat];
      lVolume[aux1]   = gVolume[nel]; 
      lGeomType[aux1] = geomType[nel];
      lFaceR[aux1]    = MAT2D(nel,aux1,faceR   ,aux2);
      lFaceL[aux1]    = MAT2D(nel,aux1,faceLd1 ,aux2);
      lDensity        = MAT2D(nel,TIME_N,density ,DENSITY_LEVEL);
      lCoefDiffD[aux1]= cDiffD[nel];
      
      for(j=0;j<ndf;j++){
        MAT2D(aux1,j,lu0   ,ndf) = MAT2D(nel,j,u0   ,ndf);
        MAT2D(aux1,j,lId   ,ndf) = MAT2D(nel,j,id   ,ndf) - 1;
      }   
     
      for(j=0;j<ndm;j++)
        MAT2D(aux1,j,lGradU0,ndm) = MAT2D(nel,j,gradU0,ndm);      

      for(i=0;i<aux1;i++)
      {
        lDcca[i]   = MAT2D(nel,i,gDcca,maxViz);
        aux2       = (maxViz+1);
        lFaceR[i]  = MAT2D(nel,i,faceR   ,aux2);
        lFaceL[i]  = MAT2D(nel,i,faceLd1 ,aux2);
        for(j=0;j<ndm;j++)
        {
          MAT2D(i,j,lXmcc,ndm)=MAT3D(nel,i,j,gXmCc,maxViz,ndm);
        }
      }

/*... propriedades por face*/  
      for(i=0;i<aux1;i++) 
      {
        idFace     = MAT2D(nel, i, cellFace, maxViz) - 1;
        cellOwner  = MAT2D(idFace, 0, fOwner, 2) - 1;
        ch         = OWNER(cellOwner, nel);
        lmKsi[i]   = fModksi[idFace];
        lfArea[i]  = fArea[idFace];
        lmvSkew[i] = fModvSkew[idFace];
        for (j = 0; j<ndm; j++) 
        {
          MAT2D(i,j,lKsi,ndm)   = ch*MAT2D(idFace,j,fKsi,ndm);
          MAT2D(i,j,lEta,ndm)   = ch*MAT2D(idFace,j,fEta,ndm);
          MAT2D(i,j,lNormal,ndm)= ch*MAT2D(idFace,j,fNormal,ndm);
          MAT2D(i,j,lXm,ndm)    = MAT2D(idFace,j,fXm,ndm);
          MAT2D(i,j,lvSkew,ndm) = MAT2D(idFace,j,fvSkew,ndm);
        }
      }

/*... loop na celulas vizinhas*/    
      for(i=0;i<aux1;i++)
      {
        vizNel        = MAT2D(nel,i,nelcon,maxViz) - 1;
        lViz[i] = vizNel;
        for(j=0;j<ndf;j++)
          MAT2D(i,j,lId ,ndf) = - 1;
/*....*/
        if( vizNel != -2)
        {
          lVolume[i]    = gVolume[vizNel]; 
          lGeomType[i]  = geomType[vizNel];
          lCoefDiffD[i] = cDiffD[vizNel];
          lMat = mat[vizNel]-1;
          for(j=0;j<ndf;j++){
            MAT2D(i,j,lu0 ,ndf)   = MAT2D(vizNel,j,u0   ,ndf);
            MAT2D(i,j,lId ,ndf)   = MAT2D(vizNel,j,id   ,ndf) - 1;
          }
          for(j=0;j<ndm;j++)
            MAT2D(i,j,lGradU0,ndm)   = MAT2D(vizNel,j,gradU0,ndm);
          }
        }  
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
      cellLibDif(loads     ,diff
                ,dModel    ,lGeomType 
                ,lViz      ,lId           
                ,lKsi      ,lmKsi
                ,lEta      ,lfArea 
                ,lNormal   ,lVolume
                ,lXm       ,lXmcc
                ,lDcca     ,&lDensity
                ,lCoefDiffD
                ,lvSkew    ,lmvSkew
                ,lA        ,lB
                ,lRcell    ,ddt
                ,lFaceR    ,lFaceL            
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
            ,a           ,ad         
            ,b  
            ,lId 
            ,lA          ,lB
            ,nEq         ,nEqNov 
            ,nAd         ,nAdR     
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
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 26/06/2018                                   *
 *-------------------------------------------------------------------* 
 * SYSTFOMTRANS : calculo do sistema de equacoes para problemas      * 
 * transporte (Ax=b)                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads   -> definicoes de cargas                                   * 
 * diff    -> tecnica da discretizacao do termo difusivo             *
 * tModel  -> configuracoes do modelo de transporte                  *
 * el      -> conetividade dos celulas                               * 
 * nelcon  -> vizinhos dos elementos                                 * 
 * nen     -> numero de nos por celulas                              * 
 * nFace   -> numero de faces por celulas                            * 
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * gVolume -> volumes das celulas                                    *
 * gDcca   -> menor distancia do centroide a faces desta celula      *
 * gXmCc   -> vetores que unem o centroide aos pontos medios das     *
 * gCc     -> centroide do elemento                                  * 
 * fModKsi -> o mudolo do vetores que unem centroide da celula       * 
 *            central aos vizinhos destas                            *
 * fKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        * 
 * fEta    -> vetores paralelos as faces das celulas                 *
 * fArea   -> modulo do vetor eta                                    *
 * fNormal -> vetores normais as faces das celulas                   *
 * fXm     -> pontos medios das faces das celulas                    *
 * fModvSkew -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face                      * 
 * fvSkew  -> vetor entre o ponto medio a intersecao que une os      *
 *            centrois compartilhado nessa face                      *
 * calType -> tipo de calculo das celulas                            * 
 * geomType-> tipo geometrico das celulas                            * 
 * prop    -> propriedades dos material                              * 
 * mat     -> material por celula                                    * 
 * density -> massa especifica com variacao temporal                 *
 * cDiff   -> coeficiente de difusao com variacao temporal           *
 * ia      -> ponteiro para as linhas da matriz esparsa              * 
 * ja      -> ponteiro para as colunas da matriz esparsa             * 
 * a       -> matriz de coeficientes esparsa                         * 
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- fora da diagonal principal        )       *
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 * b       -> vetor de forcas                                        * 
 * id      -> numera das equacoes                                    * 
 * faceR   -> restricoes por elemento                                * 
 * faceL   -> carga por elemento                                     * 
 * u0      -> solucao conhecida                                      * 
 * gradU0  -> gradiente da solucao conhecido                         * 
 * vel     -> compo de velocidade conhecido                          * 
 * rCell   -> nao definido                                           * 
 * ddt     -> discretizacao temporal                                 *
 * nEq     -> numero de equacoes                                     *
 * neqNov  -> numero de equacoes nao sobrepostas                     *
 * nAd     -> numero de termos nao nulos                             *
 * nAdR    -> numero de termos nao nulos na parte retangular         *
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
 * rCell     -> residuo por celula                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void systFormTrans(Loads *loads          
               , Advection *advT           , Diffusion *diffT
               , TransModel *tModel          
               , INT    *RESTRICT el       , INT    *RESTRICT nelcon 
               , short  *RESTRICT nen      , short  *RESTRICT nFace
               , INT *RESTRICT cellFace    , INT *RESTRICT fOwner
               , DOUBLE *RESTRICT gVolume  , DOUBLE *RESTRICT gDcca
               , DOUBLE *RESTRICT gXmCc    , DOUBLE *RESTRICT gCc
               , DOUBLE *RESTRICT fModKsi  , DOUBLE *RESTRICT fKsi
               , DOUBLE *RESTRICT fEta     , DOUBLE *RESTRICT fArea
               , DOUBLE *RESTRICT fNormal  , DOUBLE *RESTRICT fXm
               , DOUBLE *RESTRICT fModvSkew, DOUBLE *RESTRICT fvSkew
               , short  *RESTRICT geomType , DOUBLE *RESTRICT prop
               , short  *RESTRICT calType  , short  *RESTRICT mat
               , DOUBLE *RESTRICT density  , DOUBLE *RESTRICT cDiff
               , INT    *RESTRICT ia       , INT    *RESTRICT ja
               , DOUBLE *RESTRICT a        , DOUBLE *RESTRICT ad
               , DOUBLE *RESTRICT b        , INT    *RESTRICT id
               , short  *RESTRICT faceR    , short  *RESTRICT faceL
               , DOUBLE *RESTRICT u0       , DOUBLE *RESTRICT gradU0 
               , DOUBLE *RESTRICT vel                              
               , DOUBLE *RESTRICT rCell    , Temporal *ddt 
               , INT nEq                   , INT nEqNov
               , INT nAd                   , INT nAdR                     
               , short maxNo               , short maxViz
               , short ndm                 , INT numel
               , short ndf                 , short storage
               , bool forces               , bool matrix 
               , bool calRcell             , bool unsym) 
{
  short i,j;
  short nThreads = ompVar.nThreadsCell;
  INT nel,vizNel;
/*... variavel local */
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM],lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM],lVolume[MAX_NUM_FACE+1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM],lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE],lvSkew[MAX_NUM_FACE*MAX_NDM];
  short  lGeomType[MAX_NUM_FACE+1];
  short  lib;
  short  lFaceR[MAX_NUM_FACE+1];
  short  lFaceL[MAX_NUM_FACE+1];
  DOUBLE lDensity[(MAX_NUM_FACE+1)],lCoefDiffD[MAX_NUM_FACE + 1];;
  DOUBLE lA[(MAX_NUM_FACE+1)*MAX_NDF],lB[MAX_NDF];
  DOUBLE lProp[(MAX_NUM_FACE+1)*MAXPROP];
  DOUBLE lu0[(MAX_NUM_FACE+1)*MAX_NDF];
  DOUBLE lGradU0[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lVel[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lCc[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lRcell[MAX_NDF];
  INT    lId[(MAX_NUM_FACE+1)*MAX_NDF],lViz[MAX_NUM_FACE];
  INT idFace, cellOwner, ch;
  short  aux1,aux2,lMat;

/*...*/
  if (ompVar.fCell) {
/*... loop nas celulas*/
      aux2 = maxViz + 1;
#pragma omp parallel  for default(none) num_threads(nThreads)\
   private(nel,i,j,aux1,lId,lMat,lib,lVolume,lGeomType,lu0\
          ,lFaceR,lA,lB,lFaceL,lProp,lGradU0\
          ,lVel,lCc,lmKsi,lfArea,lDcca,lmvSkew,lKsi,lEta\
          ,lNormal,lXm,lXmcc,lvSkew,vizNel,lViz,lRcell,lCoefDiffD,lDensity\
          ,idFace,cellOwner,ch)\
   shared(aux2,ndm,ndf,numel,maxViz,calRcell,rCell,nFace,mat\
         ,calType,gVolume, geomType,faceR,u0,gCc,loads\
         ,faceL,density,cDiff,prop,gradU0,vel,fModKsi\
         ,fArea,gDcca,fModvSkew,fKsi,fEta,fNormal,fXm,gXmCc,fvSkew\
         ,nelcon,id,advT,diffT\
         ,ddt,nen,ia,ja,a,ad,b,nEq,nEqNov,nAd\
         ,nAdR,storage,forces,matrix,unsym,tModel\
         ,fOwner,cellFace) 
    for (nel = 0; nel<numel; nel++) 
    {
/*...*/
      if (calRcell)
        for (j = 0; j<ndf; j++)
          MAT2D(nel, j, rCell, ndf) = 0.e0;;
/*...*/
      aux1 = nFace[nel];
/*... elementos com equacoes*/
      if (MAT2D(nel, aux1, faceR, aux2) != PCCELL)
      {

/*... zerando vetores*/
        for (j = 0; j<(MAX_NUM_FACE + 1)*MAX_NDF; j++)
        {
          lId[j] = -1;
          lu0[j] = 0.e0;
        }
/*... loop na celula central*/
        lMat             = mat[nel] - 1;
        lib              = calType[lMat];
        lVolume[aux1]    = gVolume[nel];
        lGeomType[aux1]  = geomType[nel];
        lFaceR[aux1]     = MAT2D(nel, aux1, faceR, aux2);
        lFaceL[aux1]     = MAT2D(nel, aux1, faceL, aux2);
        lDensity[aux1]   = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
        lCoefDiffD[aux1] = cDiff[nel];

        for (j = 0; j<ndf; j++)
        {
          MAT2D(aux1, j, lu0, ndf) = MAT2D(nel, j, u0, ndf);
          MAT2D(aux1, j, lId, ndf) = MAT2D(nel, j, id, ndf) - 1;
        }

        for (j = 0; j<MAXPROP; j++)
          MAT2D(aux1, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);

        for (j = 0; j<ndm; j++) 
        {
          MAT2D(aux1, j, lGradU0, ndm) = MAT2D(nel, j, gradU0, ndm);
          MAT2D(aux1, j, lVel, ndm) = MAT2D(nel, j, vel, ndm);
          MAT2D(aux1, j, lCc, ndm) = MAT2D(nel, j, gCc, ndm);
        }

        for (i = 0; i<aux1; i++)
        {
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lFaceR[i] = MAT2D(nel, i, faceR, aux2);
          lFaceL[i] = MAT2D(nel, i, faceL, aux2);
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);
          }
        }
/*...................................................................*/

/*... propriedades por face*/
        for (i = 0; i<aux1; i++)
        {
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
        }
/*...................................................................*/

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++)
        {
          vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2) 
          {
            lVolume[i]    = gVolume[vizNel];
            lGeomType[i]  = geomType[vizNel];
            lCoefDiffD[i] = cDiff[vizNel];
            lDensity[i]   = MAT2D(vizNel, TIME_N, density, DENSITY_LEVEL);
            lMat          = mat[vizNel] - 1;
            for (j = 0; j<ndf; j++)
            {
              MAT2D(i, j, lu0, ndf) = MAT2D(vizNel, j, u0, ndf);
              MAT2D(i, j, lId, ndf) = MAT2D(vizNel, j, id, ndf) - 1;
            }
            for (j = 0; j<ndm; j++)
            {
              MAT2D(i, j, lGradU0, ndm) = MAT2D(vizNel, j, gradU0, ndm);
              MAT2D(i, j, lVel, ndm) = MAT2D(vizNel, j, vel, ndm);
              MAT2D(i, j, lCc, ndm) = MAT2D(vizNel, j, gCc, ndm);
            }
            for (j = 0; j<DIFPROP; j++)
              MAT2D(i, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);
          }
        }
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibTrans(loads
                   , advT     , diffT
                   , tModel
                   , lGeomType, lProp
                   , lViz     , lId
                   , lKsi     , lmKsi
                   , lEta     , lfArea
                   , lNormal  , lVolume
                   , lXm      , lXmcc
                   , lDcca
                   , lDensity , lCoefDiffD
                   , lvSkew   , lmvSkew
                   , lA       , lB
                   , lRcell   , ddt
                   , lFaceR   , lFaceL
                   , lu0      , lGradU0
                   , lVel     , lCc
                   , nen[nel] , nFace[nel]
                   , ndm      , lib
                   , nel);
/*...................................................................*/

/*... residuo da celula*/
        if (calRcell)
          for (j = 0; j<ndf; j++)
            MAT2D(nel,j,rCell,ndf) = lRcell[j];
/*...................................................................*/

/*...*/
        assbly(ia        ,ja
              ,a         ,ad
              ,b         ,lId
              ,lA        ,lB
              ,nEq       ,nEqNov
              ,nAd       ,nAdR
              ,nFace[nel],ndf
              ,storage   ,forces
              ,matrix    ,unsym);
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  else{
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
      if(MAT2D(nel,aux1,faceR ,aux2) != PCCELL){

/*... zerando vetores*/
        for(j=0;j<(MAX_NUM_FACE+1)*MAX_NDF;j++){
          lId[j] = -1;
          lu0[j] = 0.e0;    
        }
/*... loop na celula central*/    
        lMat             = mat[nel]-1;
        lib              = calType[lMat];
        lVolume[aux1]    = gVolume[nel]; 
        lGeomType[aux1]  = geomType[nel];
        lFaceR[aux1]     = MAT2D(nel,aux1,faceR   ,aux2);
        lFaceL[aux1]     = MAT2D(nel,aux1,faceL ,aux2);
        lDensity[aux1]   = MAT2D(nel,TIME_N,density ,DENSITY_LEVEL);
        lCoefDiffD[aux1] = cDiff[nel];
      
        for(j=0;j<ndf;j++)
        { 
          MAT2D(aux1,j,lu0   ,ndf) = MAT2D(nel,j,u0   ,ndf);
          MAT2D(aux1,j,lId   ,ndf) = MAT2D(nel,j,id   ,ndf) - 1;
        }
      
        for(j=0;j<MAXPROP;j++)
          MAT2D(aux1,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
      
        for(j=0;j<ndm;j++)
        {
          MAT2D(aux1,j,lGradU0,ndm) = MAT2D(nel,j,gradU0,ndm);
          MAT2D(aux1,j,lVel   ,ndm) = MAT2D(nel,j,vel   ,ndm);
          MAT2D(aux1,j,lCc    ,ndm) = MAT2D(nel,j,gCc   ,ndm);
        }

        for(i=0;i<aux1;i++)
        {
          lDcca[i]   = MAT2D(nel,i,gDcca   ,maxViz);
          lFaceR[i]  = MAT2D(nel,i,faceR,aux2);
          lFaceL[i]  = MAT2D(nel,i,faceL,aux2);
          for(j=0;j<ndm;j++)
          {
            MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,gXmCc  ,maxViz,ndm);
          }
        }

/*... propriedades por face*/
        for (i = 0; i<aux1; i++)
        {
          idFace     = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner  = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch         = OWNER(cellOwner, nel);
          lmKsi[i]   = fModKsi[idFace];
          lfArea[i]  = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lKsi, ndm)    = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm)    = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm)     = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm)  = MAT2D(idFace, j, fvSkew, ndm);
          }
        }
 
/*... loop na celulas vizinhas*/    
        for(i=0;i<aux1;i++){
          vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
          lViz[i] = vizNel;
          if( vizNel != -2) {
            lVolume[i]   = gVolume[vizNel]; 
            lGeomType[i] = geomType[vizNel];
            lCoefDiffD[i] = cDiff[vizNel];
            lDensity[i]  = MAT2D(vizNel, TIME_N,density ,DENSITY_LEVEL);
            lMat         = mat[vizNel]-1;
            for(j=0;j<ndf;j++){
              MAT2D(i,j,lu0 ,ndf)   = MAT2D(vizNel,j,u0   ,ndf);
              MAT2D(i,j,lId ,ndf)   = MAT2D(vizNel,j,id   ,ndf) - 1;
            }
            for(j=0;j<ndm;j++){
              MAT2D(i,j,lGradU0,ndm) = MAT2D(vizNel,j,gradU0,ndm);
              MAT2D(i,j,lVel   ,ndm) = MAT2D(vizNel,j,vel   ,ndm);
              MAT2D(i,j,lCc    ,ndm) = MAT2D(vizNel,j,gCc    ,ndm);
            }
            for(j=0;j<DIFPROP;j++)
              MAT2D(i,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
          }
        }  
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibTrans(loads     
                    ,advT      ,diffT
                    ,tModel    
                    ,lGeomType ,lProp 
                    ,lViz      ,lId           
                    ,lKsi      ,lmKsi
                    ,lEta      ,lfArea 
                    ,lNormal   ,lVolume
                    ,lXm       ,lXmcc
                    ,lDcca     
                    ,lDensity  ,lCoefDiffD
                    ,lvSkew    ,lmvSkew
                    ,lA        ,lB
                    ,lRcell    ,ddt
                    ,lFaceR    ,lFaceL            
                    ,lu0       ,lGradU0     
                    ,lVel      ,lCc
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
              ,a           ,ad              
              ,b  
              ,lId 
              ,lA          ,lB
              ,nEq         ,nEqNov
              ,nAd         ,nAdR     
              ,nFace[nel]  ,ndf 
              ,storage     ,forces
              ,matrix      ,unsym); 
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 23/08/2019                                   * 
 *-------------------------------------------------------------------* 
 * SYSTFOMSIMPLEVELLM: calculo do sistema de equacoes para problemas * 
 * de escomaneto de fluidos ( Vel )                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loadsVel  -> definicoes de cargas de velocidades                  * 
 * loadsPres -> definicoes de cargas de pressao                      * 
 * advVel    -> tecnica da discretizacao do termo advecao            *
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * ModelMomentum -> termos/modelos da equacao de momento linear      *
 * iCel    -> interface de elementos                                 *
 * typeSimple-> tipo do metodo simple                                *
 * el        -> conetividade dos celulas                             * 
 * nelcon    -> vizinhos dos elementos                               * 
 * nen       -> numero de nos por celulas                            * 
 * nFace     -> numero de faces por celulas                          * 
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * gVolume -> volumes das celulas                                    *
 * gDcca   -> menor distancia do centroide a faces desta celula      *
 * gXmCc   -> vetores que unem o centroide aos pontos medios das     *
 * gCc     -> centroide do elemento                                  * 
 * fModKsi -> o mudolo do vetores que unem centroide da celula       * 
 *            central aos vizinhos destas                            *
 * fKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        * 
 * fEta    -> vetores paralelos as faces das celulas                 *
 * fArea   -> modulo do vetor eta                                    *
 * fNormal -> vetores normais as faces das celulas                   *
 * fXm     -> pontos medios das faces das celulas                    *
 * fModvSkew -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face                      * 
 * fvSkew  -> vetor entre o ponto medio a intersecao que une os      *
 *            centrois compartilhado nessa face                      *
 * calType -> tipo de calculo das celulas                            * 
 * geomType-> tipo geometrico das celulas                            *  
 * mat     -> material por celula                                    * 
 * ia        -> ponteiro para as linhas da matriz esparsa            * 
 * ja        -> ponteiro para as colunas da matriz esparsa           * 
 * a         -> matriz de coeficientes esparsa                       * 
 *              ( CSR - nao utiliza                            )     *
 *              ( CSRD/CSRC- fora da diagonal principal        )     *
 * ad        -> matrix de coeficientes esparsa                       *
 *              ( CSR - matriz completa                        )     *
 *              ( CSRD/CSRC- diagonal principal                )     *
 * b         -> vetor de forcas                                      * 
 * id        -> numera das equacoes                                  * 
 * faceVelR  -> restricoes por elemento de velocidades               * 
 * faceVelL  -> carga por elemento de velocidades                    * 
 * facePresR -> restricoes por elemento de pressao                   * 
 * facePresL -> carga por elemento de pressao                        * 
 * pres      -> campo de pressao conhecido                           * 
 * gradVel   -> gradiente da solucao conhecido                       * 
 * dField    -> matriz D do metodo simple                            * 
 * underU    -> fator underrelaxtion sinple                          * 
 * vel       -> campo de velocidade conhecido                        * 
 * rCell     -> nao definido                                         *
 * stressR   -> tensor residual                                      *
 * density   -> massa especifica com variacao temporal               *  
 * dViscosity-> viscosidade dinamica com variacao temporal           *
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  *  
 * ddt       -> discretizacao temporal                               *
 * nEq       -> numero de equacoes                                   *
 * neqNov    -> numero de equacoes nao sobrepostas                   *
 * nAd       -> numero de termos nao nulos                           *
 * nAdR      -> numero de termos nao nulos na parte retangular       *
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * numel     -> numero de toral de celulas                           * 
 * ndf       -> graus de liberdade                                   * 
 * storage   -> tecnica de armazenamento da matriz esparsa           *
 * ntn       -> numero de termo do tensor de tensao                  * 
 * forces    -> mantagem no vetor de forcas                          * 
 * matrix    -> mantagem da matriz de coeficientes                   * 
 * calRcell  -> calculo do residuo de celula                         * 
 * unsym     -> matiz nao simetrica                                  * 
 * sPressure -> reconstrucao de segunda ordem para pressoes nas      *
 *              faces                                                *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ad,a      -> coeficiente da linha i     (matriz = true)           *
 * b         -> vetor de forca da linha i  (forces = true)           *
 * rCell     -> residuo por celula                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 * b     = | bx1 bx2 ... bxn by1 by2 ... byn bz1 bz2 ... bzn |       * 
 * rCell = | rx1 rx2 ... rxn ry1 ry2 ... ryn rz1 rz2 ... rzn |       *
 *                                                                   *
 * e necassirio comunicar o arranjo dField                           *   
 *********************************************************************/
void systFormSimpleVelLm(Loads *loadsVel   , Loads *loadsPres    
    , Advection *advVel                    , Diffusion *diffVel
    , Turbulence *tModel                   , MomentumModel *ModelMomentum
    , Interface *iCel                      , short typeSimple     
    , INT    *RESTRICT el                  , INT    *RESTRICT nelcon 
    , short  *RESTRICT nen                 , short  *RESTRICT nFace
    , INT *RESTRICT cellFace               , INT *RESTRICT fOwner
    , DOUBLE *RESTRICT gVolume             , DOUBLE *RESTRICT gDcca
    , DOUBLE *RESTRICT gXmCc               , DOUBLE *RESTRICT gCc
    , DOUBLE *RESTRICT fModKsi             , DOUBLE *RESTRICT fKsi
    , DOUBLE *RESTRICT fEta                , DOUBLE *RESTRICT fArea
    , DOUBLE *RESTRICT fNormal             , DOUBLE *RESTRICT fXm
    , DOUBLE *RESTRICT fModvSkew           , DOUBLE *RESTRICT fvSkew
    , short  *RESTRICT geomType            
    , short  *RESTRICT calType             , short  *RESTRICT mat
    , INT    *RESTRICT ia                  , INT    *RESTRICT ja
    , DOUBLE *RESTRICT a                   , DOUBLE *RESTRICT ad
    , DOUBLE *RESTRICT b                   , INT    *RESTRICT id
    , short  *RESTRICT faceVelR            , short  *RESTRICT faceVelL       
    , short  *RESTRICT facePresR           , short  *RESTRICT facePresL             
    , DOUBLE *RESTRICT pres                , DOUBLE *RESTRICT gradPres
    , DOUBLE *RESTRICT vel                 , DOUBLE *RESTRICT gradVel
    , DOUBLE *RESTRICT dField              , DOUBLE underU 
    , DOUBLE *RESTRICT rCell               , DOUBLE *RESTRICT stressR  
    , DOUBLE *RESTRICT density             , DOUBLE *RESTRICT dViscosity 
    , DOUBLE *RESTRICT eddyViscosity       , DOUBLE *RESTRICT wallPar
    , DOUBLE densityMed                    , Temporal *ddt                     
    , INT nEq                              , INT nEqNov
    , INT nAd                              , INT nAdR                 
    , short maxNo                          , short maxViz
    , short ndm                            , INT numel
    , short ndf                            , short storage
    , short ntn                            , bool forces      
    , bool matrix                          , bool calRcell
    , bool unsym                           , bool sPressure) 
{   
  bool fTurb       = tModel->fTurb
      ,fTurbStruct = tModel->fTurbStruct
      ,fWallModel  = tModel->fWall;
  short i,j,k;
  short nThreads = ompVar.nThreadsCell;
  INT nel,vizNel;
  
/*... variavel local */
  short  aux1,aux2,lMat,lib;
  short  lGeomType[MAX_NUM_FACE+1],
         lFaceVelR[MAX_NUM_FACE+1], lFaceVelL[MAX_NUM_FACE+1],
         lFacePresR[MAX_NUM_FACE+1], lFacePresL[MAX_NUM_FACE+1];
  INT    lId[(MAX_NUM_FACE+1)*MAX_NDF], lViz[MAX_NUM_FACE];
  INT    idFace, cellOwner, ch;
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM]    , lmKsi[MAX_NUM_FACE],
         lEta[MAX_NUM_FACE*MAX_NDM]    , lfArea[MAX_NUM_FACE],
         lNormal[MAX_NUM_FACE*MAX_NDM] , lVolume[MAX_NUM_FACE+1],
         lXm[MAX_NUM_FACE*MAX_NDM]     , lXmcc[MAX_NUM_FACE*MAX_NDM],
         lDcca[MAX_NUM_FACE]           , lmvSkew[MAX_NUM_FACE],
         lvSkew[MAX_NUM_FACE*MAX_NDM]  , lDensity[(MAX_NUM_FACE+1)],
         lViscosity[(MAX_NUM_FACE+1)*2], lA[MAX_NUM_FACE+ MAX_NDM],
         lB[MAX_NDF], 
         lPres[(MAX_NUM_FACE+1)]       , lCc[(MAX_NUM_FACE+1)*MAX_NDM],
         lGradVel[(MAX_NUM_FACE+1)*MAX_NDM*MAX_NDF],
         lGradPres[(MAX_NUM_FACE+1)*MAX_NDM], lRcell[MAX_NDF],
         lVel[(MAX_NUM_FACE+1)*MAX_NDM]     , lDfield[(MAX_NUM_FACE+1)*MAX_NDM],
         lStressR[(MAX_NUM_FACE+1)*6]       ,lWallPar[NWALLPAR];
  
/*...*/
  if(ompVar.fCell)
  {
/*... loop nas celulas*/
    aux2 = maxViz + 1;
  #pragma omp parallel  for default(none) num_threads(nThreads)\
     private(nel,i,j,k,aux1,lId,lPres,lMat,lib,lVolume,lGeomType\
          ,lFaceVelR,lA,lB,lDfield\
          ,lFaceVelL,lFacePresR,lFacePresL,lDensity,lGradPres\
          ,lVel,lCc,lGradVel,lmKsi,lfArea,lDcca,lmvSkew,lKsi,lEta\
          ,lNormal,lXm,lXmcc,lvSkew,vizNel,lViz,lRcell,lViscosity\
          ,lStressR,lWallPar,idFace,cellOwner,ch)\
     shared(aux2,ndm,ndf,numel,maxViz,calRcell,rCell,nFace,mat,ntn\
         ,calType,gVolume, geomType, faceVelR, faceVelL, facePresR\
         ,facePresL,density,gradPres,vel,gCc,gradVel,fModKsi\
         ,fArea,gDcca,fModvSkew,fKsi,fEta,fNormal,fXm,gXmCc,fvSkew\
         ,nelcon,id,loadsVel,loadsPres,advVel,diffVel,typeSimple\
         ,ddt,underU,sPressure,nen,ia,ja,a,ad,b,nEq,nEqNov,nAd\
         ,nAdR,storage,forces,matrix,unsym,pres,dField,dViscosity,eddyViscosity\
         ,stressR,tModel,ModelMomentum,wallPar,densityMed\
         ,fOwner,cellFace,fTurb,fTurbStruct,fWallModel)

/*... loop nas celulas*/
   for(nel=0;nel<numel;nel++)
    {
/*...*/
      if(calRcell)
      {
        rCell[nel]       = 0.e0;  
        rCell[numel+nel] = 0.e0;  
        if(ndf == 3 ) rCell[2*numel+nel] = 0.e0; 
      }  
/*...................................................................*/

/*...*/
      aux1    = nFace[nel];
/*... elementos com equacoes*/
      if(MAT2D(nel,aux1,faceVelR ,aux2) != PCCELL)
      {
 
/*... zerando vetores*/
        for(j=0;j<(MAX_NUM_FACE+1)*MAX_NDF;j++) 
          lId[j]  = -1;
       
        for(j=0;j<MAX_NUM_FACE+1;j++) 
          lPres[j] = 0.e0;    
       
/*... loop na celula central*/    
        lMat              = mat[nel]-1;
        lib               = calType[lMat];
        lVolume[aux1]     = gVolume[nel]; 
        lGeomType[aux1]   = geomType[nel];
        lPres[aux1]       = pres[nel];
        lFaceVelR[aux1]   = MAT2D(nel,aux1,faceVelR ,aux2);
        lFaceVelL[aux1]   = MAT2D(nel,aux1,faceVelL ,aux2);
        lFacePresR[aux1]  = MAT2D(nel,aux1,facePresR ,aux2);
        lFacePresL[aux1]  = MAT2D(nel,aux1,facePresL ,aux2);
/*... viscosidade dinamica e turbulentea*/
        MAT2D(aux1, 0, lViscosity, 2) = dViscosity[nel];
        if(fTurb)
          MAT2D(aux1, 1, lViscosity, 2) = eddyViscosity[nel];
/*...................................................................*/

/*...*/
        lDensity[aux1]    = MAT2D(nel,TIME_N,density ,DENSITY_LEVEL);
/*...................................................................*/

/*...*/
        lId[aux1]         = id[nel] - 1;
/*...................................................................*/

/*...*/
        for(j=0;j<ndm;j++)
        {
          MAT2D(aux1,j,lGradPres,ndm) = MAT2D(nel,j,gradPres,ndm);
          MAT2D(aux1,j,lVel     ,ndm) = MAT2D(nel,j,vel     ,ndm);
          MAT2D(aux1,j,lCc      ,ndm) = MAT2D(nel,j,gCc      ,ndm);
          MAT2D(aux1,j,lDfield  ,ndm) = MAT2D(nel,j,dField  ,ndm);
        }        
/*...................................................................*/

/*...*/
        if(fTurbStruct)
          for(j=0;j<ntn;j++)
            MAT2D(aux1,j,lStressR,ntn)  = MAT2D(nel,j,stressR,ntn);
/*...................................................................*/

/*...*/
        for(i=0;i<ndf;i++)
          for(j=0;j<ndm;j++)
            MAT3D(aux1,i,j,lGradVel,ndf,ndm) 
                               = MAT3D(nel,i,j,gradVel,ndf,ndm);
/*...................................................................*/

/*...*/
        for (i = 0; i<aux1; i++)
        {
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lFaceVelR[i] = MAT2D(nel, i, faceVelR, aux2);
          lFaceVelL[i] = MAT2D(nel, i, faceVelL, aux2);
          lFacePresR[i] = MAT2D(nel, i, facePresR, aux2);
          lFacePresL[i] = MAT2D(nel, i, facePresL, aux2);
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
        }
/*...................................................................*/

/*... loop na celulas vizinhas*/    
        for(i=0;i<aux1;i++)
        {
          vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
          lViz[i] = vizNel;
          if( vizNel != -2) 
          {
            lVolume[i]   = gVolume[vizNel]; 
            lGeomType[i] = geomType[vizNel];
            lMat         = mat[vizNel]-1;
            lPres[i]     = pres[vizNel];
/*... viscusidade dinamica e turbulentea*/
            MAT2D(i, 0, lViscosity, 2) = dViscosity[vizNel];
            if(fTurb)
             MAT2D(i, 1, lViscosity, 2) = eddyViscosity[vizNel];
/*...*/
            lDensity[i]  = MAT2D(vizNel, TIME_N,density ,DENSITY_LEVEL);
/*.....................................................................*/

/*...*/
            lId[i]       = id[vizNel] - 1;
/*.....................................................................*/

/*...*/
            for(j=0;j<ndm;j++)
            {
              MAT2D(i,j,lGradPres,ndm) = MAT2D(vizNel,j,gradPres,ndm);
              MAT2D(i,j,lVel     ,ndm) = MAT2D(vizNel,j,vel     ,ndm);
              MAT2D(i,j,lCc      ,ndm) = MAT2D(vizNel,j,gCc      ,ndm);
              MAT2D(i ,j, lDfield, ndm) = MAT2D(vizNel, j, dField, ndm);
            }
/*.....................................................................*/
 
/*...*/
            if(fTurbStruct)
              for(j=0;j<ntn;j++)
                MAT2D(i,j,lStressR,ntn)  = MAT2D(vizNel,j,stressR,ntn);
/*.....................................................................*/
 
/*...*/
            for(k=0;k<ndf;k++)
              for(j=0;j<ndm;j++)
                MAT3D(i,k,j,lGradVel,ndf,ndm) 
                               = MAT3D(vizNel,k,j,gradVel,ndf,ndm);
/*.....................................................................*/

          }
        }  
/*...................................................................*/

/*...*/ 
        if(fWallModel)   
          for(i=0;i<NWALLPAR;i++)
            lWallPar[i] = MAT2D(nel,i,wallPar,NWALLPAR);
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibSimpleVelLm(loadsVel , loadsPres    
                    , advVel        , diffVel  
                    , tModel        , ModelMomentum
                    , typeSimple    , lGeomType     
                    , lViz          , lId            
                    , lKsi          , lmKsi 
                    , lEta          , lfArea  
                    , lNormal       , lVolume 
                    , lXm           , lXmcc 
                    , lDcca         , lCc 
                    , lvSkew        , lmvSkew 
                    , lA            , lB 
                    , lRcell        , ddt 
                    , lFaceVelR     , lFaceVelL 
                    , lFacePresR    , lFacePresL             
                    , lPres         , lGradPres     
                    , lVel          , lGradVel 
                    , lDensity      , lViscosity 
                    , lDfield       , lStressR 
                    , lWallPar      , densityMed 
                    , underU        , sPressure 
                    , nen[nel]      , nFace[nel]  
                    , ndm           , lib    
                    , nel);    
/*...................................................................*/

/*... residuo da celula*/
        if(calRcell)
        {
          rCell[nel]       = lRcell[0];  
          rCell[numel+nel] = lRcell[1];  
          if(ndf == 3 ) rCell[2*numel+nel] = lRcell[2];  
        }
/*...................................................................*/

/*...*/
        for(j=0;j<ndm;j++)
   	      MAT2D(nel,j,dField,ndm) = MAT2D(aux1, j, lDfield, ndm);
/*...................................................................*/
      
/*...*/
        assblySimple(ia    ,ja
                    ,a           ,ad              
                    ,b           ,lId 
                    ,lA          ,lB
                    ,nEq         ,nEqNov
                    ,nAd         ,nAdR     
                    ,nFace[nel]  ,ndf 
                    ,storage     ,forces
                    ,matrix      ,unsym); 
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... sequencial*/
  else
  {
/*... loop nas celulas*/
    aux2    = maxViz+1;
    for(nel=0;nel<numel;nel++)
    {
/*...*/
      if(calRcell)
      {
        rCell[nel]       = 0.e0;  
        rCell[numel+nel] = 0.e0;  
        if(ndf == 3 ) rCell[2*numel+nel] = 0.e0; 
      }  
/*...................................................................*/

/*...*/
      aux1    = nFace[nel];
/*... elementos com equacoes*/
      if(MAT2D(nel,aux1,faceVelR ,aux2) != PCCELL)
      {
 
/*... zerando vetores*/
        for(j=0;j<(MAX_NUM_FACE+1)*MAX_NDF;j++) 
          lId[j]  = -1;
       
        for(j=0;j<MAX_NUM_FACE+1;j++) 
          lPres[j] = 0.e0;    
       
/*... loop na celula central*/    
        lMat              = mat[nel]-1;
        lib               = calType[lMat];
        lVolume[aux1]     = gVolume[nel]; 
        lGeomType[aux1]   = geomType[nel];
        lFaceVelR[aux1]   = MAT2D(nel,aux1,faceVelR ,aux2);
        lFaceVelL[aux1]   = MAT2D(nel,aux1,faceVelL ,aux2);
        lFacePresR[aux1]  = MAT2D(nel,aux1,facePresR ,aux2);
        lFacePresL[aux1]  = MAT2D(nel,aux1,facePresL ,aux2);
        lDensity[aux1]    = MAT2D(nel,TIME_N,density ,DENSITY_LEVEL);  
/*... viscosidade dinamica e turbulentea*/
        MAT2D(aux1, 0, lViscosity, 2) = dViscosity[nel];
        if(fTurb) 
          MAT2D(aux1, 1, lViscosity, 2) = eddyViscosity[nel];

        lPres[aux1]       = pres[nel];
        lId[aux1]         = id[nel] - 1;  

/*...*/
        for(j=0;j<ndm;j++)
        {
          MAT2D(aux1,j,lGradPres,ndm) = MAT2D(nel,j,gradPres,ndm);
          MAT2D(aux1,j,lVel     ,ndm) = MAT2D(nel,j,vel     ,ndm);
          MAT2D(aux1,j,lCc      ,ndm) = MAT2D(nel,j,gCc      ,ndm);
          MAT2D(aux1,j,lDfield  ,ndm) = MAT2D(nel,j,dField  ,ndm);
        }  
/*...................................................................*/

/*...*/
        if(fTurbStruct)
          for(j=0;j<ntn;j++) 
           MAT2D(aux1,j,lStressR,ntn)  = MAT2D(nel,j,stressR,ntn);
/*...................................................................*/

/*...*/
        for(i=0;i<ndf;i++)
          for(j=0;j<ndm;j++)
            MAT3D(aux1,i,j,lGradVel,ndf,ndm) 
                               = MAT3D(nel,i,j,gradVel,ndf,ndm);  
/*...................................................................*/

/*...*/
        for (i = 0; i<aux1; i++)
        {
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lFaceVelR[i] = MAT2D(nel, i, faceVelR, aux2);
          lFaceVelL[i] = MAT2D(nel, i, faceVelL, aux2);
          lFacePresR[i] = MAT2D(nel, i, facePresR, aux2);
          lFacePresL[i] = MAT2D(nel, i, facePresL, aux2);  
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);  
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
        }  
/*...................................................................*/

/*... loop na celulas vizinhas*/    
        for(i=0;i<aux1;i++)
        {
          vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
          lViz[i] = vizNel;
          if( vizNel != -2) 
          {
            lVolume[i]   = gVolume[vizNel]; 
            lGeomType[i] = geomType[vizNel];
            lDensity[i]  = MAT2D(vizNel, TIME_N,density ,DENSITY_LEVEL);  
/*... viscusidade dinamica e turbulentea*/
            MAT2D(i, 0, lViscosity, 2) = dViscosity[vizNel];
            if(fTurb)
              MAT2D(i, 1, lViscosity, 2) = eddyViscosity[vizNel];

            lMat         = mat[vizNel]-1;
            lPres[i]     = pres[vizNel];
            lId[i]       = id[vizNel] - 1;

            for(j=0;j<ndm;j++)
            {
              MAT2D(i,j,lGradPres,ndm) = MAT2D(vizNel,j,gradPres,ndm);
              MAT2D(i,j,lVel     ,ndm) = MAT2D(vizNel,j,vel     ,ndm);
              MAT2D(i,j,lCc      ,ndm) = MAT2D(vizNel,j,gCc      ,ndm);
              MAT2D(i ,j, lDfield, ndm) = MAT2D(vizNel, j, dField, ndm);
            }

            if(fTurbStruct)
              for(j=0;j<ntn;j++)
                MAT2D(i,j,lStressR,ntn)  = MAT2D(vizNel,j,stressR,ntn);

            for(k=0;k<ndf;k++)
              for(j=0;j<ndm;j++)
                MAT3D(i,k,j,lGradVel,ndf,ndm) 
                               = MAT3D(vizNel,k,j,gradVel,ndf,ndm);

          }
        }    
/*...................................................................*/

/*...*/  
        if(fWallModel)
          for(i=0;i<NWALLPAR;i++)
            lWallPar[i] = MAT2D(nel,i,wallPar,NWALLPAR);
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibSimpleVelLm(loadsVel , loadsPres    
                    , advVel        , diffVel  
                    , tModel        , ModelMomentum
                    , typeSimple 
                    , lGeomType     
                    , lViz          , lId            
                    , lKsi          , lmKsi 
                    , lEta          , lfArea  
                    , lNormal       , lVolume 
                    , lXm           , lXmcc 
                    , lDcca         , lCc 
                    , lvSkew        , lmvSkew 
                    , lA            , lB 
                    , lRcell        , ddt 
                    , lFaceVelR     , lFaceVelL 
                    , lFacePresR    , lFacePresL             
                    , lPres         , lGradPres     
                    , lVel          , lGradVel 
                    , lDensity      , lViscosity 
                    , lDfield       , lStressR 
                    , lWallPar      , densityMed 
                    , underU        , sPressure 
                    , nen[nel]      , nFace[nel]  
                    , ndm           , lib    
                    , nel);     
/*...................................................................*/

/*... residuo da celula*/
        if(calRcell)
        {
          rCell[nel]       = lRcell[0];  
          rCell[numel+nel] = lRcell[1];  
          if(ndf == 3 ) rCell[2*numel+nel] = lRcell[2];  
        }  
/*...................................................................*/

/*...*/
        for(j=0;j<ndm;j++)
   	      MAT2D(nel,j,dField,ndm) = MAT2D(aux1, j, lDfield, ndm);  
/*...................................................................*/
      
/*...*/
        assblySimple(ia    ,ja
                    ,a           ,ad              
                    ,b           ,lId 
                    ,lA          ,lB
                    ,nEq         ,nEqNov
                    ,nAd         ,nAdR     
                    ,nFace[nel]  ,ndf 
                    ,storage     ,forces
                    ,matrix      ,unsym);     
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if(mpiVar.nPrcs > 1 ) comunicateCel(iCel,dField,ndm,1);
/*...................................................................*/
}
/*********************************************************************/ 

/*********************************************************************
* Data de criacao    : 10/09/2016                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* VELEXP: calculo do sistema de equacoes para problemas             * 
* de escomaneto de fluidos ( Vel )                                  *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* loadsVel  -> definicoes de cargas de velocidades                  *
* loadsPres -> definicoes de cargas de pressao                      *
* advVel    -> tecnica da discretizacao do termo advecao            *
* diffVel   -> tecnica da discretizacao do termo difusivo           *
* el        -> conetividade dos celulas                             *
* nelcon    -> vizinhos dos elementos                               *
* nen       -> numero de nos por celulas                            *
* nFace     -> numero de faces por celulas                          *
* calType   -> tipo de calculo das celulas                          *
* geomType  -> tipo geometrico das celulas                          *
* prop      -> propriedades dos material                            *
* mat       -> material por celula                                  *
* cc        -> centroide das celulas                                *
* gKsi      -> vetores que unem centroide da celula central aos     *
*              vizinhos destas                                      *
* gmKsi     -> modulo do vetor ksi                                  *
* gEta      -> vetores paralelos as faces das celulas               *
* gfArea    -> modulo do vetor eta                                  *
* gNormal   -> vetores normais as faces das celulas                 *
* gVolume   -> volumes das celulas                                  *
* gXm       -> pontos medios das faces das celulas                  *
* gXmcc     -> vetores que unem o centroide aos pontos medios das   *
*              faces                                                *
* gvSkew    -> vetor entre o ponto medio a intersecao que une os    *
*              centrois compartilhado nessa face                    *
* gmvSkew   -> distacia entre o ponto medio a intersecao que une os *
*              centrois compartilhado nessa face                    *
* gDcca     -> menor distancia do centroide a faces desta celula    *
* density   -> massa especifica com variacao temporal               *
* id        -> numera das equacoes                                  *
* faceVelR  -> restricoes por elemento de velocidades               *
* faceVelL  -> carga por elemento de velocidades                    *
* facePresR -> restricoes por elemento de pressao                   *
* facePresL -> carga por elemento de pressao                        *
* pres      -> campo de pressao conhecido                           *
* gradVel   -> gradiente da solucao conhecido                       *
* bT        -> parte da discretizacao temporal                      *
* dField    -> matriz D do metodo simple                            *
* vel       -> campo de velocidade conhecido                        *
* velUp     -> nao definido                                         *
* ddt       -> discretizacao temporal                               *
* maxNo     -> numero de nos por celula maximo da malha             *
* maxViz    -> numero vizinhos por celula maximo da malha           *
* ndm       -> numero de dimensoes                                  *
* numel     -> numero de toral de celulas                           *
* ndf       -> graus de liberdade                                   *
* sPressure -> reconstrucao de segunda ordem para pressoes nas      *
*              faces                                                *
* fResidual -> calculo do residuo                                   *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* velUp     -> velocidaes atualizada (explicito)                    *
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
void velExp(Loads *loadsVel        ,Loads *loadsPres
             ,Advection advVel           ,Diffusion diffVel
             ,INT    *RESTRICT el        ,INT    *RESTRICT nelcon
             ,short  *RESTRICT nen       ,short  *RESTRICT nFace
             ,short  *RESTRICT geomType  ,DOUBLE *RESTRICT prop
             ,short  *RESTRICT calType   ,short  *RESTRICT mat
             ,DOUBLE *RESTRICT cc        
             ,DOUBLE *RESTRICT gKsi      ,DOUBLE *RESTRICT gmKsi
             ,DOUBLE *RESTRICT gEta      ,DOUBLE *RESTRICT gfArea
             ,DOUBLE *RESTRICT gNormal   ,DOUBLE *RESTRICT gVolume
             ,DOUBLE *RESTRICT gXm       ,DOUBLE *RESTRICT gXmcc
             ,DOUBLE *RESTRICT gvSkew    ,DOUBLE *RESTRICT gmvSkew
             ,DOUBLE *RESTRICT gDcca     ,DOUBLE *RESTRICT density
             ,short  *RESTRICT faceVelR  ,short  *RESTRICT faceVelL
             ,short  *RESTRICT facePresR ,short  *RESTRICT facePresL
             ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres
             ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT velUp
             ,DOUBLE *RESTRICT gradVel   ,DOUBLE *RESTRICT bT
             ,DOUBLE *RESTRICT dField    ,DOUBLE underU
             ,Temporal ddt
             ,short maxNo                ,short maxViz
             ,short ndm                  ,INT numel
             ,short ndf                  ,bool sPressure
             ,bool fResidual)
{
  short nThreads = ompVar.nThreadsCell;
  INT nel, vizNel;
  short i, j, k;
/*... variavel local*/
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM],lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM],lVolume[MAX_NUM_FACE + 1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM],lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE], lvSkew[MAX_NUM_FACE*MAX_NDM];
  short  lGeomType[MAX_NUM_FACE+1];
  short  lFaceVelR[MAX_NUM_FACE+1],lFacePresR[MAX_NUM_FACE+1];
  short  lFaceVelL[MAX_NUM_FACE+1],lFacePresL[MAX_NUM_FACE+1];
  DOUBLE lDensity[(MAX_NUM_FACE+1)];
  DOUBLE lProp[(MAX_NUM_FACE+1)*MAXPROP];
  DOUBLE lPres[(MAX_NUM_FACE+1)],lB[MAX_NDM],lBt[MAX_NDM];
  DOUBLE lGradVel[(MAX_NUM_FACE+1)*MAX_NDM*MAX_NDF];
  DOUBLE lGradPres[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lVel[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lCc[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lDfield[MAX_NDM];
  INT    lViz[MAX_NUM_FACE];
  short  aux1, aux2, lMat;

/*...*/
  if (ompVar.fCell) {
/*... loop nas celulas*/
    aux2 = maxViz + 1;
#pragma omp parallel  for default(none) num_threads(nThreads)\
     private(nel,i,j,k,aux1,lPres,lMat,lVolume,lGeomType\
          ,lFaceVelR,lDfield,lB,lBt\
          ,lFaceVelL,lFacePresR,lFacePresL,lDensity,lProp,lGradPres\
          ,lVel,lCc,lGradVel,lmKsi,lfArea,lDcca,lmvSkew,lKsi,lEta\
          ,lNormal,lXm,lXmcc,lvSkew,vizNel,lViz)\
     shared(aux2,ndm,ndf,numel,maxViz,nFace,mat,velUp\
         ,calType,gVolume, geomType, faceVelR, faceVelL, facePresR\
         ,facePresL,density, prop,gradPres,vel,cc,gradVel,gmKsi\
         ,gfArea,gDcca,gmvSkew,gKsi,gEta,gNormal,gXm,gXmcc,gvSkew\
         ,nelcon,loadsVel,loadsPres,advVel,diffVel,bT\
         ,ddt,sPressure,nen,pres,dField,underU,fResidual) 
    for (nel = 0; nel<numel; nel++) {
/*...*/
      aux1 = nFace[nel];
/*... elementos com equacoes*/
      if (MAT2D(nel,aux1,faceVelR,aux2) != PCCELL) {

        for(j=0;j<MAX_NUM_FACE+1;j++)
          lPres[j] = 0.e0;

/*... loop na celula central*/
        lMat             = mat[nel] - 1;
//      lib              = calType[lMat];
        lVolume[aux1]    = gVolume[nel];
        lGeomType[aux1]  = geomType[nel];
        lFaceVelR[aux1]  = MAT2D(nel, aux1, faceVelR, aux2);
        lFaceVelL[aux1]  = MAT2D(nel, aux1, faceVelL, aux2);
        lFacePresR[aux1] = MAT2D(nel, aux1, facePresR, aux2);
        lFacePresL[aux1] = MAT2D(nel, aux1, facePresL, aux2);
        lDensity[aux1]   = MAT2D(nel, 0, density, DENSITY_LEVEL);

        lPres[aux1] = pres[nel];

        for (j = 0; j<MAXPROP; j++)
          MAT2D(aux1, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);

        for (j = 0; j<ndm; j++) {
          MAT2D(aux1,j,lGradPres,ndm) = MAT2D(nel,j,gradPres,ndm);
          MAT2D(aux1,j,lVel     ,ndm) = MAT2D(nel,j,vel     ,ndm);
          MAT2D(aux1,j,lCc      ,ndm) = MAT2D(nel,j,cc      ,ndm);
          lBt[j]                      = MAT2D(nel,j,bT      ,ndm);
        }

        for (i = 0; i<ndf; i++)
          for (j = 0; j<ndm; j++)
            MAT3D(aux1,i,j,lGradVel,ndf,ndm)
            = MAT3D(nel,i,j,gradVel,ndf,ndm);

        for (i = 0; i<aux1; i++) {
          lmKsi[i]      = MAT2D(nel,i,gmKsi    ,maxViz);
          lfArea[i]     = MAT2D(nel,i,gfArea   ,maxViz);
          lDcca[i]      = MAT2D(nel,i,gDcca    ,maxViz);
          lmvSkew[i]    = MAT2D(nel,i,gmvSkew  ,maxViz);
          lFaceVelR[i]  = MAT2D(nel,i,faceVelR ,aux2);
          lFaceVelL[i]  = MAT2D(nel,i,faceVelL ,aux2);
          lFacePresR[i] = MAT2D(nel,i,facePresR,aux2);
          lFacePresL[i] = MAT2D(nel,i,facePresL,aux2);
          for (j = 0; j<ndm; j++) {
            MAT2D(i,j,lKsi   ,ndm) = MAT3D(nel,i,j,gKsi   ,maxViz,ndm);
            MAT2D(i,j,lEta   ,ndm) = MAT3D(nel,i,j,gEta   ,maxViz,ndm);
            MAT2D(i,j,lNormal,ndm) = MAT3D(nel,i,j,gNormal,maxViz,ndm);
            MAT2D(i,j,lXm    ,ndm) = MAT3D(nel,i,j,gXm    ,maxViz,ndm);
            MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,gXmcc  ,maxViz,ndm);
            MAT2D(i,j,lvSkew ,ndm) = MAT3D(nel,i,j,gvSkew ,maxViz,ndm);
          }
        }

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++) {
          vizNel = MAT2D(nel,i,nelcon,maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2) {
            lVolume[i]   = gVolume[vizNel];
            lGeomType[i] = geomType[vizNel];
            lDensity[i]  = MAT2D(vizNel,0,density,DENSITY_LEVEL);
            lMat = mat[vizNel] - 1;
            lPres[i] = pres[vizNel];

            for (j = 0; j<ndm; j++) {
              MAT2D(i,j,lGradPres,ndm) = MAT2D(vizNel,j,gradPres,ndm);
              MAT2D(i,j,lVel     ,ndm) = MAT2D(vizNel,j,vel     ,ndm);
              MAT2D(i,j,lCc      ,ndm) = MAT2D(vizNel,j,cc      ,ndm);
            }

            for (k = 0; k<ndf; k++)
              for (j = 0; j<ndm; j++)
                MAT3D(i,k,j,lGradVel,ndf,ndm)
                = MAT3D(vizNel,k,j,gradVel,ndf,ndm);

            for (j = 0; j<DIFPROP; j++)
              MAT2D(i,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
          }
        }
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibVelExp(loadsVel  ,loadsPres
                     ,advVel    ,diffVel
                     ,lGeomType ,lProp
                     ,lViz
                     ,lKsi      ,lmKsi
                     ,lEta      ,lfArea
                     ,lNormal   ,lVolume
                     ,lXm       ,lXmcc
                     ,lDcca     ,lDensity
                     ,lvSkew    ,lmvSkew
                     ,lB        ,ddt
                     ,lFaceVelR ,lFaceVelL
                     ,lFacePresR,lFacePresL
                     ,lPres     ,lGradPres
                     ,lVel      ,lGradVel
                     ,lDfield   ,lCc
                     ,lBt       ,underU
                     ,sPressure ,fResidual
                     ,nen[nel]  ,nFace[nel]
                     ,ndm       ,nel);
/*...................................................................*/

/*...*/
        for (j = 0; j<ndm; j++) {
          MAT2D(nel,j,velUp ,ndm) = lB[j];
          MAT2D(nel,j,dField,ndm) = lDfield[j];
        }
/*...................................................................*/

      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... sequencial*/
  else {
/*... loop nas celulas*/
    aux2 = maxViz + 1;
    for (nel = 0; nel<numel; nel++) {
 
/*...*/
      aux1 = nFace[nel];
      /*... elementos com equacoes*/
      if (MAT2D(nel, aux1, faceVelR, aux2) != PCCELL) {

/*... zerando vetores*/
        for (j = 0; j<MAX_NUM_FACE + 1; j++)
          lPres[j] = 0.e0;

/*... loop na celula central*/
        lMat             = mat[nel] - 1;
//      lib              = calType[lMat];
        lVolume[aux1]    = gVolume[nel];
        lGeomType[aux1]  = geomType[nel];
        lFaceVelR[aux1]  = MAT2D(nel,aux1,faceVelR ,aux2);
        lFaceVelL[aux1]  = MAT2D(nel,aux1,faceVelL ,aux2);
        lFacePresR[aux1] = MAT2D(nel,aux1,facePresR,aux2);
        lFacePresL[aux1] = MAT2D(nel,aux1,facePresL,aux2);
        lDensity[aux1]   = MAT2D(nel,0   ,density  ,DENSITY_LEVEL);

        lPres[aux1] = pres[nel];

        for (j = 0; j<MAXPROP; j++)
          MAT2D(aux1,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);

        for (j = 0; j<ndm; j++) {
          MAT2D(aux1,j,lGradPres,ndm) = MAT2D(nel,j,gradPres,ndm);
          MAT2D(aux1,j,lVel     ,ndm) = MAT2D(nel,j,vel     ,ndm);
          MAT2D(aux1,j,lCc      ,ndm) = MAT2D(nel,j,cc      ,ndm);
          lBt[j]                      = MAT2D(nel,j,bT      ,ndm);
        }

        for (i = 0; i<ndf; i++)
          for (j = 0; j<ndm; j++)
            MAT3D(aux1,i,j,lGradVel,ndf,ndm)
            = MAT3D(nel,i,j,gradVel,ndf,ndm);

        for (i = 0; i<aux1; i++) {
          lmKsi[i]      = MAT2D(nel,i,gmKsi    ,maxViz);
          lfArea[i]     = MAT2D(nel,i,gfArea   ,maxViz);
          lDcca[i]      = MAT2D(nel,i,gDcca    ,maxViz);
          lmvSkew[i]    = MAT2D(nel,i,gmvSkew  ,maxViz);
          lFaceVelR[i]  = MAT2D(nel,i,faceVelR ,aux2);
          lFaceVelL[i]  = MAT2D(nel,i,faceVelL ,aux2);
          lFacePresR[i] = MAT2D(nel,i,facePresR,aux2);
          lFacePresL[i] = MAT2D(nel,i,facePresL,aux2);
          for (j = 0; j<ndm; j++) {
            MAT2D(i,j,lKsi   ,ndm) = MAT3D(nel,i,j,gKsi   ,maxViz,ndm);
            MAT2D(i,j,lEta   ,ndm) = MAT3D(nel,i,j,gEta   ,maxViz,ndm);
            MAT2D(i,j,lNormal,ndm) = MAT3D(nel,i,j,gNormal,maxViz,ndm);
            MAT2D(i,j,lXm    ,ndm) = MAT3D(nel,i,j,gXm    ,maxViz,ndm);
            MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,gXmcc  ,maxViz,ndm);
            MAT2D(i,j,lvSkew ,ndm) = MAT3D(nel,i,j,gvSkew ,maxViz,ndm);
          }
        }

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++) {
          vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2) {
            lVolume[i]   = gVolume[vizNel];
            lGeomType[i] = geomType[vizNel];
            lDensity[i]  = MAT2D(vizNel, 0, density, DENSITY_LEVEL);
            lMat         = mat[vizNel] - 1;
            lPres[i]     = pres[vizNel];

            for (j = 0; j<ndm; j++) {
              MAT2D(i,j,lGradPres,ndm) = MAT2D(vizNel,j,gradPres,ndm);
              MAT2D(i,j,lVel     ,ndm) = MAT2D(vizNel,j,vel     ,ndm);
              MAT2D(i,j,lCc      ,ndm) = MAT2D(vizNel,j,cc      ,ndm);
            }

            for (k = 0; k<ndf; k++)
              for (j = 0; j<ndm; j++)
                MAT3D(i, k, j, lGradVel, ndf, ndm)
                = MAT3D(vizNel, k, j, gradVel, ndf, ndm);

            for (j = 0; j<DIFPROP; j++)
              MAT2D(i,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
          }
        }
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibVelExp(loadsVel  ,loadsPres
                     ,advVel    ,diffVel
                     ,lGeomType ,lProp
                     ,lViz     
                     ,lKsi      ,lmKsi
                     ,lEta      ,lfArea
                     ,lNormal   ,lVolume
                     ,lXm       ,lXmcc
                     ,lDcca     ,lDensity
                     ,lvSkew    ,lmvSkew
                     ,lB        ,ddt
                     ,lFaceVelR ,lFaceVelL
                     ,lFacePresR,lFacePresL
                     ,lPres     ,lGradPres
                     ,lVel      ,lGradVel
                     ,lDfield   ,lCc
                     ,lBt       ,underU
                     ,sPressure ,fResidual
                     ,nen[nel]  ,nFace[nel]
                     ,ndm       ,nel);
/*...................................................................*/

/*...*/
        for (j = 0; j<ndm; j++){
          MAT2D(nel,j,velUp ,ndm) = lB[j];
          MAT2D(nel,j,dField,ndm) = lDfield[j];
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 22/08/2017                                   *
 * Data de modificaco : 27/08/2019                                   *
 *-------------------------------------------------------------------*
 * SYSTFOMENERGY: calculo do sistema de equacoes para problemas      *
 * transporte de energia (Ax=b)                                      *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loads   -> definicoes de cargas                                   *
 * eModel   -> modelo da equacao de energia                          *
 * tModel  -> modelp de turbulencia                                  *
 * cModel  -> modelo de combustao                                    *
 * advT    -> tecnica da discretizacao do termo advecao              *
 * diffT   -> tecnica da discretizacao do termo difusivo             *
 * vProp     -> propedades variaveis (true|false)                    *
 * el      -> conetividade dos celulas                               *
 * nelcon  -> vizinhos dos elementos                                 *
 * nen     -> numero de nos por celulas                              *
 * nFace   -> numero de faces por celulas                            *
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * gVolume -> volumes das celulas                                    *
 * gDcca   -> menor distancia do centroide a faces desta celula      *
 * gXmCc   -> vetores que unem o centroide aos pontos medios das     *
 * gCc     -> centroide do elemento                                  * 
 * fModKsi -> o mudolo do vetores que unem centroide da celula       * 
 *            central aos vizinhos destas                            *
 * fKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        * 
 * fEta    -> vetores paralelos as faces das celulas                 *
 * fArea   -> modulo do vetor eta                                    *
 * fNormal -> vetores normais as faces das celulas                   *
 * fXm     -> pontos medios das faces das celulas                    *
 * fModvSkew -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face                      * 
 * fvSkew  -> vetor entre o ponto medio a intersecao que une os      *
 *            centrois compartilhado nessa face                      *
 * calType -> tipo de calculo das celulas                            * 
 * geomType-> tipo geometrico das celulas                            * 
 * prop    -> propriedades dos material                              * 
 * mat     -> material por celula                                    * 
 * rateHeatComb -> liberao de energia por combustao                  *
 * au      -> matriz de coeficientes esparsa                         *
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- triangular superior               )       *
 * gDcca   -> menor distancia do centroide a faces desta celula      *
 * ia      -> ponteiro para as linhas da matriz esparsa              *
 * ja      -> ponteiro para as colunas da matriz esparsa             *
 * a       -> matriz de coeficientes esparsa                         *
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- fora da diagonal principal        )       *
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 * b       -> vetor de forcas                                        *
 * id      -> numera das equacoes                                    *
 * faceR   -> restricoes por elemento                                *
 * faceLd1 -> carga por elemento                                     *
 * u0      -> solucao conhecida                                      *
 * gradU0  -> gradiente da solucao conhecido                         *
 * vel       -> campo de velocidade conhecido                        *
 * gradVel   -> gradiente rescontruido da velocidade                 *
 * pres0     -> pressao do tempo anterior                            *
 * pres      -> pressao do tempo atual                               *
 * gradPres  -> gradiente de pressao do tempo atual                  *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * rCell   -> nao definido                                           *
 * density -> massa especifica com variacao temporal                 *
 * sHeat    -> calor especifico com variacao temporal                *
 * dViscosity-> viscosidade dinamica com variacao temporal           *
 * tConductvity -> condutividade termica com variacao temporal       *
 * enthalpyk -> entalpia sensivel por especie                        *
 * gradY     -> gradiente das especies                               *
 * diffY  -> coeficiente de difusao das especie                      *
 * yFrac  -> fracao massica
 * rateHeatReComb -> taxa de liberacao de energia                    * 
 * dField    -> matriz D do metodo simple                            * 
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  *
 * ddt     -> discretizacao temporal                                 *
 * underU  -> parametro de sob relaxamento                           *
 * nEq     -> numero de equacoes                                     *
 * neqNov  -> numero de equacoes nao sobrepostas                     *
 * nAd     -> numero de termos nao nulos                             *
 * nAdR    -> numero de termos nao nulos na parte retangular         *
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
 * rCell     -> residuo por celula                                   *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void systFormEnergy(Loads *loads       , Loads *ldVel  
       , Advection *adv                , Diffusion *diff 
       , Turbulence *tModel            , EnergyModel *eModel 
       , Combustion *cModel            , PropVarFluid *vProp           
       , INT    *RESTRICT el           , INT    *RESTRICT nelcon
       , short  *RESTRICT nen          , short  *RESTRICT nFace
       , INT *RESTRICT cellFace        , INT *RESTRICT fOwner
       , DOUBLE *RESTRICT gVolume      , DOUBLE *RESTRICT gDcca
       , DOUBLE *RESTRICT gXmCc        , DOUBLE *RESTRICT gCc
       , DOUBLE *RESTRICT fModKsi      , DOUBLE *RESTRICT fKsi
       , DOUBLE *RESTRICT fEta         , DOUBLE *RESTRICT fArea
       , DOUBLE *RESTRICT fNormal      , DOUBLE *RESTRICT fXm
       , DOUBLE *RESTRICT fModvSkew    , DOUBLE *RESTRICT fvSkew
       , short  *RESTRICT geomType     
       , short  *RESTRICT calType      , short  *RESTRICT mat
       , INT    *RESTRICT ia           , INT    *RESTRICT ja
       , DOUBLE *RESTRICT a            , DOUBLE *RESTRICT ad
       , DOUBLE *RESTRICT b            , INT    *RESTRICT id
       , short  *RESTRICT faceR        , short  *RESTRICT faceL
       , short  *RESTRICT faceVelR     , short  *RESTRICT faceVelL
       , DOUBLE *RESTRICT u0           , DOUBLE *RESTRICT gradU0
       , DOUBLE *RESTRICT vel          , DOUBLE *RESTRICT gradVel
       , DOUBLE *RESTRICT pres0        , DOUBLE *RESTRICT pres 
       , DOUBLE *RESTRICT gradPres     , DOUBLE *RESTRICT rCell
       , DOUBLE *RESTRICT density      , DOUBLE *RESTRICT sHeat
       , DOUBLE *RESTRICT dViscosity   , DOUBLE *RESTRICT eddyViscosity
       , DOUBLE *RESTRICT tConductivity
       , DOUBLE *RESTRICT enthalpyk    , DOUBLE *RESTRICT gradY 
       , DOUBLE *RESTRICT diffY        , DOUBLE *RESTRICT yFrac
       , DOUBLE *RESTRICT rateHeatComb 
       , DOUBLE *RESTRICT dField       , DOUBLE *RESTRICT wallPar
       , Temporal *ddt                 , DOUBLE underU
       , INT nEq                       , INT nEqNov
       , INT nAd                       , INT nAdR
       , short maxNo                   , short maxViz
       , short ndm                     , INT numel
       , short ndf                     , short storage
       , bool forces                   , bool matrix
       , bool calRcell                 , bool unsym)
{
  bool fTurb = tModel->fTurb
      ,fWallModel  = tModel->fWall;
  short i, j, k, lib, aux1, aux2, lMat;
  short nThreads = ompVar.nThreadsCell, ns = cModel->nOfSpecies;
  INT nel, vizNel;
  
/*... variavel local */
  short  lGeomType[MAX_NUM_FACE + 1],
          lFaceR[MAX_NUM_FACE + 1], lFaceL[MAX_NUM_FACE + 1],
          lFaceVelR[MAX_NUM_FACE + 1], lFaceVelL[MAX_NUM_FACE + 1];
  INT    lId[(MAX_NUM_FACE + 1)*MAX_NDF], lViz[MAX_NUM_FACE];
  INT    idFace, cellOwner, ch;
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM], lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM], lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM], lVolume[MAX_NUM_FACE + 1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM], lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE], lvSkew[MAX_NUM_FACE*MAX_NDM];  
  DOUBLE lDensity[(MAX_NUM_FACE + 1)],lsHeat[(MAX_NUM_FACE + 1)],
         lViscosity[(MAX_NUM_FACE + 1)*2],ltConductivity[(MAX_NUM_FACE + 1)];
  DOUBLE lA[(MAX_NUM_FACE + 1)*MAX_NDF], lB[MAX_NDF];
  DOUBLE lu0[(MAX_NUM_FACE + 1)*MAX_NDF],lPres[(MAX_NUM_FACE + 1)*2];
  DOUBLE lGradVel[(MAX_NUM_FACE + 1)*MAX_NDM*MAX_NDF];
  DOUBLE lGradU0[(MAX_NUM_FACE + 1)*MAX_NDM],lDfield[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lGradPres[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lVel[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lCc[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lRcell[MAX_NDF],lWallPar[NWALLPAR];
  DOUBLE lRateHeatC,lEnthalpyk[(MAX_NUM_FACE + 1)*MAXSPECIES]
        ,lGradY[(MAX_NUM_FACE + 1)*MAX_NDM*MAXSPECIES]
        ,lDiffY[(MAX_NUM_FACE + 1)*MAXSPECIES]
        ,lYfrac[(MAX_NUM_FACE + 1)*MAXSPECIES];

/*...*/
  if (ompVar.fCell)
  { 
/*... loop nas celulas*/
    aux2 = maxViz + 1;
#pragma omp parallel  for default(none) num_threads(nThreads)\
     private(nel,i,j,k,aux1,lId,lMat,lib,lVolume,lGeomType,lu0\
          ,lFaceR,lA,lB,lFaceL,lDensity,lGradU0,lDfield\
          ,lVel,lCc,lmKsi,lfArea,lDcca,lmvSkew,lKsi,lEta,lViscosity\
          ,lNormal,lXm,lXmcc,lvSkew,vizNel,lViz,lRcell,lPres,lGradPres\
          ,lFaceVelR,lFaceVelL,lsHeat,ltConductivity,lRateHeatC\
          ,lGradVel,lWallPar,idFace,cellOwner,ch\
          ,lYfrac,lDiffY,lEnthalpyk,lGradY)\
     shared(aux2,ndm,ndf,numel,maxViz,calRcell,rCell,nFace,mat\
         ,calType,gVolume,geomType,faceR,faceVelR,faceL,faceVelL\
         ,u0,gCc,loads,ldVel,dField\
         ,density,gradU0,vel,fModKsi,dViscosity,eddyViscosity\
         ,fArea,gDcca,fModvSkew,fKsi,fEta,fNormal,fXm,gXmCc,fvSkew\
         ,nelcon,id,adv,diff,sHeat,tConductivity,pres0,pres,gradPres\
         ,ddt,nen,ia,ja,a,ad,b,nEq,nEqNov,nAd,gradVel,underU \
         ,nAdR,storage,forces,matrix,unsym,tModel,eModel,cModel,vProp,wallPar\
         ,fOwner,cellFace,enthalpyk,gradY,diffY,yFrac,rateHeatComb,ns\
         ,fTurb,fWallModel)
    for (nel = 0; nel<numel; nel++) {
/*...*/
      if (calRcell)
        for (j = 0; j<ndf; j++)
          MAT2D(nel, j, rCell, ndf) = 0.e0;;
/*...*/
      aux1 = nFace[nel];
/*... elementos com equacoes*/
      if (MAT2D(nel, aux1, faceR, aux2) != PCCELL) {

/*... zerando vetores*/
        for (j = 0; j<(MAX_NUM_FACE + 1)*MAX_NDF; j++) {
          lId[j] = -1;
          lu0[j] = 0.e0;
        }
/*... loop na celula central*/
        lMat            = mat[nel] - 1;
        lib             = calType[lMat];
        lVolume[aux1]   = gVolume[nel];
        lGeomType[aux1] = geomType[nel];
        lFaceR[aux1]    = MAT2D(nel, aux1, faceR, aux2);
        lFaceL[aux1]    = MAT2D(nel, aux1, faceL, aux2);
        lFaceVelR[aux1] = MAT2D(nel, aux1, faceVelR, aux2);
        lFaceVelL[aux1] = MAT2D(nel, aux1, faceVelL, aux2);
        lDensity[aux1]  = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
        lsHeat[aux1]    = MAT2D(nel, TIME_N, sHeat  , SHEAT_LEVEL); 
        ltConductivity[aux1] = tConductivity[nel];
/*...*/
        MAT2D(aux1, 0, lPres, 2) = pres0[nel];
        MAT2D(aux1, 1, lPres, 2) = pres[nel];
/*...................................................................*/

/*... viscosidade dinamica e turbulentea*/
        MAT2D(aux1, 0, lViscosity, 2) = dViscosity[nel];
        if(fTurb)
          MAT2D(aux1, 1, lViscosity, 2) = eddyViscosity[nel];
/*...................................................................*/

/*...*/
        for (j = 0; j<ndf; j++) {
          MAT2D(aux1, j, lu0, ndf) = MAT2D(nel, j, u0, ndf);
          MAT2D(aux1, j, lId, ndf) = MAT2D(nel, j, id, ndf) - 1;
        }
/*...................................................................*/

/*...*/
        for (j = 0; j<ndm; j++) {
          MAT2D(aux1, j, lGradU0, ndm) = MAT2D(nel, j, gradU0, ndm);
          MAT2D(aux1, j, lVel, ndm) = MAT2D(nel, j, vel, ndm);
          MAT2D(aux1, j, lCc, ndm) = MAT2D(nel, j, gCc, ndm);
          MAT2D(aux1, j, lGradPres, ndm) = MAT2D(nel, j, gradPres, ndm);
          MAT2D(aux1, j, lDfield, ndm) = MAT2D(nel, j, dField, ndm);
        }
/*...................................................................*/

/*...*/
        for(i=0;i<ndm;i++)
          for(j=0;j<ndm;j++)
            MAT3D(aux1,i,j,lGradVel,ndm,ndm) 
                               = MAT3D(nel,i,j,gradVel,ndm,ndm);
/*...................................................................*/

/*...*/
        if(cModel->fCombustion)
        {
          for(j=0;j<ns;j++)
          {
/*...*/
            MAT2D(aux1, j, lYfrac, ns) = MAT2D(nel, j, yFrac, ns); 
/*...*/
            MAT2D(aux1, j, lDiffY, ns) = MAT2D(nel, j, diffY, ns); 
/*... entalpia sensivel por especies*/
            MAT2D(aux1, j, lEnthalpyk, ns) =
                    MAT2D(nel, j, enthalpyk, ns); 
/*...*/
                for (k = 0; k<ndm; k++)
                  MAT3D(aux1, j, k, lGradY, ns, ndm) =
                         MAT3D(nel, j, k, gradY, ns, ndm); 
          }
        }  
/*...................................................................*/

/*...*/
        for (i = 0; i<aux1; i++)
        {
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lFaceR[i] = MAT2D(nel, i, faceR, aux2);
          lFaceL[i] = MAT2D(nel, i, faceL, aux2);
          lFaceVelR[i] = MAT2D(nel, i, faceVelR, aux2);
          lFaceVelL[i] = MAT2D(nel, i, faceVelL, aux2);
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
        }
/*...................................................................*/

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++) {
          vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2) {
            lVolume[i]       = gVolume[vizNel];
            lGeomType[i]     = geomType[vizNel];
            lDensity[i]      = MAT2D(vizNel, 2, density, DENSITY_LEVEL);
            lsHeat[i]        = MAT2D(vizNel, 2, sHeat  , SHEAT_LEVEL);
            ltConductivity[i]= tConductivity[vizNel];
/*...*/
            MAT2D(i, 0, lPres, 2) = pres0[vizNel];
            MAT2D(i, 1, lPres, 2) = pres[vizNel];
/*...................................................................*/

/*... viscusidade dinamica e turbulentea*/
            MAT2D(i, 0, lViscosity, 2) = dViscosity[vizNel];
            if(fTurb)
              MAT2D(i, 1, lViscosity, 2) = eddyViscosity[vizNel];

            lMat = mat[vizNel] - 1;
            for (j = 0; j<ndf; j++) {
              MAT2D(i, j, lu0, ndf) = MAT2D(vizNel, j, u0, ndf);
              MAT2D(i, j, lId, ndf) = MAT2D(vizNel, j, id, ndf) - 1;
            }
            for (j = 0; j<ndm; j++) {
              MAT2D(i, j, lGradU0, ndm) = MAT2D(vizNel, j, gradU0, ndm);
              MAT2D(i, j, lVel, ndm) = MAT2D(vizNel, j, vel, ndm);
              MAT2D(i, j, lCc, ndm) = MAT2D(vizNel, j, gCc, ndm);
              MAT2D(i, j, lGradPres, ndm) = MAT2D(vizNel, j, gradPres, ndm);
              MAT2D(i ,j, lDfield, ndm) = MAT2D(vizNel, j, dField, ndm);
            }

            for (k = 0; k<ndm; k++)
              for (j = 0; j<ndm; j++)
                MAT3D(i, k, j, lGradVel, ndm, ndm)
                = MAT3D(vizNel, k, j, gradVel, ndm, ndm);

/*...*/
            if(cModel->fCombustion)
            {
              for(j=0;j<ns;j++)
              {
/*...*/
                MAT2D(i, j, lYfrac, ns) = MAT2D(vizNel, j, yFrac, ns); 
/*...*/      
                MAT2D(i, j, lDiffY, ns) = MAT2D(vizNel, j, diffY, ns);
/*...*/
                MAT2D(i, j, lEnthalpyk, ns) =
                    MAT2D(vizNel, j, enthalpyk, ns); 
/*...*/
                for (k = 0; k<ndm; k++)
                  MAT3D(i, j, k, lGradY, ns, ndm) =
                         MAT3D(vizNel, j, k, gradY, ns, ndm); 
              }
            }  
/*...................................................................*/
          }
/*...................................................................*/
        }
/*...................................................................*/

/*...*/ 
        if(fWallModel)
          for(i=0;i<NWALLPAR;i++)
            lWallPar[i] = MAT2D(nel,i,wallPar,NWALLPAR);
/*...................................................................*/

/*...*/
        if(cModel->fCombustion)
        {
          lRateHeatC = rateHeatComb[nel];
        }  
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibEnergy(loads        , ldVel
                    , adv          , diff
                    , tModel       , eModel
                    , cModel       , vProp         
                    , lGeomType    
                    , lViz         , lId
                    , lKsi         , lmKsi
                    , lEta         , lfArea
                    , lNormal      , lVolume
                    , lXm          , lXmcc
                    , lDcca        , lCc
                    , lvSkew       , lmvSkew
                    , lA           , lB
                    , lRcell       , ddt
                    , lFaceR       , lFaceL
                    , lFaceVelR    , lFaceVelL
                    , lu0          , lGradU0
                    , lVel         , lGradVel
                    , lPres        , lGradPres
                    , lDensity     , lsHeat
                    , lViscosity   , ltConductivity
                    , lEnthalpyk   , lGradY
                    , lDiffY       , lYfrac   
                    , lRateHeatC
                    , lDfield      , lWallPar
                    , underU
                    , nen[nel]     , nFace[nel]
                    , ndm          , lib
                    , nel);
/*...................................................................*/

/*... residuo da celula*/
        if (calRcell)
          for (j = 0; j<ndf; j++)
            MAT2D(nel, j, rCell, ndf) = lRcell[j];
/*...................................................................*/

/*...*/
        assbly(ia        ,ja
              ,a         ,ad
              ,b
              ,lId
              ,lA        ,lB
              ,nEq       ,nEqNov
              ,nAd       ,nAdR
              ,nFace[nel],ndf
              ,storage   ,forces
              ,matrix    ,unsym);
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  else {
/*... loop nas celulas*/
    aux2 = maxViz + 1;
    for (nel = 0; nel<numel; nel++) {
/*...*/
      if (calRcell)
        for (j = 0; j<ndf; j++)
          MAT2D(nel, j, rCell, ndf) = 0.e0;;
/*...*/
      aux1 = nFace[nel];
/*... elementos com equacoes*/
      if (MAT2D(nel, aux1, faceR, aux2) != PCCELL) {

/*... zerando vetores*/
        for (j = 0; j<(MAX_NUM_FACE + 1)*MAX_NDF; j++) {
          lId[j] = -1;
          lu0[j] = 0.e0;
        }
/*... loop na celula central*/
        lMat            = mat[nel] - 1;
        lib             = calType[lMat];
        lVolume[aux1]   = gVolume[nel];
        lGeomType[aux1] = geomType[nel];
        lFaceR[aux1]    = MAT2D(nel, aux1, faceR, aux2);
        lFaceL[aux1]    = MAT2D(nel, aux1, faceL, aux2);
        lFaceVelR[aux1] = MAT2D(nel, aux1, faceVelR, aux2);
        lFaceVelL[aux1] = MAT2D(nel, aux1, faceVelL, aux2);
        lDensity[aux1]  = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
        lsHeat[aux1]    = MAT2D(nel, TIME_N, sHeat  , SHEAT_LEVEL); 
        ltConductivity[aux1] = tConductivity[nel];
/*...*/
        MAT2D(aux1, 0, lPres, 2) = pres0[nel];
        MAT2D(aux1, 1, lPres, 2) = pres[nel];
/*...................................................................*/

/*... viscosidade dinamica e turbulentea*/
        MAT2D(aux1, 0, lViscosity, 2) = dViscosity[nel];
        if(fTurb)
          MAT2D(aux1, 1, lViscosity, 2) = eddyViscosity[nel];
/*...................................................................*/

/*...*/
        for (j = 0; j<ndf; j++) {
          MAT2D(aux1, j, lu0, ndf) = MAT2D(nel, j, u0, ndf);
          MAT2D(aux1, j, lId, ndf) = MAT2D(nel, j, id, ndf) - 1;
        }
/*...................................................................*/

/*...*/
        for (j = 0; j<ndm; j++) {
          MAT2D(aux1, j, lGradU0, ndm) = MAT2D(nel, j, gradU0, ndm);
          MAT2D(aux1, j, lVel, ndm) = MAT2D(nel, j, vel, ndm);
          MAT2D(aux1, j, lCc, ndm) = MAT2D(nel, j, gCc, ndm);
          MAT2D(aux1, j, lGradPres, ndm) = MAT2D(nel, j, gradPres, ndm);
          MAT2D(aux1, j, lDfield, ndm) = MAT2D(nel, j, dField, ndm);
        }
/*...................................................................*/

/*...*/
        for(i=0;i<ndm;i++)
          for(j=0;j<ndm;j++)
            MAT3D(aux1,i,j,lGradVel,ndm,ndm) 
                               = MAT3D(nel,i,j,gradVel,ndm,ndm);
/*...................................................................*/

/*...*/
        if(cModel->fCombustion)
        {
          for(j=0;j<ns;j++)
          {
/*...*/
            MAT2D(aux1, j, lYfrac, ns) = MAT2D(nel, j, yFrac, ns); 
/*...*/
            MAT2D(aux1, j, lDiffY, ns) = MAT2D(nel, j, diffY, ns); 
/*... entalpia sensivel por especies*/
            MAT2D(aux1, j, lEnthalpyk, ns) =
                    MAT2D(nel, j, enthalpyk, ns); 
/*...*/
                for (k = 0; k<ndm; k++)
                  MAT3D(aux1, j, k, lGradY, ns, ndm) =
                         MAT3D(nel, j, k, gradY, ns, ndm); 
          }
        }  
/*...................................................................*/

/*...*/
        for (i = 0; i<aux1; i++)
        {
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lFaceR[i] = MAT2D(nel, i, faceR, aux2);
          lFaceL[i] = MAT2D(nel, i, faceL, aux2);
          lFaceVelR[i] = MAT2D(nel, i, faceVelR, aux2);
          lFaceVelL[i] = MAT2D(nel, i, faceVelL, aux2);
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
        }
/*...................................................................*/

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++) {
          vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2) {
            lVolume[i]       = gVolume[vizNel];
            lGeomType[i]     = geomType[vizNel];
            lDensity[i]      = MAT2D(vizNel, 2, density, DENSITY_LEVEL);
            lsHeat[i]        = MAT2D(vizNel, 2, sHeat  , SHEAT_LEVEL);
            ltConductivity[i]= tConductivity[vizNel];
/*...*/
            MAT2D(i, 0, lPres, 2) = pres0[vizNel];
            MAT2D(i, 1, lPres, 2) = pres[vizNel];
/*...................................................................*/

/*... viscusidade dinamica e turbulentea*/
            MAT2D(i, 0, lViscosity, 2) = dViscosity[vizNel];
            if(fTurb)
              MAT2D(i, 1, lViscosity, 2) = eddyViscosity[vizNel];
/*...................................................................*/
            lMat = mat[vizNel] - 1;
            for (j = 0; j<ndf; j++) {
              MAT2D(i, j, lu0, ndf) = MAT2D(vizNel, j, u0, ndf);
              MAT2D(i, j, lId, ndf) = MAT2D(vizNel, j, id, ndf) - 1;
            }
            for (j = 0; j<ndm; j++) {
              MAT2D(i, j, lGradU0, ndm) = MAT2D(vizNel, j, gradU0, ndm);
              MAT2D(i, j, lVel, ndm) = MAT2D(vizNel, j, vel, ndm);
              MAT2D(i, j, lCc, ndm) = MAT2D(vizNel, j, gCc, ndm);
              MAT2D(i, j, lGradPres, ndm) = MAT2D(vizNel, j, gradPres, ndm);
              MAT2D(i ,j, lDfield, ndm) = MAT2D(vizNel, j, dField, ndm);
            }

            for (k = 0; k<ndm; k++)
              for (j = 0; j<ndm; j++)
                MAT3D(i, k, j, lGradVel, ndm, ndm)
                = MAT3D(vizNel, k, j, gradVel, ndm, ndm);

/*...*/
            if(cModel->fCombustion)
            {
              for(j=0;j<ns;j++)
              {
/*...*/
                MAT2D(i, j, lYfrac, ns) = MAT2D(vizNel, j, yFrac, ns); 
/*...*/      
                MAT2D(i, j, lDiffY, ns) = MAT2D(vizNel, j, diffY, ns);
/*...*/
                MAT2D(i, j, lEnthalpyk, ns) =
                    MAT2D(vizNel, j, enthalpyk, ns); 
/*...*/
                for (k = 0; k<ndm; k++)
                  MAT3D(i, j, k, lGradY, ns, ndm) =
                         MAT3D(vizNel, j, k, gradY, ns, ndm); 
              }
            }  
/*...................................................................*/
          }
/*...................................................................*/
        }
/*...................................................................*/

/*...*/    
        if(fWallModel)
          for(i=0;i<NWALLPAR;i++)
            lWallPar[i] = MAT2D(nel,i,wallPar,NWALLPAR);
/*...................................................................*/

/*...*/
        if(cModel->fCombustion)
        {
          lRateHeatC = rateHeatComb[nel];
        }  
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibEnergy(loads        , ldVel
                    , adv          , diff
                    , tModel       , eModel
                    , cModel       , vProp         
                    , lGeomType   
                    , lViz         , lId
                    , lKsi         , lmKsi
                    , lEta         , lfArea
                    , lNormal      , lVolume
                    , lXm          , lXmcc
                    , lDcca        , lCc
                    , lvSkew       , lmvSkew
                    , lA           , lB
                    , lRcell       , ddt
                    , lFaceR       , lFaceL
                    , lFaceVelR    , lFaceVelL
                    , lu0          , lGradU0
                    , lVel         , lGradVel
                    , lPres        , lGradPres
                    , lDensity     , lsHeat
                    , lViscosity   , ltConductivity
                    , lEnthalpyk   , lGradY
                    , lDiffY       , lYfrac   
                    , lRateHeatC
                    , lDfield      , lWallPar
                    , underU
                    , nen[nel]     , nFace[nel]
                    , ndm          , lib
                    , nel);
/*...................................................................*/

/*... residuo da celula*/
        if (calRcell)
          for (j = 0; j<ndf; j++)
            MAT2D(nel, j, rCell, ndf) = lRcell[j];
/*...................................................................*/

/*...*/
        assbly(ia        ,ja
              ,a         ,ad
              ,b
              ,lId
              ,lA        ,lB
              ,nEq       ,nEqNov
              ,nAd       ,nAdR
              ,nFace[nel],ndf
              ,storage   ,forces
              ,matrix    ,unsym);
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 04/08/2017                                   *
 * Data de modificaco : 25/07/2019                                   *
 *-------------------------------------------------------------------*
 * SYSTFOMCOMB  : calculo do sistema de equacoes para problemas      *
 * transporte de energia (Ax=b)                                      *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loads   -> definicoes de cargas                                   *
 * tModel  -> modelo de turbulencia                                  *
 * cModel  -> modelo de combustao                                    *
 * advT    -> tecnica da discretizacao do termo advecao              *
 * diffT   -> tecnica da discretizacao do termo difusivo             *
 * vProp     -> propedades variaveis (true|false)                    *
 * el      -> conetividade dos celulas                               *
 * nelcon  -> vizinhos dos elementos                                 *
 * nen     -> numero de nos por celulas                              *
 * nFace   -> numero de faces por celulas                            *
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * gVolume -> volumes das celulas                                    *
 * gDcca   -> menor distancia do centroide a faces desta celula      *
 * gXmCc   -> vetores que unem o centroide aos pontos medios das     *
 * gCc     -> centroide do elemento                                  *
 * fModKsi -> o mudolo do vetores que unem centroide da celula       *
 *            central aos vizinhos destas                            *
 * fKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        *
 * fEta    -> vetores paralelos as faces das celulas                 *
 * fArea   -> modulo do vetor eta                                    *
 * fNormal -> vetores normais as faces das celulas                   *
 * fXm     -> pontos medios das faces das celulas                    *
 * fModvSkew -> distacia entre o ponto medio a intersecao que une os *
 *            centrois compartilhado nessa face                      *
 * fvSkew  -> vetor entre o ponto medio a intersecao que une os      *
 *            centrois compartilhado nessa face                      *
 * calType -> tipo de calculo das celulas                            *
 * geomType-> tipo geometrico das celulas                            *
 * prop    -> propriedades dos material                              *
 * mat     -> material por celula                                    *
 * au      -> matriz de coeficientes esparsa                         *
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- triangular superior               )       *
 * gDcca   -> menor distancia do centroide a faces desta celula      *
 * ia      -> ponteiro para as linhas da matriz esparsa              *
 * ja      -> ponteiro para as colunas da matriz esparsa             *
 * a       -> matriz de coeficientes esparsa                         *
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- fora da diagonal principal        )       *
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 * b       -> vetor de forcas                                        *
 * id      -> numera das equacoes                                    *
 * faceR   -> restricoes por elemento                                *
 * faceLd1 -> carga por elemento                                     *
 * u0      -> solucao conhecida                                      *
 * gradU0  -> gradiente da solucao conhecido                         *
 * wk        -> taxa de consumo massico das especies kg/(m3 s)       * 
 * vel       -> campo de velocidade conhecido                        *
 * gradVel   -> gradiente rescontruido da velocidade                 *
 * pres0     -> pressao do tempo anterior                            *
 * pres      -> pressao do tempo atual                               *
 * gradPres  -> gradiente de pressao do tempo atual                  *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * rCell   -> nao definido                                           *
 * density -> massa especifica com variacao temporal                 *
 * cDiff   -> coeficiente de difusao das especies com variacao       * 
 *            temporal                                               *
 * eddyVis -> viscosidade turbulenta com variacao temporal           *
 * dField    -> matriz D do metodo simple                            *
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  *
 * ddt     -> discretizacao temporal                                 *
 * underU  -> parametro de sob relaxamento                           *
 * nEq     -> numero de equacoes                                     *
 * neqNov  -> numero de equacoes nao sobrepostas                     *
 * nAd     -> numero de termos nao nulos                             *
 * nAdR    -> numero de termos nao nulos na parte retangular         *
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
 * rCell     -> residuo por celula                                   *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * b     = | bx1 bx2 ... bxn by1 by2 ... byn bz1 bz2 ... bzn |       * 
 * rCell = | rx1 rx2 ... rxn ry1 ry2 ... ryn rz1 rz2 ... rzn |       * 
 ********************************************************************/
void systFormComb(Loads *loads              , Loads *ldVel
                , Advection *scAdv          , Diffusion *scDiff
                , Turbulence *tModel        , Combustion *cModel
                , PropVarFluid *vProp
                , INT    *RESTRICT el       , INT    *RESTRICT nelcon
                , short  *RESTRICT nen      , short  *RESTRICT nFace
                , INT *RESTRICT cellFace    , INT *RESTRICT fOwner
                , DOUBLE *RESTRICT gVolume  , DOUBLE *RESTRICT gDcca
                , DOUBLE *RESTRICT gXmCc    , DOUBLE *RESTRICT gCc
                , DOUBLE *RESTRICT fModKsi  , DOUBLE *RESTRICT fKsi
                , DOUBLE *RESTRICT fEta     , DOUBLE *RESTRICT fArea
                , DOUBLE *RESTRICT fNormal  , DOUBLE *RESTRICT fXm
                , DOUBLE *RESTRICT fModvSkew, DOUBLE *RESTRICT fvSkew
                , short  *RESTRICT geomType 
                , short  *RESTRICT calType  , short  *RESTRICT mat
                , INT    *RESTRICT ia       , INT    *RESTRICT ja
                , DOUBLE *RESTRICT a        , DOUBLE *RESTRICT ad
                , DOUBLE *RESTRICT b        , INT    *RESTRICT id
                , short  *RESTRICT faceR    , short  *RESTRICT faceL
                , short  *RESTRICT faceVelR , short  *RESTRICT faceVelL
                , DOUBLE *RESTRICT u0       , DOUBLE *RESTRICT gradU0
                , DOUBLE *RESTRICT wk       , DOUBLE *RESTRICT vel      
                , DOUBLE *RESTRICT pres0    , DOUBLE *RESTRICT pres
                , DOUBLE *RESTRICT gradPres , DOUBLE *RESTRICT rCell
                , DOUBLE *RESTRICT density  , DOUBLE *RESTRICT diff
                , DOUBLE *RESTRICT eddyVisc , DOUBLE *RESTRICT wallPar
                , DOUBLE *RESTRICT dField   
                , Temporal *ddt             , DOUBLE underU
                , INT nEq                   , INT nEqNov
                , INT nAd                   , INT nAdR
                , short maxNo               , short maxViz
                , short ndm                 , INT numel
                , short ndf                 , short storage
                , bool forces               , bool matrix
                , bool calRcell             , bool unsym)
{
  bool fTurb      = tModel->fTurb
      ,fWallModel = tModel->fWall;
  short i, j, k, lib, aux1, aux2, lMat;
  short nThreads = ompVar.nThreadsCell
      , nReac = cModel->chem.nReac
      , ns    = cModel->nOfSpecies;
  INT nel, vizNel;

  /*... variavel local */
  short  lGeomType[MAX_NUM_FACE + 1],
    lFaceR[MAX_NUM_FACE + 1], lFaceL[MAX_NUM_FACE + 1],
    lFaceVelR[MAX_NUM_FACE + 1], lFaceVelL[MAX_NUM_FACE + 1];
  INT    lId[(MAX_NUM_FACE + 1)*MAX_NDF], lViz[MAX_NUM_FACE];
  INT    idFace, cellOwner, ch;
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM], lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM], lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM], lVolume[MAX_NUM_FACE + 1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM], lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE], lvSkew[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDensity[(MAX_NUM_FACE + 1)], lDiff[(MAX_NUM_FACE + 1)*MAX_COMB];
  DOUBLE lViscosity[(MAX_NUM_FACE + 1)*2];
  DOUBLE lA[(MAX_NUM_FACE + 1)*MAX_COMB], lB[MAX_COMB];
  DOUBLE lProp[(MAX_NUM_FACE + 1)*MAXPROP];
  DOUBLE lu0[(MAX_NUM_FACE + 1)*MAX_NDF], lPres[(MAX_NUM_FACE + 1) * 2];
  DOUBLE lGradU0[(MAX_NUM_FACE + 1)*MAX_NDM*MAX_COMB];
  DOUBLE lDfield[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lGradPres[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lVel[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lCc[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lRcell[MAX_COMB], lWallPar[NWALLPAR], lWk[MAX_COMB];

  /*...*/
  if (ompVar.fCell)
  {
    /*... loop nas celulas*/
    aux2 = maxViz + 1;
#pragma omp parallel  for default(none) num_threads(nThreads)\
    private(nel,i,j,k,aux1,lId,lMat,lib,lVolume,lGeomType,lu0\
          ,lFaceR,lA,lB,lFaceL,lDensity,lGradU0,lDfield\
          ,lVel,lCc,lmKsi,lfArea,lDcca,lmvSkew,lKsi,lEta,lViscosity\
          ,lNormal,lXm,lXmcc,lvSkew,vizNel,lViz,lRcell,lPres,lGradPres\
          ,lFaceVelR,lFaceVelL\
          ,lWallPar,idFace,cellOwner,ch\
          ,lDiff,lWk)\
     shared(aux2,ndm,ndf,numel,maxViz,calRcell,rCell,nFace,mat\
         ,calType,gVolume,geomType,faceR,faceVelR,faceL,faceVelL\
         ,u0,gCc,loads,ldVel,dField\
         ,density,gradU0,vel,fModKsi,eddyVisc\
         ,fArea,gDcca,fModvSkew,fKsi,fEta,fNormal,fXm,gXmCc,fvSkew\
         ,nelcon,id,diff,pres0,pres,gradPres\
         ,ddt,nen,ia,ja,a,ad,b,nEq,nEqNov,nAd,underU \
         ,nAdR,storage,forces,matrix,unsym,tModel,cModel,vProp,wallPar\
         ,fOwner,cellFace\
         ,ns,nReac,wk,scAdv,scDiff,fTurb,fWallModel)
    for (nel = 0; nel<numel; nel++)
    {
/*...*/
      if (calRcell)
        for (j = 0; j<ndf; j++)
          rCell[j*numel+nel] = 0.e0;  
/*...*/
      aux1 = nFace[nel];
/*... elementos com equacoes*/
      if (MAT2D(nel, aux1, faceR, aux2) != PCCELL) 
      {

/*... zerando vetores*/
        for (j = 0; j<(MAX_NUM_FACE + 1)*MAX_NDF; j++) 
        {
          lId[j] = -1;
          lu0[j] = 0.e0;
        }
/*... loop na celula central*/
        lMat                = mat[nel] - 1;
        lib                 = calType[lMat];
        lVolume[aux1]       = gVolume[nel];
        lGeomType[aux1]     = geomType[nel];
        lFaceR[aux1]        = MAT2D(nel, aux1, faceR, aux2);
        lFaceL[aux1]        = MAT2D(nel, aux1, faceL, aux2);
        lFaceVelR[aux1]     = MAT2D(nel, aux1, faceVelR, aux2);
        lFaceVelL[aux1]     = MAT2D(nel, aux1, faceVelL, aux2);
        lDensity[aux1]      = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
/*...*/
        MAT2D(aux1, 0, lPres, 2) = pres0[nel];
        MAT2D(aux1, 1, lPres, 2) = pres[nel];
/*...................................................................*/

/*... viscosidade dinamica e turbulentea*/
        MAT2D(aux1, 0, lViscosity, 2) = 0.e0;
        if(fTurb)
          MAT2D(aux1, 1, lViscosity, 2) = eddyVisc[nel];
/*...................................................................*/

/*...*/
        lId[aux1] = id[nel] - 1;        
/*...................................................................*/

/*...*/
        for (j = 0; j<ndf; j++)
          MAT2D(aux1, j, lu0, ndf) = MAT2D(nel, j, u0, ndf);       
/*...................................................................*/

/*...*/        
        for (j = 0; j<ns; j++)
          MAT2D(aux1, j, lDiff, ns) = MAT2D(nel, j, diff, ns); 
/*...................................................................*/

/*...*/
        for (j = 0; j<ndm; j++)
        {
          MAT2D(aux1, j, lVel, ndm) = MAT2D(nel, j, vel, ndm);
          MAT2D(aux1, j, lCc, ndm) = MAT2D(nel, j, gCc, ndm);
          MAT2D(aux1, j, lGradPres, ndm) = MAT2D(nel, j, gradPres, ndm);
          MAT2D(aux1, j, lDfield, ndm) = MAT2D(nel, j, dField, ndm);
        }
/*...................................................................*/

/*...*/
        for (k = 0; k < ndf; k++)
          for (j = 0; j<ndm; j++)
            MAT3D(aux1, k, j, lGradU0, ndf, ndm) 
            = MAT3D(nel, k, j, gradU0, ndf, ndm);
/*...................................................................*/

/*...*/
        for (i = 0; i<aux1; i++)
        {
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lFaceR[i] = MAT2D(nel, i, faceR, aux2);
          lFaceL[i] = MAT2D(nel, i, faceL, aux2);
          lFaceVelR[i] = MAT2D(nel, i, faceVelR, aux2);
          lFaceVelL[i] = MAT2D(nel, i, faceVelL, aux2);
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
        }
/*...................................................................*/

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++) 
        {
          vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2) 
          {
            lVolume[i] = gVolume[vizNel];
            lGeomType[i] = geomType[vizNel];
            lDensity[i] = MAT2D(vizNel, TIME_N, density, DENSITY_LEVEL);
/*...*/
            MAT2D(i, 0, lPres, 2) = pres0[vizNel];
            MAT2D(i, 1, lPres, 2) = pres[vizNel];
/*... viscusidade dinamica e turbulentea*/
            MAT2D(i, 0, lViscosity, 2) = 0.e0;
            if(fTurb)
              MAT2D(i, 1, lViscosity, 2) = eddyVisc[vizNel];
/*...................................................................*/

/*...*/
            lId[i] = id[vizNel] - 1;
/*...................................................................*/

            lMat = mat[vizNel] - 1;
/*...*/  
            for (j = 0; j<ndf; j++)
              MAT2D(i, j, lu0, ndf) = MAT2D(vizNel, j, u0, ndf);
/*...................................................................*/

/*...*/        
            for (j = 0; j<ns; j++)
              MAT2D(i, j, lDiff, ns) = MAT2D(vizNel, j, diff, ns); 
/*...................................................................*/

/*...*/  
            for (j = 0; j<ndm; j++) 
            {
              MAT2D(i, j, lVel, ndm) = MAT2D(vizNel, j, vel, ndm);
              MAT2D(i, j, lCc, ndm) = MAT2D(vizNel, j, gCc, ndm);
              MAT2D(i, j, lGradPres, ndm) = MAT2D(vizNel, j, gradPres, ndm);
              MAT2D(i, j, lDfield, ndm) = MAT2D(vizNel, j, dField, ndm);
            }
/*...................................................................*/

/*...*/  
            for (k = 0; k<ndf; k++)
              for (j = 0; j<ndm; j++)
                MAT3D(i, k, j, lGradU0, ndf, ndm) 
               = MAT3D(vizNel, k, j, gradU0, ndf,ndm);
/*...................................................................*/
          }
        }
/*...................................................................*/

/*...*/
        if(fWallModel)
          for (i = 0; i<NWALLPAR; i++)
            lWallPar[i] = MAT2D(nel, i, wallPar, NWALLPAR);
/*...................................................................*/

/*...*/
        for (i = 0; i<ndf; i++)
          lWk[i] = MAT2D(nel,i,wk,ns); 
/*...................................................................*/  

/*... chamando a biblioteca de celulas*/
        cellLibCombustion(loads        , ldVel
                        , scAdv        , scDiff
                        , tModel       , cModel
                        , vProp        
                        , lGeomType    
                        , lViz         , lId
                        , lKsi         , lmKsi
                        , lEta         , lfArea
                        , lNormal      , lVolume
                        , lXm          , lXmcc
                        , lDcca        , lCc
                        , lvSkew       , lmvSkew
                        , lA           , lB
                        , lRcell       , ddt
                        , lFaceR       , lFaceL
                        , lFaceVelR    , lFaceVelL
                        , lu0          , lGradU0
                        , lWk          , lVel      
                        , lPres        , lGradPres
                        , lDensity     , lDiff
                        , lViscosity
                        , lDfield      , lWallPar
                        , underU
                        , nen[nel]     , nFace[nel]
                        , ndm          , lib
                        , nel);      
/*...................................................................*/

/*... residuo da celula*/
        if(calRcell)
          for (j = 0; j<ndf; j++)
            rCell[j*numel+nel] = lRcell[j];  
/*...................................................................*/
  
/*...*/
          assblyBlock(ia       , ja
                   , a         , ad
                   , b         , lId
                   , lA        , lB
                   , nEq       , nEqNov
                   , nAd       , nAdR
                   , nFace[nel], ndf  
                   , ndf       , ndf    
                   , storage   , forces
                   , matrix    , unsym);          
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  else {
/*... loop nas celulas*/
    aux2 = maxViz + 1;
    for (nel = 0; nel<numel; nel++) {
/*...*/
      if (calRcell)
        for (j = 0; j<ndf; j++)
          rCell[j*numel+nel] = 0.e0;  
/*...*/
      aux1 = nFace[nel];
/*... elementos com equacoes*/
      if (MAT2D(nel, aux1, faceR, aux2) != PCCELL) 
      {

/*... zerando vetores*/
        for (j = 0; j<(MAX_NUM_FACE + 1)*MAX_NDF; j++) 
        {
          lId[j] = -1;
          lu0[j] = 0.e0;
        }
/*... loop na celula central*/
        lMat                = mat[nel] - 1;
        lib                 = calType[lMat];
        lVolume[aux1]       = gVolume[nel];
        lGeomType[aux1]     = geomType[nel];
        lFaceR[aux1]        = MAT2D(nel, aux1, faceR, aux2);
        lFaceL[aux1]        = MAT2D(nel, aux1, faceL, aux2);
        lFaceVelR[aux1]     = MAT2D(nel, aux1, faceVelR, aux2);
        lFaceVelL[aux1]     = MAT2D(nel, aux1, faceVelL, aux2);
        lDensity[aux1]      = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
/*...*/
        MAT2D(aux1, 0, lPres, 2) = pres0[nel];
        MAT2D(aux1, 1, lPres, 2) = pres[nel];
/*...................................................................*/

/*... viscosidade dinamica e turbulentea*/
        MAT2D(aux1, 0, lViscosity, 2) = 0.e0;
        if(fTurb)
          MAT2D(aux1, 1, lViscosity, 2) = eddyVisc[nel];
/*...................................................................*/

/*...*/
        lId[aux1] = id[nel] - 1;        
/*...................................................................*/

/*...*/
        for (j = 0; j<ndf; j++)
          MAT2D(aux1, j, lu0, ndf) = MAT2D(nel, j, u0, ndf);       
/*...................................................................*/

/*...*/        
        for (j = 0; j<ns; j++)
          MAT2D(aux1, j, lDiff, ns) = MAT2D(nel, j, diff, ns); 
/*...................................................................*/

/*...*/
        for (j = 0; j<ndm; j++)
        {
          MAT2D(aux1, j, lVel, ndm) = MAT2D(nel, j, vel, ndm);
          MAT2D(aux1, j, lCc, ndm) = MAT2D(nel, j, gCc, ndm);
          MAT2D(aux1, j, lGradPres, ndm) = MAT2D(nel, j, gradPres, ndm);
          MAT2D(aux1, j, lDfield, ndm) = MAT2D(nel, j, dField, ndm);
        }
/*...................................................................*/

/*...*/
        for (k = 0; k < ndf; k++)
          for (j = 0; j<ndm; j++)
            MAT3D(aux1, k, j, lGradU0, ndf, ndm) 
            = MAT3D(nel, k, j, gradU0, ndf, ndm);
/*...................................................................*/

/*...*/
        for (i = 0; i<aux1; i++)
        {
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lFaceR[i] = MAT2D(nel, i, faceR, aux2);
          lFaceL[i] = MAT2D(nel, i, faceL, aux2);
          lFaceVelR[i] = MAT2D(nel, i, faceVelR, aux2);
          lFaceVelL[i] = MAT2D(nel, i, faceVelL, aux2);
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
        }
/*...................................................................*/

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++) 
        {
          vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2) 
          {
            lVolume[i] = gVolume[vizNel];
            lGeomType[i] = geomType[vizNel];
            lDensity[i] = MAT2D(vizNel, TIME_N, density, DENSITY_LEVEL);
/*...*/
            MAT2D(i, 0, lPres, 2) = pres0[vizNel];
            MAT2D(i, 1, lPres, 2) = pres[vizNel];
/*... viscusidade dinamica e turbulentea*/
            MAT2D(i, 0, lViscosity, 2) = 0.e0;
            if(fTurb)
              MAT2D(i, 1, lViscosity, 2) = eddyVisc[vizNel];
/*...................................................................*/

/*...*/
            lId[i] = id[vizNel] - 1;
/*...................................................................*/

            lMat = mat[vizNel] - 1;
/*...*/  
            for (j = 0; j<ndf; j++)
              MAT2D(i, j, lu0, ndf) = MAT2D(vizNel, j, u0, ndf);
/*...................................................................*/

/*...*/        
            for (j = 0; j<ns; j++)
              MAT2D(i, j, lDiff, ns) = MAT2D(vizNel, j, diff, ns); 
/*...................................................................*/

/*...*/  
            for (j = 0; j<ndm; j++) 
            {
              MAT2D(i, j, lVel, ndm) = MAT2D(vizNel, j, vel, ndm);
              MAT2D(i, j, lCc, ndm) = MAT2D(vizNel, j, gCc, ndm);
              MAT2D(i, j, lGradPres, ndm) = MAT2D(vizNel, j, gradPres, ndm);
              MAT2D(i, j, lDfield, ndm) = MAT2D(vizNel, j, dField, ndm);
            }
/*...................................................................*/

/*...*/  
            for (k = 0; k<ndf; k++)
              for (j = 0; j<ndm; j++)
                MAT3D(i, k, j, lGradU0, ndf, ndm) 
               = MAT3D(vizNel, k, j, gradU0, ndf,ndm);
/*...................................................................*/
          }
        }
/*...................................................................*/

/*...*/
        if(fWallModel)
          for (i = 0; i<NWALLPAR; i++)
            lWallPar[i] = MAT2D(nel, i, wallPar, NWALLPAR);
/*...................................................................*/

/*...*/
        for (i = 0; i<ndf; i++)
          lWk[i] = MAT2D(nel,i,wk,ns); 
/*...................................................................*/    

/*... chamando a biblioteca de celulas*/
        cellLibCombustion(loads        , ldVel
                        , scAdv        , scDiff
                        , tModel       , cModel
                        , vProp        
                        , lGeomType    
                        , lViz         , lId
                        , lKsi         , lmKsi
                        , lEta         , lfArea
                        , lNormal      , lVolume
                        , lXm          , lXmcc
                        , lDcca        , lCc
                        , lvSkew       , lmvSkew
                        , lA           , lB
                        , lRcell       , ddt
                        , lFaceR       , lFaceL
                        , lFaceVelR    , lFaceVelL
                        , lu0          , lGradU0
                        , lWk          , lVel      
                        , lPres        , lGradPres
                        , lDensity     , lDiff
                        , lViscosity
                        , lDfield      , lWallPar
                        , underU
                        , nen[nel]     , nFace[nel]
                        , ndm          , lib
                        , nel);      
/*...................................................................*/

/*... residuo da celula*/
        if(calRcell)
          for (j = 0; j<ndf; j++)
            rCell[j*numel+nel] = lRcell[j];  
/*...................................................................*/
  
/*...*/
          assblyBlock(ia       , ja
                   , a         , ad
                   , b         , lId
                   , lA        , lB
                   , nEq       , nEqNov
                   , nAd       , nAdR
                   , nFace[nel], ndf  
                   , ndf       , ndf    
                   , storage   , forces
                   , matrix    , unsym);          
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 14/01/2018                                   *
 * Data de modificaco : 18/07/2018                                   *
 *-------------------------------------------------------------------*
 * systFormOneEqK: calculo do sistema de equacoes para problemas     *
 * transporte de energia cinetica turbulenta (Ax=b)                  *                    *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loads   -> definicoes de cargas K                                 *
 * ldsVel  -> definicoes de cargas para velocidades                  *
 * model   -> modelo da equacao de energia                           *
 * advT    -> tecnica da discretizacao do termo advecao              *
 * diffT   -> tecnica da discretizacao do termo difusivo             *
 * vProp     -> propedades variaveis (true|false)                    *
 * el      -> conetividade dos celulas                               *
 * nelcon  -> vizinhos dos elementos                                 *
 * nen     -> numero de nos por celulas                              *
 * nFace   -> numero de faces por celulas                            *
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * gVolume -> volumes das celulas                                    *
 * gDcca   -> menor distancia do centroide a faces desta celula      *
 * gXmCc   -> vetores que unem o centroide aos pontos medios das     *
 * gCc     -> centroide do elemento                                  *
 * fModKsi -> o mudolo do vetores que unem centroide da celula       *
 *            central aos vizinhos destas                            *
 * fKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        *
 * fEta    -> vetores paralelos as faces das celulas                 *
 * fArea   -> modulo do vetor eta                                    *
 * fNormal -> vetores normais as faces das celulas                   *
 * fXm     -> pontos medios das faces das celulas                    *
 * fModvSkew -> distacia entre o ponto medio a intersecao que une os *
 *            centrois compartilhado nessa face                      *
 * fvSkew  -> vetor entre o ponto medio a intersecao que une os      *
 *            centrois compartilhado nessa face                      *
 * calType -> tipo de calculo das celulas                            *
 * geomType-> tipo geometrico das celulas                            *
 * prop    -> propriedades dos material                              *
 * mat     -> material por celula                                    *
 * au      -> matriz de coeficientes esparsa                         *
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- triangular superior               )       *
 * ia      -> ponteiro para as linhas da matriz esparsa              *
 * ja      -> ponteiro para as colunas da matriz esparsa             *
 * a       -> matriz de coeficientes esparsa                         *
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- fora da diagonal principal        )       *
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 * b       -> vetor de forcas                                        *
 * id      -> numera das equacoes                                    *
 * faceReK -> restricoes por elemento                                *
 * faceLdK -> carga por elemento                                     *
 * faceReVel-> restricoes por elemento                               *
 * faceLdVel-> carga por elemento                                    *
 * u0      -> solucao conhecida                                      *
 * gradU0  -> gradiente da solucao conhecido                         *
 * vel       -> campo de velocidade conhecido                        *
 * gradVel   -> gradiente rescontruido da velocidade                 *
 * pres      -> pressao do tempo atual                               *
 * gradPres  -> gradiente de pressao do tempo atual                  *
 * density -> massa especifica com variacao temporal                 *
 * dViscosity-> viscosidade dinamica com variacao temporal           *
 * eddyVisc  -> viscosidade turbulenta                               *
 * dField    -> matriz D do metodo simple                            * 
 * rCell   -> nao definido                                           *
 * wallPar -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)    * 
 * cDyn    -> coeficiente dinanmicos ( Ck, Ce)                       * 
 * ddt     -> discretizacao temporal                                 *
 * nEq     -> numero de equacoes                                     *
 * neqNov  -> numero de equacoes nao sobrepostas                     *
 * nAd     -> numero de termos nao nulos                             *
 * nAdR    -> numero de termos nao nulos na parte retangular         *
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
 * rCell     -> residuo por celula                                   *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void systFormOneEqK(Loads *ldsK        ,Loads *ldsVel
       , Turbulence *tModel
       , Advection adv                 , Diffusion diff 
       , INT    *RESTRICT el           , INT    *RESTRICT nelcon
       , short  *RESTRICT nen          , short  *RESTRICT nFace
       , INT *RESTRICT cellFace        , INT *RESTRICT fOwner
       , DOUBLE *RESTRICT gVolume      , DOUBLE *RESTRICT gDcca
       , DOUBLE *RESTRICT gXmCc        , DOUBLE *RESTRICT gCc
       , DOUBLE *RESTRICT fModKsi      , DOUBLE *RESTRICT fKsi
       , DOUBLE *RESTRICT fEta         , DOUBLE *RESTRICT fArea
       , DOUBLE *RESTRICT fNormal      , DOUBLE *RESTRICT fXm
       , DOUBLE *RESTRICT fModvSkew    , DOUBLE *RESTRICT fvSkew
       , short  *RESTRICT geomType     , DOUBLE *RESTRICT prop
       , short  *RESTRICT calType      , short  *RESTRICT mat
       , INT    *RESTRICT ia           , INT    *RESTRICT ja
       , DOUBLE *RESTRICT a            , DOUBLE *RESTRICT ad
       , DOUBLE *RESTRICT b            , INT    *RESTRICT id
       , short  *RESTRICT faceReK      , short  *RESTRICT faceLdK
       , short  *RESTRICT faceReVel    , short  *RESTRICT faceLdVel
       , DOUBLE *RESTRICT u0           , DOUBLE *RESTRICT gradU0
       , DOUBLE *RESTRICT vel          , DOUBLE *RESTRICT gradVel
       , DOUBLE *RESTRICT pres         , DOUBLE *RESTRICT gradPres    
       , DOUBLE *RESTRICT density      , DOUBLE *RESTRICT dViscosity   
       , DOUBLE *RESTRICT eddyViscosity, DOUBLE *RESTRICT dField
       , DOUBLE *RESTRICT rCell        , DOUBLE *RESTRICT wallPar
       , DOUBLE *RESTRICT cDyn         , Temporal *ddt  
       , INT nEq                       , INT nEqNov
       , INT nAd                       , INT nAdR
       , short maxNo                   , short maxViz
       , short ndm                     , INT numel
       , short ndf                     , short storage
       , bool forces                   , bool matrix
       , bool calRcell                 , bool unsym)
{
  short i, j, k, aux1, aux2, lMat;;
  short nThreads = ompVar.nThreadsCell;
  INT nel, vizNel;
  
/*... variavel local */
  short  lGeomType[MAX_NUM_FACE + 1],
         lFaceReK[MAX_NUM_FACE + 1], lFaceLdK[MAX_NUM_FACE + 1],
         lFaceReVel[MAX_NUM_FACE + 1], lFaceLdVel[MAX_NUM_FACE + 1];
  INT    lId[(MAX_NUM_FACE + 1)*MAX_NDF], lViz[MAX_NUM_FACE];
  INT    idFace, cellOwner, ch;
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM], lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM], lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM], lVolume[MAX_NUM_FACE + 1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM], lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE], lvSkew[MAX_NUM_FACE*MAX_NDM];  
  DOUBLE lDensity[(MAX_NUM_FACE + 1)],lViscosity[(MAX_NUM_FACE + 1)*2];
  DOUBLE lA[(MAX_NUM_FACE + 1)*MAX_NDF], lB[MAX_NDF];
  DOUBLE lProp[(MAX_NUM_FACE + 1)*MAXPROP];
  DOUBLE lu0[(MAX_NUM_FACE + 1)*MAX_NDF],lPres[(MAX_NUM_FACE + 1)];
  DOUBLE lGradVel[(MAX_NUM_FACE + 1)*MAX_NDM*MAX_NDF];
  DOUBLE lGradU0[(MAX_NUM_FACE + 1)*MAX_NDM],lDfield[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lGradPres[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lVel[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lCc[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lRcell[MAX_NDF],lWallPar[NWALLPAR],lCdyn[2];


/*... loop nas celulas*/
/*...*/
  if (ompVar.fCell)
  {
    aux2 = maxViz + 1;
#pragma omp parallel  for default(none) num_threads(nThreads)\
    private(nel,i,j,k,aux1,lId,lMat,lVolume,lGeomType,lu0\
           ,lA,lB,lDensity,lProp,lGradU0,lDfield\
           ,lFaceReK,lFaceLdK,lFaceReVel,lFaceLdVel\
          ,lVel,lCc,lmKsi,lfArea,lDcca,lmvSkew,lKsi,lEta,lViscosity\
          ,lNormal,lXm,lXmcc,lvSkew,vizNel,lViz,lRcell,lPres,lGradPres\
          ,lGradVel,lWallPar,lCdyn,idFace,cellOwner,ch)\
    shared(aux2,ndm,ndf,numel,maxViz,calRcell,rCell,nFace,mat\
         ,calType,gVolume,geomType,faceReK,faceLdK,faceReVel,faceLdVel\
         ,u0,gCc,ldsK,ldsVel,dField\
         ,density,prop,gradU0,vel,fModKsi,dViscosity,eddyViscosity\
         ,fArea,gDcca,fModvSkew,fKsi,fEta,fNormal,fXm,gXmCc,fvSkew\
         ,nelcon,id,adv,diff,pres,gradPres\
         ,ddt,nen,ia,ja,a,ad,b,nEq,nEqNov,nAd,gradVel,wallPar\
         ,nAdR,storage,forces,matrix,unsym,tModel,cDyn\
         ,fOwner,cellFace)
/*... loop nas celulas*/
    for (nel = 0; nel<numel; nel++) {
/*...*/
      if (calRcell)
        for (j = 0; j<ndf; j++)
          MAT2D(nel, j, rCell, ndf) = 0.e0;;
/*...*/
      aux1 = nFace[nel];
/*... elementos com equacoes*/
      if (MAT2D(nel, aux1, faceReK, aux2) != PCCELL)
      {

/*... zerando vetores*/
        for (j = 0; j<(MAX_NUM_FACE + 1)*MAX_NDF; j++) 
        {
          lId[j] = -1;
          lu0[j] = 0.e0;
        }
/*... loop na celula central*/
        lMat             = mat[nel] - 1;
//      lib              = calType[lMat];
        lVolume[aux1]    = gVolume[nel];
        lGeomType[aux1]  = geomType[nel];
        lFaceReK[aux1]   = MAT2D(nel, aux1, faceReK, aux2);
        lFaceLdK[aux1]   = MAT2D(nel, aux1, faceLdK, aux2);
        lFaceReVel[aux1] = MAT2D(nel, aux1, faceReVel, aux2);
        lFaceLdVel[aux1] = MAT2D(nel, aux1, faceLdVel, aux2);
        lDensity[aux1]   = MAT2D(nel, 2, density, DENSITY_LEVEL);
        lCdyn[0]         = MAT2D(nel, 0, cDyn   , 2);
        lCdyn[1]         = MAT2D(nel, 1, cDyn   , 2);
/*...*/
        lPres[aux1] = pres[nel];
/*...................................................................*/

/*... viscosidade dinamica e turbulentea*/
        MAT2D(aux1, 0, lViscosity, 2) = dViscosity[nel];
        MAT2D(aux1, 1, lViscosity, 2) = eddyViscosity[nel];
/*...................................................................*/

/*...*/
        for (j = 0; j<ndf; j++) 
        {
          MAT2D(aux1, j, lu0, ndf) = MAT2D(nel, j, u0, ndf);
          MAT2D(aux1, j, lId, ndf) = MAT2D(nel, j, id, ndf) - 1;
        }
/*...................................................................*/

/*...*/
        for (j = 0; j<MAXPROP; j++)
          MAT2D(aux1, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);
/*...................................................................*/

/*...*/
        for (j = 0; j<ndm; j++)
        {
          MAT2D(aux1, j, lGradU0, ndm) = MAT2D(nel, j, gradU0, ndm);
          MAT2D(aux1, j, lVel, ndm) = MAT2D(nel, j, vel, ndm);
          MAT2D(aux1, j, lCc, ndm) = MAT2D(nel, j, gCc, ndm);
          MAT2D(aux1, j, lGradPres, ndm) = MAT2D(nel, j, gradPres, ndm);
          MAT2D(aux1, j, lDfield, ndm) = MAT2D(nel, j, dField, ndm);
        }
/*...................................................................*/

/*...*/
        for(i=0;i<ndm;i++)
          for(j=0;j<ndm;j++)
            MAT3D(aux1,i,j,lGradVel,ndm,ndm) 
                                = MAT3D(nel,i,j,gradVel,ndm,ndm);
/*...................................................................*/

/*...*/
        for (i = 0; i<aux1; i++) 
        {
          lDcca[i]     = MAT2D(nel, i, gDcca, maxViz);
          lFaceReK[i]   = MAT2D(nel, i, faceReK, aux2);
          lFaceLdK[i]   = MAT2D(nel, i, faceLdK, aux2);
          lFaceReVel[i] = MAT2D(nel, i, faceReVel, aux2);
          lFaceLdVel[i] = MAT2D(nel, i, faceLdVel, aux2);
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm); 
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
        }
/*...................................................................*/

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++) 
        {
          vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2) 
          {
            lVolume[i]       = gVolume[vizNel];
            lGeomType[i]     = geomType[vizNel];
            lDensity[i]      = MAT2D(vizNel, 2, density, DENSITY_LEVEL);
/*...*/
            lPres[i] = pres[vizNel];
/*...................................................................*/

/*... viscusidade dinamica e turbulentea*/
            MAT2D(i, 0, lViscosity, 2) = dViscosity[vizNel];
            MAT2D(i, 1, lViscosity, 2) = eddyViscosity[vizNel];

            lMat = mat[vizNel] - 1;
            for (j = 0; j<ndf; j++) 
            {
              MAT2D(i, j, lu0, ndf) = MAT2D(vizNel, j, u0, ndf);
              MAT2D(i, j, lId, ndf) = MAT2D(vizNel, j, id, ndf) - 1;
            }
            for (j = 0; j<ndm; j++) 
            {
                MAT2D(i, j, lGradU0, ndm) = MAT2D(vizNel, j, gradU0, ndm);
                MAT2D(i, j, lVel, ndm) = MAT2D(vizNel, j, vel, ndm);
                MAT2D(i, j, lCc, ndm) = MAT2D(vizNel, j, gCc, ndm);
                MAT2D(i, j, lGradPres, ndm) = MAT2D(vizNel, j, gradPres, ndm);
                MAT2D(i ,j, lDfield, ndm) = MAT2D(vizNel, j, dField, ndm);
            }

            for (k = 0; k<ndm; k++)
              for (j = 0; j<ndm; j++)
                MAT3D(i, k, j, lGradVel, ndm, ndm)
                       = MAT3D(vizNel, k, j, gradVel, ndm, ndm);

            for (j = 0; j<DIFPROP; j++)
              MAT2D(i, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);
          }
        }
/*...................................................................*/

/*...*/    
        for(i=0;i<NWALLPAR;i++)
          lWallPar[i] = MAT2D(nel,i,wallPar,NWALLPAR);
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibOneEqK(ldsK       , ldsVel
                    , tModel
                    , &adv       , &diff 
                    , lGeomType  , lProp
                    , lViz       , lId
                    , lKsi       , lmKsi
                    , lEta       , lfArea
                    , lNormal    , lVolume
                    , lXm        , lXmcc
                    , lDcca      , lCc
                    , lvSkew     , lmvSkew
                    , lA         , lB
                    , lRcell     , ddt
                    , lFaceReK   , lFaceLdK
                    , lFaceReVel , lFaceLdVel
                    , lu0        , lGradU0
                    , lVel       , lGradVel
                    , lPres      , lGradPres 
                    , lDensity   , lViscosity     
                    , lDfield    , lWallPar
                    , lCdyn
                    , nen[nel]   , nFace[nel]
                    , ndm        , nel);
/*...................................................................*/

/*... residuo da celula*/
        if (calRcell)
          for (j = 0; j<ndf; j++)
           MAT2D(nel, j, rCell, ndf) = lRcell[j];
/*...................................................................*/

/*...*/
        assbly(ia        ,ja
                ,a         ,ad
                ,b         ,lId
                ,lA        ,lB
                ,nEq       ,nEqNov
                ,nAd       ,nAdR
                ,nFace[nel],ndf
                ,storage   ,forces
                ,matrix    ,unsym);
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  else
  {
    aux2 = maxViz + 1;
/*... loop nas celulas*/
    for (nel = 0; nel<numel; nel++) {
/*...*/
      if (calRcell)
        for (j = 0; j<ndf; j++)
          MAT2D(nel, j, rCell, ndf) = 0.e0;;
/*...*/
      aux1 = nFace[nel];
/*... elementos com equacoes*/
      if (MAT2D(nel, aux1, faceReK, aux2) != PCCELL)
      {

/*... zerando vetores*/
        for (j = 0; j<(MAX_NUM_FACE + 1)*MAX_NDF; j++)
        {
          lId[j] = -1;
          lu0[j] = 0.e0;
        }
/*... loop na celula central*/
        lMat = mat[nel] - 1;
        //      lib              = calType[lMat];
        lVolume[aux1] = gVolume[nel];
        lGeomType[aux1] = geomType[nel];
        lFaceReK[aux1] = MAT2D(nel, aux1, faceReK, aux2);
        lFaceLdK[aux1] = MAT2D(nel, aux1, faceLdK, aux2);
        lFaceReVel[aux1] = MAT2D(nel, aux1, faceReVel, aux2);
        lFaceLdVel[aux1] = MAT2D(nel, aux1, faceLdVel, aux2);
        lDensity[aux1] = MAT2D(nel, 2, density, DENSITY_LEVEL);
        lCdyn[0] = MAT2D(nel, 0, cDyn, 2);
        lCdyn[1] = MAT2D(nel, 1, cDyn, 2);
/*...*/
        lPres[aux1] = pres[nel];
/*...................................................................*/

/*... viscosidade dinamica e turbulentea*/
        MAT2D(aux1, 0, lViscosity, 2) = dViscosity[nel];
        MAT2D(aux1, 1, lViscosity, 2) = eddyViscosity[nel];
/*...................................................................*/

/*...*/
        for (j = 0; j<ndf; j++)
        {
          MAT2D(aux1, j, lu0, ndf) = MAT2D(nel, j, u0, ndf);
          MAT2D(aux1, j, lId, ndf) = MAT2D(nel, j, id, ndf) - 1;
        }
/*...................................................................*/

/*...*/
        for (j = 0; j<MAXPROP; j++)
          MAT2D(aux1, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);
/*...................................................................*/

/*...*/
        for (j = 0; j<ndm; j++)
        {
          MAT2D(aux1, j, lGradU0, ndm) = MAT2D(nel, j, gradU0, ndm);
          MAT2D(aux1, j, lVel, ndm) = MAT2D(nel, j, vel, ndm);
          MAT2D(aux1, j, lCc, ndm) = MAT2D(nel, j, gCc, ndm);
          MAT2D(aux1, j, lGradPres, ndm) = MAT2D(nel, j, gradPres, ndm);
          MAT2D(aux1, j, lDfield, ndm) = MAT2D(nel, j, dField, ndm);
        }
/*...................................................................*/

/*...*/
        for (i = 0; i<ndm; i++)
          for (j = 0; j<ndm; j++)
            MAT3D(aux1, i, j, lGradVel, ndm, ndm)
            = MAT3D(nel, i, j, gradVel, ndm, ndm);
/*...................................................................*/

/*...*/
        for (i = 0; i<aux1; i++)
        {
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lFaceReK[i] = MAT2D(nel, i, faceReK, aux2);
          lFaceLdK[i] = MAT2D(nel, i, faceLdK, aux2);
          lFaceReVel[i] = MAT2D(nel, i, faceReVel, aux2);
          lFaceLdVel[i] = MAT2D(nel, i, faceLdVel, aux2);
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
        }
/*...................................................................*/

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++)
        {
          vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2)
          {
            lVolume[i] = gVolume[vizNel];
            lGeomType[i] = geomType[vizNel];
            lDensity[i] = MAT2D(vizNel, 2, density, DENSITY_LEVEL);
/*...*/
            lPres[i] = pres[vizNel];
/*...................................................................*/

/*... viscusidade dinamica e turbulentea*/
            MAT2D(i, 0, lViscosity, 2) = dViscosity[vizNel];
            MAT2D(i, 1, lViscosity, 2) = eddyViscosity[vizNel];

            lMat = mat[vizNel] - 1;
            for (j = 0; j<ndf; j++)
            {
              MAT2D(i, j, lu0, ndf) = MAT2D(vizNel, j, u0, ndf);
              MAT2D(i, j, lId, ndf) = MAT2D(vizNel, j, id, ndf) - 1;
            }
            for (j = 0; j<ndm; j++)
            {
              MAT2D(i, j, lGradU0, ndm) = MAT2D(vizNel, j, gradU0, ndm);
              MAT2D(i, j, lVel, ndm) = MAT2D(vizNel, j, vel, ndm);
              MAT2D(i, j, lCc, ndm) = MAT2D(vizNel, j, gCc, ndm);
              MAT2D(i, j, lGradPres, ndm) = MAT2D(vizNel, j, gradPres, ndm);
              MAT2D(i, j, lDfield, ndm) = MAT2D(vizNel, j, dField, ndm);
            }

            for (k = 0; k<ndm; k++)
              for (j = 0; j<ndm; j++)
                MAT3D(i, k, j, lGradVel, ndm, ndm)
                = MAT3D(vizNel, k, j, gradVel, ndm, ndm);

            for (j = 0; j<DIFPROP; j++)
              MAT2D(i, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);
          }
        }
/*...................................................................*/

/*...*/
        for (i = 0; i<NWALLPAR; i++)
          lWallPar[i] = MAT2D(nel, i, wallPar, NWALLPAR);
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibOneEqK(ldsK, ldsVel
          , tModel
          , &adv, &diff
          , lGeomType, lProp
          , lViz, lId
          , lKsi, lmKsi
          , lEta, lfArea
          , lNormal, lVolume
          , lXm, lXmcc
          , lDcca, lCc
          , lvSkew, lmvSkew
          , lA, lB
          , lRcell, ddt
          , lFaceReK, lFaceLdK
          , lFaceReVel, lFaceLdVel
          , lu0, lGradU0
          , lVel, lGradVel
          , lPres, lGradPres
          , lDensity, lViscosity
          , lDfield, lWallPar
          , lCdyn
          , nen[nel], nFace[nel]
          , ndm, nel);
/*...................................................................*/

/*... residuo da celula*/
        if (calRcell)
          for (j = 0; j<ndf; j++)
            MAT2D(nel, j, rCell, ndf) = lRcell[j];
/*...................................................................*/

/*...*/
        assbly(ia, ja
          , a, ad
          , b, lId
          , lA, lB
          , nEq, nEqNov
          , nAd, nAdR
          , nFace[nel], ndf
          , storage, forces
          , matrix, unsym);
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 10/09/2016                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * VELRESIDUAL: Calculo do residuo das velocidades                   *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  *
 * loadsPres -> definicoes de cargas de pressao                      *
 * advVel    -> tecnica da discretizacao do termo advecao            *
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * el        -> conetividade dos celulas                             *
 * nelcon    -> vizinhos dos elementos                               *
 * nen       -> numero de nos por celulas                            *
 * nFace     -> numero de faces por celulas                          *
 * calType   -> tipo de calculo das celulas                          *
 * geomType  -> tipo geometrico das celulas                          *
 * prop      -> propriedades dos material                            *
 * mat       -> material por celula                                  *
 * cc        -> centroide das celulas                                *
 * ad        -> coeficientes da diagonal principal                   *
 * gKsi      -> vetores que unem centroide da celula central aos     *
 *              vizinhos destas                                      *
 * gmKsi     -> modulo do vetor ksi                                  *
 * gEta      -> vetores paralelos as faces das celulas               *
 * gfArea    -> modulo do vetor eta                                  *
 * gNormal   -> vetores normais as faces das celulas                 *
 * gVolume   -> volumes das celulas                                  *
 * gXm       -> pontos medios das faces das celulas                  *
 * gXmcc     -> vetores que unem o centroide aos pontos medios das   *
 *              faces                                                *
 * gvSkew    -> vetor entre o ponto medio a intersecao que une os    *
 *              centrois compartilhado nessa face                    *
 * gmvSkew   -> distacia entre o ponto medio a intersecao que une os *
 *              centrois compartilhado nessa face                    *
 * gDcca     -> menor distancia do centroide a faces desta celula    *
 * density   -> massa especifica com variacao temporal               *
 * id        -> numera das equacoes                                  *
 * faceVelR  -> restricoes por elemento de velocidades               *
 * faceVelL  -> carga por elemento de velocidades                    *
 * facePresR -> restricoes por elemento de pressao                   *
 * facePresL -> carga por elemento de pressao                        *
 * pres      -> campo de pressao conhecido                           *
 * gradVel   -> gradiente da solucao conhecido                       *
 * dField    -> matriz D do metodo simple                            *
 * vel       -> campo de velocidade conhecido                        *
 * res       -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * maxNo     -> numero de nos por celula maximo da malha             *
 * maxViz    -> numero vizinhos por celula maximo da malha           *
 * ndm       -> numero de dimensoes                                  *
 * numel     -> numero de toral de celulas                           *
 * ndf       -> graus de liberdade                                   *
 * sPressure -> reconstrucao de segunda ordem para pressoes nas      *
 *              faces                                                *
 * fResidual -> calculo do residuo                                   *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * res       -> residuo                                              *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
*********************************************************************/
void velResidual(Loads *loadsVel            , Loads *loadsPres
              , Advection advVel           , Diffusion diffVel
              , INT    *RESTRICT el        , INT    *RESTRICT nelcon
              , short  *RESTRICT nen       , short  *RESTRICT nFace
              , short  *RESTRICT geomType  , DOUBLE *RESTRICT prop
              , short  *RESTRICT calType   , short  *RESTRICT mat
              , DOUBLE *RESTRICT cc        , DOUBLE *RESTRICT ad
              , DOUBLE *RESTRICT gKsi      , DOUBLE *RESTRICT gmKsi
              , DOUBLE *RESTRICT gEta      , DOUBLE *RESTRICT gfArea
              , DOUBLE *RESTRICT gNormal   , DOUBLE *RESTRICT gVolume
              , DOUBLE *RESTRICT gXm       , DOUBLE *RESTRICT gXmcc
              , DOUBLE *RESTRICT gvSkew    , DOUBLE *RESTRICT gmvSkew
              , DOUBLE *RESTRICT gDcca     , DOUBLE *RESTRICT density
              , short  *RESTRICT faceVelR  , short  *RESTRICT faceVelL
              , short  *RESTRICT facePresR , short  *RESTRICT facePresL
              , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres
              , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT res  
              , DOUBLE *RESTRICT gradVel   , DOUBLE *RESTRICT bT
              , DOUBLE *RESTRICT dField    , DOUBLE underU
              , Temporal ddt                 
              , short maxNo                , short maxViz
              , short ndm                  , INT numel
              , short ndf                  , bool sPressure
              , bool fResidual)
{
  short i, j, k;
  short nThreads = ompVar.nThreadsCell;
  INT nel, vizNel;
  /*... variavel local*/
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM], lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM], lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM], lVolume[MAX_NUM_FACE + 1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM], lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE], lvSkew[MAX_NUM_FACE*MAX_NDM];
  short  lGeomType[MAX_NUM_FACE + 1];
  short  lFaceVelR[MAX_NUM_FACE + 1], lFacePresR[MAX_NUM_FACE + 1];
  short  lFaceVelL[MAX_NUM_FACE + 1], lFacePresL[MAX_NUM_FACE + 1];
  DOUBLE lDensity[(MAX_NUM_FACE + 1)];
  DOUBLE lProp[(MAX_NUM_FACE + 1)*MAXPROP];
  DOUBLE lPres[(MAX_NUM_FACE + 1)], lB[2*MAX_NDM],lBt[2*MAX_NDM];
  DOUBLE lGradVel[(MAX_NUM_FACE + 1)*MAX_NDM*MAX_NDF];
  DOUBLE lGradPres[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lVel[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lCc[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lDfield[MAX_NDM];
  INT    lViz[MAX_NUM_FACE];
  short  aux1, aux2, lMat;

/*...*/
  if (ompVar.fCell) {
/*... loop nas celulas*/
    aux2 = maxViz + 1;
#pragma omp parallel  for default(none) num_threads(nThreads)\
     private(nel,i,j,k,aux1,lPres,lMat,lVolume,lGeomType\
          ,lFaceVelR,lDfield,lB,lBt\
          ,lFaceVelL,lFacePresR,lFacePresL,lDensity,lProp,lGradPres\
          ,lVel,lCc,lGradVel,lmKsi,lfArea,lDcca,lmvSkew,lKsi,lEta\
          ,lNormal,lXm,lXmcc,lvSkew,vizNel,lViz)\
     shared(aux2,ndm,ndf,numel,maxViz,nFace,mat,ad,res\
         ,calType,gVolume, geomType, faceVelR, faceVelL, facePresR\
         ,facePresL,density, prop,gradPres,vel,cc,gradVel,gmKsi\
         ,gfArea,gDcca,gmvSkew,gKsi,gEta,gNormal,gXm,gXmcc,gvSkew\
         ,nelcon,loadsVel,loadsPres,advVel,diffVel,bT\
         ,ddt,sPressure,nen,pres,dField,underU,fResidual)  
    for (nel = 0; nel<numel; nel++) {
/*...*/
      aux1 = nFace[nel];
/*... elementos com equacoes*/
      if (MAT2D(nel, aux1, faceVelR, aux2) != PCCELL) {

        for (j = 0; j<MAX_NUM_FACE + 1; j++)
          lPres[j] = 0.e0;

/*... loop na celula central*/
        lMat = mat[nel] - 1;
//      lib = calType[lMat];
        lVolume[aux1] = gVolume[nel];
        lGeomType[aux1] = geomType[nel];
        lFaceVelR[aux1] = MAT2D(nel, aux1, faceVelR, aux2);
        lFaceVelL[aux1] = MAT2D(nel, aux1, faceVelL, aux2);
        lFacePresR[aux1] = MAT2D(nel, aux1, facePresR, aux2);
        lFacePresL[aux1] = MAT2D(nel, aux1, facePresL, aux2);
        lDensity[aux1] = MAT2D(nel, 0, density, DENSITY_LEVEL);

        lPres[aux1] = pres[nel];

        for (j = 0; j<MAXPROP; j++)
          MAT2D(aux1, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);

        for (j = 0; j<ndm; j++) {
          MAT2D(aux1, j, lGradPres, ndm) = MAT2D(nel,j,gradPres,ndm);
          MAT2D(aux1, j, lVel, ndm)      = MAT2D(nel,j,vel     ,ndm);
          MAT2D(aux1, j, lCc, ndm)       = MAT2D(nel,j,cc      ,ndm);
          lBt[j]                         = MAT2D(nel,j,bT      ,ndm);
        }

        for (i = 0; i<ndf; i++)
          for (j = 0; j<ndm; j++)
            MAT3D(aux1, i, j, lGradVel, ndf, ndm)
            = MAT3D(nel, i, j, gradVel, ndf, ndm);

        for (i = 0; i<aux1; i++) {
          lmKsi[i] = MAT2D(nel, i, gmKsi, maxViz);
          lfArea[i] = MAT2D(nel, i, gfArea, maxViz);
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lmvSkew[i] = MAT2D(nel, i, gmvSkew, maxViz);
          lFaceVelR[i] = MAT2D(nel, i, faceVelR, aux2);
          lFaceVelL[i] = MAT2D(nel, i, faceVelL, aux2);
          lFacePresR[i] = MAT2D(nel, i, facePresR, aux2);
          lFacePresL[i] = MAT2D(nel, i, facePresL, aux2);
          for (j = 0; j<ndm; j++) {
            MAT2D(i, j, lKsi, ndm) = MAT3D(nel, i, j, gKsi, maxViz, ndm);
            MAT2D(i, j, lEta, ndm) = MAT3D(nel, i, j, gEta, maxViz, ndm);
            MAT2D(i, j, lNormal, ndm) = MAT3D(nel, i, j, gNormal, maxViz, ndm);
            MAT2D(i, j, lXm, ndm) = MAT3D(nel, i, j, gXm, maxViz, ndm);
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmcc, maxViz, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT3D(nel, i, j, gvSkew, maxViz, ndm);
          }
        }

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++) {
          vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2) {
            lVolume[i] = gVolume[vizNel];
            lGeomType[i] = geomType[vizNel];
            lDensity[i] = MAT2D(vizNel, 0, density, DENSITY_LEVEL);
            lMat = mat[vizNel] - 1;
            lPres[i] = pres[vizNel];

            for (j = 0; j<ndm; j++) {
              MAT2D(i, j, lGradPres, ndm) = MAT2D(vizNel, j, gradPres, ndm);
              MAT2D(i, j, lVel, ndm) = MAT2D(vizNel, j, vel, ndm);
              MAT2D(i, j, lCc, ndm) = MAT2D(vizNel, j, cc, ndm);
            }

            for (k = 0; k<ndf; k++)
              for (j = 0; j<ndm; j++)
                MAT3D(i, k, j, lGradVel, ndf, ndm)
                = MAT3D(vizNel, k, j, gradVel, ndf, ndm);

            for (j = 0; j<DIFPROP; j++)
              MAT2D(i, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);
          }
        }
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibVelExp(loadsVel  ,loadsPres
                     ,advVel    ,diffVel
                     ,lGeomType ,lProp
                     ,lViz
                     ,lKsi      ,lmKsi
                     ,lEta      ,lfArea
                     ,lNormal   ,lVolume
                     ,lXm       ,lXmcc
                     ,lDcca     ,lDensity
                     ,lvSkew    ,lmvSkew
                     ,lB        ,ddt
                     ,lFaceVelR ,lFaceVelL
                     ,lFacePresR,lFacePresL
                     ,lPres     ,lGradPres
                     ,lVel      ,lGradVel
                     ,lDfield   ,lCc
                     ,lBt       ,underU
                     ,sPressure ,fResidual
                     ,nen[nel]  ,nFace[nel]
                     ,ndm       ,nel);
/*...................................................................*/

/*...*/
        if (ndf == 2) {
          res[nel]       = lB[0];
          res[numel+nel] = lB[1];
          ad[nel]        = lB[2];
          ad[nel+1]      = lB[3];
        }
        else if (ndf == 3) {
          res[nel]         = lB[0];
          res[numel+nel]   = lB[1];
          res[2*numel+nel] = lB[2];
          ad[nel]          = lB[3];
          ad[nel+1]        = lB[4];
          ad[nel+2]        = lB[5];
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... sequencial*/
  else {
/*... loop nas celulas*/
    aux2 = maxViz + 1;
    for (nel = 0; nel<numel; nel++) {

/*...*/
      aux1 = nFace[nel];
/*... elementos com equacoes*/
      if (MAT2D(nel, aux1, faceVelR, aux2) != PCCELL) {

/*... zerando vetores*/
        for (j = 0; j<MAX_NUM_FACE + 1; j++)
          lPres[j] = 0.e0;

/*... loop na celula central*/
        lMat             = mat[nel] - 1;
//      lib              = calType[lMat];
        lVolume[aux1]    = gVolume[nel];
        lGeomType[aux1]  = geomType[nel];
        lFaceVelR[aux1]  = MAT2D(nel, aux1, faceVelR, aux2);
        lFaceVelL[aux1]  = MAT2D(nel, aux1, faceVelL, aux2);
        lFacePresR[aux1] = MAT2D(nel, aux1, facePresR, aux2);
        lFacePresL[aux1] = MAT2D(nel, aux1, facePresL, aux2);
        lDensity[aux1]   = MAT2D(nel, 0, density, DENSITY_LEVEL);

        lPres[aux1] = pres[nel];

        for (j = 0; j<MAXPROP; j++)
          MAT2D(aux1, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);

        for (j = 0; j<ndm; j++) {
          MAT2D(aux1,j,lGradPres,ndm) = MAT2D(nel,j,gradPres,ndm);
          MAT2D(aux1,j,lVel     ,ndm) = MAT2D(nel,j,vel     ,ndm);
          MAT2D(aux1,j,lCc      ,ndm) = MAT2D(nel,j,cc      ,ndm);
          lBt[j]                      = MAT2D(nel,j,bT      ,ndm);
        }

        for (i = 0; i<ndf; i++)
          for (j = 0; j<ndm; j++)
            MAT3D(aux1, i, j, lGradVel, ndf, ndm)
            = MAT3D(nel, i, j, gradVel, ndf, ndm);

        for (i = 0; i<aux1; i++) {
          lmKsi[i] = MAT2D(nel, i, gmKsi, maxViz);
          lfArea[i] = MAT2D(nel, i, gfArea, maxViz);
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lmvSkew[i] = MAT2D(nel, i, gmvSkew, maxViz);
          lFaceVelR[i] = MAT2D(nel, i, faceVelR, aux2);
          lFaceVelL[i] = MAT2D(nel, i, faceVelL, aux2);
          lFacePresR[i] = MAT2D(nel, i, facePresR, aux2);
          lFacePresL[i] = MAT2D(nel, i, facePresL, aux2);
          for (j = 0; j<ndm; j++) {
            MAT2D(i, j, lKsi, ndm) = MAT3D(nel, i, j, gKsi, maxViz, ndm);
            MAT2D(i, j, lEta, ndm) = MAT3D(nel, i, j, gEta, maxViz, ndm);
            MAT2D(i, j, lNormal, ndm) = MAT3D(nel, i, j, gNormal, maxViz, ndm);
            MAT2D(i, j, lXm, ndm) = MAT3D(nel, i, j, gXm, maxViz, ndm);
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmcc, maxViz, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT3D(nel, i, j, gvSkew, maxViz, ndm);
          }
        }

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++) {
          vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2) {
            lVolume[i] = gVolume[vizNel];
            lGeomType[i] = geomType[vizNel];
            lDensity[i] = MAT2D(vizNel, 0, density, DENSITY_LEVEL);
            lMat = mat[vizNel] - 1;
            lPres[i] = pres[vizNel];

            for (j = 0; j<ndm; j++) {
              MAT2D(i, j, lGradPres, ndm) = MAT2D(vizNel, j, gradPres, ndm);
              MAT2D(i, j, lVel, ndm) = MAT2D(vizNel, j, vel, ndm);
              MAT2D(i, j, lCc, ndm) = MAT2D(vizNel, j, cc, ndm);
            }

            for (k = 0; k<ndf; k++)
              for (j = 0; j<ndm; j++)
                MAT3D(i, k, j, lGradVel, ndf, ndm)
                = MAT3D(vizNel, k, j, gradVel, ndf, ndm);

            for (j = 0; j<DIFPROP; j++)
              MAT2D(i, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);
          }
        }
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibVelExp(loadsVel  ,loadsPres
                     ,advVel    ,diffVel
                     ,lGeomType ,lProp
                     ,lViz
                     ,lKsi      ,lmKsi
                     ,lEta      ,lfArea
                     ,lNormal   ,lVolume
                     ,lXm       ,lXmcc
                     ,lDcca     ,lDensity
                     ,lvSkew    ,lmvSkew
                     ,lB        ,ddt
                     ,lFaceVelR ,lFaceVelL
                     ,lFacePresR,lFacePresL
                     ,lPres     ,lGradPres
                     ,lVel      ,lGradVel
                     ,lDfield   ,lCc
                     ,lBt       ,underU
                     ,sPressure ,fResidual
                     ,nen[nel]  ,nFace[nel]
                     ,ndm       ,nel);
/*...................................................................*/

/*...*/
        if(ndf == 2){
          res[nel      ] = lB[0];
          res[numel+nel] = lB[1];
          ad[nel]        = lB[2];
          ad[nel+1]      = lB[3];
        }
        else if (ndf == 3) {
          res[nel]             = lB[0];
          res[numel + nel]     = lB[1];
          res[2 * numel + nel] = lB[2];
          ad[nel]              = lB[3];
          ad[nel + 1]          = lB[4];
          ad[nel + 2]          = lB[5];
        }        
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 17/09/2017                                   *
 * Data de modificaco : 27/08/2019                                   * 
 *-------------------------------------------------------------------* 
 * SYSTFOMSIMPLEPRESLM:calculo do sistema de equacoes para problemas * 
 * de escomaneto de fluidos para baixo Mach(Pres)                    * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loadsVel  -> definicoes de cargas de velocidades                  * 
 * loadsPres -> definicoes de cargas de pressao                      * 
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * eMass     -> termos/modelos da equacao de mass                    *
 * tModel    -> modelo de turbuelencia                               *
 * el      -> conetividade dos celulas                               * 
 * nelcon  -> vizinhos dos elementos                                 * 
 * nen     -> numero de nos por celulas                              * 
 * nFace   -> numero de faces por celulas                            * 
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * gVolume -> volumes das celulas                                    *
 * gDcca   -> menor distancia do centroide a faces desta celula      *
 * gXmCc   -> vetores que unem o centroide aos pontos medios das     *
 * gCc     -> centroide do elemento                                  *
 * fModKsi -> o mudolo do vetores que unem centroide da celula       *
 *            central aos vizinhos destas                            *
 * fKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        *
 * fEta    -> vetores paralelos as faces das celulas                 *
 * fArea   -> modulo do vetor eta                                    *
 * fNormal -> vetores normais as faces das celulas                   *
 * fXm     -> pontos medios das faces das celulas                    *
 * fModvSkew -> distacia entre o ponto medio a intersecao que une os *
 *            centrois compartilhado nessa face                      *
 * fvSkew  -> vetor entre o ponto medio a intersecao que une os      *
 *            centrois compartilhado nessa face                      *
 * calType -> tipo de calculo das celulas                            *
 * geomType-> tipo geometrico das celulas                            *
 * mat     -> material por celula                                    *
 * ia      -> ponteiro para as linhas da matriz esparsa              * 
 * ja      -> ponteiro para as colunas da matriz esparsa             * 
 * a       -> matriz de coeficientes esparsa                         * 
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- fora da diagonal principal        )       *
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 * b       -> vetor de forcas                                        * 
 * id      -> numera das equacoes                                    * 
 * faceVelR  -> restricoes por elemento de velocidades               * 
 * faceVelL  -> carga por elemento de velocidades                    * 
 * facePresR -> restricoes por elemento de pressao                   * 
 * facePresL -> carga por elemento de pressao                        * 
 * pres    -> campo de pressao conhecido                             * 
 * gradVel -> gradiente da solucao conhecido                         * 
 * vel     -> campo de velocidade conhecido                          * 
 * dField  -> matriz D do metodo simple ( volume/A(i,i) )            *
 * temp    -> temperatura                                            *   
 * wallPar -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)    *
 * rCell   -> nao definido                                           *  
 * density -> massa especifica com variacao temporal                 *
 * ddt     -> discretizacao temporal                                 *
 * nEq     -> numero de equacoes                                     *
 * neqNov  -> numero de equacoes nao sobrepostas                     *
 * nAd     -> numero de termos nao nulos                             *
 * nAdR    -> numero de termos nao nulos na parte retangular         *
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
 * rCell     -> residuo por celula                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void systFormSimplePresLm(Loads *loadsVel  , Loads *loadsPres 
							 ,Diffusion *diffPres        , MassEqModel *eMass 
               ,Turbulence *tModel   
               ,INT    *RESTRICT el        , INT    *RESTRICT nelcon 
               ,short  *RESTRICT nen       , short  *RESTRICT nFace
               , INT *RESTRICT cellFace    , INT *RESTRICT fOwner
               , DOUBLE *RESTRICT gVolume  , DOUBLE *RESTRICT gDcca
               , DOUBLE *RESTRICT gXmCc    , DOUBLE *RESTRICT gCc
               , DOUBLE *RESTRICT fModKsi  , DOUBLE *RESTRICT fKsi
               , DOUBLE *RESTRICT fEta     , DOUBLE *RESTRICT fArea
               , DOUBLE *RESTRICT fNormal  , DOUBLE *RESTRICT fXm
               , DOUBLE *RESTRICT fModvSkew, DOUBLE *RESTRICT fvSkew
               , short  *RESTRICT geomType 
               , short  *RESTRICT calType  , short  *RESTRICT mat
               , INT    *RESTRICT ia       , INT    *RESTRICT ja
               , DOUBLE *RESTRICT a        , DOUBLE *RESTRICT ad 
               , DOUBLE *RESTRICT b        , INT    *RESTRICT id
               , short  *RESTRICT faceVelR , short  *RESTRICT faceVelL       
               , short  *RESTRICT facePresR, short  *RESTRICT facePresL      
               , DOUBLE *RESTRICT pres     , DOUBLE *RESTRICT gradPres
               , DOUBLE *RESTRICT vel      , DOUBLE *RESTRICT dField
               , DOUBLE *RESTRICT temp     , DOUBLE *RESTRICT wallPar  
               , DOUBLE *RESTRICT rCell    , DOUBLE *RESTRICT density
               , Temporal *ddt 
               , INT nEq                   , INT nEqNov
               , INT nAd                   , INT nAdR                  
               , short maxNo               , short maxViz
               , short ndm                 , INT numel
               , short ndf                 , short storage
               , bool forces               , bool matrix 
               , bool calRcell             , bool unsym) 
{
  bool fWallModel = tModel->fWall;
  short i,j;
  short nThreads = ompVar.nThreadsCell;
  INT nel,vizNel;
/*... variavel local */
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM],lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM],lVolume[MAX_NUM_FACE+1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM],lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE],lvSkew[MAX_NUM_FACE*MAX_NDM];
  short  lGeomType[MAX_NUM_FACE+1];
  short  lib;
  short  lFaceVelR[MAX_NUM_FACE+1],lFacePresR[MAX_NUM_FACE+1];
  short  lFaceVelL[MAX_NUM_FACE+1],lFacePresL[MAX_NUM_FACE+1];
  DOUBLE lDensity[(MAX_NUM_FACE+1)*DENSITY_LEVEL];
  DOUBLE lA[(MAX_NUM_FACE+1)*MAX_NDF],lB[MAX_NDF];
  DOUBLE lPres[(MAX_NUM_FACE+1)],lTemp[(MAX_NUM_FACE + 1)];
  DOUBLE lGradPres[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lVel[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lDfield[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lRcell,lWallPar[NWALLPAR];
  INT    lId[(MAX_NUM_FACE+1)*MAX_NDF],lViz[MAX_NUM_FACE];
  INT idFace, cellOwner, ch;
  short  aux1,aux2,lMat;
/*...*/
  if(ompVar.fCell)
  {
/*... loop nas celulas*/
    aux2 = maxViz + 1;
#pragma omp parallel  for default(none) num_threads(nThreads)\
     private(nel,i,j,aux1,lId,lPres,lMat,lib,lVolume,lGeomType\
          ,lFaceVelR,lA,lB,lDfield,lTemp\
          ,lFaceVelL,lFacePresR,lFacePresL,lDensity,lGradPres\
          ,lVel,lmKsi,lfArea,lDcca,lmvSkew,lKsi,lEta\
          ,lNormal,lXm,lXmcc,lvSkew,vizNel,lViz,lRcell,lWallPar\
          ,idFace,cellOwner,ch)\
     shared(aux2,ndm,ndf,numel,maxViz,calRcell,rCell,nFace,mat\
         ,calType,gVolume, geomType, faceVelR, faceVelL, facePresR,temp\
         ,facePresL,density,gradPres,vel,fModKsi\
         ,fArea,gDcca,fModvSkew,fKsi,fEta,fNormal,fXm,gXmCc,fvSkew\
         ,nelcon,id,loadsVel,loadsPres,diffPres,eMass\
         ,ddt,nen,ia,ja,a,ad,b,nEq,nEqNov,nAd\
         ,nAdR,storage,forces,matrix,unsym,pres,dField,wallPar\
         ,fOwner,cellFace,fWallModel)
    for (nel = 0; nel<numel; nel++) {
/*...*/
      if (calRcell)
        rCell[nel] = 0.e0;;
/*...*/
      aux1 = nFace[nel];
/*... elementos com equacoes*/
      if (MAT2D(nel, aux1, facePresR, aux2) != PCCELL) {

/*... zerando vetores*/
        for (j = 0; j<(MAX_NUM_FACE + 1)*MAX_NDF; j++)
          lId[j] = -1;

        for (j = 0; j<MAX_NUM_FACE + 1; j++)
          lPres[j] = 0.e0;

/*... loop na celula central*/
        lMat = mat[nel] - 1;
        lib = calType[lMat];
        lVolume[aux1] = gVolume[nel];
        lGeomType[aux1] = geomType[nel];
        lId[aux1] = id[nel] - 1;
        lFaceVelR[aux1] = MAT2D(nel, aux1, faceVelR, aux2);
        lFaceVelL[aux1] = MAT2D(nel, aux1, faceVelL, aux2);
        lFacePresR[aux1] = MAT2D(nel, aux1, facePresR, aux2);
        lFacePresL[aux1] = MAT2D(nel, aux1, facePresL, aux2);
/*...*/
        MAT2D(aux1, 0, lDensity, DENSITY_LEVEL) = MAT2D(nel
          , TIME_N_MINUS_2, density, DENSITY_LEVEL);
        MAT2D(aux1, 1, lDensity, DENSITY_LEVEL) = MAT2D(nel
          , TIME_N_MINUS_1, density, DENSITY_LEVEL);
        MAT2D(aux1, 2, lDensity, DENSITY_LEVEL) = MAT2D(nel
          , TIME_N, density, DENSITY_LEVEL);
/*...................................................................*/

/*...*/
        lPres[aux1] = pres[nel];
        lTemp[aux1] = temp[nel];
/*...................................................................*/

/*...*/
        lDfield[aux1] = dField[nel];
        for (j = 0; j<ndm; j++)
        {
          MAT2D(aux1, j, lGradPres, ndm) = MAT2D(nel, j, gradPres, ndm);
          MAT2D(aux1, j, lVel, ndm) = MAT2D(nel, j, vel, ndm);
          MAT2D(aux1, j, lDfield, ndm) = MAT2D(nel, j, dField, ndm);
        }
/*...................................................................*/

/*...*/
        for (i = 0; i<aux1; i++)
        {
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lFaceVelR[i] = MAT2D(nel, i, faceVelR, aux2);
          lFaceVelL[i] = MAT2D(nel, i, faceVelL, aux2);
          lFacePresR[i] = MAT2D(nel, i, facePresR, aux2);
          lFacePresL[i] = MAT2D(nel, i, facePresL, aux2);
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          } 
        }
/*...................................................................*/

/*... loop na celulas vizinhas*/
        for (i = 0; i<aux1; i++)
        {
          vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
          lViz[i] = vizNel;
          if (vizNel != -2) {
            lVolume[i] = gVolume[vizNel];
            lGeomType[i] = geomType[vizNel];
            lMat = mat[vizNel] - 1;
            lId[i] = id[vizNel] - 1;
/*...*/
            MAT2D(i, 0, lDensity, DENSITY_LEVEL) = MAT2D(vizNel
              , TIME_N_MINUS_2, density, DENSITY_LEVEL);
            MAT2D(i, 1, lDensity, DENSITY_LEVEL) = MAT2D(vizNel
              , TIME_N_MINUS_1, density, DENSITY_LEVEL);
            MAT2D(i, 2, lDensity, DENSITY_LEVEL) = MAT2D(vizNel
              , TIME_N, density, DENSITY_LEVEL);
/*...................................................................*/

/*...*/
            lPres[i] = pres[vizNel];
            lTemp[i] = temp[vizNel];
/*...................................................................*/

            for (j = 0; j<ndm; j++)
            {
              MAT2D(i, j, lGradPres, ndm) = MAT2D(vizNel, j, gradPres, ndm);
              MAT2D(i, j, lVel, ndm) = MAT2D(vizNel, j, vel, ndm);
              MAT2D(i, j, lDfield, ndm) = MAT2D(vizNel, j, dField, ndm);
            }
          }
        }
/*...................................................................*/

/*...*/
        if(fWallModel)
          for (i = 0; i<NWALLPAR; i++)
            lWallPar[i] = MAT2D(nel, i, wallPar, NWALLPAR);
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibSimplePresLm(loadsVel  , loadsPres
                          , diffPres  , eMass
                          , lGeomType 
                          , lViz      , lId
                          , lKsi      , lmKsi
                          , lEta      , lfArea
                          , lNormal   , lVolume
                          , lXm       , lXmcc
                          , lDcca     , lDensity
                          , lvSkew    , lmvSkew
                          , lA, lB
                          , &lRcell   , ddt
                          , lFaceVelR , lFaceVelL
                          , lFacePresR, lFacePresL
                          , lPres     , lGradPres
                          , lVel      , lDfield
                          , lTemp     , lWallPar
                          , nen[nel]  , nFace[nel]
                          , ndm       , lib
                          , nel);
/*...................................................................*/

/*... residuo da celula*/
        if (calRcell)
          rCell[nel] = lRcell;
/*...................................................................*/

/*...*/
        assbly(ia, ja
          , a, ad
          , b, lId
          , lA, lB
          , nEq, nEqNov
          , nAd, nAdR
          , nFace[nel], 1
          , storage, forces
          , matrix, unsym);
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  else
  {
/*... loop nas celulas*/
    aux2    = maxViz+1;
    for(nel=0;nel<numel;nel++){
/*...*/
      if(calRcell)
        rCell[nel] = 0.e0;;
/*...*/
      aux1    = nFace[nel];
/*... elementos com equacoes*/
      if(MAT2D(nel,aux1,facePresR ,aux2) != PCCELL){

/*... zerando vetores*/
        for(j=0;j<(MAX_NUM_FACE+1)*MAX_NDF;j++) 
          lId[j]  = -1;
       
        for(j=0;j<MAX_NUM_FACE+1;j++) 
          lPres[j] = 0.e0;    

/*... loop na celula central*/    
        lMat            = mat[nel]-1;
        lib             = calType[lMat];
        lVolume[aux1]   = gVolume[nel]; 
        lGeomType[aux1] = geomType[nel];
        lId[aux1]       = id[nel] - 1;
        lFaceVelR[aux1] = MAT2D(nel,aux1,faceVelR ,aux2);
        lFaceVelL[aux1] = MAT2D(nel,aux1,faceVelL ,aux2);
        lFacePresR[aux1]= MAT2D(nel,aux1,facePresR ,aux2);
        lFacePresL[aux1]= MAT2D(nel,aux1,facePresL ,aux2);
/*...*/
        MAT2D(aux1,0,lDensity, DENSITY_LEVEL) = MAT2D(nel 
                             , TIME_N_MINUS_2,density, DENSITY_LEVEL);
        MAT2D(aux1,1,lDensity, DENSITY_LEVEL) = MAT2D(nel
                             , TIME_N_MINUS_1,density, DENSITY_LEVEL);
        MAT2D(aux1,2,lDensity, DENSITY_LEVEL) = MAT2D(nel
                             , TIME_N,density, DENSITY_LEVEL);
/*...................................................................*/ 

/*...*/
        lPres[aux1]     = pres[nel];
        lTemp[aux1]     = temp[nel];
/*...................................................................*/ 
    
/*...*/
        lDfield[aux1] = dField[nel];
        for(j=0;j<ndm;j++) 
        {
          MAT2D(aux1,j,lGradPres,ndm) = MAT2D(nel,j,gradPres,ndm); 
          MAT2D(aux1,j,lVel     ,ndm) = MAT2D(nel,j,vel     ,ndm);
          MAT2D(aux1,j,lDfield  ,ndm) = MAT2D(nel,j,dField  ,ndm);
        }
/*...................................................................*/

/*...*/
        for (i = 0; i<aux1; i++)
        {
          lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
          lFaceVelR[i] = MAT2D(nel, i, faceVelR, aux2);
          lFaceVelL[i] = MAT2D(nel, i, faceVelL, aux2);
          lFacePresR[i] = MAT2D(nel, i, facePresR, aux2);
          lFacePresL[i] = MAT2D(nel, i, facePresL, aux2);
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for (j = 0; j<ndm; j++)
          {
            MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
        }
/*...................................................................*/

/*... loop na celulas vizinhas*/    
        for(i=0;i<aux1;i++)
        {
          vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
          lViz[i] = vizNel;
          if( vizNel != -2)
          {
            lVolume[i]    = gVolume[vizNel]; 
            lGeomType[i]  = geomType[vizNel];
            lMat          = mat[vizNel]-1;
            lId[i]        = id[vizNel] - 1;
/*...*/
            MAT2D(i,0,lDensity, DENSITY_LEVEL) = MAT2D(vizNel
                              , TIME_N_MINUS_2,density ,DENSITY_LEVEL);
            MAT2D(i,1,lDensity, DENSITY_LEVEL) = MAT2D(vizNel
                              , TIME_N_MINUS_1,density ,DENSITY_LEVEL);
            MAT2D(i,2,lDensity, DENSITY_LEVEL) = MAT2D(vizNel
                              ,TIME_N,density ,DENSITY_LEVEL);
/*...................................................................*/
           
/*...*/
            lPres[i] = pres[vizNel];
            lTemp[i] = temp[vizNel];
/*...................................................................*/ 

            for(j=0;j<ndm;j++)
            {
              MAT2D(i,j,lGradPres,ndm) = MAT2D(vizNel,j,gradPres,ndm);
              MAT2D(i,j,lVel     ,ndm) = MAT2D(vizNel,j,vel     ,ndm);
              MAT2D(i,j,lDfield  ,ndm) = MAT2D(vizNel,j,dField  ,ndm);
            }

          }
        }  
/*...................................................................*/

/*...*/ 
        if(fWallModel)   
          for(i=0;i<NWALLPAR;i++)
            lWallPar[i] = MAT2D(nel,i,wallPar,NWALLPAR);
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibSimplePresLm(loadsVel,loadsPres 
											,diffPres   ,eMass
                      ,lGeomType  
                      ,lViz       ,lId           
                      ,lKsi       ,lmKsi
                      ,lEta       ,lfArea 
                      ,lNormal    ,lVolume
                      ,lXm        ,lXmcc
                      ,lDcca      ,lDensity
                      ,lvSkew     ,lmvSkew
                      ,lA         ,lB
                      ,&lRcell    ,ddt
                      ,lFaceVelR  ,lFaceVelL            
                      ,lFacePresR ,lFacePresL   
                      ,lPres      ,lGradPres    
                      ,lVel       ,lDfield 
                      ,lTemp      ,lWallPar
                      ,nen[nel]   ,nFace[nel] 
                      ,ndm        ,lib   
                      ,nel);      
/*...................................................................*/

/*... residuo da celula*/
        if(calRcell)
          rCell[nel] = lRcell;  
/*...................................................................*/
 
/*...*/
        assbly(ia      ,ja
            ,a           ,ad              
            ,b           ,lId 
            ,lA          ,lB
            ,nEq         ,nEqNov
            ,nAd         ,nAdR     
            ,nFace[nel]  ,1 
            ,storage     ,forces
            ,matrix      ,unsym);   
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 02/08/2016                                   *
 * Data de modificaco : 18/07/2018                                   * 
 *-------------------------------------------------------------------* 
 * SIMPLENONORTHPRES: correcao nao ortogonal para as pressoes        *
 * de correcao                                                       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * el        -> conetividade dos celulas                             * 
 * nelcon    -> vizinhos dos elementos                               * 
 * nen       -> numero de nos por celulas                            * 
 * nFace     -> numero de faces por celulas                          * 
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * gVolume -> volumes das celulas                                    *
 * gDcca   -> menor distancia do centroide a faces desta celula      *
 * gXmCc   -> vetores que unem o centroide aos pontos medios das     *
 * gCc     -> centroide do elemento                                  * 
 * fModKsi -> o mudolo do vetores que unem centroide da celula       * 
 *            central aos vizinhos destas                            *
 * fKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        * 
 * fEta    -> vetores paralelos as faces das celulas                 *
 * fArea   -> modulo do vetor eta                                    *
 * fNormal -> vetores normais as faces das celulas                   *
 * fXm     -> pontos medios das faces das celulas                    *
 * fModvSkew -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face                      * 
 * fvSkew  -> vetor entre o ponto medio a intersecao que une os      *
 *            centrois compartilhado nessa face                      *
 * calType -> tipo de calculo das celulas                            * 
 * geomType-> tipo geometrico das celulas                            * 
 * prop    -> propriedades dos material                              * 
 * mat     -> material por celula                                    * 
 * density   -> massa especifica com variacao temporal               * 
 * b         -> vetor de forcas                                      * 
 * id        -> numera das equacoes                                  * 
 * facePresR -> restricoes por elemento de pressao                   * 
 * pres      -> campo de pressao conhecido                           * 
 * gradPres  -> gradiente da pressao de correcao conhecido           * 
 * dField    -> matriz D do metodo simple ( volume/A(i,i) )          * 
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * ndf       -> graus de liberdade                                   * 
 * numel     -> numero de total de celulas                           * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * b         -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void simpleNonOrthPres(Diffusion *diffPres
               , INT    *RESTRICT el        , INT *RESTRICT nelcon 
               , short  *RESTRICT nen       , short  *RESTRICT nFace
               , INT *RESTRICT cellFace     , INT *RESTRICT fOwner
               , DOUBLE *RESTRICT gVolume   , DOUBLE *RESTRICT gDcca
               , DOUBLE *RESTRICT gXmcc     , DOUBLE *RESTRICT gCc
               , DOUBLE *RESTRICT fModKsi   , DOUBLE *RESTRICT fKsi      
               , DOUBLE *RESTRICT fEta      , DOUBLE *RESTRICT fArea 
               , DOUBLE *RESTRICT fNormal   , DOUBLE *RESTRICT fXm       
               , DOUBLE *RESTRICT fModvSkew , DOUBLE *RESTRICT fvSkew                  
               , short  *RESTRICT geomType  , DOUBLE *RESTRICT prop
               , short  *RESTRICT calType   , short  *RESTRICT mat
               , DOUBLE *RESTRICT density
               , DOUBLE *RESTRICT b         , INT    *RESTRICT id
               , short  *RESTRICT facePresR , DOUBLE *RESTRICT pres      
               , DOUBLE *RESTRICT gradPres  , DOUBLE *RESTRICT dField 
               , short maxNo                , short maxViz
               , short ndm                  , INT numel)
{ 

  short i,j;
  short  aux1, aux2, lMat;
  short  lGeomType[MAX_NUM_FACE + 1];
  INT nel, vizNel;
  INT idFace, cellOwner, ch;
  INT lViz[MAX_NUM_FACE], lNeq;
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM],lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM],lVolume[MAX_NUM_FACE+1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM],lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE],lvSkew[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDensity[(MAX_NUM_FACE+1)];
  DOUBLE lB;
  DOUBLE lPres[MAX_NUM_FACE+1];
  DOUBLE lProp[(MAX_NUM_FACE+1)*MAXPROP];
  DOUBLE lGradPres[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lDfield[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lCc[(MAX_NUM_FACE+1)*MAX_NDM];



/*... loop nas celulas*/
  aux2    = maxViz+1;
  for(nel=0;nel<numel;nel++)
  {
/*...*/
    aux1    = nFace[nel];
/*... elementos com equacoes*/
    if(MAT2D(nel,aux1,facePresR ,aux2) != PCCELL)
    {

/*... loop na celula central*/    
      lMat            = mat[nel]-1;
      lVolume[aux1]   = gVolume[nel]; 
      lGeomType[aux1] = geomType[nel];
      lDensity[aux1]  = MAT2D(nel,2   ,density ,DENSITY_LEVEL);
 
      lPres[aux1]     = pres[nel];
/*...................................................................*/
      
/*...*/
      for(j=0;j<MAXPROP;j++)
        MAT2D(aux1,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
/*...................................................................*/
      
/*...*/
      for(j=0;j<ndm;j++)
      {
        MAT2D(aux1,j,lGradPres,ndm) = MAT2D(nel,j,gradPres,ndm);
        MAT2D(aux1,j,lCc      ,ndm) = MAT2D(nel,j,gCc      ,ndm);
        MAT2D(aux1,j,lDfield  ,ndm) = MAT2D(nel,j,dField  ,ndm);
      }
/*...................................................................*/

/*...*/
      for(i=0;i<aux1;i++)
      {
        lDcca[i]      = MAT2D(nel,i,gDcca   ,maxViz);
/*... propriedades por face*/
        idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
        cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
        ch = OWNER(cellOwner, nel);
        lmKsi[i] = fModKsi[idFace];
        lfArea[i] = fArea[idFace];
        lmvSkew[i] = fModvSkew[idFace];
        for(j=0;j<ndm;j++){
          MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,gXmcc  ,maxViz,ndm);
/*... propriedades por face*/
          MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
          MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
          MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
          MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
          MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
        }
      }
/*...................................................................*/

/*... loop na celulas vizinhas*/    
      for(i=0;i<aux1;i++){
        vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
        lViz[i] = vizNel;
        if( vizNel != -2) {
          lVolume[i]    = gVolume[vizNel]; 
          lGeomType[i]  = geomType[vizNel];
          lDensity[i]   = MAT2D(vizNel,0   ,density ,DENSITY_LEVEL);
          lMat          = mat[vizNel]-1;
          lPres[i]      = pres[vizNel];

          for(j=0;j<ndm;j++){
            MAT2D(i,j,lGradPres,ndm) = MAT2D(vizNel,j,gradPres,ndm);
            MAT2D(i,j,lCc      ,ndm) = MAT2D(vizNel,j,gCc      ,ndm);
            MAT2D(i,j,lDfield  ,ndm) = MAT2D(vizNel,j,dField  ,ndm);
          }
  
          for(j=0;j<DIFPROP;j++)
            MAT2D(i,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
          
        }
      }  
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
        cellLibSimpleNonOrthPres(*diffPres
                         ,lGeomType 
                         ,lProp      ,lViz                  
                         ,lKsi       ,lmKsi
                         ,lEta       ,lfArea 
                         ,lNormal    ,lVolume
                         ,lXm        ,lXmcc
                         ,lDcca      ,lDensity
                         ,lvSkew     ,lmvSkew
                         ,&lB        
                         ,lPres      ,lGradPres 
                         ,lDfield    ,lCc
                         ,nen[nel]   ,nFace[nel] 
                         ,ndm        ,nel);    
/*...................................................................*/

/*...*/
      lNeq    = id[nel] - 1;
      b[lNeq] = lB;
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
 * loads     -> definicoes de cargas                                 * 
 * cc        -> centroide das celulas                                * 
 * faceR     -> restricoes por elmento                               * 
 * faceL     -> carga por elemento                                   * 
 * volume    -> volume das celulas                                   * 
 * id        -> numera das equacoes                                  * 
 * u         -> solucao                                              * 
 * f         -> vetor de forcas                                      * 
 * numel     -> numero de elementos                                  * 
 * ndf       -> graus de liberade                                    * 
 * ndm       -> numero de dimensao                                   * 
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
void cellPload(Loads *loads           ,DOUBLE *RESTRICT cc 
              ,short  *RESTRICT faceR ,short *RESTRICT faceL
              ,DOUBLE *RESTRICT volume,INT *RESTRICT id 
              ,DOUBLE *RESTRICT u     ,DOUBLE *RESTRICT f
              ,INT const numel        ,short const ndf
              ,short const ndm        ,short const maxViz)
{
  INT nel,lNeq;
  DOUBLE xx[3],uT[MAX_NDF];
  short carg,j;
  short col = maxViz + 1;

/*...*/  
  for(nel = 0; nel < numel;nel++){
    carg = MAT2D(nel,maxViz,faceR,col);
/*... variavel prescrita no dominio*/
    if(carg == PCCELL){
      carg= MAT2D(nel,maxViz,faceL,col)-1;
/*... valor prescrito na celula constante*/
      if( loads[carg].type == CONST)
        for(j = 0; j< ndf;j++)
          MAT2D(nel,j,u,ndf) = loads[carg].par[j];
/*...................................................................*/
    }

/*... carga na celula*/
    else if( carg == SCCELL){
      carg= MAT2D(nel,maxViz,faceL,col)-1;
/*... carga constante*/
      if( loads[carg].type == CONST)
        for(j = 0; j< ndf;j++){
          lNeq = MAT2D(nel,j,id,ndf) - 1;
          if( lNeq > -1)
            MAT2D(lNeq,j,f,ndf) 
            = volume[nel]*loads[carg].par[j];
        }
/*...................................................................*/

/*... carga senoidal*/
      else if( loads[carg].type == SINBC){
/*...*/
        if(ndm == 2){
          xx[0] = MAT2D(nel,0,cc ,ndm); 
          xx[1] = MAT2D(nel,1,cc ,ndm); 
          xx[2] = 0.0e0;                             
        }
        else{              
          xx[0] = MAT2D(nel,0,cc ,ndm); 
          xx[1] = MAT2D(nel,1,cc ,ndm); 
          xx[2] = MAT2D(nel,2,cc ,ndm); 
        }
        loadSenProd(uT,loads[carg].par,xx);
/*...................................................................*/
        for(j = 0; j< ndf;j++){
          lNeq = MAT2D(nel,j,id,ndf) - 1;
          if( lNeq > -1){
            MAT2D(lNeq,j,f,ndf) = volume[nel]*uT[j];
          }
        }
/*...................................................................*/
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
 * Data de modificaco : 04/07/2016                                   * 
 *-------------------------------------------------------------------* 
 * CELLPLOADSIMPLE : Carregamento prescrito por centro da celula no  * 
 * metodo simple                                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads     -> definicoes de cargas                                 * 
 * cc        -> centroide das celulas                                * 
 * faceR     -> restricoes de pressao por celula                     * 
 * faceL     -> carga de pressao por celula                          * 
 * volume    -> volume das celulas                                   * 
 * idVel     -> numera das equacoes Vel                              * 
 * idPres    -> numera das equacoes Vel                              * 
 * vel       -> velocidade                                           * 
 * pres      -> pressoes                                             * 
 * fVel      -> vetor de forcas das velocidades                      * 
 * fPres     -> vetor de forcas da pressoes                          * 
 * numel     -> numero de elementos                                  * 
 * ndf       -> graus de liberade                                    * 
 * ndm       -> numero de dimensao                                   * 
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
void cellPloadSimple(Loads *loadsPres       ,DOUBLE *RESTRICT cc 
                ,short  *RESTRICT faceRpres ,short *RESTRICT faceLpres
                ,DOUBLE *RESTRICT volume
                ,INT *RESTRICT idVel    ,INT *RESTRICT idPres
                ,DOUBLE *RESTRICT vel   ,DOUBLE *RESTRICT pres
                ,DOUBLE *RESTRICT fVel  ,DOUBLE *RESTRICT fPres
                ,INT const numel        ,short const ndf
                ,short const ndm        ,short const maxViz)
{
  INT nel;
  short carg;
  short col = maxViz + 1;

/*...*/  
  for(nel = 0; nel < numel;nel++){
    carg = MAT2D(nel,maxViz,faceRpres,col);
/*... variavel prescrita no dominio*/
    if(carg == PCCELL){
      carg = MAT2D(nel,maxViz,faceLpres,col)-1;
/*... valor prescrito constante para pressao na celula*/
      if( loadsPres[carg].type == CONST)
        pres[nel] = loadsPres[carg].par[0];
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}  
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------*
 * UPDATECELLVALUE:atualizacao dos valores das variaveis das celulas *
 * com os valores das respectivas equacoes                           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * u       -> variavel nas celulas                                   * 
 * x       -> solucao do sistema                                     * 
 * id      -> numera das equacoes                                    * 
 * iNeq    -> mapa de equacoes de interface                          *
 * numel   -> numero de elementos                                    * 
 * ndf     -> graus de liberdade                                     * 
 * fAdd    -> true add false sobreescreve                            * 
 * fCom    -> comunica os valores x entre as particoes               * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * u      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 void updateCellValue(DOUBLE *RESTRICT u,DOUBLE *RESTRICT x
                 ,INT *RESTRICT id      ,Interface *iNeq
                 ,INT const numel       ,short const ndf
                 ,bool const fAdd       ,bool const fCom)
{
  INT nel,lNeq;
  short jNdf;

/*... obtem os valores de x das equacoes em overlaping*/  
  if(fCom)
    comunicateNeq(iNeq,x);
/*.................................................................*/  
  
/*...*/
  if(fAdd)    
    for(nel=0;nel<numel;nel++){
      for(jNdf = 0;jNdf<ndf;jNdf++){ 
        lNeq = MAT2D(nel,jNdf,id,ndf) - 1;
        if( lNeq > -1)
          MAT2D(nel,jNdf,u,ndf) += MAT2D(lNeq,jNdf,x,ndf);
      }
    }
/*.................................................................*/  

/*...*/
  else
    for(nel=0;nel<numel;nel++){
      for(jNdf = 0;jNdf<ndf;jNdf++){ 
        lNeq = MAT2D(nel,jNdf,id,ndf) - 1;
        if( lNeq > -1)
          MAT2D(nel,jNdf,u,ndf) = MAT2D(lNeq,jNdf,x,ndf);
      }
    }
/*.................................................................*/  
  
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 20/08/2019                                   * 
 *-------------------------------------------------------------------*
 * UPDATECELLVALUE:atualizacao dos valores das variaveis das celulas *
 * com os valores das respectivas equacoes                           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * u       -> variavel nas celulas                                   * 
 * x       -> solucao do sistema                                     * 
 * id      -> numera das equacoes                                    * 
 * iNeq    -> mapa de equacoes de interface                          *
 * numel   -> numero de elementos                                    * 
 * ndf     -> graus de liberdade                                     * 
 * fAdd    -> true add false sobreescreve                            * 
 * fCom    -> comunica os valores x entre as particoes               * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * u      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 * x = [ x1 x2 ... xNeq y1 y2 ... yNeq ... z1 z2 ... zNeq            *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void updateCellValueBlock(DOUBLE *RESTRICT u    ,DOUBLE *RESTRICT x
                          ,INT *RESTRICT id     ,Interface *iNeq
                          ,INT const numel      ,INT const nEq
                          ,short const ndf                        
                          ,bool const fAdd      ,bool const fCom)
{
  INT nel,lNeq;
  short jNdf;

/*... obtem os valores de x das equacoes em overlaping*/  
  if(fCom)
  {
    for(jNdf = 0;jNdf<ndf;jNdf++)
      comunicateNeq(iNeq,&x[jNdf*nEq]);
  }
/*.................................................................*/  
  
/*...*/
  if(fAdd)    
    for(nel=0;nel<numel;nel++)
    {
      lNeq = id[nel] - 1;
      for(jNdf = 0;jNdf<ndf;jNdf++)
      {        
        if( lNeq > -1) 
          MAT2D(nel,jNdf,u,ndf) += x[lNeq+jNdf*nEq];
      }
    }
/*.................................................................*/  

/*...*/
  else
    for(nel=0;nel<numel;nel++)
    {
      lNeq = id[nel] - 1;
      for(jNdf = 0;jNdf<ndf;jNdf++)
      { 
        if( lNeq > -1)
          MAT2D(nel,jNdf,u,ndf) = x[lNeq+jNdf*nEq];
      }
    }
/*.................................................................*/  
  
/*for(nel=0;nel<numel;nel++)
  {
    fprintf(fileLogDebug," %4d ",nel);
    for(jNdf = 0;jNdf<ndf;jNdf++)
      fprintf(fileLogDebug," %e ",MAT2D(nel,jNdf,u,ndf));
    fprintf(fileLogDebug,"\n");
  }
  mpiStop();
  exit(0);
*/
  
}
/*********************************************************************/


/********************************************************************* 
 * Data de criacao    : 10/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * UPDATECELLVALUE:atualizacao dos valores das variaveis das celulas *
 * com os valores das respectivas equacoes                           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * u       -> variavel nas celulas                                   * 
 * x       -> solucao do sistema                                     * 
 * id      -> numera das equacoes                                    * 
 * iNeq    -> mapa de equacoes de interface                          *
 * numel   -> numero de elementos                                    * 
 * nEq     -> numero de equacoes                                     * 
 * ndf     -> graus de liberdade                                     * 
 * fAdd    -> true add false sobreescreve                            * 
 * fCom    -> comunica os valores x entre as particoes               * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * u      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS: f | ru(1) ru(2) ... ru(neq) |                                *
 *        | rv(1) rv(2) ... rv(neq) |                                *
 *        | rw(1) rw(2) ... bw(neq) |                                * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 void updateCellValueSimple(DOUBLE *RESTRICT u, DOUBLE *RESTRICT x
                           ,INT *RESTRICT id  , Interface *iNeq
                           ,INT const numel   , INT const nEq   
                           ,short const ndf
                           ,bool const fAdd   , bool const fCom)
{
  INT nel,lNeq;
  short jNdf;

/*... obtem os valores de x das equacoes em overlaping*/  
  if(fCom)
    comunicateNeq(iNeq,x);
/*.................................................................*/  
  
/*...*/
  if(fAdd)    
    for(jNdf = 0;jNdf<ndf;jNdf++){ 
      for(nel=0;nel<numel;nel++){
        lNeq = id[nel] - 1;
        if( lNeq > -1)
          MAT2D(jNdf,nel,u,numel) += MAT2D(jNdf,lNeq,x,nEq);
      }
    }
/*.................................................................*/  

/*...*/
  else
    for(jNdf = 0;jNdf<ndf;jNdf++){ 
      for(nel=0;nel<numel;nel++){
        lNeq = id[nel] - 1;
        if( lNeq > -1)
          MAT2D(jNdf,nel,u,numel) = MAT2D(jNdf,lNeq,x,nEq);
      }
    }
/*.................................................................*/  
  
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 11/05/2019                                   * 
 *-------------------------------------------------------------------* 
 * INTERCELLNODE: interpolacao dos valores das celulas para o no da  *
 * malha                                                             *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m       -> vetor de memoria principal                             * 
 * loads   -> definicao de cargas                                    * 
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * noU     -> nao definido                                           * 
 * elU     -> valores nas celulas                                    * 
 * el      -> conectividades das celulas                             * 
 * geomType-> tipo geometrico das celulas                            * 
 * cc      -> centroide das celulas                                  * 
 * x       -> coordenadas                                            * 
 * xm      -> pontos medios das faces das celulas                    * 
 * nen     -> numero de nos por celulas                              * 
 * faceR   -> restricoes por elmento                                 * 
 * faceL   -> carga por elemento                                     * 
 * iNo     -> mapa de nos de interface                               *
 * numelNov-> numero de elementos sem sobreposicoes                  * 
 * numel   -> numero de elementos                                    * 
 * nNodeNov-> numero de nos sem sobreposicoes                        * 
 * nNode   -> numero de nos                                          * 
 * maxNo   -> numero de nos por celula maximo da malha               * 
 * ndf1    -> graus de liberdade linha  (tensor)                     * 
 * ndf2    -> graus de liberdade coluna (tensor)                     * 
 * ndm     -> numero de dimensao                                     * 
 * fBc     -> forca condicao de controno conhecida                   * 
 * type    -> tipo de interpolacao                                   * 
 *            1 - media simples                                      * 
 *            2 - media ponderada                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * u      -> atualizado | u1 u2 u3 |                                 * 
 *                      | v1 v2 v3 |                                 *
 *                      | w1 w2 w3 |                                 *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 void interCellNode(Memoria *m             ,Loads *loads 
                   ,INT *RESTRICT cellFace, INT *RESTRICT fOwner
                   ,DOUBLE *RESTRICT noU   ,DOUBLE *RESTRICT elU
                   ,INT *RESTRICT el       ,short  *RESTRICT geomType 
                   ,DOUBLE *RESTRICT cc    ,DOUBLE *RESTRICT x
                   ,DOUBLE *RESTRICT fXm  
                   ,short *RESTRICT nen    ,short *RESTRICT nFace
                   ,short  *RESTRICT faceR ,short *RESTRICT faceL 
                   ,InterfaceNo *iNo 
                   ,INT const numelNov     ,INT const numel        
                   ,INT const nNodeNov     ,INT const nNode
                   ,short const maxNo      ,short const maxViz     
                   ,short const ndf1       ,short const ndf2
                   ,short const ndm
                   ,bool const fBc         ,short const type)

{
  int *md=NULL;
  DOUBLE *mdf=NULL;
  bool *flag=NULL;
  INT    idFace;
  short i,j,k,l,n,nodeFace,aux=maxViz+1;
  INT nel,no1,no[4];
  short  isNod[MAX_SN],nCarg,ty,typed;
  DOUBLE dist, dx;
  DOUBLE uT[MAX_NDF], xx[4], par[MAXLOADPARAMETER];;
  
  switch(type)
  {
/*... media simple*/
    case 1:
      if(ndf2 != 1)
      {
        ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__
                  ,"Opcao nao implementada",EXIT_PROG);
      }
/*...*/
      HccaAlloc(int,m,md,nNodeNov,"md",false);
      zero(md,nNodeNov,"int");
      zero(noU,ndf1*nNodeNov,DOUBLEC);
/*...................................................................*/

/*...*/
      if( ndf1 == 1 )
        for(nel = 0; nel < numelNov; nel++)
        {
          for(j = 0; j < nen[nel];j++)
          {
            no1 = MAT2D(nel,j,el,maxNo) - 1;
            noU[no1] += elU[nel];
            md[no1]++;
          }
        }
/*...................................................................*/

/*...*/
      else
        for(nel = 0; nel < numelNov; nel++)
        {
          for(j = 0; j < nen[nel];j++)
          {
            no1 = MAT2D(nel,j,el,maxNo) - 1;
            for(k = 0; k   < ndf1;k++)
              MAT2D(no1,k,noU,ndf1) += MAT2D(nel,k,elU,ndf1);
            md[no1]++;
          }
        }
/*...................................................................*/

/*... comunicacao*/
      if(mpiVar.nPrcs > 1) 
      {
        dComunicateNod(iNo,noU,ndf1,ndf2);
        iComunicateNod(iNo,md ,1   ,1);
      }
/*...................................................................*/

/*...*/
      if( ndf1 == 1 )
        for(no1 = 0; no1 < nNodeNov; no1++)
         noU[no1] /= md[no1];
/*...................................................................*/

/*...*/
      else
        for(no1 = 0; no1 < nNodeNov; no1++){
          for(k = 0; k < ndf1; k++)
            MAT2D(no1,k,noU,ndf1) /= md[no1];
        }
/*...................................................................*/
          
/*...*/
      HccaDealloc(m,md,"md",false);
/*...................................................................*/
    break;
/*...................................................................*/

/*... media ponderada*/
    case 2:
/*...*/
      HccaAlloc(DOUBLE,m,mdf,nNodeNov,"mdf",false);
      zero(mdf,nNodeNov,"double");
      zero(noU,ndf1*nNodeNov,DOUBLEC);
/*...................................................................*/

/*...*/
      if( ndf1 == 1 )
      {
/*...*/
        for(nel = 0; nel < numelNov; nel++)
        {
          for(j = 0; j < nen[nel];j++)
          {
            no1 = MAT2D(nel,j,el,maxNo) - 1;
            dist = 0.e0;
            for(k = 0; k   < ndm;k++){
              dx = MAT2D(no1,k,x,ndm) - MAT2D(nel,k,cc,ndm);
              dist += dx*dx;
            }
            dist = 1.e0/sqrt(dist);
            noU[no1]+= elU[nel]*dist;
            mdf[no1]+= dist;
          }
        }
/*...................................................................*/
      }
/*...................................................................*/

/*...*/
      else
      {
/*...*/
        if(ndf2 == 1)
        {
          for(nel = 0; nel < numelNov; nel++)
          {
            for(j = 0; j < nen[nel];j++)
            {
              no1 = MAT2D(nel,j,el,maxNo) - 1;
              dist = 0.e0;
              for(k = 0; k   < ndm;k++)
              {
                dx = MAT2D(no1,k,x,ndm) - MAT2D(nel,k,cc,ndm);
                dist += dx*dx;
              }
              dist = 1.e0/sqrt(dist);
              for(k = 0; k < ndf1;k++)
                MAT2D(no1,k,noU,ndf1) += MAT2D(nel,k,elU,ndf1)*dist;
              mdf[no1]+=dist;
            }
          }
        }
/*...................................................................*/

/*...*/
        else
        {
          for(nel = 0; nel < numelNov; nel++)
          {
            for(j = 0; j < nen[nel];j++)
            {
              no1 = MAT2D(nel,j,el,maxNo) - 1;
              dist = 0.e0;
              for(k = 0; k   < ndm;k++)
              {
                dx = MAT2D(no1,k,x,ndm) - MAT2D(nel,k,cc,ndm);
                dist += dx*dx;
              }
              dist = 1.e0/sqrt(dist);
              for(k = 0; k < ndf1;k++)
                for(l = 0; l < ndf2;l++)
                  MAT3D(no1,k,l,noU,ndf1,ndf2) 
                 += MAT3D(nel,k,l,elU,ndf1,ndf2)*dist;
              mdf[no1]+=dist;
            }
          }
        }
/*...................................................................*/
      }
/*...................................................................*/

/*... comunicacao*/
      if(mpiVar.nPrcs > 1) 
      {
        dComunicateNod(iNo,noU,ndf1,ndf2);
        dComunicateNod(iNo,mdf,1  ,1);
      }
/*...................................................................*/

/*...*/
      if( ndf1 == 1 )
        for(no1 = 0; no1 < nNodeNov; no1++)
          noU[no1] /= mdf[no1];
/*...................................................................*/

/*...*/
      else
      {
/*...*/
        if( ndf2 == 1 )
          for(no1 = 0; no1 < nNodeNov; no1++)
            for(k = 0; k < ndf1; k++)
              MAT2D(no1,k,noU,ndf1) /= mdf[no1];
/*...................................................................*/
        else
          for(no1 = 0; no1 < nNodeNov; no1++)
            for(k = 0; k < ndf1; k++)
              for(l = 0; l < ndf2;l++)
                MAT3D(no1,k,l,noU,ndf1,ndf2) /= mdf[no1];
/*...................................................................*/
      }          
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
    HccaAlloc(int ,m,md  ,nNode,"md",false);
    HccaAlloc(bool,m,flag,nNode,"flag",false);
    zero(md  ,nNode,"int");
    zero(flag,nNode,"bool");
    for(nel = 0; nel < numel; nel++)
      for(i = 0; i < nFace[nel]; i++)
      {

        idFace    = MAT2D(nel   , i, cellFace, maxViz) - 1;

        if(MAT2D(nel,i,faceR,aux)>0)
        {
/*...*/
          xx[0] = MAT2D(idFace, 0, fXm, ndm);
          xx[1] = MAT2D(idFace, 1, fXm, ndm);
          xx[2] = 0.0e0;                             
          if(ndm == 3) xx[2] = MAT2D(idFace, 2, fXm, ndm);
/*...................................................................*/
          nCarg=MAT2D(nel,i,faceL,aux)-1;
          typed = loads[nCarg].type;
/*... valor pescrito */
          if( typed == DIRICHLETBC || typed == MOVEWALL){
            ty = geomType[nel];
            nodeFace =  sn(isNod,ty,nel); 
            for(n=0;n<nodeFace ;n++){
              no[n] = MAT2D(i,n,isNod,nodeFace );
              no[n] = MAT2D(nel,no[n],el,maxNo) - 1;
              if(flag[no[n]] == false){
                flag[no[n]] = true;
                for(k = 0; k   < ndf1;k++)
                  MAT2D(no[n],k,noU,ndf1) = 0.e0;
              }
              getLoads(par,&loads[nCarg],xx);
              for(k = 0; k   < ndf1;k++) 
                MAT2D(no[n],k,noU,ndf1) += par[k];
              md[no[n]]++;
            }
          }
/*...................................................................*/

/*... valor pescrito */
          else if( typed == INLET )
          {
            ty = geomType[nel];
            nodeFace =  sn(isNod,ty,nel); 
            for(n=0;n<nodeFace ;n++){
              no[n] = MAT2D(i,n,isNod,nodeFace );
              no[n] = MAT2D(nel,no[n],el,maxNo) - 1;
              if(flag[no[n]] == false)
              {
                flag[no[n]] = true;
                for(k = 0; k   < ndf1;k++)
                  MAT2D(no[n],k,noU,ndf1) = 0.e0;
              }
              getLoads(par,&loads[nCarg],xx);
              for(k = 0; k   < ndf1;k++) 
                MAT2D(no[n],k,noU,ndf1) += par[k];
              md[no[n]]++;
            }
          }
/*...................................................................*/

/*... valor pescrito  senoidal*/
          else if( typed == SINBC)
          {
            ty = geomType[nel];
            nodeFace =  sn(isNod,ty,nel); 
            for(n=0;n<nodeFace ;n++)
            {
              no[n] = MAT2D(i,n,isNod,nodeFace );
              no[n] = MAT2D(nel,no[n],el,maxNo) - 1;
              if(flag[no[n]] == false)
              {
                flag[no[n]] = true;
                for(k = 0; k   < ndf1;k++)
                  MAT2D(no[n],k,noU,ndf1) = 0.e0;
              }
              loadSenProd(uT,loads[nCarg].par,xx);
              for(k = 0; k   < ndf1;k++) 
                MAT2D(no[n],k,noU,ndf1) += uT[k];
              md[no[n]]++;
            }
/*...................................................................*/
          }
/*...................................................................*/
        }
/*...................................................................*/

/*... condicao de parede estatica*/
        else if(MAT2D(nel,i,faceR,aux) == STATICWALL)
        {
            ty = geomType[nel];
            nodeFace =  sn(isNod,ty,nel); 
            for(n=0;n<nodeFace ;n++)
            {
              no[n] = MAT2D(i,n,isNod,nodeFace );
              no[n] = MAT2D(nel,no[n],el,maxNo) - 1;
              if(flag[no[n]] == false)
              {
                flag[no[n]] = true;
                for(k = 0; k   < ndf1;k++)
                  MAT2D(no[n],k,noU,ndf1) = 0.e0;
              }
              md[no[n]]++;
            }
        }
/*...................................................................*/
      } 
/*...*/
    for(no1 = 0; no1 < nNodeNov; no1++)
      if(flag[no1])
        for(k = 0; k < ndf1; k++)
          MAT2D(no1,k,noU,ndf1) /= md[no1];
/*...................................................................*/

/*...*/
    HccaDealloc(m,flag,"flag",false);
    HccaDealloc(m,md  ,"md"  ,false);
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 04/05/2019                                   *
 *-------------------------------------------------------------------*
 * RCGRADU: calculo do gradiente de um campo escalar ou vetorial     * 
 * conhecido.                                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m       -> vetor de memoria principal                             * 
 * loads   -> definicao de cargas                                    * 
 * el      -> conetividade dos celulas                               * 
 * cc      -> centroide das celulas                                  * 
 * x       -> cordenadas dos pontos                                  * 
 * nelcon  -> vizinhos dos elementos                                 * 
 * nen     -> numero de nos por celulas                              * 
 * nFace   -> numero de faces por celulas                            * 
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * geomType-> tipo geometrico das celulas                            * 
 * prop    -> propriedades dos material                              * 
 * mat     -> material por celula                                    * 
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
 * gXm     -> pontos medios das faces das celulas                    * 
 * gXmcc     -> vetores que unem o centroides aos pontos medios das  * 
 *            faces da celula central                                * 
 * gDcca   -> menor distancia do centroide a faces desta celula      * 
 * faceR   -> restricoes por elmento                                 * 
 * faceL   -> carga por elemento                                     * 
 * u       -> solucao conhecida por celula (atualizado)              * 
 * gradU   -> gradiente da solucao         (desatualizado)           * 
 * nU      -> solucao conhecida por no     (desatualizado)           * 
 * rcGrad  -> config do reconstrucao de gradiente                    * 
 * maxNo   -> numero de nos por celula maximo da malha               * 
 * maxViz  -> numero vizinhos por celula maximo da malha             * 
 * ndm     -> numero de dimensoes                                    * 
 * ndf     -> graus de liberdade                                     * 
 * iNo     -> interface de nos                                       * 
 * iCel    -> interface de elementos                                 * 
 * numel   -> numero de total de celulas                             * 
 * numelNov-> numero de total de celulas sem sobreposicao            * 
 * nNode   -> numero de total de nos                                 * 
 * nNodeNov-> numero de total nos em celulas sem sobreposicao        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*
 *                  |gradU1|   | dU1/dx1 dU1/dx2 dU1/dx3 |           *
 * gradU(ndf,ndm) = | ...  | = |         ...             |           *
 *                  |gradUn|   | dUn/dx1 dUn/dx2 dUn/dx3 |           *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void rcGradU(Memoria *m                , Loads *loads
           , INT    *RESTRICT el       , INT    *RESTRICT nelcon 
           , DOUBLE *RESTRICT x
           , short  *RESTRICT nen      , short  *RESTRICT nFace
           , INT *RESTRICT cellFace    , INT *RESTRICT fOwner
           , DOUBLE *RESTRICT gVolume  , DOUBLE *RESTRICT gDcca
           , DOUBLE *RESTRICT fXmCc    , DOUBLE *RESTRICT gCc
           , DOUBLE *RESTRICT fModKsi  , DOUBLE *RESTRICT fKsi
           , DOUBLE *RESTRICT fEta     , DOUBLE *RESTRICT fArea
           , DOUBLE *RESTRICT fNormal  , DOUBLE *RESTRICT fXm
           , DOUBLE *RESTRICT fModvSkew, DOUBLE *RESTRICT fvSkew
           , short  *RESTRICT geomType , DOUBLE *RESTRICT prop
           , short  *RESTRICT calType  , short  *RESTRICT mat     
           , DOUBLE *RESTRICT lSquare  , DOUBLE *RESTRICT lSquareR
           , short  *RESTRICT faceR    , short *RESTRICT faceL  
           , DOUBLE *RESTRICT u        , DOUBLE *RESTRICT gradU              
           , DOUBLE *RESTRICT nU       , RcGrad *rcGrad 
           , short maxNo               , short maxViz
           , short ndf                 , short ndm
           , InterfaceNo *iNo          , Interface *iCel
           , INT numelNov              , INT numel 
           , INT nNodeNov              , INT nNode)
{
  short i,j;
  short nThreads = ompVar.nThreadsGrad;
  INT nel,no,vizNel;
/*... variavel local */
  DOUBLE lKsi[MAX_NUM_FACE*3],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*3],lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*3],lVolume[MAX_NUM_FACE+1];
  DOUBLE lvSkew[MAX_NUM_FACE*3];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM],lXmcc[MAX_NUM_FACE*3];
  short  lFaceR[MAX_NUM_FACE+1];
  short  lFaceL[MAX_NUM_FACE+1];
  DOUBLE lu[(MAX_NUM_FACE+1)*MAX_NDF];
  DOUBLE lnU[(MAX_NUM_NODE)*MAX_NDF];
  DOUBLE lGradU[MAX_NDM*MAX_NDF];
  DOUBLE lLsquare[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lLsquareR[2*MAX_NDM];
  DOUBLE lProp[MAXPROP],lDcca[MAX_NUM_FACE];
  INT    lViz[MAX_NUM_FACE];
  short  aux1,aux2,lMat;
  short  isNod[MAX_SN],ty;
  INT idFace, cellOwner, ch;

/*... reconstrucao de gradiente Green-Gauss nodal*/
  if(rcGrad->type ==  RCGRADGAUSSN){
    interCellNode(m         ,loads  
                 ,cellFace  ,fOwner
                 ,nU        ,u
                 ,el        ,geomType
                 ,gCc       ,x
                 ,fXm
                 ,nen       ,nFace
                 ,faceR     ,faceL
                 ,iNo       
                 ,numelNov  ,numel     
                 ,nNodeNov  ,nNode  
                 ,maxNo     ,maxViz
                 ,ndf       ,1  
                 ,ndm    
                 ,true      ,2);                
  }
/*.....................................................................*/

/*...*/
  if (ompVar.fGrad) {
/*... */
    aux2 = maxViz + 1;
#pragma omp parallel  for default(none) num_threads(nThreads)\
     private(nel,i,j,aux1,lMat,lVolume,lmKsi,lfArea,lDcca,lKsi,lEta\
          ,lNormal,lXm,lXmcc,lvSkew,vizNel,lViz,lFaceR,lFaceL,lu,lnU\
          ,lProp,lLsquare,lLsquareR,lGradU,ty,isNod,no\
          ,idFace,cellOwner,ch)\
     shared(aux2,ndm,ndf,numel,maxViz,nFace,mat,u,nU,lSquare,lSquareR\
         ,gVolume, geomType,prop,fModKsi,faceL,faceR,gradU,numelNov,rcGrad\
         ,fArea,gDcca,fKsi,fEta,fNormal,fXm,fXmCc,fvSkew,maxNo\
         ,nelcon,loadsVel,loadsPres,nen,loads,el\
         ,fOwner,cellFace)
    for (nel = 0; nel<numelNov; nel++) {
      aux1 = nFace[nel];

/*... loop na celula central*/
      lMat = mat[nel] - 1;
      lVolume[aux1] = gVolume[nel];
      lFaceR[aux1] = MAT2D(nel, aux1, faceR, aux2);
      lFaceL[aux1] = MAT2D(nel, aux1, faceL, aux2);

      for (i = 0; i<ndf; i++)
        MAT2D(aux1, i, lu, ndf) = MAT2D(nel, i, u, ndf);

      for (j = 0; j<MAXPROP; j++)
        lProp[j] = MAT2D(lMat, j, prop, MAXPROP);

/*... valor da funcao nodal nodias*/
      if (rcGrad->type == RCGRADGAUSSN)
      {
        for (i = 0; i<nen[nel]; i++)
        {
          no = MAT2D(nel, i, el, maxNo) - 1;
          for (j = 0; j<ndf; j++)
            MAT2D(i, j, lnU, ndf) = MAT2D(no, j, nU, ndf);
        }
      }

/*... leastSquare*/
      if (rcGrad->type == RCLSQUARE)
        for (i = 0; i<ndm; i++)
          for (j = 0; j<aux1; j++)
            MAT2D(i, j, lLsquare, aux1)
            = MAT3D(nel, i, j, lSquare, ndm, maxViz);
/*...................................................................*/

/*... leastSquare-QR*/
      if (rcGrad->type == RCLSQUAREQR)
      {
        for (i = 0; i<ndm; i++)
          for (j = 0; j<aux1; j++)
            MAT2D(i, j, lLsquare, aux1)
            = MAT3D(nel, i, j, lSquare, ndm, maxViz);
/*... R*/
        if (ndm == 2)
        {
          lLsquareR[0] = MAT2D(nel, 0, lSquareR, 3);
          lLsquareR[1] = MAT2D(nel, 1, lSquareR, 3);
          lLsquareR[2] = MAT2D(nel, 2, lSquareR, 3);
        }
        else if (ndm == 3)
        {
          lLsquareR[0] = MAT2D(nel, 0, lSquareR, 6);
          lLsquareR[1] = MAT2D(nel, 1, lSquareR, 6);
          lLsquareR[2] = MAT2D(nel, 2, lSquareR, 6);
          lLsquareR[3] = MAT2D(nel, 3, lSquareR, 6);
          lLsquareR[4] = MAT2D(nel, 4, lSquareR, 6);
          lLsquareR[5] = MAT2D(nel, 5, lSquareR, 6);
        }
      }
/*...................................................................*/

/*...*/
      for (i = 0; i<ndf; i++)
        for (j = 0; j<ndm; j++)
          MAT2D(i, j, lGradU, ndm) = MAT3D(nel, i, j, gradU, ndf, ndm);
/*...................................................................*/

/*...*/
      for (i = 0; i<aux1; i++)
      {
        lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
        lFaceR[i] = MAT2D(nel, i, faceR, aux2);
        lFaceL[i] = MAT2D(nel, i, faceL, aux2);

        for (j = 0; j<ndm; j++)
          MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, fXmCc, maxViz, ndm);

        vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
        lViz[i] = vizNel;
        if (vizNel != -2)
        {
          lVolume[i] = gVolume[vizNel];
          for (j = 0; j<ndf; j++)
            MAT2D(i, j, lu, ndf) = MAT2D(vizNel, j, u, ndf);
        }
      }
/*...................................................................*/

/*... propriedades por face*/
      for (i = 0; i<aux1; i++)
      {
        idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
        cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
        ch = OWNER(cellOwner, nel);
        lmKsi[i] = fModKsi[idFace];
        lfArea[i] = fArea[idFace];
        for (j = 0; j<ndm; j++)
        {
          MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
          MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
          MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
          MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
          MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
        }
      }
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
      ty = geomType[nel];
      sn(isNod, ty, nel);
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
      cellLibRcGrad(loads       , rcGrad
                  , lViz        , lProp
                  , lLsquare    , lLsquareR
                  , lKsi        , lmKsi
                  , lEta        , lfArea
                  , lNormal     , lVolume
                  , lvSkew      
                  , lXm         , lXmcc
                  , lDcca       
                  , lFaceR      , lFaceL
                  , lu          , lGradU
                  , lnU         , ty
                  , nFace[nel]  , ndm
                  , ndf
                  , isNod       , nel);
/*...................................................................*/

/*...*/
      for (i = 0; i<ndf; i++)
        for (j = 0; j<ndm; j++)
          MAT3D(nel, i, j, gradU, ndf, ndm) = MAT2D(i, j, lGradU, ndm);
/*...................................................................*/
    }
/*...................................................................*/
  }
/*.....................................................................*/

/*... */
  else
  {
/*... */
    aux2    = maxViz+1;
    for(nel=0;nel<numelNov;nel++){
      aux1    = nFace[nel];

/*... loop na celula central*/    
      lMat            = mat[nel]-1;
      lVolume[aux1]   = gVolume[nel]; 
      lFaceR[aux1]    = MAT2D(nel,aux1,faceR ,aux2);
      lFaceL[aux1]    = MAT2D(nel,aux1,faceL ,aux2);
      
      for(i=0;i<ndf;i++)
        MAT2D(aux1,i,lu,ndf) = MAT2D(nel,i,u,ndf);
     
      for(j=0;j<MAXPROP;j++)
        lProp[j]= MAT2D(lMat,j,prop,MAXPROP);

/*... valor da funcao nodal nodias*/    
      if(rcGrad->type ==  RCGRADGAUSSN)
      {
        for(i=0;i<nen[nel];i++)
        {
          no = MAT2D(nel,i,el,maxNo)-1;
          for(j=0;j<ndf;j++)
            MAT2D(i,j,lnU,ndf) = MAT2D(no,j,nU,ndf);
        }
      }

/*... leastSquare*/
      if(rcGrad->type ==  RCLSQUARE)
        for(i=0;i<ndm;i++)
          for(j=0;j<aux1;j++)
            MAT2D(i,j,lLsquare,aux1) 
            = MAT3D(nel,i,j,lSquare,ndm,maxViz); 
/*...................................................................*/

/*... leastSquare-QR*/
      if(rcGrad->type ==  RCLSQUAREQR)
      {
        for(i=0;i<ndm;i++)
          for(j=0;j<aux1;j++)
            MAT2D(i,j,lLsquare,aux1) 
            = MAT3D(nel,i,j,lSquare,ndm,maxViz);
/*... R*/
        if(ndm == 2)
        { 
          lLsquareR[0] = MAT2D(nel,0,lSquareR,3); 
          lLsquareR[1] = MAT2D(nel,1,lSquareR,3); 
          lLsquareR[2] = MAT2D(nel,2,lSquareR,3); 
        }
        else if(ndm == 3)
        { 
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
          MAT2D(i,j,lGradU,ndm) =  MAT3D(nel,i,j,gradU,ndf,ndm);
/*...................................................................*/
      
/*...*/      
      for(i=0;i<aux1;i++)
      {
        lDcca[i]  = MAT2D(nel,i,gDcca ,maxViz);
        lFaceR[i] = MAT2D(nel,i,faceR ,aux2);
        lFaceL[i] = MAT2D(nel,i,faceL ,aux2);
        
        for(j=0;j<ndm;j++)
          MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,fXmCc  ,maxViz,ndm);
        
        vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
        lViz[i] = vizNel;
        if (vizNel != -2)
        {
          lVolume[i] = gVolume[vizNel];
          for (j = 0; j<ndf; j++)
            MAT2D(i, j, lu, ndf) = MAT2D(vizNel, j, u, ndf);
        }
      }
/*...................................................................*/

/*... propriedades por face*/
      
      for (i = 0; i<aux1; i++)
      {
        idFace     = MAT2D(nel, i, cellFace, maxViz) - 1;
        cellOwner  = MAT2D(idFace, 0, fOwner, 2) - 1;
        ch         = OWNER(cellOwner, nel);
        lmKsi[i]   = fModKsi[idFace];
        lfArea[i]  = fArea[idFace];
        for (j = 0; j<ndm; j++)
        {
          MAT2D(i, j, lKsi, ndm)    = ch * MAT2D(idFace, j, fKsi, ndm);
          MAT2D(i, j, lEta, ndm)    = ch * MAT2D(idFace, j, fEta, ndm);
          MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
          MAT2D(i, j, lXm, ndm)     = MAT2D(idFace, j, fXm, ndm);
          MAT2D(i, j, lvSkew, ndm)  = MAT2D(idFace, j, fvSkew, ndm);
        }
      }
/*...................................................................*/
    
/*... chamando a biblioteca de celulas*/
      ty = geomType[nel];
      sn(isNod,ty,nel); 
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
      cellLibRcGrad(loads       ,rcGrad
                   ,lViz        ,lProp
                   ,lLsquare    ,lLsquareR
                   ,lKsi        ,lmKsi
                   ,lEta        ,lfArea 
                   ,lNormal     ,lVolume
                   ,lvSkew      
                   ,lXm         ,lXmcc 
                   ,lDcca       
                   ,lFaceR      ,lFaceL
                   ,lu          ,lGradU
                   ,lnU         ,ty     
                   ,nFace[nel]  ,ndm 
                   ,ndf
                   ,isNod       ,nel);    
/*...................................................................*/

/*...*/
      for(i=0;i<ndf;i++)
        for(j=0;j<ndm;j++)
          MAT3D(nel,i,j,gradU,ndf,ndm)  = MAT2D(i,j,lGradU,ndm);
/*...................................................................*/
    }
/*...................................................................*/
  }  
/*...................................................................*/

/*...*/
  if(mpiVar.nPrcs > 1 ) comunicateCel(iCel,gradU,ndm,ndf);
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 18/07/2018                                   *
 *-------------------------------------------------------------------*
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
void rcLeastSquare(INT *RESTRICT cellFace, INT *RESTRICT fOwner
                  ,DOUBLE *RESTRICT fModKsi, DOUBLE *RESTRICT fKsi
                  ,DOUBLE *RESTRICT lSquare ,DOUBLE *RESTRICT lSquareR
                  ,short *RESTRICT nFace       
                  ,INT const numel          ,short const maxViz
                  ,short const type         ,short const ndm){

  short lnFace, i, j;
  INT nEl;
  INT idFace, cellOwner, ch;
  DOUBLE lLsquare[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lLsquareR[2*MAX_NDM];
  DOUBLE lKsi[MAX_NUM_FACE*3],lmKsi[MAX_NUM_FACE];
  
  for(nEl=0;nEl<numel;nEl++)
  {
    lnFace  = nFace[nEl];

    for (i = 0; i<lnFace; i++)
    {
      idFace = MAT2D(nEl, i, cellFace, maxViz) - 1;
      cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
      ch = OWNER(cellOwner, nEl);
      lmKsi[i] = fModKsi[idFace];
      for (j = 0; j<ndm; j++)
        MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
    }

    leastSquareMatrix(lKsi    ,lmKsi
                     ,lLsquare,lLsquareR
                     ,type
                     ,lnFace  ,ndm);

    for(i=0;i<ndm;i++)
      for(j=0;j<lnFace;j++) 
        MAT3D(nEl,i,j,lSquare,ndm,maxViz) = MAT2D(i,j,lLsquare,lnFace);

    if(type == RCLSQUAREQR)
    {
      if(ndm == 2)
      {
        MAT2D(nEl,0,lSquareR,3) = lLsquareR[0];
        MAT2D(nEl,1,lSquareR,3) = lLsquareR[1];
        MAT2D(nEl,2,lSquareR,3) = lLsquareR[2];
      }
      else if(ndm == 3)
      {
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
  void convTempForKelvin(DOUBLE *RESTRICT u,INT const n
                        ,bool const fKelvin){
    int i;
/*...*/
    if(fKelvin){
      for(i=0;i<n;i++){
        u[i] = CELSIUS_FOR_KELVIN(u[i]); 
      }
    }
/*...................................................................*/ 

/*...*/
    else{
      for(i=0;i<n;i++){
        u[i] = KELVIN_FOR_CELSIUS(u[i]); 
      }
    }
/*...................................................................*/    
  }
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 18/07/2018                                   *
 *-------------------------------------------------------------------*
 * MESHQUALITY: calulo das propriedades da malha                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * nFace     -> numero de faces por celulas                          *
 * gVolume -> volumes das celulas                                    * 
 * fKsi    -> vetores que unem centroide da celula central aos       *
 *            vizinhos destas                                        *
 * fNormal -> vetores normais as faces das celulas                   *
 * fModvSkew -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face                      * 
 * gDcca   -> menor distancia do centroide a faces desta celula      *
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
               , INT *RESTRICT cellFace    , INT *RESTRICT fOwner
               , short  *RESTRICT nFace    , DOUBLE *RESTRICT volume
               , DOUBLE *RESTRICT fKsi     , DOUBLE *RESTRICT fNormal
               , DOUBLE *RESTRICT fModvSkew, DOUBLE *RESTRICT gDcca
               , short const maxViz        , short const ndm
               , INT const numel)        
{
  INT nEl,k=0;
  INT idFace;
  DOUBLE volumeTotal = 0.e0,nk,nkMin=1.e0,nkMed=0.e0;
  DOUBLE skewMax=0.0e0, skewMed=0.0e0,teta;
  DOUBLE lenth,lMax,lMin,aspectRaMax=0.0,aspectRaMin=1.e+32;
  short nf,j;

  for(nEl=0;nEl<numel;nEl++){
/*... Volume total da malha*/
    volumeTotal += volume[nEl];
/*...*/
    lMax=0.0;
    lMin=1.0e+32;
/*...*/
    for(nf=0;nf<nFace[nEl];nf++)
    {      
      idFace    = MAT2D(nEl, nf, cellFace, maxViz) - 1;
/*... k * normal*/
      for(j=0,nk= 0.e0;j<ndm;j++)
        nk += MAT2D(idFace, j, fKsi, ndm)
             *MAT2D(idFace, j, fNormal, ndm);
/*...................................................................*/ 
  
/*... nao-ortoganilidade*/   
      nkMed += nk;
      nkMin = min(nkMin,nk); 
      k++;   
      nk       = fModvSkew[idFace];
      skewMed += nk;
      skewMax  = max(skewMax,nk);
/*...................................................................*/

/*... aspcto ratio*/
      lenth = MAT2D(nEl,nf,gDcca,maxViz);
      lMax  = max(lMax,lenth);
      lMin  = min(lMin,lenth);
/*...................................................................*/
    }
/*...................................................................*/ 

/*... aspectRa*/    
    aspectRaMax = max(aspectRaMax,lMax/lMin);
    aspectRaMin = min( aspectRaMin,lMax/lMin);
/*...................................................................*/

  }
/*...................................................................*/ 

  skewMed    /= k;
  nkMed      /= k;

  teta            = acos(nkMed);
  teta            = radToDeg(teta);
  mq->nonOrthMed  = teta; 

  teta            = acos(nkMin);
  teta            = radToDeg(teta);
  mq->nonOrthMax  = teta; 
  
  mq->volume      = volumeTotal; 
  mq->skewMed     = skewMed; 
  mq->skewMax     = skewMax; 

  mq->aspectRaMax =  aspectRaMax;
  mq->aspectRaMin =  aspectRaMin;
 
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 04/07/2016                                   *
 * Data de modificaco : 24/11/2017                                   * 
 *-------------------------------------------------------------------* 
 * WALLFLUID : identifica as paredes impermeveis estaticas           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * faceR   -> condicoes de contorno das velocidades                  * 
 * nelcon  -> vizinhos dos elementos                                 * 
 * nFace   -> numero de faces por celulas                            * 
 * nEl     -> numero de toral de celulas                             * 
 * maxViz  -> numero vizinhos por celula maximo da malha             * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void wallFluid(short *RESTRICT faceR ,INT *RESTRICT nelcon
              ,short *RESTRICT nFace     
              ,INT const nEl         ,short const maxViz){
  
  INT i,j,vizNel,aux2;                   


  aux2 = maxViz+1;
  for(i=0;i<nEl;i++)
    for(j=0;j<nFace[i];j++){
      vizNel = MAT2D(i,j,nelcon,maxViz);
/*... contorno*/
      if( vizNel == -1){
/*... parede estatic*/
        if(MAT2D(i,j,faceR,aux2) == 0){
          MAT2D(i,j,faceR,aux2) = STATICWALL;
        }
/*....................................................................*/

/*... parede estatic*/
        else if(MAT2D(i,j,faceR,aux2) == SLIP){
          MAT2D(i,j,faceR,aux2) = 0;
        }
/*....................................................................*/
      }
/*....................................................................*/
    }
}   
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 28/08/2017                                   *
 * Data de modificaco : 20/08/2019                                   *
 *-------------------------------------------------------------------*
 * PARAMETERCELLLM:                                                  *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * vel       - velocidade                                            *
 * prop      - propriedade definidas por materiaç                    * 
 * density   - densidade por celula                                  *
 * sHeat     - calor especifico por celula                           *
 * tCond     - condutividade termica                                 *
 * viscosity -  viscosidade dinamica                                 *
 * volume    - volume                                                *
 * mat        - material por celula                                  *
 * cfl        - cfl                                                  *
 * reynolds- numero de reynolds                                      *
 * peclet  - peclet                                                  *
 * fPrameter - parametro a serem calculados                          *
 *              0 - cfl                                              *
 *              1 - reynolds                                         *
 *              2 - peclet                                           *
 *              3 - massa total                                      *
 * dt         - passo de tempo                                       *
 * nEl        - numero de celulas                                    *
 * ndm        - numero de dimensoes                                  *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * cfl      - numero de peclet nesse passo de tempo                  *
 * reynolds - numero de reynolds                                     *
 * peclet   - numero de peclet                                       *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
*********************************************************************/
void parameterCellLm(DOUBLE *RESTRICT vel    , DOUBLE *RESTRICT prop
                , DOUBLE *RESTRICT density   , DOUBLE *RESTRICT sHeat
                , DOUBLE *RESTRICT tCond     , DOUBLE *RESTRICT dViscosity
                , DOUBLE *RESTRICT volume    , short  *RESTRICT mat      
                , DOUBLE *cfl                , DOUBLE *reynolds
                , DOUBLE *peclet             , DOUBLE *mass   
                , bool *fParameter           , DOUBLE const dt
                , INT const nEl              , short const ndm)

{
  INT i;
  DOUBLE modVel,lc,den,sHeat0,viscosity,coefDif,tmp,v[3]
        ,dm=0.e0;
  DOUBLE cflMax=0.0e0,reynoldsMax=0.e0,pecletMax=0.e0;
#ifdef _MPI_
  DOUBLE gg;
#endif

  for(i=0;i<nEl;i++)
  {

/*... numero de Courant–Friedrichs–Lewy*/
    if(fParameter[0])
    {
/*... modulo das velocidades*/
      if(ndm == 2){
        v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
        v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
        modVel = sqrt( v[0] + v[1] );
      }
      else{
        v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
        v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
        v[2] = MAT2D(i, 2, vel, ndm)*MAT2D(i, 2, vel, ndm);
        modVel = sqrt(v[0] + v[1] + v[2] );
      }
/*..................................................................*/

/*... tamanho caracteristico*/
      lc = sizeCar(volume[i],ndm);
/*..................................................................*/

/*...*/
      cflMax = max(modVel*dt/lc,cflMax);  
/*..................................................................*/
    }
/*..................................................................*/  

/*... numero de Reynolds*/
    if(fParameter[1])
    {
/*... modulo das velocidades*/
      if (ndm == 2)
      {
        v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
        v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
        modVel = sqrt(v[0] + v[1]);
      }
      else
      {
        v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
        v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
        v[2] = MAT2D(i, 2, vel, ndm)*MAT2D(i, 2, vel, ndm);
        modVel = sqrt(v[0] + v[1] + v[2]);
      }
/*..................................................................*/

/*... tamanho caracteristico*/
      lc = sizeCar(volume[i],ndm);
/*..................................................................*/

/*...*/
      den       = MAT2D(i,2 ,density ,DENSITY_LEVEL);
      viscosity = dViscosity[i];
/*..................................................................*/


/*...*/
      tmp         = den*modVel*lc/viscosity;
      reynoldsMax = max(tmp,reynoldsMax);  
/*..................................................................*/
    }
/*..................................................................*/

/*... numero de peclet*/
      if (fParameter[2])
      {
/*... modulo das velocidades*/
        if (ndm == 2)
        {
          v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
          v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
          modVel = sqrt(v[0] + v[1]);
        }
        else
        {
          v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
          v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
          v[2] = MAT2D(i, 2, vel, ndm)*MAT2D(i, 2, vel, ndm);
          modVel = sqrt(v[0] + v[1] + v[2]);
        }
/*..................................................................*/

/*... tamanho caracteristico*/
        lc = sizeCar(volume[i], ndm);
/*..................................................................*/

/*...*/
        den     = MAT2D(i, 2, density, DENSITY_LEVEL);
        sHeat0  = MAT2D(i, 2, sHeat  , SHEAT_LEVEL);
        coefDif = tCond[i];
/*..................................................................*/

/*...*/
        tmp = den*sHeat0*modVel*lc / coefDif;
        pecletMax = max(tmp, pecletMax);
/*..................................................................*/
      }
/*..................................................................*/

/*... massa total*/
      if (fParameter[3])
      {
/*...*/
        den = MAT2D(i, 2, density, DENSITY_LEVEL);
        dm += den*volume[i];
/*..................................................................*/
      }
/*..................................................................*/

  }
/*..................................................................*/  

/*....*/
#ifdef _MPI_
  if(mpiVar.nPrcs>1)
  { 
    tm.overHeadMiscMpi = getTimeC() - tm.overHeadMiscMpi;
    if(fParameter[0])
    {      
      MPI_Allreduce(&cflMax,&gg ,1,MPI_DOUBLE,MPI_MAX,mpiVar.comm);
      cflMax = gg;
    } 
    if(fParameter[1])
    {
      MPI_Allreduce(&reynoldsMax,&gg ,1,MPI_DOUBLE,MPI_MAX,mpiVar.comm);
      reynoldsMax = gg;
    }
    if(fParameter[2])
    {
      MPI_Allreduce(&pecletMax,&gg ,1,MPI_DOUBLE,MPI_MAX,mpiVar.comm);
      pecletMax = gg;
    }
    if(fParameter[3])
    {
      MPI_Allreduce(&dm,&gg ,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
      dm = gg;
    }
    tm.overHeadMiscMpi = getTimeC() - tm.overHeadMiscMpi;
  }
#endif
/*...................................................................*/

/*...*/
  *cfl      = cflMax;
  *reynolds = reynoldsMax;
  *peclet   = pecletMax;
  *mass     = dm;
/*...................................................................*/
}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 26/08/2017                                   *
 *-------------------------------------------------------------------*
 * PARAMETERCELL:                                                    *
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
void parameterCell(DOUBLE *RESTRICT vel     , DOUBLE *RESTRICT prop
                 , DOUBLE *RESTRICT density , DOUBLE *RESTRICT volume
                 , short  *RESTRICT mat       
                 , DOUBLE *cfl              , DOUBLE *reynolds
                 , bool *fParameter         , DOUBLE const dt
                 , INT const nEl            , short const ndm)

{

  DOUBLE modVel, lc, den, viscosity, tmp, v[3];
  DOUBLE cflMax = 0.0e0, reynoldsMax = 0.e0;
  INT i;
  short lMat;

/*...*/
  for (i = 0; i<nEl; i++)
  {
 /*... numero de Courant–Friedrichs–Lewy*/
    if (fParameter[0])
    {
 /*... modulo das velocidades*/
      if (ndm == 2)
      {
        v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
        v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
        modVel = sqrt(v[0] + v[1]);
      }
      else
      {
        v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
        v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
        v[2] = MAT2D(i, 2, vel, ndm)*MAT2D(i, 2, vel, ndm);
        modVel = sqrt(v[0] + v[1] + v[2]);
      }
 /*..................................................................*/

/*... tamanho caracteristico*/
      lc = sizeCar(volume[i], ndm);
/*..................................................................*/

/*...*/
      cflMax = max(modVel*dt / lc, cflMax);
/*..................................................................*/
    }
/*..................................................................*/

/*... numero de Reynolds*/
    if (fParameter[1]) 
    {
/*... modulo das velocidades*/
      if (ndm == 2)
      {
        v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
        v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
        modVel = sqrt(v[0] + v[1]);
      }
      else
      {
        v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
        v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
        v[2] = MAT2D(i, 2, vel, ndm)*MAT2D(i, 2, vel, ndm);
        modVel = sqrt(v[0] + v[1] + v[2]);
      }
/*..................................................................*/

/*... tamanho caracteristico*/
      lc = sizeCar(volume[i], ndm);
/*..................................................................*/

/*...*/
      lMat = mat[i] - 1;
      den = MAT2D(i, 0, density, DENSITY_LEVEL);
      viscosity = MAT2D(lMat, DYNAMICVISCOSITY, prop, MAXPROP);
/*..................................................................*/


/*...*/
      tmp = den*modVel*lc / viscosity;
      reynoldsMax = max(tmp, reynoldsMax);
/*..................................................................*/
    }
/*..................................................................*/

 
  }
/*..................................................................*/

  *cfl = cflMax;
  *reynolds = reynoldsMax;
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 15/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * OPENDOMAIN: verifica se o dominio e aberto ou fechado             *
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
bool openDomain(Loads *loadVel
              , short  *RESTRICT faceVelLoad, short  *RESTRICT nFace
              , INT const numel             , short const maxViz  ) {

  bool fOpen = false;
  short  nCarg, j, type, aux1, aux2;
  INT nel;

  aux2 = maxViz + 1;
  for (nel = 0; nel<numel; nel++) {
    aux1 = nFace[nel];
/*... elementos com equacoes*/
    for(j=0;j<aux1;j++){
      nCarg = MAT2D(nel, j, faceVelLoad, aux2);
      if(nCarg){
        nCarg--;
        type  = loadVel[nCarg].type;
        if (  type == INLET || type == OUTLET
            || type == OPEN || type == INLETSTAICTPRES
            || type == INLETTOTALPRES) {
          fOpen = true;
          break; 
        }     
      }
    }
  }

  return fOpen;

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 01/10/2017                                   *
 * Data de modificaco : 23/05/2019                                   *
 *-------------------------------------------------------------------*
 * MASSFLUXOPENDOMAIN: calcula da massa entrando e saindo            *
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
void massFluxOpenDomain(Loads *loadVel    , Temporal const ddt
              , INT *RESTRICT cellFace      , INT *RESTRICT fOwner
              , short  *RESTRICT faceVelLoad, short  *RESTRICT nFace
              , DOUBLE *RESTRICT fArea      , DOUBLE *RESTRICT fNormal
              , DOUBLE *RESTRICT fXm    
              , DOUBLE *RESTRICT density    , DOUBLE *RESTRICT vel
              , DOUBLE *massInOut           , DOUBLE *deltaMass
              , INT const numel             , short const ndm  
              , short const maxViz  )
{

  short  nCarg, j, type, aux1, aux2;
  INT nel;
  INT    idFace, cellOwner, ch;
  DOUBLE lDensity,aFace,mOut,mIn,n[3],v[3],wfn,dt;
  DOUBLE par[MAXLOADPARAMETER],xx[4];
#ifdef _MPI_
  DOUBLE gg;
#endif
  mOut = mIn = 0.e0;
  aux2 = maxViz + 1;
  for (nel = 0; nel<numel; nel++) {
    aux1 = nFace[nel];

/*... elementos com equacoes*/
    for(j=0;j<aux1;j++)
    {

      idFace    = MAT2D(nel, j, cellFace, maxViz) - 1;
      cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
      ch        = OWNER(cellOwner, nel);

      nCarg = MAT2D(nel, j, faceVelLoad, aux2);
      if(nCarg){
        nCarg--;
        type  = loadVel[nCarg].type;
/*...*/
        if (type == INLET) {
/*...*/
          xx[0] = MAT2D(idFace, j, fXm, ndm);
          xx[1] = MAT2D(idFace, j, fXm, ndm);
          xx[2] = 0.0e0;                             
          if(ndm == 3) xx[2] = MAT2D(idFace, j, fXm, ndm);
          xx[3] = ddt.t;
/*...................................................................*/
          getLoads(par,&loadVel[nCarg],xx); 
          aFace = fArea[idFace];
          n[0]  = ch * MAT2D(idFace, 0, fNormal, ndm);
          n[1]  = ch * MAT2D(idFace, 1, fNormal, ndm);
          v[0]  = par[0];
          v[1]  = par[1];
          wfn   = v[0]*n[0] + v[1]*n[1];
          if ( ndm == 3 ) {
            n[2] = ch * MAT2D(idFace, 2, fNormal, ndm);
            v[2] = par[2];
            wfn +=  v[2]*n[2];
          }
          lDensity = par[ndm];
      
          mIn     += lDensity*wfn*aFace;          
        }
/*...................................................................*/

/*...*/
        else if ( type == OUTLET ){
          lDensity  = MAT2D(nel, 2, density ,DENSITY_LEVEL);
          aFace     = fArea[idFace];
          n[0]      = ch * MAT2D(idFace, 0, fNormal, ndm);
          n[1]      = ch * MAT2D(idFace, 1, fNormal, ndm);
          v[0]      = MAT2D(nel, 0, vel, ndm);
          v[1]      = MAT2D(nel, 1, vel, ndm);
          wfn       = v[0]*n[0] + v[1]*n[1];
          if ( ndm == 3 ) {
            n[2] = ch * MAT2D(idFace, 2, fNormal, ndm);
            v[2] = MAT2D(nel, 2, vel, ndm);
            wfn +=  v[2]*n[2];
          }
          mOut       += lDensity*wfn*aFace;  
        }
/*...................................................................*/

/*...*/
        else if( type == OPEN || type == INLETSTAICTPRES || type == INLETTOTALPRES) {
          lDensity  = MAT2D(nel, 2, density ,DENSITY_LEVEL);
          aFace     = fArea[idFace];
          n[0]      = ch * MAT2D(idFace, 0, fNormal, ndm);
          n[1]      = ch * MAT2D(idFace, 1, fNormal, ndm);
          v[0]      = MAT2D(nel, 0, vel, ndm);
          v[1]      = MAT2D(nel, 1, vel, ndm);
          wfn = v[0]*n[0] + v[1]*n[1];
          if (ndm == 3) {
            n[2] = ch * MAT2D(idFace, 2, fNormal, ndm);
            v[2] = MAT2D(nel, 2, vel, ndm);
            wfn +=  v[2]*n[2];
          } 
          if(wfn > 0.e0)
            mOut += lDensity*wfn*aFace;
          else{
            lDensity = loadVel[nCarg].par[0];
            mIn   += lDensity*wfn*aFace;  
          }  
/*...................................................................*/          
        }     
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*....*/
#ifdef _MPI_
  if(mpiVar.nPrcs>1)
  { 
    tm.overHeadMiscMpi = getTimeC() - tm.overHeadMiscMpi;
    MPI_Allreduce(&mIn,&gg ,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    mIn = gg;
    MPI_Allreduce(&mOut,&gg ,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    mOut = gg;
    tm.overHeadMiscMpi = getTimeC() - tm.overHeadMiscMpi;
  }
#endif
/*...................................................................*/

/*...*/
  massInOut[0] = mIn;
  massInOut[1] = mOut;
  dt = ddt.dt[1];
  *deltaMass = -(mOut + mIn) * dt;
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 01/10/2017                                   *
 * Data de modificaco : 16/08/2019                                   *
 *-------------------------------------------------------------------*
 * TOTALMASS ; calculo da massa total do sistema                     *
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
DOUBLE totalMass(DOUBLE *RESTRICT density  , DOUBLE *RESTRICT volume
                ,INT const nEl){

  INT i;
  DOUBLE dm, den, gDm;

  dm = 0.e0;
  for(i=0;i<nEl;i++)
  {
    den = MAT2D(i, 2, density, DENSITY_LEVEL);
    dm += den*volume[i];
   }
/*..................................................................*/  

/*....*/
#ifdef _MPI_
  if(mpiVar.nPrcs>1)
  { 
    tm.overHeadMiscMpi = getTimeC() - tm.overHeadMiscMpi;
    MPI_Allreduce(&dm,&gDm,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.overHeadMiscMpi = getTimeC() - tm.overHeadMiscMpi;
  }
  else  
    gDm = dm;
#else
    gDm = dm;
#endif
/*...................................................................*/

  return gDm;

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/02/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * hPres : pressao com parcela hidroestatica                         *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * pres0   -> pressao do tempo (n-1)                                 *
 * pres    -> pressao do tempo (n)                                   *
 * dFluid  -> densidade do Fluid                                     *
 * gravity -> vetor campo de gravidade                               *
 * xRef    -> pontos de referencia                                   *
 * nEl     -> numero de elementos                                    *
 * ndm     -> numero de dimensoes                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * pres0   -> p= p - rho*g*h                                         *
 * pres    -> p= p - rho*g*h                                         *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
*********************************************************************/
void hPres(DOUBLE *RESTRICT pres0  , DOUBLE *RESTRICT pres
         , DOUBLE *RESTRICT dFluid , DOUBLE *RESTRICT cc
         , DOUBLE *RESTRICT gravity, DOUBLE *RESTRICT xRef
         , INT const nEl           , short const ndm) 
{
  short j,nD=DENSITY_LEVEL;
  INT i;
  DOUBLE xc[3],tmp,dm;
  
  xRef[0] =  0.0e0;
  xRef[1] =  0.0e0;
  xRef[2] =  0.0e0;

  for(i=0; i < nEl; i++){
    
    for(tmp = 0.e0, j=0; j < ndm; j++){
      xc[j] = MAT2D(i,j,cc,ndm);
      tmp  += (xc[j]-xRef[j])*gravity[j];
    }
    dm       = MAT2D(i,0,dFluid,nD);
    pres0[i] += dm*tmp;
    pres[i]   = pres0[i];
  }  
}
/*********************************************************************/

void systFormDifOld(Loads *loads, Diffusion *diff
                  , INT    *RESTRICT el, INT    *RESTRICT nelcon
                  , short  *RESTRICT nen, short  *RESTRICT nFace
                  , DOUBLE *RESTRICT gVolume, DOUBLE *RESTRICT gDcca
                  , DOUBLE *RESTRICT gXmCc
                  , DOUBLE *RESTRICT gModksi, DOUBLE *RESTRICT gKsi
                  , DOUBLE *RESTRICT gEta, DOUBLE *RESTRICT gfArea
                  , DOUBLE *RESTRICT gNormal, DOUBLE *RESTRICT gXm
                  , DOUBLE *RESTRICT gModvSkew, DOUBLE *RESTRICT gvSkew
                  , short  *RESTRICT geomType, DOUBLE *RESTRICT prop
                  , short  *RESTRICT calType, short  *RESTRICT mat
                  , DOUBLE *RESTRICT density
                  , INT    *RESTRICT ia, INT    *RESTRICT ja
                  , DOUBLE *RESTRICT a, DOUBLE *RESTRICT ad
                  , DOUBLE *RESTRICT b, INT    *RESTRICT id
                  , short  *RESTRICT faceR, short  *RESTRICT faceL
                  , DOUBLE *RESTRICT u0, DOUBLE *RESTRICT gradU0
                  , DOUBLE *RESTRICT rCell, Temporal *ddt
                  , INT nEq, INT nEqNov
                  , INT nAd, INT nAdR
                  , short maxNo, short maxViz
                  , short ndm, INT numel
                  , short ndf, short storage
                  , bool forces, bool matrix
                  , bool calRcell, bool unsym)
{
  short i, j;
  INT nel, vizNel;
  /*... variavel local */
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM], lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM], lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM], lVolume[MAX_NUM_FACE + 1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM], lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE], lvSkew[MAX_NUM_FACE*MAX_NDM];
  short  lGeomType[MAX_NUM_FACE + 1];
  short  lib;
  short  lFaceR[MAX_NUM_FACE + 1];
  short  lFaceL[MAX_NUM_FACE + 1];
  DOUBLE lDensity;
  DOUBLE lA[(MAX_NUM_FACE + 1)*MAX_NDF], lB[MAX_NDF];
  DOUBLE lProp[(MAX_NUM_FACE + 1)*MAXPROP];
  DOUBLE lu0[(MAX_NUM_FACE + 1)*MAX_NDF];
  DOUBLE lGradU0[(MAX_NUM_FACE + 1)*MAX_NDM];
  DOUBLE lRcell[MAX_NDF];
  INT    lId[(MAX_NUM_FACE + 1)*MAX_NDF], lViz[MAX_NUM_FACE];
  short  aux1, aux2, lMat;

/*... loop nas celulas*/
  aux2 = maxViz + 1;
  for (nel = 0; nel<numel; nel++) {
/*...*/
    if (calRcell)
      for (j = 0; j<ndf; j++)
        MAT2D(nel, j, rCell, ndf) = 0.e0;;
    /*...*/
    aux1 = nFace[nel];
/*... elementos com equacoes*/
    if (MAT2D(nel, aux1, faceR, aux2) != PCCELL) {

/*... zerando vetores*/
      for (j = 0; j<(MAX_NUM_FACE + 1)*MAX_NDF; j++) {
        lId[j] = -1;
        lu0[j] = 0.e0;
      }

/*... loop na celula central*/
      lMat = mat[nel] - 1;
      lib = calType[lMat];
      lVolume[aux1] = gVolume[nel];
      lGeomType[aux1] = geomType[nel];
      lFaceR[aux1] = MAT2D(nel, aux1, faceR, aux2);
      lFaceL[aux1] = MAT2D(nel, aux1, faceL, aux2);
      lDensity = MAT2D(nel, 0, density, DENSITY_LEVEL);

      for (j = 0; j<ndf; j++) {
        MAT2D(aux1, j, lu0, ndf) = MAT2D(nel, j, u0, ndf);
        MAT2D(aux1, j, lId, ndf) = MAT2D(nel, j, id, ndf) - 1;
      }

      for (j = 0; j<MAXPROP; j++)
        MAT2D(aux1, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);

      for (j = 0; j<ndm; j++)
        MAT2D(aux1, j, lGradU0, ndm) = MAT2D(nel, j, gradU0, ndm);

      for (i = 0; i<aux1; i++)
      {
        lDcca[i] = MAT2D(nel, i, gDcca, maxViz);
        aux2 = (maxViz + 1);
        lFaceR[i] = MAT2D(nel, i, faceR, aux2);
        lFaceL[i] = MAT2D(nel, i, faceL, aux2);
        for (j = 0; j<ndm; j++)
        {
          MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);
        }
      }

      for (i = 0; i<aux1; i++) 
      {
        lmKsi[i]   = MAT2D(nel, i, gModksi, maxViz);
        lfArea[i]  = MAT2D(nel, i, gfArea, maxViz);
        lDcca[i]   = MAT2D(nel, i, gDcca, maxViz);
        lmvSkew[i] = MAT2D(nel, i, gModvSkew, maxViz);
        aux2 = (maxViz + 1);
        lFaceR[i] = MAT2D(nel, i, faceR, aux2);
        lFaceL[i] = MAT2D(nel, i, faceL, aux2);
        for (j = 0; j<ndm; j++) 
        {
          MAT2D(i, j, lKsi, ndm) = MAT3D(nel, i, j, gKsi, maxViz, ndm);
          MAT2D(i, j, lEta, ndm) = MAT3D(nel, i, j, gEta, maxViz, ndm);
          MAT2D(i, j, lNormal, ndm) = MAT3D(nel, i, j, gNormal, maxViz, ndm);
          MAT2D(i, j, lXm, ndm) = MAT3D(nel, i, j, gXm, maxViz, ndm);
          MAT2D(i, j, lXmcc, ndm) = MAT3D(nel, i, j, gXmCc, maxViz, ndm);
          MAT2D(i, j, lvSkew, ndm) = MAT3D(nel, i, j, gvSkew, maxViz, ndm);
        }
      }

/*... loop na celulas vizinhas*/
      for (i = 0; i<aux1; i++)
      {
        vizNel = MAT2D(nel, i, nelcon, maxViz) - 1;
        lViz[i] = vizNel;
        if (vizNel != -2)
        {
          lVolume[i] = gVolume[vizNel];
          lGeomType[i] = geomType[vizNel];
          lMat = mat[vizNel] - 1;
          for (j = 0; j<ndf; j++) {
            MAT2D(i, j, lu0, ndf) = MAT2D(vizNel, j, u0, ndf);
            MAT2D(i, j, lId, ndf) = MAT2D(vizNel, j, id, ndf) - 1;
          }
          for (j = 0; j<ndm; j++)
            MAT2D(i, j, lGradU0, ndm) = MAT2D(vizNel, j, gradU0, ndm);
          for (j = 0; j<DIFPROP; j++)
            MAT2D(i, j, lProp, MAXPROP) = MAT2D(lMat, j, prop, MAXPROP);
        }
      }
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
/*    cellLibDif(loads    , diff
               , lGeomType, lProp
               , lViz     , lId
               , lKsi     , lmKsi
               , lEta     , lfArea
               , lNormal  , lVolume
               , lXm      , lXmcc
               , lDcca    , &lDensity
               , lvSkew   , lmvSkew
               , lA       , lB
               , lRcell   , ddt
               , lFaceR   , lFaceL
               , lu0      , lGradU0
               , nen[nel] , nFace[nel]
               , ndm      , lib
               , nel);*/
/*...................................................................*/

/*... residuo da celula*/
      if (calRcell)
        for (j = 0; j<ndf; j++)
          MAT2D(nel, j, rCell, ndf) = lRcell[j];
/*...................................................................*/

/*...*/
      assbly(ia, ja
        , a, ad
        , b
        , lId
        , lA, lB
        , nEq, nEqNov
        , nAd, nAdR
        , nFace[nel], ndf
        , storage, forces
        , matrix, unsym);
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/
