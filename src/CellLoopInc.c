#include<CellLoop.h>

/********************************************************************* 
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 06/10/2019                                   * 
 *-------------------------------------------------------------------* 
 * SYSTFOMSIMPLEVEL: calculo do sistema de equacoes para problemas   * 
 * de escomaneto de fluidos ( Vel )                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loadsVel  -> definicoes de cargas de velocidades                  * 
 * loadsPres -> definicoes de cargas de pressao                      * 
 * advVel    -> tecnica da discretizacao do termo advecao            *
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * tModel    -> modelo de turbulencia                                *
 * momentumModel -> termos/modelos da equacao de momento linear      *
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
 * prop    -> propriedades dos material                              *
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
 * facePresR -> restricoes por elemento de pressao                   * 
 * pres      -> campo de pressao conhecido                           * 
 * grdPres   -> campo de gradiente de pressao conhecido              * 
 * gradVel   -> gradiente da solucao conhecido                       * 
 * dField    -> matriz D do metodo simple                            * 
 * underU    -> fator underrelaxtion sinple                          * 
 * vel       -> campo de velocidade conhecido                        * 
 * rCell     -> nao definido                                         * 
 * stressR   -> tensao do modelo turbulento estrutural               *
 * eddyVisc  -> viscosidade turbulneta                               *
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
 * ntn       -> numero de termo do tensor de tensao                  *
 * storage   -> tecnica de armazenamento da matriz esparsa           * 
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
 *********************************************************************/
  void systFormSimpleVel(Loads *loadsVel   , Loads *loadsPres
               , Advection *advVel         , Diffusion *diffVel
               , Turbulence *tModel        , MomentumModel *momentumModel    
               , Interface *iCel           , short typeSimple          
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
               , INT    *RESTRICT ia       , INT    *RESTRICT ja
               , DOUBLE *RESTRICT a        , DOUBLE *RESTRICT ad
               , DOUBLE *RESTRICT b        , INT    *RESTRICT id
               , short  *RESTRICT faceVelR , short  *RESTRICT facePresR
               , DOUBLE *RESTRICT pres     , DOUBLE *RESTRICT gradPres
               , DOUBLE *RESTRICT vel      , DOUBLE *RESTRICT gradVel
               , DOUBLE *RESTRICT dField   , DOUBLE underU
               , DOUBLE *RESTRICT rCell    , DOUBLE *RESTRICT stressR
               , DOUBLE *RESTRICT eddyVisc , DOUBLE *RESTRICT wallPar 
               , Temporal *ddt
               , INT nEq                   , INT nEqNov
               , INT nAd                   , INT nAdR
               , short maxNo               , short maxViz
               , short ndm                 , INT numel
               , short ndf                 , short ntn    
               , short storage
               , bool forces               , bool matrix
               , bool calRcell             , bool unsym)
{
  bool fTurb       = tModel->fTurb
      ,fTurbStruct = tModel->fTurbStruct
      ,fWallModel  = tModel->fWall;
  short i,j,k;
  short nThreads = ompVar.nThreadsCell;
  INT nel,vizNel;

/*... variavel local */
  short  aux1,aux2,lMat;
  short  lGeomType[MAX_NUM_FACE+1];
  short  lib;
  short  lFaceVelR[MAX_NUM_FACE+1],lFacePresR[MAX_NUM_FACE+1];
  INT    idFace, cellOwner, ch;
  INT    lId[(MAX_NUM_FACE+1)*MAX_NDF],lViz[MAX_NUM_FACE];
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM],lmKsi[MAX_NUM_FACE];
  DOUBLE lEta[MAX_NUM_FACE*MAX_NDM],lfArea[MAX_NUM_FACE];
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM],lVolume[MAX_NUM_FACE+1];
  DOUBLE lXm[MAX_NUM_FACE*MAX_NDM],lXmcc[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lDcca[MAX_NUM_FACE];
  DOUBLE lmvSkew[MAX_NUM_FACE],lvSkew[MAX_NUM_FACE*MAX_NDM];
  DOUBLE lA[MAX_NUM_FACE+ MAX_NDM],lB[MAX_NDF];
  DOUBLE lProp[(MAX_NUM_FACE+1)*MAXPROP];
  DOUBLE lPres[(MAX_NUM_FACE+1)];
  DOUBLE lGradVel[(MAX_NUM_FACE+1)*MAX_NDM*MAX_NDF];
  DOUBLE lGradPres[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lVel[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lCc[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lRcell[MAX_NDF],lDfield[(MAX_NUM_FACE+1)*MAX_NDM];
  DOUBLE lStressR[(MAX_NUM_FACE+1)*6],lEddyVisc[MAX_NUM_FACE+1];
  DOUBLE lWallPar[NWALLPAR];

/*...*/
  if(ompVar.fCell)
  {
/*... loop nas celulas*/
    aux2 = maxViz + 1;
  #pragma omp parallel  for default(none) num_threads(nThreads)\
     private(nel,i,j,k,aux1,lId,lPres,lMat,lib,lVolume,lGeomType\
          ,lA,lB,lDfield,lProp,lEddyVisc\
          ,lFaceVelR,lFacePresR,lGradPres\
          ,lVel,lCc,lGradVel,lmKsi,lfArea,lDcca,lmvSkew,lKsi,lEta\
          ,lNormal,lXm,lXmcc,lvSkew,vizNel,lViz,lRcell\
          ,lStressR,lWallPar,idFace,cellOwner,ch)\
     shared(aux2,ndm,ndf,numel,maxViz,calRcell,rCell,nFace,mat,ntn\
         ,calType,gVolume, geomType, faceVelR, facePresR\
         ,gradPres,vel,gCc,gradVel,fModKsi,prop\
         ,fArea,gDcca,fModvSkew,fKsi,fEta,fNormal,fXm,gXmCc,fvSkew\
         ,nelcon,id,loadsVel,loadsPres,advVel,diffVel,typeSimple\
         ,ddt,underU,nen,ia,ja,a,ad,b,nEq,nEqNov,nAd\
         ,nAdR,storage,forces,matrix,unsym,pres,dField,eddyVisc\
         ,stressR,tModel,momentumModel,wallPar\
         ,fOwner,cellFace,fTurb,fTurbStruct,fWallModel)
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
      if(MAT2D(nel,aux1,faceVelR ,aux2) < 1)
      {

/*... loop na celula central*/    
        lMat              = mat[nel]-1;
        lib               = calType[lMat];
        lVolume[aux1]     = gVolume[nel]; 
        lGeomType[aux1]   = geomType[nel];
        lPres[aux1]       = pres[nel];
        lFaceVelR[aux1]   = MAT2D(nel,aux1,faceVelR ,aux2);
        lFacePresR[aux1]  = MAT2D(nel,aux1,facePresR ,aux2);
/*... viscosidade dinamica e turbulentea*/
        if(fTurb)
          lEddyVisc[aux1] = eddyVisc[nel];
/*...................................................................*/

/*...*/      
        for(j=0;j<MAXPROP;j++)
          MAT2D(aux1,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
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
        for(i=0;i<aux1;i++)
        {
          lDcca[i]      = MAT2D(nel,i,gDcca   ,maxViz);
          lFaceVelR[i]  = MAT2D(nel,i,faceVelR ,aux2);
          lFacePresR[i] = MAT2D(nel,i,facePresR ,aux2);
/*... propriedades por face*/
          idFace     = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner  = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch         = OWNER(cellOwner, nel);
          lmKsi[i]   = fModKsi[idFace];
          lfArea[i]  = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for(j=0;j<ndm;j++)
          {
            MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,gXmCc  ,maxViz,ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
/*...................................................................*/

/*... loop na celulas vizinhas*/  
          vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
          lViz[i] = vizNel;
          lId[i]  = - 1;
          if( vizNel != -2) 
          {
            lVolume[i]   = gVolume[vizNel]; 
            lGeomType[i] = geomType[vizNel];
            lMat         = mat[vizNel]-1;
            lPres[i]     = pres[vizNel];
/*... viscosidade dinamica e turbulentea*/
            if(fTurb)
              lEddyVisc[i] = eddyVisc[vizNel];
/*...................................................................*/

/*...*/      
            for(j=0;j<MAXPROP;j++)
              MAT2D(i,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
/*...................................................................*/

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
        cellLibSimpleVel(loadsVel ,loadsPres   
                      ,advVel     ,diffVel
                      ,tModel     ,momentumModel     
                      ,typeSimple
                      ,lGeomType  ,lProp 
                      ,lViz       ,lId           
                      ,lKsi       ,lmKsi
                      ,lEta       ,lfArea 
                      ,lNormal    ,lVolume
                      ,lXm        ,lXmcc
                      ,lDcca      ,lCc
                      ,lvSkew     ,lmvSkew
                      ,lA         ,lB
                      ,lRcell     ,ddt
                      ,lFaceVelR  ,lFacePresR           
                      ,lPres      ,lGradPres    
                      ,lVel       ,lGradVel
                      ,lDfield    ,lStressR 
                      ,lEddyVisc  ,lWallPar
                      ,underU     
                      ,nen[nel]   ,nFace[nel] 
                      ,ndm        ,lib   
                      ,nel);    
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
   	      MAT2D(nel,j,dField,ndm) = lDfield[j];
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
      if(MAT2D(nel,aux1,faceVelR ,aux2) < 1)
      {

/*... loop na celula central*/    
        lMat              = mat[nel]-1;
        lib               = calType[lMat];
        lVolume[aux1]     = gVolume[nel]; 
        lGeomType[aux1]   = geomType[nel];
        lPres[aux1]       = pres[nel];
        lFaceVelR[aux1]   = MAT2D(nel,aux1,faceVelR ,aux2);
        lFacePresR[aux1]  = MAT2D(nel,aux1,facePresR ,aux2);     
/*... viscosidade dinamica e turbulentea*/
        if(fTurb)
          lEddyVisc[aux1] = eddyVisc[nel];
/*...................................................................*/

/*...*/      
        for(j=0;j<MAXPROP;j++)
          MAT2D(aux1,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
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
        for(i=0;i<aux1;i++)
        {
          lDcca[i]      = MAT2D(nel,i,gDcca   ,maxViz);
          lFaceVelR[i]  = MAT2D(nel,i,faceVelR ,aux2);
          lFacePresR[i] = MAT2D(nel,i,facePresR ,aux2);
/*... propriedades por face*/
          idFace     = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner  = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch         = OWNER(cellOwner, nel);
          lmKsi[i]   = fModKsi[idFace];
          lfArea[i]  = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];
          for(j=0;j<ndm;j++)
          {
            MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,gXmCc  ,maxViz,ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
/*...................................................................*/

/*... loop na celulas vizinhas*/  
          vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
          lViz[i] = vizNel;
          lId[i]  = - 1;
          if( vizNel != -2) 
          {
            lVolume[i]   = gVolume[vizNel]; 
            lGeomType[i] = geomType[vizNel];
            lMat         = mat[vizNel]-1;
            lPres[i]     = pres[vizNel];
/*... viscosidade dinamica e turbulentea*/
            if(fTurb)
              lEddyVisc[i] = eddyVisc[vizNel];
/*...................................................................*/

/*...*/      
            for(j=0;j<MAXPROP;j++)
              MAT2D(i,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
/*...................................................................*/

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
        cellLibSimpleVel(loadsVel ,loadsPres   
                      ,advVel     ,diffVel
                      ,tModel     ,momentumModel     
                      ,typeSimple
                      ,lGeomType  ,lProp 
                      ,lViz       ,lId           
                      ,lKsi       ,lmKsi
                      ,lEta       ,lfArea 
                      ,lNormal    ,lVolume
                      ,lXm        ,lXmcc
                      ,lDcca      ,lCc
                      ,lvSkew     ,lmvSkew
                      ,lA         ,lB
                      ,lRcell     ,ddt
                      ,lFaceVelR  ,lFacePresR                                 
                      ,lPres      ,lGradPres    
                      ,lVel       ,lGradVel
                      ,lDfield    ,lStressR 
                      ,lEddyVisc  ,lWallPar
                      ,underU     
                      ,nen[nel]   ,nFace[nel] 
                      ,ndm        ,lib   
                      ,nel);    
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
   	      MAT2D(nel,j,dField,ndm) = lDfield[j];
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
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 01/07/2016                                   *
 * Data de modificaco : 06/10/2019                                   * 
 *-------------------------------------------------------------------* 
 * SYSTFOMSIMPLEPRES:calculo do sistema de equacoes para problemas   * 
 * de escomaneto de fluidos (Pres)                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loadsVel  -> definicoes de cargas de velocidades                  * 
 * loadsPres -> definicoes de cargas de pressao                      * 
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
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
 * prop    -> propriedades dos material                              *
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
 * facePresR -> restricoes por elemento de pressao                   * 
 * pres    -> campo de pressao conhecido                             * 
 * gradVel -> gradiente da solucao conhecido                         * 
 * dField  -> matriz D do metodo simple ( volume/A(i,i) )            * 
 * wallPar -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)    *
 * vel     -> campo de velocidade conhecido                          * 
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
void systFormSimplePres(Loads *loadsVel    , Loads *loadsPres 
							 , Diffusion *diffPres       ,Turbulence *tModel      
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
               , INT    *RESTRICT ia       , INT    *RESTRICT ja
               , DOUBLE *RESTRICT a        , DOUBLE *RESTRICT ad 
               , DOUBLE *RESTRICT b        , INT    *RESTRICT id
               , short  *RESTRICT faceVelR , short  *RESTRICT facePresR      
               , DOUBLE *RESTRICT pres     , DOUBLE *RESTRICT gradPres
               , DOUBLE *RESTRICT vel      , DOUBLE *RESTRICT dField    
               , DOUBLE *RESTRICT wallPar  , DOUBLE *RESTRICT rCell
               , Temporal *ddt              
               , INT nEq                   , INT  nEqNov
               , INT nAd                   , INT nAdR                  
               , short maxNo               , short  maxViz
               , short ndm                 , INT  numel
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
  DOUBLE lA[(MAX_NUM_FACE+1)*MAX_NDF],lB[MAX_NDF];
  DOUBLE lProp[(MAX_NUM_FACE+1)*MAXPROP];
  DOUBLE lPres[(MAX_NUM_FACE+1)];
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
          ,lA,lB,lDfield,lProp\
          ,lFaceVelR,lFacePresR,lGradPres\
          ,lVel,lmKsi,lfArea,lDcca,lmvSkew,lKsi,lEta\
          ,lNormal,lXm,lXmcc,lvSkew,vizNel,lViz,lRcell,lWallPar\
          ,idFace,cellOwner,ch)\
     shared(aux2,ndm,ndf,numel,maxViz,calRcell,rCell,nFace,mat\
         ,calType,gVolume,prop, geomType, faceVelR, facePresR\
         ,gradPres,vel,fModKsi\
         ,fArea,gDcca,fModvSkew,fKsi,fEta,fNormal,fXm,gXmCc,fvSkew\
         ,nelcon,id,loadsVel,loadsPres,diffPres\
         ,ddt,nen,ia,ja,a,ad,b,nEq,nEqNov,nAd\
         ,nAdR,storage,forces,matrix,unsym,pres,dField,wallPar\
         ,fOwner,cellFace,fWallModel)

    for(nel=0;nel<numel;nel++)
    {
/*...*/
      if(calRcell)
        rCell[nel] = 0.e0;;
/*...*/
      aux1    = nFace[nel];
/*... elementos com equacoes*/
      if(MAT2D(nel,aux1,facePresR ,aux2) < 1)
      {

/*... loop na celula central*/    
        lMat            = mat[nel]-1;
        lib             = calType[lMat];
        lVolume[aux1]   = gVolume[nel]; 
        lGeomType[aux1] = geomType[nel];
        lFaceVelR[aux1] = MAT2D(nel,aux1,faceVelR ,aux2);
        lFacePresR[aux1]= MAT2D(nel,aux1,facePresR ,aux2);
 
        lPres[aux1]     = pres[nel];
        lId[aux1]       = id[nel] - 1;
/*...................................................................*/
      
/*...*/
        for(j=0;j<MAXPROP;j++)
          MAT2D(aux1,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
/*...................................................................*/
      
/*...*/
        lDfield[aux1] = dField[nel];
        for(j=0;j<ndm;j++){
          MAT2D(aux1,j,lGradPres,ndm) = MAT2D(nel,j,gradPres,ndm);
          MAT2D(aux1,j,lVel     ,ndm) = MAT2D(nel,j,vel     ,ndm);
          MAT2D(aux1,j,lDfield  ,ndm) = MAT2D(nel,j,dField  ,ndm);
        }
/*...................................................................*/

/*... propriedades por face*/
        for(i=0;i<aux1;i++)
        {
          lDcca[i]      = MAT2D(nel,i,gDcca   ,maxViz);
          lFaceVelR[i]  = MAT2D(nel,i,faceVelR ,aux2);
          lFacePresR[i] = MAT2D(nel,i,facePresR ,aux2);
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];

          for(j=0;j<ndm;j++)
          {
            MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,gXmCc  ,maxViz,ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
/*...................................................................*/

/*... loop na celulas vizinhas*/  
          vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
          lViz[i] = vizNel;
          lId[i]  = - 1;
          if( vizNel != -2) 
          {
            lVolume[i]    = gVolume[vizNel]; 
            lGeomType[i]  = geomType[vizNel];
            lMat          = mat[vizNel]-1;
            lPres[i]      = pres[vizNel];
            lId[i]        = id[vizNel] - 1;

            for(j=0;j<ndm;j++)
            {
              MAT2D(i,j,lGradPres,ndm) = MAT2D(vizNel,j,gradPres,ndm);
              MAT2D(i,j,lVel     ,ndm) = MAT2D(vizNel,j,vel     ,ndm);
              MAT2D(i,j,lDfield  ,ndm) = MAT2D(vizNel,j,dField  ,ndm);
            }

            for(j=0;j<MAXPROP;j++)
              MAT2D(i,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
          }
        }  
/*...................................................................*/

/*...*/
        if(fWallModel)
          for (i = 0; i<NWALLPAR; i++)
            lWallPar[i] = MAT2D(nel, i, wallPar, NWALLPAR);
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
          cellLibSimplePres(loadsVel,loadsPres 
											,diffPres   
                      ,lGeomType  ,lProp 
                      ,lViz       ,lId           
                      ,lKsi       ,lmKsi
                      ,lEta       ,lfArea 
                      ,lNormal    ,lVolume
                      ,lXm        ,lXmcc
                      ,lDcca      
                      ,lvSkew     ,lmvSkew
                      ,lA         ,lB
                      ,&lRcell    ,ddt
                      ,lFaceVelR  ,lFacePresR    
                      ,lPres      ,lGradPres    
                      ,lVel       
                      ,lDfield    ,lWallPar        
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

/*...*/
  else
  {
/*... loop nas celulas*/
    aux2    = maxViz+1;
    for(nel=0;nel<numel;nel++)
    {
/*...*/
      if(calRcell)
        rCell[nel] = 0.e0;;
/*...*/
      aux1    = nFace[nel];
/*... elementos com equacoes*/
      if(MAT2D(nel,aux1,facePresR ,aux2) < 1)
      {

/*... loop na celula central*/    
        lMat            = mat[nel]-1;
        lib             = calType[lMat];
        lVolume[aux1]   = gVolume[nel]; 
        lGeomType[aux1] = geomType[nel];
        lFaceVelR[aux1] = MAT2D(nel,aux1,faceVelR ,aux2);
        lFacePresR[aux1]= MAT2D(nel,aux1,facePresR ,aux2);
 
        lPres[aux1]     = pres[nel];
        lId[aux1]       = id[nel] - 1;
/*...................................................................*/
      
/*...*/
        for(j=0;j<MAXPROP;j++)
          MAT2D(aux1,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
/*...................................................................*/
      
/*...*/
        lDfield[aux1] = dField[nel];
        for(j=0;j<ndm;j++){
          MAT2D(aux1,j,lGradPres,ndm) = MAT2D(nel,j,gradPres,ndm);
          MAT2D(aux1,j,lVel     ,ndm) = MAT2D(nel,j,vel     ,ndm);
          MAT2D(aux1,j,lDfield  ,ndm) = MAT2D(nel,j,dField  ,ndm);
        }
/*...................................................................*/

/*... propriedades por face*/
        for(i=0;i<aux1;i++)
        {
          lDcca[i]      = MAT2D(nel,i,gDcca   ,maxViz);
          lFaceVelR[i]  = MAT2D(nel,i,faceVelR ,aux2);
          lFacePresR[i] = MAT2D(nel,i,facePresR ,aux2);
/*... propriedades por face*/
          idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
          cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
          ch = OWNER(cellOwner, nel);
          lmKsi[i] = fModKsi[idFace];
          lfArea[i] = fArea[idFace];
          lmvSkew[i] = fModvSkew[idFace];

          for(j=0;j<ndm;j++)
          {
            MAT2D(i,j,lXmcc  ,ndm) = MAT3D(nel,i,j,gXmCc  ,maxViz,ndm);
/*... propriedades por face*/
            MAT2D(i, j, lKsi, ndm) = ch * MAT2D(idFace, j, fKsi, ndm);
            MAT2D(i, j, lEta, ndm) = ch * MAT2D(idFace, j, fEta, ndm);
            MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
            MAT2D(i, j, lXm, ndm) = MAT2D(idFace, j, fXm, ndm);
            MAT2D(i, j, lvSkew, ndm) = MAT2D(idFace, j, fvSkew, ndm);
          }
/*...................................................................*/

/*... loop na celulas vizinhas*/  
          vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
          lViz[i] = vizNel;
          lId[i]  = - 1;
          if( vizNel != -2) 
          {
            lVolume[i]    = gVolume[vizNel]; 
            lGeomType[i]  = geomType[vizNel];
            lMat          = mat[vizNel]-1;
            lPres[i]      = pres[vizNel];
            lId[i]        = id[vizNel] - 1;

            for(j=0;j<ndm;j++)
            {
              MAT2D(i,j,lGradPres,ndm) = MAT2D(vizNel,j,gradPres,ndm);
              MAT2D(i,j,lVel     ,ndm) = MAT2D(vizNel,j,vel     ,ndm);
              MAT2D(i,j,lDfield  ,ndm) = MAT2D(vizNel,j,dField  ,ndm);
            }

            for(j=0;j<MAXPROP;j++)
              MAT2D(i,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
          }
        }  
/*...................................................................*/

/*...*/
        if(fWallModel)
          for (i = 0; i<NWALLPAR; i++)
            lWallPar[i] = MAT2D(nel, i, wallPar, NWALLPAR);
/*...................................................................*/

/*... chamando a biblioteca de celulas*/
          cellLibSimplePres(loadsVel,loadsPres 
											,diffPres   
                      ,lGeomType  ,lProp 
                      ,lViz       ,lId           
                      ,lKsi       ,lmKsi
                      ,lEta       ,lfArea 
                      ,lNormal    ,lVolume
                      ,lXm        ,lXmcc
                      ,lDcca      
                      ,lvSkew     ,lmvSkew
                      ,lA         ,lB
                      ,&lRcell    ,ddt
                      ,lFaceVelR  ,lFacePresR  
                      ,lPres      ,lGradPres    
                      ,lVel       
                      ,lDfield    ,lWallPar        
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