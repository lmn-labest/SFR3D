#include<Turbulence.h>

/********************************************************************* 
 * Data de criacao    : 19/01/2018                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * kTurbPositve: mante o valor sempre positivos                      *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * kTurb    -> enercia cinetica turbulheta                           *
 * eddyVisc -> nao definico                                          *
 * density  -> densidade                                             *
 * volume   -> volume da celula                                      *
 * nEl      -> numero de elementos                                   *
 * ndm      -> numera de dimensao                                    *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * eddyVis -> viscosidade turbulenta                                 *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
static void eddyViscKtEq(Turbulence *tModel
                  , DOUBLE *RESTRICT kTurb  ,DOUBLE *RESTRICT eddyVisc
                  , DOUBLE *RESTRICT density,DOUBLE *RESTRICT volume
                  , DOUBLE *RESTRICT cDyn   
                  , INT const nEl           ,short const ndm) {

  bool fDynamic = tModel->fDynamic;
  INT i;
  DOUBLE lDensity,delta,ck=tModel->eK.ck;

/*...*/
  for (i = 0; i < nEl; i++) {

    if(kTurb[i] < 0.0)
      kTurb[i] = 0.e0;
    
    lDensity = MAT2D(i, 0, density, DENSITY_LEVEL);

    if(fDynamic)
      ck = MAT2D(i, 0, cDyn, 2);

/*...*/
    if(ndm == 3)
      delta = pow(volume[i],D1DIV3);
    else
      delta = sqrt(volume[i]);
/*...................................................................*/

    eddyVisc[i] = lDensity*ck*sqrt(kTurb[i])*delta;
//  printf("%d %e %e\n",i,ck,eddyVisc[i]);
  }
/*...................................................................*/

}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 30/11/2017                                   *
 * Data de modificaco : 19/07/2018                                   * 
 *-------------------------------------------------------------------* 
 * TURBULENCE: Calculo da viscosidae turbulenta                      *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m      -> memoria principal                                       *
 * pMesh  -> particionamento da malha                                * 
 * tModel -> turbulencia                                             * 
 * mesh   -> malha                                                   *
 * ndf    -> graus de liberdade                                      * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * eddyViscosity-> viscosidade turbulenta atualizada                 *
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void turbulence(Memoria *m            
              , Loads *lKturb         , Loads *lVel
              , PartMesh *pMesh       , Turbulence *tModel
              , Mesh *mesh            , Scheme *sc
              , Simple *sp            , SistEq *sistEq    
              , Solv *solv            , short const ndf)                  
{   
 
  short type = tModel->type,typeLes = tModel->typeLes;
  DOUBLE *nVel=NULL,*nDen=NULL;
  
/*...*/
  wallParametersCellLoop(lVel    , tModel             
      , mesh->elm.cellFace       , mesh->face.owner
      , mesh->elm.adj.nelcon     , mesh->elm.adj.nViz 
      , mesh->face.normal        , mesh->elm.geom.volume
      , mesh->elm.geom.dcca  
      , mesh->elm.faceRvel       , mesh->elm.faceLoadVel               
      , mesh->elm.vel            , mesh->elm.densityFluid
      , mesh->elm.dViscosity     , mesh->elm.wallParameters  
      , mesh->maxNo              , mesh->maxViz
      , mesh->ndm                , mesh->numelNov     
      , ndf); 
/*...................................................................*/

/*...*/
  switch (type)
  { 
/*...*/ 
    case LES:
/*...*/
      if(typeLes == LESFUNCMODEL )
      {
/*...*/ 
        if(tModel->fDynamic)
          lesDynamicMean(m                      , tModel 
                        , mesh->elm.adj.nelcon  , mesh->elm.nen    
                        , mesh->elm.adj.nViz    , mesh->elm.geom.volume  
                        , mesh->elm.vel         , mesh->elm.gradVel 
                        , mesh->elm.densityFluid, mesh->elm.dViscosity
                        , mesh->elm.cd
                        , mesh->maxNo           , mesh->maxViz         
                        , mesh->ndm             , mesh->numelNov 
                        , ndf                   , tModel->typeDynamic); 
/*...................................................................*/ 

/*...*/ 
        turbulenceCellLoop(lVel             , tModel   
                   , &pMesh->iEl              
                   , mesh->elm.node         , mesh->elm.adj.nelcon 
                   , mesh->elm.nen          , mesh->elm.adj.nViz 
                   , mesh->elm.cellFace     , mesh->face.owner
                   , mesh->elm.geom.volume  , mesh->elm.geom.dcca
                   , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
                   , mesh->face.mksi        , mesh->face.ksi
                   , mesh->face.eta         , mesh->face.area
                   , mesh->face.normal      , mesh->face.xm
                   , mesh->face.mvSkew      , mesh->face.vSkew
                   , mesh->elm.geomType     , mesh->elm.material.prop
                   , mesh->elm.material.type, mesh->elm.mat
                   , mesh->elm.faceRvel     , mesh->elm.faceLoadVel     
                   , mesh->elm.vel          , mesh->elm.gradVel      
                   , mesh->elm.densityFluid , mesh->elm.dViscosity
                   , mesh->elm.eddyViscosity, mesh->elm.wallParameters
                   , mesh->elm.cd   
                   , mesh->maxNo            , mesh->maxViz         
                   , mesh->ndm              , mesh->numelNov 
                   , ndf);   
/*...................................................................*/ 
      }
/*...................................................................*/ 

/*...*/
      else if (typeLes == LESMIXEDMODEL)
      {  
        HccaAlloc(DOUBLE,m,nVel, mesh->nnode*mesh->ndm,"nU"  ,_AD_); 
        HccaAlloc(DOUBLE,m,nDen, mesh->nnode*3  ,"nDen",_AD_); 

/*... modelo funcional*/ 
        turbulenceCellLoop(lVel            , tModel
                  , &pMesh->iEl  
                  , mesh->elm.node         , mesh->elm.adj.nelcon
                  , mesh->elm.nen          , mesh->elm.adj.nViz
                  , mesh->elm.cellFace     , mesh->face.owner
                  , mesh->elm.geom.volume  , mesh->elm.geom.dcca
                  , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
                  , mesh->face.mksi        , mesh->face.ksi
                  , mesh->face.eta         , mesh->face.area
                  , mesh->face.normal      , mesh->face.xm
                  , mesh->face.mvSkew      , mesh->face.vSkew
                  , mesh->elm.geomType     , mesh->elm.material.prop
                  , mesh->elm.material.type, mesh->elm.mat
                  , mesh->elm.faceRvel     , mesh->elm.faceLoadVel
                  , mesh->elm.vel          , mesh->elm.gradVel
                  , mesh->elm.densityFluid , mesh->elm.dViscosity
                  , mesh->elm.eddyViscosity, mesh->elm.wallParameters
                  , mesh->elm.cd
                  , mesh->maxNo            , mesh->maxViz
                  , mesh->ndm              , mesh->numelNov
                  , ndf);
/*...................................................................*/ 

/*... interpolacao das variaveis da celulas para pos nos (vel)*/
        if( tModel->typeMixed[ESTMODEL] == BARDINA)
        {
          interCellNode(m                , lVel
                     , mesh->elm.cellFace, mesh->face.owner
                     , nVel              , mesh->elm.vel        
                     , mesh->elm.node    , mesh->elm.geomType                      
                     , mesh->elm.geom.cc , mesh->node.x               
                     , mesh->face.xm           
                     , mesh->elm.nen     , mesh->elm.adj.nViz               
                     , mesh->elm.faceRvel, mesh->elm.faceLoadVel                
                     , &pMesh->iNo     
                     , mesh->numelNov    , mesh->numelNov              
                     , mesh->nnode       , mesh->nnode        
                     , mesh->maxNo       , mesh->maxViz           
                     , mesh->ndm         , 1
                     , mesh->ndm     
                     , false             , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (vel)*/
          interCellNode(m                , lVel
                     , mesh->elm.cellFace, mesh->face.owner
                     , nDen              , mesh->elm.densityFluid        
                     , mesh->elm.node    , mesh->elm.geomType                      
                     , mesh->elm.geom.cc , mesh->node.x               
                     , mesh->face.xm           
                     , mesh->elm.nen      , mesh->elm.adj.nViz               
                     , mesh->elm.faceRvel , mesh->elm.faceLoadVel                
                     , &pMesh->iNo      
                     , mesh->numelNov     , mesh->numelNov              
                     , mesh->nnode        , mesh->nnode        
                     , mesh->maxNo        , mesh->maxViz            
                     , 3                  , 1
                     , mesh->ndm     
                     , false              , 2);
/*...................................................................*/
        }
/*... modelo estrutual*/
        sLesModCellLoop(tModel  
                     , mesh->node.x        , mesh->elm.node
                     , mesh->elm.adj.nelcon, mesh->elm.nen 
                     , mesh->elm.adj.nViz  , mesh->elm.geom.volume    
                     , nVel                , mesh->elm.vel           
                     , nDen                , mesh->elm.densityFluid
                     , mesh->elm.gradVel   ,  mesh->elm.stressR 
                     , mesh->elm.cd
                     , mesh->maxNo         , mesh->maxViz         
                     , mesh->ndm           , mesh->numelNov 
                     , ndf                 , true); 
/*...................................................................*/ 

/*...*/
        HccaDealloc(m,nDen,"nDen",_AD_); 
        HccaDealloc(m,nVel,"nU",_AD_); 
 /*...................................................................*/
      }
/*...................................................................*/ 

/*...*/
      else if (typeLes == LESMIXEDTWOMODEL)
      {
        
        lesDynamicMean(m                      , tModel 
                      , mesh->elm.adj.nelcon  , mesh->elm.nen    
                      , mesh->elm.adj.nViz    , mesh->elm.geom.volume  
                      , mesh->elm.vel         , mesh->elm.gradVel 
                      , mesh->elm.densityFluid, mesh->elm.dViscosity
                      , mesh->elm.cd
                      , mesh->maxNo           , mesh->maxViz         
                      , mesh->ndm             , mesh->numelNov 
                      , ndf                   , tModel->typeDynamic); 

/*...*/ 
        turbulenceCellLoop(lVel             , tModel  
                   , &pMesh->iEl                      
                   , mesh->elm.node         , mesh->elm.adj.nelcon 
                   , mesh->elm.nen          , mesh->elm.adj.nViz 
                   , mesh->elm.cellFace     , mesh->face.owner
                   , mesh->elm.geom.volume  , mesh->elm.geom.dcca
                   , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
                   , mesh->face.mksi        , mesh->face.ksi
                   , mesh->face.eta         , mesh->face.area
                   , mesh->face.normal      , mesh->face.xm
                   , mesh->face.mvSkew      , mesh->face.vSkew
                   , mesh->elm.geomType     , mesh->elm.material.prop
                   , mesh->elm.material.type, mesh->elm.mat
                   , mesh->elm.faceRvel     , mesh->elm.faceLoadVel     
                   , mesh->elm.vel          , mesh->elm.gradVel      
                   , mesh->elm.densityFluid , mesh->elm.dViscosity
                   , mesh->elm.eddyViscosity, mesh->elm.wallParameters
                   , mesh->elm.cd   
                   , mesh->maxNo            , mesh->maxViz         
                   , mesh->ndm             , mesh->numelNov 
                   , ndf);   
/*...................................................................*/ 

/*... modelo estrutual*/
        sLesModCellLoop(tModel  
                     , mesh->node.x        , mesh->elm.node
                     , mesh->elm.adj.nelcon, mesh->elm.nen 
                     , mesh->elm.adj.nViz  , mesh->elm.geom.volume    
                     , nVel                , mesh->elm.vel           
                     , nDen                , mesh->elm.densityFluid
                     , mesh->elm.gradVel   , mesh->elm.stressR 
                     , mesh->elm.cd
                     , mesh->maxNo         , mesh->maxViz         
                     , mesh->ndm           , mesh->numelNov 
                     , ndf                 , false); 
/*...................................................................*/ 

      }
/*...................................................................*/   

/*...*/
      else if (typeLes == LESFUNCMODELONEEQK)
      {
        oneEquationK(m        
                   , lKturb  , lVel
                   , pMesh   , tModel
                   , mesh    , sc
                   , sp  
                   , sistEq  , solv); 
      }
/*...................................................................*/   

/*...*/
      else      
        ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__,
                  "LES Model: Opcao invalida !! ",EXIT_PROG);
/*...................................................................*/ 
    break;
/*...................................................................*/ 

/*...*/
    default:       
      ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__,
                "Turb Model: Opcao invalida !! ",EXIT_PROG);
/*...................................................................*/ 
  }

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 11/09/2017                                   *
 * Data de modificaco : 20/11/2017                                   * 
 *-------------------------------------------------------------------* 
 * oneEquationK: modelo de turbulencia de uma equacao                *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * eddyViscosity-> viscosidade turbulenta atualizada                 *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
  *********************************************************************/
void oneEquationK(Memoria *m        
                , Loads *lKturb     , Loads *lVel
                , PartMesh *pMesh   , Turbulence *tModel
                , Mesh *mesh        , Scheme *sc
                , Simple *sp  
                , SistEq *sistEq    , Solv *solv)                          
{   

  bool fRes = tModel->eK.fRes;
  short ndfVel;
  INT numel = mesh->numel;
  INT i,maxIt = tModel->eK.maxIt;
  DOUBLE rCell,rCell0,conv,tol=tModel->eK.tol;

  ndfVel = max(mesh->ndfF-1, mesh->ndfFt-2);

/*... loop nao-linear*/   
  for(i=0;i<maxIt;i++){

/*... dinamico oneKtEq*/ 
    if (tModel->fDynamic)
      lesDynamicOneKeq(m                    , tModel
                  , mesh->elm.adj.nelcon   , mesh->elm.nen    
                  , mesh->elm.adj.nViz     , mesh->elm.geom.volume 
                  , mesh->elm.vel          , mesh->elm.gradVel
                  , mesh->elm.densityFluid , mesh->elm.dViscosity
                  , mesh->elm.eddyViscosity, mesh->elm.cd                                      
                  , mesh->maxNo            , mesh->maxViz
                  , mesh->ndm              , mesh->numelNov   
                  , ndfVel                 , 1); 
/*...................................................................*/

/*... reconstruindo do gradiente*/
    rcGradU(m                       , lKturb
           , mesh->elm.node         , mesh->elm.adj.nelcon
           , mesh->node.x
           , mesh->elm.nen          , mesh->elm.adj.nViz
           , mesh->elm.cellFace     , mesh->face.owner
           , mesh->elm.geom.volume  , mesh->elm.geom.dcca
           , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
           , mesh->face.mksi        , mesh->face.ksi
           , mesh->face.eta         , mesh->face.area
           , mesh->face.normal      , mesh->face.xm
           , mesh->face.mvSkew      , mesh->face.vSkew
           , mesh->elm.geomType     , mesh->elm.material.prop
           , mesh->elm.material.type, mesh->elm.mat
           , mesh->elm.leastSquare  , mesh->elm.leastSquareR
           , mesh->elm.faceReKturb  , mesh->elm.faceLoadKturb    
           , mesh->elm.kTurb        , mesh->elm.gradKturb               
           , mesh->node.kTurb       , &sc->rcGrad
           , mesh->maxNo            , mesh->maxViz
           , 1                      , mesh->ndm
           , &pMesh->iNo            , &pMesh->iEl  
           , mesh->numelNov         , mesh->numel        
           , mesh->nnodeNov         , mesh->nnode); 
/*...................................................................*/           

/*... calculo de: A(i),b(i)*/
    systFormOneEqK(lKturb                , lVel
               , tModel
               , sc->advKturb            , sc->diffKturb
               , mesh->elm.node          , mesh->elm.adj.nelcon
               , mesh->elm.nen           , mesh->elm.adj.nViz
               , mesh->elm.cellFace      , mesh->face.owner
               , mesh->elm.geom.volume   , mesh->elm.geom.dcca
               , mesh->elm.geom.xmcc     , mesh->elm.geom.cc
               , mesh->face.mksi         , mesh->face.ksi
               , mesh->face.eta          , mesh->face.area
               , mesh->face.normal       , mesh->face.xm
               , mesh->face.mvSkew       , mesh->face.vSkew
               , mesh->elm.geomType      , mesh->elm.material.prop
               , mesh->elm.material.type , mesh->elm.mat
               , sistEq->ia              , sistEq->ja      
               , sistEq->al              , sistEq->ad       
               , sistEq->b               , sistEq->id       
               , mesh->elm.faceReKturb   , mesh->elm.faceLoadKturb  
               , mesh->elm.faceRvel      , mesh->elm.faceLoadVel  
               , mesh->elm.kTurb         , mesh->elm.gradKturb             
               , mesh->elm.vel           , mesh->elm.gradVel   
               , mesh->elm.pressure      , mesh->elm.gradPres
               , mesh->elm.densityFluid  , mesh->elm.dViscosity
               , mesh->elm.eddyViscosity , sp->d                                    
               , mesh->elm.rCellKturb    , mesh->elm.wallParameters
               , mesh->elm.cd            , &sc->ddt
               , sistEq->neq             , sistEq->neqNov      
               , sistEq->nad             , sistEq->nadr      
               , mesh->maxNo             , mesh->maxViz
               , mesh->ndm               , mesh->numelNov
               , 1                       , sistEq->storage
               , true                    , true   
               , true                    , sistEq->unsym);   
/*...................................................................*/

/*... soma o vetor b(i) = b(i) + b0(i)*/
    addVector(1.0e0           , sistEq->b
             ,1.0e0           , sistEq->b0
             ,sistEq->neqNov  , sistEq->b);
/*...................................................................*/

/*... soma o vetor R(i) = R(i) + b0(i)*/ 
    updateCellValue(mesh->elm.rCellKturb  , sistEq->b0
                   ,sistEq->id            , &sistEq->iNeq
                   ,mesh->numelNov        , 1
                   ,true                  , false);
/*...................................................................*/

/*...*/ 
    if( i == 0 ){
      rCell  = rCell0 = sqrt(dot(mesh->elm.rCellKturb
                                ,mesh->elm.rCellKturb 
                                ,mesh->numelNov));
      conv   = rCell0*tol;
    }
    else
      rCell  = sqrt(dot(mesh->elm.rCellKturb 
                       ,mesh->elm.rCellKturb 
                       ,mesh->numelNov));
        
    if(rCell < conv) break;
//  printf("%d %e %e\n",i,rCell,rCell/rCell0);
/*...................................................................*/

/*...*/
    solverC(m               
           ,sistEq->neq       ,sistEq->neqNov  
           ,sistEq->nad       ,sistEq->nadr
           ,sistEq->ia        ,sistEq->ja  
           ,sistEq->al        ,sistEq->ad,sistEq->au
           ,sistEq->b         ,sistEq->x
           ,&sistEq->iNeq     ,&sistEq->omp 
           ,solv->tol         ,solv->maxIt     
           ,sistEq->storage   ,solv->solver
           ,solv->fileSolv    ,solv->log  
           ,true              ,sistEq->unsym);  
/*...................................................................*/

/*... x -> kTur*/
    updateCellValue(mesh->elm.kTurb,sistEq->x
                   ,sistEq->id     ,&sistEq->iNeq
                   ,mesh->numel    ,1
                   ,fRes           ,true);
/*...................................................................*/

/*...*/
     eddyViscKtEq(tModel
                , mesh->elm.kTurb       , mesh->elm.eddyViscosity
                , mesh->elm.densityFluid, mesh->elm.geom.volume
                , mesh->elm.cd         
                , numel                 , mesh->ndm);
/*...................................................................*/

   }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 11/09/2017                                   *
 * Data de modificaco : 25/08/2019                                   * 
 *-------------------------------------------------------------------* 
 * turbulenceCellLoop: Calculo da viscosidae turbulenta              *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * loads   -> definicao de cargas de velocidade                      * 
 * tModel  -> modelo de turbulencia                                  *
 * iCel    -> interface de elementos                                 * 
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
 * gradVel   -> gradiente da solucao conhecido                       * 
 * density   -> massa especifica com variacao temporal               *  
 * eddyViscosity-> viscosidade turbulenta                            *
 * cDyn      -> coeficiente dinamicamente determinado                *
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * numel     -> numero de toral de celulas                           * 
 * ndf       -> graus de liberdade                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * eddyViscosity-> viscosidade turbulenta atualizada                 *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
  *********************************************************************/
void turbulenceCellLoop(Loads *lVel         , Turbulence *tModel  
      , Interface *iCel           
      , INT    *RESTRICT el                 , INT    *RESTRICT nelcon 
      , short  *RESTRICT nen                , short  *RESTRICT nFace 
      , INT *RESTRICT cellFace              , INT *RESTRICT fOwner
      , DOUBLE *RESTRICT gVolume            , DOUBLE *RESTRICT gDcca
      , DOUBLE *RESTRICT gXmCc              , DOUBLE *RESTRICT gCc
      , DOUBLE *RESTRICT fModKsi            , DOUBLE *RESTRICT fKsi
      , DOUBLE *RESTRICT fEta               , DOUBLE *RESTRICT fArea
      , DOUBLE *RESTRICT fNormal            , DOUBLE *RESTRICT fXm
      , DOUBLE *RESTRICT fModvSkew          , DOUBLE *RESTRICT fvSkew 
      , short  *RESTRICT geomType           , DOUBLE *RESTRICT prop
      , short  *RESTRICT calType            , short  *RESTRICT mat
      , short  *RESTRICT faceVelR           , short *RESTRICT faceVelL               
      , DOUBLE *RESTRICT vel                , DOUBLE *RESTRICT gradVel
      , DOUBLE *RESTRICT density            , DOUBLE *RESTRICT dViscosity        
      , DOUBLE *RESTRICT eddyViscosity      , DOUBLE *RESTRICT wallPar  
      , DOUBLE *RESTRICT cDyn                 
      , short const maxNo                   , short const maxViz
      , short const ndm                     , INT const numel     
      , short const ndf)                      
{   
  short i,j,k;
  bool fWall    = tModel->fWall,
       fDynamic = tModel->fDynamic; 
  INT nel,vizNel;
  
/*... variavel local */
  short  aux1,aux2,lMat,lib;
  short  lGeomType[MAX_NUM_FACE+1],  
         lFaceVelR[MAX_NUM_FACE+1],lFaceVelL[MAX_NUM_FACE+1];       
  INT    lViz[MAX_NUM_FACE];
  INT    idFace, cellOwner, ch;
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM], lmKsi[MAX_NUM_FACE],
         lEta[MAX_NUM_FACE*MAX_NDM], lfArea[MAX_NUM_FACE],
         lNormal[MAX_NUM_FACE*MAX_NDM], lVolume[MAX_NUM_FACE+1],
         lXm[MAX_NUM_FACE*MAX_NDM], lXmcc[MAX_NUM_FACE*MAX_NDM],
         lDcca[MAX_NUM_FACE], lmvSkew[MAX_NUM_FACE],
         lvSkew[MAX_NUM_FACE*MAX_NDM], lDensity[(MAX_NUM_FACE+1)],
         lProp[(MAX_NUM_FACE+1)*MAXPROP],
         lCc[(MAX_NUM_FACE+1)*MAX_NDM],
         lGradVel[(MAX_NUM_FACE+1)*MAX_NDM*MAX_NDF],
         lVel0[(MAX_NUM_FACE+1)*MAX_NDM],
         lEddyViscosity,lDviscosity,lWallPar[4],lcDyn=0.e0;
           
/*... loop nas celulas*/
  aux2    = maxViz+1;
  for(nel=0;nel<numel;nel++)
  {
/*...*/
    aux1    = nFace[nel];
/*... loop na celula central*/    
    lMat                 = mat[nel]-1;
    lib                  = calType[lMat];
    lVolume[aux1]        = gVolume[nel]; 
    lGeomType[aux1]      = geomType[nel];
    lFaceVelR[aux1]      = MAT2D(nel,aux1,faceVelR ,aux2);
    lFaceVelL[aux1]      = MAT2D(nel,aux1,faceVelL ,aux2);
    lDensity[aux1]       = MAT2D(nel,2   ,density ,DENSITY_LEVEL);
    
/*...*/
    lDviscosity       = dViscosity[nel];
/*...*/      
    for(j=0;j<MAXPROP;j++)
      MAT2D(aux1,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
/*...................................................................*/

/*...*/
    for(j=0;j<ndm;j++)
    {
       MAT2D(aux1,j,lVel0,ndm)  = MAT2D(nel, j, vel, ndm);
       MAT2D(aux1, j, lCc, ndm) = MAT2D(nel, j, gCc, ndm);
    }
/*...................................................................*/

/*...*/
    for(i=0;i<ndf;i++)
      for(j=0;j<ndm;j++)
        MAT3D(aux1, i, j, lGradVel, ndf, ndm) 
                = MAT3D(nel, i, j, gradVel, ndf, ndm);
/*...................................................................*/

/*...*/
    for(i=0;i<aux1;i++)
    {
      lDcca[i]      = MAT2D(nel,i,gDcca   ,maxViz);
      lFaceVelR[i]  = MAT2D(nel,i,faceVelR,aux2);
      lFaceVelL[i]  = MAT2D(nel,i,faceVelL,aux2);
/*... propriedades por face*/
      idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
      cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
      ch = OWNER(cellOwner, nel);
      lmKsi[i] = fModKsi[idFace];
      lfArea[i] = fArea[idFace];
      lmvSkew[i] = fModvSkew[idFace];
      for(j=0;j<ndm;j++)
      {
        MAT2D(i,j,lXmcc  ,ndm)    = MAT3D(nel,i,j,gXmCc  ,maxViz,ndm);
/*... propriedades por face*/
        MAT2D(i, j, lKsi, ndm)    = ch * MAT2D(idFace, j, fKsi, ndm);
        MAT2D(i, j, lEta, ndm)    = ch * MAT2D(idFace, j, fEta, ndm);
        MAT2D(i, j, lNormal, ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
        MAT2D(i, j, lXm, ndm)     = MAT2D(idFace, j, fXm, ndm);
        MAT2D(i, j, lvSkew, ndm)  = MAT2D(idFace, j, fvSkew, ndm);
      }
    }

/*... loop na celulas vizinhas*/    
    for(i=0;i<aux1;i++)
    {
      vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
      lViz[i] = vizNel;
      if( vizNel != -2)
      {
        lVolume[i]   = gVolume[vizNel]; 
        lGeomType[i] = geomType[vizNel];
        lDensity[i]  = MAT2D(vizNel,2   ,density ,DENSITY_LEVEL);

        lMat         = mat[vizNel]-1;
 
        for(j=0;j<ndm;j++)
        {
          MAT2D(i,j,lCc      ,ndm)  = MAT2D(vizNel, j, gCc, ndm);
          MAT2D(i,j,lVel0     ,ndm) = MAT2D(vizNel, j, vel, ndm);
        }
      
        for(k=0;k<ndf;k++)
          for(j=0;j<ndm;j++)
            MAT3D(i,k,j,lGradVel,ndf,ndm) 
                           = MAT3D(vizNel,k,j,gradVel,ndf,ndm);

        for(j=0;j<DIFPROP;j++)
          MAT2D(i,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
      }
    }  
/*...................................................................*/
     if(fDynamic)
       lcDyn = MAT2D(nel,0,cDyn,2);
/*...*/
     if(fWall)
       for(j=0;j<NWALLPAR;j++)
         lWallPar[j] = MAT2D(nel,j,wallPar,NWALLPAR);


/*... chamando a biblioteca de celulas*/
     cellLibTurbulence(lVel            , *tModel        
                      , lGeomType      , lProp
                      , lViz           , lKsi
                      , lmKsi          , lEta  
                      , lfArea         , lNormal        
                      , lVolume        , lXm 
                      , lXmcc          , lDcca  
                      , lCc            , lvSkew         
                      , lmvSkew        , lFaceVelR 
                      , lFaceVelL      , lVel0   
                      , lGradVel       , lDensity
                      , lDviscosity    , &lEddyViscosity
                      , lWallPar       , lcDyn
                      , nen[nel]       , nFace[nel] 
                      , ndm            , lib   
                      , nel);    
/*...................................................................*/
 
/*...*/
    eddyViscosity[nel]  = lEddyViscosity;
  }
/*...................................................................*/

/*...*/
  if(mpiVar.nPrcs > 1 ) comunicateCel(iCel,eddyViscosity,1,1);
/*...................................................................*/

}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 17/01/2018                                   *
 * Data de modificaco : 19/07/2018                                   * 
 *-------------------------------------------------------------------* 
 * wallParametersCellLoop: Calculo dos paramentros de parede         *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*  
 * nelcon    -> vizinhos dos elementos                               * 
 * nFace     -> numero de faces por celulas                          * 
 * cellFace-> faces que formam a celulas                             *
 * owner   -> elementos que compartilham a face(0- o dono,1 - viz)   *
 * fNormal   -> vetores normais as faces das celulas                 * 
 * gVolume   -> volumes das celulas                                  * 
 * gDcca     -> menor distancia do centroide a faces desta celula    * 
 * density   -> massa especifica com variacao temporal               *  
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * numel     -> numero de toral de celulas                           * 
 * ndf       -> graus de liberdade                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
  *********************************************************************/
void wallParametersCellLoop(Loads *lVel, Turbulence *tModel 
      , INT *RESTRICT cellFace         , INT *RESTRICT fOwner
      , INT    *RESTRICT nelcon        , short  *RESTRICT nFace 
      , DOUBLE *RESTRICT fNormal       , DOUBLE *RESTRICT gVolume
      , DOUBLE *RESTRICT gDcca  
      , short  *RESTRICT faceVelR      , short *RESTRICT faceVelL               
      , DOUBLE *RESTRICT vel           , DOUBLE *RESTRICT density
      , DOUBLE *RESTRICT dViscosity    , DOUBLE *RESTRICT wallPar  
      , short const maxNo              , short const maxViz
      , short const ndm                , INT const numel     
      , short const ndf)                      
{   
  bool fWall = false; 
  short i,j,wallType;
  INT nel,vizNel;
  
/*... variavel local */
  short  aux1,aux2;
  short  lFaceVelR[MAX_NUM_FACE+1],lFaceVelL[MAX_NUM_FACE+1];       
  INT    lViz[MAX_NUM_FACE];
  INT    idFace, cellOwner, ch;
  DOUBLE lNormal[MAX_NUM_FACE*MAX_NDM],lDcca[MAX_NUM_FACE],
         lDensityC,lVel0[MAX_NDM],lDviscosityC,lWallPar[4],dMin;
 
  wallType   = tModel->wallType;  
  dMin = 0.0e0;        
/*... loop nas celulas*/
  aux2    = maxViz+1;
  for(nel=0;nel<numel;nel++)
  {
/*...*/
    aux1    = nFace[nel];
/*... loop na celula central*/    
    lFaceVelR[aux1] = MAT2D(nel,aux1,faceVelR ,aux2);
    lFaceVelL[aux1] = MAT2D(nel,aux1,faceVelL ,aux2);
    lDensityC       = MAT2D(nel,2   ,density ,DENSITY_LEVEL);
    
/*...*/
    lDviscosityC = dViscosity[nel];

/*...*/
    for(j=0;j<ndm;j++)
      lVel0[j] = MAT2D(nel,j,vel,ndm);
/*...................................................................*/

/*...*/
    for(i=0;i<aux1;i++)
    {
      lDcca[i]      = MAT2D(nel,i,gDcca   ,maxViz);
      lFaceVelR[i]  = MAT2D(nel,i,faceVelR,aux2);
      lFaceVelL[i]  = MAT2D(nel,i,faceVelL,aux2);

      idFace    = MAT2D(nel, i, cellFace, maxViz) - 1;
      cellOwner = MAT2D(idFace, 0, fOwner, 2) - 1;
      ch        = OWNER(cellOwner, nel);

      for(j=0;j<ndm;j++)
        MAT2D(i,j,lNormal,ndm) = ch * MAT2D(idFace, j, fNormal, ndm);
    }

/*... loop na celulas vizinhas*/    
    for(i=0;i<aux1;i++)
    {
      vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
      lViz[i] = vizNel;
    }  
/*...................................................................*/

/*...*/    
    for(i=0;i<NWALLPAR;i++)
      lWallPar[i] = 0.e0;
/*...................................................................*/

/*...*/
    fWall = wallParameters(lVel               
                         ,lViz         ,lVel0
                         ,lNormal      ,lDcca
                         ,lFaceVelR    ,lFaceVelL
                         ,lDviscosityC ,lDensityC
                         ,lWallPar     ,&dMin
                         ,wallType     ,nFace[nel]);
/*...................................................................*/
 
/*...*/
    if(fWall)
    {
      for(i=0;i<NWALLPAR;i++)
        MAT2D(nel,i,wallPar,NWALLPAR) = lWallPar[i];
    }
  }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 30/11/2017                                   *
 * Data de modificaco : 13/01/2018                                   * 
 *-------------------------------------------------------------------* 
 * lesDynamicMean : cacula os  M:L e M:M separados para fazer a media*
 * do LES dinamico                                                   *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m         -> memoria principal                                    *
 * tModel    -> turbulencia                                          * 
 * nelcon    -> vizinhos dos elementos                               * 
 * nen       -> numero de nos por celulas                            * 
 * nFace     -> numero de faces por celulas                          * 
 * gradVel   -> gradiente da solucao conhecido                       * 
 * density   -> massa especifica com variacao temporal               *
 * dViscosity-> viscosidade dinamica                                 * 
 * cDyn      -> nao definido                                         * 
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * numel     -> numero de toral de celulas                           * 
 * ndf       -> graus de liberdade                                   * 
 * iCod      -> tipo de coeficiente dinamico                         *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * cDyn     -> coeficiente calculados dinamicamente (cf,cs)          *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 * 1 paramentro:                                                     *
 * dynamic = ( L:M, M:M)                                             *
 * 2 paramentro:                                                     *
 * dynamic = ( L:Ms, L:Mf, Ms:Ms, Ms:Mf, Mf:Ms, Mf:Mf )              *
 *********************************************************************/
void lesDynamicMean(Memoria *m               , Turbulence *tModel
                  , INT    *RESTRICT nelcon  , short  *RESTRICT nen    
                  , short  *RESTRICT nFace   , DOUBLE *RESTRICT gVolume 
                  , DOUBLE *RESTRICT vel     , DOUBLE *RESTRICT gradVel 
                  , DOUBLE *RESTRICT density , DOUBLE *RESTRICT dViscosity
                  , DOUBLE *RESTRICT cDyn                                      
                  , short const maxNo        , short const maxViz
                  , short const ndm          , INT const numel     
                  , short const ndf          , short const iCod)                      
{   
  short i,j,k,nPar=1;
  INT nel,vizNel;
  
/*... variavel local */
  short  aux1;
  INT    lViz[MAX_NUM_FACE];
  DOUBLE lVolume[MAX_NUM_FACE+1],lDensity[(MAX_NUM_FACE+1)],gmm[4],gCoef,
         lGradVel[(MAX_NUM_FACE+1)*MAX_NDM*MAX_NDF],lDyn[2],
         lVel[(MAX_NUM_FACE+1)*MAX_NDM],lDynamic[(MAX_NUM_FACE+1)*6],
         lViscosity;   
  DOUBLE *dynamic=NULL;
   
/*...*/
  if(iCod == LDYNAMIC || iCod == GDYNAMIC )
    nPar = 3;
  else if(iCod == GDYNAMICMOD )  
    nPar = 4;
  else if(iCod == TWOPARDYNAMIC)
    nPar = 6;
/*...................................................................*/
    
  HccaAlloc(DOUBLE,m,dynamic, numel*nPar,"dynamic",_AD_);
    
/*... loop nas celulas*/
  for(nel=0;nel<numel;nel++){
/*...*/
    aux1    = nFace[nel];

/*... loop na celula central*/    
    lVolume[aux1]     = gVolume[nel]; 
    lDensity[aux1]    = MAT2D(nel,2   ,density ,DENSITY_LEVEL);

/*...*/
    for(j=0;j<ndm;j++){
       MAT2D(aux1,j,lVel,ndm)  = MAT2D(nel,j,vel     ,ndm);
    }
/*...................................................................*/

/*...*/
    for(i=0;i<ndf;i++)
      for(j=0;j<ndm;j++)
        MAT3D(aux1, i, j, lGradVel, ndf, ndm) 
                = MAT3D(nel, i, j, gradVel, ndf, ndm);
/*...................................................................*/

/*... loop na celulas vizinhas*/    
    for(i=0;i<aux1;i++){
      vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
      lViz[i] = vizNel;
      if( vizNel != -2) {
        lVolume[i]   = gVolume[vizNel]; 
        lDensity[i]  = MAT2D(vizNel,2   ,density ,DENSITY_LEVEL);

        for(j=0;j<ndm;j++){
          MAT2D(i,j,lVel      ,ndm) = MAT2D(vizNel,j,vel     ,ndm);
        }
      
        for(k=0;k<ndf;k++)
          for(j=0;j<ndm;j++)
            MAT3D(i,k,j,lGradVel,ndf,ndm) 
                           = MAT3D(vizNel,k,j,gradVel,ndf,ndm);
      }
    }  
/*...................................................................*/

/*...*/
    if(iCod == LDYNAMIC || iCod == GDYNAMIC)
      lesDynOnePar(lViz                       , lVolume
                 , lDensity                   , lVel
                 , lGradVel                   , lDynamic 
                 , tModel->typeMixed[FUNMODEL], nFace[nel]);
    else if(iCod == GDYNAMICMOD)
     lesDynOneParG(lViz                       , lVolume
                 , lGradVel                   , lDynamic 
                 , tModel->typeMixed[FUNMODEL], nFace[nel]);
    else if(iCod == TWOPARDYNAMIC)
      lesDynTwoPar(lViz       , lVolume
                 , lDensity   , lVel
                 , lGradVel   , lDynamic
                 , nFace[nel]  );
/*...................................................................*/

/// if (nel == 3099 || nel == 62299) {
//    printf("%d %e %e %lf %e\n",nel,lDynamic[0],lDynamic[1]
//                               ,lDynamic[2],lDynamic[0]/lDynamic[1]);
//  }

/*...*/
    for(k=0;k<nPar;k++)
      MAT2D(nel,k,dynamic,nPar) = lDynamic[k];
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  gmm[0] = gmm[1] = gmm[2] = gmm[3] = 0.e0;
  for (nel = 0; nel < numel; nel++) {
/*...*/
    aux1    = nFace[nel];
    for(k=0;k<nPar;k++)
      MAT2D(aux1,k,lDynamic,nPar) = MAT2D(nel,k,dynamic,nPar); 

/*... loop na celulas vizinhas*/    
    for(i=0;i<aux1;i++){
      vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1; 
      lViz[i] = vizNel;
      if( vizNel != -2){
        for(k=0;k<nPar;k++)
          MAT2D(i,k,lDynamic,nPar) = MAT2D(vizNel,k,dynamic,nPar);
      }
    }  
/*...................................................................*/

/*...*/
    if(iCod == LDYNAMIC)
      MAT2D(nel,0,cDyn,2) = oneParLes(tModel 
                                    , lViz      , lDynamic                                    
                                    , nFace[nel], true);
/*...................................................................*/

/*...*/
    else if (iCod == GDYNAMIC) {
/*... <L:M>global:*/
      gmm[0] += lDynamic[0];
/*... <M:M>global:*/
      gmm[1] += lDynamic[1];
    }
/*...................................................................*/

/*...*/
    else if (iCod == GDYNAMICMOD) {
/*... <tesFilter(a:a)>global:*/
      gmm[0] += lDynamic[0];
/*... <tesFilter(a):tesFilter(a)>global:*/
      gmm[1] += lDynamic[1];
/*... <op*delta*tesFilter(S:S)>global:*/
      gmm[2] += lDynamic[2];
/*... <opT*deltaT*tesFilter(S):tesFilter(S)>global:*/
      gmm[3] += lDynamic[3];
    }
/*...................................................................*/

/*...*/
    else if(iCod == TWOPARDYNAMIC){
      twoParLes( tModel
               , lViz       , lDynamic
               , lDyn
               , nFace[nel], true);
      MAT2D(nel,0,cDyn,2) = lDyn[0];
      MAT2D(nel,1,cDyn,2) = lDyn[1];
    }
 /*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if (iCod == GDYNAMIC){
    gCoef = oneGoParLes(tModel, gmm, true);
    for (nel = 0; nel < numel; nel++)
      MAT2D(nel,0,cDyn,2) = gCoef;
  }
  else if(iCod == GDYNAMICMOD ){
    for (nel = 0; nel < numel; nel++){
      lDensity[0] = MAT2D(nel,2   ,density ,DENSITY_LEVEL); 
      lViscosity = dViscosity[nel];
      MAT2D(nel,0,cDyn,2) = oneGoParLesMod(tModel     , gmm
                                         , lDensity[0], lViscosity
                                         , false);
    }
  }
/*...................................................................*/

  HccaDealloc(m,dynamic,"dynamic",_AD_);
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 30/11/2017                                   *
 * Data de modificaco : 13/01/2018                                   * 
 *-------------------------------------------------------------------* 
 * lesDynamicOneKeq : calculo dos coeficientes dinamicos para o LES  *
 * de uma equacao K                                                  *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m         -> memoria principal                                    *
 * tModel    -> turbulencia                                          * 
 * nelcon    -> vizinhos dos elementos                               * 
 * nen       -> numero de nos por celulas                            * 
 * nFace     -> numero de faces por celulas                          * 
 * gradVel   -> gradiente da solucao conhecido                       * 
 * density   -> massa especifica com variacao temporal               *
 * dViscosity-> viscosidade dinamica                                 * 
 * eddyVisc  -> viscosidade turbulenta                               * 
 * cDyn      -> nao definido                                         * 
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * numel     -> numero de toral de celulas                           * 
 * ndf       -> graus de liberdade                                   * 
 * iCod      -> tipo de coeficiente dinamico                         *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * cDyn     -> coeficiente calculados dinamicamente (ck,ce)          *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 * 1 paramentro:                                                     *
 * dynamic = ( L:L, M:M, Ce)                                         *
 *********************************************************************/
void lesDynamicOneKeq(Memoria *m             , Turbulence *tModel
                  , INT    *RESTRICT nelcon  , short  *RESTRICT nen    
                  , short  *RESTRICT nFace   , DOUBLE *RESTRICT gVolume 
                  , DOUBLE *RESTRICT vel     , DOUBLE *RESTRICT gradVel 
                  , DOUBLE *RESTRICT density , DOUBLE *RESTRICT dViscosity
                  , DOUBLE *RESTRICT eddyVis , DOUBLE *RESTRICT cDyn      
                  , short const maxNo        , short const maxViz
                  , short const ndm          , INT const numel     
                  , short const ndf          , short const iCod)                      
{   
  short i,j,k,nPar;
  INT nel,vizNel;
  
/*... variavel local */
  short  aux1;
  INT    lViz[MAX_NUM_FACE];
  DOUBLE lVolume[MAX_NUM_FACE+1],lDensity[(MAX_NUM_FACE+1)],
         lGradVel[(MAX_NUM_FACE+1)*MAX_NDM*MAX_NDF],lDyn[2],
         lVel[(MAX_NUM_FACE+1)*MAX_NDM],lDynamic[(MAX_NUM_FACE+1)*6];   
  DOUBLE *dynamic=NULL;
   
/*...*/
  nPar = 4;
/*...................................................................*/
    
  HccaAlloc(DOUBLE,m,dynamic, numel*nPar,"dynamic",_AD_);
    
/*... loop nas celulas*/
  for(nel=0;nel<numel;nel++){
/*...*/
    aux1    = nFace[nel];

/*... loop na celula central*/    
    lVolume[aux1]  = gVolume[nel]; 
    lDensity[aux1] = MAT2D(nel,2   ,density ,DENSITY_LEVEL);

/*...*/
    for(j=0;j<ndm;j++){
       MAT2D(aux1,j,lVel,ndm) = MAT2D(nel,j,vel     ,ndm);
    }
/*...................................................................*/

/*...*/
    for(i=0;i<ndf;i++)
      for(j=0;j<ndm;j++)
        MAT3D(aux1, i, j, lGradVel, ndf, ndm) 
                = MAT3D(nel, i, j, gradVel, ndf, ndm);
/*...................................................................*/

/*... loop na celulas vizinhas*/    
    for(i=0;i<aux1;i++){
      vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
      lViz[i] = vizNel;
      if( vizNel != -2) {
        lVolume[i]   = gVolume[vizNel]; 
        lDensity[i]  = MAT2D(vizNel,2   ,density ,DENSITY_LEVEL);

        for(j=0;j<ndm;j++){
          MAT2D(i,j,lVel      ,ndm) = MAT2D(vizNel,j,vel     ,ndm);
        }
      
        for(k=0;k<ndf;k++)
          for(j=0;j<ndm;j++)
            MAT3D(i,k,j,lGradVel,ndf,ndm) 
                           = MAT3D(vizNel,k,j,gradVel,ndf,ndm);
      }
    }  
/*...................................................................*/

/*...*/
    lesOneKeqLoc(lViz           , lVolume
               , lDensity       , lVel
               , lGradVel       , lDynamic
               , dViscosity[nel], eddyVis[nel]
               , nel            , nFace[nel]);
/*...................................................................*/

/*...*/
    for(k=0;k<nPar;k++)
      MAT2D(nel,k,dynamic,nPar) = lDynamic[k];
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  for (nel = 0; nel < numel; nel++) {
/*...*/
    aux1    = nFace[nel];
    for(k=0;k<nPar;k++)
      MAT2D(aux1,k,lDynamic,nPar) = MAT2D(nel,k,dynamic,nPar); 
      lVolume[aux1]  = gVolume[nel];
/*... loop na celulas vizinhas*/    
    for(i=0;i<aux1;i++){
      vizNel     = MAT2D(nel,i,nelcon,maxViz) - 1; 
      lVolume[i] = gVolume[vizNel]; 
      lViz[i] = vizNel;
      if( vizNel != -2){
        for(k=0;k<nPar;k++)
          MAT2D(i,k,lDynamic,nPar) = MAT2D(vizNel,k,dynamic,nPar);
      }
    }  
/*...................................................................*/

/*...*/
    lesOneKeqPar( tModel   , lViz 
               , lDynamic  , lDyn
               , lVolume   , nFace[nel]
               , nel       , nPar
               , true);
    MAT2D(nel,0,cDyn,2) = lDyn[0];
    MAT2D(nel,1,cDyn,2) = lDyn[1];
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  HccaDealloc(m,dynamic,"dynamic",_AD_);
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 04/12/2017                                   *
 * Data de modificaco : 12/12/2017                                   * 
 *-------------------------------------------------------------------* 
 * sLesModCellLoop : cacula do tensor residual por metodos           *
 * estrutarais                                                       *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * tModel    -> modelo de turbulencia                                * 
 * x         -> cordenadas dos pontos                                * 
 * el      -> conetividade dos celulas                               * 
 * nelcon    -> vizinhos dos elementos                               * 
 * nen       -> numero de nos por celulas                            * 
 * nFace     -> numero de faces por celulas                          * 
 * nVel      -> velocidade nodal                                     *
 * eVel      -> velocidade por elemento                              *
 * nDensity  -> massa especifica com variacao temporal (nodal)       *
 * eDensity  -> massa especifica com variacao temporal (celula)      * 
 * gradVel   -> gradiente da solucao conhecido                       * 
 * stressR   -> nao definido                                         *
 * cd        -> coeficientes calculados dinamicamente                *   
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * numel     -> numero de toral de celulas                           * 
 * ndf       -> graus de liberdade                                   *
 * fNode     -> valores nodais                                       *  
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * stressR  -> tensor residual                                       *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 * stressR(s11,s22,s33,s12,s23,s13)                                  *
 * stressR(sxx,syy,szz,sxy,syz,sxz)                                  *
 *                                                                   *
 *      | x   y  z |  no1                                            *
 * lx = | x   y  z |  no2                                            *
 *      | x   y  z |  no3                                            *
 *                                                                   *
 *       | u   v  w |  no1                                           *
 * vel = | u   v  w |  no2                                           *
 *       | u   v  w |  no3                                           *
 *                                                                   *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void sLesModCellLoop(Turbulence *tModel      
                   , DOUBLE *RESTRICT x       , INT *RESTRICT el       
                   , INT *RESTRICT nelcon     , short  *RESTRICT nen    
                   , short *RESTRICT nFace    , DOUBLE *RESTRICT gVolume 
                   , DOUBLE *RESTRICT nVel    , DOUBLE *RESTRICT eVel 
                   , DOUBLE *RESTRICT nDensity, DOUBLE *RESTRICT eDensity  
                   , DOUBLE *RESTRICT gradVel , DOUBLE *RESTRICT stressR   
                   , DOUBLE *RESTRICT cDyn                               
                   , short const maxNo        , short const maxViz
                   , short const ndm          , INT const numel     
                   , short const ndf          , bool const fNode)                      
{   
  short i,j,k;
  INT nel,vizNel;
  
/*... variavel local */
  short  aux1;
  INT    no,lViz[MAX_NUM_FACE];
  DOUBLE lVolume[MAX_NUM_FACE+1],lDensity[(MAX_NUM_FACE+1)],
         lGradVel[(MAX_NUM_FACE+1)*MAX_NDM*MAX_NDF],
         lVel[(MAX_NUM_FACE+1)*MAX_NDM],lStress[6],
         lx[(MAX_NUM_NODE)*MAX_NDM],lnVel[(MAX_NUM_NODE)*MAX_NDM],
         lnDensity[MAX_NUM_NODE],lcDyn;
           
/*... loop nas celulas*/

  for(nel=0;nel<numel;nel++){
/*...*/
    aux1    = nFace[nel];

/*... coordenadas dos pontos dos verticeis da celula central*/
    if(fNode)
      for(j=0;j<nen[nel];j++){
/*...*/
        no = MAT2D(nel,j,el,maxNo) - 1;

        lnDensity[j] = MAT2D(no,2,nDensity ,DENSITY_LEVEL);
        for(k=0;k<ndm;k++){
          MAT2D(j,k,lx,ndm)    = MAT2D(no,k,x,ndm);
          MAT2D(j,k,lnVel,ndm) = MAT2D(no,k,nVel,ndm);
        }     
      }
/*...................................................................*/

/*... loop na celula central*/    
    lVolume[aux1]     = gVolume[nel]; 
    lDensity[aux1]    = MAT2D(nel,2   ,eDensity ,DENSITY_LEVEL);

/*...*/
    for(j=0;j<ndm;j++){
       MAT2D(aux1,j,lVel,ndm)  = MAT2D(nel,j,eVel     ,ndm);
    }
/*...................................................................*/

/*...*/
    for(i=0;i<ndf;i++)
      for(j=0;j<ndm;j++)
        MAT3D(aux1, i, j, lGradVel, ndf, ndm) 
                = MAT3D(nel, i, j, gradVel, ndf, ndm);
/*...................................................................*/

/*... loop na celulas vizinhas*/    
    for(i=0;i<aux1;i++){
      vizNel  = MAT2D(nel,i,nelcon,maxViz) - 1;
      lViz[i] = vizNel;
      if( vizNel != -2) {
        lVolume[i]   = gVolume[vizNel]; 
        lDensity[i]  = MAT2D(vizNel,2   ,eDensity ,DENSITY_LEVEL);

        for(j=0;j<ndm;j++){
          MAT2D(i,j,lVel,ndm) = MAT2D(vizNel,j,eVel,ndm);
        }
      
        for(k=0;k<ndf;k++)
          for(j=0;j<ndm;j++)
            MAT3D(i,k,j,lGradVel,ndf,ndm) 
                           = MAT3D(vizNel,k,j,gradVel,ndf,ndm);
      }
    }  
/*...................................................................*/
  
    lcDyn = MAT2D(nel,1,cDyn,2);

/*...*/
    structuralStress(*tModel      , lx
                   , lViz         , lStress  
                   , lnVel        , lVel   
                   , lnDensity    , lDensity     
                   , lGradVel     , lVolume
                   , lcDyn
                   , ndm          , nFace[nel]
                   , nel);
/*...................................................................*/

/*...*/
    if(ndm == 3){
      for(k=0;k<6;k++)
        MAT2D(nel,k,stressR,6) = lStress[k];
    }  
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 11/09/2017                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* CELLLES: Calculo da viscosidae turbulenta para o LES              *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* lnFace    -> numero de faces da celula central e seus vizinhos    *
* lGeomType -> tipo geometrico da celula central e seus vizinhos    *
* lprop     -> propriedade fisicas das celulas                      *
* lViz      -> viznhos da celula central                            *
* lId       -> numeracoes das equacoes das celulas                  *
* Ksi       -> vetores que unem centroide da celula central aos     *
*              vizinhos desta                                       *
* mKsi      -> modulo do vetor ksi                                  *
* eta       -> vetores paralelos as faces das celulas               *
* mEta      -> modulo do vetor eta                                  *
* normal    -> vetores normais as faces das celulas                 *
* area      -> area da celula central                               *
* xm        -> pontos medios das faces da celula central            *
* xmcc      -> vetores que unem o centroide aos pontos medios das   *
*              faces da celula central                              *
* vSkew     -> vetor entre o ponto medio a intersecao que une os    *
*            centrois compartilhado nessa face da celula central    *
* mvSkew    -> distacia entre o ponto medio a intersecao que une os *
*              centrois compartilhado nessa face da celula central  *
* dcca      -> menor distancia do centroide central a faces desta   *
*              celula                                               *
* cc        -> centroides da celula centra e seus vizinhos          *
* lDensity  -> massa especifica sem variacao temporal               *
* gradVel   -> gradiente rescontruido das velocidades               *
* lDensity  -> massa especifica sem variacao temporal               *
* lDviscosity-> viscosidade dinamica com variacao temporal          *
* nEn       -> numero de nos da celula central                      *
* nFace     -> numero de faces da celula central                    *
* ndm       -> numero de dimensoes                                  *
* nel       -> numero da celula                                     *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* eddyViscosity                                                     *
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*                                                                   *
* gradVel(nFace+1,ndf,ndm) -> gradVel(nFace+1,2,2)                  *
*                                                                   *
*                  | du1dx1 du1dx2 |                                *
* grad(viz1,*,*) = |               |                                *
*                  | du2dx1 du2dx2 |                                *
*                                                                   *
*                  | du1dx1 du1dx2 |                                *
* grad(viz2,*,*) = |               |                                *
*                  | du2dx1 du2dx2 |                                *
*                                                                   *
*********************************************************************/
void cellLes(Loads *lVel               , Turbulence tModel           
          , short *RESTRICT lGeomType  , DOUBLE *RESTRICT prop 
          , INT *RESTRICT lViz         , DOUBLE *RESTRICT ksi 
          , DOUBLE *RESTRICT mKsi
          , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT mEta 
          , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT area 
          , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc 
          , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT cc 
          , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew 
          , short *RESTRICT lFaceVelR  , short *RESTRICT lFaceVelL 
          , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT gradVel       
          , DOUBLE *RESTRICT lDensity  , DOUBLE const dViscosity
          , DOUBLE *viscosity          
          , const short nEn            , short const nFace 
          , const short ndm            , INT const nel) 
{
/*...*/
  bool fWall = false;
  short nAresta, nCarg, idCell = nFace,wallType;
  INT vizNel;
  DOUBLE modS, tmp, densityC, viscosityC, cs, s[3], gradVelC[2][2],delta;
  DOUBLE wt,velC[2],vParallel[2],lNormal[2],lMin,dMin,xx[4]
        ,yPlus,uPlus,velB[2],yPlusMax,par[MAXLOADPARAMETER];

/*...*/
  cs         = tModel.cs;
  wallType   = tModel.wallType;
  densityC   = lDensity[idCell];
  viscosityC = dViscosity;
/*...................................................................*/

/*...*/
  velC[0] = MAT2D(idCell, 0, vel, ndm);
  velC[1] = MAT2D(idCell, 1, vel, ndm);
/*...................................................................*/

/*... | du1/dx1 du1/dx2*/
  gradVelC[0][0] = MAT3D(idCell, 0, 0, gradVel, 2, ndm);
  gradVelC[0][1] = MAT3D(idCell, 0, 1, gradVel, 2, ndm);
/*... | du2/dx1 du2/dx2*/
  gradVelC[1][0] = MAT3D(idCell, 1, 0, gradVel, 2, ndm);
  gradVelC[1][1] = MAT3D(idCell, 1, 1, gradVel, 2, ndm);
/*...................................................................*/
  
/*... wall model*/
  yPlusMax = 0.e0;
  dMin     = 1.e+16;
  for (nAresta = 0; nAresta<nFace; nAresta++) {
    vizNel         = lViz[nAresta]; 
    yPlus = 0.e0; 
/*... contorno*/
    if (vizNel  == -2) {
      if (lFaceVelR[nAresta] > 0) {
        nCarg = lFaceVelL[nAresta] - 1;
        if( lVel[nCarg].type == MOVEWALL){
          fWall = true;
/*... velocidade da parede*/
          nCarg     = lFaceVelL[nAresta] - 1;
          getLoads(par,&loadsVel[nCarg],xx);
          velB[0]   = par[0];
          velB[1]   = par[1];
/*...*/
          lNormal[0] = MAT2D(nAresta, 0, normal, ndm);
          lNormal[1] = MAT2D(nAresta, 1, normal, ndm);   
/*... calculo da velociade paralela a face*/
          wt = velC[0] * lNormal[0] + velC[1] * lNormal[1];
          vParallel[0] = velC[0] - wt * lNormal[0] - velB[0];
          vParallel[1] = velC[1] - wt * lNormal[1] - velB[1];
          wt = vParallel[0]*vParallel[0]+ vParallel[1]*vParallel[1];
          wt = sqrt(wt);
/*...*/
          if (wt > 0.e0)
            wallModel(wt      , viscosityC
                            , densityC, dcca[nAresta]
                            , &yPlus  , &uPlus
                            , wallType); 
/*...................................................................*/
          yPlusMax = max(yPlus,yPlusMax);
          dMin     = min(dcca[nAresta],dMin);
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...*/
      else if (lFaceVelR[nAresta] == STATICWALL){
        fWall = true;
        lNormal[0] = MAT2D(nAresta, 0, normal, ndm);
        lNormal[1] = MAT2D(nAresta, 1, normal, ndm);   
/*... calculo da velociade paralela a face*/
        wt = velC[0] * lNormal[0] + velC[1] * lNormal[1];
        vParallel[0] = velC[0] - wt * lNormal[0];
        vParallel[1] = velC[1] - wt * lNormal[1];
        wt = vParallel[0]*vParallel[0]+ vParallel[1]*vParallel[1];
        wt = sqrt(wt);
/*...*/
        if (wt > 0.e0)
          wallModel(wt      , viscosityC
                  , densityC, dcca[nAresta]
                  , &yPlus  , &uPlus
                  , wallType); 
/*...................................................................*/
        yPlusMax = max(yPlus,yPlusMax);
        dMin     = min(dcca[nAresta],dMin);
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... |S|*/
  tmp  = gradVelC[0][0] + gradVelC[1][1];
  s[0] = gradVelC[0][0] - D1DIV3*tmp;
  s[1] = gradVelC[1][1] - D1DIV3*tmp;
  s[2] = 0.5e0*(gradVelC[0][1] + gradVelC[1][0]);

  modS = s[0]*s[0] + s[1]*s[1] + 2.e0*( s[2]*s[2] );
  modS = sqrt(2.e0*modS);
/*...................................................................*/
 
/*...*/
  delta = sqrt(area[idCell]);
  lMin = cs*delta;
  if (fWall) {
    tmp  = 1.e0-exp(-yPlusMax/VANDRIEST); 
    lMin = min(VONKARMAN*dMin,tmp*lMin); 
  }
/*...................................................................*/
  *viscosity = densityC*lMin*lMin*modS;  
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 03/10/2017                                   *
 * Data de modificaco : 01/04/2018                                   *
 *-------------------------------------------------------------------*
 * eddyViscosity3D: Calculo da viscosidae turbulenta para o LES      *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * lVel      -> definicoes de cargas de velocidades                  *
 * tModel    -> config do modelo de turbulencia                      *
 * geomType  -> tipo geometrico das celulas                          * 
 * prop      -> propriedades dos material                            *  
 * lViz      -> viznhos da celula central                            *
 * fArea     -> area da face                                         *
 * normal    -> vetores normais as faces das celulas                 *
 * volume    -> volume da celula                                     *
 * dcca      -> menor distancia do centroide central a faces desta   *
 *              celula                                               *
 * faceVelR  -> restricoes por elemento de velocidades               *
 * faceVelL  -> carga por elemento de velocidades                    * 
 * vel       -> velocidades fica sem variacao temporal               *
 * gradVel   -> gradiente de velocidade                              *
 * density   -> massa especifica                                     *
 * dViscosity-> viscosidade dinamica                                 *
 * eddyViscosity - nao definido                                      *
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * eddyViscosity -> viscosidade turbulenta                           * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *                                                                   *  
 * gradVel(nFace+1,ndf,ndm) -> gradVel(nFace+1,2,2)                  * 
 *                                                                   *  
 *                  | du1dx1 du1dx2 du1dx3 |                         *
 * grad(viz1,*,*) = | du2dx1 du2dx2 du2dx3 |                         *
 *                  | du3dx1 du3dx2 du3dx3 |                         *
 *                                                                   *  
 *                  | du1dx1 du1dx2 du1dx3 |                         *
 * grad(viz2,*,*) = | du2dx1 du2dx2 du2dx3 |                         *
 *                  | du3dx1 du3dx2 du3dx3 |                         *
 *                                                                   *  
 *                  | du1dx1 du1dx2 du1dx3 |                         *
 * grad(P   ,*,*) = | du2dx1 du2dx2 du2dx3 |                         *
 *                  | du3dx1 du3dx2 du3dx3 |                         *
 *                                                                   *  
 *********************************************************************/
void eddyViscosity3D(Loads *lVel             , Turbulence tModel           
          , short *RESTRICT lGeomType  , DOUBLE *RESTRICT prop 
          , INT *RESTRICT lViz         , DOUBLE *RESTRICT fArea
          , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume 
          , DOUBLE *RESTRICT dcca      , short *RESTRICT lFaceVelR  
          , short *RESTRICT lFaceVelL  , DOUBLE *RESTRICT vel       
          , DOUBLE *RESTRICT gradVel   , DOUBLE *RESTRICT lDensity  
          , DOUBLE const dViscosity    , DOUBLE *viscosity          
          , DOUBLE *RESTRICT wallPar   , DOUBLE const cDyn
          , const short nEn            , short const nFace 
          , const short ndm            , INT const nel) 
{
/*...*/
  bool fDynamic;
  short idCell,type,typeLes;
  DOUBLE op, tmp, density, cf2, s[6], delta,*iGradVel;
  DOUBLE lMin, yPlusMax, c ;

  idCell = nFace;
 
  type     = tModel.typeMixed[FUNMODEL];
  typeLes  = tModel.typeLes;
  fDynamic = tModel.fDynamic;
  c        = tModel.c; 

/*...*/
  yPlusMax = wallPar[0];
/*...................................................................*/

/*...*/
  if (fDynamic)
    cf2 = cDyn;
  else     
    cf2 = tModel.cf*tModel.cf;
/*...................................................................*/

/*...*/
  delta    = pow(volume[idCell],D2DIV3);
  density  = lDensity[idCell];
/*...................................................................*/

/*... (cs*delta)^2*/
  lMin= cf2*delta;
/*...................................................................*/

/*...*/
  switch (type) {
    case SMAGORINSKY:
/*.. calculo Sij*/
      iGradVel = &MAT3D(idCell,0,0,gradVel,3,3);
      tensorS(s,iGradVel,false);
/*... |S| = sqrt(2S:S)*/
      op = sqrt(2.e0*doubleDotSym(s));
/*...................................................................*/
 
/*...*/
      if (!fDynamic) {
        tmp  = 1.e0-exp(-yPlusMax/VANDRIEST);
        lMin = tmp*lMin; 
        lMin = lMin*lMin;
      }
/*...................................................................*/
      break;
/*...................................................................*/

/*...*/
    case WALEMODEL: 
/*... operador(grad(u))*/
      iGradVel = &MAT3D(idCell,0,0,gradVel,3,3);
      op = waleModel(iGradVel);
/*...................................................................*/     
      break;
/*...................................................................*/

/*...*/
    case VREMAN: 
/*... operador(grad(u))*/
      iGradVel = &MAT3D(idCell,0,0,gradVel,3,3);
      op = vremanModel(iGradVel); 
/*...................................................................*/
      break;
/*...................................................................*/

/*...*/
    case SIGMAMODEL: 
/*... operador(grad(u))*/
      iGradVel = &MAT3D(idCell,0,0,gradVel,3,3);
      op = sigmaModel(iGradVel,nFace,ndm);
/*...................................................................*/
      break;
/*...................................................................*/

/*...*/
    default: 
      ERRO_OP(__FILE__,__func__,type);
/*...................................................................*/
  }
/*...................................................................*/

/*...*/   
  *viscosity = density*lMin*op; 
/*...................................................................*/

/*... metodo misto com coeficentes constantes*/
  if (typeLes == LESMIXEDMODEL) 
    *viscosity *= c;
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 24/10/2017                                   *
 * Data de modificaco : 19/12/2017                                   *
 *-------------------------------------------------------------------*
 * lesDynOnePar: Metodo dinamico para o calculo da constente do LES  *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * lViz       -> vizinhos da celula central                          *
 * volume     -> volume da celula central                            *
 * lDensity   -> massa especifica sem variacao temporal              *
 * vel        -> campo de velocidades                                *
 * gradVel    -> gradiente rescontruido das velocidades              *
 * lDynamic   -> nao definido                                        *
 * typeLesFunc-> metodo funcional                                    *
 * nFace      -> numero de faces da celula central                   * 
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lDynamic  -> retorna o produto M:L e M:M da celula central        *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * delta a largura do filtro que é igual a raiz cubica do volume da  *
 * celula central                                                    *
 *-------------------------------------------------------------------*
 *********************************************************************/
void lesDynOnePar(INT *RESTRICT lViz       , DOUBLE *RESTRICT volume
                , DOUBLE *RESTRICT lDensity, DOUBLE *RESTRICT vel
                , DOUBLE *RESTRICT gradVel , DOUBLE *RESTRICT lDynamic
                , short const typeLesFunc  , short const nFace  ) {

  short i,j,k,idCell;
  INT vizNel;
  DOUBLE tmp,s[6],ss[6];
  DOUBLE v[3],m[6],l[6],deltaT,delta,volW,volTotal,density,g[3][3];
  DOUBLE tFilterDenVv[6],rFiltergVel[3][3];
  DOUBLE tFilterDen,tFilterDenV[3],op=0.e0;
  DOUBLE tFilterFf[6];
  DOUBLE *iGradVel=NULL;

  idCell = nFace;

  tFilterDenV[0] = tFilterDenV[1] = tFilterDenV[2] = 0.e0;
  tFilterDen     = volTotal = 0.e0;
  
  iGradVel = rFiltergVel[0];
  for(i = 0; i< 9 ; i++ )
    iGradVel[i] = 0.e0;

  for(i = 0; i<6 ; i++ ){
    tFilterFf[i]       = 0.e0;  
    tFilterDenVv[i]    = 0.e0;
  }

/*...*/
  for(i = 0; i<nFace+1;i++ ){
    vizNel = 0;
    if( i < nFace )
      vizNel = lViz[i];
/*... celula central e elementos vizinhos*/        
    if( vizNel  > -1 ){
      volW      = volume[i];
      volTotal += volW;

/*...*/
      v[0] = MAT2D(i, 0, vel, 3);
      v[1] = MAT2D(i, 1, vel, 3);
      v[2] = MAT2D(i, 2, vel, 3);
/*...................................................................*/

/*...*/
      density = lDensity[i];
/*...................................................................*/

/*... | du1/dx1 du1/dx2 du1/dx3*/
      g[0][0] = MAT3D(i,0,0,gradVel,3,3);
      g[0][1] = MAT3D(i,0,1,gradVel,3,3);
      g[0][2] = MAT3D(i,0,2,gradVel,3,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
      g[1][0] = MAT3D(i,1,0,gradVel,3,3);
      g[1][1] = MAT3D(i,1,1,gradVel,3,3);
      g[1][2] = MAT3D(i,1,2,gradVel,3,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
      g[2][0] = MAT3D(i,2,0,gradVel,3,3);
      g[2][1] = MAT3D(i,2,1,gradVel,3,3);
      g[2][2] = MAT3D(i,2,2,gradVel,3,3);
/*...................................................................*/

/*... tensor L*/

/*... tesFilter( gradVel ) -> vol*gradVel*/      
      for(j=0;j<3;j++)
        for(k=0;k<3;k++)
          rFiltergVel[j][k] += volW*g[j][k];

/*... tesFilter(desinty*vv) -> vol*desinty*vv */
      tmp              = volW*density;
      tFilterDenVv[0] += tmp*v[0]*v[0]; /*v1v1*/
      tFilterDenVv[1] += tmp*v[1]*v[1]; /*v2v2*/
      tFilterDenVv[2] += tmp*v[2]*v[2]; /*v3v3*/
      tFilterDenVv[3] += tmp*v[0]*v[1]; /*v1v2*/
      tFilterDenVv[4] += tmp*v[0]*v[2]; /*v1v3*/
      tFilterDenVv[5] += tmp*v[1]*v[2]; /*v2v3*/
/*... tesFilter(desinty*v) -> vol*desinty*v */
      tFilterDenV[0] += tmp*v[0]; 
      tFilterDenV[1] += tmp*v[1]; 
      tFilterDenV[2] += tmp*v[2];
/*... tesFilter(desinty) -> vol*desinty */
      tFilterDen += tmp;
/*...................................................................*/

/*... smagorinsky*/
/*.. calculo Sij*/
      tensorS(s,g[0],true);
/*...*/
      if( typeLesFunc == SMAGORINSKY ){
/*... |S| = |2S:S|*/
        tensorS(ss,g[0],false);
        op = sqrt(2.e0*doubleDotSym(ss));
      }
/*...................................................................*/

/*... Sigmamodel*/
      else if (typeLesFunc == SIGMAMODEL) {
/*... OP = sg3*(sg1 -sg2)*(sg2 -sg3)/sg1^2 )*/
        iGradVel = &MAT3D(i,0,0,gradVel,3,3);
        op = sigmaModel(iGradVel,nFace ,3);
      }
/*...................................................................*/

/*... wale*/
      else if (typeLesFunc == WALEMODEL) {
/*... OP = ( (Sd:Sd)^1.5 ) / ( (S:S)^2.5 + (Sd:Sd)^1.25 )*/
        iGradVel = &MAT3D(i,0,0,gradVel,3,3);
        op = waleModel(iGradVel);
      }
/*...................................................................*/

/*... vreman*/
      else if (typeLesFunc == VREMAN) {
/*... OP = sqrt(B/g:g)*/
        iGradVel = &MAT3D(i,0,0,gradVel,3,3);
        op = vremanModel(iGradVel);
      }
/*...................................................................*/

      tmp = volW*op*density;
      for(j=0;j<6;j++)
/*... tesFilter( desinty*OP*s) -> vol*desinty*OP*s */
        tFilterFf[j] += tmp*s[j];  
/*...................................................................*/
    }
/*...................................................................*/
  } 
/*...................................................................*/

/*... (filtros)^2*/
  delta    = pow(volume[idCell],D2DIV3);
  deltaT   = pow(volTotal      ,D2DIV3);
/*...................................................................*/

/*... tensor L*/

/*... tesFilter( gradVel ) */
  iGradVel = rFiltergVel[0];
  for(i = 0; i< 9 ; i++ )
    iGradVel[i] /= volTotal;

/*... tesFilter(density*vv)*/
  for(i=0;i<6;i++)
    tFilterDenVv[i] /= volTotal; /*v1v1*/
 
/*... tesFilter(density*v)*/
  tFilterDenV[0] /= volTotal; 
  tFilterDenV[1] /= volTotal; 
  tFilterDenV[2] /= volTotal;
/*... tesFilter(desinty)*/
  tFilterDen /= volTotal;
      
/*... L*/
  l[0] = tFilterDenVv[0] - tFilterDenV[0]*tFilterDenV[0]/tFilterDen;  /*l11*/
  l[1] = tFilterDenVv[1] - tFilterDenV[1]*tFilterDenV[1]/tFilterDen;  /*l22*/
  l[2] = tFilterDenVv[2] - tFilterDenV[2]*tFilterDenV[2]/tFilterDen;  /*l33*/
  l[3] = tFilterDenVv[3] - tFilterDenV[0]*tFilterDenV[1]/tFilterDen;  /*l1v2*/
  l[4] = tFilterDenVv[4] - tFilterDenV[0]*tFilterDenV[2]/tFilterDen;  /*l1v3*/
  l[5] = tFilterDenVv[5] - tFilterDenV[1]*tFilterDenV[2]/tFilterDen;  /*l2v3*/
/*... parte desviadora de L*/
  tmp   = D1DIV3*( l[0] + l[1] + l[2] );
  l[0] -=tmp;
  l[1] -=tmp;
  l[2] -=tmp;

/*.... tesFilter(s)*/
/*... Sij - (1/3)*Skk*deltaij*/
  tensorS(s,rFiltergVel[0],true);
/*... samgorisnky*/
  if( typeLesFunc == SMAGORINSKY ){ 
/*... |S| = |2S:S|*/
    tensorS(ss,rFiltergVel[0],false);
    op = sqrt(2.e0*doubleDotSym(ss));
  }
/*...................................................................*/

/*... Sigmamodel*/
  else if (typeLesFunc == SIGMAMODEL) {
/*... OP=sg3*(sg1 -sg2)*(sg2 -sg3)/sg1^2 )*/
    iGradVel = rFiltergVel[0];
    op = sigmaModel(iGradVel, nFace  , 3);
  }
/*...................................................................*/

/*... wale*/
  else if (typeLesFunc == WALEMODEL) {
/*... OP= ( (Sd:Sd)^1.5 ) / ( (S:S)^2.5 + (Sd:Sd)^1.25 )*/
    iGradVel = rFiltergVel[0];
    op = waleModel(iGradVel);
  }
/*...................................................................*/
  
/*... vreman*/
  else if (typeLesFunc == VREMAN) {
/*... OP = sqrt(B/g:g)*/
    iGradVel = rFiltergVel[0];
    op = vremanModel(iGradVel);
  }
/*...................................................................*/

/*...*/
  tmp = deltaT*tFilterDen*op;
  for(i=0;i<6;i++){
/*.... tesFilter( density*|s|*s ) */
    tFilterFf[i] /= volTotal;
/*...*/
    m[i] = 2.e0*(delta*tFilterFf[i] - tmp*s[i]);  
  }
/*...................................................................*/

/*... LijMij*/
  lDynamic[0] = doubleDotSym2(l,m);
/*... MijMij*/
  lDynamic[1] = doubleDotSym(m);      
/*...................................................................*/

/*... SVS*/
  lDynamic[2] = 0.0e0; 
  if( typeLesFunc == WALEMODEL ){
    iGradVel = &MAT3D(idCell,0,0,gradVel,3,3);    
    lDynamic[2] = svs(iGradVel);
  }
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 21/01/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * lesDynOneKeq: Metodo dinamico para o calculo da constente do LES  *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * lViz       -> vizinhos da celula central                          *
 * volume     -> volume da celula central                            *
 * lDensity   -> massa especifica sem variacao temporal              *
 * vel        -> campo de velocidades                                *
 * gradVel    -> gradiente rescontruido das velocidades              *
 * lDynamic   -> nao definido                                        *
 * typeLesFunc-> metodo funcional                                    *
 * nFace      -> numero de faces da celula central                   * 
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lDynamic  -> retorna o produto M:L e M:M da celula central        *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * delta a largura do filtro que é igual a raiz cubica do volume da  *
 * celula central                                                    *
 *-------------------------------------------------------------------*
 *********************************************************************/
void lesOneKeqLoc(INT *RESTRICT lViz       , DOUBLE *RESTRICT volume
                , DOUBLE *RESTRICT lDensity, DOUBLE *RESTRICT vel
                , DOUBLE *RESTRICT gradVel , DOUBLE *RESTRICT lDynamic
                , DOUBLE const dViscosity  , DOUBLE const eddyViscosity
                , INT const nel            , short const nFace  ) {

  short i,j,k;
  INT vizNel;
  DOUBLE tmp,tmp1,s[6];
  DOUBLE v[3],m[6],l[6],deltaT,volW,volTotal,density,g[3][3];
  DOUBLE tFilterDenVv[6],rFiltergVel[3][3];
  DOUBLE tFilterDen,tFilterDenV[3],kTest;
  DOUBLE tFilterVv[3],tFilterV[3],tFilterGvGv[3][3];
  DOUBLE *p1=NULL,*p2=NULL;

  tFilterDenV[0] = tFilterDenV[1] = tFilterDenV[2] = 0.e0;
  tFilterVv[0] = tFilterVv[1] = tFilterVv[2] = 0.e0;
  tFilterV[0]   = tFilterV[1] = tFilterV[2] = 0.e0;
  tFilterDen     = volTotal = 0.e0;
  
  p1 = rFiltergVel[0];
  p2 = tFilterGvGv[0];
  for(i = 0; i< 9 ; i++ ){
    p1[i] = 0.e0;
    p2[i] = 0.e0;
  }

  for(i = 0; i<6 ; i++ ){
    tFilterDenVv[i] = 0.e0;    
  }

/*...*/
  for(i = 0; i<nFace+1;i++ ){
    vizNel = 0;
    if( i < nFace )
      vizNel = lViz[i];
/*... celula central e elementos vizinhos*/        
    if( vizNel  > -1 ){
      volW      = volume[i];
      volTotal += volW;

/*...*/
      v[0] = MAT2D(i, 0, vel, 3);
      v[1] = MAT2D(i, 1, vel, 3);
      v[2] = MAT2D(i, 2, vel, 3);
/*...................................................................*/

/*...*/
      density = lDensity[i];
/*...................................................................*/

/*... | du1/dx1 du1/dx2 du1/dx3*/
      g[0][0] = MAT3D(i,0,0,gradVel,3,3);
      g[0][1] = MAT3D(i,0,1,gradVel,3,3);
      g[0][2] = MAT3D(i,0,2,gradVel,3,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
      g[1][0] = MAT3D(i,1,0,gradVel,3,3);
      g[1][1] = MAT3D(i,1,1,gradVel,3,3);
      g[1][2] = MAT3D(i,1,2,gradVel,3,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
      g[2][0] = MAT3D(i,2,0,gradVel,3,3);
      g[2][1] = MAT3D(i,2,1,gradVel,3,3);
      g[2][2] = MAT3D(i,2,2,gradVel,3,3);
/*...................................................................*/

/*... tesFilter( gradVel ) -> vol*gradVel*/      
      for(j=0;j<3;j++)
        for(k=0;k<3;k++)
          rFiltergVel[j][k] += volW*g[j][k];

/*... tesFilter(desinty*vv) -> vol*desinty*vv */
      tmp              = volW*density;
      tFilterDenVv[0] += tmp*v[0]*v[0]; /*v1v1*/
      tFilterDenVv[1] += tmp*v[1]*v[1]; /*v2v2*/
      tFilterDenVv[2] += tmp*v[2]*v[2]; /*v3v3*/
      tFilterDenVv[3] += tmp*v[0]*v[1]; /*v1v2*/
      tFilterDenVv[4] += tmp*v[0]*v[2]; /*v1v3*/
      tFilterDenVv[5] += tmp*v[1]*v[2]; /*v2v3*/
/*... tesFilter(desinty*v) -> vol*desinty*v */
      tFilterDenV[0] += tmp*v[0]; 
      tFilterDenV[1] += tmp*v[1]; 
      tFilterDenV[2] += tmp*v[2];
/*... tesFilter(desinty) -> vol*desinty */
      tFilterDen += tmp;
/*...................................................................*/

/*... tesFilter(vv) -> vol*vv */
      tFilterVv[0] += volW*v[0]*v[0]; /*v1v1*/
      tFilterVv[1] += volW*v[1]*v[1]; /*v2v2*/
      tFilterVv[2] += volW*v[2]*v[2]; /*v3v3*/ 
/*... tesFilter(v) -> vol*v */
      tFilterV[0] += volW*v[0]; 
      tFilterV[1] += volW*v[1]; 
      tFilterV[2] += volW*v[2];
/*...................................................................*/

/*... tesFilter(gVgV) -> vol*GvGv */
      for(j=0;j<3;j++)
        for(k=0;k<3;k++)
          tFilterGvGv[j][k] += volW*g[j][k]*g[j][k];
    }
/*...................................................................*/
  } 
/*...................................................................*/

/*... (filtros)^2*/
  deltaT   = pow(volTotal      ,D1DIV3);
/*...................................................................*/

/*...*/
  p1 = rFiltergVel[0];
  p2 = tFilterGvGv[0];
  for(i = 0; i< 9 ; i++ ){
/*... tesFilter( gradVel ) */
    p1[i] /= volTotal;
/*... tesFilter( gVgV ) */
    p2[i] /= volTotal;
  }
 
/*... tesFilter(density*vv)*/
  for(i=0;i<6;i++){
    tFilterDenVv[i] /= volTotal; 
  }

  for(i=0;i<3;i++){
/*... tesFilter(density*v)*/
    tFilterDenV[i] /= volTotal; 
/*... tesFilter(Vv)*/
    tFilterVv[i] /= volTotal;
/*... tesFilter(v)*/ 
    tFilterV[i] /= volTotal;
  }
/*... tesFilter(desinty)*/
  tFilterDen /= volTotal;
      
/*... L*/
  l[0] = tFilterDenVv[0] - tFilterDenV[0]*tFilterDenV[0]/tFilterDen;  /*l11*/
  l[1] = tFilterDenVv[1] - tFilterDenV[1]*tFilterDenV[1]/tFilterDen;  /*l22*/
  l[2] = tFilterDenVv[2] - tFilterDenV[2]*tFilterDenV[2]/tFilterDen;  /*l33*/
  l[3] = tFilterDenVv[3] - tFilterDenV[0]*tFilterDenV[1]/tFilterDen;  /*l1v2*/
  l[4] = tFilterDenVv[4] - tFilterDenV[0]*tFilterDenV[2]/tFilterDen;  /*l1v3*/
  l[5] = tFilterDenVv[5] - tFilterDenV[1]*tFilterDenV[2]/tFilterDen;  /*l2v3*/
/*... parte desviadora de L*/
  tmp   = D1DIV3*( l[0] + l[1] + l[2] );
  l[0] -=tmp;
  l[1] -=tmp;
  l[2] -=tmp;

/*.... tesFilter(s)*/
/*... Sij - (1/3)*Skk*deltaij*/
  tensorS(s,rFiltergVel[0],true);
 /*... */
  kTest = 0.e0;
  for(i=0;i<3;i++)    
    kTest += tFilterVv[i] - tFilterV[i]*tFilterV[i];
  kTest *= 0.5e0;

  if(kTest <= 0.e0)
    kTest = 1.e-14;
                            
/*...*/
  tmp = 2.e0*deltaT*sqrt(kTest);
  for(i=0;i<6;i++){
    m[i] = -tmp*s[i];  
  }
/*...................................................................*/

/*... LijLij*/
  lDynamic[0] = doubleDotSym2(l,m);
/*... MijMij*/
  lDynamic[1] = doubleDotSym(m);      
/*...................................................................*/

/*... Ce*/
  p1 = rFiltergVel[0];
  p2 = tFilterGvGv[0];
  tmp1 = 0.e0;
  for(i = 0; i< 9 ; i++ )
    tmp1 += p2[i] - p1[i]*p1[i];

  lDynamic[2] = (dViscosity + eddyViscosity)*tmp1; 
  lDynamic[3] = pow(kTest,1.5)/deltaT;
/*...................................................................*/

}
/*********************************************************************/


/*********************************************************************
 * Data de criacao    : 13/01/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * lesDynOneParG: Metodo dinamico global para o calculo da constente *
 * do LES                                                            *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * lViz       -> vizinhos da celula central                          *
 * volume     -> volume da celula central                            *
 * gradVel    -> gradiente rescontruido das velocidades              *
 * lDynamic   -> nao definido                                        *
 * typeLesFunc-> metodo funcional                                    *
 * nFace      -> numero de faces da celula central                   * 
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lDynamic  -> retorna o produto:                                   *
 *             testFilter(Al:Al)                                     *
 *             testFilter(Al):testFilter(Al)                         *
 *             op*delta*tesFilter(S:S)                               *
 *             opT*deltaT*tesFilter(S):tesFilter(S)                  *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * delta a largura do filtro que é igual a raiz cubica do volume da  *
 * celula central                                                    *
 *-------------------------------------------------------------------*
 *********************************************************************/
void lesDynOneParG(INT *RESTRICT lViz      , DOUBLE *RESTRICT volume
                , DOUBLE *RESTRICT gradVel , DOUBLE *RESTRICT lDynamic
                , short const typeLesFunc  , short const nFace  ) {

  short i,j,k,idCell;
  INT vizNel;
  DOUBLE deltaT,delta,volW,volTotal,op,opt;
  DOUBLE rFilterAl[3][3],rFilterAlAl[3][3],g[3][3];
  DOUBLE tFilterS[6],tFilterSs[6],s[6];
  DOUBLE *p1=NULL,*p2=NULL;

  idCell = nFace;
  volTotal = 0.e0;
  op = opt = 0.e0;
  
  p1 = rFilterAl[0];
  p2 = rFilterAlAl[0];
  for(i = 0; i< 9 ; i++ ){
    p1[i] = 0.e0;
    p2[i] = 0.e0;
  }

  for(i = 0; i<6 ; i++ ){
    tFilterS[i]  = 0.e0;  
    tFilterSs[i] = 0.e0; 
  }

/*...*/
  for(i = 0; i<nFace+1;i++ ){
    vizNel = 0;
    if( i < nFace )
      vizNel = lViz[i];
/*... celula central e elementos vizinhos*/        
    if( vizNel  > -1 ){
      volW      = volume[i];
      volTotal += volW;

/*... | du1/dx1 du1/dx2 du1/dx3*/
      g[0][0] = MAT3D(i,0,0,gradVel,3,3);
      g[0][1] = MAT3D(i,0,1,gradVel,3,3);
      g[0][2] = MAT3D(i,0,2,gradVel,3,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
      g[1][0] = MAT3D(i,1,0,gradVel,3,3);
      g[1][1] = MAT3D(i,1,1,gradVel,3,3);
      g[1][2] = MAT3D(i,1,2,gradVel,3,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
      g[2][0] = MAT3D(i,2,0,gradVel,3,3);
      g[2][1] = MAT3D(i,2,1,gradVel,3,3);
      g[2][2] = MAT3D(i,2,2,gradVel,3,3);
/*...................................................................*/

/*...*/      
      for(j=0;j<3;j++){
        for(k=0;k<3;k++){
/*... tesFilter( gradVel ) -> vol*gradVel*/      
          rFilterAl[j][k]   += volW*g[j][k];
/*... tesFilter( gradVel * gradVel) -> vol*gradVel*gradVel*/
          rFilterAlAl[j][k] += volW*g[j][k]*g[j][k];
        }  
      }
/*...................................................................*/

/*.. calculo Sij*/
      tensorS(s,g[0],true);
      for(j=0;j<6;j++){
/*... tesFilter( s ) -> vol*s */
        tFilterS[j] += volW*s[j];  
/*... tesFilter( sijsij ) -> vol*sij*sij */
        tFilterSs[j] += volW*s[j]*s[j];
      }
/*...................................................................*/
    }
/*...................................................................*/
  } 
/*...................................................................*/

/*... (filtros)^2*/
  delta    = pow(volume[idCell],D2DIV3);
  deltaT   = pow(volTotal      ,D2DIV3);
/*...................................................................*/

/*...*/
  p1 = rFilterAl[0];
  p2 = rFilterAlAl[0];
  for(i = 0; i< 9 ; i++ ){
/*... tesFilter( alpha )*/
    p1[i] /= volTotal;
/*... tesFilter( alpha*alpha )*/
    p2[i] /= volTotal;
  }
/*...................................................................*/
    
/*...*/
  for(i = 0; i< 6 ; i++ ){
/*... tesFilter( alpha )*/
    tFilterS[i]  /= volTotal;
/*... tesFilter( alpha*alpha )*/
    tFilterSs[i] /= volTotal;
  }
/*...................................................................*/

/*... | du1/dx1 du1/dx2 du1/dx3*/
  g[0][0] = MAT3D(idCell,0,0,gradVel,3,3);
  g[0][1] = MAT3D(idCell,0,1,gradVel,3,3);
  g[0][2] = MAT3D(idCell,0,2,gradVel,3,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
  g[1][0] = MAT3D(idCell,1,0,gradVel,3,3);
  g[1][1] = MAT3D(idCell,1,1,gradVel,3,3);
  g[1][2] = MAT3D(idCell,1,2,gradVel,3,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
  g[2][0] = MAT3D(idCell,2,0,gradVel,3,3);
  g[2][1] = MAT3D(idCell,2,1,gradVel,3,3);
  g[2][2] = MAT3D(idCell,2,2,gradVel,3,3);
/*...................................................................*/
 
/*....*/
/*... Sij - (1/3)*Skk*deltaij*/
  tensorS(s,rFilterAl[0],true);
/*... Sigmamodel*/
  if (typeLesFunc == SIGMAMODEL) {
/*... OP=sg3*(sg1 -sg2)*(sg2 -sg3)/sg1^2 )*/
    p1 = rFilterAl[0];
    p2 = g[0];
    opt = sigmaModel(p1, nFace  , 3);
    op  = sigmaModel(p2, nFace  , 3);
  }
/*...................................................................*/

/*... wale*/
  else if (typeLesFunc == WALEMODEL) {
/*... OP= ( (Sd:Sd)^1.5 ) / ( (S:S)^2.5 + (Sd:Sd)^1.25 )*/
    p1 = rFilterAl[0];
    p2 = g[0];
    opt = waleModel(p1);
    op  = waleModel(p2);
  }
/*...................................................................*/
  
/*... vreman*/
  else if (typeLesFunc == VREMAN) {
/*... OP = sqrt(B/g:g)*/
    p1 = rFilterAl[0];
    p2 = g[0];
    opt = vremanModel(p1);
    op  = vremanModel(p2);
  }
/*...................................................................*/

/*... tesFilter(Al:Al)*/
  lDynamic[0] = rFilterAlAl[0][0] + rFilterAlAl[0][1] + rFilterAlAl[0][2]
              + rFilterAlAl[1][0] + rFilterAlAl[1][1] + rFilterAlAl[1][2]
              + rFilterAlAl[2][0] + rFilterAlAl[2][1] + rFilterAlAl[2][2];
/*... tesFilter(Al):tesFilter(Al)*/
  lDynamic[1] = doubleDot(rFilterAl[0]);      
/*... deltaT*opt*tesFilter(S:S)*/
  lDynamic[2] = tFilterSs[0] + tFilterSs[1] + tFilterSs[2] 
              + 2.e0*(tFilterSs[3] + tFilterSs[4] + tFilterSs[5]);  
  lDynamic[2] = delta*op*lDynamic[2];
/*... delta*op*tesFilter(S):tesFilter(S)*/
  lDynamic[3] = deltaT*opt*doubleDotSym(tFilterS);      
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 24/10/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * lesDynTwoPar: Metodo dinamico para o calculo da metodos mistos de *
 * dois parametros                                                   *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * lViz      -> vizinhos da celula central                           *
 * volume    -> volume da celula central                             *
 * lDensity  -> massa especifica sem variacao temporal               *
 * vel       -> campo de velocidades                                 *
 * gradVel   -> gradiente rescontruido das velocidades               *
 * lDynamic  -> nao definido                                         *
 * nFace     -> numero de faces da celula central                    * 
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lDynamic  -> retorna o produto M:L e M:M da celula central        *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * delta a largura do filtro que é igual a raiz cubica do volume da  *
 * celula central                                                    *
 *-------------------------------------------------------------------*
 *********************************************************************/
void lesDynTwoPar(INT *RESTRICT lViz       , DOUBLE *RESTRICT volume
                , DOUBLE *RESTRICT lDensity, DOUBLE *RESTRICT vel
                , DOUBLE *RESTRICT gradVel , DOUBLE *RESTRICT lDynamic
                , short const nFace  ) {

  short i,j,k,idCell;
  INT vizNel;
  DOUBLE tmp,s[6],modS;
  DOUBLE v[3],ms[6],mf[6],l[6],deltaT,delta,volW,volTotal,density;
  DOUBLE g[3][3];
  DOUBLE fs[6],tFilterFs[6],ff[6],tFilterFf[6],rFiltergVel[3][3];
  DOUBLE tFilterModS,tFilterDen,tFilterDenV[3],tFilterS[6],tFilterDenVv[6];

  idCell = nFace;

  tFilterDenV[0] = tFilterDenV[1] = tFilterDenV[2] = 0.e0;
  tFilterDen     = tFilterModS = volTotal = 0.e0;


  for(i = 0; i< 3 ; i++ )
    for(j = 0; j< 3 ; j++ )
      rFiltergVel[i][j] = 0.e0;

  for(i = 0; i< 6 ; i++ ){
    fs[i]              = 0.e0;
    ff[i]              = 0.e0;
    tFilterFs[i]       = 0.e0;
    tFilterFf[i]       = 0.e0;
    tFilterS[i]        = 0.e0;
    tFilterDenVv[i]    = 0.e0;
  }

/*...*/
  for(i = 0; i<nFace+1;i++ ){
    vizNel = 0;
    if( i < nFace )
      vizNel = lViz[i];
/*... celula central e elementos vizinhos*/        
    if( vizNel  > -1 ){
      volW      = volume[i];
      volTotal += volW;

/*...*/
      v[0] = MAT2D(i, 0, vel, 3);
      v[1] = MAT2D(i, 1, vel, 3);
      v[2] = MAT2D(i, 2, vel, 3);
/*...................................................................*/

/*...*/
      density = lDensity[i];
/*...................................................................*/

/*... | du1/dx1 du1/dx2 du1/dx3*/
      g[0][0] = MAT3D(i,0,0,gradVel,3,3);
      g[0][1] = MAT3D(i,0,1,gradVel,3,3);
      g[0][2] = MAT3D(i,0,2,gradVel,3,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
      g[1][0] = MAT3D(i,1,0,gradVel,3,3);
      g[1][1] = MAT3D(i,1,1,gradVel,3,3);
      g[1][2] = MAT3D(i,1,2,gradVel,3,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
      g[2][0] = MAT3D(i,2,0,gradVel,3,3);
      g[2][1] = MAT3D(i,2,1,gradVel,3,3);
      g[2][2] = MAT3D(i,2,2,gradVel,3,3);
/*...................................................................*/

/*... clark model*/

/*... tesFilter( gradVel ) -> vol*gradVel*/      
      for(k=0;k<3;k++)
        for(j=0;j<3;j++)
          rFiltergVel[k][j] += volW*g[k][j];

/*... tesFilter( StressClarkDesviador ) -> vol*gradVel*/  
      tmp = D1DIV3*doubleDot(g[0]);

      for(k=0;k<3;k++)
        tFilterFs[k] += volW*(g[0][k]*g[0][k] 
                            + g[1][k]*g[1][k] 
                            + g[2][k]*g[2][k] - tmp);

      tFilterFs[3] += volW*(g[0][0]*g[0][1] 
                          + g[1][0]*g[1][1] 
                          + g[2][0]*g[2][1]);

      tFilterFs[4] += volW*(g[0][1]*g[0][2] 
                          + g[1][1]*g[1][2] 
                          + g[2][1]*g[2][2]);

      tFilterFs[5] += volW*(g[0][0]*g[0][2] 
                          + g[1][0]*g[1][2] 
                          + g[2][0]*g[2][2]);

/*... smagorinsky*/

/*.. calculo Sij*/
      tmp  = D1DIV3*(g[0][0] + g[1][1] + g[2][2]);
      s[0] = g[0][0] - tmp;
      s[1] = g[1][1] - tmp;
      s[2] = g[2][2] - tmp;
      s[3] = 0.5e0*(g[0][1] + g[1][0]); /*s12*/
      s[4] = 0.5e0*(g[0][2] + g[2][0]); /*s13*/
      s[5] = 0.5e0*(g[1][2] + g[2][1]); /*s23*/
/*... |S| = |2S:S|*/
      modS = doubleDotSym(s);
      modS = sqrt(2.e0*modS);
/*... tesFilter( desinty*|s|s) -> vol*desinty*modS*s */
      tmp                 = volW*modS*density;
      for(k=0;k<6;k++)
        tFilterFf[k] += tmp*s[k];

/*... tesFilter( |s|) -> vol*modS */
      tFilterModS += volW*modS;
/*... tesFilter(s) -> vol*s */
      for(k=0;k<6;k++)
        tFilterS[k] += volW*s[k];

/*... tensor L*/

/*... tesFilter(desinty*vv) -> vol*desinty*vv */
      tmp              = volW*density;
      tFilterDenVv[0] += tmp*v[0]*v[0]; /*v1v1*/
      tFilterDenVv[1] += tmp*v[1]*v[1]; /*v2v2*/
      tFilterDenVv[2] += tmp*v[2]*v[2]; /*v3v3*/
      tFilterDenVv[3] += tmp*v[0]*v[1]; /*v1v2*/
      tFilterDenVv[4] += tmp*v[0]*v[2]; /*v1v3*/
      tFilterDenVv[5] += tmp*v[1]*v[2]; /*v2v3*/
/*... tesFilter(desinty*v) -> vol*desinty*v */
      tFilterDenV[0] += tmp*v[0]; 
      tFilterDenV[1] += tmp*v[1]; 
      tFilterDenV[2] += tmp*v[2];
/*... tesFilter(desinty) -> vol*desinty */
      tFilterDen += tmp; 
    }
/*...................................................................*/
  } 
/*...................................................................*/

/*... (filtros)^2*/
  delta    = pow(volume[idCell],D2DIV3);
  deltaT   = pow(volTotal      ,D2DIV3);

/*... tensor L*/
/*... tesFilter(density*vv)*/
  for(i=0;i<6;i++)
    tFilterDenVv[i] /= volTotal; 

/*... tesFilter(density*v)*/
  tFilterDenV[0] /= volTotal; 
  tFilterDenV[1] /= volTotal; 
  tFilterDenV[2] /= volTotal;
/*... tesFilter(desinty)*/
  tFilterDen     /= volTotal;

/*... L*/
  l[0] = tFilterDenVv[0] - tFilterDenV[0]*tFilterDenV[0]/tFilterDen;  /*l11*/
  l[1] = tFilterDenVv[1] - tFilterDenV[1]*tFilterDenV[1]/tFilterDen;  /*l22*/
  l[2] = tFilterDenVv[2] - tFilterDenV[2]*tFilterDenV[2]/tFilterDen;  /*l33*/
  l[3] = tFilterDenVv[3] - tFilterDenV[0]*tFilterDenV[1]/tFilterDen;  /*l1v2*/
  l[4] = tFilterDenVv[4] - tFilterDenV[0]*tFilterDenV[2]/tFilterDen;  /*l1v3*/
  l[5] = tFilterDenVv[5] - tFilterDenV[1]*tFilterDenV[2]/tFilterDen;  /*l2v3*/

/*... parte desviadora de L*/
  tmp   = D1DIV3*( l[0] + l[1] + l[2] );
  l[0] -=tmp;
  l[1] -=tmp;
  l[2] -=tmp;
/*...................................................................*/

/*... clark model L*/

/*tesFilter( gradVel ) */
  for(k=0;k<3;k++)
    for(j=0;j<3;j++)
      rFiltergVel[k][j] /= volTotal;

/*tesFilter( (dui/dxk)*(duj/dxk) - 1/3 (dul/dxk)*(dul/dxk) ) */
  for(i=0;i<6;i++)
    tFilterFs[i] /= volTotal;

  tmp = D1DIV3*doubleDot(rFiltergVel[0]);

/*fs*/
  for(k=0;k<3;k++)
    fs[k] = deltaT*(rFiltergVel[0][k]*rFiltergVel[0][k] 
                  + rFiltergVel[1][k]*rFiltergVel[1][k] 
                  + rFiltergVel[2][k]*rFiltergVel[2][k] - tmp);

  fs[3] = deltaT*(rFiltergVel[0][0]*rFiltergVel[0][1] 
                + rFiltergVel[1][0]*rFiltergVel[1][1] 
                + rFiltergVel[2][0]*rFiltergVel[2][1]);

  fs[4] = deltaT*(rFiltergVel[0][1]*rFiltergVel[0][2] 
                + rFiltergVel[1][1]*rFiltergVel[1][2] 
                + rFiltergVel[2][1]*rFiltergVel[2][2]);

  fs[5] = deltaT*(rFiltergVel[0][0]*rFiltergVel[0][2] 
                + rFiltergVel[1][0]*rFiltergVel[1][2] 
                + rFiltergVel[2][0]*rFiltergVel[2][2]);

  for(i=0;i<6;i++)
    ms[i] = (fs[i] - delta*tFilterFs[i])/12.e0;  /*l11*/

/*... samgorisnky*/
 
/*... tesFilter( |s|) -> vol*modS */
  tFilterModS /= volTotal;

  tmp = deltaT*tFilterDen*tFilterModS;
  for(i=0;i<6;i++){
/*tesFilter( density*|s|*s ) */
    tFilterFf[i] /= volTotal;
/*... tesFilter(s) -> vol*s */
    tFilterS[i] /= volTotal;

    ff[i] = tmp*tFilterS[i];
    mf[i] = 2.0e0*(delta*tFilterFf[i] - ff[i]); 
//  mf[i] = -2.0e0*(ff[i] - tFilterFf[i]);         
  }

/*... Ld:Ms*/
  lDynamic[0] = doubleDotSym2(l,ms);
/*... Ld:Mf*/
  lDynamic[1] = doubleDotSym2(l,mf);
/*... Ms:Ms*/
  lDynamic[2] = doubleDotSym(ms);
/*... Ms:Mf*/
  lDynamic[3] = doubleDotSym2(ms,mf);
/*... Mf:Ms*/
  lDynamic[4] = lDynamic[3];
/*... Mf:Mf*/
  lDynamic[5] = doubleDotSym(mf);
  
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 21/09/2017                                   *
 * Data de modificaco : 11/01/2017                                   *
 *-------------------------------------------------------------------*
 * WALLMODEL : modelo de parade                                      *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * vt        -> velocidade tangencial a parede da celula             *
 * viscosity -> viscosidade dinamica molecular                       *
 * density   -> densidade                                            *
 * dWall     -> distancia da celcula a parede                        *
 * iCod      -> tipo de funcao de parede                             *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * y+ -> retorna o y+                                                *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void wallModel(DOUBLE const vt     , DOUBLE const viscosity
             , DOUBLE const density, DOUBLE const dWall 
             , DOUBLE *yP          , DOUBLE *uP   
             , short const iCod) {
  
  unsigned short i, maxIt=10000;
  DOUBLE stressW,f,df,tol = 1.e-11, nu = viscosity/density;
  DOUBLE yPlus,uPlus,fu,fu0, conv,onePlusB,yPlus3,yPlus4,uPlusL,uPlusT,
         gama,gamai,temp1;
  DOUBLE invKarman = 1.e0/VONKARMAN;

  yPlus = uPlus = 0.e0;

  switch (iCod) {

/*...*/
    case STANDARDWALL:
/* ... wall shear stress (viscosidade)*/
      stressW = viscosity*vt/dWall;
/* ... friction velocity*/
      fu     = sqrt(stressW/density);
      conv   = fu*tol;
      for(i = 0; i < maxIt; i++){
        fu0      = fu;
/*...*/
        yPlus    = fu*dWall/nu;
        uPlus    = vt/fu;
/*...................................................................*/
        
/*...*/
        if( yPlus < 11.81e0){
          f   = yPlus - uPlus;
          df  = (yPlus + uPlus)/nu;
          fu -= f/df;
        }
        else {
          f  = invKarman*log(E_WALLMODEL*yPlus) - uPlus;
          df = (1.e0 + VONKARMAN*uPlus)/(VONKARMAN*fu);
          fu -= f/df;
        }
        if(fabs(fu-fu0) < conv) break; 
      }
/*...................................................................*/ 
      yPlus    = fu*dWall/nu;      
      if ( yPlus < 0.e0 )
        ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__,"WallModel: y+ < 0"
                  ,EXIT_PROG);   
/*...................................................................*/        
      break;
/*...................................................................*/

/*...*/
    case ENHANCEDWALL:
/* ... wall shear stress (viscosidade)*/
      stressW = viscosity*vt/dWall;
/* ... friction velocity*/
      fu     = sqrt(stressW/density);
      conv   = fu*tol;
      for(i = 0; i < maxIt; i++){
        fu0      = fu;
/*...*/
        yPlus    = fu*dWall/nu;
        uPlus    = vt/fu;
/*...................................................................*/
        
/*...*/
        onePlusB = 1.0e0 + BKADER*yPlus; 
        yPlus3   = yPlus*yPlus*yPlus;
        yPlus4   = yPlus3*yPlus;
/*... u+ laminar (linear)*/
        uPlusL   = yPlus;
/* ... u+ trubulento (log)*/
        uPlusT   = invKarman*log(E_WALLMODEL*yPlus);
/* ... parametro T = -(a(y+)**4)/(1+by+)*/
        gama  = - (AKADER*yPlus4) / onePlusB;
        gamai = 1.0e0/gama;
/*... derivada dT/dy+*/
        temp1 = AKADER*BKADER*yPlus4-4.0e0*AKADER*yPlus3*onePlusB;
        temp1 = temp1/(onePlusB*onePlusB);
/* ... f = e(T)*(y+) + e(1/T)*(1/k)ln(Ey+) - u+*/
        f = exp(gama)*uPlusL + exp(gamai)*uPlusT - uPlus;
/* ... derivada de f*/
        df = (temp1*yPlus+1.0e0)*exp(gama)*yPlus   
           + invKarman*exp(gamai)*( 1.0e0 
           - (temp1/(gama*gama))*yPlus*log(E_WALLMODEL*yPlus))
           + uPlus;
/* .....................................................................*/

/*...*/
//      fu -= f/df;
        fu = fu*(1.e0 - f/df);
        if(fabs(fu-fu0) < conv) break; 
      }
/*...................................................................*/ 
      yPlus    = fu*dWall/nu;      
      if ( yPlus < 0.e0 )
        ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__
                  ,"WallModel: y+ < 0",EXIT_PROG);   
/*...................................................................*/        
      break;
/*...................................................................*/


/*...*/
      default: 
        ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__,
                  "Wall Model: Opcao invalida !! ",EXIT_PROG);
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  *yP = yPlus;
  *uP = uPlus;
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 23/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * WALLMODELHEAT : modelo de parade par conducao termica             *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * yPLus     -> distancia normalizada a parede                       *
 * prM       -> numero de Prandlt molecular                          *
 * prT       -> numero de Prandlt turbulento                         * 
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * y+ -> retorna o y+                                                *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE wallModelHeat(DOUBLE const yPlus,DOUBLE const prM
                   , DOUBLE const prT) {
  
  DOUBLE tempPlus,tmp,bt;  
  
  if (yPlus < 11.81e0) 
    tempPlus = prM*yPlus;
  else{
    tmp = 3.85e0*pow(prM,D1DIV3);
    bt = tmp*tmp + 2.12*log (prM);
    tempPlus = (prT/VONKARMAN)*log(yPlus) + bt;
  }
  
  return tempPlus;
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/11/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * wallParameters : calcula a distancia a parede e a distancia       *
 * normalizada                                                       *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  *
 * lViz      -> viznhos da celula central                            *
 * vel       -> velocidades fica sem variacao temporal               *
 * normal    -> vetores normais as faces das celulas                 *
 * dcca      -> menor distancia do centroide central a faces desta   *
 *              celula                                               *
 * faceVelR  -> restricoes por elemento de velocidades               *
 * faceVelL  -> carga por elemento de velocidades                    *
 * dViscosity-> viscosidade dinamica                                 *
 * density   -> massa especifica                                     *
 * wallPar   -> parametros de parede ( yPlus, uPlus, uFri,sTressW)   *            * 
 * dWall     -> nao definido                                         *
 * nFace     -> numero de faces da celula central                    * 
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 * dWall - distancia a parede                                        *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
bool wallParameters(Loads *lVel               
             , INT *RESTRICT lViz       , DOUBLE *RESTRICT v
             , DOUBLE *RESTRICT normal  , DOUBLE *RESTRICT dcca
             , short *RESTRICT lFaceVelR, short *RESTRICT lFaceVelL
             , DOUBLE const viscosity   , DOUBLE const density
             , DOUBLE *RESTRICT wallPar , DOUBLE *dWall     
             , short const wallType     , short const nFace )
{ 
  bool fWall = false;
  short nf, nCarg;
  INT vizNel;
  DOUBLE wt, yPlusMax, uPlusMax,uFricMax,sWallMax,par[MAXLOADPARAMETER];
  DOUBLE dMin, vParallel[3], lNormal[3],yPlus,uPlus,uFric,sW,velB[3],xx[4];

  sWallMax = uFricMax = uPlusMax = yPlusMax = 0.e0;
  dMin     = 1.e+16;
  for (nf = 0; nf<nFace; nf++) {
    vizNel         = lViz[nf]; 
    yPlus = 0.e0;
    uPlus = 1.e0; 
/*... contorno*/
    if (vizNel  == -2) {
      if (lFaceVelR[nf] > 0) {
        nCarg = lFaceVelL[nf] - 1;
        if( lVel[nCarg].type == MOVEWALL){
          fWall = true;
/*... velocidade da parede*/
          nCarg     = lFaceVelL[nf] - 1;
          getLoads(par,&lVel[nCarg],xx);
          velB[0]   = par[0];
          velB[1]   = par[1];
          velB[2]   = par[2];
/*...*/
          lNormal[0] = MAT2D(nf, 0, normal, 3);
          lNormal[1] = MAT2D(nf, 1, normal, 3);   
          lNormal[2] = MAT2D(nf, 2, normal, 3);   
/*... calculo da velociade paralela a face*/
          wt = v[0] * lNormal[0] 
             + v[1] * lNormal[1] 
             + v[2] * lNormal[2];
          vParallel[0] = v[0] - wt * lNormal[0] - velB[0];
          vParallel[1] = v[1] - wt * lNormal[1] - velB[1];
          vParallel[2] = v[2] - wt * lNormal[2] - velB[2];
          wt = vParallel[0]*vParallel[0] 
             + vParallel[1]*vParallel[1]
             + vParallel[1]*vParallel[1];
          wt = sqrt(wt);
/*...*/
          if (wt > 0.e0)
            wallModel(wt      , viscosity
                     , density , dcca[nf]
                     , &yPlus  , &uPlus
                     , wallType); 
/*...................................................................*/
          dMin     = min(dcca[nf],dMin);
          yPlusMax = max(yPlus   ,yPlusMax);
          uPlusMax = max(uPlus   ,uPlusMax);
          uFric    = wt/uPlus;
          sW       = density*uFric*uFric;
          uFricMax = max(uFric,uFricMax);
          sWallMax = max(sW   ,sWallMax);
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...*/
      else if (lFaceVelR[nf] == STATICWALL){
        fWall = true;
        lNormal[0] = MAT2D(nf, 0, normal, 3);
        lNormal[1] = MAT2D(nf, 1, normal, 3);   
        lNormal[2] = MAT2D(nf, 2, normal, 3);   
/*... calculo da velociade paralel a a face*/
        wt = v[0] * lNormal[0] 
           + v[1] * lNormal[1]
           + v[2] * lNormal[2];
        vParallel[0] = v[0] - wt * lNormal[0];
        vParallel[1] = v[1] - wt * lNormal[1];
        vParallel[2] = v[2] - wt * lNormal[2];
        wt = vParallel[0]*vParallel[0] 
           + vParallel[1]*vParallel[1]
           + vParallel[2]*vParallel[2];
        wt = sqrt(wt);
/*...*/
        if (wt > 0.e0)
          wallModel(wt     , viscosity
                 , density , dcca[nf]
                 , &yPlus  , &uPlus
                 , wallType); 
/*...................................................................*/
        dMin     = min(dcca[nf],dMin);
        yPlusMax = max(yPlus   ,yPlusMax);
        uPlusMax = max(uPlus   ,uPlusMax);
        uFric    = wt/uPlus;
        sW       = density*uFric*uFric;
        uFricMax = max(uFric,uFricMax);
        sWallMax = max(sW   ,sWallMax);
      }
/*...................................................................*/
    }
/*...................................................................*/
  }

  wallPar[0] = yPlusMax;
  wallPar[1] = uPlusMax;
  wallPar[2] = uFricMax;
  wallPar[3] = sWallMax;
  *dWall     = dMin;

  return  fWall; 
/*...................................................................*/

}

/*********************************************************************
 * Data de criacao    : 27/11/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * sigmaModel: Metodo dinamico para o calculo da constente do LES    *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * gradVel   -> gradiente rescontruido das velocidades               *
 * nFace     -> numero de faces da celula central                    * 
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * retorna sg3*(sg1 -sg2)*(sg2 -sg3)/sg1^2                           *                                     *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * delta a largura do filtro que é igual a raiz cubica do volume da  *
 * celula central                                                    *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE sigmaModel(DOUBLE *RESTRICT gradVel 
                , short const nFace       , short const ndm) {

  INT ip[3];
  DOUBLE g[9],x[9],sg[3],op;
 
/*... g11*/
  g[0] = MAT2D(0,0,gradVel,3)*MAT2D(0,0,gradVel,3)
       + MAT2D(0,1,gradVel,3)*MAT2D(0,1,gradVel,3)
       + MAT2D(0,2,gradVel,3)*MAT2D(0,2,gradVel,3);
/*... g12*/
  g[1] = MAT2D(0,0,gradVel,3)*MAT2D(1,0,gradVel,3)
       + MAT2D(0,1,gradVel,3)*MAT2D(1,1,gradVel,3)
       + MAT2D(0,2,gradVel,3)*MAT2D(1,2,gradVel,3);
/*... g13*/
  g[2] = MAT2D(0,0,gradVel,3)*MAT2D(2,0,gradVel,3)
       + MAT2D(0,1,gradVel,3)*MAT2D(2,1,gradVel,3)
       + MAT2D(0,2,gradVel,3)*MAT2D(2,2,gradVel,3);
/*... g21*/
  g[3] = g[1];
/*... g22*/
  g[4] = MAT2D(1,0,gradVel,3)*MAT2D(1,0,gradVel,3)
       + MAT2D(1,1,gradVel,3)*MAT2D(1,1,gradVel,3)
       + MAT2D(1,2,gradVel,3)*MAT2D(1,2,gradVel,3);
/*... g23*/
  g[5] = MAT2D(1,0,gradVel,3)*MAT2D(2,0,gradVel,3)
       + MAT2D(1,1,gradVel,3)*MAT2D(2,1,gradVel,3)
       + MAT2D(1,2,gradVel,3)*MAT2D(2,2,gradVel,3);
/*... g31*/
  g[6] = g[2];   
/*... g32*/              
  g[7] = g[5];
/*... g33*/
  g[8] = MAT2D(2,0,gradVel,3)*MAT2D(2,0,gradVel,3)
       + MAT2D(2,1,gradVel,3)*MAT2D(2,1,gradVel,3)
       + MAT2D(2,2,gradVel,3)*MAT2D(2,2,gradVel,3);
/*...................................................................*/

/*... autovalores*/
  cyclic_jacobi(g    ,x
               ,sg   ,ip
               ,3    ,1.e-11
               ,10000,false);  
/*...................................................................*/

/*... valores singulares de g*/
  sg[0] = sqrt(fabs(sg[0]));
  sg[1] = sqrt(fabs(sg[1]));
  sg[2] = sqrt(fabs(sg[2]));
/*...................................................................*/

/*...*/  
   op = sg[2]*(sg[0]-sg[1])*(sg[1]-sg[2])/(sg[0]*sg[0]);  
// printf("%e %e %e %e\n",s[0],s[1],s[2],op);
/*...................................................................*/
   return op;
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 04/12/2017                                   *
 * Data de modificaco : 13/12/2017                                   * 
 *-------------------------------------------------------------------* 
 * structuralStress : calculo do tensor residual por metodos         *
 * estrutarais                                                       *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * tModel    -> modelo de turbulencia                                *
 * lViz      -> vizinhos dos elementos                               * 
 * stressR   -> nao definido                                         *
 * nVel      -> velocidade nodal                                     *
 * eVel      -> velocidade por elemento                              *
 * nDensity  -> massa especifica com variacao temporal (nodal)       *
 * eDensity  -> massa especifica com variacao temporal (celula)      * 
 * gradVel   -> gradiente da solucao conhecido                       * 
 * vol       -> volume das celulas                                   *
 * cs        -> constantde dinamica     
 * ndm       -> numero de dimensoes                                  * 
 * nFace     -> numero de faces por celulas                          * 
 * nel       -> numero do elemento                                   *
 * ----------------------------------------------------------------- * 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * s        -> tensor residual                                       *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 * s(s11,s22,s33,s12,s23,s13)                                        *
 * s(sxx,syy,szz,sxy,syz,sxz)                                        *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void structuralStress(Turbulence tModel   , DOUBLE *RESTRICT xl
                , INT *RESTRICT lViz      , DOUBLE *RESTRICT stressR  
                , DOUBLE *RESTRICT nVel   , DOUBLE *RESTRICT eVel
                , DOUBLE *RESTRICT nDen   , DOUBLE *RESTRICT eDen
                , DOUBLE *RESTRICT gradVel, DOUBLE *RESTRICT vol   
                , DOUBLE const cd  
                , short const ndm         , short const nFace
                , INT const nEl ) {

  bool dynamic = tModel.fDynamic;
  short i, typeLes,type, idCell = nFace;
  DOUBLE c,cs,g[3][3],tmp,delta;

/*...*/
  type    = tModel.typeMixed[ESTMODEL];
  typeLes = tModel.typeLes;
  c       = 1.e0-tModel.type;
/*...................................................................*/

/*...*/
  if (dynamic)
    cs = cd;
  else
    cs = tModel.cs;
/*...................................................................*/

/*...*/
  switch (type) {

/*...*/
    case CLARK:

/*... | du1/dx1 du1/dx2 du1/dx3*/
      g[0][0] = MAT3D(idCell,0,0,gradVel,3,3);
      g[0][1] = MAT3D(idCell,0,1,gradVel,3,3);
      g[0][2] = MAT3D(idCell,0,2,gradVel,3,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
      g[1][0] = MAT3D(idCell,1,0,gradVel,3,3);
      g[1][1] = MAT3D(idCell,1,1,gradVel,3,3);
      g[1][2] = MAT3D(idCell,1,2,gradVel,3,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
      g[2][0] = MAT3D(idCell,2,0,gradVel,3,3);
      g[2][1] = MAT3D(idCell,2,1,gradVel,3,3);
      g[2][2] = MAT3D(idCell,2,2,gradVel,3,3);
/*...................................................................*/

/*... delta*delta*/
      delta    = pow(vol[idCell],D2DIV3);
/*...................................................................*/

/*...*/
      tmp = cs*eDen[idCell]*delta/12.e0;
/*...s00 = sxx*/
      stressR[0] = tmp*( g[0][0]*g[0][0] 
                       + g[0][1]*g[0][1] 
                       + g[0][2]*g[0][2]  );
/*...s11 = syy*/
      stressR[1] = tmp*( g[1][0]*g[1][0] 
                       + g[1][1]*g[1][1] 
                       + g[1][2]*g[1][2]  );
/*...s22 = szz*/
      stressR[2] = tmp*( g[2][0]*g[2][0] 
                       + g[2][1]*g[2][1] 
                       + g[2][2]*g[2][2]  );
/*...s01 = sxy*/     
      stressR[3] = tmp*( g[0][0]*g[1][0] 
                       + g[0][1]*g[1][1] 
                       + g[0][2]*g[1][2]  );
/*...s12 = syz*/     
      stressR[4] = tmp*( g[1][0]*g[2][0] 
                       + g[1][1]*g[2][1] 
                       + g[1][2]*g[2][2]  );

/*...s02 = sxz*/     
      stressR[5] = tmp*( g[0][0]*g[2][0] 
                       + g[0][1]*g[2][1] 
                       + g[0][2]*g[2][2]  );  

/*...................................................................*/
    break;
/*...................................................................*/

/*...*/
    case BARDINA:
      bardinaModel(xl       , stressR  
                 , nVel     , nDen
                 , vol[idCell]
                 , cs       , lViz[idCell]
                 , ndm      , nFace);
    break;
/*...................................................................*/

/*...*/
    case BARDINAMOD:
      bardinaModelMod(lViz     , stressR
                     ,eVel     , gradVel
                     ,vol      , eDen
                     ,cs           
                     ,ndm      , nFace);
    break;
/*...................................................................*/

/*...*/
    default: 
      ERRO_OP(__FILE__,__func__,type);
/*...................................................................*/ 
  }
/*...................................................................*/

/*... metodo misto com coeficentes constantes*/
  if (typeLes == LESMIXEDMODEL) {
    for(i=0;i<6;i++)
      stressR[i] *= c; 
  }
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 05/12/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * bardinaModelMod : tensor residual baseado no tensor Bardina       *
 * (Liu-Meneveau-Katz)                                               *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lViz      -> vizinhos dos elementos                               * 
 * stressR   -> nao definido                                         *
 * vel       -> campo de velocidades                                 *
 * gradVel   -> gradiente de velocidades                             * 
 * vol       -> volume da  celula central                            *
 * density   -> massa especifica com variacao temporal               *   
 * ndm       -> numero de dimensoes                                  * 
 * nFace     -> numero de faces por celulas                          * 
 * ----------------------------------------------------------------- * 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * stressR        -> tensor residual                                 *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 * stressR(s11,s22,s33,s12,s23,s13)                                  *
 * stressR(sxx,syy,szz,sxy,syz,sxz)                                  *
 *       | v1 v2 v3 | no1                                            *
 * vel = | v1 v2 v3 | no2                                            *
         | v1 v2 v3 | no3                                            *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void bardinaModel(DOUBLE *RESTRICT xl   , DOUBLE *RESTRICT stressR  
                , DOUBLE *RESTRICT nVel , DOUBLE *RESTRICT nDensity
                , DOUBLE const vol      , DOUBLE const cs       
                , INT const nEl
                , short const ndm       , short const nFace) {

  short i,j,l,k,nint;
  bool af,dev,ninter;
  DOUBLE hx[8],hy[8],hz[8],N[8],de[8],dn[8],dz[8];
  DOUBLE alf,bet,teta,eps,nn,ze,det,wt,denI,velI[3];
  DOUBLE denVel[3],denVelVel[6],den;

/*...*/
  den          = 0.e0;
  denVelVel[0] = 0.e0;
  denVelVel[1] = 0.e0;
  denVelVel[2] = 0.e0;
  denVelVel[3] = 0.e0;
  denVelVel[4] = 0.e0;
  denVelVel[5] = 0.e0;
  denVel[0]    = 0.e0;
  denVel[1]    = 0.e0;
  denVel[2]    = 0.e0;  
/*...................................................................*/

/*... Matriz de rigidez - 2 pontos de integração*/
  nint  =     3;
  af    = true;/* calculo do jacobiano*/
  ninter= true;/* funcoes de interpolacao*/
  dev   = true;/* derivadas das funcoes de interpolacao*/
  for(i=0;i<nint;i++){      
    eps = intponto(nint,i);
    alf = peso(nint,i);
    for(j=0;j<nint;j++){
      nn = intponto(nint,j);
      bet= peso(nint,j);
      for(l=0;l<nint;l++){
        ze   = intponto(nint,l);
	      teta = peso(nint,l);

/*...*/
        sfHexa8(eps,nn,ze,N,de,dn,dz,ninter,dev);
        det  = jacob3d(xl,de,dn,dz,hx,hy,hz,8,af,dev,nEl);
/*...................................................................*/

/*...*/
        denI    = 0.e0;      
        velI[0] = 0.e0;
        velI[1] = 0.e0;
        velI[2] = 0.e0;
/*...................................................................*/

/*...*/
        for (k = 0; k < 8; k++) {
          denI    += nDensity[k]*N[k];
          velI[0] += MAT2D(k,0,nVel,3)*N[k];
          velI[1] += MAT2D(k,1,nVel,3)*N[k];
          velI[2] += MAT2D(k,2,nVel,3)*N[k]; 
        }
/*...................................................................*/

/*...*/         
        wt            = det*alf*bet*teta;
        den          += denI*wt;

        denVelVel[0] += (denI*velI[0]*velI[0])*wt;
        denVelVel[1] += (denI*velI[1]*velI[1])*wt;
        denVelVel[2] += (denI*velI[2]*velI[2])*wt;
        denVelVel[3] += (denI*velI[0]*velI[1])*wt;
        denVelVel[4] += (denI*velI[1]*velI[2])*wt;
        denVelVel[5] += (denI*velI[0]*velI[2])*wt;

        denVel[0]    += denI*velI[0]*wt;
        denVel[1]    += denI*velI[1]*wt;
        denVel[2]    += denI*velI[2]*wt;
/*...................................................................*/
       }
    }	
  }       

/*... filtro(den)*/
  den       /= vol;
/*... filtro(den*v)*/
  denVel[0] /= vol;
  denVel[1] /= vol;
  denVel[2] /= vol;
/*... filtro(den*vv)*/
  denVelVel[0] /= vol;
  denVelVel[1] /= vol;
  denVelVel[2] /= vol;
  denVelVel[3] /= vol;
  denVelVel[4] /= vol;
  denVelVel[5] /= vol;

/*...*/
  stressR[0] = cs*(denVelVel[0] - denVel[0]*denVel[0]/den);
  stressR[1] = cs*(denVelVel[1] - denVel[1]*denVel[1]/den);
  stressR[2] = cs*(denVelVel[2] - denVel[2]*denVel[2]/den);
  stressR[3] = cs*(denVelVel[3] - denVel[0]*denVel[1]/den);
  stressR[4] = cs*(denVelVel[4] - denVel[1]*denVel[2]/den);
  stressR[5] = cs*(denVelVel[5] - denVel[0]*denVel[2]/den);
//  printf("%e %e %e %e %e %e %e %e\n",stressR[0],stressR[1],stressR[2]
//                             ,stressR[3],stressR[4],stressR[5],vol,den);
/*...................................................................*/

}
/********************************************************************/

/********************************************************************* 
 * Data de criacao    : 05/12/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * bardinaModelMod : tensor residual baseado no tensor Bardina       *
 * (Liu-Meneveau-Katz)                                               *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lViz      -> vizinhos dos elementos                               * 
 * stressR   -> nao definido                                         *
 * vel       -> campo de velocidades                                 *
 * gradVel   -> gradiente de velocidades                             * 
 * vol       -> volume das celulas                                   *
 * density   -> massa especifica com variacao temporal               *   
 * ndm       -> numero de dimensoes                                  * 
 * nFace     -> numero de faces por celulas                          * 
 * ----------------------------------------------------------------- * 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * stressR        -> tensor residual                                 *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 * stressR(s11,s22,s33,s12,s23,s13)                                  *
 * stressR(sxx,syy,szz,sxy,syz,sxz)                                  *
 * fonte:                                                            *
 * Large Eddy Simulation for incompressible Flows - Pierre Sagaut    *
 * 3° Ed                                                             * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void bardinaModelMod(INT *RESTRICT lViz  , DOUBLE *RESTRICT stressR  
                   , DOUBLE *RESTRICT vel, DOUBLE *RESTRICT gradVel
                   , DOUBLE *RESTRICT vol, DOUBLE *RESTRICT lDensity
                   , DOUBLE const cs 
                   , short const ndm     , short const nFace) {

  short i, idCell = nFace;
  INT vizNel;
  DOUBLE volW,volTotal,tmp,g[3][3]; 
  DOUBLE v[3],filterVv[6],filterV[3],l[6],s[6],ll,lm,ss,iLs,fs;

/*... | du1/dx1 du1/dx2 du1/dx3*/
  g[0][0] = MAT3D(idCell,0,0,gradVel,3,3);
  g[0][1] = MAT3D(idCell,0,1,gradVel,3,3);
  g[0][2] = MAT3D(idCell,0,2,gradVel,3,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
  g[1][0] = MAT3D(idCell,1,0,gradVel,3,3);
  g[1][1] = MAT3D(idCell,1,1,gradVel,3,3);
  g[1][2] = MAT3D(idCell,1,2,gradVel,3,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
  g[2][0] = MAT3D(idCell,2,0,gradVel,3,3);
  g[2][1] = MAT3D(idCell,2,1,gradVel,3,3);
  g[2][2] = MAT3D(idCell,2,2,gradVel,3,3);
/*...................................................................*/

  volTotal = 0.e0;
  filterV[0] = filterV[1] = filterV[2] = 0.e0;
  for(i = 0; i<6 ; i++ )
    filterVv[i]    = 0.e0;

/*...*/
  for(i = 0; i<nFace+1;i++ ){
    vizNel = 0;
    if( i < nFace )
      vizNel = lViz[i];
/*... celula central e elementos vizinhos*/        
    if( vizNel  > -1 ){
      volW      = vol[i];
      volTotal += volW;

/*...*/
      v[0] = MAT2D(i, 0, vel, 3);
      v[1] = MAT2D(i, 1, vel, 3);
      v[2] = MAT2D(i, 2, vel, 3);
/*...................................................................*/

/*... tesFilter(vv) -> vol**vv */
      tmp          = volW;
      filterVv[0] += tmp*v[0]*v[0]; /*v1v1*/
      filterVv[1] += tmp*v[1]*v[1]; /*v2v2*/
      filterVv[2] += tmp*v[2]*v[2]; /*v3v3*/
      filterVv[3] += tmp*v[0]*v[1]; /*v1v2*/
      filterVv[4] += tmp*v[0]*v[2]; /*v1v3*/
      filterVv[5] += tmp*v[1]*v[2]; /*v2v3*/
/*... tesFilter(desinty*v) -> vol*desinty*v */
      filterV[0] += tmp*v[0]; 
      filterV[1] += tmp*v[1]; 
      filterV[2] += tmp*v[2];
    }
/*...................................................................*/
  } 
/*...................................................................*/

/*... tesFilter(vv)*/
  filterVv[0] /= volTotal; /*v1v1*/
  filterVv[1] /= volTotal; /*v2v2*/
  filterVv[2] /= volTotal; /*v3v3*/
  filterVv[3] /= volTotal; /*v1v2*/
  filterVv[4] /= volTotal; /*v1v3*/
  filterVv[5] /= volTotal; /*v2v3*/
/*... tesFilter(v)*/
  filterV[0] /= volTotal; 
  filterV[1] /= volTotal; 
  filterV[2] /= volTotal;

/*... bardina*/
  l[0] = filterVv[0] - filterV[0]*filterV[0];
  l[1] = filterVv[1] - filterV[1]*filterV[1];
  l[2] = filterVv[2] - filterV[2]*filterV[2];
  l[3] = filterVv[3] - filterV[0]*filterV[1];
  l[4] = filterVv[4] - filterV[1]*filterV[2];
  l[5] = filterVv[5] - filterV[0]*filterV[2]; 
/*...................................................................*/

/*.. calculo Sij*/
  s[0] = g[0][0];
  s[1] = g[1][1];
  s[2] = g[2][2];
  s[3] = 0.5e0*(g[0][1] + g[1][0]); /*s12*/
  s[4] = 0.5e0*(g[0][2] + g[2][0]); /*s13*/
  s[5] = 0.5e0*(g[1][2] + g[2][1]); /*s23*/
/*...................................................................*/

/*... LijSij*/
  lm = l[0]*s[0] + l[1]*s[1] + l[2]*s[2]
     + 2.e0*( l[3]*s[3] + l[4]*s[4] + l[5]*s[5]);
/*... LijLij*/
  ll = doubleDotSym(l);
/*... LijLij*/
  ss = doubleDotSym(s);  

/*...*/
  iLs = lm/(sqrt(ll)*sqrt(ss));
  if(iLs >= 0.0)
    fs = 1.e0 - exp(-10.e0*iLs*iLs);
  else
    fs = 0.e0;
/*...................................................................*/

//printf("%e\n",fs);
/*...*/    
  tmp = cs*fs*lDensity[idCell];
  stressR[0] = tmp*fs*l[0];
  stressR[1] = tmp*fs*l[1];
  stressR[2] = tmp*fs*l[2];
  stressR[3] = tmp*fs*l[3];
  stressR[4] = tmp*fs*l[4];
  stressR[5] = tmp*fs*l[5];
/*...................................................................*/

}
/********************************************************************/

/********************************************************************* 
 * Data de criacao    : 13/01/2018                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * vremanModel : calculo do operador do modelo VREMAN                * 
 *-------------------------------------------------------------------* 
 * gradVel  -> gradiente de velociade                                *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*                                                                      *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE vremanModel(DOUBLE *RESTRICT gradVel) {

  short i,j;
  DOUBLE bb,aa,g[3][3],b[3][3];

/*... | du1/dx1 du1/dx2 du1/dx3*/
  g[0][0] = MAT2D(0,0,gradVel,3);
  g[0][1] = MAT2D(0,1,gradVel,3);
  g[0][2] = MAT2D(0,2,gradVel,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
  g[1][0] = MAT2D(1,0,gradVel,3);
  g[1][1] = MAT2D(1,1,gradVel,3);
  g[1][2] = MAT2D(1,2,gradVel,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
  g[2][0] = MAT2D(2,0,gradVel,3);
  g[2][1] = MAT2D(2,1,gradVel,3);
  g[2][2] = MAT2D(2,2,gradVel,3);
/*...................................................................*/

/*... beta*/
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      b[i][j] = ( g[0][i]*g[0][j] 
                + g[1][i]*g[1][j]
                + g[2][i]*g[2][j]);
/*...*/
  bb = b[0][0]*b[1][1] - b[0][1]*b[0][1] + b[0][0]*b[2][2]
     - b[0][2]*b[0][2] + b[1][1]*b[2][2] - b[1][2]*b[1][2];
  if(bb < 0.e0) bb = 0.0e0;
/*... |a:a|*/
  aa = doubleDot(g[0]);
/*...................................................................*/

/*...*/
  return sqrt(bb/aa);
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 17/12/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * waleModel : calculo do operador do modelo WALE                    * 
 *-------------------------------------------------------------------* 
 * gradVel  -> gradiente de velociade                                *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*                                                                      *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE waleModel(DOUBLE *RESTRICT gradVel) {

  DOUBLE lm,mm,modS,modSd,tmp,g[3][3],b[3][3],m[6],s[6],ce;

/*... | du1/dx1 du1/dx2 du1/dx3*/
  g[0][0] = MAT2D(0,0,gradVel,3);
  g[0][1] = MAT2D(0,1,gradVel,3);
  g[0][2] = MAT2D(0,2,gradVel,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
  g[1][0] = MAT2D(1,0,gradVel,3);
  g[1][1] = MAT2D(1,1,gradVel,3);
  g[1][2] = MAT2D(1,2,gradVel,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
  g[2][0] = MAT2D(2,0,gradVel,3);
  g[2][1] = MAT2D(2,1,gradVel,3);
  g[2][2] = MAT2D(2,2,gradVel,3);
/*...................................................................*/

/*... Sd = Sd:Sd*/
/*... g00 = g00*g00 + g01*g10 + g02*g20
      g01 = g00*g01 + g01*g11 + g02*g21
      g02 = g00*g02 + g01*g12 + g02*g22*/ 
  b[0][0] = g[0][0]*g[0][0] + g[0][1]*g[1][0] + g[0][2]*g[2][0];
  b[0][1] = g[0][0]*g[0][1] + g[0][1]*g[1][1] + g[0][2]*g[2][1]; 
  b[0][2] = g[0][0]*g[0][2] + g[0][1]*g[1][2] + g[0][2]*g[2][2]; 
/*... g10 = g10*g00 + g11*g10 + g12*g20
      g11 = g10*g01 + g11*g11 + g12*g21
      g12 = g10*g02 + g11*g12 + g12*g22*/  
  b[1][0] = g[1][0]*g[0][0] + g[1][1]*g[1][0] + g[1][2]*g[2][0]; 
  b[1][1] = g[1][0]*g[0][1] + g[1][1]*g[1][1] + g[1][2]*g[2][1]; 
  b[1][2] = g[1][0]*g[0][2] + g[1][1]*g[1][2] + g[1][2]*g[2][2]; 
/*... g20 = g20*g00 + g21*g10 + g22*g20
      g21 = g20*g01 + g21*g11 + g22*g21
      g22 = g20*g02 + g21*g12 + g22*g22*/ 
  b[2][0] = g[2][0]*g[0][0] + g[2][1]*g[1][0] + g[2][2]*g[2][0]; 
  b[2][1] = g[2][0]*g[0][1] + g[2][1]*g[1][1] + g[2][2]*g[2][1]; 
  b[2][2] = g[2][0]*g[0][2] + g[2][1]*g[1][2] + g[2][2]*g[2][2]; 
/*...*/
  tmp  = D1DIV3*(b[0][0] + b[1][1] + b[2][2]);
  m[0] = b[0][0] - tmp;
  m[1] = b[1][1] - tmp;
  m[2] = b[2][2] - tmp;     
  m[3] = 0.5e0*(b[0][1] + b[1][0]); /*g12*/
  m[4] = 0.5e0*(b[0][2] + b[2][0]); /*g13*/
  m[5] = 0.5e0*(b[1][2] + b[2][1]); /*g23*/
   
  modSd = doubleDotSym(m);

/*... S = S:S*/
  s[0] = g[0][0];
  s[1] = g[1][1];
  s[2] = g[2][2];
  s[3] = 0.5e0*(g[0][1] + g[1][0]); /*s12*/
  s[4] = 0.5e0*(g[0][2] + g[2][0]); /*s13*/
  s[5] = 0.5e0*(g[1][2] + g[2][1]); /*s23*/

  modS = doubleDotSym(s);
/*...................................................................*/
      
/*...*/
  mm = pow(modS,2.5e0) + pow(modSd,1.25e0);
  lm = pow(modSd,1.5e0);
  ce = lm/mm;
  if (mm == 0.0e0) ce = 0.e0;
/*...................................................................*/

/*...*/
  return ce;
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 17/12/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * waleModel : calculo do operador do modelo WALE                    * 
 *-------------------------------------------------------------------* 
 * gradVel  -> gradiente de velociade                                *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*                                                                      *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE svs(DOUBLE *RESTRICT gradVel) {

  DOUBLE lm,mm,modS,modSd,tmp,g[3][3],b[3][3],m[6],s[6];

/*... | du1/dx1 du1/dx2 du1/dx3*/
  g[0][0] = MAT2D(0,0,gradVel,3);
  g[0][1] = MAT2D(0,1,gradVel,3);
  g[0][2] = MAT2D(0,2,gradVel,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
  g[1][0] = MAT2D(1,0,gradVel,3);
  g[1][1] = MAT2D(1,1,gradVel,3);
  g[1][2] = MAT2D(1,2,gradVel,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
  g[2][0] = MAT2D(2,0,gradVel,3);
  g[2][1] = MAT2D(2,1,gradVel,3);
  g[2][2] = MAT2D(2,2,gradVel,3);
/*...................................................................*/

/*... Sd = Sd:Sd*/
/*... g00 = g00*g00 + g01*g10 + g02*g20
      g01 = g00*g01 + g01*g11 + g02*g21
      g02 = g00*g02 + g01*g12 + g02*g22*/ 
  b[0][0] = g[0][0]*g[0][0] + g[0][1]*g[1][0] + g[0][2]*g[2][0];
  b[0][1] = g[0][0]*g[0][1] + g[0][1]*g[1][1] + g[0][2]*g[2][1]; 
  b[0][2] = g[0][0]*g[0][2] + g[0][1]*g[1][2] + g[0][2]*g[2][2]; 
/*... g10 = g10*g00 + g11*g10 + g12*g20
      g11 = g10*g01 + g11*g11 + g12*g21
      g12 = g10*g02 + g11*g12 + g12*g22*/  
  b[1][0] = g[1][0]*g[0][0] + g[1][1]*g[1][0] + g[1][2]*g[2][0]; 
  b[1][1] = g[1][0]*g[0][1] + g[1][1]*g[1][1] + g[1][2]*g[2][1]; 
  b[1][2] = g[1][0]*g[0][2] + g[1][1]*g[1][2] + g[1][2]*g[2][2]; 
/*... g20 = g20*g00 + g21*g10 + g22*g20
      g21 = g20*g01 + g21*g11 + g22*g21
      g22 = g20*g02 + g21*g12 + g22*g22*/ 
  b[2][0] = g[2][0]*g[0][0] + g[2][1]*g[1][0] + g[2][2]*g[2][0]; 
  b[2][1] = g[2][0]*g[0][1] + g[2][1]*g[1][1] + g[2][2]*g[2][1]; 
  b[2][2] = g[2][0]*g[0][2] + g[2][1]*g[1][2] + g[2][2]*g[2][2]; 
/*...*/
  tmp  = D1DIV3*(b[0][0] + b[1][1] + b[2][2]);
  m[0] = b[0][0] - tmp;
  m[1] = b[1][1] - tmp;
  m[2] = b[2][2] - tmp;     
  m[3] = 0.5e0*(b[0][1] + b[1][0]); /*g12*/
  m[4] = 0.5e0*(b[0][2] + b[2][0]); /*g13*/
  m[5] = 0.5e0*(b[1][2] + b[2][1]); /*g23*/
   
  modSd = doubleDotSym(m);

/*... S = S:S*/
  s[0] = g[0][0];
  s[1] = g[1][1];
  s[2] = g[2][2];
  s[3] = 0.5e0*(g[0][1] + g[1][0]); /*s12*/
  s[4] = 0.5e0*(g[0][2] + g[2][0]); /*s13*/
  s[5] = 0.5e0*(g[1][2] + g[2][1]); /*s23*/

  modS = doubleDotSym(s);
/*...................................................................*/
      
/*...*/
  mm = lm = 1.0e-64;
  mm = pow(modSd,1.5e0) + pow(modS,3.e0);
  lm = pow(modSd,1.5e0);
/*...................................................................*/

/*...*/
  return lm/mm;
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 16/12/2017                                   *
 * Data de modificaco : 21/01/2018                                   * 
 *-------------------------------------------------------------------* 
 * tensorS : calculo de Sij                                          * 
 *-------------------------------------------------------------------* 
 * s        -> nao definido                                          * 
 * gradVel  -> gradiente de velociade                                *
 * flag     -> true|false                                            *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * flag     -> true -  Sij - (1/3)*Skk*deltaij                       * 
 *             false-  Sij                                           *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 * s =(s11,s22,s33,s12,s23,s13)                                      *
 *********************************************************************/
void tensorS( DOUBLE *RESTRICT s,  DOUBLE *RESTRICT gradVel
            , bool const flag ) {

  DOUBLE tmp,g[3][3];

/*... | du1/dx1 du1/dx2 du1/dx3*/
  g[0][0] = MAT2D(0,0,gradVel,3);
  g[0][1] = MAT2D(0,1,gradVel,3);
  g[0][2] = MAT2D(0,2,gradVel,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
  g[1][0] = MAT2D(1,0,gradVel,3);
  g[1][1] = MAT2D(1,1,gradVel,3);
  g[1][2] = MAT2D(1,2,gradVel,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
  g[2][0] = MAT2D(2,0,gradVel,3);
  g[2][1] = MAT2D(2,1,gradVel,3);
  g[2][2] = MAT2D(2,2,gradVel,3);
/*...................................................................*/

  s[0] = g[0][0];
  s[1] = g[1][1];
  s[2] = g[2][2];
  s[3] = 0.5e0*(g[0][1] + g[1][0]); /*s12*/
  s[4] = 0.5e0*(g[0][2] + g[2][0]); /*s13*/
  s[5] = 0.5e0*(g[1][2] + g[2][1]); /*s23*/ 

  if (flag) {
    tmp  = D1DIV3*(s[0] + s[1] + s[2]);
    s[0] -= tmp;
    s[1] -= tmp;
    s[2] -= tmp;
  }

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 21/01/2018                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * lesOneKeqPar : calculo do coefciente do les de um paramento       * 
 *-------------------------------------------------------------------*
 * tModel   -> modelo de turbulencia                                 *
 * lViz     -> vizinhos                                              * 
 * dynamic  -> L:L, M:M e Ce  nas celulas central e vizinhas         *
 * cDyn     -> nao definido                                          *
 * nFace    -> numero de faces da celula                             *
 * cap      -> true limita o range de variacao                       *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * cDyn     ->(Ck,Ce)                                                *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 * dynamic(*,3) = ( L:L, M:M, Ce )                                   *
 *********************************************************************/
void lesOneKeqPar(Turbulence *tModel, INT *RESTRICT lViz
         , DOUBLE *RESTRICT dynamic , DOUBLE *RESTRICT cDyn
         , DOUBLE *RESTRICT volume  , short const nFace
         , INT const nel            , short const nPar 
         , bool const fCap) {

  short i;
  INT vizNel;
  DOUBLE volW,volTotal;
  DOUBLE capCkFmin,capCkFmax,capCeFmin,capCeFmax;
  DOUBLE fC[4];

/*...*/  
  capCkFmin = 0.0e0;
  capCeFmin = 0.e0;
  capCkFmax = tModel->eK.ck;
  capCeFmax = tModel->eK.ce;   
/*...................................................................*/

/*...*/
  fC[0] = fC[1] = fC[2] = fC[3] = volTotal = 0.e0;
  for(i = 0; i<nFace+1;i++ ){
    vizNel = 0;
    if( i < nFace )
      vizNel = lViz[i];
/*... celula central e elementos vizinhos*/        
    if( vizNel  > -1 ){
      volW      = volume[i];
      volTotal += volW;
/*...*/      
      fC[0] += volW*MAT2D(i,0,dynamic,nPar);
      fC[1] += volW*MAT2D(i,1,dynamic,nPar);
      fC[2] += volW*MAT2D(i,2,dynamic,nPar);
      fC[3] += volW*MAT2D(i,3,dynamic,nPar);
/*......................................................*/
    }
/*...................................................................*/
  } 
/*...................................................................*/

/*...*/ 
  fC[0] /= volTotal;
  fC[1] /= volTotal;
  fC[2] /= volTotal;
  fC[3] /= volTotal;
/*...................................................................*/

/*...*/      
  cDyn[0] = fC[0]/fC[1];
/*...................................................................*/

/*...*/      
  cDyn[1] = fC[2]/fC[3];
/*...................................................................*/

/*... 0.0 < sqrt(c) < cap*/
  if(fCap){
    if(cDyn[0]  <= capCkFmin)
      cDyn[0] = capCkFmin;
    else if (cDyn[0] > capCkFmax)
      cDyn[0] = capCkFmax;

   if(cDyn[1]  <= capCeFmin)
      cDyn[1] = capCeFmin;
    else if (cDyn[1] > capCeFmax)
      cDyn[1] = capCeFmax;
  }
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 12/12/2017                                   *
 * Data de modificaco : 15/01/2018                                   * 
 *-------------------------------------------------------------------* 
 * oneParLes : calculo do coefciente do les de um paramento          * 
 *-------------------------------------------------------------------*
 * tModel   -> modelo de turbulencia                                 *
 * lViz     -> vizinhos                                              * 
 * dynamicl -> L:M, M:M e SVS nas celulas central e vizinhas         *
 * nFace    -> numero de faces da celula                             *
 * cap      -> true limita o range de variacao                       *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 * dynamic(*,3) = ( L:M, M:M, SVS)                                   *
 *********************************************************************/
DOUBLE oneParLes(Turbulence *tModel
               , INT *RESTRICT lViz, DOUBLE *RESTRICT dynamic
               , short const nFace , bool const fCap) {

  short i, idCell, type;
  INT vizNel;
  DOUBLE lm,mm,cf,capFmin,capFmax;

  type    = tModel->typeMixed[FUNMODEL];
  capFmin = 0.e0;
  capFmax = tModel->cf;

  idCell = nFace;
/*... cs =<M:L>/<M:M>*/
  lm = mm = 1.e-64;
  for (i = 0; i < nFace; i++) {
    vizNel = lViz[i];
    if(vizNel  > -1){
      lm += MAT2D(i,0,dynamic,3);
      mm += MAT2D(i,1,dynamic,3);
    }
  }
  lm += MAT2D(idCell,0,dynamic,3);
  mm += MAT2D(idCell,1,dynamic,3);      
/*...................................................................*/

/*...*/      
  cf = lm/mm;
/*...................................................................*/

/*...*/ 
  if (type == WALEMODEL) { 
    if( MAT2D(idCell,2,dynamic,3) < 0.09e0 )
      cf = capFmax*capFmax; 
  }
/*...................................................................*/

/*... 0.0 < sqrt(c) < cap*/
  if(fCap){
    if(cf <= capFmin)
      cf = capFmin;
    else if (sqrt(cf) > capFmax)
      cf = capFmax*capFmax;
  }
/*...................................................................*/

  return  cf;

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 13/01/2018                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * oneParLes : calculo do coefciente do les de um paramento (global) * 
 *-------------------------------------------------------------------*
 * tModel   -> modelo de turbulencia                                 *
 * gmm      -> L:M e M:M medias globais                              * 
 * cap      -> true limita o range de variacao                       *
 * ------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 * gm(2) = (L:M, M:M)                                                *
 *********************************************************************/
DOUBLE oneGoParLes(Turbulence *tModel , DOUBLE *gmm
                 , bool const fCap) {

  DOUBLE cf,capFmin,capFmax;

/*...*/      
  cf = gmm[0]/gmm[1];
/*...................................................................*/

/*... 0.0 < sqrt(c) < cap*/
  if(fCap){
    capFmin = 0.e0;
    capFmax = tModel->cf;
    if(cf <= capFmin)
      cf = capFmin; 
    else if (sqrt(cf) > capFmax)
      cf = capFmax*capFmax;
  }
/*...................................................................*/

  return  cf;

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 13/01/2018                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * oneParLesMod : calculo do coefciente do les de um paramento       *
 * (global)                                                          * 
 *-------------------------------------------------------------------*
 * tModel   -> modelo de turbulencia                                 *
 * gmm      -> L:M e M:M medias globais                              * 
 * density  -> densidade                                             * 
 * viscosity-> viscosidade dinamica                                  * 
 * cap      -> true limita o range de variacao                       *
 * ------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 * gm(4) = testFilter(Al:Al)                                         *
 *          testFilter(Al):testFilter(Al)                            *
 *          op*delta*tesFilter(S:S)                                  *
 *          opT*deltaT*tesFilter(S):tesFilter(S)                     *
 *********************************************************************/
DOUBLE oneGoParLesMod(Turbulence *tModel  , DOUBLE *gmm
                    , DOUBLE const density, DOUBLE const viscosity
                    , bool const fCap) {

  DOUBLE cf,capFmin,capFmax;

/*...*/      
  cf =-0.5e0*(viscosity/density)*(gmm[0] - gmm[1])/(gmm[2]-gmm[3]);
/*...................................................................*/

/*... 0.0 < sqrt(c) < cap*/
  if(fCap){
    capFmin = 0.e0;
    capFmax = tModel->cf;
    if(cf <= capFmin)
      cf = capFmin; 
    else if (sqrt(cf) > capFmax)
      cf = capFmax*capFmax;
  }
/*...................................................................*/

  return  cf;

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 12/12/2017                                   *
 * Data de modificaco : 15/01/2018                                   * 
 *-------------------------------------------------------------------* 
 * twoParLes : calculo do coefciente do les de um paramento          * 
 *-------------------------------------------------------------------* 
 * tModel   -> modelo de turbulencia                                 *
 * lViz     -> vizinhos                                              * 
 * dynamicl -> L:M e M:M nas celulas central e vizinhas              *
 * capFmin  -> valor minimo (Funcional)                              *        
 * capFmax  -> valor maximo (Funcional)                              *
 * capSmin  -> valor minimo (Estrutural)                             *
 * capSmax  -> valor maximo (Estrutural)                             *
 * nFace    -> numero de faces da celula                             *
 * cap      -> true limita o range de variacao                       *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * cDyn -> coeficiente dinamicos (cf,cs)                             *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 * dynamic = ( L:Ms, L:Mf, Ms:Ms, Ms:Mf, Mf:Ms, Mf:Mf )              *
 *                  | d  -b |   | a   b |                            *
 * inv(A) = 1/det(A)|       | ; |       |                            *
 *                  |-c   a |   | c   d |                            *
 *********************************************************************/
void twoParLes(Turbulence *tModel      , INT *RESTRICT lViz
             , DOUBLE *RESTRICT dynamic, DOUBLE *RESTRICT cDyn
             , short const nFace       , bool const fCap) {

  bool fUncoupled = false; 
  short i, idCell,count=1;
  INT vizNel;
  DOUBLE b[2],a[2][2],cf,cs,idet
        ,capFmin,capSmin,capFmax,capSmax;

  idCell = nFace;

  b[0] = b[1] = 0.e0;
  a[0][0] = a[0][1] = a[1][0] = a[1][1] = 0.e0;
/*... */
  for (i = 0; i < nFace; i++) {
    vizNel = lViz[i];
    if(vizNel  > -1){
      count++;
      b[0]    += MAT2D(i,0,dynamic,6);
      b[1]    += MAT2D(i,1,dynamic,6);
      a[0][0] += MAT2D(i,2,dynamic,6);
      a[0][1] += MAT2D(i,3,dynamic,6);
      a[1][0] += MAT2D(i,4,dynamic,6);
      a[1][1] += MAT2D(i,5,dynamic,6);
    }
  }
  b[0]    += MAT2D(idCell,0,dynamic,6);
  b[1]    += MAT2D(idCell,1,dynamic,6);
  a[0][0] += MAT2D(idCell,2,dynamic,6);
  a[0][1] += MAT2D(idCell,3,dynamic,6);
  a[1][0] += MAT2D(idCell,4,dynamic,6);
  a[1][1] += MAT2D(idCell,5,dynamic,6);      

  b[0]    /= count;
  b[1]    /= count;
  a[0][0] /= count;
  a[0][1] /= count;
  a[1][0] /= count;
  a[1][1] /= count; 
  if(fUncoupled) a[1][0]  = 0.e0;
/*...................................................................*/

/*... */ 
  idet = a[0][0]*a[1][1] - a[0][1]*a[1][0];
  idet = 1.0/idet;     
  cs = idet*(a[1][1]*b[0] - a[0][1]*b[1]);
  cf = idet*(a[0][0]*b[1] - a[1][0]*b[0]);
/*...................................................................*/
  
/*... */
  if(fCap){
    capFmin = capSmin = 0.0;
    capFmax = tModel->cf;
    capSmax = tModel->cs;
/*... capFmin < sqrt(c) < capFmax*/
    if(cf < capFmin)
      cf = capFmin;
    else if (sqrt(cf) > capFmax)
      cf = capFmax*capFmax;
/*...................................................................*/

/*... capSmin < sqrt(c) < capSmax*/      
    if(cs < capSmin)
      cs =  capSmin;
    else if (cs > capSmax)
      cs =  capSmax;
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  cDyn[0] = cf;
  cDyn[1] = cs;
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 19/10/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * DoubleDot : produto duplo para tensores                           * 
 *-------------------------------------------------------------------* 
 * t         -> tensor                                               * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 * | t11 t12 t13 |                                                   *  
 * | t21 t22 t23 | = (t11,t12,t13,t21,t22,t23,t31,t32,t33            *                                   *
 * | t31 t32 t33 |                                                   *
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE doubleDot(DOUBLE *t) {

  return  t[0]*t[0] + t[1]*t[1] + t[2]*t[2]
        + t[3]*t[3] + t[4]*t[4] + t[5]*t[5]
        + t[6]*t[6] + t[7]*t[7] + t[8]*t[8];

}
/*********************************************************************/


/********************************************************************* 
 * Data de criacao    : 19/10/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * DoubleDotSym : produto duplo para tensores simentricos            * 
 *-------------------------------------------------------------------* 
 * t         -> tensor (t11,t22,t33,t12,t13,t23)                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE doubleDotSym(DOUBLE *t) {

  return  t[0]*t[0] + t[1]*t[1] + t[2]*t[2]
         + 2.e0*( t[3]*t[3] + t[4]*t[4] + t[5]*t[5]);

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 12/12/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * DoubleDotSym2 : produto duplo para tensores simentricos           * 
 *-------------------------------------------------------------------* 
 * t         -> tensor (t11,t22,t33,t12,t13,t23)                     *
 * q         -> tensor (q11,q22,q33,q12,q13,q23)                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE doubleDotSym2(DOUBLE *RESTRICT t,DOUBLE *RESTRICT q) {

  return  t[0]*q[0] + t[1]*q[1] + t[2]*q[2]
         + 2.e0*( t[3]*q[3] + t[4]*q[4] + t[5]*q[5]);

}
/*********************************************************************/


/**********************************************************************
 * Data de criacao    : 07/12/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 * ------------------------------------------------------------------ *
 * sfHexa8 : funcao de interpolacao do hexaedro                       * 
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * eps - coordenada do ponto de integração                            *  
 * nn  - coordenada do ponto de integração                            * 
 * ze  - coordenada do ponto de integração                            * 
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * N  - Funçao de interpolação (ninter = true)                        * 
 * dn - derivada de N em relação a n                                  *        
 * de - derivada de N em relação a e                                  * 
 * -------------------------------------------------------------------* 
 *  OBS:                                                              *
 * no  ( r, s, t)                                                     *
 * no1 ( 1, 1, 1)                                                     * 
 * no2 (-1, 1, 1)                                                     * 
 * no3 (-1,-1, 1)                                                     * 
 * no4 ( 1,-1, 1)                                                     * 
 * no5 ( 1, 1,-1)                                                     * 
 * no6 (-1, 1,-1)                                                     * 
 * no7 (-1,-1,-1)                                                     * 
 * no8 ( 1,-1,-1)                                                     * 
 * -------------------------------------------------------------------* 
 **********************************************************************/
void sfHexa8(DOUBLE const eps ,DOUBLE const nn
           , DOUBLE const ze  
           , DOUBLE *RESTRICT N        ,DOUBLE *RESTRICT de
           , DOUBLE *RESTRICT dn       ,DOUBLE *RESTRICT dz
           , bool const ninter         ,bool const dev){


/*... Função de interpolação */
  if(ninter){
    N[4]  = ( 1.0+eps ) * ( 1.0+nn ) * (1.0+ze) * 0.125;
    N[5]  = ( 1.0-eps ) * ( 1.0+nn ) * (1.0+ze) * 0.125;
    N[6]  = ( 1.0-eps ) * ( 1.0-nn ) * (1.0+ze) * 0.125;
    N[7]  = ( 1.0+eps ) * ( 1.0-nn ) * (1.0+ze) * 0.125;
    N[0]  = ( 1.0+eps ) * ( 1.0+nn ) * (1.0-ze) * 0.125;
    N[1]  = ( 1.0-eps ) * ( 1.0+nn ) * (1.0-ze) * 0.125;
    N[2]  = ( 1.0-eps ) * ( 1.0-nn ) * (1.0-ze) * 0.125;
    N[3]  = ( 1.0+eps ) * ( 1.0-nn ) * (1.0-ze) * 0.125;
  }
/*.................................................................*/

/*... derivadas*/
  if(dev) {
/*... derivadas em relacao a e :*/

    de[4]  =   (1.0+nn) * (1.0+ze) * 0.125;
    de[5]  = - (1.0+nn) * (1.0+ze) * 0.125;
    de[6]  =   (nn-1.0) * (1.0+ze) * 0.125;
    de[7]  =   (1.0-nn) * (1.0+ze) * 0.125;
    de[0]  =   (1.0+nn) * (1.0-ze) * 0.125;
    de[1]  =   (1.0+nn) * (ze-1.0) * 0.125;
    de[2]  =   (1.0-nn) * (ze-1.0) * 0.125;
    de[3]  =   (1.0-nn) * (1.0-ze) * 0.125;

/*... derivadas em relacao a n :*/

    dn[4]  =  (1.0+eps) * (1.0+ze)* 0.125;
    dn[5]  =  (1.0-eps) * (1.0+ze)* 0.125;
    dn[6]  =  (eps-1.0) * (1.0+ze)* 0.125;
    dn[7]  = -(1.0+eps) * (1.0+ze)* 0.125;
    dn[0]  =  (1.0+eps) * (1.0-ze)* 0.125;
    dn[1]  =  (1.0-eps) * (1.0-ze)* 0.125;
    dn[2]  =  (1.0-eps) * (ze-1.0)* 0.125;
    dn[3]  =  (1.0+eps) * (ze-1.0)* 0.125;

/*... derivadas em relacao a z   :*/

    dz[4]  =  ( 1.0+eps ) * ( 1.0+nn )* 0.125;
    dz[5]  =  ( 1.0-eps ) * ( 1.0+nn )* 0.125;
    dz[6]  =  ( 1.0-eps ) * ( 1.0-nn )* 0.125;
    dz[7]  =  ( 1.0+eps ) * ( 1.0-nn )* 0.125;
    dz[0]  = -( 1.0+eps ) * ( 1.0+nn )* 0.125;
    dz[1]  =  (eps -1.0 ) * ( 1.0+nn )* 0.125;
    dz[2]  =  ( 1.0-eps ) * ( nn -1.0)* 0.125;
    dz[3]  =  ( 1.0+eps ) * ( nn -1.0)* 0.125;
  }
/*..................................................................*/  
} 
/********************************************************************/ 


/*********************************************************************
 * Data de criacao    : 07/12/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ----------------------------------------------------------------- *
 * JACAB3D: Calcula a matriz jacobiana inversa e o det do            *   
 * ponto (e,n)                                                       * 
 * e as derivada da em função de intepolação em relação a x,y,z      *   
 * ----------------------------------------------------------------- *
 * parametros de entrada:                                            * 
 * ----------------------------------------------------------------- *
 * xl  - coordenadas                                                 *
 * de - derivadas das funcoes de interpolacao em relacao a e         *
 * dn - derivadas das funcaem de interpolacao em relacao a n         *  
 * dz - derivadas das funcoes de interpolacao em relacao a e         *
 * nel - numero do elemneto                                          *
 * nen - numero de pontos                                            *  
 * alf - calculo do jacobiano                                        * 
 * dev - derivadas das funcoes de interplocao                        *
 * nel - numero do elemento                                          *
 * ----------------------------------------------------------------- *
 * parametros de saida:                                              * 
 * ----------------------------------------------------------------- *
 * hx - derivada da funcao de interpolacap em relacao a x            *
 * hy - derivada da funcao de interpolacap em relacao a y            *
 * hz - derivada da funcao de interpolacap em relacao a y            *
 * det - determinante da matriz jacobiana normal                     *
 * -------------------------------------------------------------------* 
 *  OBS:                                                              *
 * -------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE jacob3d(DOUBLE *RESTRICT xl, DOUBLE *RESTRICT de
           , DOUBLE *RESTRICT dn  , DOUBLE *RESTRICT dz 
           , DOUBLE *RESTRICT hx  , DOUBLE *RESTRICT hy
           , DOUBLE *RESTRICT hz   
 	         , short const nen      , bool const afl
           , bool const dev       , INT const nel ){

  double jac[3][3],jaci[3][3],DET,xx;
  long i,j;

  DET = 1.e0;
 
/*... Calculo da matriz jacobiana */
  if(afl){
    for(i=0;i<3;i++){
      jac[0][i] = 0.0;
      jac[1][i] = 0.0;
      jac[2][i] = 0.0;
      for(j=0;j<nen;j++){
        xx         = MAT2D(j,i,xl,3); 
        jac[0][i] += de[j] * xx;
        jac[1][i] += dn[j] * xx;
        jac[2][i] += dz[j] * xx;
      }
    }
/*....................................................................*/ 

/*... determinante da matriz jacobiana:*/ 
    DET = jac[0][0]*jac[1][1]*jac[2][2] + jac[0][1]*jac[1][2]*jac[2][0] 
        + jac[0][2]*jac[1][0]*jac[2][1] - jac[2][0]*jac[1][1]*jac[0][2] 
        - jac[0][1]*jac[1][0]*jac[2][2] - jac[0][0]*jac[2][1]*jac[1][2];
 
    if (DET <= 0.e0)   {
      printf("Determinante nulo ou negativo do elemento ->  %d ",nel); 
      exit(EXIT_FAILURE); 
    }
/*...................................................................*/
  }  
/*...................................................................*/

/*... Inversa da matriz Jacobiana:*/  
  if(dev){
    jaci[0][0] =  ( jac[1][1] * jac[2][2] 
               -    jac[1][2] * jac[2][1] ) / DET;
    jaci[1][0] = -( jac[1][0] * jac[2][2] 
               -    jac[1][2] * jac[2][0] ) / DET;
    jaci[2][0] =  ( jac[1][0] * jac[2][1] 
               -     jac[1][1] * jac[2][0] ) / DET;
    jaci[0][1] = -( jac[0][1] * jac[2][2] 
               -     jac[0][2] * jac[2][1] ) / DET;
    jaci[1][1] =  ( jac[0][0] * jac[2][2] 
               -    jac[0][2] * jac[2][0] ) / DET;
    jaci[2][1] = -( jac[0][0] * jac[2][1] 
               -    jac[0][1] * jac[2][0] ) / DET;
    jaci[0][2] =  ( jac[0][1] * jac[1][2] 
               -    jac[0][2] * jac[1][1] ) / DET;
    jaci[1][2] = -( jac[0][0] * jac[1][2] 
               -    jac[0][2] * jac[1][0] ) / DET;
    jaci[2][2] =  ( jac[0][0] * jac[1][1] 
               -    jac[0][1] * jac[1][0] ) / DET;
/*...................................................................*/
 

/*... Derivadas(x,y,z) das funcoes de interpolacao:*/
    for(i=0;i<nen;i++){   
      hx[i] = jaci[0][0]*de[i] + jaci[1][0]*dn[i] 
            + jaci[2][0]*dz[i];
      hy[i] = jaci[0][1]*de[i] + jaci[1][1]*dn[i] 
            + jaci[2][1]*dz[i];
      hz[i] = jaci[0][2]*de[i] + jaci[1][2]*dn[i] 
            + jaci[2][2]*dz[i];
    }
  }
/*....................................................................*/
  return DET;
}
/*********************************************************************/   

/********************************************************************* 
 * Data de criacao    : 27/11/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * eigenValue3x3 : autovalores de uma matrix 3x3 simetrica           * 
 *-------------------------------------------------------------------* 
 * g -> matriz 3x3                                                   *
 * s -> nao definido                                                 *  
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * s -> autovalores s1>s2>s3                                         *  
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void eigenValue3x3(DOUBLE *RESTRICT g, DOUBLE *RESTRICT s) {

  DOUBLE i1,i2,i3,a1,a2,a3,tmp;

/*... primeiro invariante*/
  i1 =  MAT2D(0,0,g,3) + MAT2D(1,1,g,3) + MAT2D(2,2,g,3); 
/*... segundo invariante*/
  tmp = MAT2D(0,0,g,3)*MAT2D(0,0,g,3)   
      + MAT2D(0,1,g,3)*MAT2D(1,0,g,3) 
      + MAT2D(0,2,g,3)*MAT2D(2,0,g,3) 
/*...*/
      + MAT2D(1,0,g,3)*MAT2D(0,1,g,3) 
      + MAT2D(1,1,g,3)*MAT2D(1,1,g,3) 
      + MAT2D(1,2,g,3)*MAT2D(2,1,g,3) 
/*...*/
      + MAT2D(2,0,g,3)*MAT2D(0,2,g,3) 
      + MAT2D(2,1,g,3)*MAT2D(1,2,g,3) 
      + MAT2D(2,2,g,3)*MAT2D(2,2,g,3);
  i2 = 0.5*(i1*i1 - tmp); 
/*... terceiro incariante*/
  i3 = MAT2D(0,0,g,3)*MAT2D(1,1,g,3)*MAT2D(2,2,g,3)  
     + MAT2D(0,1,g,3)*MAT2D(1,2,g,3)*MAT2D(2,0,g,3)  
     + MAT2D(0,2,g,3)*MAT2D(1,0,g,3)*MAT2D(2,1,g,3) 
/*...*/ 
     - MAT2D(2,0,g,3)*MAT2D(1,1,g,3)*MAT2D(0,2,g,3)  
     - MAT2D(2,1,g,3)*MAT2D(1,2,g,3)*MAT2D(0,0,g,3) 
     - MAT2D(2,2,g,3)*MAT2D(1,0,g,3)*MAT2D(0,1,g,3);  
/*...*/
  a1 = i1*i1 - 3.e0*i2;
  a2 = 2.e0*i1*i1*i1 - 9.e0*i1*i2 + 27.e0*i3;
/* ... problema de precisao que gera valore ligeiramente superiores a 1*/ 
  tmp = 0.5*a2/pow(a1,1.5);
  if(tmp > 1.0) 
    tmp = 1.e0;
  else if( tmp < -1.0)
    tmp = -1.e0;

  a3=D1DIV3*acos(tmp);
/* ..................................................................*/

/*...*/
  s[0] = D1DIV3*i1 + D2DIV3*sqrt(a1)*cos(a3);
  s[1] = D1DIV3*i1 + D2DIV3*sqrt(a1)*cos(a3-D2DIV3*PI);
  s[2] = D1DIV3*i1 + D2DIV3*sqrt(a1)*cos(a3-2.e0*D2DIV3*PI);

  printf("Matrix:\n");
  printf("%e %e %e\n%e %e %e\n%e %e %e\n"
        ,MAT2D(0,0,g,3),MAT2D(0,1,g,3),MAT2D(0,2,g,3)
        ,MAT2D(1,0,g,3),MAT2D(1,1,g,3),MAT2D(1,2,g,3)
        ,MAT2D(2,0,g,3),MAT2D(2,1,g,3),MAT2D(2,2,g,3));

  printf("%e %e %e %e\n",s[0],s[1],s[2],a3);
}
/*********************************************************************/ 
