#include<Turbulence.h>

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
static DOUBLE doubleDotSym(DOUBLE *t) {

  return  t[0]*t[0] + t[1]*t[1] + t[2]*t[2]
         + 2.e0*( t[3]*t[3] + t[4]*t[4] + t[5]*t[5]);

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 30/11/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * TURBULENCE: Calculo da viscosidae turbulenta                      *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m     -> memoria principal                                        *
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
 * mKsi     -> modulo do vetor ksi                                   * 
 * eta      -> vetores paralelos as faces das celulas                * 
 * fArea    -> modulo do vetor eta                                   * 
 * normal   -> vetores normais as faces das celulas                  * 
 * volume   -> volumes das celulas                                   * 
 * xm       -> pontos medios das faces das celulas                   * 
 * xmcc     -> vetores que unem o centroide aos pontos medios das    * 
 *              faces                                                * 
 * vSkew    -> vetor entre o ponto medio a intersecao que une os     * 
 *              centrois compartilhado nessa face                    * 
 * mvSkew   -> distacia entre o ponto medio a intersecao que une os  * 
 *              centrois compartilhado nessa face                    * 
 * dcca     -> menor distancia do centroide a faces desta celula     * 
 * gradVel   -> gradiente da solucao conhecido                       * 
 * density   -> massa especifica com variacao temporal               *  
 * eddyViscosity-> viscosidade turbulenta                            *
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
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void turbulence(Memoria *m
      , Loads *lVel                 , Turbulence tModel             
      , INT    *RESTRICT el                 , INT    *RESTRICT nelcon 
      , short  *RESTRICT nen                , short  *RESTRICT nFace 
      , short  *RESTRICT geomType           , DOUBLE *RESTRICT prop  
      , short  *RESTRICT calType            , short  *RESTRICT mat        
      , DOUBLE *RESTRICT cc                 , DOUBLE *RESTRICT ksi                      
      , DOUBLE *RESTRICT mKsi               , DOUBLE *RESTRICT eta  
      , DOUBLE *RESTRICT fArea              , DOUBLE *RESTRICT normal  
      , DOUBLE *RESTRICT volume             , DOUBLE *RESTRICT xm  
      , DOUBLE *RESTRICT xmcc               , DOUBLE *RESTRICT vSkew    
      , DOUBLE *RESTRICT mvSkew             , DOUBLE *RESTRICT dcca  
      , short  *RESTRICT faceRvel           , short *RESTRICT faceVelL               
      , DOUBLE *RESTRICT vel                , DOUBLE *RESTRICT gradVel
      , DOUBLE *RESTRICT densityFluid       , DOUBLE *RESTRICT dViscosity        
      , DOUBLE *RESTRICT eddyViscosity      , DOUBLE *RESTRICT wallPar                        
      , short const maxNo                   , short const maxViz
      , short const ndm                     , INT const numel     
      , short const ndf)                      
{   
  bool flag = false;
  short type = tModel.type;
  DOUBLE *dynamic=NULL;

  HccaAlloc(DOUBLE,m,dynamic, numel*2,"dynamic"  ,_AD_); 
/*...*/
  switch (type) {  
    case DYNAMIC:
/*...*/ 
      lesDynamicMean(nelcon,nen          , nFace 
                    , volume      , vel   
                    , gradVel     , densityFluid 
                    , dynamic
                    , maxNo        , maxViz         
                    , ndm          , numel 
                    , ndf); 
/*...................................................................*/ 

    case SMAGORINSKY:
    case WALEMODEL:
    case VREMAN: 
    case SIGMAMODEL: 
/*...*/ 
      turbulenceCellLoop(lVel         , tModel               
                   , el           , nelcon 
                   , nen          , nFace 
                   , geomType     , prop 
                   , calType      , mat 
                   , cc           , ksi 
                   , mKsi         , eta     
                   , fArea        , normal 
                   , volume       , xm     
                   , xmcc         , vSkew    
                   , mvSkew       , dcca 
                   , faceRvel     , faceVelL     
                   , vel          , gradVel      
                   , densityFluid , dViscosity
                   , eddyViscosity, wallPar
                   , dynamic   
                   , maxNo        , maxViz         
                   , ndm          , numel 
                   , ndf); 
/*...................................................................*/ 
      break;
/*...................................................................*/ 

/*...*/
    default: 
      ERRO_OP(__FILE__,__func__,type);
/*...................................................................*/ 
  }

/*...*/
  HccaDealloc(m,dynamic,"dynamic",_AD_);
/*...................................................................*/ 

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 11/09/2017                                   *
 * Data de modificaco : 20/11/2017                                   * 
 *-------------------------------------------------------------------* 
 * TURBULENCE: Calculo da viscosidae turbulenta                      *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
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
 * gradVel   -> gradiente da solucao conhecido                       * 
 * density   -> massa especifica com variacao temporal               *  
 * eddyViscosity-> viscosidade turbulenta                            *
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
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
  *********************************************************************/
void turbulenceCellLoop(Loads *lVel         , Turbulence tModel             
      , INT    *RESTRICT el                 , INT    *RESTRICT nelcon 
      , short  *RESTRICT nen                , short  *RESTRICT nFace 
      , short  *RESTRICT geomType           , DOUBLE *RESTRICT prop  
      , short  *RESTRICT calType            , short  *RESTRICT mat        
      , DOUBLE *RESTRICT cc                 , DOUBLE *RESTRICT gKsi                      
      , DOUBLE *RESTRICT gmKsi              , DOUBLE *RESTRICT gEta  
      , DOUBLE *RESTRICT gfArea             , DOUBLE *RESTRICT gNormal  
      , DOUBLE *RESTRICT gVolume            , DOUBLE *RESTRICT gXm  
      , DOUBLE *RESTRICT gXmcc              , DOUBLE *RESTRICT gvSkew    
      , DOUBLE *RESTRICT gmvSkew            , DOUBLE *RESTRICT gDcca  
      , short  *RESTRICT faceVelR           , short *RESTRICT faceVelL               
      , DOUBLE *RESTRICT vel                , DOUBLE *RESTRICT gradVel
      , DOUBLE *RESTRICT density            , DOUBLE *RESTRICT dViscosity        
      , DOUBLE *RESTRICT eddyViscosity      , DOUBLE *RESTRICT wallPar  
      , DOUBLE *RESTRICT dynamic                 
      , short const maxNo                   , short const maxViz
      , short const ndm                     , INT const numel     
      , short const ndf)                      
{   
  short i,j,k;
  INT nel,vizNel;
  
/*... variavel local */
  short  aux1,aux2,lMat,lib;
  short  lGeomType[MAX_NUM_FACE+1],  
         lFaceVelR[MAX_NUM_FACE+1],lFaceVelL[MAX_NUM_FACE+1];       
  INT    lViz[MAX_NUM_FACE];
  DOUBLE lKsi[MAX_NUM_FACE*MAX_NDM], lmKsi[MAX_NUM_FACE],
         lEta[MAX_NUM_FACE*MAX_NDM], lfArea[MAX_NUM_FACE],
         lNormal[MAX_NUM_FACE*MAX_NDM], lVolume[MAX_NUM_FACE+1],
         lXm[MAX_NUM_FACE*MAX_NDM], lXmcc[MAX_NUM_FACE*MAX_NDM],
         lDcca[MAX_NUM_FACE], lmvSkew[MAX_NUM_FACE],
         lvSkew[MAX_NUM_FACE*MAX_NDM], lDensity[(MAX_NUM_FACE+1)],
         lProp[(MAX_NUM_FACE+1)*MAXPROP],
         lCc[(MAX_NUM_FACE+1)*MAX_NDM],
         lGradVel[(MAX_NUM_FACE+1)*MAX_NDM*MAX_NDF],
         lVel0[(MAX_NUM_FACE+1)*MAX_NDM],lDyn[(MAX_NUM_FACE+1)*2],
         lEddyViscosity,lDviscosity,lWallPar[4];
           
/*... loop nas celulas*/
  aux2    = maxViz+1;
  for(nel=0;nel<numel;nel++){
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
    MAT2D(aux1,0,lDyn,2) = MAT2D(nel,0,dynamic,2); 
    MAT2D(aux1,1,lDyn,2) = MAT2D(nel,1,dynamic,2); 
/*...*/
    lDviscosity       = dViscosity[nel];
/*...*/      
    for(j=0;j<MAXPROP;j++)
      MAT2D(aux1,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
/*...................................................................*/

/*...*/
    for(j=0;j<ndm;j++){
       MAT2D(aux1,j,lVel0,ndm)  = MAT2D(nel,j,vel     ,ndm);
       MAT2D(aux1, j, lCc, ndm) = MAT2D(nel, j, cc, ndm);
    }
/*...................................................................*/

/*...*/
    for(i=0;i<ndf;i++)
      for(j=0;j<ndm;j++)
        MAT3D(aux1, i, j, lGradVel, ndf, ndm) 
                = MAT3D(nel, i, j, gradVel, ndf, ndm);
/*...................................................................*/

/*...*/
    for(i=0;i<aux1;i++){
      lmKsi[i]      = MAT2D(nel,i,gmKsi   ,maxViz);
      lfArea[i]     = MAT2D(nel,i,gfArea  ,maxViz);
      lDcca[i]      = MAT2D(nel,i,gDcca   ,maxViz);
      lmvSkew[i]    = MAT2D(nel,i,gmvSkew ,maxViz);
      lFaceVelR[i]  = MAT2D(nel,i,faceVelR,aux2);
      lFaceVelL[i]  = MAT2D(nel,i,faceVelL,aux2);
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
        lDensity[i]  = MAT2D(vizNel,2   ,density ,DENSITY_LEVEL);
        MAT2D(i,0,lDyn,2) = MAT2D(vizNel,0,dynamic,2); 
        MAT2D(i,1,lDyn,2) = MAT2D(vizNel,1,dynamic,2); 

        lMat         = mat[vizNel]-1;
 
        for(j=0;j<ndm;j++){
          MAT2D(i,j,lCc      ,ndm)  = MAT2D(vizNel,j,cc      ,ndm);
          MAT2D(i,j,lVel0     ,ndm) = MAT2D(vizNel,j,vel     ,ndm);
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

/*... chamando a biblioteca de celulas*/
     cellLibTurbulence(lVel            , tModel        
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
                      , lWallPar       , lDyn
                      , nen[nel]       , nFace[nel] 
                      , ndm            , lib   
                      , nel);    
/*...................................................................*/
 
/*...*/
    eddyViscosity[nel]  = lEddyViscosity;
    MAT2D(nel,0,wallPar,4) = lWallPar[0];
    MAT2D(nel,1,wallPar,4) = lWallPar[1];
    MAT2D(nel,2,wallPar,4) = lWallPar[2];
    MAT2D(nel,3,wallPar,4) = lWallPar[3];

  }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 30/11/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * lesDynamicMean : cacula do M:L e M:M separados para fazer a media *
 * do LES dinamico                                                   *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * nelcon    -> vizinhos dos elementos                               * 
 * nen       -> numero de nos por celulas                            * 
 * nFace     -> numero de faces por celulas                          * 
 * gradVel   -> gradiente da solucao conhecido                       * 
 * density   -> massa especifica com variacao temporal               * 
 * dynamic   -> nao definido                                         * 
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * numel     -> numero de toral de celulas                           * 
 * ndf       -> graus de liberdade                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * dynamic  -> retorna o produto M:L e M:M                           *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void lesDynamicMean(INT    *RESTRICT nelcon  , short  *RESTRICT nen    
                  , short  *RESTRICT nFace   , DOUBLE *RESTRICT gVolume 
                  , DOUBLE *RESTRICT vel     , DOUBLE *RESTRICT gradVel 
                  , DOUBLE *RESTRICT density , DOUBLE *RESTRICT dynamic                                  
                  , short const maxNo        , short const maxViz
                  , short const ndm          , INT const numel     
                  , short const ndf)                      
{   
  short i,j,k;
  INT nel,vizNel;
  
/*... variavel local */
  short  aux1,aux2,lMat,lib;
  INT    lViz[MAX_NUM_FACE];
  DOUBLE lVolume[MAX_NUM_FACE+1],lDensity[(MAX_NUM_FACE+1)],
         lGradVel[(MAX_NUM_FACE+1)*MAX_NDM*MAX_NDF],
         lVel[(MAX_NUM_FACE+1)*MAX_NDM],lDynamic[2];
           
/*... loop nas celulas*/
  aux2    = maxViz+1;
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
    lesDynamic(lViz     , lVolume
              , lDensity , lVel
              , lGradVel , lDynamic 
              , nFace[nel]);
/*...................................................................*/

/*...*/
    MAT2D(nel,0,dynamic,2) = lDynamic[0];
    MAT2D(nel,1,dynamic,2) = lDynamic[1];
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
  DOUBLE wt,velC[2],vParallel[2],lNormal[2],lMin,dMin
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
          getLoads(par,loadsVel[nCarg]);
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
 * Data de modificaco : 30/11/2017                                   *
 *-------------------------------------------------------------------*
 * celLes3D: Calculo da viscosidae turbulenta para o LES             *
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
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
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
void cellLes3D(Loads *lVel             , Turbulence tModel           
          , short *RESTRICT lGeomType  , DOUBLE *RESTRICT prop 
          , INT *RESTRICT lViz         , DOUBLE *RESTRICT fArea
          , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume 
          , DOUBLE *RESTRICT dcca      , short *RESTRICT lFaceVelR  
          , short *RESTRICT lFaceVelL  , DOUBLE *RESTRICT vel       
          , DOUBLE *RESTRICT gradVel   , DOUBLE *RESTRICT lDensity  
          , DOUBLE const dViscosity    , DOUBLE *viscosity          
          , DOUBLE *RESTRICT wallPar   , DOUBLE *RESTRICT dynamic
          , const short nEn            , short const nFace 
          , const short ndm            , INT const nel) 
{
/*...*/
  bool fWall = false;
  short i, j, nf, nCarg, idCell,type,wallType;
  INT vizNel;
  DOUBLE modS, tmp, density, viscosityC, cs, s[6],delta,*iGradVel;
  DOUBLE wt,v[3],vParallel[3],lNormal[3],lMin,dMin, g[3][3]
        ,yPlus  ,uPlus,velB[3],yPlusMax, modSd, b[3][3];
  DOUBLE m[6],l[6],mm,lm,volTotal;

  idCell = nFace;
  type=tModel.type;

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

/*... wall model*/
  wallType   = tModel.wallType;
  viscosityC = dViscosity;
  density    = lDensity[idCell];
/*...*/
  v[0] = MAT2D(idCell, 0, vel, 3);
  v[1] = MAT2D(idCell, 1, vel, 3);
  v[2] = MAT2D(idCell, 2, vel, 3);
/*...................................................................*/

  fWall = wallDist(lVel               
                  ,lViz       ,v
                  ,normal     ,dcca
                  ,lFaceVelR  ,lFaceVelL
                  ,viscosityC ,density
                  ,wallPar    ,&dMin
                  ,wallType   ,nFace );
   
  yPlusMax = wallPar[0];

/*...*/
  switch (type) {
    case SMAGORINSKY:
/*...*/
      cs         = tModel.cs;
      density    = lDensity[idCell];
/*...................................................................*/

/*.. calculo Sij*/
      s[0] = g[0][0];
      s[1] = g[1][1];
      s[2] = g[2][2];
      s[3] = 0.5e0*(g[0][1] + g[1][0]); /*s12*/
      s[4] = 0.5e0*(g[0][2] + g[2][0]); /*s13*/
      s[5] = 0.5e0*(g[1][2] + g[2][1]); /*s23*/
/*... |S| = sqrt(2S:S)*/
      modS = doubleDotSym(s);
      modS = sqrt(2.e0*modS);
/*...................................................................*/
 
/*...*/
      delta = pow(volume[idCell],D1DIV3);
      lMin = cs*delta;      
      if (fWall) {
        tmp  = 1.e0-exp(-yPlusMax/VANDRIEST); 
        lMin = min(VONKARMAN*dMin,tmp*lMin); 
      }
/*...................................................................*/
      *viscosity = density*lMin*lMin*modS;  
/*...................................................................*/
      break;
/*...................................................................*/

/*...*/
    case WALEMODEL: 

/*...*/
      cs         = tModel.cs;
      density    = lDensity[idCell];
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
      tmp = D1DIV3*(b[0][0] + b[1][1] + b[2][2]);
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
      mm  = pow(modS,2.5e0) +  pow(modSd,1.25e0);
      lm  = pow(modSd,1.5e0);
      tmp = lm/mm; 
      delta = pow(volume[idCell],D1DIV3);
/*...*/
      if(fWall)
        lMin = min(VONKARMAN*dMin,cs*delta);
      else
        lMin  = cs*delta;
/*...................................................................*/
      *viscosity = density*lMin*lMin*tmp;  
/*...................................................................*/
        
      break;
/*...................................................................*/

/*...*/
    case VREMAN: 

/*...*/
      cs       = 2.5e0*tModel.cs*tModel.cs;
      density  = lDensity[idCell];
      delta    = pow(volume[idCell],D1DIV3);
/*...................................................................*/

/*... beta*/
      tmp = delta*delta;
      for(i=0;i<3;i++)
        for(j=0;j<3;j++)
          b[i][j] = tmp*( g[0][i]*g[0][j] 
                        + g[1][i]*g[1][j]
                        + g[2][i]*g[2][j]);
/*...*/
      tmp = b[0][0]*b[1][1] - b[0][1]*b[0][1] + b[0][0]*b[2][2]
          - b[0][2]*b[0][2] + b[1][1]*b[2][2] - b[1][2]*b[1][2];
/*... |S| = |2S:S|*/
      modS = 1.e-32;
      for(i=0;i<3;i++)
        for(j=0;j<3;j++)
          modS += g[i][j]*g[i][j];
/*...................................................................*/
      
/*...*/  
      if ( tmp < 0.e0 ) tmp = 0.e0;
      tmp = sqrt(tmp/modS);      
      *viscosity = density*cs*tmp;     
/*...................................................................*/
        
      break;
/*...................................................................*/

/*...*/
    case DYNAMIC: 

/*... cs =(1/2)<M:L>/<M:M>*/
      lm = mm = 1.e-64;
      for (i = 0; i < nFace; i++) {
        vizNel = lViz[i];
        if(vizNel  > -1){
          lm += MAT2D(i,0,dynamic,2);
          mm += MAT2D(i,1,dynamic,2);
        }
      }
      lm += MAT2D(idCell,0,dynamic,2);
      mm += MAT2D(idCell,1,dynamic,2);      
/*...................................................................*/

/*... 0.0 < sqrt(c) < 0.23*/      
      cs = 0.5e0*lm/mm;
      if(cs < 0.0e0)
        cs = 0.e0;
      else if (sqrt(cs) > DYNAMIC_CLIP)
        cs = DYNAMIC_CLIP*DYNAMIC_CLIP;
/*...................................................................*/
      
/*...*/
      density  = lDensity[idCell];
      delta    = pow(volume[idCell],D2DIV3);
/*...................................................................*/

/*... (cs*delta)^2*/
      lMin= cs*delta;
/*...................................................................*/

/*.. calculo Sij*/
      s[0] = g[0][0];
      s[1] = g[1][1];
      s[2] = g[2][2];
      s[3] = 0.5e0*(g[0][1] + g[1][0]); /*s12*/
      s[4] = 0.5e0*(g[0][2] + g[2][0]); /*s13*/
      s[5] = 0.5e0*(g[1][2] + g[2][1]); /*s23*/
/*... |S| = |2S:S|*/
      modS = doubleDotSym(s);
      modS = sqrt(2.e0*modS);
/*...................................................................*/

/*...*/   
      tmp      = density*lMin*modS;
       *viscosity = tmp; 
/*...................................................................*/
      break;
/*...................................................................*/

/*...*/
    case SIGMAMODEL: 
/*...*/
      cs       = tModel.cs*tModel.cs;
      density  = lDensity[idCell];
      delta    = pow(volume[idCell],D2DIV3);
/*...................................................................*/

/*... (cs*delta)^2*/
      lMin = cs*delta;
/*...................................................................*/

/*... operador(grad(u))*/
      iGradVel = &MAT3D(idCell,0,0,gradVel,3,3);
      tmp = sigmaModel(s    ,iGradVel 
                      ,nFace,ndm);
/*...................................................................*/

/*...*/   
       *viscosity = density*lMin*tmp; 
/*...................................................................*/
      break;
/*...................................................................*/

/*...*/
    default: 
      ERRO_OP(__FILE__,__func__,type);
/*...................................................................*/
  }
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 21/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
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
  DOUBLE stressW,f,df,tol = 1.e-7, nu = viscosity/density;
  DOUBLE yPlus,uPlus,fu,fu0, conv;
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
      yPlus    = fu*dWall/nu;      
      if ( yPlus < 0.e0 )
        ERRO_GERAL(__FILE__,__func__,__LINE__,"WallModel: y+ < 0");   
/*...................................................................*/        
      break;
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
 * Data de criacao    : 24/10/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * LESDYNAMIC: Metodo dinamico para o calculo da constente do LES    *
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
void lesDynamic(INT *RESTRICT lViz       , DOUBLE *RESTRICT volume
              , DOUBLE *RESTRICT lDensity, DOUBLE *restrict vel
              , DOUBLE *restrict gradVel , DOUBLE *RESTRICT lDynamic
              , short const nFace  ) {

  short i,idCell;
  INT vizNel;
  DOUBLE tmp,s[6],modS,cs;
  DOUBLE v[3],m[6],l[6],mm,lm,deltaT,delta,volW,volTotal,density,g[3][3];
  DOUBLE tFilterDenModSs[6],tFilterMods[6],tFilterDenVv[6];
  DOUBLE tFilterModS,tFilterDen,tFilterDenV[3],tFilterS[6];

  idCell = nFace;

  tFilterDenV[0] = tFilterDenV[1] = tFilterDenV[2] = 0.e0;
  tFilterDen     = tFilterModS = volTotal = 0.e0;

  for(i = 0; i<6 ; i++ ){
    tFilterDenModSs[i] = 0.e0;
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
      tFilterDenModSs[0] += tmp*s[0];  /*modS*s11*/
      tFilterDenModSs[1] += tmp*s[1];  /*modS*s22*/
      tFilterDenModSs[2] += tmp*s[2];  /*modS*s33*/
      tFilterDenModSs[3] += tmp*s[3];  /*modS*s12*/
      tFilterDenModSs[4] += tmp*s[4];  /*modS*s13*/
      tFilterDenModSs[5] += tmp*s[5];  /*modS*s23*/
/*... tesFilter( |s|) -> vol*modS */
      tFilterModS += volW*modS;
/*... tesFilter(s) -> vol*s */
      tFilterS[0] += volW*s[0];     /*s11*/
      tFilterS[1] += volW*s[1];     /*s22*/
      tFilterS[2] += volW*s[2];     /*s33*/
      tFilterS[3] += volW*s[3];     /*s12*/
      tFilterS[4] += volW*s[4];     /*s13*/
      tFilterS[5] += volW*s[5];     /*s23*/
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

/*tesFilter( density*|s|*s ) */
  tFilterDenModSs[0] /= volTotal;
  tFilterDenModSs[1] /= volTotal;
  tFilterDenModSs[2] /= volTotal;
  tFilterDenModSs[3] /= volTotal;
  tFilterDenModSs[4] /= volTotal;
  tFilterDenModSs[5] /= volTotal;
/*... tesFilter( |s|) -> vol*modS */
  tFilterModS /= volTotal;
/*... tesFilter(s) -> vol*s */
  tFilterS[0] /= volTotal;
  tFilterS[1] /= volTotal;
  tFilterS[2] /= volTotal;
  tFilterS[3] /= volTotal;
  tFilterS[4] /= volTotal;
  tFilterS[5] /= volTotal;
/*... tesFilter(density*vv)*/
  tFilterDenVv[0] /= volTotal; /*v1v1*/
  tFilterDenVv[1] /= volTotal; /*v2v2*/
  tFilterDenVv[2] /= volTotal; /*v3v3*/
  tFilterDenVv[3] /= volTotal; /*v1v2*/
  tFilterDenVv[4] /= volTotal; /*v1v3*/
  tFilterDenVv[5] /= volTotal; /*v2v3*/
/*... tesFilter(density*v)*/
  tFilterDenV[0] /= volTotal; 
  tFilterDenV[1] /= volTotal; 
  tFilterDenV[2] /= volTotal;
/*... tesFilter(desinty)*/
  tFilterDen     /= volTotal;
      
/*... (filtros)^2*/
  delta    = pow(volume[idCell],D2DIV3);
  deltaT   = pow(volTotal      ,D2DIV3);

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
    
/*... m*/
  tmp  = deltaT*tFilterDen*tFilterModS;
  m[0] = delta*tFilterDenModSs[0] - tmp*tFilterS[0];  /*l11*/
  m[1] = delta*tFilterDenModSs[1] - tmp*tFilterS[1];  /*l22*/
  m[2] = delta*tFilterDenModSs[2] - tmp*tFilterS[2];  /*l33*/
  m[3] = delta*tFilterDenModSs[3] - tmp*tFilterS[3];  /*l1v2*/
  m[4] = delta*tFilterDenModSs[4] - tmp*tFilterS[4];  /*l1v3*/
  m[5] = delta*tFilterDenModSs[5] - tmp*tFilterS[5];  /*l2v3*/

/*... LijMij*/
  lm = l[0]*m[0] + l[1]*m[1] + l[2]*m[2]
     + 2.e0*( l[3]*m[3] + l[4]*m[4] + l[5]*m[5]);
/*... MijMij*/
  mm = doubleDotSym(m);      

/*...*/
  lDynamic[0] = lm;
  lDynamic[1] = mm;
/*...................................................................*/

}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/11/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * wallDist  : calcula a distancia a parede e a distancia normalizada*
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
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri)          * 
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
bool  wallDist(Loads *lVel               
             , INT *RESTRICT lViz       , DOUBLE *RESTRICT v
             , DOUBLE *RESTRICT normal  , DOUBLE *RESTRICT dcca
             , short *RESTRICT lFaceVelR, short *RESTRICT lFaceVelL
             , DOUBLE const viscosity   , DOUBLE const density
             , DOUBLE *RESTRICT wallPar , DOUBLE *dWall     
             , short const wallType     , short const nFace )
{ 
  bool fWall = false;
  short i, j, nf, nCarg, idCell,type;
  INT vizNel;
  DOUBLE wt, yPlusMax, uPlusMax,uFricMax,sWallMax,par[MAXLOADPARAMETER];
  DOUBLE dMin, vParallel[3], lNormal[3],yPlus,uPlus,uFric,sW,velB[3];

  sWallMax = uFricMax = uPlusMax = yPlusMax = 0.e0;
  dMin     = 1.e+16;
  for (nf = 0; nf<nFace; nf++) {
    vizNel         = lViz[nf]; 
    yPlus = 0.e0; 
/*... contorno*/
    if (vizNel  == -2) {
      if (lFaceVelR[nf] > 0) {
        nCarg = lFaceVelL[nf] - 1;
        if( lVel[nCarg].type == MOVEWALL){
          fWall = true;
/*... velocidade da parede*/
          nCarg     = lFaceVelL[nf] - 1;
          getLoads(par,loadsVel[nCarg]);
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
 * LESDYNAMIC: Metodo dinamico para o calculo da constente do LES    *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * lViz      -> vizinhos da celula central                           *
 * volume    -> volume da celula central                             *
 * lDensity  -> massa especifica sem variacao temporal               *
 * vel       -> campo de velocidades                                 *
 * gradVel   -> gradiente rescontruido das velocidades               *
 * nFace     -> numero de faces da celula central                    * 
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * retorna cs*delta                                                  *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * delta a largura do filtro que é igual a raiz cubica do volume da  *
 * celula central                                                    *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE sigmaModel(DOUBLE *RESTRICT s, DOUBLE *restrict gradVel 
                , short const nFace , short const ndm) {

  short i,idCell;
  INT vizNel,ip[3];
  DOUBLE g[9],x[9],op;
 
/*... tgg*/
//g[0] = MAT2D(0,0,gradVel,3);
//g[1] = MAT2D(0,1,gradVel,3);
//g[2] = MAT2D(0,2,gradVel,3);
//g[3] = MAT2D(1,0,gradVel,3);
//g[4] = MAT2D(1,1,gradVel,3);
//g[5] = MAT2D(1,2,gradVel,3);
//g[6] = MAT2D(2,0,gradVel,3);
//g[7] = MAT2D(2,1,gradVel,3);  
//g[8] = MAT2D(2,2,gradVel,3);
/*...*/
//gt[0] = MAT2D(0,0,gradVel,3);
//gt[1] = MAT2D(1,0,gradVel,3);
//gt[2] = MAT2D(2,0,gradVel,3);
//gt[3] = MAT2D(0,1,gradVel,3);
//gt[4] = MAT2D(1,1,gradVel,3);
//gt[5] = MAT2D(2,1,gradVel,3);
//gt[6] = MAT2D(0,2,gradVel,3);
//gt[7] = MAT2D(1,2,gradVel,3);  
//gt[8] = MAT2D(2,2,gradVel,3);

//dgemm(3,3,3,g,gt,gtg);
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
               ,s    ,ip
               ,3    ,1.e-14
               ,10000,false);  
/*...................................................................*/

/*... valores singulares de g*/
  s[0] = sqrt(fabs(s[0]));
  s[1] = sqrt(fabs(s[1]));
  s[2] = sqrt(fabs(s[2]));
/*...................................................................*/

/*...*/  
   op = s[2]*(s[0]-s[1])*(s[1]-s[2])/(s[0]*s[0]);  
// printf("%e %e %e %e\n",s[0],s[1],s[2],op);
/*...................................................................*/
   return op;
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

  a3=D1DIV3*acos(tmp);;
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