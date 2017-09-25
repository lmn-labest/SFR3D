#include<Turbulence.h>
/********************************************************************* 
 * Data de criacao    : 11/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
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
void turbulence(Loads *loadsVel             , Turbulence tModel             
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
      , DOUBLE *RESTRICT density
      , DOUBLE *RESTRICT dViscosity         , DOUBLE *RESTRICT eddyViscosity                          
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
         lVel[(MAX_NUM_FACE+1)*MAX_NDM],
         lEddyViscosity,lDviscosity; 
           

/*... loop nas celulas*/
  aux2    = maxViz+1;
  for(nel=0;nel<numel;nel++){
/*...*/
    aux1    = nFace[nel];

/*... loop na celula central*/    
    lMat              = mat[nel]-1;
    lib               = calType[lMat];
    lVolume[aux1]     = gVolume[nel]; 
    lGeomType[aux1]   = geomType[nel];
    lFaceVelR[aux1]   = MAT2D(nel,aux1,faceVelR ,aux2);
    lFaceVelL[aux1]   = MAT2D(nel,aux1,faceVelL ,aux2);
    lDensity[aux1]    = MAT2D(nel,2   ,density ,DENSITY_LEVEL);
/*...*/
    lDviscosity       = dViscosity[nel];
/*...*/      
    for(j=0;j<MAXPROP;j++)
      MAT2D(aux1,j,lProp,MAXPROP) = MAT2D(lMat,j,prop,MAXPROP);
/*...................................................................*/

/*...*/
    for(j=0;j<ndm;j++){
       MAT2D(aux1,j,lVel,ndm)  = MAT2D(nel,j,vel     ,ndm);
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

        lMat         = mat[vizNel]-1;
 
        for(j=0;j<ndm;j++){
          MAT2D(i,j,lCc      ,ndm) = MAT2D(vizNel,j,cc      ,ndm);
          MAT2D(i,j,lVel     ,ndm) = MAT2D(vizNel,j,vel     ,ndm);
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
     cellLibTurbulence(loadsVel        , tModel        
                      , lGeomType      , lProp
                      , lViz           , lKsi
                      , lmKsi 
                      , lEta           , lfArea 
                      , lNormal        , lVolume
                      , lXm            , lXmcc
                      , lDcca          , lCc
                      , lvSkew         , lmvSkew
                      , lFaceVelR      , lFaceVelL 
                      , lVel           , lGradVel
                      , lDensity       , lDviscosity 
                      , &lEddyViscosity
                      , nen[nel]       , nFace[nel] 
                      , ndm            , lib   
                      , nel);    
/*...................................................................*/
 
/*...*/
    eddyViscosity[nel]  = lEddyViscosity;
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
void cellLes(Loads *loadsVel           , Turbulence tModel           
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
  short nAresta, vizNel, nCarg, idCell = nFace,wallType;
  DOUBLE modS, tmp, densityC, viscosityC, cs, s[3], gradVelC[2][2],delta;
  DOUBLE wt,velC[2],vParallel[2],lNormal[2],lMin,dMin
        ,yPlus,uPlus,velB[2],yPlusMax;

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
        if( loadsVel[nCarg].type == MOVEWALL){
          fWall = true;
/*... velocidade da parede*/
          nCarg     = lFaceVelL[nAresta] - 1;
          velB[0]   = loadsVel[nCarg].par[0];
          velB[1]   = loadsVel[nCarg].par[1];
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

  switch (iCod) {

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