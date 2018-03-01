#include<CellLib.h>
/*********************************************************************
 * Biblioteca de celulas                                             *
 *-------------------------------------------------------------------*
 * Celulas 2D                                                        *
 *-------------------------------------------------------------------*
 * cellTrans2D - Celula de transporte geral                          *
 * cellEnergy2D- Celula de transporte da equacao de energia          *
 *-------------------------------------------------------------------*
 * Celulas 3D                                                        *
 *-------------------------------------------------------------------*
 * cellTrans3D - Celula de transporte geral                          *
 * cellEnergy3D- Celula de transporte da equacao de energia          *
 ********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 07/02/2018                                   *
 *-------------------------------------------------------------------*
 * CELLTRANS2D: Celula 2D para transporte                            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads     -> definicoes de cargas                                 * 
 * advT      -> tecnica da discretizacao do termo advecao            * 
 * diffT     -> tecnica da discretizacao do termo difusivo           *
 * lnFace    -> numero de faces da celula central e seus vizinhos    * 
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lprop     -> propriedade fisicas das celulas                      * 
 * lViz      -> viznhos da celula central                            * 
 * lId       -> numeracoes das equacoes das celulas                  * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * mEta      -> modulo do vetor eta                                  * 
 * normal    -> vetores normais as faces das celulas                 * 
 * area      -> area da celula central                               * 
 * xm        -> pontos medios das faces da celula central            * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * vSkew  -> vetor entre o ponto medio a intersecao que une os       * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew -> distacia entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca      -> menor distancia do centroide central a faces desta   *
 *              celula                                               * 
 * lDensity  -> massa especifica com variacao temporal               * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceR     -> restricoes por elmento                               * 
 * faceL     -> carga por elemento                                   * 
 * u0        -> solucao conhecida                                    * 
 * gradU0    -> gradiente rescontruido da solucao conhecida          * 
 * vel       -> campo de velocidade conhecido                        * 
 * cc        -> centroides da celula centra e seus vizinhos          *
 * nEn       -> numero de nos da celula central                      * 
 * nFace     -> numero de faces da celula central                    * 
 * ndm       -> numero de dimensoes                                  * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 * lRcell    -> residuo por celula                                   *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellTrans2D(Loads *loads           
              ,Advection advT           ,Diffusion diffT
              ,short *RESTRICT lGeomType,DOUBLE *RESTRICT prop
              ,INT *RESTRICT lViz       ,INT *RESTRICT lId  
              ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
              ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT mEta
              ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
              ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
              ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
              ,DOUBLE *RESTRICT lRcell  ,Temporal const ddt
              ,short  *RESTRICT lFaceR  ,short *RESTRICT lFaceL
              ,DOUBLE *RESTRICT u0      ,DOUBLE *RESTRICT gradU0
              ,DOUBLE *RESTRICT vel     ,DOUBLE *RESTRICT cc
              ,const short nEn          ,short const nFace    
              ,const short ndm          ,INT const nel)
{ 

  DOUBLE coefDifC,coefDif,coefDifV,rCell,dt,dt0;
  DOUBLE densityC,densityF,densityM;
  DOUBLE p,sP,dfd,gfKsi,lvSkew[2];
  DOUBLE v[2],gradUcomp[2],lKsi[2],lNormal[2],gf[2];
  DOUBLE dPviz,lModKsi,lModEta,du,duDksi,lXmcc[2],lXm[2];
  DOUBLE gradUp[2],gradUv[2],ccV[2];
  DOUBLE alpha,alphaMenosUm;
  DOUBLE tA,coef,tmp;
/*... nonOrtogonal*/
  DOUBLE e[2], t[2], modE, dfdc;
  DOUBLE xx[3];
/*... */
  DOUBLE wfn,wf[2],velC[2],velF[2],cv,cvc;
  short iCodAdv1 = advT.iCod1;
  short iCodAdv2 = advT.iCod2;
  short iCodDif  = diffT.iCod;
/*...*/
  short idCell = nFace;
  short nAresta,nCarg,typeTime;
  INT vizNel;
  bool fTime;
/*...*/
  DOUBLE pAdv[NPADV];

/*...*/
  dt       = ddt.dt[0];
  dt0      = ddt.dt[1];
  typeTime = ddt.type;
  fTime    = ddt.flag;
  densityC = lDensity[idCell];
/*...................................................................*/

/*... propriedades da celula*/
  coefDifC = MAT2D(idCell,COEFDIF,prop,MAXPROP);
/*...................................................................*/

/*...*/
  gradUp[0] = MAT2D(idCell,0,gradU0,ndm);
  gradUp[1] = MAT2D(idCell,1,gradU0,ndm);
  velC[0]   = MAT2D(idCell,0,vel,ndm);
  velC[1]   = MAT2D(idCell,1,vel,ndm);
/*...................................................................*/

  p          = 0.0e0;
  sP         = 0.0e0;
  for(nAresta=0;nAresta<nFace;nAresta++){
    vizNel     = lViz[nAresta];
    lNormal[0] = MAT2D(nAresta,0,normal,ndm);
    lNormal[1] = MAT2D(nAresta,1,normal,ndm);
    lModEta    = mEta[nAresta];
    lXmcc[0]   = MAT2D(nAresta,0,xmcc,ndm);
    lXmcc[1]   = MAT2D(nAresta,1,xmcc,ndm);
/*... dominio*/
    if( vizNel  > -1 ){
/*...*/
      velF[0] = MAT2D(nAresta, 0, vel, ndm);
      velF[1] = MAT2D(nAresta, 1, vel, ndm);
      densityF = lDensity[nAresta];
/*...*/
      lKsi[0]    = MAT2D(nAresta,0,ksi,ndm);
      lKsi[1]    = MAT2D(nAresta,1,ksi,ndm);
      lModKsi    = mKsi[nAresta];
 /*...*/
      lvSkew[0]  = MAT2D(nAresta,0,vSkew,ndm);
      lvSkew[1]  = MAT2D(nAresta,1,vSkew,ndm);
 /*...*/
      duDksi     = (u0[nAresta] - u0[idCell]) / lModKsi;
 /*...*/
      gradUv[0]  = MAT2D(nAresta,0,gradU0,ndm);
      gradUv[1]  = MAT2D(nAresta,1,gradU0,ndm);

      ccV[0] = MAT2D(nAresta, 0, cc, ndm);
      ccV[1] = MAT2D(nAresta, 1, cc, ndm);
    
      lXm[0] = MAT2D(nAresta,0,xm,ndm);
      lXm[1] = MAT2D(nAresta,1,xm,ndm);
/*...................................................................*/

/*... termo difusivo
grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      difusionScheme(lNormal,lKsi
                    ,lModEta,lModKsi
                    ,e      ,t
                    ,ndm    ,iCodDif);
/*...................................................................*/

/*...*/
      v[0]         = lvSkew[0] + lXmcc[0];
      v[1]         = lvSkew[1] + lXmcc[1];
      dPviz        = sqrt(v[0]*v[0] + v[1]*v[1]);
      alpha        = dPviz/lModKsi;
      alphaMenosUm = 1.0e0 - alpha; 
/*...................................................................*/

/*... media harmonica*/
      coefDifV = MAT2D(nAresta,COEFDIF,prop,MAXPROP); 
      coefDif  = alpha/coefDifC + alphaMenosUm/coefDifV;
      coefDif  = 1.0e0/coefDif;
/*...................................................................*/

/*... difusao direta*/
      coef = coefDif;
      modE = sqrt(e[0]*e[0] + e[1]*e[1]);
      dfd = coef*modE / lModKsi;
/*...................................................................*/
      
/*...*/
      gf[0]    = alphaMenosUm*gradUp[0] + alpha*gradUv[0];
      gf[1]    = alphaMenosUm*gradUp[1] + alpha*gradUv[1];
      wf[0]    =   alphaMenosUm*velC[0] + alpha*velF[0];
      wf[1]    =   alphaMenosUm*velC[1] + alpha*velF[1];
      densityM = alphaMenosUm*densityC  + alpha*densityF;
/*...................................................................*/

/*... velocidade normal a face*/
      wfn = wf[0]*lNormal[0] + wf[1]*lNormal[1];
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi = gf[0] *lKsi[0]    + gf[1] *lKsi[1];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
      du           = duDksi  - gfKsi;
      gradUcomp[0] = gf[0] + du*lKsi[0];
      gradUcomp[1] = gf[1] + du*lKsi[1];
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi = gradUcomp[0]*t[0] + gradUcomp[1]*t[1];
/*...................................................................*/

/*... correcao nao-ortogonal*/
      dfdc = coefDif*gfKsi;    
/*...................................................................*/

/*... fluxo convectivo upwind de primeira ordem*/
      cv   = densityM*wfn*lModEta;
/*...................................................................*/

/*...*/
      v[0] = lXm[0] - ccV[0];
      v[1] = lXm[1] - ccV[1];
      advectiveSchemeScalar(u0[idCell]   ,u0[nAresta]
                           ,gradUp       ,gradUv
                           ,gradUcomp    ,lvSkew
                           ,lXmcc        ,v    
                           ,lKsi         ,lModKsi
                           ,cv           ,&cvc
                           ,alphaMenosUm ,alpha    
                           ,pAdv         ,ndm
                           ,iCodAdv1     ,iCodAdv2);
/*...................................................................*/

/*...*/
      lA[nAresta] = dfd - min(cv,0e0);
      sP         += cv;    
/*... correcao nao ortogonal e do fluxo advectivo*/        
      p          += dfdc - cv*cvc;
/*...................................................................*/
    }
/*... contorno*/
    else{
      lA[nAresta] = 0.e0;
      if(lFaceR[nAresta]){
        wfn = velC[0]*lNormal[0] + velC[1]*lNormal[1];
/*...cargas*/
        nCarg=lFaceL[nAresta]-1;
        xx[0] = MAT2D(nAresta,0,xm,2);
        xx[1] = MAT2D(nAresta,1,xm,2);
        xx[2] = 0.e0;                    
        pLoad(&sP           ,&p
             ,&tA
             ,coefDifC      ,densityC
             ,wfn           ,xx 
             ,lModEta       ,dcca[nAresta]
             ,loads[nCarg]  ,true);
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... distretizacao temporal*/
  if (fTime) {
/*... EULER*/
    if (typeTime == EULER)
      sP += densityC*volume[idCell] / dt;
/*...BACKWARD*/
    else if (typeTime == BACKWARD) {
      tmp = 1.e0 / dt + 1.e0 / (dt + dt0);
      sP += tmp*densityC*volume[idCell];
    }
  }
/*...................................................................*/

/*...*/
  if(nFace == 3){
    lA[idCell] = sP + lA[0] + lA[1] + lA[2];
  }
  else if(nFace == 4){
    lA[idCell] = sP + lA[0] + lA[1] + lA[2] + lA[3];
  }
/*...................................................................*/

/*...*/
  rCell = 0.0e0;
  for(nAresta=0;nAresta<nFace;nAresta++){
    if( lViz[nAresta] > -1){
/*... pasando os valoeres conhecidos para o lado direito*/
      if(lId[nAresta] == -2)
        p += lA[nAresta]*u0[nAresta]; 
      else
/*residuo (R = F-KvizUviz ) e valores prescritos por elemento*/
        rCell += lA[nAresta]*u0[nAresta]; 
    }
  }
/*... residuo: R = F - KpUp*/ 
  rCell += p -lA[idCell]*u0[idCell];   
/*...................................................................*/

/*...*/
  if(nFace == 3){
    lA[0] *= -1.e0;
    lA[1] *= -1.e0;
    lA[2] *= -1.e0;
  }
  else if(nFace == 4){
    lA[0] *= -1.e0;
    lA[1] *= -1.e0;
    lA[2] *= -1.e0;
    lA[3] *= -1.e0;
  }
/*...................................................................*/

/*...*/
  lB[0]     = p;
  lRcell[0] = rCell;
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 20/08/2017                                   *
 * Data de modificaco : 24/09/2017                                   *
 *-------------------------------------------------------------------*
 * CELLENERGY2D: Celula 2D para transporte                           *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loads     -> definicoes de cargas                                 *
 * model     -> modelo da equacao de energia                         *
 * advT      -> tecnica da discretizacao do termo advecao            *
 * diffT     -> tecnica da discretizacao do termo difusivo           *
 * tModel    -> modelo de turbulencia                                *
 * vProp     -> propedades variaveis (true|false)                    *
 * lnFace    -> numero de faces da celula central e seus vizinhos    *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lId       -> numeracoes das equacoes das celulas                  *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * mEta      -> modulo do vetor eta                                  *
 * normal    -> vetores normais as faces das celulas                 *
 * area      -> area da celula central                               *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * vSkew  -> vetor entre o ponto medio a intersecao que une os       *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew -> distacia entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * dcca      -> menor distancia do centroide central a faces desta   *
 *              celula                                               *
 * lDensity  -> massa especifica com variacao temporal               *
 * lSheat    -> calor especifico com variacao temporal               *
 * lDviscosity-> viscosidade dinamica com variacao temporal          *
 * lTconductivity-> condutividade termica com variacao temporal      *
 * dField    -> matriz D do metodo simple                            * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceR     -> restricoes por elmento                               *
 * faceL     -> carga por elemento                                   *
 * u0        -> solucao conhecida                                    *
 * gradU0    -> gradiente rescontruido da solucao conhecida          *
 * vel       -> campo de velocidade conhecido                        *
 * gradVel   -> gradiente rescontruido da velocidade                 *
 * pres      -> pressao do tempo atual e do tempo anterior           *
 * gradPres  -> gradiente de pressao do tempo atual                  *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * underU    -> parametro de sob relaxamento                         *
 * fSheat    -> calor especifico com variacao com a Temperatura      *
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 * lRcell    -> residuo por celula                                   *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellEnergy2D(Loads *loads            , Loads *loadsVel
            , Advection advT              , Diffusion diffT
            , Turbulence tModel           , EnergyModel model
            , PropVar vProp
            , short *RESTRICT lGeomType   , DOUBLE *RESTRICT prop
            , INT *RESTRICT lViz          , INT *RESTRICT lId
            , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
            , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT mEta
            , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume
            , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
            , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
            , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
            , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
            , DOUBLE *RESTRICT lRcell     , Temporal const ddt
            , short  *RESTRICT lFaceR     , short *RESTRICT lFaceL
            , short  *RESTRICT lFaceVelR  , short *RESTRICT lFaceVelL
            , DOUBLE *RESTRICT u0         , DOUBLE *RESTRICT gradU0
            , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
            , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres  
            , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lSheat
            , DOUBLE *RESTRICT lViscosity , DOUBLE *RESTRICT lTconductivity
            , DOUBLE *RESTRICT dField
            , DOUBLE const underU           
            , const short nEn             , short const nFace
            , const short ndm             , INT const nel)
{
  bool fTime, fDisp, fRes, fPresWork, fTemp, fTurb, fWallModel, fKelvin
     , fSheat;
  short iCodAdv1, iCodAdv2, iCodDif, wallType, idCell, nAresta, nCarg1
        , nCarg2, typeTime;
/*...*/
  INT vizNel;
/*...*/
  DOUBLE thermCoefC, thermCoefV, densityC, densityV, densityM,
         diffEffC, diffEffV, diffEff, sHeatC, sHeatV, sHeatM, 
         eddyViscosityC, eddyViscosityV, viscosityC, tA, coef,
         tmp, tmp1, prTwall, prTsgs;
  DOUBLE p, sP, dfd, gfKsi, lvSkew[2], alpha, alphaMenosUm;
  DOUBLE v[2], gradUcomp[2], lKsi[2], lNormal[2], gf[2];
  DOUBLE dPviz, lModKsi, lModEta, du, duDksi, lXmcc[2], lXm[2];
  DOUBLE gradUp[2], gradUv[2], ccV[2], gradVelC[2][2], rCell, dt, dt0; 
  DOUBLE uC, uV;
/*... nonOrtogonal*/
  DOUBLE e[2], t[2], modE, dfdc, xx[3];
/*... */
  DOUBLE presC, presC0, presV, gradPresC[2], gradPresV[2], wfn
        , velC[2], velV[2], dFieldC[2], dFieldV[2], dFieldF[2], cv, cvc;       

/*...*/
  DOUBLE phi,psi,lambda;
/*...*/
  DOUBLE pAdv[NPADV];

/*...*/
  idCell   = nFace;
  iCodAdv1 = advT.iCod1;
  iCodAdv2 = advT.iCod2;
  pAdv[0]  = advT.par[0];
  iCodDif  = diffT.iCod;
/*...................................................................*/

/*...*/  
  dt         = ddt.dt[0];
  dt0        = ddt.dt[1];
  typeTime   = ddt.type;
  fTime      = ddt.flag;  
  fDisp      = model.fDissipation;
  fRes       = model.fRes;
  fPresWork  = model.fPresWork;
  fTemp      = model.fTemperature;
  fKelvin    = model.fKelvin;
  fTurb      = tModel.fTurb;
  prTwall    = tModel.PrandltTwall;
  prTsgs     = tModel.PrandltTsgs;
  fWallModel = tModel.fWall;
  wallType   = tModel.wallType;
  fSheat     = vProp.fSpecificHeat;
/*...................................................................*/

/*... propriedades da celula*/
  sHeatM      =  eddyViscosityC = eddyViscosityV = 0.e0;
  densityC    = lDensity[idCell];
  sHeatC      = lSheat[idCell];
  thermCoefC  = lTconductivity[idCell];  
  viscosityC  = MAT2D(idCell, 0, lViscosity, 2);
  if(fTurb) eddyViscosityC = MAT2D(idCell, 1, lViscosity, 2);
  lambda      = -D2DIV3*(viscosityC + eddyViscosityC);  
/*...................................................................*/

/*...*/
  if(fTemp)
    diffEffC = thermCoefC + sHeatC*viscosityC/prTsgs;
  else
    diffEffC = thermCoefC/sHeatC + viscosityC/prTsgs;
/*...................................................................*/

/*...*/
  gradUp[0]    = MAT2D(idCell, 0, gradU0, ndm);
  gradUp[1]    = MAT2D(idCell, 1, gradU0, ndm);
  velC[0]      = MAT2D(idCell, 0, vel, ndm);
  velC[1]      = MAT2D(idCell, 1, vel, ndm);
  uC           = u0[idCell];
/*... p(n)*/
  presC0 = MAT2D(idCell, 0, pres, 2); 
/*... p(n+1)*/
  presC  = MAT2D(idCell, 1, pres, 2);
/*...*/
  gradPresC[0] = MAT2D(idCell, 0, gradPres, ndm);
  gradPresC[1] = MAT2D(idCell, 1, gradPres, ndm);
  dFieldC[0]   = MAT2D(idCell, 0, dField, 2); 
  dFieldC[1]   = MAT2D(idCell, 1, dField, 2);
/*...................................................................*/

  p = 0.0e0;
  sP = 0.0e0;
  for (nAresta = 0; nAresta<nFace; nAresta++) {
    vizNel = lViz[nAresta];
    lNormal[0] = MAT2D(nAresta, 0, normal, ndm);
    lNormal[1] = MAT2D(nAresta, 1, normal, ndm);
    lModEta = mEta[nAresta];
    lXmcc[0] = MAT2D(nAresta, 0, xmcc, ndm);
    lXmcc[1] = MAT2D(nAresta, 1, xmcc, ndm);
/*... dominio*/
    if (vizNel  > -1) {
/*...*/
      densityV = lDensity[nAresta]; 
      sHeatV   = lSheat[nAresta];    
      thermCoefV = lTconductivity[nAresta];  
      if(fTurb) eddyViscosityV = MAT2D(nAresta, 1, lViscosity, 2);
/*... p(n+1)*/
      presV     = MAT2D(nAresta, 1, pres, 2);
      gradPresV[0] = MAT2D(nAresta, 0, gradPres, ndm);
      gradPresV[1] = MAT2D(nAresta, 1, gradPres, ndm);
      dFieldV[0]   = MAT2D(nAresta, 0, dField, 2);
      dFieldV[1]   = MAT2D(nAresta, 1, dField, 2);
/*...*/
      uV      = u0[nAresta];
      velV[0] = MAT2D(nAresta, 0, vel, ndm);
      velV[1] = MAT2D(nAresta, 1, vel, ndm);      
/*...*/
      lKsi[0] = MAT2D(nAresta, 0, ksi, ndm);
      lKsi[1] = MAT2D(nAresta, 1, ksi, ndm);
      lModKsi = mKsi[nAresta];
/*...*/
      lvSkew[0] = MAT2D(nAresta, 0, vSkew, ndm);
      lvSkew[1] = MAT2D(nAresta, 1, vSkew, ndm);
/*...*/
      duDksi = (uV - uC) / lModKsi;
/*...*/
      gradUv[0] = MAT2D(nAresta, 0, gradU0, ndm);
      gradUv[1] = MAT2D(nAresta, 1, gradU0, ndm);

      ccV[0] = MAT2D(nAresta, 0, cc, ndm);
      ccV[1] = MAT2D(nAresta, 1, cc, ndm);

      lXm[0] = MAT2D(nAresta, 0, xm, ndm);
      lXm[1] = MAT2D(nAresta, 1, xm, ndm);
/*...................................................................*/

/*... termo difusivo
      grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      difusionScheme(lNormal ,lKsi
                    ,lModEta ,lModKsi
                    ,e       ,t
                    ,ndm     ,iCodDif);
/*...................................................................*/

/*...*/
      v[0] = lvSkew[0] + lXmcc[0];
      v[1] = lvSkew[1] + lXmcc[1];
      dPviz = sqrt(v[0] * v[0] + v[1] * v[1]);
      alpha = dPviz / lModKsi;
      alphaMenosUm = 1.0e0 - alpha;
/*...................................................................*/

/*... media harmonica*/ 
      if(fTemp){
        diffEffV = thermCoefV + sHeatV*eddyViscosityV/prTsgs;
        diffEff = alpha / diffEffC + alphaMenosUm / diffEffV;
        diffEff = 1.0e0 /  diffEff;
      }
      else {        
        diffEffV = thermCoefV/sHeatV + eddyViscosityV/prTsgs;
        diffEff = alpha / diffEffC + alphaMenosUm / diffEffV;
        diffEff = 1.0e0 / diffEff;
      }
/*...................................................................*/

/*... difusao direta*/
      coef = diffEff;
      modE = sqrt(e[0] * e[0] + e[1] * e[1]);
      dfd = coef*modE / lModKsi;
/*...................................................................*/

/*...*/
      gf[0] = alphaMenosUm*gradUp[0] + alpha*gradUv[0];
      gf[1] = alphaMenosUm*gradUp[1] + alpha*gradUv[1];
      dFieldF[0] = alphaMenosUm*dFieldC[0] + alpha*dFieldV[0];
      dFieldF[1] = alphaMenosUm*dFieldC[1] + alpha*dFieldV[1];
      densityM = alphaMenosUm*densityC + alpha*densityV;
      if(fTemp)
        sHeatM   = alphaMenosUm*sHeatC + alpha*sHeatV;
/*...................................................................*/

/*... velocidade normal a face*/
      wfn = interpolFaceVel(velC         ,velV
                           ,presC        ,presV
                           ,gradPresC    ,gradPresV
                           ,lNormal      ,lKsi
                           ,lModKsi      ,dFieldF
                           ,alphaMenosUm ,alpha
                           ,ndm);    
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi = gf[0] * lKsi[0] + gf[1] * lKsi[1];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
      du = duDksi - gfKsi;
      gradUcomp[0] = gf[0] + du*lKsi[0];
      gradUcomp[1] = gf[1] + du*lKsi[1];
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi = gradUcomp[0] * t[0] + gradUcomp[1] * t[1];
/*...................................................................*/

/*... correcao nao-ortogonal*/
      dfdc = coef*gfKsi;
/*...................................................................*/

/*... fluxo convectivo upwind de primeira ordem*/
      if(fTemp)
        cv = sHeatM*densityM*wfn*lModEta;
      else
        cv = densityM*wfn*lModEta;
/*...................................................................*/
 
/*...*/
      v[0] = lXm[0] - ccV[0];
      v[1] = lXm[1] - ccV[1];
      advectiveSchemeScalar(uC        ,uV
                           ,gradUp    ,gradUv
                           ,gradUcomp ,lvSkew
                           ,lXmcc     ,v
                           ,lKsi      ,lModKsi
                           ,cv        ,&cvc
                           ,alphaMenosUm ,alpha
                           ,pAdv      ,ndm   
                           ,iCodAdv1  ,iCodAdv2);
/*...................................................................*/

/*...*/
      lA[nAresta] = dfd - min(cv, 0e0);
      sP += cv;
/*... correcao nao ortogonal e do fluxo advectivo*/
      p += dfdc - cv*cvc;
/*...................................................................*/
    }
/*... contorno*/
    else {
      lA[nAresta] = 0.e0;
      if (lFaceR[nAresta]) {
        wfn = velC[0] * lNormal[0] + velC[1] * lNormal[1];
/*...cargas*/
        nCarg1 = lFaceL[nAresta] - 1;
        nCarg2 = lFaceVelL[nAresta] - 1;
        xx[0] = MAT2D(nAresta, 0, xm, 2);
        xx[1] = MAT2D(nAresta, 1, xm, 2);
        xx[2] = 0.e0;
/*      pLoadEnergy(&sP          , &p
                  , &tA          , velC
                  , uC           , lNormal  
                  , thermCoefC   , densityC
                  , viscosityC   , sHeatC
                  , prTwall      , xx                   
                  , lModEta      , dcca[nAresta]
                  , loads[nCarg1], loadsVel[nCarg2]
                  , ndm
                  , true         , fTemp
                  , fKelvin      , fSheat
                  , fWallModel   , wallType);  */
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... | du1/dx1 du1/dx2*/
  gradVelC[0][0] = MAT3D(idCell, 0, 0, gradVel, 2, ndm);
  gradVelC[0][1] = MAT3D(idCell, 0, 1, gradVel, 2, ndm);
/*... | du2/dx1 du2/dx2*/
  gradVelC[1][0] = MAT3D(idCell, 1, 0, gradVel, 2, ndm);
  gradVelC[1][1] = MAT3D(idCell, 1, 1, gradVel, 2, ndm);
/*...................................................................*/

/*...*/
  if(fDisp){
    tmp = gradVelC[0][0] + gradVelC[1][1];
    psi = tmp*tmp;
/*...*/
    tmp  = gradVelC[0][0]*gradVelC[0][0] + gradVelC[1][1]*gradVelC[1][1];
    tmp1 = gradVelC[0][1] + gradVelC[1][0];
    phi = 2.0*tmp + tmp1*tmp1;    
    tmp = ( lambda*psi + (viscosityC + eddyViscosityC)*phi );
/*...*/
    p+= tmp*volume[idCell];
  }
/*.....................................................................*/ 

/*...*/
  if(fPresWork){
/*... derivada materia da pressao*/
    tmp = (pres[1] - pres[0])/dt 
        + velC[0]*gradPresC[0] + velC[1]*gradPresC[1];
/*...*/
    p+= tmp*volume[idCell];
/*.....................................................................*/
  }
/*.....................................................................*/

/*... distretizacao temporal*/
  if(fTemp)
    tmp1 = sHeatC*densityC;
  else
    tmp1 = densityC;

  if (fTime) {
/*... EULER*/
    if (typeTime == EULER)
      sP += tmp1*volume[idCell] / dt;
/*...BACKWARD*/
    else if (typeTime == BACKWARD) {
      tmp = 1.e0 / dt + 1.e0 / (dt + dt0);
      sP += tmp*tmp1*volume[idCell];
    }
  }
/*...................................................................*/

/*...*/
  if (nFace == 3) {
    lA[idCell] = sP + lA[0] + lA[1] + lA[2];
  }
  else if (nFace == 4) {
    lA[idCell] = sP + lA[0] + lA[1] + lA[2] + lA[3];
  }
/*...................................................................*/

/*...*/
  rCell = 0.0e0;
  for (nAresta = 0; nAresta<nFace; nAresta++) {
    if (lViz[nAresta] > -1) {
/*... pasando os valoeres conhecidos para o lado direito*/
      if (lId[nAresta] == -2) 
        p += lA[nAresta] * u0[nAresta];
      else
/*residuo (R = F-KvizUviz ) e valores prescritos por elemento*/
        rCell += lA[nAresta] * u0[nAresta];
    }
  }
/*... residuo: R = F - KpUp*/
  rCell += p - lA[idCell] * u0[idCell];
/*...................................................................*/

/*... under-relaxation(simple)*/
  lA[idCell] = lA[idCell] / underU;
  if(!fRes)
    p += (1.e0 - underU)*lA[idCell] * u0[idCell];
/*...................................................................*/

/*...*/
  if (nFace == 3) {
    lA[0] *= -1.e0;
    lA[1] *= -1.e0;
    lA[2] *= -1.e0;
  }
  else if (nFace == 4) {
    lA[0] *= -1.e0;
    lA[1] *= -1.e0;
    lA[2] *= -1.e0;
    lA[3] *= -1.e0;
  }
/*...................................................................*/

/*...*/
  if(fRes)
    lB[0]     = rCell;
  else
    lB[0]     = p;
/*...*/
  lRcell[0] = rCell;
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 22/08/2016                                   *
 *-------------------------------------------------------------------*
 * CELLTRANS3D: Celula 3D para transporte                            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads     -> definicoes de cargas                                 * 
 * advT      -> tecnica da discretizacao do termo advecao            *
 * diffT   -> tecnica da discretizacao do termo difusivo             *
 * lnFace    -> numero de faces da celula central e seus vizinhos    * 
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lprop     -> propriedade fisicas das celulas                      * 
 * lViz      -> viznhos da celula central                            * 
 * lId       -> numeracoes das equacoes das celulas                  * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area das faces da celula central                     * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> volume da celula central                             * 
 * xm        -> pontos medios das faces da celula central            * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * vSkew  -> vetor entre o ponto medio a intersecao que une os       * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew -> distacia entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca      -> menor distancia do centroide central a faces desta   *
 *              celula                                               * 
 * lDensity  -> massa especifica com variacao temporal               * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceR     -> restricoes por elemento                              * 
 * faceL     -> carga por elemento                                   * 
 * u0        -> fsolucao conhecida                                   * 
 * gradU0    -> gradiente rescontruido da solucao conhecida          * 
 * vel       -> campo de velocidade conhecido                        * 
 * cc        -> centroides da celula centra e seus vizinhos          *
 * nEn       -> numero de nos da celula central                      * 
 * nFace     -> numero de faces da celula central                    * 
 * ndm       -> numero de dimensoes                                  * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 * lRcell    -> residuo por celula                                   *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellTrans3D(Loads *loads           
              ,Advection advT           ,Diffusion diffT
              ,short *RESTRICT lGeomType,DOUBLE *RESTRICT prop
              ,INT *RESTRICT lViz       ,INT *RESTRICT lId  
              ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
              ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea
              ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
              ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
              ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
              ,DOUBLE *RESTRICT lRcell  ,Temporal const ddt             
              ,short  *RESTRICT lFaceR  ,short  *RESTRICT lFaceL  
              ,DOUBLE *RESTRICT u0      ,DOUBLE *RESTRICT gradU0
              ,DOUBLE *RESTRICT vel     ,DOUBLE *RESTRICT cc
              ,const short nEn          ,short const nFace    
              ,const short ndm          ,INT const nel)
{ 

  DOUBLE coefDifC,coefDif,coefDifV,rCell,dt,dt0;
  DOUBLE densityC,densityF,densityM;
  DOUBLE p,sP,dfd,gfKsi,lvSkew[3];
  DOUBLE v[3],gradUcomp[3],lKsi[3],lNormal[3],gf[3];
  DOUBLE dPviz,lModKsi,lfArea,du,duDksi,lXmcc[3],lXm[3];
  DOUBLE gradUp[3],gradUv[3],ccV[3];
  DOUBLE alpha,alphaMenosUm;
  DOUBLE tA,coef,tmp;
/*... nonOrtogonal*/
  DOUBLE e[3], t[3], modE, dfdc;
  DOUBLE xx[3];
/*...*/
  DOUBLE wfn,wf[3],velC[3],velF[3],cv,cvc;
  short iCodAdv1=advT.iCod1;
  short iCodAdv2=advT.iCod2;
  short iCodDif = diffT.iCod;
/*...*/
  short idCell = nFace;
  short nf,nCarg,typeTime;
  INT vizNel;
  bool fTime;
/*...*/
  DOUBLE pAdv[NPADV];

/*...*/
  dt       = ddt.dt[0];
  dt0      = ddt.dt[1];
  typeTime = ddt.type;
  fTime    = ddt.flag;
  densityC = lDensity[idCell];
/*...................................................................*/
  
/*... propriedades da celula*/
  coefDifC = MAT2D(idCell,COEFDIF,prop,MAXPROP);
/*...................................................................*/

/*...*/
  gradUp[0] = MAT2D(idCell,0,gradU0,ndm);
  gradUp[1] = MAT2D(idCell,1,gradU0,ndm);
  gradUp[2] = MAT2D(idCell,2,gradU0,ndm);
  velC[0]   = MAT2D(idCell,0,vel,ndm);
  velC[1]   = MAT2D(idCell,1,vel,ndm);
  velC[2]   = MAT2D(idCell,2,vel,ndm);
/*...................................................................*/
      
  p          = 0.0e0;
  sP         = 0.0e0;
  for(nf=0;nf<nFace;nf++){
    vizNel     = lViz[nf];
    lNormal[0] = MAT2D(nf,0,normal,ndm);
    lNormal[1] = MAT2D(nf,1,normal,ndm);
    lNormal[2] = MAT2D(nf,2,normal,ndm);
    lfArea     = fArea[nf];
    lXmcc[0] = MAT2D(nf,0,xmcc,ndm);
    lXmcc[1] = MAT2D(nf,1,xmcc,ndm);
    lXmcc[2] = MAT2D(nf,2,xmcc,ndm);
/*... dominio*/
    if( vizNel  > -1 ){
/*...*/
      velF[0] = MAT2D(nf, 0, vel, ndm);
      velF[1] = MAT2D(nf, 1, vel, ndm);
      velF[2] = MAT2D(nf, 2, vel, ndm);
      densityF = lDensity[nf];
 /*...*/
      lKsi[0]    = MAT2D(nf,0,ksi,ndm);
      lKsi[1]    = MAT2D(nf,1,ksi,ndm);
      lKsi[2]    = MAT2D(nf,2,ksi,ndm);
      lModKsi    = mKsi[nf];
 /*...*/
      lvSkew[0]  = MAT2D(nf,0,vSkew,ndm);
      lvSkew[1]  = MAT2D(nf,1,vSkew,ndm);
      lvSkew[2]  = MAT2D(nf,2,vSkew,ndm);
 /*...*/
      duDksi     = (u0[nf] - u0[idCell]) / lModKsi;
 /*...*/   
      gradUv[0]  = MAT2D(nf,0,gradU0,ndm);
      gradUv[1]  = MAT2D(nf,1,gradU0,ndm);
      gradUv[2]  = MAT2D(nf,2,gradU0,ndm);

      ccV[0] = MAT2D(nf,0,cc,ndm);
      ccV[1] = MAT2D(nf,1,cc,ndm);
      ccV[2] = MAT2D(nf,2,cc,ndm);

      lXm[0] = MAT2D(nf,0,xm,ndm);
      lXm[1] = MAT2D(nf,1,xm,ndm);
      lXm[2] = MAT2D(nf,2,xm,ndm);
/*...................................................................*/

/*... termo difusivo
grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      difusionScheme(lNormal,lKsi
                    ,lfArea ,lModKsi
                    ,e      ,t
                    ,ndm    ,iCodDif);
/*...................................................................*/

/*...*/
      v[0]  = lvSkew[0] + lXmcc[0];
      v[1]  = lvSkew[1] + lXmcc[1];
      v[2]  = lvSkew[2] + lXmcc[2];
      dPviz = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
      alpha        = dPviz/lModKsi;
      alphaMenosUm = 1.0e0 - alpha; 
/*...................................................................*/

/*... media harmonica*/
      coefDifV = MAT2D(nf,COEFDIF,prop,MAXPROP); 
      coefDif  = alpha/coefDifC + alphaMenosUm/coefDifV;
      coefDif  = 1.0e0/coefDif;
/*...................................................................*/

/*... difusao direta*/
      coef = coefDif;
      modE = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
      dfd = coef*modE / lModKsi;
/*...................................................................*/
      
/*...*/
      gf[0]    = alphaMenosUm*gradUp[0] + alpha*gradUv[0];
      gf[1]    = alphaMenosUm*gradUp[1] + alpha*gradUv[1];
      gf[2]    = alphaMenosUm*gradUp[2] + alpha*gradUv[2];
      wf[0]    =   alphaMenosUm*velC[0] + alpha*velF[0];
      wf[1]    =   alphaMenosUm*velC[1] + alpha*velF[1];
      wf[2]    =   alphaMenosUm*velC[2] + alpha*velF[2];
      densityM = alphaMenosUm*densityC  + alpha*densityF;
/*...................................................................*/

/*... velocidade normal a face*/
      wfn = wf[0]*lNormal[0] + wf[1]*lNormal[1] + wf[2]*lNormal[2];
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi        = gf[0] *lKsi[0]  
                   + gf[1] *lKsi[1] 
                   + gf[2] *lKsi[2];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
      du           = duDksi  - gfKsi;
      gradUcomp[0] = gf[0] + du*lKsi[0];
      gradUcomp[1] = gf[1] + du*lKsi[1];
      gradUcomp[2] = gf[2] + du*lKsi[2];
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi = gradUcomp[0]*t[0] 
            + gradUcomp[1]*t[1] 
            + gradUcomp[2]*t[2];
/*...................................................................*/

/*... correcao nao-ortogonal*/
      dfdc = coefDif*gfKsi;    
/*...................................................................*/

/*... fluxo convectivo upwind de primeira ordem*/
      cv   = densityM*wfn*lfArea;
/*...................................................................*/

/*...*/
      v[0] = lXm[0] - ccV[0];
      v[1] = lXm[1] - ccV[1];
      v[2] = lXm[2] - ccV[2];
      advectiveSchemeScalar(u0[idCell]  ,u0[nf]
                           ,gradUp      ,gradUv
                           ,gradUcomp   ,lvSkew
                           ,lXmcc       ,v
                           ,lKsi        ,lModKsi
                           ,wfn         ,&cvc
                           ,alphaMenosUm,alpha     
                           ,pAdv        ,ndm
                           ,iCodAdv1    ,iCodAdv2);
/*...................................................................*/

/*...*/
      lA[nf] = dfd - min(cv,0.e0);
      sP    += cv;                    
/*... correcao nao ortogonal e do fluxo advectivo*/        
      p     += dfdc - cv*cvc;
/*...................................................................*/
    }
/*... contorno*/
    else{
      lA[nf] = 0.0e0;
      if(lFaceR[nf]){
        wfn = velC[0]*lNormal[0] 
            + velC[1]*lNormal[1] 
            + velC[2]*lNormal[2];
/*...cargas*/
        nCarg=lFaceL[nf]-1;
        xx[0] = MAT2D(nf,0,xm,3);
        xx[1] = MAT2D(nf,1,xm,3);
        xx[2] = MAT2D(nf,2,xm,3);        
        pLoad(&sP           ,&p
             ,&tA
             ,coefDifC      ,densityC
             ,wfn           ,xx 
             ,fArea[nf]     ,dcca[nf]
             ,loads[nCarg]  ,true);
      }
/*...................................................................*/
    }
/*...................................................................*/
  }

/*... distretizacao temporal*/
  if (fTime) {
    /*... EULER*/
    if (typeTime == EULER)
      sP += densityC*volume[idCell] / dt;
    /*...BACKWARD*/
    else if (typeTime == BACKWARD) {
      tmp = 1.e0 / dt + 1.e0 / (dt + dt0);
      sP += tmp*densityC*volume[idCell];
    }
  }
/*...................................................................*/

/*...*/ 
  if(nf == 4){
    lA[idCell] = sP + lA[0] + lA[1] + lA[2] + lA[3];
  }
  else if(nf == 6){
    lA[idCell] = sP + lA[0] + lA[1] + lA[2] + lA[3] + lA[4] + lA[5];
  }
/*...................................................................*/

/*...*/
  rCell = 0.0e0;
  for(nf=0;nf<nFace;nf++){
    if( lViz[nf] > -1){
/*... pasando os valoeres conhecidos para o lado direito*/
      if(lId[nf] == -2)
        p += lA[nf]*u0[nf]; 
      else
/*residuo (R = F-KvizUviz ) e valores prescritos por elemento*/
        rCell += lA[nf]*u0[nf]; 
    }
  }
/*... residuo: R = F - KpUp*/ 
  rCell += p -lA[idCell]*u0[idCell]; 
/*...................................................................*/

/*...*/  
/*
  for(nf=0;nf<nFace;nf++){
   lA[nf] *= -1.e0;
  }
*/
  if(nf == 4){
    lA[0] *= -1.e0;
    lA[1] *= -1.e0;
    lA[2] *= -1.e0;
    lA[3] *= -1.e0;
  }
  else if(nf == 6){
    lA[0] *= -1.e0;
    lA[1] *= -1.e0;
    lA[2] *= -1.e0;
    lA[3] *= -1.e0;
    lA[4] *= -1.e0;
    lA[5] *= -1.e0;
  }
/*...................................................................*/

/*...*/
  lB[0]     = p;
  lRcell[0] = rCell;
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 03/10/2017                                   *
 * Data de modificaco : 18/01/2018                                   *
 *-------------------------------------------------------------------*
 * CELLENERGY3D: Celula 3D para transporte                           *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loads     -> definicoes de cargas                                 *
 * model     -> modelo da equacao de energia                         *
 * advT      -> tecnica da discretizacao do termo advecao            *
 * diffT     -> tecnica da discretizacao do termo difusivo           *
 * tModel    -> modelo de turbulencia                                *
 * vProp     -> propedades variaveis (true|false)                    *
 * lnFace    -> numero de faces da celula central e seus vizinhos    *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lId       -> numeracoes das equacoes das celulas                  *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * fArea     -> area das faces                                       *
 * normal    -> vetores normais as faces das celulas                 *
 * area      -> area da celula central                               *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * vSkew  -> vetor entre o ponto medio a intersecao que une os       *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew -> distacia entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * dcca      -> menor distancia do centroide central a faces desta   *
 *              celula                                               *
 * lDensity  -> massa especifica com variacao temporal               *
 * lSheat    -> calor especifico com variacao temporal               *
 * lDviscosity-> viscosidade dinamica com variacao temporal          *
 * lTconductivity-> condutividade termica com variacao temporal      *
 * dField    -> matriz D do metodo simple                            * 
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceR     -> restricoes por elmento                               *
 * faceL     -> carga por elemento                                   *
 * u0        -> solucao conhecida                                    *
 * gradU0    -> gradiente rescontruido da solucao conhecida          *
 * vel       -> campo de velocidade conhecido                        *
 * gradVel   -> gradiente rescontruido da velocidade                 *
 * pres      -> pressao do tempo atual e do tempo anterior           *
 * gradPres  -> gradiente de pressao do tempo atual                  *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * underU    -> parametro de sob relaxamento                         *
 * fSheat    -> calor especifico com variacao com a Temperatura      *
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 * lRcell    -> residuo por celula                                   *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellEnergy3D(Loads *loads            , Loads *lVel
            , Advection advT              , Diffusion diffT
            , Turbulence tModel           , EnergyModel model
            , PropVar vProp
            , short *RESTRICT lGeomType   , DOUBLE *RESTRICT prop
            , INT *RESTRICT lViz          , INT *RESTRICT lId
            , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
            , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT fArea
            , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume
            , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
            , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
            , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
            , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
            , DOUBLE *RESTRICT lRcell     , Temporal const ddt
            , short  *RESTRICT lFaceR     , short *RESTRICT lFaceL
            , short  *RESTRICT lFaceVelR  , short *RESTRICT lFaceVelL
            , DOUBLE *RESTRICT u0         , DOUBLE *RESTRICT gradU0
            , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
            , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres  
            , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lSheat
            , DOUBLE *RESTRICT lViscosity , DOUBLE *RESTRICT lTconductivity
            , DOUBLE *RESTRICT dField     , DOUBLE *RESTRICT wallPar
            , DOUBLE const underU           
            , const short nEn             , short const nFace
            , const short ndm             , INT const nel)
{
  bool fTime, fDisp, fRes, fPresWork, fTemp, fTurb, fWallModel, fKelvin;
  short iCodAdv1, iCodAdv2, iCodDif, wallType, idCell, nf, nCarg1
        , nCarg2, typeTime, fSheat, iCodPolFace;
/*...*/
  INT vizNel;
/*...*/
  DOUBLE thermCoefC, thermCoefV, densityC, densityV, densityM,
         diffEffC, diffEffV, diffEff, sHeatC, sHeatV, sHeatM, 
         eddyViscosityC, eddyViscosityV, viscosityC, tA, coef,
         tmp, tmp1, tmp2, tmp3, prTwall, prTsgs;
  DOUBLE p, sP, dfd, gfKsi, lvSkew[3], alpha, alphaMenosUm;
  DOUBLE v[3], gradUcomp[3], lKsi[3], lNormal[3], gf[3];
  DOUBLE dPviz, lModKsi, lFarea, du, duDksi, lXmcc[3], lXm[3];
  DOUBLE gradUp[3], gradUv[3], ccV[3], gradVelC[3][3], rCell, dt, dt0; 
  DOUBLE uC, uV;
/*... nonOrtogonal*/
  DOUBLE e[3], t[3], s[3], modE, dfdc, xx[3];
/*... */
  DOUBLE presC, presC0, presV, gradPresC[3], gradPresV[3], wfn
        , velC[3], velV[3], dFieldC[3], dFieldV[3], dFieldF[3], cv, cvc;       
/*...*/
  DOUBLE phi,psi,lambda;
  DOUBLE tW,tC;
/*...*/
  DOUBLE pAdv[NPADV];

/*...*/
  idCell      = nFace;
  iCodAdv1    = advT.iCod1;
  iCodAdv2    = advT.iCod2;
  pAdv[0]     = advT.par[0];
  iCodDif     = diffT.iCod;
  iCodPolFace = INTPOLFACELINEAR;
/*...................................................................*/

/*...*/  
  dt         = ddt.dt[0];
  dt0        = ddt.dt[1];
  typeTime   = ddt.type;
  fTime      = ddt.flag;  
  fDisp      = model.fDissipation;

  fRes       = model.fRes;
  fPresWork  = model.fPresWork;
  fTemp      = model.fTemperature;
  fKelvin    = model.fKelvin;
  fTurb      = tModel.fTurb;
  prTwall    = tModel.PrandltTwall;
  prTsgs     = tModel.PrandltTsgs;
  fWallModel = tModel.fWall;
  wallType   = tModel.wallType;
  fSheat     = vProp.fSpecificHeat;
/*...................................................................*/

/*... propriedades da celula*/
  sHeatM      = eddyViscosityC = eddyViscosityV = 0.e0;
  densityC    = lDensity[idCell];
  sHeatC      = lSheat[idCell];
  thermCoefC  = lTconductivity[idCell];  
  viscosityC  = MAT2D(idCell, 0, lViscosity, 2);
  if(fTurb) eddyViscosityC = MAT2D(idCell, 1, lViscosity, 2);
  lambda      = -D2DIV3*(viscosityC + eddyViscosityC);  
/*...................................................................*/

/*...*/
  if(fTemp)
    diffEffC = thermCoefC + sHeatC*eddyViscosityC/prTsgs;
  else
    diffEffC = thermCoefC/sHeatC + eddyViscosityC/prTsgs;
/*...................................................................*/

/*...*/
  gradUp[0]    = MAT2D(idCell, 0, gradU0, 3);
  gradUp[1]    = MAT2D(idCell, 1, gradU0, 3);
  gradUp[2]    = MAT2D(idCell, 2, gradU0, 3);
  velC[0]      = MAT2D(idCell, 0, vel, 3);
  velC[1]      = MAT2D(idCell, 1, vel, 3);
  velC[2]      = MAT2D(idCell, 2, vel, 3);
  uC           = u0[idCell];
/*... p(n)*/
  presC0 = MAT2D(idCell, 0, pres, 2); 
/*... p(n+1)*/
  presC  = MAT2D(idCell, 1, pres, 2);
/*...*/
  gradPresC[0] = MAT2D(idCell, 0, gradPres, 3);
  gradPresC[1] = MAT2D(idCell, 1, gradPres, 3);
  gradPresC[2] = MAT2D(idCell, 2, gradPres, 3);
  dFieldC[0]   = MAT2D(idCell, 0, dField, 3); 
  dFieldC[1]   = MAT2D(idCell, 1, dField, 3);
  dFieldC[2]   = MAT2D(idCell, 2, dField, 3);
/*...................................................................*/

  p = 0.0e0;
  sP = 0.0e0;
  for (nf = 0; nf<nFace; nf++) {
    vizNel = lViz[nf];
    lFarea = fArea[nf];
    lNormal[0] = MAT2D(nf, 0, normal, 3);
    lNormal[1] = MAT2D(nf, 1, normal, 3);
    lNormal[2] = MAT2D(nf, 2, normal, 3);    
    lXmcc[0] = MAT2D(nf, 0, xmcc, 3);
    lXmcc[1] = MAT2D(nf, 1, xmcc, 3);
    lXmcc[2] = MAT2D(nf, 2, xmcc, 3);
/*... dominio*/
    if (vizNel  > -1) {
/*...*/
      densityV = lDensity[nf]; 
      sHeatV   = lSheat[nf];    
      thermCoefV = lTconductivity[nf];  
      if(fTurb) eddyViscosityV = MAT2D(nf, 1, lViscosity, 2);
/*... p(n+1)*/
      presV     = MAT2D(nf, 1, pres, 2);
      gradPresV[0] = MAT2D(nf, 0, gradPres, 3);
      gradPresV[1] = MAT2D(nf, 1, gradPres, 3);
      gradPresV[2] = MAT2D(nf, 2, gradPres, 3);
      dFieldV[0]   = MAT2D(nf, 0, dField, 3);
      dFieldV[1]   = MAT2D(nf, 1, dField, 3);
      dFieldV[2]   = MAT2D(nf, 2, dField, 3);
/*...*/
      uV      = u0[nf];
      velV[0] = MAT2D(nf, 0, vel, 3);
      velV[1] = MAT2D(nf, 1, vel, 3);      
      velV[2] = MAT2D(nf, 2, vel, 3);      
/*...*/
      lKsi[0] = MAT2D(nf, 0, ksi, 3);
      lKsi[1] = MAT2D(nf, 1, ksi, 3);
      lKsi[2] = MAT2D(nf, 2, ksi, 3);
      lModKsi = mKsi[nf];
/*...*/
      lvSkew[0] = MAT2D(nf, 0, vSkew, 3);
      lvSkew[1] = MAT2D(nf, 1, vSkew, 3);
      lvSkew[2] = MAT2D(nf, 2, vSkew, 3);
/*...*/
      duDksi = (uV - uC) / lModKsi;
/*...*/
      gradUv[0] = MAT2D(nf, 0, gradU0, 3);
      gradUv[1] = MAT2D(nf, 1, gradU0, 3);
      gradUv[2] = MAT2D(nf, 2, gradU0, 3);

      ccV[0] = MAT2D(nf, 0, cc, 3);
      ccV[1] = MAT2D(nf, 1, cc, 3);
      ccV[2] = MAT2D(nf, 2, cc, 3);

      lXm[0] = MAT2D(nf, 0, xm, 3);
      lXm[1] = MAT2D(nf, 1, xm, 3);
      lXm[2] = MAT2D(nf, 2, xm, 3);
/*...................................................................*/

/*... termo difusivo
      grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      s[0] = lFarea*lNormal[0];
      s[1] = lFarea*lNormal[1];
      s[2] = lFarea*lNormal[2];
/*...*/
			difusionSchemeNew(s  , lKsi
			                , e  , t
				              , ndm, iCodDif);
/*...................................................................*/

/*...*/
      alpha = interpolFace(lvSkew           ,lXmcc
                          ,volume[idCell]   ,volume[nf]
                          ,lModKsi          ,ndm
                          ,iCodPolFace);
      alphaMenosUm = 1.0e0 - alpha;
/*...................................................................*/

/*... media harmonica*/ 
      if(fTemp){
        diffEffV = thermCoefV + sHeatV*eddyViscosityV/prTsgs;
        diffEff = alpha / diffEffC + alphaMenosUm / diffEffV;
        diffEff = 1.0e0 /  diffEff;
      }
      else {        
        diffEffV = thermCoefV/sHeatV + eddyViscosityV/prTsgs;
        diffEff = alpha / diffEffC + alphaMenosUm / diffEffV;
        diffEff = 1.0e0 / diffEff;
      }
/*...................................................................*/

/*... difusao direta*/
      coef = diffEff;
      modE = sqrt(e[0] * e[0] + e[1] * e[1]  + e[2] * e[2]);
      dfd = coef*modE / lModKsi;
/*...................................................................*/

/*...*/
      gf[0] = alphaMenosUm*gradUp[0] + alpha*gradUv[0];
      gf[1] = alphaMenosUm*gradUp[1] + alpha*gradUv[1];
      gf[2] = alphaMenosUm*gradUp[2] + alpha*gradUv[2];
      dFieldF[0] = alphaMenosUm*dFieldC[0] + alpha*dFieldV[0];
      dFieldF[1] = alphaMenosUm*dFieldC[1] + alpha*dFieldV[1];
      dFieldF[2] = alphaMenosUm*dFieldC[2] + alpha*dFieldV[2];
      densityM = alphaMenosUm*densityC + alpha*densityV;
      if(fTemp)
        sHeatM   = alphaMenosUm*sHeatC + alpha*sHeatV;
/*...................................................................*/

/*... velocidade normal a face*/
      wfn = interpolFaceVel(velC         ,velV
                           ,presC        ,presV
                           ,gradPresC    ,gradPresV
                           ,lNormal      ,lKsi
                           ,lModKsi      ,dFieldF
                           ,alphaMenosUm ,alpha
                           ,ndm);          
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi = gf[0] * lKsi[0] + gf[1] * lKsi[1] + gf[2] * lKsi[2];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
      du = duDksi - gfKsi;
      gradUcomp[0] = gf[0] + du*lKsi[0];
      gradUcomp[1] = gf[1] + du*lKsi[1];
      gradUcomp[2] = gf[2] + du*lKsi[2];
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi = gradUcomp[0] * t[0] 
            + gradUcomp[1] * t[1]
            + gradUcomp[2] * t[2];
/*...................................................................*/

/*... correcao nao-ortogonal*/
      dfdc = coef*gfKsi;
/*...................................................................*/

/*... fluxo convectivo upwind de primeira ordem*/
      if(fTemp)
        cv = sHeatM*densityM*wfn*lFarea;
      else
        cv = densityM*wfn*lFarea;
/*...................................................................*/
 
/*...*/
      v[0] = lXm[0] - ccV[0];
      v[1] = lXm[1] - ccV[1];
      v[2] = lXm[2] - ccV[2];
      advectiveSchemeScalar(uC           ,uV
                           ,gradUp       ,gradUv
                           ,gradUcomp    ,lvSkew
                           ,lXmcc        ,v
                           ,lKsi         ,lModKsi
                           ,wfn          ,&cvc
                           ,alphaMenosUm ,alpha
                           ,pAdv         ,ndm   
                           ,iCodAdv1     ,iCodAdv2);
/*...................................................................*/

/*...*/
      lA[nf] = dfd - min(cv, 0e0);
      sP += cv;
/*... correcao nao ortogonal e do fluxo advectivo*/
      p += dfdc - cv*cvc;
/*...................................................................*/
    }
/*... contorno*/
    else {
      lA[nf] = 0.e0;
      if (lFaceR[nf]) {
        wfn = velC[0] * lNormal[0] 
            + velC[1] * lNormal[1]
            + velC[2] * lNormal[2];
/*...cargas*/
        nCarg1 = lFaceL[nf] - 1;
        nCarg2 = lFaceVelL[nf] - 1;
        xx[0] = MAT2D(nf, 0, xm, 3);
        xx[1] = MAT2D(nf, 1, xm, 3);
        xx[2] = MAT2D(nf, 2, xm, 3);
        pLoadEnergy(&sP          , &p
                  , &tA          , velC
                  , uC           , lNormal  
                  , thermCoefC   , densityC
                  , viscosityC   , sHeatC
                  , prTwall      , xx                   
                  , lFarea       , dcca[nf]
                  , &loads[nCarg1], &loadsVel[nCarg2]
                  , wallPar      , ndm           
                  , true         , fTemp
                  , fKelvin      , fSheat
                  , fWallModel   , wallType);  
/*...................................................................*/
      }
/*...................................................................*/

/*...*/
      else {
/*... inertial sub-layer*/     
         if ( wallPar[0] > 11.81e0 ){
/*... energia na forma da temperatura*/
          if (fTemp) 
            tW = uC;      
/*...................................................................*/
      
/*... energia na forma da entalpia*/
          else{
            tC = specificEnthalpyForTemp(uC ,sHeatC,fSheat,fKelvin); 
            tW = tC;
            tW = tempForSpecificEnthalpy(tW,sHeatC,fSheat,fKelvin); 
          }
/*...................................................................*/
      
/*...*/  
          tmp   = thermCoefC*lFarea /dcca[nf];
          if(!fTemp) tmp /=  sHeatC;
          sP += tmp;
          p  += tmp*tW;      
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... | du1/dx1 du1/dx2 du1/dx3*/
  gradVelC[0][0] = MAT3D(idCell,0,0,gradVel,3,3);
  gradVelC[0][1] = MAT3D(idCell,0,1,gradVel,3,3);
  gradVelC[0][2] = MAT3D(idCell,0,2,gradVel,3,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
  gradVelC[1][0] = MAT3D(idCell,1,0,gradVel,3,3);
  gradVelC[1][1] = MAT3D(idCell,1,1,gradVel,3,3);
  gradVelC[1][2] = MAT3D(idCell,1,2,gradVel,3,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
  gradVelC[2][0] = MAT3D(idCell,2,0,gradVel,3,3);
  gradVelC[2][1] = MAT3D(idCell,2,1,gradVel,3,3);
  gradVelC[2][2] = MAT3D(idCell,2,2,gradVel,3,3);
/*...................................................................*/

/*...*/
  if(fDisp){
    tmp = gradVelC[0][0] + gradVelC[1][1] + gradVelC[2][2];
    psi = tmp*tmp;
/*...*/
    tmp  = gradVelC[0][0]*gradVelC[0][0] 
         + gradVelC[1][1]*gradVelC[1][1]
         + gradVelC[2][2]*gradVelC[2][2];
    tmp1 = gradVelC[0][1] + gradVelC[1][0];
    tmp2 = gradVelC[0][2] + gradVelC[2][0];
    tmp3 = gradVelC[1][2] + gradVelC[2][1];
    phi = 2.0*tmp + tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;    
    tmp = ( lambda*psi + (viscosityC + eddyViscosityC)*phi );
/*...*/
    p+= tmp*volume[idCell];
  }
/*.....................................................................*/ 

/*...*/
  if(fPresWork){
/*... derivada materia da pressao*/
    tmp = (pres[1] - pres[0])/dt 
        + velC[0]*gradPresC[0] 
        + velC[1]*gradPresC[1]
        + velC[2]*gradPresC[2];
/*...*/
    p+= tmp*volume[idCell];
/*.....................................................................*/
  }
/*.....................................................................*/

/*... distretizacao temporal*/
  if(fTemp)
    tmp1 = sHeatC*densityC;
  else
    tmp1 = densityC;

  if (fTime) {
/*... EULER*/
    if (typeTime == EULER)
      sP += tmp1*volume[idCell] / dt;
/*...BACKWARD*/
    else if (typeTime == BACKWARD) {
      tmp = 1.e0 / dt + 1.e0 / (dt + dt0);
      sP += tmp*tmp1*volume[idCell];
    }
  }
/*...................................................................*/

/*...*/ 
  if(nFace == 4){
    lA[idCell] = sP + lA[0] + lA[1] + lA[2] + lA[3];
  }
  else if(nFace == 6){
    lA[idCell] = sP + lA[0] + lA[1] + lA[2] + lA[3] + lA[4] + lA[5];
  }
/*...................................................................*/

/*...*/
  rCell = 0.0e0;
  for (nf = 0; nf<nFace; nf++) {
    if (lViz[nf] > -1) {
/*... pasando os valoeres conhecidos para o lado direito*/
      if (lId[nf] == -2) 
        p += lA[nf] * u0[nf];
      else
/*residuo (R = F-KvizUviz ) e valores prescritos por elemento*/
        rCell += lA[nf] * u0[nf];
    }
  }
/*... residuo: R = F - KpUp*/
  rCell += p - lA[idCell] * u0[idCell];
/*...................................................................*/

/*... under-relaxation(simple)*/
  lA[idCell] = lA[idCell] / underU;
  if(!fRes)
    p += (1.e0 - underU)*lA[idCell] * u0[idCell];
/*...................................................................*/

/*...*/
  if(nFace == 4){
    lA[0] *= -1.e0;
    lA[1] *= -1.e0;
    lA[2] *= -1.e0;
    lA[3] *= -1.e0;
  }
  else if(nFace == 6){
    lA[0] *= -1.e0;
    lA[1] *= -1.e0;
    lA[2] *= -1.e0;
    lA[3] *= -1.e0;
    lA[4] *= -1.e0;
    lA[5] *= -1.e0;
  }
/*...................................................................*/

/*...*/
  if(fRes)
    lB[0]     = rCell;
  else
    lB[0]     = p;
/*...*/
  lRcell[0] = rCell;
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/01/2018                                   *
 * Data de modificaco : 22/01/2018                                   *
 *-------------------------------------------------------------------*
 * cellKinectTurb3D: Celula 3D para transporte da energia kinetica   *
 * turbulenta                                                        * 
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * ldsK    -> definicoes de cargas K                                 *
 * ldsVel  -> definicoes de cargas para velocidades                  *
 * model     -> modelo da equacao de energia                         *
 * advT      -> tecnica da discretizacao do termo advecao            *
 * diffT     -> tecnica da discretizacao do termo difusivo           *
 * tModel    -> modelo de turbulencia                                *
 * vProp     -> propedades variaveis (true|false)                    *
 * lnFace    -> numero de faces da celula central e seus vizinhos    *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lId       -> numeracoes das equacoes das celulas                  *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * fArea     -> area das faces                                       *
 * normal    -> vetores normais as faces das celulas                 *
 * area      -> area da celula central                               *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * vSkew  -> vetor entre o ponto medio a intersecao que une os       *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew -> distacia entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * dcca      -> menor distancia do centroide central a faces desta   *
 *              celula                                               *
 * lDensity  -> massa especifica com variacao temporal               *
 * lSheat    -> calor especifico com variacao temporal               *
 * lDviscosity-> viscosidade dinamica com variacao temporal          *
 * lTconductivity-> condutividade termica com variacao temporal      *
 * dField    -> matriz D do metodo simple                            * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceReK -> restricoes por elemento                                *
 * faceLdK -> carga por elemento                                     *
 * faceReVel-> restricoes por elemento                               *
 * faceLdVel-> carga por elemento                                    *
 * u0        -> solucao conhecida                                    *
 * gradU0    -> gradiente rescontruido da solucao conhecida          *
 * vel       -> campo de velocidade conhecido                        *
 * gradVel   -> gradiente rescontruido da velocidade                 *
 * pres      -> pressao do tempo atual e do tempo anterior           *
 * gradPres  -> gradiente de pressao do tempo atual                  *    
 * density -> massa especifica com variacao temporal                 *
 * dViscosity-> viscosidade dinamica com variacao temporal           *
 * dField    -> matriz D do metodo simple                            * 
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  *
 * cDyn    -> coeficiente dinanmicos ( Ck, Ce)                       * 
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 * lRcell    -> residuo por celula                                   *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellKinectTurb3D(Loads *ldsK         , Loads *ldsVel
            , Turbulence *tModel
            , Advection *advT             , Diffusion *diffT
            , short *RESTRICT lGeomType   , DOUBLE *RESTRICT prop
            , INT *RESTRICT lViz          , INT *RESTRICT lId
            , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
            , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT fArea
            , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume
            , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
            , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
            , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
            , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
            , DOUBLE *RESTRICT lRcell     , Temporal *ddt
            , short  *RESTRICT lFaceReK   , short *RESTRICT lFaceLdK
            , short  *RESTRICT lFaceReVel , short *RESTRICT lFaceLdVel
            , DOUBLE *RESTRICT u0         , DOUBLE *RESTRICT gradU0
            , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
            , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres  
            , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lViscosity 
            , DOUBLE *RESTRICT dField     , DOUBLE *RESTRICT wallPar      
            , DOUBLE *RESTRICT cDyn
            , const short nEn             , short const nFace
            , const short ndm             , INT const nel)
{
  bool fTime, fRes, fTurb, fWallModel, flag, fDynamic;
  short iCodAdv1, iCodAdv2, iCodDif, wallType, idCell, nf, nCarg1
        , nCarg2, typeTime;
/*...*/
  INT vizNel;
/*...*/
  DOUBLE ck,ce,sk, densityC, densityV, densityM;
  DOUBLE diffEffC, diffEffV, diffEff; 
  DOUBLE eddyViscosityC, eddyViscosityV, viscosityC, viscosityV;
  DOUBLE tmp, tmp1, prTwall, prTsgs, tA, coef;
  DOUBLE p, sP, dfd, gfKsi, lvSkew[3], alpha, alphaMenosUm;
  DOUBLE v[3], gradUcomp[3], lKsi[3], lNormal[3], gf[3];
  DOUBLE dPviz, lModKsi, lFarea, du, duDksi, lXmcc[3], lXm[3];
  DOUBLE gradUp[3], gradUv[3], ccV[3], gVelC[3][3], rCell, dt, dt0; 
  DOUBLE uC, uV, ksgsQC;
/*... nonOrtogonal*/
  DOUBLE e[3], t[3], modE, dfdc, xx[3];
/*... */
  DOUBLE presC, presC0, presV, gradPresC[3], gradPresV[3], wfn
        , velC[3], velV[3], dFieldC[3], dFieldV[3], dFieldF[3], cv, cvc
        , s[6], sd[6], delta;       
/*...*/
  DOUBLE pAdv[NPADV];

/*...*/
  idCell   = nFace;
  iCodAdv1 = advT->iCod1;
  iCodAdv2 = advT->iCod2;
  pAdv[0]  = advT->par[0];
  iCodDif  = diffT->iCod;
/*...................................................................*/

/*...*/  
  dt         = ddt->dt[0];
  dt0        = ddt->dt[1];
  typeTime   = ddt->type;
  fTime      = ddt->flag;  

  fDynamic   = tModel->fDynamic;
  flag       = false;
  fRes       = tModel->eK.fRes;
  fTurb      = tModel->fTurb;
  prTwall    = tModel->PrandltTwall;
  prTsgs     = tModel->PrandltTsgs;
  fWallModel = tModel->fWall;
  wallType   = tModel->wallType;
  
  if(fDynamic){
    ck = cDyn[0];
    ce = cDyn[1];
    sk = tModel->eK.sk;
  }
  else {
    ck = tModel->eK.ck;
    ce = tModel->eK.ce;
    sk = tModel->eK.sk;
  }
/*...................................................................*/

/*... propriedades da celula*/
  densityC       = lDensity[idCell];
  viscosityC     = MAT2D(idCell, 0, lViscosity, 2);
  eddyViscosityC = MAT2D(idCell, 1, lViscosity, 2);
/*...................................................................*/

/*...*/
  diffEffC = viscosityC + eddyViscosityC/sk;
/*...................................................................*/

/*...*/
  gradUp[0]    = MAT2D(idCell, 0, gradU0, 3);
  gradUp[1]    = MAT2D(idCell, 1, gradU0, 3);
  gradUp[2]    = MAT2D(idCell, 2, gradU0, 3);
  velC[0]      = MAT2D(idCell, 0, vel, 3);
  velC[1]      = MAT2D(idCell, 1, vel, 3);
  velC[2]      = MAT2D(idCell, 2, vel, 3);
  uC           = u0[idCell];
/*... p(n+1)*/
  presC  = pres[idCell];
/*...*/
  gradPresC[0] = MAT2D(idCell, 0, gradPres, 3);
  gradPresC[1] = MAT2D(idCell, 1, gradPres, 3);
  gradPresC[2] = MAT2D(idCell, 2, gradPres, 3);
/*...*/
  dFieldC[0]   = MAT2D(idCell, 0, dField, 3); 
  dFieldC[1]   = MAT2D(idCell, 1, dField, 3);
  dFieldC[2]   = MAT2D(idCell, 2, dField, 3);
/*...*/ 
  ksgsQC = sqrt(uC);
  delta  = pow(volume[idCell],D1DIV3);
/*...................................................................*/

  p = 0.0e0;
  sP = 0.0e0;
  for (nf = 0; nf<nFace; nf++) {
    vizNel = lViz[nf];
    lFarea = fArea[nf];
    lNormal[0] = MAT2D(nf, 0, normal, 3);
    lNormal[1] = MAT2D(nf, 1, normal, 3);
    lNormal[2] = MAT2D(nf, 2, normal, 3);    
    lXmcc[0] = MAT2D(nf, 0, xmcc, 3);
    lXmcc[1] = MAT2D(nf, 1, xmcc, 3);
    lXmcc[2] = MAT2D(nf, 2, xmcc, 3);
/*... dominio*/
    if (vizNel  > -1) {
/*...*/
      densityV = lDensity[nf]; 
      viscosityV  = MAT2D(nf, 0, lViscosity, 2);
      eddyViscosityV = MAT2D(nf, 1, lViscosity, 2);
/*... p(n+1)*/
      presV        = pres[nf];
      gradPresV[0] = MAT2D(nf, 0, gradPres, 3);
      gradPresV[1] = MAT2D(nf, 1, gradPres, 3);
      gradPresV[2] = MAT2D(nf, 2, gradPres, 3);
      dFieldV[0]   = MAT2D(nf, 0, dField, 3);
      dFieldV[1]   = MAT2D(nf, 1, dField, 3);
      dFieldV[2]   = MAT2D(nf, 2, dField, 3);
/*...*/
      uV      = u0[nf];
      velV[0] = MAT2D(nf, 0, vel, 3);
      velV[1] = MAT2D(nf, 1, vel, 3);      
      velV[2] = MAT2D(nf, 2, vel, 3);      
/*...*/
      lKsi[0] = MAT2D(nf, 0, ksi, 3);
      lKsi[1] = MAT2D(nf, 1, ksi, 3);
      lKsi[2] = MAT2D(nf, 2, ksi, 3);
      lModKsi = mKsi[nf];
/*...*/
      lvSkew[0] = MAT2D(nf, 0, vSkew, 3);
      lvSkew[1] = MAT2D(nf, 1, vSkew, 3);
      lvSkew[2] = MAT2D(nf, 2, vSkew, 3);
/*...*/
      duDksi = (uV - uC) / lModKsi;
/*...*/
      gradUv[0] = MAT2D(nf, 0, gradU0, 3);
      gradUv[1] = MAT2D(nf, 1, gradU0, 3);
      gradUv[2] = MAT2D(nf, 2, gradU0, 3);

      ccV[0] = MAT2D(nf, 0, cc, 3);
      ccV[1] = MAT2D(nf, 1, cc, 3);
      ccV[2] = MAT2D(nf, 2, cc, 3);

      lXm[0] = MAT2D(nf, 0, xm, 3);
      lXm[1] = MAT2D(nf, 1, xm, 3);
      lXm[2] = MAT2D(nf, 2, xm, 3);
/*...................................................................*/

/*... termo difusivo
      grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      difusionScheme(lNormal ,lKsi
                    ,lFarea  ,lModKsi
                    ,e       ,t
                    ,ndm     ,iCodDif);
/*...................................................................*/

/*...*/
      v[0] = lvSkew[0] + lXmcc[0];
      v[1] = lvSkew[1] + lXmcc[1];
      v[2] = lvSkew[2] + lXmcc[2];
      dPviz = sqrt(v[0] * v[0] + v[1] * v[1] + v[2]*v[2]);
      alpha = dPviz / lModKsi;
      alphaMenosUm = 1.0e0 - alpha;
/*...................................................................*/

/*... media harmonica*/ 
      diffEffV = viscosityV + eddyViscosityV/sk;
      diffEff = alpha / diffEffC + alphaMenosUm / diffEffV;
      diffEff = 1.0e0 / diffEff;
/*...................................................................*/

/*... difusao direta*/
      coef = diffEff;
      modE = sqrt(e[0] * e[0] + e[1] * e[1]  + e[2] * e[2]);
      dfd = coef*modE / lModKsi;
/*...................................................................*/

/*...*/
      gf[0] = alphaMenosUm*gradUp[0] + alpha*gradUv[0];
      gf[1] = alphaMenosUm*gradUp[1] + alpha*gradUv[1];
      gf[2] = alphaMenosUm*gradUp[2] + alpha*gradUv[2];
      dFieldF[0] = alphaMenosUm*dFieldC[0] + alpha*dFieldV[0];
      dFieldF[1] = alphaMenosUm*dFieldC[1] + alpha*dFieldV[1];
      dFieldF[2] = alphaMenosUm*dFieldC[2] + alpha*dFieldV[2];
      densityM = alphaMenosUm*densityC + alpha*densityV;
/*...................................................................*/

/*... velocidade normal a face*/
      wfn = interpolFaceVel(velC         ,velV
                           ,presC        ,presV
                           ,gradPresC    ,gradPresV
                           ,lNormal      ,lKsi
                           ,lModKsi      ,dFieldF
                           ,alphaMenosUm ,alpha
                           ,ndm);    
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi = gf[0] * lKsi[0] + gf[1] * lKsi[1] + gf[2] * lKsi[2];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
      du = duDksi - gfKsi;
      gradUcomp[0] = gf[0] + du*lKsi[0];
      gradUcomp[1] = gf[1] + du*lKsi[1];
      gradUcomp[2] = gf[2] + du*lKsi[2];
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi = gradUcomp[0] * t[0] 
            + gradUcomp[1] * t[1]
            + gradUcomp[2] * t[2];
/*...................................................................*/

/*... correcao nao-ortogonal*/
      dfdc = coef*gfKsi;
/*...................................................................*/

/*... fluxo convectivo upwind de primeira ordem*/
      cv = densityM*wfn*lFarea;
/*...................................................................*/
 
/*...*/
      v[0] = lXm[0] - ccV[0];
      v[1] = lXm[1] - ccV[1];
      v[2] = lXm[2] - ccV[2];
      advectiveSchemeScalar(uC           ,uV
                           ,gradUp       ,gradUv
                           ,gradUcomp    ,lvSkew
                           ,lXmcc        ,v
                           ,lKsi         ,lModKsi
                           ,wfn          ,&cvc
                           ,alphaMenosUm ,alpha
                           ,pAdv         ,ndm   
                           ,iCodAdv1     ,iCodAdv2);
/*...................................................................*/

/*...*/
      lA[nf] = dfd - min(cv, 0e0);
      sP += cv;
/*... correcao nao ortogonal e do fluxo advectivo*/
      p += dfdc - cv*cvc;
/*...................................................................*/
    }
/*... contorno*/
    else {
      lA[nf] = 0.e0;
      if (lFaceReK[nf]) {
        wfn = velC[0] * lNormal[0] 
            + velC[1] * lNormal[1]
            + velC[2] * lNormal[2];
/*...cargas*/
        nCarg1 = lFaceLdK[nf] - 1;
        nCarg2 = lFaceLdVel[nf] - 1;
        xx[0] = MAT2D(nf, 0, xm, 3);
        xx[1] = MAT2D(nf, 1, xm, 3);
        xx[2] = MAT2D(nf, 2, xm, 3);
        pLoadOneEqK(&sP          , &p
                  , &tA          , velC
                  , uC           , lNormal  
                  , densityC     , viscosityC 
                  , prTwall      , xx                   
                  , lFarea       , dcca[nf]
                  , &ldsK[nCarg1], &ldsVel[nCarg2]
                  , lFaceReK[nf] , wallPar   
                  , ndm          , true 
                  , fWallModel   , wallType);   
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... | du1/dx1 du1/dx2 du1/dx3*/
  gVelC[0][0] = MAT3D(idCell,0,0,gradVel,3,3);
  gVelC[0][1] = MAT3D(idCell,0,1,gradVel,3,3);
  gVelC[0][2] = MAT3D(idCell,0,2,gradVel,3,3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
  gVelC[1][0] = MAT3D(idCell,1,0,gradVel,3,3);
  gVelC[1][1] = MAT3D(idCell,1,1,gradVel,3,3);
  gVelC[1][2] = MAT3D(idCell,1,2,gradVel,3,3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
  gVelC[2][0] = MAT3D(idCell,2,0,gradVel,3,3);
  gVelC[2][1] = MAT3D(idCell,2,1,gradVel,3,3);
  gVelC[2][2] = MAT3D(idCell,2,2,gradVel,3,3);
/*...................................................................*/

/*...*/
  if(flag){
    tmp = D2DIV3*uC*(gVelC[0][0] + gVelC[1][1] + gVelC[2][2]);
/*...................................................................*/

/*... tensao residual*/
    tensorS(sd,gVelC[0],true);
    tmp1 = -2.e0*eddyViscosityC/densityC;
    sd[0] = tmp1*sd[0];       /*s11*/
    sd[1] = tmp1*sd[1];       /*s22*/
    sd[2] = tmp1*sd[2];       /*s33*/
    sd[3] = tmp1*sd[3];       /*s12*/
    sd[4] = tmp1*sd[4];       /*s13*/ 
    sd[5] = tmp1*sd[5];       /*s23*/
/*....................................................................*/

/*... density*(s:gradVel)*/
    tmp1= sd[0]*gVelC[0][0] 
        + sd[1]*gVelC[1][1] 
        + sd[2]*gVelC[2][2]
        + sd[3]*(gVelC[0][1] + gVelC[1][0])
        + sd[4]*(gVelC[0][2] + gVelC[2][0])
        + sd[5]*(gVelC[1][2] + gVelC[2][1]);
/*....................................................................*/
    p -= densityC*(tmp+tmp1)*volume[idCell];
  }
/*....................................................................*/

/*...*/
  else{
    tensorS(s ,gVelC[0],false);
    tensorS(sd,gVelC[0],true);
    tmp  = D2DIV3*uC*(s[0] + s[1] + s[2]);
/*...................................................................*/

/*... tensao residual*/    
    tmp1 = -2.e0*eddyViscosityC/densityC;
    sd[0] = tmp1*sd[0];   /*s11*/
    sd[1] = tmp1*sd[1];   /*s22*/
    sd[2] = tmp1*sd[2];   /*s33*/
    sd[3] = tmp1*sd[3];   /*s12*/
    sd[4] = tmp1*sd[4];   /*s13*/ 
    sd[5] = tmp1*sd[5];   /*s23*/
/*....................................................................*/

/*... density*(s:sd)*/
    tmp1 = doubleDotSym2(s,sd); 
/*....................................................................*/
    p -= densityC*(tmp+tmp1)*volume[idCell];
  }
/*....................................................................*/

/*... termo nao linear ce*ksg(3/2)/delta*/
  sP += ce*densityC*ksgsQC*volume[idCell]/delta;
/*....................................................................*/

/*... distretizacao temporal*/
  tmp1 = densityC;
  if (fTime) {
/*... EULER*/
    if (typeTime == EULER)
      sP += tmp1*volume[idCell] / dt;
/*...BACKWARD*/
    else if (typeTime == BACKWARD) {
      tmp = 1.e0 / dt + 1.e0 / (dt + dt0);
      sP += tmp*tmp1*volume[idCell];
    }
  }
/*...................................................................*/

/*...*/ 
  if(nFace == 4){
    lA[idCell] = sP + lA[0] + lA[1] + lA[2] + lA[3];
  }
  else if(nFace == 6){
    lA[idCell] = sP + lA[0] + lA[1] + lA[2] + lA[3] + lA[4] + lA[5];
  }
/*...................................................................*/

/*...*/
  rCell = 0.0e0;
  for (nf = 0; nf<nFace; nf++) {
    if (lViz[nf] > -1) {
/*... pasando os valoeres conhecidos para o lado direito*/
      if (lId[nf] == -2) 
        p += lA[nf] * u0[nf];
      else
/*residuo (R = F-KvizUviz ) e valores prescritos por elemento*/
        rCell += lA[nf] * u0[nf];
    }
  }
/*... residuo: R = F - KpUp*/
  rCell += p - lA[idCell] * u0[idCell];
/*...................................................................*/

/*...*/
  if(nFace == 4){
    lA[0] *= -1.e0;
    lA[1] *= -1.e0;
    lA[2] *= -1.e0;
    lA[3] *= -1.e0;
  }
  else if(nFace == 6){
    lA[0] *= -1.e0;
    lA[1] *= -1.e0;
    lA[2] *= -1.e0;
    lA[3] *= -1.e0;
    lA[4] *= -1.e0;
    lA[5] *= -1.e0;
  }
/*...................................................................*/

/*...*/
  if(fRes)
    lB[0]     = rCell;
  else
    lB[0]     = p;
/*...*/
  lRcell[0] = rCell;
/*...................................................................*/
}
/*********************************************************************/
