#include<CellLib.h>

/********************************************************************
* Data de criacao    : 05/08/2018                                   *
* Data de modificaco : 12/08/2018                                   *
*-------------------------------------------------------------------*
* cellCombustion3D: Celula 3D para combustao                        *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* loads     -> definicoes de cargas                                 *
* model     -> modelo da equacao de energia                         *
* advT      -> tecnica da discretizacao do termo advecao            *
* diffT     -> tecnica da discretizacao do termo difusivo           *
* tModel    -> modelo de turbulencia                                *
* cModel    -> modelo de combustao                                  *
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
* lDiff     -> calor especifico com variacao temporal               *
* lEddyVisc -> Viscosidade turbulenta                               * 
* dField    -> matriz D do metodo simple                            *
* wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  *
* lRcell    -> nao definido                                         *
* ddt       -> discretizacao temporal                               *
* faceR     -> restricoes por elmento                               *
* faceL     -> carga por elemento                                   *
* u0        -> solucao conhecida                                    *
* gradU0    -> gradiente rescontruido da solucao conhecida          *
* rateFuel  -> taxa de consumo de combustivel                       * 
* vel       -> campo de velocidade conhecido                        *
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
void cellCombustion3D(Loads *loads              , Loads *lVel
                    , Advection *advT           , Diffusion *diffT
                    , Turbulence *tModel        , Combustion *cModel
                    , PropVarFluid *vProp
                    , short *RESTRICT lGeomType , DOUBLE *RESTRICT prop
                    , INT *RESTRICT lViz        , INT *RESTRICT lId
                    , DOUBLE *RESTRICT ksi      , DOUBLE *RESTRICT mKsi
                    , DOUBLE *RESTRICT eta      , DOUBLE *RESTRICT fArea
                    , DOUBLE *RESTRICT normal   , DOUBLE *RESTRICT volume
                    , DOUBLE *RESTRICT xm       , DOUBLE *RESTRICT xmcc
                    , DOUBLE *RESTRICT dcca     , DOUBLE *RESTRICT cc
                    , DOUBLE *RESTRICT vSkew    , DOUBLE *RESTRICT mvSkew
                    , DOUBLE *RESTRICT lA       , DOUBLE *RESTRICT lB
                    , DOUBLE *RESTRICT lRcell   , Temporal const ddt
                    , short  *RESTRICT lFaceR   , short *RESTRICT lFaceL
                    , short  *RESTRICT lFaceVelR, short *RESTRICT lFaceVelL
                    , DOUBLE *RESTRICT u0       , DOUBLE *RESTRICT gradU0
                    , DOUBLE const rateFuel     , DOUBLE *RESTRICT vel
                    , DOUBLE *RESTRICT pres     , DOUBLE *RESTRICT gradPres
                    , DOUBLE *RESTRICT lDensity , DOUBLE *RESTRICT lDiff 
                    , DOUBLE *RESTRICT lEddyVisc
                    , DOUBLE *RESTRICT dField   , DOUBLE *RESTRICT wallPar
                    , DOUBLE const underU
                    , const short nEn           , short const nFace
                    , const short ndm           , INT const nel)
{
  bool fTime, fRes, fTurb, fWallModel;
  short iCodAdv1, iCodAdv2, iCodDif, wallType, idCell, nf, nCarg1
    , nCarg2, typeTime, iCodPolFace, nComb, nst;
/*...*/
  INT vizNel;
/*...*/
  DOUBLE densityC, densityV, densityM, diffCeofC[3], diffCeofV[3], diffCeof[3],
    diffEffC[3], diffEffV[3], diffEff[3],
    eddyViscosityC, eddyViscosityV, viscosityC, tA[3], coef[3],
    tmp, tmp1, tmp2, tmp3, prTwall, prTsgs;
  DOUBLE p[3], sP, sPc[3], dfd[3], gfKsi[3], lvSkew[3], alpha, alphaMenosUm;
  DOUBLE v[3], gradUComp[3][3], lKsi[3], lNormal[3], gf[3][3];
  DOUBLE lModKsi, lFarea, du[3], duDksi[3], lXmcc[3], lXm[3];
  DOUBLE gradUp[3][3], gradUv[3][3], ccV[3], rCell[3], dt, dt0;
  DOUBLE uC[3], uV[3],wf[3];
/*... nonOrtogonal*/
  DOUBLE e[3], t[3], s[3], modE, dfdc[3], xx[3];
/*... */
  DOUBLE presC, presC0, presV, gradPresC[3], gradPresV[3], wfn
    , velC[3], velV[3], dFieldC[3], dFieldV[3], dFieldF[3], cv, cvc[3];
/*...*/
  DOUBLE pAdv[NPADV];

/*...*/
  idCell = nFace;
  nst = nFace + 1;
  iCodAdv1 = advT->iCod1;
  iCodAdv2 = advT->iCod2;
  pAdv[0] = advT->par[0];
  iCodDif = diffT->iCod;
  iCodPolFace = INTPOLFACELINEAR;
/*...................................................................*/

/*...*/
  nComb = cModel->nComb;
/*...................................................................*/

/*...*/
  dt = ddt.dt[0];
  dt0 = ddt.dt[1];
  typeTime = ddt.type;
  fTime = ddt.flag;
  fRes = cModel->fRes;
  fTurb = tModel->fTurb;
  prTwall = tModel->PrandltTwall;
  prTsgs = tModel->PrandltTsgs;
  fWallModel = tModel->fWall;
  wallType = tModel->wallType;
/*...................................................................*/

/*... propriedades da celula*/
  eddyViscosityC = eddyViscosityV = 0.e0;
  densityC = lDensity[idCell];
  diffCeofC[0] = densityC*MAT2D(idCell, 0, lDiff, 3);
  diffCeofC[1] = densityC*MAT2D(idCell, 1, lDiff, 3);
  diffCeofC[2] = densityC*MAT2D(idCell, 2, lDiff, 3);

//if (fTurb) eddyViscosityC = MAT2D(idCell, 1, lViscosity, 2);
  diffEffC[0] = diffCeofC[0] + eddyViscosityC / prTsgs;
  diffEffC[1] = diffCeofC[1] + eddyViscosityC / prTsgs;
  diffEffC[2] = diffCeofC[2] + eddyViscosityC / prTsgs;
/*...................................................................*/

/*... | du1/dx1 du1/dx2 du1/dx3*/
  gradUp[0][0] = MAT3D(idCell, 0, 0, gradU0, 3, 3);
  gradUp[0][1] = MAT3D(idCell, 0, 1, gradU0, 3, 3);
  gradUp[0][2] = MAT3D(idCell, 0, 2, gradU0, 3, 3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
  gradUp[1][0] = MAT3D(idCell, 1, 0, gradU0, 3, 3);
  gradUp[1][1] = MAT3D(idCell, 1, 1, gradU0, 3, 3);
  gradUp[1][2] = MAT3D(idCell, 1, 2, gradU0, 3, 3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
  gradUp[2][0] = MAT3D(idCell, 2, 0, gradU0, 3, 3);
  gradUp[2][1] = MAT3D(idCell, 2, 1, gradU0, 3, 3);
  gradUp[2][2] = MAT3D(idCell, 2, 2, gradU0, 3, 3);
/*...................................................................*/

/*...*/
  uC[0]   = MAT2D(idCell, 0, u0, 3);
  uC[1]   = MAT2D(idCell, 1, u0, 3);
  uC[2]   = MAT2D(idCell, 2, u0, 3);
  velC[0] = MAT2D(idCell, 0, vel, 3);
  velC[1] = MAT2D(idCell, 1, vel, 3);
  velC[2] = MAT2D(idCell, 2, vel, 3);
/*... p(n)*/
  presC0 = MAT2D(idCell, 0, pres, 2);
/*... p(n+1)*/
  presC = MAT2D(idCell, 1, pres, 2);
/*...*/
  gradPresC[0] = MAT2D(idCell, 0, gradPres, 3);
  gradPresC[1] = MAT2D(idCell, 1, gradPres, 3);
  gradPresC[2] = MAT2D(idCell, 2, gradPres, 3);
  dFieldC[0] = MAT2D(idCell, 0, dField, 3);
  dFieldC[1] = MAT2D(idCell, 1, dField, 3);
  dFieldC[2] = MAT2D(idCell, 2, dField, 3);
/*...................................................................*/

/*...*/
  p[0] = p[1] = p[2] = 0.e0;
  sPc[0] = sPc[1] = sPc[2] = 0.0e0;
  sP = 0.0e0;
  for (nf = 0; nf<nFace; nf++) 
  {
    vizNel = lViz[nf];
    lFarea = fArea[nf];
    lNormal[0] = MAT2D(nf, 0, normal, 3);
    lNormal[1] = MAT2D(nf, 1, normal, 3);
    lNormal[2] = MAT2D(nf, 2, normal, 3);
    lXmcc[0] = MAT2D(nf, 0, xmcc, 3);
    lXmcc[1] = MAT2D(nf, 1, xmcc, 3);
    lXmcc[2] = MAT2D(nf, 2, xmcc, 3);
/*... dominio*/
    if (vizNel  > -1) 
    {
/*...*/
      densityV   = lDensity[nf];
      diffCeofV[0] = densityV*MAT2D(nf, 0, lDiff, 3);
      diffCeofV[1] = densityV*MAT2D(nf, 1, lDiff, 3);
      diffCeofV[2] = densityV*MAT2D(nf, 2, lDiff, 3);
//    if (fTurb) eddyViscosityV = MAT2D(nf, 1, lViscosity, 2);
/*... p(n+1)*/
      presV = MAT2D(nf, 1, pres, 2);
      gradPresV[0] = MAT2D(nf, 0, gradPres, 3);
      gradPresV[1] = MAT2D(nf, 1, gradPres, 3);
      gradPresV[2] = MAT2D(nf, 2, gradPres, 3);
      dFieldV[0] = MAT2D(nf, 0, dField, 3);
      dFieldV[1] = MAT2D(nf, 1, dField, 3);
      dFieldV[2] = MAT2D(nf, 2, dField, 3);
/*...*/
      uV[0]   = MAT2D(nf, 0, u0, 3);
      uV[1]   = MAT2D(nf, 1, u0, 3);
      uV[2]   = MAT2D(nf, 2, u0, 3);
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
      duDksi[0] = (uV[0] - uC[0]) / lModKsi;
      duDksi[1] = (uV[1] - uC[1]) / lModKsi;
      duDksi[2] = (uV[2] - uC[2]) / lModKsi;
/*... | du1/dx1 du1/dx2 du1/dx3*/
      gradUv[0][0] = MAT3D(nf, 0, 0, gradU0, 3, 3);
      gradUv[0][1] = MAT3D(nf, 0, 1, gradU0, 3, 3);
      gradUv[0][2] = MAT3D(nf, 0, 2, gradU0, 3, 3);
/*... | du2/dx1 du2/dx2 du2/dx3*/
      gradUv[1][0] = MAT3D(nf, 1, 0, gradU0, 3, 3);
      gradUv[1][1] = MAT3D(nf, 1, 1, gradU0, 3, 3);
      gradUv[1][2] = MAT3D(nf, 1, 2, gradU0, 3, 3);
/*... | du3/dx1 du3/dx2 du3/dx3*/
      gradUv[2][0] = MAT3D(nf, 2, 0, gradU0, 3, 3);
      gradUv[2][1] = MAT3D(nf, 2, 1, gradU0, 3, 3);
      gradUv[2][2] = MAT3D(nf, 2, 2, gradU0, 3, 3);
/*...................................................................*/

      ccV[0] = MAT2D(nf, 0, cc, 3);
      ccV[1] = MAT2D(nf, 1, cc, 3);
      ccV[2] = MAT2D(nf, 2, cc, 3);

      lXm[0] = MAT2D(nf, 0, xm, 3);
      lXm[1] = MAT2D(nf, 1, xm, 3);
      lXm[2] = MAT2D(nf, 2, xm, 3);
/*...................................................................*/

/*... termo difusivo
      grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      s[0] = lFarea * lNormal[0];
      s[1] = lFarea * lNormal[1];
      s[2] = lFarea * lNormal[2];
/*...*/
      difusionSchemeNew(s, lKsi, e, t, ndm, iCodDif);
/*...................................................................*/

/*...*/
      alpha = interpolFace(lvSkew        , lXmcc
                         , volume[idCell], volume[nf]
                         , lModKsi       , ndm
                         , iCodPolFace);
      alphaMenosUm = 1.0e0 - alpha;
/*...................................................................*/

/*... media harmonica*/
      diffEffV[0] = diffCeofV[0] + eddyViscosityV / prTsgs;
      diffEffV[1] = diffCeofV[1] + eddyViscosityV / prTsgs;
      diffEffV[2] = diffCeofV[2] + eddyViscosityV / prTsgs;
      diffEff[0] = alpha / diffEffC[0] + alphaMenosUm / diffEffV[0];
      diffEff[1] = alpha / diffEffC[1] + alphaMenosUm / diffEffV[1];
      diffEff[2] = alpha / diffEffC[2] + alphaMenosUm / diffEffV[2];
      diffEff[0] = 1.0e0 / diffEff[0];
      diffEff[1] = 1.0e0 / diffEff[1];
      diffEff[2] = 1.0e0 / diffEff[2];
/*...................................................................*/

/*... difusao direta*/
      coef[0] = diffEff[0];
      coef[1] = diffEff[1];
      coef[2] = diffEff[2];
      modE = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
      tmp = modE / lModKsi;
      dfd[0] = coef[0] * tmp;
      dfd[1] = coef[1] * tmp;
      dfd[2] = coef[2] * tmp;
/*...................................................................*/

/*...*/
      gf[0][0] = alphaMenosUm * gradUp[0][0] + alpha * gradUv[0][0];
      gf[0][1] = alphaMenosUm * gradUp[0][1] + alpha * gradUv[0][1];
      gf[0][2] = alphaMenosUm * gradUp[0][2] + alpha * gradUv[0][2];
/*...*/
      gf[1][0] = alphaMenosUm * gradUp[1][0] + alpha * gradUv[1][0];
      gf[1][1] = alphaMenosUm * gradUp[1][1] + alpha * gradUv[1][1];
      gf[1][2] = alphaMenosUm * gradUp[1][2] + alpha * gradUv[1][2];
/*...*/
      gf[2][0] = alphaMenosUm * gradUp[2][0] + alpha * gradUv[2][0];
      gf[2][1] = alphaMenosUm * gradUp[2][1] + alpha * gradUv[2][1];
      gf[2][2] = alphaMenosUm * gradUp[2][2] + alpha * gradUv[2][2];
/*...*/
      dFieldF[0] = alphaMenosUm * dFieldC[0] + alpha * dFieldV[0];
      dFieldF[1] = alphaMenosUm * dFieldC[1] + alpha * dFieldV[1];
      dFieldF[2] = alphaMenosUm * dFieldC[2] + alpha * dFieldV[2];
/*...*/
      densityM = alphaMenosUm * densityC + alpha * densityV;
/*...................................................................*/

/*... velocidade normal a face*/
      wfn = interpolFaceVel(velC       , velV
                         , presC       , presV
                         , gradPresC   , gradPresV
                         , lNormal     , lKsi
                         , lModKsi     , dFieldF
                         , alphaMenosUm, alpha
                         , ndm);
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi[0] = gf[0][0] * lKsi[0]
               + gf[0][1] * lKsi[1]
               + gf[0][2] * lKsi[2];
/*...*/
      gfKsi[1] = gf[1][0] * lKsi[0]
               + gf[1][1] * lKsi[1]
               + gf[1][2] * lKsi[2];
/*...*/
      gfKsi[2] = gf[2][0] * lKsi[0]
               + gf[2][1] * lKsi[1]
               + gf[2][2] * lKsi[2];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
      du[0] = duDksi[0] - gfKsi[0];
      du[1] = duDksi[1] - gfKsi[1];
      du[2] = duDksi[2] - gfKsi[2];
/*...*/
      gradUComp[0][0] = gf[0][0] + du[0] * lKsi[0];
      gradUComp[0][1] = gf[0][1] + du[0] * lKsi[1];
      gradUComp[0][2] = gf[0][2] + du[0] * lKsi[2];
/*...*/
      gradUComp[1][0] = gf[1][0] + du[1] * lKsi[0];
      gradUComp[1][1] = gf[1][1] + du[1] * lKsi[1];
      gradUComp[1][2] = gf[1][2] + du[1] * lKsi[2];
/*...*/
      gradUComp[2][0] = gf[2][0] + du[2] * lKsi[0];
      gradUComp[2][1] = gf[2][1] + du[2] * lKsi[1];
      gradUComp[2][2] = gf[2][2] + du[2] * lKsi[2];
/*...................................................................*/

//*... derivadas direcionais*/
/*...*/
      gfKsi[0] = gradUComp[0][0] * t[0]
               + gradUComp[0][1] * t[1]
               + gradUComp[0][2] * t[2];
/*...*/
      gfKsi[1] = gradUComp[1][0] * t[0]
               + gradUComp[1][1] * t[1]
               + gradUComp[1][2] * t[2];
/*...*/
      gfKsi[2] = gradUComp[2][0] * t[0]
               + gradUComp[2][1] * t[1]
               + gradUComp[2][2] * t[2];
/*...................................................................*/

/*... correcao nao-ortogonal*/
      dfdc[0] = coef[0] * gfKsi[0];
      dfdc[1] = coef[1] * gfKsi[1];
      dfdc[2] = coef[2] * gfKsi[2];
/*...................................................................*/

/*... fluxo convectivo upwind de primeira ordem*/
      cv = densityM*wfn*lFarea;
/*...................................................................*/

/*...*/
      v[0] = lXm[0] - ccV[0];
      v[1] = lXm[1] - ccV[1];
      v[2] = lXm[2] - ccV[2];
      advectiveScheme(uC          , uV
                    , gradUp[0]   , gradUv[0]
                    , gradUComp[0], lvSkew
                    , lXmcc       , v
                    , lKsi        , lModKsi
                    , cv          , cvc
                    , alphaMenosUm, alpha
                    , pAdv        , ndm
                    , iCodAdv1    , iCodAdv2);
/*...................................................................*/

/*...*/
      tmp = min(cv, 0e0);
      MAT2D( 0, nf, lA, nst) = dfd[0] - tmp;
      MAT2D( 1, nf, lA, nst) = dfd[1] - tmp;
      MAT2D( 2, nf, lA, nst) = dfd[2] - tmp;
      sP += cv;
/*... correcao nao ortogonal e do fluxo advectivo*/
      p[0] += dfdc[0] - cv * cvc[0];
      p[1] += dfdc[1] - cv * cvc[1];
      p[2] += dfdc[2] - cv * cvc[2];
/*...................................................................*/
    }
/*... contorno*/
    else 
    {
      MAT2D( 0, nf, lA, nst) = 0.e0;
      MAT2D( 1, nf, lA, nst) = 0.e0;
      MAT2D( 2, nf, lA, nst) = 0.e0;
      if (lFaceR[nf]) 
      {
        wfn = velC[0] * lNormal[0]
            + velC[1] * lNormal[1]
            + velC[2] * lNormal[2];
/*...cargas*/
        nCarg1 = lFaceL[nf] - 1;
        nCarg2 = lFaceVelL[nf] - 1;
        xx[0] = MAT2D(nf, 0, xm, 3);
        xx[1] = MAT2D(nf, 1, xm, 3);
        xx[2] = MAT2D(nf, 2, xm, 3);
        pLoadCombustion(vProp
                      , sPc           , p
                      , tA            , velC
                      , uC            , lNormal
                      , diffCeofC     , densityC
                      , densityC
                      , prTwall       , xx
                      , lFarea        , dcca[nf]
                      , &loads[nCarg1], &loadsVel[nCarg2]
                      , wallPar       , ndm
                      , true          , fWallModel
                      , nComb         , wallType);    
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... distretizacao temporal*/
  if (fTime) 
  {
/*... EULER*/
    if (typeTime == EULER)
      sP += densityC * volume[idCell] / dt;
/*...BACKWARD*/
    else if (typeTime == BACKWARD) 
    {
      tmp = 1.e0 / dt + 1.e0 / (dt + dt0);
      sP += tmp * densityC*volume[idCell];
    }
  }
/*...................................................................*/

/*... reacao*/
  tmp  = cModel->sMass;
  tmp1 = rateFuel*volume[idCell];
  p[0] -= tmp*tmp1;
  p[1] -= tmp1;
  p[2] += (1.e0+tmp)*tmp1;
/*...................................................................*/

/*...*/
  MAT2D(0, idCell, lA, nst) = sP + sPc[0];
  MAT2D(1, idCell, lA, nst) = sP + sPc[1];
  MAT2D(2, idCell, lA, nst) = sP + sPc[2];
  for (nf = 0; nf<nFace; nf++)
  {
    MAT2D(0, idCell, lA, nst) += MAT2D(0, nf, lA, nst);
    MAT2D(1, idCell, lA, nst) += MAT2D(1, nf, lA, nst);
    MAT2D(2, idCell, lA, nst) += MAT2D(2, nf, lA, nst);
  }
/*...................................................................*/

/*...*/
  rCell[0] = 0.0e0;
  rCell[1] = 0.0e0;
  rCell[2] = 0.0e0;
  for (nf = 0; nf<nFace; nf++) 
  {
    if (lViz[nf] > -1) 
    {
/*... passando os valores conhecidos para o lado direito*/
      if (lId[nf] == -2)
      {
        p[0] += MAT2D(0, nf,  lA, nst) * MAT2D(nf, 0, u0, 3);
        p[1] += MAT2D(1, nf,  lA, nst) * MAT2D(nf, 1, u0, 3);
        p[2] += MAT2D(2, nf,  lA, nst) * MAT2D(nf, 2, u0, 3);
      }  
      else
      {
/*... residuo (R = F-KvizUviz ) e valores prescritos por elemento*/
        rCell[0] += MAT2D( 0, nf, lA, nst) * MAT2D(nf, 0, u0, 3);
        rCell[1] += MAT2D( 1, nf, lA, nst) * MAT2D(nf, 1, u0, 3);
        rCell[2] += MAT2D( 2, nf, lA, nst) * MAT2D(nf, 2, u0, 3);
      }
    }
  }
/*... residuo: R = F - KpUp*/
  rCell[0] += p[0] - MAT2D( 0, idCell, lA, nst) * MAT2D(idCell, 0, u0, 3);
  rCell[1] += p[1] - MAT2D( 1, idCell, lA, nst) * MAT2D(idCell, 1, u0, 3);
  rCell[2] += p[2] - MAT2D( 2, idCell, lA, nst) * MAT2D(idCell, 2, u0, 3);
/*...................................................................*/

/*... under-relaxation(simple)*/
  MAT2D( 0, idCell, lA, nst) /= underU; 
  MAT2D( 1, idCell, lA, nst) /= underU;
  MAT2D( 2, idCell, lA, nst) /= underU;
  if (!fRes)
  {
    tmp = (1.e0 - underU);
    p[0] += tmp*MAT2D( 0, idCell,  lA, nst)*MAT2D(nf,0,u0,3);
    p[1] += tmp*MAT2D( 1, idCell,  lA, nst)*MAT2D(nf,1,u0,3);
    p[2] += tmp*MAT2D( 2, idCell,  lA, nst)*MAT2D(nf,2,u0,3);
  }
/*...................................................................*/

/*...*/
  for (nf = 0; nf<nFace; nf++)
  {
    MAT2D( 0, nf, lA, nst) *= -1.e0;
    MAT2D( 1, nf, lA, nst) *= -1.e0;
    MAT2D( 2, nf, lA, nst) *= -1.e0;
  }   
/*...................................................................*/

/*...*/
  if (fRes) 
  {
    lB[0] = rCell[0];
    lB[1] = rCell[1];
    lB[2] = rCell[2];
  }
  else 
  {
    lB[0] = p[0];
    lB[1] = p[1];
    lB[2] = p[2];
  }
  lRcell[0] = rCell[0];
  lRcell[1] = rCell[1];
  lRcell[2] = rCell[2];
/*...................................................................*/
}
/*********************************************************************/
