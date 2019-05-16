#include<CellLib.h>
/********************************************************************
* Data de criacao    : 05/08/2018                                   *
* Data de modificaco : 03/05/2019                                   *
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
* fLump - true                                                      *
* Fuel + Air -> Prod                                                *
* z(0) - Air                                                        *
* z(1) - Fuel                                                       *
* z(2) - Prod                                                       *
* fLump - false                                                     *
* CH4 + O2 + N2 -> CO2 + H2O + N2                                   *
* z(0) - CH4                                                        *
* z(1) - O2                                                         *
* z(2) - N2                                                         *
* z(3) - CO2                                                        *
* z(4) - H2O                                                        *
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
  bool fTime, fRes, fTurb, fWallModel, fLump, fDiffCoor;
  short iCodAdv1, iCodAdv2, iCodDif, wallType, idCell, nf, nCarg1
    , nCarg2, typeTime, iCodPolFace, nComb, nst, i, nSp, nSpLump;
/*...*/
  INT vizNel;
/*...*/
  DOUBLE densityC, densityV, densityM, diffCeofC[MAX_COMB]
        , diffCeofV[MAX_COMB]
        , diffEffC[MAX_COMB], diffEffV[MAX_COMB], diffEff[MAX_COMB]
        , eddyViscosityC, eddyViscosityV, viscosityC
        , tA[MAX_COMB], coef[MAX_COMB]
        , tmp, tmp1, prTwall, prTsgs;
  DOUBLE p[MAX_COMB], sP, sPc[MAX_COMB], dfd[MAX_COMB], gfKsi[MAX_COMB];
  DOUBLE lvSkew[3], alpha, alphaMenosUm;
  DOUBLE v[3], gradUComp[MAX_COMB][3], lKsi[3], lNormal[3], gf[MAX_COMB][3];
  DOUBLE lModKsi, lFarea, du[MAX_COMB], duDksi[MAX_COMB], lXmcc[3], lXm[3];
  DOUBLE gradUp[MAX_COMB][3], gradUv[MAX_COMB][3], ccV[3];
  DOUBLE rCell[MAX_COMB], dt, dt0;
  DOUBLE uC[MAX_COMB], uV[MAX_COMB];
/*... nonOrtogonal*/
  DOUBLE e[3], t[3], s[3], modE, dfdc[MAX_COMB], xx[3];
/*... */
  DOUBLE presC, presC0, presV, gradPresC[3], gradPresV[3], wfn
    , velC[3], velV[3], dFieldC[3], dFieldV[3], dFieldF[3], cv, cvc[MAX_COMB];
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
  fLump      = cModel->fLump;
  nComb      = cModel->nComb;
  nSp        = cModel->nOfSpecies;
  nSpLump    = cModel->nOfSpeciesLump;
  fDiffCoor  = cModel->fCorrectVel;  
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
//if (fTurb) eddyViscosityC = MAT2D(idCell, 1, lViscosity, 2);
/*...................................................................*/

/*...*/
  tmp1 = eddyViscosityC / prTsgs;
  for(i=0;i<nComb;i++)
  {
/*... propriedades da celula*/
    diffCeofC[i] = densityC*MAT2D(idCell, i, lDiff, nComb);

    diffEffC[i] = diffCeofC[i] + tmp1;
/*... | du1/dx1 du1/dx2 du1/dx3 |
      | du2/dx1 du2/dx2 du2/dx3 |
      | du3/dx1 du3/dx2 du3/dx3 |
      | du4/dx1 du4/dx2 du4/dx3 |
*/
    gradUp[i][0] = MAT3D(idCell, i, 0, gradU0, nComb, 3);
    gradUp[i][1] = MAT3D(idCell, i, 1, gradU0, nComb, 3);
    gradUp[i][2] = MAT3D(idCell, i, 2, gradU0, nComb, 3);
/*...*/
    uC[i] = MAT2D(idCell, i, u0, nComb);
/*...*/
    p[i]   = 0.e0;
    sPc[i] = 0.0e0;
  }
/*...................................................................*/

/*...*/
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

/*... correcao dos coefciente de difusao para o transporte de 
      especies*/
   if(fDiffCoor)
     velCorrectCombustion(diffEffC       ,gradUp[0]
                         , velC          ,ndm
                         , nComb );
/*...................................................................*/

/*...*/
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
      for(i=0;i<nComb;i++)
      {
        diffCeofV[i] = densityV*MAT2D(nf, i, lDiff, nComb);
        uV[i]   = MAT2D(nf, i, u0, nComb);
      }
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
      for(i=0;i<nComb;i++)
      {
        duDksi[i] = (uV[i] - uC[i]) / lModKsi;
/*... | du1/dx1 du1/dx2 du1/dx3*/
        gradUv[i][0] = MAT3D(nf, 0, 0, gradU0, nComb, 3);
        gradUv[i][1] = MAT3D(nf, i, 1, gradU0, nComb, 3);
        gradUv[i][2] = MAT3D(nf, i, 2, gradU0, nComb, 3);
      }
/*...................................................................*/

/*...*/
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
      modE = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
      tmp = modE / lModKsi;
      tmp1 = eddyViscosityV / prTsgs;
      for(i=0;i<nComb;i++)
      {
        diffEffV[i] = diffCeofV[i] + tmp1;
        diffEff[i] = alpha / diffEffC[i] + alphaMenosUm / diffEffV[i];
        diffEff[i] = 1.0e0 / diffEff[i];
/*... difusao direta*/
        coef[i] = diffEff[i];
        dfd[i] = coef[i] * tmp;
/*...................................................................*/

/*...*/
        gf[i][0] = alphaMenosUm * gradUp[i][0] + alpha * gradUv[i][0];
        gf[i][1] = alphaMenosUm * gradUp[i][1] + alpha * gradUv[i][1];
        gf[i][2] = alphaMenosUm * gradUp[i][2] + alpha * gradUv[i][2];
      }
/*...................................................................*/

/*...*/
      dFieldF[0] = alphaMenosUm * dFieldC[0] + alpha * dFieldV[0];
      dFieldF[1] = alphaMenosUm * dFieldC[1] + alpha * dFieldV[1];
      dFieldF[2] = alphaMenosUm * dFieldC[2] + alpha * dFieldV[2];
/*...*/
      densityM = alphaMenosUm * densityC + alpha * densityV;
/*...................................................................*/

/*... correcao dos coefciente de difusao para o transporte de 
      especies*/
      if(fDiffCoor)
        velCorrectCombustion(diffEffV       ,gradUv[0]
                            , velV          ,ndm
                            , nComb );
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
      for(i=0;i<nComb;i++)
      {
        gfKsi[i] = gf[i][0] * lKsi[0]
                 + gf[i][1] * lKsi[1]
                 + gf[i][2] * lKsi[2];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
        du[i] = duDksi[i] - gfKsi[i];
/*...*/
        gradUComp[i][0] = gf[i][0] + du[i] * lKsi[0];
        gradUComp[i][1] = gf[i][1] + du[i] * lKsi[1];
        gradUComp[i][2] = gf[i][2] + du[i] * lKsi[2];
/*...................................................................*/

/*... derivadas direcionais*/
        gfKsi[i] = gradUComp[i][0] * t[0]
                 + gradUComp[i][1] * t[1]
                 + gradUComp[i][2] * t[2];
/*...................................................................*/

/*... correcao nao-ortogonal*/
        dfdc[i] = coef[i] * gfKsi[i];
      }
/*...................................................................*/

/*... fluxo convectivo upwind de primeira ordem*/
      cv = densityM*wfn*lFarea;
/*...................................................................*/

/*...*/
      v[0] = lXm[0] - ccV[0];
      v[1] = lXm[1] - ccV[1];
      v[2] = lXm[2] - ccV[2];
      advectiveSchemeNdim(uC      , uV
                    , gradUp[0]   , gradUv[0]
                    , gradUComp[0], lvSkew
                    , lXmcc       , v
                    , lKsi        , lModKsi
                    , cv          , cvc
                    , alphaMenosUm, alpha
                    , pAdv        , ndm
                    , nComb
                    , iCodAdv1    , iCodAdv2);
/*...................................................................*/

/*...*/
      tmp = min(cv, 0e0);
      for(i=0;i<nComb;i++)
      {
        MAT2D( i, nf, lA, nst) = dfd[i] - tmp;
/*... correcao nao ortogonal e do fluxo advectivo*/
        p[i] += dfdc[i] - cv * cvc[i];
      }
      sP += cv;
/*...................................................................*/
    }
/*... contorno*/
    else 
    {
      for(i=0;i<nComb;i++)
        MAT2D( i, nf, lA, nst) = 0.e0;
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
  if(fLump)
  {
    tmp1 = rateFuel*volume[idCell];
    p[SL_FUEL] -= tmp1;
    if( nComb == nSpLump) p[SL_AIR ] -= cModel->sMassAir*tmp1;
    p[SL_PROD] += (1.e0+cModel->sMassAir)*tmp1;
  }
  else
  {
    tmp1 = rateFuel*volume[idCell];
    p[SP_FUEL] -= tmp1;
    p[SP_O2  ] -= cModel->sMassO2*tmp1;
    if( nComb == nSp) p[SP_N2]   += 0.0e0;
    p[SP_CO2]  += cModel->sMassCO2p*tmp1; 
    p[SP_H2O]  += cModel->sMassH2Op*tmp1;
  }
/*...................................................................*/

/*...*/
  for(i=0;i<nComb;i++)
  {
    MAT2D(i, idCell, lA, nst) = sP + sPc[i];
    for (nf = 0; nf<nFace; nf++)
    {
      MAT2D(i, idCell, lA, nst) += MAT2D(i, nf, lA, nst);
    }
  
/*...*/
    rCell[i] = 0.0e0;
    for (nf = 0; nf<nFace; nf++) 
    {
      if (lViz[nf] > -1) 
      {
/*... passando os valores conhecidos para o lado direito*/
        if (lId[nf] == -2)
        {
          p[i] += MAT2D(i, nf,  lA, nst) * MAT2D(nf, i, u0, nComb);
        }  
        else
        {
/*... residuo (R = F-KvizUviz ) e valores prescritos por elemento*/
          rCell[i] += MAT2D( i, nf, lA, nst) * MAT2D(nf, i, u0, nComb);
        }
      }
    }
/*... residuo: R = F - KpUp*/
    rCell[i] += p[i] 
          - MAT2D( i, idCell, lA, nst) * MAT2D(idCell, i, u0, nComb);
/*...................................................................*/

/*... under-relaxation(simple)*/
    MAT2D( i, idCell, lA, nst) /= underU; 
    if (!fRes)
    {
      tmp = (1.e0 - underU);
      p[i] += tmp*MAT2D( i, idCell,  lA, nst)*MAT2D(nf,i,u0,nComb);
    }
/*...................................................................*/

/*...*/
    for (nf = 0; nf<nFace; nf++)
      MAT2D( i, nf, lA, nst) *= -1.e0;
/*...................................................................*/

/*...*/
    if (fRes) 
      lB[i] = rCell[i];
    else 
      lB[i] = p[i];
    lRcell[i] = rCell[i];
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/
