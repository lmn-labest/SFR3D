#include<CellLib.h>
/*********************************************************************
 * Biblioteca de celulas                                             *
 *-------------------------------------------------------------------*
 * Celulas 2D                                                        *
 *-------------------------------------------------------------------*
 * ------------------ INCOMPRESSIVEl ------------------------------- *
 *-------------------------------------------------------------------*
 * Celulas 3D                                                        *
 *-------------------------------------------------------------------*
 * ------------------ INCOMPRESSIVEl ------------------------------- *
 *                                                                   *
 * CellSimpleVel3D: Celula 3D para equacao de correcao de pressao    *
 * metodo simple em escoamento imcompressivel                        *
 *                                                                   *
 * ------------------ LEVEMENTE COMPRESSIVEl ------------------------*
 *                                                                   *
 * CellSimpleVel3DLm: Celula 3D para equacao de correcao de pressao  *
 * metodo simple em escoamento para baixo mach                       *
 *                                                                   *
 *-------------------------------------------------------------------*
*********************************************************************/
/*********************************************************************
 * Data de criacao    : 11/07/2016                                   *
 * Data de modificaco : 28/09/2019                                   * 
 *-------------------------------------------------------------------* 
 * CELLSIMPLEVE3D: Celula 3D para velocidade do metodo simple        * 
 * em escoamento imcompressivel                                      * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loadsVel  -> definicoes de cargas de velocidades                  * 
 * loadsPres -> definicoes de cargas de pressao                      * 
 * advVel    -> tecnica da discretizacao do termo advecao            *
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * tModel    -> modelo de turbulencia                                *
 * momentumModel -> termos/modelos da equacao de momento linear      *
 * typeSimple-> tipo do metodo simple                                *
 * lnFace    -> numero de faces da celula central e seus vizinhos    * 
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lprop     -> propriedade fisicas das celulas                      * 
 * lViz      -> viznhos da celula central                            * 
 * lId       -> numeracoes das equacoes das celulas                  * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *              vizinhos desta                                       * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area da face                                         * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> volume da celula central                             * 
 * xm        -> pontos medios das faces da celula central            * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *              faces da celula central                              * 
 * dcca      -> menor distancia do centroide central a faces desta   *
 *              celula                                               *
 * cc        -> centroides da celula centra e seus vizinhos          * 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os * 
 *              centrois compartilhado nessa face da celula central  * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceVelR  -> restricoes por elemento de velocidades               * 
 * facePresR -> restricoes por elemento de pressao                   * 
 * pres      -> campo de pressao conhecido                           * 
 * gradPes   -> gradiente reconstruido da pressao                    * 
 * vel       -> campo de velocidade conhecido                        * 
 * gradVel   -> gradiente rescontruido das velocidades (Transposto)  * 
 * dField    -> matriz D do metodo simple                            * 
 * stressR   -> tensao residual estrutural                           *
 * eddyVisc  -> viscosidade turbulenta                               *
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  *
 * nEn       -> numero de nos da celula central                      *
 * underU    -> parametro de sobre relaxamento                       * 
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
 *                                                                   *  
 * gradVel(nFace+1,ndf,ndm) -> gradVel(nFace+1,2,2)                  * 
 *                                                                   *  
 *                  | du1dx1 du1dx2 du1dx3 |                         *
 * grad(viz1,*,*) = | du2dx1 du2dx2 du2dx3 | = djui                  *
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
void cellSimpleVel3D(Loads *lVel         ,Loads *lPres 
            ,Advection *advVel           ,Diffusion *diffVel
            ,Turbulence *tModel          ,MomentumModel *momentumModel
            ,short const typeSimple 
            ,short *RESTRICT lGeomType   ,DOUBLE *RESTRICT prop
            ,INT *RESTRICT lViz          ,INT *RESTRICT lId  
            ,DOUBLE *RESTRICT ksi        ,DOUBLE *RESTRICT mKsi
            ,DOUBLE *RESTRICT eta        ,DOUBLE *RESTRICT fArea
            ,DOUBLE *RESTRICT normal     ,DOUBLE *RESTRICT volume 
            ,DOUBLE *RESTRICT xm         ,DOUBLE *RESTRICT xmcc
            ,DOUBLE *RESTRICT dcca       ,DOUBLE *RESTRICT cc
            ,DOUBLE *RESTRICT vSkew      ,DOUBLE *RESTRICT mvSkew
            ,DOUBLE *RESTRICT lA         ,DOUBLE *RESTRICT lB
            ,DOUBLE *RESTRICT lRcell     ,Temporal *ddt
            ,short  *RESTRICT lFaceVelR  ,short  *RESTRICT lFacePresR
            ,DOUBLE *RESTRICT pres       ,DOUBLE *RESTRICT gradPres 
            ,DOUBLE *RESTRICT vel        ,DOUBLE *RESTRICT gradVel
            ,DOUBLE *RESTRICT dField     ,DOUBLE *RESTRICT stressR
            ,DOUBLE *RESTRICT eddyVisc   ,DOUBLE *RESTRICT wallPar 
            ,DOUBLE const underU         
            ,const short nEn             ,short const nFace    
            ,const short ndm             ,INT const nel)
{ 
/*...*/
  bool fTime, fRes, fTurb, fRhieInt, fWallModel, fStruc;
  bool fViscosity,fSoPressure;
/*...*/
  short iCodAdv1, iCodAdv2, iCodDif, wallType, idCell, nf
       ,nCarg   , typeTime, iCodPolFace, i, typeFacePres;
/*...*/
  INT vizNel;
/*...*/
  DOUBLE viscosityC    ,viscosityV     ,viscosity
       , eddyViscosityC, eddyViscosityV, effViscosityC, effViscosityV
       , densityC      ,densityV       ,density;
/*...*/
  DOUBLE rCell[3], dt, dt0, p[3], sP, sPc[3], tA[3]
        ,dFieldC[3], dFieldV[3], dFieldF[3];
/*...*/
  DOUBLE v[3],lKsi[3],lNormal[3],lXmcc[3],ccV[3];
  DOUBLE lModKsi,lFarea,du[3],duDksi[3],lXm[3];
  DOUBLE coef,lAn,tmp;
/*...*/
  DOUBLE dfd,cv,cvc[3],lvSkew[3];
/*... nonOrtogonal*/
	DOUBLE e[3], t[3], modE, dfdc[3];
/*... interpolacao linear*/
  DOUBLE alpha,alphaMenosUm;
/*... */
  DOUBLE pFace,pf[3],presC,presV,gradPresC[3],gradPresV[3];
/*... */
  DOUBLE gf[3][3],gfKsi[3];
  DOUBLE wfn,velC[3],velV[3],wf[3];
  DOUBLE gradVelC[3][3],gradVelV[3][3],gradVelComp[3][3];
/*...*/
  DOUBLE stressRc[6],stressRv[6],s[6];
/*...*/
  DOUBLE xx[4],ts;
/*...*/
  DOUBLE pAdv[NPADV];

/*...*/
  idCell      = nFace;
  iCodAdv1    = advVel->iCod1;
  iCodAdv2    = advVel->iCod2;
  pAdv[0]     = advVel->par[0];
  iCodDif     = diffVel->iCod;
  iCodPolFace = INTPOLFACELINEAR;
/*...................................................................*/

/*...*/  
  fWallModel       = tModel->fWall;
  wallType         = tModel->wallType;
  fTurb            = tModel->fTurb;
  fStruc           = tModel->fTurbStruct;
/*...*/  
  fSoPressure = momentumModel->fSoPressure;
  fRhieInt    = momentumModel->fRhieChowInt;
  fRes        = momentumModel->fRes;
  fViscosity  = momentumModel->fViscosity;
//typeFacePres     = INT_BUOYANT_FORCE;
  if(fSoPressure)
    typeFacePres   = INT_PRESSURE_SO;
  else
    typeFacePres   = INT_PRESSURE_FO;  
/*...................................................................*/

/*...*/
  ts       = ddt->t; 
  dt       = ddt->dt[0];
  dt0      = ddt->dt[1];
  typeTime = ddt->type;
  fTime    = ddt->flag;
/*...................................................................*/

/*... propriedades da celula*/
  densityC   = MAT2D(idCell, DENSITY         ,prop,MAXPROP);
  viscosityC = MAT2D(idCell, DYNAMICVISCOSITY,prop,MAXPROP);
/*...................................................................*/

/*...*/
  eddyViscosityC = eddyViscosityV = 0.0e0;
  if(fTurb) eddyViscosityC= eddyVisc[idCell];
  effViscosityC = viscosityC + eddyViscosityC;
/*...................................................................*/

/*...*/
  sP    = 0.e0;
  presC = pres[idCell];
  for (i = 0; i < 3; i++)
  {
/*... | dui/dx1 dui/dx2 dui/dx3*/
    gradVelC[i][0] = MAT3D(idCell,i,0,gradVel,3,3);
    gradVelC[i][1] = MAT3D(idCell,i,1,gradVel,3,3);
    gradVelC[i][2] = MAT3D(idCell,i,2,gradVel,3,3);
/*...*/
    gradPresC[i]  = MAT2D(idCell, i, gradPres, 3);
    dFieldC[i]    = MAT2D(idCell, i, dField  , 3); 
/*...*/
    velC[i]  = MAT2D(idCell,i,vel,3);
/*...................................................................*/

/*...*/
    p[i]   = 0.e0;
    pf[i]  = 0.0e0;
    sPc[i] = 0.0e0;
/*...................................................................*/
  }
/*...................................................................*/

/*... tensor residual*/
  if(fStruc)
    for (i = 0; i < 6; i++)
      stressRc[i] = MAT2D(idCell,i,stressR,6);
/*...................................................................*/

/*....*/
  for(nf=0;nf<nFace;nf++)
  {
    vizNel  = lViz[nf];
    lFarea  = fArea[nf];
    lModKsi = mKsi[nf];
/*..*/
    for (i = 0; i < 3; i++)
    {
      lNormal[i] = MAT2D(nf,i,normal,ndm);
      lXmcc[i]   = MAT2D(nf,i,xmcc,ndm);
      lKsi[i] = MAT2D(nf,i,ksi,3);
    }
/*...................................................................*/

/*... dominio*/
    if( vizNel  > -1 )
    {
/*...*/
      densityV       = MAT2D(nf, DENSITY         ,prop,MAXPROP);
      viscosityV     = MAT2D(nf, DYNAMICVISCOSITY,prop,MAXPROP);
      if (fTurb) eddyViscosityV = eddyVisc[nf];
/*...*/
      presV        = pres[nf]; 
      for(i=0;i<3;i++)
      {
        velV[i]        = MAT2D(nf,i,vel,ndm);
/*...*/
        gradPresV[i] = MAT2D(nf, i, gradPres, 3);
        dFieldV[i]   = MAT2D(nf, i, dField, 3);
/*..*/
        lvSkew[i]      = MAT2D(nf,i,vSkew,3);
/*..*/
        duDksi[i]      = (velV[i] - velC[i]) / lModKsi;
/*... | dui/dx1 dui/dx2 dui/dx3|*/
        gradVelV[i][0] = MAT3D(nf,i,0,gradVel,3,3);
        gradVelV[i][1] = MAT3D(nf,i,1,gradVel,3,3);
        gradVelV[i][2] = MAT3D(nf,i,2,gradVel,3,3);
/*...*/
        ccV[i] = MAT2D(nf, i, cc, 3);
        lXm[i] = MAT2D(nf, i, xm, 3);
/*...................................................................*/
      }
/*...................................................................*/

/*... tensor residual*/
      if(fStruc)
        for(i=0;i<6;i++)
          stressRv[i] = MAT2D(nf,i,stressR,6);
/*...................................................................*/
 
/*... termo difusivo
      grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
/*...*/
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
      effViscosityV = viscosityV + eddyViscosityV;
      viscosity = mediaHarmonica(effViscosityV,effViscosityC
                                ,alphaMenosUm ,alpha        );
/*...................................................................*/

/*... difusao direta*/
      coef = viscosity;
      modE = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
      dfd  = coef*modE/lModKsi;
/*...................................................................*/
      
/*...*/
      for(i=0;i<3;i++)
      {
        gf[i][0] = alphaMenosUm*gradVelC[i][0] + alpha*gradVelV[i][0];
        gf[i][1] = alphaMenosUm*gradVelC[i][1] + alpha*gradVelV[i][1];
        gf[i][2] = alphaMenosUm*gradVelC[i][2] + alpha*gradVelV[i][2];
/*...*/
        wf[i]    = alphaMenosUm*velC[i] + alpha*velV[i];
      }
/*...................................................................*/
      density  = alphaMenosUm*densityC  + alpha*densityV;
/*...................................................................*/

/*... velocidade normal a face*/
      if(fRhieInt) 
      {
        dFieldF[0] = alphaMenosUm*dFieldC[0] + alpha*dFieldV[0];
        dFieldF[1] = alphaMenosUm*dFieldC[1] + alpha*dFieldV[1];
        dFieldF[2] = alphaMenosUm*dFieldC[2] + alpha*dFieldV[2];
        wfn = interpolFaceVel(velC         ,velV
                             ,presC        ,presV
                             ,gradPresC    ,gradPresV
                             ,lNormal      ,lKsi
                             ,lModKsi      ,dFieldF
                             ,alphaMenosUm ,alpha
                             ,ndm);
      }
      else
        wfn   = wf[0]*lNormal[0] + wf[1]*lNormal[1] + wf[2]*lNormal[2];    
/*...................................................................*/

/*...*/
      for(i=0;i<3;i++)
      {
/*... derivadas direcionais*/
        gfKsi[i] = gf[i][0]*lKsi[0] 
                 + gf[i][1]*lKsi[1] 
                 + gf[i][2]*lKsi[2];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
        du[i]             = duDksi[i] - gfKsi[i];
/*...*/
        gradVelComp[i][0] = gf[i][0] + du[i]*lKsi[0];
        gradVelComp[i][1] = gf[i][1] + du[i]*lKsi[1];
        gradVelComp[i][2] = gf[i][2] + du[i]*lKsi[2];
/*...................................................................*/

/*... derivadas direcionais*/
/*...*/
        gfKsi[i] = gradVelComp[i][0]*t[0] 
                 + gradVelComp[i][1]*t[1] 
                 + gradVelComp[i][2]*t[2];
/*...................................................................*/

/*... correcao nao-ortogonal*/
        dfdc[i] = coef*gfKsi[i];
      }
/*...................................................................*/

/*... fluxo convectivo upwind de primeira ordem*/
      cv   = density*wfn*lFarea;
/*...................................................................*/

/*... correcao do fluxo advectivo*/
      v[0] = lXm[0] - ccV[0];
      v[1] = lXm[1] - ccV[1];
      v[2] = lXm[2] - ccV[2];
      advectiveScheme(velC          ,velV
                     ,gradVelC[0]   ,gradVelV[0]
                     ,gradVelComp[0],lvSkew
                     ,lXmcc         ,v  
                     ,lKsi          ,lModKsi
                     ,cv            ,cvc
                     ,alphaMenosUm  ,alpha
                     ,pAdv          ,ndm
                     ,iCodAdv1      ,iCodAdv2);  
/*...................................................................*/

/*...*/
      lA[nf] = dfd - min(cv,0.0e0);
      sP         += cv;
/*... correcao nao ortogonal e do fluxo advectivo*/
      p[0]       += dfdc[0] - cv*cvc[0];
      p[1]       += dfdc[1] - cv*cvc[1];
      p[2]       += dfdc[2] - cv*cvc[2];
/*...................................................................*/

/*... gradiente da pressao com resconstrucao de segunda ordem
      (forma conservativa)*/
      facePressure(gradPresC       ,gradPresV
                 ,lXmcc            ,lXm 
                 ,ccV              ,pf
                 ,lNormal          ,NULL
                 ,0.e0             ,0.e0
                 ,lFarea
                 ,presC            ,presV
                 ,alphaMenosUm     ,alpha    
                 ,false            ,typeFacePres);
/*...................................................................*/

/*... termos viscosos explicitos*/      
      if(fViscosity)
        viscosityPartExp(p          ,gradVelComp[0]
                        ,lNormal
                        ,viscosityC ,lFarea
                        ,false);
/*...................................................................*/

/*... termos viscosos explicitos*/
      if(fStruc)
        turbStructIntegral(stressRc    ,stressRv
                          ,lNormal     ,p
                          ,lFarea
                          ,alphaMenosUm,alpha    
                          ,false);
/*...................................................................*/

    }
/*... contorno*/
    else{
      lA[nf] = 0.0e0;
      wfn = velC[0]*lNormal[0] 
          + velC[1]*lNormal[1] 
          + velC[2]*lNormal[2];

/*...*/
      s[0] = dFieldC[0]*lFarea*lNormal[0];
      s[1] = dFieldC[1]*lFarea*lNormal[1];
      s[2] = dFieldC[2]*lFarea*lNormal[2];
/*...*/
      pFace = presC;
/*... pressao prescrita*/
      if(lFacePresR[nf])
      {
/*...cargas*/
        nCarg = lFacePresR[nf]-1;
        pLoadSimplePres( &sP             , p
                       , tA              , lXmcc
                       , presC           , gradPresC   
                       , s               , e         
                       , t               , lNormal
                       , NULL            , NULL
                       , velC            , 0.e0   
                       , densityC        , densityC
                       , lFarea          , dcca[nf]
                       , &lPres[nCarg]   , ndm
                       , 0               , false     
                       , false); 
      } 
/*...................................................................*/

/*... gradiente da pressao com resconstrucao de segunda ordem*/
      facePressure(gradPresC       ,gradPresV
                 ,lXmcc            ,lXm
                 ,ccV              ,pf
                 ,lNormal          ,NULL
                 ,pFace            ,0.e0    
                 ,lFarea
                 ,presC            ,0.e0
                 ,0.0              ,0.e0     
                 ,true             ,typeFacePres);
/*...................................................................*/

/*... termos viscosos explicitos*/
      if(fStruc)
        turbStructIntegral(stressRc    ,stressRv
                          ,lNormal     ,p
                          ,lFarea
                          ,0.e0        ,0.e0     
                          ,true);
/*...................................................................*/

/*...*/
      s[0] = lFarea*lNormal[0];
      s[1] = lFarea*lNormal[1];
      s[2] = lFarea*lNormal[2];
/*...*/
      difusionSchemeNew(s   , lKsi
                       , e  , t
                       , ndm, iCodDif);
/*...................................................................*/

/*... velocidades*/
      if(lFaceVelR[nf])
      {
        xx[0] = MAT2D(nf, 0, xm, 3);
        xx[1] = MAT2D(nf, 1, xm, 3);
        xx[2] = MAT2D(nf, 2, xm, 3);
        xx[3] = ts;
/*...cargas*/
        nCarg = lFaceVelR[nf]-1;
        pLoadSimple(sPc         , p
                  , tA          , lXmcc
                  , velC        , gradVelC[0] 
                  , presC       , gradPresC
                  , viscosityC  , effViscosityC
                  , xx  
                  , s           , e         
                  , t           , lNormal
                  , densityC    , wallPar
                  , lFarea      , dcca[nf]
                  , &lVel[nCarg], ndm
                  , true        , false
                  , false       , fWallModel 
                  , wallType     , nel );      
      }  
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... distretização temporal*/
  if(fTime)
  {
/*... EULER*/
    if(typeTime == EULER) 
      sP     += densityC*volume[idCell]/dt;
/*...BACKWARD*/
    else if(typeTime == BACKWARD)
    {
      tmp = 1.e0/dt + 1.e0/(dt+dt0);
      sP += tmp*densityC*volume[idCell];
    }
  }
/*...................................................................*/

/*...*/
  lAn = 0.e0;
  for (nf = 0; nf<nFace; nf++)
    lAn += lA[nf];
  lA[idCell  ] = sP + sPc[0] + lAn;
  lA[idCell+1] = sP + sPc[1] + lAn;
  lA[idCell+2] = sP + sPc[2] + lAn;
/*...................................................................*/

/*... under-relaxation(simple)*/
  lA[idCell]   = lA[idCell  ]/underU;
  lA[idCell+1] = lA[idCell+1]/underU;
  lA[idCell+2] = lA[idCell+2]/underU;
  p[0]        += (1-underU)*lA[idCell  ]*velC[0];
  p[1]        += (1-underU)*lA[idCell+1]*velC[1];
  p[2]        += (1-underU)*lA[idCell+2]*velC[2];
/*...................................................................*/

/*... GradP*/
  p[0] -= pf[0];
  p[1] -= pf[1];
  p[2] -= pf[2];
/*...................................................................*/

/*...*/
  rCell[0] = 0.0e0;
  rCell[1] = 0.0e0;
  rCell[2] = 0.0e0;
  for(nf=0;nf<nFace;nf++){
    if( lViz[nf] > -1){
/*... pasando os valores conhecidos para o lado direito*/
      if(lId[nf] == -2){
        p[0] += lA[nf]*MAT2D(nf,0,vel,ndm);
        p[1] += lA[nf]*MAT2D(nf,1,vel,ndm);
        p[2] += lA[nf]*MAT2D(nf,2,vel,ndm);
      }
      else{
/*residuo (R = F-KvizUviz ) e valores prescritos por elemento*/
        rCell[0] += lA[nf]*MAT2D(nf,0,vel,ndm);
        rCell[1] += lA[nf]*MAT2D(nf,1,vel,ndm);
        rCell[2] += lA[nf]*MAT2D(nf,2,vel,ndm);
      } 
    }
  }
/*... residuo: R = F - KpUp*/ 
   rCell[0] += p[0] -lA[idCell  ]*velC[0];
   rCell[1] += p[1] -lA[idCell+1]*velC[1];
   rCell[2] += p[2] -lA[idCell+2]*velC[2];
/*...................................................................*/

/*...*/
   if(typeSimple == SIMPLEC){ 
    dField[0] = volume[idCell]/(lA[idCell  ] - lAn);
    dField[1] = volume[idCell]/(lA[idCell+1] - lAn);
    dField[2] = volume[idCell]/(lA[idCell+2] - lAn);
  }
  else if(typeSimple == SIMPLE){ 
    dField[0] = volume[idCell]/lA[idCell];
    dField[1] = volume[idCell]/lA[idCell+1];
    dField[2] = volume[idCell]/lA[idCell+2];
  } 
/*...................................................................*/

/*...*/
  for (nf = 0; nf<nFace; nf++)
    lA[nf] *= -1.e0;
/*...................................................................*/

/*...*/
  if(fRes)
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

/*********************************************************************
 * Data de criacao    : 03/10/2017                                   *
 * Data de modificaco : 12/10/2019                                   * 
 *-------------------------------------------------------------------* 
 * CELLSIMPLEVE3DLM: Celula 3D para velocidade do metodo simple      * 
 * em escoamento levemento compressivel (Low Mach)                   * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loadsVel  -> definicoes de cargas de velocidades                  *
 * loadsPres -> definicoes de cargas de pressao                      *
 * advVel    -> tecnica da discretizacao do termo advecao            *
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * tModel    -> modelo de turbulencia                                *
 * momentumModel -> termos/modelos da equacao de momento linear      *
 * typeSimple-> tipo do metodo simple                                *
 * lnFace    -> numero de faces da celula central e seus vizinhos    *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lViz      -> viznhos da celula central                            *
 * lId       -> numeracoes das equacoes das celulas                  *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *              vizinhos desta                                       *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * fArea     -> area da face                                         *
 * normal    -> vetores normais as faces das celulas                 *
 * volume    -> volume da celula central                             *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *              faces da celula central                              *
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os *
 *              centrois compartilhado nessa face da celula central  *
 * dcca      -> menor distancia do centroide central a faces desta   *
 *              celula                                               *
 * cc        -> centroides da celula centra e seus vizinhos          *                                             *
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceVelR  -> restricoes por elemento de velocidades               *
 * facePresR -> restricoes por elemento de pressao                   *
 * pres      -> campo de pressao conhecido                           *
 * gradPes   -> gradiente reconstruido da pressao                    *
 * vel       -> campo de velocidade conhecido                        *
 * gradVel   -> gradiente rescontruido das velocidades               *
 * gradRho   -> gradiente rescontruido da densidade                  *
 * lDensity  -> massa especifica sem variacao temporal               *
 * lDviscosity-> viscosidade dinamica com variacao temporal          *
 * dField    -> matriz D do metodo simple                            * 
 * stressR   -> tensao residual estrutural                           *
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  *
 * densityMed-> densidade media do meio                              *
 * underU    -> parametro de sobre relaxamento                       *
 * nEn       -> numero de nos da celula central                      *
 * sPressure -> reconstrucao de segunda ordem para pressoes nas      *
 *              faces                                                *
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
void cellSimpleVel3DLm(Loads *lVel        , Loads *lPres 
            , Advection *advVel           , Diffusion *diffVel
            , Turbulence *tModel          , MomentumModel *momentumModel            
            , short const typeSimple 
            , short *RESTRICT lGeomType   
            , INT *RESTRICT lViz          , INT *RESTRICT lId  
            , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
            , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT fArea
            , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume 
            , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
            , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
            , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
            , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
            , DOUBLE *RESTRICT lRcell     , Temporal *ddt
            , short  *RESTRICT lFaceVelR  , short *RESTRICT lFacePresR
            , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres 
            , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
            , DOUBLE *RESTRICT gradRho
            , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lViscosity 
            , DOUBLE *RESTRICT dField     , DOUBLE *RESTRICT stressR
            , DOUBLE *RESTRICT wallPar    , DOUBLE const densityMed
            , DOUBLE const underU         , const bool sPressure
            , const short nEn             , short const nFace    
            , const short ndm             , INT const nel)
{ 
/*...*/
  bool fTime, fRes, fTurb, fRhieInt, fWallModel, fStruc;
  bool fDiv,fViscosity,fSoPressure,fFacePres;
/*...*/
  short iCodAdv1, iCodAdv2, iCodDif, wallType, idCell, nf
      , nCarg, typeTime, iCodBuoyant, iCodPolFace, i,typeFacePres;
/*...*/
  INT vizNel;
/*...*/
  DOUBLE viscosityC, viscosityV, viscosity, 
         eddyViscosityC, eddyViscosityV, effViscosityC, effViscosityV,
         densityC, densityV, density, densityBd;
/*...*/
  DOUBLE rCell[3], dt, dt0, p[3], sP, sPc[3], tA[3]
        ,dFieldC[3], dFieldV[3], dFieldF[3];
/*...*/
  DOUBLE v[3], lKsi[3], lNormal[3], lXmcc[3], ccV[3];
  DOUBLE lModKsi, lFarea, du[3], duDksi[3], lXm[3];
  DOUBLE coef, lAn, tmp;
/*...*/
  DOUBLE dfd, cv, cvc[3], lvSkew[3];
/*... nonOrtogonal*/
  DOUBLE e[3], t[3], modE, dfdc[3];
/*... interpolacao linear*/
  DOUBLE alpha, alphaMenosUm;
/*... */
  DOUBLE pFace, pf[3], presC, presV, gradPresC[3], gradPresV[3];
/*... */
  DOUBLE gf[3][3], gfKsi[3];
  DOUBLE wfn, velC[3], velV[3], wf[3];
  DOUBLE gradVelC[3][3], gradVelV[3][3], gradVelComp[3][3];
/*...*/
  DOUBLE stressRc[6],stressRv[6],s[6];
/*...*/
  DOUBLE xx[4],ts;
/*... */
  DOUBLE g[3],gh;
/*...*/
  DOUBLE pAdv[NPADV];
  
/*...*/
  idCell      = nFace;
  iCodAdv1    = advVel->iCod1;
  iCodAdv2    = advVel->iCod2;
  pAdv[0]     = advVel->par[0];
  iCodDif     = diffVel->iCod;
  iCodPolFace = INTPOLFACELINEAR;
  fFacePres   = true;
/*...................................................................*/

/*...*/
  ts       = ddt->t; 
  dt       = ddt->dt[0];  
  dt0      = ddt->dt[1];
  typeTime = ddt->type;
  fTime    = ddt->flag;
/*...................................................................*/

/*...*/  
  fWallModel       = tModel->fWall;
  wallType         = tModel->wallType;
  fTurb            = tModel->fTurb;
  fStruc           = tModel->fTurbStruct;
/*...*/
  fSoPressure      = momentumModel->fSoPressure;
  fRhieInt         = momentumModel->fRhieChowInt;
  fRes             = momentumModel->fRes;
  fViscosity       = momentumModel->fViscosity;
  fDiv             = momentumModel->fDiv;
  iCodBuoyant      = momentumModel->iCodBuoyant;
//typeFacePres     = INT_BUOYANT_FORCE;
  if(fSoPressure)
    typeFacePres   = INT_PRESSURE_SO;
  else
    typeFacePres   = INT_PRESSURE_FO;  
/*...................................................................*/

/*...*/
  ts       = ddt->t;
  dt       = ddt->dt[0];
  dt0      = ddt->dt[1];
  typeTime = ddt->type;
  fTime    = ddt->flag;
/*...................................................................*/

/*...*/
  densityC      = lDensity[idCell];
  if(iCodBuoyant == BUOYANT_PRGH)
    densityBd = densityC;
  else if (iCodBuoyant == BUOYANT_RHOREF)
    densityBd = densityMed;
/*...................................................................*/

/*...*/
  eddyViscosityC = eddyViscosityV = 0.0e0;
  viscosityC     = MAT2D(idCell, 0, lViscosity, 2);
  if(fTurb) eddyViscosityC= MAT2D(idCell, 1, lViscosity, 2);
  effViscosityC = viscosityC + eddyViscosityC;
/*...................................................................*/

/*...*/
  sP = gh = 0.e0;
  presC = pres[idCell];
  for (i = 0; i < 3; i++)
  {
/*... | dui/dx1 dui/dx2 dui/dx3*/
    gradVelC[i][0] = MAT3D(idCell, i, 0, gradVel, 3, 3);
    gradVelC[i][1] = MAT3D(idCell, i, 1, gradVel, 3, 3);
    gradVelC[i][2] = MAT3D(idCell, i, 2, gradVel, 3, 3);
/*...*/
    gradPresC[i] = MAT2D(idCell, i, gradPres, 3);
    dFieldC[i] = MAT2D(idCell, i, dField, 3);
/*...*/
    velC[i] = MAT2D(idCell, i, vel, 3);
/*...................................................................*/

/*...*/
    p[i] = pf[i] = sPc[i] =  0.e0;
/*...................................................................*/

/*...*/
    g[i]  = gravity[i];
    gh   += (MAT2D(idCell,i,cc  ,3) - xRef[i])*g[i];
/*...................................................................*/

  }
/*...................................................................*/

/*... tensor residual*/
  if (fStruc)
    for (i = 0; i < 6; i++)
      stressRc[i] = MAT2D(idCell, i, stressR, 6);
/*...................................................................*/

/*...*/
  for(nf=0;nf<nFace;nf++)
  {
    vizNel = lViz[nf];
    lFarea = fArea[nf];
    lModKsi = mKsi[nf];
    /*..*/
    for (i = 0; i < 3; i++)
    {
      lNormal[i] = MAT2D(nf, i, normal, ndm);
      lXmcc[i] = MAT2D(nf, i, xmcc, ndm);
      lKsi[i] = MAT2D(nf, i, ksi, 3);
    }
/*...................................................................*/

/*... dominio*/
    if( vizNel  > -1 )
    {
/*...*/
      densityV    = lDensity[nf];
      viscosityV  = MAT2D(nf, 0, lViscosity, 2);
      if(fTurb) eddyViscosityV = MAT2D(nf, 1, lViscosity, 2);
/*...*/
      presV = pres[nf];
      for (i = 0; i<3; i++)
      {
        velV[i] = MAT2D(nf, i, vel, ndm);
/*...*/
        gradPresV[i] = MAT2D(nf, i, gradPres, 3);
        dFieldV[i] = MAT2D(nf, i, dField, 3);
/*..*/
        lvSkew[i] = MAT2D(nf, i, vSkew, 3);
/*..*/
        duDksi[i] = (velV[i] - velC[i]) / lModKsi;
/*... | dui/dx1 dui/dx2 dui/dx3*/
        gradVelV[i][0] = MAT3D(nf, i, 0, gradVel, 3, 3);
        gradVelV[i][1] = MAT3D(nf, i, 1, gradVel, 3, 3);
        gradVelV[i][2] = MAT3D(nf, i, 2, gradVel, 3, 3);
/*...*/
        ccV[i] = MAT2D(nf, i, cc, 3);
        lXm[i] = MAT2D(nf, i, xm, 3);
/*...................................................................*/

      }
/*...................................................................*/

/*... tensor residual*/
      if (fStruc)
        for (i = 0; i<6; i++)
          stressRv[i] = MAT2D(nf, i, stressR, 6);
/*...................................................................*/

/*... termo difusivo
grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
/*...*/
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
      effViscosityV = viscosityV + eddyViscosityV;
      viscosity = mediaHarmonica(effViscosityV,effViscosityC
                                ,alphaMenosUm ,alpha        );
/*...................................................................*/

/*... difusao direta*/
      coef = viscosity;
      modE = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
      dfd  = coef*modE/lModKsi;
/*...................................................................*/
      
/*...*/
      for (i = 0; i<3; i++)
      {
        gf[i][0] = alphaMenosUm * gradVelC[i][0] + alpha * gradVelV[i][0];
        gf[i][1] = alphaMenosUm * gradVelC[i][1] + alpha * gradVelV[i][1];
        gf[i][2] = alphaMenosUm * gradVelC[i][2] + alpha * gradVelV[i][2];
/*...*/
        wf[i] = alphaMenosUm * velC[i] + alpha * velV[i];
      }
/*...................................................................*/
      density  = alphaMenosUm*densityC  + alpha*densityV;
/*...................................................................*/

/*... velocidade normal a face*/
      if(fRhieInt)
      {
        dFieldF[0] = alphaMenosUm*dFieldC[0] + alpha*dFieldV[0];
        dFieldF[1] = alphaMenosUm*dFieldC[1] + alpha*dFieldV[1];
        dFieldF[2] = alphaMenosUm*dFieldC[2] + alpha*dFieldV[2];
        wfn = interpolFaceVel(velC         ,velV
                             ,presC        ,presV
                             ,gradPresC    ,gradPresV
                             ,lNormal      ,lKsi
                             ,lModKsi      ,dFieldF
                             ,alphaMenosUm ,alpha
                             ,ndm);
      }
      else 
        wfn   = wf[0]*lNormal[0] + wf[1]*lNormal[1] + wf[2]*lNormal[2]; 
/*...................................................................*/

/*...*/
      for (i = 0; i<3; i++)
      {
/*... derivadas direcionais*/
        gfKsi[i] = gf[i][0] * lKsi[0]
          + gf[i][1] * lKsi[1]
          + gf[i][2] * lKsi[2];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
        du[i] = duDksi[i] - gfKsi[i];
/*...*/
        gradVelComp[i][0] = gf[i][0] + du[i] * lKsi[0];
        gradVelComp[i][1] = gf[i][1] + du[i] * lKsi[1];
        gradVelComp[i][2] = gf[i][2] + du[i] * lKsi[2];
/*...................................................................*/

/*... derivadas direcionais*/
/*...*/
        gfKsi[i] = gradVelComp[i][0] * t[0]
          + gradVelComp[i][1] * t[1]
          + gradVelComp[i][2] * t[2];
/*...................................................................*/

/*... correcao nao-ortogonal*/
        dfdc[i] = coef * gfKsi[i];
      }
/*...................................................................*/

/*... fluxo convectivo upwind de primeira ordem*/
      cv   = density*wfn*lFarea;
/*...................................................................*/

/*... correcao do fluxo advectivo*/
      v[0] = lXm[0] - ccV[0];
      v[1] = lXm[1] - ccV[1];
      v[2] = lXm[2] - ccV[2];
      advectiveScheme(velC          ,velV
                     ,gradVelC[0]   ,gradVelV[0]
                     ,gradVelComp[0],lvSkew
                     ,lXmcc         ,v  
                     ,lKsi          ,lModKsi
                     ,cv            ,cvc
                     ,alphaMenosUm  ,alpha
                     ,pAdv          ,ndm
                     ,iCodAdv1      ,iCodAdv2);
/*...................................................................*/

/*...*/
      lA[nf] = dfd - min(cv,0.0e0);
      sP         += cv;
/*... correcao nao ortogonal e do fluxo advectivo*/
      p[0]       += dfdc[0] - cv*cvc[0];
      p[1]       += dfdc[1] - cv*cvc[1];
      p[2]       += dfdc[2] - cv*cvc[2];  
/*...................................................................*/

/*... gradiente da pressao com resconstrucao de segunda ordem
(forma conservativa)*/
      facePressure(gradPresC   , gradPresV
                 , lXmcc       , lXm
                 , ccV         , pf
                 , lNormal     , g
                 , 0.e0        , densityBd
                 , lFarea
                 , presC       , presV
                 , alphaMenosUm, alpha
                 , false       , typeFacePres);
/*...................................................................*/

/*... termos viscosos explicitos*/
      if (fViscosity)
        viscosityPartExp(p         , gradVelComp[0]
                       , lNormal
                       , viscosityC, lFarea
                       , fDiv);
/*...................................................................*/

/*... termos viscosos explicitos*/
      if (fStruc)
        turbStructIntegral(stressRc, stressRv
          , lNormal, p
          , lFarea
          , alphaMenosUm, alpha
          , false);
/*...................................................................*/

    }
/*... contorno*/
    else{
      lA[nf] = 0.0e0;
      wfn = velC[0]*lNormal[0] 
          + velC[1]*lNormal[1] 
          + velC[2]*lNormal[2];
/*...*/
      s[0] = dFieldC[0]*lFarea*lNormal[0];
      s[1] = dFieldC[1]*lFarea*lNormal[1];
      s[2] = dFieldC[2]*lFarea*lNormal[2];
/*...*/
//    difusionSchemeAnisotropic(s, lKsi, e, t, ndm, iCodDif);
/*...................................................................*/

/*...*/
      pFace = presC;
/*... pressao prescrita*/
      if(lFacePresR[nf])
      {
/*...cargas*/
        nCarg = lFacePresR[nf]-1;
        pLoadSimplePres(&sP             , p
                      , &pFace          , lXmcc
                      , presC           , gradPresC   
                      , s               , e         
                      , t               , lNormal
                      , g               , gradRho
                      , velC            , gh    
                      , densityC        , densityMed  
                      , lFarea          , dcca[nf]
                      , &lPres[nCarg]   , ndm
                      , iCodBuoyant     , true      
                      , false); 
      } 
/*...................................................................*/

/*... gradiente da pressao com resconstrucao de segunda ordem*/
      facePressure(gradPresC  , gradPresV
                 , lXmcc      , lXm
                 , ccV        , pf
                 , lNormal    , g
                 , pFace      , densityBd
                 , lFarea
                 , presC      , 0.e0
                 , 0.0        , 0.e0
                 , true       , typeFacePres);
/*...................................................................*/

/*... termos viscosos explicitos*/
      if (fStruc)
        turbStructIntegral(stressRc, stressRv
                        , lNormal, p
                        , lFarea
                        , 0.e0, 0.e0
                        , true);
/*...................................................................*/

/*...*/
      s[0] = lFarea*lNormal[0];
      s[1] = lFarea*lNormal[1];
      s[2] = lFarea*lNormal[2];
/*...*/
      difusionSchemeNew(s   , lKsi
                       , e  , t
                       , ndm, iCodDif);
/*...................................................................*/

/*... velocidades*/
      if(lFaceVelR[nf])
      {
        xx[0] = MAT2D(nf, 0, xm, 3);
        xx[1] = MAT2D(nf, 1, xm, 3);
        xx[2] = MAT2D(nf, 2, xm, 3);
        xx[3] = ts;
/*...cargas*/
        nCarg = lFaceVelR[nf]-1;
        pLoadSimple(sPc         , p
                  , tA          , lXmcc
                  , velC        , gradVelC[0] 
                  , presC       , gradPresC
                  , viscosityC  , effViscosityC
                  , xx  
                  , s           , e         
                  , t           , lNormal
                  , densityC    , wallPar
                  , lFarea      , dcca[nf]
                  , &lVel[nCarg], ndm
                  , true        , false
                  , fDiv        , fWallModel 
                  , wallType    , nel );    
      }  
/*...................................................................*/

    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if(iCodBuoyant == BUOYANT_HYDROSTATIC)
  {
    tmp   = densityC*volume[idCell];
    p[0] += tmp*g[0];
    p[1] += tmp*g[1];
    p[2] += tmp*g[2];
  }
  else if (iCodBuoyant == BUOYANT_PRGH) 
  {
    tmp   = gh*volume[idCell];
    p[0] -= tmp*gradRho[0];
    p[1] -= tmp*gradRho[1];
    p[2] -= tmp*gradRho[2]; 
  }
  else if (iCodBuoyant == BUOYANT_RHOREF)
  {
    tmp   = (densityC-densityMed)*volume[idCell];
    p[0] += tmp*g[0];
    p[1] += tmp*g[1];
    p[2] += tmp*g[2];
  }
/*...................................................................*/
  
/*... distretização temporal*/
  if(fTime)
  {
/*... EULER*/
    if(typeTime == EULER) 
      sP     += densityC*volume[idCell]/dt;
/*...BACKWARD*/
    else if(typeTime == BACKWARD)
    {
      tmp = 1.e0/dt + 1.e0/(dt+dt0);
      sP += tmp*densityC*volume[idCell];
    }
  }
/*...................................................................*/

/*...*/
  lAn = 0.e0;
  for (nf = 0; nf<nFace; nf++)
    lAn += lA[nf];
  lA[idCell  ] = sP + sPc[0] + lAn;
  lA[idCell+1] = sP + sPc[1] + lAn;
  lA[idCell+2] = sP + sPc[2] + lAn;
/*...................................................................*/

/*... under-relaxation(simple)*/
  lA[idCell]   = lA[idCell  ]/underU;
  lA[idCell+1] = lA[idCell+1]/underU;
  lA[idCell+2] = lA[idCell+2]/underU;
  p[0]        += (1-underU)*lA[idCell  ]*velC[0];
  p[1]        += (1-underU)*lA[idCell+1]*velC[1];
  p[2]        += (1-underU)*lA[idCell+2]*velC[2];
/*...................................................................*/

/*... GradP*/
  if(fFacePres)
  {
    p[0] -= pf[0];
    p[1] -= pf[1];
    p[2] -= pf[2];  
  }
  else
  {
    p[0] -= gradPresC[0]*volume[idCell];
    p[1] -= gradPresC[1]*volume[idCell];
    p[2] -= gradPresC[2]*volume[idCell];
  }
/*...................................................................*/

/*...*/
  rCell[0] = 0.0e0;
  rCell[1] = 0.0e0;
  rCell[2] = 0.0e0;
  for(nf=0;nf<nFace;nf++){
    if( lViz[nf] > -1){
/*... pasando os valores conhecidos para o lado direito*/
      if(lId[nf] == -2){
        p[0] += lA[nf]*MAT2D(nf,0,vel,ndm);
        p[1] += lA[nf]*MAT2D(nf,1,vel,ndm);
        p[2] += lA[nf]*MAT2D(nf,2,vel,ndm);
      }
      else{
/*residuo (R = F-KvizUviz ) e valores prescritos por elemento*/
        rCell[0] += lA[nf]*MAT2D(nf,0,vel,ndm);
        rCell[1] += lA[nf]*MAT2D(nf,1,vel,ndm);
        rCell[2] += lA[nf]*MAT2D(nf,2,vel,ndm);
      } 
    }
  }
/*... residuo: R = F - KpUp*/
   rCell[0] += p[0] -lA[idCell  ]*velC[0];
   rCell[1] += p[1] -lA[idCell+1]*velC[1];
   rCell[2] += p[2] -lA[idCell+2]*velC[2];
/*...................................................................*/

/*...*/
   if(typeSimple == SIMPLEC){ 
     MAT2D(idCell,0,dField,3) = volume[idCell] / (lA[idCell]     - lAn);
     MAT2D(idCell,1,dField,3) = volume[idCell] / (lA[idCell + 1] - lAn);
     MAT2D(idCell,2,dField,3) = volume[idCell] / (lA[idCell + 2] - lAn);   
  }
  else if(typeSimple == SIMPLE){ 
    MAT2D(idCell,0,dField,3) = volume[idCell] / lA[idCell];
    MAT2D(idCell,1,dField,3) = volume[idCell] / lA[idCell + 1];
    MAT2D(idCell,2,dField,3) = volume[idCell] / lA[idCell + 2];  
  } 
/*...................................................................*/

/*...*/
  for (nf = 0; nf<nFace; nf++)
    lA[nf] *= -1.e0;
/*...................................................................*/

/*...*/
  if(fRes){
    lB[0] = rCell[0];
    lB[1] = rCell[1];
    lB[2] = rCell[2];
  }
  else{
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