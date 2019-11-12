#include<CellLib.h>
/*********************************************************************
 * Biblioteca de celulas                                             *
 *-------------------------------------------------------------------*
 * Celulas 2D                                                        *
 *-------------------------------------------------------------------*
 * ------------------ INCOMPRESSIVEl ------------------------------- *
 *                                                                   *
 * cellSimplePres2D: Celula 2D para equacao de correcao de pressoa   *
 * metodo simple em escoamento imcompressivel                        *
 *                                                                   *
 * cellSimpleNonOrthPres2D: correcao nao ortogonal da Celula 2D      *
 * para equacao de correcao de pressao                               *
 *                                                                   *
 * ------------------ LEVEMENTE COMPRESSIVEl ------------------------*
 *                                                                   *
 * cellSimplePres2DLm: Celula 2D para equacao de correcao de pressoa *
 * metodo simple em escoamento compressivel para baixo mach          *
 *                                                                   *
 *-------------------------------------------------------------------*
 * Celulas 3D                                                        *
 *-------------------------------------------------------------------*
 * ------------------ INCOMPRESSIVEl ------------------------------- *
 *                                                                   *
 * cellSimplePres3D: Celula 3D para equacao de correcao de pressoa   *
 * metodo simple em escoamento imcompressivel                        *
 *                                                                   *
 * cellSimpleNonOrthPres3D: correcao nao ortogonal da Celula 2D      *
 * para equacao de correcao de pressao                               *
 *                                                                   *
 * ------------------ LEVEMENTE COMPRESSIVEl ------------------------*
 *                                                                   *
 * cellSimplePres3DLm: Celula 3D para equacao de correcao de pressoa *
 * metodo simple em escoamento compressivel para baixo mach          *
 *                                                                   *
 *-------------------------------------------------------------------*
*********************************************************************/

/*********************************************************************
 * Data de criacao    : 01/07/2016                                   *
 * Data de modificaco : 19/07/2018                                   * 
 *-------------------------------------------------------------------* 
 * CELLSIMPLEPRES2D: Celula 2D para equacao de correcao de pressoa   *
 * metodo simple em escoamento imcompressivel                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  * 
 * loadsPres -> definicoes de cargas de pressao                      * 
 * diffPres  -> tecnica da discretizacao do termo difusivo           *
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
 * lDensity  -> massa especifica sem variacao temporal               * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * faceVelR  -> restricoes por elemento de velocidades               * 
 * faceVelL  -> carga por elemento de velocidades                    * 
 * facePresR -> restricoes por elemento de pressao                   * 
 * facePresL -> carga por elemento de pressao                        * 
 * pres      -> campo de pressao conhecido                           * 
 * gradPes   -> gradiente reconstruido da pressao                    * 
 * vel       -> campo de velocidade conhecido                        * 
 * gradVel   -> gradiente rescontruido das velocidades               * 
 '   
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
void cellSimplePres2D(Loads *loadsVel     ,Loads *loadsPres 
							,Diffusion *diffPres	
              ,short *RESTRICT lGeomType,DOUBLE *RESTRICT prop
              ,INT *RESTRICT lViz       ,INT *RESTRICT lId  
              ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
              ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT mEta
              ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT area   
              ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
              ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
              ,DOUBLE *RESTRICT lRcell  
              ,short  *RESTRICT lFaceVelR  ,short *RESTRICT lFaceVelL
              ,short  *RESTRICT lFacePresR ,short *RESTRICT lFacePresL
              ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres 
              ,DOUBLE *RESTRICT vel     ,DOUBLE *RESTRICT dField  
              ,const short nEn          ,short const nFace    
              ,const short ndm          ,INT const nel)
{ 

  DOUBLE densityC,densityV ,density;
  DOUBLE dFieldC[2],dFieldV[2],dFieldF[2];
/*...*/
  DOUBLE rCell;
  DOUBLE p,sP;
/*...*/
  DOUBLE v[2],lKsi[2],lNormal[2];
  DOUBLE dPviz,lModKsi,lModEta;
/*...*/
  DOUBLE gradPresC[2],gradPresV[2];
/*...*/
  DOUBLE dfd,coef,lvSkew[2];
/*... nonOrtogonal*/
	DOUBLE e[2],t[2],s[2],modE;
/*... interpolacao linear*/
  DOUBLE alpha,alphaMenosUm;
/*... */
  DOUBLE wfn,velC[2],velF[2],presC,presV;
/*...*/
	short iCodDif = diffPres->iCod;
/*...*/
  short idCell = nFace;
  short nAresta;
  INT vizNel;

/*...*/
  densityC = lDensity[idCell];
/*...................................................................*/

/*...*/
  velC[0]      = MAT2D(idCell,0,vel,ndm);
  velC[1]      = MAT2D(idCell,1,vel,ndm);
  
  gradPresC[0] = MAT2D(idCell,0,gradPres,ndm);
  gradPresC[1] = MAT2D(idCell,1,gradPres,ndm);

  presC     = pres[idCell];

  dFieldC[0] = MAT2D(idCell,0,dField,2); 
  dFieldC[1] = MAT2D(idCell,1,dField,2);
/*...................................................................*/

  p          = 0.0e0;
  sP         = 0.0e0;
  rCell      = 0.0e0;
  for(nAresta=0;nAresta<nFace;nAresta++){
    vizNel     = lViz[nAresta];
    lNormal[0] = MAT2D(nAresta,0,normal,ndm);
    lNormal[1] = MAT2D(nAresta,1,normal,ndm);
    lModEta    = mEta[nAresta];
    lKsi[0] = MAT2D(nAresta, 0, ksi, ndm);
    lKsi[1] = MAT2D(nAresta, 1, ksi, ndm);
    lModKsi = mKsi[nAresta];
/*... dominio*/
    if( vizNel  > -1 ){
/*...*/
      velF[0]       = MAT2D(nAresta,0,vel,ndm);
      velF[1]       = MAT2D(nAresta,1,vel,ndm);
      lvSkew[0]     = MAT2D(nAresta,0,vSkew,ndm);
      lvSkew[1]     = MAT2D(nAresta,1,vSkew,ndm);
      presV         = pres[nAresta];
      densityV      = lDensity[nAresta];
      dFieldV[0]    = MAT2D(nAresta,0,dField,2);
      dFieldV[1]    = MAT2D(nAresta,1,dField,2);
      gradPresV[0]  = MAT2D(nAresta,0,gradPres,ndm);
      gradPresV[1]  = MAT2D(nAresta,1,gradPres,ndm);
/*...................................................................*/

/*...*/
      v[0]         = lvSkew[0] + MAT2D(nAresta,0,xmcc,ndm);
      v[1]         = lvSkew[1] + MAT2D(nAresta,1,xmcc,ndm);
      dPviz        = sqrt(v[0]*v[0] + v[1]*v[1]);
      alpha        = dPviz/lModKsi;
      alphaMenosUm = 1.0e0 - alpha; 
/*...................................................................*/
      
/*... interpolacao das propriedades*/
      density    = alphaMenosUm*densityC   + alpha*densityV;
      dFieldF[0] = alphaMenosUm*dFieldC[0] + alpha*dFieldV[0];
      dFieldF[1] = alphaMenosUm*dFieldC[1] + alpha*dFieldV[1];
/*...................................................................*/

/*... termo difusivo
grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      s[0] = dFieldF[0]*lModEta*lNormal[0];
      s[1] = dFieldF[1]*lModEta*lNormal[1];
      difusionSchemeAnisotropic(s,lKsi,e,t,ndm,iCodDif);
/*...................................................................*/

/*... difusao direta*/
			coef = density;
			modE = sqrt(e[0]*e[0]+e[1]*e[1]);
  		dfd  = coef*modE/lModKsi;
/*...................................................................*/

/*...*/
      lA[nAresta] = dfd;
/*...................................................................*/

/*...*/
      wfn = interpolFaceVel(velC         ,velF
                           ,presC        ,presV
                           ,gradPresC    ,gradPresV
                           ,lNormal      ,lKsi
                           ,lModKsi      ,dFieldF
                           ,alphaMenosUm ,alpha
                           ,ndm);
/*...................................................................*/
      p   -= density*wfn*lModEta;
/*...................................................................*/
    }
/*... contorno*/
    else{
      lA[nAresta] = 0.0e0;
      wfn = velC[0]*lNormal[0] + velC[1]*lNormal[1];
      if(lFacePresR[nAresta]){
/*... termo difusivo
        grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
        s[0] = dFieldC[0]*lModEta*lNormal[0];
        s[1] = dFieldC[1]*lModEta*lNormal[1];
/*...*/
        difusionSchemeAnisotropic(s,lKsi,e,t,ndm,iCodDif);
/*...................................................................*/

/*...*/
        modE = sqrt(e[0]*e[0] + e[1]*e[1]);
/*...................................................................*/

/*...cargas*/
//      nCarg = lFacePresL[nAresta]-1;
/*      pLoadSimplePres(&sP             , &p
                      , tA              , lKsi
                      , presC           , gradPresC    
                      , s               , e
                      , t               , lNormal  
                      , densityC        , velC                      
                      , lModEta         , lModKsi         
                      , &loadsPres[nCarg], ndm
                      , true);*/
/*...................................................................*/
      }
/*...................................................................*/

/*... velocidades*/
      if(lFaceVelR[nAresta] > 0){
//      nCarg = lFaceVelL[nAresta]-1;
/*      pLoadSimple(&sP            , &p
                  , tA             , &ddum
                  , velC           , &ddum
                  , ddum           , &ddum
                  , ddum           , ddum           
                  , &ddum          , &ddum   
                  , &ddum          , lNormal    
                  , densityC       
                  , lModEta        , dcca[nAresta]
                  , loadsVel[nCarg], ndm
                  , false          , true
                  , false          , 0);   */
      } 
/*...................................................................*/
    }
/*...................................................................*/
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

/*... residuo de massa por celula*/
  rCell = p;
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
 * Data de criacao    : 17/09/2017                                   *
 * Data de modificaco : 17/07/2018                                   * 
 *-------------------------------------------------------------------* 
 * CELLSIMPLEPRES2DLM: Celula 2D para equacao de correcao de pressoa *
 * metodo simple em escoamento levemene compressivel                 *                      * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  * 
 * loadsPres -> definicoes de cargas de pressao                      * 
 * diffPres  -> tecnica da discretizacao do termo difusivo           *
 * lnFace    -> numero de faces da celula central e seus vizinhos    * 
 * eMass     -> termos/modelos da equacao de mass                    *
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
 * lDensity  -> massa especifica sem variacao temporal               * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceVelR  -> restricoes por elemento de velocidades               * 
 * faceVelL  -> carga por elemento de velocidades                    * 
 * facePresR -> restricoes por elemento de pressao                   * 
 * facePresL -> carga por elemento de pressao                        * 
 * pres      -> campo de pressao conhecido                           * 
 * gradPes   -> gradiente reconstruido da pressao                    * 
 * vel       -> campo de velocidade conhecido                        * 
 * gradVel   -> gradiente rescontruido das velocidades               *   
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
void cellSimplePres2DLm(Loads *lVel      , Loads *lPres 
					 , Diffusion *diffPres	       , MassEqModel *eMass
           , short *RESTRICT lGeomType   
           , INT *RESTRICT lViz          , INT *RESTRICT lId  
           , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
           , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT mEta
           , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT area   
           , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
           , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT lDensity
           , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
           , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
           , DOUBLE *RESTRICT lRcell     , Temporal const ddt 
           , short  *RESTRICT lFaceVelR  , short *RESTRICT lFaceVelL
           , short  *RESTRICT lFacePresR , short *RESTRICT lFacePresL
           , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres 
           , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT dField  
           , DOUBLE *RESTRICT temp 
           , const short nEn             , short const nFace    
           , const short ndm             , INT const nel)
{ 
/*...*/
	short iCodDif = diffPres->iCod;
/*...*/
  short idCell = nFace;
  short nAresta, nCarg, typeTime;
  bool  fTime,fLhsDensity,fRhsDensity;
  INT vizNel;
/*...*/
  DOUBLE densityC,densityC0,densityC00,densityV ,density;
  DOUBLE dFieldC[2],dFieldV[2],dFieldF[2];
/*...*/
  DOUBLE rCell,p,sP, dt, dt0, tmp, tmp0, tmp00, tmp1;
/*...*/
  DOUBLE v[2],lKsi[2],lNormal[2];
  DOUBLE dPviz,lModKsi,lModEta;
/*...*/
  DOUBLE gradPresC[2],gradPresV[2];
/*...*/
  DOUBLE dfd,coef,lvSkew[2];
/*... nonOrtogonal*/
	DOUBLE e[2],t[2],s[2],modE;
/*... interpolacao linear*/
  DOUBLE alpha,alphaMenosUm;
/*... */
  DOUBLE wfn,velC[2],velF[2],presC,presV,tempC;

/*...*/
  dt        = ddt.dt[0];  
  dt0       = ddt.dt[1];
  typeTime  = ddt.type;
  fTime     = ddt.flag;
  fLhsDensity = eMass->LhsDensity;
  fRhsDensity = eMass->RhsDensity; 
/*...................................................................*/

/*...*/
  densityC00 = MAT2D(idCell,0,lDensity,3);
  densityC0  = MAT2D(idCell,1,lDensity,3);
  densityC   = MAT2D(idCell,2,lDensity,3);
/*...................................................................*/

/*...*/
  velC[0]      = MAT2D(idCell,0,vel,ndm);
  velC[1]      = MAT2D(idCell,1,vel,ndm);
  
  gradPresC[0] = MAT2D(idCell,0,gradPres,ndm);
  gradPresC[1] = MAT2D(idCell,1,gradPres,ndm);

  presC     = pres[idCell];
  tempC      = temp[idCell];
  
  dFieldC[0] = MAT2D(idCell,0,dField,2); 
  dFieldC[1] = MAT2D(idCell,1,dField,2);
/*...................................................................*/

  p          = 0.0e0;
  sP         = 0.0e0;
  rCell      = 0.0e0;
  for(nAresta=0;nAresta<nFace;nAresta++){
    vizNel     = lViz[nAresta];
    lNormal[0] = MAT2D(nAresta,0,normal,ndm);
    lNormal[1] = MAT2D(nAresta,1,normal,ndm);
    lModEta    = mEta[nAresta];
    lKsi[0] = MAT2D(nAresta, 0, ksi, ndm);
    lKsi[1] = MAT2D(nAresta, 1, ksi, ndm);
    lModKsi = mKsi[nAresta];
/*... dominio*/
    if( vizNel  > -1 ){
/*...*/
      velF[0]       = MAT2D(nAresta,0,vel,ndm);
      velF[1]       = MAT2D(nAresta,1,vel,ndm);
      lvSkew[0]     = MAT2D(nAresta,0,vSkew,ndm);
      lvSkew[1]     = MAT2D(nAresta,1,vSkew,ndm);
      presV         = pres[nAresta];
      densityV      = MAT2D(nAresta,2,lDensity,3);
      dFieldV[0]    = MAT2D(nAresta,0,dField,2);
      dFieldV[1]    = MAT2D(nAresta,1,dField,2);
      gradPresV[0]  = MAT2D(nAresta,0,gradPres,ndm);
      gradPresV[1]  = MAT2D(nAresta,1,gradPres,ndm);
/*...................................................................*/

/*...*/
      v[0]         = lvSkew[0] + MAT2D(nAresta,0,xmcc,ndm);
      v[1]         = lvSkew[1] + MAT2D(nAresta,1,xmcc,ndm);
      dPviz        = sqrt(v[0]*v[0] + v[1]*v[1]);
      alpha        = dPviz/lModKsi;
      alphaMenosUm = 1.0e0 - alpha; 
/*...................................................................*/
      
/*... interpolacao das propriedades*/
      density    = alphaMenosUm*densityC   + alpha*densityV;
      dFieldF[0] = alphaMenosUm*dFieldC[0] + alpha*dFieldV[0];
      dFieldF[1] = alphaMenosUm*dFieldC[1] + alpha*dFieldV[1];
/*...................................................................*/

/*... termo difusivo
grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      s[0] = dFieldF[0]*lModEta*lNormal[0];
      s[1] = dFieldF[1]*lModEta*lNormal[1];
      difusionSchemeAnisotropic(s,lKsi,e,t,ndm,iCodDif);
/*...................................................................*/

/*... difusao direta*/
			coef = density;
			modE = sqrt(e[0]*e[0]+e[1]*e[1]);
  		dfd  = coef*modE/lModKsi;
/*...................................................................*/

/*...*/
      lA[nAresta] = dfd;
/*...................................................................*/

/*...*/
      wfn = interpolFaceVel(velC         ,velF
                           ,presC        ,presV
                           ,gradPresC    ,gradPresV
                           ,lNormal      ,lKsi
                           ,lModKsi      ,dFieldF
                           ,alphaMenosUm ,alpha
                           ,ndm);
/*...................................................................*/
      p   -= density*wfn*lModEta;
/*...................................................................*/
    }
/*... contorno*/
    else{
      lA[nAresta] = 0.0e0;
      wfn = velC[0]*lNormal[0] + velC[1]*lNormal[1];
      if(lFacePresR[nAresta]){
/*... termo difusivo
        grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
        s[0] = dFieldC[0]*lModEta*lNormal[0];
        s[1] = dFieldC[1]*lModEta*lNormal[1];
/*...*/
        difusionSchemeAnisotropic(s,lKsi,e,t,ndm,iCodDif);
/*...................................................................*/

/*...*/
        modE = sqrt(e[0]*e[0] + e[1]*e[1]);
/*...................................................................*/

/*...cargas*/
        nCarg = lFacePresL[nAresta]-1;
/*      pLoadSimplePres(&sP             , &p
                      , tA              , lKsi
                      , presC           , gradPresC      
                      , s               , e
                      , t               , lNormal  
                      , densityC        , velC                      
                      , lModEta         , lModKsi         
                      , &lPres[nCarg]    , ndm
                      , true);*/
/*...................................................................*/
      }
/*...................................................................*/

/*... velocidades*/
      if(lFaceVelR[nAresta] > 0){
        nCarg = lFaceVelL[nAresta]-1;
/*      pLoadSimple(&sP            , &p
                  , tA             , &ddum
                  , velC           , &ddum
                  , ddum           , &ddum
                  , ddum           , ddum           
                  , &ddum          , &ddum   
                  , &ddum          , lNormal    
                  , densityC
                  , lModEta        , dcca[nAresta]
                  , loadsVel[nCarg], ndm
                  , false          , true
                  , false          , 0);     */
      } 
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if(fTime){
/*... EULER*/
    if (typeTime == EULER){    
      tmp = (densityC - densityC0)/dt;
      if( fRhsDensity) p  -= tmp*area[idCell];
/*...*/
      tempC = CELSIUS_FOR_KELVIN(tempC);
      tmp1 = 1.e0/(IDEALGASR*tempC);
      if( fLhsDensity) sP += tmp1*area[idCell]/dt;
    } 
/*...BACKWARD*/
    else if (typeTime == BACKWARD) {
      tmp1  = dt + dt0 ;
      tmp   = 1.e0 / dt + 1.e0 / tmp1;
      tmp0  = -tmp1/(dt*dt0);
      tmp00 = dt/( dt0*tmp1 );
      tmp1 = tmp*densityC + tmp0*densityC0 + tmp00*densityC00; 
      if( fRhsDensity) p -= tmp1*area[idCell];
/*...*/
      tempC = CELSIUS_FOR_KELVIN(tempC);
      tmp1 = 1.e0/(IDEALGASR*tempC);
      if( fLhsDensity) sP += tmp*area[idCell];
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

/*... residuo de massa por celula*/
  rCell = p;
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
 * Data de criacao    : 07/08/2016                                   *
 * Data de modificaco : 15/08/2016                                   *
 *-------------------------------------------------------------------*
 * CELLSIMPLENONORTHPRES2D: correcao nao ortogonal da Celula 2D      *
 * para equacao de correcao de pressao                               *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
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
 * lDensity  -> massa especifica sem variacao temporal               *
 * lB        -> nao definido                                         *
 * pres      -> campo de pressao conhecido                           *
 * gradPes   -> gradiente reconstruido da pressao                    *
 * dField    -> matriz D do metodo simple                            *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * Fonte: Ferziger-Precic                                            *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellSimpleNonOrthPres2D(Diffusion diffPres
              ,short *RESTRICT lGeomType
              ,DOUBLE *RESTRICT prop    ,INT *RESTRICT lViz
              ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
              ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT mEta
              ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT area
              ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
              ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT lB      
              ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres
              ,DOUBLE *RESTRICT dField  ,DOUBLE *RESTRICT cc
              ,const short nEn          ,short const nFace
              ,const short ndm          ,INT const nel)
{ 

  DOUBLE densityC,densityV ,density;
  DOUBLE dFieldC[2],dFieldV[2],dFieldF[2];
/*...*/
  DOUBLE p,coef;
/*...*/
  DOUBLE v[2],lKsi[2],lNormal[2],presC,presV;
//DOUBLE ccV[2],ccC[2],lXm[2],lXmcc[2];
  DOUBLE dPviz,lModKsi,lModEta,du,duDksi;
/*...*/
  DOUBLE gradPresC[2],gradPresV[2],gradPresComp[2],gf[2],gfKsi;
/*...*/
  DOUBLE lvSkew[2];
/*... nonOrtogonal*/
	DOUBLE e[2],t[2],s[2];
/*... interpolacao linear*/
  DOUBLE alpha,alphaMenosUm;
/*...*/
	short iCodDif = diffPres.iCod;
/*...*/
  short idCell = nFace;
  short nAresta;
  INT vizNel;

/*...*/
  densityC = lDensity[idCell];
/*...................................................................*/

/*...*/
  
  gradPresC[0] = MAT2D(idCell,0,gradPres,ndm);
  gradPresC[1] = MAT2D(idCell,1,gradPres,ndm);
      
//ccC[0]     =  MAT2D(idCell,0,cc,ndm);
//ccC[1]     =  MAT2D(idCell,1,cc,ndm);

  presC      =  pres[idCell];

  dFieldC[0] = MAT2D(idCell, 0, dField, 2);
  dFieldC[1] = MAT2D(idCell, 1, dField, 2);
/*...................................................................*/

  p          = 0.0e0;
  for(nAresta=0;nAresta<nFace;nAresta++){
    vizNel     = lViz[nAresta];
/*... dominio*/
    if( vizNel  > -1 ){
/*...*/
      presV      = pres[nAresta];
/*...*/
      lNormal[0] = MAT2D(nAresta,0,normal,ndm);
      lNormal[1] = MAT2D(nAresta,1,normal,ndm);
      
      lModEta    = mEta[nAresta];
      
      dFieldV[0] = MAT2D(nAresta, 0, dField, 2);
      dFieldV[1] = MAT2D(nAresta, 1, dField, 2);
      
      lKsi[0]    = MAT2D(nAresta,0,ksi,ndm);
      lKsi[1]    = MAT2D(nAresta,1,ksi,ndm);

      lModKsi    = mKsi[nAresta];

      lvSkew[0]  = MAT2D(nAresta,0,vSkew,ndm);
      lvSkew[1]  = MAT2D(nAresta,1,vSkew,ndm);

      duDksi     = (presV - presC) / lModKsi;

      densityV   = lDensity[nAresta];

/*    ccV[0]     =  MAT2D(nAresta,0,cc,ndm);
      ccV[1]     =  MAT2D(nAresta,1,cc,ndm);

      lXm[0]     =  MAT2D(nAresta,0,xm,ndm);
      lXm[1]     =  MAT2D(nAresta,1,xm,ndm);

      lXmcc[0]   = MAT2D(nAresta,0,xmcc,ndm);
      lXmcc[1]   = MAT2D(nAresta,1,xmcc,ndm);*/
      
      gradPresV[0] =  MAT2D(nAresta,0,gradPres,ndm);
      gradPresV[1] =  MAT2D(nAresta,1,gradPres,ndm);
/*...................................................................*/

/*...*/
      v[0]         = lvSkew[0] + MAT2D(nAresta,0,xmcc,ndm);
      v[1]         = lvSkew[1] + MAT2D(nAresta,1,xmcc,ndm);
      dPviz        = sqrt(v[0]*v[0] + v[1]*v[1]);
      alpha        = dPviz/lModKsi;
      alphaMenosUm = 1.0e0 - alpha; 
/*...................................................................*/

/*... interpolacao das propriedades*/
      density  = alphaMenosUm*densityC  + alpha*densityV;
      dFieldF[0] = alphaMenosUm*dFieldC[0] + alpha*dFieldV[0];
      dFieldF[1] = alphaMenosUm*dFieldC[1] + alpha*dFieldV[1];
/*...................................................................*/

/*... termo difusivo
grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      s[0] = dFieldF[0] * lModEta*lNormal[0];
      s[1] = dFieldF[1] * lModEta*lNormal[1];
      difusionSchemeAnisotropic(s, lKsi, e, t, ndm, iCodDif);
/*...................................................................*/

/*... difusao direta*/
      coef = density;
/*...................................................................*/

/*...*/
      gf[0] = alphaMenosUm*gradPresC[0] + alpha*gradPresV[0];
      gf[1] = alphaMenosUm*gradPresC[1] + alpha*gradPresV[1];
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi = gf[0]*lKsi[0] + gf[1]*lKsi[1];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
      du              = duDksi  - gfKsi;
      gradPresComp[0] = gf[0] + du*lKsi[0];   
      gradPresComp[1] = gf[1] + du*lKsi[1];
/*...................................................................*/

/*...*/
      gfKsi = gradPresComp[0]*t[0]
            + gradPresComp[1]*t[1];
/*...................................................................*/

/*...*/
      p += coef*gfKsi;
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  *lB = p;
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 11/07/2016                                   *
 * Data de modificaco : 06/10/2019                                   *
 *-------------------------------------------------------------------*
 * CELLSIMPLEPRES3D: Celula 3D para equacao de correcao de pressao   *
 * metodo simple em escoamento imcompressivel                        *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  *
 * loadsPres -> definicoes de cargas de pressao                      *
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
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
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * faceVelR  -> restricoes por elemento de velocidades               *
 * facePresR -> restricoes por elemento de pressao                   *
 * pres      -> campo de pressao conhecido                           *
 * gradPes   -> gradiente reconstruido da pressao                    *
 * vel       -> campo de velocidade conhecido                        *
 * dField    -> matriz D do metodo simple                            *
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  *
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
 * Fonte: Ferziger-Precic                                            *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellSimplePres3D(Loads *lVel          ,Loads *lPres 
							,Diffusion *diffPres         
              ,short *RESTRICT lGeomType   ,DOUBLE *RESTRICT prop
              ,INT *RESTRICT lViz          ,INT *RESTRICT lId  
              ,DOUBLE *RESTRICT ksi        ,DOUBLE *RESTRICT mKsi
              ,DOUBLE *RESTRICT eta        ,DOUBLE *RESTRICT fArea
              ,DOUBLE *RESTRICT normal     ,DOUBLE *RESTRICT volume
              ,DOUBLE *RESTRICT xm         ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca       
              ,DOUBLE *RESTRICT vSkew      ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT lA         ,DOUBLE *RESTRICT lB
              ,DOUBLE *RESTRICT lRcell
              ,short  *RESTRICT lFaceVelR  ,short  *RESTRICT lFacePresR 
              ,DOUBLE *RESTRICT pres       ,DOUBLE *RESTRICT gradPres
              ,DOUBLE *RESTRICT vel            
              ,DOUBLE *RESTRICT dField     ,DOUBLE *RESTRICT wallPar
              ,const short nEn             ,short const nFace
              ,const short ndm             ,INT const nel)
{ 
/*...*/
	short iCodDif = diffPres->iCod,iCodPolFace;
/*...*/
  short idCell = nFace;
  short i, nf, nCarg;
/*...*/
  INT vizNel;
/*...*/
  DOUBLE densityC,densityV ,density;
  DOUBLE dFieldC[3],dFieldV[3],dFieldF[3];
/*...*/
  DOUBLE rCell,p,sP;
  DOUBLE tA[3],ddum=0.e0;
/*...*/
  DOUBLE lKsi[3],lXmcc[3],lNormal[3],lModKsi,lFarea;
/*...*/
  DOUBLE gradPresC[3],gradPresV[3],presC,presV;
/*...*/
  DOUBLE dfd,coef,lvSkew[3];
/*... nonOrtogonal*/
	DOUBLE e[3],t[3],s[3],modE;
/*... interpolacao linear*/
  DOUBLE alpha,alphaMenosUm;
/*... */
  DOUBLE wfn,velC[3],velF[3];
/*...*/
  DOUBLE ts,xx[4]; 

/*...*/
  densityC     =  MAT2D(idCell, DENSITY         ,prop,MAXPROP);
/*...................................................................*/

/*...*/
  iCodPolFace = INTPOLFACELINEAR;
/*...................................................................*/

/*...*/
  ts = 0.e0;
/*...................................................................*/

/*...*/
  presC     = pres[idCell];
  for (i = 0; i < 3; i++)
  {
    velC[i]      = MAT2D(idCell,i,vel,3);
  
    gradPresC[i] = MAT2D(idCell,i,gradPres,3);

    dFieldC[i] = MAT2D(idCell,i,dField,3);
  }
/*...................................................................*/

/*...*/
  p = sP = 0.0e0;
  for(nf=0;nf<nFace;nf++){
    vizNel     = lViz[nf];
    
    lNormal[0] = MAT2D(nf,0,normal,ndm);
    lNormal[1] = MAT2D(nf,1,normal,ndm);
    lNormal[2] = MAT2D(nf,2,normal,ndm);
    lFarea     = fArea[nf];

    lKsi[0] = MAT2D(nf, 0, ksi, ndm);
    lKsi[1] = MAT2D(nf, 1, ksi, ndm);
    lKsi[2] = MAT2D(nf, 2, ksi, ndm);
    lModKsi = mKsi[nf];
    lXmcc[0] = MAT2D(nf,0,xmcc,3);
    lXmcc[1] = MAT2D(nf,1,xmcc,3);
    lXmcc[2] = MAT2D(nf,2,xmcc,3);

/*... dominio*/
    if( vizNel  > -1 ){
/*...*/
      velF[0]  = MAT2D(nf,0,vel,3);
      velF[1] = MAT2D(nf,1,vel,3);
      velF[2] = MAT2D(nf,2,vel,3);

      lvSkew[0] = MAT2D(nf,0,vSkew,3);
      lvSkew[1] = MAT2D(nf,1,vSkew,3);
      lvSkew[2] = MAT2D(nf,2,vSkew,3);

      presV      = pres[nf];
      densityV   =  MAT2D(nf, DENSITY         ,prop,MAXPROP);

      dFieldV[0] = MAT2D(nf, 0, dField, 3);
      dFieldV[1] = MAT2D(nf, 1, dField, 3);
      dFieldV[2] = MAT2D(nf, 2, dField, 3);

      gradPresV[0]  = MAT2D(nf,0,gradPres,3);
      gradPresV[1]  = MAT2D(nf,1,gradPres,3);
      gradPresV[2]  = MAT2D(nf,2,gradPres,3);
/*...................................................................*/

/*...*/
      alpha = interpolFace(lvSkew           ,lXmcc
                          ,volume[idCell]   ,volume[nf]
                          ,lModKsi          ,ndm
                          ,iCodPolFace);
      alphaMenosUm = 1.0e0 - alpha;
/*...................................................................*/

/*... interpolacao das propriedades*/
      density    = alphaMenosUm*densityC   + alpha*densityV;
      dFieldF[0] = alphaMenosUm*dFieldC[0] + alpha*dFieldV[0];
      dFieldF[1] = alphaMenosUm*dFieldC[1] + alpha*dFieldV[1];
      dFieldF[2] = alphaMenosUm*dFieldC[2] + alpha*dFieldV[2];
/*...................................................................*/
 
/*... termo difusivo
        grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      s[0] = dFieldF[0]*lFarea*lNormal[0];
      s[1] = dFieldF[1]*lFarea*lNormal[1];
      s[2] = dFieldF[2]*lFarea*lNormal[2];
      difusionSchemeAnisotropic(s,lKsi,e,t,ndm,iCodDif);
/*...................................................................*/

/*... difusao direta*/
			coef = density;
			modE = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
			dfd  = coef*modE/lModKsi;
/*...................................................................*/

/*...*/
      lA[nf] = dfd;
/*...................................................................*/

/*...*/
      wfn = interpolFaceVel(velC         ,velF
                           ,presC        ,presV
                           ,gradPresC    ,gradPresV
                           ,lNormal      ,lKsi
                           ,lModKsi      ,dFieldF
                           ,alphaMenosUm ,alpha
                           ,ndm);
/*...................................................................*/

/*...*/
      p   -= density*wfn*lFarea;
/*...................................................................*/
    }
/*...................................................................*/

/*... contorno*/
    else
    {
      xx[0] = MAT2D(nf, 0, xm, 3);
      xx[1] = MAT2D(nf, 1, xm, 3);
      xx[2] = MAT2D(nf, 2, xm, 3);
      xx[3] = ts;
      lA[nf] = 0.0e0;
      wfn = velC[0]*lNormal[0] 
          + velC[1]*lNormal[1] 
          + velC[2]*lNormal[2];
      if(lFacePresR[nf])
      {
/*... termo difusivo
      grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
        s[0] = dFieldC[0]*lFarea*lNormal[0];
        s[1] = dFieldC[1]*lFarea*lNormal[1];
        s[2] = dFieldC[2]*lFarea*lNormal[2];
/*...*/
        difusionSchemeAnisotropic(s, lKsi, e, t, ndm, iCodDif);
//      difusionSchemeAnisotropic(s, lKsi, e, t, ndm, 1);
/*...................................................................*/

/*...cargas*/
        nCarg = lFacePresR[nf] - 1;
        pLoadSimplePres(&sP             , &p
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
                      , true); 
/*...................................................................*/
      }
/*...................................................................*/

/*... velocidades*/
      if(lFaceVelR[nf])
      {
        nCarg = lFaceVelR[nf]-1;
        pLoadSimple(&sP             , &p
                  , tA              , lXmcc
                  , velC            , &ddum       
                  , presC           , gradPresC
                  , ddum            , ddum  
                  , xx         
                  , s               , e       
                  , t               , lNormal    
                  , densityC        , wallPar
                  , lFarea          , lModKsi
                  , &loadsVel[nCarg], ndm
                  , false           , true
                  , false           , false       
                  , 0               , nel);     
      } 
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  lA[idCell] = sP;
  for (nf = 0; nf<nFace; nf++)
    lA[idCell] += lA[nf];
/*...................................................................*/

/*... residuo de massa por celula*/
  rCell = p;
/*...................................................................*/

/*...*/
  for (nf = 0; nf<nFace; nf++)
    lA[nf] *= -1.e0;
/*...................................................................*/

/*...*/
  lB[0]     = p;
  lRcell[0] = rCell;
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 03/10/2017                                   *
 * Data de modificaco : 13/10/2019                                   *
 *-------------------------------------------------------------------*
 * CELLSIMPLEPRES3D: Celula 3D para equacao de correcao de pressao   *
 * metodo simple em escoamento levemene compressivel                 *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  *
 * loadsPres -> definicoes de cargas de pressao                      *
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * momentumModel -> termos/modelos da equacao de momento linear      *
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
 * lDensity  -> massa especifica sem variacao temporal               *
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * faceVelR  -> restricoes por elemento de velocidades               *
 * faceVelL  -> carga por elemento de velocidades                    *
 * facePresR -> restricoes por elemento de pressao                   *
 * facePresL -> carga por elemento de pressao                        *
 * pres      -> campo de pressao conhecido                           *
 * gradPes   -> gradiente reconstruido da pressao                    *
 * gradRho   -> gradiente rescontruido das densidade                 *
 * mMolar    -> massa molar da mistura                               *
 * vel       -> campo de velocidade conhecido                        *
 * gradVel   -> gradiente rescontruido das velocidades               *
 * dField    -> matriz D do metodo simple                            *
 * temp      -> temperatura                                          *  
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  *  
 * densityMed-> densidade media do meio                              *
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
 * Fonte: Ferziger-Precic                                            *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellSimplePres3DLm(Loads *lVel        , Loads *lPres 
							, Diffusion *diffPres        
              , MassEqModel *eMass         , MomentumModel *momentumModel
              , short *RESTRICT lGeomType  
              , INT *RESTRICT lViz         , INT *RESTRICT lId  
              , DOUBLE *RESTRICT ksi       , DOUBLE *RESTRICT mKsi
              , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT fArea
              , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume
              , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc
              , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT cc   
              , DOUBLE *RESTRICT lDensity
              , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew
              , DOUBLE *RESTRICT lA        , DOUBLE *RESTRICT lB
              , DOUBLE *RESTRICT lRcell    , Temporal *ddt 
              , short  *RESTRICT lFaceVelR , short  *RESTRICT lFacePresR
              , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres
              , DOUBLE *RESTRICT gradRho   , DOUBLE *RESTRICT mMolar  
              , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT dField
              , DOUBLE *RESTRICT temp      , DOUBLE *RESTRICT wallPar
              , DOUBLE const densityMed
              , const short nEn            , short const nFace
              , const short ndm            , INT const nel)
{ 

/*...*/
	short iCodDif = diffPres->iCod,iCodPolFace,iCodBuoyant;
/*...*/
  short idCell = nFace;
  short nf, nCarg, typeTime;
  bool  fTime,fLhsDensity,fRhsDensity;
  INT vizNel;
/*...*/
  DOUBLE densityC,densityC0,densityC00,densityV ,density;
  DOUBLE dFieldC[3],dFieldV[3],dFieldF[3];
/*...*/
  DOUBLE rCell,p,sP, dt, dt0, tmp, tmp0, tmp00, tmp1;
  DOUBLE tA[3],ddum=0.e0;
/*...*/
  DOUBLE lKsi[3],lXmcc[3],lNormal[3],lModKsi,lFarea;
/*...*/
  DOUBLE gradPresC[3],gradPresV[3],presC,presV;
/*...*/
  DOUBLE dfd,coef,lvSkew[3];
/*... nonOrtogonal*/
  DOUBLE e[3],t[3],s[3],modE;
/*... interpolacao linear*/
  DOUBLE alpha,alphaMenosUm;
/*... */
  DOUBLE wfn,velC[3],velF[3],tempC;
/*...*/
  DOUBLE ts,xx[4];
/*...*/
  DOUBLE g[3],gh;
/*...*/
//DOUBLE mW,mW0,mW00;
  DOUBLE mW;

/*...*/
  ts          = ddt->t;
  dt          = ddt->dt[0];  
  dt0         = ddt->dt[1];
  typeTime    = ddt->type;
  fTime       = ddt->flag;
  fLhsDensity = eMass->LhsDensity;
  fRhsDensity = eMass->RhsDensity; 
  iCodPolFace = INTPOLFACELINEAR;
/*...................................................................*/

/*...*/

  iCodBuoyant = momentumModel->iCodBuoyant;
/*...................................................................*/

/*...*/
  densityC00 = MAT2D(idCell,TIME_N_MINUS_2,lDensity,DENSITY_LEVEL);
  densityC0  = MAT2D(idCell,TIME_N_MINUS_1,lDensity,DENSITY_LEVEL);
  densityC   = MAT2D(idCell,TIME_N        ,lDensity,DENSITY_LEVEL);
/*...................................................................*/

/*...*/
/*mW00       = mMolar[TIME_N_MINUS_2];
  mW0        = mMolar[TIME_N_MINUS_1];*/
  mW         = mMolar[TIME_N];
/*...................................................................*/

/*...*/
  velC[0] = MAT2D(idCell,0,vel,3);
  velC[1] = MAT2D(idCell,1,vel,3);
  velC[2] = MAT2D(idCell,2,vel,3);
  
  gradPresC[0] = MAT2D(idCell,0,gradPres,3);
  gradPresC[1] = MAT2D(idCell,1,gradPres,3);
  gradPresC[2] = MAT2D(idCell,2,gradPres,3);

  presC = pres[idCell];
  tempC = temp[idCell];
  
  dFieldC[0] = MAT2D(idCell,0,dField,3); 
  dFieldC[1] = MAT2D(idCell,1,dField,3);
  dFieldC[2] = MAT2D(idCell,2,dField,3);

  g[0] = gravity[0];
  g[1] = gravity[1];
  g[2] = gravity[2];

  gh  = (MAT2D(idCell,0,cc,3) - xRef[0])*g[0]
      + (MAT2D(idCell,1,cc,3) - xRef[1])*g[1]
      + (MAT2D(idCell,2,cc,3) - xRef[2])*g[2];
/*...................................................................*/

/*...*/
  p = sP = 0.0e0;
  for(nf=0;nf<nFace;nf++){
    vizNel     = lViz[nf];
    
    lNormal[0] = MAT2D(nf,0,normal,3);
    lNormal[1] = MAT2D(nf,1,normal,3);
    lNormal[2] = MAT2D(nf,2,normal,3);
    
    lFarea     = fArea[nf];

    lKsi[0] = MAT2D(nf,0,ksi,3);
    lKsi[1] = MAT2D(nf,1,ksi,3);
    lKsi[2] = MAT2D(nf,2,ksi,3);
    lModKsi = mKsi[nf];

    lXmcc[0] = MAT2D(nf,0,xmcc,3);
    lXmcc[1] = MAT2D(nf,1,xmcc,3);
    lXmcc[2] = MAT2D(nf,2,xmcc,3);

/*... dominio*/
    if( vizNel  > -1 ){
/*...*/
      velF[0] = MAT2D(nf,0,vel,3);
      velF[1] = MAT2D(nf,1,vel,3);
      velF[2] = MAT2D(nf,2,vel,3);

/*..*/
      lvSkew[0] = MAT2D(nf,0,vSkew,3);
      lvSkew[1] = MAT2D(nf,1,vSkew,3);
      lvSkew[2] = MAT2D(nf,2,vSkew,3);

      presV    = pres[nf];
      densityV = MAT2D(nf,TIME_N,lDensity,3);

      dFieldV[0] = MAT2D(nf, 0, dField, 3);
      dFieldV[1] = MAT2D(nf, 1, dField, 3);
      dFieldV[2] = MAT2D(nf, 2, dField, 3);

      gradPresV[0]  = MAT2D(nf,0,gradPres,3);
      gradPresV[1]  = MAT2D(nf,1,gradPres,3);
      gradPresV[2]  = MAT2D(nf,2,gradPres,3);
/*...................................................................*/

/*...*/
      alpha = interpolFace(lvSkew           ,lXmcc
                          ,volume[idCell]   ,volume[nf]
                          ,lModKsi          ,ndm
                          ,iCodPolFace);
      alphaMenosUm = 1.0e0 - alpha;
/*...................................................................*/

/*... interpolacao das propriedades*/
      density    = alphaMenosUm*densityC   + alpha*densityV;
      dFieldF[0] = alphaMenosUm*dFieldC[0] + alpha*dFieldV[0];
      dFieldF[1] = alphaMenosUm*dFieldC[1] + alpha*dFieldV[1];
      dFieldF[2] = alphaMenosUm*dFieldC[2] + alpha*dFieldV[2];
/*...................................................................*/
 
/*... termo difusivo
        grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      s[0] = dFieldF[0]*lFarea*lNormal[0];
      s[1] = dFieldF[1]*lFarea*lNormal[1];
      s[2] = dFieldF[2]*lFarea*lNormal[2];
      difusionSchemeAnisotropic(s,lKsi,e,t,ndm,iCodDif);
/*...................................................................*/

/*... difusao direta*/
			coef = density;
			modE = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
			dfd  = coef*modE/lModKsi;
/*...................................................................*/

/*...*/
      lA[nf] = dfd;
/*...................................................................*/

/*...*/
      wfn = interpolFaceVel(velC         ,velF
                           ,presC        ,presV
                           ,gradPresC    ,gradPresV
                           ,lNormal      ,lKsi
                           ,lModKsi      ,dFieldF
                           ,alphaMenosUm ,alpha
                           ,ndm);
/*...................................................................*/

/*...*/
      p   -= density*wfn*lFarea;
/*...................................................................*/
    }
/*...................................................................*/

/*... contorno*/
    else
    {
      xx[0] = MAT2D(nf, 0, xm, 3);
      xx[1] = MAT2D(nf, 1, xm, 3);
      xx[2] = MAT2D(nf, 2, xm, 3);
      xx[3] = ts;
      lA[nf] = 0.0e0;
      wfn = velC[0]*lNormal[0] 
          + velC[1]*lNormal[1] 
          + velC[2]*lNormal[2];
      if(lFacePresR[nf])
      {
/*... termo difusivo
      grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
        s[0] = dFieldC[0]*lFarea*lNormal[0];
        s[1] = dFieldC[1]*lFarea*lNormal[1];
        s[2] = dFieldC[2]*lFarea*lNormal[2];
/*...*/
        difusionSchemeAnisotropic(s, lKsi, e, t, ndm, iCodDif);
/*...................................................................*/

/*...cargas*/
        nCarg = lFacePresR[nf] - 1;
        pLoadSimplePres(&sP             , &p
                      , tA              , lXmcc
                      , presC           , gradPresC   
                      , s               , e         
                      , t               , lNormal
                      , g               , gradRho
                      , velC            , gh  
                      , densityC        , densityMed
                      , lFarea          , dcca[nf]
                      , &lPres[nCarg]   , ndm
                      , iCodBuoyant     , true      
                      , true); 
/*...................................................................*/
      }
/*...................................................................*/

/*... velocidades*/
      if(lFaceVelR[nf])
      {
        nCarg = lFaceVelR[nf]-1;
        pLoadSimple(&sP             , &p
                  , tA              , lXmcc
                  , velC            , &ddum
                  , presC           , gradPresC
                  , ddum            , ddum  
                  , xx         
                  , s               , e       
                  , t               , lNormal    
                  , densityC        , wallPar
                  , lFarea          , lModKsi
                  , &loadsVel[nCarg], ndm
                  , false           , true
                  , false           , false       
                  , 0               , nel);      
      } 
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if(fTime)
  {
/*... EULER*/
    if (typeTime == EULER){    
      tmp = (densityC - densityC0)/dt;
      if( fRhsDensity) p -= tmp*volume[idCell];
/*...*/
      tempC = CELSIUS_FOR_KELVIN(tempC);
//    tmp1 = 1.e0/(IDEALGASR*tempC);
      tmp1 = fKsi(mW,tempC,IDEALGASR);
      if( fLhsDensity) sP += tmp1*volume[idCell]/dt;
    } 
/*...BACKWARD*/
    else if (typeTime == BACKWARD) 
    {
      tmp1  = dt + dt0 ;
      tmp   = 1.e0 / dt + 1.e0 / tmp1;
      tmp0  = -tmp1/(dt*dt0);
      tmp00 = dt/( dt0*tmp1 );
      tmp1 = tmp*densityC + tmp0*densityC0 + tmp00*densityC00; 
      if(fRhsDensity) p -= tmp1*volume[idCell];
/*...*/
      tempC = CELSIUS_FOR_KELVIN(tempC);
      tmp1 = fKsi(mW,tempC,IDEALGASR);
      if(fLhsDensity) sP += tmp1*tmp*volume[idCell];
    }
  }
/*...................................................................*/

/*...*/
  lA[idCell] = sP;
  for (nf = 0; nf<nFace; nf++)
    lA[idCell] += lA[nf];
/*...................................................................*/

/*... residuo de massa por celula*/
  rCell = p;
/*...................................................................*/

/*...*/
  for (nf = 0; nf<nFace; nf++)
    lA[nf] *= -1.e0;
/*...................................................................*/

/*...*/
  lB[0]     = p;
  lRcell[0] = rCell;
/*...................................................................*/
}
/*********************************************************************/


/*********************************************************************
 * Data de criacao    : 02/08/2016                                   *
 * Data de modificaco : 09/08/2016                                   *
 *-------------------------------------------------------------------*
 * CELLSIMPLENONORTHPRES3D: correcao nao ortogonal da Celula 3D      *
 * para equacao de correcao de pressao                               *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * lnFace    -> numero de faces da celula central e seus vizinhos    *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lId       -> numeracoes das equacoes das celulas                  *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * mEta      -> area das faces                                       *
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
 * lDensity  -> massa especifica sem variacao temporal               *
 * lB        -> nao definido                                         *
 * pres      -> campo de pressao conhecido                           *
 * gradPes   -> gradiente reconstruido da pressao                    *
 * dField    -> matriz D do metodo simple                            *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * Fonte: Ferziger-Precic                                            *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellSimpleNonOrthPres3D(Diffusion diffPres
              ,short *RESTRICT lGeomType
              ,DOUBLE *RESTRICT prop    ,INT *RESTRICT lViz
              ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
              ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea
              ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
              ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
              ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT lB      
              ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres
              ,DOUBLE *RESTRICT dField  ,DOUBLE *RESTRICT cc
              ,const short nEn          ,short const nFace
              ,const short ndm          ,INT const nel)
{ 

  DOUBLE densityC,densityV ,density;
  DOUBLE dFieldC[3],dFieldV[3],dFieldF[3];
/*...*/
  DOUBLE p,coef;
/*...*/
  DOUBLE v[3],lKsi[3],lNormal[3],presC,presV;
/*DOUBLE ccV[3],ccC[3],lXm[3],lXmcc[3];*/
  DOUBLE dPviz,lModKsi,lFarea,du,duDksi;
/*...*/
  DOUBLE gradPresC[3],gradPresV[3],gradPresComp[3],gf[3],gfKsi;
/*... nonOrtogonal*/
	DOUBLE e[3], t[3],s[4];
/*...*/
  DOUBLE lvSkew[3];
/*... interpolacao linear*/
  DOUBLE alpha,alphaMenosUm;
/*...*/
	short iCodDif = diffPres.iCod;
/*...*/
  short idCell = nFace;
  short nAresta;
  INT vizNel;

/*...*/
  densityC = lDensity[idCell];
/*...................................................................*/

/*...*/
  
  gradPresC[0] = MAT2D(idCell,0,gradPres,ndm);
  gradPresC[1] = MAT2D(idCell,1,gradPres,ndm);
  gradPresC[2] = MAT2D(idCell,2,gradPres,ndm);
      
/*ccC[0]     =  MAT2D(idCell,0,cc,ndm);
  ccC[1]     =  MAT2D(idCell,1,cc,ndm);
  ccC[2]     =  MAT2D(idCell,2,cc,ndm);*/

  presC      =  pres[idCell];


  dFieldC[0] = MAT2D(idCell, 0, dField, 3);
  dFieldC[1] = MAT2D(idCell, 1, dField, 3);
  dFieldC[2] = MAT2D(idCell, 2, dField, 3);
/*...................................................................*/

  p          = 0.0e0;
  for(nAresta=0;nAresta<nFace;nAresta++){
    vizNel     = lViz[nAresta];
/*... dominio*/
    if( vizNel  > -1 ){
/*...*/
      presV      = pres[nAresta];
/*...*/
      lNormal[0] = MAT2D(nAresta,0,normal,ndm);
      lNormal[1] = MAT2D(nAresta,1,normal,ndm);
      lNormal[2] = MAT2D(nAresta,2,normal,ndm);
      
      lFarea     = fArea[nAresta];
      
      dFieldV[0] = MAT2D(nAresta, 0, dField, 3);
      dFieldV[1] = MAT2D(nAresta, 1, dField, 3);
      dFieldV[2] = MAT2D(nAresta, 2, dField, 3);

      lKsi[0]    = MAT2D(nAresta,0,ksi,ndm);
      lKsi[1]    = MAT2D(nAresta,1,ksi,ndm);
      lKsi[2]    = MAT2D(nAresta,2,ksi,ndm);

      lModKsi    = mKsi[nAresta];

      lvSkew[0]  = MAT2D(nAresta,0,vSkew,ndm);
      lvSkew[1]  = MAT2D(nAresta,1,vSkew,ndm);
      lvSkew[2]  = MAT2D(nAresta,2,vSkew,ndm);

      duDksi     = (presV - presC) / lModKsi;

      densityV   = lDensity[nAresta];

/*    ccV[0]     =  MAT2D(nAresta,0,cc,ndm);
      ccV[1]     =  MAT2D(nAresta,1,cc,ndm);
      ccV[2]     =  MAT2D(nAresta,2,cc,ndm);

      lXm[0]     =  MAT2D(nAresta,0,xm,ndm);
      lXm[1]     =  MAT2D(nAresta,1,xm,ndm);
      lXm[2]     =  MAT2D(nAresta,2,xm,ndm);

      lXmcc[0]   = MAT2D(nAresta,0,xmcc,ndm);
      lXmcc[1]   = MAT2D(nAresta,1,xmcc,ndm);
      lXmcc[2]   = MAT2D(nAresta,2,xmcc,ndm);*/
      
      gradPresV[0] =  MAT2D(nAresta,0,gradPres,ndm);
      gradPresV[1] =  MAT2D(nAresta,1,gradPres,ndm);
      gradPresV[2] =  MAT2D(nAresta,2,gradPres,ndm);
/*...................................................................*/

/*...*/
      v[0]         = lvSkew[0] + MAT2D(nAresta,0,xmcc,ndm);
      v[1]         = lvSkew[1] + MAT2D(nAresta,1,xmcc,ndm);
      v[2]         = lvSkew[2] + MAT2D(nAresta,2,xmcc,ndm);
      dPviz        = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
      alpha        = dPviz/lModKsi;
      alphaMenosUm = 1.0e0 - alpha; 
/*...................................................................*/

/*... interpolacao das propriedades*/
      density    = alphaMenosUm*densityC   + alpha*densityV;
      dFieldF[0] = alphaMenosUm*dFieldC[0] + alpha*dFieldV[0];
      dFieldF[1] = alphaMenosUm*dFieldC[1] + alpha*dFieldV[1];
      dFieldF[2] = alphaMenosUm*dFieldC[2] + alpha*dFieldV[2];
/*...................................................................*/

/*... termo difusivo
      grad(phi)*S = (grad(phi)*E)Imp + (grad(phi)*T)Exp*/
      s[0] = dFieldF[0]*lFarea*lNormal[0];
      s[1] = dFieldF[1]*lFarea*lNormal[1];
      s[2] = dFieldF[2]*lFarea*lNormal[2];
      difusionSchemeAnisotropic(s, lKsi, e, t, ndm, iCodDif);
/*...................................................................*/

/*... difusao direta*/
      coef = density;
/*...................................................................*/

/*...*/
      gf[0] = alphaMenosUm*gradPresC[0] + alpha*gradPresV[0];
      gf[1] = alphaMenosUm*gradPresC[1] + alpha*gradPresV[1];
      gf[2] = alphaMenosUm*gradPresC[2] + alpha*gradPresV[2];
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi = gf[0]*lKsi[0] + gf[1]*lKsi[1] + gf[2]*lKsi[2];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
      du              = duDksi  - gfKsi;
      gradPresComp[0] = gf[0] + du*lKsi[0];
      gradPresComp[1] = gf[1] + du*lKsi[1];
      gradPresComp[2] = gf[2] + du*lKsi[2];
/*...................................................................*/

/*...*/
      gfKsi = gradPresComp[0]*t[0]
            + gradPresComp[1]*t[1]
            + gradPresComp[2]*t[2];
/*...................................................................*/

/*...*/
      p += coef*gfKsi;
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  *lB = p;
/*...................................................................*/
}
/*********************************************************************/
