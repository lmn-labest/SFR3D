#include<CellLoop.h>
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 22/08/2016                                   *
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
              ,short *restrict lGeomType,DOUBLE *restrict prop
              ,INT *restrict lViz       ,INT *restrict lId  
              ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
              ,DOUBLE *restrict eta     ,DOUBLE *restrict mEta
              ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
              ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
              ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
              ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
              ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
              ,DOUBLE *restrict lRcell  ,Temporal const ddt
              ,short  *restrict lFaceR  ,short *restrict lFaceL
              ,DOUBLE *restrict u0      ,DOUBLE *restrict gradU0
              ,DOUBLE *restrict vel     ,DOUBLE *restrict cc
              ,const short nEn          ,short const nFace    
              ,const short ndm          ,INT const nel)
{ 

  DOUBLE coefDifC,coefDif,coefDifV,rCell,dt,dt0;
  DOUBLE densityC,densityF,densityM;
  DOUBLE p,sP,nk,dfd,gfKsi,lvSkew[2];
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
      advectiveSchemeScalar(u0[idCell],u0[nAresta]
                           ,gradUp    ,gradUv
                           ,gradUcomp ,lvSkew
                           ,lXmcc     ,v    
                           ,lKsi      ,lModKsi
                           ,cv        ,&cvc
                           ,ndm
                           ,iCodAdv1  ,iCodAdv2);
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
              ,short *restrict lGeomType,DOUBLE *restrict prop
              ,INT *restrict lViz       ,INT *restrict lId  
              ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
              ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea
              ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
              ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
              ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
              ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
              ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
              ,DOUBLE *restrict lRcell  ,Temporal const ddt             
              ,short  *restrict lFaceR  ,short  *restrict lFaceL  
              ,DOUBLE *restrict u0      ,DOUBLE *restrict gradU0
              ,DOUBLE *restrict vel     ,DOUBLE *restrict cc
              ,const short nEn          ,short const nFace    
              ,const short ndm          ,INT const nel)
{ 

  DOUBLE coefDifC,coefDif,coefDifV,rCell,dt,dt0;
  DOUBLE densityC,densityF,densityM;
  DOUBLE p,sP,nk,dfd,gfKsi,lvSkew[3];
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
      advectiveSchemeScalar(u0[idCell],u0[nf]
                           ,gradUp    ,gradUv
                           ,gradUcomp ,lvSkew
                           ,lXmcc     ,v
                           ,lKsi      ,lModKsi
                           ,cv        ,&cvc
                           ,ndm
                           ,iCodAdv1, iCodAdv2);
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

/*
  lA[idCell] = sP;
  for(nf=0;nf<nFace;nf++){
    lA[idCell] += lA[nf];
  }
*/
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

