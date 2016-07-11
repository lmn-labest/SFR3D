#include<CellLoop.h>
/*********************************************************************
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 03/07/2016                                   * 
 *-------------------------------------------------------------------* 
 * CELLSIMPLEVELINC2D: Celula 2D para velocidade do metodo simple    * 
 * em escoamento imcompressivel                                      * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loadsVel  -> definicoes de cargas de velocidades                  * 
 * loadsPres -> definicoes de cargas de pressao                      * 
 * advVel    -> tecnica da discretizacao do termo advecao            * 
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
 * dField    -> matriz D do metodo simple                            * 
 * nEn       -> numero de nos da celula central                      *
 * cc        -> centroides da celula centra e seus vizinhos          * 
 * underU    -> parametro de sobre relaxamento                       * 
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
 *                  | du1dx1 du1dx2 |                                *
 * grad(viz1,*,*) = |               |                                *
 *                  | du2dx1 du2dx2 |                                *
 *                                                                   *  
 *                  | du1dx1 du1dx2 |                                *
 * grad(viz2,*,*) = |               |                                *
 *                  | du2dx1 du2dx2 |                                *
 *                                                                   *  
 *                  | du1dx1 du1dx2 |                                *
 * grad(P   ,*,*) = |               |                                *
 *                  | du2dx1 du2dx2 |                                *
 *********************************************************************/
void cellSimpleVel2D(Loads *loadsVel     ,Loads *loadsPres 
            ,Advection advVel            ,short const typeSimple 
            ,short *restrict lGeomType   ,DOUBLE *restrict prop
            ,INT *restrict lViz          ,INT *restrict lId  
            ,DOUBLE *restrict ksi        ,DOUBLE *restrict mKsi
            ,DOUBLE *restrict eta        ,DOUBLE *restrict mEta
            ,DOUBLE *restrict normal     ,DOUBLE *restrict area   
            ,DOUBLE *restrict xm         ,DOUBLE *restrict xmcc
            ,DOUBLE *restrict dcca       ,DOUBLE *restrict lDensity
            ,DOUBLE *restrict vSkew      ,DOUBLE *restrict mvSkew
            ,DOUBLE *restrict lA         ,DOUBLE *restrict lB
            ,DOUBLE *restrict lRcell     ,Temporal const ddt
            ,short  *restrict lFaceVelR  ,short *restrict lFaceVelL
            ,short  *restrict lFacePresR ,short *restrict lFacePresL
            ,DOUBLE *restrict pres       ,DOUBLE *restrict gradPres 
            ,DOUBLE *restrict vel        ,DOUBLE *restrict gradVel
            ,DOUBLE *restrict dField     ,DOUBLE *restrict cc
            ,DOUBLE const underU         ,const bool sPressure
            ,const short nEn             ,short const nFace    
            ,const short ndm             ,INT const nel)
{ 
  DOUBLE viscosityC ,viscosityV,viscosity;
  DOUBLE densityC   ,densityV  ,density;
/*...*/
  DOUBLE rCell[2],dt;
  DOUBLE p[2],sP;
/*...*/
  DOUBLE v[2],lKsi[2],lNormal[2],lXmcc[2],wf[2],ccV[2];
  DOUBLE dPviz,lModKsi,lModEta,du[2],duDksi[2],lXm[2];
  DOUBLE coef,lAn;
/*...*/
  DOUBLE nk,dfd,dfdc[2],cv,cvc[2],modE,lvSkew[2],nMinusKsi[2];
/*... interpolacao linear*/
  DOUBLE alpha,alphaMenosUm;
  DOUBLE aP,tA[2];
/*... */
  DOUBLE pFace,pf[2],p1,p2;
/*... */
  DOUBLE wfn,velC[2],velV[2],presC,presF;
  DOUBLE gradPresC[2],gradPresV[2];
  DOUBLE gradVelC[2][2],gradVelV[2][2],gf[2][2],gfKsi[2];
  DOUBLE gradVelComp[2][2];
/*...*/
  short iCod=advVel.iCod;
/*...*/
  short idCell = nFace;
  short nAresta,nCarg,typeTime;
  INT vizNel;
  bool fTime;

/*...*/
  dt       = ddt.dt;
  typeTime = ddt.type;
  fTime    = ddt.flag;
  densityC = lDensity[idCell];
/*...................................................................*/

/*... propriedades da celula*/
  viscosityC = MAT2D(idCell,VISCOSITY,prop,MAXPROP);
/*...................................................................*/

/*... | du1/dx1 du1/dx2*/
  gradVelC[0][0] = MAT3D(idCell,0,0,gradVel,2,ndm);
  gradVelC[0][1] = MAT3D(idCell,0,1,gradVel,2,ndm);
/*... | du2/dx1 du2/dx2*/
  gradVelC[1][0] = MAT3D(idCell,1,0,gradVel,2,ndm);
  gradVelC[1][1] = MAT3D(idCell,1,1,gradVel,2,ndm);
/*...................................................................*/

  velC[0]   = MAT2D(idCell,0,vel,ndm);
  velC[1]   = MAT2D(idCell,1,vel,ndm);

  presC         = pres[idCell];
  gradPresC[0]  = MAT2D(idCell,0,gradPres,ndm);
  gradPresC[1]  = MAT2D(idCell,1,gradPres,ndm);
/*...................................................................*/

  p[0]       = 0.0e0;
  p[1]       = 0.0e0;
  pf[0]      = 0.0e0;
  pf[1]      = 0.0e0;
  sP         = 0.0e0;
  lA[idCell] = 0.0e0;
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
      velV[0]        = MAT2D(nAresta,0,vel,ndm);
      velV[1]        = MAT2D(nAresta,1,vel,ndm);
      presF          = pres[nAresta];
      densityV       = lDensity[nAresta];
/*..*/
      lKsi[0]        = MAT2D(nAresta,0,ksi,ndm);
      lKsi[1]        = MAT2D(nAresta,1,ksi,ndm);
      lModKsi        = mKsi[nAresta];
/*..*/
      lvSkew[0]      = MAT2D(nAresta,0,vSkew,ndm);
      lvSkew[1]      = MAT2D(nAresta,1,vSkew,ndm);
/*..*/
      duDksi[0]      = (velV[0] - velC[0]) / lModKsi;
      duDksi[1]      = (velV[1] - velC[1]) / lModKsi;
/*... | du1/dx1 du1/dx2*/
      gradVelV[0][0] = MAT3D(nAresta,0,0,gradVel,2,ndm);
      gradVelV[0][1] = MAT3D(nAresta,0,1,gradVel,2,ndm);
/*... | du2/dx1 du2/dx2*/
      gradVelV[1][0] = MAT3D(nAresta,1,0,gradVel,2,ndm);
      gradVelV[1][1] = MAT3D(nAresta,1,1,gradVel,2,ndm);
/*...................................................................*/
      gradPresV[0] =  MAT2D(nAresta,0,gradPres,ndm);
      gradPresV[1] =  MAT2D(nAresta,1,gradPres,ndm);
      ccV[0]       =  MAT2D(nAresta,0,cc,ndm);
      ccV[1]       =  MAT2D(nAresta,1,cc,ndm);
      lXm[0]       =  MAT2D(nAresta,0,xm,ndm);
      lXm[1]       =  MAT2D(nAresta,1,xm,ndm);
/*...................................................................*/

/*... produtos interno*/
      nk  =    lKsi[0]*lNormal[0] + lKsi[1]*lNormal[1];
/*...................................................................*/
      
/*... correcao sobre-relaxada*/
      modE       = 1.0e0/nk;
/*...................................................................*/

/*...*/
      v[0]         = lvSkew[0] + MAT2D(nAresta,0,xmcc,ndm);
      v[1]         = lvSkew[1] + MAT2D(nAresta,1,xmcc,ndm);
      dPviz        = sqrt(v[0]*v[0] + v[1]*v[1]);
      alpha        = dPviz/lModKsi;
      alphaMenosUm = 1.0e0 - alpha; 
/*...................................................................*/

/*... media harmonica*/
      viscosityV = MAT2D(nAresta,VISCOSITY,prop,MAXPROP); 
      viscosity  = alpha/viscosityC + alphaMenosUm/viscosityV;
      viscosity  = 1.0e0/viscosity;
/*...................................................................*/

/*... difusao direta*/
      coef = viscosity*lModEta;
      dfd  =  coef*modE/lModKsi;
/*...................................................................*/
      
/*...*/
/*... | du1/dx1 du1/dx2 |*/
      gf[0][0] = alphaMenosUm*gradVelC[0][0] + alpha*gradVelV[0][0];
      gf[0][1] = alphaMenosUm*gradVelC[0][1] + alpha*gradVelV[0][1];
/*... | du2/dx1 du2/dx2* |*/
      gf[1][0] = alphaMenosUm*gradVelC[1][0] + alpha*gradVelV[1][0];
      gf[1][1] = alphaMenosUm*gradVelC[1][1] + alpha*gradVelV[1][1];
/*...*/
      wf[0]    =   alphaMenosUm*velC[0] + alpha*velV[0];
      wf[1]    =   alphaMenosUm*velC[1] + alpha*velV[1];
      density  = alphaMenosUm*densityC  + alpha*densityV;
/*...................................................................*/

/*... velocidade normal a face*/
      wfn = wf[0]*lNormal[0] + wf[1]*lNormal[1];
/*...................................................................*/

/*... derivadas direcionais*/
      gfKsi[0]      = gf[0][0]*lKsi[0] + gf[0][1] *lKsi[1];
/*...*/
      gfKsi[1]      = gf[1][0]*lKsi[0] + gf[1][1] *lKsi[1];
/*...................................................................*/

/*... gradiente compacto (Darwish e Moukalled)*/
      du[0]             = duDksi[0]  - gfKsi[0];
      du[1]             = duDksi[0]  - gfKsi[0];
/*...*/
      gradVelComp[0][0] = gf[0][0] + du[0]*lKsi[0];
      gradVelComp[0][1] = gf[0][1] + du[0]*lKsi[1];
/*...*/
      gradVelComp[1][0] = gf[1][0] + du[1]*lKsi[0];
      gradVelComp[1][1] = gf[1][1] + du[1]*lKsi[1];
/*...................................................................*/

/*... derivadas direcionais*/
      nMinusKsi[0] = lNormal[0] - modE*lKsi[0];
      nMinusKsi[1] = lNormal[1] - modE*lKsi[1];
/*...*/
      gfKsi[0] = gradVelComp[0][0]*nMinusKsi[0] 
               + gradVelComp[0][1]*nMinusKsi[1];
/*...*/
      gfKsi[1] = gradVelComp[1][0]*nMinusKsi[0] 
               + gradVelComp[1][1]*nMinusKsi[1];
/*...................................................................*/

/*... correcao nao-ortogonal*/
      dfdc[0] = coef*gfKsi[0];    
      dfdc[1] = coef*gfKsi[1];    
/*...................................................................*/

/*... fluxo convectivo upwind de primeira ordem*/
      cv   = density*wfn*lModEta;
/*...................................................................*/

/*... correcao do fluxo advectivo*/
      if(FOUP == iCod){
        cvc[0] = 0.e0;
        cvc[1] = 0.e0;
      }
/*...................................................................*/

/*...*/
      else{
        cvc[0] = faceBaseTvdV1(velC[0]    ,velV[0]
                              ,gradVelC[0],gradVelV[0]
                              ,lKsi       ,lModKsi   
                              ,cv
                              ,iCod       ,ndm);
/*...*/
        cvc[1] = faceBaseTvdV1(velC[1]    ,velV[1]
                              ,gradVelC[1],gradVelV[1]
                              ,lKsi       ,lModKsi   
                              ,cv
                              ,iCod       ,ndm);
     }
/*...................................................................*/

/*...*/
      lA[nAresta] = dfd - min(cv,0.0e0);
      sP         += cv;    
/*... correcao nao ortogonal e do fluxo advectivo*/        
      p[0]       += dfdc[0] - cv*cvc[0];
      p[1]       += dfdc[1] - cv*cvc[1];
/*...................................................................*/

/*... gradiente da pressao com resconstrucao de segunda ordem
      (forma conservativa)*/
      if(sPressure){
/*.. vetor que une o centroide da celula viznha ao ponto medio 
     da face*/
        v[0]         = lXm[0] - ccV[0];
        v[1]         = lXm[1] - ccV[1];

        p1 = presC + gradPresC[0]*lXmcc[0] + gradPresC[1]*lXmcc[1];
        p2 = presF + gradPresV[0]*v[0]     + gradPresV[1]*v[1];
     
        pFace = 0.5e0*(p1+p2);
      }
/*...................................................................*/

/*... interpolacao linear*/
      else
        pFace = alphaMenosUm*presC + alpha*presF;
/*...................................................................*/
      pf[0]+= pFace*lModEta*lNormal[0];
      pf[1]+= pFace*lModEta*lNormal[1];
/*...................................................................*/

      
    }
/*... contorno*/
    else{
      lA[nAresta] = 0.0e0;
      wfn = velC[0]*lNormal[0] + velC[1]*lNormal[1];
/*...*/
      pFace = presC;
/*... pressao prescrita*/
      if(lFacePresR[nAresta]){
/*...cargas*/
        nCarg = lFacePresL[nAresta]-1;
        pLoadSimplePres(&sP        ,p
                       ,&pFace
                       ,viscosityC      ,densityC
                       ,wfn                 
                       ,lModEta         ,dcca[nAresta]
                       ,loadsPres[nCarg],false); 
      } 
/*...................................................................*/

/*... gradiente da pressao com resconstrucao de segunda ordem*/
      if(sPressure)  
        pFace += gradPresC[0]*lXmcc[0] + gradPresC[1]*lXmcc[1];
      pf[0]+= pFace*lModEta*lNormal[0];
      pf[1]+= pFace*lModEta*lNormal[1];
/*...................................................................*/


/*... velocidades*/
      if(lFaceVelR[nAresta] > 0){
        tA[0] = velC[0];
        tA[1] = velC[1];
/*...cargas*/
        nCarg = lFaceVelL[nAresta]-1;
        pLoadSimple(&sP            ,p
                   ,tA             ,lNormal  
                   ,viscosityC     ,densityC
                   ,lModEta        ,dcca[nAresta]
                   ,loadsVel[nCarg]
                   ,true           ,false);
      }  
/*...................................................................*/

/*... parede impermevavel*/
      else if(lFaceVelR[nAresta] == STATICWALL){
        aP   = viscosityC*lModEta/dcca[nAresta];
        sP  += aP;
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
  
/*... distretização temporal*/
  if(fTime){
/*... EULER*/
    if(typeTime == EULER) 
      sP     += densityC*area[idCell]/dt;
/*...BACKWARD*/
    else if(typeTime == BACKWARD) 
      sP     += 1.5e0*densityC*area[idCell]/dt;
  }
/*...................................................................*/


/*...*/
  lAn = 0.e0;
  if(nFace == 3){
    lAn        = lA[0] + lA[1] + lA[2];
  }
  else if(nFace == 4){
    lAn        = lA[0] + lA[1] + lA[2] + lA[3];
  }
  lA[idCell] = sP + lAn;                   
/*...................................................................*/

/*... under-relaxation(simple)*/
  lA[idCell] = lA[idCell]/underU;
  p[0]      += (1-underU)*lA[idCell]*velC[0];
  p[1]      += (1-underU)*lA[idCell]*velC[1];
/*...................................................................*/

/*... GradP*/
    p[0] -= pf[0];
    p[1] -= pf[1];
/*...................................................................*/

/*...*/
  rCell[0] = 0.0e0;
  rCell[1] = 0.0e0;
  for(nAresta=0;nAresta<nFace;nAresta++){
    if( lViz[nAresta] > -1){
/*... pasando os valores conhecidos para o lado direito*/
      if(lId[nAresta] == -2){
        p[0] += lA[nAresta]*MAT2D(nAresta,0,vel,ndm);
        p[1] += lA[nAresta]*MAT2D(nAresta,0,vel,ndm);
      }
      else{
/*residuo (R = F-KvizUviz ) e valores prescritos por elemento*/
        rCell[0] += lA[nAresta]*MAT2D(nAresta,0,vel,ndm);
        rCell[1] += lA[nAresta]*MAT2D(nAresta,1,vel,ndm);
      } 
    }
  }
/*... residuo: R = F - KpUp*/ 
   rCell[0] += p[0] -lA[idCell]*velC[0];   
   rCell[1] += p[1] -lA[idCell]*velC[1];   
/*...................................................................*/

/*...*/
   if(typeSimple == SIMPLEC) 
    *dField = area[idCell]/(lA[idCell] - lAn);
  else if(typeSimple == SIMPLE) 
    *dField = area[idCell]/lA[idCell]; 
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
  lB[0]     = p[0];
  lB[1]     = p[1];
  lRcell[0] = rCell[0];
  lRcell[1] = rCell[1];
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 01/07/2016                                   *
 * Data de modificaco : 09/07/2016                                   * 
 *-------------------------------------------------------------------* 
 * CELLSIMPLEPRES2D: Celula 2D para equacao de correcao de pressoa   *
 * metodo simple em escoamento imcompressivel                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  * 
 * loadsPres -> definicoes de cargas de pressao                      *  
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
 * dField    -> matriz D do metodo simple                            * 
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
              ,short *restrict lGeomType,DOUBLE *restrict prop
              ,INT *restrict lViz       ,INT *restrict lId  
              ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
              ,DOUBLE *restrict eta     ,DOUBLE *restrict mEta
              ,DOUBLE *restrict normal  ,DOUBLE *restrict area   
              ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
              ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
              ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
              ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
              ,DOUBLE *restrict lRcell  
              ,short  *restrict lFaceVelR  ,short *restrict lFaceVelL
              ,short  *restrict lFacePresR ,short *restrict lFacePresL
              ,DOUBLE *restrict pres    ,DOUBLE *restrict gradPres 
              ,DOUBLE *restrict vel     ,DOUBLE *restrict dField  
              ,const short nEn          ,short const nFace    
              ,const short ndm          ,INT const nel)
{ 

  DOUBLE densityC,densityV ,density;
  DOUBLE dFieldC ,dFieldV  ,dFieldF;
/*...*/
  DOUBLE rCell;
  DOUBLE p,sP;
/*...*/
  DOUBLE v[2],lKsi[2],lNormal[2],gf[2],wf[2],gfp[2];
  DOUBLE dPviz,lModKsi,lModEta;
/*...*/
  DOUBLE gradPresC[2],gradPresV[2];
/*...*/
  DOUBLE nk,dfd,lvSkew[2];
/*... interpolacao linear*/
  DOUBLE alpha,alphaMenosUm;
  DOUBLE tA;
/*... */
  DOUBLE wfn,velC[2],velF[2],presC,presV;
/*...*/
  short idCell = nFace;
  short nAresta,nCarg;
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

  dFieldC   = dField[idCell]; 
  
/*...................................................................*/

  p          = 0.0e0;
  sP         = 0.0e0;
  lA[idCell] = 0.0e0;
  for(nAresta=0;nAresta<nFace;nAresta++){
    vizNel     = lViz[nAresta];
    lNormal[0] = MAT2D(nAresta,0,normal,ndm);
    lNormal[1] = MAT2D(nAresta,1,normal,ndm);
    lModEta    = mEta[nAresta];
    dFieldV    = dField[nAresta];
/*... dominio*/
    if( vizNel  > -1 ){
/*...*/
      velF[0]       = MAT2D(nAresta,0,vel,ndm);
      velF[1]       = MAT2D(nAresta,1,vel,ndm);
      lKsi[0]       = MAT2D(nAresta,0,ksi,ndm);
      lKsi[1]       = MAT2D(nAresta,1,ksi,ndm);
      lModKsi       = mKsi[nAresta];
      lvSkew[0]     = MAT2D(nAresta,0,vSkew,ndm);
      lvSkew[1]     = MAT2D(nAresta,1,vSkew,ndm);
      presV         = pres[nAresta];
      densityV      = lDensity[nAresta];
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
      density  = alphaMenosUm*densityC  + alpha*densityV;
      dFieldF  = alphaMenosUm*dFieldC   + alpha*dFieldV;
/*...................................................................*/

/*... produtos interno*/
      nk  =    lKsi[0]*lNormal[0] + lKsi[1]*lNormal[1];
/*...................................................................*/

/*... difusao direta*/
      dfd = (density*lModEta*dFieldF)/(nk*lModKsi);
/*...................................................................*/

/*... Ferziger-Precic*/
      wf[0]    =   alphaMenosUm*velC[0] + alpha*velF[0];
      wf[1]    =   alphaMenosUm*velC[1] + alpha*velF[1];
/*...................................................................*/

/*... velocidade normal a face*/
      wfn = wf[0]*lNormal[0] + wf[1]*lNormal[1];
/*...................................................................*/
      
/*... interpolacao linear dos gradientes das pressoes*/
      gf[0]    = alphaMenosUm*gradPresC[0] + alpha*gradPresV[0];
      gf[1]    = alphaMenosUm*gradPresC[1] + alpha*gradPresV[1];
      gfp[0]   = presV - presC;
/*...................................................................*/

/*...*/
      gfp[1] = (gf[0]*lKsi[0] + gf[1]*lKsi[1])*lModKsi;  
      gfp[1] = (gfp[1] - gfp[0])/(nk*lModKsi);  
/*...................................................................*/

/*...*/
      wfn += dFieldF*gfp[1];
      wfn *= lModEta;  
/*...................................................................*/

/*...*/
      lA[nAresta] = dfd; 
      p          -=   density*wfn;
/*...................................................................*/

    }
/*... contorno*/
    else{
      lA[nAresta] = 0.0e0;
      wfn = velC[0]*lNormal[0] + velC[1]*lNormal[1];
      if(lFacePresR[nAresta]){
/*...cargas*/
        nCarg = lFacePresL[nAresta]-1;
        pLoadSimplePres(&sP             ,&p
                       ,&tA   
                       ,dFieldV         ,densityC
                       ,wfn                      
                       ,lModEta         ,dcca[nAresta]
                       ,loadsPres[nCarg],true); 
      }
/*... velocidades*/
      if(lFaceVelR[nAresta] > 0){
        nCarg = lFaceVelL[nAresta]-1;
        pLoadSimple(&sP            ,&p
                   ,&tA            ,lNormal
                   ,dFieldV        ,densityC
                   ,lModEta        ,dcca[nAresta]
                   ,loadsVel[nCarg]
                   ,false          ,true);   
      } 
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
//if( nel == 0 ) sP += 1.0e60;
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
