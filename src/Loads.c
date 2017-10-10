#include<Loads.h>

/********************************************************************* 
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 30/09/2017                                   * 
 *-------------------------------------------------------------------* 
 * GETLOADS:                                                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * tA     -> nao definido                                            * 
 * par    -> parametros                                              * 
 * xm     -> ponto (x1,x2,x3)                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * tA     -> valor calculado da funcao                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void getLoads(DOUBLE *par, Loads ld) {


  if (ld.type == DIRICHLETBC) {
    par[0] = ld.par[0];
  }
  
  else if( ld.type == ROBINBC){
    par[0] = ld.par[0];
    par[1] = ld.par[1];
  }
  
  else if (ld.type == NEUMANNBC) {
    par[0] = ld.par[0];
  }

  else if (ld.type == MOVEWALL) {
    par[0]   = ld.par[0];
    par[1]   = ld.par[1];
    par[2]   = ld.par[2];
  }
  else if (ld.type == INLETSTAICTPRES) {

  }
  else if (ld.type == INLET) {
    par[0] = ld.par[0];
    par[1] = ld.par[1];
    par[2] = ld.par[2];
    par[3] = ld.par[3];
  }

  else if (ld.type == OPEN) {
    par[0]   = ld.par[0];
    par[1]   = ld.par[1];
  }
}
/********************************************************************/

/********************************************************************* 
 * LOADSENDPROD :                                                    * 
 * ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * tA     -> nao definido                                            * 
 * par    -> parametros                                              * 
 * xm     -> ponto (x1,x2,x3)                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * tA     -> valor calculado da funcao                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void loadSenProd(DOUBLE *tA,DOUBLE *par,DOUBLE *xm){

  DOUBLE a1,a2,a3,w1,w2,w3,c1,c2,c3;

/*... a1*/
  a1 = par[0];
/*... w1*/
  w1 = par[1];
/*... c1*/
  c1 = par[2];
/*... a2*/
  a2 = par[3];
/*... w2*/
  w2 = par[4];
/*... c2*/
  c2 = par[5];
/*... a3*/
  a3 = par[6];
/*... w3*/
  w3 = par[7];
/*... c3*/
  c3 = par[8];
/*...*/
  tA[0] = a1*sin(w1*xm[0]+c1)
         *a2*sin(w2*xm[1]+c2) 
         *a3*sin(w3*xm[2]+c3);
/*...................................................................*/
} 
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 30/09/2017                                   * 
 *-------------------------------------------------------------------* 
 * pLoadSimple : condicao de contorno para velocidades               *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sP         -> termo da diagonal                                   * 
 * p          -> forca local                                         * 
 * tA         -> velocidade na face                                  * 
 * VelC       -> velocidade na centro do elemento no passo (n-1)     * 
 * n          -> normal                                              * 
 * gradVel    -> dradiente da velocidade                             * 
 * xmcc       -> vetores que unem o centroide aos pontos medios das  * 
 *              faces da celula central                              * 
 * viscosityC -> coeficiente de viscosidade dinamica                 * 
 * effViscosityC -> coeficiente de viscosidade effeyiva              * 
 * densityC   -> massa especifica                                    * 
 * fArea      -> area da face                                        * 
 * dcca       -> menor distancia do centroide central a face desta   *
 *               celula                                              * 
 * ld         -> definicao da carga                                  * 
 * ndm        -> numero de dimensoes                                 * 
 * fCalVel    -> true - atualizada sP e p pela equacao de velocidades* 
 *               false- nao atualizada sP e p                        * 
 * fCalPres   -> true - atualizada sP e p pela da pressao            * 
 *               false- nao atualizada sP e p                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal atualizado                           * 
 * p       -> forcas locais atualizada                               * 
 * tA      -> valor pescrito na face                                 * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void pLoadSimple(DOUBLE *RESTRICT sP  ,DOUBLE *RESTRICT p
          ,DOUBLE *RESTRICT tA        ,DOUBLE *RESTRICT velC
          ,DOUBLE *RESTRICT n       
          ,DOUBLE *RESTRICT gradVel   ,DOUBLE *RESTRICT xmcc
          ,DOUBLE const viscosityC    ,DOUBLE const effViscosityC  
          ,DOUBLE const densityC
          ,DOUBLE const fArea         ,DOUBLE const dcca
          ,Loads ld                   ,short  const ndm 
          ,bool const fCalVel         ,bool const fCalPres
          ,bool const fWallModel      ,short const wallType){

  DOUBLE aP,wfn,wf[3],m,tmp[4],gradVelFace[9],modVel,yPlus,uPlus;
  DOUBLE viscosityWall,densityEnv;

  tmp[0] = tmp[1] = tmp[2] = tmp[3] = 0.e0;
/*... parade impermeavel movel*/
  if( ld.type == MOVEWALL){
    tA[0]   = ld.par[0];
    tA[1]   = ld.par[1];
    if( ndm == 3 )  tA[2]   = ld.par[2];
    if(!fCalVel)return;
/*...*/
    yPlus = 0.e0;
    uPlus = 0.e0;
    if (fWallModel) {
/*... calculo da velociade paralela a face*/
      wfn   = velC[0] * n[0] + velC[1] * n[1];
      wf[0] = velC[0] - wfn * n[0] - tA[0];
      wf[1] = velC[1] - wfn * n[1] - tA[1];
      if( ndm == 3 ) wf[2] = velC[2] - wfn*n[2] - tA[2];
      wfn = wf[0]*wf[0] + wf[1]*wf[1];
      if( ndm == 3 ) wfn += wf[2]*wf[2];
      wfn   = sqrt(wfn);
/*...*/
      if (wfn > 0.e0)
        wallModel(wfn     , viscosityC
                 ,densityC, dcca
                 ,&yPlus  , &uPlus     
                 ,wallType);
/*...................................................................*/
      if( yPlus > 11.81e0)
        viscosityWall = viscosityC*yPlus/uPlus;
      else
        viscosityWall = viscosityC;
    }
/*...................................................................*/ 

/*...*/
    else
      viscosityWall = effViscosityC;
/*...................................................................*/

/*...*/
    aP     = viscosityWall*fArea/dcca;
    if( ndm == 2) {
/*...*/
      sP[0] += aP*(1.e0 - n[0]*n[0]);
      sP[1] += aP*(1.e0 - n[1]*n[1]);
/*...................................................................*/

/*... x*/
      p[0]  += aP*(tA[0]*(1.0-n[0]*n[0]) 
            + (velC[1]-tA[1])*n[1]*n[0]);
/*... y*/
      p[1]  += aP*(tA[1]*(1.0-n[1]*n[1]) 
            + (velC[0]-tA[0])*n[0]*n[1]);
    }
/*...................................................................*/

/*...*/
    else if( ndm == 3 ){
/*...*/
      sP[0] += aP*(1.e0 - n[0] * n[0]);
      sP[1] += aP*(1.e0 - n[1] * n[1]);
      sP[2] += aP*(1.e0 - n[2] * n[2]);
/*...................................................................*/

/*... x*/
      p[0]  += aP*(tA[0]*(1.0-n[0]*n[0]) 
            + (velC[1]-tA[1])*n[1]*n[0]
            + (velC[2]-tA[2])*n[2]*n[0]);
/*... y*/
      p[1]  += aP*(tA[1]*(1.0-n[1]*n[1]) 
            + (velC[0]-tA[0])*n[0]*n[1]
            + (velC[2]-tA[2])*n[2]*n[1]);
/*... z*/
      p[2]  += aP*(tA[2]*(1.0-n[2]*n[2]) 
            + (velC[0]-tA[0])*n[0]*n[2]
            + (velC[1]-tA[1])*n[1]*n[2]);
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... entrada de massa(Pressao estatica)*/
  else if (ld.type == INLETSTAICTPRES) {
    if (ndm == 2) {
/*... direcao*/
      tmp[0] = ld.par[0];
      tmp[1] = ld.par[1];
/*... velocidade normal*/
      wfn    = velC[0]*n[0] + velC[1]*n[1];
      tmp[3] = wfn/(tmp[0]*n[0] + tmp[1]*n[1]);
    }
    else {
/*... direcao*/
      tmp[0] = ld.par[0];
      tmp[1] = ld.par[1];
      tmp[2] = ld.par[2];
/*... velocidade normal*/
      wfn = velC[0]*n[0] + velC[1]*n[1] + velC[2]*n[2];
      tmp[3] = wfn / (tmp[0]*n[0] + tmp[1]*n[1] + tmp[2]*n[2]);
    }
/*...*/
    m = densityC*wfn*fArea;
/*... eq momentum*/
    if (fCalVel) {
      p[0] -= m*tmp[0]*tmp[3];
      p[1] -= m*tmp[1]*tmp[3];
      if (ndm == 3) p[2] -= m*tmp[2]*tmp[3];
    }

/*... eq pressao*/
    if (fCalPres) p[0] -= m;
  }
/*...................................................................*/

/*... entrada de massa ( Velocidade)*/
  else if( ld.type ==  INLET){
    tA[0] = ld.par[1];
    tA[1] = ld.par[2];
    wfn   = tA[0]*n[0] + tA[1]*n[1];
    if( ndm == 3){
      tA[2] = ld.par[3];
      wfn  += tA[2]*n[2];
    }
 /*... eq momentum*/
    densityEnv = ld.par[0];
    m  = densityEnv*wfn*fArea;
    if(fCalVel){
      p[0] -= m*tA[0];
      p[1] -= m*tA[1];
      if(ndm == 3) p[2] -= m*tA[2];
    }
/*... eq pressao*/
    if(fCalPres) p[0] -= m;
  } 
/*...................................................................*/

/*... saida (derivada nula, pressa estatica)*/
  else if( ld.type ==  OUTLET){
    
    wfn = velC[0]*n[0] + velC[1]*n[1];
    if ( ndm == 3 ) wfn += velC[2]*n[2];

    m    = densityC*wfn*fArea;
/*... vb = vc + Grad(Vc)*r*/
    if(fCalVel){
/*...*/
      if( ndm == 2){
        gradVelFace[0] = 0.0e0;
        gradVelFace[1] = 0.0e0;
        gradVelFace[2] = 0.0e0;
        gradVelFace[3] = 0.0e0;
        
        gradFaceNull(gradVelFace,gradVel,xmcc,ndm);
        
        tmp[0] = MAT2D(0,0,gradVelFace,2)*xmcc[0] 
               + MAT2D(0,1,gradVelFace,2)*xmcc[1];      
         
        tmp[1] = MAT2D(1,0,gradVelFace,2)*xmcc[0] 
               + MAT2D(1,1,gradVelFace,2)*xmcc[1];      

        p[0] -= m*tmp[0];
        p[1] -= m*tmp[1];
/*...................................................................*/

/*...*/
        sP[0] += m;
        sP[1] += m;
/*...................................................................*/
      }
/*...................................................................*/

/*...*/
      else{
        gradVelFace[0] = 0.0e0;
        gradVelFace[1] = 0.0e0;
        gradVelFace[2] = 0.0e0;
        gradVelFace[3] = 0.0e0;
        gradVelFace[4] = 0.0e0;
        gradVelFace[5] = 0.0e0;
        gradVelFace[6] = 0.0e0;
        gradVelFace[7] = 0.0e0;
        gradVelFace[8] = 0.0e0;
        
        gradFaceNull(gradVelFace,gradVel,xmcc,ndm);
        
        tmp[0] = MAT2D(0,0,gradVelFace,3)*xmcc[0] 
               + MAT2D(0,1,gradVelFace,3)*xmcc[1]       
               + MAT2D(0,2,gradVelFace,3)*xmcc[2];      
        
        tmp[1] = MAT2D(1,0,gradVelFace,3)*xmcc[0] 
               + MAT2D(1,1,gradVelFace,3)*xmcc[1]       
               + MAT2D(1,2,gradVelFace,3)*xmcc[2];      
        
        tmp[2] = MAT2D(2,0,gradVelFace,3)*xmcc[0] 
               + MAT2D(2,1,gradVelFace,3)*xmcc[1]       
               + MAT2D(2,2,gradVelFace,3)*xmcc[2]; 
     
        p[0] -= m*tmp[0];
        p[1] -= m*tmp[1];
        p[2] -= m*tmp[2];
/*...................................................................*/

/*...*/
        sP[0] += m;
        sP[1] += m;
        sP[2] += m;
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...*/
    if(fCalPres) p[0] -=  m;
/*...................................................................*/
  }   
/*... fronteira aberta(Pressao estatica)*/
  else if (ld.type == OPEN) {

/*...*/
    wfn = velC[0]*n[0] + velC[1]*n[1];
    if( ndm == 3 ) wfn += velC[2]*n[2];     
/*...................................................................*/

/*... saida de massa*/
    if(wfn >= 0.e0){
      m = densityC*wfn*fArea;
/*... vb = vc + Grad(Vc)*r*/
      if(fCalVel){
/*...*/
        if(ndm == 2){
          gradVelFace[0] = 0.0e0; 
          gradVelFace[1] = 0.0e0;
          gradVelFace[2] = 0.0e0;
          gradVelFace[3] = 0.0e0;
        
          gradFaceNull(gradVelFace,gradVel,xmcc,ndm);
        
          tmp[0] = MAT2D(0,0,gradVelFace,2)*xmcc[0] 
                 + MAT2D(0,1,gradVelFace,2)*xmcc[1];      
         
          tmp[1] = MAT2D(1,0,gradVelFace,2)*xmcc[0] 
                 + MAT2D(1,1,gradVelFace,2)*xmcc[1];      

          p[0] -= m*tmp[0];
          p[1] -= m*tmp[1];
/*...................................................................*/

/*...*/
          sP[0] += m;
          sP[1] += m;
/*...................................................................*/
        }
/*...................................................................*/

/*...*/
        else{
          gradVelFace[0] = 0.0e0;
          gradVelFace[1] = 0.0e0;
          gradVelFace[2] = 0.0e0;
          gradVelFace[3] = 0.0e0;
          gradVelFace[4] = 0.0e0;
          gradVelFace[5] = 0.0e0;
          gradVelFace[6] = 0.0e0;
          gradVelFace[7] = 0.0e0;
          gradVelFace[8] = 0.0e0;
        
          gradFaceNull(gradVelFace,gradVel,xmcc,ndm);
        
          tmp[0] = MAT2D(0,0,gradVelFace,3)*xmcc[0] 
                 + MAT2D(0,1,gradVelFace,3)*xmcc[1]       
                 + MAT2D(0,2,gradVelFace,3)*xmcc[2];      
        
          tmp[1] = MAT2D(1,0,gradVelFace,3)*xmcc[0] 
                 + MAT2D(1,1,gradVelFace,3)*xmcc[1]       
                 + MAT2D(1,2,gradVelFace,3)*xmcc[2];      
        
          tmp[2] = MAT2D(2,0,gradVelFace,3)*xmcc[0] 
                 + MAT2D(2,1,gradVelFace,3)*xmcc[1]       
                 + MAT2D(2,2,gradVelFace,3)*xmcc[2]; 
     
          p[0] -= m*tmp[0];
          p[1] -= m*tmp[1];
          p[2] -= m*tmp[2];
/*...................................................................*/

/*...*/
          sP[0] += m;
          sP[1] += m;
          sP[2] += m;
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... entrada*/
    else{
/*...*/
      densityEnv = ld.par[0];
/*...................................................................*/
      m = densityEnv*wfn*fArea;
/*... eq momentum*/
      if (fCalVel) {
/*... direacao definida pelo proprio campo de velocidae*/
        modVel = velC[0]*velC[0] + velC[1]*velC[1];
        if(ndm == 3) modVel += velC[2]*velC[2];
        modVel = 1.e0/sqrt(modVel);
        aP    = m*wfn*modVel;
        p[0] -= aP*velC[0];
        p[1] -= aP*velC[1];
        if (ndm == 3) p[2] -= aP*velC[2];  
/*...................................................................*/

/*... direcao normal*/
/*      aP = m*wfn;
        p[0] -= aP*n[0];
        p[1] -= aP*n[1];
        if (ndm == 3) p[2] -=  aP*n[2];*/
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... eq pressao*/
    if (fCalPres) p[0] -= m;
/*...................................................................*/
  }
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 01/07/2016                                   *
 * Data de modificaco : 20/08/2016                                   * 
 *-------------------------------------------------------------------* 
 * pLoadSimplePres : condicao de contorno para pressao               *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sP         -> termo da diagonal                                   * 
 * p          -> forca local                                         * 
 * tA         -> nao definido                                        * 
 * dField     -> coeficiente de D na face                            * 
 * densityC   -> massa especifica                                    * 
 * wfn        -> velocidade normal a face                            * 
 * xm         -> coordenada do ponto medio da face                   * 
 * fArea      -> area da face                                        * 
 * dd         -> distancia ate o ponto medio da face                 * 
 * ld         -> definicao da carga                                  * 
 * fCal       -> true - atualizada sP e p                            * 
 *               false- atualizada sP e p                            * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal atualizado                           * 
 * p       -> forcas locais atualizada                               * 
 * tA      -> valor pescrito na face                                 * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void pLoadSimplePres(DOUBLE *RESTRICT sP,DOUBLE *RESTRICT p
          ,DOUBLE *RESTRICT tA          ,DOUBLE const df         
          ,DOUBLE const densityC        ,DOUBLE const wfn                                              
          ,DOUBLE const fArea           ,DOUBLE const dd
          ,Loads ld                     ,bool const fCal){

/*... pressao prescrita*/
  if( ld.type == DIRICHLETBC){
    tA[0]   = ld.par[0];
    if(fCal){
      *sP += densityC*df/dd;
    }
  }
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * pLoad : cargas                                                    *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal                                      * 
 * p       -> forca local                                            * 
 * tA      -> nao definido                                           * 
 * coefDifC-> coeficiente de difusao                                 * 
 * densityC-> densidade                                              * 
 * wfn     -> velocidade normal a face                               * 
 * xm      -> coordenada do ponto medio da face                      * 
 * fArea   -> area da face                                           * 
 * dcca    -> menor distancia do centroide central a face desta      *
 *            celula                                                 * 
 * ld      -> definicao da carga                                     * 
 * fCal    -> true - atualizada sP e p                               * 
 *            false- atualizada sP e p                               * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal atualizado                           * 
 * p       -> forcas locais atualizada                               * 
 * tA      -> valor pescrito na face                                 * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void pLoad(DOUBLE *RESTRICT sP  ,DOUBLE *RESTRICT p
          ,DOUBLE *RESTRICT tA
          ,DOUBLE const coefDifC,DOUBLE const densityC
          ,DOUBLE const wfn     ,DOUBLE *RESTRICT xm                   
          ,DOUBLE const fArea   ,DOUBLE const dcca
          ,Loads ld             ,bool const fCal){

  DOUBLE aP,h;


/*... potencial prescrito*/
  if( ld.type == DIRICHLETBC){
    tA[0]   = ld.par[0];
/*...*/
    if(fCal){
      aP   = coefDifC*fArea/dcca;
      *sP += aP;
      *p  += aP*tA[0];
    } 
/*...................................................................*/
  }
/*...................................................................*/

/*... lei de resfriamento de newton*/
  else if( ld.type == ROBINBC){
    h     = ld.par[0];
    tA[0] = ld.par[1];
/*...*/
    if(fCal){
      aP  = ((coefDifC*h)/(coefDifC+h*dcca))*fArea;
      *sP += aP;
      *p  += aP*tA[0];
    }    
/*...................................................................*/
  }
/*...................................................................*/

/*... fluxo prestrito diferente de zero*/
   else if( ld.type == NEUMANNBC){
     tA[0]   = ld.par[0];
/*...*/
     if(fCal)
       *p += fArea*tA[0];
/*...................................................................*/
   }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
   else if( ld.type == SINBC){
     loadSenProd(tA,ld.par,xm); 
/*...*/
     if(fCal){ 
       aP   = coefDifC*fArea/dcca;
       *sP += aP;
       *p  += aP*tA[0];
     }  
/*...................................................................*/
   }
/*...................................................................*/

/*... potencial prescrito (entra)*/
   else if( ld.type == INLET){
     tA[0]   = ld.par[0];
/*...*/
     if(fCal)
       *p -= wfn*densityC*fArea*tA[0];
   }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
   else if( ld.type == OUTLET){
/*...*/
     if(fCal)
       *sP += wfn*densityC*fArea;
/*...................................................................*/

   }
/*...................................................................*/
    
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 21/09/2017                                   *
 * Data de modificaco : 30/09/2017                                   * 
 *-------------------------------------------------------------------* 
 * pLoadEnergy : cargas da equacao de energia                        *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal                                      * 
 * p       -> forca local                                            * 
 * tA      -> nao definido                                           * 
 * coefDifC-> coeficiente de difusao                                 * 
 * densityC-> densidade                                              * 
 * wfn     -> velocidade normal a face                               * 
 * xm      -> coordenada do ponto medio da face                      * 
 * fArea   -> area da face                                           * 
 * dcca    -> menor distancia do centroide central a face desta      *
 *            celula                                                 * 
 * ld      -> definicao da carga                                     * 
 * ldVel   -> definicao da carga para velociades                     * 
 * fCal    -> true - atualizada sP e p                               * 
 *            false- atualizada sP e p                               * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal atualizado                           * 
 * p       -> forcas locais atualizada                               * 
 * tA      -> valor pescrito na face                                 * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void pLoadEnergy(DOUBLE *RESTRICT sP     , DOUBLE *RESTRICT p
               , DOUBLE *RESTRICT tA     , DOUBLE *RESTRICT velC
               , DOUBLE const uC         , DOUBLE *RESTRICT n  
               , DOUBLE const thermCoef  , DOUBLE const densityC
               , DOUBLE const viscosityC , DOUBLE const sHeatC
               , DOUBLE const prT        , DOUBLE *RESTRICT xm                   
               , DOUBLE const fArea      , DOUBLE const dcca
               , Loads ld                , Loads ldVel 
               , short  const ndm
               , bool const fCal         , bool const fTemp
               , bool const iKelvin      , bool const fSheat
               , bool const fWallModel   , short const wallType){

  DOUBLE aP,h,wfn,wf[3],tempPlus,yPlus,uPlus,tC,tW,vW[3],densityEnv;
  DOUBLE prM = viscosityC*sHeatC/thermCoef;

/*...*/
  tempPlus = uPlus = yPlus = 0.e0;
  vW[0] = vW[1] = vW[2] = 0.e0;
  wf[0] = wf[1] = wf[2] = 0.e0;
  if (fWallModel && fCal) {
    if (ldVel.type == MOVEWALL) {
      vW[0]   = ldVel.par[0];
      vW[1]   = ldVel.par[1];
      if( ndm == 3 )  vW[2]   = ldVel.par[2]; 
    }
/*... calculo da velocidade paralela a face*/
    
    wfn   = velC[0] * n[0] + velC[1] * n[1];
    wf[0] = velC[0] - wfn * n[0] - vW[0];
    wf[1] = velC[1] - wfn * n[1] - vW[1];
    if( ndm == 3 ) wf[2] = velC[2] - wfn*n[2] - vW[2];
    wfn = wf[0]*wf[0] + wf[1]*wf[1];
    if( ndm == 3 ) wfn += wf[2]*wf[2];
    wfn   = sqrt(wfn);
/*...*/
    if (wfn > 0.e0){
      wallModel(wfn      , viscosityC
              , densityC , dcca
              , &yPlus   , &uPlus         
              , wallType);
/*...*/ 
      tempPlus = wallModelHeat(yPlus,prM,prT);
    }
/*...................................................................*/ 
  }
/*...................................................................*/ 

/*... potencial prescrito (Parede)*/
  if( ld.type == DIRICHLETBC){
    tA[0]   = ld.par[0];
    if(!fCal) return;
/*... inertial sub-layer*/
    if ( yPlus > 11.81e0 ){
/*...*/
      if (fTemp) {
        tW = tA[0];
        tC = uC;  
      }
      else{
        tW =specificEnthalpyForTemp(tA[0],sHeatC,fSheat,iKelvin);
        tC =specificEnthalpyForTemp(uC   ,sHeatC,fSheat,iKelvin); 
      }
/*...................................................................*/ 

/*...*/
      aP = sHeatC*viscosityC*yPlus*(tW - tC)/(tempPlus*dcca);
      *p += aP*fArea;
/*...................................................................*/ 
    }
/*... viscosity sub-layer*/
    else{
    
/*...*/
      aP   = thermCoef*fArea/dcca;
      if(!fTemp) aP /=  sHeatC;
      *sP += aP;
      *p  += aP*tA[0];
    }
/*...................................................................*/ 
  }
/*...................................................................*/

/*... lei de resfriamento de newton*/
  else if( ld.type == CONVECTIONHEAT){
    tA[0] = ld.par[0];
    if(!fCal) return;

/*...*/
    if (fTemp) {
      tW = tA[0];
      tC = uC;  
    }
    else{
      tW = tA[0];
      tC = specificEnthalpyForTemp(uC   ,sHeatC,fSheat,iKelvin); 
    }
/*...................................................................*/ 

/*... inertial sub-layer*/  
    if ( yPlus > 11.81e0 )
      h  = sHeatC*viscosityC*tempPlus/(tempPlus*dcca);   
    else
      h     = thermCoef/dcca;
/*...................................................................*/
  
    aP  = h*(tW-tC)*fArea;
    *p  += aP;
/*...................................................................*/
  }
/*...................................................................*/

/*... condicao de ROBIN*/
  else if( ld.type == ROBINBC){
    h     = ld.par[0];
    tA[0] = ld.par[1];
  
    aP  = ((thermCoef*h)/(thermCoef+h*dcca))*fArea;
    *sP += aP;
    *p  += aP*tA[0];
/*...................................................................*/
  }
/*...................................................................*/

/*... fluxo prestrito diferente de zero*/
   else if( ld.type == NEUMANNBC){
     tA[0]   = ld.par[0];
     if(!fCal) return;
/*... inertial sub-layer*/     
     if ( yPlus > 11.81e0 ){
/*... energia na forma da temperatura*/
       if (fTemp) 
        tW = uC + tA[0]*dcca*tempPlus/(sHeatC*viscosityC*yPlus);      
/*...................................................................*/
      
/*... energia na forma da entalpia*/
      else{
        tC = specificEnthalpyForTemp(uC ,sHeatC,fSheat,iKelvin); 
        tW = tC + tA[0]*dcca*tempPlus/(sHeatC*viscosityC*yPlus);
        tW = tempForSpecificEnthalpy(tW,sHeatC,fSheat,iKelvin); 
       }
/*...................................................................*/
      
/*...*/
      aP   = thermCoef*fArea/dcca;
      if(!fTemp) aP /=  sHeatC;
      *sP += aP;
      *p  += aP*tW;      
/*...................................................................*/
    }
/*...................................................................*/

/*... viscosity sub-layer*/
    else
      *p += fArea*tA[0];
/*...................................................................*/
   }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
   else if( ld.type == SINBC){
//   loadSenProd(tA,ld.par,xm); 
/*...*/
//   if(fCal){ 
//     aP   = thermCoef*fArea/dcca;
//     *sP += aP;
//     *p  += aP*tA[0];
//   }  
/*...................................................................*/
   }
/*...................................................................*/

/*... potencial prescrito (entra)*/
   else if( ld.type == INLET){
     tA[0]   = ld.par[1];
/*...*/
     if(fCal){
        densityEnv = ld.par[0];
        wfn   = velC[0] * n[0] + velC[1] * n[1];
       *p -= wfn*densityEnv*fArea*tA[0];
     } 
   }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
   else if( ld.type == OUTLET){
/*...*/
     if(fCal){
       wfn   = velC[0] * n[0] + velC[1] * n[1];
       *sP += wfn*densityC*fArea;
     } 
/*...................................................................*/
  }
/*...................................................................*/

/*... fronteira aberta(Pressao estatica)*/
  else if (ld.type == OPEN) {
    tA[0]   = ld.par[1];
    if(!fCal) return;
/*...*/
    wfn  = velC[0]*n[0] + velC[1]*n[1];
    if (ndm == 3) wfn += velC[2]*n[2];
/*... saida de massa*/
    if(wfn >= 0.e0)   
      *sP += wfn*densityC*fArea;
/*...................................................................*/

/*... entrada*/
    else{
      densityEnv = ld.par[0];
      densityEnv = densityC;
      *p -= wfn*densityEnv*fArea*tA[0]; 
    } 
/*...................................................................*/
   }
/*...................................................................*/
}
/*********************************************************************/
