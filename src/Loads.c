#include<Loads.h>
//#include<MyLibs.h>

/********************************************************************* 
 * Data de criacao    : 27/01/2017                                   *
 * Data de modificaco : 27/01/2018                                   * 
 *-------------------------------------------------------------------* 
 * InletFunction:                                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * r      -> posicao                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * tA     -> valor calculado da funcao                               * 
 *-------------------------------------------------------------------* 
 * u(r) = vMax*(1-(r-d)/R)^(n)                                       * 
 * n = 2 laminar                                                     * 
 * n = 7 turbulento                                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE InletFunction(DOUBLE const x) {

  DOUBLE vMax,h,h2,f,r;


  vMax = 1.0;
  h  = 0.1;
  h2 = h*0.5e0;
  r  = (x - h2)/h2;

  f = vMax*(1.0 - pow(fabs(r),7.0));

  return  f;

}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 14/08/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * sen:                                                              * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE functionFace(DOUBLE const x1,DOUBLE const x2, DOUBLE const x3)
{
  return  sin(PI*x1); 
}

/********************************************************************* 
 * Data de criacao    : 14/08/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * sen:                                                              * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 *********************************************************************/
static void functionGeneral(DOUBLE *x,DOUBLE *c, short const iCod)
{
//DOUBLE x1,x2,x3,t,r,R,Vel;

//x1 = x[0];
//x2 = x[1];
//x3 = x[2];
//t  = x[3];
/* 
  c[0] = 0.e0;
  if(x1< 0.e0 && x1>-0.5 && x2 < 0.5)
    c[0] = 1.e0;

  c[1] = 2.e0*x2*(1-x1*x1);
  c[2] =-2.e0*x1*(1-x2*x2);
  c[3] = 0.e0;

  c[4] = 1.e0;
*/

//r = x1*x1 + x2*x2;
//R = 1.e0;
//Vel = 1.e0;

//c[0] = c[1] = 0.e0;
//c[2] = 2.e0*Vel*(1.e0 - r/(R*R));
//c[3] = 1.e0;
}
/********************************************************************/

/********************************************************************* 
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 14/09/2019                                   * 
 *-------------------------------------------------------------------* 
 * GETLOADS:                                                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * tA     -> nao definido                                            * 
 * par    -> parametros                                              * 
 * xm     -> ponto (x1,x2,x3,t)                                      * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * tA     -> valor calculado da funcao                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void getLoads(DOUBLE *par, Loads *ld, DOUBLE *xx)
{

  short i;

  if (ld->type == DIRICHLETBC)
  {
    if (ld->nTypeVar == LVARCONST)
      for(i = 0; i< ld->np; par[i] = ld->par[i], i++);
    else if(ld->nTypeVar == LFUNC)
       functionGeneral(xx,par,0);
  } 
  else if( ld->type == ROBINBC)
  {
    par[0] = ld->par[0];
    par[1] = ld->par[1];
  }
  
  else if (ld->type == NEUMANNBC)
    for(i = 0; i< ld->np; par[i] = ld->par[i], i++);

  else if (ld->type == MOVEWALL) 
  {
    for (i = 0; i< 4; par[i] = ld->par[i], i++);
  }

  else if (ld->type == INLETSTAICTPRES) 
  {
  }

  else if (ld->type == INLETTOTALPRES) 
  {
    for (i = 0; i< 5; par[i] = ld->par[i], i++);
  }

  else if (ld->type == INLET)
  {
    if (ld->nTypeVar == LVARCONST)
      for(i = 0; i< ld->np; par[i] = ld->par[i], i++);
    else if(ld->nTypeVar == LFUNC)
      functionGeneral(xx,par,0);
  }

  else if (ld->type == OPEN)
  {
    for(i = 0; i< ld->np; par[i] = ld->par[i], i++);
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
 * Data de modificaco : 06/10/2019                                   * 
 *-------------------------------------------------------------------* 
 * pLoadSimple : condicao de contorno para velocidades               *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sP         -> termo da diagonal                                   * 
 * p          -> forca local                                         * 
 * tA         -> velocidade na face                                  * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * VelC       -> velocidade na centro do elemento no passo (n-1)     * 
 * gradVel    -> dradiente da velocidade                             * 
 * presC      -> pressao no centro da celula                         * 
 * gradPresC  -> gradiente no centro da celulca                      * 
 * viscosityC -> coeficiente de viscosidade dinamica                 * 
 * effViscosityC -> coeficiente de viscosidade effeyiva              * 
 * xx         -> pontos medios da face e tempo (x1,x2,x3,t)
 * s          -> vetor s                                             * 
 * e          -> vetor e                                             *
 * t          -> vetor t                                             *  
 * n          -> vetorno normal a face                               * 
 * densityC   -> massa especifica                                    * 
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  *  
 * fArea      -> area da face                                        * 
 * dcca       -> menor distancia do centroide central a face desta   *
 *               celula                                              * 
 * ld         -> definicao da carga                                  * 
 * ndm        -> numero de dimensoes                                 * 
 * fCalVel    -> true - atualizada sP e p pela equacao de velocidades* 
 *               false- nao atualizada sP e p                        * 
 * fCalPres   -> true - atualizada sP e p pela da pressao            * 
 *               false- nao atualizada sP e p                        * 
 * nEl        -> numero do elemento                                  *
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
void pLoadSimple(DOUBLE *RESTRICT sP, DOUBLE *RESTRICT p
          , DOUBLE *RESTRICT tA     , DOUBLE *RESTRICT xmcc               
          , DOUBLE *RESTRICT velC   , DOUBLE *RESTRICT gradVel
          , DOUBLE const presC      , DOUBLE *RESTRICT gradPresC
          , DOUBLE const viscosityC , DOUBLE const effViscosityC  
          , DOUBLE *RESTRICT xx
          , DOUBLE *RESTRICT s      , DOUBLE *RESTRICT e        
          , DOUBLE *RESTRICT t      , DOUBLE *RESTRICT n  
          , DOUBLE const densityC   , DOUBLE *RESTRICT wallPar   
          , DOUBLE const fArea      , DOUBLE const dcca
          , Loads *ld               , short  const ndm     
          , bool const fCalVel      , bool const fCalPres
          , bool const fDiv         , bool const fWallModel
          , short const wallType    , INT const nEl )
{

  DOUBLE aP,wfn,m,tmp[5],gradVelFace[9],modVel,yPlus,uPlus,lambda;
  DOUBLE viscosityWall,densityEnv,par[MAXLOADPARAMETER],ev[3],ss[6];

  lambda = 0.0e+00;
  if(fDiv)
    lambda = -2.e0/3.e0*viscosityC;
  tmp[0] = tmp[1] = tmp[2] = tmp[3] = 0.e0;
    
/*... parade impermeavel movel*/
  if( ld->type == MOVEWALL)
  {
    getLoads(par,ld,xx);
    tA[0] = par[0];
    tA[1] = par[1];
    if(ndm == 3) tA[2] = par[2];
    if(!fCalVel)return;
/*...*/
    yPlus = 0.e0;
    uPlus = 1.e0;
    if (fWallModel) {    
      yPlus = wallPar[0];
      uPlus = wallPar[1];
/*...................................................................*/
      if(yPlus > 0.01e0)          
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
    moveWall(velC ,tA 
           ,n
           ,sP    ,p
           ,dcca  ,viscosityWall
           ,fArea);
/*...................................................................*/
  }
/*...................................................................*/

/*... parade impermeavel movel*/
  else if( ld->type == STATICWALL)
  {
    tA[0] = 0.e0;
    tA[1] = 0.e0;
    if(ndm == 3) tA[2] = 0.e0;
    if(!fCalVel)return;
/*...*/
    yPlus = 0.e0;
    uPlus = 1.e0;
    if (fWallModel) {    
      yPlus = wallPar[0];
      uPlus = wallPar[1];
/*...................................................................*/
      if(yPlus > 0.01e0)          
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
    moveWall(velC ,tA 
           ,n
           ,sP    ,p
           ,dcca  ,viscosityWall
           ,fArea);
/*...................................................................*/
  }
/*...................................................................*/

/*... entrada de massa(Pressao estatica)*/
  else if (ld->type == INLETSTAICTPRES)
  {
    getLoads(par,ld,xx);
    densityEnv = par[ndm];
/*... direcao*/
    tmp[0] = par[0];
    tmp[1] = par[1];
/*... velocidade normal*/
    if (ndm == 3) {
      tmp[2] = par[2];
      wfn = velC[0]*n[0] + velC[1]*n[1] + velC[2]*n[2];
      tmp[3] = wfn / (tmp[0]*n[0] + tmp[1]*n[1] + tmp[2]*n[2]);
    }
    else {
      wfn    = velC[0]*n[0] + velC[1]*n[1];
      tmp[3] = wfn/(tmp[0]*n[0] + tmp[1]*n[1]);
    }
/*...................................................................*/    

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
  else if( ld->type ==  INLET)
  {
    getLoads(par,ld, xx);
    tA[0] = par[0];
    tA[1] = par[1];
    wfn   = tA[0]*n[0] + tA[1]*n[1];
    if( ndm == 3){
      tA[2] = par[2];
      wfn  += tA[2]*n[2];
    }
    
 /*... eq momentum*/
    densityEnv = par[ld->np-1];
    m  = densityEnv*wfn*fArea;
    if(fCalVel)
    {
/*....*/
      stress(ss        ,gradVel
            ,viscosityC,lambda
            ,3);
      p[0] += ss[0]*s[0] + ss[3]*s[1] + ss[5]*s[2]; 
      p[1] += ss[3]*s[0] + ss[1]*s[1] + ss[4]*s[2]; 
      p[2] += ss[5]*s[0] + ss[4]*s[1] + ss[2]*s[2]; 
/*....................................................................*/
      p[0] -= m*tA[0];
      p[1] -= m*tA[1];
      if(ndm == 3) p[2] -= m*tA[2];
    }
/*... eq pressao*/
    if(fCalPres) p[0] -= m;
  } 
/*...................................................................*/

/*... saida (derivada nula, pressa estatica)*/
  else if( ld->type ==  OUTLET)
  {    
    wfn = velC[0]*n[0] + velC[1]*n[1];
    if ( ndm == 3 ) wfn += velC[2]*n[2];

    m    = densityC*wfn*fArea;
/*... vb = vc + Grad(Vc)*r*/
    if(fCalVel)
    {
/*....*/
      stress(ss        ,gradVel
            ,viscosityC,lambda
            ,3);
      p[0] += ss[0]*s[0] + ss[3]*s[1] + ss[5]*s[2];  
      p[1] += ss[3]*s[0] + ss[1]*s[1] + ss[4]*s[2]; 
      p[2] += ss[5]*s[0] + ss[4]*s[1] + ss[2]*s[2];   
/*....................................................................*/

/*...*/
      if( ndm == 2)
      {
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
      else
      {
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
  else if (ld->type == OPEN) 
  {
/*...*/
    wfn = velC[0]*n[0] + velC[1]*n[1];
    if( ndm == 3 ) wfn += velC[2]*n[2];     
/*...................................................................*/

/*... saida de massa*/
    if(wfn >= 0.e0)
    {
      m = densityC*wfn*fArea;
/*... vb = vc + Grad(Vc)*r*/
      if(fCalVel)
      {
/*...*/
        if(ndm == 2)
        {
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
        else
        {
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
    else
    {
/*...*/
      densityEnv = ld->density;
/*...................................................................*/
      m = densityEnv*wfn*fArea;
/*... eq momentum*/
      if (fCalVel) 
      {
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
//      aP = m*wfn;
//      p[0] -= aP*n[0];
//      p[1] -= aP*n[1];
//      if (ndm == 3) p[2] -=  aP*n[2];  
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

/*... entrada de massa(Pressao total)*/
  else if (ld->type == INLETTOTALPRES) 
  {
    getLoads(par,ld,xx);
    densityEnv = par[0];   
    ev[0]      = par[2];
    ev[1]      = par[3];
    ev[2]      = par[4];

/*...*/
    if(ndm == 3)
    {
      tmp[0] = velC[0]*n[0] + velC[1]*n[1] + velC[2]*n[2];
      tmp[1] = (ev[0]*n[0] + ev[1]*n[1] + ev[2]*n[2])*fArea;
   }
/*...................................................................*/

/*...*/
    else 
    {
      tmp[0] = velC[0]*n[0] + velC[1]*n[1];
      tmp[1] = (ev[0]*n[0] + ev[1]*n[1])*fArea;
    }
/*...................................................................*/
    
/*...*/
    m      = densityEnv*tmp[0]*fArea;  
    modVel = m/(densityEnv*tmp[1]);
/*...................................................................*/

/*... eq momentum*/
    if (fCalVel) 
    {
      tmp[0] = m*modVel;
      p[0] -= tmp[0]*ev[0];
      p[1] -= tmp[0]*ev[1];
      if (ndm == 3) p[2] -= tmp[0]*ev[2];
    }
/*... eq pressao*/
    if (fCalPres) p[0] -= m;
/*...................................................................*/
  }
/*...................................................................*/


}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 01/07/2016                                   *
 * Data de modificaco : 16/01/2020                                   * 
 *-------------------------------------------------------------------* 
 * pLoadSimplePres : condicao de contorno para pressao               *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sP         -> termo da diagonal                                   * 
 * p          -> forca local                                         * 
 * tA         -> nao definido                                        *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                * 
 * presC      -> pressao no centro da celula                         * 
 * gradPresC  -> gradiente no centro da celulca                      * 
 * sl         -> vetor s                                             * 
 * e          -> vetor e                                             *
 * t          -> vetor t                                             *  
 * n          -> vetorno normal a face                               * 
 * densityC   -> massa especifica                                    * 
 * velC       -> velocidade da celucla central                       * 
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
 * termo difusivo anisotropico                                       *
 * sl = e + t, sl = D*S                                              *
 *********************************************************************/
void pLoadSimplePres(DOUBLE *RESTRICT sP, DOUBLE *RESTRICT p
          , DOUBLE *RESTRICT tA         , DOUBLE *RESTRICT xmcc
          , DOUBLE const presC          , DOUBLE *RESTRICT gradPresC 
          , DOUBLE *RESTRICT s          , DOUBLE *RESTRICT e        
          , DOUBLE *RESTRICT t          , DOUBLE *RESTRICT n
          , DOUBLE *RESTRICT g          , DOUBLE *RESTRICT gradRho      
          , DOUBLE *RESTRICT velC       , DOUBLE const gh
          , DOUBLE const densityC       , DOUBLE const densityRef
          , DOUBLE const fArea          , DOUBLE const dd
          , Loads *ld                   , short  const ndm 
          , short const iCodPres        , bool const fBuoyant
          , bool const fCal)
{              

  DOUBLE modVel,par[MAXLOADPARAMETER],ev[3],dc,densityEnv,m,presT;
  DOUBLE tmp[5],modE,xx[4],gradPb[3];

/*...*/
  tmp[0] = e[0]*e[0] + e[1]*e[1];
  if(ndm == 3) tmp[0] += e[2]*e[2];
  modE = sqrt(tmp[0]);
/*...................................................................*/

/*... pressao prescrita*/
  if( ld->type == DIRICHLETBC)
  {
    tA[0]   = ld->par[0];
    if(fCal){
      *sP += densityC*modE/dd;
    }
  }
/*...................................................................*/

/*... pressao prescrita*/
  else if( ld->type == FLUXPRES)
  {
    if(fBuoyant)
    {
      gradPbBuoyant(gradPb  , gradRho
                  , g       , gh  
                  , densityC, densityRef
                  , iCodPres);
/*...................................................................*/

/*...*/
       tA[0] += (gradPb[0]*n[0]+gradPb[1]*n[1]+gradPb[2]*n[2])*dd;
/*...................................................................*/

    }
/*...................................................................*/
  }
/*...................................................................*/

/*... pressao total prescrita*/
  else if( ld->type == INLETTOTALPRES){
    getLoads(par,ld,xx);
    densityEnv = par[0];   
    presT      = par[1];
    ev[0]      = par[2];
    ev[1]      = par[3];
    ev[2]      = par[4];

    dc = modE/dd;

/*...*/
    if(ndm == 3){
/*...*/
      tmp[0] = velC[0]*n[0] + velC[1]*n[1] + velC[2]*n[2];
      tmp[1] = (ev[0]*n[0] + ev[1]*n[1] + ev[2]*n[2])*fArea;
    }
/*...................................................................*/

/*...*/
    else {
      tmp[0] = velC[0]*n[0] + velC[1]*n[1];
      tmp[1] = (ev[0]*n[0] + ev[1]*n[1])*fArea;
    }
/*...................................................................*/

    m      = densityEnv*tmp[0]*fArea;
    modVel = m/(densityEnv*tmp[1]);
/*... pressa estatica*/
    tA[0] = presT - 0.5*densityEnv*modVel*modVel;

    if(fCal){    
      tmp[0] = densityEnv*modVel;
      tmp[1] = m - dc*tmp[0]*tmp[0];
      *sP += densityEnv*m*dc/tmp[1];
//    *sP += densityC*modE/dd;
    }
  }
/*...................................................................*/

/*... pressao prescrita*/
  else if( ld->type == OPEN)
  {
    getLoads(par,ld,xx);
    tmp[0] = velC[0]*n[0] + velC[1]*n[1];
    if(ndm==3)  tmp[0] += velC[2]*n[2]; 

/*... saida*/
    if(tmp[0] > 0.e0)
    {
      tA[0]   = par[0];
      if(fCal)
      {
        *sP += densityC*modE/dd;
      }
      if(fCal){
      *sP += densityC*modE/dd;
      }
    }
/*...................................................................*/

  }
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 02/06/2018                                   * 
 *-------------------------------------------------------------------* 
 * pLoad : cargas                                                    *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal                                      * 
 * p       -> forca local                                            * 
 * tA      -> nao definido                                           * 
 * velC    -> Velocidade da celula central                           *
 * n       -> normal externa                                         *
 * coefDifC-> coeficiente de difusao                                 * 
 * densityC-> densidade                                              * 
 * xm      -> coordenada do ponto medio da face                      * 
 * fArea   -> area da face                                           * 
 * dcca    -> menor distancia do centroide central a face desta      *
 *            celula                                                 * 
 * ld      -> definicao da carga                                     *
 * ndm     -> dimensao espacial                                      *
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
          ,DOUBLE *RESTRICT tA  ,DOUBLE *RESTRICT velC
          ,DOUBLE *RESTRICT n
          ,DOUBLE const coefDifC,DOUBLE const densityC
          ,DOUBLE *RESTRICT xm  ,DOUBLE const fArea   
          ,DOUBLE const dcca    ,Loads *ld    
          ,short const ndm      ,bool const fCal){

  DOUBLE aP,h,par[MAXLOADPARAMETER], velB[3], densityEnv, wfn;

  velB[0] = velB[1] = velB[2] = 0.e0;

/*... potencial prescrito*/
  if( ld->type == DIRICHLETBC){
    getLoads(par, ld, xm);
    tA[0]   = par[0];
/*...*/
    if(fCal)
    {
      aP   = coefDifC*fArea/dcca;
      *sP += aP;
      *p  += aP*tA[0];
    } 
/*...................................................................*/
  }
/*...................................................................*/

/*... lei de resfriamento de newton*/
  else if( ld->type == ROBINBC){
    getLoads(par, ld, xm);
    h     = par[0];
    tA[0] = par[1];
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
   else if( ld->type == NEUMANNBC){
     getLoads(par, ld, xm);
     tA[0] = par[0];
/*...*/
     if(fCal)
       *p += fArea*tA[0];
/*...................................................................*/
   }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
   else if( ld->type == SINBC){
     getLoads(par, ld, xm);
     loadSenProd(tA,par,xm); 
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
   else if( ld->type == INLET)
   {
     getLoads(par, ld, xm);
     tA[0] = par[0];

/*...*/
     if(fCal)
     {
       densityEnv = par[ndm+1];
       velB[0]    = par[1];
       velB[1]    = par[2];
       wfn = velB[0] * n[0] + velB[1] * n[1];
       if (ndm == 3)
       {
         velB[2] = par[3];
         wfn += velB[2] * n[2];
       }
/*... termo advectivo*/
       *p -= wfn*densityEnv*fArea*tA[0];
/*... termo difusivo*/
        aP   = coefDifC*fArea/dcca;
        *sP += aP;
        *p  += aP*tA[0];
/*...................................................................*/

     }
   }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
   else if( ld->type == OUTLET){
/*...*/
     if(fCal){
       wfn = velC[0] * n[0] + velC[1] * n[1];
       if (ndm == 3) wfn += velC[2] * n[2];
       *sP += wfn * densityC*fArea;
     }
/*...................................................................*/
   }
/*...................................................................*/
    
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 21/09/2017                                   *
 * Data de modificaco : 15/07/2018                                   * 
 *-------------------------------------------------------------------* 
 * pLoadEnergy : cargas da equacao de energia                        *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * vProp   -> estrutura que guarda a propriedades do fluido          *
 * sP      -> termo da diagonal                                      * 
 * p       -> forca local                                            * 
 * tA      -> nao definido                                           * 
 * n       -> normal externa                                         *
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
void pLoadEnergy(PropVarFluid *vProp
               , DOUBLE *RESTRICT sP     , DOUBLE *RESTRICT p
               , DOUBLE *RESTRICT tA     , DOUBLE *RESTRICT velC
               , DOUBLE const uC         , DOUBLE *RESTRICT n  
               , DOUBLE const thermCoef  , DOUBLE const densityC
               , DOUBLE const viscosityC , DOUBLE const sHeatC
               , DOUBLE const prT        , DOUBLE *RESTRICT xx                   
               , DOUBLE const fArea      , DOUBLE const dcca
               , Loads *ld               
               , DOUBLE *RESTRICT wallPar, short  const ndm          
               , bool const fCal         , bool const fTemp
               , bool const iKelvin      , bool const fSheat
               , bool const fWallModel   , short const wallType){

  DOUBLE aP,h,wfn,wf[3],tempPlus,yPlus,uPlus,tC,tW,densityEnv;
  DOUBLE prM = viscosityC*sHeatC/thermCoef,diff;
  DOUBLE par[MAXLOADPARAMETER],velB[3];
  
  velB[0] = velB[1] = velB[2] = 0.e0;

/*...*/
  tempPlus = uPlus = yPlus = 0.e0;
  wf[0] = wf[1] = wf[2] = 0.e0;
  if (fWallModel && fCal) {
/*... calculo da velocidade paralela a face*/
      yPlus = wallPar[0];
/*...*/ 
      tempPlus = wallModelHeat(yPlus,prM,prT);
/*...................................................................*/ 
  }
/*...................................................................*/ 

/*... potencial prescrito (Parede)*/
  if( ld->type == DIRICHLETBC)
  {
    getLoads(par,ld,xx);
    tA[0]   = par[0];
    if(!fCal) return;
/*... inertial sub-layer*/
    if ( yPlus > 11.81e0 )
    {
/*...*/
      if (fTemp)
      {
        tW = tA[0];
        tC = uC;  
      }
      else
      {
        tW =specificEnthalpyForTemp(&vProp->sHeat ,25.0
                                   ,tA[0]         ,sHeatC
                                   ,fSheat        ,iKelvin);

        tC =specificEnthalpyForTemp(&vProp->sHeat ,25.0
                                   ,uC            ,sHeatC
                                   ,fSheat        ,iKelvin); 
      }
/*...................................................................*/ 

/*...*/
      aP = sHeatC*viscosityC*yPlus*(tW - tC)/(tempPlus*dcca);
      *p += aP*fArea;
/*...................................................................*/ 
    }
/*... viscosity sub-layer*/
    else
    {    
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
  else if( ld->type == CONVECTIONHEAT){
    getLoads(par,ld,xx);
    tA[0] = par[0];    
    if(!fCal) return;

/*...*/
    if (fTemp) {
      tW   = tA[0];
      tC   = uC; 
      diff =  thermCoef; 
    }
    else{
//    tW = tA[0];
      tW = tempForSpecificEnthalpy(&vProp->sHeat
                                  ,tA[0],sHeatC,fSheat,iKelvin);
//    tC = specificEnthalpyForTemp(uC   ,sHeatC,fSheat,iKelvin); 
      diff =  thermCoef / sHeatC; 
    }
/*...................................................................*/ 

/*... inertial sub-layer*/  
    if ( yPlus > 11.81e0 )
      h  = sHeatC*viscosityC*tempPlus/(tempPlus*dcca);   
    else
      h  = thermCoef/dcca;
/*...................................................................*/
  
//  aP  = h*(tW-tC)*fArea;
//  *p  += aP;
    aP  = ((diff*h)/(diff+h*dcca))*fArea;
    *sP += aP;
    *p  += aP*tW;
/*...................................................................*/
  }
/*...................................................................*/

/*... condicao de ROBIN*/
  else if( ld->type == ROBINBC){
    getLoads(par,ld,xx); 
    h     = par[0];
    tA[0] = par[1];
    aP  = ((thermCoef*h)/(thermCoef+h*dcca))*fArea;
    *sP += aP;
    *p  += aP*tA[0];
/*...................................................................*/
  }
/*...................................................................*/

/*... fluxo prestrito diferente de zero*/
   else if( ld->type == NEUMANNBC){
     getLoads(par,ld,xx); 
     tA[0] = par[0];
     if(!fCal) return;
/*... inertial sub-layer*/     
     if ( yPlus > 11.81e0 ){
/*... energia na forma da temperatura*/
       if (fTemp) 
        tW = uC + tA[0]*dcca*tempPlus/(sHeatC*viscosityC*yPlus);      
/*...................................................................*/
      
/*... energia na forma da entalpia*/
       else
       {
         tC = specificEnthalpyForTemp(&vProp->sHeat,298.15
                                    ,uC           ,sHeatC
                                    ,fSheat       ,iKelvin); 
         tW = tC + tA[0]*dcca*tempPlus/(sHeatC*viscosityC*yPlus);
         tW = tempForSpecificEnthalpy(&vProp->sHeat
                                    ,tW           ,sHeatC
                                    ,fSheat       ,iKelvin); 
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
   else if( ld->type == SINBC){
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
   else if( ld->type == INLET){
     getLoads(par,ld,xx); 
     tA[0]   = par[0];
/*...*/
     if(fCal){
        densityEnv = ld->density;
        velB[0]    = ld->vel[0];
        velB[1]    = ld->vel[1];
        
        wfn   = velB[0] * n[0] + velB[1] * n[1];
        if (ndm == 3) {
          velB[2] = ld->vel[2];  
          wfn += velB[2]*n[2];
        }
       *p -= wfn*densityEnv*fArea*tA[0];
     } 
   }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
   else if( ld->type == OUTLET){
/*...*/
     if(fCal){
       wfn   = velC[0] * n[0] + velC[1] * n[1];
       if (ndm == 3) wfn += velC[2]*n[2];
       *sP += wfn*densityC*fArea;
     } 
/*...................................................................*/
  }
/*...................................................................*/

/*... fronteira aberta(Pressao estatica)*/
  else if (ld->type == OPEN)
  {
    getLoads(par,ld,xx); 
    if(!fCal) return;
/*...*/
    wfn  = velC[0]*n[0] + velC[1]*n[1];
    if (ndm == 3) wfn += velC[2]*n[2];
/*... saida de massa*/
    if(wfn >= 0.e0)   
      *sP += wfn*densityC*fArea;
/*...................................................................*/

/*... entrada*/
    else
    {
      tA[0]      = par[0];
      densityEnv = ld->density;
      *p -= wfn*densityEnv*fArea*tA[0]; 
    } 
/*...................................................................*/
   }
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 21/09/2017                                   *
 * Data de modificaco : 10/05/2019                                   * 
 *-------------------------------------------------------------------* 
 * pLoadCombustion : cargas das equacoes de combustao                *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * vProp   -> estrutura que guarda a propriedades do fluido          *
 * sP      -> termo da diagonal                                      * 
 * p       -> forca local                                            * 
 * tA      -> nao definido                                           * 
 * n       -> normal externa                                         *
 * coefDif -> coeficiente de difusao                                 * 
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
void pLoadCombustion(PropVarFluid *vProp
               , DOUBLE *RESTRICT sP      , DOUBLE *RESTRICT p
               , DOUBLE *RESTRICT tA      , DOUBLE *RESTRICT velC
               , DOUBLE *RESTRICT uC      , DOUBLE *RESTRICT n  
               , DOUBLE *RESTRICT diffCoef, DOUBLE const densityC
               , DOUBLE const viscosityC 
               , DOUBLE const prT         , DOUBLE *RESTRICT xx                   
               , DOUBLE const fArea       , DOUBLE const dcca
               , Loads *ld                
               , DOUBLE *RESTRICT wallPar , short  const ndm          
               , bool const fCal          , bool const fWallModel  
               , short const nComb        , short const wallType)
{

  short i;
  DOUBLE aP,wfn,wf[3],yPlus,uPlus,densityEnv;
  DOUBLE par[MAXLOADPARAMETER],velB[3];
  
  velB[0] = velB[1] = velB[2] = 0.e0;

/*...*/
//tempPlus = uPlus = yPlus = 0.e0;
  wf[0] = wf[1] = wf[2] = 0.e0;
//if (fWallModel && fCal) {
/*... calculo da velocidade paralela a face*/
//    yPlus = wallPar[0];
/*...*/ 
//    tempPlus = wallModelHeat(yPlus,prM,prT);
/*...................................................................*/ 
//}
/*...................................................................*/ 

/*... potencial prescrito (Parede)*/
  if( ld->type == DIRICHLETBC)
  {
    getLoads(par,ld,xx);
      
    for(i=0;i<nComb;i++)
      tA[i]   = par[i];
    
    if(!fCal) return;

    for(i=0;i<nComb;i++)
    {
      aP     = diffCoef[i]*fArea/dcca;
      sP[i] += aP;
      p[i]  += aP*tA[i];
    }
/*...................................................................*/ 
  }
/*...................................................................*/

/*... lei de resfriamento de newton*/
  else if( ld->type == CONVECTIONHEAT);
/*...................................................................*/

/*... condicao de ROBIN*/
  else if( ld->type == ROBINBC);
/*...................................................................*/

/*... fluxo prestrito diferente de zero*/
   else if( ld->type == NEUMANNBC)
   {
     getLoads(par,ld,xx); 

    for(i=0;i<nComb;i++)
      tA[i]   = par[i];

     if(!fCal) return;

     for(i=0;i<nComb;i++)
       p[i] += fArea*tA[i];
/*...................................................................*/
   }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
   else if( ld->type == SINBC);
/*...................................................................*/

/*... potencial prescrito (entra)*/
   else if( ld->type == INLET){
     getLoads(par,ld,xx); 
     for(i=0;i<nComb;i++)
       tA[i]   = par[i];

/*...*/
     if(fCal){
        densityEnv = ld->density;
        velB[0]    = ld->vel[0];
        velB[1]    = ld->vel[1];
        
        wfn   = velB[0] * n[0] + velB[1] * n[1];
        if (ndm == 3) {
          velB[2] = ld->vel[2];
          wfn += velB[2]*n[2];
        }
        for(i=0;i<nComb;i++)
          p[i] -= wfn*densityEnv*fArea*tA[i];
     } 
   }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
   else if( ld->type == OUTLET){
/*...*/
     if(fCal){
       wfn   = velC[0] * n[0] + velC[1] * n[1];
       if (ndm == 3) wfn += velC[2]*n[2];
       for(i=0;i<nComb;i++)
         sP[i] += wfn*densityC*fArea;
     } 
/*...................................................................*/
  }
/*...................................................................*/

/*... fronteira aberta(Pressao estatica)*/
  else if (ld->type == OPEN)
  {
    getLoads(par,ld,xx); 
/*...*/
    wfn  = velC[0]*n[0] + velC[1]*n[1];
    if (ndm == 3) wfn += velC[2]*n[2];
/*... saida de massa*/
    if(wfn >= 0.e0)
    {
      if(!fCal) return;
      for(i=0;i<nComb;i++)
        sP[i] += wfn*densityC*fArea;
    }
/*...................................................................*/

/*... entrada*/
    else
    {
      densityEnv = ld->density;
      for(i=0;i<nComb;i++)
      {
        tA[i] = par[i];
        if(fCal) 
          p[i] -= wfn*densityEnv*fArea*tA[i];
      }
      
    } 
/*...................................................................*/
   }
/*...................................................................*/
/*...................................................................*/
}
/*********************************************************************/


/********************************************************************* 
 * Data de criacao    : 21/09/2017                                   *
 * Data de modificaco : 11/05/2019                                   * 
 *-------------------------------------------------------------------* 
 * pLoadOneEqK : cargas da equacao K                                 *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sP      -> termo da diagonal                                      * 
 * p       -> forca local                                            * 
 * tA      -> nao definido                                           * 
 * velC    -> Velocidade da celula central                           *
 * n       -> normal externa                                         *
 * coefDifC-> coeficiente de difusao                                 * 
 * densityC-> densidade                                              * 
 * xm      -> coordenada do ponto medio da face                      * 
 * fArea   -> area da face                                           * 
 * dd      -> menor distancia do centroide central a face desta      *
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
void pLoadOneEqK(DOUBLE *RESTRICT sP     , DOUBLE *RESTRICT p
               , DOUBLE *RESTRICT tA     , DOUBLE *RESTRICT velC
               , DOUBLE const uC         , DOUBLE *RESTRICT n  
               , DOUBLE const densityC   , DOUBLE const viscosityC
               , DOUBLE const prT        , DOUBLE *RESTRICT xm                   
               , DOUBLE const fArea      , DOUBLE const dd
               , Loads *ld               , Loads *ldVel
               , short const lFaceReK  
               , DOUBLE *RESTRICT wallPar, short  const ndm          
               , bool const fCal         , bool const fWallModel   
               , short const wallType){

  short i;
  DOUBLE aP,wfn,densityEnv,in,vv,vB[3];
  DOUBLE par[MAXLOADPARAMETER],xx[4];


  densityEnv = wfn = 0.e0;
  for(i=0;i<MAXLOADPARAMETER;i++)
    par[i] = 0.e0;

//yPlus = 0.e0;
/*...*/
//if (fWallModel && fCal) {
/*... calculo da velocidade paralela a face*/
//  yPlus = wallPar[0];
//  uF    = wallPar[2];
/*...*/ 
//    tempPlus = wallModelHeat(yPlus,prM,prT);
/*...................................................................*/ 
//}
/*...................................................................*/ 

/*...*/
  if (lFaceReK == STATICWALL || ld->type == MOVEWALL) {
    tA[0]   = 0.0;
    if(!fCal) return;
    aP   = viscosityC*fArea/dd;
//  if(yPlus <= 20.0)
    *sP += aP;
 // *p += aP*uF*uF/0.3;
  }
/*...................................................................*/ 

/*... potencial prescrito (entra)*/
   else if( ld->type == INLET){
     getLoads(par,ldVel,xx);
     densityEnv = par[4];
     vB[0]      = par[0];
     vB[1]      = par[1];
     vv         = vB[0]*vB[0] + vB[1]*vB[1]; 
     wfn        = vB[0]*n[0] + vB[1]*n[1];
     if( ndm == 3){
       vB[2] = par[2];
       wfn += vB[2]*n[2];
       vv  += vB[2]*vB[2];
     }
     getLoads(par,ld,xx);
     in    = par[0];
     tA[0] = 1.5e0*in*in*vv;
/*...*/
     if(fCal)
       *p -= wfn*densityEnv*fArea*tA[0]; 
   }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
   else if( ld->type == OUTLET){
/*...*/
     if(fCal){
       wfn   = velC[0] * n[0] + velC[1] * n[1];
       if (ndm == 3) wfn += velC[2]*n[2];
       *sP += wfn*densityC*fArea;
     } 
/*...................................................................*/
  }
/*...................................................................*/

}
/*********************************************************************/

/**********************************************************************
 * Data de criacao    : 10/05/2019                                    *
 * Data de modificaco : 00/00/0000                                    *
 *--------------------------------------------------------------------*
 * initPropCD : inicializacao                                         *
 *--------------------------------------------------------------------*
 * Parametros de entrada:                                             *
 *--------------------------------------------------------------------*
 *--------------------------------------------------------------------*
 * Parametros de saida:                                               *
 *--------------------------------------------------------------------*
 *--------------------------------------------------------------------*
 * OBS:                                                               *
 *--------------------------------------------------------------------*
 **********************************************************************/
void initLoads()
{
  short i;

  for(i=0;i<MAXLOAD;i++)
  {
    loadsD1[i].type      = 0;
    loadsT1[i].type      = 0;
    loadsPres[i].type    = 0;
    loadsPresC[i].type   = 0;
    loadsEnergy[i].type  = 0;
    loadsTemp[i].type    = 0;
    loadsKturb[i].type   = 0;
    loadsZcomb[i].type   = 0;

    loadsD1[i].fUse      = false;
    loadsT1[i].fUse      = false;
    loadsPres[i].fUse    = false;
    loadsPresC[i].fUse   = false;
    loadsEnergy[i].fUse  = false;
    loadsTemp[i].fUse    = false;
    loadsKturb[i].fUse   = false;
    loadsZcomb[i].fUse   = false;

  }
}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 13/10/2019                                    *
 * Data de modificaco : 00/00/0000                                    *
 *--------------------------------------------------------------------*
 * gradPbBuoyant :                                                    *
 *--------------------------------------------------------------------*
 * Parametros de entrada:                                             *
 *--------------------------------------------------------------------*
 *--------------------------------------------------------------------*
 * Parametros de saida:                                               *
 *--------------------------------------------------------------------*
 *--------------------------------------------------------------------*
 * OBS:                                                               *
 *--------------------------------------------------------------------*
 **********************************************************************/
void gradPbBuoyant(DOUBLE *RESTRICT gradPb, DOUBLE *RESTRICT gradRho
                 , DOUBLE *RESTRICT g     , DOUBLE const gh  
                 , DOUBLE const densityC  , DOUBLE const densityRef
                 , short const iCod)
{

  DOUBLE tmp;

  
  if(iCod == BUOYANT_RHOREF)
  {
    tmp       = densityC - densityRef;
    gradPb[0] = tmp*g[0];
    gradPb[1] = tmp*g[1];
    gradPb[2] = tmp*g[2];    
  }
  else if(iCod == BUOYANT_PRGH)
  {
    gradPb[0] =-gh*gradRho[0];
    gradPb[1] =-gh*gradRho[1];
    gradPb[2] =-gh*gradRho[2];
  }
  else
    gradPb[0] = gradPb[1] = gradPb[2] = 0.e0;
/*...................................................................*/



}