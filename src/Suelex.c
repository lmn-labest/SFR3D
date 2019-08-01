#include<Edo.h>

#define KMAXX 12

/*********************************************************************
 * Data de criacao    : 16/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * polYextr: polinomio de interpolacao                               *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*
 * *******************************************************************/
static void polYextr(DOUBLE *RESTRICT t   ,DOUBLE *RESTRICT c
                   , DOUBLE *RESTRICT y   ,short const nEdo
                   , short const k        )
{
  
  INT i,j,l,ct,cc=KMAXX + 1;

  l = ct = nEdo;

  for (j=k-1; j>0; j--)
  {
    for (i=0; i<l; i++)
      MAT2D(j-1,i,t,ct) = MAT2D(j,i,t,ct) 
     + MAT2D(k,j,c,cc)*(MAT2D(j,i,t,ct) -MAT2D(j-1,i,t,ct) );

    for (i=0; i<l; i++)
      y[i] = MAT2D(0,i,t,ct) 
       + MAT2D(k,0,c,cc)*(MAT2D(0,i,t,ct) - y[i] );

  }

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 16/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * dy:                                                               *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*
 * *******************************************************************/
static bool dy(DOUBLE *RESTRICT a       , DOUBLE *RESTRICT w   
       , DOUBLE *RESTRICT dfdy   , DOUBLE *RESTRICT del 
       , DOUBLE *RESTRICT y      , DOUBLE *RESTRICT yt
       , DOUBLE *RESTRICT yEnd   , DOUBLE *RESTRICT dyt  
       , DOUBLE *RESTRICT scale  , void  **pt
       , INT *RESTRICT p         , short *RESTRICT nSeq 
       , DOUBLE const hTot       , DOUBLE const x
       , short const nEdo        , short const k
       , void (*rhs)()           ) 
{
  short i,j,nn,nStep=nSeq[k];
  DOUBLE xNew,del1,del2,theta,h;

  h = hTot/nStep;

/*... A = (1/h - dydx)*/
  for (i = 0; i<nEdo; i++)
  { 
      for (j = 0; j<nEdo; j++) 
        MAT2D(i,j,a,nEdo) = -MAT2D(i, j, dfdy, nEdo);
      MAT2D(i, i, a, nEdo) += 1.e0 / h;
  }
/*...................................................................*/

/*... decomposicao de a*/
  fatLUpp(a,p, nEdo, LUKIJPP);
/*...................................................................*/

/*...*/
  xNew = x + h;       /*Special step for nonautonomous system.*/
  rhs(xNew, y, del, pt);
/*...*/
  for (i = 0; i<nEdo; i++)
    yt[i] = y[i];        
/*...................................................................*/

/*...*/
  solvLUpp(a, del, w, p, nEdo);
/*...................................................................*/

/*...*/
  for(nn=1;nn<nStep;nn++)
  {
    for (i=0;i<nEdo;i++)
      yt[i] += del[i];
    
    xNew += h;
    rhs(xNew, yt, yEnd, pt);
/*... Stability test and test for recomputing Jac*/
    if (nn ==1 && k<=1)
    { 

      for (i=0,del1=0.0;i<nEdo;i++)
        del1 += sqr(del[i]/scale[i]);
      del1=sqrt(del1);

      rhs(x+h,yt,dyt,pt);
      for (i=0;i<nEdo;i++)
        del[i]=dyt[i]-del[i]/h;

      solvLUpp(a, del, w, p, nEdo);
    
      for (i=0,del2=0.0;i<nEdo;i++)
        del2 += sqr(del[i]/scale[i]);
      del2=sqrt(del2);
      
      theta=del2/max(1.0,del1);
     
      if (theta > 1.0)
        return false;
    } 
/*...................................................................*/

/*...*/
     solvLUpp(a, yEnd, w, p, nEdo);
     for(i=0;i<nEdo;i++)
       del[i] = yEnd[i];
/*...................................................................*/
    
  }
/*...................................................................*/

/*... Last Step*/
  for (i=0;i<nEdo;i++) 
    yEnd[i]=yt[i]+del[i];
/*...................................................................*/


  return true;
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 16/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * StepperSie :                                                      *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*
 * y     -> valores inicial                                          *
 * param -> paramatros utilizado no calculo da funcao f              *
 * x1    -> limite inferio do intervalo de integracao                *
 * x2    -> limite superior do intervalo de integraçao               *
 * hInit -> passo inicial                                            *
 * maxIt -> numero maximo de iteração no processo de integracao      *
 * nEdo  -> numero de edos                                           *
 * maxInt-> numero maximo de passos de integracao                    * 
 * fStopIt -> true para quando maxInt é atingido                     *
 *           false nao a limite de passo de initegra                 *
 * outPut  -> cod de escrita                                         *
 *            1 - tela                                               *
 *            2 - arquivo                                            *
 * fName   -> nomde do aquivo de saida                               *
 * rhs     -> f                                                      *
 * jacY    -> parte do jacabiano da f dependo explicidamente de y    *
 * mass    -> matriz de massa                                        *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * y       -> resultado                                              *
 * it      -> numero total passos dados                              *
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*
 * *******************************************************************/
int StepperSie(DOUBLE *y                  , void **pt
               , DOUBLE const x1          , DOUBLE const x2
               , DOUBLE const hInit       , DOUBLE const hMax 
               , DOUBLE const aTol        , DOUBLE const rTol        
               , short const nEdo         , DOUBLE *xF
               , unsigned INT const maxInt, bool const fStopIt
               , short const  outCod      , const char *const fName
               , void(*rhs)()             , void(*jacY)()) 
{
   int const iMaxx = KMAXX + 1;
   const DOUBLE costJac = 5.0e0
                    , costLu = 1.0e0
                    , costFunc = 1.0e0
                    , costSolve = 1.0e0;
  
   const DOUBLE STEPFAC1=0.6e0
               ,STEPFAC2=0.93e0
               ,STEPFAC3=0.1e0
               ,STEPFAC4=4.0e0
               ,STEPFAC5=0.5e0
               ,KFAC1=0.7e0
               ,KFAC2=0.9e0;
  FILE *file = NULL;
  bool firstStep=true
      ,lastStep=false      
      ,reject=false
      ,prevReject=false
      ,forward
      ,firstK
      ,sucess  = true 
      ,calcJac = false;
  short i,j,k,nSeq[KMAXX + 1],maxChangeSteps, iChangeSteps;
  long INT fac[KMAXX + 1];
  INT kTarg,kOpt,kMax=0;
  DOUBLE logFact,ratio,coeff[KMAXX + 1][KMAXX + 1];
  DOUBLE theta,jacRedo,hNext,hNew,hDid,h,x,xOld,err,errOld,eps;
  DOUBLE cost[KMAXX+1],work[KMAXX+1],hOpt[KMAXX+1],expo;
  DOUBLE facMin,fa,ySafe,delt,uround;  
  unsigned INT it,IntegralStepMax;
  INT p[MAXSPECIES+1];
  DOUBLE scale[MAXSPECIES+1]
        ,dfdy[(MAXSPECIES+1)*(MAXSPECIES+1)]
        ,ySeq[MAXSPECIES+1]
        ,table[(MAXSPECIES+1)*KMAXX]
        ,ySav[MAXSPECIES+1]
        ,a[(MAXSPECIES+1)*(MAXSPECIES+1)]
        ,w[MAXSPECIES+1]
        ,del[MAXSPECIES+1]
        ,yt[MAXSPECIES+1]
        ,dyt[MAXSPECIES+1];

/*...*/
  it      = 0;
  hNext   = hInit;
  x       = x1;
  eps     = 1.e-15;
  uround  = 1.e-06;  
/*...................................................................*/

/*...*/  
  maxChangeSteps  = 100;
  IntegralStepMax = maxInt;
/*...................................................................*/

/*...*/
  if (outCod == FILE_OUT)
    file = openFile(fName, "w");
/*...................................................................*/

/*...  
        y0    y1      yn   
       T00   T00 ... T00    
       T10   T10 ... T10 
       T11   T11 ... T11
       T20   T20 ... T20
       T21   T21 ... T21
       T22   T22 ... T22
       T30   T30 ... T30
       T31   T31 ... T31 
       T32   T32 ... T32
       T33   T33 ... T33
       ...   ... ... ...
       TK1   TK1 ... TK1
       TK2   TK2 ... TK2
       TK3   TK3 ... TK3
       ...   ... ... ...
       TKK   TKK ... TKK
 */
/*...*/
  jacRedo=min(1.0e-4,rTol);
  theta=2.0*jacRedo; /*garante que o jacobiano e calcula do primeiro passo*/
/*...................................................................*/

/*... sequencia*/
  nSeq[0] = 2;
  nSeq[1] = 3;
  for(i=2;i<iMaxx;i++)
    nSeq[i] = 2*nSeq[i-2];
/*.....................................................................*/

/*...*/
  cost[0] = costJac + costLu + nSeq[0]*(costFunc+costSolve);
  for(i=0;i<KMAXX;i++)
    cost[i+1]=cost[i]+(nSeq[i+1]-1)*(costFunc+costSolve)+costLu;
/*.....................................................................*/
  
//hNext=-1.0e99;

  logFact=-log10(rTol+aTol)*0.6e0+0.5e0;
  kTarg  =max(2,min(KMAXX-1,(INT) logFact)); 
  for (i=0; i<iMaxx; i++)
  { 
    for(j=0; j<i; j++) 
    { 
      ratio= (DOUBLE) nSeq[i]/nSeq[j];
      coeff[i][j]=1.e0/(ratio-1.e0);
    }
  }
/*... fatorial*/
  fac[0]=1;
  for (i=0; i<iMaxx-1; i++)
    fac[i+1]=(i+1)*fac[i]; 
/*.....................................................................*/

/*...*/
  do
  {    
    iChangeSteps    = 0; 
    work[0]=1.e30;
    h = min(hNext,hMax);
    
    if(fabs(x+h)>fabs(x2))
      h = x2 - x;
     
    hOpt[0] = 0.1e0*fabs(h);
    forward = h > 0 ? true : false;
/*... Guardando os valores inicias*/
    for (i=0;i<nEdo;i++) 
      ySav[i] = y[i]; 

    if(h!=hNext && !firstStep)
      lastStep = true;

/*... ultimo h rejetado */
    if(reject)
    {
      prevReject = true;
      lastStep   = false;
      theta      = 2.0*jacRedo; /* garante o caluclo do jac*/
    }
/*.....................................................................*/

/*...*/
    for(i=0;i<nEdo;i++)
      scale[i] = aTol + rTol*fabs(y[i]);
/*.....................................................................*/

/*...*/
    reject = false;
    firstK = true;
    hNew   = fabs(h);
/*... Evaluate Jacobian*/
    compute_jac:
    if (theta > jacRedo && !calcJac)
    { 
/*... diferenca finita*/
      if (jacY == NULL)
      {
        rhs(x, y,dyt,pt);
        for(i=0;i<nEdo;i++)
        {
          ySafe = y[i];
          delt  = sqrt(uround*max(1.e-5,fabs(ySafe)));
          y[i]  = ySafe + delt;
          rhs(x, y,yt,pt);
          for(j=0;j<nEdo;j++)
            MAT2D(j,i,dfdy,nEdo) = (yt[j] - dyt[j])/delt;
          y[i]   = ySafe;
        }
      }

/*... analitico*/
      else
        jacY(x, y, dfdy ,pt);    
/*... derivada df/dy*/

      calcJac=true;
    }
/*...................................................................*/

/*... repetir ate um h aceitavel*/
    while(firstK || reject)
    {
      h = forward ? hNew : -hNew;
      firstK=false;
      reject=false;

/*...*/
      if (x == x + h)
      {
        printf("stepsize underflow in suelex!!\n");
        exit(EXIT_FAILURE);
      }
/*...................................................................*/
      
      errOld = 0.e0;

/*... loop da tabela*/
      for(k=0;k<=kTarg+1;k++)
      {
        sucess =  dy(a     , w
                  , dfdy , del
                  , ySav , yt
                  , ySeq , dyt
                  , scale, pt 
                  , p    , nSeq
                  , h    , x
                  , nEdo , k
                  , rhs  );
/*... reduce stepsize*/
        if(!sucess)
        {
          reject=true;
          hNew=fabs(h)*STEPFAC5;
          break; /*sai do loop da tabela*/
        }
/*...................................................................*/

/*... primeira entrada da tabela*/
        if(k==0)
        {
          for (i = 0; i<nEdo; i++)
            y[i] = ySeq[i];
        }
        else
        {
          for (i = 0; i<nEdo; i++)
            MAT2D(k-1,i,table,nEdo) = ySeq[i];
        }
/*...................................................................*/

/*...*/
        if(k != 0)
        {
/*...extrapolacao*/
          polYextr(table,coeff[0],y,nEdo,k); 
/*...*/
          for (i=0,err=0.0;i<nEdo;i++)
          {
//          scale[i]=aTol+rTol*fabs(ySav[i]);
            err+=sqr((y[i]-MAT2D(0,i,table,nEdo))/scale[i]);
          }
          err=sqrt(err/nEdo);  
/*...................................................................*/

/*... reduce stepsize*/
          if(err > 1.0/eps || ( k > 1 && err >= errOld))
          {
            reject = true;
            hNew   = fabs(h)*STEPFAC5;
            break;/*sai do loop da tabela*/ 
          }
          errOld = max(4.0*err,1.e0);
/*...................................................................*/

          expo=1.e0/(k+1.e0);
/*Compute optimal stepsize for this order. Note k instead of 2k in exponent.*/
          facMin=pow(STEPFAC3,expo);
          if (err == 0.0)
            fa=1.0/facMin;
          else 
          {
            fa=STEPFAC2/pow(err/STEPFAC1,expo); /*(S2/errk)^(1/k+1)*/
            fa=max(facMin/STEPFAC4,min(1.e0/facMin,fa));
          }
          hOpt[k] = fabs(h*fa);                /*... Hk = HS1(S2/errk)^(1/k+1)*/
          work[k] = cost[k]/hOpt[k];
          if ((firstStep || lastStep) && err <= 1.0)
            break; /*sai do loop da tabela*/ 
/*...*/
          if (k == kTarg-1 && !prevReject && !firstStep && !lastStep)
          {
            if (err <= 1.0) /*Converged within order window.*/
              break;/*sai do loop da tabela*/ 
            else if (err>nSeq[kTarg]*nSeq[kTarg+1]*4.0) 
            {
              reject=true;        /*No convergence expected by k_targ+1.*/
              kTarg=k;
              if (kTarg>1 && work[k-1]<KFAC1*work[k])
                kTarg--;
              hNew=hOpt[kTarg];
              break;/*sai do loop da tabela*/ 
            }
          }
/*...................................................................*/

/*...*/
          if (k == kTarg) 
          {
            if (err <= 1.0) /*Converged within order window.*/
              break;/*sai do loop da tabela*/ 
            else if (err>nSeq[k+1]*2.0) 
            {
              reject=true; /*No convergence expected by k_targ+1.*/
              if (kTarg>1 && work[k-1]<KFAC1*work[k])
                kTarg--;
              hNew=hOpt[kTarg];
              break;/*sai do loop da tabela*/ 
            }
          }
/*...................................................................*/

/*...*/
          if (k == kTarg+1)
          {
            if (err > 1.0)
            {
              reject=true;
              if (kTarg>1 && work[kTarg-1]<KFAC1*work[kTarg])
                kTarg--;
              hNew=hOpt[kTarg];
            }
            break;/*sai do loop da tabela*/ 
          }
/*...................................................................*/
      
/*... Go back and try next k.*/
          if (reject) /*Arrive here from any break in for loop.*/
          {            
            prevReject=true;
            if (!calcJac)
            {
              theta=2.0*jacRedo;
              goto compute_jac;
            }
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/   

/*...*/
      iChangeSteps++;
      if(iChangeSteps == maxChangeSteps)
      {
        printf("Numero maximo de mudanca de passo alcancada!!\n");
        exit(EXIT_FAILURE);
      }
/*...................................................................*/  
    }
/*...................................................................*/

/*...*/
    calcJac = false;  /*Successful step. Allow Jacobian to be recomputed
                                                   if theta too big.*/  
    xOld     = x;
    x       += h;
    hDid     = h;
    firstStep=false;
    kMax     = max(kMax,kTarg);

/*... Determine optimal order for next step*/
    if(k==1)
    {
      kOpt = max(3,KMAXX-1);
      if(reject)
        kOpt = 2;
    }  
/*...................................................................*/

/*...*/
    else if (k <= kTarg)
    {
      kOpt = k;
      if (work[k-1] < KFAC1*work[k])
        kOpt=k-1;
      else if (work[k] < KFAC2*work[k-1])
        kOpt=min(k+1,KMAXX-1);
    }
/*...................................................................*/

/*...*/
    else
    {
      kOpt=k-1;
      if (k > 2 && (work[k-2] < KFAC1*work[k-1]))
        kOpt=k-2;
      if (work[k] < KFAC2*work[kOpt])
        kOpt=min(k,KMAXX-1);
    }  
/*...................................................................*/

/*...After a rejected step neither order nor size should increase*/
    if (prevReject)
    {  
      kTarg = min(kOpt,k);
      hNew=min(fabs(h),hOpt[kTarg]);
      prevReject=false;
    }  
/*... Stepsize control for next step.*/
    else 
    { 
      if (kOpt <= k)
        hNew=hOpt[kOpt];
      else
      {
        if (k<kTarg && work[k]<KFAC2*work[k-1])
          hNew=hOpt[k]*cost[kOpt+1]/cost[k];
        else
          hNew=hOpt[k]*cost[kOpt]/cost[k];
      }
      kTarg=kOpt;
    }
/*...................................................................*/

/*...*/
    if (forward)
      hNext=hNew;
    else
      hNext=-hNew;
/*...................................................................*/

/*...*/
     writeOutput(it            , x
                ,h             , y
                ,nEdo          , outCod
                ,file);
/*.....................................................................*/

/*...*/
    it++;
/*...*/
    if(it > IntegralStepMax && fStopIt)
    {
       printf("Numero maximo de passos para integracao alcancado !!\n");
      exit(EXIT_FAILURE);
    }
/*.....................................................................*/

  }while(x < x2);
/*.....................................................................*/

  *xF = x;

  return it;
}
/***********************************************************************/