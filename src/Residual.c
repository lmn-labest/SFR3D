#include<Residual.h>
/********************************************************************* 
 * Data de criacao    : 12/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * rscaled:                                                          *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * u     -> campo de u                                               *
 * ad    -> diagonal das equacoes                                    *
 * rCell -> residuo das equacoes por celulas                         *
 * id    -> numero da equacao da celula i                            *
 * nEl   -> numero de elementos                                      *
 * nEq   -> numero de equacoes                                       *
 * ndf   -> graus de liberdade                                       *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * rU       -> residuo                                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void rScaled(DOUBLE *RESTRICT u    ,DOUBLE *RESTRICT ad
            ,DOUBLE *RESTRICT rCell,DOUBLE *RESTRICT rU
            ,INT  *RESTRICT id
            ,INT const nEl         ,INT const nEq
            ,short const ndf) 
{
  INT i,j,lNeq;
  DOUBLE mod,max[MAXSPECIES],v,rScale;


  for(j=0;j<ndf;j++)
    max[j] = rU[j] = 0.e0;

/*... max(Ap*uP) */
  for(i=0;i<nEl;i++)
  {
    lNeq = id[i] - 1;
    if(lNeq > -1)
    { 
      for(j=0;j<ndf;j++)
      {
        v       = MAT2D(i,j,u,ndf);
        lNeq   += j*nEq;  
        mod     = fabs(ad[lNeq]*v);
        max[j] = max(max[j],mod);
      }
    }
  }
/*...................................................................*/
      
/*... max ( | F - Ax |P / max(Ap*uP) )*/
  for(j=0;j<ndf;j++)
  {
    for(i=0;i<nEl;i++)
    {
      mod    = fabs(MAT2D(j,i,rCell,nEl));
      if(max[j] != 0.e0)
        rScale = mod/max[j];
      else
        rScale = 0.e0;
      rU[j]  = max(rU[j],rScale);
    }
  }  
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 12/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * rscaled:                                                          *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * u     -> campo de u                                               *
 * ad    -> diagonal das equacoes                                    *
 * rCell -> residuo das equacoes por celulas                         *
 * id    -> numero da equacao da celula i                            *
 * nEl   -> numero de elementos                                      *
 * nEq   -> numero de equacoes                                       *
 * ndf   -> graus de liberdade                                       *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * rU       -> residuo                                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void rScaledSum(DOUBLE *RESTRICT u    ,DOUBLE *RESTRICT ad
               ,DOUBLE *RESTRICT rCell,DOUBLE *RESTRICT rU
               ,INT  *RESTRICT id
               ,INT const nEl         ,INT const nEq
               ,short const ndf) 
{
  INT i,j,lNeq;
  DOUBLE mod,sum[MAXSPECIES],v,rScale;


  for(j=0;j<ndf;j++)
    sum[j] = rU[j] = 0.e0;

/*... max(Ap*uP) */
  for(i=0;i<nEl;i++)
  {
    lNeq = id[i] - 1;
    if(lNeq > -1)
    { 
      for(j=0;j<ndf;j++)
      {
        v       = MAT2D(i,j,u,ndf);
        lNeq   += j*nEq;  
        mod     = fabs(ad[lNeq]*v);
        sum[j] += mod;
      }
    }
  }
/*...................................................................*/
      
/*... max ( | F - Ax |P / max(Ap*uP) )*/
  for(j=0;j<ndf;j++)
  {
    for(i=0;i<nEl;i++)
    {
      mod    = fabs(MAT2D(j,i,rCell,nEl));
      rU[j] += mod;
    }
    if( sum[j] > rU[j]*SZERO)
      rU[j]  /=  sum[j];  
  }  
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 12/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * rscaled:                                                          *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * u     -> campo de u                                               *
 * ad    -> diagonal das equacoes                                    *
 * rCell -> residuo das equacoes por celulas                         *
 * id    -> numero da equacao da celula i                            *
 * nEl   -> numero de elementos                                      *
 * nEq   -> numero de equacoes                                       *
 * ndf   -> graus de liberdade                                       *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * rU       -> residuo                                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void rScaledSumMax(DOUBLE *RESTRICT u    ,DOUBLE *RESTRICT ad
                  ,DOUBLE *RESTRICT rCell,DOUBLE *RESTRICT rU
                  ,INT  *RESTRICT id
                  ,INT const nEl         ,INT const nEq
                  ,short const ndf) 
{
  INT i,j,lNeq;
  DOUBLE mod,sum[MAXSPECIES],sumMax,v,rScale;


  for(j=0;j<ndf;j++)
    sum[j] = rU[j] = 0.e0;

/*... max(Ap*uP) */
  for(i=0;i<nEl;i++)
  {
    lNeq = id[i] - 1;
    if(lNeq > -1)
    { 
      for(j=0;j<ndf;j++)
      {
        v       = MAT2D(i,j,u,ndf);
        lNeq   += j*nEq;  
        mod     = fabs(ad[lNeq]*v);
        sum[j] += mod;
      }
    }
  }
/*...................................................................*/
      
/*...*/
  sumMax = sum[0];
  for (j = 1; j < ndf; j++) 
    sumMax = max(sumMax,sum[j]);
/*...................................................................*/

/*... sum ( | F - Ax |P / sum( |Ap*velP| ) )*/
  for(j=0;j<ndf;j++)
  {
    for(i=0;i<nEl;i++)
    {
      mod    = fabs(MAT2D(j,i,rCell,nEl));
      rU[j] += mod;
    }
    rU[j]  /= sumMax;  
  }  
/*...................................................................*/
}
/*********************************************************************/