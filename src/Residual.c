#include<Residual.h>

/********************************************************************* 
 * Data de criacao    : 13/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * RESIDUALCOMBUSTION : calculo dos residuos no metodo simple        *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * vel      -  > campo de velocidade                                 *
 * energy     -> campo de energia                                    *  
 * zComb      -> campo com a fracao massica das especies             *
 * rCellVel   -> residuo das equacoes das velocidade por celulas     * 
 * rCellMass  -> residuo de massa da equacoes de pressao por celulas * 
 * rCellEnergy-> residuo de massa da equacoes de energia por celulas * 
 * adVel    -> diagonal da equacoes de velocidades                   *
 * adEnergy -> diagonal da equacoes de energia                       *  
 * rU       -> nao definido                                          * 
 * rMass    -> nao definido                                          * 
 * rEnergy  -> nao definido                                          * 
 * rComb    -> nao definido                                          * 
 * idVel    -> numero da equacao da celula i                         * 
 * idEnergy -> numero da equacao da celula i                         * 
 * idComb   -> numero da equacao da celula i                         * 
 * nEl      -> numero de elementos                                   * 
 * nEqVel   -> numero de equacoes de velocidade                      * 
 * ndm      -> numero de dimensoes                                   * 
 * nComb    -> numero de especies explicitamente resolvidas          * 
 * iCod     -> tipo de residuo                                       * 
 *          RSCALED - residuo com escala de grandeza                 * 
 *          RSQRT   - norma p-2 ( norma euclidiana)                  * 
 *          RSCALEDM- residuo com escala de grandeza                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * rU       -> residuo das velocidades                               * 
 * rMass    -> residuo de mass                                       * 
 * rEnergy  -> residuo de energia                                    * 
 * rComb    -> resido das especies                                   *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *            | rvx1 rvx2 ... rvxnEl |                               *
 * rCellVel = | rvy1 rvy2 ... rvynEl |                               * 
 *            | rvz1 rvz2 ... rvznEl |                               * 
 *                                                                   *
 *         | adx1 adx2 ... adxnEq |                                  *
 * adVel = | ady1 ady2 ... adynEq |                                  *
 *         | adz1 adz2 ... adznEq |                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void residualCombustion(DOUBLE *RESTRICT vel ,DOUBLE *RESTRICT energy
            ,DOUBLE *RESTRICT zComb         
            ,DOUBLE *RESTRICT rCellVel    ,DOUBLE *RESTRICT rCellMass
            ,DOUBLE *RESTRICT rCellEnergy ,DOUBLE *RESTRICT rCellComb      
            ,DOUBLE *RESTRICT adVel       ,DOUBLE *RESTRICT adEnergy 
            ,DOUBLE *RESTRICT adComb   
            ,DOUBLE *RESTRICT rU          ,DOUBLE *rMass
            ,DOUBLE *rEnergy              ,DOUBLE *rComb    
            ,INT  *RESTRICT idVel         ,INT  *RESTRICT idEnergy
            ,INT  *RESTRICT idComb    
            ,INT const nEl                ,INT const nEqVel
            ,INT const nEqComb
            ,short const ndm              ,short const nComb
            ,short *iCod)

{

/*... velocidade*/
  residualType(vel    ,rCellVel
              ,adVel  ,rU
              ,idVel  ,nEl   
              ,nEqVel ,ndm 
              ,iCod[0]);
/*...................................................................*/

/*... Energia*/
  residualType(energy   ,rCellEnergy
              ,adEnergy ,rEnergy
              ,idEnergy ,nEl   
              ,0        ,1   
              ,iCod[1]);
/*...................................................................*/

/*... Conservao das especies*/
  residualType(zComb    , rCellComb  
              ,adComb   , rComb  
              ,idComb   , nEl   
              ,nEqComb  , nComb   
              ,iCod[2]);
/*...................................................................*/

/*... Conservao de massa*/
  *rMass = residualMass(rCellMass,nEl,iCod[3]);
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 13/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * residualType   : calculo dos residuos                             *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * u        -  > campo de velocidade                                 *
 * rCell      -> residuo das equacoes das velocidade por celulas     * 
 * ad       -> diagonal da equacoes de velocidades                   *  
 * rU       -> nao definido                                          * 
 * id       -> numero da equacao da celula i                         * 
 * nEl      -> numero de elementos                                   * 
 * nEq      -> numero de equacoes de velocidade                      * 
 * ndf      -> graus de liberdade                                    *   
 * iCod     -> tipo de residuo                                       * 
 *          RSCALED - residuo com escala de grandeza                 * 
 *          RSQRT   - norma p-2 ( norma euclidiana)                  * 
 *          RSCALEDM- residuo com escala de grandeza                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * rU       -> residuo das velocidades                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *            | rvx1 rvx2 ... rvxnEl |                               *
 * rCellVel = | rvy1 rvy2 ... rvynEl |                               * 
 *            | rvz1 rvz2 ... rvznEl |                               * 
 *                                                                   *
 *         | adx1 adx2 ... adxnEq |                                  *
 * adVel = | ady1 ady2 ... adynEq |                                  *
 *         | adz1 adz2 ... adznEq |                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void residualType(DOUBLE *RESTRICT u      ,DOUBLE *RESTRICT rCellU
                 ,DOUBLE *RESTRICT adU    ,DOUBLE *RESTRICT rU
                 ,INT  *RESTRICT idU      ,INT const nEl   
                 ,INT const nEqU          ,short const ndf 
                 ,short iCod)

{
  short j;
  DOUBLE *p=NULL;
  
/*... Velocidade*/
  switch(iCod)
  {

/*... scaled*/
    case RSCALED:
      rScaled(u                     ,adU  
             ,rCellU                ,rU
             ,idU   
             ,nEl                   ,nEqU   
             ,ndf);
    break;
/*...................................................................*/

/*... norma euclidiana*/
    case RSQRT:
/*...*/
      for(j=0;j<ndf;j++)
      {
        p     = &rCellU[j*nEl]; 
        rU[j] = sqrt(dot(p,p,nEl));
      }
/*...................................................................*/
    break;
/*...................................................................*/

/*... scaled*/
    case RSCALEDSUM:
      rScaledSum(u     , adU   
                ,rCellU, rU
                ,idU
                ,nEl   ,nEqU
                ,ndf); 

    break;
/*...................................................................*/

/*... scaled*/
    case RSCALEDSUMMAX:
     rScaledSumMax(u     ,adU
                  ,rCellU,rU
                  ,idU
                  ,nEl   ,nEqU
                  ,ndf); 
    break;
/*...................................................................*/

/*... */
    default:
      ERRO_OP(__FILE__,__func__,iCod);
    break;
/*...................................................................*/
  }

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 12/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * residualMass:                                                     *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * rCellMass  -> residuo de massa da equacoes de pressao por celulas * 
 * nEl   -> numero de elementos                                      *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * rU       -> residuo                                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE residualMass(DOUBLE *RESTRICT rCellMass,INT const nEl
                  , short const iCod)
{

  INT i;
  DOUBLE tmp,v;

/*...*/
  if( iCod == RSQRT)
    tmp = sqrt(dot(rCellMass,rCellMass,nEl));
/*...................................................................*/

/*...*/
  else
  {
    tmp = 0.e0;
    for(i=0;i<nEl;i++)
    {
      v    = fabs(rCellMass[i]);
      tmp += v;
    }
  } 
/*...................................................................*/

  return tmp;

}
/********************************************************************/


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
 * rscaledSum:                                                       *
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
 * rscaledSumMax:                                                    *
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