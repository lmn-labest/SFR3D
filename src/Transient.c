#include<Transient.h>  
/********************************************************************* 
 * Data de criacao    : 00/00/2016                                   *
 * Data de modificaco : 27/08/2016                                   *
 *-------------------------------------------------------------------*
 * SETTRANSIENTSCHEME : set a discretizacao temporal                 * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * word -> str com a discretizacao                                   * 
 * type -> nao definido                                              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * type -> tipo de discretizacao                                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void setTransientScheme(char *word,short *type){

  if(!strcmp(word,"EULER"))
   *type = EULER;
  else if(!strcmp(word,"BACKWARD"))
   *type =  BACKWARD;

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/2016                                   *
 * Data de modificaco : 27/08/2016                                   *
 *-------------------------------------------------------------------*
 * CELLTRANSIENT : discretizacao temporal                            *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * volume    -> volume das celulas                                   *
 * id        -> numera das equacoes                                  *
 * u0        -> solucao (n-1)                                        *
 * u         -> solucao (n)                                          *
 * density   -> massa especifica no tempo   (n)                      *
 * f         -> vetor de forcas                                      *
 * ddt       -> delta t                                              *
 * numel     -> numero de elementos                                  *
 * ndf       -> graus de liberade                                    *
 * type      -> tipo de discretizacao temporal                       *
 * fAdd      -> acumula o valor em f                                 *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * f -> atualizado com a discreticao temporal ( fAdd = true)         *
 *      sobreescreve com a discreticao temporal (fAdd = false)       *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
*********************************************************************/
void cellTransient(DOUBLE *RESTRICT volume ,INT *RESTRICT id
                  ,DOUBLE *RESTRICT u0     ,DOUBLE *RESTRICT u
                  ,DOUBLE *RESTRICT density,DOUBLE *RESTRICT f
                  ,Temporal const ddt      ,INT const numel
                  ,short const ndf         ,bool const fAdd)
{
  INT nel, lNeq;
  DOUBLE t0, t00, tmp1;
  DOUBLE dt = ddt.dt[0], dt0 = ddt.dt[1];
  short type = ddt.type;
  short j, nD = DENSITY_LEVEL;

/*...*/
  switch (type) 
  {
/*... EULER de primeira ordem*/
    case EULER:
/*... acumula em f*/
      if (fAdd) 
      {
/*... ndf = 1*/
        if (ndf == 1) 
        {
          for (nel = 0; nel < numel; nel++) 
          {
            lNeq = id[nel] - 1;
            if (lNeq > -1) 
            {
              t0 = MAT2D(nel, TIME_N_MINUS_1, density, nD)*u[nel];
              f[lNeq] += volume[nel]*t0/dt;
            }
          }
        }
/*...................................................................*/

/*... ndf > 1*/
        else 
        {
          for (nel = 0; nel < numel; nel++) 
          {
            for (j = 0; j< ndf; j++) 
            {
              lNeq = MAT2D(nel, j, id, ndf) - 1;
              if (lNeq > -1) 
              {
                t0 = MAT2D(nel,TIME_N_MINUS_1,density,nD)
                    *MAT2D(nel,j,u,ndf);
                MAT2D(lNeq, j, f, ndf) += volume[nel] * t0 / dt;
              }
            }
          }
        }
/*...................................................................*/
      }
/*...................................................................*/

/*... sobrecreve f*/
      else 
      {
/*... ndf = 1*/
        if (ndf == 1) 
        {
          for (nel = 0; nel < numel; nel++) 
          {
            lNeq = id[nel] - 1;
            if (lNeq > -1)
            {
              t0 = MAT2D(nel, TIME_N_MINUS_1, density, nD)*u[nel];
              f[lNeq] = volume[nel]*t0/dt;
            }
          }
        }
/*...................................................................*/

/*... ndf > 1*/
        else 
        {
          for (nel = 0; nel < numel; nel++) 
          {
            for (j = 0; j< ndf; j++) 
            {
              lNeq = MAT2D(nel, j, id, ndf) - 1;
              if (lNeq > -1) 
              {
                t0 = MAT2D(nel,TIME_N_MINUS_1,density,nD)
                    *MAT2D(nel,j,u,ndf);
                MAT2D(lNeq,j,f,ndf) = volume[nel]*t0/dt;
              }
            }
          }
        }
/*...................................................................*/
      }
/*...................................................................*/
    break;
/*...................................................................*/

/*... BACKWARD de segunda ordem*/
  case BACKWARD:
/*... acumula em f*/
    if (fAdd) 
    {
/*... ndf = 1*/
      if (ndf == 1) 
      {
        for (nel = 0; nel < numel; nel++) 
        {
          lNeq = id[nel] - 1;
          if (lNeq > -1) 
          {
            tmp1 = dt + dt0;
/*...*/
            t0 = MAT2D(nel, TIME_N_MINUS_1, density, nD)*u[nel];
            t0 *= (tmp1/(dt*dt0));
/*...................................................................*/

/*...*/
            t00= MAT2D(nel, TIME_N_MINUS_2, density, nD)*u0[nel];
            t00*= (dt/(dt0*tmp1));
/*...................................................................*/
            f[lNeq] += volume[nel]*(t0-t00);
          }
        }
      }
/*...................................................................*/

/*... ndf > 1*/
      else 
      {
        for (nel = 0; nel < numel; nel++) 
        {
          for (j = 0; j< ndf; j++) 
          {
            lNeq = MAT2D(nel, j, id, ndf) - 1;
            if (lNeq > -1) 
            {
              tmp1 = dt + dt0;
/*...*/
              t0 = MAT2D(nel,TIME_N_MINUS_1,density,nD)
                  *MAT2D(nel, j, u, ndf);
              t0 *= (tmp1 / (dt*dt0));
/*...................................................................*/

/*...*/
              t00= MAT2D(nel, TIME_N_MINUS_2, density, nD)
                  *MAT2D(nel, j, u0, ndf);
              t00*= (dt / (dt0*tmp1));
/*...................................................................*/
              MAT2D(lNeq,j,f,ndf) += volume[nel] * (t0 - t00);
            }
          }
        }
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... sobrecreve f*/
    else 
    {
/*... ndf = 1*/
      if (ndf == 1) 
      {
        for (nel = 0; nel < numel; nel++) 
        {
          lNeq = id[nel] - 1;
          if (lNeq > -1) 
          {
            tmp1 = dt + dt0;
/*...*/
            t0 = MAT2D(nel, TIME_N_MINUS_1, density, nD)*u[nel];
            t0 *= (tmp1 / (dt*dt0));
/*...................................................................*/

/*...*/
            t00= MAT2D(nel, TIME_N_MINUS_2, density, nD)*u0[nel];
            t00*= (dt / (dt0*tmp1));
/*...................................................................*/
            f[lNeq] = volume[nel] *(t0-t00);
          }
        }
      }
/*...................................................................*/

/*... ndf > 1*/
      else 
      {
        for (nel = 0; nel < numel; nel++) 
        {
          for (j = 0; j< ndf; j++) 
          {
            lNeq = MAT2D(nel, j, id, ndf) - 1;
            if (lNeq > -1) 
            {
              tmp1 = dt + dt0;
/*...*/
              t0 = MAT2D(nel,TIME_N_MINUS_1,density,nD)*MAT2D(nel,j,u,ndf);
              t0 *= (tmp1 / (dt*dt0));
/*...................................................................*/

/*...*/
              t00= MAT2D(nel,TIME_N_MINUS_2,density,nD)*MAT2D(nel,j,u0,ndf);
              t00*= (dt / (dt0*tmp1));
/*...................................................................*/
              MAT2D(lNeq,j,f,ndf) = volume[nel]*(t0-t00);
            }
          }
        }
      }
/*...................................................................*/
    }
/*...................................................................*/
    break;
  }
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 00/00/2016                                   *
* Data de modificaco : 26/08/2017                                   *
*-------------------------------------------------------------------*
* CELLTRANSIENTEnergy : discretizacao temporal                      *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* volume    -> volume das celulas                                   *
* id        -> numera das equacoes                                  *
* u0        -> solucao (n-1)                                        *
* u         -> solucao (n)                                          *
* density   -> massa especifica no tempo   (n,n-1)                  *
* sHeat     -> calor especifico no tempo   (n,n-1)                  *
* f         -> vetor de forcas                                      *
* ddt       -> delta t                                              *
* numel     -> numero de elementos                                  *
* type      -> tipo de discretizacao temporal                       *
* fAdd      -> acumula o valor em f                                 *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* f -> atualizado com a discreticao temporal ( fAdd = true)         *
*      sobreescreve com a discreticao temporal (fAdd = false)       *
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
void cellTransientEnergy(DOUBLE *RESTRICT volume ,INT *RESTRICT id
                        ,DOUBLE *RESTRICT u0     ,DOUBLE *RESTRICT u
                        ,DOUBLE *RESTRICT density, DOUBLE *RESTRICT sHeat
                        ,DOUBLE *RESTRICT f
                        ,Temporal const ddt      ,INT const numel
                        ,bool const fAdd)
{
  INT nel, lNeq;
  DOUBLE t0, t00, tmp1,density0,sHeat0;
  DOUBLE dt = ddt.dt[0], dt0 = ddt.dt[1];
  short type = ddt.type;
  short nD = DENSITY_LEVEL,nH = SHEAT_LEVEL;

/*...*/
  switch (type) {
/*... EULER de primeira ordem*/
  case EULER:
/*... acumula em f*/
    if (fAdd) {
/*... ndf = 1*/
      for (nel = 0; nel < numel; nel++) {
        lNeq = id[nel] - 1;
        if (lNeq > -1) {
/*...*/
          density0 = MAT2D(nel, TIME_N_MINUS_1, density, nD);
          sHeat0   = MAT2D(nel, TIME_N_MINUS_1, sHeat  , nH);
/*...................................................................*/
          f[lNeq] += volume[nel] * sHeat0*density0*u[nel] / dt;
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... sobrecreve f*/
    else {
/*... */
      for (nel = 0; nel < numel; nel++) {
        lNeq = id[nel] - 1;
        if (lNeq > -1){
/*...*/
          density0 = MAT2D(nel, TIME_N_MINUS_1, density, nD);
          sHeat0   = MAT2D(nel, TIME_N_MINUS_1, sHeat, nH);
/*...................................................................*/
          f[lNeq] = volume[nel]*sHeat0*density0*u[nel] / dt;
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
    break;
/*...................................................................*/

/*... BACKWARD de segunda ordem*/
  case BACKWARD:
/*... acumula em f*/
    if (fAdd) {
      for (nel = 0; nel < numel; nel++) {
        lNeq = id[nel] - 1;
        if (lNeq > -1) {
          tmp1 = dt + dt0;
 /*... t(n-1)*/
          density0 = MAT2D(nel, TIME_N_MINUS_1, density, nD);
          sHeat0   = MAT2D(nel, TIME_N_MINUS_1, sHeat, nH);
 /*...................................................................*/

 /*...*/
          t0 = sHeat0*density0*u[nel];
          t0 *= (tmp1 / (dt*dt0));
 /*...................................................................*/

 /*...t(n-2)*/
          density0 = MAT2D(nel, TIME_N_MINUS_2, density, nD);
          sHeat0   = MAT2D(nel, TIME_N_MINUS_2, sHeat, nH);
 /*...................................................................*/

 /*...*/
          t00 = sHeat0*density0*u0[nel];
          t00 *= (dt / (dt0*tmp1));
 /*...................................................................*/
          f[lNeq] += volume[nel] * (t0 - t00);
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... sobrecreve f*/
    else {
      for (nel = 0; nel < numel; nel++) {
        lNeq = id[nel] - 1;
        if (lNeq > -1) {
          tmp1 = dt + dt0;
 /*... t(n-1)*/
          density0 = MAT2D(nel, TIME_N_MINUS_1, density, nD);
          sHeat0   = MAT2D(nel, TIME_N_MINUS_1, sHeat, nH);
/*...................................................................*/

/*... t(n-2)*/
          t0  = sHeat0*density0*u[nel];
          t0 *= (tmp1 / (dt*dt0));
/*...................................................................*/

/*...*/ 
          density0 = MAT2D(nel, TIME_N_MINUS_2, density, nD);
          sHeat0   = MAT2D(nel, TIME_N_MINUS_2, sHeat, nH);
/*...................................................................*/

/*...*/
          t00= sHeat0*density0*u0[nel];
          t00*= (dt / (dt0*tmp1));
/*...................................................................*/
          f[lNeq] = volume[nel] * (t0 - t00);
        }
/*...................................................................*/
      } 
/*...................................................................*/
    }
/*...................................................................*/

    break;
  }
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 00/00/2016                                   *
* Data de modificaco : 27/08/2016                                   *
*-------------------------------------------------------------------*
* CELLTRANSIENTSIMPLE : discretizacao temporal                      *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* volume    -> volume das celulas                                   *
* id        -> numera das equacoes                                  *
* u0        -> solucao (n-1)                                        *
* u         -> solucao (n)                                          *
* density   -> massa especifica no tempo   (n)                      *
* f         -> vetor de forcas                                      *
* ddt       -> delta t                                              *
* nEq       -> numero de equacoes                                   *
* numel     -> numero de elementos                                  *
* ndf       -> graus de liberade                                    *
* fAdd      -> acumula o valor em f                                 *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* f -> atualizado com a discreticao temporal ( fAdd = true)         *
*      sobreescreve com a discreticao temporal (fAdd = false)       *
*-------------------------------------------------------------------*
* OBS: f | bu(1) bu(2) ... bu(neq) |                                *
*        | bv(1) bv(2) ... bv(neq) |                                *
*        | bw(1) bw(2) ... bw(neq) |                                *
*-------------------------------------------------------------------*
*********************************************************************/
void cellTransientSimple(DOUBLE *RESTRICT volume ,INT *RESTRICT id
                        ,DOUBLE *RESTRICT u0     ,DOUBLE *RESTRICT u
                        ,DOUBLE *RESTRICT density,DOUBLE *RESTRICT f
                        ,Temporal const ddt      ,INT const nEq
                        ,INT const numel         ,short const ndf
                        ,bool const fAdd)
{
  INT nel,lNeq;
  DOUBLE t0,t00,tmp1;
  DOUBLE dt=ddt.dt[0],dt0 = ddt.dt[1];
  short type = ddt.type;
  short j,nD = DENSITY_LEVEL;

  /*...*/
  switch (type) {
/*... EULER de primeira ordem*/
    case EULER:
/*... acumula em f*/
      if (fAdd) {
        for (j = 0; j< ndf; j++) {
          for (nel = 0; nel < numel; nel++) {
            lNeq = id[nel] - 1;
            if (lNeq > -1) {
              t0 = MAT2D(nel,TIME_N_MINUS_1,density,nD)*MAT2D(nel,j,u,ndf);
              MAT2D(j, lNeq, f, nEq) += volume[nel]*t0/dt;
            }
          }
        } 
/*...................................................................*/
      }
/*...................................................................*/

/*... sobrecreve f*/
      else {
        for (j = 0; j< ndf; j++) {
          for (nel = 0; nel < numel; nel++) {
            lNeq = id[nel] - 1;
            if (lNeq > -1) {
              t0 = MAT2D(nel,TIME_N_MINUS_1,density,nD)*MAT2D(nel,j,u,ndf);
              MAT2D(j, lNeq, f, nEq) = volume[nel]*t0/dt;
            }
          }
        }
/*...................................................................*/
      }
/*...................................................................*/
    break;
/*...................................................................*/

/*... BACKWARD de segunda ordem*/
    case BACKWARD:
/*... acumula em f*/
      if (fAdd) {
        for (j = 0; j< ndf; j++) {
          for (nel = 0; nel < numel; nel++) {
            lNeq = id[nel] - 1;
            if (lNeq > -1) {
              tmp1 = dt + dt0;
/*...*/
              t0 = MAT2D(nel,TIME_N_MINUS_1,density,nD)*MAT2D(nel,j,u,ndf);
              t0 *= (tmp1/(dt*dt0));
/*...................................................................*/

/*...*/
              t00 = MAT2D(nel,TIME_N_MINUS_2,density,nD)*MAT2D(nel,j,u0,ndf);
              t00*= (dt/(dt0*tmp1));
/*...................................................................*/
              MAT2D(j,lNeq,f,nEq) += volume[nel]*(t0-t00);
            }
          }
        }
      }
/*...................................................................*/

/*... sobrecreve f*/
      else {
        for (j = 0; j< ndf; j++) {
          for (nel = 0; nel < numel; nel++) {
            lNeq = id[nel] - 1;
            if (lNeq > -1) {
              tmp1 = dt + dt0;
/*..*/
              t0 = MAT2D(nel, TIME_N_MINUS_1, density, nD)*MAT2D(nel, j, u, ndf);
              t0 *= (tmp1 / (dt*dt0));
/*...................................................................*/

/*..*/
              t00= MAT2D(nel, TIME_N_MINUS_1, density, nD)*MAT2D(nel, j, u0, ndf);
              t00*= (dt / (dt0*tmp1));
/*...................................................................*/
              MAT2D(j, lNeq, f, nEq) = volume[nel] * (t0 - t00);
            }
          }
        }
      }
/*...................................................................*/
    break;
  }
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 16/09/2016                                   *
 * Data de modificaco : 02/09/2017                                   *
 *-------------------------------------------------------------------*
 * CELLTRANSIENTPRIME:discretizacao temporal                         *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * volume    -> volume das celulas                                   *
 * u0        -> solucao (n-1)                                        *
 * u         -> solucao (n)                                          *
 * density   -> massa especifica no tempo   (n)                      *
 * f         -> vetor de forcas                                      *
 * ddt       -> delta t                                              *
 * nEq       -> numero de equacoes                                   *
 * numel     -> numero de elementos                                  *
 * ndf       -> graus de liberade                                    *
 * fAdd      -> acumula o valor em f                                 *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * f -> atualizado com a discreticao temporal ( fAdd = true)         *
 *      sobreescreve com a discreticao temporal (fAdd = false)       *
 *-------------------------------------------------------------------*
 * OBS: f | bu(1) bv(1) bw(1) |                                      *
 *        |       ...         |                                      *
 *        | bu(n) bv(n) bw(n) |                                      *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellTransientPrime(DOUBLE *RESTRICT volume 
                       ,DOUBLE *RESTRICT u0     ,DOUBLE *RESTRICT u
                       ,DOUBLE *RESTRICT density,DOUBLE *RESTRICT f
                       ,Temporal const ddt      
                       ,INT const numel         ,short const ndf
                       ,bool const fAdd) 
{
  INT nel;
  DOUBLE t0,t00,tmp1;
  DOUBLE dt=ddt.dt[0],dt0 = ddt.dt[1];
  short type = ddt.type;
  short j,nD = DENSITY_LEVEL;

  /*...*/
  switch (type) {
/*... EULER de primeira ordem*/
    case EULER:
/*... acumula em f*/
      if (fAdd) {
        for (nel = 0; nel < numel; nel++) { 
          for (j = 0; j< ndf; j++) {
            t0 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
            MAT2D(nel,j,f,ndf) += volume[nel]*t0/dt;
          }
        }
      } 
/*...................................................................*/

/*... sobrecreve f*/
      else {
       for (nel = 0; nel < numel; nel++) { 
          for (j = 0; j< ndf; j++) {
            t0 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
            MAT2D(nel,j,f,ndf) = volume[nel]*t0/dt;
          }
        }
      }
/*...................................................................*/
    break;
/*...................................................................*/

/*... BACKWARD de segunda ordem*/
    case BACKWARD:
/*... acumula em f*/
      if (fAdd) {
        for (nel = 0; nel < numel; nel++) {
          for (j = 0; j< ndf; j++) {
            tmp1 = dt + dt0;
/*...*/
            t0 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
            t0 *= (tmp1/(dt*dt0));
/*...................................................................*/

/*...*/
            t00 = MAT2D(nel,1,density,nD)*MAT2D(nel,j,u0,ndf);
            t00*= (dt/(dt0*tmp1));
/*...................................................................*/
            MAT2D(nel,j,f,ndf) += volume[nel]*(t0-t00);
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/

/*... sobrecreve f*/
      else {
        for (nel = 0; nel < numel; nel++) {
          for (j = 0; j< ndf; j++) {
            tmp1 = dt + dt0;
/*..*/
            t0 = MAT2D(nel, 0, density, nD)*MAT2D(nel, j, u, ndf);
            t0 *= (tmp1 / (dt*dt0));
/*...................................................................*/

/*..*/
            t00= MAT2D(nel, 1, density, nD)*MAT2D(nel, j, u0, ndf);
            t00*= (dt / (dt0*tmp1));
/*...................................................................*/
            MAT2D(nel,j,f,ndf) = volume[nel]*(t0 - t00);
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/
    break;
  }
}
/*********************************************************************/


