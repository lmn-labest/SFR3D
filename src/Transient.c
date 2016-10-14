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
void cellTransient(DOUBLE *restrict volume  ,INT *restrict id
                  ,DOUBLE *restrict u0     ,DOUBLE *restrict u
                  ,DOUBLE *restrict density,DOUBLE *restrict f
                  ,Temporal const ddt      ,INT const numel
                  ,short const ndf         ,bool const fAdd)
{
  INT nel, lNeq;
  DOUBLE t1, t2, tmp1;
  DOUBLE dt = ddt.dt[0], dt0 = ddt.dt[1];
  short type = ddt.type;
  short j, nD = DENSITY_LEVEL;

/*...*/
  switch (type) {
/*... EULER de primeira ordem*/
    case EULER:
/*... acumula em f*/
      if (fAdd) {
/*... ndf = 1*/
        if (ndf == 1) {
          for (nel = 0; nel < numel; nel++) {
            lNeq = id[nel] - 1;
            if (lNeq > -1) {
              f[lNeq] += volume[nel]*MAT2D(nel,0,density,nD)*u[nel]/dt;
            }
          }
        }
/*...................................................................*/

/*... ndf > 1*/
        else {
          for (nel = 0; nel < numel; nel++) {
            for (j = 0; j< ndf; j++) {
              lNeq = MAT2D(nel, j, id, ndf) - 1;
              if (lNeq > -1) {
                t1 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
                MAT2D(lNeq, j, f, ndf) += volume[nel] * t1 / dt;
              }
            }
          }
        }
/*...................................................................*/
      }
/*...................................................................*/

/*... sobrecreve f*/
      else {
/*... ndf = 1*/
        if (ndf == 1) {
          for (nel = 0; nel < numel; nel++) {
            lNeq = id[nel] - 1;
            if (lNeq > -1)
              f[lNeq] = volume[nel]*MAT2D(nel,0,density,nD)*u[nel]/dt;
          }
        }
/*...................................................................*/

/*... ndf > 1*/
        else {
          for (nel = 0; nel < numel; nel++) {
            for (j = 0; j< ndf; j++) {
              lNeq = MAT2D(nel, j, id, ndf) - 1;
              if (lNeq > -1) {
                t1 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
                MAT2D(lNeq,j,f,ndf) = volume[nel]*t1/dt;
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
    if (fAdd) {
/*... ndf = 1*/
      if (ndf == 1) {
        for (nel = 0; nel < numel; nel++) {
          lNeq = id[nel] - 1;
          if (lNeq > -1) {
            tmp1 = dt + dt0;
/*...*/
            t1 = MAT2D(nel, 0, density, nD)*u[nel];
            t1 *= (tmp1/(dt*dt0));
/*...................................................................*/

/*...*/
            t2 = MAT2D(nel, 1, density, nD)*u0[nel];
            t2 *= (dt/(dt0*tmp1));
/*...................................................................*/
            f[lNeq] += volume[nel]*(t1-t2);
          }
        }
      }
/*...................................................................*/

/*... ndf > 1*/
      else {
        for (nel = 0; nel < numel; nel++) {
          for (j = 0; j< ndf; j++) {
            lNeq = MAT2D(nel, j, id, ndf) - 1;
            if (lNeq > -1) {
              tmp1 = dt + dt0;
/*...*/
              t1 = MAT2D(nel, 0, density, nD)*MAT2D(nel, j, u, ndf);
              t1 *= (tmp1 / (dt*dt0));
/*...................................................................*/

/*...*/
              t2 = MAT2D(nel, 1, density, nD)*MAT2D(nel, j, u0, ndf);
              t2 *= (dt / (dt0*tmp1));
/*...................................................................*/
              MAT2D(lNeq,j,f,ndf) += volume[nel] * (t1 - t2);
            }
          }
        }
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... sobrecreve f*/
    else {
/*... ndf = 1*/
      if (ndf == 1) {
        for (nel = 0; nel < numel; nel++) {
          lNeq = id[nel] - 1;
          if (lNeq > -1) {
/*...*/
            t1 = MAT2D(nel, 0, density, nD)*u[nel];
            t1 *= (tmp1 / (dt*dt0));
/*...................................................................*/

/*...*/
            t2 = MAT2D(nel, 1, density, nD)*u0[nel];
            t2 *= (dt / (dt0*tmp1));
/*...................................................................*/
            f[lNeq] = volume[nel] *(t1-t2);
          }
        }
      }
/*...................................................................*/

/*... ndf > 1*/
      else {
        for (nel = 0; nel < numel; nel++) {
          for (j = 0; j< ndf; j++) {
            lNeq = MAT2D(nel, j, id, ndf) - 1;
            if (lNeq > -1) {
              tmp1 = dt + dt0;
/*...*/
              t1 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
              t1 *= (tmp1 / (dt*dt0));
/*...................................................................*/

/*...*/
              t2 = MAT2D(nel,1,density,nD)*MAT2D(nel,j,u0,ndf);
              t2 *= (dt / (dt0*tmp1));
/*...................................................................*/
              MAT2D(lNeq,j,f,ndf) = volume[nel]*(t1-t2);
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
void cellTransientSimple(DOUBLE *restrict volume ,INT *restrict id
                        ,DOUBLE *restrict u0     ,DOUBLE *restrict u
                        ,DOUBLE *restrict density,DOUBLE *restrict f
                        ,Temporal const ddt      ,INT const nEq
                        ,INT const numel         ,short const ndf
                        ,bool const fAdd)
{
  INT nel,lNeq;
  DOUBLE t1,t2,tmp1;
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
              t1 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
              MAT2D(j, lNeq, f, nEq) += volume[nel]*t1/dt;
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
              t1 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
              MAT2D(j, lNeq, f, nEq) = volume[nel]*t1/dt;
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
              t1 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
              t1 *= (tmp1/(dt*dt0));
/*...................................................................*/

/*...*/
              t2  = MAT2D(nel,1,density,nD)*MAT2D(nel,j,u0,ndf);
              t2 *= (dt/(dt0*tmp1));
/*...................................................................*/
              MAT2D(j,lNeq,f,nEq) += volume[nel]*(t1-t2);
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
              t1 = MAT2D(nel, 0, density, nD)*MAT2D(nel, j, u, ndf);
              t1 *= (tmp1 / (dt*dt0));
/*...................................................................*/

/*..*/
              t2 = MAT2D(nel, 1, density, nD)*MAT2D(nel, j, u0, ndf);
              t2 *= (dt / (dt0*tmp1));
/*...................................................................*/
              MAT2D(j, lNeq, f, nEq) = volume[nel] * (t1 - t2);
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
 * Data de modificaco : 00/00/0000                                   *
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
void cellTransientPrime(DOUBLE *restrict volume 
                       ,DOUBLE *restrict u0     ,DOUBLE *restrict u
                       ,DOUBLE *restrict density,DOUBLE *restrict f
                       ,Temporal const ddt      
                       ,INT const numel         ,short const ndf
                       ,bool const fAdd) 
{
  INT nel;
  DOUBLE t1,t2,tmp1;
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
            t1 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
            MAT2D(nel,j,f,ndf) += volume[nel]*t1/dt;
          }
        }
      } 
/*...................................................................*/

/*... sobrecreve f*/
      else {
       for (nel = 0; nel < numel; nel++) { 
          for (j = 0; j< ndf; j++) {
            t1 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
            MAT2D(nel,j,f,ndf) = volume[nel]*t1/dt;
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
            t1 = MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
            t1 *= (tmp1/(dt*dt0));
/*...................................................................*/

/*...*/
            t2  = MAT2D(nel,1,density,nD)*MAT2D(nel,j,u0,ndf);
            t2 *= (dt/(dt0*tmp1));
/*...................................................................*/
            MAT2D(nel,j,f,ndf) += volume[nel]*(t1-t2);
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
            t1 = MAT2D(nel, 0, density, nD)*MAT2D(nel, j, u, ndf);
            t1 *= (tmp1 / (dt*dt0));
/*...................................................................*/

/*..*/
            t2 = MAT2D(nel, 1, density, nD)*MAT2D(nel, j, u0, ndf);
            t2 *= (dt / (dt0*tmp1));
/*...................................................................*/
            MAT2D(nel,j,f,ndf) = volume[nel]*(t1 - t2);
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


