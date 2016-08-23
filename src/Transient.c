#include<Transient.h>  
/********************************************************************* 
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
* dt        -> delta t                                              *
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
void cellTransient(DOUBLE *restrict volume, INT *restrict id
                   , DOUBLE *restrict u0, DOUBLE *restrict u
                   , DOUBLE *restrict density, DOUBLE *restrict f
                   , DOUBLE const dt
                   , INT const numel, short const ndf
                   , short const type, bool const fAdd)
{
  INT nel, lNeq;
  DOUBLE t1, t2;
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
            t1 = 2.0e0*MAT2D(nel, 0, density, nD)*u[nel];
            t2 = 0.5e0*MAT2D(nel, 1, density, nD)*u0[nel];
            f[lNeq] += volume[nel] * (t1 - t2) / dt;
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
              t1 = 2.0e0*MAT2D(nel, 0, density, nD)*MAT2D(nel, j, u, ndf);
              t2 = 0.5e0*MAT2D(nel, 1, density, nD)*MAT2D(nel, j, u0, ndf);
              MAT2D(lNeq, j, f, ndf) += volume[nel] * (t1 - t2) / dt;
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
            t1 = 2.0e0*MAT2D(nel, 0, density, nD)*u[nel];
            t2 = 0.5e0*MAT2D(nel, 1, density, nD)*u0[nel];
            f[lNeq] = volume[nel] * (t1 - t2) / dt;
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
              t1 = 2.0e0*MAT2D(nel, 0, density, nD)*MAT2D(nel, j, u, ndf);
              t2 = 0.5e0*MAT2D(nel, 1, density, nD)*MAT2D(nel, j, u0, ndf);
              MAT2D(lNeq, j, f, ndf) = volume[nel] * (t1 - t2) / dt;
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
* dt        -> delta t                                              *
* nEq       -> numero de equacoes                                   *
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
* OBS: f | bu(1) bu(2) ... bu(neq) |                                *
*        | bv(1) bv(2) ... bv(neq) |                                *
*        | bw(1) bw(2) ... bw(neq) |                                *
*-------------------------------------------------------------------*
*********************************************************************/
void cellTransientSimple(DOUBLE *restrict volume, INT *restrict id
                         , DOUBLE *restrict u0, DOUBLE *restrict u
                         , DOUBLE *restrict density, DOUBLE *restrict f
                         , DOUBLE const dt, INT const nEq
                         , INT const numel, short const ndf
                         , short const type, bool const fAdd)
{
  INT nel, lNeq;
  DOUBLE t1, t2;
  short j, nD = DENSITY_LEVEL;

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
            t1 = 2.0e0*MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
            t2 = 0.5e0*MAT2D(nel,1,density,nD)*MAT2D(nel,j,u0,ndf);
            MAT2D(j,lNeq,f,nEq) += volume[nel]*(t1-t2)/dt;
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
            t1 = 2.0e0*MAT2D(nel,0,density,nD)*MAT2D(nel,j,u,ndf);
            t2 = 0.5e0*MAT2D(nel,1,density,nD)*MAT2D(nel,j,u0,ndf);
            MAT2D(j,lNeq,f,nEq) = volume[nel]*(t1-t2)/dt;
          }
        }
      }
    }
/*...................................................................*/
    break;
  }
}
/*********************************************************************/



