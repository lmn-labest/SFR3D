#include<TimeInterpol.h>

/*********************************************************************
 * Data de criacao    : 14/11/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * interPolTime: interpolação linear entre dois passos de tempo      *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void interPolTime(DOUBLE *RESTRICT ui ,DOUBLE *RESTRICT u
                 ,DOUBLE *RESTRICT u0 ,DOUBLE const ts
                 ,DOUBLE const t1     ,DOUBLE const t0
                 ,INT const nl        ,short const nc
                 ,bool const fTimePlot,bool const fStepPlot)
{
  short j;
  INT i; 
  DOUBLE f1,f2;

  if(fTimePlot)
  {
    f1 = (ts-t0)/(t1-t0);
    f2 = 1.e0 -f1;
    for(i=0;i<nl;i++)
      for(j=0;j<nc;j++)
       MAT2D(i,j,ui,nc) = f1*MAT2D(i,j,u,nc) + f2*MAT2D(i,j,u0,nc);
  }

  else if(fStepPlot)
  {
    for(i=0;i<nl;i++)
      for(j=0;j<nc;j++)
       MAT2D(i,j,ui,nc) = MAT2D(i,j,u,nc);
  }

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 14/11/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * initTimeStruct:                                                   *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void initTimeStruct(Memoria *m        ,TimeInterpol *ti
                   ,Mesh *mesh0       ,Mesh *mesh
                   ,Combustion *cModel,FileOpt *opt)
{
  short ndm = mesh->ndm
       ,ns  = cModel->nOfSpecies 
       ,nD  = DENSITY_LEVEL;
  INT nel = mesh->numelNov
     ,nelG = mesh->numel;

/*... vel*/ 
  if(opt->vel)
  {

    if(opt->fTimePlot)
    {
      ti->vel0 = mesh->elm.vel0;
      HccaAlloc(DOUBLE,m,ti->vel0,nel*ndm,"ivel0"  ,_AD_);
      ti->vel  = mesh->elm.vel;
    }  
/*...................................................................*/

/*...*/
    else if(opt->fStepPlot)
    {
      ti->vel0 = ti->vel = NULL;
    }
    HccaAlloc(DOUBLE,m,ti->veli,nel*ndm,"iveli"  ,_AD_);
    ti->velG = ti->veli;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->velG,nelG*ndm,"ivelG"  ,_AD_);
  } 
/*...................................................................*/

/*... pres*/
  if(opt->pres)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->p0,nel,"p0"  ,_AD_);
      ti->p0 = mesh->elm.pressure0;
      ti->p  = mesh->elm.pressure;
    }  
    else if(opt->fStepPlot)
    {
      ti->p0 = ti->p = NULL;
    }
    HccaAlloc(DOUBLE,m,ti->pi,nel,"pi"  ,_AD_);
    ti->pG = ti->pi;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->pG,nelG,"pG"  ,_AD_);
  } 
/*...................................................................*/

/*... temp*/
  if(opt->temp)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->temp0,nel,"itemp0"  ,_AD_);
      ti->temp  = mesh->elm.temp;
    }  
    else if(opt->fStepPlot)
    {
      ti->temp0 = ti->temp = NULL;
    }
    HccaAlloc(DOUBLE,m,ti->tempi,nel,"tempi"  ,_AD_);
    ti->tempG = ti->tempi;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->tempG,nelG,"tempG"  ,_AD_);
  } 
/*...................................................................*/

/*... yFrac*/
  if(opt->yFrac)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->y0,nel*ns,"y0"  ,_AD_);
      ti->y  = mesh->elm.yFrac;
    }  
    else if(opt->fStepPlot)
    {
      ti->y0 = ti->y = NULL;
    }
    HccaAlloc(DOUBLE,m,ti->yi,nel*ns,"yi"  ,_AD_);
    ti->yG = ti->yi;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->yG,nelG*ns,"yG"  ,_AD_);
  } 
/*...................................................................*/

/*... wT*/
  if(opt->rateHeatComb)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->wT0,nel,"tIwT0"  ,_AD_);
      ti->wT  = mesh->elm.rateHeatReComb;
    }  
    else if(opt->fStepPlot)
    {
      ti->wT0 = ti->wT = NULL;
    }
    HccaAlloc(DOUBLE,m,ti->wTi,nel,"tIwTi",_AD_);
    ti->wTG = ti->wTi;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->wTG,nelG,"tIwTG"  ,_AD_);
  } 
/*...................................................................*/

/*... densityFluid*/
  //if(opt->densityFluid)
 // {
//  if(opt->fTimePlot)
//  {
//    HccaAlloc(DOUBLE,m,ti->rho,nel,"irho0"  ,_AD_);
//    ti->wT  = &mesh->elm.densityFluid[nel*nD];
//  }  
//  else if(opt->fStepPlot)
 // {
 //   ti->rho0 = ti->rho = NULL;
 // }
 // HccaAlloc(DOUBLE,m,ti->wTi,nel,"irho0",_AD_);
 // ti->rhoG = ti->rhoi;
 // if(!mpiVar.myId && mpiVar.nPrcs>1)
 //   HccaAlloc(DOUBLE,m,ti->rhoG,nelG,"irhoG"  ,_AD_);
//} 
/*...................................................................*/

/*... dVisc*/
  if(opt->dViscosity)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->dVisc0,nel,"dVisc0"  ,_AD_);
      ti->dVisc = mesh->elm.dViscosity;
    }  
    else if(opt->fStepPlot)
    {
      ti->dVisc0 = ti->dVisc = NULL;
    }
    HccaAlloc(DOUBLE,m,ti->dVisci,nel,"dVisci",_AD_);
    ti->dViscG = ti->dVisci;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->dViscG,nelG,"dViscG"  ,_AD_);
  } 
/*...................................................................*/

/*... tCond*/
  if(opt->tConductivity)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->tCond0,nel,"tCond0"  ,_AD_);
      ti->tCond = mesh->elm.tConductivity;
    }  
    else if(opt->fStepPlot)
    {
      ti->tCond0 = ti->tCond = NULL;
    }
    HccaAlloc(DOUBLE,m,ti->tCondi,nel,"tCondi",_AD_);
    ti->tCondG = ti->tCondi;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->tCondG,nelG,"tCondG"  ,_AD_);
  } 
/*...................................................................*/

/*... cDiff*/
  if(opt->coefDiffSp)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->cDiff0,nel*ns,"cDiff0"  ,_AD_);
      ti->cDiff = mesh->elm.cDiffComb;
    }  
    else if(opt->fStepPlot)
    {
      ti->cDiff0 = ti->cDiff = NULL;
    }
    HccaAlloc(DOUBLE,m,ti->cDiffi,nel*ns,"cDiffi",_AD_);
    ti->cDiffG = ti->cDiffi;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->cDiffG,nelG*ns,"cDiffG"  ,_AD_);
  } 
/*...................................................................*/


}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 14/11/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * initTimeStruct:                                                   *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void updateTimeStruct(Memoria *m        ,TimeInterpol *ti
                     ,Mesh *mesh
                     ,Combustion *cModel,FileOpt *opt)
{
  short j
       ,ndm = mesh->ndm
       ,ns  = cModel->nOfSpecies;
  INT nel = mesh->numelNov
     ,i;

  if(opt->fStepPlot) 
    return;

/*... vel*/
  if(opt->vel)
  {
    for(i=0;i<nel;i++)
      for(j=0;j<ndm;j++)
        MAT2D(i,j,ti->vel0,ndm) = MAT2D(i,j,mesh->elm.vel,ndm);
  } 
/*...................................................................*/

/*... pres*/
  if(opt->pres)
  {
    for(i=0;i<nel;i++)
      ti->p0[i] = mesh->elm.pressure[i];
  } 
/*...................................................................*/

/*... temp*/
  if(opt->temp)
  {
    for(i=0;i<nel;i++)
      ti->temp0[i] = mesh->elm.temp[i];
  } 
/*...................................................................*/

/*... yFrac*/
  if(opt->yFrac)
  {
    for(i=0;i<nel;i++)
      for(j=0;j<ns;j++)
        MAT2D(i,j,ti->y0,ns) = MAT2D(i,j,mesh->elm.yFrac,ns);
  } 
/*...................................................................*/

/*... wT*/
  if(opt->rateHeatComb)
  {
    for(i=0;i<nel;i++)
      ti->wT0[i] = mesh->elm.rateHeatReComb[i];
  } 
/*...................................................................*/

/*... dVisc*/
  if(opt->dViscosity)
  {
    for(i=0;i<nel;i++)
      ti->dVisc0[i] = mesh->elm.dViscosity[i];
  } 
/*...................................................................*/

/*... tCond*/
  if(opt->tConductivity)
  {
    for(i=0;i<nel;i++)
      ti->tCond0[i] = mesh->elm.tConductivity[i];
  } 
/*...................................................................*/

/*... cDiff*/
  if(opt->coefDiffSp)
  {
    for(i=0;i<nel;i++)
      for(j=0;j<ns;j++)
        MAT2D(i,j,ti->cDiff0,ns) = MAT2D(i,j,mesh->elm.cDiffComb,ns);
  } 
/*...................................................................*/

}
/*********************************************************************/