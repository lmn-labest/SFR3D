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
       ,ns  = cModel->nOfSpecies; 
  INT nel =  mesh->numel
     ,nelG  = mesh0->numel
     ,nno  = mesh->nnode
     ,nnoG = mesh0->nnode;

/*... vel*/ 
  if(opt->vel || opt->gradVel)
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
      ti->vel0 = ti->vel = mesh->elm.vel;
    }
    
    HccaAlloc(DOUBLE,m,ti->velI,nel*ndm,"iveli"  ,_AD_);
    ti->velG = ti->velI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)   
      HccaAlloc(DOUBLE,m,ti->velG,nelG*ndm,"ivelG"  ,_AD_);      
    
    HccaAlloc(DOUBLE,m,ti->nVelI,nno*ndm,"nveli"  ,_AD_);
    ti->nVelG = ti->nVelI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)   
      HccaAlloc(DOUBLE,m,ti->nVelG,nnoG*ndm,"nvelG"  ,_AD_);

  } 
/*...................................................................*/

/*... pres*/
  if(opt->pres || opt->gradPres )
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->p0,nel,"p0"  ,_AD_);
      ti->p  = mesh->elm.pressure;
    }  
    else if(opt->fStepPlot)
    {
      ti->p0 = ti->p = mesh->elm.pressure;
    }
    HccaAlloc(DOUBLE,m,ti->pI,nel,"pi"  ,_AD_);
    ti->pG = ti->pI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->pG,nelG,"pG"  ,_AD_);

    HccaAlloc(DOUBLE,m,ti->nPresI,nno,"npi"  ,_AD_);
    ti->nPresG = ti->nPresI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)   
      HccaAlloc(DOUBLE,m,ti->nPresG,nnoG,"npG"  ,_AD_);

  } 
/*...................................................................*/

/*... temp*/
  if(opt->temp || opt->gradTemp)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->temp0,nel,"itemp0"  ,_AD_);
      ti->temp  = mesh->elm.temp;
    }  
    else if(opt->fStepPlot)
    {
      ti->temp0 = ti->temp = mesh->elm.temp;
    }
    HccaAlloc(DOUBLE,m,ti->tempI,nel,"tempi"  ,_AD_);
    ti->tempG = ti->tempI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->tempG,nelG,"tempG"  ,_AD_);
    
    HccaAlloc(DOUBLE,m,ti->nTempI,nno,"nTempI"  ,_AD_);
    ti->nTempG = ti->nTempI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)   
      HccaAlloc(DOUBLE,m,ti->nTempG,nnoG,"nTempG"  ,_AD_);

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
      ti->y0 = ti->y = mesh->elm.yFrac;
    }
    HccaAlloc(DOUBLE,m,ti->yI,nel*ns,"yi"  ,_AD_);
    ti->yG = ti->yI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->yG,nelG*ns,"yG"  ,_AD_);

    HccaAlloc(DOUBLE,m,ti->nYI,nno*ns,"nYI"  ,_AD_);
    ti->nYG = ti->nYI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)   
      HccaAlloc(DOUBLE,m,ti->nYG,nnoG*ns,"nYG"  ,_AD_);

  } 
/*...................................................................*/

/*... wT*/
  if(opt->wT)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->wT0,nel,"tIwT0"  ,_AD_);
      ti->wT  = mesh->elm.wT;
    }  
    else if(opt->fStepPlot)
    {
      ti->wT0 = ti->wT = mesh->elm.wT;
    }
    HccaAlloc(DOUBLE,m,ti->wTI,nel,"tIwTi",_AD_);
    ti->wTG = ti->wTI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->wTG,nelG,"tIwTG"  ,_AD_);

    HccaAlloc(DOUBLE,m,ti->nWTI,nno,"nWGI",_AD_);
    ti->nWTG = ti->nWTI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->nWTG,nnoG,"nWTG"  ,_AD_);
  } 
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
      ti->dVisc0 = ti->dVisc = mesh->elm.dViscosity;
    }
    HccaAlloc(DOUBLE,m,ti->dViscI,nel,"dVisci",_AD_);
    ti->dViscG = ti->dViscI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->dViscG,nelG,"dViscG"  ,_AD_);

    HccaAlloc(DOUBLE,m,ti->ndViscI,nno,"ndViscI",_AD_);
    ti->ndViscG = ti->ndViscI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->ndViscG,nnoG,"ndViscG"  ,_AD_);

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
      ti->tCond0 = ti->tCond = mesh->elm.tConductivity;
    }
    HccaAlloc(DOUBLE,m,ti->tCondI,nel,"tCondi",_AD_);
    ti->tCondG = ti->tCondI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->tCondG,nelG,"tCondG"  ,_AD_);

    HccaAlloc(DOUBLE,m,ti->ntCondI,nno,"ntCondI",_AD_);
    ti->ntCondG = ti->ntCondI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->ntCondG,nnoG,"ntCondG"  ,_AD_);

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
      ti->cDiff0 = ti->cDiff = mesh->elm.cDiffComb;
    }
    HccaAlloc(DOUBLE,m,ti->cDiffI,nel*ns,"cDiffi",_AD_);
    ti->cDiffG = ti->cDiffI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->cDiffG,nelG*ns,"cDiffG"  ,_AD_);

    HccaAlloc(DOUBLE,m,ti->ncDiffI,nno*ns,"ncDiffIi",_AD_);
    ti->ncDiffG = ti->ncDiffI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->ncDiffG,nnoG*ns,"ncDiffG"  ,_AD_);

  } 
/*...................................................................*/

/*... sHeat*/
  if(opt->specificHeat)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->sHeat0,nel,"iSheat0"  ,_AD_);
      ti->sHeat = mesh->elm.specificHeat.t;
    }  
    else if(opt->fStepPlot)
    {
      ti->sHeat0 = ti->sHeat = mesh->elm.specificHeat.t;
    }
    HccaAlloc(DOUBLE,m,ti->sHeatI,nel,"iSheati",_AD_);
    ti->sHeatG = ti->sHeatI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->sHeatG,nelG,"iSheatG"  ,_AD_);

     HccaAlloc(DOUBLE,m,ti->nsHeatI,nno,"nsHeatI",_AD_);
    ti->nsHeatG = ti->nsHeatI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->nsHeatG,nnoG,"nsHeatG"  ,_AD_);

  } 
/*...................................................................*/

/*... density*/
  if(opt->densityFluid)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->rho0,nel,"iRhot0"  ,_AD_);
      ti->rho = mesh->elm.densityFluid.t;
    }  
    else if(opt->fStepPlot)
    {
      ti->rho0 = ti->rho = mesh->elm.densityFluid.t;
    }

    HccaAlloc(DOUBLE,m,ti->rhoI,nel,"iRhoi",_AD_);
    ti->rhoG = ti->rhoI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->rhoG,nelG,"iRhoG"  ,_AD_);

    HccaAlloc(DOUBLE,m,ti->nRhoI,nno,"nRhoI",_AD_);
    ti->nRhoG = ti->nRhoI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->nRhoG,nnoG,"nRhoG"  ,_AD_);


  } 
/*...................................................................*/

/*... molarMass*/
  if(opt->mMolar)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->mMolar0,nel,"iMmolar0"  ,_AD_);
      ti->mMolar = mesh->elm.mMolar.t;
    }  
    else if(opt->fStepPlot)
    {
      ti->mMolar0 = ti->mMolar = mesh->elm.mMolar.t;
    }
    HccaAlloc(DOUBLE,m,ti->mMolarI,nel,"iMmolari",_AD_);
    ti->mMolarG = ti->mMolarI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->mMolarG,nelG,"iMmolarG"  ,_AD_);

    HccaAlloc(DOUBLE,m,ti->nmMolarI,nno,"nMmolarI",_AD_);
    ti->nmMolarG = ti->nmMolarI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->nmMolarG,nnoG,"nMmolarG"  ,_AD_);

  } 
/*...................................................................*/

/*... gradTemp*/
  if(opt->gradTemp)
  {
    ti->gradTempI = mesh->elm.gradTemp;
    ti->gradTempG = ti->gradTempI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->gradTempG,nelG*ndm,"iGradTG"  ,_AD_);

    HccaAlloc(DOUBLE,m,ti->nGradTempI,nno*ndm,"nGradTempI"  ,_AD_);
    ti->nGradTempG = ti->nGradTempI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)   
      HccaAlloc(DOUBLE,m,ti->nGradTempG,nnoG*ndm,"nGradTempG"  ,_AD_);

  } 
/*...................................................................*/

/*... gradPres*/
  if(opt->gradPres)
  {
    ti->gradPresI = mesh->elm.gradPres;
    ti->gradPresG = ti->gradPresI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->gradPresG,nelG*ndm,"iGradPresG"  ,_AD_);

    HccaAlloc(DOUBLE,m,ti->nGradPresI,nno*ndm,"nGradPresI"  ,_AD_);
    ti->nGradPresG = ti->nGradPresI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)   
      HccaAlloc(DOUBLE,m,ti->nGradPresG,nnoG*ndm,"nGradPresG"  ,_AD_);

  } 
/*...................................................................*/

/*... gradVel*/
  if(opt->gradVel)
  {
    ti->gradVelI = mesh->elm.gradVel;
    ti->gradVelG = ti->gradVelI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->gradVelG,nelG*ndm*ndm,"iGradVelG"  ,_AD_);
    
    HccaAlloc(DOUBLE,m,ti->nGradVelI,nno*ndm*ndm,"nGradVeli"  ,_AD_);
    ti->nGradVelG = ti->nGradVelI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)   
      HccaAlloc(DOUBLE,m,ti->nGradVelG,nnoG*ndm*ndm,"nGradVelG"  ,_AD_);  
  } 
/*...................................................................*/


/*... gradY*/
  if(opt->gradY)
  {
    HccaAlloc(DOUBLE,m,ti->gradYI,nel*ndm*ns,"iGradYI"  ,_AD_);
    ti->gradYG = ti->gradYI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->gradYG,nelG*ndm*ns,"iGradYG"  ,_AD_);
    
    HccaAlloc(DOUBLE,m,ti->nGradYI,nno*ndm*ns,"nGradYI"  ,_AD_);
    ti->nGradYG = ti->nGradYI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)   
      HccaAlloc(DOUBLE,m,ti->nGradVelG,nnoG*ndm*ns,"nGradYG"  ,_AD_);  
  } 
/*...................................................................*/


/*... eddyViscosity*/
  if(opt->eddyViscosity)
  {
    if(opt->fTimePlot)
    {
      HccaAlloc(DOUBLE,m,ti->eddyVisc0,nel,"eddyVisc0"  ,_AD_);
      ti->eddyVisc = mesh->elm.eddyViscosity;
    }  
    else if(opt->fStepPlot)
    {
      ti->eddyVisc = ti->dVisc = mesh->elm.eddyViscosity;
    }
    HccaAlloc(DOUBLE,m,ti->eddyViscI,nel,"eddyVisci",_AD_);
    ti->eddyViscG = ti->eddyViscI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->eddyViscG,nelG,"eddyViscG"  ,_AD_);

    HccaAlloc(DOUBLE,m,ti->neddyViscI,nno,"neddydViscI",_AD_);
    ti->neddyViscG = ti->neddyViscI;
    if(!mpiVar.myId && mpiVar.nPrcs>1)
      HccaAlloc(DOUBLE,m,ti->neddyViscG,nnoG,"neddyViscG"  ,_AD_);

  } 
/*...................................................................*/

/*... eddyViscosity*/
  if(opt->tReactor)
  {
    fprintf(fileLogDebug,"Print: Nao implementado tReactor!!\n");
    exit(-1);
  } 
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 14/11/2019                                   *
 * Data de modificaco : 17/01/2019                                   *
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
  INT nel = mesh->numel
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
  if(opt->wT)
  {
    for(i=0;i<nel;i++)
      ti->wT0[i] = mesh->elm.wT[i];
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

/*... sHeat*/
  if(opt->specificHeat)
  {
    for(i=0;i<nel;i++)
      ti->sHeat0[i] = mesh->elm.specificHeat.t[i];
  } 
/*...................................................................*/

/*... density*/
  if(opt->densityFluid)
  {
    for(i=0;i<nel;i++)
      ti->rho0[i] = mesh->elm.densityFluid.t[i];
  } 
/*...................................................................*/

/*... mMolar*/
  if(opt->mMolar)
  {
    for(i=0;i<nel;i++)
      ti->mMolar0[i] = mesh->elm.mMolar.t[i];
  } 
/*...................................................................*/

/*... eddyVisc*/
  if(opt->eddyViscosity)
  {
    for(i=0;i<nel;i++)
      ti->eddyVisc0[i] = mesh->elm.eddyViscosity[i];
  } 
/*...................................................................*/


}
/*********************************************************************/