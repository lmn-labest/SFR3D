#include<Reaction.h>

/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * arrheniusLaw: Lei de arrhenius                                    *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * T     -> temperatura em kelvin                                    *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*  
 *********************************************************************/
static DOUBLE arrheniusLaw(Arrhenius *arr,DOUBLE const T)
{
  DOUBLE A = arr->A, beta = arr->beta, Ta = arr->Ta;   

  return A*pow(T,beta)*exp(-Ta/T);    

}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 21/07/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * molaRateReaction: taxa de reacao molar                            *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * chem  -> concentracao molar das especies                          *
 * Q     -> taxe de reacao molar o reacao (kmol/m3s)                 *
 * w     -> nao definido                                             *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * w     -> taxa de reacao das especies em kg/m3 s                   *
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*  
 *********************************************************************/
void massRateReaction(Chemical *chem  ,DOUBLE *RESTRICT Q
                     ,DOUBLE *RESTRICT w)
{
  unsigned short i, j, nReac =chem->nReac, nSp = chem->nSp;
  DOUBLE de;
  
  for(j=0;j<nSp;j++)
    w[j] = 0.e0;

  for(i=0;i<nReac;i++)
  { 
    for(j=0;j<nSp;j++)
    {
/*... */
      de = chem->reac[i].stch[2][j];
/*... kg/m3s*/
      w[j] += de*Q[i]*chem->sp[j].mW;
    }
  }
}
/********************************************************************/


/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * massActionMass:                                                   *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * c     -> concentracao molar das especies                          *
 * w     -> nao definido                                             *
 * T     -> temperatura em kelvin                                    *
 * nSp   -> numero de especies                                       *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS: mol/(cm3 s)                                                  *
 * ------------------------------------------------------------------*  
 *********************************************************************/
DOUBLE massActionMass(Reaction *reac     ,Prop *sHeatPol 
                     ,DOUBLE  *RESTRICT c
                     ,DOUBLE const T     ,unsigned short const nSp)
{
  unsigned short i;
  DOUBLE kf,kr,qf,qr;  

/*... reacao direta*/
  kf = arrheniusLaw(&reac->ArrF, T);

  for(i=0,qf=kf;i<nSp;i++)
    qf *= pow(c[i],reac->exp[0][i]);

/*... reacao inversao*/
  qr = 0.e0;
  if(reac->reverse)
  {
/*  kr = reverseK(reac       ,sHeatPol 
                 ,T          ,R
                 ,kf         ,1.0     
                 ,nSp);*/
    kr = arrheniusLaw(&reac->ArrR, T);

    for(i=0,qr=kr;i<nSp;i++)
      qr *= pow(c[i],reac->exp[1][i]);
  }

  return qf-qr;
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 17/08/2019                                   *
 *-------------------------------------------------------------------*
 * rateReaction: calculo da taxa de consumo do combustivel           *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel        -> modelo de combustao                              *
 * zComb         -> fracao massica                                   *
 * temp          -> temperatura                                      *
 * rate          -> nao definido                                     * 
 * density       -> densidade                                        *
 * eddyViscosity -> viscosidade turbulenta                           *
 * dViscosity    -> viscosidade dinamica                             *
 * tReactor      -> paramentro do reator                             *
 * numel         -> numero de elementos                              *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * rate    -> taxa de consumo massico das especies                   * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void rateReaction(Combustion *cModel         , Turbulence *tModel
             , PropVarFluid *pFluid
             , DOUBLE *RESTRICT zComb        , DOUBLE *RESTRICT temp        
             , DOUBLE *RESTRICT rate         , DOUBLE *RESTRICT density 
             , DOUBLE *RESTRICT gradVel      , DOUBLE *RESTRICT eddyViscosity
             , DOUBLE *RESTRICT dViscosity   , DOUBLE *RESTRICT tReactor 
             , DOUBLE const dt               , DOUBLE const Pth 
             , short const ndm               , INT const numel
             , bool const fKelvin            , bool const fOmp           
             , short const nThreads )  
{
  short nComb = cModel->nComb
       , iCod = cModel->reactionKinetic
       , nReac= cModel->chem.nReac
       , nSp  = cModel->nOfSpecies;
  short i;
  INT nel,it=0;
  DOUBLE *pz;
  DOUBLE omega, densityC, eddyC;
  DOUBLE tK, y[MAXSPECIES],cM[MAXSPECIES],Q[MAXREAC],w[MAXSPECIES];

/*...*/
  switch(iCod)
  {
/*...*/
    case ARRHENIUS:
      for(nel = 0; nel < numel; nel++)
      {
/*...*/
        pz = &MAT2D(nel,0,zComb,nComb);
        densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
        getSpeciesPrimitivesCc(cModel,y,pz);
        concetracionOfSpecies(cModel            ,y
                             ,cM                ,densityC);

        TEMP(tK,temp[nel],fKelvin);
/*...................................................................*/

/*... */
        for(i=0;i<nReac;i++)
        {
/*... mol/(cm3 s)*/
          omega = massActionMass(&cModel->chem.reac[i],&pFluid->sHeat 
                                ,cM
                                ,tK                   ,nComb);
/*... kmol/(m3 s)*/
          Q[i] = 1.e+03*omega;
        }
/*...................................................................*/

      massRateReaction(&cModel->chem,Q,w);
/*... */     
      for(i=0;i<nSp;i++)
        MAT2D(nel,i,rate,nSp) = w[i];
/*...................................................................*/
    } 
/*...................................................................*/    
    break;
/*...................................................................*/

/*...*/
    case EDC:
/*... omp*/
      if(fOmp)
      {
/*...*/
#pragma omp parallel  for default(none) num_threads(nThreads)\
        private(nel,pz,i,it,densityC,eddyC,y,w)\
        shared(cModel,zComb,nComb,density,nSp,rate,pFluid,temp\
              ,tReactor,eddyViscosity,dViscosity,thDynamic\
              ,tModel)
        for(nel = 0; nel < numel; nel++)
        {
          pz = &MAT2D(nel,0,zComb,nComb);
          densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
          getSpeciesPrimitivesCc(cModel,y,pz);     
          if(tModel->fTurb)
            eddyC = eddyViscosity[nel];
/*...*/
          it = edc(cModel                         ,pFluid
           ,y                                     ,w
           ,&MAT2D(nel,0,tReactor,N_TERMS_REACTOR),densityC
           ,dt                                    ,temp[nel]
           ,eddyC                                 ,dViscosity[nel]
           ,thDynamic.pTh[2]                      ,fKelvin
           ,nel );
/*...................................................................*/
          for(i=0;i<nSp;i++)
            MAT2D(nel,i,rate,nSp) = w[i];
        }     
/*...................................................................*/
      }
/*...................................................................*/ 

/*...*/
      else
      {
/*...*/
        for(nel = 0; nel < numel; nel++)
        {
          pz = &MAT2D(nel,0,zComb,nComb);
          densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
          getSpeciesPrimitivesCc(cModel,y,pz);
          if(tModel->fTurb)
            eddyC = eddyViscosity[nel];
/*...*/
          it = edc(cModel                         ,pFluid
           ,y                                     ,w
           ,&MAT2D(nel,0,tReactor,N_TERMS_REACTOR),densityC
           ,dt                                    ,temp[nel]
           ,eddyC                                 ,dViscosity[nel]
           ,thDynamic.pTh[2]                      ,fKelvin
           ,nel );
/*...................................................................*/
          for(i=0;i<nSp;i++)
            MAT2D(nel,i,rate,nSp) = w[i];
        }     
/*...................................................................*/
      }
/*...................................................................*/
      break;
/*...................................................................*/

/*...*/
    case EDM:
/*...*/
/*... omp*/
      if(fOmp)
      {
/*...*/
#pragma omp parallel  for default(none) num_threads(nThreads)\
        private(nel,pz,i,it,densityC,y,w)\
        shared(cModel,zComb,nComb,density,nSp,rate,tReactor)
        for(nel = 0; nel < numel; nel++)
        { 
          pz = &MAT2D(nel,0,zComb,nComb);
          densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
          getSpeciesPrimitivesCc(cModel,y,pz);
/*...*/
          edm(cModel,y,w,densityC,MAT2D(nel,0,tReactor,N_TERMS_REACTOR));
/*...................................................................*/
          for(i=0;i<nSp;i++)
            MAT2D(nel,i,rate,nSp) = w[i];
        }
/*...................................................................*/
      }
/*...................................................................*/

/*...*/
      else
      {
/*...*/
        for(nel = 0; nel < numel; nel++)
        { 
          pz = &MAT2D(nel,0,zComb,nComb);
          densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
          getSpeciesPrimitivesCc(cModel,y,pz);
/*...*/
          edm(cModel,y,w,densityC,MAT2D(nel,0,tReactor,N_TERMS_REACTOR));
/*...................................................................*/
          for(i=0;i<nSp;i++)
            MAT2D(nel,i,rate,nSp) = w[i];

        }
/*...................................................................*/
      }
/*...................................................................*/
      break;
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 01/11/2019                                   *
 *-------------------------------------------------------------------*
 * timeChemical: calculo da taxa de consumo do combustivel           *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> modelo de combustao                                    *
 * zFrac   -> fracao massica                                         *
 * temp    -> temperatura                                            *
 * density -> densidade                                              *
 * eddyVis -> viscosidade turbulenta                                 * 
 * cDiffY  -> coeficiente de difusaõ massico                         *
 * dVisc   -> viscosidade dinamica                                   *
 * tReactor-> escala de tempo do reator                              * 
 * numel   -> numero de elementos                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * rate    -> taxa de consumo do combustivel                         * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void timeChemical(Combustion *cModel      , Turbulence *tModel
             , PropVarFluid *pFluid
             , DOUBLE *RESTRICT zComb        , DOUBLE *RESTRICT temp  
             , DOUBLE *RESTRICT density      , DOUBLE *RESTRICT gradVel 
             , DOUBLE *RESTRICT eddyViscosity, DOUBLE *RESTRICT cDiffY
             , DOUBLE *RESTRICT volume  
             , DOUBLE *RESTRICT dViscosity   , DOUBLE *RESTRICT tReactor
             , short const ndm               , INT const numel   
             , bool const fKelvin )
{
  short  nComb = cModel->nComb    
       , nReac=cModel->chem.nReac
       , nSp  = cModel->nOfSpecies;
  short i;
  INT nel;
  DOUBLE sT[6],*iGradVel,*pz;
  DOUBLE omega, densityC,tmp,modS,tMix,dF;
  DOUBLE tTurb,tDiff;
  DOUBLE tK,y[MAXSPECIES],cM[MAXSPECIES],Q[MAXREAC],w[MAXSPECIES],tChemical; 

/*...*/
  for(nel = 0; nel < numel; nel++)
  {
    pz = &MAT2D(nel,0,zComb,nComb);
    densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
    getSpeciesPrimitivesCc(cModel,y,pz);
    concetracionOfSpecies(cModel            ,y
                         ,cM                ,densityC);

    TEMP(tK,temp[nel],fKelvin);

/*.. calculo Sij*/
    iGradVel = &MAT3D(nel,0,0,gradVel,ndm,ndm);
    tensorS(sT,iGradVel,false);
/*... |S| = sqrt(2S:S)*/
    modS = sqrt(2.e0*doubleDotSym(sT)) + 1.0e-32;
/*...................................................................*/

/*...*/
    for(i=0;i<cModel->chem.nReac;i++)
    {
/*... mol/(cm3 s)*/
      omega = massActionMass(&cModel->chem.reac[i],&pFluid->sHeat 
                            ,cM
                            ,tK                   ,nComb);
/*...................................................................*/

/*... kmol/(m3 s)*/
      Q[i] = 1.e+03*omega;
    }
/*...................................................................*/ 

    massRateReaction(&cModel->chem,Q ,w);
    dF = MAT2D(nel,0,cDiffY,nSp);
    for(i=0,tChemical = 1.e+32;i<nSp;i++)
    {
      if(fabs(w[i]) > 1.e-16)
      {
        tmp       = fabs(w[i]/densityC);
        tmp       = y[i]/tmp;
        tChemical = max(tmp,tChemical); 
      }
      dF = min(dF,MAT2D(nel,i,cDiffY,nSp));
    }
/*...................................................................*/ 
      
/*...*/ 
    tDiff = pow(volume[nel], D2DIV3)/dF;
    tTurb = 1.e0/modS;
    tMix  = min(tTurb,tDiff);
/*...................................................................*/

/*...*/  
    MAT2D(nel,0,tReactor,N_TERMS_REACTOR) = tMix;
    MAT2D(nel,1,tReactor,N_TERMS_REACTOR) = tChemical; 
    MAT2D(nel,2,tReactor,N_TERMS_REACTOR) = tMix/tChemical; 
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 12/05/2019                                   *
 *-------------------------------------------------------------------*
 * initLumpedMatrix : inicializa a matriz de especies                *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> modelo de combustao                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void initLumpedMatrix(Combustion *cModel)
{
  
  short nl = cModel->nOfSpeciesLump
      , ns =  cModel->nOfSpecies 
      , posN2 = cModel->chem.sN2;
  DOUBLE mO2,mN2,mCO2p,mH2Op,mN2p,mAir,mProd;

  mO2 = mN2 = mCO2p = mH2Op = mN2p = 0.e0; 

  mO2   = cModel->chem.reac[0].stch[cModel->chem.sO2][0]
         *cModel->chem.sp[cModel->chem.sO2].mW;
  mN2   = cModel->chem.reac[0].stch[posN2][0]
         *cModel->chem.sp[posN2].mW;
/*...*/
  mCO2p = cModel->chem.reac[0].stch[cModel->chem.sCO2][1]
         *cModel->chem.sp[cModel->chem.sCO2].mW;
  mH2Op = cModel->chem.reac[0].stch[cModel->chem.sH2O][1]
         *cModel->chem.sp[cModel->chem.sH2O].mW;
  mN2p  = cModel->chem.reac[0].stch[posN2][1]*cModel->chem.sp[posN2].mW;

  mAir  = mO2 + mN2;
  mProd = mCO2p + mH2Op + mN2p;
/*... Fuel*/
  MAT2D(cModel->chem.sCH4,0,cModel->lumpedMatrix,nl) = 1.e0; 
  MAT2D(cModel->chem.sCH4,1,cModel->lumpedMatrix,nl) = 0.e0; 
  MAT2D(cModel->chem.sCH4,2,cModel->lumpedMatrix,nl) = 0.e0; 
/*...................................................................*/

/*... O2*/
  MAT2D(cModel->chem.sO2,0,cModel->lumpedMatrix,nl) = 0.e0;
  MAT2D(cModel->chem.sO2,1,cModel->lumpedMatrix,nl) = mO2/mAir; 
  MAT2D(cModel->chem.sO2,2,cModel->lumpedMatrix,nl) = 0.0e0; 
/*...................................................................*/

/*... CO2*/
  MAT2D(cModel->chem.sCO2,0,cModel->lumpedMatrix,nl) = 0.e0; 
  MAT2D(cModel->chem.sCO2,1,cModel->lumpedMatrix,nl) = 0.0e0; 
  MAT2D(cModel->chem.sCO2,2,cModel->lumpedMatrix,nl) = mCO2p/mProd;  
/*...................................................................*/

/*... H2O*/
  MAT2D(cModel->chem.sH2O,0,cModel->lumpedMatrix,nl) = 0.e0; 
  MAT2D(cModel->chem.sH2O,1,cModel->lumpedMatrix,nl) = 0.e0; 
  MAT2D(cModel->chem.sH2O,2,cModel->lumpedMatrix,nl) =  mH2Op/mProd; 
/*...................................................................*/

/*... N2*/
  MAT2D(cModel->chem.sN2,0,cModel->lumpedMatrix,nl) = 0.e0; 
  MAT2D(cModel->chem.sN2,1,cModel->lumpedMatrix,nl) = mN2/mAir; 
  MAT2D(cModel->chem.sN2,2,cModel->lumpedMatrix,nl) = mN2p/mProd;  
/*...................................................................*/


}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * initLumpedMatrix : inicializa a matriz de especies                *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * y       -> nao definido                                           *
 * a       -> matriz das especies agrupadas                          *
 * z       -> especies agrupadas                                     *
 * ns      -> numero de especies primitivas                          *
 * nl      -> numero de especies agrupadas                           *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * y       -> especies primitivas                                    *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * c[0] - CH4                                                        *
 * c[1] - O2                                                         *
 * c[2] - CO2                                                        *
 * c[3] - H2O                                                        *
 * c[4] - N2                                                         *
 *********************************************************************/
void yLumpedMatrixZ(DOUBLE *RESTRICT y, DOUBLE *RESTRICT a
                  , DOUBLE *RESTRICT z
                  , short const ns    , short const nl)
{
  short i,j;

/*...*/
  for(i=0;i<ns;i++){
    y[i] = 0.e0;
    for(j=0;j<nl;j++)
      y[i] += MAT2D(i,j,a,nl)*z[j];  
  }
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * initMolarMass: inicializa a massa molar das especies              *
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
void initMolarMass(Combustion *cModel)
{
  
  short i,
        eO  = cModel->chem.eO,
        eN  = cModel->chem.eN,
        eC  = cModel->chem.eC,
        eH  = cModel->chem.eH;

  for(i=0;i<cModel->chem.nSp;i++)
    cModel->chem.sp[i].mW =  cModel->chem.sp[i].nO*cModel->chem.mE[eO]
                          +  cModel->chem.sp[i].nN*cModel->chem.mE[eN]
                          +  cModel->chem.sp[i].nC*cModel->chem.mE[eC]
                          +  cModel->chem.sp[i].nH*cModel->chem.mE[eH];

}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 23/08/2018                                   *
 * Data de modificaco : 21/05/2019                                   *
 *-------------------------------------------------------------------*
 * stoichiometricCoeff:                                              *
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
void stoichiometricCoeff(Combustion *cModel)
{
  short i;
  for(i=0;i<cModel->chem.nReac;i++)
  {
    fprintf(fileLogExc,"%d)\n",i);
    globalReac(cModel, i);
  }
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/08/2018                                   *
 * Data de modificaco : 26/07/2019                                   *
 *-------------------------------------------------------------------*
 * initMolarMass: inicializa a massa molar das especies              *
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
void initEntalpyOfFormation(Combustion *cModel, Prop *sHeatPol)
{
  short i;

/*...*/
  for(i=0;i<cModel->chem.nSp;i++)
  {
    cModel->chem.sp[i].entalphyOfForm = 
    cModel->chem.sp[i].entalphyOfFormMolar = polNasaH(&sHeatPol->nasa[i]
                                                     , 298.15e0, true);
    cModel->chem.sp[i].entalphyOfForm /= cModel->chem.sp[i].mW;
  }
/*...................................................................*/
  
/*...*/
  if(!mpiVar.myId)
  {
    fprintf(fileLogExc, "%-20s:\n","entalphy of Form (kj/kmol)");
    for(i=0;i<cModel->chem.nSp;i++)
      fprintf(fileLogExc, "%-20s: %+.6e\n",cModel->chem.sp[i].name
                              ,cModel->chem.sp[i].entalphyOfFormMolar);
    fprintf(fileLogExc, "%-20s:\n","entalphy of Form (kj/kg)");
    for(i=0;i<cModel->chem.nSp;i++)
      fprintf(fileLogExc, "%-20s: %+.6e\n",cModel->chem.sp[i].name
                              ,cModel->chem.sp[i].entalphyOfForm);  
  }
/*..................................................................*/
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/08/2018                                   *
 * Data de modificaco : 05/05/2019                                   *
 *-------------------------------------------------------------------*
 * initMolarMass: inicializa a massa molar das especies              *
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
void initEntalpyOfCombustion(Combustion *cModel)
{
  
  DOUBLE pCO2,pH2O,hFuel,hCO2,hH2O;

  pCO2 = cModel->CO2InProd;
  pH2O = cModel->H2OInProd;

/*...  Fuel - KJ/kMol*/
  hFuel = -74870.e0; 
/*... CO2 - KJ/kMol*/
  hCO2 = -393520.e0;
/*... H2O - KJ/kMol*/
  hH2O = -241830.e0;

/*... KJ/Kmol de fuel*/
  cModel->entalphyOfCombustion = (hCO2 + 2.0*hH2O) - hFuel;

  fprintf(fileLogExc,"Primitives species:\n");

  fprintf(fileLogExc,"Entalphy of combustion (KJ/Kmol) = %lf\n"  
                    ,cModel->entalphyOfCombustion);


/*... KJ/KG de fuel */
  cModel->entalphyOfCombustion  /= cModel->chem.sp[cModel->chem.sCH4].mW;

  fprintf(fileLogExc,"Entalphy of combustion (KJ/KG)   = %lf\n\n"  
                    ,cModel->entalphyOfCombustion);
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 21/07/2019                                   *
 *-------------------------------------------------------------------*
 * concetracionOfSpecies:                                            *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> nao definido                                           *
 * y       -> especies primitivas   (cModel->fLump = true)           *
 * c       -> concetracao das especies primitivas                    *
 * density -> densidade da mistura                                   *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * c       -> concentracao das especies                              *
 *-------------------------------------------------------------------*
 * OBS: mol/cm^3                                                             *
 *-------------------------------------------------------------------*
 *********************************************************************/
void concetracionOfSpecies(Combustion *cModel,DOUBLE *RESTRICT y
                          ,DOUBLE *RESTRICT c,DOUBLE const density)
{
  short i, ns = cModel->nOfSpecies;
  DOUBLE tmp;
  for(i=0;i<ns;i++)
  {
    tmp  = max(y[i],0.e0);
    c[i] =1.e-03*tmp*density/cModel->chem.sp[i].mW;
  }
 

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 20/07/2019                                   *
 *-------------------------------------------------------------------*
 * mixtureMolarMass: massa molar da mistura                          *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> nao definido                                           *
 * y       -> especies primitivas                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE mixtureMolarMass(Combustion *cModel,DOUBLE *RESTRICT y) 
{                    

  short i, ns = cModel->nOfSpecies;
  DOUBLE tmp;

  for(i=0,tmp=0.e0;i<ns;i++)
    tmp += y[i]/cModel->chem.sp[i].mW;

  return 1.e0/tmp;

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/08/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * mixtureMolarMassMed: massa molar media da mistura                 *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> nao definido                                           *
 * y       -> especies primitivas                                    *
 * volume  -> volume da celula                                       *
 * numel   -> numero de elementos                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE mixtureMolarMassMed(Combustion *cModel     ,DOUBLE *RESTRICT y
                          ,DOUBLE *RESTRICT volume,DOUBLE const numel) 
{                    

  short nSp = cModel->nOfSpecies;
  INT nel;
  DOUBLE vm,mm;

  for(nel=0,mm=0.e0,vm=0.e0;nel<numel;nel++)
  {
    vm += volume[nel];
    mm += mixtureMolarMass(cModel,&MAT2D(nel,0,y,nSp))*volume[nel];
  }

  return mm/vm;

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/08/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * mixtureMolarMassMedMpi:  massa molar media da mistura             *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> nao definido                                           *
 * y       -> especies primitivas                                    *
 * volume  -> volume da celula
 * numel   -> numero de elementos                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE mixtureMolarMassMedMpi(Combustion *cModel    ,DOUBLE *RESTRICT y
                            ,DOUBLE *RESTRICT volume,DOUBLE const numel) 
{                    

  short nSp = cModel->nOfSpecies;
  INT nel;
  DOUBLE vm,mm,mW;
#ifdef _MPI_
  DOUBLE gVm,gMm;
#endif

/*...*/
  for(nel=0,mm=0.e0,vm=0.e0;nel<numel;nel++)
  {
    vm += volume[nel];
    mm += mixtureMolarMass(cModel,&MAT2D(nel,0,y,nSp))*volume[nel]; 
  }
/*...................................................................*/

/*....*/
#ifdef _MPI_
  if(mpiVar.nPrcs>1)
  { 
    tm.overHeadMiscMpi = getTimeC() - tm.overHeadMiscMpi;
    MPI_Allreduce(&vm ,&gVm ,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    MPI_Allreduce(&mm ,&gMm ,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.overHeadMiscMpi = getTimeC() - tm.overHeadMiscMpi;
    mW = gMm/gVm; 
  }
  else  
    mW = mm/vm;
#else
    mW = mm/vm;
#endif
/*...................................................................*/

  return mW;

}
/*********************************************************************/



void plugFlowReactor(DOUBLE const t    ,DOUBLE *RESTRICT y
                    ,DOUBLE *RESTRICT w,void **pt)
{
  Combustion *cModel = pt[0];
  Prop *den    = pt[1],
       *sHeat  = pt[2];
  DOUBLE *T    = pt[3];
  DOUBLE *pres = pt[4];
  DOUBLE rho,cp,molarMassGas,c[MAXSPECIES],Q[MAXREAC];
  short i,
        nSp =cModel->nOfSpecies;

/*... massa especifoca*/
  molarMassGas = mixtureMolarMass(cModel,y);
  rho          = mixtureSpeciesDensity(den  ,molarMassGas
                                     ,*T    ,*pres
                                     ,*pres  ,true);
/*...............................................................*/

/*... calor especifico*/
  cp = mixtureSpecifiHeat(sHeat,y,*T,nSp,true);
/*...............................................................*/

/*...*/
  concetracionOfSpecies(cModel,y,c,rho);
/*...............................................................*/

/*...*/
  for(i=0;i<cModel->chem.nReac;i++)
    Q[i] = massActionMass(&cModel->chem.reac[i],sHeat 
                  ,c
                  ,*T                 ,nSp);  
  massRateReaction(&cModel->chem  ,Q,w);
/*...............................................................*/

/*...*/
  for(i=0;i<nSp;i++)
    w[i] /= rho;
/*...............................................................*/
}