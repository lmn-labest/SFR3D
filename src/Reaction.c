#include<Reaction.h>


void reactor( Combustion *cModel      , Prop *sHeatPol
            , DOUBLE *RESTRICT yFrac
            , DOUBLE *RESTRICT dy     , DOUBLE const density
            , DOUBLE const tMix       , DOUBLE const temp 
            , bool const fKelvin)
{
 
}

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
 * Data de modificaco : 26/07/2019                                   *
 *-------------------------------------------------------------------*
 * rateReaction: calculo da taxa de consumo do combustivel           *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> modelo de combustao                                    *
 * zComb   -> fracao massica                                         *
 * temp    -> temperatura                                            *
 * rate    -> nao definido                                           * 
 * density -> densidade                                              *
 * numel   -> numero de elementos                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * rate    -> taxa de consumo massico das especies                   * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * fLump - true                                                      *
 * z(nel,0) -> fuel(comburente)                                      *
 * z(nel,1) -> air (oxidante)                                        *
 * z(nel,2) -> produto                                               *
 * z(nel,0) -> fuel(comburente)                                      *
 * fLump - false                                                     *
 * z(nel,1) -> O2  (oxidante)                                        * 
 * z(nel,2) -> CO2                                                   *
 * z(nel,3) -> H2O                                                   *
 * z(nel,4) -> N2                                                    *
 *********************************************************************/
void rateReaction(Combustion *cModel      , Turbulence *tModel
             , Prop *sHeatPol
             , DOUBLE *RESTRICT zComb        , DOUBLE *RESTRICT diffComb
             , DOUBLE *RESTRICT temp         , DOUBLE *RESTRICT rate
             , DOUBLE *RESTRICT density      , DOUBLE *RESTRICT gradVel
             , DOUBLE *RESTRICT eddyViscosity, DOUBLE *RESTRICT dViscosity
             , DOUBLE *RESTRICT volume
             , short const ndm               , INT const numel
             , bool const fKelvin )
{
  short nComb = cModel->nComb
       , iCod = cModel->reactionKinetic
       , nReac= cModel->chem.nReac
       , nSp  = cModel->nOfSpecies;
  short i,j;
  INT nel;
  DOUBLE s,tMix,eddy,sT[6],*iGradVel,df,*pz;
  DOUBLE omega, densityC, alpha, coefA;
  DOUBLE modS, c[3], e1, e2;
  DOUBLE tK, ru, y[MAXSPECIES],cM[MAXSPECIES],Q[MAXREAC],w[MAXSPECIES];


/*...*/
  switch(iCod)
  {
/*...*/
    case ARRHENIUS:
      for(nel = 0; nel < numel; nel++)
      {
        for(i=0;i<cModel->chem.nReac;i++)
        {
/*...*/
          pz = &MAT2D(nel,0,zComb,nComb);
          densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
          getSpeciesPrimitivesCc(cModel,y,pz);
          concetracionOfSpecies(cModel            ,y
                               ,cM                ,densityC);

          TEMP(tK,temp[nel],fKelvin);
/*... mol/(cm3 s)*/
          omega = massActionMass(&cModel->chem.reac[i],sHeatPol 
                                ,cM
                                ,tK                   ,nComb);
/*... kmol/(m3 s)*/
          Q[i] = 1.e+03*omega;
        }
/*...................................................................*/

      massRateReaction(&cModel->chem,Q,w);
      
      for(i=0;i<nSp;i++)
        MAT2D(nel,i,rate,nSp) = w[i];
 
    } 
/*...................................................................*/    
    break;
/*...................................................................*/

/*...*/
    case EDC:
//    for(i=0;i<cModel->chem.nReac;i++)
//    {
/*...*/
//      for(nel = 0; nel < numel; nel++)
//      {
//        pz = &MAT2D(nel,0,zComb,nComb);
//        densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
//        getSpeciesPrimitivesCc(cModel,y,pz);
/*.. calculo Sij*/
//        iGradVel = &MAT3D(nel,0,0,gradVel,ndm,ndm);
//        tensorS(sT,iGradVel,false);
/*... |S| = sqrt(2S:S)*/
//        modS = sqrt(2.e0*doubleDotSym(sT));
/*...................................................................*/
//        df    = MAT2D(nel,cModel->sp_fuel[i],diffComb,nComb);
/*...*/
/*        omega  = edc(y              ,iComb
                      ,iOx            ,iProd 
                      ,cModel->nSpeciesPart[i][1]
                      ,s              ,densityC
                      ,volume[nel]    ,eddyViscosity[nel]
                      ,c              ,modS
                      ,dViscosity[nel],df
                      ,tMix           ,tChemical
                      ,cModel->edc.type);*/
//        omega = edc2(cModel         ,sHeatPol 
//                    ,y              ,iComb
//                    ,iOx            ,iProd 
//                    ,cModel->nSpeciesPart[i][1]
//                    ,s               ,densityC
//                    ,volume[nel]     ,eddyViscosity[nel]
//                    ,c               ,modS
//                    ,dViscosity[nel] ,df
//                    ,tMix            ,tChemical
//                    ,temp[nel]
//                    ,cModel->edc.type,fKelvin);
/*...................................................................*/
//        MAT2D(nel,i,rate,nReac) = omega; 
//      }
//    }
/*...................................................................*/
      break;
/*...................................................................*/

/*...*/
    case EDM:
/*...*/
      for(nel = 0; nel < numel; nel++)
      {
        pz = &MAT2D(nel,0,zComb,nComb);
        densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
        getSpeciesPrimitivesCc(cModel,y,pz);
/*.. calculo Sij*/
        iGradVel = &MAT3D(nel,0,0,gradVel,ndm,ndm);
        tensorS(sT,iGradVel,false);
/*... |S| = sqrt(2S:S)*/
        modS = sqrt(2.e0*doubleDotSym(sT));
/*...................................................................*/
        edm(cModel,y,w,densityC,modS);
/*...................................................................*/
        for(i=0;i<nSp;i++)
          MAT2D(nel,i,rate,nSp) = w[i];
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
 * Data de modificaco : 24/05/2019                                   *
 *-------------------------------------------------------------------*
 * rateReaction: calculo da taxa de consumo do combustivel        *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> modelo de combustao                                    *
 * zFrac   -> fracao massica                                         *
 * temp    -> temperatura                                            *
 * density -> densidade                                              *
 * eddyVis -> viscosidade turbulenta                                 * 
 * dVisc   -> viscosidade dinamica                                   *
 * tReacot -> escala de tempo do reator                              * 
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
             , DOUBLE *RESTRICT zComb     , DOUBLE *RESTRICT temp  
             , DOUBLE *RESTRICT density   , DOUBLE *RESTRICT gradVel 
             , DOUBLE *RESTRICT eddyViscosity
             , DOUBLE *RESTRICT dViscosity, DOUBLE *RESTRICT tReactor
             , short const ndm            , INT const numel   
             , bool const fKelvin )
{
  short nc = cModel->nComb       
       , nReac=cModel->chem.nReac;
  short iComb,iOx,iProd[3],i,j;
  INT nel;
  DOUBLE s,tMix,eddy,sT[6],*iGradVel,*pz;
  DOUBLE omega, densityC, alpha, coefA,tmp;
  DOUBLE modS, e1, e2;
  DOUBLE mWfuel,mWox, tempA, ru, y[MAXSPECIES],tChemical; 

  tChemical = 0.e0;

  ru = IDEALGASRC;
/*...*/
  for(nel = 0; nel < numel; nel++)
  {
    pz = &MAT2D(nel,0,zComb,nc);
    densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
    getSpeciesPrimitivesCc(cModel,y,pz);

/*.. calculo Sij*/
    iGradVel = &MAT3D(nel,0,0,gradVel,ndm,ndm);
    tensorS(sT,iGradVel,false);
/*... |S| = sqrt(2S:S)*/
    modS = sqrt(2.e0*doubleDotSym(sT));
/*...................................................................*/

/*...*/
//  for(i=0;i<cModel->nReac;i++)
//  {
//    iComb  = cModel->sp_fuel[i];
//    alpha  = cModel->arrhenius[i].alpha;
//    mWfuel = cModel->mW[iComb];

//    tempA  = cModel->arrhenius[i].energyAtivation/ru;
//    coefA  = cModel->arrhenius[i].a;
//    e1     = cModel->arrhenius[i].e1; 
//    e2     = cModel->arrhenius[i].e2; 
//    s      = cModel->stoich[i][0][cModel->sp_O2];
//    mWox   = cModel->mW[cModel->sp_O2];    
/*...................................................................*/

/*...*/
//    omega = arrhenius(y[iComb]      ,y[cModel->sp_O2]
//                     ,e1       ,e2   
//                     ,mWfuel   ,mWox
//                     ,temp[nel] ,alpha
//                     ,densityC  ,tempA 
//                     ,coefA     ,fKelvin);
/*...................................................................*/
//    tmp       = fabs(omega/densityC);
//    tmp       = y[iComb]/tmp;
//    tChemical = max(tmp,tChemical);
//  }
/*...................................................................*/ 

/*...*/    
//  MAT2D(nel,0,tReactor,2) = 1.e0/modS;
//  MAT2D(nel,1,tReactor,2) = tChemical; 
//  if( nel == 5599)
//    fprintf(fileLogDebug,"%e %e %e\n", MAT2D(nel,0,tReactor,2)
//                                     , MAT2D(nel,1,tReactor,2)
//                                     , temp[nel]);
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

/*...  CH4*/
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
      fprintf(fileLogExc, "%-20s: %lf\n",cModel->chem.sp[i].name
                              ,cModel->chem.sp[i].entalphyOfFormMolar);
    fprintf(fileLogExc, "%-20s:\n","entalphy of Form (kj/kg)");
    for(i=0;i<cModel->chem.nSp;i++)
      fprintf(fileLogExc, "%-20s: %lf\n",cModel->chem.sp[i].name
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