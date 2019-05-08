#ifndef _PROPERTIES_H_
  #define  _PROPERTIES_H_
/*...*/
  #include<math.h>
  #include<stdio.h>
  #include<stdlib.h>
  #include<string.h>
/*...................................................................*/

/*...*/
  #include<Combustion.h>
  #include<Define.h>
  #include<HccaStdBool.h>
  #include<Erro.h>
  #include<File.h>
  #include<Mesh.h>
/*...................................................................*/

/*...agua*/
  DOUBLE waterDensity(DOUBLE const t);
  DOUBLE waterSpecifiHeat(DOUBLE const t);
  DOUBLE waterDynamicViscosity(DOUBLE const t);
/*...................................................................*/

/*... gas ideal incompressivel(Ar)*/
  DOUBLE airDensity(PropPol *den
                   ,DOUBLE const t, DOUBLE const presRef
                   ,DOUBLE const p, bool const fKelvin);
  DOUBLE airSpecifiHeat(PropPol *sHeatPol
                       ,DOUBLE const t,bool const fKelvin);
  DOUBLE airDynamicViscosity(PropPol *dVisc,DOUBLE const t
                            ,bool const fKelvin);
  DOUBLE airThermalConductvity(PropPol *thC,DOUBLE const t
                              ,bool const fKelvin);
  DOUBLE specificEnthalpyForTemp(PropPol *sHeatPol
                               , DOUBLE const hs, DOUBLE const sHeatRef
                               , bool const fSheat, bool const fKelvin); 
/*...*/
  DOUBLE tempForSpecificEnthalpy(PropPol *sHeatPol
                               , DOUBLE const t, DOUBLE const sHeatRef
                               , bool const fSheat, bool const fKelvin);
/*...................................................................*/

/*...*/
  void updateDensity(PropPol *pDen
                    ,DOUBLE *RESTRICT temp   , DOUBLE *RESTRICT pressure
                    ,DOUBLE *RESTRICT density                 
                    ,DOUBLE const alpha        ,bool const iKelvin 
                    ,INT const nEl             ,char  const iCod);

  void updateSpecificHeat(PropPol *sHeatPol
                         ,DOUBLE *RESTRICT temp,DOUBLE *RESTRICT sHeat
                         ,bool const iKelvin  
                         ,INT const nEl        ,char  const iCod);
  void updateDynamicViscosity(PropPol *dVisc
                            ,DOUBLE *RESTRICT temp,DOUBLE *RESTRICT visc    
                            ,bool const iKelvin   
                            ,INT const nEl);
  void updateThermalconductivity(PropPol *thC
                              ,DOUBLE *RESTRICT t,DOUBLE *RESTRICT thCond    
                              ,bool const iKelvin       
                              ,INT const nEl);
  
  void updateDensityCD(PropPol *pol, DOUBLE *RESTRICT u
                     , DOUBLE *RESTRICT density, INT nEl
                     , char  iCod);

  void updateProp(PropPol *pol, DOUBLE *RESTRICT u
                , DOUBLE *RESTRICT coef, INT nEl);


  void initPropTemp(PropVarFluid *propFluid
                  ,DOUBLE *RESTRICT prop    ,DOUBLE *RESTRICT t 
                  ,DOUBLE *RESTRICT pressure,DOUBLE *RESTRICT propMat
                  ,short *RESTRICT mat
                  ,short const np           ,INT const nCell 
                  ,bool const iKelvin       ,short const iProp);

  void initPropCD(PropPol *pol, DOUBLE *RESTRICT prop
                 , DOUBLE *RESTRICT u, DOUBLE *RESTRICT propMat
                 , short *RESTRICT mat
                 , short np, INT    nCell
                 , short iProp);

  void initPropStructCD(PropVarCD *propVar, short const n);
/*...................................................................*/

/*...*/
  void getTempForEnergy(PropPol *sHeatPol
                     ,DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                     ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                     ,INT const nCell      ,bool const fTemp
                     ,bool const fSheat    ,bool const fKelvin
                     ,bool const fOmp      ,short const nThreads );

  void getEnergyForTemp(PropPol *sHeatPol
                     ,DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                     ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                     ,INT const nCell     
                     ,bool const fSheat    ,bool const fKelvin
                     ,bool const fOmp      ,short const nThreads );
/*...................................................................*/

/*...*/
  DOUBLE diffProp(PropPol *pol, DOUBLE u);
/*...................................................................*/

/*...*/
  void specificMassRef(DOUBLE *RESTRICT density, DOUBLE *RESTRICT volume                  
                  , DOUBLE *RESTRICT prop      , short  *RESTRICT mat
                  , INT const nCell);
  void presRef(DOUBLE *RESTRICT temp0 , DOUBLE *RESTRICT temp  
             , DOUBLE *RESTRICT volume  , DOUBLE *pTh                  
             , INT const nCell          , bool const fKelvin);
  void initPresRef(DOUBLE *RESTRICT temp  , DOUBLE *RESTRICT volume
               , DOUBLE *pTh            , DOUBLE *RESTRICT prop  
               , short  *RESTRICT mat   , DOUBLE const molarMass                  
               , INT const nCell        , bool const fKelvin);
/*...................................................................*/

/*...*/
  void initDensityPol(PropPol *prop, char *s, FILE *file);
  void initSheatPol(PropPol *prop, char *s, FILE *file);
  void initDviscosityPol(PropPol *prop, char *s, FILE *file);
  void initThCondPol(PropPol *prop, char *s, FILE *file);
/*...................................................................*/

/*... mistura gasosa*/
  void initMixtureSpeciesfiHeat(PropPol *prop, char *s, FILE *file);
  void updateMixSpecificHeat(PropPol *sHeatPol
                         , DOUBLE *RESTRICT temp  , DOUBLE *RESTRICT yFrac  
                         , DOUBLE *RESTRICT sHeat , short const nOfPrSp
                         , bool const iKelvin
                         , INT const nEl          , char  const iCod);
  void updateMixDensity(PropPol *pDen         , Combustion *cModel
                 , DOUBLE *RESTRICT temp    , DOUBLE *RESTRICT pressure
                 , DOUBLE *RESTRICT density , DOUBLE *RESTRICT zComb
                 , DOUBLE const alpha       , bool const iKelvin    
                 , INT const nEl            , char  const iCod);
  DOUBLE mixtureSpecifiHeat(PropPol *sHeat   , DOUBLE *yFrac
                         , DOUBLE const t    , short const nOfPrSp
                         , bool const fKelvin); 
  DOUBLE specieSpecifiHeat(PropPol *sHeat     , short const kSpecie
                        , DOUBLE const t      , bool const fKelvin); 
  DOUBLE mixtureSpeciesDensity(PropPol *den        ,DOUBLE const malorMassMix
                            ,DOUBLE const t      ,DOUBLE const p
                            ,DOUBLE const presRef,bool const fKelvin);
  DOUBLE specificEnthalpyForTempOfMix(PropPol *sHeatPol
                             , DOUBLE const hs        , DOUBLE *yFrac
                             , DOUBLE const sHeatRef  , short const nOfPrSp
                             , bool const fSheat      , bool const fKelvin
                             , INT const nel ); 

  DOUBLE tempForSpecificEnthalpyMix(PropPol *sHeat    , DOUBLE *yFrac
                                , DOUBLE const t    , DOUBLE const sHeatRef
                                , short const nOfPrSp
                                , bool const fSheat , bool const fKelvin); 
  DOUBLE tempForSpecificEnthalpySpecies(PropPol *sHeat, short const kSpecie
                               , DOUBLE const t   
                               , bool const fSheat , bool const fKelvin);

  void getEnergyForTempMix(PropPol *sHeatPol  ,DOUBLE *RESTRICT yFrac 
                        ,DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                        ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                        ,INT const nCell      ,short const nOfPrSp
                        ,bool const fSheat    ,bool const fKelvin
                        ,bool const fOmp      ,short const nThreads );


  void getTempForEnergyMix(PropPol *sHeatPol    ,DOUBLE *RESTRICT yFrac
                        ,DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                        ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                        ,INT const nCell      ,short const nOfPrSp 
                        ,bool const fTemp     ,bool const fSheat    
                        ,bool const fKelvin
                        ,bool const fOmp      ,short const nThreads );

  void initPropTempMix(PropVarFluid *propFluid, Combustion *cModel
                 ,DOUBLE *RESTRICT prop 
                 ,DOUBLE *RESTRICT t        ,DOUBLE *RESTRICT pressure
                 ,DOUBLE *RESTRICT yFrac    ,DOUBLE *RESTRICT propMat
                 ,short *RESTRICT mat       ,short const nOfPrSp 
                 ,short const np            ,INT    const nCell 
                 ,bool const iKelvin        ,short const iProp);
/*...................................................................*/

/*...*/
  void initCdPol(PropPol *prop, char *s, FILE *file);
/*...................................................................*/

/*...*/
//DOUBLE mixtureMolarMass(Combustion *cModel,DOUBLE *RESTRICT z);
/*...................................................................*/

  int readFileLineSimple(DOUBLE *x, FILE *file);

#endif /*_PROPERTIES_H_*/