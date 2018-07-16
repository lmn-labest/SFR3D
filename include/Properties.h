#ifndef _PROPERTIES_H_
  #define  _PROPERTIES_H_
/*...*/
  #include<math.h>
  #include<stdio.h>
  #include<stdlib.h>
  #include<string.h>
/*...................................................................*/

/*...*/
  #include<Define.h>
  #include<HccaStdBool.h>
  #include<Erro.h>
  #include<File.h>
  #include<Mesh.h>
/*...................................................................*/

/*...*/
  typedef struct{
    bool fDensityRef;
    bool fPresTh;
    DOUBLE pTh[3];
  }ThermoDynamic;
   ThermoDynamic thDynamic;
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
  void initPresRef(DOUBLE *RESTRICT temp  
                 , DOUBLE *RESTRICT volume, DOUBLE *pTh   
                 , DOUBLE *RESTRICT prop  , short  *RESTRICT mat                     
                 , INT const nCell        , bool const fKelvin);
/*...................................................................*/

/*...*/
  void initDensityPol(PropPol *prop, char *s, FILE *file);
  void initSheatPol(PropPol *prop, char *s, FILE *file);
  void initDviscosityPol(PropPol *prop, char *s, FILE *file);
  void initThCondPol(PropPol *prop, char *s, FILE *file);
/*...................................................................*/

/*...*/
  void initCdPol(PropPol *prop, char *s, FILE *file);
/*...................................................................*/

  int readFileLineSimple(DOUBLE *x, FILE *file);

#endif /*_PROPERTIES_H_*/