#ifndef _PROPERTIES_H_
  #define  _PROPERTIES_H_
/*...*/
  #include<math.h>
  #include<stdio.h>
  #include<stdlib.h>
/*...................................................................*/

/*...*/
  #include<Define.h>
  #include<HccaStdBool.h>
  #include<Erro.h>
/*...................................................................*/

/*...*/
  typedef struct{
    short nPol;
    DOUBLE a[10];   
  }PropPol;
/*...................................................................*/

/*...*/
  #define TREF      288.15e+00    /*Kelvin         */
  #define PREREF    1.01325e-01   /*MPa            */
  #define IDEALGASR 8.3144598e-06 /*MJ/(mol.kelvin)*/
  #define MMOLARAR  2.896e-2      /*kg/mol         */
/*...................................................................*/

/*...*/
  #define PROP_UPDATE_SIMPLE_LOOP 0
  #define PROP_UPDATE_OLD_TIME    1
/*...................................................................*/

/*...agua*/
  DOUBLE waterDensity(DOUBLE const t);
  DOUBLE waterSpecifiHeat(DOUBLE const t);
  DOUBLE waterDynamicViscosity(DOUBLE const t);
/*...................................................................*/

/*... gas ideal incompressivel(Ar)*/
  DOUBLE airDensity(DOUBLE const t,bool fKelvin);
  DOUBLE airSpecifiHeat(DOUBLE const t,bool fKelvin);
  DOUBLE airDynamicViscosity(DOUBLE const t,bool fKelvin);
  DOUBLE airThermalConductvity(DOUBLE const t,bool fKelvin);
  DOUBLE specificEnthalpyForTemp(DOUBLE const hs,bool fKelvin);
  DOUBLE tempForSpecificEnthalpy(DOUBLE const t,bool fKelvin);
/*...................................................................*/

/*...*/
  void updateDensity(DOUBLE *RESTRICT temp     ,DOUBLE *RESTRICT density
                    ,INT const nEl             ,char  const iCod);
  void updateSpecificHeat(DOUBLE *RESTRICT temp,DOUBLE *RESTRICT sHeat
                         ,INT const nEl        ,char  const iCod);
  void updateDynamicViscosity(DOUBLE *RESTRICT temp,DOUBLE *RESTRICT visc    
                            ,INT const nEl        );
  void updateThermalCondutivty(DOUBLE *RESTRICT t,DOUBLE *RESTRICT thCond    
                              ,INT const nEl);
  
 void initPropTemp(DOUBLE *RESTRICT prop,DOUBLE *RESTRICT t 
                  ,DOUBLE *RESTRICT propMat,short *RESTRICT mat
                  ,short const np          ,INT const nCell 
                  ,short const iProp);
/*...................................................................*/

/*...*/
  void initSheatPol(void);
/*...................................................................*/

/*...*/
  bool iKelvin;
  PropPol sHeat;
/*...................................................................*/
#endif /*_PROPERTIES_H_*/