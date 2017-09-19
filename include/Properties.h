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
/*...................................................................*/

/*...*/
  #define POL        1
  #define SUTHERLAND 2
/*...................................................................*/

/*...*/
  typedef struct{
    bool fDensityRef;
    bool fPresTh;
    DOUBLE pTh[3];
  }ThermoDynamic;
   ThermoDynamic thDynamic;
/*...................................................................*/  

/*...*/
  typedef struct{
    unsigned char type;
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
  #define TEMP_FOR_ENTHALPY(cp,t,tr) ( cp*(t-tr) )
  #define ENTHALPY_FOR_TEMP(cp,hs,tr) ( hs/cp+tr )
/*...................................................................*/

/*...*/
  #define PRESREF(dRef,R,T,Mg) (dRef*R*T/Mg) 
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
  DOUBLE airDensity(DOUBLE const t, DOUBLE const presRef
                   ,bool const fKelvin);
  DOUBLE airSpecifiHeat(DOUBLE const t,bool const fKelvin);
  DOUBLE airDynamicViscosity(DOUBLE const t,bool const fKelvin);
  DOUBLE airThermalConductvity(DOUBLE const t,bool const fKelvin);
  DOUBLE specificEnthalpyForTemp(DOUBLE const hs,bool const fKelvin);
  DOUBLE tempForSpecificEnthalpy(DOUBLE const t,bool const fKelvin);
/*...................................................................*/

/*...*/
  void updateDensity(DOUBLE *RESTRICT temp     ,DOUBLE *RESTRICT density
                    ,DOUBLE const alpha        ,bool const iKelvin 
                    ,INT const nEl             ,char  const iCod);
  void updateSpecificHeat(DOUBLE *RESTRICT temp,DOUBLE *RESTRICT sHeat
                         ,bool const iKelvin  
                         ,INT const nEl        ,char  const iCod);
  void updateDynamicViscosity(DOUBLE *RESTRICT temp,DOUBLE *RESTRICT visc    
                            ,bool const iKelvin   
                            ,INT const nEl);
  void updateThermalCondutivty(DOUBLE *RESTRICT t,DOUBLE *RESTRICT thCond    
                              ,bool const iKelvin       
                              ,INT const nEl);
  
 void initPropTemp(DOUBLE *RESTRICT prop   ,DOUBLE *RESTRICT t 
                  ,DOUBLE *RESTRICT propMat,short *RESTRICT mat
                  ,short const np          ,INT const nCell 
                  ,bool const iKelvin      ,short const iProp);
/*...................................................................*/

/*...*/
 void getTempForEnergy(DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                     ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                     ,INT const nCell      ,bool const fTemp
                     ,bool const fSheat    ,bool const fKelvin);

  void getEnergyForTemp(DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                       ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                       ,INT const nCell     
                       ,bool const fHeat     ,bool const fKelvin);
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
  void initSheatPol(void);
  void initDviscosityPol(char *s);
  void initThCondPol(char *s);
/*...................................................................*/

/*...*/
  PropPol sHeat,dVisc,thCond;
/*...................................................................*/
#endif /*_PROPERTIES_H_*/