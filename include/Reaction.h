#ifndef _REACTION_H_
  #define _REACTION_H_
/*...*/
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<HccaStdBool.h>
  #include<Define.h>
  #include<StructDefine.h>
  #include<CellLoop.h>
/*...................................................................*/

/*...*/
  DOUBLE arrhenius(DOUBLE const y1   ,DOUBLE const y2
                ,DOUBLE const e1     ,DOUBLE const e2
                ,DOUBLE const mW1    ,DOUBLE const mW2
                ,DOUBLE const t      ,DOUBLE const alpha
                ,DOUBLE const density,DOUBLE const tA    
                ,DOUBLE const coefA  ,bool const fKelvin);
/*...................................................................*/

/*...*/
  DOUBLE massActionMass(Reaction *reac     ,Prop *sHeatPol 
                     ,DOUBLE  *RESTRICT c
                     ,DOUBLE const T     ,unsigned short const nSp);
/*...................................................................*/

/*...*/
  void massRateReaction(Chemical *chem  ,DOUBLE *RESTRICT Q
                     ,DOUBLE *RESTRICT w);
/*...................................................................*/

/*...*/
  void rateReaction(Combustion *cModel    , Turbulence *tModel
             , Prop *sHeatPol 
             , DOUBLE *RESTRICT zComb        , DOUBLE *RESTRICT diffComb
             , DOUBLE *RESTRICT temp         , DOUBLE *RESTRICT rate
             , DOUBLE *RESTRICT density      , DOUBLE *RESTRICT gradVel
             , DOUBLE *RESTRICT eddyViscosity, DOUBLE *RESTRICT dViscosity
             , DOUBLE *RESTRICT volume
             , short const ndm               , INT const numel
             , bool const fKelvin  );
/*...................................................................*/

/*...*/
 void timeChemical(Combustion *cModel      , Turbulence *tModel
             , DOUBLE *RESTRICT zFrac     , DOUBLE *RESTRICT temp  
             , DOUBLE *RESTRICT density    
             , DOUBLE *RESTRICT gradVel   , DOUBLE *RESTRICT eddyViscosity
             , DOUBLE *RESTRICT dViscosity, DOUBLE *RESTRICT tReactor
             , short const ndm            , INT const numel   
             , bool const fKelvin );
/*...................................................................*/

/*...*/
  void initLumpedMatrix(Combustion *cModel);
  void yLumpedMatrixZ(DOUBLE *RESTRICT y, DOUBLE *RESTRICT a
                  , DOUBLE *RESTRICT z
                  , short const ni    , short const nj);
  void initMolarMass(Combustion *cModel);
  void initEntalpyOfFormation(Combustion *cModel, Prop *sHeatPol);
  void stoichiometricCoeff(Combustion *cModel);
  void globalReac(Combustion *c, short const iReac);
  void initEntalpyOfCombustion(Combustion *cModel);

  void concetracionOfSpecies(Combustion *cModel,DOUBLE *RESTRICT z
                          ,DOUBLE *RESTRICT c,DOUBLE const density);
/*...................................................................*/

/*...*/
  DOUBLE mixtureMolarMass(Combustion *cModel,DOUBLE *RESTRICT z);
/*...................................................................*/

#endif/*_REACTION_H_*/