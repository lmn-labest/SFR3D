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
  #include<Edo.h>
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
  void rateReaction(Combustion *cModel         , Turbulence *tModel
             , PropVarFluid *pFluid
             , DOUBLE *RESTRICT zComb        , DOUBLE *RESTRICT temp        
             , DOUBLE *RESTRICT rate         , DOUBLE *RESTRICT density 
             , DOUBLE *RESTRICT gradVel      , DOUBLE *RESTRICT eddyViscosity
             , DOUBLE *RESTRICT dViscosity   , DOUBLE *RESTRICT volume
             , DOUBLE const dt               , DOUBLE const Pth 
             , short const ndm               , INT const numel
             , bool const fKelvin            , bool const fOmp           
             , short const nThreads);
/*...................................................................*/

/*...*/
  void timeChemical(Combustion *cModel      , Turbulence *tModel
             , PropVarFluid *pFluid
             , DOUBLE *RESTRICT zComb        , DOUBLE *RESTRICT temp  
             , DOUBLE *RESTRICT density      , DOUBLE *RESTRICT gradVel 
             , DOUBLE *RESTRICT eddyViscosity, DOUBLE *RESTRICT sHeat
             , DOUBLE *RESTRICT tCond        , DOUBLE *RESTRICT volume  
             , DOUBLE *RESTRICT dViscosity   , DOUBLE *RESTRICT tReactor
             , short const ndm               , INT const numel   
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
  DOUBLE mixtureMolarMass(Combustion *cModel,DOUBLE *RESTRICT y);
  DOUBLE mixtureMolarMassMed(Combustion *cModel  ,DOUBLE *RESTRICT y
                            ,DOUBLE *RESTRICT volume,DOUBLE const numel);
  DOUBLE mixtureMolarMassMedMpi(Combustion *cModel  ,DOUBLE *RESTRICT y
                            ,DOUBLE *RESTRICT volume,DOUBLE const numel);
/*...................................................................*/

/*... reatores*/
  void plugFlowReactor(DOUBLE const t    ,DOUBLE *RESTRICT y
                    ,DOUBLE *RESTRICT w,void **pt);
/*...................................................................*/

#endif/*_REACTION_H_*/