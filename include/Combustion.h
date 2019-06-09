#ifndef _COMBUSTION_H_
  #define _COMBUSTION_H_
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<Memoria.h>
  #include<Mesh.h>
  #include<HccaStdBool.h>
  #include<CellLoop.h>
  #include<Define.h>
//#include<Properties.h>
  #include<Sisteq.h> 
  #include<Solv.h>
/*...................................................................*/

/*...*/
void combustionSolver(Memoria *m          , PropVarFluid *propF
                    , Loads *loadsVel     , Loads *loadsPres
                    , Loads *loadsEnergy  , Loads *loadsKturb
                    , Loads *loadsComb
                    , EnergyModel *eModel , Combustion *cModel
                    , MassEqModel *eMass  , MomentumModel *ModelMomentum
                    , Turbulence *tModel  , ThermoDynamic *thDynamic
                    , Mesh *mesh0         , Mesh *mesh
                    , SistEq *sistEqVel   , SistEq *sistEqPres
                    , SistEq *sistEqEnergy, SistEq *sistEqKturb
                    , SistEq *sistEqComb
                    , Solv *solvVel       , Solv *solvPres
                    , Solv *solvEnergy    , Solv *solvKturb
                    , Solv *solvComb
                    , Simple *sp          , Scheme *sc
                    , PartMesh *pMesh     , Mean *media
                    , FileOpt *opt        , char *preName
                    , char *nameOut       , FILE *fileOut);
/*...................................................................*/

/*...*/
  void combustionModel(Memoria *m       , PropVarFluid *prop
                   , Loads *loadsComb   , Loads *loadsVel
                   , Turbulence *tModel , Combustion *cModel
                   , EnergyModel *eModel, Mesh *mesh
                   , SistEq *sistEqComb , Solv *solvComb  
                   , Simple *sp         , Scheme *sc        
                   , PartMesh *pMesh    , FileOpt *opt
                   , bool *fComb        , short itSimple    ); 
/*...................................................................*/

/*...*/
  void rateFuelConsume(Combustion *cModel    , Turbulence *tModel
             , DOUBLE *RESTRICT zComb        , DOUBLE *RESTRICT diffComb
             , DOUBLE *RESTRICT temp         , DOUBLE *RESTRICT rate
             , DOUBLE *RESTRICT density      , DOUBLE *RESTRICT gradVel
             , DOUBLE *RESTRICT eddyViscosity, DOUBLE *RESTRICT dViscosity
             , DOUBLE *RESTRICT volume
             , short const ndm               , INT const numel
             , bool const fKelvin  );
/*...................................................................*/

/*...*/
  void rateHeatRealeseCombustion(Combustion *cModel,Prop *sHeat   
                   , DOUBLE *RESTRICT q      , DOUBLE *RESTRICT temp
                   , DOUBLE *RESTRICT zComb0 , DOUBLE *RESTRICT zComb
                   , DOUBLE *RESTRICT density, DOUBLE *RESTRICT rateFuel 
                   , DOUBLE *RESTRICT prop   , short  *RESTRICT mat
                   , DOUBLE const dt         , INT const numel
                   , bool const fsHeat       , bool const fKelvin);
/*...................................................................*/

/*...*/
  void getSpeciesPrimitives(Combustion *cModel 
                        , DOUBLE *RESTRICT y,DOUBLE *RESTRICT z
                        , INT const numel);
  void getSpeciesPrimitivesCc(Combustion *cModel 
                           , DOUBLE *RESTRICT y,DOUBLE *RESTRICT z);
/*...................................................................*/

/*...*/
  void getEnthalpySpecies(Combustion *cModel      ,  PropVarFluid *propF
                      , DOUBLE *RESTRICT enthalpyk, DOUBLE *RESTRICT temp 
                      , INT const numel           , bool const fKelvin);
  void  getGradSpecies(Combustion *cModel   
                   , DOUBLE *RESTRICT gradZ, DOUBLE *RESTRICT gradY
                   , INT const numel       , short const ndm);
/*...*/
  DOUBLE maxArray(DOUBLE *RESTRICT x,INT const n);
  DOUBLE getVolumeMed(DOUBLE *RESTRICT x,DOUBLE *RESTRICT vol
                   ,INT const n);
  void regularZ(DOUBLE *RESTRICT z,INT const numel, short const nComb);
  void initLumpedMatrix(Combustion *cModel);
  void yLumpedMatrixZ(DOUBLE *RESTRICT y, DOUBLE *RESTRICT a
                  , DOUBLE *RESTRICT z
                  , short const ni    , short const nj);
  void initMolarMass(Combustion *cModel);
  void initEntalpyOfFormation(Combustion *cModel);
  void stoichiometricCoeff(Combustion *cModel);
  void globalReac(Combustion *c, short const iReac);
  void initEntalpyOfCombustion(Combustion *cModel);

  void sumFracZ(DOUBLE *z      ,DOUBLE *zComb 
              ,INT const n     ,short const nComb);

  DOUBLE totalHeatRealeseComb(DOUBLE *RESTRICT q, DOUBLE *RESTRICT vol  
                          , DOUBLE const dt     , INT const numel);

/*...*/
  DOUBLE mixtureMolarMass(Combustion *cModel,DOUBLE *RESTRICT z);
/*...................................................................*/

/*...*/
  DOUBLE edc(DOUBLE *y           ,short const iYf
          ,short const iYox   ,short *iProd
          ,short const nProd
          ,DOUBLE const s     ,DOUBLE const density
          ,DOUBLE const vol   ,DOUBLE const eddyVisc
          ,DOUBLE *c          ,DOUBLE const modS 
          ,DOUBLE const dVisc ,DOUBLE const df 
          ,DOUBLE const tMix  ,short const iCod);
  DOUBLE arrhenius(DOUBLE const y1   ,DOUBLE const y2
                ,DOUBLE const a1     ,DOUBLE const a2
                ,DOUBLE const mW1    ,DOUBLE const mW2
                ,DOUBLE const t      ,DOUBLE const alpha
                ,DOUBLE const density,DOUBLE const tA    
                ,DOUBLE const coefA  ,bool const fKelvin);
/*...................................................................*/

#endif/*_COMBUSTION_H_*/