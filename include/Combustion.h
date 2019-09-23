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
  #include<Reaction.h>
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
  void rateHeatRealesedReaction(Combustion *cModel,Prop *sHeat   
                   , DOUBLE *RESTRICT q      , DOUBLE *RESTRICT temp
                   , DOUBLE *RESTRICT zComb0 , DOUBLE *RESTRICT zComb
                   , DOUBLE *RESTRICT density, DOUBLE *RESTRICT Q 
                   , DOUBLE *RESTRICT prop   , short  *RESTRICT mat
                   , DOUBLE const dt         , INT const numel
                   , bool const fsHeat       , bool const fKelvin
                   , bool const fOmp         , short const nThreads);
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
                      , INT numel                 , bool fKelvin
                      , bool fOmp                 , short nThreads );
  void  getGradSpecies(Combustion *cModel   
                   , DOUBLE *RESTRICT gradZ, DOUBLE *RESTRICT gradY
                   , INT const numel       , short const ndm);
/*...*/
  DOUBLE maxArray(DOUBLE *RESTRICT x,INT const n);
  DOUBLE getVolumeMed(DOUBLE *RESTRICT x,DOUBLE *RESTRICT vol
                   ,INT const n);
  void regularZ(DOUBLE *RESTRICT z,INT const numel, short const nComb);


  void sumFracZ(DOUBLE *z      ,DOUBLE *zComb 
              ,INT const n     ,short const nComb);

  DOUBLE totalHeatRealeseComb(DOUBLE *RESTRICT q, DOUBLE *RESTRICT vol  
                          , DOUBLE const dt     , INT const numel);

/*...*/
  INT edc(Combustion *c           ,PropVarFluid *pFluid
        ,DOUBLE *RESTRICT y       ,DOUBLE *RESTRICT w
        ,DOUBLE *RESTRICT tReactor,DOUBLE const density   
        ,DOUBLE const dt          ,DOUBLE const temp 
        ,DOUBLE const eddyVisc    ,DOUBLE const dVisc
        ,DOUBLE const Pth         ,bool const fKelvin
        ,INT const nel);
  
  void edm(Combustion *c       ,DOUBLE *RESTRICT y  
          ,DOUBLE *RESTRICT w  ,DOUBLE const density
          ,DOUBLE const itMix);
/*...................................................................*/

#endif/*_COMBUSTION_H_*/