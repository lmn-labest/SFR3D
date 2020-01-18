#ifndef _SIMPLE_H_
  #define _SMIPLE_H_
/*...*/ 
  #include<stdlib.h>
  #include<stdio.h>
 
/*...*/
  #include<File.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<CellLoop.h>
  #include<Combustion.h>
  #include<Mesh.h>
  #include<WriteVtk.h>
  #include<Sisteq.h>
  #include<Solv.h>
  #include<Transient.h>
  #include<Turbulence.h>
  #include<Properties.h>
  #include<Energy.h>
  #include<Media.h>
  #include<Residual.h>
/*...................................................................*/

  void simpleSolver(Memoria *m        
                   ,Loads *loadsVel   ,Loads *loadsPres 
                   ,MassEqModel *eMass, MomentumModel *momentumModel
                   ,Turbulence *tModel 
                   ,Mesh *mesh0       ,Mesh *mesh       
                   ,SistEq *sistEqVel ,SistEq *sistEqPres
                   ,Solv *solvVel     ,Solv *solvPres 
                   ,Simple *sp
                   ,Scheme *sc        ,PartMesh *pMesh 
                   ,FileOpt *opt       ,char *preName 
                   ,char *nameOut     ,FILE *fileOut);
  
  void simpleSolverLm(Memoria *m          , PropVarFluid *prop
                    , Loads *loadsVel     , Loads *loadsPres 
                    , Loads *loadsEnergy  , Loads *loadsKturb  
                    , EnergyModel *eModel , Combustion *cModel
                    , MassEqModel *eMass  , MomentumModel *momentumModel
                    , Turbulence *tModel  , ThermoDynamic *thDynamic
                    , Mesh *mesh0         , Mesh *mesh
                    , SistEq *sistEqVel   , SistEq *sistEqPres
                    , SistEq *sistEqEnergy, SistEq *sistEqKturb 
                    , Solv *solvVel       , Solv *solvPres 
                    , Solv *solvEnergy    , Solv *solvKturb      
                    , Simple *sp          , Scheme *sc          
                    , PartMesh *pMesh     , Mean *media
                    , FileOpt *opt        , char *preName
                    , char *nameOut       , FILE *fileOut);

/*...*/
  void updateCellSimpleVel(DOUBLE *RESTRICT w
                ,DOUBLE *RESTRICT u1 ,DOUBLE *RESTRICT u2
                ,INT *RESTRICT id    ,INT const nEl
                ,short const ndm);
  
  void updateCellSimpleVel3D(DOUBLE *RESTRICT w
                ,DOUBLE *RESTRICT u1 ,DOUBLE *RESTRICT u2
                ,DOUBLE *RESTRICT u3     
                ,INT *RESTRICT id    ,INT const nEl
                ,short const ndm);

  void updateCellSimpleVelR(DOUBLE  *RESTRICT w  ,DOUBLE  *RESTRICT u1
                         ,DOUBLE  *RESTRICT u2 ,DOUBLE  *RESTRICT u3
                         ,INT  *RESTRICT id    ,Interface *iNeq
                         ,INT const nEl        ,short const ndm  
                         ,bool const fRes      ,bool const fCom);


  void updateCellSimplePres(DOUBLE  *RESTRICT presC,DOUBLE  *RESTRICT xp   
                         ,INT  *RESTRICT id      ,Interface *iNeq
                         ,INT const nEl          ,bool const fCom);

  void simpleUpdate(DOUBLE *RESTRICT w     ,DOUBLE *RESTRICT pressure
                 ,DOUBLE *RESTRICT PresC ,DOUBLE *RESTRICT GradPresC
                 ,DOUBLE *RESTRICT dField         
                 ,INT const nEl          ,short const ndm
                 ,DOUBLE const alphaPres);
 
  void residualCombustionOld(DOUBLE *RESTRICT vel ,DOUBLE *RESTRICT energy
            ,DOUBLE *RESTRICT zComb
            ,DOUBLE *RESTRICT rCellVel     ,DOUBLE *RESTRICT rCellMass
            ,DOUBLE *RESTRICT rCellEnergy ,DOUBLE *RESTRICT rCellComb      
            ,DOUBLE *RESTRICT adVel       ,DOUBLE *RESTRICT adEnergy 
            ,DOUBLE *RESTRICT adComb   
            ,DOUBLE *RESTRICT rU          ,DOUBLE *rMass
            ,DOUBLE *rEnergy              ,DOUBLE *rComb    
            ,INT  *RESTRICT idVel         ,INT  *RESTRICT idEnergy
            ,INT  *RESTRICT idComb    
            ,INT const nEl                ,INT const nEqVel
            ,INT const nEqComb
            ,short const ndm              ,short const nComb
            ,short iCod);
/*...................................................................*/

/*...*/
  void dynamicDeltat(DOUBLE *RESTRICT vel   , DOUBLE *RESTRICT volume
                  , DOUBLE *RESTRICT density, DOUBLE *RESTRICT sHeat
                  , DOUBLE *RESTRICT tCond  , DOUBLE *RESTRICT dViscosity
                  , Temporal *ddt           , INT const nEl
                  , short const ndm         , short const iCod); 

  void dynamicDeltatChe(DOUBLE *RESTRICT vel    , DOUBLE *RESTRICT volume
                    , DOUBLE *RESTRICT density, DOUBLE *RESTRICT sHeat
                    , DOUBLE *RESTRICT tCond  , DOUBLE *RESTRICT dViscosity
                    , DOUBLE *RESTRICT wk     , Temporal *ddt        
                    , INT const nEl           , short const ns
                    , short const ndm         , short const iCod); 
/*...................................................................*/

/*...*/
  void setSimpleCombustionScheme(char *word , short const ndm
                              , Simple *sp, FILE *fileIn); 
/*...................................................................*/

/*...*/
  void velPresCoupling(Memoria *m         , PropVarFluid *propF
                      , Loads *loadsVel   , Loads *loadsPres
                      , MassEqModel *eMass, MomentumModel *momentumModel
                      , Turbulence *tModel, Mesh *mesh                     
                      , SistEq *sistEqVel , SistEq *sistEqPres
                      , Solv *solvVel     , Solv *solvPres
                      , Simple *sp        , Scheme *sc
                      , PartMesh *pMesh   , DOUBLE *rCellPc
                      , bool *xMomentum   , bool *yMomentum
                      , bool *zMomentum   , bool *pCor
                      , bool fPrint       , short itSimple);
/*...................................................................*/

#endif/*_SIMPLE_H_*/
