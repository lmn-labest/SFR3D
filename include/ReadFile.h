#ifndef _READ_FILE_
  #define _READ_FILE_
  #define NPARAMETROS   6
  #define NCONFIG       4
  #ifdef NMACROS 
    #undef NMACROS
  #endif  
  #define NMACROS       39
/*...................................................................*/  
  #include<stdio.h>
  #include<stdlib.h>
/*...*/
  #include<Mesh.h>
  #include<CellLoop.h>
  #include<File.h>
  #include<HccaStdBool.h>
  #include<HccaBlas.h>
  #include<Memoria.h>
  #include<ParallelMpi.h>
  #include<Properties.h>
  #include<Simple.h>
  #include<Prime.h>
  #include<Print.h>
/*...................................................................*/  

  void parametros(INT   *nnode,INT *nel    
                 ,short *maxNo,short *maxViz
                 ,short *ndm  ,short *numat
                 ,FILE  *file);
  
  void readFileFvMesh( Memoria *m              , Mesh *mesh
                   , PropVarFluid *propF           
                   , PropVarCD *propD        , PropVarCD *propT
                   , EnergyModel *energyModel, Turbulence *tModel     
                   , Combustion *cModel      , Mean *media
                   , FILE* file);

  void readVfMat(DOUBLE *prop,short *type,short numat,FILE *file);

  void readFileMat(DOUBLE *prop, short *type, short numat, FILE *file);

  void readVfCoor(DOUBLE *x,INT nn, short ndm,FILE *file);
  void readVfElmt(INT *el    ,short *mat ,short *nen,short *nFace
                 ,short *ty  ,INT  nel   ,short maxno
                 ,FILE *file);
  void readVfRes(short *id,INT numel,short maxno,char *str,FILE *file);
  void readVfSource(DOUBLE *f    ,INT numel,short const maxCarg
                   ,char *str,FILE *file);
  void readVfInitial(DOUBLE *f    ,INT numel,short const ndf
                    ,char *str,FILE *file);
  void readVfLoads(Loads *loads,char *str,FILE* file);

  void config(FileOpt *opt ,Reord *reord,FILE* f);
  
  void readEdo(Mesh *mesh,FILE *file);  
  void readPropVar(PropVarFluid *pf, PropVarCD *pd, PropVarCD *pt,Combustion *cModel
                 , FILE *file);
  void readPropVarDiff(PropVarCD *p, FILE *file);
  void readPropVarTrans(PropVarCD *p, FILE *file);
  void readPropVarFluid(PropVarFluid *p, FILE *file);
  void readPropVarMixture(PropVarFluid *p,Combustion *cModel,FILE *file);
  void readGravity(DOUBLE *gravity,FILE *file);
  void readModel(EnergyModel *e    , Turbulence *t
               , MassEqModel *eMass, MomentumModel *ModelMomentum
               , DiffModel *dModel , TransModel *tModel
               , Combustion *cModel
               , FILE *file);
  void readMean(Memoria *m, FILE *fileIn
              , Mesh *mesh, Mean *media);
  void setMixedModelLes(Turbulence *t  , FILE *file);
  void setDynamicModelLes(Turbulence *t, FILE *file);
  void setEdc(Edc *e       , FILE *file);
  void setEdm(Edm *e       , FILE *file); 
  void setPrint(FileOpt *opt,FILE *file);
  void setReGrad(short *rcGrad, FILE *fileIn);

  void initProp(DOUBLE *RESTRICT prop 
             ,DOUBLE *RESTRICT propMat,short *RESTRICT mat
             ,short const np          ,INT const nCell 
             ,short const iProp);
  
   void uniformField(DOUBLE *field, INT const n, short const ndf
                    , FILE* file);
   void readAdvectionScheme(FILE *fileIn, Scheme *sc);
   void readDiffusionScheme(FILE *fileIn, Scheme *sc);
   void readSetSimple(Memoria *m    , FILE *fileIn
                 , Mesh *mesh0   , Mesh *mesh
                 , Simple *simple, bool *fSolvSimple);
   void readSetSimpleComb(Memoria *m    , FILE *fileIn
                     , Mesh *mesh0   , Mesh *mesh
                    , Simple *simple, bool *fSolvComb);
/*...*/
   void readSolvFluid(Memoria *m, Mesh *mesh
                    , Reord *reordMesh
                    , Solv *solvVel, SistEq* sistEqVel, bool *fSolvVel
                    , Solv *solvPres, SistEq* sistEqPres, bool *fSolvPres
                    , Solv *solvEnergy, SistEq* sistEqEnergy, bool *fSolvEnergy
                    , Solv *solvKturb, SistEq* sistEqKturb, bool *fSolvKturb
                    , char* auxName, char* preName, char* nameOut
                    , FILE *fileIn, FileOpt *opt);
/*...*/
  void readSolvDiff (Memoria *m   , Mesh *mesh      , Reord *reordMesh
                   , Solv *solvD1 , SistEq* sistEqD1, bool *fSolvD1
                   , char* auxName, char* preName   , char* nameOut
                   , FILE *fileIn , FileOpt *opt);

  void readSolvTrans(Memoria *m   , Mesh *mesh      , Reord *reordMesh
                   , Solv *solvT1 , SistEq* sistEqT1, bool *fSolvT1
                   , char* auxName, char* preName   , char* nameOut
                   , FILE *fileIn , FileOpt *opt);

  void readSolvComb(Memoria *m      , Mesh *mesh          , Reord *reordMesh
                , Combustion *cModel 
                , Solv *solvVel   , SistEq* sistEqVel   , bool *fSolvVel
                , Solv *solvPres  , SistEq* sistEqPres  , bool *fSolvPres
                , Solv *solvEnergy, SistEq* sistEqEnergy, bool *fSolvEnergy
                , Solv *solvKturb , SistEq* sistEqKturb , bool *fSolvKturb
                , Solv *solvComb  , SistEq* sistEqComb  , bool *fSolvComb
                , PartMesh *pMesh 
                , char* auxName   , char* preName       , char* nameOut
                , FILE *fileIn    , FileOpt *opt);

   void readChemical(Combustion *c, FILE *file);
   void readResidual(Simple *sc, FILE *file);
   void readNlIt(Scheme *sc, FILE *fileIn);
   
   void readSetPrime(Memoria *m   , FILE *fileIn
                   , Mesh *mesh0  , Mesh *mesh
                   , Prime  *prime, bool *fSolvPrime);

   bool readBool(char *word);
   void typeResidual(char *word,Residual *re);
/*...................................................................*/
   void help(FILE *f);
#endif  /*_READ_FILE_*/
