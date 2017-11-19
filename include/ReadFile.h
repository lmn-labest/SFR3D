#ifndef _READ_FILE_
  #define _READ_FILE_
  #define NPARAMETROS   6
  #define NCONFIG       6
  #ifdef NMACROS 
    #undef NMACROS
  #endif  
  #define NMACROS       39
/*...................................................................*/  
  #include<stdio.h>
  #include<stdlib.h>
  #include<ctype.h>
/*...*/
  #include<Mesh.h>
  #include<CellLoop.h>
  #include<string.h>
  #include<File.h>
  #include<HccaStdBool.h>
  #include<HccaBlas.h>
  #include<Memoria.h>
  #include<ParallelMpi.h>
  #include<Properties.h>
/*...................................................................*/  

  void parametros(INT   *nnode,INT *nel    
                 ,short *maxNo,short *maxViz
                 ,short *ndm  ,short *numat
                 ,FILE  *file);
  void readFileFvMesh(Memoria *m  ,Mesh *mesh
                     ,PropVar prop,EnergyModel eModel
                     ,FILE *file);

  void readVfMat(DOUBLE *prop,short *type,short numat,FILE *file);
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

  void config(FileOpt *opt ,Reord *reord
             ,short *rcGrad
             ,FILE* f);
  
  void readEdo(Mesh *mesh,FILE *file);
  void readPropVar(PropVar *p,FILE *file);
  void readGravity(DOUBLE *gravity,FILE *file);
  void readModel(EnergyModel *e    , Turbulence *t
               , MassEqModel *eMass, MomentumModel *eMomentum
               , FILE *file);
  void setPrintFluid(FileOpt *opt,FILE *file);

  void initProp(DOUBLE *RESTRICT prop 
             ,DOUBLE *RESTRICT propMat,short *RESTRICT mat
             ,short const np          ,INT const nCell 
             ,short const iProp);
  
   void uniformField(DOUBLE *field, INT const n, short const ndf
                    , FILE* file);
 
   void convStringLower(char *s);
   void help(FILE *f);
#endif  /*_READ_FILE_*/
