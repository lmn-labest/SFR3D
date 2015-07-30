#ifndef _READ_FILE_
  #define _READ_FILE_
  #define NPARAMETROS   7
  #define NCONFIG       6
  #ifdef NMACROS 
    #undef NMACROS
  #endif  
  #define NMACROS       24
/*...................................................................*/  
  #include<stdio.h>
  #include<stdlib.h>
  #include<Mesh.h>
  #include<string.h>
  #include<File.h>
  #include<HccaStdBool.h>
  #include<Memoria.h>
  void parametros(INT   *nnode,INT *nel    ,short *maxNo,short *maxViz
                 ,short *ndm  ,short *numat,short  *ndf ,FILE  *file);
  void readFileFvMesh(Memoria *m,Mesh *mesh, FILE *file);
  void readVfMat(double *prop,short *type,short numat,FILE *file);
  void readVfCoor(double *x,INT nn, short ndm,FILE *file);
  void readVfElmt(INT *el    ,short *mat ,short *nen,short *nFace
                 ,short *ty  ,INT  nel   ,short maxno
                 ,FILE *file);
  void readVfRes(short *id,INT numel,short maxno,char *str,FILE *file);
  void readVfSource(double *f    ,INT numel,short maxno
                   ,char *str,FILE *file);
  void readVfLoads(Loads *loads,char *str,FILE* file);

  void config(FileOpt *opt ,Reord *reord
             ,short *rcGrad
             ,FILE* f);

#endif
