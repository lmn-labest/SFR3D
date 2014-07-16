#ifndef _READ_FILE_
  #define _READ_FILE_
  #define NPARAMETROS   6
  #define NCONFIG       3
  #ifdef NMACROS 
    #undef NMACROS
  #endif  
  #define NMACROS       7
/*...................................................................*/  
  #include<stdio.h>
  #include<stdlib.h>
  #include<Mesh.h>
  #include<string.h>
  #include<File.h>
  #include<Mystdbool.h>
  #include<Memoria.h>
  void parametros(INT   *nnode,INT *nel   ,short *maxno,short *ndm
                 ,short *numat,short  *ndf,FILE  *file);
  void readFileFvMesh(Memoria *m,Mesh *mesh, FILE *file);
  void readVfCoor(double *x,INT nn, short ndm,FILE *file);
  void readVfElmt(INT *el    ,short *mat ,short *nen,short *nFace
                 ,short *ty  ,INT  nel   ,short maxno
                 ,FILE *file);
  void readVfRes(short *id,INT numel,short maxno,char *str,FILE *file);
  void readVfSource(double *f ,INT numel,short maxno
                   ,char *str ,FILE *file);
  void config(bool *bvtk,Reord *reord,FILE* f);
#endif
