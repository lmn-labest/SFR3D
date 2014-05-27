#ifndef _READ_FILE_
  #define _READ_FILE_
  #define NPARAMETROS   6
  #define NCONFIG       1
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
  void parametros(long int*nnode,long int *nel,short *maxno,short *ndm
                 ,short   *numat,short    *ndf,FILE  *file);
  void readFileFv(Memoria *m,Mesh *mesh, FILE *file);
  void readVfCoor(double *x,long int nn, short int ndm,FILE *file);
  void readVfElmt(long int *el    ,short int *mat ,short int *nen
                 ,short int *ty  ,long int nel   ,short int maxno
                 ,FILE *file);
  void readVfRes(short int *id,long int numel,short int maxno
                ,char *str    ,FILE *file);
  void readVfSource(double *f ,long int numel,short int maxno
                   ,char *str ,FILE *file);
  void config(bool *bvtk,FILE* f);
#endif
