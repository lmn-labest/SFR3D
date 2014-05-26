#ifndef _WRITE_VTK_H
  #define _WRITE_VTK_H
  #define INTEGER 1
  #define DOUBLE  2
  #include<stdio.h>
  #include<stdlib.h>
  #include<Mystdbool.h>
  #include<Vtk.h>
  #include<Mesh.h>
  #include<Memoria.h>
  #include<File.h>
  void wGeoVtk(Memoria *m        ,double *x      ,int *el 
              ,short int *mat    ,short int *nen ,short int *type
              ,short int *faceRt1,double *faceSt1,long int nnode 
              ,long int numel    ,long int ndm   ,short int maxno
              ,short int *ndfT   ,char *nameOut  ,bool iws      
              ,FILE *f); 
  
#endif/*_WRITE_VTK_H*/
