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
  void wGeoVtk(Memoria *m    ,double *x      ,INT *el 
              ,short *mat    ,short int *nen ,short *type
              ,short *faceRt1,double *faceSt1,INT nnode 
              ,INT numel     ,INT ndm        ,short maxno
              ,short *ndfT   ,char *nameOut  ,bool iws      
              ,FILE *f); 
  
#endif/*_WRITE_VTK_H*/
