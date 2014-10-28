#ifndef _WRITE_VTK_H
  #define _WRITE_VTK_H
  #define INTEGER 1
  #define DOUBLEV 2
  #include<stdio.h>
  #include<stdlib.h>
  #include<HccaStdBool.h>
  #include<Vtk.h>
  #include<Mesh.h>
  #include<Memoria.h>
  #include<File.h>
  
  void wGeoVtk(Memoria *m     ,double *x      
              ,INT *el        ,short *mat    
              ,short *nen     ,short *typeGeom
              ,double *prop   ,short *typeCal
              ,short *faceRt1 ,double *faceSt1
              ,INT nnode      ,INT numel    
              ,short ndm        ,short maxno
              ,short numat    ,short *ndfT   
              ,char *nameOut  ,bool iws      
              ,FILE *f); 
  
#endif/*_WRITE_VTK_H*/
