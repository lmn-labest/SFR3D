#ifndef _PRINT_H_
  #define _PRINT_H_
  
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
/*...*/
  #include<Define.h>
  #include<HccaStdBool.h>
  #include<CellLoop.h>
  #include<Mesh.h>
  #include<Memoria.h>
  #include<File.h>
  #include<WriteVtk.h>
  #include<ParallelMpi.h>
  #include<Media.h>

  void printFluid(Memoria *m        
                , Turbulence *turbModel, EnergyModel *eModel
                , PartMesh *pMesh      , Scheme sc
                , Loads *loadsVel      , Loads *loadsPres 
                , Loads *loadsTemp     , FileOpt opt
                , Mesh *mesh0          , Mesh *mesh  
                , Mean *media          
                , char *preName        , char *nameOut);
  
  void printDiff(Memoria *m
               , PartMesh *pMesh, Scheme *sc
               , Loads *loadsD1 , FileOpt *opt
               , Mesh *mesh0    , Mesh *mesh
               , char *preName  , char *nameOut);

  void reScaleMesh(DOUBLE *x, INT nnode, short ndm, FILE *fileIn);


#endif/*_WRITE_VTK_H_*/