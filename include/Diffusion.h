#ifndef _DIFFUSION_H_
  #define _DIFFUSION_H_
  
  #include<File.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<CellLoop.h>
  #include<Mesh.h>
  #include<WriteVtk.h>
  #include<Sisteq.h>
  #include<Solv.h>
  #include<Transient.h>


  void diffusion(Memoria *m ,Loads *loadsDif
                ,Mesh *mesh0,Mesh *mesh     ,SistEq *sistEqD
                ,Solv *solvD,Scheme sc      ,PartMesh *pMesh
                ,FileOpt opt,char *preName  ,char *nameOut
                ,FILE *fileOut);
#endif/*_DIFFUSION_H_*/
