#ifndef _TRANSPORT_H
  #define _TRANSPORT_H
  
  #include<File.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<CellLoop.h>
  #include<Mesh.h>
  #include<WriteVtk.h>
  #include<Sisteq.h>
  #include<Solv.h>
  #include<Transient.h>


  void transport(Memoria *m ,Loads *loadsTrans
                ,Mesh *mesh0,Mesh *mesh       ,SistEq *sistEqT
                ,Solv *solvT,Scheme sc        ,PartMesh *pMesh
                ,FileOpt opt,char *preName    ,char *nameOut
                ,FILE *fileOut);
#endif/*_DIFFUSION*/
