#ifndef _ADJCENCY_H
  #define _ADJCENCY_H_
  #include<Memoria.h>
  #include<Mystdbool.h>
  #include<Mesh.h>
  #include<Define.h>
  void viz(Memoria *m ,INT *el   ,INT *nelcon
          ,short *nViz,short *nen ,INT nnode
          ,INT numel ,short maxNo,short maxViz);
  void adj2d(INT *el   ,INT *nodcon ,INT *nelcon
            ,short *nen ,INT numel   , INT nnode 
            ,short maxNo, short maxViz,INT *nEdeg);

#endif/*_MESH_*/

