#ifndef _ADJCENCY_H
  #define _ADJCENCY_H_
  #include<Memoria.h>
  #include<Mystdbool.h>
  #include<Mesh.h>
  void viz(Memoria *m ,long *el   ,long *nelcon
          ,short *nViz,short *nen ,long nnode
          ,long numel ,short maxNo,short maxViz);
  void adj2d(long *el   ,long *nodcon ,long *nelcon
            ,short *nen ,long numel   , long nnode 
            ,short maxNo, short maxViz,long *nEdeg);

#endif/*_MESH_*/

