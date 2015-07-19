#ifndef _ADJCENCY_H
  #define _ADJCENCY_H_
/*...*/
  #include<Memoria.h>
  #include<HccaStdBool.h>
  #include<Mesh.h>
  #include<Define.h>
/*...................................................................*/
  void viz(Memoria *m       ,INT *el           ,INT *nelcon
          ,short *nen       ,INT const nnode   ,INT const numel  
          ,short const maxNo,short const maxViz,short const ndm);
  void adj2d(INT *el          ,INT *nodcon       ,INT *nelcon
            ,short *nen       ,INT const numel   , INT const nnode 
            ,short const maxNo, short const maxViz,INT *nEdeg);
  
  void adjHex8(INT *el         ,INT *nodcon         ,INT *nelcon
            ,INT const nnode    , INT const numel 
            ,short const maxNo,short const maxViz);

  void hexa8fNod(INT const nEl     ,short const face
             ,INT *restrict el   ,INT *restrict nodeFace
             ,short const maxNo);
       
  short hexa8face(INT const nEl         ,INT *restrict el
                 ,INT *restrict nodeFace ,short const maxNo);


#endif/*_MESH_*/

