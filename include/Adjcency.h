#ifndef _ADJCENCY_H_
  #define _ADJCENCY_H_
/*...*/
  #include<Memoria.h>
  #include<HccaStdBool.h>
  #include<Define.h>
/*...................................................................*/

/*...*/
  void neighbors(Memoria *m         , INT *RESTRICT el
              , INT *RESTRICT nelcon, short *RESTRICT nViz
              , short *RESTRICT nen , INT *nFace
              , INT nnode           , INT numel
              , short maxNo         , short maxViz
              , short ndm);
/*...................................................................*/

/*...*/
  void makeFaces(Memoria *m             , INT *RESTRICT el
                , INT *RESTRICT nelcon  , short *RESTRICT nViz
                , INT *RESTRICT cellFace
                , INT *RESTRICT fNode   , INT *RESTRICT owner
                , INT nnode             , INT numel
                , INT nFaces
                , short maxNo           , short maxViz
                , short ndm);
/*...................................................................*/

/*... malha mista de triangulos e quadrilateros*/
  void adj2d(INT *el          ,INT *nodcon       ,INT *nelcon
            ,short *nen       ,INT const numel   ,INT const nnode 
            ,short const maxNo,short const maxViz,INT *nEdeg);
/*...................................................................*/

/*... malha mista de tetraedros, piramide e hexaedros*/
  void adj3d(INT *el          , INT *nodcon       , INT *nelcon
           , short *nViz      , INT const nnode   , INT const numel
           , short const maxNo, short const maxViz, INT *nFace);
/*...................................................................*/

/*... malha de tetraedros*/
  void adjTetra4(INT *RESTRICT el         ,INT *RESTRICT nodcon
            ,INT *RESTRICT nelcon
            ,INT const nnode    , INT const numel 
            ,short const maxNo,short const maxViz);

  void tetra4fNod(INT const nEl     ,short const face
             ,INT *RESTRICT el   ,INT *RESTRICT nodeFace
             ,short const maxNo);
       
  short tetra4face(INT const nEl         ,INT *RESTRICT el
                 ,INT *RESTRICT nodeFace ,short const maxNo);
/*...................................................................*/
  
/*... malha de hexaedros*/
  void adjHex8(INT *el         ,INT *nodcon         ,INT *nelcon
            ,INT const nnode   , INT const numel 
            ,short const maxNo,short const maxViz   ,INT *nFace);

  void hexa8fNod(INT const nEl     ,short const face
             ,INT *RESTRICT el   ,INT *RESTRICT nodeFace
             ,short const maxNo);
       
  short hexa8face(INT const nEl         ,INT *RESTRICT el
                 ,INT *RESTRICT nodeFace ,short const maxNo);
/*...................................................................*/

/*... malha de piramides*/
  void pira5fNod(INT const nEl   , short const face
               , INT *RESTRICT el, INT *RESTRICT nodeFace
               , short const maxNo);

  short pira5face(INT const nEl         , INT *RESTRICT el
                , INT *RESTRICT nodeFace, short const maxNo);
  void pira5Aux1(short idFace, INT *nodcon, INT *node, INT nel1
               , bool *imiss);
  void pira5Aux2(short k     , short j
               , short maxNo , short maxViz
               , INT *el     , INT *nodcon, INT *nelcon
               , INT *node   , INT nel1   , INT nel2
               , INT *nFace  , bool *imiss);
  void pira5Aux3(short j    , short maxViz
               , INT *nodcon, INT *nelcon , INT *node
               , INT nel1   , INT nel2    , INT *nFace
               , bool *imiss);
/*...................................................................*/
#endif/*_ADJCENCY_H_*/

