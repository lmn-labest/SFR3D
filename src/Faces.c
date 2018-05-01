#include<FaceStruct.h>


/*********************************************************************
* Data de criacao    : 08/04/2018                                   *
* Data de modificaco : 00/00/0000                                   *
* ------------------------------------------------------------------*
* faceStruct: estrutura das faces                                   *
* ------------------------------------------------------------------*
* Parametros de entrada:                                            *
* ------------------------------------------------------------------*
* m     -> memoria principal                                        *
* ------------------------------------------------------------------*
* Paramanetros de saida:                                            *
* ------------------------------------------------------------------*
* ------------------------------------------------------------------*
* ------------------------------------------------------------------*
* OBS:                                                              *
* ------------------------------------------------------------------*
*********************************************************************/
void faceStruct(Memoria *m, Mesh *mesh0) {


/*...*/
  HccaAlloc(INT, m, mesh0->face.node
          , mesh0->nFaces*mesh0->maxNo, "faceNode", _AD_);
  HccaAlloc(INT, m, mesh0->face.owner
          , mesh0->nFaces * 2, "fOwner", _AD_);
  /*...................................................................*/
  
/*... gerando as faces*/ 
  makeFaces(m                     , mesh0->elm.node
          , mesh0->elm.adj.nelcon , mesh0->elm.adj.nViz
          , mesh0->elm.cellFace
          , mesh0->face.node      , mesh0->face.owner
          , mesh0->nnode          , mesh0->numel
          , mesh0->nFaces         , mesh0->maxNo   
          , mesh0->maxViz         , mesh0->ndm);  
/*...................................................................*/

/*... alocando */
  HccaAlloc(DOUBLE, m, mesh0->face.ksi
           , mesh0->nFaces*mesh0->ndm, "faceKsi", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh0->face.mksi
           , mesh0->nFaces, "faceModKsi", _AD_);

  HccaAlloc(DOUBLE, m, mesh0->face.eta
           , mesh0->nFaces*mesh0->ndm, "faceEta", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh0->face.area
           , mesh0->nFaces, "faceArea", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh0->face.normal
           , mesh0->nFaces*mesh0->ndm, "faceNormal", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh0->face.xm
           , mesh0->nFaces*mesh0->ndm, "faceXm", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh0->face.vSkew
           , mesh0->nFaces*mesh0->ndm, "faceVskew", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh0->face.mvSkew
           , mesh0->nFaces           , "faceMvSkew", _AD_);
/*...................................................................*/
}
/*********************************************************************/