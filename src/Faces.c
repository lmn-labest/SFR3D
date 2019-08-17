#include<FaceStruct.h>

/*********************************************************************
 * Data de criacao    : 08/04/2018                                   *
 * Data de modificaco : 12/08/2019                                   *
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
 * Nao e preciso calcular as faces das celulas em sobreposicao       *
 * ------------------------------------------------------------------*
 *********************************************************************/
void faceStruct(Memoria *m, Mesh *mesh) {


/*...*/
  HccaAlloc(INT, m, mesh->face.node
          , mesh->nFaces*mesh->maxNo, "faceNode", _AD_);
  HccaAlloc(INT, m, mesh->face.owner
          , mesh->nFaces * 2, "fOwner", _AD_);
/*...................................................................*/

/*... gerando as faces*/ 
   makeFaces(m                    , mesh->elm.node
          , mesh->elm.adj.nelcon , mesh->elm.adj.nViz
          , mesh->elm.cellFace
          , mesh->face.node      , mesh->face.owner
          , mesh->nnode          , mesh->numelNov
          , mesh->nFaces         , mesh->maxNo   
          , mesh->maxViz         , mesh->ndm);  
/*...................................................................*/

/*... alocando */
  HccaAlloc(DOUBLE, m, mesh->face.ksi
           , mesh->nFaces*mesh->ndm, "faceKsi", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh->face.mksi
           , mesh->nFaces, "faceModKsi", _AD_);

  HccaAlloc(DOUBLE, m, mesh->face.eta
           , mesh->nFaces*mesh->ndm, "faceEta", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh->face.area
           , mesh->nFaces, "faceArea", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh->face.normal
           , mesh->nFaces*mesh->ndm, "faceNormal", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh->face.xm
           , mesh->nFaces*mesh->ndm, "faceXm", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh->face.vSkew
           , mesh->nFaces*mesh->ndm, "faceVskew", _AD_);
  
  HccaAlloc(DOUBLE, m, mesh->face.mvSkew
           , mesh->nFaces           , "faceMvSkew", _AD_);
/*...................................................................*/
}
/*********************************************************************/