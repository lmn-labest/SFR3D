#include<ElIncid.h>

/*********************************************************************
 * MKELINCID: gera a conectividade de elmento da malha               *
 * ----------------------------------------------------------------- *
 * Parametros de entrada :                                           *
 * ----------------------------------------------------------------- *
 * m       -> vetor de memoria                                       *
 * mesh    -> malha                                                  *
 * pn      -> nao definido                                           *
 * ----------------------------------------------------------------- *
 * Parametros de saida :                                             *
 * ----------------------------------------------------------------- *
 * pn       - conectividade de elmentos da malha                     *
 *********************************************************************/
void mkElIncid(Memoria *m,Pnode *pn,Mesh *mesh){

  INT maxGrade;

/*...*/
  HccaAlloc(INT,m,pn->nincid,mesh->nnode
           ,"nIncid"             ,false);
  zero(pn->nincid,mesh->nnode,INTC);

  nodeGrade(mesh->elm.node,pn->nincid
           ,mesh->elm.nen  ,&maxGrade
           ,mesh->nnode    ,mesh->numel
           ,mesh->maxNo);
/*...................................................................*/
 
/*...*/   
  HccaAlloc(INT,m,pn->incid,mesh->nnode*maxGrade
           ,"incid"        ,false);
  zero(pn->incid,mesh->nnode*maxGrade,INTC);
  elmIncid(mesh->elm.node,pn->nincid
          ,pn->incid     ,mesh->elm.nen
          ,mesh->nnode   ,mesh->numel
          ,maxGrade      ,mesh->maxNo);
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * NODEGRADE: Função do arquivo Colormesh do George. Determina o nume*
 * ro maximo de elementos compartilhando um mesmo no                 *
 * ----------------------------------------------------------------- *
 * Parametros de entrada :                                           *
 * ----------------------------------------------------------------- *
 * el       - conectivida da malha(sem material)                     *
 * nincid   - nao definido                                           *
 * nen      - nos por elementos                                      *
 * maxgrade - nao definido                                           *
 * nnode    - numero de nos                                          *
 * numel    - numero de elementos                                    *
 * maxNo    - numero maximo de nos por elemento                      *
 * ----------------------------------------------------------------- *
 * Parametros de saida :                                             *
 * ----------------------------------------------------------------- *
 * nincid   - numero de incidencias por nos                          *
 * maxgrade - numero de incidencias maximo na malha                  *
 *********************************************************************/
void nodeGrade(INT *restrict el       ,INT *restrict nincid
              ,short *restrict nen    ,INT *maxGrade
              ,INT const nNode        ,INT const numel
              ,short const maxNo){
  INT i,j;
  INT node,grade;

/*...*/
  *maxGrade = 0;
   grade    = 0;
/*... zera o vetor*/
  for(i=0;i<nNode;i++)
    nincid[i] = 0;
/*...................................................................*/

/*...*/
  for(i=0;i<numel;i++){
    for(j=0;j<nen[i];j++){
      node  = MAT2D(i,j,el,maxNo) - 1;
      nincid[node] += 1;
      grade  = nincid[node];
      if( grade > *maxGrade)
        *maxGrade = grade;
    }
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * ELMINCID : Função do arquivo Colormesh do George. Determina       *
 * a concetividade dos elementos                                     *
 * ----------------------------------------------------------------- *
 * Parametros de entrada :                                           *
 * ----------------------------------------------------------------- *
 * el       - conectivida da malha(sem material)                     *
 * incid    - nao definido                                           *
 * nincid   - nao definido                                           *
 * nen      - numero de nos por elementos                            *
 * nnode    - numero de nos                                          *
 * numel    - numero de elementos                                    *
 * maxgrade - numero de incidencias maximo na malha                  *
 * maxNo    - numero maximo de nos por elemento                      *
 * ----------------------------------------------------------------- *
 * Parametros de saida :                                             *
 * ----------------------------------------------------------------- *
 * nincid   - numero de incidencias por nos                          *
 * incid    - incidencia dos elementos por no                        *
 *********************************************************************/
void elmIncid(INT *restrict el    ,INT *restrict incid
               ,INT *restrict nincid,short *restrict nen  
               ,INT const nnode     ,INT const numel
               ,INT const maxGrade  ,short const maxNo){
  INT node,i,j;

/*... zera o vetor*/
  for(i=0;i<nnode;i++){
    nincid[i] = 0;
  }
  for(i=0;i<nnode*maxGrade;i++){
    incid[i] = 0;
  }
/*...................................................................*/

/*...*/
  for(i=0;i<numel;i++){
    for(j=0;j<nen[i];j++){
      node  = MAT2D(i,j,el,maxNo) - 1;
      MAT2D(node,nincid[node],incid,maxGrade) = i; 
      nincid[node] += 1;
    }
  }
/*...................................................................*/
}
/*********************************************************************/

