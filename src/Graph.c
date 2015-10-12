#include<Graph.h>
/********************************************************************* 
 * CONVGRAPH: grafo da malha no formato CSR                          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * xAdj     -> nao definido                                          * 
 * adjncy   -> nao definido                                          * 
 * adj      -> vininhacas das celulas                                * 
 * nViz     -> numero de vizinhos                                    * 
 * maxViz   -> numero maximo de vizinhos                             * 
 * numel    -> numero de celulas/elementos                           * 
 * xAdjFlag -> monta xAdj                                            * 
 * adjFlag  -> monta adjncy                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * xAdj   -> ponteiro do celula i (CSR ia)                           * 
 * adjncy -> numeracao das celulas vizinhas a celula i  (CSR ja)     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/ 
void convGraph(INT *restrict xAdj     ,INT *restrict adjncy 
              ,INT const *adj         ,short const *nViz     
              ,short const maxViz     ,INT const numel
              ,bool const xAdjFlag    ,bool const adjFlag){

  INT nel,kk=0,viz,n;
  int aux = 0;
  
  if(xAdj) xAdj[0] = 0;
  for(nel=0;nel<numel;nel++){
    for(viz=0;viz<nViz[nel];viz++){
      n = MAT2D(nel,viz,adj,maxViz);
      if( n != -1){
        aux++;
        if(adjFlag) 
          adjncy[kk++] = n - 1;  
      }
    } 
    if(xAdjFlag) xAdj[nel+1] = xAdj[nel] + aux;
    aux = 0;
  }

}
/*********************************************************************/ 

/********************************************************************* 
 * CONVGRAPHPART : grafo da malha particionada no formato CSR        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * xAdj     -> nao definido                                          * 
 * adjncy   -> nao definido                                          * 
 * adj      -> vininhacas das celulas                                * 
 * nViz     -> numero de vizinhos                                    * 
 * maxViz   -> numero maximo de vizinhos                             * 
 * numelNov -> numero de celulas/elementos sem sobreposicao          * 
 * xAdjFlag -> monta xAdj                                            * 
 * adjFlag  -> monta adjncy                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * xAdj   -> ponteiro do celula i (CSR ia)                           * 
 * adjncy -> numeracao das celulas vizinhas a celula i  (CSR ja)     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/ 
void convGraphPart(INT *restrict xAdj     ,INT *restrict adjncy 
                  ,INT const *adj         ,short const *nViz     
                  ,short const maxViz   
                  ,INT const numelNov
                  ,bool const xAdjFlag    ,bool const adjFlag){

  INT nel,kk=0,viz,n;
  int aux = 0;
  
  if(xAdj) xAdj[0] = 0;
  for(nel=0;nel<numelNov;nel++){
    for(viz=0;viz<nViz[nel];viz++){
      n = MAT2D(nel,viz,adj,maxViz);
      if( n != -1 && n < numelNov){
        aux++;
        if(adjFlag) 
          adjncy[kk++] = n - 1;  
      }
    } 
    if(xAdjFlag) xAdj[nel+1] = xAdj[nel] + aux;
    aux = 0;
  }

}

/*********************************************************************/ 

/********************************************************************* 
 * SORTGRAPHCSR: ordena o gafo no formato CSR                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia     -> arranjo CSR/CSRC/CSRD                                   * 
 * ja     -> colunas do CSR/CSRC/CSRD                                * 
 * n      -> numera de linhas                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ja     -> colunas do CSR                                          * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/ 
void sortGraphCsr(INT *restrict ia,INT *restrict ja,INT const n){

  INT i,nl;
  
  for(i=0;i<n;i++){
    nl = ia[i+1] - ia[i];
    if(nl>1) bubblesort(&ja[ia[i]],nl);
  }
}
/*********************************************************************/ 

