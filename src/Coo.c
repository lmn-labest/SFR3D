#include<Coo.h>
/********************************************************************* 
 * COONNZ: conta o numero total de nao zeros na marte restangular da * 
 * matriz                                                            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * id     -> numeracao das equacoes por elemento                     * 
 * num    -> renumeracao dos elementos                               * 
 * adj    -> adjacencia dos elementos                                * 
 * nViz   -> numero de vizinhos por elemento                         * 
 * numel  -> numero de elementos                                     * 
 * neq    -> numero de equacoes                                      * 
 * ndf    -> numero de graus de liberade                             * 
 * maxViz -> numero maximo de vizinho da malha                       * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ia     -> ponteiro do CSRC-CSR                                    * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
INT cooNnzR(INT *restrict id   
           ,INT *restrict num    ,INT *restrict adj
           ,short *restrict nViz
           ,INT const numel      ,INT const nEq
           ,short const maxViz   ,short  const ndf){
  
  INT  i,nel1,neq,neqi,neqj,viz1,nnz=0;
  short jNdf,kNdf,j;
/*... gerando arranjo ia*/  
  
  neq = nEq-1;
  for(i=0;i<numel;i++){
    nel1= num[i]-1;
    for(jNdf=0;jNdf<ndf;jNdf++){
      neqi= MAT2D(nel1,jNdf,id,ndf)-1;
      if(neqi != -2){
/*... equacoes ligadas as celulas vizinhas*/
        for(j=0;j<nViz[nel1];j++){
          viz1 = MAT2D(nel1,j,adj,maxViz) - 1;
          if( viz1 != -2) {
            for(kNdf=0;kNdf<ndf;kNdf++){
              neqj   = MAT2D(viz1,kNdf,id,ndf)-1;
              if( neqj != -2 && neqj > neq) nnz++;
            }
/*...................................................................*/
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
  return nnz;
}
/*********************************************************************/ 

/********************************************************************* 
 * COOIAJAR: gera os vetores ia e ja da parte retangular da matriz no* 
 * formato COO e um vertor auxiliar ia0 semelhante au vetor ia do    *
 * CSRD                                                              *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia0    -> nao definido                                            * 
 * ia     -> nao definido                                            * 
 * ja     -> nao definido                                            * 
 * id     -> numeracao das equacoes por elemento                     * 
 * num    -> renumeracao dos elementos                               * 
 * adj    -> adjacencia dos elementos                                * 
 * nViz   -> numero de vizinhos por elemento                         * 
 * numel  -> numero de elementos                                     * 
 * neq    -> numero de equacoes                                      * 
 * ndf    -> numero de graus de liberade                             * 
 * maxViz -> numero maximo de vizinho da malha                       * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ia0    -> ponteiro para linas do COO no formato CSRD              * 
 * ia     -> linhas  COO                                             * 
 * ja     -> colunas COO                                             * 
 *-------------------------------------------------------------------* 
 * OBS: ia0 e usado apenas para auxiliar a montagem da matriz global * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cooIaJaR(INT *restrict ia0
             ,INT *restrict ia     ,INT *restrict ja
             ,INT *restrict id     ,INT *restrict num   
             ,INT *restrict adj    ,short *restrict nViz
             ,INT const numel      ,INT const nEq
             ,short const maxViz   ,short  const ndf){
  
  INT  i,nel1,neq,viz1,neqi,neqj,aux=0,aux1=0;
  short jNdf,kNdf,j;

  for(i=0;i<nEq+1;i++)
    ia0[i] = 0;
 
/*... gerando arranjo ia*/  
  neq = nEq-1;
  for(i=0;i<numel;i++){
    nel1= num[i]-1;
    for(jNdf=0;jNdf<ndf;jNdf++){
      aux1 = 0;
      neqi= MAT2D(nel1,jNdf,id,ndf)-1;
      if(neqi != -2){
/*... equacoes ligadas ao proprio elemento
        (mais de um grau de liberdade por celula)*/
        for(kNdf=0;kNdf<ndf;kNdf++){
          neqj = MAT2D(nel1,kNdf,id,ndf)-1;
          if( neqj != -2 && neqj > neq) { 
            ia[aux] = neqi;
            ja[aux] = neqj;
            aux++;
            aux1++;
          }
        }
/*...................................................................*/

/*... equacoes ligadas as celulas vizinhas*/
        for(j=0;j<nViz[nel1];j++){
          viz1 = MAT2D(nel1,j,adj,maxViz) - 1;
          if( viz1 != -2) {
            for(kNdf=0;kNdf<ndf;kNdf++){
              neqj  = MAT2D(viz1,kNdf,id,ndf)-1;
              if( neqj!= -2 && neqj > neq){
                ia[aux] = neqi;
                ja[aux] = neqj;
                aux++;
                aux1++;
              }  
            }
/*...................................................................*/
          }
/*...................................................................*/
        }
/*...................................................................*/
        ia0[neqi+1]=ia0[neqi] + aux1;
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

