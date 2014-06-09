#include<Csr.h>
/********************************************************************* 
 * CSRIA:                                                            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia     -> indefinido                                              * 
 * id     -> numeracao das equacoes por elemento                     * 
 * num    -> renumeracao dos elementos                               * 
 * adj    -> adjacencia dos elementos                                * 
 * nViz   -> numero de vizinhos por elemento                         * 
 * numel  -> numero de elementos                                     * 
 * neq    -> numero de equacoes                                      * 
 * maxViz -> numero maximo de vizinho da malha                       * 
 * upper  -> armazenamento da parte superior da matriz (CSR/CSRC)    * 
 * diag   -> armazenamento da diagonal (CSR/CSRC)                    * 
 * lower  -> armazenamenro da parte inferior (CSR/CSRC)              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ia     -> ponteiro do CSR                                         * 
 *-------------------------------------------------------------------* 
 * OBS: a funcao retorna o numero do termos nao nulor no CSR/CSRC    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
long csrIa(long *ia  ,long *id    ,long *num   ,long  *adj, short *nViz
          ,long numel,long neq    ,short maxViz, bool upper
          ,bool diag , bool lower){
  
  long  i,nel1,neq1,viz1,col,aux;
  short j;
/*... gerando arranjo ia*/  
  ia[0] = 0;
  for(i=0;i<numel;i++){
    nel1= num[i]-1;
    aux = 0;
    neq1= id[nel1]-1;
    if(neq1 != -2){
/*...*/
      for(j=0;j<nViz[nel1];j++){
        viz1 = NELCON(nel1,j,adj,maxViz) - 1;
        if( viz1 != -2) {
          col   = id[viz1]-1;
          if( col != -2){
/*... parte superior*/
            if(lower && col < neq1) 
              aux++;
/*...................................................................*/

/*... parte inferior*/            
            else if(upper && col > neq1)
              aux++;
          }    
/*...................................................................*/
        }
      }
/*...................................................................*/

/*... diagonal princial*/      
      if(diag) aux++;
      ia[neq1+1] = ia[neq1] + aux;
/*...................................................................*/
    }
  }
/*...................................................................*/
  return ia[neq] - ia[0];
}
/*********************************************************************/ 


/********************************************************************* 
 * CSRJA:                                                            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia     -> arranjo CSR/CSRC                                        * 
 * ja     -> indefinido                                              * 
 * id     -> numeracao das equacoes por elemento                     * 
 * num    -> renumeracao dos elementos                               * 
 * adj    -> adjacencia dos elementos                                * 
 * nViz   -> numero de vizinhos por elemento                         * 
 * numel  -> numero de elementos                                     * 
 * neq    -> numero de equacoes                                      * 
 * maxViz -> numero maximo de vizinho da malha                       * 
 * upper  -> armazenamento da parte superior da matriz (CSR/CSRC)    * 
 * diag   -> armazenamento da diagonal (CSR/CSRC)                    * 
 * lower  -> armazenamenro da parte inferior (CSR/CSRC)              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ja     -> ponteiro do CSR                                         * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void csrJa(long *ia    ,long *ja 
          ,long *id    ,long *num   ,long  *adj, short *nViz
          ,long numel,long neq    ,short maxViz, bool upper
          ,bool diag , bool lower){
  
  long  i,nel1,neq1,viz1,col,aux,ipont;
  short j;

/*... gerando arranjo ia*/  
  for(i=0;i<numel;i++){
    nel1= num[i]-1;
    aux = 0;
    neq1= id[nel1]-1;
    ipont = ia[neq1];
    if(neq1 != -2){
/*...*/
      for(j=0;j<nViz[nel1];j++){
        viz1 = NELCON(nel1,j,adj,maxViz) - 1;
        if( viz1 != -2) {
          col   = id[viz1]-1;
          if( col != -2){
/*... parte superior*/
            if(lower && col < neq1){
              ja[ipont+aux] = col; 
              aux++;
            }
/*...................................................................*/

/*... parte inferior*/            
            else if(upper && col > neq1){
              ja[ipont+aux] = col; 
              aux++;
            }
/*...................................................................*/
          }    
        }
      }
/*...................................................................*/

/*... diagonal princial*/      
      if(diag) {
        ja[ipont+aux] = neq1;  
        aux++;
      }
/*...................................................................*/
    }
  }
/*...................................................................*/
}
/*********************************************************************/ 


/********************************************************************* 
 * SORTGRAPHCSR: ordena o gafo no formato CSR                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia     -> arranjo CSR/CSRC                                        * 
 * ja     -> indefinido                                              * 
 * n      -> numera de linhas                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ja     -> ponteiro do CSR                                         * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/ 
void sortGraphCsr(long *ia,long *ja,long n){

  long i,nl;
  
  for(i=0;i<n;i++){
    nl = ia[i+1] - ia[i];
    if(nl!=0) bubblesort(&ja[ia[i]],nl);
  }
}

