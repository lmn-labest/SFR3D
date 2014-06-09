#include<Sisteq.h>
/*********************************************************************
 * numEq : numeracao da equacoes                                     *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * m     -> memoria principal                                        * 
 * id    -> indefinido                                               *
 * num   -> renumeracao dos elementos                                *
 * rt    -> restricoes                                               * 
 * nen   -> numero de aresta por elemento                            * 
 * numel -> numero de elementos                                      * 
 * maxViz-> numero maximo de vizinhos                                * 
 * str   -> nome do vetor                                            * 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * id    -> numerocao das equacoes por celula                        *
 * ------------------------------------------------------------------*
 * *******************************************************************/
long numeq(Memoria *m , long *id   ,long *num, short *rt
          , short *nen,long numel,short maxViz, char *str){

  long i,neq,nel;

/*...*/  
  neq = 0;
  for(i=0;i<numel;i++){
    nel = num[i] -1; 
    if( VET(nel,nen[nel],rt,maxViz) )
      id[nel] = -1; 
    else
      id[nel] = ++neq;
  }
/*...................................................................*/

  return neq;
}
