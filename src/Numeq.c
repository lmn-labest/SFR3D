#include<Sisteq.h>
/*********************************************************************
 * numEq : numeracao das equacoes                                    *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * m     -> memoria principal                                        * 
 * id    -> indefinido                                               *
 * num   -> renumeracao dos elementos                                *
 * rt    -> restricoes                                               * 
 * nFace -> numero de faces por elemento                             * 
 * numel -> numero de elementos                                      * 
 * maxViz-> numero maximo de vizinhos                                * 
 * str   -> nome do vetor                                            * 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * id    -> numerocao das equacoes por celula                        *
 * ------------------------------------------------------------------*
 * *******************************************************************/
INT numeq(INT *restrict id       ,INT *restrict num
         ,short *restrict rt     ,short *restrict nFace 
         ,INT const numel        ,short const maxViz  
         ,short const ndf){

  INT i,neq,nel;
  short maxRes = (maxViz + 1)*ndf;
  short j,aux;

/*...*/  
  neq = 0;
  for(i=0;i<numel;i++){
    nel = num[i] -1;
    for(j=0;j<ndf;j++){
      aux = nFace[nel]*ndf+j;
      if( MAT2D(nel,aux,rt,maxRes) == 1)
        MAT2D(nel,j,id,ndf) = -1; 
      else
        MAT2D(nel,j,id,ndf) = ++neq;
    } 
  }
/*...................................................................*/

  return neq;
}
