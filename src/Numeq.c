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
INT numeq(Memoria *m , INT *id   ,INT *num, short *rt
          , short *nen,INT numel  ,short maxViz, short ndf){

  INT i,neq,nel;
  short maxRes = (maxViz + 1)*ndf;
  short j,aux;

/*...*/  
  neq = 0;
  for(i=0;i<numel;i++){
    nel = num[i] -1;
    for(j=0;j<ndf;j++){
      aux = nen[nel]*ndf+j;
      if( MAT2D(nel,aux,rt,maxRes) )
        MAT2D(nel,j,id,ndf) = -1; 
      else
        MAT2D(nel,j,id,ndf) = ++neq;
    } 
  }
/*...................................................................*/

  return neq;
}
