#include<Reord.h>
/********************************************************************* 
 * REORD : reordanacao do grafo dos elementos                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * num      -> indefinido                                            * 
 * graph    -> grafo dos elementos                                   * 
 * numel    -> numero de elementos total                             * 
 * numelNov -> numero de elementos sem sobreposicao                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * num   -> elementos renumerados                                    * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *********************************************************************/
void reord(Memoria *m            ,INT *num           ,INT const *adj
          ,short const*nViz      ,short const maxViz  
          ,INT const numel       ,INT const numelNov            
          ,bool const flag       ,unsigned  short const nPrcs){

  INT *xAdj=NULL,*adjncy=NULL,*RESTRICT perm=NULL;
  INT i,nDeg;

/*...*/  
  if(flag){
/*...*/
/*...................................................................*/

/*... armazena a malha no formato CSR*/
/*... calculo do vetor xAdj*/    
    HccaAlloc(INT,m,xAdj  ,(numelNov+1)      ,"xAdj"  ,_AD_);
/*...*/
    if(nPrcs > 1)
      convGraphPart(xAdj   ,adjncy
                   ,adj    ,nViz
                   ,maxViz,numelNov
                   ,true   ,false);
/*...................................................................*/
    
/*...*/
    else
      convGraph(xAdj ,adjncy
               ,adj   ,nViz
               ,maxViz,numelNov
               ,true  ,false);
/*...................................................................*/
    
/*... calculo do vetor adjncy*/    
    nDeg = xAdj[numelNov] -xAdj[0];
    HccaAlloc(INT,m,adjncy,nDeg,"adjncy",_AD_);
/*...*/
    if( nPrcs > 1)
      convGraphPart(xAdj   ,adjncy
                   ,adj    ,nViz
                   ,maxViz,numelNov
                   ,false  ,true);
/*...................................................................*/
    
/*...*/
    else
      convGraph(xAdj ,adjncy
               ,adj   ,nViz
               ,maxViz,numelNov
               ,false ,true);
/*...................................................................*/

/*... ordena o grafo CSR em ordem crescente*/
    sortGraphCsr(xAdj,adjncy,numelNov);
/*... soma + 1 em todas as posicoes dos vetores*/
    vectorPlusOne(xAdj,(numelNov+1),i);
    vectorPlusOne(adjncy,nDeg,i);  
/*...................................................................*/

/*...*/    
    HccaAlloc(INT,m,perm  ,numelNov        ,"perm"  ,_AD_);
/*...................................................................*/

/*...*/    
    genrcm (numelNov,xAdj,adjncy,perm);
/*...................................................................*/
 
    
    for(i=0;i<numelNov;i++){
      num[i] = perm[i];
    }
    
    for(;i<numel;i++){
      num[i] = i+1;
    }
    
/*...................................................................*/
   
    HccaDealloc(m,perm  ,"perm"  ,false);
    HccaDealloc(m,adjncy,"adjncy",false);
    HccaDealloc(m,xAdj  ,"xAdj"  ,false);
    
  } 
/*....................................................................*/  

/*... numeracao inicial*/  
  else
    for(i=0;i<numel;i++){
      num[i] = i+1;
    }
  
/*....................................................................*/  
//  for(i=0;i<numel;i++)
//    printf("%3d %3d\n",i+1,num[i]); 

}
/*********************************************************************/ 

