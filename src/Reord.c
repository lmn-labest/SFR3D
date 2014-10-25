#include<Reord.h>
/********************************************************************* 
 * REORD : reordanacao do grapho dos elementos                       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * num   -> indefinido                                               * 
 * graph  -> grafo dos elementos                                     * 
 * numel -> numero de elementos                                      * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * num   -> elementos renumerados                                    * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *********************************************************************/
void reord(Memoria *m,INT *num,INT const *adj ,short const*nViz
          ,short const maxViz  ,INT numel, bool flag){

  INT *xAdj=NULL,*adjncy=NULL,*restrict perm=NULL;
  INT i,nDeg;

/*...*/  
  if(flag){
/*...*/
/*...................................................................*/

/*... armazena a malha no formato CSR*/
/*... calculo do vetor xAdj*/    
    HccaAlloc(INT,m,xAdj  ,(numel+1)      ,"xAdj"  ,_AD_);
    convGraph(xAdj,adjncy,adj,nViz,maxViz,numel,true,false);
/*... calculo do vetor adjncy*/    
    nDeg = xAdj[numel] -xAdj[0];
    HccaAlloc(INT,m,adjncy,nDeg,"adjncy",_AD_);
    convGraph(xAdj,adjncy,adj,nViz,maxViz,numel,false,true);
/*... ordena o grafo CSR em ordem crescente*/
    sortGraphCsr(xAdj,adjncy,numel);
/*... soma + 1 em todas as posicoes dos vetores*/
     vectorPlusOne(xAdj,(numel+1),i);
     vectorPlusOne(adjncy,nDeg,i);  
/*...................................................................*/

/*...*/    
    HccaAlloc(INT,m,perm  ,numel        ,"perm"  ,_AD_);
/*...................................................................*/

/*...*/    
    genrcm (numel,xAdj,adjncy,perm);
/*...................................................................*/
 
    
    for(i=0;i<numel;i++){
      num[i] = perm[i];
    }
    
/*...................................................................*/
//  for(i=0;i<numel;i++)
//    printf("%3ld %3ld %3ld\n",i+1,perm[i],num[i]); 
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

}
/*********************************************************************/ 

