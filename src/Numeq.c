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
INT numeq(INT *RESTRICT id       ,INT *RESTRICT num
         ,short *RESTRICT rt     ,short *RESTRICT nFace 
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
      if( MAT2D(nel,aux,rt,maxRes) == PCCELL)
        MAT2D(nel,j,id,ndf) = -1; 
      else
        MAT2D(nel,j,id,ndf) = ++neq;
    } 
  }
/*...................................................................*/

  return neq;
}
/*********************************************************************/ 

/*********************************************************************
 * numEqV1: numeracao das equacoes para as velocidade descopladas    *
 * do metodo simple                                                  *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * m     -> memoria principal                                        * 
 * id    -> indefinido                                               *
 * num   -> renumeracao dos elementos                                *
 * numel -> numero de elementos                                      * 
 * str   -> nome do vetor                                            * 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * id    -> numeracao das equacoes por celula                        *
 * ------------------------------------------------------------------*
 * OBS: a numera das equações u, v, e w são iguais e todas as celulas*
 * tem sempre uma equacao.                                           *
 * *******************************************************************/
INT numEqV1(INT *RESTRICT id       ,INT *RESTRICT num
           ,INT const numel        ){

  INT i,neq,nel;

/*...*/  
  neq = 0;
  for(i=0;i<numel;i++){
    nel = num[i] -1;
    id[nel] = ++neq;
  }
/*...................................................................*/

  return neq;
}
/*********************************************************************/ 

/*********************************************************************
 * numEqV2: numeracao das equacoes das pressoes do metodo simple     *
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
 * id    -> numeracao das equacoes por celula                        *
 * ------------------------------------------------------------------*
 * *******************************************************************/
INT numEqV2(INT *RESTRICT id       ,INT *RESTRICT num
           ,short *RESTRICT rt     ,short *RESTRICT nFace 
           ,INT const numel        ,short const maxViz){

  INT i,neq,nel;
  short maxRes = (maxViz + 1);
  short aux;

/*...*/  
  neq = 0;
  for(i=0;i<numel;i++){
    nel = num[i] -1;
    aux = nFace[nel];
    if( MAT2D(nel,aux,rt,maxRes) == PCCELL)
      id[nel] = -1; 
    else
      id[nel] = ++neq;
  }
/*...................................................................*/

  return neq;
}
/*********************************************************************/ 

/*********************************************************************
 * NUMEQNOV : conta as equacoes sem sobre posicao                    *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * m     -> memoria principal                                        * 
 * num   -> renumeracao dos elementos                                *
 * rt    -> restricoes                                               * 
 * nFace -> numero de faces por elemento                             * 
 * numel -> numero de elementos                                      * 
 * maxViz-> numero maximo de vizinhos                                * 
 * str   -> nome do vetor                                            * 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * *******************************************************************/
INT countEq(INT *RESTRICT num
           ,short *RESTRICT rt     ,short *RESTRICT nFace 
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
      if( MAT2D(nel,aux,rt,maxRes) != 1)
        ++neq;
    } 
  }
/*...................................................................*/

  return neq;
}
/*********************************************************************/ 

/*********************************************************************
 * FRONT : mapa de equacoes na interface de comunicao                *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * m     -> memoria principal                                        * 
 * num   -> renumeracao dos elementos                                *
 * rt    -> restricoes                                               * 
 * nFace -> numero de faces por elemento                             * 
 * numel -> numero de elementos                                      * 
 * maxViz-> numero maximo de vizinhos                                * 
 * str   -> nome do vetor                                            * 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * *******************************************************************/
void front(Memoria *m
          ,PartMesh *pMesh, SistEq *sistEq, short const ndf){

#ifdef _MPICH_
  INT sNeq = 0;
  INT rNeq = 0;
  INT k,kk,j,nel,nEq1,aux=0;
  short jj;
  unsigned short nPart,nVizParts = pMesh->iEl.nVizPart; 
  

/*... buffer de recebimentos*/
  for(nPart = 0; nPart < nVizParts; nPart++){
    k      = pMesh->iEl.iaRcvs[nPart];
    kk     = pMesh->iEl.iaRcvs[nPart+1];
    for(j=k;j<kk;j++){
      nel  =  pMesh->iEl.fMap[j];
      for(jj=0;jj<ndf;jj++){
        nEq1 =  MAT2D(nel,jj,sistEq->id,ndf); 
        if( nEq1 != -1 )
          rNeq++;
      } 
    }
  }
/*...................................................................*/

/*... buffer de envios*/
  for(nPart = 0; nPart < nVizParts; nPart++){
    k      = pMesh->iEl.iaSends[nPart];
    kk     = pMesh->iEl.iaSends[nPart+1];
    for(j=k;j<kk;j++){
      nel  =  pMesh->iEl.fMap[j];
      for(jj=0;jj<ndf;jj++){
        nEq1 =  MAT2D(nel,jj,sistEq->id,ndf); 
        if( nEq1 != -1 )
          sNeq++;
      } 
    }
  }
/*...................................................................*/

/*... alocando memoria*/

  kk = rNeq+sNeq;
  sistEq->iNeq.nRcvs  = rNeq; 
  sistEq->iNeq.nSends = sNeq;
   
/*... fMap(buffer de comunicacao)*/ 
  HccaAlloc(INT      ,m,sistEq->iNeq.fMap,kk   
           ,"fMapNeq",false);
  zero(sistEq->iNeq.fMap,kk,INTC);

/*... alocando memoria*/
  kk = nVizParts + 1;
  
  HccaAlloc(INT         ,m    ,sistEq->iNeq.iaSends , kk        
           ,"iaSendsNeq",false);
  zero(sistEq->iNeq.iaSends, kk,INTC);
        
  HccaAlloc(INT        ,m   ,sistEq->iNeq.iaRcvs , kk       
           ,"iaRcvsNeq",false);
  zero(sistEq->iNeq.iaRcvs, kk,INTC);
/*...................................................................*/


/*... buffer de comunicacao*/
  kk = rNeq+sNeq;
  HccaAlloc(DOUBLE,m,sistEq->iNeq.xb 
           ,kk    ,"xBufferMpiNeq"              ,false);
  zero(sistEq->iNeq.xb, kk,DOUBLEC);
/*...................................................................*/


/*... buffer de recebimentos*/
  for(nPart = 0; nPart < nVizParts; nPart++){
    k      = pMesh->iEl.iaRcvs[nPart];
    kk     = pMesh->iEl.iaRcvs[nPart+1];
    sistEq->iNeq.iaRcvs[nPart] = aux;
    for(j=k;j<kk;j++){
      nel  =  pMesh->iEl.fMap[j];
      for(jj=0;jj<ndf;jj++){
        nEq1 =  MAT2D(nel,jj,sistEq->id,ndf); 
        if( nEq1 != -1 ){
          sistEq->iNeq.fMap[aux++] = nEq1-1;
        }
      } 
    }
  }
  sistEq->iNeq.iaRcvs[nVizParts] = aux;
/*...................................................................*/

/*... buffer de envios*/
  for(nPart = 0; nPart < nVizParts; nPart++){
    k      = pMesh->iEl.iaSends[nPart];
    kk     = pMesh->iEl.iaSends[nPart+1];
    sistEq->iNeq.iaSends[nPart] = aux;
    for(j=k;j<kk;j++){
      nel  =  pMesh->iEl.fMap[j];
      for(jj=0;jj<ndf;jj++){
        nEq1 =  MAT2D(nel,jj,sistEq->id,ndf); 
        if( nEq1 != -1 ){
          sistEq->iNeq.fMap[aux++] = nEq1-1;
        }
      } 
    }
  }
  sistEq->iNeq.iaSends[nVizParts] = aux;
/*...................................................................*/

/*...*/
  sistEq->iNeq.nVizPart = nVizParts;
  HccaAlloc(short,m,sistEq->iNeq.vizPart , nVizParts   
           ,"nVizPartNeq"                ,false);
  zero(sistEq->iNeq.vizPart              ,nVizParts,"short");
  for(nPart = 0; nPart < nVizParts; nPart++){
    sistEq->iNeq.vizPart[nPart] = pMesh->iEl.vizPart[nPart];
  }

/*...................................................................*/

/*  
  printf("r %d s %d\n",sistEq->iNeq.nRcvs,sistEq->iNeq.nSends);
   
  for(nPart = 0; nPart < nVizParts; nPart++){
    printf("viz %d\n",sistEq->iNeq.vizPart[nPart]);
  }
 
  for(nPart = 0; nPart < nVizParts; nPart++){
    printf("iaR %d\n",sistEq->iNeq.iaRcvs[nPart]);
    printf("iaS %d\n",sistEq->iNeq.iaSends[nPart]);
  }
  
  
  for(nPart = 0; nPart <  rNeq+sNeq; nPart++)
    printf("fMap %d\n",sistEq->iNeq.fMap[nPart]+1);
*/  
#endif
}  
/*********************************************************************/ 
