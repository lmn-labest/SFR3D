#ifndef _CSR_H  
  #define _CSR_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
/*...................................................................*/

/*...*/
  #include<Define.h>
  #include<HccaStdBool.h>
  #include<ParallelMpi.h>
/*...................................................................*/

/*...*/
  INT csrIa(INT *restrict ia     ,INT *restrict id 
           ,INT *restrict num    ,INT *restrict adj
           ,short *restrict nViz
           ,INT const numel      ,INT const neq  
           ,short const maxViz   ,short const ndf 
           ,bool const upper     ,bool const diag 
           ,bool const lower     );

  INT csrIaR(INT *restrict ia     ,INT *restrict id   
            ,INT *restrict num    ,INT *restrict adj
            ,short *restrict nViz
            ,INT const numel      ,INT const neq
            ,short const maxViz   ,short  const ndf);

  void csrJa(INT *restrict ia   ,INT *restrict ja 
           ,INT *restrict id    ,INT *restrict num
           ,INT *restrict adj   ,short *restrict nViz
           ,INT const numel     ,INT const neq 
           ,short const maxViz  ,short const ndf
           ,bool const upper    ,bool const diag
           ,bool const lower);
  
  void csrJaR(INT *restrict ia   ,INT *restrict ja 
             ,INT *restrict id    ,INT *restrict num
             ,INT *restrict adj   ,short *restrict nViz
             ,INT const numel     ,INT const neq 
             ,short const maxViz  ,short const ndf);
/*...................................................................*/

/*...*/  
  void csr(INT    *restrict  ia,INT *restrict ja 
          ,DOUBLE *restrict a  ,DOUBLE *restrict ad
          ,DOUBLE *restrict b  ,INT *restrict lId                       
          ,DOUBLE *restrict lA ,DOUBLE *restrict lB 
          ,INT const nEq       ,INT const nEqNov 
          ,INT const nAd       ,INT const nAdR 
          ,short const nFace   ,short const ndf  
          ,short const storade ,bool  const forces
          ,bool const matrix   ,bool  const  unsym);
/*...................................................................*/

/*...*/  
  void csrSimple(INT    *restrict  ia,INT *restrict ja 
          ,DOUBLE *restrict a  ,DOUBLE *restrict ad
          ,DOUBLE *restrict b  ,INT *restrict lId                       
          ,DOUBLE *restrict lA ,DOUBLE *restrict lB 
          ,INT const nEq       ,INT const nEqNov 
          ,INT const nAd       ,INT const nAdR 
          ,short const nFace   ,short const ndf  
          ,short const storade ,bool  const forces
          ,bool const matrix   ,bool  const  unsym);
/*...................................................................*/


/*...*/
  INT bandCsr(INT *ia,INT *ja,INT  neq,short type);
  INT bandCsrC(INT *ia,INT *ja,INT  neq,short type);
/*...................................................................*/

/*...*/
  void sortGraphCsr(INT *restrict ia,INT *restrict ja,INT const n);
/*...................................................................*/

/*... divisao da matriz para o openmp*/
  void partitionCsrByNonzeros(INT *restrict ia      ,INT *restrictja
                 ,INT const nEq
                 ,int nThreads          ,INT *restrict thBegin
                 ,INT *restrict thEnd   ,INT *restrict thSize
                 ,INT *restrict thHeight,short type);
  void computeEffectiveWork(INT *restrict ia        ,INT *restrict ja
                 ,INT const nEq
                 ,INT *restrict thBegin,INT *restrict thEnd
                 ,INT *restrict thSize ,INT *restrict thHeight);
/*...................................................................*/

#endif/*_CSR_H*/
