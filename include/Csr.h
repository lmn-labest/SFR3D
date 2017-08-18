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
  INT csrIa(INT *RESTRICT ia     ,INT *RESTRICT id 
           ,INT *RESTRICT num    ,INT *RESTRICT adj
           ,short *RESTRICT nViz
           ,INT const numel      ,INT const neq  
           ,short const maxViz   ,short const ndf 
           ,bool const upper     ,bool const diag 
           ,bool const lower     );

  INT csrIaR(INT *RESTRICT ia     ,INT *RESTRICT id   
            ,INT *RESTRICT num    ,INT *RESTRICT adj
            ,short *RESTRICT nViz
            ,INT const numel      ,INT const neq
            ,short const maxViz   ,short  const ndf);

  void csrJa(INT *RESTRICT ia   ,INT *RESTRICT ja 
           ,INT *RESTRICT id    ,INT *RESTRICT num
           ,INT *RESTRICT adj   ,short *RESTRICT nViz
           ,INT const numel     ,INT const neq 
           ,short const maxViz  ,short const ndf
           ,bool const upper    ,bool const diag
           ,bool const lower);
  
  void csrJaR(INT *RESTRICT ia   ,INT *RESTRICT ja 
             ,INT *RESTRICT id    ,INT *RESTRICT num
             ,INT *RESTRICT adj   ,short *RESTRICT nViz
             ,INT const numel     ,INT const neq 
             ,short const maxViz  ,short const ndf);
/*...................................................................*/

/*...*/  
  void csr(INT    *RESTRICT  ia,INT *RESTRICT ja 
          ,DOUBLE *RESTRICT a  ,DOUBLE *RESTRICT ad
          ,DOUBLE *RESTRICT b  ,INT *RESTRICT lId                       
          ,DOUBLE *RESTRICT lA ,DOUBLE *RESTRICT lB 
          ,INT const nEq       ,INT const nEqNov 
          ,INT const nAd       ,INT const nAdR 
          ,short const nFace   ,short const ndf  
          ,short const storade ,bool  const forces
          ,bool const matrix   ,bool  const  unsym);
/*...................................................................*/

/*...*/  
  void csrSimple(INT    *RESTRICT  ia,INT *RESTRICT ja 
          ,DOUBLE *RESTRICT a  ,DOUBLE *RESTRICT ad
          ,DOUBLE *RESTRICT b  ,INT *RESTRICT lId                       
          ,DOUBLE *RESTRICT lA ,DOUBLE *RESTRICT lB 
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
  void sortGraphCsr(INT *RESTRICT ia,INT *RESTRICT ja,INT const n);
/*...................................................................*/

/*... divisao da matriz para o openmp*/
  void partitionCsrByNonzeros(INT *RESTRICT ia      ,INT *RESTRICT ja
                 ,INT const nEq
                 ,int nThreads          ,INT *RESTRICT thBegin
                 ,INT *RESTRICT thEnd   ,INT *RESTRICT thSize
                 ,INT *RESTRICT thHeight,short type);
  void computeEffectiveWork(INT *RESTRICT ia        ,INT *RESTRICT ja
                 ,INT const nEq
                 ,INT *RESTRICT thBegin,INT *RESTRICT thEnd
                 ,INT *RESTRICT thSize ,INT *RESTRICT thHeight);
/*...................................................................*/

#endif/*_CSR_H*/
