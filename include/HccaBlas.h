#ifndef _HCCABLAS_H
  #define _HCCABLAS_H

/*...*/
  #include<Mesh.h>
  #include<ParallelMpi.h>
  #include<Erro.h>  
  #include<HccaStdBool.h>
  #include<Define.h>
  #include<HccaTime.h>
  #include<OpenMp.h>
/*...................................................................*/

/*...*/
  long flopDot(INT nDim);
  void prodVet(DOUBLE *restrict a
              ,DOUBLE *restrict b
              ,DOUBLE *restrict c);
/*...................................................................*/


/*======================== level 1 =================================*/
/*... soma de vetores*/
  void addVector(DOUBLE const alpha,DOUBLE *restrict a
                ,DOUBLE const beta ,DOUBLE *restrict b
                ,INT const nDim    ,DOUBLE *restrict c);
/*... multiplicacao de vetor por um escalar*/  
  void  alphaProdVector(DOUBLE const alpha,DOUBLE *restrict a
                       ,INT const nDim    ,DOUBLE *restrict c);
/*...................................................................*/

/*... produto enterno (MPI)*/
  DOUBLE dot(DOUBLE *restrict a,DOUBLE *restrict b,INT const nDim);
  DOUBLE dotO2I2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotI2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotI4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotI6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotI8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotO2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotO4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotO6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotO8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
/*...................................................................*/

/*... produto enterno (OPENMP)*/
  DOUBLE dotOmp(DOUBLE *restrict a,DOUBLE *restrict b,INT const nDim);
  DOUBLE dotOmpO2I2(DOUBLE *restrict x1,DOUBLE *restrict x2
                   ,INT const n);
  DOUBLE dotOmpO2I4(DOUBLE *restrict x1,DOUBLE *restrict x2
                   ,INT const n);
  DOUBLE dotOmpI2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotOmpI4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotOmpI6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotOmpI8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotOmpO2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotOmpO4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotOmpO6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
  DOUBLE dotOmpO8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n);
/*...................................................................*/
/*==================================================================*/

/*======================== level 2 =================================*/
/*... CsrDSym*/
  void matVecCsrDSym(INT const neq           
                    ,INT *restrict ia   ,INT *restrict ja
                    ,DOUBLE *restrict al,DOUBLE *restrict ad
                    ,DOUBLE *restrict x ,DOUBLE *restrict y);
/*...................................................................*/

/*... CsrDSymOmp*/
  void matVecCsrDsymOmp(INT const nEq
                       ,INT *restrict ia, INT *restrict ja
                       ,DOUBLE *restrict al, DOUBLE *restrict ad
                       ,DOUBLE *restrict x, DOUBLE *restrict y
                       ,INT  *restrict thBegin, INT *restrict thEnd
                       ,INT  *restrict thHeight
                       ,DOUBLE *restrict thY, int const nThreads);
/*...................................................................*/

/*... CsrD*/ 
  void     matVecCsrD(INT const neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDI2(INT const neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDI4(INT const neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDI6(INT const neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDO2(INT const neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDO4(INT const neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void   matVecCsrDO6(INT const neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
  void matVecCsrDO2I2(INT const neq           
                     ,INT *restrict ia  ,INT *restrict ja
                     ,DOUBLE *restrict a,DOUBLE *restrict ad
                     ,DOUBLE *restrict x,DOUBLE *restrict y);
/*...................................................................*/

/*... OPENMP - CsrD*/
  void matVecCsrDomp(INT const neq
                    ,INT *restrict ia      ,INT *restrict ja
                    ,DOUBLE *restrict a    ,DOUBLE *restrict ad
                    ,DOUBLE *restrict x    ,DOUBLE *restrict y
                    ,INT  *restrict thBegin,INT *restrict thEnd
                    ,INT  *restrict thHeight
                    ,DOUBLE *ddum          ,int const nThreads);
/*...................................................................*/

/*... MPI - CsrD */ 
/*... Simetrico*/
/*... CSRD+CSR*/
  void mpiMatVecCsrDSym(INT const nEq    ,INT const *nAd 
                     ,INT *restrict ia   ,INT *restrict ja
                     ,DOUBLE *restrict al,DOUBLE *restrict ad
                     ,DOUBLE *restrict x ,DOUBLE *restrict y
                     ,Interface *iNeq);
/*... CSRD+COO*/
  void mpiMatVecCsrDcooSym(INT const nEq      ,INT const *nAd       
                          ,INT *restrict ia   ,INT *restrict ja
                          ,DOUBLE *restrict al,DOUBLE *restrict ad
                          ,DOUBLE *restrict x ,DOUBLE *restrict y
                          ,Interface *iNeq);
/*... Geral*/
  void mpiMatVecCsrD(INT const neq   ,INT const *nAd 
                  ,INT *restrict ia  ,INT *restrict ja
                  ,DOUBLE *restrict a,DOUBLE *restrict ad
                  ,DOUBLE *restrict x,DOUBLE *restrict y
                  ,Interface *iNeq);
  
  
/*...................................................................*/

/*... CsrC */ 
/*... Simetrico*/
/*... ESTRUTURALMENTE SIMETRICA*/
  void matVecCsrC(INT const neq                       
                 ,INT *restrict ia  ,INT *restrict ja
                 ,DOUBLE *restrict a,DOUBLE *restrict ad
                 ,DOUBLE *restrict x,DOUBLE *restrict y);

/*... CsrCsymOmp*/
  void matVecCsrComp(INT const nEq
                    ,INT *restrict ia      ,INT *restrict ja
                    ,DOUBLE *restrict a    ,DOUBLE *restrict ad
                    ,DOUBLE *restrict x    ,DOUBLE *restrict y
                    ,INT  *restrict thBegin,INT *restrict thEnd
                    ,INT  *restrict thHeight
                    ,DOUBLE *restrict thY, int const nThreads);
/*...................................................................*/

/*... CSRD+CSR*/
  void mpiMatVecCsrC(INT const nEq      ,INT const *nAd      
                    ,INT *restrict ia   ,INT *restrict ja
                    ,DOUBLE *restrict a ,DOUBLE *restrict ad
                    ,DOUBLE *restrict x ,DOUBLE *restrict y
                    ,Interface *iNeq);

/*... CSRC+COO*/
  void mpiMatVecCsrCcoo(INT const nEq      ,INT const *nAd      
                       ,INT *restrict ia   ,INT *restrict ja
                       ,DOUBLE *restrict a ,DOUBLE *restrict ad
                       ,DOUBLE *restrict x ,DOUBLE *restrict y
                       ,Interface *iNeq);
  
/*...................................................................*/


/*... EllPack*/ 
  void matVecEllPack(INT const nEq           
                  ,INT *restrict ia  ,INT *restrict ja
                  ,DOUBLE *restrict a,DOUBLE *restrict ad
                  ,DOUBLE *restrict x,DOUBLE *restrict y);
  void matVecEllPackO2(INT const nEq           
                      ,INT *restrict ia  ,INT *restrict ja
                      ,DOUBLE *restrict a,DOUBLE *restrict ad
                      ,DOUBLE *restrict x,DOUBLE *restrict y);
  void matVecEllPackO4(INT const nEq           
                      ,INT *restrict ia  ,INT *restrict ja
                      ,DOUBLE *restrict a,DOUBLE *restrict ad
                      ,DOUBLE *restrict x,DOUBLE *restrict y);
/*... MPI*/
  void mpiMatVecEllPack(INT const nEq  ,INT const *nAd 
                      ,INT *restrict ia  ,INT *restrict ja
                      ,DOUBLE *restrict a,DOUBLE *restrict ad
                      ,DOUBLE *restrict x,DOUBLE *restrict y
                      ,Interface *iNeq);
/*...................................................................*/


/*==================================================================*/

/*======================== level 3 =================================*/
/*==================================================================*/

#endif
