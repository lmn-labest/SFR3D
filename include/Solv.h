#ifndef _SOLV_H
  #define _SOLV_H
/*...*/
  #include<stdlib.h>
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<preCond.h>
  #include<HccaBlas.h>
  #include<Define.h>
  #include<Memoria.h>
/*...................................................................*/
  
/*....*/
  typedef struct{
    DOUBLE            tol;
    unsigned int    maxIt;
    short          solver;
    FILE        *fileSolv;
    bool              log;
    bool             flag;
  }Solv;
/*...................................................................*/

/*....*/
  void solverC(Memoria *m          ,INT const neq     ,INT const nad
              ,INT *ia             ,INT *ja   
              ,DOUBLE *al          ,DOUBLE *ad        ,DOUBLE *au
              ,DOUBLE *b           ,DOUBLE *x
              ,DOUBLE tol          ,unsigned int maxit
              ,short const storage ,short const solver
              ,FILE* fileSolvLog   ,bool const fLog
              ,bool const newX     ,bool const openmp  
              ,bool const unsym    ,bool const loopwise);
/*...................................................................*/

/*========================= Iterativos ==============================*/
/*... gradiente conjugado precondicionado*/
  void pcg(INT const neq      ,INT const nad  
          ,INT *restrict ia   ,INT *restrict ja
          ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
          ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict x
          ,DOUBLE *restrict z ,DOUBLE *restrict r ,DOUBLE const tol
          ,unsigned int maxit ,bool const newX          
          ,FILE* fileSolvLog  ,bool const log
          ,bool const fPrint
          ,void(*matvec)()    ,DOUBLE(*dot)());
/*...................................................................*/
/*===================================================================*/

#endif/*_SOLV_H*/
