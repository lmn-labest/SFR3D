#ifndef _SOLV_H
  #define _SOLV_H
/*...*/
  #include<stdlib.h>
  #include<math.h>
  #include<string.h>
/*...................................................................*/

/*...*/
  #include<PreCond.h>
  #include<HccaBlas.h>
  #include<Sisteq.h>
  #include<Define.h>
  #include<Memoria.h>
  #include<Mesh.h>
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

/*... funcao de apoio*/
  void setSolver(char * word,short *solver);
  DOUBLE smachn();
/*...................................................................*/

/*....*/
  void solverC(Memoria *m     
              ,INT const nEq       ,INT const nEqNov
              ,INT const nAd       ,INT const nAdR
              ,INT *ia             ,INT *ja   
              ,DOUBLE *al          ,DOUBLE *ad        ,DOUBLE *au
              ,DOUBLE *b           ,DOUBLE *x
              ,Interface *iNeq
              ,DOUBLE const tol    ,unsigned int maxit
              ,short const storage ,short const solver
              ,FILE* fileSolvLog   ,bool const fLog
              ,bool const newX     ,bool const openmp
              ,bool const unsym    ,bool const loopwise);
/*...................................................................*/

/*========================= Iterativos ==============================*/
/*... gradiente conjugado precondicionado*/
  void pcg(INT const nEq        
          ,INT const nad  
          ,INT *restrict ia   ,INT *restrict ja
          ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
          ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict x
          ,DOUBLE *restrict z ,DOUBLE *restrict r ,DOUBLE const tol
          ,unsigned int maxit ,bool const newX          
          ,FILE* fileSolvLog  ,bool const log
          ,bool const fPrint
          ,void(*matvec)()    ,DOUBLE(*dot)());
  
  void mpiPcg(INT const nEq  ,INT const nEqNov  
          ,INT const nAd      ,INT const nAdr 
          ,INT *restrict ia   ,INT *restrict ja
          ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
          ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict x
          ,DOUBLE *restrict z ,DOUBLE *restrict r ,DOUBLE const tol
          ,unsigned int maxit ,bool const newX          
          ,FILE* fileSolvLog  ,bool const log
          ,bool const fPrint
          ,Interface *iNeq
          ,void(*matvec)()    ,DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado*/
  void pbicgstab(INT const nEq    ,INT const nAd
       ,INT *restrict ia   ,INT *restrict ja
       ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
       ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict x
       ,DOUBLE *restrict t ,DOUBLE *restrict v ,DOUBLE *restrict r
       ,DOUBLE *restrict p ,DOUBLE *restrict z ,DOUBLE const tol
       ,unsigned int maxIt ,bool const newX          
       ,FILE* fLog         ,bool const log
       ,bool const fPrint
       ,void(*matvec)()    ,DOUBLE(*dot)());
 
  void mpiPbicgstab(INT const nEq,INT const nEqNov      
      ,INT const nAd         ,INT const nAdR
      ,INT *restrict ia      ,INT *restrict ja
      ,DOUBLE *restrict al   ,DOUBLE *restrict ad,DOUBLE *restrict au
      ,DOUBLE *restrict m    ,DOUBLE *restrict b ,DOUBLE *restrict x
      ,DOUBLE *restrict t    ,DOUBLE *restrict v ,DOUBLE *restrict r
      ,DOUBLE *restrict p    ,DOUBLE *restrict z ,DOUBLE const tol
      ,unsigned int maxIt    ,bool const newX          
      ,FILE* fLog            ,bool const log
      ,bool const fPrint
      ,Interface *iNeq    
      ,void(*matvec)()     ,DOUBLE(*dot)());
/*...................................................................*/

/*===================================================================*/

#endif/*_SOLV_H*/
