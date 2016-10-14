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
  #include<OpenMp.h>
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
  void setSolverConfig(char *word,Solv *solv,FILE *fileIn);
  DOUBLE smachn();
	void setDot(DOUBLE(**dotC)(), short const iCod);
	void setMatVec(void(**matVecC)(),short const storage
                ,bool const unSym ,bool const openMp);
/*...................................................................*/

/*....*/
  void solverC(Memoria *m     
              ,INT const nEq       ,INT const nEqNov
              ,INT const nAd       ,INT const nAdR
              ,INT *ia             ,INT *ja   
              ,DOUBLE *al          ,DOUBLE *ad        ,DOUBLE *au
              ,DOUBLE *b           ,DOUBLE *x
              ,Interface *iNeq     ,BufferOmp *bOmp
              ,DOUBLE const tol    ,unsigned int maxit
              ,short const storage ,short const solver
              ,FILE* fileSolvLog   ,bool const fLog
              ,bool const newX     ,bool const unsym    );
/*...................................................................*/

/*========================= Iterativos ==============================*/
  void callCg(INT const nEq      ,INT const nEqNov
              ,INT const nAd      ,INT const nAdr
              ,INT *restrict ia   ,INT *restrict ja
              ,DOUBLE *restrict al,DOUBLE *restrict ad
              ,DOUBLE *restrict m ,DOUBLE *restrict b
              ,DOUBLE *restrict x ,DOUBLE *restrict z 
              ,DOUBLE *restrict r ,DOUBLE *restrict p
              ,DOUBLE const tol   ,unsigned int maxit 
              ,bool const newX    ,FILE* fileSolvLog  
              ,bool const log     ,bool const fPrint  
              ,Interface *iNeq    ,BufferOmp *bOmp
              ,void(*matvec)()    ,DOUBLE(*dot)());

/*... gradiente conjugado precondicionado*/
	void pcg(INT const nEq,INT const nAd
		,INT *restrict ia   ,INT *restrict ja
		,DOUBLE *restrict al,DOUBLE *restrict ad
		,DOUBLE *restrict m ,DOUBLE *restrict b
    ,DOUBLE *restrict x ,DOUBLE *restrict z
    ,DOUBLE *restrict r ,DOUBLE *restrict p
		,DOUBLE const tol   ,unsigned int maxIt
    ,bool const newX  	,FILE* fLog
    ,FILE *fileHistLog  ,bool const log 
    ,bool const fHistLog,bool const fPrint
		,void(*matvec)()    ,DOUBLE(*dot)());

/*...  gradiente conjugado precondicionado (OPENMP)*/  
  void pcgOmp(INT const nEq      ,INT const nAd
             ,INT *restrict ia   ,INT *restrict ja
             ,DOUBLE *restrict a ,DOUBLE *restrict ad
             ,DOUBLE *restrict m ,DOUBLE *restrict b
             ,DOUBLE *restrict x ,DOUBLE *restrict z
             ,DOUBLE *restrict r ,DOUBLE *restrict p
             ,DOUBLE const tol   ,unsigned int maxIt
             ,bool const newX    ,FILE* fLog
             ,FILE *fileHistLog  ,bool const log
             ,bool const fHistLog, bool const fPrint
             ,BufferOmp *bOmp
             ,void(*matvec)(), DOUBLE(*dot)());

/*...  gradiente conjugado precondicionado (MPI)*/
  void mpiPcg(INT const nEq  ,INT const nEqNov  
          ,INT const nAd      ,INT const nAdr 
          ,INT *restrict ia   ,INT *restrict ja
          ,DOUBLE *restrict al,DOUBLE *restrict ad
          ,DOUBLE *restrict m ,DOUBLE *restrict b 
          ,DOUBLE *restrict x ,DOUBLE *restrict z
          ,DOUBLE *restrict r 
          ,DOUBLE const tol   ,unsigned int maxit 
          ,bool const newX    ,FILE* fileSolvLog  
          ,bool const log     ,bool const fPrint
          ,Interface *iNeq
          ,void(*matvec)()    ,DOUBLE(*dot)());
/*...................................................................*/

/*...*/
  void callBicgStab(INT const nEq     ,INT const nEqNov
                   ,INT const nAd     ,INT const nAdR
                   ,INT *restrict ia  ,INT *restrict ja
                   ,DOUBLE *restrict a,DOUBLE *restrict ad
                   ,DOUBLE *restrict m,DOUBLE *restrict b
                   ,DOUBLE *restrict x,DOUBLE *restrict t
                   ,DOUBLE *restrict v,DOUBLE *restrict r
                   ,DOUBLE *restrict p,DOUBLE *restrict z
                   ,DOUBLE *restrict h
                   ,DOUBLE const tol, unsigned int maxIt
                   ,bool const newX, FILE* fSolvLog
                   ,bool const fLog, bool const fPrint
                   ,Interface *iNeq, BufferOmp *bOmp
                   ,void(*matVec)(), DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado*/
  void pbicgstab(INT const nEq    ,INT const nAd
       ,INT *restrict ia   ,INT *restrict ja
       ,DOUBLE *restrict al,DOUBLE *restrict ad
       ,DOUBLE *restrict m ,DOUBLE *restrict b 
       ,DOUBLE *restrict x ,DOUBLE *restrict t 
       ,DOUBLE *restrict v ,DOUBLE *restrict r
       ,DOUBLE *restrict p ,DOUBLE *restrict z 
       ,DOUBLE *restrict r0,DOUBLE const tol
       ,unsigned int maxIt ,bool const newX          
	     ,FILE* fLog         ,FILE *fileHistLog
	     ,bool const log     , bool const fHistLog
       ,bool const fPrint
       ,void(*matvec)()    ,DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado (OMP)*/
  void pbicgstabOmp(INT const nEq, INT const nAd
                   ,INT *restrict ia, INT *restrict ja
                   ,DOUBLE *restrict a, DOUBLE *restrict ad
                   ,DOUBLE *restrict m, DOUBLE *restrict b
                   ,DOUBLE *restrict x, DOUBLE *restrict t
                   ,DOUBLE *restrict v, DOUBLE *restrict r
                   ,DOUBLE *restrict p, DOUBLE *restrict z
                   ,DOUBLE *restrict r0
                   ,DOUBLE const tol, unsigned int maxIt
                   ,bool const newX, FILE* fLog
                   ,FILE *fileHistLog, bool const log
                   ,bool const fHistLog, bool const fPrint
                   ,BufferOmp *bOmp
                   ,void(*matvec)(), DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado (MPI)*/
  void mpiPbicgstab(INT const nEq,INT const nEqNov      
      ,INT const nAd         ,INT const nAdR
      ,INT *restrict ia      ,INT *restrict ja
      ,DOUBLE *restrict al   ,DOUBLE *restrict ad
      ,DOUBLE *restrict m    ,DOUBLE *restrict b 
      ,DOUBLE *restrict x    ,DOUBLE *restrict t 
      ,DOUBLE *restrict v    ,DOUBLE *restrict r
      ,DOUBLE *restrict p    ,DOUBLE *restrict z 
      ,DOUBLE const tol
      ,unsigned int maxIt    ,bool const newX          
      ,FILE* fLog            ,bool const log
      ,bool const fPrint     ,Interface *iNeq    
      ,void(*matvec)()       ,DOUBLE(*dot)());
/*...................................................................*/

/*...*/
  void callBicgStabl2(INT const nEq, INT const nEqNov
                      , INT const nAd, INT const nAdR
                      , INT *restrict ia, INT *restrict ja
                      , DOUBLE *restrict a, DOUBLE *restrict ad
                      , DOUBLE *restrict m, DOUBLE *restrict b
                      , DOUBLE *restrict x, DOUBLE *restrict t
                      , DOUBLE *restrict v, DOUBLE *restrict r
                      , DOUBLE *restrict u, DOUBLE *restrict r0
                      , DOUBLE *restrict w, DOUBLE *restrict s
                      , DOUBLE *restrict p, DOUBLE *restrict h
                      , DOUBLE *restrict z
                      , DOUBLE const tol, unsigned int maxIt
                      , bool const newX, FILE* fSolvLog
                      , bool const fLog, bool const fPrint
                      , Interface *iNeq, BufferOmp *bOmp
                      , void(*matVec)(), DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado l=2*/
  void pbicgstabl2(INT const nEq, INT const nAd
                   , INT *restrict ia, INT *restrict ja
                   , DOUBLE *restrict a, DOUBLE *restrict ad
                   , DOUBLE *restrict m, DOUBLE *restrict b
                   , DOUBLE *restrict x, DOUBLE *restrict t
                   , DOUBLE *restrict v, DOUBLE *restrict r
                   , DOUBLE *restrict u, DOUBLE *restrict r0
                   , DOUBLE *restrict w, DOUBLE *restrict s
                   , DOUBLE *restrict p, DOUBLE *restrict h
                   , DOUBLE *restrict z
                   , DOUBLE const tol, unsigned int maxIt
                   , bool const newX, FILE* fLog
                   , FILE *fileHistLog, bool const log
                   , bool const fHistLog, bool const fPrint
                   , void(*matvec)(), DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado (OMP)*/
  void pbicgstabl2Omp(INT const nEq, INT const nAd
                      , INT *restrict ia, INT *restrict ja
                      , DOUBLE *restrict a, DOUBLE *restrict ad
                      , DOUBLE *restrict m, DOUBLE *restrict b
                      , DOUBLE *restrict x, DOUBLE *restrict t
                      , DOUBLE *restrict v, DOUBLE *restrict r
                      , DOUBLE *restrict u, DOUBLE *restrict r0
                      , DOUBLE *restrict w, DOUBLE *restrict s
                      , DOUBLE *restrict p, DOUBLE *restrict h
                      , DOUBLE *restrict z
                      , DOUBLE const tol, unsigned int maxIt
                      , bool const newX, FILE* fLog
                      , FILE *fileHistLog, bool const log
                      , bool const fHistLog, bool const fPrint
                      , BufferOmp *bOmp
                      , void(*matvec)(), DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado (MPI)*/
/*void mpiPbicgstab(INT const nEq, INT const nEqNov
                    , INT const nAd, INT const nAdR
                    , INT *restrict ia, INT *restrict ja
                    , DOUBLE *restrict al, DOUBLE *restrict ad
                    , DOUBLE *restrict m, DOUBLE *restrict b
                    , DOUBLE *restrict x, DOUBLE *restrict t
                    , DOUBLE *restrict v, DOUBLE *restrict r
                    , DOUBLE *restrict p, DOUBLE *restrict z
                    , DOUBLE const tol
                    , unsigned int maxIt, bool const newX
                    , FILE* fLog, bool const log
                    , bool const fPrint, Interface *iNeq
                    , void(*matvec)(), DOUBLE(*dot)());*/
/*...................................................................*/

/*...*/
  void callGmres(INT const nEq, INT const nEqNov
                 , INT const nAd, INT const nAdR
                 , INT *restrict ia, INT *restrict ja
                 , DOUBLE *restrict a, DOUBLE *restrict ad
                 , DOUBLE *restrict m, DOUBLE *restrict b
                 , DOUBLE *restrict x, DOUBLE *restrict g
                 , DOUBLE *restrict h, DOUBLE *restrict y
                 , DOUBLE *restrict c, DOUBLE *restrict s
                 , DOUBLE *restrict e, short const nKrylov
                 , DOUBLE const tol  , unsigned int maxIt
                 , bool const newX   , FILE* fSolvLog
                 , bool const fLog   , bool const fPrint
                 , Interface *iNeq  , BufferOmp *bOmp
                 , void(*matVec)()  , DOUBLE(*dot)());
/*...................................................................*/

/*... GMRES*/
  void gmres(INT const nEq       ,INT const nAd
            ,INT *restrict ia    ,INT *restrict ja
            ,DOUBLE *restrict a  ,DOUBLE *restrict ad
            ,DOUBLE *restrict m  ,DOUBLE *restrict b
            ,DOUBLE *restrict x  ,DOUBLE *restrict g
            ,DOUBLE *restrict h  ,DOUBLE *restrict y
            ,DOUBLE *restrict c  ,DOUBLE *restrict s
            ,DOUBLE *restrict e
            ,short const nKrylov
            ,DOUBLE const tol    ,unsigned int nCycles
            ,bool const newX     ,FILE* fLog
            ,FILE *fileHistLog   ,bool const log
            ,bool const fHistLog ,bool const fPrint
            ,void(*matvec)()     ,DOUBLE(*dot)());
/*...................................................................*/

/*... Omp*/
  void gmresOmp(INT const nEq, INT const nAd
                , INT *restrict ia, INT *restrict ja
                , DOUBLE *restrict a, DOUBLE *restrict ad
                , DOUBLE *restrict m, DOUBLE *restrict b
                , DOUBLE *restrict x, DOUBLE *restrict g
                , DOUBLE *restrict h, DOUBLE *restrict y
                , DOUBLE *restrict c, DOUBLE *restrict s
                , DOUBLE *restrict e, short const nKrylov
                , DOUBLE const tol, unsigned int nCycles
                , bool const newX, FILE* fLog
                , FILE *fileHistLog, bool const log
                , bool const fHistLog, bool const fPrint
                , BufferOmp *bOmp
                , void(*matvec)(), DOUBLE(*dot)());
/*...................................................................*/
/*===================================================================*/

#endif/*_SOLV_H*/
