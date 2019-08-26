#ifndef _SOLV_H_
  #define _SOLV_H_
/*...*/
  #include<stdlib.h>
  #include<math.h>
  #include<string.h>
/*...................................................................*/

/*...*/
  #if _MKL_
    #if defined(MKL_ILP64)
      #define MKL_INT long long
    #else
      #define MKL_INT int
    #endif  
    extern MKL_INT pardiso
    (void *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
    double *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
    MKL_INT *, double *, double *, MKL_INT *);
    extern void mkl_set_num_threads(MKL_INT *nt);
  #endif
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
  
/*...*/
  #define SING1(x) x > 0.0 ? 1.e0:-1.e0 
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
              ,INT *RESTRICT ia   ,INT *RESTRICT ja
              ,DOUBLE *RESTRICT al,DOUBLE *RESTRICT ad
              ,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b
              ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT z 
              ,DOUBLE *RESTRICT r ,DOUBLE *RESTRICT p
              ,DOUBLE const tol   ,unsigned int maxit 
              ,bool const newX    ,FILE* fileSolvLog  
              ,bool const log     ,bool const fPrint  
              ,Interface *iNeq    ,BufferOmp *bOmp
              ,void(*matvec)()    ,DOUBLE(*dot)());

/*... gradiente conjugado precondicionado*/
	void pcg(INT const nEq,INT const nAd
		,INT *RESTRICT ia   ,INT *RESTRICT ja
		,DOUBLE *RESTRICT al,DOUBLE *RESTRICT ad
		,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b
    ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT z
    ,DOUBLE *RESTRICT r ,DOUBLE *RESTRICT p
		,DOUBLE const tol   ,unsigned int maxIt
    ,bool const newX  	,FILE* fLog
    ,FILE *fileHistLog  ,bool const log 
    ,bool const fHistLog,bool const fPrint
		,void(*matvec)()    ,DOUBLE(*dot)());

/*...  gradiente conjugado precondicionado (OPENMP)*/  
  void pcgOmp(INT const nEq      ,INT const nAd
             ,INT *RESTRICT ia   ,INT *RESTRICT ja
             ,DOUBLE *RESTRICT a ,DOUBLE *RESTRICT ad
             ,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b
             ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT z
             ,DOUBLE *RESTRICT r ,DOUBLE *RESTRICT p
             ,DOUBLE const tol   ,unsigned int maxIt
             ,bool const newX    ,FILE* fLog
             ,FILE *fileHistLog  ,bool const log
             ,bool const fHistLog, bool const fPrint
             ,BufferOmp *bOmp
             ,void(*matvec)(), DOUBLE(*dot)());

/*...  gradiente conjugado precondicionado (MPI)*/
  void mpiPcg(INT const nEq   ,INT const nEqNov
        ,INT const nAd      ,INT const nAdR  
        ,INT *RESTRICT ia   ,INT *RESTRICT ja
        ,DOUBLE *RESTRICT al,DOUBLE *RESTRICT ad
        ,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b 
        ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT z 
        ,DOUBLE *RESTRICT r ,DOUBLE *RESTRICT p
        ,DOUBLE const tol   ,unsigned int maxIt
        ,bool const newX    ,FILE* fileLog  
        ,FILE *fileHistLog  ,bool const log     
        ,bool const fHistLog,bool const fPrint
        ,Interface *iNeq                      
        ,void(*matvec)()    ,DOUBLE(*dot)());
/*...................................................................*/

/*...*/
  void callBicgStab(INT const nEq     ,INT const nEqNov
                   ,INT const nAd     ,INT const nAdR
                   ,INT *RESTRICT ia  ,INT *RESTRICT ja
                   ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                   ,DOUBLE *RESTRICT m,DOUBLE *RESTRICT b
                   ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT t
                   ,DOUBLE *RESTRICT v,DOUBLE *RESTRICT r
                   ,DOUBLE *RESTRICT p,DOUBLE *RESTRICT z
                   ,DOUBLE *RESTRICT h
                   ,DOUBLE const tol, unsigned int maxIt
                   ,bool const newX, FILE* fSolvLog
                   ,bool const fLog, bool const fPrint
                   ,Interface *iNeq, BufferOmp *bOmp
                   ,void(*matVec)(), DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado*/
  void pbicgstab(INT const nEq    ,INT const nAd
       ,INT *RESTRICT ia   ,INT *RESTRICT ja
       ,DOUBLE *RESTRICT al,DOUBLE *RESTRICT ad
       ,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b 
       ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT t 
       ,DOUBLE *RESTRICT v ,DOUBLE *RESTRICT r
       ,DOUBLE *RESTRICT p ,DOUBLE *RESTRICT z 
       ,DOUBLE *RESTRICT r0,DOUBLE const tol
       ,unsigned int maxIt ,bool const newX          
	     ,FILE* fLog         ,FILE *fileHistLog
	     ,bool const log     , bool const fHistLog
       ,bool const fPrint
       ,void(*matvec)()    ,DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado (OMP)*/
  void pbicgstabOmp(INT const nEq, INT const nAd
                   ,INT *RESTRICT ia, INT *RESTRICT ja
                   ,DOUBLE *RESTRICT a, DOUBLE *RESTRICT ad
                   ,DOUBLE *RESTRICT m, DOUBLE *RESTRICT b
                   ,DOUBLE *RESTRICT x, DOUBLE *RESTRICT t
                   ,DOUBLE *RESTRICT v, DOUBLE *RESTRICT r
                   ,DOUBLE *RESTRICT p, DOUBLE *RESTRICT z
                   ,DOUBLE *RESTRICT r0
                   ,DOUBLE const tol, unsigned int maxIt
                   ,bool const newX, FILE* fLog
                   ,FILE *fileHistLog, bool const log
                   ,bool const fHistLog, bool const fPrint
                   ,BufferOmp *bOmp
                   ,void(*matvec)(), DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado (MPI)*/
  void mpiPbicgstab(INT const nEq,INT const nEqNov      
          ,INT const nAd       ,INT const nAdR
          ,INT *RESTRICT ia    ,INT *RESTRICT ja
          ,DOUBLE *RESTRICT al ,DOUBLE *RESTRICT ad
          ,DOUBLE *RESTRICT m  ,DOUBLE *RESTRICT b 
          ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT t
          ,DOUBLE *RESTRICT v  ,DOUBLE *RESTRICT r
          ,DOUBLE *RESTRICT p  ,DOUBLE *RESTRICT z 
          ,DOUBLE *RESTRICT r0
          ,DOUBLE const tol    ,unsigned int maxIt
          ,bool const newX     ,FILE* fileLog 
          ,FILE *fileHistLog   ,bool const log
          ,bool const fHistLog ,bool const fPrint 
          ,Interface *iNeq    
          ,void(*matvec)()     ,DOUBLE(*dot)());
/*...................................................................*/

/*...*/
  void callBicgStabl2(INT const nEq, INT const nEqNov
                      , INT const nAd, INT const nAdR
                      , INT *RESTRICT ia, INT *RESTRICT ja
                      , DOUBLE *RESTRICT a, DOUBLE *RESTRICT ad
                      , DOUBLE *RESTRICT m, DOUBLE *RESTRICT b
                      , DOUBLE *RESTRICT x, DOUBLE *RESTRICT t
                      , DOUBLE *RESTRICT v, DOUBLE *RESTRICT r
                      , DOUBLE *RESTRICT u, DOUBLE *RESTRICT r0
                      , DOUBLE *RESTRICT w, DOUBLE *RESTRICT s
                      , DOUBLE *RESTRICT p, DOUBLE *RESTRICT h
                      , DOUBLE *RESTRICT z
                      , DOUBLE const tol, unsigned int maxIt
                      , bool const newX, FILE* fSolvLog
                      , bool const fLog, bool const fPrint
                      , Interface *iNeq, BufferOmp *bOmp
                      , void(*matVec)(), DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado l=2*/
  void pbicgstabl2(INT const nEq, INT const nAd
                   , INT *RESTRICT ia, INT *RESTRICT ja
                   , DOUBLE *RESTRICT a, DOUBLE *RESTRICT ad
                   , DOUBLE *RESTRICT m, DOUBLE *RESTRICT b
                   , DOUBLE *RESTRICT x, DOUBLE *RESTRICT t
                   , DOUBLE *RESTRICT v, DOUBLE *RESTRICT r
                   , DOUBLE *RESTRICT u, DOUBLE *RESTRICT r0
                   , DOUBLE *RESTRICT w, DOUBLE *RESTRICT s
                   , DOUBLE *RESTRICT p, DOUBLE *RESTRICT h
                   , DOUBLE *RESTRICT z
                   , DOUBLE const tol, unsigned int maxIt
                   , bool const newX, FILE* fLog
                   , FILE *fileHistLog, bool const log
                   , bool const fHistLog, bool const fPrint
                   , void(*matvec)(), DOUBLE(*dot)());
/*...................................................................*/

/*... gradiente conjugado bi-ortoganilizados precondicionado (OMP)*/
  void pbicgstabl2Omp(INT const nEq, INT const nAd
                      , INT *RESTRICT ia, INT *RESTRICT ja
                      , DOUBLE *RESTRICT a, DOUBLE *RESTRICT ad
                      , DOUBLE *RESTRICT m, DOUBLE *RESTRICT b
                      , DOUBLE *RESTRICT x, DOUBLE *RESTRICT t
                      , DOUBLE *RESTRICT v, DOUBLE *RESTRICT r
                      , DOUBLE *RESTRICT u, DOUBLE *RESTRICT r0
                      , DOUBLE *RESTRICT w, DOUBLE *RESTRICT s
                      , DOUBLE *RESTRICT p, DOUBLE *RESTRICT h
                      , DOUBLE *RESTRICT z
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
                    , INT *RESTRICT ia, INT *RESTRICT ja
                    , DOUBLE *RESTRICT al, DOUBLE *RESTRICT ad
                    , DOUBLE *RESTRICT m, DOUBLE *RESTRICT b
                    , DOUBLE *RESTRICT x, DOUBLE *RESTRICT t
                    , DOUBLE *RESTRICT v, DOUBLE *RESTRICT r
                    , DOUBLE *RESTRICT p, DOUBLE *RESTRICT z
                    , DOUBLE const tol
                    , unsigned int maxIt, bool const newX
                    , FILE* fLog, bool const log
                    , bool const fPrint, Interface *iNeq
                    , void(*matvec)(), DOUBLE(*dot)());*/
/*...................................................................*/

/*...*/
  void callMinRes(INT const nEq     , INT const nEqNov
               , INT const nAd      , INT const nAdR
               , INT *RESTRICT ia   , INT *RESTRICT ja
               , DOUBLE *RESTRICT al, DOUBLE *RESTRICT ad
	             , DOUBLE *RESTRICT m , DOUBLE *RESTRICT b
               , DOUBLE *RESTRICT x , DOUBLE *RESTRICT v0
               , DOUBLE *RESTRICT v , DOUBLE *RESTRICT w
               , DOUBLE *RESTRICT w0, DOUBLE *RESTRICT w00
               , DOUBLE *RESTRICT z , DOUBLE *RESTRICT z0
               , DOUBLE *RESTRICT p
               , DOUBLE const tol   , unsigned int maxIt
               , bool const newX    , FILE* fSolvLog
               , bool const fLog    , bool const fPrint
               , Interface *iNeq    , BufferOmp *bOmp
               , void(*matVec)()    , DOUBLE(*dot)());
/*...*/
  void minres(INT const nEq        , INT const nAd
	        , INT *RESTRICT ia    , INT *RESTRICT ja
	        , DOUBLE *RESTRICT al , DOUBLE *RESTRICT ad
	        , DOUBLE *RESTRICT b  , DOUBLE *RESTRICT x 
          , DOUBLE *RESTRICT v0 , DOUBLE *RESTRICT v
          , DOUBLE *RESTRICT w  , DOUBLE *RESTRICT w0
          , DOUBLE *RESTRICT w00, DOUBLE *RESTRICT z 
          , DOUBLE *RESTRICT p
	        , DOUBLE const tol   , unsigned int maxIt
          , bool const newX    , FILE* fileLog   
          , FILE *fileHistLog	, bool const log 
          , bool const fHistLog, bool const fPrint
	        , void(*matvec)()    , DOUBLE(*dot)());

/*...*/
  void pminres(INT const nEq    , INT const nAd
	         , INT *RESTRICT ia   , INT *RESTRICT ja
	         , DOUBLE *RESTRICT al, DOUBLE *RESTRICT ad
	         , DOUBLE *RESTRICT m , DOUBLE *RESTRICT b
           , DOUBLE *RESTRICT x , DOUBLE *RESTRICT v0
           , DOUBLE *RESTRICT v , DOUBLE *RESTRICT w
           , DOUBLE *RESTRICT w0, DOUBLE *RESTRICT w00
           , DOUBLE *RESTRICT z , DOUBLE *RESTRICT z0
           , DOUBLE *RESTRICT p
	         , DOUBLE const tol   , unsigned int maxIt
           , bool const newX    , FILE* fileLog   
           , FILE *fileHistLog	, bool const log 
           , bool const fHistLog, bool const fPrint
	         , void(*matvec)()    , DOUBLE(*dot)());
/*...................................................................*/

/*...*/
  void callGmres(INT const nEq, INT const nEqNov
                 , INT const nAd, INT const nAdR
                 , INT *RESTRICT ia, INT *RESTRICT ja
                 , DOUBLE *RESTRICT a, DOUBLE *RESTRICT ad
                 , DOUBLE *RESTRICT m, DOUBLE *RESTRICT b
                 , DOUBLE *RESTRICT x, DOUBLE *RESTRICT g
                 , DOUBLE *RESTRICT h, DOUBLE *RESTRICT y
                 , DOUBLE *RESTRICT c, DOUBLE *RESTRICT s
                 , DOUBLE *RESTRICT e, short const nKrylov
                 , DOUBLE const tol  , unsigned int maxIt
                 , bool const newX   , FILE* fSolvLog
                 , bool const fLog   , bool const fPrint
                 , Interface *iNeq  , BufferOmp *bOmp
                 , void(*matVec)()  , DOUBLE(*dot)());
/*...................................................................*/

/*... GMRES*/
  void gmres(INT const nEq       ,INT const nAd
            ,INT *RESTRICT ia    ,INT *RESTRICT ja
            ,DOUBLE *RESTRICT a  ,DOUBLE *RESTRICT ad
            ,DOUBLE *RESTRICT m  ,DOUBLE *RESTRICT b
            ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT g
            ,DOUBLE *RESTRICT h  ,DOUBLE *RESTRICT y
            ,DOUBLE *RESTRICT c  ,DOUBLE *RESTRICT s
            ,DOUBLE *RESTRICT e
            ,short const nKrylov
            ,DOUBLE const tol    ,unsigned int nCycles
            ,bool const newX     ,FILE* fLog
            ,FILE *fileHistLog   ,bool const log
            ,bool const fHistLog ,bool const fPrint
            ,void(*matvec)()     ,DOUBLE(*dot)());
/*...................................................................*/

/*... Omp*/
  void gmresOmp(INT const nEq, INT const nAd
                , INT *RESTRICT ia, INT *RESTRICT ja
                , DOUBLE *RESTRICT a, DOUBLE *RESTRICT ad
                , DOUBLE *RESTRICT m, DOUBLE *RESTRICT b
                , DOUBLE *RESTRICT x, DOUBLE *RESTRICT g
                , DOUBLE *RESTRICT h, DOUBLE *RESTRICT y
                , DOUBLE *RESTRICT c, DOUBLE *RESTRICT s
                , DOUBLE *RESTRICT e, short const nKrylov
                , DOUBLE const tol, unsigned int nCycles
                , bool const newX, FILE* fLog
                , FILE *fileHistLog, bool const log
                , bool const fHistLog, bool const fPrint
                , BufferOmp *bOmp
                , void(*matvec)(), DOUBLE(*dot)());
/*...................................................................*/

/*...*/
  void callMklPardiso(INT  nEq          , INT  mtype
                    , INT   *RESTRICT ia, INT   *RESTRICT ja
                    , DOUBLE *RESTRICT a, DOUBLE *RESTRICT b
                    , DOUBLE *RESTRICT x, DOUBLE *RESTRICT z  
                    , DOUBLE *RESTRICT r                        
                    , bool const fPrint);
/*...................................................................*/
/*===================================================================*/

#endif/*_SOLV_H_*/
