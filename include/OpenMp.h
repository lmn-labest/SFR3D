#ifndef _OPENMP_H_
    #define _OPENMP_H_

/*...*/
    #include<omp.h>
/*...................................................................*/

/*...*/
    #include<Define.h>
    #include<HccaStdBool.h>
    #include<Memoria.h>
    #include<Sisteq.h>
/*...................................................................*/

    extern DOUBLE tmpDotOmp;
    extern DOUBLE tmpDotOmp1;
    extern DOUBLE tmpDotOmp2;
    extern DOUBLE tmpDotOmp3;
    extern DOUBLE tmpDotOmp4;
    extern DOUBLE tmpDotOmp5;
    extern DOUBLE tmpDotOmp6;
    extern DOUBLE tmpDotOmp7;
    extern DOUBLE tmpDotOmp8;

/*... OPENMP*/
  typedef struct {
    bool           flag,fSolver,fCell,fUpdate,fGrad,fReaction;
    unsigned short nThreadsSolver;
    unsigned short nThreadsCell;
    unsigned short nThreadsUpdate;
    unsigned short nThreadsGrad;
    unsigned short nThreadsReaction;
    INT            nEqMax;
    DOUBLE         *buffer;
  }Omp;
/*...................................................................*/

  extern Omp ompVar;

  void openMpCheck(bool omp);

  void pMatrixSolverOmp(Memoria *m,SistEq *eq
                       ,char *s1  ,char *s2
                       ,char *s3  ,char *s4);

  void openMpSet(FILE *fileIn, Omp *ompVar);
#endif /*_OPENMP_H_*/
