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
  
    DOUBLE tmpDotOmp;
    DOUBLE tmpDotOmp1;
    DOUBLE tmpDotOmp2;
    DOUBLE tmpDotOmp3;
    DOUBLE tmpDotOmp4;
    DOUBLE tmpDotOmp5;
    DOUBLE tmpDotOmp6;
    DOUBLE tmpDotOmp7;
    DOUBLE tmpDotOmp8;
  
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
  Omp ompVar;
  
  void openMpCheck(bool omp);

  void pMatrixSolverOmp(Memoria *m,SistEq *eq
                       ,char *s1  ,char *s2
                       ,char *s3  ,char *s4);

  void openMpSet(FILE *fileIn, Omp *ompVar);  

#endif /*_OPENMP_H_*/
