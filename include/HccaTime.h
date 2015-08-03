#ifndef _MYTIME_H_
  #define _MYTIME_H_
  #if _OPENMP
    #include<Openmp.h>
  #elif _WIN32 
    #include<time.h>  
  #else  
    #include<sys/time.h>
  #endif
  #include<stdio.h>
  typedef struct{
    double adjcency;
    double geom;
    double leastSquareMatrix;
    double reord;
/*... D1*/
    double solvD1;
    double solvEdoD1;
    double numeqD1;
    double dataStructD1;
    double CellPloadD1;
    double CellTransientD1;
    double systFormD1;
    double rcGradD1;
/*... blas*/
    double matVecSparse;
    double dot;
/*... iterativos*/
    double pcg;
/*... precondicionador*/
    double precondDiag;

/*...*/
    double total;
  }Time;
  Time tm;
  double getTimeC(void);
#endif/*_MYTIME_H_*/
