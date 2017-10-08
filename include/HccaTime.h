#ifndef _MYTIME_H_
  #define _MYTIME_H_
  #if _OPENMP
    #include<omp.h>
  #elif _WIN32 
    #include<time.h> 
  #elif _MPICH_ 
    #include<mpi.h>
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
    double solvEdpD1;
    double numeqD1;
    double dataStructD1;
    double CellPloadD1;
    double CellTransientD1;
    double systFormD1;
    double rcGradD1;
/*... T1*/
    double solvT1;
    double solvEdpT1;
    double numeqT1;
    double dataStructT1;
    double CellPloadT1;
    double CellTransientT1;
    double systFormT1;
    double rcGradT1;
/*... fluid*/
    double solvPres;
    double solvVel;
    double solvEnergy;
    double numeqPres;
    double numeqVel;
    double numeqEnergy;
    double dataStructVel;
    double dataStructPres;
    double dataStructEnergy;
    double solvEdpFluid;
    double cellPloadSimple;
    double cellTransientSimple;
    double systFormPres;                                 
    double systFormVel;
    double systFormEnergy;
    double velExp;
    double rcGradPres;
    double rcGradVel;
    double rcGradEnergy;
/*... blas*/
    double matVecOverHeadMpi;
    double matVecSparse;
    double dot;
    double dotOverHeadMpi;
/*... iterativos*/
    double pcg;
    double pbicgstab;
    double gmres;
    double minres;
/*... precondicionador*/
    double precondDiag;
/*... particionamento*/
    double partdMesh;
    double partdMeshCom;
/*... comunicacao entre as particoes*/
    double overHeadCelMpi;
    double overHeadNodMpi;
    double overHeadNeqMpi;
    double overHeadGNodMpi;
    double overHeadGCelMpi;
    double overHeadTotalMpi;
/*...*/
    double tempForEnergy;
/*...*/
    double turbulence;
/*...*/
    double total;
  }Time;
  Time tm;
  double getTimeC(void);
#endif/*_MYTIME_H_*/
