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
    double solv  ;
    double dot   ;
    double matvec;
    double pform;
    double tform;
    double color;
    double total;
  }Time;
  Time tm;
  double getTimeC(void);
#endif/*_MYTIME_H_*/
