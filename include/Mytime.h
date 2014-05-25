#ifndef _MYTIME_H_
  #define _MYTIME_H_
  #if _OPENMP
    #include<Openmp.h>
  #else  
    #include <time.h>
  #endif  
  typedef struct Time{
    double solv  ;
    double dot   ;
    double matvec;
    double pform;
    double tform;
    double color;
    double total;
  }Time;
  Time tm;
  double gettimec(void);
#endif/*_MYTIME_H_*/
