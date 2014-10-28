#ifndef _HCCABLAS_H
  #define _HCCABLAS_H

/*...*/  
  #include<HccaStdBool.h>
  #include<Define.h>
/*...................................................................*/

/*...*/
  void prodVet(double *restrict a
              ,double *restrict b
              ,double *restrict c);
/*...................................................................*/

/*... level 1*/
  void addVector(double const alpha,double *restrict a
                ,double const beta ,double *restrict b
                ,INT const nDim    ,double *restrict c);
/*...................................................................*/
#endif
