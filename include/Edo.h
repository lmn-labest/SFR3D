#ifndef _EDO_H_
  #define _EDO_H_
/*...*/
  #include<stdlib.h>
  #include<stdio.h>
  #include<math.h>
/*...................................................................*/

/*...*/ 
  #include<Define.h> 
  #include<Gauss.h>
/*...................................................................*/

  #define SCREEN_OUT 1
  #define FILE_OUT   2  
 
  int StepperSie(DOUBLE *y                , void **pt
               , DOUBLE const x1          , DOUBLE const x2
               , DOUBLE const hInit       , DOUBLE const hMax
               , DOUBLE const aTol        , DOUBLE const rTol    
               , short const nEdo         , DOUBLE *xF
               , unsigned INT const maxInt, bool const fStopIt
               , short const  outCod      , const char *const fName
               , void(*rhs)()             , void(*jacY)());

  void writeOutput(int const it     , double const x
                , double const h  , double *y
                , short const nEdo, short const outCod
                , FILE *file );

  double edoError(DOUBLE *y        ,DOUBLE *yOut
               ,DOUBLE *yErr     ,DOUBLE const atol
               ,DOUBLE const rtol,short const nEdo);

  DOUBLE sqr(DOUBLE const x); 

#endif
/*_EDO_H_*/