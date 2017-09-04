#ifndef _REORD_H_  
  #define _REORD_H_
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false

/*...*/
  #include<HccaStdBool.h>
  #include<Memoria.h>
  #include<Graph.h>
  #include<Rcm.h>
/*...................................................................*/

/*...*/
  typedef struct{
    INT *num;
    INT flag;
  }Reord;
/*...................................................................*/
  
/*...*/
  void reord(Memoria *m         ,INT *num          ,INT const *adj
            ,short const *nViz  ,short const maxViz
            ,INT const numel    ,INT const numelNov
            ,bool const flag    ,unsigned  short const nPrcs);
/*...................................................................*/


#endif/*_REORD_H_*/
