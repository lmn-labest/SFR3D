#ifndef _REORD_H  
  #define _REORD_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
  #include<Mystdbool.h>
  #include<Memoria.h>
  #include<Graph.h>
  #include<Rcm.h>

  typedef struct Reord{
    long *num;
    long flag;
  }Reord;
  
  void reord(Memoria *m  ,long *num  ,long const *adj, short const *nViz
            ,short const maxViz,long numel , bool flag);


#endif/*_REORD_H*/
