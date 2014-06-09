#ifndef _REORD_H  
  #define _REORD_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
  #include<Mystdbool.h>

  typedef struct Reord{
    long *num;
    long flag;
  }Reord;
  
  void reord(long *num  ,long *graph, long numel, bool flag);

/*...*/
  void bubblesort(long *ja,long n);

#endif/*_REORD_H*/
