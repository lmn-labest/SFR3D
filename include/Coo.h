#ifndef _COO_H  
  #define _COO_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
  #include<mmio/mmio.h>
  #include<Memoria.h>
  #include<stdio.h>
  #include<stdlib.h>
/*...*/
  void writeCoo(Memoria *m,long *ia   ,long *ja,long neq
               ,long nad,short type
               ,bool unsym,char* name);

/*...*/
  void csrToCoo(int  *linha,int  *col,double *val,long *ia,long *ja
               ,long neq ,long nad);


#endif/*_COO_H*/
