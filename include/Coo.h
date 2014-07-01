#ifndef _COO_H  
  #define _COO_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
  #ifdef _MMIO_
    #include<mmio/mmio.h>
  #endif
  #include<Memoria.h>
  #include<stdio.h>
  #include<stdlib.h>
  #include<Define.h>
/*...*/
  void writeCoo(Memoria *m,INT *ia   ,INT *ja,INT neq
               ,INT nad,short type
               ,bool unsym,char* name);

/*...*/
  void csrToCoo(int  *linha,int  *col,double *val,INT *ia,INT *ja
               ,INT neq ,INT nad);


#endif/*_COO_H*/
