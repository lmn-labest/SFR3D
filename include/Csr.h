#ifndef _CSR_H  
  #define _CSR_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
  #include<Mystdbool.h>
  #include<Reord.h>
  #include<Mesh.h>
/*...*/
  long csrIa(long *ia  ,long *id    ,long *num   ,long  *adj , short *nViz
            ,long numel,long neq    ,short maxViz, bool upper
            ,bool diag , bool lower);
/*...*/
  void csrJa(long *ia     ,long *ja 
           ,long *id    ,long *num   ,long  *adj, short *nViz
           ,long numel,long neq    ,short maxViz, bool upper
           ,bool diag , bool lower);
/*...*/
  void sortGraphCsr(long *ia,long *ja,long n);

#endif/*_CSR_H*/
