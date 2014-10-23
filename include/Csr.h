#ifndef _CSR_H  
  #define _CSR_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
  #include<stdlib.h>
  #include<Mystdbool.h>
  #include<Define.h>
/*...*/
  INT csrIa(INT *ia  ,INT *id    ,INT *num   ,INT  *adj ,short *nViz
            ,INT numel,INT neq    ,short maxViz,short ndf  ,bool upper
            ,bool diag , bool lower);
/*...*/
  void csrJa(INT *ia     ,INT *ja 
           ,INT *id  ,INT *num ,INT  *adj, short *nViz
           ,INT numel,INT neq ,short maxViz,short ndf
           ,bool upper,bool diag,bool lower);

/*...*/
  INT bandCsr(INT *ia,INT *ja,INT  neq,short type);

/*...*/
  void sortGraphCsr(INT *ia,INT *ja,INT n);

#endif/*_CSR_H*/
