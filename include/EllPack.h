#ifndef _ELLPACK_H  
  #define _ELLPACK_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
  #include<string.h>
/*...................................................................*/

/*...*/
  #include<Define.h>
  #include<Erro.h>
  #include<HccaStdBool.h>
/*...................................................................*/

/*...*/
  INT ellPackJa(INT *restrict ifEllPack,INT *restrict ja
               ,INT *restrict id       ,INT *restrict num
               ,INT *restrict adj      ,short *restrict nViz
               ,INT const numel        ,INT const nEq    
               ,short const maxLineNzr ,const short ndf);
/*...................................................................*/

/*...*/
 void ellPack(INT    *restrict  ifEllPack
        ,DOUBLE *restrict ad            ,DOUBLE *restrict a
        ,DOUBLE *restrict b
        ,INT *restrict lId                       
        ,DOUBLE *restrict lA            ,DOUBLE *restrict lB 
        ,short const nFace              ,short const ndf  
        ,short const storage            ,bool  const forces
        ,bool const matrix              );
/*...................................................................*/

/*...*/
  INT bandEllPack(INT *ifEllPack,INT *ja,INT  neq,short type);
/*...................................................................*/


#endif/*_ELLPACK_H*/
