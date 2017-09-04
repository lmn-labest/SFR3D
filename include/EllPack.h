#ifndef _ELLPACK_H_  
  #define _ELLPACK_H_
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
  INT ellPackJa(INT *RESTRICT ifEllPack,INT *RESTRICT ja
               ,INT *RESTRICT id       ,INT *RESTRICT num
               ,INT *RESTRICT adj      ,short *RESTRICT nViz
               ,INT const numel        ,INT const nEq    
               ,short const maxLineNzr ,const short ndf);
/*...................................................................*/

/*...*/
 void ellPack(INT    *RESTRICT  ifEllPack
        ,DOUBLE *RESTRICT ad            ,DOUBLE *RESTRICT a
        ,DOUBLE *RESTRICT b
        ,INT *RESTRICT lId                       
        ,DOUBLE *RESTRICT lA            ,DOUBLE *RESTRICT lB 
        ,short const nFace              ,short const ndf  
        ,short const storage            ,bool  const forces
        ,bool const matrix              );
/*...................................................................*/

/*...*/
  INT bandEllPack(INT *ifEllPack,INT *ja,INT  neq,short type);
/*...................................................................*/


#endif/*_ELLPACK_H_*/
