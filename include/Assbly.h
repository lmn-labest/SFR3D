#ifndef _ASSBLY_H
  #define _ASSBLY_H
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
/*...................................................................*/

/*...*/
  #include<Csr.h>
  #include<EllPack.h>
  #include<Erro.h>
  #include<HccaStdBool.h>
  #include<Define.h>
/*...................................................................*/

/*...*/
  void assbly(INT    *restrict  ia,INT *restrict ja 
             ,DOUBLE *restrict au ,DOUBLE *restrict ad
             ,DOUBLE *restrict al ,DOUBLE *restrict b
             ,INT *restrict lId
             ,DOUBLE *restrict lA ,DOUBLE *restrict lB
             ,INT const nEq       ,INT const nAd  
             ,short const nFace   ,short const ndf 
             ,short const storage ,bool  const forces  
             ,bool const matrix   ,bool  const  unsym);
/*...................................................................*/

#endif/*_CSR_H*/
