#ifndef _ASSBLY_H
  #define _ASSBLY_H
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
/*...................................................................*/

/*...*/
  #include<Csr.h>
  #include<HccaStdBool.h>
  #include<Define.h>
/*...................................................................*/

/*...*/
  void assbly(INT    *restrict  ia,INT *restrict ja 
             ,double *restrict au ,double *restrict ad
             ,double *restrict al ,double *restrict b
             ,INT *restrict lId
             ,double *restrict lA ,double *restrict lB
             ,short const nFace   ,short const ndf 
             ,short const storage ,bool  const forces  
             ,bool const matrix   ,bool  const  unsym);
/*...................................................................*/

#endif/*_CSR_H*/
