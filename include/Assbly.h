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
  void assbly(INT    *RESTRICT  ia,INT *RESTRICT ja 
             ,DOUBLE *RESTRICT al ,DOUBLE *RESTRICT ad
             ,DOUBLE *RESTRICT b  ,INT *RESTRICT lId
             ,DOUBLE *RESTRICT lA ,DOUBLE *RESTRICT lB
             ,INT const nEq       ,INT const nEqNov 
             ,INT const nAd       ,INT const nAdr  
             ,short const nFace   ,short const ndf 
             ,short const storage ,bool  const forces  
             ,bool const matrix   ,bool  const  unsym);
/*...................................................................*/

/*...*/
  void assblySimple(INT    *RESTRICT  ia,INT *RESTRICT ja 
             ,DOUBLE *RESTRICT al ,DOUBLE *RESTRICT ad
             ,DOUBLE *RESTRICT b  ,INT *RESTRICT lId
             ,DOUBLE *RESTRICT lA ,DOUBLE *RESTRICT lB
             ,INT const nEq       ,INT const nEqNov 
             ,INT const nAd       ,INT const nAdr  
             ,short const nFace   ,short const ndf 
             ,short const storage ,bool  const forces  
             ,bool const matrix   ,bool  const  unsym);
/*...................................................................*/

#endif/*_CSR_H*/
