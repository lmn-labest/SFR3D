#ifndef _ASSBLY_H_
  #define _ASSBLY_H_
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

/*...*/
  void assblyBlock(INT    *RESTRICT  ia, INT *RESTRICT ja
                 , DOUBLE *RESTRICT a  , DOUBLE *RESTRICT ad
                 , DOUBLE *RESTRICT b  , INT *RESTRICT lId
                 , DOUBLE *RESTRICT lA , DOUBLE *RESTRICT lB
                 , INT const nEq       , INT const nEqNov
                 , INT const nAd       , INT const nAdR
                 , short const nFace   , short const nForceBlock
                 , short const nAdBlock, short const nAlBlock
                 , short const storage , bool  const forces
                 , bool const matrix   , bool  const  unsym);
/*...................................................................*/
#endif/*_ASSBLY_H_*/
