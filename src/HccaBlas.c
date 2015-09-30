#include<HccaBlas.h>
/********************************************************************** 
 * HCCABLAS : biblioteca de algebra linear                            * 
 *------------------------------------------------------------------- * 
 * nome das funcoes                                                   * 
 *------------------------------------------------------------------- * 
 * flopMatVecCsr  -> FLOP para operacao matriz*vector, onde a matriz  *
 * esta armazenada em um formato da familia CSR                       *
 * flopDot        -> FLOP para operacao produto interno               *
 * pordVet        -> produto vetorial                                 * 
 * ------------------------- VECTOR --------------------------------- * 
 * alphaProdVector-> produto de vetor por um escalar                  *
 * addVector      -> adicao de vetores                                *
 * ------------------------- DOT ------------------------------------ * 
 * dot            -> produto interno entre dois vetores               * 
 * dotO2L2        -> produto interno entre dois vetores               * 
 * dotL2          -> produto interno entre dois vetores               * 
 * dotL4          -> produto interno entre dois vetores               * 
 * dotL6          -> produto interno entre dois vetores               * 
 * dotL8          -> produto interno entre dois vetores               * 
 * dotO2          -> produto interno entre dois vetores               * 
 * dotO4          -> produto interno entre dois vetores               * 
 * dotO6          -> produto interno entre dois vetores               * 
 * dotO8          -> produto interno entre dois vetores               * 
 * -------------------------- CSRD ---------------------------------- *
 * ....................... SIMETRICA ................................ *
 * mpiMatVecCsrDSym -> matizt vetor para matriz simetrica no formato  *           
 * CSRD                                                               *
 * ................... SIMETRICA - MPI .............................. *
 * mpiMatVecCsrDSym -> matizt vetor para matriz simetrica no formato  *           
 * CSRD                                                               *           
 * ...........................GERAL ................................. *
 * matVecCsrD     -> matriz vetor para matriz geral no formato CSRD   * 
 * matVecCsrDI2   -> matriz vetor para matriz geral no formato CSRD   * 
 * matVecCsrDI4   -> matriz vetor para matriz geral no formato CSRD   * 
 * matVecCsrDI6   -> matriz vetor para matriz geral no formato CSRD   * 
 * matVecCsrDO2   -> matriz vetor para matriz geral no formato CSRD   * 
 * matVecCsrDO4   -> matriz vetor para matriz geral no formato CSRD   * 
 * matVecCsrDO6   -> matriz vetor para matriz geral no formato CSRD   * 
 * matVecCsrDO2I2 -> matriz vetor para matriz geral no formato CSRD   * 
 * ------------------------ ELLPACK --------------------------------- * 
 * ...........................GERAL ................................. *
 * matVecEllPack  -> matriz vetor para matriz geral no formato EllPack* 
 * matVecEllPackO2-> matriz vetor para matriz geral no formato EllPack* 
 * matVecEllPackO4-> matriz vetor para matriz geral no formato EllPack* 
 *********************************************************************/

/********************************************************************* 
 * FLOPDOT: Calculo do numero de flops do produto interno            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * nDim -> numero de dimensoes                                       * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * Numero de flops                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
long flopDot(INT nDim){
  return (2*nDim-1);
}
/*********************************************************************/ 

/********************************************************************* 
 * FLOPMATVECCSR  : Calculo do numero de flops da operacao de matiz  * 
 * por um vetor, onde a matriz e uma matriz esparsa no formato CSR   * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq  -> numero de linhas                                          * 
 * nad  -> numero total de elementos nao nulos                       * 
 * ty   -> tipo do formato csr                                       * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * Numero de flops                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
INT flopMatVecCsr(INT neq,INT nad,short ty){
  
  INT flop;
  switch(ty){
/*... armazenamento CSR(a)*/
    case CSR:
      flop = 2*nad;
    break;
/*... armazenamento CSRD(a,ad)*/
    case CSRD:
      flop = 2*nad + neq;
    break;
/*... armazenamento CSRC(ad,au,al)*/
    case CSRC:
      flop = 4*nad + neq;
    break;
/*...*/
    default:
      ERRO_OP(__FILE__,__func__,ty);
    break;
  }
/*...................................................................*/
  return flop;
}
/*********************************************************************/ 


/********************************************************************* 
 * PRODVET: produto vetorial                                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a  - vetor a                                                      * 
 * b  - vetor b                                                      * 
 * c  - indefinido                                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * c  - produto vetoria entre a e b                                  * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void prodVet(DOUBLE *restrict a,DOUBLE *restrict b, DOUBLE *restrict c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}
/********************************************************************/ 

/*======================== level 1 =================================*/

/********************************************************************* 
 * ALPAHAPRODVECTOR : produto de vetor por um escalar                * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * alpha -> escalar que multiplica o vetor a                         * 
 * a     -> vetor a                                                  * 
 * nDim  -> dimensao dos vetores                                     * 
 * c  - indefinido                                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * c  - produto vetoria entre alpha*a                                * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void alphaProdVector(DOUBLE const alpha,DOUBLE *restrict a
                    ,INT const nDim    ,DOUBLE *restrict c) 

{
  INT i;

  if(alpha == 1.0e0){
    for(i=0;i<nDim;i++)
      c[i] = a[i];
  }
  else 
    for(i=0;i<nDim;i++)
      c[i] = alpha*a[i];

}
/*********************************************************************/ 

/********************************************************************* 
 * ADDVECTOR : adicao de vetores                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * alpha -> escalar que multiplica o vetor a                         * 
 * a     -> vetor a                                                  * 
 * beto  -> escalar que multiplica o vetor b                         * 
 * b     -> vetor b                                                  * 
 * nDim  -> dimensao dos vetores                                     * 
 * c  - indefinido                                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * c  - produto vetoria entre alpha*a e beta*b                       * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void addVector(DOUBLE const alpha,DOUBLE *restrict a
              ,DOUBLE const beta ,DOUBLE *restrict b
              ,INT const nDim    ,DOUBLE *restrict c) 

{
  INT i;
   
  if(alpha == 1.e0)
    for(i=0;i<nDim;i++)
      c[i] = a[i] + beta*b[i];

  else if(beta == 1.e0)
    for(i=0;i<nDim;i++)
      c[i] = alpha*a[i] + b[i];

  else
    for(i=0;i<nDim;i++)
      c[i] = alpha*a[i] + beta*b[i];

}
/*********************************************************************/ 


/********************************************************************* 
 * DOT: produto interno entre dois vetores                           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x1 - vetor x1                                                     * 
 * x2 - vetor x2                                                     * 
 * n  - numero de dimensoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE dot(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n)
{
  INT i;
  DOUBLE tmp=0.e0,dot=0.e0;

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/

/*...*/
  for(i=0;i<n;i++){
    tmp += x1[i]*x2[i];
  }
/*...................................................................*/

/*...*/
#ifdef _MPICH_
  if(mpiVar.nPrcs>1){ 
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
    MPI_Allreduce(&tmp,&dot,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
  }
  else  
    dot = tmp;
#else
    dot = tmp;
#endif
/*...................................................................*/

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/

  return dot;
} 
/*********************************************************************/ 

/********************************************************************* 
 * DOTO2L2: produto interno entre dois vetores                       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x1 - vetor x1                                                     * 
 * x2 - vetor x2                                                     * 
 * n  - numero de dimensoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE dotO2L2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n)
{
  INT i,resto;
  DOUBLE tmp1=0.e0,tmp2=0.e0;
  DOUBLE tmp,dot=0.e0;
  
/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/

/*...*/ 
  resto = n%4;
  if(resto == 1)
     tmp1 = x1[0]*x2[0]; 
  else if(resto == 2)
     tmp1 = x1[0]*x2[0] + x1[1]*x2[1]; 
  else if(resto == 3)
     tmp1 = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]; 
    
  for(i=resto;i<n;i+=4){
    tmp1 += x1[i  ]*x2[i  ] + x1[i+1]*x2[i+1]; 
    tmp2 += x1[i+2]*x2[i+2] + x1[i+3]*x2[i+3]; 
  }
  tmp = tmp1 + tmp2;
/*...................................................................*/

  
/*...*/
#ifdef _MPICH_
  if(mpiVar.nPrcs>1){ 
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
    MPI_Allreduce(&tmp,&dot,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
  }
  else  
    dot = tmp;
#else
    dot = tmp;
#endif
/*...................................................................*/

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/
  return dot;
} 
/*********************************************************************/ 

/********************************************************************* 
 * DOTL2: produto interno entre dois vetores                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x1 - vetor x1                                                     * 
 * x2 - vetor x2                                                     * 
 * n  - numero de dimensoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE dotL2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n)
{
  INT i,resto;
  DOUBLE tmp1=0.e0,tmp2=0.e0;
  DOUBLE tmp,dot=0.e0;
  resto = n%2;

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/

/*...*/ 
  if(resto)
    tmp1 = x1[0]*x2[0];
    
  for(i=resto;i<n;i+=2){
    tmp1 += x1[i  ]*x2[i  ]; 
    tmp2 += x1[i+1]*x2[i+1];
  }
  tmp = tmp1 + tmp2;
/*...................................................................*/

/*...*/
#ifdef _MPICH_
  if(mpiVar.nPrcs>1){ 
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
    MPI_Allreduce(&tmp,&dot,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
  }
  else  
    dot = tmp;
#else
    dot = tmp;
#endif
/*...................................................................*/

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/

  return dot;
} 
/*********************************************************************/ 

/********************************************************************* 
 * DOTL4: produto interno entre dois vetores                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x1 - vetor x1                                                     * 
 * x2 - vetor x2                                                     * 
 * n  - numero de dimensoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE dotL4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n)
{
  INT i,resto;
  DOUBLE tmp1=0.e0,tmp2=0.e0,tmp3=0.e0,tmp4=0.e0;
  DOUBLE tmp,dot=0.e0;
  
/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/
  
/*...*/ 
  resto = n%4;
  if(resto == 1)
     tmp1 = x1[0]*x2[0]; 
  else if(resto == 2)
     tmp1 = x1[0]*x2[0] + x1[1]*x2[1]; 
  else if(resto == 3)
     tmp1 = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]; 
    
  
  for(i=resto;i<n;i+=4){
    tmp1 += x1[i  ]*x2[i  ];
    tmp2 += x1[i+1]*x2[i+1];
    tmp3 += x1[i+2]*x2[i+2];
    tmp4 += x1[i+3]*x2[i+3];
  }
  tmp = tmp1 + tmp2 + tmp3 + tmp4;
/*...................................................................*/

/*...*/
#ifdef _MPICH_
  if(mpiVar.nPrcs>1){ 
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
    MPI_Allreduce(&tmp,&dot,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
  }
  else  
    dot = tmp;
#else
    dot = tmp;
#endif
/*...................................................................*/

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/

  return dot;
} 
/*********************************************************************/ 

/********************************************************************* 
 * DOTL6: produto interno entre dois vetores                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x1 - vetor x1                                                     * 
 * x2 - vetor x2                                                     * 
 * n  - numero de dimensoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE dotL6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n)
{
  INT i,resto;
  DOUBLE tmp1=0.e0,tmp2=0.e0,tmp3=0.e0,tmp4=0.e0;
  DOUBLE tmp5=0.e0,tmp6=0.e0;
  DOUBLE tmp,dot=0.e0;

  
/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/
  
/*...*/ 
  resto = n%6;
  if(resto == 1)
     tmp1= x1[0]*x2[0]; 
  else if(resto == 2)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1]; 
  else if(resto == 3)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]; 
  else if(resto == 3)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]; 
  else if(resto == 4)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2] + x1[3]*x2[3]; 
  else if(resto == 5)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1] 
         + x1[2]*x2[2] + x1[3]*x2[3]  
         + x1[4]*x2[4];               
  
  for(i=resto;i<n;i+=6){
    tmp1 += x1[i  ]*x2[i  ];
    tmp2 += x1[i+1]*x2[i+1];
    tmp3 += x1[i+2]*x2[i+2];
    tmp4 += x1[i+3]*x2[i+3];
    tmp5 += x1[i+4]*x2[i+4];
    tmp6 += x1[i+5]*x2[i+5];
  }
  tmp = tmp1+tmp2+tmp3+tmp4+tmp5+tmp6;
/*...................................................................*/

/*...*/
#ifdef _MPICH_
  if(mpiVar.nPrcs>1){ 
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
    MPI_Allreduce(&tmp,&dot,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
  }
  else  
    dot = tmp;
#else
    dot = tmp;
#endif
/*...................................................................*/

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/

  return dot;
} 
/*********************************************************************/ 

/********************************************************************* 
 * DOTL8: produto interno entre dois vetores                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x1 - vetor x1                                                     * 
 * x2 - vetor x2                                                     * 
 * n  - numero de dimensoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE dotL8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n)
{
  INT i,resto;
  DOUBLE tmp1=0.e0,tmp2=0.e0,tmp3=0.e0,tmp4=0.e0;
  DOUBLE tmp5=0.e0,tmp6=0.e0,tmp7=0.e0,tmp8=0.e0;
  DOUBLE tmp,dot=0.e0;

  
/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/
  
/*...*/ 
  resto = n%8;
  if(resto == 1)
     tmp1= x1[0]*x2[0]; 
  else if(resto == 2)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1]; 
  else if(resto == 3)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]; 
  else if(resto == 3)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]; 
  else if(resto == 4)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2] + x1[3]*x2[3]; 
  else if(resto == 5)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1] 
         + x1[2]*x2[2] + x1[3]*x2[3]  
         + x1[4]*x2[4];               
  else if(resto == 6)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1] 
         + x1[2]*x2[2] + x1[3]*x2[3]  
         + x1[4]*x2[4] + x1[5]*x2[5]; 
  else if(resto == 7)
     tmp1= x1[0]*x2[0] + x1[1]*x2[1] 
         + x1[2]*x2[2] + x1[3]*x2[3]  
         + x1[4]*x2[4] + x1[5]*x2[5]  
         + x1[6]*x2[6];
  
  for(i=resto;i<n;i+=8){
    tmp1 += x1[i  ]*x2[i  ];
    tmp2 += x1[i+1]*x2[i+1];
    tmp3 += x1[i+2]*x2[i+2];
    tmp4 += x1[i+3]*x2[i+3];
    tmp5 += x1[i+4]*x2[i+4];
    tmp6 += x1[i+5]*x2[i+5];
    tmp7 += x1[i+6]*x2[i+6];
    tmp8 += x1[i+7]*x2[i+7];
  }
  dot = tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8;
/*...................................................................*/

/*...*/
#ifdef _MPICH_
  if(mpiVar.nPrcs>1){ 
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
    MPI_Allreduce(&tmp,&dot,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
  }
  else  
    dot = tmp;
#else
    dot = tmp;
#endif
/*...................................................................*/

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/

  return dot;
} 
/*********************************************************************/ 

/********************************************************************* 
 * DOTO2: produto interno entre dois vetores                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x1 - vetor x1                                                     * 
 * x2 - vetor x2                                                     * 
 * n  - numero de dimensoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE dotO2(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n)
{
  INT i,resto;
  DOUBLE tmp=0.e0;
  resto = n%2;
  DOUBLE dot=0.e0;

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/
  
/*...*/ 
  if(resto)
    tmp = x1[0]*x2[0];
    
  for(i=resto;i<n;i+=2)
    tmp += x1[i]*x2[i] + x1[i+1]*x2[i+1];
/*...................................................................*/

/*...*/
#ifdef _MPICH_
  if(mpiVar.nPrcs>1){ 
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
    MPI_Allreduce(&tmp,&dot,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
  }
  else  
    dot = tmp;
#else
    dot = tmp;
#endif
/*...................................................................*/

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/

  return dot;
} 
/*********************************************************************/ 

/********************************************************************* 
 * DOTO4: produto interno entre dois vetores                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x1 - vetor x1                                                     * 
 * x2 - vetor x2                                                     * 
 * n  - numero de dimensoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE dotO4(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n)
{
  INT i,resto;
  DOUBLE tmp=0.e0;
  DOUBLE dot=0.e0;
 
/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/
 
/*...*/ 
  resto = n%4;
  if(resto == 1)
     tmp = x1[0]*x2[0]; 
  else if(resto == 2)
     tmp = x1[0]*x2[0] + x1[1]*x2[1]; 
  else if(resto == 3)
     tmp = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]; 
    
  
  for(i=resto;i<n;i+=4)
    tmp += x1[i  ]*x2[i  ] + x1[i+1]*x2[i+1] 
         + x1[i+2]*x2[i+2] + x1[i+3]*x2[i+3];
/*...................................................................*/

/*...*/
#ifdef _MPICH_
  if(mpiVar.nPrcs>1){ 
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
    MPI_Allreduce(&tmp,&dot,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
  }
  else  
    dot = tmp;
#else
    dot = tmp;
#endif
/*...................................................................*/

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/

  return dot;
} 
/*********************************************************************/ 

/********************************************************************* 
 * DOTO6: produto interno entre dois vetores                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x1 - vetor x1                                                     * 
 * x2 - vetor x2                                                     * 
 * n  - numero de dimensoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE dotO6(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n)
{
  INT i,resto;
  DOUBLE tmp=0.e0;
  DOUBLE dot=0.e0;


/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/
  
/*...*/ 
  resto = n%6;
  if(resto == 1)
     tmp = x1[0]*x2[0]; 
  else if(resto == 2)
     tmp = x1[0]*x2[0] + x1[1]*x2[1]; 
  else if(resto == 3)
     tmp = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]; 
  else if(resto == 4)
     tmp = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2] + x1[3]*x2[3];
  else if(resto == 5)
     tmp = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2] 
         + x1[3]*x2[3] + x1[4]*x2[4];
    
  
  for(i=resto;i<n;i+=6)
    tmp += x1[i  ]*x2[i  ] + x1[i+1]*x2[i+1] 
         + x1[i+2]*x2[i+2] + x1[i+3]*x2[i+3] 
         + x1[i+4]*x2[i+4] + x1[i+5]*x2[i+5];
/*...................................................................*/
  
/*...*/
#ifdef _MPICH_
  if(mpiVar.nPrcs>1){ 
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
    MPI_Allreduce(&tmp,&dot,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
  }
  else  
    dot = tmp;
#else
    dot = tmp;
#endif
/*...................................................................*/

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/


  return dot;
} 
/*********************************************************************/ 

/********************************************************************* 
 * DOTO8: produto interno entre dois vetores                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x1 - vetor x1                                                     * 
 * x2 - vetor x2                                                     * 
 * n  - numero de dimensoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE dotO8(DOUBLE *restrict x1,DOUBLE *restrict x2,INT const n)
{
  INT i,resto;
  DOUBLE tmp=0.e0;
  DOUBLE dot=0.e0;

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/
  
/*...*/ 
  resto = n%8;
  if(resto == 1)
     tmp = x1[0]*x2[0]; 
  else if(resto == 2)
     tmp = x1[0]*x2[0] + x1[1]*x2[1]; 
  else if(resto == 3)
     tmp = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]; 
  else if(resto == 3)
     tmp = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]; 
  else if(resto == 4)
     tmp = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2] + x1[3]*x2[3]; 
  else if(resto == 5)
     tmp = x1[0]*x2[0] + x1[1]*x2[1] 
         + x1[2]*x2[2] + x1[3]*x2[3]  
         + x1[4]*x2[4];               
  else if(resto == 6)
     tmp = x1[0]*x2[0] + x1[1]*x2[1] 
         + x1[2]*x2[2] + x1[3]*x2[3]  
         + x1[4]*x2[4] + x1[5]*x2[5]; 
  else if(resto == 7)
     tmp = x1[0]*x2[0] + x1[1]*x2[1] 
         + x1[2]*x2[2] + x1[3]*x2[3]  
         + x1[4]*x2[4] + x1[5]*x2[5]  
         + x1[6]*x2[6];

  for(i=resto;i<n;i+=8)
    tmp += x1[i  ]*x2[i  ] + x1[i+1]*x2[i+1] 
         + x1[i+2]*x2[i+2] + x1[i+3]*x2[i+3]
         + x1[i+4]*x2[i+4] + x1[i+5]*x2[i+5]
         + x1[i+6]*x2[i+6] + x1[i+7]*x2[i+7];
/*...................................................................*/
 
/*...*/
#ifdef _MPICH_
  if(mpiVar.nPrcs>1){ 
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
    MPI_Allreduce(&tmp,&dot,1,MPI_DOUBLE,MPI_SUM,mpiVar.comm);
    tm.dotOverHeadMpi = getTimeC() - tm.dotOverHeadMpi;
  }
  else  
    dot = tmp;
#else
    dot = tmp;
#endif
/*...................................................................*/

/*...*/ 
  tm.dot = getTimeC() - tm.dot;
/*...................................................................*/

  return dot;
} 
/*********************************************************************/ 

/*==================================================================*/

/*======================== level 2 =================================*/

/*********************** CSRD SIMETRICO *****************************/

/********************************************************************* 
 * MATVECCSRDSYM:produto matriz vetor para uma matriz no formato CSRD* 
 * (y=Ax, A uma matriz simentrica)                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * al  -> vetor com os valores inferiores da matriz                  * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS: ja guarda os indiceis da para inferior da matriz             * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecCsrDSym(INT const neq           
                  ,INT *restrict ia   ,INT *restrict ja
                  ,DOUBLE *restrict al,DOUBLE *restrict ad
                  ,DOUBLE *restrict x ,DOUBLE *restrict y)
{
  INT    i,k,jak;
  DOUBLE xi,tmp,sAu;

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
 
  y[0] = ad[0]*x[0]; 
  for(i=1;i<neq;i++){
    xi  = x[i];
    tmp = ad[i]*xi;
    for(k=ia[i];k<ia[i+1];k++){
      jak = ja[k];
      sAu = al[k];
/*... produto da linha i pelo vetor x*/
      tmp    += sAu*x[jak];
/*... produto dos coef. da parte superior da matriz por x(i)*/
      y[jak] += sAu*xi;
    }
/*... armazena o resultado em y(i) */
    y[i] = tmp;
  }
  
/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
} 
/*********************************************************************/ 

/********************************************************************* 
 * MPIMATVECCSRDSYM:produto matriz vetor para uma matriz no formato  *
 * CSRD+CSR (y=Ax, A uma matriz simentrica)                          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * al  -> vetor com os valores inferiores da matriz                  * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS: ja guarda os indiceis da para inferior da matriz             * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void mpiMatVecCsrDSym(INT const nEq           
                     ,INT *restrict ia   ,INT *restrict ja
                     ,DOUBLE *restrict al,DOUBLE *restrict ad
                     ,DOUBLE *restrict x ,DOUBLE *restrict y
                     ,Interface *iNeq)
{
#ifdef _MPICH_
  INT    i,k,jak,nAd;
  DOUBLE xi,tmp,sAu;
  INT *iar,*jar;
  DOUBLE *ar;

  nAd = ia[nEq] - ia[0];
  iar = &ia[nEq+1];
  jar = &ja[nAd];
  ar  = &al[nAd];

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  
/*...*/
  tm.matVecOverHeadMpi        = getTimeC() - tm.matVecOverHeadMpi;
  comunicateNeq(iNeq,x);
  tm.matVecOverHeadMpi        = getTimeC() - tm.matVecOverHeadMpi;
/*...................................................................*/
      
/*...*/
  for(i=0;i<nEq;i++){
    xi  = x[i];
    tmp = ad[i]*xi;
    for(k=ia[i];k<ia[i+1];k++){
      jak = ja[k];
      sAu = al[k];
/*... produto da linha i pelo vetor x*/
      tmp    += sAu*x[jak];
/*... produto dos coef. da parte superior da matriz por x(i)*/
      y[jak] += sAu*xi;
    }
    
/*... parte retangular*/
    for(k=iar[i];k<iar[i+1];k++){
      jak  = jar[k];
      tmp += ar[k]*x[jak];
    } 
/*... armazena o resultado em y(i) */
    y[i] = tmp;
  }
/*...................................................................*/
 
/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
#endif
} 
/*********************************************************************/ 

/*********************** CSRD GERAL **********************************/

/********************************************************************* 
 * MATVECCSRD :produto matriz vetor para uma matriz no formato CSRD  * 
 * (y=Ax, A uma matriz nao simentrica)                               * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * a   -> vetor com os valores da matriz fora diagonal               * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecCsrD(INT const neq
               ,INT *restrict ia  ,INT *restrict ja
               ,DOUBLE *restrict a,DOUBLE *restrict ad
               ,DOUBLE *restrict x,DOUBLE *restrict y)
{

  INT i,j;
  DOUBLE tmp;

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/

  for(i=0;i<neq;i++){
    tmp = ad[i]*x[i];
    for(j=ia[i];j<ia[i+1];j++)
      tmp += a[j]*x[ja[j]];
    y[i] = tmp;
  }

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * MATVECCSRDI2:  produto matriz vetor para uma matriz generica no   *
 * formato csr com a diagonal principal retirada                     *
 * (y=Ax, A uma matriz geral)                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * a   -> vetor com os valores da matriz                             * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecCsrDI2(INT const neq           
                 ,INT *restrict ia  ,INT *restrict ja
                 ,DOUBLE *restrict a,DOUBLE *restrict ad
                 ,DOUBLE *restrict x,DOUBLE *restrict y)
{
  INT i;
  INT ia1,ia2,ja1;
  int resto,n;
  DOUBLE tmp;

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/

  for(i=0;i<neq;i++){
    tmp   = ad[i]*x[i];
    ia1   = ia[i  ];
    ia2   = ia[i+1];
    
    n     = ia2 - ia1;
    resto = n%2;
    
    if(resto)
      tmp += a[ia1]*x[ja[ia1]];

    for(ja1=ia1+resto;ja1<ia2;ja1+=2)
      tmp += a[  ja1]*x[ja[  ja1]]
           + a[ja1+1]*x[ja[ja1+1]];

    y[i] = tmp;
  }

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  
} 
/*********************************************************************/ 

/********************************************************************* 
 * MATVECCSRDI4:  produto matriz vetor para uma matriz generica no   *
 * formato csr com a diagonal principal retirada                     *
 * (y=Ax, A uma matriz geral)                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * a   -> vetor com os valores da matriz                             * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecCsrDI4(INT const neq           
                 ,INT *restrict ia  ,INT *restrict ja
                 ,DOUBLE *restrict a,DOUBLE *restrict ad
                 ,DOUBLE *restrict x,DOUBLE *restrict y)
{
  INT i;
  INT ia1,ia2,ja1;
  int resto,n;
  DOUBLE tmp;

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/

  for(i=0;i<neq;i++){
    tmp   = ad[i]*x[i];
    ia1   = ia[i  ];
    ia2   = ia[i+1];
    
    n     = ia2 - ia1;
    resto = n%4;
    
    if(resto == 3)
      tmp +=   a[ia1]*x[  ja[ia1]] 
           + a[ia1+1]*x[ja[ia1+1]] 
           + a[ia1+2]*x[ja[ia1+2]]; 
    else if(resto == 2)
      tmp +=   a[ia1]*x[  ja[ia1]] 
           + a[ia1+1]*x[ja[ia1+1]];
    else if(resto == 1)
      tmp += a[ia1]*x[ja[ia1]]; 

    for(ja1=ia1+resto;ja1<ia2;ja1+=4)
      tmp += a[  ja1]*x[ja[  ja1]]
           + a[ja1+1]*x[ja[ja1+1]] 
           + a[ja1+2]*x[ja[ja1+2]] 
           + a[ja1+3]*x[ja[ja1+3]];

    y[i] = tmp;
  }

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  
} 
/*********************************************************************/ 

/********************************************************************* 
 * MATVECCSRDI6:  produto matriz vetor para uma matriz generica no   *
 * formato csr com a diagonal principal retirada                     *
 * (y=Ax, A uma matriz geral)                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * a   -> vetor com os valores da matriz                             * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecCsrDI6(INT const neq           
                 ,INT *restrict ia  ,INT *restrict ja
                 ,DOUBLE *restrict a,DOUBLE *restrict ad
                 ,DOUBLE *restrict x,DOUBLE *restrict y)
{
  INT i;
  INT ia1,ia2,ja1;
  int resto,n;
  DOUBLE tmp;

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  for(i=0;i<neq;i++){
    tmp   = ad[i]*x[i];
    ia1   = ia[  i];
    ia2   = ia[i+1];
    
    n     = ia2 - ia1;
    resto = n%6;
    
    if(resto == 5)
      tmp +=   a[ia1]*x[  ja[ia1]] 
           + a[ia1+1]*x[ja[ia1+1]] 
           + a[ia1+2]*x[ja[ia1+2]]  
           + a[ia1+3]*x[ja[ia1+3]]  
           + a[ia1+4]*x[ja[ia1+4]]; 
    else if(resto == 4)
      tmp +=   a[ia1]*x[  ja[ia1]] 
           + a[ia1+1]*x[ja[ia1+1]] 
           + a[ia1+2]*x[ja[ia1+2]]  
           + a[ia1+3]*x[ja[ia1+3]]; 
    else if(resto == 3)
      tmp +=   a[ia1]*x[  ja[ia1]] 
           + a[ia1+1]*x[ja[ia1+1]] 
           + a[ia1+2]*x[ja[ia1+2]]; 
    else if(resto == 2)
      tmp +=   a[ia1]*x[  ja[ia1]] 
           + a[ia1+1]*x[ja[ia1+1]];
    else if(resto == 1)
      tmp += a[ia1]*x[ja[ia1]]; 

    for(ja1=ia1+resto;ja1<ia2;ja1+=6)
      tmp += a[  ja1]*x[ja[  ja1]]
           + a[ja1+1]*x[ja[ja1+1]] 
           + a[ja1+2]*x[ja[ja1+2]] 
           + a[ja1+3]*x[ja[ja1+3]] 
           + a[ja1+4]*x[ja[ja1+4]] 
           + a[ja1+5]*x[ja[ja1+5]];

    y[i] = tmp;
  }

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  
} 
/*********************************************************************/ 
 
/********************************************************************* 
 * MATVECCSRDO2:  produto matriz vetor para uma matriz generica no   *
 * formato csr com a diagonal principal retirada                     *
 * (y=Ax, A uma matriz geral)                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * a   -> vetor com os valores da matriz                             * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecCsrDO2(INT const neq           
                 ,INT *restrict ia  ,INT *restrict ja
                 ,DOUBLE *restrict a,DOUBLE *restrict ad
                 ,DOUBLE *restrict x,DOUBLE *restrict y)
{
  INT i;
  INT ia1,ia2,ia3,ja1;
  int resto;
  DOUBLE tmp1,tmp2;

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  resto = neq%2;
  
  if(resto){
    tmp1  = ad[0]*x[0];
    for(ja1=0;ja1<ia[1];ja1++)
      tmp1 += a[ja1]*x[ja[ja1]];
    y[0] = tmp1;
  }
  
  for(i=resto;i<neq;i+=2){
    tmp1  = ad[i  ]*x[  i];
    tmp2  = ad[i+1]*x[i+1];
    ia1   = ia[i  ];
    ia2   = ia[i+1];
    ia3   = ia[i+2];
    if(ia1 == ia2) goto linhai1;

/*linha i*/
    for(ja1=ia1;ja1<ia2;ja1++)
      tmp1 += a[ja1]*x[ja[ja1]];
/*...................................................................*/

/*linha i+1*/
linhai1:
    y[i] = tmp1;
    if(ia2 == ia3) goto linhai2;

    for(ja1=ia2;ja1<ia3;ja1++)
      tmp2 += a[ja1]*x[ja[ja1]];
/*...................................................................*/
linhai2:
    y[i+1] = tmp2;
  }

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  
} 
/*********************************************************************/ 

/********************************************************************* 
 * MATVECCRSDO4:  produto matriz vetor para uma matriz generica no   *
 * formato csr com a diagonal principal retirada                     *
 * (y=Ax, A uma matriz geral)                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * a   -> vetor com os valores da matriz                             * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecCsrDO4(INT const neq           
                 ,INT *restrict ia  ,INT *restrict ja
                 ,DOUBLE *restrict a,DOUBLE *restrict ad
                 ,DOUBLE *restrict x,DOUBLE *restrict y)
{
  INT i,j;
  INT ia1,ia2,ia3,ia4,ia5;
  int resto;
  DOUBLE tmp1,tmp2,tmp3,tmp4;

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  resto = neq %4;
  
  if(resto == 3){
    tmp1  = ad[0]*x[0];
    tmp2  = ad[1]*x[1];
    tmp3  = ad[2]*x[2];
    for(j=ia[0];j<ia[1];j++)
      tmp1 += a[j]*x[ja[j]];
    for(j=ia[1];j<ia[2];j++)
      tmp2 += a[j]*x[ja[j]];
    for(j=ia[2];j<ia[3];j++)
      tmp3 += a[j]*x[ja[j]];
    y[0] = tmp1;
    y[1] = tmp2;
    y[2] = tmp3;
  }

  else if(resto == 2){
    tmp1  = ad[0]*x[0];
    tmp2  = ad[1]*x[1];
    for(j=ia[0];j<ia[1];j++)
      tmp1 += a[j]*x[ja[j]];
    for(j=ia[1];j<ia[2];j++)
      tmp2 += a[j]*x[ja[j]];
    y[0] = tmp1;
    y[1] = tmp2;
  }
  
  else if(resto == 1){
    tmp1  = ad[0]*x[0];
    for(j=ia[0];j<ia[1];j++)
      tmp1 += a[j]*x[ja[j]];
    y[0] = tmp1;
  }
  
  for(i=resto;i<neq;i+=4){
    tmp1  = ad[i  ]*x[  i];
    tmp2  = ad[i+1]*x[i+1];
    tmp3  = ad[i+2]*x[i+2];
    tmp4  = ad[i+3]*x[i+3];
    ia1   = ia[i  ];
    ia2   = ia[i+1];
    ia3   = ia[i+2];
    ia4   = ia[i+3];
    ia5   = ia[i+4];
    if(ia1 == ia2) goto linhai1;

/*linha i*/
    for(j=ia1;j<ia2;j++)
      tmp1 += a[j]*x[ja[j]];
/*...................................................................*/

/*linha i+1*/
linhai1:
    y[i] = tmp1;
    if(ia2 == ia3) goto linhai2;

    for(j=ia2;j<ia3;j++)
      tmp2 += a[j]*x[ja[j]];
/*...................................................................*/

/*linha i+2*/
linhai2:
    y[i+1] = tmp2;
    if(ia4 == ia3) goto linhai3;

    for(j=ia3;j<ia4;j++)
      tmp3 += a[j]*x[ja[j]];
/*...................................................................*/

/*linha i+3*/
linhai3:
    y[i+2] = tmp3;
    if(ia5 == ia4) goto linhai4;

    for(j=ia4;j<ia5;j++)
      tmp4 += a[j]*x[ja[j]];
/*...................................................................*/
linhai4:
    y[i+3] = tmp4;
  }

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  
} 
/*********************************************************************/ 

/********************************************************************* 
 * MATVECCRSDO6:  produto matriz vetor para uma matriz generica no   *
 * formato csr com a diagonal principal retirada                     *
 * (y=Ax, A uma matriz geral)                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * a   -> vetor com os valores da matriz                             * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecCsrDO6(INT const neq           
                 ,INT *restrict ia  ,INT *restrict ja
                 ,DOUBLE *restrict a,DOUBLE *restrict ad
                 ,DOUBLE *restrict x,DOUBLE *restrict y)
{
  INT i,j;
  INT ia1,ia2,ia3,ia4,ia5,ia6,ia7;
  int resto;
  DOUBLE tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  resto = neq%6;
  
  if(resto == 5){
    tmp1  = ad[0]*x[0];
    tmp2  = ad[1]*x[1];
    tmp3  = ad[2]*x[2];
    tmp4  = ad[3]*x[3];
    tmp5  = ad[4]*x[4];
    for(j=ia[0];j<ia[1];j++)
      tmp1 += a[j]*x[ja[j]];
    for(j=ia[1];j<ia[2];j++)
      tmp2 += a[j]*x[ja[j]];
    for(j=ia[2];j<ia[3];j++)
      tmp3 += a[j]*x[ja[j]];
    for(j=ia[3];j<ia[4];j++)
      tmp4 += a[j]*x[ja[j]];
    for(j=ia[4];j<ia[5];j++)
      tmp5 += a[j]*x[ja[j]];
    y[0] = tmp1;
    y[1] = tmp2;
    y[2] = tmp3;
    y[3] = tmp4;
    y[4] = tmp5;
  }
  if(resto == 4){
    tmp1  = ad[0]*x[0];
    tmp2  = ad[1]*x[1];
    tmp3  = ad[2]*x[2];
    tmp4  = ad[3]*x[3];
    for(j=ia[0];j<ia[1];j++)
      tmp1 += a[j]*x[ja[j]];
    for(j=ia[1];j<ia[2];j++)
      tmp2 += a[j]*x[ja[j]];
    for(j=ia[2];j<ia[3];j++)
      tmp3 += a[j]*x[ja[j]];
    for(j=ia[3];j<ia[4];j++)
      tmp4 += a[j]*x[ja[j]];
    y[0] = tmp1;
    y[1] = tmp2;
    y[2] = tmp3;
    y[3] = tmp4;
  }
  if(resto == 3){
    tmp1  = ad[0]*x[0];
    tmp2  = ad[1]*x[1];
    tmp3  = ad[2]*x[2];
    for(j=ia[0];j<ia[1];j++)
      tmp1 += a[j]*x[ja[j]];
    for(j=ia[1];j<ia[2];j++)
      tmp2 += a[j]*x[ja[j]];
    for(j=ia[2];j<ia[3];j++)
      tmp3 += a[j]*x[ja[j]];
    y[0] = tmp1;
    y[1] = tmp2;
    y[2] = tmp3;
  }
  
  else if(resto == 2){
    tmp1  = ad[0]*x[0];
    tmp2  = ad[1]*x[1];
    for(j=ia[0];j<ia[1];j++)
      tmp1 += a[j]*x[ja[j]];
    for(j=ia[1];j<ia[2];j++)
      tmp2 += a[j]*x[ja[j]];
    y[0] = tmp1;
    y[1] = tmp2;
  }
  
  else if(resto == 1){
    tmp1  = ad[0]*x[0];
    for(j=ia[0];j<ia[1];j++)
      tmp1 += a[j]*x[ja[j]];
    y[0] = tmp1;
  }
  
for(i=resto;i<neq;i+=6){
    tmp1  = ad[i  ]*x[  i];
    tmp2  = ad[i+1]*x[i+1];
    tmp3  = ad[i+2]*x[i+2];
    tmp4  = ad[i+3]*x[i+3];
    tmp5  = ad[i+4]*x[i+4];
    tmp6  = ad[i+5]*x[i+5];
    ia1   = ia[i  ];
    ia2   = ia[i+1];
    ia3   = ia[i+2];
    ia4   = ia[i+3];
    ia5   = ia[i+4];
    ia6   = ia[i+5];
    ia7   = ia[i+6];
    if(ia1 == ia2) goto linhai1;

/*linha i*/
    for(j=ia1;j<ia2;j++)
      tmp1 += a[j]*x[ja[j]];
/*...................................................................*/

linhai1:
    y[i] = tmp1;

/*linha i+1*/
    if(ia2 == ia3) goto linhai2;

    for(j=ia2;j<ia3;j++)
      tmp2 += a[j]*x[ja[j]];
/*...................................................................*/

linhai2:
    y[i+1] = tmp2;
/*linha i+2*/
    if(ia4 == ia3) goto linhai3;

    for(j=ia3;j<ia4;j++)
      tmp3 += a[j]*x[ja[j]];
/*...................................................................*/

linhai3:
    y[i+2] = tmp3;
/*linha i+3*/
    if(ia5 == ia4) goto linhai4;

    for(j=ia4;j<ia5;j++)
      tmp4 += a[j]*x[ja[j]];
/*...................................................................*/

linhai4:
    y[i+3] = tmp4;
/*linha i+4*/
    if(ia6 == ia5) goto linhai5;

    for(j=ia5;j<ia6;j++)
      tmp5 += a[j]*x[ja[j]];
/*...................................................................*/

linhai5:
    y[i+4] = tmp5;
/*linha i+5*/
    if(ia7 == ia6) goto linhai6;

    for(j=ia6;j<ia7;j++)
      tmp6 += a[j]*x[ja[j]];
/*...................................................................*/
linhai6:
    y[i+5] = tmp6;
  }

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  
} 
/*********************************************************************/ 

/********************************************************************* 
 * MATVECCRSDO2I2: produto matriz vetor para uma matriz generica no  *
 * formato csr com a diagonal principal retirada                     *
 * (y=Ax, A uma matriz geral)                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * a   -> vetor com os valores da matriz                             * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecCsrDO2I2(INT const neq           
                   ,INT *restrict ia  ,INT *restrict ja
                   ,DOUBLE *restrict a,DOUBLE *restrict ad
                   ,DOUBLE *restrict x,DOUBLE *restrict y)
{
  INT i,j;
  INT ia1,ia2,ia3;
  int resto1,resto2,n;
  DOUBLE tmp1,tmp2;
  
/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  resto1 = neq % 2;
  
  if(resto1){
    tmp1  = ad[0]*x[0];
    for(j=ia[0];j<ia[1];j++)
      tmp1 += a[j]*x[ja[j]];
    y[0] = tmp1;
  }
  
  for(i=resto1;i<neq;i+=2){
    tmp1  = ad[i  ]*x[  i];
    tmp2  = ad[i+1]*x[i+1];
    ia1   = ia[i  ];
    ia2   = ia[i+1];
    ia3   = ia[i+2];
    if(ia1 == ia2) goto linhai1;
    n      = ia2 - ia1;
    resto2 = n%2;
    if(resto2)
      tmp1 += a[ia1]*x[ja[ia1]];

/*linha i*/
    for(j=ia1+resto2;j<ia2;j+=2)
      tmp1 +=a[  j]*x[ja[  j]] 
           + a[j+1]*x[ja[j+1]];
/*...................................................................*/

linhai1:
    y[i] = tmp1;
/*linha i+1*/
    if(ia2 == ia3) goto linhai2;
    n      = ia3 - ia2;
    resto2 = n%2;
    if(resto2)
      tmp2 += a[ia2]*x[ja[ia2]];

    for(j=ia2+resto2;j<ia3;j+=2)
      tmp2 += a[  j]*x[  ja[j]]
            + a[j+1]*x[ja[j+1]];
/*...................................................................*/
linhai2:
    y[i+1] = tmp2;
  }

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  
} 
/*********************************************************************/ 

/************************ ELLPACK GERAL ******************************/

/********************************************************************* 
 * MATVECELLPACK : produto matriz vetor para uma matriz generica no  *
 * formato ELLPACK                                                   *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * a   -> vetor com os valores da matriz                             * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecEllPack(INT const nEq           
                  ,INT *restrict ia  ,INT *restrict ja
                  ,DOUBLE *restrict a,DOUBLE *restrict ad
                  ,DOUBLE *restrict x,DOUBLE *restrict y)
{
 
  INT i,maxLineNzr=ia[0];
  DOUBLE ti;

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
 
  if( maxLineNzr == 6){
    for(i = 0;i<nEq;i++){
      ti   = ad[i]*x[i]; 
      ti  += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti  += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti  += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      ti  += MAT2D(i,3,a,maxLineNzr)*x[MAT2D(i,3,ja,maxLineNzr)];
      ti  += MAT2D(i,4,a,maxLineNzr)*x[MAT2D(i,4,ja,maxLineNzr)];
      ti  += MAT2D(i,5,a,maxLineNzr)*x[MAT2D(i,5,ja,maxLineNzr)];
      y[i] = ti;   
    }
  }
  
  else if( maxLineNzr == 4){
    for(i = 0;i<nEq;i++){
      ti    = ad[i]*x[i]; 
      ti   += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti   += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti   += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      ti   += MAT2D(i,3,a,maxLineNzr)*x[MAT2D(i,3,ja,maxLineNzr)];
      y[i] = ti;                 
    }
  }
  
  else if( maxLineNzr == 3){
    for(i = 0;i<nEq;i++){
      ti   = ad[i]*x[i]; 
      ti  += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti  += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti  += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      y[i] = ti;                 
    }
  
  }
  
/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  
} 
/*********************************************************************/ 

/********************************************************************* 
 * MATVECELLPACKO2 : produto matriz vetor para uma matriz generica no*
 * formato ELLPACK                                                   *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * a   -> vetor com os valores da matriz                             * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecEllPackO2(INT const nEq           
                    ,INT *restrict ia  ,INT *restrict ja
                    ,DOUBLE *restrict a,DOUBLE *restrict ad
                    ,DOUBLE *restrict x,DOUBLE *restrict y)
{
 
  INT i,ii,maxLineNzr=ia[0];
  int resto;
  DOUBLE ti,tii;

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
   
  if( maxLineNzr == 6){
    resto = nEq%2;
    if(resto){
      ti   = ad[0]*x[0]; 
      ti  += MAT2D(0,0,a,maxLineNzr)*x[MAT2D(0,0,ja,maxLineNzr)];
      ti  += MAT2D(0,1,a,maxLineNzr)*x[MAT2D(0,1,ja,maxLineNzr)];
      ti  += MAT2D(0,2,a,maxLineNzr)*x[MAT2D(0,2,ja,maxLineNzr)];
      ti  += MAT2D(0,3,a,maxLineNzr)*x[MAT2D(0,3,ja,maxLineNzr)];
      ti  += MAT2D(0,4,a,maxLineNzr)*x[MAT2D(0,4,ja,maxLineNzr)];
      ti  += MAT2D(0,5,a,maxLineNzr)*x[MAT2D(0,5,ja,maxLineNzr)];
      y[0] = ti;                                        
    }
    for(i = resto;i<nEq;i+=2){
      ii   = i + 1;
/*...*/
      ti   = ad[i]*x[i]; 
      ti  += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti  += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti  += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      ti  += MAT2D(i,3,a,maxLineNzr)*x[MAT2D(i,3,ja,maxLineNzr)];
      ti  += MAT2D(i,4,a,maxLineNzr)*x[MAT2D(i,4,ja,maxLineNzr)];
      ti  += MAT2D(i,5,a,maxLineNzr)*x[MAT2D(i,5,ja,maxLineNzr)];
      y[i] = ti;                                        
/*...*/
      tii   = ad[ii]*x[ii]; 
      tii  += MAT2D(ii,0,a,maxLineNzr)*x[MAT2D(ii,0,ja,maxLineNzr)];
      tii  += MAT2D(ii,1,a,maxLineNzr)*x[MAT2D(ii,1,ja,maxLineNzr)];
      tii  += MAT2D(ii,2,a,maxLineNzr)*x[MAT2D(ii,2,ja,maxLineNzr)];
      tii  += MAT2D(ii,3,a,maxLineNzr)*x[MAT2D(ii,3,ja,maxLineNzr)];
      tii  += MAT2D(ii,4,a,maxLineNzr)*x[MAT2D(ii,4,ja,maxLineNzr)];
      tii  += MAT2D(ii,5,a,maxLineNzr)*x[MAT2D(ii,5,ja,maxLineNzr)];
      y[ii] = tii;                                        
    }
  }
  else if( maxLineNzr == 4){
    resto = nEq%2;
    if(resto){
      ti   = ad[0]*x[0]; 
      ti  += MAT2D(0,0,a,maxLineNzr)*x[MAT2D(0,0,ja,maxLineNzr)];
      ti  += MAT2D(0,1,a,maxLineNzr)*x[MAT2D(0,1,ja,maxLineNzr)];
      ti  += MAT2D(0,2,a,maxLineNzr)*x[MAT2D(0,2,ja,maxLineNzr)];
      ti  += MAT2D(0,3,a,maxLineNzr)*x[MAT2D(0,3,ja,maxLineNzr)];
      y[0] = ti;                                        
    }
    for(i = resto;i<nEq;i+=2){
      ii   = i + 1;
/*...*/
      ti   = ad[i]*x[i]; 
      ti  += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti  += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti  += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      ti  += MAT2D(i,3,a,maxLineNzr)*x[MAT2D(i,3,ja,maxLineNzr)];
      y[i] = ti;                                        
/*...*/
      tii   = ad[ii]*x[ii]; 
      tii  += MAT2D(ii,0,a,maxLineNzr)*x[MAT2D(ii,0,ja,maxLineNzr)];
      tii  += MAT2D(ii,1,a,maxLineNzr)*x[MAT2D(ii,1,ja,maxLineNzr)];
      tii  += MAT2D(ii,2,a,maxLineNzr)*x[MAT2D(ii,2,ja,maxLineNzr)];
      tii  += MAT2D(ii,3,a,maxLineNzr)*x[MAT2D(ii,3,ja,maxLineNzr)];
      y[ii] = tii;                                        
    }
  }
  
  else if( maxLineNzr == 3){
    resto = nEq%2;
    if(resto){
      ti   = ad[0]*x[0]; 
      ti  += MAT2D(0,0,a,maxLineNzr)*x[MAT2D(0,0,ja,maxLineNzr)];
      ti  += MAT2D(0,1,a,maxLineNzr)*x[MAT2D(0,1,ja,maxLineNzr)];
      ti  += MAT2D(0,2,a,maxLineNzr)*x[MAT2D(0,2,ja,maxLineNzr)];
      y[0] = ti;                                        
    }
    for(i = resto;i<nEq;i+=2){
      ii   = i + 1;
/*...*/
      ti   = ad[i]*x[i]; 
      ti  += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti  += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti  += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      y[i] = ti;                                        
/*...*/
      tii   = ad[ii]*x[ii]; 
      tii  += MAT2D(ii,0,a,maxLineNzr)*x[MAT2D(ii,0,ja,maxLineNzr)];
      tii  += MAT2D(ii,1,a,maxLineNzr)*x[MAT2D(ii,1,ja,maxLineNzr)];
      tii  += MAT2D(ii,2,a,maxLineNzr)*x[MAT2D(ii,2,ja,maxLineNzr)];
      y[ii] = tii;                                        
    }
  }
/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  
} 
/*********************************************************************/ 

/********************************************************************* 
 * MATVECELLPACKO4 : produto matriz vetor para uma matriz generica no*
 * formato ELLPACK                                                   *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * neq -> numero de equacoes                                         * 
 * ia  -> vetor csr                                                  * 
 * ja  -> vetor csr                                                  * 
 * a   -> vetor com os valores da matriz                             * 
 * ad  -> vetor com os valores da diagonal principal da matriz       * 
 * x   -> vetor a ser multiplicado                                   * 
 * y   -> indefinido                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * y   -> vetor com o resultado da multiplicacao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void matVecEllPackO4(INT const nEq           
                    ,INT *restrict ia  ,INT *restrict ja
                    ,DOUBLE *restrict a,DOUBLE *restrict ad
                    ,DOUBLE *restrict x,DOUBLE *restrict y)
{
 
  INT i,ii,iii,iv,maxLineNzr=ia[0];
  int resto;
  DOUBLE ti,tii,tiii,tiv;

/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
   
  if( maxLineNzr == 6){
    resto = nEq%4;
    for(i = 0;i<resto;i++){
      ti   = ad[i]*x[i]; 
      ti  += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti  += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti  += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      ti  += MAT2D(i,3,a,maxLineNzr)*x[MAT2D(i,3,ja,maxLineNzr)];
      ti  += MAT2D(i,4,a,maxLineNzr)*x[MAT2D(i,4,ja,maxLineNzr)];
      ti  += MAT2D(i,5,a,maxLineNzr)*x[MAT2D(i,5,ja,maxLineNzr)];
      y[i] = ti;                                        
    }
    for(i = resto;i<nEq;i+=4){
      ii   = i + 1;
      iii  = i + 2;
      iv   = i + 3;
/*...*/
      ti   = ad[i]*x[i]; 
      ti  += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti  += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti  += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      ti  += MAT2D(i,3,a,maxLineNzr)*x[MAT2D(i,3,ja,maxLineNzr)];
      ti  += MAT2D(i,4,a,maxLineNzr)*x[MAT2D(i,4,ja,maxLineNzr)];
      ti  += MAT2D(i,5,a,maxLineNzr)*x[MAT2D(i,5,ja,maxLineNzr)];
      y[i] = ti;                                        
/*...*/
      tii   = ad[ii]*x[ii]; 
      tii  += MAT2D(ii,0,a,maxLineNzr)*x[MAT2D(ii,0,ja,maxLineNzr)];
      tii  += MAT2D(ii,1,a,maxLineNzr)*x[MAT2D(ii,1,ja,maxLineNzr)];
      tii  += MAT2D(ii,2,a,maxLineNzr)*x[MAT2D(ii,2,ja,maxLineNzr)];
      tii  += MAT2D(ii,3,a,maxLineNzr)*x[MAT2D(ii,3,ja,maxLineNzr)];
      tii  += MAT2D(ii,4,a,maxLineNzr)*x[MAT2D(ii,4,ja,maxLineNzr)];
      tii  += MAT2D(ii,5,a,maxLineNzr)*x[MAT2D(ii,5,ja,maxLineNzr)];
      y[ii] = tii;                                        
/*...*/
      tiii  = ad[iii]*x[iii]; 
      tiii += MAT2D(iii,0,a,maxLineNzr)*x[MAT2D(iii,0,ja,maxLineNzr)];
      tiii += MAT2D(iii,1,a,maxLineNzr)*x[MAT2D(iii,1,ja,maxLineNzr)];
      tiii += MAT2D(iii,2,a,maxLineNzr)*x[MAT2D(iii,2,ja,maxLineNzr)];
      tiii += MAT2D(iii,3,a,maxLineNzr)*x[MAT2D(iii,3,ja,maxLineNzr)];
      tiii += MAT2D(iii,4,a,maxLineNzr)*x[MAT2D(iii,4,ja,maxLineNzr)];
      tiii += MAT2D(iii,5,a,maxLineNzr)*x[MAT2D(iii,5,ja,maxLineNzr)];
      y[iii] = tiii;                                        
/*...*/
      tiv   = ad[iv]*x[iv]; 
      tiv  += MAT2D(iv,0,a,maxLineNzr)*x[MAT2D(iv,0,ja,maxLineNzr)];
      tiv  += MAT2D(iv,1,a,maxLineNzr)*x[MAT2D(iv,1,ja,maxLineNzr)];
      tiv  += MAT2D(iv,2,a,maxLineNzr)*x[MAT2D(iv,2,ja,maxLineNzr)];
      tiv  += MAT2D(iv,3,a,maxLineNzr)*x[MAT2D(iv,3,ja,maxLineNzr)];
      tiv  += MAT2D(iv,4,a,maxLineNzr)*x[MAT2D(iv,4,ja,maxLineNzr)];
      tiv  += MAT2D(iv,5,a,maxLineNzr)*x[MAT2D(iv,5,ja,maxLineNzr)];
      y[iv] = tiv;                                        
    }
  }
/*... */
  else if( maxLineNzr == 4){
    resto = nEq%4;
    for(i = 0;i<resto;i++){
      ti   = ad[i]*x[i]; 
      ti  += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti  += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti  += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      ti  += MAT2D(i,3,a,maxLineNzr)*x[MAT2D(i,3,ja,maxLineNzr)];
      y[i] = ti;                                        
    }
    for(i = resto;i<nEq;i+=4){
      ii   = i + 1;
      iii  = i + 2;
      iv   = i + 3;
/*...*/
      ti   = ad[i]*x[i]; 
      ti  += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti  += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti  += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      ti  += MAT2D(i,3,a,maxLineNzr)*x[MAT2D(i,3,ja,maxLineNzr)];
      y[i] = ti;                                        
/*...*/
      tii   = ad[ii]*x[ii]; 
      tii  += MAT2D(ii,0,a,maxLineNzr)*x[MAT2D(ii,0,ja,maxLineNzr)];
      tii  += MAT2D(ii,1,a,maxLineNzr)*x[MAT2D(ii,1,ja,maxLineNzr)];
      tii  += MAT2D(ii,2,a,maxLineNzr)*x[MAT2D(ii,2,ja,maxLineNzr)];
      tii  += MAT2D(ii,3,a,maxLineNzr)*x[MAT2D(ii,3,ja,maxLineNzr)];
      y[ii] = tii;                                        
/*...*/
      tiii  = ad[iii]*x[iii]; 
      tiii += MAT2D(iii,0,a,maxLineNzr)*x[MAT2D(iii,0,ja,maxLineNzr)];
      tiii += MAT2D(iii,1,a,maxLineNzr)*x[MAT2D(iii,1,ja,maxLineNzr)];
      tiii += MAT2D(iii,2,a,maxLineNzr)*x[MAT2D(iii,2,ja,maxLineNzr)];
      tiii += MAT2D(iii,3,a,maxLineNzr)*x[MAT2D(iii,3,ja,maxLineNzr)];
      y[iii] = tiii;                                        
/*...*/
      tiv   = ad[iv]*x[iv]; 
      tiv  += MAT2D(iv,0,a,maxLineNzr)*x[MAT2D(iv,0,ja,maxLineNzr)];
      tiv  += MAT2D(iv,1,a,maxLineNzr)*x[MAT2D(iv,1,ja,maxLineNzr)];
      tiv  += MAT2D(iv,2,a,maxLineNzr)*x[MAT2D(iv,2,ja,maxLineNzr)];
      tiv  += MAT2D(iv,3,a,maxLineNzr)*x[MAT2D(iv,3,ja,maxLineNzr)];
      y[iv] = tiv;                                        
    }
  }
/*... */
  else if( maxLineNzr == 3){
    resto = nEq%4;
    for(i = 0;i<resto;i++){
      ti   = ad[i]*x[i]; 
      ti  += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti  += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti  += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      y[i] = ti;                                        
    }
    for(i = resto;i<nEq;i+=4){
      ii   = i + 1;
      iii  = i + 2;
      iv   = i + 3;
/*...*/
      ti   = ad[i]*x[i]; 
      ti  += MAT2D(i,0,a,maxLineNzr)*x[MAT2D(i,0,ja,maxLineNzr)];
      ti  += MAT2D(i,1,a,maxLineNzr)*x[MAT2D(i,1,ja,maxLineNzr)];
      ti  += MAT2D(i,2,a,maxLineNzr)*x[MAT2D(i,2,ja,maxLineNzr)];
      y[i] = ti;                                        
/*...*/
      tii   = ad[ii]*x[ii]; 
      tii  += MAT2D(ii,0,a,maxLineNzr)*x[MAT2D(ii,0,ja,maxLineNzr)];
      tii  += MAT2D(ii,1,a,maxLineNzr)*x[MAT2D(ii,1,ja,maxLineNzr)];
      tii  += MAT2D(ii,2,a,maxLineNzr)*x[MAT2D(ii,2,ja,maxLineNzr)];
      y[ii] = tii;                                        
/*...*/
      tiii  = ad[iii]*x[iii]; 
      tiii += MAT2D(iii,0,a,maxLineNzr)*x[MAT2D(iii,0,ja,maxLineNzr)];
      tiii += MAT2D(iii,1,a,maxLineNzr)*x[MAT2D(iii,1,ja,maxLineNzr)];
      tiii += MAT2D(iii,2,a,maxLineNzr)*x[MAT2D(iii,2,ja,maxLineNzr)];
      y[iii] = tiii;                                        
/*...*/
      tiv   = ad[iv]*x[iv]; 
      tiv  += MAT2D(iv,0,a,maxLineNzr)*x[MAT2D(iv,0,ja,maxLineNzr)];
      tiv  += MAT2D(iv,1,a,maxLineNzr)*x[MAT2D(iv,1,ja,maxLineNzr)];
      tiv  += MAT2D(iv,2,a,maxLineNzr)*x[MAT2D(iv,2,ja,maxLineNzr)];
      y[iv] = tiv;                                        
    }
  }
/*...*/ 
  tm.matVecSparse             = getTimeC() - tm.matVecSparse;
/*...................................................................*/
  
} 
/*********************************************************************/ 

 
/*==================================================================*/

/*======================== level 3 =================================*/
/*==================================================================*/
