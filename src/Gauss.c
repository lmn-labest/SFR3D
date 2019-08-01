#include<Gauss.h>

/********************************************************************* 
 * fatLDMT : fatora da matriz A em L, D e Mt                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa                       * 
 * d(neq)      -> nao definido                                       * 
 * r(neq)      -> vetor auxiliar (nao definido)                      * 
 * w(neq)      -> vetor auxiliar (nao definido)                      * 
 * nEq         -> numero de equacoes                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * a           -> atualizado com L e Mt                              * 
 * d           -> atualizado com D                                   * 
 *-------------------------------------------------------------------* 
 * OBS: matriz cheia                                                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void fatLDMt(DOUBLE *RESTRICT a,DOUBLE *RESTRICT d
            ,DOUBLE *RESTRICT r,DOUBLE *RESTRICT w
            ,INT const nEq){
  
  int row;
  int i,p;
  DOUBLE zero = 1.0e0-15;
  DOUBLE tmp1,tmp2;
  

  for(row=0;row<nEq;row++){
    for(p=0;p<row;p++){
      r[p] = d[p]*MAT2D(p,row,a,nEq);
      w[p] = d[p]*MAT2D(row,p,a,nEq);
    } 
    
    tmp1=0.0e0;
    for(p=0;p<row;p++){
      tmp1 += MAT2D(row,p,a,nEq)*r[p];
    }
    
    d[row] = MAT2D(row,row,a,nEq) - tmp1;
    
    if( fabs(d[row]) < zero  ){
      printf("Erro: termo nulo na diagonal da matriz D\n");
      exit(EXIT_FAILURE);
      break;
    }
    for(i=row+1;i<nEq;i++){
      tmp1=0.0e0;
      tmp2=0.0e0;
      for(p=0;p<row;p++){
        tmp1 += MAT2D(i,p,a,nEq)*r[p];
        tmp2 += MAT2D(p,i,a,nEq)*w[p];
      }
      MAT2D(i,row,a,nEq) = (MAT2D(i,row,a,nEq) - tmp1)/d[row];
      MAT2D(row,i,a,nEq) = (MAT2D(row,i,a,nEq) - tmp2)/d[row];
    }
     
  }
  
  for(p=0;p<nEq;p++)
    MAT2D(p,p,a,nEq) = 0.0e0;;


}
/*********************************************************************/ 

/********************************************************************* 
 * fatLDLT : fatora da matriz A em L, D e Lt ( A - simetrica)        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa                       * 
 * d(neq)      -> nao definido                                       * 
 * r(neq)      -> vetor auxiliar (nao definido)                      * 
 * nEq         -> numero de equacoes                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * a           -> atualizado com L                                   * 
 * d           -> atualizado com D                                   * 
 *-------------------------------------------------------------------* 
 * OBS: Parte superio da matriz A fica inalterada                    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void fatLDLt(DOUBLE *RESTRICT a,DOUBLE *RESTRICT d
            ,DOUBLE *RESTRICT r,INT const nEq){
  
  int row;
  int i,p;
  DOUBLE zero = 1.0e0-15;
  DOUBLE tmp;

  for(row=0;row<nEq;row++){
    for(p=0;p<row;p++)
      r[p] = d[p]*MAT2D(row,p,a,nEq);
    
    tmp=0.0e0;
    for(p=0;p<row;p++){
      tmp += MAT2D(row,p,a,nEq)*r[p];
    }
    d[row] = MAT2D(row,row,a,nEq) - tmp;
     
    if( fabs(d[row]) < zero  ){
      printf("Erro: termo nulo na diagonal da matriz D\n");
      exit(EXIT_FAILURE);
      break;
    }
    
    for(i=row+1;i<nEq;i++){
      tmp=0.0e0;
      for(p=0;p<row;p++)
        tmp += MAT2D(i,p,a,nEq)*r[p];
      MAT2D(i,row,a,nEq) = (MAT2D(i,row,a,nEq) - tmp)/d[row];
    } 
  }
  
  for(p=0;p<nEq;p++)
    MAT2D(p,p,a,nEq) = 0.0e0;;
  
}
/*********************************************************************/ 

/********************************************************************* 
 * fatGGT:   fatora da matriz A em G e Gt (Cholesky) (A - simetrica  * 
 * definida positiva - SDP)                                          *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa                       * 
 * nEq         -> numero de equacoes                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * a           -> atualizado com G                                   * 
 *-------------------------------------------------------------------* 
 * OBS: Parte superio da matriz A fica inalterada                    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void fatGGt(DOUBLE *RESTRICT a,INT const nEq){
  
  int row;
  int i,p;
  DOUBLE tmp;

  for(row=0;row<nEq;row++){
    
    tmp=0.0e0;
    for(p=0;p<row;p++){
      tmp += MAT2D(row,p,a,nEq)*MAT2D(row,p,a,nEq);
    }
    MAT2D(row,row,a,nEq) = sqrt(MAT2D(row,row,a,nEq) - tmp);
     
    for(i=row+1;i<nEq;i++){
      tmp=0.0e0;
      for(p=0;p<row;p++)
        tmp += MAT2D(i,p,a,nEq)*MAT2D(row,p,a,nEq);
      MAT2D(i,row,a,nEq) =  (MAT2D(i,row,a,nEq) - tmp)
                            /MAT2D(row,row,a,nEq);
    } 
  }
  
}
/*********************************************************************/ 

/********************************************************************* 
 * iZeroFatGGT:  fatoracao incompleta de (Cholesky) (A - simetrica   * 
 * definida positiva - SDP)                                          *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa                       * 
 * nEq         -> numero de equacoes                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * a           -> atualizado com G                                   * 
 *-------------------------------------------------------------------* 
 * OBS: Parte superio da matriz A fica inalterada                    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void ic0(DOUBLE *RESTRICT a,INT const nEq){
  
  int j,k,i;
  DOUBLE tmp1,tmp2,tmp3;

  for(j=0;j<nEq;j++){
    
    MAT2D(j,j,a,nEq) = sqrt(MAT2D(j,j,a,nEq));
    
    for(k=0;k<j;k++){
      tmp1 = MAT2D(j,k,a,nEq);
      if( tmp1 != 0.0) 
        for(i=j+1;i<nEq;i++){
          tmp2 = MAT2D(i,k,a,nEq);
          tmp3 = MAT2D(i,j,a,nEq);
          if(tmp2 != 0.0 && tmp3 != 0.0)
            MAT2D(i,j,a,nEq) = tmp3 - tmp2*tmp1;
        }
    }

    for(i=j+1;i<nEq;i++){
      tmp1   = MAT2D(i,j,a,nEq);
      if( tmp1 != 0.0) {
        MAT2D(i,j,a,nEq)  = tmp1 / MAT2D(j,j,a,nEq);
        MAT2D(i,i,a,nEq) -= MAT2D(i,j,a,nEq)*MAT2D(i,j,a,nEq);
      }
    }
    
  }  
}
/*********************************************************************/ 

/********************************************************************* 
 * fatLU :   fatora da matriz A em L e U                             * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa                       * 
 * nEq         -> numero de equacoes                                 * 
 * code        ->  4 - LUKIJ                                         * 
 *             ->  5 - LUIKJ                                         * 
 *             ->  6 - LUKJI                                         * 
 *             ->  7 - LUJKI                                         * 
 *             ->  8 - LUIJK                                         * 
 *             ->  9 - DOOLITTLELC                                   *  
 *             -> 10 - DOOLITTLECL                                   *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * a           -> atualizado com L e U                               * 
 *-------------------------------------------------------------------* 
 * OBS: matriz cheia                                                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void fatLU(DOUBLE *RESTRICT a,INT const nEq,short const code){
  
  int i;
  int j;
  int k;
  int ii;
  DOUBLE tmp,conv,diag;
  DOUBLE zero = 1.0e-15;

  tmp = 0.0e0;
  for(i=0;i<nEq;i++)
    tmp += MAT2D(i,i,a,nEq)*MAT2D(i,i,a,nEq);
  conv = sqrt(tmp)*zero;

  switch(code){
/*... versao KIJ - (escalar vector)*/
    case LUKIJ:
      for(k=0;k<nEq-1;k++){
        diag  = MAT2D(k,k,a,nEq);
        if( fabs(diag) <= conv ){ 
          printf("Erro: termo nulo na diagonal da matriz D\n");
          exit(EXIT_FAILURE);
        }
        ii = k + 1;
        for(i=ii;i<nEq;i++){
          tmp = MAT2D(i,k,a,nEq) / diag; 
          MAT2D(i,k,a,nEq) = tmp;
          addVector(1.0e0  ,&MAT2D(i,ii,a,nEq)
                   ,-tmp   ,&MAT2D(k,ii,a,nEq)
                   ,nEq-ii ,&MAT2D(i,ii,a,nEq));
/*        for(j=ii;j<nEq;j++)
              MAT2D(i,j,a,nEq) -= tmp*MAT2D(k,j,a,nEq);*/
        }
      }
    break;
/*...................................................................*/

/*... versao IKJ - (escalar vector)*/
    case LUIKJ:
      for(i=1;i<nEq;i++)
        for(k=0;k<i;k++){
          diag  = MAT2D(k,k,a,nEq);
          if( fabs(diag) <= conv ){ 
            printf("Erro: termo nulo na diagonal da matriz D\n");
            exit(EXIT_FAILURE);
          }
          tmp = MAT2D(i,k,a,nEq) / diag; 
          MAT2D(i,k,a,nEq) = tmp;
          ii = k + 1;
          addVector(1.0e0  ,&MAT2D(i,ii,a,nEq)
                   ,-tmp   ,&MAT2D(k,ii,a,nEq)
                   ,nEq-ii ,&MAT2D(i,ii,a,nEq));
/*        for(j=ii;j<nEq;j++)
            MAT2D(i,j,a,nEq) -= tmp*MAT2D(k,j,a,nEq);*/
        }
      break;
/*...................................................................*/

/*... versao KJI - (escalar vector)*/
    case LUKJI:
      for(k=0;k<nEq-1;k++){
        diag  = MAT2D(k,k,a,nEq);
        if( fabs(diag) <= conv ){ 
          printf("Erro: termo nulo na diagonal da matriz D\n");
          exit(EXIT_FAILURE);
        }
        ii = k + 1;
        for(i=ii;i<nEq;i++)
          MAT2D(i,k,a,nEq) /= diag; 
        for(j=ii;j<nEq;j++){
          for(i=ii;i<nEq;i++)
            MAT2D(i,j,a,nEq) -= MAT2D(i,k,a,nEq)*MAT2D(k,j,a,nEq);
        }
      }
      break;
/*...................................................................*/

/*... versao JKI */
    case LUJKI:
      for(j=0;j<nEq;j++){
        for(k=0;k<j;k++){
          for(i=k+1;i<nEq;i++)
            MAT2D(i,j,a,nEq) -= MAT2D(i,k,a,nEq)*MAT2D(k,j,a,nEq);
        }
        for(i=j+1;i<nEq;i++){
          diag  = MAT2D(j,j,a,nEq);
          if( fabs(diag) <= conv ){ 
            printf("Erro: termo nulo na diagonal da matriz D\n");
            exit(EXIT_FAILURE);
          }
          MAT2D(i,j,a,nEq) /= diag; 
        }
      }
      break;
/*...................................................................*/

/*... versao IJK (dot - version - DooLittte - compact form)*/
/*   varredura                                                     */
/* | i------f | i - inicio                                         */
/* | i------f | f - fim                                            */
/* | i------f |                                                    */
/* | i------f |                                                    */
/*                                                                 */
    case DOOLITTLEL:
      for(i=0;i<nEq;i++){
/*... parte superior*/
        for(j=i;j<nEq;j++){
          for(k=0;k<i;k++)
            MAT2D(i,j,a,nEq) -= MAT2D(i,k,a,nEq)*MAT2D(k,j,a,nEq);
        }
/*.....................................................................*/
  
/*... parte inferior*/
        for(j=i+1;j<nEq;j++){ 
          for(k=0;k<i;k++)
            MAT2D(j,i,a,nEq) -= MAT2D(j,k,a,nEq)*MAT2D(k,i,a,nEq);
          diag = MAT2D(i,i,a,nEq);
          if( fabs(diag) <= conv ){ 
            printf("Erro: termo nulo na diagonal da matriz D\n");
            exit(EXIT_FAILURE);
          }
          MAT2D(j,i,a,nEq) /= diag;
        } 
/*.....................................................................*/
      }
      break;
/*...................................................................*/

/*... versao IJK (dot - version - DooLittte - compact form)*/
/*   varredura                                                     */
/* | i-------f| i - inicio                                         */
/* | | i ---f | f - fim                                            */
/* | | |      |                                                    */
/* | f f      |                                                    */
/*                                                                 */
    case DOOLITTLELC:
      for(i=1;i<nEq;i++){
/*... */
        for(j=1;j<=i;j++){
          ii   = j - 1;
          diag = MAT2D(ii,ii,a,nEq);
          if( fabs(diag) <= conv ){ 
            printf("Erro: termo nulo na diagonal da matriz D\n");
            exit(EXIT_FAILURE);
          }
          MAT2D(i,ii,a,nEq) /= diag;
          for(k=0;k<j;k++)
            MAT2D(i,j,a,nEq) -= MAT2D(i,k,a,nEq)*MAT2D(k,j,a,nEq);
        }
/*.....................................................................*/
  
/*... */
        for(j=i+1;j<nEq;j++){ 
          for(k=0;k<i;k++)
            MAT2D(i,j,a,nEq) -= MAT2D(i,k,a,nEq)*MAT2D(k,j,a,nEq);
        } 
/*.....................................................................*/
      }
      break;
/*...................................................................*/

/*... versao IJK (dot - version - DooLittte - compact form)*/
/*   varredura                                                     */
/* | i i i i  | i - inicio                                         */
/* | i f | |  | f - fim                                            */
/* | i - f |  |                                                    */
/* | i - - f  |                                                    */
/*                                                                 */
    case DOOLITTLECL:
      for(i=1;i<nEq;i++){
/*... inferior*/
        for(j=0;j<i;j++){
/*... inferior*/
          for(k=0;k<j;k++)
            MAT2D(i,j,a,nEq) -= MAT2D(i,k,a,nEq)*MAT2D(k,j,a,nEq);
          diag = MAT2D(j,j,a,nEq);
          if( fabs(diag) <= conv ){ 
            printf("Erro: termo nulo na diagonal da matriz D\n");
            exit(EXIT_FAILURE);
          }
          MAT2D(i,j,a,nEq) /= diag;
/*....................................................................*/

/*... superior*/
          ii = j + 1;
          for(k=0;k<ii;k++)
            MAT2D(ii,i,a,nEq) -= MAT2D(ii,k,a,nEq)*MAT2D(k,i,a,nEq);
/*.....................................................................*/
        }
/*.....................................................................*/
      }
      break;
/*...................................................................*/

/*...*/
    default :
      ERRO_OP_NEW(__FILE__,__func__,__LINE__,"fatLU",code);
/*..................................................................*/
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * fatLUpp :  fatora da matriz A em L e U com pivatiamento parcial   *
 * dinamico                                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa                       * 
 * p(neq)      -> nao definido                                       * 
 * nEq         -> numero de equacoes                                 * 
 * code        ->  4 - LUKIJ                                         * 
 *             ->  5 - LUIKJ                                         * 
 *             ->  6 - LUKJI                                         * 
 *             ->  7 - LUJKI                                         * 
 *             ->  8 - LUIJK                                         * 
 *             ->  9 - DOOLITTLELC                                   *  
 *             -> 10 - DOOLITTLECL                                   *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * a           -> atualizado com L e U                               * 
 * p           -> matriz de permitacao                               * 
 *-------------------------------------------------------------------* 
 * OBS: matriz cheia                                                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void fatLUpp(DOUBLE *RESTRICT a,INT    *RESTRICT p
            ,INT const nEq     ,short const code){
  
  int i,j,k;
  int ii,jj;
  DOUBLE itroca;
  DOUBLE tmp1,tmp2;
  DOUBLE conv,diag;
  DOUBLE zero = 1.0e-15;

  tmp1 = 0.0e0;
  for(i=0;i<nEq;i++)
    tmp1 += MAT2D(i,i,a,nEq)*MAT2D(i,i,a,nEq);
  conv = sqrt(tmp1)*zero;

  switch(code){
/*... versao KIJ - (escalar vector)*/
    case LUKIJPP:
      for(k=0;k<nEq;k++)
        p[k] = k;
      for(k=0;k<nEq-1;k++){
/*...busca na coluna k com o maior pivo em valor absoluto*/
        tmp1   = fabs(MAT2D(k,k,a,nEq));
        itroca = 0;
        jj     = k;
        for(j=k+1  ;j<nEq;j++){
          tmp2 = fabs(MAT2D(j,k,a,nEq));
          if( tmp1 < tmp2){
            tmp1   = tmp2;
            jj     = j;
            itroca = 1;
          } 
        }
        if(itroca){
          ii       = p[k];
          p[k]     = p[jj];
          p[jj]    = ii; 
          for(j=0;j<nEq;j++){
            tmp1                = MAT2D(k   ,j,a,nEq);
            MAT2D(k   ,j,a,nEq) = MAT2D(jj  ,j,a,nEq);
            MAT2D(jj  ,j,a,nEq) = tmp1;
          }
        } 
/*...................................................................*/
          
          
/*... fatoracao*/        
        diag  = MAT2D(k,k,a,nEq);
        if( fabs(diag) <= conv ){ 
          printf("Erro: termo nulo na diagonal da matriz D\n");
          exit(EXIT_FAILURE);
        }
        ii = k + 1;
        for(i=ii;i<nEq;i++){
          tmp1 = MAT2D(i,k,a,nEq) / diag; 
          MAT2D(i,  k,a,nEq) = tmp1;
          addVector(1.0e0  ,&MAT2D(i   ,ii,a,nEq)
                   ,-tmp1  ,&MAT2D(k   ,ii,a,nEq)
                   ,nEq-ii ,&MAT2D(i   ,ii,a,nEq));
//          for(j=ii;j<nEq;j++)
//              MAT2D(i,j,a,nEq) -= tmp1*MAT2D(k,j,a,nEq);  
        }
      }
    break;
/*...................................................................*/

/*... versao IJK (dot - version - DooLittte - compact form)*/
/*   varredura                                                     */
/* | i i i i  | i - inicio                                         */
/* | i f | |  | f - fim                                            */
/* | i - f |  |                                                    */
/* | i - - f  |                                                    */
/*                                                                 */
    case DOOLITTLECLPP:
/*...*/
      for(k=0;k<nEq;k++)
          p[k] = k;
      for(i=0;i<nEq;i++){
/*...busca na coluna k com o maior pivo em valor absoluto*/
        tmp1   = fabs(MAT2D(i,i,a,nEq));
        itroca = 0;
        jj     = i;
        for(j=i+1  ;j<nEq;j++){
          tmp2 = fabs(MAT2D(j,i,a,nEq));
          if( tmp1 < tmp2){
            tmp1   = tmp2;
            jj     = j;
            itroca = 1;
          } 
        }
        if(itroca){
          ii       = p[i];
          p[i]     = p[jj];
          p[jj]    = ii; 
          for(j=0;j<nEq;j++){
            tmp1                = MAT2D(i   ,j,a,nEq);
            MAT2D(i   ,j,a,nEq) = MAT2D(jj  ,j,a,nEq);
            MAT2D(jj  ,j,a,nEq) = tmp1;
          }
        }
/*...................................................................*/

/*... inferior*/
        for(j=0;j<i;j++){
/*... inferior*/
          for(k=0;k<j;k++)
            MAT2D(i,j,a,nEq) -= MAT2D(i,k,a,nEq)*MAT2D(k,j,a,nEq);
          diag = MAT2D(j,j,a,nEq);
          if( fabs(diag) <= conv ){ 
            printf("Erro: termo nulo na diagonal da matriz D\n");
            exit(EXIT_FAILURE);
          }
          MAT2D(i,j,a,nEq) /= diag;
/*....................................................................*/

/*... superior*/
          ii = j + 1;
          for(k=0;k<ii;k++)
            MAT2D(ii,i,a,nEq) -= MAT2D(ii,k,a,nEq)*MAT2D(k,i,a,nEq);
/*.....................................................................*/
        }
      }
      break;
/*...................................................................*/

/*...*/
    default :
      ERRO_OP_NEW(__FILE__,__func__,__LINE__,"fatLU",code);
/*..................................................................*/
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * solvTri : reslove os sistemas Ly=b, Dz=y e Mtx=z (A=LDMt)         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa L e U                 * 
 * d(neq)      -> matriz D                                           * 
 * b(neq)      -> vetor de forcas                                    * 
 * nEq         -> numero de equacoes                                 * 
 * code        -> 1 - LDMt                                           * 
 *             -> 2 - LDLt                                           * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * b           -> atualizado com o relustado x                       * 
 *-------------------------------------------------------------------* 
 * OBS: matriz cheia                                                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void solvTri(DOUBLE *RESTRICT a  ,DOUBLE *RESTRICT d
            ,DOUBLE *RESTRICT b  ,INT const nEq
            ,short const code){

  int i,j;

/*... Ly = b*/
  for(i=1;i<nEq;i++)
    for(j=0;j<i;j++) 
      b[i]-= MAT2D(i,j,a,nEq)*b[j]; 
/*...................................................................*/

/*...Dz = y */
  for(i=0;i<nEq;i++)
    b[i] /= d[i];
/*...................................................................*/

/*... Mtx = z*/
  if(code == LDMT)
    for(i=nEq-2;i>-1;i--)
      for(j=i+1;j<nEq;j++) 
        b[i] -= MAT2D(i,j,a,nEq)*b[j]; 
/*... Ltx = z*/
  else if(code == LDLT)
    for(i=nEq-2;i>-1;i--)
      for(j=i+1;j<nEq;j++) 
        b[i] -= MAT2D(j,i,a,nEq)*b[j]; 
/*...................................................................*/

}
/*********************************************************************/ 

/********************************************************************* 
 * solvCholesky : reslove os sistemas Gy=b e Gtx=y (A=GGt) (Cholesky)*
 * (A - simetrica definida positiva - SDP)                           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa G e Gt                * 
 * b(neq)      -> vetor de forcas                                    * 
 * nEq         -> numero de equacoes                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * b           -> atualizado com o relustado x                       * 
 *-------------------------------------------------------------------* 
 * OBS: matriz cheia                                                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void solvCholesky(DOUBLE *RESTRICT a,DOUBLE *RESTRICT b,INT const nEq){

  int i,j;

/*... Gy = b*/
  b[0] /= MAT2D(0,0,a,nEq);
  for(i=1;i<nEq;i++){
    for(j=0;j<i;j++) 
      b[i]-= MAT2D(i,j,a,nEq)*b[j]; 
    b[i] /= MAT2D(i,i,a,nEq);
  }
/*...................................................................*/

/*... Gtx = y*/
  i = nEq - 1;
  b[i] /= MAT2D(i,i,a,nEq);
  for(i=nEq-2;i>-1;i--){
    for(j=i+1;j<nEq;j++) 
      b[i] -= MAT2D(j,i,a,nEq)*b[j];
    b[i] /= MAT2D(i,i,a,nEq);
  }
/*...................................................................*/
  
}
/*********************************************************************/ 

/********************************************************************* 
 * solvLU : reslove os sistemas Ly=b e Ux=y (A=LU)                   *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa L e U                 * 
 * b(neq)      -> vetor de forcas                                    * 
 * nEq         -> numero de equacoes                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * b           -> atualizado com o relustado x                       * 
 *-------------------------------------------------------------------* 
 * OBS: matriz cheia                                                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void solvLU(DOUBLE *RESTRICT a,DOUBLE *RESTRICT b,INT const nEq){

  INT i,j;

/*... Ly = b*/
  for(i=1;i<nEq;i++){
    for(j=0;j<i;j++) 
      b[i]-= MAT2D(i,j,a,nEq)*b[j]; 
  }
/*...................................................................*/

/*... Ux = y*/
  i     = nEq - 1;
  b[i] /= MAT2D(i,i,a,nEq);
  for(i=nEq-2;i>-1;i--){
    for(j=i+1;j<nEq;j++) 
      b[i] -= MAT2D(i,j,a,nEq)*b[j];
    b[i] /= MAT2D(i,i,a,nEq);
  }
/*...................................................................*/
  
}
/*********************************************************************/ 

/********************************************************************* 
 * solvLUpp : reslove os sistemas com pivotiamento parcial           *
 * Ly=Pb e Ux=y (A=LU)                                               *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa L e U                 * 
 * b(neq)      -> vetor de forcas                                    * 
 * w(neq)      -> nao definido                                       * 
 * p(neq)      -> matriz de permutacao                               * 
 * nEq         -> numero de equacoes                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * b           -> atualizado com o relustado x                       * 
 *-------------------------------------------------------------------* 
 * OBS: matriz cheia                                                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void solvLUpp(DOUBLE *RESTRICT a,DOUBLE *RESTRICT b
             ,DOUBLE *RESTRICT w,INT *RESTRICT p
             ,int const nEq){

  int i,j;
/*... Pb*/
  for(i=0;i<nEq;i++)
    w[i] = b[p[i]]; 
/*...................................................................*/
  
/*... Ly = Pb*/
  for(i=1;i<nEq;i++){
    for(j=0;j<i;j++) 
      w[i] -= MAT2D(i,j,a,nEq)*w[j]; 
  }
/*...................................................................*/

/*... Ux = y*/
  i     = nEq - 1;
  w[i] /= MAT2D(i,i,a,nEq);
  for(i=nEq-2;i>-1;i--){
    for(j=i+1;j<nEq;j++){ 
      w[i] -= MAT2D(i,j,a,nEq)*w[j];
    }
    w[i] /= MAT2D(i,i,a,nEq);
  }
/*...................................................................*/

/*... */
  for(i=0;i<nEq;i++)
    b[i] = w[i]; 
/*...................................................................*/
  
}
/*********************************************************************/ 

/********************************************************************* 
 * SOLVERD: solver direto                                            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa                       * 
 * b(neq)      -> vetor de forcas                                    * 
 * nEq         -> numero de equacoes                                 * 
 * code        ->  1 - LDMt                                          * 
 *             ->  2 - LDLt                                          * 
 *             ->  3 - GGt                                           * 
 *             ->  4 - LUKIJ                                         * 
 *             ->  5 - LUIKJ                                         * 
 *             ->  6 - LUKJI                                         * 
 *             ->  7 - LUJKI                                         * 
 *             ->  8 - LUIJK                                         * 
 *             ->  9 - DOOLITTLELC                                   *  
 *             -> 10 - DOOLITTLECL                                   *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * a           -> modificado                                         * 
 * b           -> atualizado com a solucao x                         * 
 *-------------------------------------------------------------------* 
 * OBS: matriz cheia                                                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void solverD(DOUBLE *RESTRICT a,DOUBLE *RESTRICT b
            ,int const nEq
            ,short const code  ,bool fscaling){

  DOUBLE  *d=NULL,*r=NULL,*w=NULL,*scaling=NULL,*yr=NULL; 
  LDOUBLE *itr=NULL,*br=NULL,*ar=NULL,*xr=NULL; 
  int    *p=NULL; 
  int sizeDouble  = sizeof(DOUBLE);
  int sizeInt     = sizeof(int);
  int sizeLDouble = sizeof(LDOUBLE);

  switch(code){
/*...*/    
    case LDMT:
/*... rescalando a matriz*/
      if(fscaling){
        scaling = (DOUBLE*) calloc(nEq,sizeDouble);
        ERRO_MALLOC(scaling,"scaling"
                   ,__LINE__,__FILE__,__func__);
/*...*/
        diagScaling(a,scaling,nEq);
        scalingSystem(a,b,scaling,nEq);
/*...................................................................*/

/*...*/
        free(scaling);
      }
/*...................................................................*/

/*... alloc*/
      d = (DOUBLE*) calloc(nEq,sizeDouble);
      w = (DOUBLE*) calloc(nEq,sizeDouble);
      r = (DOUBLE*) calloc(nEq,sizeDouble);
      ERRO_MALLOC(d,"d",__LINE__,__FILE__,__func__);
      ERRO_MALLOC(w,"w",__LINE__,__FILE__,__func__);
      ERRO_MALLOC(r,"r",__LINE__,__FILE__,__func__);
/*...................................................................*/

/*... fatoracao LDMt*/
      fatLDMt(a,d,r,w,nEq);
      printFat(a,d,nEq,LDMT);
      solvTri(a,d,b,nEq,LDMT);
/*...................................................................*/

/*... dealloc*/
      free(r);
      free(w);
      free(d);
/*...................................................................*/
    break;
/*...................................................................*/

/*...*/    
    case LDLT:
/*... alloc*/
      d = (DOUBLE*) calloc(nEq,sizeDouble);
      ERRO_MALLOC(d,"d",__LINE__,__FILE__,__func__);
      r = (DOUBLE*) calloc(nEq,sizeDouble);
      ERRO_MALLOC(r,"r",__LINE__,__FILE__,__func__);
/*...................................................................*/

/*... fatoracao LDLt*/
      fatLDLt(a,d,r,nEq);
      printFat(a,d,nEq,LDLT);
      solvTri(a,d,b,nEq,LDLT);
/*...................................................................*/

/*... dealloc*/
      free(r);
      free(d);
/*...................................................................*/
    break;
/*...................................................................*/

/*...*/    
    case GGT:
/*... fatoracao GGt*/
      ic0(a,nEq);
//    fatGGt(a,nEq);
      printFat(a,d,nEq,GGT);
      solvCholesky(a,b,nEq);
    break;
/*...................................................................*/

/*...*/    
    case LUKIJ:
/*... rescalando a matriz*/
      if(fscaling){
        scaling = (DOUBLE*) calloc(nEq,sizeDouble);
        ERRO_MALLOC(scaling,"scaling"
                   ,__LINE__,__FILE__,__func__);
/*...*/
        diagScaling(a,scaling,nEq);
        scalingSystem(a,b,scaling,nEq);
/*...................................................................*/

/*...*/
        free(scaling);
      }
/*...................................................................*/

/*... fatoracao LU*/
      fatLU(a,nEq,LUKIJ);
      printFat(a,d,nEq,LUKIJ);
      solvLU(a,b,nEq);
/*...................................................................*/
    break;
/*...................................................................*/

/*...*/    
    case LUIKJ:
/*... rescalando a matriz*/
      if(fscaling){
        scaling = (DOUBLE*) calloc(nEq,sizeDouble);
        ERRO_MALLOC(scaling,"scaling"
                   ,__LINE__,__FILE__,__func__);
/*...*/
        diagScaling(a,scaling,nEq);
        scalingSystem(a,b,scaling,nEq);
/*...................................................................*/

/*...*/
        free(scaling);
      }
/*...................................................................*/

/*... fatoracao LU*/
      fatLU(a,nEq,LUIKJ);
      printFat(a,d,nEq,LUIKJ);
      solvLU(a,b,nEq);
/*...................................................................*/
    break;
/*...................................................................*/

/*...*/    
    case LUKJI:
/*... rescalando a matriz*/
      if(fscaling){
        scaling = (DOUBLE*) calloc(nEq,sizeDouble);
        ERRO_MALLOC(scaling,"scaling"
                   ,__LINE__,__FILE__,__func__);
/*...*/
        diagScaling(a,scaling,nEq);
        scalingSystem(a,b,scaling,nEq);
/*...................................................................*/

/*...*/
        free(scaling);
      }
/*...................................................................*/

/*... fatoracao LU*/
      fatLU(a,nEq,LUKJI);
      printFat(a,d,nEq,LUKJI);
      solvLU(a,b,nEq);
/*...................................................................*/
    break;
/*...................................................................*/

/*...*/    
    case LUJKI:
/*... rescalando a matriz*/
      if(fscaling){
        scaling = (DOUBLE*) calloc(nEq,sizeDouble);
        ERRO_MALLOC(scaling,"scaling"
                   ,__LINE__,__FILE__,__func__);
/*...*/
        diagScaling(a,scaling,nEq);
        scalingSystem(a,b,scaling,nEq);
/*...................................................................*/

/*...*/
        free(scaling);
      }
/*...................................................................*/

/*... fatoracao LU*/
      fatLU(a,nEq,LUJKI);
      printFat(a,d,nEq,LUJKI);
      solvLU(a,b,nEq);
/*...................................................................*/

/*...*/    
    case DOOLITTLEL:
/*... rescalando a matriz*/
      if(fscaling){
        scaling = (DOUBLE*) calloc(nEq,sizeDouble);
        ERRO_MALLOC(scaling,"scaling"
                   ,__LINE__,__FILE__,__func__);
/*...*/
        diagScaling(a,scaling,nEq);
        scalingSystem(a,b,scaling,nEq);
/*...................................................................*/

/*...*/
        free(scaling);
      }
/*...................................................................*/

/*... fatoracao LU*/
      fatLU(a,nEq,DOOLITTLEL);
      printFat(a,d,nEq,DOOLITTLEL);
      solvLU(a,b,nEq);
/*...................................................................*/
    break;
/*...................................................................*/

/*...*/    
    case DOOLITTLELC:
/*... rescalando a matriz*/
      if(fscaling){
        scaling = (DOUBLE*) calloc(nEq,sizeDouble);
        ERRO_MALLOC(scaling,"scaling"
                   ,__LINE__,__FILE__,__func__);
/*...*/
        diagScaling(a,scaling,nEq);
        scalingSystem(a,b,scaling,nEq);
/*...................................................................*/

/*...*/
        free(scaling);
      }
/*...................................................................*/

/*... fatoracao LU*/
      fatLU(a,nEq,DOOLITTLELC);
      printFat(a,d,nEq,DOOLITTLELC);
      solvLU(a,b,nEq);
/*...................................................................*/
    break;
/*...................................................................*/

/*...*/    
    case DOOLITTLECL:
/*...................................................................*/

/*... rescalando a matriz*/
      if(fscaling){
        scaling = (DOUBLE*) calloc(nEq,sizeDouble);
        ERRO_MALLOC(scaling,"scaling"
                   ,__LINE__,__FILE__,__func__);
/*...*/
        diagScaling(a,scaling,nEq);
        scalingSystem(a,b,scaling,nEq);
/*...................................................................*/

/*...*/
        free(scaling);
      }
/*...................................................................*/

/*... fatoracao LU*/
      fatLU(a,nEq,DOOLITTLECL);
      printFat(a,d,nEq,DOOLITTLECL);
      solvLU(a,b,nEq);
/*...................................................................*/
    break;
/*...................................................................*/

/*... fatoracao com pivotiamento parcial*/    
    case LUKIJPP:
/*...*/
      w       = (DOUBLE*) calloc(nEq,sizeDouble);
      ERRO_MALLOC(w,"w",__LINE__,__FILE__,__func__);
      p       = (int*) calloc(nEq,sizeInt);
      ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
/*...................................................................*/

/*... rescalando a matriz*/
      if(fscaling){
        scaling = (DOUBLE*) calloc(nEq,sizeDouble);
        ERRO_MALLOC(scaling,"scaling"
                   ,__LINE__,__FILE__,__func__);
/*...*/
        diagScaling(a,scaling,nEq);
        scalingSystem(a,b,scaling,nEq);
/*...................................................................*/

/*...*/
        free(scaling);
      }
/*...................................................................*/

/*... fatoracao LU*/
      fatLUpp(a,p,nEq,LUKIJPP);
      printFat(a,d,nEq,LUKIJPP);
      solvLUpp(a,b,w,p,nEq);
/*...................................................................*/

/*... dealloc*/
      free(p);
      free(w);
/*...................................................................*/
    break;
/*...................................................................*/

/*... fatoracao com pivotiamento parcial*/    
    case DOOLITTLECLPP:
/*... rescalando a matriz*/
      if(fscaling){
        scaling = (DOUBLE*) calloc(nEq,sizeDouble);
        ERRO_MALLOC(scaling,"scaling"
                   ,__LINE__,__FILE__,__func__);
/*...*/
        diagScaling(a,scaling,nEq);
        scalingSystem(a,b,scaling,nEq);
/*...................................................................*/

/*...*/
        free(scaling);
      }
      printf("A %lf %lf %lf %lf %lf\n",a[0] ,a[1],a[2],a[3],a[4]);
      printf("A %lf %lf %lf %lf %lf\n",a[5] ,a[6],a[7],a[8],a[9]);
      printf("A %lf %lf %lf %lf %lf\n",a[10],a[11],a[12],a[13],a[14]);
      printf("A %lf %lf %lf %lf %lf\n",a[15],a[16],a[17],a[18],a[19]);
      printf("A %lf %lf %lf %lf %lf\n",a[20],a[21],a[22],a[23],a[24]);
/*...................................................................*/

/*... refinamento iterativo*/
      itr     = (LDOUBLE*) calloc(nEq,sizeLDouble);
      ERRO_MALLOC(itr,"itr",__LINE__,__FILE__,__func__);
      br      = (LDOUBLE*) calloc(nEq,sizeLDouble);
      ERRO_MALLOC(itr,"itr",__LINE__,__FILE__,__func__);
      xr      = (LDOUBLE*) calloc(nEq,sizeLDouble);
      ERRO_MALLOC(xr ,"xr" ,__LINE__,__FILE__,__func__);
      ar      = (LDOUBLE*) calloc(nEq*nEq,sizeLDouble);
      ERRO_MALLOC(itr,"itr",__LINE__,__FILE__,__func__);
      yr      =  (DOUBLE*) calloc(nEq    ,sizeLDouble);
      ERRO_MALLOC(yr,"yr",__LINE__,__FILE__,__func__);
      
/*...*/
      initIterImprovmentData(ar,br,a,b,nEq);
/*...................................................................*/

/*...*/
      w       = (DOUBLE*) calloc(nEq,sizeDouble);
      ERRO_MALLOC(w,"w",__LINE__,__FILE__,__func__);
      p       = (int*) calloc(nEq,sizeInt);
      ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
/*...................................................................*/

      printf("%lf %lf %lf %lf\n",a[0],a[1],a[2],a[3]);
      printf("%lf %lf\n",b[0],b[1]);
/*... fatoracao LU*/
      fatLUpp(a,p,nEq ,DOOLITTLECLPP);
      printf("%d %d\n",p[0],p[1]);
      printFat(a,d,nEq,DOOLITTLECLPP);
      solvLUpp(a,b,w,p,nEq);
/*..................................................................*/

/*... refinamento iterativo*/
      iterImprovement(ar,br,itr,xr,yr,a,b,w,p,nEq);
      free(yr );
      free(ar );
      free(xr );
      free(br );
      free(itr);
/*..................................................................*/

/*... dealloc*/
      free(p);
      free(w);
/*..................................................................*/
    break; 
/*..................................................................*/

/*...*/
    default :
      ERRO_OP_NEW(__FILE__,__func__,__LINE__,"solverD",code);
/*..................................................................*/
  }
}
/*********************************************************************/ 
  
/********************************************************************* 
 * DIAGSACLING: calculo da matrix diagonal de escala                 * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa                       * 
 * scaling     -> nao definido                                       * 
 * nEq         -> numero de equacoes                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * scaling     -> matrix diagonal de escala                          * 
 *-------------------------------------------------------------------* 
 * OBS: norma infinita das linhas                                    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void diagScaling(DOUBLE *RESTRICT a,DOUBLE *RESTRICT scaling
                ,INT const nEq)
{
  int i,j;
  DOUBLE tmp;
  for(i=0;i<nEq;i++){
    tmp = 0.l;
    for(j=0;j<nEq;j++){
      tmp = max(tmp,fabs(MAT2D(i,j,a,nEq)));  
    }
    scaling[i] = tmp;
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * SACLINGSYSTEM: calculo da sistema equacoes equivalente            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa                       * 
 * b(neq)      -> vetor de forcas                                    * 
 * scaling(neq)-> matrix diagonal de escala                          * 
 * nEq         -> numero de equacoes                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes atualizada com D-1A         * 
 * b(neq)      -> vetor de forcas atualizado com D-1A                * 
 *-------------------------------------------------------------------* 
 * OBS: norma infinita das linhas                                    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void scalingSystem(DOUBLE *RESTRICT a      ,DOUBLE *RESTRICT b
                  ,DOUBLE *RESTRICT scaling
                  ,int const nEq)
{
  int i,j;
  for(i=0;i<nEq;i++){
    for(j=0;j<nEq;j++){
      MAT2D(i,j,a,nEq)/=scaling[i];  
    }
    b[i]/=scaling[i];
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * INITITERIMPROVEMNTDATA: inicializacao da estrutura de dados do    * 
 * processo de refinamento iterativo                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ar(neq,neq) -> nao definido                                       * 
 * br(neq)     -> nao definido                                       * 
 *  a(neq,neq) -> matriz de coeficientes                             * 
 *  b(neq)     -> vetor de forcas                                    * 
 * nEq         -> numero de equacoes                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ar(neq,neq) -> matriz de coeficientes armazenado em precisao      *
 * extendida                                                         *
 * br(neq)     -> vetor de forcas aramzendo em precisao extendida    * 
 *-------------------------------------------------------------------* 
 * OBS: norma infinita das linhas                                    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void initIterImprovmentData(LDOUBLE *RESTRICT ar,LDOUBLE *RESTRICT br
                       ,DOUBLE  *RESTRICT a ,DOUBLE *RESTRICT b
                       ,INT const nEq)
{
  INT i,j;
  for(i=0;i<nEq;i++){
    br[i] = b[i];
    for(j=0;j<nEq;j++)
      MAT2D(i,j,ar,nEq) = MAT2D(i,j,a,nEq);
  }

}
/*********************************************************************/ 

/********************************************************************* 
 * INITIMPROVEMNTDATA: refininamento iterativo da solucao            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ar(neq,neq) -> matriz de coeficientes armazenado em precisao      *
 * extendida                                                         *
 * br(neq)     -> vetor de forcas aramzendo em precisao extendida    * 
 * itr(neq)    -> vetor auxiliar em precisao extendida               * 
 * xr(neq)     -> vetor auxiliar em precisao extendia                * 
 * x (neq)     -> vetor de forcas aramzendo em precisao ex tendida   * 
 *  a(neq,neq) -> matriz de coeficientes fatorada                    * 
 *  b(neq)     -> resultado do sistema de equacoes                   * 
 *  w(neq)     -> vetor auxiliar                                     * 
 *  p(neq)     -> matriz de permutacao                               * 
 * nEq         -> numero de equacoes                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *  b(neq)     -> resultado refinado                                 * 
 *-------------------------------------------------------------------* 
 * OBS: norma infinita das linhas                                    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void iterImprovement(LDOUBLE *RESTRICT ar ,LDOUBLE *RESTRICT br
                    ,LDOUBLE *RESTRICT itr,LDOUBLE *RESTRICT xr
                    ,DOUBLE  *RESTRICT x
                    ,DOUBLE  *RESTRICT a  ,DOUBLE *RESTRICT b
                    ,DOUBLE  *RESTRICT w                     
                    ,INT *RESTRICT p      ,INT const nEq)
{
  INT i;
  DOUBLE  tmp;
  short it;
  LDOUBLE normR,normR0=1.0e0;

  for(it=0;it<5;it++){
    for(i=0;i<nEq;i++){
      tmp   = b[i];
      xr[i] = tmp;
      x[i]  = tmp;
     }
/*... Ax*/
    lmatVecFull(ar,xr,itr,nEq,nEq);
/*..................................................................*/

/*...r = b -Ax*/
    for(i=0;i<nEq;i++)
      itr[i] = br[i] - itr[i];
    if(it == 0)  normR0 = lnormInf(itr,nEq,1);
    normR = lnormInf(itr,nEq,1);
    printf("it: %d Residuo do refinamento: %Le\n",it,normR/normR0);
/*..................................................................*/

/*...*/
    for(i=0;i<nEq;i++)
      b[i] = itr[i];
/*..................................................................*/

/*...Az=r*/
    solvLUpp(a,b,w,p,nEq);
/*..................................................................*/
  
/*...x = x + z*/
    for(i=0;i<nEq;i++)
      b[i] += x[i];
/*..................................................................*/
  }
/*..................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * PRINTFAT : Imprime a matriz A fatorada                            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a(neq,neq)  -> matriz de coeficientes densa L e U                 * 
 * d(neq)      -> matriz D                                           * 
 * n           -> numero de linhas                                   * 
 * code        ->  1 - LDMt                                          * 
 *             ->  2 - LDLt                                          * 
 *             ->  3 - GGt                                           * 
 *              >  4 - LUKIJ                                         * 
 *             ->  5 - LUIKJ                                         * 
 *             ->  6 - LUKJI                                         * 
 *             ->  7 - LUJKI                                         * 
 *             ->  8 - LUIJK                                         * 
 *             ->  9 - DOOLITTLELC                                   *  
 *             -> 10 - DOOLITTLECL                                   *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void printFat(DOUBLE *a,DOUBLE *d,int const n,short code){
  
  int i,j;

  switch(code){
/*... LDMt*/
    case LDMT:
      printf("matrix L:\n");
      for(i=0;i<n;i++){
        for(j=0;j<n;j++)
          if( i >= j) printf("%lf ",MAT2D(i,j,a,n));
          else        printf("%lf ",0.0e0         );
        printf("\n");
      }  
  
      printf("matrix D:\n");
      for(i=0;i<n ;i++){
        printf("%lf ",d[i]);
        printf("\n");
      }
  
      printf("matrix MT:\n");
      for(i=0;i<n;i++){
        for(j=0;j<n;j++)
          if( i <= j) printf("%lf ",MAT2D(i,j,a,n));
          else       printf("%lf ",0.0e0         );
        printf("\n");
      }
    break;  
/*...................................................................*/

/*... LDLt*/
    case LDLT:
      printf("matrix L:\n");
      for(i=0;i<n;i++){
        for(j=0;j<n;j++)
          if( i >= j) printf("%lf ",MAT2D(i,j,a,n));
          else        printf("%lf ",0.0e0         );
        printf("\n");
      }  
  
      printf("matrix D:\n");
      for(i=0;i<n ;i++){
        printf("%lf ",d[i]);
        printf("\n");
      }
  
      printf("matrix LT:\n");
      for(i=0;i<n;i++){
        for(j=0;j<n;j++)
          if( i <= j) printf("%lf ",MAT2D(j,i,a,n));
          else       printf("%lf ",0.0e0         );
        printf("\n");
      }
    break;
/*...................................................................*/

/*... GGt*/
    case GGT:
      printf("matrix G:\n");
      for(i=0;i<n;i++){
        for(j=0;j<n;j++)
          if( i >= j) printf("%lf ",MAT2D(i,j,a,n));
          else       printf("%lf ",0.0e0         );
        printf("\n");
      }  
  
      printf("matrix GT:\n");
      for(i=0;i<n;i++){
        for(j=0;j<n;j++)
          if( i <= j) printf("%lf ",MAT2D(j,i,a,n));
          else       printf("%lf ",0.0e0         );
        printf("\n");
      }
    break;
/*...................................................................*/

/*... LU*/
    case LUKIJ:
    case LUIKJ:
    case LUKJI:
    case LUJKI:
    case DOOLITTLEL:
    case DOOLITTLELC: 
    case DOOLITTLECL: 
    case LUKIJPP: 
    case DOOLITTLECLPP: 
      printf("matrix L:\n");
      for(i=0;i<n;i++){
        for(j=0;j<n;j++)
          if( i > j) printf("%lf ",MAT2D(i,j,a,n));
          else if( i == j) printf("%lf ",1.0e0);
          else       printf("%lf ",0.0e0         );
        printf("\n");
      }  
  
      printf("matrix U:\n");
      for(i=0;i<n;i++){
        for(j=0;j<n;j++)
          if( i <= j) printf("%lf ",MAT2D(i,j,a,n));
          else       printf("%lf ",0.0e0         );
        printf("\n");
      }
    break;
/*...................................................................*/
  }
}
/*********************************************************************/ 
