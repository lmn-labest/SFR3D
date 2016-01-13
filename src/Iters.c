#include<Solv.h>
/**********************************************************************
 * PCG  : metodo do gradiente conjugado com precondiconador diagonal  *
 * (M-1Ax=M-1b) (matriz simentrica)                                   *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 *  neq -> numero de equacoes                                         *
 *  nad -> numero de elementos nao nulos fora da diagonal             *
 *  ia  -> estrutura de dados para matriz esparsa A                   *
 *  ja  -> estrutura de dados para matriz esparsa A                   *
 *  al  -> parte inferior da matriz A                                 *
 *  ad  -> diagnal da matriz A                                        *
 *  au  -> parte superior da matriz A                                 *
 *   p  -> precondiconador diagonal                                   *
 *   b  -> vetor b (Ax=b)                                             *
 *   x  -> vetor de solucao                                           *
 *   z  -> vetor auxiliar                                             *
 *   r  -> vetor auxiliar                                             *
 * newX -> vetor inicial iniciado com zero                            *
 * fLog -> arquivo de log do solver                                   *
 * log  -> log de arquivo (true|false)                                *
 * tol  -> tolerancia do solver                                       *
 * maxIt-> numero maximo de iteracoes                                 *
 *  newX-> true zero o vetor inicial                                  *
 *  log -> escreve o log do solver                                    *
 *fPrint-> saida de informacao na tela                                *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *  
 * b[]-> alterado                                                     *
 * ad,al,au-> inalterado                                              *
 * -------------------------------------------------------------------*
*********************************************************************/
void pcg(INT const nEq      ,INT const nAd  
        ,INT *restrict ia   ,INT *restrict ja
        ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
        ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict x
        ,DOUBLE *restrict z ,DOUBLE *restrict r ,DOUBLE const tol
        ,unsigned int maxIt ,bool const newX          
        ,FILE* fLog         ,bool const log
        ,bool const fPrint
        ,void(*matvec)()    ,DOUBLE(*dot)())
{
  INT i,j;
  DOUBLE alpha,beta,d,conv,energy;
  DOUBLE timei,timef;
  timei = getTimeC();

/* chute inicial*/
  if(newX)  
    for(i = 0; i < nEq; i++)  
      x[i] = 0.e0;
  
  matvec(nEq,ia,ja,al,ad,x,z);
  
  for(i = 0; i < nEq; i++)   {
    r[i] = b[i] - z[i];
    z[i] = r[i] * m[i];
    b[i] = z[i];
  }
  
  d    = dot(r,z,nEq);
  conv = tol * sqrt(fabs(d));
/*--------------------------------------------------------*/   
  for(j = 0; j < maxIt; j++)   {
    matvec(nEq,ia,ja,al,ad,b,z);
    alpha = d / dot(b,z,nEq);
    for(i = 0; i < nEq; i++)   {
      x[i] +=  alpha * b[i];
      r[i] -=  alpha * z[i];
      z[i]  =   r[i] * m[i];
    }
    beta = dot(r,z,nEq)/d;
    for(i = 0; i < nEq; i++)   {
      b[i] = z[i] + beta * b[i];
    }
    d = beta * d;
    if (sqrt(fabs(d)) < conv) break;
  }
/* -------------------------------------------------------*/
  matvec(nEq,ia,ja,al,ad,x,z);
/*norma de energia = xT*A*x */
  energy = dot(x,z,nEq);
/* -------------------------------------------------------*/
  timef = getTimeC() - timei;

  if(fPrint){ 
    printf("\tnad         :      %20d\n"  ,nAd);
    printf("\tSolver tol  :      %20.2e\n",tol);
    printf(" (PCG) solver:\n"
           "\tEquations   =      %20d\n"
           "\tIterarions  =      %20d\n"
	         "\tEnergy norm =      %20.12e\n"
	         "\tCPU time(s) =      %20.5lf\n" 
	         ,nEq,j+1,energy,timef);
  }
  
  if(j == maxIt){ 
    printf(" *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
    exit(EXIT_FAILURE);
  }
  if(log)
    fprintf(fLog          
           ,"PCG: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);
}
/**********************************************************************/

/**********************************************************************
 * MPIPCG : metodo do gradiente conjugado com precondiconador diagonal*
 * (M-1Ax=M-1b) (matriz simentrica)                                   *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 *  neq -> numero de equacoes                                         *
 *neqNov-> numero de equacoes nao sobrepostas                         *
 *  nad -> numero de elementos nao nulos fora da diagonal             *
 *  nadR-> numero de elementos nao nulos na parte retangular          *
 *  ia  -> estrutura de dados para matriz esparsa A                   *
 *  ja  -> estrutura de dados para matriz esparsa A                   *
 *  al  -> parte inferior da matriz A                                 *
 *  ad  -> diagnal da matriz A                                        *
 *  au  -> parte superior da matriz A                                 *
 *   p  -> precondiconador diagonal                                   *
 *   b  -> vetor b (Ax=b)                                             *
 *   x  -> vetor de solucao                                           *
 *   z  -> vetor auxiliar                                             *
 *   r  -> vetor auxiliar                                             *
 * newX -> vetor inicial iniciado com zero                            *
 * fLog -> arquivo de log do solver                                   *
 * log  -> log de arquivo (true|false)                                *
 * tol  -> tolerancia do solver                                       *
 * maxIt-> numero maximo de iteracoes                                 *
 *  newX-> true zero o vetor inicial                                  *
 *  log -> escreve o log do solver                                    *
 *fPrint-> saida de informacao na tela                                *
 * iNeq -> mapa de interface de equacoes                              *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *  
 * b[]-> alterado                                                     *
 * ad,al,au-> inalterado                                              *
 * -------------------------------------------------------------------*
*********************************************************************/
void mpiPcg(INT const nEq   ,INT const nEqNov
        ,INT const nAd      ,INT const nAdR  
        ,INT *restrict ia   ,INT *restrict ja
        ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
        ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict x
        ,DOUBLE *restrict z ,DOUBLE *restrict r ,DOUBLE const tol
        ,unsigned int maxIt ,bool const newX          
        ,FILE* fLog         ,bool const log
        ,bool const fPrint
        ,Interface *iNeq                      
        ,void(*matvec)()    ,DOUBLE(*dot)())
{
  INT i,j;
  DOUBLE alpha,beta,d,conv,energy;
  DOUBLE timei,timef;
  INT param[2];
  short st=0;

  timei = getTimeC();

  param[0] = nAd;
  param[1] = nAdR;

/* chute inicial*/
  if(newX)  
    for(i = 0; i < nEq; i++)  
      x[i] = 0.e0;
  
  matvec(nEqNov,param,ia,ja,al,ad,x,z,iNeq);
  
  for(i = 0; i < nEqNov; i++)   {
    r[i] = b[i] - z[i];
    z[i] = r[i] * m[i];
    b[i] = z[i];
  }
  
  d    = dot(r,z,nEqNov);
  conv = tol * sqrt(fabs(d));
/*--------------------------------------------------------*/   
  for(j = 0; j < maxIt; j++)   {
    matvec(nEqNov,param,ia,ja,al,ad,b,z,iNeq);
    
    alpha = d / dot(b,z,nEqNov);
    for(i = 0; i < nEqNov; i++)   {
      x[i] +=  alpha * b[i];
      r[i] -=  alpha * z[i];
      z[i]  =   r[i] * m[i];
    }
    beta = dot(r,z,nEqNov)/d;
    for(i = 0; i < nEqNov; i++)   {
      b[i] = z[i] + beta * b[i];
    }
    d = beta * d;
    if (sqrt(fabs(d)) < conv) break;
  }
/* -------------------------------------------------------*/
  matvec(nEqNov,param,ia,ja,al,ad,x,z,iNeq);
/*norma de energia = xT*A*x */
  energy = dot(x,z,nEqNov);
/* -------------------------------------------------------*/
  timef = getTimeC() - timei;

  if(fPrint){ 
    printf("\tnad         :      %20d\n"  ,nAd);
    printf("\tnadR        :      %20d\n"  ,nAdR);
    printf("\tSolver tol  :      %20.2e\n",tol);
    printf(" (MPIPCG) solver:\n"
           "\tEquations   =      %20d\n"
           "\tEquations   =      %20d\n"
           "\tIterarions  =      %20d\n"
	         "\tEnergy norm =      %20.12e\n"
	         "\tCPU time(s) =      %20.5lf\n" 
	         ,nEqNov,nEq,j+1,energy,timef);
  }
  
  if(j == maxIt) st = -1; 
  
  MPI_Bcast(&st,1,MPI_SHORT,0,mpiVar.comm);
  if( st == -1){
    mpiStop();
    printf(" *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
    exit(EXIT_SUCCESS); 
  }

  if(!mpiVar.myId && log)
    fprintf(fLog          
           ,"MPIPCG: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);
}
/**********************************************************************/

/**********************************************************************
 * PBICGSTAB : metodo do gradiente conjugado bi-ortaganalizado        *
 * com precondiconador diagonal (M-1Ax=M-1b) (matriz geral)           *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 *  neq -> numero de equacoes                                         *
 *  nad -> numero de elementos nao nulos fora da diagonal             *
 *  ia  -> estrutura de dados para matriz esparsa A                   *
 *  ja  -> estrutura de dados para matriz esparsa A                   *
 *  al  -> parte inferior da matriz A                                 *
 *  ad  -> diagnal da matriz A                                        *
 *  au  -> parte superior da matriz A                                 *
 *   p  -> precondiconador diagonal                                   *
 *   b  -> vetor b (Ax=b)                                             *
 *   x  -> vetor de solucao                                           *
 *   t  -> vetor auxiliar                                             *
 *   v  -> vetor auxiliar                                             *
 *   r  -> vetor auxiliar                                             *
 *   p  -> vetor auxiliar                                             *
 *   z  -> vetor auxiliar                                             *
 * newX -> vetor inicial iniciado com zero                            *
 * fLog -> arquivo de log do solver                                   *
 * log  -> log de arquivo (true|false)                                *
 * tol  -> tolerancia do solver                                       *
 * maxIt-> numero maximo de iteracoes                                 *
 *  newX-> true zero o vetor inicial                                  *
 *  log -> escreve o log do solver                                    *
 *fPrint-> saida de informacao na tela                                *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *  
 * b[]-> alterado                                                     *
 * ad,al,au-> inalterado                                              *
 * -------------------------------------------------------------------*
*********************************************************************/
void pbicgstab(INT const nEq  ,INT const nAd
          ,INT *restrict ia   ,INT *restrict ja
          ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
          ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict x
          ,DOUBLE *restrict t ,DOUBLE *restrict v ,DOUBLE *restrict r
          ,DOUBLE *restrict p ,DOUBLE *restrict z ,DOUBLE const tol
          ,unsigned int maxIt ,bool const newX          
          ,FILE* fLog         ,bool const log
          ,bool const fPrint
          ,void(*matvec)()    ,DOUBLE(*dot)())
{
  INT i,j;
  DOUBLE alpha,beta,d,conv,energy,w,rr0;
  DOUBLE timei,timef;

  timei = getTimeC();

/* chute inicial*/
  if(newX)  
    for (i = 0; i < nEq; i++)  
      x[i] = 0.e0;
      
/* Residuo inicial*/  
  matvec(nEq,ia,ja,al,ad,x,z);
  
  for(i = 0; i < nEq; i++)   {
    r[i] = b[i] - z[i];
    p[i] = r[i];
    b[i] = p[i];
    z[i] = p[i]*m[i];
  }
  d = dot(r,z,nEq);
  conv = tol * sqrt(fabs(d));
/*--------------------------------------------------------*/   
  for(j = 0; j < maxIt; j++)   {
    matvec(nEq,ia,ja,al,ad,z,v);
    rr0   = dot(b,r,nEq);
    alpha = rr0/dot(v,r,nEq);
    for(i = 0; i < nEq; i++)   {
      x[i] +=  alpha * z[i];
      b[i] -=  alpha * v[i];
      z[i]  = b[i] * m[i];
    }
    matvec(nEq,ia,ja,al,ad,z,t);
    w = dot(t,b,nEq)/ dot(t,t,nEq);
    for(i = 0; i < nEq; i++)   {
      x[i] += w*z[i];
      b[i] -= w*t[i];
    }

    d = dot(b,z,nEq);
    if(sqrt(fabs(d)) < conv) break;
    beta = (dot(r,b,nEq)/rr0)*(alpha/w);
    for(i = 0; i < nEq; i++)   {
      p[i]  = b[i] + beta*(p[i]-w*v[i]);
      z[i]  = p[i]*m[i];
    }

  }
/* -------------------------------------------------------*/
  matvec(nEq,ia,ja,al,ad,x,z);
/*norma de energia = xT*A*x */
  energy = dot(x,z,nEq);
/* -------------------------------------------------------*/
  timef = getTimeC() - timei;

  if(fPrint){ 
    printf("\tnad         :      %20d\n"  ,nAd);
    printf("\tSolver tol  :      %20.2e\n",tol);
    printf(" (PBICGSTAB) solver :\n"
           "\tEquations   =      %20d\n"
           "\tIterarions  =      %20d\n"
	         "\tEnergy norm =      %20.12e\n"
	         "\tCPU time(s) =      %20.5lf\n" 
	         ,nEq,j+1,energy,timef);
  }
  
  if(j == maxIt){ 
    printf(" *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
    exit(EXIT_FAILURE);
  }
  if(log)
    fprintf(fLog          
           ,"PBICGSTAB: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);


}
/**********************************************************************/

/**********************************************************************
 * MPIPBICGSTAB : metodo do gradiente conjugado bi-ortaganalizado     *
 * com precondiconador diagonal (M-1Ax=M-1b) (matriz geral)           *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 *  neq -> numero de equacoes                                         *
 *neqNov-> numero de equacoes nao sobrepostas                         *
 *  nad -> numero de elementos nao nulos fora da diagonal             *
 *  nadR-> numero de elementos nao nulos na parte retangular          *
 *  ia  -> estrutura de dados para matriz esparsa A                   *
 *  ja  -> estrutura de dados para matriz esparsa A                   *
 *  al  -> parte inferior da matriz A                                 *
 *  ad  -> diagnal da matriz A                                        *
 *  au  -> parte superior da matriz A                                 *
 *   p  -> precondiconador diagonal                                   *
 *   b  -> vetor b (Ax=b)                                             *
 *   x  -> vetor de solucao                                           *
 *   t  -> vetor auxiliar                                             *
 *   v  -> vetor auxiliar                                             *
 *   r  -> vetor auxiliar                                             *
 *   p  -> vetor auxiliar                                             *
 *   z  -> vetor auxiliar                                             *
 * newX -> vetor inicial iniciado com zero                            *
 * fLog -> arquivo de log do solver                                   *
 * log  -> log de arquivo (true|false)                                *
 * tol  -> tolerancia do solver                                       *
 * maxIt-> numero maximo de iteracoes                                 *
 *  newX-> true zero o vetor inicial                                  *
 *  log -> escreve o log do solver                                    *
 *fPrint-> saida de informacao na tela                                *
 * iNeq -> mapa de interface de equacoes                              *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *  
 * b[]-> alterado                                                     *
 * ad,al,au-> inalterado                                              *
 * -------------------------------------------------------------------*
*********************************************************************/
void mpiPbicgstab(INT const nEq,INT const nEqNov      
          ,INT const nAd       ,INT const nAdR
          ,INT *restrict ia    ,INT *restrict ja
          ,DOUBLE *restrict al ,DOUBLE *restrict ad,DOUBLE *restrict au
          ,DOUBLE *restrict m  ,DOUBLE *restrict b ,DOUBLE *restrict x
          ,DOUBLE *restrict t  ,DOUBLE *restrict v ,DOUBLE *restrict r
          ,DOUBLE *restrict p  ,DOUBLE *restrict z ,DOUBLE const tol
          ,unsigned int maxIt  ,bool const newX          
          ,FILE* fLog          ,bool const log
          ,bool const fPrint
          ,Interface *iNeq    
          ,void(*matvec)()     ,DOUBLE(*dot)())
{
  INT i,j;
  DOUBLE alpha,beta,d,conv,energy,w,rr0;
  DOUBLE timei,timef;
  INT param[2];

  timei = getTimeC();
  
  param[0] = nAd;
  param[1] = nAdR;

/* chute inicial*/
  if(newX)  
    for (i = 0; i < nEq; i++)  
      x[i] = 0.e0;
      
/* Residuo inicial*/  
  matvec(nEqNov,param,ia,ja,al,ad,x,z,iNeq);
  for(i = 0; i < nEqNov; i++)   {
    r[i] = b[i] - z[i];
    p[i] = r[i];
    b[i] = p[i];
    z[i] = p[i]*m[i];
  }
  d = dot(r,z,nEqNov);
  conv = tol * sqrt(fabs(d));
/*--------------------------------------------------------*/   
  for(j = 0; j < maxIt; j++)   {
    matvec(nEqNov,param,ia,ja,al,ad,z,v,iNeq);

    rr0   = dot(b,r,nEqNov);
    alpha = rr0/dot(v,r,nEqNov);
    for(i = 0; i < nEqNov; i++)   {
      x[i] +=  alpha * z[i];
      b[i] -=  alpha * v[i];
      z[i]  = b[i] * m[i];
    }
    matvec(nEqNov,param,ia,ja,al,ad,z,t,iNeq);
    w = dot(t,b,nEqNov)/ dot(t,t,nEqNov);
    for(i = 0; i < nEqNov; i++)   {
      x[i] += w*z[i];
      b[i] -= w*t[i];
    }

    d = dot(b,z,nEqNov);
    if(sqrt(fabs(d)) < conv) break;
    beta = (dot(r,b,nEqNov)/rr0)*(alpha/w);
    for(i = 0; i < nEqNov; i++)   {
      p[i]  = b[i] + beta*(p[i]-w*v[i]);
      z[i]  = p[i]*m[i];
    }

  }
/* -------------------------------------------------------*/
  matvec(nEqNov,param,ia,ja,al,ad,x,z,iNeq);
/*norma de energia = xT*A*x */
  energy = dot(x,z,nEqNov);
/* -------------------------------------------------------*/
  timef = getTimeC() - timei;

  if(fPrint){ 
    printf("\tnad         :      %20d\n"  ,nAd);
    printf("\tnad         :      %20d\n"  ,nAdR);
    printf("\tSolver tol  :      %20.2e\n",tol);
    printf(" (PBICGSTAB) solver :\n"
           "\tEquations   =      %20d\n"
           "\tEquations   =      %20d\n"
           "\tIterarions  =      %20d\n"
	         "\tEnergy norm =      %20.12e\n"
	         "\tCPU time(s) =      %20.5lf\n" 
	         ,nEqNov,nEq,j+1,energy,timef);
  }
  
  if(j == maxIt){ 
    printf(" *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
    exit(EXIT_FAILURE);
  }
  
  if(!mpiVar.myId && log)
    fprintf(fLog          
           ,"MPIPBICGSTAB: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);

}
/**********************************************************************/



