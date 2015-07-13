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
void pcg(INT const neq      ,INT const nad  
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
    for(i = 0; i < neq; i++)  
      x[i] = 0.l;
  
      
  matvec(neq,ia,ja,al,ad,x,z);
  
  for(i = 0; i < neq; i++)   {
    r[i] = b[i] - z[i];
    z[i] = r[i] / m[i];
    b[i] = z[i];
  }
  
  d    = dot(r,z,neq);
  conv = tol * sqrt(fabs(d));
/*--------------------------------------------------------*/   
  for(j = 0; j < maxIt; j++)   {
    matvec(neq,ia,ja,al,ad,b,z);
    alpha = d / dot(b,z,neq);
    for(i = 0; i < neq; i++)   {
      x[i] +=  alpha * b[i];
      r[i] -=  alpha * z[i];
      z[i]  = r[i] / m[i];
    }
    beta = dot(r,z,neq)/d;
    for(i = 0; i < neq; i++)   {
      b[i] = z[i] + beta * b[i];
    }
    d = beta * d;
    if (sqrt(fabs(d)) < conv) break;
  }
/* -------------------------------------------------------*/
  matvec(neq,ia,ja,al,ad,x,z);
/*norma de energia = xT*A*x */
  energy = dot(x,z,neq);
/* -------------------------------------------------------*/
  timef = getTimeC() - timei;

  if(fPrint){ 
    printf("\tnad         :      %20d\n"  ,nad);
    printf("\tSolver tol  :      %20.2e\n",tol);
    printf(" (PCG) solver:\n"
           "\tEquations   =      %20d\n"
           "\tIterarions  =      %20d\n"
	         "\tEnergy norm =      %20.12e\n"
	         "\tCPU time(s) =      %20.5lf\n" 
	         ,neq,j+1,energy,timef);
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


