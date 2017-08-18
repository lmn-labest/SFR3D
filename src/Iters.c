#include<Solv.h>
/**********************************************************************
c * Metodos iterativos para solucao de sistemas lineares              *
c * ----------------------------------------------------------------  *
c * simetricos:                                                       *
c * ----------------------------------------------------------------  *
c *                                                                   *
c * PCG - gradiente conjugados com precondicionador diagonal          *
c *                                                                   *
c * ----------------------------------------------------------------  *
c * nao - simetricos:                                                 *
c * ----------------------------------------------------------------  *
c *                                                                   *
c * pbicgstab - gradiente bi - conjugados estabilizados  com          *
c * precondicionador diagonal                                         *
c *                                                                   *
c * pbicgstabl2- gradiente bi-conjugados estabilizados (l=2) com      *
c * precondicionador diagonal                                         *
c *                                                                   *
c *********************************************************************/

/**********************************************************************
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 19/07/2016                                    *
 * -------------------------------------------------------------------*
 * PCG  : metodo do gradiente conjugado com precondiconador diagonal  *
 * (M-1Ax=M-1b) (matriz simentrica)                                   *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 *  neq				-> numero de equacoes                                   *
 *  nad				-> numero de elementos nao nulos fora da diagonal       *
 *  ia				-> estrutura de dados para matriz esparsa A             *
 *  ja				-> estrutura de dados para matriz esparsa A             *
 *  al				-> parte inferior da matriz A                           *
 *  ad				-> diagnal da matriz A                                  *
 *   p				-> precondiconador diagonal                             *
 *   b				-> vetor b (Ax=b)                                       *
 *   x				-> vetor de solucao                                     *
 *   z				-> vetor auxiliar                                       *
 *   r				-> vetor auxiliar                                       *
 *   p				-> vetor auxiliar                                       *
 * tol				-> tolerancia de convergencia                           *
 * newX				-> vetor inicial iniciado com zero                      *
 * fileLog				-> arquivo de log do solver                         *
 * log				-> log de arquivo (true|false)                          *
 * tol				-> tolerancia do solver                                 *
 * maxIt			-> numero maximo de iteracoes                           *
 *  newX			-> true zero o vetor inicial                            *
 * fileLog				-> arquivo de saida do log                          *
 * fileHistLog-> arquivo de log da iteracoes                          *
 *  log       -> escreve o log do solver                              *
 * fPrint			-> saida de informacao na tela                          *
 * fHistLog		-> log das iteracoes                                    *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *  
 * b[]-> inalterado                                                   *
 * ad,al-> inalterado                                                 *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 * ------------------------------------------------------------------ *
*********************************************************************/
void pcg(INT const nEq,INT const nAd
	,INT *RESTRICT ia   ,INT *RESTRICT ja
	,DOUBLE *RESTRICT al,DOUBLE *RESTRICT ad
	,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b
  ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT z
  ,DOUBLE *RESTRICT r ,DOUBLE *RESTRICT p
	,DOUBLE const tol   ,unsigned int maxIt
  ,bool const newX    ,FILE* fileLog   
  ,FILE *fileHistLog	,bool const log 
  ,bool const fHistLog,bool const fPrint
	,void(*matvec)()    ,DOUBLE(*dot)())
{
	unsigned int j, jj;
	INT i;
	DOUBLE alpha, beta, d, conv, xKx, norm_b, norm, norm_r_m, norm_r, di;
	DOUBLE timei, timef;
	timei = getTimeC();

/*... chute inicial*/
	if (newX)
		for (i = 0; i < nEq; i++)
			x[i] = 0.e0;
/*...................................................................*/

/*... conv = tol * |(M-1)b|m = tol(b,M(-1)b) */
	for (i = 0; i < nEq; i++)
		z[i] = b[i] * m[i];

	d = dot(b, z, nEq);
	norm_b = sqrt(fabs(d));
	conv = tol * norm_b;
/*...................................................................*/

/*... Ax0*/
	matvec(nEq, ia, ja, al, ad, x, z);
/*...................................................................*/

/*...*/
	for (i = 0; i < nEq; i++) {
/* ... r0 = b - Ax0*/
		r[i] = b[i] - z[i];
/* ... z0 = (M-1)r0*/
		z[i] = r[i] * m[i];
/* ... p0 = r0*/
		p[i] = z[i];
	}
/*...................................................................*/

/*... (r(0), z(0)) = (r(0), (M-1)r0)*/
	d = dot(r, z, nEq);
/*...................................................................*/

/*...*/
	jj = 1;
	for (j = 0; j < maxIt; j++) {
/*... z = Ap(j)*/
		matvec(nEq, ia, ja, al, ad, p, z);
/*...................................................................*/

/*... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))*/
    alpha = d / dot(z, p, nEq);
/*...................................................................*/

/*...*/
		for (i = 0; i < nEq; i++) {
/*... x(j + 1) = x(j) + alpha*p*/
			x[i] += alpha * p[i];
/*... r(j + 1) = r(j) - alpha*Ap*/
			r[i] -= alpha * z[i];
/*... z = (M-1)r0*/
			z[i] = r[i] * m[i];
		}
/*...................................................................*/

/* ... (r(j + 1), (M - 1)r(j + 1)) = (r(j + 1), z)*/
		di = dot(r, z, nEq);
/*... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) */
		beta = di / d;
/*...................................................................*/

/* ... p(j + 1) = (M - 1)r(j + 1) + beta*p(j) = z + beta*p(j)*/
 		for (i = 0; i < nEq; i++)
		  p[i] = z[i] + beta * p[i];
/*...................................................................*/

/*...*/
		if (fHistLog)
		  fprintf(fileHistLog, "%d %20.9e\n", j, sqrt(fabs(d)) / norm_b);
/*...................................................................*/

/*...*/
		d = di;
		if (sqrt(fabs(d)) < conv) break;
/*...................................................................*/

/*...*/
		if (jj == 2000) {
  		jj = 0;
			printf("PCG: %d %20.9e %20.9e\n", j+1, sqrt(fabs(d)), conv);
		}
		jj++;
/*...................................................................*/
	}
/*...................................................................*/

/*... Energy norm:  x*Kx*/
	matvec(nEq, ia, ja, al, ad, x, z);
/*norma de energia = xT*A*x */
	xKx    = dot(x, z, nEq);
/*...................................................................*/

/*... norm - 2 = || x ||*/
	norm = sqrt(dot(x, x, nEq));
/*...................................................................*/

/*... r = M(-1)(b - Ax) (calculo do residuo explicito)*/
	for (i = 0; i < nEq; i++) {
		r[i] = b[i] - z[i];
		z[i] = r[i] * m[i];
	}
	norm_r_m = dot(r, z, nEq);
	norm_r_m = sqrt(fabs(norm_r_m));
	norm_r = dot(r, r, nEq);
	norm_r = sqrt(norm_r);
	if(fPrint && norm_r_m > 3.16e0*conv)
	  printf("PCG: %20.9e > %20.9e!!\n", norm_r_m, conv);
/*...................................................................*/

  timef = getTimeC() - timei;

  if(fPrint){ 
		printf(" (PCG) solver:\n"
					"\tEquations      = %20d\n"
          "\tnad            = %20d\n"        
          "\tSolver tol     = %20.2e\n"      
          "\tIterarions     = %20d\n"
					"\tx * Kx         = %20.2e\n"
				  "\t|| x ||        = %20.2e\n"
					"\t|| b - Ax ||   = %20.2e\n"
					"\t|| b - Ax || m = %20.2e\n"
	        "\tCPU time(s)    = %20.2lf\n" 
	         ,nEq,nAd,tol,j+1,xKx,norm,norm_r,norm_r_m,timef);
  }
  
  if(j == maxIt){ 
    printf("PCG: *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
    exit(EXIT_FAILURE);
  }
  if(log)
    fprintf(fileLog          
           ,"PCG: tol %20.2e " 
            "iteration %d " 
						"xKx %20.12e "
		        "norma(x*x) %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1, xKx,norm,timef);
}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 28 / 08 / 2016                                *
 * Data de modificaco : 00 / 00 / 0000                                *
 * ------------------------------------------------------------------ *
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
 *   p  -> precondiconador diagonal                                   *
 *   b  -> vetor b (Ax=b)                                             *
 *   x  -> vetor de solucao                                           *
 *   z  -> vetor auxiliar                                             *
 *   r  -> vetor auxiliar                                             *
 * newX -> vetor inicial iniciado com zero                            *
 * fileLog -> arquivo de log do solver                                *
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
 * ad,al-> inalterado                                                 *
 * -------------------------------------------------------------------*
*********************************************************************/
void mpiPcg(INT const nEq   ,INT const nEqNov
        ,INT const nAd      ,INT const nAdR  
        ,INT *RESTRICT ia   ,INT *RESTRICT ja
        ,DOUBLE *RESTRICT al,DOUBLE *RESTRICT ad
        ,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b 
        ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT z 
        ,DOUBLE *RESTRICT r 
        ,DOUBLE const tol   ,unsigned int maxIt
        ,bool const newX    ,FILE* fileLog  
        ,bool const log     ,bool const fPrint
        ,Interface *iNeq                      
        ,void(*matvec)()    ,DOUBLE(*dot)())
{
#ifdef _MPICH_
  unsigned int j;
  INT i;
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
	         "\tCPU time(s) =      %20.2lf\n" 
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
    fprintf(fileLog          
           ,"MPIPCG: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);
#endif
}
/**********************************************************************/

/**********************************************************************
* Data de criacao    : 28/08/2016                                    *
* Data de modificaco : 00/00/0000                                    *
* -------------------------------------------------------------------*
* PCGOMP: metodo do gradiente conjugado com precondiconador diagonal *
* (M-1Ax=M-1b) (matriz simentrica)                                   *
* -------------------------------------------------------------------*
* Parametros de Entrada:                                             *
* -------------------------------------------------------------------*
*  neq			  -> numero de equacoes                                  *
*  nad			  -> numero de elementos nao nulos fora da diagonal      *
*  ia				  -> estrutura de dados para matriz esparsa A            *
*  ja				  -> estrutura de dados para matriz esparsa A            *
*  a				  -> parte inferior da matriz A                          *
*  ad				  -> diagnal da matriz A                                 *
*   p				  -> precondiconador diagonal                            *
*   b				  -> vetor b (Ax=b)                                      *
*   x				  -> vetor de solucao                                    *
*   z				  -> vetor auxiliar                                      *
*   r				  -> vetor auxiliar                                      *
*   p				  -> vetor auxiliar                                      *
* tol				  -> tolerancia de convergencia                          *
* newX			  -> vetor inicial iniciado com zero                     *
* fileLog			  -> arquivo de log do solver                          *
* log				  -> log de arquivo (true|false)                         *
* tol				  -> tolerancia do solver                                *
* maxIt			  -> numero maximo de iteracoes                          *
*  newX			  -> true zero o vetor inicial                           *
* fileLog				-> arquivo de saida do log                           *
* fileHistLog -> arquivo de log da iteracoes                         *
*  log        -> escreve o log do solver                             *
* fPrint			-> saida de informacao na tela                         *
* fHistLog		-> log das iteracoes                                   *
* bOmp        -> Openmp                                              *
* -------------------------------------------------------------------*
* Parametros de Saida:                                               *
* -------------------------------------------------------------------*
* x[]-> atualizado                                                   *
* b[]-> inalterado                                                   *
* ad,al-> inalterado                                                 *
* ------------------------------------------------------------------ *
* OBS:                                                               *
* ------------------------------------------------------------------ *
*********************************************************************/
void pcgOmp(INT const nEq, INT const nAd
            , INT *RESTRICT ia, INT *RESTRICT ja
            , DOUBLE *RESTRICT a, DOUBLE *RESTRICT ad
            , DOUBLE *RESTRICT m, DOUBLE *RESTRICT b
            , DOUBLE *RESTRICT x, DOUBLE *RESTRICT z
            , DOUBLE *RESTRICT r, DOUBLE *RESTRICT p
            , DOUBLE const tol, unsigned int maxIt
            , bool const newX, FILE* fileLog
            , FILE *fileHistLog, bool const log
            , bool const fHistLog, bool const fPrint
            , BufferOmp *bOmp
            , void(*matvec)(), DOUBLE(*dot)())
{
  short nThreads = ompVar.nThreadsSolver;
  unsigned int j, jj;
  INT i;
  DOUBLE alpha, beta, d, conv, xKx, norm_b, norm, norm_r_m, norm_r, di,tmp;
  DOUBLE timei, timef;
  timei = getTimeC();

#pragma omp parallel default(none) num_threads(nThreads)\
   private(i,j,jj,conv,norm_b,tmp,norm,d,di,alpha,beta)\
   shared(ia,ja,a,ad,m,b,x,z,r,p,tol,maxIt,newX,fileLog,xKx,norm_r_m,norm_r\
         ,fileHistLog,log,fHistLog,fPrint,matvec,dot,bOmp,nThreads)
  {
/*... chute inicial*/
    if (newX)
#pragma omp for
      for (i = 0; i < nEq; i++)
        x[i] = 0.e0;
/*...................................................................*/

/*... conv = tol * |(M-1)b|m = tol(b,M(-1)b) */
#pragma omp for
    for (i = 0; i < nEq; i++)
      z[i] = b[i] * m[i];

    d = dot(b, z, nEq);
    norm_b = sqrt(fabs(d));
    conv = tol * norm_b;
/*...................................................................*/

/*... Ax0*/
    matvec(nEq
           , ia, ja
           , a, ad
           , x, z
           , bOmp->thBegin, bOmp->thEnd
           , bOmp->thHeight, bOmp->thY
           , nThreads);
/*...................................................................*/

/*...*/
#pragma omp for
    for (i = 0; i < nEq; i++) {
/* ... r0 = b - Ax0*/
      r[i] = b[i] - z[i];
/* ... z0 = (M-1)r0*/
      z[i] = r[i] * m[i];
/* ... p0 = r0*/
      p[i] = z[i];
    }
/*...................................................................*/

/*... (r(0), z(0)) = (r(0), (M-1)r0)*/
    d = dot(r, z, nEq);
/*...................................................................*/

/*...*/
    jj = 1;
    for (j = 0; j < maxIt; j++) {
/*... z = Ap(j)*/
      matvec(nEq
             , ia, ja
             , a, ad
             , p, z
             , bOmp->thBegin, bOmp->thEnd
             , bOmp->thHeight, bOmp->thY
             , nThreads);
/*...................................................................*/

/*... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))*/
      alpha = d / dot(z, p, nEq);
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i < nEq; i++) {
/*... x(j + 1) = x(j) + alpha*p*/
        x[i] += alpha * p[i];
/*... r(j + 1) = r(j) - alpha*Ap*/
        r[i] -= alpha * z[i];
/*... z = (M-1)r0*/
        z[i] = r[i] * m[i];
      }
/*...................................................................*/

/* ... (r(j + 1), (M - 1)r(j + 1)) = (r(j + 1), z)*/
      di = dot(r, z, nEq);
/*... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) */
      beta = di / d;
/*...................................................................*/

/* ... p(j + 1) = (M - 1)r(j + 1) + beta*p(j) = z + beta*p(j)*/
#pragma omp for
      for (i = 0; i < nEq; i++)
        p[i] = z[i] + beta * p[i];
/*...................................................................*/

/*...*/
#pragma omp master
      if (fHistLog)
        fprintf(fileHistLog, "%d %20.9e\n", j, sqrt(fabs(d)) / norm_b);
/*...................................................................*/

/*...*/
      d = di;
      if (sqrt(fabs(d)) < conv) break;
/*...................................................................*/

/*...*/
      if (jj == 2000) {
        jj = 0;
#pragma omp master
        printf("PCGOMP: %d %20.9e %20.9e\n", j + 1, sqrt(fabs(d)), conv);
      }
      jj++;
/*...................................................................*/
    }
/*...................................................................*/

/*... Energy norm:  x*Kx*/
    matvec(nEq
           , ia, ja
           , a, ad
           , x, z
           , bOmp->thBegin, bOmp->thEnd
           , bOmp->thHeight, bOmp->thY
           , nThreads);
/*norma de energia = xT*A*x */
    tmp = dot(x, z, nEq);
#pragma omp single
    xKx = tmp;
/*...................................................................*/

/*... norm - 2 = || x ||*/
    tmp  = sqrt(dot(x, x, nEq));
#pragma omp single
    norm = tmp;
/*...................................................................*/

/*... r = M(-1)(b - Ax) (calculo do residuo explicito)*/
#pragma omp for
    for (i = 0; i < nEq; i++) {
      r[i] = b[i] - z[i];
      z[i] = r[i] * m[i];
    }

    tmp = dot(r, z, nEq);
#pragma omp single
    norm_r_m = sqrt(fabs(tmp));
    
    tmp = dot(r, r, nEq);
#pragma omp single
    norm_r = sqrt(tmp);

#pragma omp single
    if (fPrint && norm_r_m > 3.16e0*conv)
      printf("PCGOMP: %20.9e > %20.9e!!\n", norm_r_m, conv);
/*...................................................................*/
  }
/*...................................................................*/

  timef = getTimeC() - timei;

  if (fPrint) {
    printf(" (PCGOMP) solver:\n"
           "\tEquations      = %20d\n"
           "\tnad            = %20d\n"
           "\tSolver tol     = %20.2e\n"
           "\tIterarions     = %20d\n"
           "\tx * Kx         = %20.2e\n"
           "\t|| x ||        = %20.2e\n"
           "\t|| b - Ax ||   = %20.2e\n"
           "\t|| b - Ax || m = %20.2e\n"
           "\tCPU time(s)    = %20.2lf\n"
           , nEq, nAd, tol, j + 1, xKx, norm, norm_r, norm_r_m, timef);
  }

  if (j == maxIt) {
    printf("PCGOMP: *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n", maxIt);
    exit(EXIT_FAILURE);
  }
  if (log)
    fprintf(fileLog
            , "PCGOMP: tol %20.2e "
            "iteration %d "
            "xKx %20.12e "
            "norma(x*x) %20.12e "
            "time %20.5lf\n"
            , tol, j + 1, xKx, norm, timef);
}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 19/07/2016                                    *
 * -------------------------------------------------------------------*
 * PBICGSTAB : metodo do gradiente conjugado bi-ortaganalizado        *
 * com precondiconador diagonal (matriz geral)                        *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 *  neq       -> numero de equacoes                                   *
 *  nad       -> numero de elementos nao nulos fora da diagonal       *
 *  ia        -> estrutura de dados para matriz esparsa A             *
 *  ja        -> estrutura de dados para matriz esparsa A             *
 *  a         -> coef fora da diagonal principal                      *
 *  ad        -> diagnal da matriz A                                  *
 *   p        -> precondiconador diagonal                             *
 *   b        -> vetor b (Ax=b)                                       *
 *   x        -> vetor de solucao                                     *
 *   t        -> vetor auxiliar                                       *
 *   v        -> vetor auxiliar                                       *
 *   r        -> vetor auxiliar                                       *
 *   p        -> vetor auxiliar                                       *
 *   z        -> vetor auxiliar                                       *
 * r0          > vetor auxiliar                                       *
 * newX       -> vetor inicial iniciado com zero                      *
 * fileLog       -> arquivo de log do solver                          *
 * log        -> log de arquivo (true|false)                          *
 * tol        -> tolerancia do solver                                 *
 * maxIt			-> numero maximo de iteracoes                           *
 *  newX			-> true zero o vetor inicial                            *
 * fileLog				-> arquivo de saida do log                          *
 * fileHistLog-> arquivo de log da iteracoes                          *
 *  log       -> escreve o log do solver                              *
 * fPrint			-> saida de informacao na tela                          *
 * fHistLog		-> log das iteracoes                                    *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *  
 * b[]-> inalterado                                                   *
 * ad,a-> inalterado                                                  *
 * -------------------------------------------------------------------*
 * OBS:                                                               *
 * A(M-1)y=b precondicionador a direita                               *
 * -------------------------------------------------------------------*
***********************************************************************/
void pbicgstab(INT const nEq  ,INT const nAd
          ,INT *RESTRICT ia   ,INT *RESTRICT ja
          ,DOUBLE *RESTRICT al,DOUBLE *RESTRICT ad
          ,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b 
          ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT t
          ,DOUBLE *RESTRICT v ,DOUBLE *RESTRICT r
          ,DOUBLE *RESTRICT p ,DOUBLE *RESTRICT z
          ,DOUBLE *RESTRICT r0
					,DOUBLE const tol   ,unsigned int maxIt
          ,bool const newX    ,FILE* fileLog 
          ,FILE *fileHistLog  ,bool const log     
          ,bool const fHistLog,bool const fPrint
          ,void(*matvec)()    ,DOUBLE(*dot)())
{
  unsigned int j,jj;
  INT i;
  DOUBLE alpha,beta,d,conv,xKx,norm_r,norm,norm_b,w,rr0;
  DOUBLE timei,timef;

  timei = getTimeC();
  
/*... chute inicial*/
  if(newX)  
    for (i = 0; i < nEq; i++)  
      x[i] = 0.e0;
/*...................................................................*/

/*... conv = tol * |b|*/	
  d        = dot(b,b,nEq);
	norm_b   = sqrt(d);
	conv     = tol*norm_b;
//breaktol = btol*sqrt(d);
/*...................................................................*/

/*... Ax0*/  
  matvec(nEq,ia,ja,al,ad,x,z);
/*...................................................................*/

/*...*/
  for(i = 0; i < nEq; i++)   {
/*... r0 = b - Ax0*/
    r0[i] = b[i] - z[i];
/*... p = r0*/
    p[i] = r0[i];
/*... r = r0*/
    r[i] = r0[i];
/*... z = M(-1)p*/
    z[i] = p[i]*m[i];
  }
/*...................................................................*/

/*...*/
	jj = 1;
  for(j = 0; j < maxIt; j++)   {
/* ... v = Az(j) = AM(-1)p(j)*/
    matvec(nEq,ia,ja,al,ad,z,v);
/*...................................................................*/

/*... alpha = (r(j), r0) / (AM(-1)p(j), r0))*/
    rr0   = dot(r,r0,nEq);
    alpha = rr0/dot(v,r0,nEq);
/*...................................................................*/

/*...*/
    for(i = 0; i < nEq; i++)   {
/*... x(j + 1) = x(j) + alpha*M(-1)p*/
      x[i] += alpha * z[i];
/*... s(j) = r(j) - alpha*AM(-1)p(j)*/
      r[i] -= alpha * v[i];
/*... z = M(-1)s*/
      z[i]  = r[i] * m[i];
    }
/*...................................................................*/

/*... (s, s)*/
		d = dot(r, r, nEq);
/*...*/	
		if (sqrt(d) < conv) break;
/*...................................................................*/

/*... t = Az = AM(-1)s(j)*/
    matvec(nEq,ia,ja,al,ad,z,t);
/*...................................................................*/

/*... w = (AM(-1)s(j), s(j)) / (AM(-1)s(j), AM(-1)s(j))*/
    w = dot(t,r,nEq)/ dot(t,t,nEq);
/*...................................................................*/

/*...*/
    for(i = 0; i < nEq; i++)   {
/* ... x(j + 1) = x(j) + w*M(-1)s*/
      x[i] += w*z[i];
/* ... r(j + 1) = s(j) - w*AM(-1)s*/
      r[i] -= w*t[i];
    }
/*...................................................................*/

/*... (r, r)*/
    d = dot(r,r,nEq);
    if(sqrt(d) < conv) break;
/*...................................................................*/

/*...*/
		if (fHistLog)
			fprintf(fileHistLog, "%d %20.9e\n", j, sqrt(fabs(d)) / norm_b);
/*...................................................................*/

/*... beta = (r(j + 1), r0) / (r(j), r0)) * (alpha / w)*/
    beta = (dot(r,r0,nEq)/rr0)*(alpha/w);
/*...................................................................*/

/*...*/
    for(i = 0; i < nEq; i++)   {
/*... p(j + 1) = r(i) + beta*(p(j) - w*v(i))*/
      p[i]  = r[i] + beta*(p[i]-w*v[i]);
/*... z = M(-1)p*/
      z[i]  = p[i]*m[i];
    }
/*...................................................................*/

/*...*/
		if (jj == 2000) {
			jj = 0;
			printf("BICGSATB: %d %20.9e %20.9e\n"
            ,j + 1, sqrt(fabs(d)), conv);
		}
		jj++;
/*...................................................................*/
  }
/*...................................................................*/

/*... Energy norm:  x*Kx*/
	matvec(nEq, ia, ja, al, ad, x, z);
/*norma de energia = xT*A*x */
	xKx = dot(x, z, nEq);
/*...................................................................*/

/*... norm - 2 = || x ||*/
	norm = sqrt(dot(x, x, nEq));
/*...................................................................*/

/*... r = (b - Ax) (calculo do residuo explicito)*/
	for (i = 0; i < nEq; i++) {
		r[i] = b[i] - z[i];
	}
	norm_r = dot(r, r, nEq);
	norm_r = sqrt(norm_r);
	if(fPrint && norm_r > 3.16e0*conv)
		printf("BICGSATB: %20.9e > %20.9e!!\n", norm_r, conv);
/*...................................................................*/
  timef = getTimeC() - timei;

  if(fPrint){ 
		printf(" (BICGSATB) solver:\n"
			"\tEquations      = %20d\n"
			"\tnad            = %20d\n"
			"\tSolver tol     = %20.2e\n"
			"\tIterarions     = %20d\n"
			"\tx * Kx         = %20.2e\n"
			"\t|| x ||        = %20.2e\n"
			"\t|| b - Ax ||   = %20.2e\n"
			"\tCPU time(s)    = %20.2lf\n"
			, nEq, nAd, tol, j + 1, xKx, norm, norm_r, timef);
  }
  
  if(j == maxIt){ 
    printf(" PBICGSTAB *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
    exit(EXIT_FAILURE);
  }
  if(log)
    fprintf(fileLog          
           ,"PBICGSTAB: tol %20.2e " 
						"iteration %d "
						"xKx %20.12e "
						"norma(x*x) %20.12e "
						"time %20.5lf\n"
						, tol, j + 1, xKx, norm, timef);


}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 28 / 08 / 2016                                *
 * Data de modificaco : 00 / 00 / 0000                                *
 * ------------------------------------------------------------------ *
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
 *  a         -> coef fora da diagonal principal                      *
 *  ad  -> diagnal da matriz A                                        *
 *   p  -> precondiconador diagonal                                   *
 *   b  -> vetor b (Ax=b)                                             *
 *   x  -> vetor de solucao                                           *
 *   t  -> vetor auxiliar                                             *
 *   v  -> vetor auxiliar                                             *
 *   r  -> vetor auxiliar                                             *
 *   p  -> vetor auxiliar                                             *
 *   z  -> vetor auxiliar                                             *
 * newX -> vetor inicial iniciado com zero                            *
 * fileLog -> arquivo de log do solver                                *
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
          ,INT *RESTRICT ia    ,INT *RESTRICT ja
          ,DOUBLE *RESTRICT al ,DOUBLE *RESTRICT ad
          ,DOUBLE *RESTRICT m  ,DOUBLE *RESTRICT b 
          ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT t
          ,DOUBLE *RESTRICT v  ,DOUBLE *RESTRICT r
          ,DOUBLE *RESTRICT p  ,DOUBLE *RESTRICT z 
          ,DOUBLE const tol
          ,unsigned int maxIt  ,bool const newX          
          ,FILE* fileLog          ,bool const log
          ,bool const fPrint   ,Interface *iNeq    
          ,void(*matvec)()     ,DOUBLE(*dot)())
{
	unsigned int j;
	INT i;
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
	         "\tCPU time(s) =      %20.2lf\n" 
	         ,nEqNov,nEq,j+1,energy,timef);
  }
  
  if(j == maxIt){ 
    printf(" *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
    exit(EXIT_FAILURE);
  }
  
  if(!mpiVar.myId && log)
    fprintf(fileLog          
           ,"MPIPBICGSTAB: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);

}
/**********************************************************************/

/**********************************************************************
* Data de criacao    : 28 / 08 / 2016                                 *
* Data de modificaco : 00 / 00 / 0000                                 *
* ------------------------------------------------------------------- *
* PBICGSTABOMP : metodo do gradiente conjugado bi - ortaganalizado    *
* com precondiconador diagonal(matriz geral)                          *
* ------------------------------------------------------------------- *
* Parametros de Entrada :                                             *
* ------------------------------------------------------------------- *
* neq         ->numero de equacoes                                    *
* nad         ->numero de elementos nao nulos fora da diagonal        *
* ia          ->estrutura de dados para matriz esparsa A              *
* ja          ->estrutura de dados para matriz esparsa A              *
* a           ->parte fora da diagonal principal                      *
* ad          ->diagnal da matriz A                                   *
* p           ->precondiconador diagonal                              *
* b           ->vetor b(Ax = b)                                       *
* x           ->vetor de solucao                                      *
* t           ->vetor auxiliar                                        *
* v           ->vetor auxiliar                                        *
* r           ->vetor auxiliar                                        *
* p           ->vetor auxiliar                                        *
* z           ->vetor auxiliar                                        *
* r0          ->vetor auxiliar                                        *
* newX        ->vetor inicial iniciado com zero                       *
* fileLog        ->arquivo de log do solver                           *
* log         ->log de arquivo(true | false)                          *
* tol         ->tolerancia do solver                                  *
* maxIt       ->numero maximo de iteracoes                            *
* newX			  -> true zero o vetor inicial                            *
* fileLog        ->arquivo de saida do log                            *
* fileHistLog ->arquivo de log da iteracoes                           *
* log         ->escreve o log do solver                               *
* fPrint      ->saida de informacao na tela                           *
* fHistLog    ->log das iteracoes                                     *
* bOmp        -> Openmp                                               *
* ------------------------------------------------------------------- *
* Parametros de Saida :                                               *
* ------------------------------------------------------------------- *
* x[]->atualizado                                                     *
* b[]->inalterado                                                     *
* ad, a->inalterado                                                   *
* ------------------------------------------------------------------- *
* OBS:                                                                *
* A(M - 1)y = b precondicionador a direita                            *
* ------------------------------------------------------------------- *
***********************************************************************/
void pbicgstabOmp(INT const nEq    ,INT const nAd
               ,INT *RESTRICT ia   ,INT *RESTRICT ja
               ,DOUBLE *RESTRICT a ,DOUBLE *RESTRICT ad
               ,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b
               ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT t
               ,DOUBLE *RESTRICT v ,DOUBLE *RESTRICT r
               ,DOUBLE *RESTRICT p ,DOUBLE *RESTRICT z
               ,DOUBLE *RESTRICT r0
               ,DOUBLE const tol   ,unsigned int maxIt
               ,bool const newX    ,FILE* fileLog
               ,FILE *fileHistLog  ,bool const log
               ,bool const fHistLog,bool const fPrint
               ,BufferOmp *bOmp
               ,void(*matvec)()    ,DOUBLE(*dot)())
{
  short nThreads = ompVar.nThreadsSolver;
  unsigned int j, jj,js;
  INT i;
  DOUBLE alpha,beta,d,conv,xKx,norm_r,norm,norm_b,w,rr0,tmp;
  DOUBLE timei, timef;

  timei = getTimeC();

#pragma omp parallel default(none) num_threads(nThreads)\
   private(i,j,jj,conv,norm_b,tmp,d,alpha,beta,rr0,w)\
   shared(js,ia,ja,a,ad,m,b,x,z,r,p,r0,v,t,tol,maxIt,newX,fileLog\
         ,xKx,norm,norm_r\
         ,fileHistLog,log,fHistLog,fPrint,matvec,dot,bOmp,nThreads)
  {
/*... chute inicial*/
    if (newX)
#pragma omp for
      for (i = 0; i < nEq; i++)
        x[i] = 0.e0;
/*...................................................................*/

/*... conv = tol * |b|*/
    d      = dot(b, b, nEq);
    norm_b = sqrt(d);
    conv   = tol*norm_b;
  //breaktol = btol*sqrt(d);
/*...................................................................*/

/*... Ax0*/
    matvec(nEq
          ,ia            ,ja 
          ,a             ,a
          ,x             ,z
          ,bOmp->thBegin ,bOmp->thEnd
          ,bOmp->thHeight,bOmp->thY
          ,nThreads);
/*...................................................................*/

/*...*/
#pragma omp for
    for (i = 0; i < nEq; i++) {
/*... r0 = b - Ax0*/
      r0[i] = b[i] - z[i];
/*... p = r0*/
      p[i] = r0[i];
/*... r = r0*/
      r[i] = r0[i];
/*... z = M(-1)p*/
      z[i] = p[i]*m[i];
    }
/*...................................................................*/

/*...*/
    jj = 1;
    for (j = 0; j < maxIt; j++) {
/* ... v = Az(j) = AM(-1)p(j)*/
      matvec(nEq
            ,ia            ,ja
            ,a             ,ad
            ,z             ,v
            ,bOmp->thBegin ,bOmp->thEnd
            ,bOmp->thHeight,bOmp->thY
            ,nThreads);
/*...................................................................*/

/*... alpha = (r(j), r0) / (AM(-1)p(j), r0))*/
      rr0   = dot(r, r0, nEq);
      alpha = rr0/ dot(v, r0, nEq);
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i < nEq; i++) {
/*... x(j + 1) = x(j) + alpha*M(-1)p*/
        x[i] += alpha * z[i];
/*... s(j) = r(j) - alpha*AM(-1)p(j)*/
        r[i] -= alpha * v[i];
/*... z = M(-1)s*/
        z[i] = r[i]*m[i];
      }
/*...................................................................*/

/*... (s, s)*/
      d = dot(r, r, nEq);
/*...*/
      if (sqrt(d) < conv) break;
/*...................................................................*/

/*... t = Az = AM(-1)s(j)*/
      matvec(nEq
            ,ia            ,ja
            ,a             ,ad
            ,z             ,t
            ,bOmp->thBegin ,bOmp->thEnd
            ,bOmp->thHeight,bOmp->thY
            ,nThreads);
/*...................................................................*/

/*... w = (AM(-1)s(j), s(j)) / (AM(-1)s(j), AM(-1)s(j))*/
      w = dot(t,r,nEq)/dot(t,t,nEq);
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i < nEq; i++) {
/* ... x(j + 1) = x(j) + w*M(-1)s*/
        x[i] += w*z[i];
/* ... r(j + 1) = s(j) - w*AM(-1)s*/
        r[i] -= w*t[i];
      }
/*...................................................................*/

/*... (r, r)*/
      d = dot(r,r,nEq);
      if (sqrt(d) < conv) break;
/*...................................................................*/

/*...*/
#pragma omp master
      if (fHistLog)
        fprintf(fileHistLog, "%d %20.9e\n", j, sqrt(fabs(d)) / norm_b);
/*...................................................................*/

/*... beta = (r(j + 1), r0) / (r(j), r0)) * (alpha / w)*/
      beta = (dot(r,r0,nEq)/rr0)*(alpha/w);
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i < nEq; i++) {
/*... p(j + 1) = r(i) + beta*(p(j) - w*v(i))*/
        p[i] = r[i] + beta*(p[i] - w*v[i]);
/*... z = M(-1)p*/
        z[i] = p[i]*m[i];
      }
/*...................................................................*/

/*...*/
      if (jj == 2000) {
        jj = 0;
#pragma omp master
        printf("BICGSATBOMP: %d %20.9e %20.9e\n"
               , j + 1, sqrt(fabs(d)), conv);
      }
      jj++;
/*...................................................................*/
    }
/*...................................................................*/

#pragma omp single
    js = j;

/*... Energy norm:  x*Kx*/
    matvec(nEq
          ,ia            ,ja
          ,a             ,ad
          ,x             ,z
          ,bOmp->thBegin ,bOmp->thEnd
          ,bOmp->thHeight,bOmp->thY
          ,nThreads);
/*norma de energia = xT*A*x */
    tmp = dot(x, z, nEq);
#pragma omp master
    xKx = tmp; 
/*...................................................................*/

/*... norm - 2 = || x ||*/
    tmp = sqrt(dot(x, x, nEq));
#pragma omp master
    norm = tmp;
/*...................................................................*/

/*... r = (b - Ax) (calculo do residuo explicito)*/
#pragma omp for
    for (i = 0; i < nEq; i++) {
      r[i] = b[i] - z[i];
    }

    tmp = dot(r, r, nEq);

#pragma omp single
    {
      norm_r = sqrt(tmp);
      if(fPrint && norm_r > 3.16e0*conv)
      printf("BICGSATB: %20.9e > %20.9e!!\n", norm_r, conv);
    }
/*...................................................................*/
  }
/*...................................................................*/

  timef = getTimeC() - timei;

  if (fPrint) {
    printf(" (BICGSATBOMP) solver:\n"
           "\tEquations      = %20d\n"
           "\tnad            = %20d\n"
           "\tSolver tol     = %20.2e\n"
           "\tIterarions     = %20d\n"
           "\tx * Kx         = %20.2e\n"
           "\t|| x ||        = %20.2e\n"
           "\t|| b - Ax ||   = %20.2e\n"
           "\tCPU time(s)    = %20.2lf\n"
           , nEq, nAd, tol, js + 1, xKx, norm, norm_r, timef);
  }

  if (js == maxIt ) {
    printf(" PBICGSTABOMP *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n", maxIt);
    exit(EXIT_FAILURE);
  }
  if (log)
    fprintf(fileLog
            , "PBICGSTABOMP: tol %20.2e "
            "iteration %d "
            "xKx %20.12e "
            "norma(x*x) %20.12e "
            "time %20.5lf\n"
            , tol, j + 1, xKx, norm, timef);


}
/**********************************************************************/

/**********************************************************************
* Data de criacao    : 01/09/2016                                    *
* Data de modificaco : 00/00/0000                                    *
* -------------------------------------------------------------------*
* PBICGSTABl2 : metodo do gradiente conjugado bi-ortaganalizado      *
* com precondiconador diagonal (matriz geral)                        *
* -------------------------------------------------------------------*
* Parametros de Entrada:                                             *
* -------------------------------------------------------------------*
*  neq       -> numero de equacoes                                   *
*  nad       -> numero de elementos nao nulos fora da diagonal       *
*  ia        -> estrutura de dados para matriz esparsa A             *
*  ja        -> estrutura de dados para matriz esparsa A             *
*  a         -> coef fora da diagonal principal                      *
*  ad        -> diagnal da matriz A                                  *
*   p        -> precondiconador diagonal                             *
*   b        -> vetor b (Ax=b)                                       *
*   x        -> vetor de solucao                                     *
*   t        -> vetor auxiliar                                       *
*   v        -> vetor auxiliar                                       *
*   r        -> vetor auxiliar                                       *
*   u        -> vetor auxiliar                                       *
*   r0       -> vetor auxiliar                                       *
*   w        -> vetor auxiliar                                       *
*   s        -> vetor auxiliar                                       *
*   p        -> vetor auxiliar                                       *
*   h        -> vetor auxiliar                                       *
*  z          > vetor auxiliar                                       *
* newX       -> vetor inicial iniciado com zero                      *
* fileLog       -> arquivo de log do solver                          *
* log        -> log de arquivo (true|false)                          *
* tol        -> tolerancia do solver                                 *
* maxIt			-> numero maximo de iteracoes                            *
*  newX			-> true zero o vetor inicial                             *
* fileLog				-> arquivo de saida do log                           *
* fileHistLog-> arquivo de log da iteracoes                          *
*  log       -> escreve o log do solver                              *
* fPrint			-> saida de informacao na tela                         *
* fHistLog		-> log das iteracoes                                   *
* -------------------------------------------------------------------*
* Parametros de Saida:                                               *
* -------------------------------------------------------------------*
* x[]-> atualizado                                                   *
* b[]-> inalterado                                                   *
* ad,a-> inalterado                                                  *
* -------------------------------------------------------------------*
* OBS:                                                               *
c * Versão do livro Iterative krylov Method for large linear Systems *
c * A(M-1)y=b precondicionador a direita                             *
* -------------------------------------------------------------------*
***********************************************************************/
void pbicgstabl2(INT const nEq      ,INT const nAd
               ,INT *RESTRICT ia   ,INT *RESTRICT ja
               ,DOUBLE *RESTRICT a ,DOUBLE *RESTRICT ad
               ,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b
               ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT t
               ,DOUBLE *RESTRICT v ,DOUBLE *RESTRICT r
               ,DOUBLE *RESTRICT u ,DOUBLE *RESTRICT r0
               ,DOUBLE *RESTRICT w ,DOUBLE *RESTRICT s
               ,DOUBLE *RESTRICT p ,DOUBLE *RESTRICT h
               ,DOUBLE *RESTRICT z 
               ,DOUBLE const tol   ,unsigned int maxIt
               ,bool const newX    ,FILE* fileLog
               ,FILE *fileHistLog  ,bool const log
               ,bool const fHistLog,bool const fPrint
               ,void(*matvec)()    ,DOUBLE(*dot)())
{
  unsigned int j, jj;
  INT i;
  DOUBLE alpha,beta,d,conv,xKx,norm_r,norm,norm_b,rr0,rr1
        ,omega1,omega2,gamma,mi,ni,tau;
  DOUBLE timei, timef;

  timei = getTimeC();

/*... chute inicial*/
  if (newX)
    for (i = 0; i < nEq; i++)
      x[i] = 0.e0;
/*...................................................................*/

/*... conv = tol * |b|*/
  d       = dot(b, b, nEq);
  norm_b = sqrt(d);
  conv   = tol*norm_b;
  //breaktol = btol*sqrt(d);
/*...................................................................*/

/*... Ax0*/
  matvec(nEq, ia, ja, a, ad, x, t);
/*...................................................................*/

/*...*/
  for (i = 0; i < nEq; i++) {
/*... r0 = b - Ax0*/
    r0[i] = b[i] - t[i];
/*... r = r0*/
    r[i] = r0[i];
 /*... z = M(-1)p*/
    u[i] = 0.e0;
  }
/*...................................................................*/

/*...*/
  rr0    = 1.e0;
  alpha  = 0.e0;
  omega2 = 1.e0; 
/*...................................................................*/

/*...*/
  jj = 1;
  for (j = 0; j < maxIt; j++) {
/*... rr0 = -w2*rr0*/
    rr0 = -omega2*rr0;
/*...................................................................*/

/*... even BiCG step*/
/*... rr1 = (r,r0)*/
      rr1 = dot(r,r0,nEq);
/*... beta = alpha * rr1 / rr0*/
      beta = alpha*rr1/rr0;
/*... rr0 = rr1*/
      rr0 = rr1;
/*...................................................................*/

/*...*/
      for(i=0;i<nEq;i++){
/*... u(j) = r(j) - beta*u(i)*/
        u[i] = r[i] - beta*u[i];
/*... p = M(-1)u*/
        z[i] = u[i]*m[i];
      }
/*...................................................................*/

/*... v = Ap(j) = AM(-1)u*/
    matvec(nEq,ia,ja,a,ad,z,v);
/*...................................................................*/

/*... alpha = (v, r0) / (AM(-1)u(j), r0))*/
    gamma = dot(v, r0, nEq);
    alpha = rr0 / gamma;
/*...................................................................*/

/*...*/
    for(i=0;i<nEq;i++) {
/*... s(j) = r(j) - alpha*AM(-1)u(j)*/
      r[i] -= alpha * v[i];
/*... p = M(-1)r*/
      p[i] = r[i]*m[i];
    }
/*...................................................................*/

/*... s = Ap(j) = AM(-1)r*/
    matvec(nEq,ia,ja,a,ad,p,s);
/*...................................................................*/

/*...*/
    for(i=0;i<nEq;i++){
/*... x(j) = x(j) + alpha*(M-1)u*/
      x[i] += alpha*z[i];
    }
/*...................................................................*/

/*... odd BiCG step*/
/*... rr1 = (s, r0)*/
    rr1 = dot(s,r0,nEq);
/*... beta = alpha*rr1 / rr0*/
    beta = alpha*rr1/rr0;
/*... rr0 = rr1*/
    rr0 = rr1;
/*...................................................................*/
        
/*...*/
    for(i=0;i<nEq;i++){
/*... v(j) = s(j) - beta*v*/
      v[i] = s[i] - beta*v[i];
/*... p = M(-1)v*/
      p[i] = v[i]*m[i];
    }
/*...................................................................*/

/*... s = Ap(j) = AM(-1)v*/
    matvec(nEq,ia,ja,a,ad,p,w);
/*...................................................................*/

/*... gamma = (w,r0)*/
    gamma = dot(w,r0,nEq);
    alpha = rr0 / gamma;
/*...................................................................*/

/*...*/
      for(i=0;i<nEq;i++){
/*... u(j) = r(j) - beta*u(i)*/
        u[i] = r[i] - beta*u[i];
/*... r(j) = r(j) - alpha*v(i)*/
        r[i] -= alpha*v[i];
/*... s(j) = s(j) - alpha*w(i)*/
        s[i] -= alpha * w[i];
/*... p = M(-1)s*/
        p[i] = s[i]*m[i];
/*... z = M(-1)u*/
        z[i] = u[i]*m[i];
/*... h = M(-1)r*/
        h[i] = r[i]*m[i];
      }
/*...................................................................*/

/*... t = Ap = AM(-1)s(j)*/
      matvec(nEq,ia,ja,a,ad,p,t);
/*...................................................................*/

/*... GCR(2) - part*/
/*... w1 = (r, s)*/
      omega1 = dot(r,s,nEq);
/*... mi = (s, s)*/
      mi = dot(s,s,nEq);
/*... ni = (s, t)*/
      ni = dot(s,t,nEq);
/*... tau = (t, t)*/
      tau = dot(t,t,nEq);
/*... w2*/
      omega2 = dot(r,t,nEq);
/*... tau = tau - ni*ni / mi*/
      tau = tau - ni*ni / mi;
/*... w2 = (w2 - ni*w1 / mi) / tau*/
      omega2 = (omega2 - ni*omega1 / mi) / tau;
/*... w1 = (w1 - ni*w2) / mi*/
      omega1 = (omega1 - ni*omega2) / mi;
/*...................................................................*/

/*...*/
    for (i = 0; i < nEq; i++) {
/*... x(j + 2) = x(i) + w1 * r(i) + w2 * s(i) + alpha * u(i)*/
        x[i] += omega1*h[i] + omega2 *p[i] + alpha*z[i];
/*... r(j + 2) = r(j) - w1 * s(i) - w2 * t(i)*/
        r[i] -= (omega1 * s[i] + omega2 * t[i]);
    }
/*...................................................................*/

/*... (r, r)*/
    d = dot(r, r, nEq);
    if (sqrt(d) < conv) break;
/*...................................................................*/

/*...*/
    if (fHistLog)
      fprintf(fileHistLog, "%d %20.9e\n", j, sqrt(fabs(d)) / norm_b);
/*...................................................................*/

/*...*/
    for (i = 0; i < nEq; i++) {
/*... u(j + 1) = u(j) - w1 * v(i) - w2 * w(i)*/
      u[i] -= (omega1*v[i] + omega2*w[i]);
    }
/*...................................................................*/

/*...*/
    if (fHistLog)
      fprintf(fileHistLog, "%d %20.9e\n", j, sqrt(fabs(d)) / norm_b);
/*...................................................................*/

/*...*/
    if (jj == 2000) {
      jj = 0;
      printf("BICGSATB(2): %d %20.9e %20.9e\n"
             , j + 1, sqrt(fabs(d)), conv);
    }
    jj++;
/*...................................................................*/
  }
/*...................................................................*/

/*... Energy norm:  x*Kx*/
  matvec(nEq,ia,ja,a,ad,x,z);
  /*norma de energia = xT*A*x */
  xKx = dot(x, z, nEq);
/*...................................................................*/

/*... norm - 2 = || x ||*/
  norm = sqrt(dot(x, x, nEq));
/*...................................................................*/

/*... r = (b - Ax) (calculo do residuo explicito)*/
  for (i = 0; i < nEq; i++) {
    r[i] = b[i] - z[i];
  }
  norm_r = dot(r, r, nEq);
  norm_r = sqrt(norm_r);
  if (fPrint && norm_r > 3.16e0*conv)
    printf("BICGSATB(2): %20.9e > %20.9e!!\n", norm_r, conv);
/*...................................................................*/
  timef = getTimeC() - timei;

  if (fPrint) {
    printf(" (BICGSATB(2) solver:\n"
           "\tEquations      = %20d\n"
           "\tnad            = %20d\n"
           "\tSolver tol     = %20.2e\n"
           "\tIterarions     = %20d\n"
           "\tx * Kx         = %20.2e\n"
           "\t|| x ||        = %20.2e\n"
           "\t|| b - Ax ||   = %20.2e\n"
           "\tCPU time(s)    = %20.2lf\n"
           , nEq, nAd, tol, j + 1, xKx, norm, norm_r, timef);
  }

  if (j == maxIt) {
    printf(" PBICGSTAB(2) *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n", maxIt);
    exit(EXIT_FAILURE);
  }
  if (log)
    fprintf(fileLog
            , "PBICGSTAB(2): tol %20.2e "
            "iteration %d "
            "xKx %20.12e "
            "norma(x*x) %20.12e "
            "time %20.5lf\n"
            , tol, j + 1, xKx, norm, timef);


}
/**********************************************************************/

/**********************************************************************
* Data de criacao    : 01/09/2016                                    *
* Data de modificaco : 00/00/0000                                    *
* -------------------------------------------------------------------*
* PBICGSTABl2OMP : metodo do gradiente conjugado bi-ortaganalizado   *
* com precondiconador diagonal (matriz geral)                        *
* -------------------------------------------------------------------*
* Parametros de Entrada:                                             *
* -------------------------------------------------------------------*
*  neq       -> numero de equacoes                                   *
*  nad       -> numero de elementos nao nulos fora da diagonal       *
*  ia        -> estrutura de dados para matriz esparsa A             *
*  ja        -> estrutura de dados para matriz esparsa A             *
*  a         -> coef fora da diagonal principal                      *
*  ad        -> diagnal da matriz A                                  *
*   p        -> precondiconador diagonal                             *
*   b        -> vetor b (Ax=b)                                       *
*   x        -> vetor de solucao                                     *
*   t        -> vetor auxiliar                                       *
*   v        -> vetor auxiliar                                       *
*   r        -> vetor auxiliar                                       *
*   u        -> vetor auxiliar                                       *
*   r0       -> vetor auxiliar                                       *
*   w        -> vetor auxiliar                                       *
*   s        -> vetor auxiliar                                       *
*   p        -> vetor auxiliar                                       *
*   h        -> vetor auxiliar                                       *
*  z          > vetor auxiliar                                       *
* newX       -> vetor inicial iniciado com zero                      *
* fileLog       -> arquivo de log do solver                          *
* log        -> log de arquivo (true|false)                          *
* tol        -> tolerancia do solver                                 *
* maxIt			-> numero maximo de iteracoes                            *
*  newX			-> true zero o vetor inicial                             *
* fileLog				-> arquivo de saida do log                           *
* fileHistLog-> arquivo de log da iteracoes                          *
*  log       -> escreve o log do solver                              *
* fPrint			-> saida de informacao na tela                         *
* fHistLog		-> log das iteracoes                                   *
* -------------------------------------------------------------------*
* Parametros de Saida:                                               *
* -------------------------------------------------------------------*
* x[]-> atualizado                                                   *
* b[]-> inalterado                                                   *
* ad,a-> inalterado                                                  *
* -------------------------------------------------------------------*
* OBS:                                                               *
c * Versão do livro Iterative krylov Method for large linear Systems *
c * A(M-1)y=b precondicionador a direita                             *
* -------------------------------------------------------------------*
***********************************************************************/
void pbicgstabl2Omp(INT const nEq, INT const nAd
                 , INT *RESTRICT ia, INT *RESTRICT ja
                 , DOUBLE *RESTRICT a, DOUBLE *RESTRICT ad
                 , DOUBLE *RESTRICT m, DOUBLE *RESTRICT b
                 , DOUBLE *RESTRICT x, DOUBLE *RESTRICT t
                 , DOUBLE *RESTRICT v, DOUBLE *RESTRICT r
                 , DOUBLE *RESTRICT u, DOUBLE *RESTRICT r0
                 , DOUBLE *RESTRICT w, DOUBLE *RESTRICT s
                 , DOUBLE *RESTRICT p, DOUBLE *RESTRICT h
                 , DOUBLE *RESTRICT z
                 , DOUBLE const tol, unsigned int maxIt
                 , bool const newX, FILE* fileLog
                 , FILE *fileHistLog, bool const log
                 , bool const fHistLog, bool const fPrint
                 , BufferOmp *bOmp
                 , void(*matvec)(), DOUBLE(*dot)())
{
  short nThreads = ompVar.nThreadsSolver;
  unsigned int j, jj,js;
  INT i;
  DOUBLE alpha, beta, d, conv, xKx, norm_r, norm, norm_b, rr0, rr1
    , omega1, omega2, gamma, mi, ni, tau,tmp;
  DOUBLE timei, timef;

  timei = getTimeC();

#pragma omp parallel default(none) num_threads(nThreads)\
   private(i,j,jj,conv,norm_b,tmp,d,alpha,beta,rr0,rr1,omega1,omega2\
          ,gamma,mi,ni,tau)\
   shared(js,ia,ja,a,ad,m,b,x,t,v,r,u,r0,w,s,p,h,z,tol,maxIt,newX,fileLog\
         ,xKx,norm,norm_r\
         ,fileHistLog,log,fHistLog,fPrint,matvec,dot,bOmp,nThreads)
  {
/*... chute inicial*/
    if (newX)
#pragma omp for
      for (i = 0; i < nEq; i++)
        x[i] = 0.e0;
/*...................................................................*/

/*... conv = tol * |b|*/
    d      = dot(b,b,nEq);
    norm_b = sqrt(d);
    conv   = tol*norm_b;
/*...................................................................*/

/*... Ax0*/
    matvec(nEq
          ,ia            ,ja
          ,a             ,ad
          ,x             ,t
          ,bOmp->thBegin ,bOmp->thEnd
          ,bOmp->thHeight,bOmp->thY
          ,nThreads);
/*...................................................................*/

/*...*/
#pragma omp for
    for (i = 0; i < nEq; i++) {
    /*... r0 = b - Ax0*/
      r0[i] = b[i] - t[i];
    /*... r = r0*/
      r[i] = r0[i];
    /*... z = M(-1)p*/
      u[i] = 0.e0;
    }
/*...................................................................*/

/*...*/
    rr0 = 1.e0;
    alpha = 0.e0;
    omega2 = 1.e0;
/*...................................................................*/

/*...*/
    jj = 1;
    for (j = 0; j < maxIt; j++) {
/*... rr0 = -w2*rr0*/
      rr0 = -omega2*rr0;
/*...................................................................*/

/*... even BiCG step*/
/*... rr1 = (r,r0)*/
      rr1 = dot(r, r0, nEq);
/*... beta = alpha * rr1 / rr0*/
      beta = alpha*rr1 / rr0;
/*... rr0 = rr1*/
      rr0 = rr1;
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i<nEq; i++) {
/*... u(j) = r(j) - beta*u(i)*/
        u[i] = r[i] - beta*u[i];
/*... p = M(-1)u*/
        z[i] = u[i] * m[i];
      }
/*...................................................................*/

/*... v = Ap(j) = AM(-1)u*/
      matvec(nEq
            ,ia            ,ja
            ,a             ,ad
            ,z             ,v
            ,bOmp->thBegin ,bOmp->thEnd
            ,bOmp->thHeight,bOmp->thY
            ,nThreads);
/*...................................................................*/

/*... alpha = (v, r0) / (AM(-1)u(j), r0))*/
      gamma = dot(v, r0, nEq);
      alpha = rr0 / gamma;
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i<nEq; i++) {
/*... s(j) = r(j) - alpha*AM(-1)u(j)*/
        r[i] -= alpha * v[i];
/*... p = M(-1)r*/
        p[i] = r[i] * m[i];
      }
/*...................................................................*/

/*... s = Ap(j) = AM(-1)r*/
      matvec(nEq
          ,ia            ,ja
          ,a             ,ad
          ,p             ,s
          ,bOmp->thBegin ,bOmp->thEnd
          ,bOmp->thHeight,bOmp->thY
          ,nThreads);
/*...................................................................*/
     
/*...*/
#pragma omp for
      for (i = 0; i<nEq; i++) {
/*... x(j) = x(j) + alpha*(M-1)u*/
        x[i] += alpha*z[i];
      }
/*...................................................................*/

/*... odd BiCG step*/
/*... rr1 = (s, r0)*/
      rr1 = dot(s, r0, nEq);
/*... beta = alpha*rr1 / rr0*/
      beta = alpha*rr1 / rr0;
/*... rr0 = rr1*/
      rr0 = rr1;
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i<nEq; i++) {
/*... v(j) = s(j) - beta*v*/
        v[i] = s[i] - beta*v[i];
/*... p = M(-1)v*/
        p[i] = v[i] * m[i];
      }
/*...................................................................*/

/*... s = Ap(j) = AM(-1)v*/
      matvec(nEq
            ,ia            ,ja
            ,a             ,ad 
            ,p             ,w
            ,bOmp->thBegin ,bOmp->thEnd
            ,bOmp->thHeight,bOmp->thY
            ,nThreads);
/*...................................................................*/

/*... gamma = (w,r0)*/
      gamma = dot(w, r0, nEq);
      alpha = rr0 / gamma;
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i<nEq; i++) {
/*... u(j) = r(j) - beta*u(i)*/
        u[i] = r[i] - beta*u[i];
/*... r(j) = r(j) - alpha*v(i)*/
        r[i] -= alpha*v[i];
/*... s(j) = s(j) - alpha*w(i)*/
        s[i] -= alpha * w[i];
/*... p = M(-1)s*/
        p[i] = s[i] * m[i];
/*... z = M(-1)u*/
        z[i] = u[i] * m[i];
/*... h = M(-1)r*/
        h[i] = r[i] * m[i];
      }
/*...................................................................*/

/*... t = Ap = AM(-1)s(j)*/
      matvec(nEq
            ,ia            ,ja
            ,a             ,ad
            ,p             ,t
            ,bOmp->thBegin ,bOmp->thEnd
            ,bOmp->thHeight,bOmp->thY
            ,nThreads);
/*...................................................................*/

/*... GCR(2) - part*/
/*... w1 = (r, s)*/
      omega1 = dot(r, s, nEq);
/*... mi = (s, s)*/
      mi = dot(s, s, nEq);
/*... ni = (s, t)*/
      ni = dot(s, t, nEq);
/*... tau = (t, t)*/
      tau = dot(t, t, nEq);
/*... w2*/
      omega2 = dot(r, t, nEq);
/*... tau = tau - ni*ni / mi*/
      tau = tau - ni*ni / mi;
/*... w2 = (w2 - ni*w1 / mi) / tau*/
      omega2 = (omega2 - ni*omega1 / mi) / tau;
/*... w1 = (w1 - ni*w2) / mi*/
      omega1 = (omega1 - ni*omega2) / mi;
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i < nEq; i++) {
/*... x(j + 2) = x(i) + w1 * r(i) + w2 * s(i) + alpha * u(i)*/
        x[i] += omega1*h[i] + omega2 *p[i] + alpha*z[i];
/*... r(j + 2) = r(j) - w1 * s(i) - w2 * t(i)*/
        r[i] -= (omega1 * s[i] + omega2 * t[i]);
      }
/*...................................................................*/

/*... (r, r)*/
      d = dot(r, r, nEq);
      if (sqrt(d) < conv) break;
/*...................................................................*/

/*...*/
#pragma omp master
      if (fHistLog)
        fprintf(fileHistLog, "%d %20.9e\n", j, sqrt(fabs(d)) / norm_b);
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i < nEq; i++)
      /*... u(j + 1) = u(j) - w1 * v(i) - w2 * w(i)*/
        u[i] -= (omega1*v[i] + omega2*w[i]);
/*...................................................................*/

/*...*/
      if (jj == 2000) {
        jj = 0;
#pragma omp master
        printf("BICGSATBOMP(2): %d %20.9e %20.9e\n"
               , j + 1, sqrt(d), conv);
      }
      jj++;
/*...................................................................*/
    }
/*...................................................................*/

#pragma omp single
    js = j;

/*... Energy norm:  x*Kx*/
    matvec(nEq
          ,ia            ,ja
          ,a             ,ad
          ,x             ,z
          ,bOmp->thBegin ,bOmp->thEnd
          ,bOmp->thHeight,bOmp->thY
          ,nThreads);
/*... norma de energia = xT*A*x*/
    tmp = dot(x,z,nEq);
#pragma omp master
    xKx = tmp;
/*...................................................................*/

/*... norm - 2 = || x ||*/
    tmp = sqrt(dot(x, x, nEq));
#pragma omp master
    norm = tmp;
/*...................................................................*/

/*... r = (b - Ax) (calculo do residuo explicito)*/
#pragma omp for
    for (i = 0; i < nEq; i++) {
      r[i] = b[i] - z[i];
    }
    tmp = dot(r, r, nEq);
#pragma omp single
    {  
      norm_r = sqrt(tmp);
      if (fPrint && norm_r > 3.16e0*conv)
      printf("BICGSATBOMP(2): %20.9e > %20.9e!!\n", norm_r, conv);
    }
/*...................................................................*/
  }
/*...................................................................*/

  timef = getTimeC() - timei;

  if (fPrint) {
    printf(" (BICGSATBOMP(2) solver:\n"
           "\tEquations      = %20d\n"
           "\tnad            = %20d\n"
           "\tSolver tol     = %20.2e\n"
           "\tIterarions     = %20d\n"
           "\tx * Kx         = %20.2e\n"
           "\t|| x ||        = %20.2e\n"
           "\t|| b - Ax ||   = %20.2e\n"
           "\tCPU time(s)    = %20.2lf\n"
           , nEq, nAd, tol, j + 1, xKx, norm, norm_r, timef);
  }

  if (js == maxIt) {
    printf(" PBICGSTABOMP(2) *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n", maxIt);
    exit(EXIT_FAILURE);
  }
  if (log)
    fprintf(fileLog
            , "PBICGSTABOMP(2): tol %20.2e "
            "iteration %d "
            "xKx %20.12e "
            "norma(x*x) %20.12e "
            "time %20.5lf\n"
            , tol, j + 1, xKx, norm, timef);


}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 04/09/2016                                    *
 * Data de modificaco : 00/00/0000                                    *
 * -------------------------------------------------------------------*
 * GMRES: Solucao iterativa de sistemas simetricos e nao-simetricos   *
 *        pelo metodo GMRES com precondicionador diagonal.            *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 *  neq       -> numero de equacoes                                   *
 *  nad       -> numero de elementos nao nulos fora da diagonal       *
 *  ia        -> estrutura de dados para matriz esparsa A             *
 *  ja        -> estrutura de dados para matriz esparsa A             *
 *  a         -> coef fora da diagonal principal                      *
 *  ad        -> diagnal da matriz A                                  *
 *   p        -> precondiconador diagonal                             *
 *   b        -> vetor b (Ax=b)                                       *
 *   x        -> vetor de solucao                                     *
 *   g        -> vetor auxiliar                                       *
 *   h        -> vetor auxiliar                                       *
 *   y        -> vetor auxiliar                                       *
 *   c        -> vetor auxiliar                                       *
 *   s        -> vetor auxiliar                                       *
 *   e        -> vetor auxiliar                                       *
 * newX       -> vetor inicial iniciado com zero                      *
 * fileLog       -> arquivo de log do solver                          *
 * log        -> log de arquivo (true|false)                          *
 * tol        -> tolerancia do solver                                 *
 * maxIt			-> numero maximo de iteracoes                           *
 *  newX			-> true zero o vetor inicial                            *
 * fileLog				-> arquivo de saida do log                          *
 * fileHistLog-> arquivo de log da iteracoes                          *
 *  log       -> escreve o log do solver                              *
 * fPrint			-> saida de informacao na tela                          *
 * fHistLog		-> log das iteracoes                                    *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *
 * b[]-> inalterado                                                   *
 * ad,a-> inalterado                                                  *
 * -------------------------------------------------------------------*
 * OBS:                                                               *
 * Versão do livro Iterative krylov Method for large linear Systems   *
 *                                                                    *
 * g(nKrylov+1,neq)                                                   *
 * h(nKrylov+1,nKrylov)                                               *
 * y(nKrylov)                                                         *
 * c(nKrylov)                                                         *
 * s(nKrylov)                                                         *
 * e(nKrylov+1)                                                       *
 **********************************************************************/
void gmres(INT const nEq       ,INT const nAd
          ,INT *RESTRICT ia    ,INT *RESTRICT ja
          ,DOUBLE *RESTRICT a  ,DOUBLE *RESTRICT ad
          ,DOUBLE *RESTRICT m  ,DOUBLE *RESTRICT b
          ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT g
          ,DOUBLE *RESTRICT h  ,DOUBLE *RESTRICT y
          ,DOUBLE *RESTRICT c  ,DOUBLE *RESTRICT s 
          ,DOUBLE *RESTRICT e  ,short const nKrylov
          ,DOUBLE const tol    ,unsigned int nCycles
          ,bool const newX     ,FILE* fileLog
          ,FILE *fileHistLog   ,bool const log
          ,bool const fHistLog ,bool const fPrint
          ,void(*matvec)()     ,DOUBLE(*dot)())
{
  int i,j,jj,ni,nIt;
  unsigned short l,nCol = nKrylov;
  DOUBLE *g1,*g2;
  DOUBLE tmp,norm,norm_r,eConv,beta,h1,h2,aux1,aux2,r,xKx;
  INT iLong;
  DOUBLE timei, timef;

  timei = getTimeC();

/*... chute inicial*/
  if (newX)
    for (i = 0; i < nEq; i++)
      x[i] = 0.e0;
/*...................................................................*/

/*...g(1,i) = (M-1)*b*/
  for (iLong = 0; iLong < nEq; iLong++)
    g[iLong] = b[iLong] * m[iLong];
/*...................................................................*/

/*... Limite de convergencia*/
  norm = dot(g, g, nEq);
  norm = sqrt(norm);
  eConv = tol*norm;
/*...................................................................*/

/*... Ciclos Gmres*/
  nIt = 0;
  jj = 0;
  for (l = 0; l < nCycles; l++) {
/*... Residuo g1 = b - Ax*/
    matvec(nEq, ia, ja, a, ad, x, g);
/*...................................................................*/

/*...g1 = (M-1)*r*/
    for (iLong = 0; iLong < nEq; iLong++)
      g[iLong] = (b[iLong] - g[iLong]) * m[iLong];
/*...................................................................*/

/*... Norma do residuo*/
    e[0] = sqrt(dot(g, g, nEq));
    tmp = 1.e0 / e[0];
/*...................................................................*/

/*...g(1,i) = g(1,i)/|g(1,i)|*/
    for (iLong = 0; iLong < nEq; iLong++)
      g[iLong] *= tmp;
/*...................................................................*/

/*... Iteracoe Gmres*/
     for (ni = 0; ni < nKrylov; ni++, nIt++) {
/*... g1 = g(ni)*/
      g1 = &g[ni*nEq];
/*... g2 = g(ni+1)*/
      g2 = &g[(ni+1)*nEq];
/*... Residuo  g(ni+1) = A*g(ni)*/
      matvec(nEq, ia, ja, a, ad, g1, g2);
/*...................................................................*/

/*... g2 = (M-1)g2*/
      for (iLong = 0; iLong < nEq; iLong++)
        g2[iLong] *= m[iLong];
/*...................................................................*/

/*... Ortogonalizacao (Gram-Schmidt modificado)*/
      for (j = 0; j<=ni; j++) {
/*... g1 = g(j)*/
//      g1 = &MAT2D(j, 0, g, nEq);
        g1   = &g[j*nEq];
        beta = dot(g2, g1, nEq);
        for (iLong = 0; iLong < nEq; iLong++)
          g2[iLong] -= beta*g1[iLong];
        MAT2D(j, ni, h, nCol) = beta;
      }
/*...................................................................*/

/*... g(i+1) = |g(i+1)|*/
      norm = sqrt(dot(g2, g2, nEq));
      MAT2D(ni+1,ni,h,nCol) = norm;
/*...................................................................*/

/*...g(1,i+1) = g(1,i+1)/|g(1,i+1)|*/
      tmp = 1.e0 / norm;
      for (iLong = 0; iLong < nEq; iLong++)
        g2[iLong] *= tmp;
/*...................................................................*/
     
/*...*/
      for(j=0; j<ni ; j++) {
        h1                   = MAT2D(j  ,ni,h,nCol);
        h2                   = MAT2D(j+1,ni,h,nCol);
        aux1                 =  c[j]*h1 + s[j]*h2;
        aux2                 = -s[j]*h1 + c[j]*h2;
        MAT2D(j  ,ni,h,nCol) = aux1;
        MAT2D(j+1,ni,h,nCol) = aux2;
      }
      h1                    = MAT2D(ni  ,ni,h,nCol);
      h2                    = MAT2D(ni+1,ni,h,nCol);
      r                     = sqrt(h1*h1 + h2*h2);
      tmp                   = 1.e0 / r;
      c[ni]                 = h1*tmp;
      s[ni]                 = h2*tmp;
      MAT2D(ni  ,ni,h,nCol) = r;
      MAT2D(ni+1,ni,h,nCol) = 0.e0;
      e[ni+1]               = -s[ni]*e[ni];
      e[ni]                 =  c[ni]*e[ni];
      if (fabs(e[ni+1]) <= eConv);
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    if( ni == nKrylov) ni--;
/*.....................................................................*/

/*... h y = e*/
     y[ni] = e[ni] / MAT2D(ni, ni, h, nCol);
     for (i = ni - 1; i>=0; i--) {
       y[i] = 0.e0;
       for (j = i + 1; j<ni; j++)
         y[i] -= MAT2D(i, j, h, nCol)*y[j];
       y[i] = (y[i] + e[i]) / MAT2D(i, i, h, nCol);
     }
/*...................................................................*/

/*... x = Vy = (Vt)y*/
     for (j = 0; j<ni; j++){
        tmp = y[j];
        for (iLong = 0; iLong<nEq; iLong++)    
          x[iLong] += MAT2D(j, iLong, g, nEq)*tmp;
     }
/*...................................................................*/

/*... Verifica a convergencia:*/
    if (fabs(e[ni+1]) <= eConv) break;
/*...................................................................*/
  }

/*... Energy norm:  x*Kx*/
  matvec(nEq, ia, ja, a, ad, x, g);
/*norma de energia = xT*A*x */
  xKx = dot(x, g, nEq);
/*...................................................................*/

/*... norm - 2 = || x ||*/
  norm = sqrt(dot(x, x, nEq));
/*...................................................................*/

/*... r = (b - Ax) (calculo do residuo explicito)*/
  g1 = g;
  g2 = &g[nEq];
  for (iLong = 0; iLong < nEq; iLong++)
    g2[iLong] = b[iLong] - g1[iLong];
  norm_r = dot(g2,g2,nEq);
  norm_r = sqrt(norm_r);
  if (fPrint && norm_r > 3.16e0*eConv)
    printf("GMRES: %20.9e > %20.9e!!\n", norm_r, eConv);
/*...................................................................*/
  timef = getTimeC() - timei;

  if (fPrint) {
    printf(" (GMRES) solver:\n"
           "\tEquations      = %20d\n"
           "\tnad            = %20d\n"
           "\tnKrylov        = %20d\n"
           "\tSolver tol     = %20.2e\n"
           "\tCycles         = %20d\n"
           "\tIterarions     = %20d\n"
           "\tx * Kx         = %20.2e\n"
           "\t|| x ||        = %20.2e\n"
           "\t|| b - Ax ||   = %20.2e\n"
           "\tCPU time(s)    = %20.2lf\n"
           , nEq, nAd, nKrylov,tol,l+1,nIt+1,xKx, norm, norm_r, timef);
  }

  if (l == nCycles) {
    printf(" GMRES *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n", nCycles);
    exit(EXIT_FAILURE);
  }

  if (log)
    fprintf(fileLog
            , "GMRES: tol %20.2e "
            " Cycles %d"
            " iteration %d "
            " nKrylov %d "
            " xKx %20.12e "
            " norma(x*x) %20.12e "
            " time %20.5lf\n"
            , tol,l+1,nIt+1,nKrylov,xKx,norm,timef);

}
/**********************************************************************/

/**********************************************************************
* Data de criacao    : 04/09/2016                                    *
* Data de modificaco : 00/00/0000                                    *
* -------------------------------------------------------------------*
* GMRESOMP: Solucao iterativa de sistemas simetricos e nao-simetricos*
*        pelo metodo GMRES com precondicionador diagonal.            *
* -------------------------------------------------------------------*
* Parametros de Entrada:                                             *
* -------------------------------------------------------------------*
*  neq       -> numero de equacoes                                   *
*  nad       -> numero de elementos nao nulos fora da diagonal       *
*  ia        -> estrutura de dados para matriz esparsa A             *
*  ja        -> estrutura de dados para matriz esparsa A             *
*  a         -> coef fora da diagonal principal                      *
*  ad        -> diagnal da matriz A                                  *
*   p        -> precondiconador diagonal                             *
*   b        -> vetor b (Ax=b)                                       *
*   x        -> vetor de solucao                                     *
*   g        -> vetor auxiliar                                       *
*   h        -> vetor auxiliar                                       *
*   y        -> vetor auxiliar                                       *
*   c        -> vetor auxiliar                                       *
*   s        -> vetor auxiliar                                       *
*   e        -> vetor auxiliar                                       *
* newX       -> vetor inicial iniciado com zero                      *
* fileLog       -> arquivo de log do solver                          *
* log        -> log de arquivo (true|false)                          *
* tol        -> tolerancia do solver                                 *
* maxIt			-> numero maximo de iteracoes                            *
*  newX			-> true zero o vetor inicial                             *
* fileLog				-> arquivo de saida do log                           *
* fileHistLog-> arquivo de log da iteracoes                          *
*  log       -> escreve o log do solver                              *
* fPrint			-> saida de informacao na tela                         *
* fHistLog		-> log das iteracoes                                   *
* -------------------------------------------------------------------*
* Parametros de Saida:                                               *
* -------------------------------------------------------------------*
* x[]-> atualizado                                                   *
* b[]-> inalterado                                                   *
* ad,a-> inalterado                                                  *
* -------------------------------------------------------------------*
* OBS:                                                               *
* Versão do livro Iterative krylov Method for large linear Systems   *
*                                                                    *
* g(nKrylov+1,neq)                                                   *
* h(nKrylov+1,nKrylov)                                               *
* y(nKrylov)                                                         *
* c(nKrylov)                                                         *
* s(nKrylov)                                                         *
* e(nKrylov+1)                                                       *
**********************************************************************/
void gmresOmp(INT const nEq      ,INT const nAd
             ,INT *RESTRICT ia   ,INT *RESTRICT ja
             ,DOUBLE *RESTRICT a ,DOUBLE *RESTRICT ad
             ,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b
             ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT g
             ,DOUBLE *RESTRICT h ,DOUBLE *RESTRICT y
             ,DOUBLE *RESTRICT c ,DOUBLE *RESTRICT s
             ,DOUBLE *RESTRICT e ,short const nKrylov
             ,DOUBLE const tol   ,unsigned int nCycles
             ,bool const newX    ,FILE* fileLog
             ,FILE *fileHistLog  ,bool const log
             ,bool const fHistLog,bool const fPrint
             ,BufferOmp *bOmp
             ,void(*matvec)()    ,DOUBLE(*dot)())
{
  short nThreads = ompVar.nThreadsSolver;
  int i, j, jj, ni, nIt;
  unsigned short l,nCol = nKrylov;
  DOUBLE *g1, *g2;
  DOUBLE tmp, norm, norm_r, eConv, beta, h1, h2, aux1, aux2, r, xKx;
  INT iLong;
  DOUBLE timei, timef;

  timei = getTimeC();

  omp_set_num_threads(nThreads);

/*... chute inicial*/
  if (newX)
    for (i = 0; i < nEq; i++)
      x[i] = 0.e0;
/*...................................................................*/

/*...g(1,i) = (M-1)*b*/
#pragma omp parallel default(none) private(tmp,iLong)\
        shared(nEq,g,b,m,norm,dot)
  {
  #pragma omp for
    for (iLong = 0; iLong < nEq; iLong++)
      g[iLong] = b[iLong] * m[iLong];
/*...................................................................*/

/*... Limite de convergencia*/
    tmp = dot(g, g, nEq);
  #pragma omp single
    norm = tmp;
  }
/*...................................................................*/

  norm = sqrt(norm);
  eConv = tol*norm;
/*...................................................................*/

/*... Ciclos Gmres*/
  nIt = 0;
  jj = 0;
  for (l = 0; l < nCycles; l++) {
/*... Residuo g1 = b - Ax*/
#pragma omp parallel default(none) private(tmp)\
        shared(nEq,ia,ja,a,ad,x,g,b,m,e,bOmp,nThreads,matvec,dot)
  {
    matvec(nEq
          ,ia            ,ja
          ,a             ,ad
          ,x             ,g
          ,bOmp->thBegin ,bOmp->thEnd
          ,bOmp->thHeight,bOmp->thY
          ,nThreads);
/*...................................................................*/

/*...g1 = (M-1)*r*/
  #pragma omp for
    for (iLong = 0; iLong < nEq; iLong++)
      g[iLong] = (b[iLong] - g[iLong]) * m[iLong];
/*...................................................................*/

/*... Norma do residuo*/
    tmp  = sqrt(dot(g, g, nEq));
  #pragma omp single
    e[0] = tmp;
    tmp = 1.e0/tmp;
/*...................................................................*/

/*...g(1,i) = g(1,i)/|g(1,i)|*/
  #pragma omp for
    for (iLong = 0; iLong < nEq; iLong++)
      g[iLong] *= tmp;
/*...................................................................*/
  }
/*...................................................................*/

/*... Iteracoe Gmres*/
    for (ni = 0; ni < nKrylov; ni++, nIt++) {
/*... g1 = g(ni)*/
      g1 = &g[ni*nEq];
/*... g2 = g(ni+1)*/
      g2 = &g[(ni + 1)*nEq];
#pragma omp parallel default(none) private(tmp)\
        shared(nEq,ia,ja,a,ad,g1,g2,m,bOmp,nThreads,matvec,dot)
  {
/*... Residuo  g(ni+1) = A*g(ni)*/
      matvec(nEq
            ,ia            ,ja
            ,a             ,ad
            ,g1            ,g2
            ,bOmp->thBegin ,bOmp->thEnd
            ,bOmp->thHeight,bOmp->thY
            ,nThreads);
/*...................................................................*/

/*... g2 = (M-1)g2*/
    #pragma omp for
      for (iLong = 0; iLong < nEq; iLong++)
        g2[iLong] *= m[iLong];
/*...................................................................*/
  }
/*...................................................................*/

/*... Ortogonalizacao (Gram-Schmidt modificado)*/
      for (j = 0; j <= ni; j++) {
/*... g1 = g(j)*/
        g1   = &g[j*nEq];
#pragma omp  parallel default(none) private(tmp,iLong)\
       shared(nEq,g1,g2,beta,dot) 
  {
        tmp = dot(g2, g1, nEq);
    #pragma omp single
        beta = tmp;  
/*...................................................................*/

/*...*/
    #pragma omp for
        for (iLong = 0; iLong < nEq; iLong++)
          g2[iLong] -= beta*g1[iLong];
  }
/*...................................................................*/
        MAT2D(j, ni, h, nCol) = beta;
      }
/*...................................................................*/

#pragma omp parallel default(none) private(tmp,iLong)\
      shared(nEq,norm,h,g1,g2,beta,dot,ni,nCol)
  { 
/*... g(i+1) = |g(i+1)|*/
      tmp  = sqrt(dot(g2, g2, nEq));
    #pragma omp single
    {
      norm = tmp;
      MAT2D(ni + 1, ni, h, nCol) = norm;
    }
/*...................................................................*/

/*...g(1,i+1) = g(1,i+1)/|g(1,i+1)|*/
      tmp = 1.e0 / tmp;
#pragma omp for
      for (iLong = 0; iLong < nEq; iLong++)
        g2[iLong] *= tmp;
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
      for (j = 0; j<ni; j++) {
        h1 = MAT2D(j, ni, h, nCol);
        h2 = MAT2D(j + 1, ni, h, nCol);
        aux1 = c[j] * h1 + s[j] * h2;
        aux2 = -s[j] * h1 + c[j] * h2;
        MAT2D(j, ni, h, nCol) = aux1;
        MAT2D(j + 1, ni, h, nCol) = aux2;
      }
      h1 = MAT2D(ni, ni, h, nCol);
      h2 = MAT2D(ni + 1, ni, h, nCol);
      r = sqrt(h1*h1 + h2*h2);
      tmp = 1.e0 / r;
      c[ni] = h1*tmp;
      s[ni] = h2*tmp;
      MAT2D(ni, ni, h, nCol) = r;
      MAT2D(ni + 1, ni, h, nCol) = 0.e0;
      e[ni + 1] = -s[ni] * e[ni];
      e[ni] = c[ni] * e[ni];
      if (fabs(e[ni + 1]) <= eConv);
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    if (ni == nKrylov) ni--;
/*.....................................................................*/

/*... h y = e*/
    y[ni] = e[ni] / MAT2D(ni, ni, h, nCol);
    for (i = ni - 1; i >= 0; i--) {
      tmp = 0.e0;
#pragma omp parallel for reduction(-:tmp)
      for (j = i + 1; j<ni; j++)
        tmp -= MAT2D(i, j, h, nCol)*y[j];
      y[i] = (tmp + e[i]) / MAT2D(i, i, h, nCol);
    }
/*...................................................................*/

/*... x = Vy = (Vt)y*/
    for (j = 0; j<ni; j++) {
      tmp = y[j];
#pragma omp parallel for default (none) private(iLong) shared(tmp,x,g,j,nEq)
      for (iLong = 0; iLong<nEq; iLong++)
        x[iLong] += MAT2D(j, iLong, g, nEq)*tmp;
    }
/*...................................................................*/

/*... Verifica a convergencia:*/
    if (fabs(e[ni+1]) <= eConv) break;
/*...................................................................*/
  }

/*... Energy norm:  x*Kx*/
#pragma omp parallel default(none) private(tmp)\
        shared(nEq,ia,ja,a,ad,x,g,bOmp,nThreads,xKx,norm,matvec,dot)
  {
    matvec(nEq
          ,ia            ,ja
          ,a             ,ad
          ,x             ,g
          ,bOmp->thBegin ,bOmp->thEnd
          ,bOmp->thHeight,bOmp->thY
          ,nThreads);
/*norma de energia = xT*A*x */
    tmp = dot(x, g, nEq);
  #pragma omp single
    xKx = tmp;
/*...................................................................*/

/*... norm - 2 = || x ||*/
    tmp = sqrt(dot(x, x, nEq));
#pragma omp single
    norm = tmp;
/*...................................................................*/
  }
/*...................................................................*/

/*... r = (b - Ax) (calculo do residuo explicito)*/
  g1 = g;
  g2 = &g[nEq];
  for (iLong = 0; iLong < nEq; iLong++)
    g2[iLong] = b[iLong] - g1[iLong];
  norm_r = dot(g2, g2, nEq);
  norm_r = sqrt(norm_r);
  if (fPrint && norm_r > 3.16e0*eConv)
    printf("GMRES: %20.9e > %20.9e!!\n", norm_r, eConv);
/*...................................................................*/
  timef = getTimeC() - timei;

  if (fPrint) {
    printf(" (GMRESOMP) solver:\n"
           "\tEquations      = %20d\n"
           "\tnad            = %20d\n"
           "\tnKrylov        = %20d\n"
           "\tSolver tol     = %20.2e\n"
           "\tCycles         = %20d\n"
           "\tIterarions     = %20d\n"
           "\tx * Kx         = %20.2e\n"
           "\t|| x ||        = %20.2e\n"
           "\t|| b - Ax ||   = %20.2e\n"
           "\tCPU time(s)    = %20.2lf\n"
           , nEq, nAd, nKrylov, tol, l + 1, nIt + 1, xKx, norm, norm_r, timef);
  }

  if (l == nCycles) {
    printf(" GMRESOMP *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n", nCycles);
    exit(EXIT_FAILURE);
  }

  if (log)
    fprintf(fileLog
            , "GMRESOMP: tol %20.2e "
            " Cycles %d"
            " iteration %d "
            " nKrylov %d "
            " xKx %20.12e "
            " norma(x*x) %20.12e "
            " time %20.5lf\n"
            , tol, l + 1, nIt + 1, nKrylov, xKx, norm, timef);

}
/**********************************************************************/