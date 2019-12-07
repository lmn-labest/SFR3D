#include<Solv.h>

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
void pcg(INT const nEq      ,INT const nAd
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

	d      = dot(b, z, nEq);
	norm_b = sqrt(fabs(d));
	conv   = tol * norm_b;
/*...................................................................*/

/*... Ax0*/
	matvec(nEq, ia, ja, al, ad, x, z);
/*...................................................................*/

/*...*/
	for (i = 0; i < nEq; i++)
  {
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
	for (j = 0; j < maxIt; j++)
  {
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
		if (jj == 5000) {
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

  if(fPrint)
  { 
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
 * Data de modificaco : 23 / 08 / 2019                                *
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
        ,DOUBLE *RESTRICT r ,DOUBLE *RESTRICT p
        ,DOUBLE const tol   ,unsigned int maxIt
        ,bool const newX    ,FILE* fileLog  
        ,FILE *fileHistLog  ,bool const log     
        ,bool const fHistLog,bool const fPrint
        ,Interface *iNeq                      
        ,void(*matvec)()    ,DOUBLE(*dot)())
{
#ifdef _MPI_
  unsigned int j, jj;
	INT i;
	DOUBLE alpha, beta, d, conv, xKx, norm_b, norm, norm_r_m, norm_r, di;
  DOUBLE timei,timef;
  INT param[2];

  timei = getTimeC();

  param[0] = nAd;
  param[1] = nAdR;

/* chute inicial*/
  if(newX)  
    for(i = 0; i < nEq; i++)  
      x[i] = 0.e0;
/*...................................................................*/
  
/*... conv = tol * |(M-1)b|m = tol(b,M(-1)b) */
	for (i = 0; i < nEqNov; i++)
		z[i] = b[i] * m[i];

	d = dot(b, z, nEqNov);
	norm_b = sqrt(fabs(d));
	conv = tol * norm_b;
/*...................................................................*/

/*... Ax0*/
  matvec(nEqNov,param,ia,ja,al,ad,x,z,iNeq);
/*...................................................................*/ 

/*...*/
  for(i = 0; i < nEqNov; i++)
  {
/* ... r0 = b - Ax0*/
		r[i] = b[i] - z[i];
/* ... z0 = (M-1)r0*/
		z[i] = r[i] * m[i];
/* ... p0 = r0*/
		p[i] = z[i];
	}
/*...................................................................*/
  
/*... (r(0), z(0)) = (r(0), (M-1)r0)*/
	d = dot(r, z, nEqNov);
/*...................................................................*/

/*...*/
	jj = 1;
  for(j = 0; j < maxIt; j++)
  {
/*... z = Ap(j)*/
    matvec(nEqNov,param,ia,ja,al,ad,p,z,iNeq);
/*...................................................................*/

/*... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))*/    
    alpha = d / dot(z,p,nEqNov);
/*...................................................................*/

/*...*/
    for(i = 0; i < nEqNov; i++)
    {
/*... x(j + 1) = x(j) + alpha*p*/
			x[i] += alpha * p[i];
/*... r(j + 1) = r(j) - alpha*Ap*/
			r[i] -= alpha * z[i];
/*... z = (M-1)r0*/
			z[i] = r[i] * m[i];
		}
/*...................................................................*/

/* ... (r(j + 1), (M - 1)r(j + 1)) = (r(j + 1), z)*/
		di = dot(r, z, nEqNov);
/*... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) */
		beta = di / d;
/*...................................................................*/

/* ... p(j + 1) = (M - 1)r(j + 1) + beta*p(j) = z + beta*p(j)*/
 		for (i = 0; i < nEqNov; i++)
		  p[i] = z[i] + beta * p[i];
/*...................................................................*/

/*...*/
		if (fHistLog && !mpiVar.myId)
		  fprintf(fileHistLog, "%d %20.9e\n", j, sqrt(fabs(d)) / norm_b);
/*...................................................................*/

/*...*/
		d = di;
		if (sqrt(fabs(d)) < conv) break;
/*...................................................................*/

/*...*/
		if (jj == 5000)
    {
  		jj = 0;
			printf("MPIPCG: %d %20.9e %20.9e\n", j+1, sqrt(fabs(d)), conv);
		}
		jj++;
/*...................................................................*/
	}
/*...................................................................*/

/*... Energy norm:  x*Kx*/
  matvec(nEqNov,param,ia,ja,al,ad,x,z,iNeq);
/*norma de energia = xT*A*x */
	xKx    = dot(x, z, nEqNov);
/*...................................................................*/

/*... norm - 2 = || x ||*/
	norm = sqrt(dot(x, x, nEqNov));
/*...................................................................*/

/*... r = M(-1)(b - Ax) (calculo do residuo explicito)*/
	for (i = 0; i < nEqNov; i++) {
		r[i] = b[i] - z[i];
		z[i] = r[i] * m[i];
	}
	norm_r_m = dot(r, z, nEqNov);
	norm_r_m = sqrt(fabs(norm_r_m));
	norm_r = dot(r, r, nEqNov);
	norm_r = sqrt(norm_r);
	if(fPrint && norm_r_m > 3.16e0*conv && !mpiVar.myId)
	  printf("MPIPCG: %20.9e > %20.9e!!\n", norm_r_m, conv);
/*...................................................................*/

  timef = getTimeC() - timei;

/*...*/
  if(!mpiVar.myId && fPrint)
  { 
    printf(" (MPIPCG) solver:\n"
			"\tEquations      = %20d\n"
      "\tEquations      = %20d\n"
			"\tnad            = %20d\n"
      "\tnadR           = %20d\n"
			"\tSolver tol     = %20.2e\n"
			"\tIterarions     = %20d\n"
			"\tx * Kx         = %20.2e\n"
			"\t|| x ||        = %20.2e\n"
			"\t|| b - Ax ||   = %20.2e\n"
			"\tCPU time(s)    = %20.2lf\n"
	    ,nEq,nEqNov,nAd,nAdR,tol,j+1,xKx, norm, norm_r,timef);
  }
/*..................................................................*/

/*...*/  
  if(j == maxIt)
  { 
    if(!mpiVar.myId)
    {
      printf(" MPIPCG *** WARNING: no convergende reached !\n");
      printf("MAXIT = %d \n",maxIt);
    }
    mpiStop();
    exit(EXIT_SOLVER);
  }
/*..................................................................*/ 

/*...*/
  if(!mpiVar.myId && log)
    fprintf(fileLog          
           ,"MPIPCG: tol %20.2e " 
            "iteration %d " 
						"xKx %20.12e "
		        "norma(x*x) %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1, xKx,norm,timef);
 /*..................................................................*/
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
  unsigned int j, jj, jG;
  INT i;
  DOUBLE alpha, beta, d, conv, xKx, norm_b, norm, norm_r_m, norm_r;
  DOUBLE di,tmp;
  DOUBLE timei, timef;
  timei = getTimeC();

  norm = 0.e0;

#pragma omp parallel default(none) num_threads(nThreads)\
  private(i,j,jj,conv,norm_b,tmp,norm,d,di,alpha,beta)\
  shared(ia,ja,a,ad,m,b,x,z,r,p,maxIt,fileLog,xKx,norm_r_m,norm_r\
         ,fileHistLog,matvec,dot,bOmp,nThreads,jG)
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
    for (i = 0; i < nEq; i++)
    {
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
    for (j = 0; j < maxIt; j++)
    {
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
      for (i = 0; i < nEq; i++) 
      {
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
      if (jj == 5000) 
      {
        jj = 0;
#pragma omp master
        printf("PCGOMP: %d %20.9e %20.9e\n", j + 1, sqrt(fabs(d)), conv);
      }
      jj++;
/*...................................................................*/
    }
/*...................................................................*/

#pragma omp master
    jG = j;

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
#pragma omp master
    xKx = tmp;
/*...................................................................*/

/*... norm - 2 = || x ||*/
    tmp  = sqrt(dot(x, x, nEq));
#pragma omp master
    norm = tmp;
/*...................................................................*/

/*... r = M(-1)(b - Ax) (calculo do residuo explicito)*/
#pragma omp for
    for (i = 0; i < nEq; i++) 
    {
      r[i] = b[i] - z[i];
      z[i] = r[i] * m[i];
    }

    tmp = dot(r, z, nEq);
#pragma omp master
    norm_r_m = sqrt(fabs(tmp));
    
    tmp = dot(r, r, nEq);
#pragma omp master
    norm_r = sqrt(tmp);

#pragma omp master
    if (fPrint && norm_r_m > 3.16e0*conv)
      printf("PCGOMP: %20.9e > %20.9e!!\n", norm_r_m, conv);
/*...................................................................*/
  }
/*...................................................................*/

  timef = getTimeC() - timei;

/*...*/
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
           , nEq, nAd, tol, jG + 1, xKx, norm, norm_r, norm_r_m, timef);
  }
/*..................................................................*/

/*...*/  
  if (jG == maxIt) 
  {
    printf("PCGOMP: *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n", maxIt);
    exit(EXIT_SOLVER);
  }
/*..................................................................*/ 


/*...*/
  if (log)
    fprintf(fileLog
            , "PCGOMP: tol %20.2e "
            "iteration %d "
            "xKx %20.12e "
            "norma(x*x) %20.12e "
            "time %20.5lf\n"
            , tol, jG + 1, xKx, norm, timef);
/*...................................................................*/
}
/**********************************************************************/

/**********************************************************************
* Data de criacao    : 27/08/2019                                    *
* Data de modificaco : 00/00/0000                                    *
* -------------------------------------------------------------------*
* MPIPCGOMP: metodo do gradiente conjugado com precondiconador       *
* diagonal (M-1Ax=M-1b) (matriz simentrica)                          *
* -------------------------------------------------------------------*
* Parametros de Entrada:                                             *
* -------------------------------------------------------------------*
*  neq -> numero de equacoes                                         *
*neqNov-> numero de equacoes nao sobrepostas                         *
*  nad -> numero de elementos nao nulos fora da diagonal             *
*  nadR-> numero de elementos nao nulos na parte retangular          *
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
 * iNeq -> mapa de interface de equacoes                              *
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
void mpiPcgOmp(INT const nEq     , INT const nEqNov
            ,INT const nAd       , INT const nAdR
            , INT *RESTRICT ia   , INT *RESTRICT ja
            , DOUBLE *RESTRICT a , DOUBLE *RESTRICT ad
            , DOUBLE *RESTRICT m , DOUBLE *RESTRICT b
            , DOUBLE *RESTRICT x , DOUBLE *RESTRICT z
            , DOUBLE *RESTRICT r , DOUBLE *RESTRICT p
            , DOUBLE const tol   , unsigned int maxIt
            , bool const newX    , FILE* fileLog
            , FILE *fileHistLog  , bool const log
            , bool const fHistLog, bool const fPrint
            , BufferOmp *bOmp    ,Interface *iNeq
            , void(*matvec)()    , DOUBLE(*dot)())
{
  short nThreads = ompVar.nThreadsSolver;
  unsigned int j, jj, jG;
  INT i;
  DOUBLE alpha, beta, d, conv, xKx, norm_b, norm, norm_r_m, norm_r;
  DOUBLE di,tmp;
  DOUBLE timei, timef;
  INT param[2];
  timei = getTimeC();

  param[0] = nAd;
  param[1] = nAdR;

  norm = 0.e0;

#pragma omp parallel default(none) num_threads(nThreads)\
  private(i,j,jj,conv,norm_b,tmp,norm,d,di,alpha,beta)\
  shared(ia,ja,a,ad,m,b,x,z,r,p,maxIt,fileLog,xKx,norm_r_m,norm_r\
         ,fileHistLog,matvec,dot,bOmp,nThreads,jG,param,iNeq,mpiVar,fileLogDebug)
  {
/*... chute inicial*/
    if (newX)
#pragma omp for
      for (i = 0; i < nEqNov; i++)
        x[i] = 0.e0;
/*...................................................................*/

/*... conv = tol * |(M-1)b|m = tol(b,M(-1)b) */
#pragma omp for
    for (i = 0; i < nEqNov; i++)
      z[i] = b[i] * m[i];

    d = dot(b, z, nEqNov);
    norm_b = sqrt(fabs(d));
    conv = tol * norm_b;
/*...................................................................*/

/*... Ax0*/
    matvec(nEqNov        , param
         , ia            , ja
         , a             , ad
         , x             , z
         , bOmp->thBegin , bOmp->thEnd
         , bOmp->thHeight, bOmp->thY
         , nThreads      , iNeq);
/*...................................................................*/

/*...*/
#pragma omp for
    for (i = 0; i < nEqNov; i++)
    {
/* ... r0 = b - Ax0*/
      r[i] = b[i] - z[i];
/* ... z0 = (M-1)r0*/
      z[i] = r[i] * m[i];
/* ... p0 = r0*/
      p[i] = z[i];
    }
/*...................................................................*/

/*... (r(0), z(0)) = (r(0), (M-1)r0)*/
    d = dot(r, z, nEqNov);
/*...................................................................*/

/*...*/
    jj = 1;
    for (j = 0; j < maxIt; j++)
    {
/*... z = Ap(j)*/
      matvec(nEqNov         , param
            , ia            , ja
            , a             , ad
            , p             , z
            , bOmp->thBegin , bOmp->thEnd
            , bOmp->thHeight, bOmp->thY
            , nThreads      , iNeq);
/*...................................................................*/

/*... alpha = ( r(j),z(j) ) / ( Ap(j), p(j) ))*/
      alpha = d / dot(z, p, nEqNov);
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i < nEqNov; i++) 
      {
/*... x(j + 1) = x(j) + alpha*p*/
        x[i] += alpha * p[i];
/*... r(j + 1) = r(j) - alpha*Ap*/
        r[i] -= alpha * z[i];
/*... z = (M-1)r0*/
        z[i] = r[i] * m[i];
      }
/*...................................................................*/

/* ... (r(j + 1), (M - 1)r(j + 1)) = (r(j + 1), z)*/
      di = dot(r, z, nEqNov);
/*... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) */
      beta = di / d;
/*...................................................................*/

/* ... p(j + 1) = (M - 1)r(j + 1) + beta*p(j) = z + beta*p(j)*/
#pragma omp for
      for (i = 0; i < nEqNov; i++)
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
      if (jj == 5000) 
      {
        jj = 0;
#pragma omp master
        printf("MPIPCGOMP: %d %20.9e %20.9e\n", j + 1
                                              , sqrt(fabs(d)), conv);
      }
      jj++;
/*...................................................................*/
    }
/*...................................................................*/

#pragma omp master
    jG = j;

/*... Energy norm:  x*Kx*/
    matvec(nEqNov          , param 
           , ia            , ja
           , a             , ad
           , x             , z
           , bOmp->thBegin , bOmp->thEnd
           , bOmp->thHeight, bOmp->thY
           , nThreads      , iNeq);
/*norma de energia = xT*A*x */
    tmp = dot(x, z, nEqNov);
#pragma omp master
    xKx = tmp;
/*...................................................................*/

/*... norm - 2 = || x ||*/
    tmp  = sqrt(dot(x, x, nEqNov));
#pragma omp master
    norm = tmp;
/*...................................................................*/

/*... r = M(-1)(b - Ax) (calculo do residuo explicito)*/
#pragma omp for
    for (i = 0; i < nEqNov; i++) 
    {
      r[i] = b[i] - z[i];
      z[i] = r[i] * m[i];
    }

    tmp = dot(r, z, nEqNov);
#pragma omp master
    norm_r_m = sqrt(fabs(tmp));
    
    tmp = dot(r, r, nEqNov);
#pragma omp master
    norm_r = sqrt(tmp);

#pragma omp master
    if (fPrint && norm_r_m > 3.16e0*conv && !mpiVar.myId)
      printf("MPIPCGOMP: %20.9e > %20.9e!!\n", norm_r_m, conv);
/*...................................................................*/
  }
/*...................................................................*/

  timef = getTimeC() - timei;

/*...*/
  if(!mpiVar.myId && fPrint)
  { 
    printf(" (MPIPCG) solver:\n"
			"\tEquations      = %20d\n"
      "\tEquations      = %20d\n"
			"\tnad            = %20d\n"
      "\tnadR           = %20d\n"
			"\tSolver tol     = %20.2e\n"
			"\tIterarions     = %20d\n"
			"\tx * Kx         = %20.2e\n"
			"\t|| x ||        = %20.2e\n"
			"\t|| b - Ax ||   = %20.2e\n"
			"\tCPU time(s)    = %20.2lf\n"
	    ,nEq,nEqNov,nAd,nAdR,tol, jG +1,xKx, norm, norm_r,timef);
  }
/*..................................................................*/

/*...*/  
  if (jG == maxIt) 
  {
    if(!mpiVar.myId)
    {
      printf("MPIPCGOMP: *** WARNING: no convergende reached !\n");
      printf("MAXIT = %d \n", maxIt);
    }
    mpiStop();
    exit(EXIT_SOLVER);
  }
/*..................................................................*/ 


/*...*/
  if (!mpiVar.myId && log)
    fprintf(fileLog
            , "MPIPCGOMP: tol %20.2e "
            "iteration %d "
            "xKx %20.12e "
            "norma(x*x) %20.12e "
            "time %20.5lf\n"
            , tol, jG + 1, xKx, norm, timef);
/*...................................................................*/
}
/**********************************************************************/