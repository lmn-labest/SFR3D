#include<Solv.h>
/*********************************************************************
c * Metodos iterativos para solucao de sistemas lineares              *
c * ----------------------------------------------------------------  *
c * simetricos:                                                       *
c * ----------------------------------------------------------------  *
c * CG - gradiente conjugados                                         *
c *                                                                   *
c * PCG - gradiente conjugados com precondicionador diagonal          *
c *                                                                   *
c * ----------------------------------------------------------------  *
c * nao - simetricos:                                                 *
c * ----------------------------------------------------------------  *
c * bicgstab - gradiente bi - conjugados estabilizados                *
c *                                                                   *
c * pbicgstab - gradiente bi - conjugados estabilizados  com          *
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
 *  au				-> parte superior da matriz A                           *
 *   p				-> precondiconador diagonal                             *
 *   b				-> vetor b (Ax=b)                                       *
 *   x				-> vetor de solucao                                     *
 *   z				-> vetor auxiliar                                       *
 *   r				-> vetor auxiliar                                       *
 *   p				-> vetor auxiliar                                       *
 * tol				-> tolerancia de convergencia                           *
 * newX				-> vetor inicial iniciado com zero                      *
 * fLog				-> arquivo de log do solver                             *
 * log				-> log de arquivo (true|false)                          *
 * tol				-> tolerancia do solver                                 *
 * maxIt			-> numero maximo de iteracoes                           *
 *  newX			-> true zero o vetor inicial                            *
 * fLog				-> arquivo de saida do log                              *
 * fileHistLog-> arquivo de log da iteracoes                          *
 *  log       -> escreve o log do solver                              *
 * fPrint			-> saida de informacao na tela                          *
 * fHistLog		-> log das iteracoes                                    *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *  
 * b[]-> inalterado                                                   *
 * ad,al,au-> inalterado                                              *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 * ------------------------------------------------------------------ *
*********************************************************************/
void pcg(INT const nEq, INT const nAd
	, INT *restrict ia, INT *restrict ja
	, DOUBLE *restrict al, DOUBLE *restrict ad, DOUBLE *restrict au
	, DOUBLE *restrict m, DOUBLE *restrict b, DOUBLE *restrict x
	, DOUBLE *restrict z, DOUBLE *restrict r, DOUBLE *restrict p
	, DOUBLE const tol
	, unsigned int maxIt,bool const newX
	, FILE* fLog        ,FILE *fileHistLog
	, bool const log    ,bool const fHistLog, bool const fPrint
	, void(*matvec)(), DOUBLE(*dot)())
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
    fprintf(fLog          
           ,"PCG: tol %20.2e " 
            "iteration %d " 
						"xKx %20.12e "
		        "norma(x*x) %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1, xKx,norm,timef);
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
    fprintf(fLog          
           ,"MPIPCG: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);
#endif
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
 *  al        -> parte inferior da matriz A                           *
 *  ad        -> diagnal da matriz A                                  *
 *  au        -> parte superior da matriz A                           *
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
 * fLog       -> arquivo de log do solver                             *
 * log        -> log de arquivo (true|false)                          *
 * tol        -> tolerancia do solver                                 *
 * maxIt			-> numero maximo de iteracoes                           *
 *  newX			-> true zero o vetor inicial                            *
 * fLog				-> arquivo de saida do log                              *
 * fileHistLog-> arquivo de log da iteracoes                          *
 *  log       -> escreve o log do solver                              *
 * fPrint			-> saida de informacao na tela                          *
 * fHistLog		-> log das iteracoes                                    *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *  
 * b[]-> inalterado                                                   *
 * ad,al,au-> inalterado                                              *
 * -------------------------------------------------------------------*
 * OBS:                                                               *
 * A(M-1)y=b precondicionador a direita                               *
 * -------------------------------------------------------------------*
***********************************************************************/
void pbicgstab(INT const nEq  ,INT const nAd
          ,INT *restrict ia   ,INT *restrict ja
          ,DOUBLE *restrict al,DOUBLE *restrict ad,DOUBLE *restrict au
          ,DOUBLE *restrict m ,DOUBLE *restrict b ,DOUBLE *restrict x
          ,DOUBLE *restrict t ,DOUBLE *restrict v ,DOUBLE *restrict r
          ,DOUBLE *restrict p ,DOUBLE *restrict z ,DOUBLE *restrict r0
					,DOUBLE const tol
          ,unsigned int maxIt ,bool const newX          
          ,FILE* fLog         ,FILE *fileHistLog
          ,bool const log     ,bool const fHistLog, bool const fPrint
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
    fprintf(fLog          
           ,"PBICGSTAB: tol %20.2e " 
						"iteration %d "
						"xKx %20.12e "
						"norma(x*x) %20.12e "
						"time %20.5lf\n"
						, tol, j + 1, xKx, norm, timef);


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
    fprintf(fLog          
           ,"MPIPBICGSTAB: tol %20.2e " 
            "iteration %d " 
		        "norma %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,energy,timef);

}
/**********************************************************************/



