#include<Solv.h>

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
  for(i = 0; i < nEq; i++)
  {
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
  for(j = 0; j < maxIt; j++)
  {
/* ... v = Az(j) = AM(-1)p(j)*/
    matvec(nEq,ia,ja,al,ad,z,v);
/*...................................................................*/

/*... alpha = (r(j), r0) / (AM(-1)p(j), r0))*/
    rr0   = dot(r,r0,nEq);
    alpha = rr0/dot(v,r0,nEq);
/*...................................................................*/

/*...*/
    for(i = 0; i < nEq; i++)
    {
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
    for(i = 0; i < nEq; i++) 
    {
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
    for(i = 0; i < nEq; i++)
    {
/*... p(j + 1) = r(i) + beta*(p(j) - w*v(i))*/
      p[i]  = r[i] + beta*(p[i]-w*v[i]);
/*... z = M(-1)p*/
      z[i]  = p[i]*m[i];
    }
/*...................................................................*/

/*...*/
		if (jj == 2000)
    {
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
/*...*/
  if(fPrint)
  { 
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
/*..................................................................*/  
 
/*...*/
  if(j == maxIt){ 
    printf(" PBICGSTAB *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
    exit(EXIT_SOLVER);
  }
/*..................................................................*/

/*...*/ 
  if(log)
    fprintf(fileLog          
           ,"PBICGSTAB: tol %20.2e " 
						"iteration %d "
						"xKx %20.12e "
						"||x*x|| %20.12e "
						"time %20.5lf\n"
						, tol, j + 1, xKx, norm, timef);
/*..................................................................*/

}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 28 / 08 / 2016                                *
 * Data de modificaco : 18 / 08 / 2019                                *
 * ------------------------------------------------------------------ *
 * MPIPBICGSTAB : metodo do gradiente conjugado bi-ortaganalizado     *
 * com precondiconador diagonal (M-1Ax=M-1b) (matriz geral)           *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 * neq -> numero de equacoes                                          *
 * neqNov-> numero de equacoes nao sobrepostas                        *
 *  nad -> numero de elementos nao nulos fora da diagonal             *
 *  nadR-> numero de elementos nao nulos na parte retangular          *
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
 * iNeq -> mapa de interface de equacoes                              *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *  
 * b[]-> alterado                                                     *
 * ad,al,au-> inalterado                                              *
 * -------------------------------------------------------------------*
 * OBS:                                                               *
 * A(M-1)y=b precondicionador a direita                               *
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
          ,DOUBLE *RESTRICT r0
          ,DOUBLE const tol    ,unsigned int maxIt
          ,bool const newX     ,FILE* fileLog 
          ,FILE *fileHistLog   ,bool const log
          ,bool const fHistLog ,bool const fPrint 
          ,Interface *iNeq    
          ,void(*matvec)()     ,DOUBLE(*dot)())
{
  unsigned int j,jj;
  INT i;
  DOUBLE alpha,beta,d,conv,xKx,norm_r,norm,norm_b,w,rr0;
  DOUBLE timei,timef;
  INT param[2];

  timei = getTimeC();
  
  param[0] = nAd;
  param[1] = nAdR;
/* chute inicial*/
  if(newX)  
    for (i = 0; i < nEq; i++)  
      x[i] = 0.e0;
/*...................................................................*/

/*... conv = tol * |b|*/	
  d        = dot(b,b,nEqNov);
	norm_b   = sqrt(d);
	conv     = tol*norm_b;
//breaktol = btol*sqrt(d);
/*...................................................................*/

/*... Ax0*/ 
  matvec(nEqNov,param,ia,ja,al,ad,x,z,iNeq);
/*...................................................................*/

/*...*/  
  for(i = 0; i < nEqNov; i++)
  {
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
  for(j = 0; j < maxIt; j++) 
  {
/* ... v = Az(j) = AM(-1)p(j)*/
    matvec(nEqNov,param,ia,ja,al,ad,z,v,iNeq);
/*...................................................................*/

/*... alpha = (r(j), r0) / (AM(-1)p(j), r0))*/
    rr0   = dot(r,r0,nEqNov);
    alpha = rr0/dot(v,r0,nEqNov);
/*...................................................................*/
  
/*...*/
    for(i = 0; i < nEqNov; i++)
    {
/*... x(j + 1) = x(j) + alpha*M(-1)p*/
      x[i] += alpha * z[i];
/*... s(j) = r(j) - alpha*AM(-1)p(j)*/
      r[i] -= alpha * v[i];
/*... z = M(-1)s*/
      z[i]  = r[i] * m[i];
    }
/*...................................................................*/

/*... (s, s)*/
		d = dot(r, r, nEqNov);
/*...*/	
		if (sqrt(d) < conv) break;
/*...................................................................*/

/*... t = Az = AM(-1)s(j)*/
    matvec(nEqNov,param,ia,ja,al,ad,z,t,iNeq);
/*... w = (AM(-1)s(j), s(j)) / (AM(-1)s(j), AM(-1)s(j))*/
    w = dot(t,r,nEqNov)/ dot(t,t,nEqNov);
/*...................................................................*/

/*...*/    
    for(i = 0; i < nEqNov; i++)
    {
/* ... x(j + 1) = x(j) + w*M(-1)s*/
      x[i] += w*z[i];
/* ... r(j + 1) = s(j) - w*AM(-1)s*/
      r[i] -= w*t[i];
    }
/*...................................................................*/

/*... (r, r)*/
    d = dot(r,r,nEqNov);
    if(sqrt(fabs(d)) < conv) break;
/*...................................................................*/

/*...*/
		if(!mpiVar.myId  && fHistLog)
			fprintf(fileHistLog, "%d %20.9e\n", j, sqrt(fabs(d)) / norm_b);
/*...................................................................*/

/*... beta = (r(j + 1), r0) / (r(j), r0)) * (alpha / w)*/
    beta = (dot(r,r0,nEqNov)/rr0)*(alpha/w);
/*...................................................................*/

/*...*/
    for(i = 0; i < nEqNov; i++)
    {
/*... p(j + 1) = r(i) + beta*(p(j) - w*v(i))*/
      p[i]  = r[i] + beta*(p[i]-w*v[i]);
/*... z = M(-1)p*/
      z[i]  = p[i]*m[i];
    }
/*...................................................................*/
 
/*...*/
		if (jj == 2000) 
    {
			jj = 0;
      if(!mpiVar.myId)
			  printf("MPIBICGSATB: %d %20.9e %20.9e\n"
              ,j + 1, sqrt(fabs(d)), conv);
		}
		jj++;
/*...................................................................*/
  }
/*...................................................................*/

/*... Energy norm:  x*Kx*/
  matvec(nEqNov,param,ia,ja,al,ad,x,z,iNeq);
/*norma de energia = xT*A*x */
  xKx  = dot(x,z,nEqNov);
/*...................................................................*/

/*... norm - 2 = || x ||*/
	norm = sqrt(dot(x, x, nEqNov));
/*...................................................................*/

/*... r = (b - Ax) (calculo do residuo explicito)*/
	for (i = 0; i < nEq; i++) {
		r[i] = b[i] - z[i];
	}
	norm_r = dot(r, r, nEq);
	norm_r = sqrt(norm_r);
	if(fPrint && norm_r > 3.16e0*conv)
		if(!mpiVar.myId) printf("MPIBICGSATB: %20.9e > %20.9e!!\n", norm_r, conv);
/*...................................................................*/
  timef = getTimeC() - timei;

/*...*/
  if(!mpiVar.myId && fPrint)
  { 
    printf(" (MPIBICGSATB) solver:\n"
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
      printf(" MPIPBICGSTAB *** WARNING: no convergende reached !\n");
      printf("MAXIT = %d \n",maxIt);
    }
    mpiStop();
    exit(EXIT_SOLVER);
  }
/*..................................................................*/ 

/*...*/
  if(!mpiVar.myId && log)
    fprintf(fileLog          
           ,"MPIPBICGSTAB: tol %20.2e " 
            "iteration %d " 
		        "xKx %20.12e "
            "||x*x|| %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1,xKx,norm,timef);
 /*..................................................................*/
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
  shared(js,ia,ja,a,ad,m,b,x,z,r,p,r0,v,t,maxIt,xKx,norm,norm_r\
         ,fileHistLog,matvec,dot,bOmp,nThreads)
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
    for (i = 0; i < nEq; i++)
    {
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
    for (j = 0; j < maxIt; j++)
    {
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
      for (i = 0; i < nEq; i++)
      {
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
      for (i = 0; i < nEq; i++)
      {
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
      for (i = 0; i < nEq; i++)
      {
/*... p(j + 1) = r(i) + beta*(p(j) - w*v(i))*/
        p[i] = r[i] + beta*(p[i] - w*v[i]);
/*... z = M(-1)p*/
        z[i] = p[i]*m[i];
      }
/*...................................................................*/

/*...*/
      if (jj == 2000)
      {
        jj = 0;
#pragma omp master
        printf("BICGSATBOMP: %d %20.9e %20.9e\n"
               , j + 1, sqrt(fabs(d)), conv);
      }
      jj++;
/*...................................................................*/
    }
/*...................................................................*/

#pragma omp master
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

#pragma omp master
    {
      norm_r = sqrt(tmp);
      if(fPrint && norm_r > 3.16e0*conv)
      printf("BICGSATB: %20.9e > %20.9e!!\n", norm_r, conv);
    }
/*...................................................................*/
  }
/*...................................................................*/

  timef = getTimeC() - timei;
/*...*/
  if (fPrint)
  {
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
/*..................................................................*/  

/*...*/
  if (js == maxIt )
  {
    printf(" PBICGSTABOMP *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n", maxIt);
    exit(EXIT_FAILURE);
  }
/*..................................................................*/  

/*...*/
  if (log)
    fprintf(fileLog
            , "PBICGSTABOMP: tol %20.2e "
            "iteration %d "
            "xKx %20.12e "
            "norma(x*x) %20.12e "
            "time %20.5lf\n"
            , tol, js + 1, xKx, norm, timef);
/*..................................................................*/  

}
/**********************************************************************/

/**********************************************************************
* Data de criacao    : 12 / 09 / 2019                                 *
* Data de modificaco : 00 / 00 / 0000                                 *
* ------------------------------------------------------------------- *
* PBICGSTABOMP : metodo do gradiente conjugado bi - ortaganalizado    *
* com precondiconador diagonal(matriz geral)                          *
* ------------------------------------------------------------------- *
* Parametros de Entrada :                                             *
* ------------------------------------------------------------------- *
* neq  -> numero de equacoes                                          *
* neqNov -> numero de equacoes nao sobrepostas                        *
* nad   -> numero de elementos nao nulos fora da diagonal             *
* nadR  -> numero de elementos nao nulos na parte retangular          *
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
* iNeq         -> mapa de interface de equacoes                       *
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
void mpiPbicgstabOmp(INT const nEq,INT const nEqNov   
               ,INT const nAd      ,INT const nAdR
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
               ,BufferOmp *bOmp    ,Interface *iNeq
               ,void(*matvec)()    ,DOUBLE(*dot)())
{
  short nThreads = ompVar.nThreadsSolver;
  unsigned int j, jj,js;
  INT i;
  DOUBLE alpha,beta,d,conv,xKx,norm_r,norm,norm_b,w,rr0,tmp;
  DOUBLE timei, timef;
  INT param[2];

  timei = getTimeC();

  param[0] = nAd;
  param[1] = nAdR;

#pragma omp parallel default(none) num_threads(nThreads)\
  private(i,j,jj,conv,norm_b,tmp,d,alpha,beta,rr0,w)\
  shared(js,ia,ja,a,ad,m,b,x,z,r,p,r0,v,t,maxIt,xKx,norm,norm_r\
         ,fileHistLog,matvec,dot,bOmp,nThreads,param,iNeq)
  {
/*... chute inicial*/
    if (newX)
#pragma omp for
      for (i = 0; i < nEq; i++)
        x[i] = 0.e0;
/*...................................................................*/

/*... conv = tol * |b|*/
    d      = dot(b, b, nEqNov);
    norm_b = sqrt(d);
    conv   = tol*norm_b;
  //breaktol = btol*sqrt(d);
/*...................................................................*/

/*... Ax0*/
    matvec(nEqNov        ,param  
          ,ia            ,ja 
          ,a             ,a
          ,x             ,z
          ,bOmp->thBegin ,bOmp->thEnd
          ,bOmp->thHeight,bOmp->thY
          ,nThreads      ,iNeq);
/*...................................................................*/

/*...*/
#pragma omp for
    for (i = 0; i < nEqNov; i++)
    {
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
    for (j = 0; j < maxIt; j++)
    {
/* ... v = Az(j) = AM(-1)p(j)*/
      matvec(nEqNov        ,param
            ,ia            ,ja
            ,a             ,ad
            ,z             ,v
            ,bOmp->thBegin ,bOmp->thEnd
            ,bOmp->thHeight,bOmp->thY
            ,nThreads      ,iNeq);
/*...................................................................*/

/*... alpha = (r(j), r0) / (AM(-1)p(j), r0))*/
      rr0   = dot(r, r0, nEqNov);
      alpha = rr0/ dot(v, r0, nEqNov);
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i < nEqNov; i++)
      {
/*... x(j + 1) = x(j) + alpha*M(-1)p*/
        x[i] += alpha * z[i];
/*... s(j) = r(j) - alpha*AM(-1)p(j)*/
        r[i] -= alpha * v[i];
/*... z = M(-1)s*/
        z[i] = r[i]*m[i];
      }
/*...................................................................*/

/*... (s, s)*/
      d = dot(r, r, nEqNov);
/*...*/
      if (sqrt(d) < conv) break;
/*...................................................................*/

/*... t = Az = AM(-1)s(j)*/
      matvec(nEqNov        ,param
            ,ia            ,ja
            ,a             ,ad
            ,z             ,t
            ,bOmp->thBegin ,bOmp->thEnd
            ,bOmp->thHeight,bOmp->thY
            ,nThreads      ,iNeq);
/*...................................................................*/

/*... w = (AM(-1)s(j), s(j)) / (AM(-1)s(j), AM(-1)s(j))*/
      w = dot(t,r,nEqNov)/dot(t,t,nEqNov);
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i < nEqNov; i++)
      {
/* ... x(j + 1) = x(j) + w*M(-1)s*/
        x[i] += w*z[i];
/* ... r(j + 1) = s(j) - w*AM(-1)s*/
        r[i] -= w*t[i];
      }
/*...................................................................*/

/*... (r, r)*/
      d = dot(r,r,nEqNov);
      if (sqrt(d) < conv) break;
/*...................................................................*/

/*...*/
#pragma omp master
      if (fHistLog)
        fprintf(fileHistLog, "%d %20.9e\n", j, sqrt(fabs(d)) / norm_b);
/*...................................................................*/

/*... beta = (r(j + 1), r0) / (r(j), r0)) * (alpha / w)*/
      beta = (dot(r,r0,nEqNov)/rr0)*(alpha/w);
/*...................................................................*/

/*...*/
#pragma omp for
      for (i = 0; i < nEqNov; i++)
      {
/*... p(j + 1) = r(i) + beta*(p(j) - w*v(i))*/
        p[i] = r[i] + beta*(p[i] - w*v[i]);
/*... z = M(-1)p*/
        z[i] = p[i]*m[i];
      }
/*...................................................................*/

/*...*/
      if (jj == 2000)
      {
        jj = 0;
#pragma omp master
        printf("MPIBICGSATBOMP: %d %20.9e %20.9e\n"
               , j + 1, sqrt(fabs(d)), conv);
      }
      jj++;
/*...................................................................*/
    }
/*...................................................................*/

#pragma omp master
    js = j;

/*... Energy norm:  x*Kx*/
    matvec(nEqNov        ,param
          ,ia            ,ja
          ,a             ,ad
          ,x             ,z
          ,bOmp->thBegin ,bOmp->thEnd
          ,bOmp->thHeight,bOmp->thY
          ,nThreads);
/*norma de energia = xT*A*x */
    tmp = dot(x, z, nEqNov);
#pragma omp master
    xKx = tmp; 
/*...................................................................*/

/*... norm - 2 = || x ||*/
    tmp = sqrt(dot(x, x, nEqNov));
#pragma omp master
    norm = tmp;
/*...................................................................*/

/*... r = (b - Ax) (calculo do residuo explicito)*/
#pragma omp for
    for (i = 0; i < nEqNov; i++) {
      r[i] = b[i] - z[i];
    }

    tmp = dot(r, r, nEqNov);

#pragma omp master
    {
      norm_r = sqrt(tmp);
      if(fPrint && norm_r > 3.16e0*conv)
      printf("MPIBICGSATOMP: %20.9e > %20.9e!!\n", norm_r, conv);
    }
/*...................................................................*/
  }
/*...................................................................*/

  timef = getTimeC() - timei;
/*...*/
  if (fPrint)
  {
    printf(" (MPIBICGSATBOMP) solver:\n"
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
	    ,nEq,nEqNov,nAd,nAdR,tol,js+1,xKx, norm, norm_r,timef);
  }
/*..................................................................*/  

/*...*/
  if (js == maxIt )
  {
    printf(" MPIPBICGSTABOMP *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n", maxIt);
    mpiStop();
    exit(EXIT_SOLVER);
  }
/*..................................................................*/  

/*...*/
  if (!mpiVar.myId && log)
    fprintf( fileLog
            , "MPIPBICGSTABOMP: tol %20.2e "
            "iteration %d "
            "xKx %20.12e "
            "norma(x*x) %20.12e "
            "time %20.5lf\n"
            , tol, js + 1, xKx, norm, timef);
/*..................................................................*/  

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

/*...*/
  if (fPrint)
  {
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
/*..................................................................*/  

/*...*/
  if (j == maxIt) {
    printf(" PBICGSTAB(2) *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n", maxIt);
    exit(EXIT_FAILURE);
  }
/*..................................................................*/

/*...*/
  if (log)
    fprintf(fileLog
            , "PBICGSTAB(2): tol %20.2e "
            "iteration %d "
            "xKx %20.12e "
            "norma(x*x) %20.12e "
            "time %20.5lf\n"
            , tol, j + 1, xKx, norm, timef);
/*..................................................................*/

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
  shared(js,ia,ja,a,ad,m,b,x,t,v,r,u,r0,w,s,p,h,z,maxIt,fileLog\
         ,xKx,norm,norm_r,fileHistLog,matvec,dot,bOmp,nThreads)
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
      for (i = 0; i<nEq; i++)
      {
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
      for (i = 0; i<nEq; i++)
      {
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
      for (i = 0; i<nEq; i++) 
      {
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
      for (i = 0; i<nEq; i++) 
      {
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
      for (i = 0; i < nEq; i++) 
      {
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
      if (jj == 2000) 
      {
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
    for (i = 0; i < nEq; i++) 
    {
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
/*...*/
  if (fPrint) 
  {
    printf(" (BICGSATBOMP(2) solver:\n"
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
/*..................................................................*/  

/*...*/
  if (js == maxIt)
  {
    printf(" PBICGSTABOMP(2) *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n", maxIt);
    exit(EXIT_FAILURE);
  }
/*..................................................................*/

/*...*/
  if (log)
    fprintf(fileLog
            , "PBICGSTABOMP(2): tol %20.2e "
            "iteration %d "
            "xKx %20.12e "
            "norma(x*x) %20.12e "
            "time %20.5lf\n"
            , tol, js + 1, xKx, norm, timef);
 /*..................................................................*/
}
/**********************************************************************/