#include<Solv.h>
void symOrtho(DOUBLE const a, DOUBLE const b
                   , DOUBLE *c    , DOUBLE *s
                   , DOUBLE *r);

/**********************************************************************
 * Data de criacao    : 18/09/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 * -------------------------------------------------------------------*
 * MINRES : Solucao de sistemas de equacoes pelo metodo MINRES        *
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
 * v0(neq)  - arranjo local de trabalho                               *
 * w(neq)   - arranjo local de trabalho                               *
 * w0(neq)  - arranjo local de trabalho                               *
 * w00(neq) - arranjo local de trabalho                               *
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
 * fonte: Iterative Krylov Methods for Large Linear Systems           *
 * Henk A. van de Vorst                                               *
 * ------------------------------------------------------------------ * 
*********************************************************************/
void minres(INT const nEq        , INT const nAd
	        , INT *RESTRICT ia    , INT *RESTRICT ja
	        , DOUBLE *RESTRICT al , DOUBLE *RESTRICT ad
	        , DOUBLE *RESTRICT b  , DOUBLE *RESTRICT x 
          , DOUBLE *RESTRICT v0 , DOUBLE *RESTRICT v
          , DOUBLE *RESTRICT w  , DOUBLE *RESTRICT w0
          , DOUBLE *RESTRICT w00, DOUBLE *RESTRICT z 
          , DOUBLE *RESTRICT p
	        , DOUBLE const tol   , unsigned int maxIt
          , bool const newX    , FILE* fileLog   
          , FILE *fileHistLog	, bool const log 
          , bool const fHistLog, bool const fPrint
	        , void(*matvec)()    , DOUBLE(*dot)())
{
  unsigned int j, jj, it;
  INT i;
  DOUBLE tmp1, tmp2, conv, beta, beta_old, neta
       , c, c_old, s, s_old, ro1, ro2, ro3, alpha, delta, xKx
       , norm_r, norm_b, norm;
       
  DOUBLE timei, timef;
	timei = getTimeC();
  it       = 0;
  ro1 = ro2 = ro3 = 1.e0;
/*... chute inicial*/
	if (newX)
		for (i = 0; i < nEq; i++)
			x[i] = 0.e0;
/*...................................................................*/

/*...tol |b|*/
  norm_b = sqrt(dot(b, b, nEq));
  conv   = tol * sqrt(norm_b);
/*...................................................................*/

/*... Ax0*/
	matvec(nEq, ia, ja, al, ad, x, z);
/*...................................................................*/

/*...*/
	for (i = 0; i < nEq; i++) {
/* ... v1 = b - Ax0*/
		v[i] = b[i] - z[i];
/* ...v0 = w0 = w1 = 0.0*/
		v0[i] = w0[i] = w00[i] = 0.e0;
	}
/*...................................................................*/

/*... beta = sqrt (v, v) = (v, v)*/
	beta = sqrt(fabs(dot(v, v, nEq)));
  norm_r = beta;
/*...................................................................*/

/*...*/
  neta     = beta;
  c_old    = 1.e0;
  c        = 1.e0;
  s_old    = 0.e0;
  s        = 0.e0;
/*...................................................................*/

/*...*/
	jj = 1;
	for (j = 0; j < maxIt; j++) {
    it ++;
/*... The Lancos recurrence:*/
    tmp1 = 1.e0/beta;
    for (i = 0; i < nEq; i++){
/*... v(j) = v(j)/beta*/
      v[i] = v[i]*tmp1;
    }
/*.....................................................................*/

/*... z = Av(j)*/
		matvec(nEq, ia, ja, al, ad, v, z);
/*...................................................................*/

/*... alpha = ( v(j), z  ) = ( v(j), Av(j))*/
    alpha = dot(z, v, nEq)*tmp1;
/*...................................................................*/

/*...*/
		for (i = 0; i < nEq; i++) {
/* ... v(j+1) = Az(j) - alpha(j)*v(j) - beta(j)*z(j-1)*/
      z[i] = z[i] - alpha * v[i] - beta * v0[i];
		}
/*...................................................................*/

/* ... beta(j+1) = ||v(j+1)||*/
    beta_old  = beta;
    beta      = sqrt(dot(z,z,nEq));
/* ..................................................................*/

/*... QR part:*/
/*... delta = c(j)*alfa(j) - c(j-1)*s(j)*beta(j)*/
    delta = c*alpha - c_old*s*beta_old;
/*... ro2 = s(j)*alpha(j) + c(j-1)*c(j)*beta(j)*/
    ro2 = s*alpha + c_old*c*beta_old;
/*... ro3 = s(j-1)*beta(j)*/
    ro3 = s_old*beta_old;
/*...................................................................*/  
 
/*...*/
    c_old = c;
    s_old = s;
/* ..................................................................*/

/* ... New Givens rotation for subdiag element:*/
    symOrtho(delta,beta,&c,&s,&ro1);
/* ..................................................................*/

/*... Update of solution (W = VR^-1)*/ 
    tmp1 = 1.e0/ro1;
    tmp2 = c*neta;
    for (i = 0; i < nEq; i++) {
/*... w(j) = (v(j) - ro3*w(j-2) - ro2*w(j-1))/ro1*/
      w[i] = (v[i] - ro3*w00[i] - ro2*w0[i])*tmp1;  
/*... x(j) = x(j-1) + c(j+1)*neta*w(i)*/
      x[i] += tmp2*w[i];

      v0[i] = v[i];
      v[i]  = p[i];

      w00[i] = w0[i];
      w0[i]  = w[i];  
    }                
/*...................................................................*/  

/*... ||r(j) || = |c(j+1)| || r(j-1) ||*/
    norm_r = fabs(s)*norm_r; 
    if ( norm_r < conv) break;  
/*...................................................................*/  

/*...*/
		if (fHistLog)
		  fprintf(fileHistLog, "%d %20.9e\n", j, norm_r / norm_b);
/*...................................................................*/
 
/*... neta = - s(j+1) * neta*/
    neta = - s*neta;
/*...................................................................*/  

/*...*/
		if (jj == 2000) {
  		jj = 0;
			printf("PMINRES: %d %20.9e %20.9e\n", j+1, norm_r, conv);
		}
		jj++;
/*...................................................................*/
	}
/*...................................................................*/
/*
  if( nRestart != 0 ){
    nRestart--;
    goto RESTART;
  }
*/
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
		v[i] = b[i] - z[i];
	}
	norm_r   = sqrt(dot(v, v, nEq));
	if(fPrint && norm_r > 3.16e0*conv)
	  printf("PMINRES: %20.9e > %20.9e!!\n", norm_r, conv);
/*...................................................................*/

  timef = getTimeC() - timei;

  if(fPrint){ 
		printf(" (PMINRES) solver:\n"
					"\tEquations            = %20d\n"
          "\tnad                  = %20d\n"        
          "\tSolver tol           = %20.2e\n"
          "\ttol * || b ||        = %20.2e\n"     
          "\tIterarions           = %20d\n"
					"\tx * Kx               = %20.2e\n"
				  "\t|| x ||              = %20.2e\n"
					"\t|| b - Ax ||         = %20.2e\n"
	        "\tCPU time(s)          = %20.2lf\n" 
	         ,nEq,nAd,tol,conv,j+1,xKx,norm,norm_r,timef);
  }
  
  if(j == maxIt){ 
    printf("PMINRES: *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
//  exit(EXIT_FAILURE);
  }
  if(log)
    fprintf(fileLog          
           ,"PMINRES: tol %20.2e " 
            "iteration %d " 
						"xKx %20.12e "
		        "norma(x*x) %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1, xKx,norm,timef);
}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 18/09/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 * -------------------------------------------------------------------*
 * PMINRES :                                                          *
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
 * v0(neq)  - arranjo local de trabalho                               *
 * w(neq)   - arranjo local de trabalho                               *
 * w0(neq)  - arranjo local de trabalho                               *
 * w00(neq) - arranjo local de trabalho                               *
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
 * fonte:                                                             *
 * 1. Iterative Krylov Methods for Large Linear Systems               *
 * Henk A. van de Vorst                                               *
 * 2. Precondicionador Iterative methods for singular linear          *  
 * equations and lest-squares problems                                * 
 *( Sou-Cheng (Terrya) Choi - 2006                                    *
 * D(-1/2)AD(-1/2)x = D(-1/2)b                                        *
 * ------------------------------------------------------------------ * 
*********************************************************************/
void pminres(INT const nEq      , INT const nAd
	         , INT *RESTRICT ia   , INT *RESTRICT ja
	         , DOUBLE *RESTRICT al, DOUBLE *RESTRICT ad
	         , DOUBLE *RESTRICT m , DOUBLE *RESTRICT b
           , DOUBLE *RESTRICT x , DOUBLE *RESTRICT v0
           , DOUBLE *RESTRICT v , DOUBLE *RESTRICT w
           , DOUBLE *RESTRICT w0, DOUBLE *RESTRICT w00
           , DOUBLE *RESTRICT z , DOUBLE *RESTRICT z0
           , DOUBLE *RESTRICT p
	         , DOUBLE const tol   , unsigned int maxIt
           , bool const newX    , FILE* fileLog   
           , FILE *fileHistLog	, bool const log 
           , bool const fHistLog, bool const fPrint
	         , void(*matvec)()    , DOUBLE(*dot)())
{
  unsigned int j, jj, it;
  INT i;
  DOUBLE tmp1, tmp2, tmp3, tmp4, conv, beta, beta_old, neta
       , c, c_old, s, s_old, ro1, ro2, ro3, alpha, delta, xKx
       , norm_r, norm_b, norm, norm_r_m;
       
  DOUBLE timei, timef;
	timei = getTimeC();
  it       = 0;
  ro1 = ro2 = ro3 = 1.e0;
/*... chute inicial*/
	if (newX)
		for (i = 0; i < nEq; i++)
			x[i] = 0.e0;
/*...................................................................*/

/*...*/
  norm_b = sqrt(dot(b, b, nEq));
/*...................................................................*/

/*... conv = tol * |M(-1/2)b| */
	for (i = 0; i < nEq; i++)
		z[i] = b[i] * sqrt(m[i]);

	tmp1   = dot(z, z, nEq);
	conv   = tol * sqrt(tmp1);
/*...................................................................*/

/*... Ax0*/
	matvec(nEq, ia, ja, al, ad, x, z);
/*...................................................................*/

/*...*/
	for (i = 0; i < nEq; i++) {
/* ... v1 = b - Ax0*/
		v[i] = b[i] - z[i];
/* ...v0 = w0 = w1 = 0.0*/
		v0[i] = w0[i] = w00[i] = 0.e0;
/* ...  z  =M(-1)v*/
		z0[i] = v[i]*m[i];
	}
/*...................................................................*/

/*... beta = sqrt (v, z) = (v, (M-1)v)*/
	beta = sqrt(fabs(dot(v, z0, nEq)));
  norm_r = beta;
/*...................................................................*/

/*...*/
  neta     = beta;
  beta_old = beta;
  c_old    = 1.e0;
  c        = 1.e0;
  s_old    = 0.e0;
  s        = 0.e0;
/*...................................................................*/

/*...*/
	jj = 1;
	for (j = 0; j < maxIt; j++) {
    it ++;
/*... p = Az(j)*/
		matvec(nEq, ia, ja, al, ad, z0, p);
/*...................................................................*/

/*... alpha = ( z(j),p(j) ) / beta^2*/
    tmp1 = 1.e0/(beta*beta);
    alpha = dot(z0, p, nEq)*tmp1;
/*...................................................................*/

/*...*/
    tmp1 =  1.0e0/beta;
    tmp2 = -alpha/beta;
    tmp3 = -beta/beta_old;
		for (i = 0; i < nEq; i++) {
/* ... v(j+1) = Az(j) - alpha(j)*v(j) - beta(j)*z(j-1)*/
      p[i] = tmp1*p[i] + tmp2* v[i] + tmp3 * v0[i];
/* ... z  =M(-1)v*/
      z[i] = p[i]*m[i];
		}
/*...................................................................*/

/* ... beta(j+1) = raiz( (zt(j+1),v(j+1)) )*/
    beta_old  = beta;
    tmp1      = dot(z,p,nEq);
    beta      = sqrt(tmp1);
/* ..................................................................*/

/*... QR part:*/
/*... delta = c(j)*alfa(j) - c(j-1)*s(j)*beta(j)*/
    delta = c*alpha - c_old*s*beta_old;
/*... ro2 = s(j)*alpha(j) + c(j-1)*c(j)*beta(j)*/
    ro2 = s*alpha + c_old*c*beta_old;
/*... ro3 = s(j-1)*beta(j)*/
    ro3 = s_old*beta_old;
/*...................................................................*/  
 
/*...*/
    c_old = c;
    s_old = s;
/* ..................................................................*/

/* ... New Givens rotation for subdiag element:*/
    symOrtho(delta,beta,&c,&s,&ro1);
//  ro1 = sqrt(delta*delta + beta*beta);
//... c(j+1) = delta*ro1
//  c = delta/ro1;
//... s(j+1) = beta(j+1)/ro1
//  s = beta/ro1;
/* ..................................................................*/

/*... Update of solution (W = VR^-1)*/ 
    tmp1 = 1.e0/(beta_old*ro1);
    tmp2 = -ro2/ro1;       
    tmp3 = -ro3/ro1;
    tmp4 = c*neta;
    for (i = 0; i < nEq; i++) {
/*    w(j) = (v(j)/beta(j) - ro2*w(j-1) - ro3*w(j-2))/ro1*/
      w[i] = tmp1*z0[i] + tmp2*w0[i] + tmp3*w00[i];  
/*... x(j) = x(j-1) + c(j+1)*neta*w(i)*/
      x[i] += tmp4*w[i];

      v0[i] = v[i];
      v[i]  = p[i];

      z0[i] = z[i];

      w00[i] = w0[i];
      w0[i]  = w[i];  
    }                
/*...................................................................*/  

/*... ||r(j) || = |c(j+1)| || r(j-1) ||*/
    norm_r = fabs(s)*norm_r; 
    if ( norm_r < conv) break;  
/*...................................................................*/  

/*...*/
		if (fHistLog)
		  fprintf(fileHistLog, "%d %20.9e\n", j, norm_r / norm_b);
/*...................................................................*/
 
/*... neta = - s(j+1) * neta*/
    neta = - s*neta;
/*...................................................................*/  

/*...*/
		if (jj == 10000) {
  		jj = 0;
			printf("PMINRES: %d %20.9e %20.9e\n", j+1
                                             , norm_r, conv);
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
		v[i] = b[i] - z[i];
		z[i] = v[i] * sqrt(fabs(m[i]));
	}
	norm_r   = sqrt(dot(v, v, nEq));
	norm_r_m = sqrt(dot(z, z, nEq));
	if(fPrint && norm_r_m > 3.16e0*conv)
	  printf("PMINRES: %20.9e > %20.9e!!\n", norm_r_m, conv);
/*...................................................................*/

  timef = getTimeC() - timei;

  if(fPrint){ 
		printf(" (PMINRES) solver:\n"
					"\tEquations            = %20d\n"
          "\tnad                  = %20d\n"        
          "\tSolver tol           = %20.2e\n"
          "\ttol * || M(-1/2)b || = %20.2e\n"     
          "\tIterarions           = %20d\n"
					"\tx * Kx               = %20.2e\n"
				  "\t|| x ||              = %20.2e\n"
					"\t|| b - Ax ||         = %20.2e\n"
					"\t|| b - Ax || m       = %20.2e\n"
	        "\tCPU time(s)          = %20.2lf\n" 
	         ,nEq,nAd,tol,conv,j+1,xKx,norm,norm_r,norm_r_m,timef);
  }
  
  if(j == maxIt){ 
    printf("PMINRES: *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
//  exit(EXIT_FAILURE);
  }
  if(log)
    fprintf(fileLog          
           ,"PMINRES: tol %20.2e " 
            "iteration %d " 
						"xKx %20.12e "
		        "norma(x*x) %20.12e "
		        "time %20.5lf\n"
           ,tol,j+1, xKx,norm,timef);
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

/*...*/
     if( jj == 1000){
        jj = 0;
        printf("GMRES(%d): %d %d %20.9e %20.9e\n",nKrylov,l,nIt
             ,fabs(e[ni+1]),eConv);
      }
      jj++;
/*...................................................................*/

/*... Verifica a convergencia:*/
    if (fabs(e[ni+1]) < eConv) break;
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
  int i, j, ni, nIt;
  unsigned short l,nCol = nKrylov;
  DOUBLE *g1, *g2;
  DOUBLE tmp, norm, norm_r, eConv, beta, h1, h2, aux1, aux2, r, xKx;
  INT iLong;
  DOUBLE timei, timef;

  timei = getTimeC();

#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#endif

/*... chute inicial*/
  if (newX)
    for (i = 0; i < nEq; i++)
      x[i] = 0.e0;
/*...................................................................*/

/*...g(1,i) = (M-1)*b*/
#pragma omp parallel default(none) private(tmp,iLong)\
        shared(g,b,m,norm,dot)
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
  for (l = 0; l < nCycles; l++) {
/*... Residuo g1 = b - Ax*/
#pragma omp parallel default(none) private(tmp)\
        shared(ia,ja,a,ad,x,g,b,m,e,bOmp,nThreads,matvec,dot)
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
        shared(ia,ja,a,ad,g1,g2,m,bOmp,nThreads,matvec,dot)
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
       shared(g1,g2,beta,dot) 
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
      shared(norm,h,g1,g2,beta,dot,ni,nCol)
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
#pragma omp parallel for default (none) private(iLong) shared(tmp,x,g,j)
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
        shared(ia,ja,a,ad,x,g,bOmp,nThreads,xKx,norm,matvec,dot)
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

/*********************************************************************** 
  * Data de criacao    : 18/09/2017                                    *
  * Data de modificaco : 00/00/0000                                    * 
  * ------------------------------------------------------------------ *
  * SYM_ORTHO: Givens rotation                                         * 
  * (versa com melhor comportamento numerico)                          *      
  * ------------------------------------------------------------------ * 
  * Parametros de entrada:                                             *
  * ------------------------------------------------------------------ *
  * a   - paramentro                                                   *
  * b   - paramentro s da diagonal principal                           *
  * c   - nao definido                                                 *
  * s   - nao definido                                                 *
  * r   - nao definido                                                 *
  * ------------------------------------------------------------------ * 
  * Parametros de saida:                                               *
  * ------------------------------------------------------------------ * 
  * c   - cos(teta)                                                    *
  * s   - seno(teta)                                                   *
  * r   - raiz(a^2 + b^2)                                              *
  * ------------------------------------------------------------------ * 
  * OBS:                                                               *
  * | c  s | | a |    | raiz(a^2 + b^2) |   | r |                      *
    |      | |   |  = |                 | = |   |                      *
  * |-s  c | | b |    |        0        |   | 0 |                      *
  * ------------------------------------------------------------------ *
  *********************************************************************/
      void symOrtho(DOUBLE const a, DOUBLE const b
                   , DOUBLE *c    , DOUBLE *s
                   , DOUBLE *r)
{
      DOUBLE ma, mb, sa, sb, t;
/*...*/      
      ma = fabs(a);
      mb = fabs(b);
      sa = SING1(a);
      sb = SING1(b);
/*.....................................................................*/

/*...*/     
     if(b == 0.e0){
        *s = 0.e0;
        *r = ma;
        if ( a == 0.e0)
          *c = 1.0e0;
        else
          *c = sa;
     } 
     else if( a == 0.e0){
       *c = 0.e0;
       *s = sb;
       *r = mb;
     } 
     else if( mb > a){
       t = a/b;
       *s = sb/sqrt(1.e0 + t*t);
       *c = *s*t;
       *r = b / *s;
     }
     else if( ma > mb ){
       t = b / a;
       *c = sa/sqrt(1.e0 + t*t);
       *s = *c * t;
       *r = a/ *c;
     }
/*.....................................................................*/
 }
/**********************************************************************/ 
