#include<Solv.h>
/**********************************************************************
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 22/07/2016                                    *
 * -------------------------------------------------------------------*
 * SOLVERC: resolucao do sistema linear(Ax=b)                         *
 * -------------------------------------------------------------------*
 * Parametro de entrada                                               *
 * -------------------------------------------------------------------*
 * m       -> vetor de memoria                                        *
 * neq     -> numero de equacoes                                      *
 * neqNov  -> numero de equacoes nao sobrepostas                      *
 * nad     -> numero de elemetos nao nulos fora da diagonal principal *
 * nAdR    -> numero de termos nao nulos na parte retangular         *
 * ia      -> vetor das linhas da estrutura de dados                  *
 * ja      -> vetor das colunas da estrutura de dados                 *
 * al      -> parte inferior da matriz a                              *
 * ad      -> diagonal principal da matriz a                          *
 * au      -> parte superior da matriz a                              *
 *  b      -> vetor de forcas                                         *
 * iNeq    -> mapa de interface de equacoes                           *
 *  x      -> vetor da solucao                                        *
 * tol     -> tolerancia do solver                                    *
 * maxit   -> numero maximo de iteracao do solver iterativo           *
 * storage -> tecnica de armazenamento da matriz esparsa              * 
 * solver  -> PCG                                                     *
 * fSolvLog-> arquivo de log para solver                              *
 * fLog    -> log de arquivo (true|false)                             *
 * newX    -> vetor inicial iniciado com zero                         *
 * openmp  -> flag do openmp true|false                               *
 * unsym   -> matriz nao simetrica                                    *
 * loopwise-                                                          *
 * -------------------------------------------------------------------*
 * Parametro de saida :                                               * 
 * -------------------------------------------------------------------*
 * x      - atualizado                                                *
 * b      - modificado                                                *
 * -------------------------------------------------------------------*
 * OBS: csrD e csrC são iguais para matrizes simetricas               *
**********************************************************************/
void solverC(Memoria *m    
            ,INT const nEq       ,INT const nEqNov  
            ,INT const nAd       ,INT const nAdR
            ,INT *ia             ,INT *ja   
            ,DOUBLE *al          ,DOUBLE *ad,DOUBLE *au
            ,DOUBLE *b           ,DOUBLE *x
            ,Interface *iNeq       
            ,DOUBLE const tol    ,unsigned int maxIt
            ,short const storage ,short const solver
            ,FILE* fSolvLog      ,bool const fLog
            ,bool const newX     ,bool const openMp   
            ,bool const unSym    ,bool const loopWise)
{
  DOUBLE *z=NULL,*r=NULL,*pc=NULL,*t=NULL,*v=NULL,*p=NULL,*h=NULL;
  void   (*matVecC)();
  DOUBLE (*dotC)();
/*...*/
	dotC    = NULL;
	matVecC = NULL;
/*...................................................................*/

  switch(solver){
/*... gradientes conjugados com precondicionador diagonal*/
    case PCG:
/*... precondiconador diagonal*/
      HccaAlloc(DOUBLE,m,pc,nEqNov,"pc",false);
      zero(pc,nEqNov,DOUBLEC);
      tm.precondDiag = getTimeC() - tm.precondDiag;
      preCondDiag(pc,ad,nEqNov);
      tm.precondDiag = getTimeC() - tm.precondDiag;
/*...................................................................*/

/*... arranjos auxiliares do pcg*/
      HccaAlloc(DOUBLE,m,z,nEq,"z",false);
      HccaAlloc(DOUBLE,m,r,nEq,"r",false);
			HccaAlloc(DOUBLE,m,p,nEq,"p", false);
      zero(z,nEq,DOUBLEC);
      zero(r,nEq,DOUBLEC);
			zero(p,nEq, DOUBLEC);
/*...................................................................*/
      
/*...*/
			setMatVec(&matVecC,storage,unSym);
/*...................................................................*/

/*...*/
			setDot(&dotC, DOT);
/*...................................................................*/

/*... gradientes conjugados*/
      tm.pcg = getTimeC() - tm.pcg;
/*... PCG-MPI*/
      if( mpiVar.nPrcs > 1)
        mpiPcg(nEq      ,nEqNov    
              ,nAd      ,nAdR
              ,ia       ,ja
              ,al       ,ad   ,au
              ,pc       ,b    ,x
              ,z        ,r    ,tol
              ,maxIt    ,true 
              ,fSolvLog ,fLog
              ,false 
              ,iNeq           
              ,matVecC ,dotC);   
/*...................................................................*/
      
/*... PCG*/
      else
        pcg(nEq     ,nAd
           ,ia      ,ja
           ,al      ,ad   ,au
           ,pc      ,b    ,x
           ,z       ,r    ,p
           ,tol
           ,maxIt   ,newX 
           ,fSolvLog,NULL
           ,fLog    ,false   ,false
           ,matVecC ,dotC);   
/*...................................................................*/
      tm.pcg = getTimeC() - tm.pcg;
/*...................................................................*/
      
/*... liberando arranjos auxiliares do pcg*/
			HccaDealloc(m,p,"p",false);
      HccaDealloc(m,r,"r",false);
      HccaDealloc(m,z,"z",false);
/*...................................................................*/

/*... liberando arranjos do precondicionador*/
      HccaDealloc(m,pc,"pc",false);
/*...................................................................*/
    break;
/*...................................................................*/

/*... gradientes conjugados bi-ortoganilizado com precondicionador
     diagonal*/
    case PBICGSTAB:
/*... precondiconador diagonal*/
      HccaAlloc(DOUBLE,m,pc,nEqNov,"pc",false);
      zero(pc,nEqNov,DOUBLEC);
      tm.precondDiag = getTimeC() - tm.precondDiag;
      preCondDiag(pc,ad,nEqNov);
      tm.precondDiag = getTimeC() - tm.precondDiag;
/*...................................................................*/

/*... arranjos auxiliares do pbicgstab*/
      HccaAlloc(DOUBLE,m,z,nEq,"z",false);
      HccaAlloc(DOUBLE,m,r,nEq,"r",false);
      HccaAlloc(DOUBLE,m,t,nEq,"tt",false);
      HccaAlloc(DOUBLE,m,v,nEq,"vv",false);
      HccaAlloc(DOUBLE,m,p,nEq,"pp",false);
			HccaAlloc(DOUBLE,m,h,nEq,"hh", false);
      zero(z,nEq,DOUBLEC);
      zero(r,nEq,DOUBLEC);
      zero(t,nEq,DOUBLEC);
      zero(v,nEq,DOUBLEC);
      zero(p,nEq,DOUBLEC);
			zero(h,nEq,DOUBLEC);
/*...................................................................*/

/*...*/
			setMatVec(&matVecC, storage, unSym);
/*...................................................................*/
      
/*...*/
			setDot(&dotC,DOT);
/*...................................................................*/

/*... gradientes conjugados bi-ortoganilizado*/
      tm.pbicgstab = getTimeC() - tm.pbicgstab;
/*... PBICGSTAB-MPI*/
      if( mpiVar.nPrcs > 1)  
        mpiPbicgstab(nEq   ,nEqNov
                 ,nAd      ,nAdR
                 ,ia       ,ja
                 ,al       ,ad    ,au
                 ,pc       ,b     ,x
                 ,t        ,v     ,r
                 ,p        ,z     ,tol
                 ,maxIt    ,true          
                 ,fSolvLog ,fLog
                 ,false
                 ,iNeq           
                 ,matVecC  ,dotC);  
/*...................................................................*/
      
/*... PBICGSTAB*/
      else
        pbicgstab(nEq      ,nAd
                 ,ia       ,ja
                 ,al       ,ad    ,au
                 ,pc       ,b     ,x
                 ,t        ,v     ,r
                 ,p        ,z     ,h
                 ,tol
                 ,maxIt    ,newX          
                 ,fSolvLog ,NULL
                 ,fLog     ,false,false    
                 ,matVecC  ,dotC);
/*...................................................................*/
      tm.pbicgstab = getTimeC() - tm.pbicgstab;
/*...................................................................*/
      
/*... liberando arranjos auxiliares do pbicgstab*/
			HccaDealloc(m,h,"hh",false);
      HccaDealloc(m,p,"pp",false);
      HccaDealloc(m,v,"vv",false);
      HccaDealloc(m,t,"tt",false);
      HccaDealloc(m,r,"r" ,false);
      HccaDealloc(m,z,"z" ,false);
/*...................................................................*/

/*... liberando arranjos do precondicionador*/
      HccaDealloc(m,pc,"pc",false);
/*...................................................................*/
    break;
/*...................................................................*/

/*...*/
    default:
      ERRO_OP(__FILE__,__func__,solver);
    break;
/*...................................................................*/
 
  }
/*...................................................................*/
}
/*********************************************************************/      

/**********************************************************************
 *SETSOLVER : escolhe o solver                                        *
 **********************************************************************/
void setSolver(char *word,short *solver)
{

  if(!strcmp(word,"PCG"))
   *solver = PCG;
  else if(!strcmp(word,"PBICGSTAB"))
   *solver = PBICGSTAB;

} 
/*********************************************************************/      


/**********************************************************************
 *SMACHN: calcula a precisao da maquina para double                   *
 **********************************************************************/
DOUBLE smachn(){
  
  DOUBLE mq = 1.e0;

  while(mq+1.e0 != 1.e0)
    mq/=2.e0;

  return mq;
}
/*********************************************************************/      

/**********************************************************************
* Data de criacao :    22 / 07 / 2016                                 *
* Data de modificaco : 00 / 00 / 0000  															  *
* ------------------------------------------------------------------- *
* SETDOT :escolhe o produto interno                                   *
* ------------------------------------------------------------------- *
* Parametro de entrada                                                *
* ------------------------------------------------------------------- *
* dotC -> -> nao definido                                             *
* iCod -> versao do produto interno  																	*
* ------------------------------------------------------------------- *
* Parametro de saida :																							  *
* ------------------------------------------------------------------- *
* dotC -> funcao escolhida                                            *
* ------------------------------------------------------------------- *
* OBS :                                                               *
**********************************************************************/
void setDot(DOUBLE(**dotC)(),short const iCod) {

/*...*/
	switch(iCod){
		case DOT:
			*dotC = dot;
			break;
		case DOTL2:
	    *dotC    = dotL2;
			break;
		case DOTL4:
	    *dotC    = dotL4;
			break;
		case DOTL6:
	    *dotC    = dotL6;
			break;
		case DOTL8:
	    *dotC    = dotL8;
			break;
		case DOTO2:
	    *dotC    = dotO2;
			break;
		case DOTO4:
	    *dotC    = dotO4;
			break;
		case  DOTO6:
	    *dotC    = dotO6;
			break;
		case  DOTO8:
	    *dotC    = dotO8;
			break;
		case  DOTO2L2:
	    *dotC    = dotO2L2;
		break;
    default:
			ERRO_OP(__FILE__, __func__, iCod);
		break;	
	}
/*...................................................................*/
}
/*********************************************************************/

/**********************************************************************
* Data de criacao :    22 / 07 / 2016                                 *
* Data de modificaco : 00 / 00 / 0000  															  *
* ------------------------------------------------------------------- *
* SETMATVEC: escolhe o produto matriz vetor desejado                  *
* ------------------------------------------------------------------- *
* Parametro de entrada                                                *
* ------------------------------------------------------------------- *
* matVecC -> nao definido                                             *
* storage -> tecnica de armazenamento																	*
* unSym   -> matriz nao-simetrica (true:false 												*
* ------------------------------------------------------------------- *
* Parametro de saida :																							  *
* ------------------------------------------------------------------- *
* matVecC -> funcao escolhida                                         *
* ------------------------------------------------------------------- *
* OBS :                                                               *
**********************************************************************/
void setMatVec(void (**matVecC)(), short const storage
              ,bool const unSym) {
	
	bool fCoo = false;
/*... estrutura de dados de armazenamentos da matriz esparcas*/
	switch (storage) {
/*... armazenamento CSR(a)*/
	case CSR:
		*matVecC = NULL;
		break;
/*...................................................................*/

/*... armazenamento CSRD(a,ad)*/
/*...CSRD+COO(symetric)*/
	case CSRDCOO:
		fCoo = true;
	case CSRD:
/*... nao - simetrica*/
		if (unSym) {
/*... mpi*/
			if (mpiVar.nPrcs > 1) {
				*matVecC = mpiMatVecCsrD;
			}
/*...................................................................*/

/*... sequencial*/
			else {
				*matVecC = matVecCsrD;
			}
/*...................................................................*/
		}
/*...................................................................*/

/*... simetrica*/
		else {
/*... mpi*/
			if (mpiVar.nPrcs > 1) {
/*... CSRD+COO*/
				if (fCoo) {
					*matVecC = mpiMatVecCsrDcooSym;
				}
/*...................................................................*/

/*... CSRD+CSR*/
				else {
					*matVecC = mpiMatVecCsrDSym;
/*...................................................................*/
				}
			}
/*..................................................................*/

/*... sequencial*/
			else {
				*matVecC = matVecCsrDSym;
			}
/*..................................................................*/
		}
/*..................................................................*/
		break;
/*...................................................................*/

/*... armazenamento CSRC(ad,au,al)*/
/*...CSRC+COO*/
	case CSRCCOO:
		fCoo = true;
	case CSRC:
/*... nao - simetrica*/
		if (unSym) {
/*... mpi*/
			if (mpiVar.nPrcs > 1) {
/*... CSRC+COO*/
				if (fCoo) {
					*matVecC = mpiMatVecCsrCcoo;
				}
/*...................................................................*/

/*... CSRC+CSR*/
				else {
					*matVecC = mpiMatVecCsrC;
/*...................................................................*/
				}
			}
/*..................................................................*/
		
/*... sequencial*/
			else {
				*matVecC = matVecCsrC;
			}
/*..................................................................*/
    }		
/*..................................................................*/

/*... simetrica*/
		else {
/*... mpi*/
			if (mpiVar.nPrcs > 1) {
/*... CSRD+COO*/
				if (fCoo) {
					*matVecC = mpiMatVecCsrDcooSym;
				}
/*...................................................................*/

/*... CSRD+CSR*/
				else {
					*matVecC = mpiMatVecCsrDSym;
/*...................................................................*/
				}
			}
/*..................................................................*/

/*... sequencial*/
			else {
				*matVecC = matVecCsrDSym;
			}
/*..................................................................*/
		}
/*..................................................................*/
  break;
/*...................................................................*/

/*... armazenamento ELLPACK(ad,a)*/
	case ELLPACK:
		/*... mpil*/
		if (mpiVar.nPrcs > 1) {
			*matVecC = mpiMatVecEllPack;
		}
/*...................................................................*/

/*... sequencial*/
		else {
			*matVecC = matVecEllPack;
			//          matVecC = matVecEllPackO2;
			//          matVecC = matVecEllPackO4;
		}
/*...................................................................*/
		break;
/*...................................................................*/

/*...*/
	default:
		ERRO_OP(__FILE__, __func__, storage);
		break;
/*...................................................................*/
	}
/*...................................................................*/

}
/*********************************************************************/
