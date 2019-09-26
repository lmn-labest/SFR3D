#include<Solv.h>
/**********************************************************************
c * Metodos iterativos para solucao de sistemas lineares              *
c * ----------------------------------------------------------------  *
c * simetricos:                                                       *
c * ----------------------------------------------------------------  *
c *                                                                   *
c * PCG - gradiente conjugados com precondicionador diagonal          *
c *                                                                   * 
c * MINRES - MINRES com precondicionador diagonal M=D(1/2)D(1/2)      *
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
c * gmres      - gmres com precondicionador diagonal                  *
c *                                                                   *
c *********************************************************************/

/**********************************************************************
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 27/08/2016                                    *
 * -------------------------------------------------------------------*
 * SOLVERC: resolucao do sistema linear(Ax=b)                         *
 * -------------------------------------------------------------------*
 * Parametro de entrada                                               *
 * -------------------------------------------------------------------*
 * m       -> vetor de memoria                                        *
 * neq     -> numero de equacoes                                      *
 * neqNov  -> numero de equacoes nao sobrepostas                      *
 * nad     -> numero de elemetos nao nulos fora da diagonal principal *
 * nAdR    -> numero de termos nao nulos na parte retangular          *
 * ia      -> vetor das linhas da estrutura de dados                  *
 * ja      -> vetor das colunas da estrutura de dados                 *
 * al      -> parte inferior da matriz a                              *
 * ad      -> diagonal principal da matriz a                          *
 * au      -> parte superior da matriz a                              *
 *  b      -> vetor de forcas                                         *
 * iNeq    -> mapa de interface de equacoes                           *
 * bOmp    -> Openmp                                                  *
 *  x      -> vetor da solucao                                        *
 * tol     -> tolerancia do solver                                    *
 * maxit   -> numero maximo de iteracao do solver iterativo           *
 * storage -> tecnica de armazenamento da matriz esparsa              * 
 * solver  -> PCG                                                     *
 * fSolvLog-> arquivo de log para solver                              *
 * fLog    -> log de arquivo (true|false)                             *
 * newX    -> vetor inicial iniciado com zero                         *
 * unsym   -> matriz nao simetrica                                    *
 * -------------------------------------------------------------------*
 * Parametro de saida :                                               * 
 * -------------------------------------------------------------------*
 * x      - atualizado                                                *
 * b      - modificado                                                *
 * -------------------------------------------------------------------*
 * OBS: csrD e csrC sao iguais para matrizes simetricas               *
 **********************************************************************/
void solverC(Memoria *m    
            ,INT const nEq       ,INT const nEqNov  
            ,INT const nAd       ,INT const nAdR
            ,INT *ia             ,INT *ja   
            ,DOUBLE *al          ,DOUBLE *ad,DOUBLE *au
            ,DOUBLE *b           ,DOUBLE *x
            ,Interface *iNeq     ,BufferOmp *bOmp
            ,DOUBLE const tol    ,unsigned int maxIt
            ,short const storage ,short const solver
            ,FILE* fSolvLog      ,bool const fLog
            ,bool const newX     ,bool const unSym    )
{
  DOUBLE *z=NULL,*r=NULL,*pc=NULL,*t=NULL,*v=NULL,*p=NULL,*h=NULL
        ,*i=NULL,*o=NULL,*s=NULL,*d=NULL;
  void   (*matVecC)();
  DOUBLE (*dotC)();
  INT nn;
  unsigned short nKrylov=15;
  short mType;
  bool fPrint = false,openMp = ompVar.fSolver;
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
			HccaAlloc(DOUBLE,m,p,nEq,"p",false);
      zero(z,nEq,DOUBLEC);
      zero(r,nEq,DOUBLEC);
			zero(p,nEq, DOUBLEC);
/*...................................................................*/
      
/*...*/
			setMatVec(&matVecC,storage,unSym,openMp);
/*...................................................................*/

/*...*/
			setDot(&dotC, DOT);
/*...................................................................*/

/*... gradientes conjugados*/
      tm.pcg = getTimeC() - tm.pcg;
      callCg(nEq        ,nEqNov
             ,nAd        ,nAdR
             ,ia         ,ja
             ,al         ,ad
             ,pc         ,b   
             ,x          ,z      
             ,r          ,p
             ,tol
             ,maxIt      ,newX
             ,fSolvLog   ,fLog
             ,fPrint     
             ,iNeq       ,bOmp
             ,matVecC    ,dotC);
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
			setMatVec(&matVecC, storage, unSym,openMp);
/*...................................................................*/
      
/*...*/
			setDot(&dotC,DOT);
/*...................................................................*/

/*... gradientes conjugados bi-ortoganilizado*/
      tm.pbicgstab = getTimeC() - tm.pbicgstab;
/*...*/
      callBicgStab(nEq   ,nEqNov
                 ,nAd    ,nAdR
                 ,ia     ,ja
                 ,al     ,ad    
                 ,pc     ,b     
                 ,x      ,t  
                 ,v      ,r
                 ,p      ,z
                 ,h 
                 ,tol    ,maxIt 
                 ,newX   ,fSolvLog
                 ,fLog   ,fPrint
                 ,iNeq   ,bOmp           
                 ,matVecC,dotC);  
/*...................................................................*/

/*...*/
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

/*... gradientes conjugados bi-ortoganilizado com precondicionador
      diagonal*/
    case PBICGSTABL2:
/*... precondiconador diagonal*/
      HccaAlloc(DOUBLE, m, pc, nEqNov, "pc", false);
      zero(pc, nEqNov, DOUBLEC);
      tm.precondDiag = getTimeC() - tm.precondDiag;
      preCondDiag(pc, ad, nEqNov);
      tm.precondDiag = getTimeC() - tm.precondDiag;
/*...................................................................*/

/*... arranjos auxiliares do pbicgstab*/
      HccaAlloc(DOUBLE, m, z, nEq, "z", false);
      HccaAlloc(DOUBLE, m, r, nEq, "r", false);
      HccaAlloc(DOUBLE, m, t, nEq, "tt", false);
      HccaAlloc(DOUBLE, m, v, nEq, "vv", false);
      HccaAlloc(DOUBLE, m, p, nEq, "pp", false);
      HccaAlloc(DOUBLE, m, h, nEq, "hh", false);
      HccaAlloc(DOUBLE, m, i, nEq, "ii", false);
      HccaAlloc(DOUBLE, m, o, nEq, "oo", false);
      HccaAlloc(DOUBLE, m, s, nEq, "ss", false);
      HccaAlloc(DOUBLE, m, d, nEq, "dd", false);
      zero(z, nEq, DOUBLEC);
      zero(r, nEq, DOUBLEC);
      zero(t, nEq, DOUBLEC);
      zero(v, nEq, DOUBLEC);
      zero(p, nEq, DOUBLEC);
      zero(h, nEq, DOUBLEC);
      zero(i, nEq, DOUBLEC);
      zero(o, nEq, DOUBLEC);
      zero(s, nEq, DOUBLEC);
      zero(d, nEq, DOUBLEC);
/*...................................................................*/

/*...*/
      setMatVec(&matVecC, storage, unSym, openMp);
/*...................................................................*/

/*...*/
      setDot(&dotC, DOT);
/*...................................................................*/

/*... gradientes conjugados bi-ortoganilizado*/
      tm.pbicgstab = getTimeC() - tm.pbicgstab;
/*...*/
      callBicgStabl2(nEq  , nEqNov
                    ,nAd ,nAdR
                    ,ia  ,ja
                    ,al  ,ad
                    ,pc  ,b
                    ,x   ,t
                    ,v   ,r
                    ,p   ,z
                    ,h   ,i
                    ,o   ,s
                    ,d     
                    ,tol , maxIt
                    ,newX, fSolvLog
                    ,fLog, fPrint
                    ,iNeq, bOmp
                    ,matVecC, dotC);
/*...................................................................*/

/*...*/
      tm.pbicgstab = getTimeC() - tm.pbicgstab;
/*...................................................................*/

/*... liberando arranjos auxiliares do pbicgstab*/
      HccaDealloc(m, d, "dd", false);
      HccaDealloc(m, s, "ss", false);
      HccaDealloc(m, o, "oo", false);
      HccaDealloc(m, i, "ii", false);
      HccaDealloc(m, h, "hh", false);
      HccaDealloc(m, p, "pp", false);
      HccaDealloc(m, v, "vv", false);
      HccaDealloc(m, t, "tt", false);
      HccaDealloc(m, r, "r", false);
      HccaDealloc(m, z, "z", false);
/*...................................................................*/

/*... liberando arranjos do precondicionador*/
      HccaDealloc(m, pc, "pc", false);
/*...................................................................*/
    break;
/*...................................................................*/

/*... GMRES*/
    case GMRES:
/*... precondiconador diagonal*/
      HccaAlloc(DOUBLE, m, pc, nEqNov, "pc", false);
      zero(pc, nEqNov, DOUBLEC);
/*...*/
      tm.precondDiag = getTimeC() - tm.precondDiag;
      preCondDiag(pc, ad, nEqNov);
      tm.precondDiag = getTimeC() - tm.precondDiag;
/*...................................................................*/

/*... arranjos auxiliares do pbicgstab*/
      nn = nEq*(nKrylov+1);
      HccaAlloc(DOUBLE, m, z, nn, "gg", false);
      zero(z, nn, DOUBLEC);
      
      nn = nKrylov*(nKrylov + 1);
      HccaAlloc(DOUBLE, m, h, nn,"hh", false);
      zero(h, nn, DOUBLEC);
      
      nn = nKrylov+1;
      HccaAlloc(DOUBLE, m, p, nn, "ee", false);
      zero(p, nn, DOUBLEC);
      
      nn = nKrylov;
      HccaAlloc(DOUBLE,m,t,nn,"cc", false);
      HccaAlloc(DOUBLE,m,v,nn,"ss", false);
      HccaAlloc(DOUBLE,m,d,nn,"yy", false);
      zero(t,nn,DOUBLEC);
      zero(v,nn,DOUBLEC);
      zero(d,nn,DOUBLEC);
/*...................................................................*/

/*...*/
      setMatVec(&matVecC, storage, unSym, openMp);
/*...................................................................*/

/*...*/
      setDot(&dotC, DOT);
/*...................................................................*/

/*... Gmres*/
      tm.gmres = getTimeC() - tm.gmres;
/*...*/
      callGmres(nEq    ,nEqNov
               ,nAd    ,nAdR
               ,ia     ,ja
               ,al     ,ad
               ,pc     ,b
               ,x      ,z
               ,h      ,d
               ,t      ,v
               ,p      ,nKrylov
               ,tol    ,maxIt
               ,newX   ,fSolvLog
               ,fLog   ,fPrint
               ,iNeq   ,bOmp
               ,matVecC,dotC);
/*...................................................................*/

/*...*/
      tm.gmres = getTimeC() - tm.gmres;
/*...................................................................*/

/*... liberando arranjos auxiliares do pbicgstab*/
       HccaDealloc(m, d, "yy", false);
       HccaDealloc(m, v, "ss", false);
       HccaDealloc(m, t, "cc", false);
       HccaDealloc(m, p, "ee", false);
       HccaDealloc(m, h, "hh", false);
       HccaDealloc(m, z, "gg", false);
/*...................................................................*/

/*... liberando arranjos do precondicionador*/
      HccaDealloc(m, pc, "pc", false);
/*...................................................................*/
    break;
/*...................................................................*/ 

/*... MINRES com precondicionador diagonal*/
    case MINRES:
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
      HccaAlloc(DOUBLE,m,o,nEq,"oo", false);
      HccaAlloc(DOUBLE,m,d,nEq,"dd", false);
      zero(z,nEq,DOUBLEC);
      zero(r,nEq,DOUBLEC);
      zero(t,nEq,DOUBLEC);
      zero(v,nEq,DOUBLEC);
      zero(p,nEq,DOUBLEC);
			zero(h,nEq,DOUBLEC);
      zero(o,nEq,DOUBLEC);
      zero(d,nEq,DOUBLEC);
/*...................................................................*/

/*...*/
			setMatVec(&matVecC, storage, unSym,openMp);
/*...................................................................*/
      
/*...*/
			setDot(&dotC,DOT);
/*...................................................................*/

/*... gradientes conjugados bi-ortoganilizado*/
      tm.minres = getTimeC() - tm.minres;
/*...*/
      callMinRes(nEq    , nEqNov
                ,nAd    , nAdR
                ,ia     , ja
                ,al     , ad    
                ,pc     , b     
                ,x      , t  
                ,v      , r
                ,p      , z
                ,h      , o
                ,d
                ,tol    , maxIt 
                ,newX   , fSolvLog
                ,fLog   , fPrint
                ,iNeq   , bOmp           
                ,matVecC, dotC);  
/*...................................................................*/

/*...*/
      tm.minres = getTimeC() - tm.minres;
/*...................................................................*/
      
/*... liberando arranjos auxiliares do pbicgstab*/
      HccaDealloc(m,d,"dd",false);
      HccaDealloc(m,o,"oo",false);
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

/*... PARDISO MKL*/
    case PARDISO:

/*... arranjos auxiliares*/
      HccaAlloc(DOUBLE,m,z,nEq,"z",false);
      HccaAlloc(DOUBLE,m,r,nEq,"r",false);
      zero(z,nEq,DOUBLEC);
      zero(r,nEq,DOUBLEC);
/*...................................................................*/

/*...*/
      if(unSym)
        mType = 3;
      else
        mType = 2;
/*...................................................................*/

/*... pardiso*/
      tm.pardiso = getTimeC() - tm.pardiso;
      callMklPardiso(nEq    , mType
                   , ia     , ja
                   , ad     , b
                   , x      , z
                   , r
                   , false);  
      tm.pardiso = getTimeC() - tm.pardiso;
/*...................................................................*/
      
/*... liberando arranjos auxiliares*/
      HccaDealloc(m,r,"r" ,false);
      HccaDealloc(m,z,"z" ,false);
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
 *SetSolver : escolhe o solver                                        *
 **********************************************************************/
void setSolver(char *word,short *solver)
{

  if(!strcmp(word,"cg"))
    *solver = PCG;
  else if(!strcmp(word,"bicgstab"))
    *solver = PBICGSTAB;
  else if (!strcmp(word, "bicgstabl1"))
    *solver = PBICGSTABL2;
  else if (!strcmp(word, "gmres"))
    *solver = GMRES;
 else if (!strcmp(word, "minres"))
    *solver = MINRES;
  else if (!strcmp(word, "pardiso"))
    *solver = PARDISO;
  else
  {
    ERRO_OP_WORD(fileLogDebug,__FILE__,__func__,__LINE__
                ,"Solver config", word, EXIT_SOLVER_CONFIG);
  }

} 
/*********************************************************************/      

/********************************************************************* 
 * Data de criacao    : 25/09/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * SetSolverConfig : escolhe o solver                                *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * f - arquivo                                                       *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void setSolverConfig(char *fileName,Solv *solv,short *storage)
{
  FILE *fileIn=NULL;
  short j=0;
  char str[] = {"end"};
  char word[WORD_SIZE];
  char macro[][WORD_SIZE] = { "name" ,"tol"             /*0,1*/
                             ,"maxit","storage"         /*2,3*/   
                             ,""     ,""                /*4,5*/
                            };   
  
  fileIn = openFile(fileName,"r");

/*...*/
  readMacroV2(fileIn, word, false, true);
  do
    {
      j = 0;
/*... name*/
      if (!strcmp(word, macro[j++]))
      {
        readMacroV2(fileIn, word, false, true);
        setSolver(word   ,&solv->solver);
        if(!mpiVar.myId )
          fprintf(fileLogExc,"%-20s : %s\n","Solver",word);
      }
/*.....................................................................*/

/*... tol*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn,"%lf",&solv->tol);
        if(solv->tol == 0.e0) 
          solv->tol = smachn();
        if(!mpiVar.myId )
           fprintf(fileLogExc,"%-20s : %e\n","tol",solv->tol);
      }
/*.....................................................................*/

/*... maxIt*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn,"%u" ,&solv->maxIt);
        if(!mpiVar.myId )
           fprintf(fileLogExc,"%-20s : %u\n","MaxIt",solv->maxIt);
      }
/*.....................................................................*/

/*... storage*/
      else if (!strcmp(word, macro[j++]))
      {
        readMacroV2(fileIn, word, false, true);
        setDataStruct(word, storage);
      }
/*.....................................................................*/

      readMacroV2(fileIn, word, false, true);
    }while(strcmp(word,str));
/*.....................................................................*/
    fclose(fileIn);

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
      if(ompVar.fSolver)
        *dotC = dotOmp;
      else
			  *dotC = dot;
			break;
		case DOTI2:
      if (ompVar.fSolver)
        *dotC = dotOmpI2;
      else
	      *dotC = dotI2;
			break;
		case DOTI4:
      if (ompVar.fSolver)
        *dotC = dotOmpI4;
      else
	      *dotC = dotI4;
			break;
		case DOTI6:
      if (ompVar.fSolver)
        *dotC = dotOmpI6;
      else
	      *dotC = dotI6;
			break;
		case DOTI8:
      if (ompVar.fSolver)
        *dotC = dotOmpI8;
      else
	      *dotC = dotI8;
			break;
		case DOTO2:
      if (ompVar.fSolver)
        *dotC = dotOmpO2;
      else
	      *dotC = dotO2;
			break;
		case DOTO4:
      if (ompVar.fSolver)
        *dotC = dotOmpO4;
      else
	      *dotC = dotO4;
			break;
		case  DOTO6:
      if (ompVar.fSolver)
        *dotC = dotOmpO6;
      else
	      *dotC = dotO6;
			break;
		case  DOTO8:
      if (ompVar.fSolver)
        *dotC = dotOmpO8;
      else
	      *dotC = dotO8;
			break;
		case  DOTO2I2:
      if (ompVar.fSolver)
        *dotC = dotOmpO2I2;
      else
	      *dotC = dotO2I2;
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
* Data de modificaco : 27 / 08 / 2019  															  *
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
void setMatVec(void (**matVecC)(),short const storage
              ,bool const unSym  ,bool const openMp) {
	
	bool fCoo = false;
/*... estrutura de dados de armazenamentos da matriz esparcas*/
	switch (storage) {
/*... armazenamento CSR(a)*/
	case CSR:
    if(unSym)
      *matVecC = matVecCsr;
    else
		  *matVecC = matVecCsrSym;
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
/*... Omp*/
        if (ompVar.fSolver)
          *matVecC = matVecCsrDomp;
/*... sequencial*/
        else
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
/*... Omp*/
        if (ompVar.fSolver)
          *matVecC = matVecCsrDsymOmp;
/*... sequencial*/
        else
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
		if (unSym) 
    {
/*... mpi*/
			if (mpiVar.nPrcs > 1)
      {
/*... CSRC+COO*/
				if (fCoo) 
        {
					*matVecC = mpiMatVecCsrCcoo;
				}
/*...................................................................*/

/*... CSRC+CSR*/
				else
        {
/*... Omp*/
          if (ompVar.fSolver)
            *matVecC = mpiMatVecCsrComp;
/*... sequencial*/
          else
  					*matVecC = mpiMatVecCsrC;
/*...................................................................*/
				}
			}
/*..................................................................*/
		
/*... sem MPI*/
			else {
/*... Omp*/
        if (ompVar.fSolver)
          *matVecC = matVecCsrComp;
/*... sequencial*/
        else
				  *matVecC = matVecCsrC;
			}
/*..................................................................*/
    }		
/*..................................................................*/

/*... simetrica*/
		else {
/*... mpi*/
			if (mpiVar.nPrcs > 1)
      {
/*... CSRD+COO*/
				if (fCoo)
        {
					*matVecC = mpiMatVecCsrDcooSym;
				}
/*...................................................................*/

/*... CSRD+CSR*/
				else
        {
          if(ompVar.fSolver)
            *matVecC = mpiMatVecCsrDSymOmp;
          else
				    *matVecC = mpiMatVecCsrDSym;
/*...................................................................*/
				}
			}
/*..................................................................*/

/*... sequencial*/
			else 
      {
        if(ompVar.fSolver)
          *matVecC = matVecCsrDsymOmp;
        else
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

/**********************************************************************
* Data de criacao    : 27/08/2016                                    *
* Data de modificaco : 27/08/2019                                    *
* -------------------------------------------------------------------*
* CALLCG : chama o gradiente conjugados escolhido                    *
* -------------------------------------------------------------------*
* Parametro de entrada                                               *
* -------------------------------------------------------------------*
* -------------------------------------------------------------------*
* Parametro de saida :                                               *
* -------------------------------------------------------------------*
* -------------------------------------------------------------------*
* OBS:                                                               *
**********************************************************************/
void callCg(INT const nEq      ,INT const nEqNov
            ,INT const nAd      ,INT const nAdR
            ,INT *RESTRICT ia   ,INT *RESTRICT ja
            ,DOUBLE *RESTRICT al,DOUBLE *RESTRICT ad
            ,DOUBLE *RESTRICT m ,DOUBLE *RESTRICT b 
            ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT z
            ,DOUBLE *RESTRICT r ,DOUBLE *RESTRICT p
            ,DOUBLE const tol   ,unsigned int maxIt 
            ,bool const newX    ,FILE* fSolvLog     
            ,bool const fLog    ,bool const fPrint  
            ,Interface *iNeq    ,BufferOmp *bOmp
            ,void(*matVec)()    ,DOUBLE(*dot)())
{

/*... PCG-MPI*/
  if (mpiVar.nPrcs > 1)
  {
/*... OpenMp*/
    if(ompVar.fSolver)
    {
/*    mpiPcg(nEq     ,nEqNov
          ,nAd     ,nAdR
          ,ia      ,ja
          ,al      ,ad 
          ,m       ,b  
          ,x       ,z
          ,r       ,p
          ,tol     ,maxIt
          ,newX    ,fSolvLog
          ,NULL    ,fLog
          ,false   ,fPrint
          ,iNeq
          ,matVec  ,dot);*/
    }
/*...................................................................*/

/*... sequencial*/
    else
    {
      mpiPcg(nEq     ,nEqNov
          ,nAd     ,nAdR
          ,ia      ,ja
          ,al      ,ad 
          ,m       ,b  
          ,x       ,z
          ,r       ,p
          ,tol     ,maxIt
          ,newX    ,fSolvLog
          ,NULL    ,fLog
          ,false   ,fPrint
          ,iNeq
          ,matVec  ,dot);
    }
/*...................................................................*/
  }
/*... PCG*/
  else
  {
/*... OpenMp*/
    if(ompVar.fSolver)
    {
      pcgOmp(nEq     ,nAd
            ,ia      ,ja
            ,al      ,ad
            ,m       ,b
            ,x       ,z
            ,r       ,p
            ,tol     ,maxIt 
            ,newX    ,fSolvLog
            ,NULL    ,fLog    
            ,false   ,fPrint
            ,bOmp
            ,matVec  ,dot);
    }
/*...................................................................*/

/*... sequencial*/
    else
    {  
      pcg(nEq     ,nAd
         ,ia      ,ja
         ,al      ,ad  
         ,m       ,b    
         ,x       ,z   
         ,r       ,p
         ,tol     ,maxIt  
         ,newX    ,fSolvLog
         ,NULL    ,fLog 
         ,false   ,fPrint
         ,matVec  ,dot);
    }  
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/**********************************************************************
* Data de criacao    : 27/08/2016                                    *
* Data de modificaco : 27/08/2019                                    *
* -------------------------------------------------------------------*
* CALLBICGCSTAB : chama o gradiente biconjugados estabilizados       *                   *
* -------------------------------------------------------------------*
* Parametro de entrada                                               *
* -------------------------------------------------------------------*
* -------------------------------------------------------------------*
* Parametro de saida :                                               *
* -------------------------------------------------------------------*
* -------------------------------------------------------------------*
* OBS:                                                               *
**********************************************************************/
void callBicgStab(INT const nEq     ,INT const nEqNov
                 ,INT const nAd     ,INT const nAdR
                 ,INT *RESTRICT ia  ,INT *RESTRICT ja
                 ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                 ,DOUBLE *RESTRICT m,DOUBLE *RESTRICT b
                 ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT t
                 ,DOUBLE *RESTRICT v,DOUBLE *RESTRICT r
                 ,DOUBLE *RESTRICT p,DOUBLE *RESTRICT z
                 ,DOUBLE *RESTRICT h
                 ,DOUBLE const tol  ,unsigned int maxIt
                 ,bool const newX   ,FILE* fSolvLog
                 ,bool const fLog   ,bool const fPrint
                 ,Interface *iNeq   ,BufferOmp *bOmp
                 ,void(*matVec)()   ,DOUBLE(*dot)())
{

/*... MPI*/
  if (mpiVar.nPrcs > 1)
/*... OpenMp*/
    if (ompVar.fSolver)
       mpiPbicgstabOmp(nEq    ,nEqNov
                  ,nAd    ,nAdR
                  ,ia    ,ja
                  ,a     ,ad
                  ,m     ,b
                  ,x     ,t
                  ,v     ,r
                  ,p     ,z
                  ,h
                  ,tol   ,maxIt
                  ,newX  ,fSolvLog
                  ,NULL  ,fLog
                  ,false ,fPrint
                  ,bOmp  ,iNeq
                  ,matVec,dot);
/*...................................................................*/

/*... sequencial*/
    else
      mpiPbicgstab(nEq    ,nEqNov
                ,nAd    ,nAdR
                ,ia     ,ja
                 ,a     ,ad
                 ,m     ,b 
                 ,x     ,t   
                 ,v     ,r
                 ,p     ,z 
                 ,h
                 ,tol   ,maxIt  
                 ,newX  ,fSolvLog
                 ,NULL  ,fLog  
                 ,false ,fPrint
                 ,iNeq  
                 ,matVec,dot);
/*...................................................................*/

/*... */
  else
  {
/*... OpenMp*/
    if (ompVar.fSolver)
    {
      pbicgstabOmp(nEq   ,nAd
                  ,ia    ,ja
                  ,a     ,ad
                  ,m     ,b
                  ,x     ,t
                  ,v     ,r
                  ,p     ,z
                  ,h
                  ,tol   ,maxIt
                  ,newX  ,fSolvLog
                  ,NULL  ,fLog
                  ,false ,fPrint
                  ,bOmp
                  ,matVec,dot);
    }
/*...................................................................*/

/*... sequencial*/
    else
    {
      pbicgstab(nEq   ,nAd
               ,ia    ,ja
               ,a     ,ad
               ,m     ,b 
               ,x     ,t   
               ,v     ,r
               ,p     ,z 
               ,h
               ,tol   ,maxIt  
               ,newX  ,fSolvLog
               ,NULL  ,fLog  
               ,false ,fPrint
               ,matVec,dot);
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/**********************************************************************
* Data de criacao    : 01/09/2016                                    *
* Data de modificaco : 00/00/0000                                    *
* -------------------------------------------------------------------*
* CALLBICGCSTAB(2) : chama o gradiente biconjugados estabilizados    *                   
* -------------------------------------------------------------------*
* Parametro de entrada                                               *
* -------------------------------------------------------------------*
* -------------------------------------------------------------------*
* Parametro de saida :                                               *
* -------------------------------------------------------------------*
* -------------------------------------------------------------------*
* OBS:                                                               *
**********************************************************************/
void callBicgStabl2(INT const nEq     ,INT const nEqNov
                   ,INT const nAd     ,INT const nAdR
                   ,INT *RESTRICT ia  ,INT *RESTRICT ja
                   ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                   ,DOUBLE *RESTRICT m,DOUBLE *RESTRICT b
                   ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT t
                   ,DOUBLE *RESTRICT v,DOUBLE *RESTRICT r
                   ,DOUBLE *RESTRICT u,DOUBLE *RESTRICT r0
                   ,DOUBLE *RESTRICT w,DOUBLE *RESTRICT s
                   ,DOUBLE *RESTRICT p,DOUBLE *RESTRICT h
                   ,DOUBLE *RESTRICT z
                   ,DOUBLE const tol  ,unsigned int maxIt
                   ,bool const newX   ,FILE* fSolvLog
                   ,bool const fLog   ,bool const fPrint
                   ,Interface *iNeq   ,BufferOmp *bOmp
                   ,void(*matVec)()   ,DOUBLE(*dot)())
{

/*... MPI*/
  if (mpiVar.nPrcs > 1)
  {
    printf("MpiBicgStabl2: Nao implementado\n");
    mpiStop();
    exit(EXIT_FAILURE);
  }
/*...................................................................*/

/*... */
  else {
/*... OpenMp*/
    if (ompVar.fSolver) 
    {
      pbicgstabl2Omp(nEq   ,nAd
                    ,ia    ,ja
                    ,a     ,ad
                    ,m     ,b
                    ,x     ,t
                    ,v     ,r
                    ,u     ,r0
                    ,w     ,s
                    ,p     ,h
                    ,z
                    ,tol   ,maxIt
                    ,newX  ,fSolvLog
                    ,NULL  ,fLog
                    ,false ,fPrint
                    ,bOmp 
                    ,matVec,dot);
    }
/*...................................................................*/

/*... sequencial*/
    else 
    {
      pbicgstabl2(nEq   ,nAd
               ,ia    ,ja
               ,a     ,ad
               ,m     ,b
               ,x     ,t
               ,v     ,r
               ,u     ,r0
               ,w     ,s
               ,p     ,h
               ,z
               ,tol   ,maxIt
               ,newX  ,fSolvLog
               ,NULL  ,fLog
               ,false ,fPrint
               ,matVec,dot);
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/**********************************************************************
* Data de criacao    : 04/09/2016                                    *
* Data de modificaco : 00/00/0000                                    *
* -------------------------------------------------------------------*
* CALLGMRES : chama o GMRES                                          *
* -------------------------------------------------------------------*
* Parametro de entrada                                               *
* -------------------------------------------------------------------*
* -------------------------------------------------------------------*
* Parametro de saida :                                               *
* -------------------------------------------------------------------*
* -------------------------------------------------------------------*
* OBS:                                                               *
**********************************************************************/
void callGmres(INT const nEq     ,INT const nEqNov
              ,INT const nAd     ,INT const nAdR
              ,INT *RESTRICT ia  ,INT *RESTRICT ja
              ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
              ,DOUBLE *RESTRICT m,DOUBLE *RESTRICT b
              ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT g
              ,DOUBLE *RESTRICT h,DOUBLE *RESTRICT y
              ,DOUBLE *RESTRICT c,DOUBLE *RESTRICT s
              ,DOUBLE *RESTRICT e,short const nKrylov
              ,DOUBLE const tol  ,unsigned int maxIt
              ,bool const newX   ,FILE* fSolvLog
              ,bool const fLog   ,bool const fPrint
              ,Interface *iNeq   ,BufferOmp *bOmp
              ,void(*matVec)()   ,DOUBLE(*dot)())
{

/*... MPI*/
  if (mpiVar.nPrcs > 1);
/*...................................................................*/

/*... */
  else {
/*... OpenMp*/
    if (ompVar.fSolver)
      gmresOmp(nEq  ,nAd
                   ,ia   ,ja
                   ,a    ,ad
                   ,m    ,b
                   ,x    ,g
                   ,h    ,y
                   ,c    ,s
                   ,e    ,nKrylov
                   ,tol  ,maxIt
                   ,newX ,fSolvLog
                   ,NULL ,fLog
                   ,fLog ,fPrint
                   ,bOmp
                   ,matVec,dot);
/*...................................................................*/

/*... sequencial*/
    else {
      gmres(nEq    ,nAd
           ,ia     ,ja
           ,a      ,ad
           ,m      ,b
           ,x      ,g
           ,h      ,y
           ,c      ,s
           ,e      ,nKrylov
           ,tol    ,maxIt
           ,newX   ,fSolvLog
           ,NULL   ,fLog
           ,fLog   ,fPrint
           ,matVec ,dot);
    }
/*...................................................................*/
  }
/*...................................................................*/

}
/*********************************************************************/


/**********************************************************************
* Data de criacao    : 18/09/2017                                    *
* Data de modificaco : 00/00/0000                                    *
* -------------------------------------------------------------------*
* CALLMINRES : chama o minres                                        *                   
* -------------------------------------------------------------------*
* Parametro de entrada                                               *
* -------------------------------------------------------------------*
* -------------------------------------------------------------------*
* Parametro de saida :                                               *
* -------------------------------------------------------------------*
* -------------------------------------------------------------------*
* OBS:                                                               *
**********************************************************************/
void callMinRes(INT const nEq     ,INT const nEqNov
               ,INT const nAd     ,INT const nAdR
               ,INT *RESTRICT ia  ,INT *RESTRICT ja
               , DOUBLE *RESTRICT al, DOUBLE *RESTRICT ad
	             , DOUBLE *RESTRICT m , DOUBLE *RESTRICT b
               , DOUBLE *RESTRICT x , DOUBLE *RESTRICT v0
               , DOUBLE *RESTRICT v , DOUBLE *RESTRICT w
               , DOUBLE *RESTRICT w0, DOUBLE *RESTRICT w00
               , DOUBLE *RESTRICT z , DOUBLE *RESTRICT z0
               , DOUBLE *RESTRICT p
               ,DOUBLE const tol  ,unsigned int maxIt
               ,bool const newX   ,FILE* fSolvLog
               ,bool const fLog   ,bool const fPrint
               ,Interface *iNeq   ,BufferOmp *bOmp
               ,void(*matVec)()   ,DOUBLE(*dot)())
{

/*... MPI*/
  if (mpiVar.nPrcs > 1);
/*  mpiPbicgstab(nEq, nEqNov
                 , nAd, nAdR
                 , ia, ja
                 , a, ad
                 , m, b
                 , x, t
                 , v, r
                 , p, z
                 , tol
                 , maxIt, newX
                 , fSolvLog, fLog
                 , fPrint, iNeq
                 , matVec, dot);  */
/*...................................................................*/

/*... */
  else {
/*... OpenMp*/
    if (ompVar.fSolver) {
/*    pbicgstabl2Omp(nEq   ,nAd
                    ,ia    ,ja
                    ,a     ,ad
                    ,m     ,b
                    ,x     ,t
                    ,v     ,r
                    ,u     ,r0
                    ,w     ,s
                    ,p     ,h
                    ,z
                    ,tol   ,maxIt
                    ,newX  ,fSolvLog
                    ,NULL  ,fLog
                    ,false ,fPrint
                    ,bOmp 
                    ,matVec,dot);*/
    }
/*...................................................................*/

/*... sequencial*/
    else {
      pminres(nEq   ,nAd
             ,ia    ,ja
             ,al    ,ad
             ,m     ,b
             ,x     ,v0
             ,v     ,w
             ,w0    ,w00
             ,z     ,z0
             ,p 
             ,tol   ,maxIt
             ,newX  ,fSolvLog
             ,NULL  ,fLog
             ,false ,fPrint
             ,matVec,dot);  
     }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/***********************************************************************
 * Data de criacao    : 30/04/2016                                    *
 * Data de modificaco : 00/00/0000                                    * 
 * ------------------------------------------------------------------ *   
 * CALL_PARDISO : chama o sover pardiso mkl                           *    
 * ------------------------------------------------------------------ * 
 * Parametros de entrada:                                             *
 * ------------------------------------------------------------------ * 
 * nEq      - numero de equacoes                                      *
 * mtype    - -2 simetrico indefinido                                 *
               2 simetrico positivo definido                          *
 * ia(*)    - ponteiro das colunas no formato CSR                     *
 * ja(*)    - ponteiro das colunas no formato CSR                     *
 * a(*)     - matriz A                                                *
 * b(neq)   - vetor de forcas                                         *
 * x(neq)   - nao definido                                            *   
 * z(neq)   - vetor auxiliar                                          *
 * r(neq)   - vetor auxiliar                                          *
 * ------------------------------------------------------------------ * 
 * Parametros de saida:                                               *
 * ------------------------------------------------------------------ *
 * x(neq)      - vetor solucao                                        *
  * ----------------------------------------------------------------- * 
 * OBS:                                                               *
 * ------------------------------------------------------------------ *
 * A matriz simetrica guarda apenas a parte superior no formato CSR   *
 * parada (ia,ja,a)                                                   *
 **********************************************************************/ 
void callMklPardiso(INT nEq           , INT mtype
                  , INT   *RESTRICT ia, INT   *RESTRICT ja
                  , DOUBLE *RESTRICT a, DOUBLE *RESTRICT b
                  , DOUBLE *RESTRICT x, DOUBLE *RESTRICT z 
                  , DOUBLE *RESTRICT r                        
                  , bool const fPrint)
{
  
#if _MKL_
  int nThreads = ompVar.nThreadsCell;
  DOUBLE norm,normR,xKx,time,ddum;
/*... variavel interna do mkl( 64 btis)*/
  void *pt[64];
/*...*/
  INT i,iparm[64],msglvl,error,idum;
/*...*/
  INT maxfct,mnum,phase,nrhs;
/*...*/
  error  = 0;
  maxfct = 1;
  mnum   = 1;
  nrhs   = 1;
/*...................................................................*/

/*...................................................................*/
  for (i = 0; i < 64; i++) {
    pt[i]    = 0;
    iparm[i] = 0;
  }
/*...................................................................*/

/*... simetrico indefinido*/
  if( mtype == -2){
    iparm[0]  = 1;/*no solver default*/
    iparm[1]  = 2;/*fill-in reordering from METIS*/
    iparm[6]  = 2;/*numbers of iterative refinement steps*/
    iparm[ 9] = 8;/*perturbe the pivot elements with 1E-08*/
    iparm[20] = 1;/*Pivoting for symmetric indefinite matrices.*/                
    iparm[23] = 0;/*Parallel factorization control.*/
    iparm[34] = 1; /*C-style indexing*/
  }
/*...................................................................*/

/*... simetrico definido positivo*/
   else if( mtype == 2){
     iparm[0]  = 1;/*no solver default*/
     iparm[1]  = 2;/*fill-in reordering from METIS*/
     iparm[6]  = 2;/*numbers of iterative refinement steps*/
//   iparm[16] =-1;
//   iparm[17] =-1;
     iparm[23] = 0; /*Parallel factorization control.*/
     iparm[34] = 1; /*C-style indexing*/
   }             
/*...................................................................*/

/*... simetrico definido positivo*/
   else if( mtype == 11){
     iparm[0]  = 1;/*no solver default*/
     iparm[1]  = 2;/*fill-in reordering from METIS*/
     iparm[6]  = 2;/*numbers of iterative refinement steps*/
//   iparm[16] =-1;
//   iparm[17] =-1;
     iparm[23] = 0; /*Parallel factorization control.*/
     iparm[34] = 1; /*C-style indexing*/
   }             
/*...................................................................*/

/*...*/
   time = getTimeC(); 
/*...................................................................*/
 
   for (i = 0; i < nEq; i++) 
     r[i] = b[i] - z[i];
/*...*/
   phase = 13;       
   msglvl = 1;
   mkl_set_num_threads(&nThreads);
   pardiso (pt   , &maxfct, &mnum, &mtype, &phase
          , &nEq , a      , ia   , ja
          , &idum, &nrhs  , iparm, &msglvl, b, x, &error);
   time = getTimeC() - time;
/*...................................................................*/
 
/*...*/
/* mem = ( max(iparm[14], iparm[15] + iparm[16] ) )/1024e0;*/
/*...................................................................*/
  
/*... produto:  x*Kx*/
   matVecCsrSym(nEq,ia,ja,a,a,x,z);
   xKx = dot(x,z,nEq);
/*...................................................................*/  
 
/*... norm-2 = || x ||*/
   norm = sqrt(dot(x,x,nEq));
/*....................................................................*/

/*... Termination and release of memory*/
   phase = -1;/*release internal memory*/
   msglvl = 0;
   pardiso(pt    , &maxfct, &mnum, &mtype , &phase
          , &nEq , &ddum, &idum  , &idum
          , &idum, &nrhs  , iparm, &msglvl
          , &ddum, &ddum  , &error);
/*....................................................................*/

/*... r = (b - Ax) (calculo do residuo explicito)*/
	for (i = 0; i < nEq; i++) 
    r[i] = b[i] - z[i];

  normR = dot(r, r, nEq);
/*...*/
  if(fPrint)
	  printf(" (PARDISO) solver:\n"
		  			"\tEquations      = %20d\n"
 			  		"\tx * Kx         = %20.2e\n"
				    "\t|| x ||        = %20.2e\n"
					  "\t|| b - Ax ||   = %20.2e\n"
	          "\tCPU time(s)    = %20.2lf\n" 
	          ,nEq,xKx,norm,normR,time);
/*....................................................................*/
#endif
}
/**********************************************************************/
