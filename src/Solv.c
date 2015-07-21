#include<Solv.h>
/**********************************************************************
 * SOLVERC: resolucao do sistema linear(Ax=b)                         *
 * -------------------------------------------------------------------*
 * Parametro de entrada                                               *
 * -------------------------------------------------------------------*
 * m       -> vetor de memoria                                        *
 * neq     -> numero de equacoes                                      *
 * nad     -> numero de elemetos nao nulos fora da diagonal principal *
 * ia      -> vetor das linhas da estrutura de dados                  *
 * ja      -> vetor das colunas da estrutura de dados                 *
 * al      -> parte inferior da matriz a                              *
 * ad      -> diagonal principal da matriz a                          *
 * au      -> parte superior da matriz a                              *
 *  b      -> vetor de forcas                                         *
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
**********************************************************************/
void solverC(Memoria *m    ,INT neq   ,INT nad
            ,INT *ia       ,INT *ja   
            ,DOUBLE *al    ,DOUBLE *ad,DOUBLE *au
            ,DOUBLE *b     ,DOUBLE *x
            ,DOUBLE tol    ,unsigned int maxIt
            ,short const storage ,short const solver
            ,FILE* fSolvLog      ,bool const fLog
            ,bool const newX     ,bool const openMp   
            ,bool const unSym    ,bool const loopWise)
{
  DOUBLE *z=NULL,*r=NULL,*pc=NULL;
  void   (*matVecC)();
  DOUBLE (*dotC)();

  switch(solver){
/*gradientes conjugados com precondicionador diagonal*/
    case PCG:
/*... precondiconador diagonal*/
      HccaAlloc(DOUBLE,m,pc,neq,"pc",false);
      tm.precondDiag = getTimeC() - tm.precondDiag;
      preCondDiag(pc,ad,neq);
      tm.precondDiag = getTimeC() - tm.precondDiag;
/*...................................................................*/

/*... arranjos auxiliares do pcg*/
      HccaAlloc(DOUBLE,m,z,neq,"z",false);
      zero(z,neq,DOUBLEC);
      HccaAlloc(DOUBLE,m,r,neq,"r",false);
      zero(r,neq,DOUBLEC);
/*...................................................................*/
      
/*... estrutura de dados de armazenamentos da matriz esparcas*/
      switch(storage){
/*... armazenamento CSR(a)*/
        case CSR:
            matVecC = NULL;
        break;
/*...................................................................*/

/*... armazenamento CSRD(a,ad)*/
        case CSRD:
          if(unSym)
            matVecC = matVecCsrD;
          else
            matVecC = matVecCsrDSym;
        break;
/*...................................................................*/

/*... armazenamento CSRC(ad,au,al)*/
        case CSRC:
            matVecC = NULL;
        break;
/*...................................................................*/

/*...*/
        default:
          ERRO_OP(__FILE__,__func__,storage);
        break;
/*...................................................................*/
      }
/*...................................................................*/
      
/*...*/
      dotC    = dot;
/*...................................................................*/

/*... gradientes conjugados*/
      tm.pcg = getTimeC() - tm.pcg;
      pcg(neq     ,nad
         ,ia      ,ja
         ,al      ,ad   ,au
         ,pc      ,b    ,x
         ,z       ,r    ,tol
         ,maxIt   ,true 
         ,fSolvLog,fLog
         ,false 
         ,matVecC ,dotC);  
      tm.pcg = getTimeC() - tm.pcg;
/*...................................................................*/
      
/*... liberando arranjos auxiliares do pcg*/
      HccaDealloc(m,r,"r",false);
      HccaDealloc(m,z,"z",false);
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

  if(!strcmp(word,"PCG"));
   *solver = PCG;

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

