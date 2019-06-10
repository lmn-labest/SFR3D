#include<OpenMp.h>

/**********************************************************************
 * Data de criacao    : 28/08/2016                                    *
 * Data de modificaco : 00/00/0000                                    *
 * -------------------------------------------------------------------*
 * PMATRIXSOLVEROMP: balanciamento do treabalho no openmp no solver   *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 * m  -> vetor de memoria                                             *
 * eq -> sistema de equacao                                           *
 * s1 -> nome do vetor                                                *
 * s2 -> nome do vetor                                                *
 * s3 -> nome do vetor                                                *
 * s4 -> nome do vetor                                                *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 * ------------------------------------------------------------------ *
*********************************************************************/
void pMatrixSolverOmp(Memoria *m,SistEq *eq
                     ,char *s1  ,char *s2
                     ,char *s3  ,char *s4){
#if _OPENMP
  short nth = ompVar.nThreadsSolver;

  HccaAlloc(INT, m, eq->omp.thBegin , nth, s1, false);
  HccaAlloc(INT, m, eq->omp.thEnd   , nth, s2, false);
  HccaAlloc(INT, m, eq->omp.thSize  , nth, s3, false);
  HccaAlloc(INT, m, eq->omp.thHeight, nth, s4, false);

  partitionCsrByNonzeros(eq->ia          ,eq->ja
                        ,eq->neq
                        ,nth             ,eq->omp.thBegin
                        ,eq->omp.thEnd   ,eq->omp.thSize
                        ,eq->omp.thHeight,eq->storage);
#endif
}
/*********************************************************************/

/**********************************************************************
 * Data de criacao    : 04/10/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 * -------------------------------------------------------------------*
 * openMpCheck : verifica se o codigo foi compilado com openmp        *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 * omp -> run-time openmp                                             *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 * ------------------------------------------------------------------ *
*********************************************************************/
void openMpCheck(bool omp) {

 char str[] = "Opnemp enable in run-time but disable in compile-time ";

  if(omp){
#ifndef _OPENMP
    ERRO_GERAL(__FILE__,__func__,__LINE__,str);
#endif
  }

}
/*********************************************************************/

/**********************************************************************
 * Data de criacao    : 28/01/2018                                    *
 * Data de modificaco : 00/00/0000                                    *
 * -------------------------------------------------------------------*
 * openMpSet: configuracao do openMP                                  *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 * omp -> run-time openmp                                             *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 * ------------------------------------------------------------------ *
*********************************************************************/
void openMpSet(FILE *fileIn, Omp *ompVar){
  char word[WORD_SIZE];
  unsigned short nOmp;
  
  readMacro(fileIn, word, false);
  fprintf(fileLogExc,"%-20s:\n","OpenMp:");
  nOmp = (short)atol(word);
  ompVar->flag = true;
  do {
    readMacro(fileIn, word, false);
/*... solver*/
    if (!strcmp(word, "solver") || !strcmp(word, "Solver")) {
      readMacro(fileIn, word, false);
/*... codigo da da funcao limitadora de fluxo*/
      ompVar->nThreadsSolver = (short)atol(word);
      ompVar->fSolver        = true;
/*...................................................................*/

/*...*/    
      fprintf(fileLogExc,"%-20s: %d\n","Solver nThreads"
                        ,ompVar->nThreadsSolver);
/*...................................................................*/
      nOmp--;
   }
/*...................................................................*/

/*... cell*/
   else if (!strcmp(word, "Cell") || !strcmp(word, "cell")) {
     readMacro(fileIn, word, false);
/*...*/
     ompVar->nThreadsCell = (short)atol(word);
     ompVar->fCell = true;
/*...................................................................*/

/*...*/       
     fprintf(fileLogExc,"%-20s: %d\n","Cell nThreads"
                       , ompVar->nThreadsCell);
/*...................................................................*/
     nOmp--;
   }
/*...................................................................*/

/*... update*/
    else if (!strcmp(word, "Update") || !strcmp(word, "update")) {
      readMacro(fileIn, word, false);
/*...*/
      ompVar->nThreadsUpdate = (short)atol(word);
      ompVar->fUpdate = true;
/*...................................................................*/

/*...*/       
      fprintf(fileLogExc,"%-20s: %d\n","Update nThreads"
                        , ompVar->nThreadsUpdate);
/*...................................................................*/
      nOmp--;
    }
/*...................................................................*/

/*... grad*/
    else if (!strcmp(word, "Grad") || !strcmp(word, "grad")) {
      readMacro(fileIn, word, false);
/*...*/
      ompVar->nThreadsGrad = (short)atol(word);
      ompVar->fGrad = true;
/*...................................................................*/

/*...*/       
      fprintf(fileLogExc,"%-20s: %d\n","Update nThreads"
                        , ompVar->nThreadsGrad);
/*...................................................................*/
      nOmp--;
    }
/*...................................................................*/

  } while (nOmp);

  openMpCheck(ompVar->flag);

}
/*********************************************************************/