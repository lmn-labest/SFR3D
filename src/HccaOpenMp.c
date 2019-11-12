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

#ifndef _OPENMP
 char str[] = "Opnemp enable in run-time but disable in compile-time ";

  if(omp)
  {
    ERRO_GERAL(__FILE__,__func__,__LINE__,str);
  }
#endif
}
/*********************************************************************/

/**********************************************************************
 * Data de criacao    : 28/01/2018                                    *
 * Data de modificaco : 05/08/2019                                    *
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
  
  readMacroV2(fileIn, word, false, true);
  fprintf(fileLogExc,"%-20s:\n","OpenMp:");
  nOmp = (short)atol(word);
  ompVar->flag = true;
  do {
    readMacroV2(fileIn, word, false, true);
/*... solver*/
    if (!strcmp(word, "solver")) {
      readMacroV2(fileIn, word, false, true);
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
   else if (!strcmp(word, "cell")) {
     readMacroV2(fileIn, word, false, true);
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
    else if (!strcmp(word, "update")) {
      readMacroV2(fileIn, word, false, true);
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
    else if (!strcmp(word, "grad")) {
      readMacroV2(fileIn, word, false, true);
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

/*... reaction*/
    else if (!strcmp(word, "reaction")) {
      readMacroV2(fileIn, word, false, true);
/*...*/
      ompVar->nThreadsReaction = (short)atol(word);
      ompVar->fReaction = true;
/*...................................................................*/

/*...*/       
      fprintf(fileLogExc,"%-20s: %d\n","Update nThreads"
                        , ompVar->nThreadsReaction);
/*...................................................................*/
      nOmp--;
    }
/*...................................................................*/


  } while (nOmp);

  openMpCheck(ompVar->flag);

}
/*********************************************************************/