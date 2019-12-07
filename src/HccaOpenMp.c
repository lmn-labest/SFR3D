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
 * Data de modificaco : 07/12/2019                                    *
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
  int nthMax,nth;
  #if _OPENMP
    nthMax = omp_get_max_threads();
  #endif
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
      ompVar->fSolver        = true;
      nth = atol(word);
      if(nth<0)
        ompVar->nThreadsSolver = (short)nthMax;
      else
        ompVar->nThreadsSolver = (short) nth;
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
      ompVar->fCell = true;
      nth = atol(word);
      if(nth<0)
        ompVar->nThreadsCell = (short) nthMax;
      else
        ompVar->nThreadsCell = (short) nth;
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
      ompVar->fUpdate = true;
      nth = atol(word);
      if(nth<0)
        ompVar->fUpdate  = (short) nthMax;
      else
        ompVar->fUpdate  = (short) nth;
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
      ompVar->fGrad = true;
      nth = atol(word);
      if(nth<0)
        ompVar->fGrad = (short) nthMax;
      else
        ompVar->fGrad = (short) nth;
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
      ompVar->fReaction = true;
      nth = atol(word);
      if(nth<0)
        ompVar->nThreadsReaction = (short) nthMax;
      else
        ompVar->nThreadsReaction = (short) nth;
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