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
