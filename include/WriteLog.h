#ifndef _WRITELOG_H_
  #define _WRITELOG_H_
/*...*/
  #include<Mesh.h>
  #include<Sisteq.h>
  #include<Solv.h>
  #include<HccaTime.h>
  #include<Define.h>
/*...*/
  #include<stdio.h>

  void writeLog(Mesh *mesh          ,Scheme *sc
             ,Solv *solvD1          ,SistEq *sistEqD1
             ,Solv *solvT1          ,SistEq *sistEqT1
             ,Solv *solvVel         ,SistEq *sistEqVel
             ,Solv *solvPres        ,SistEq *sistEqPres
             ,Solv *solvEnergy      ,SistEq *sistEqEnergy
             ,Solv *solvComb        ,SistEq *sistEqComb
             ,Time *t
             ,bool const fSolvD1    ,bool const fSolvT1
             ,bool const fSolvVel   ,bool const fSolvPres
             ,bool const fEnergy    ,bool const fTurb
             ,bool const fCombustion
             ,Omp omp
             ,char *nameIn          ,FILE *file);

  void writeLogMeanTime(Mesh *mesh         ,Scheme *sc
             ,Solv *solvD1                 ,SistEq *sistEqD1
             ,Solv *solvT1                 ,SistEq *sistEqT1
             ,Solv *solvVel                ,SistEq *sistEqVel
             ,Solv *solvPres               ,SistEq *sistEqPres
             ,Solv *solvEnergy             ,SistEq *sistEqEnergy
             ,Solv *solvComb               ,SistEq *sistEqComb
             ,Time *t                      ,Omp *omp
             ,bool const fSolvD1           ,bool const fSolvT1
             ,bool const fSolvVel          ,bool const fSolvPres
             ,bool const fEnergy           ,bool const fTurb
             ,bool const fCombustion
             ,char *nameIn                 ,FILE *file);

#endif/*_WRITELOG_H_*/
