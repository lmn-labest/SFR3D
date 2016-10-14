#ifndef _PRIME_H
  #define _PRIME_H

  #include<File.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<CellLoop.h>
  #include<Mesh.h>
  #include<WriteVtk.h>
  #include<Sisteq.h>
  #include<Solv.h>
  #include<Transient.h>

/*...*/
  void primeSolver(Memoria *m
                  ,Loads *loadsVel   ,Loads *loadsPres
                  ,Mesh *mesh0       ,Mesh *mesh
                  ,SistEq *sistEqPres,Solv *solvPres
                  ,Prime *pr         ,Scheme sc    
                  ,PartMesh *pMesh   ,FileOpt opt      
                  ,char *preName     ,char *nameOut   
                  ,FILE *fileOut     ,short const ndf);
/*...................................................................*/

/*...*/  
  void updateCellPrimePres(DOUBLE *restrict presC
                          ,DOUBLE *restrict xp
                          ,INT *restrict id, INT const nEl);
/*...................................................................*/

/*...*/
  void primeUpdate(DOUBLE *restrict w        ,DOUBLE *restrict wUp
                  ,DOUBLE *restrict pressure ,DOUBLE *restrict presC 
                  ,DOUBLE *restrict gradPresC,DOUBLE *restrict dField
                  ,INT const nEl             ,short const ndm
                  ,DOUBLE const alphaPres);

  void residualPrime(DOUBLE *restrict vel  ,DOUBLE *restrict adVel
                 ,DOUBLE *restrict rCellVel,DOUBLE *restrict rCellMass
                 ,DOUBLE *restrict rU      ,DOUBLE *restrict rMass
                 ,INT const nEl            ,short const ndm
                 ,short const iCod );

  void setPrimeScheme(char *word, Prime *sp, FILE *fileIn);
/*...................................................................*/

#endif/*_PRIME_H*/
