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
                  ,FILE *fileOut     );
/*...................................................................*/

/*...*/  
  void updateCellPrimePres(DOUBLE *RESTRICT presC
                          ,DOUBLE *RESTRICT xp
                          ,INT *RESTRICT id, INT const nEl);
/*...................................................................*/

/*...*/
  void primeUpdate(DOUBLE *RESTRICT w        ,DOUBLE *RESTRICT wUp
                  ,DOUBLE *RESTRICT pressure ,DOUBLE *RESTRICT presC 
                  ,DOUBLE *RESTRICT gradPresC,DOUBLE *RESTRICT dField
                  ,INT const nEl             ,short const ndm
                  ,DOUBLE const alphaPres);

  void residualPrime(DOUBLE *RESTRICT vel  ,DOUBLE *RESTRICT adVel
                 ,DOUBLE *RESTRICT rCellVel,DOUBLE *RESTRICT rCellMass
                 ,DOUBLE *RESTRICT rU      ,DOUBLE *RESTRICT rMass
                 ,INT const nEl            ,short const ndm
                 ,short const iCod );

  void setPrimeScheme(char *word, Prime *sp, FILE *fileIn);
/*...................................................................*/

#endif/*_PRIME_H*/
