#ifndef _SIMPLE_H
  #define _SMIPLE_H
  
  #include<File.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<CellLoop.h>
  #include<Mesh.h>
  #include<WriteVtk.h>
  #include<Sisteq.h>
  #include<Solv.h>
  #include<Transient.h>

  void simpleSolver(Memoria *m        
                   ,Loads *loadsVel   ,Loads *loadsPres 
                   ,Mesh *mesh0       ,Mesh *mesh       
                   ,SistEq *sistEqVel ,SistEq *sistEqPres
                   ,Solv *solvVel     ,Solv *solvPres 
                   ,Simple *sp
                   ,Scheme sc         ,PartMesh *pMesh 
                   ,FileOpt opt       ,char *preName  
                   ,char *nameOut     ,FILE *fileOut);

/*...*/
  void updateCellSimpleVel(DOUBLE *restrict w
                ,DOUBLE *restrict u1 ,DOUBLE *restrict u2
                ,INT *restrict id    ,INT const nEl
                ,short const ndm);

  void updateCellSimplePres(DOUBLE *restrict presC,DOUBLE *restrict xp
                           ,INT *restrict id,INT const nEl);

  void simpleUpdate(DOUBLE *restrict w     ,DOUBLE *restrict pressure
                 ,DOUBLE *restrict PresC ,DOUBLE *restrict GradPresC
                 ,DOUBLE *restrict dField         
                 ,INT const nEl          ,short const ndm
                 ,DOUBLE const alphaPres);
/*...................................................................*/

/*...*/
  void setSimpleScheme(char *word,Simple *sp,FILE *fileIn);
/*...................................................................*/

#endif/*_SIMPLE_H*/
