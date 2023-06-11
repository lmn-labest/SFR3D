#ifndef _GAUSS_H_
  #define _GAUSS_H_
/*...*/
  #include<math.h>
  #include<stdlib.h>
  #include<stdio.h>
/*...................................................................*/

/*... */
  #include<Define.h>
  #include<HccaBlas.h>
  #include<HccaStdBool.h>
/*...................................................................*/


/*...*/
  void printFat(DOUBLE *a,DOUBLE *d,INT const n,short code);
/*...................................................................*/

/*...*/
  void initIterImprovmentData(LDOUBLE *RESTRICT ar,LDOUBLE *RESTRICT br
                             ,DOUBLE  *RESTRICT a ,DOUBLE *RESTRICT b
                             ,INT const nEq);
  void iterImprovement(LDOUBLE *RESTRICT ar ,LDOUBLE *RESTRICT br
                      ,LDOUBLE *RESTRICT itr,LDOUBLE *RESTRICT xr
                      ,DOUBLE  *RESTRICT x
                      ,DOUBLE  *RESTRICT a  ,DOUBLE *RESTRICT b
                      ,DOUBLE  *RESTRICT w
                      ,INT *RESTRICT p      ,INT const nEq);
/*...................................................................*/

/*... decomposicoes*/
  void fatLDMt(DOUBLE *RESTRICT a,DOUBLE *RESTRICT d,DOUBLE *RESTRICT r
              ,DOUBLE *RESTRICT w,INT const nEq);
  void fatLDLt(DOUBLE *RESTRICT a,DOUBLE *RESTRICT d,DOUBLE *RESTRICT r
              ,INT const nEq);
  void  fatGGt(DOUBLE *RESTRICT a,INT const nEq);
  void   fatLU(DOUBLE *RESTRICT a,INT const nEq,short const code);
  void fatLUpp(DOUBLE *RESTRICT a,INT    *RESTRICT p
              ,INT const nEq     ,short const code);
/*...................................................................*/

/*... decomposicoes imcomplretas*/
  void  ic0   (DOUBLE *RESTRICT a,INT const nEq);

/*... solver*/
  void solvTri(DOUBLE *RESTRICT a,DOUBLE *RESTRICT d
              ,DOUBLE *RESTRICT b,INT const nEq
              ,short const code);
  void solvCholesky(DOUBLE *RESTRICT a,DOUBLE *RESTRICT b,INT const nEq);
  void solvLU(DOUBLE *RESTRICT a      ,DOUBLE *RESTRICT b,INT const nEq);
  void solvLUpp(DOUBLE *RESTRICT a    ,DOUBLE *RESTRICT b
               ,DOUBLE *RESTRICT r    ,INT *RESTRICT p
               ,INT const nEq);
/*...................................................................*/

/*...*/
  void   diagScaling(DOUBLE *RESTRICT a,DOUBLE *RESTRICT scaling
                    ,INT const nEq);
  void scalingSystem(DOUBLE *RESTRICT a      ,DOUBLE *RESTRICT b
                    ,DOUBLE *RESTRICT scaling
                    ,INT const nEq);
/*...................................................................*/

/*...*/
  void solverD(DOUBLE *RESTRICT a,DOUBLE *RESTRICT b,INT const nEq
              ,short const code   ,bool fscaling);
/*...................................................................*/


#endif /*GAUSS_H*/
