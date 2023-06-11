#ifndef _JACOBI_H_
  #define _JACOBI_H_

/*...*/
  #include<stdio.h>
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<Define.h>
  #include<HccaStdBool.h>
/*...................................................................*/

/*...*/
  void cyclic_jacobi(DOUBLE *RESTRICT a   , DOUBLE *RESTRICT x
                 , DOUBLE *RESTRICT w   , INT *RESTRICT ip
                 , INT const        n   , DOUBLE const tol
                 , short const maxciclos, bool const cflag);
/*...................................................................*/

/*... funcoe de auxiliares*/
  void indenty(DOUBLE *RESTRICT x, INT const n);
  void xax(DOUBLE *RESTRICT a, DOUBLE const c, DOUBLE const s
                 , INT const n       , INT const p   , INT const q);
  void xx(DOUBLE *RESTRICT x, DOUBLE const c, DOUBLE const s
        , INT const n       , INT const p   , INT const q);
  void symSchur2(DOUBLE const aqq, DOUBLE const app
               , DOUBLE const apq
               , DOUBLE *c       , DOUBLE *s);
/*...................................................................*/

/*...*/
  void ordena(DOUBLE *RESTRICT a, DOUBLE *RESTRICT x
            , DOUBLE *RESTRICT w, INT *RESTRICT ip
            , INT const n       , bool const fDes );
  DOUBLE offForbeniusNormSym(DOUBLE *RESTRICT a,INT const n);
  DOUBLE forbeniusNormSym(DOUBLE *RESTRICT a,INT const n);
/*...................................................................*/
#endif /*_JACOBI_H*/
