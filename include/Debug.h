#ifndef _DEBUG_H
   #define _DEBUG_H

/*...*/
   #include<stdio.h>
/*...................................................................*/

/*...*/
   #include<Define.h>
   #include<HccaStdBool.h>
/*...................................................................*/

/*...*/
   void testeGeom(double *cc
                 ,double *ksi   ,double *mksi 
                 ,double *eta   ,double *meta 
                 ,double *normal,double *volume
                 ,double *xm    ,double *xmcc   
                 ,double *mkm   ,double *dcca                 
                 ,INT numel     ,short ndm
                 ,short maxViz);
/*...................................................................*/

/*...*/
   void testeSist(INT *ia      ,INT *ja
              ,double *au   ,double *ad
              ,double *al   ,double *b
              ,INT const neq,bool const unsym);
/*...................................................................*/
#endif
