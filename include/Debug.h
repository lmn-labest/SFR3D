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
   void testeGeom(DOUBLE *cc
                 ,DOUBLE *ksi   ,DOUBLE *mksi 
                 ,DOUBLE *eta   ,DOUBLE *meta 
                 ,DOUBLE *normal,DOUBLE *volume
                 ,DOUBLE *xm    ,DOUBLE *xmcc   
                 ,DOUBLE *vSkew ,DOUBLE *mvSkew   
                 ,DOUBLE *dcca                 
                 ,INT numel     ,short ndm
                 ,short maxViz);
/*...................................................................*/

/*...*/
   void testeSist(INT *ia      ,INT *ja
                 ,DOUBLE *au   ,DOUBLE *ad
                 ,DOUBLE *al   ,DOUBLE *b
                 ,INT const neq,bool const unsym);
/*...................................................................*/
#endif
