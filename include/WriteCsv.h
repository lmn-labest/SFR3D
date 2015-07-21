#ifndef _WRITECSV_
  #define _WRITECSV_
  #include<Define.h>
  #include<stdlib.h>

/*... Cell*/
  void writeCsvCell(DOUBLE *u      ,DOUBLE *gradU
                 ,DOUBLE *cc 
                 ,INT const numel,short const ndf
                 ,short const ndm,FILE *file);
  
/*... Node*/
  void writeCsvNode(DOUBLE *u      ,DOUBLE *gradU
                   ,DOUBLE *x  
                   ,INT const nNode,short const ndf
                   ,short const ndm,FILE *file);
  
#endif/*_WRITECSV_*/
