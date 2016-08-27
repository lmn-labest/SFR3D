#ifndef _TRANSIENT_H_
  #define _TRASIENT_H_
/*...*/
  #include<string.h>
/*...*/
  #include<Mesh.h>
  #include<Define.h>
  #include<HccaStdBool.h>
/*...................................................................*/
  void setTransientScheme(char *word,short *type);
/*...*/
  void cellTransient(DOUBLE *restrict volume ,INT *restrict id
                    ,DOUBLE *restrict u0     ,DOUBLE *restrict u
                    ,DOUBLE *restrict density,DOUBLE *restrict f
                    ,Temporal const ddt      ,INT const numel
                    ,short const ndf         ,bool const fAdd);
/*...................................................................*/

/*...*/
  void cellTransientSimple(DOUBLE *restrict volume ,INT *restrict id
                          ,DOUBLE *restrict u0     ,DOUBLE *restrict u
                          ,DOUBLE *restrict density,DOUBLE *restrict f
                          ,Temporal const ddt      ,INT const nEq
                          ,INT const numel         ,short const ndf
                          ,bool const fAdd);
/*...................................................................*/

#endif /*_TRANSIENT_H_*/
