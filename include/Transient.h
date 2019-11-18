#ifndef _TRANSIENT_H_
  #define _TRASIENT_H_
/*...*/
  #include<string.h>
/*...*/
  #include<Mesh.h>
  #include<Define.h>
  #include<HccaStdBool.h>
/*...................................................................*/
  void  changeSchemeTemporal(Temporal *ddt);
  void setTransientScheme(char *word,Temporal *ddt);
/*...*/
  void cellTransient(DOUBLE *RESTRICT volume ,INT *RESTRICT id
                    ,DOUBLE *RESTRICT u0     ,DOUBLE *RESTRICT u
                    ,LevelTime density       ,DOUBLE *RESTRICT f
                    ,Temporal const ddt      ,INT const numel
                    ,short const ndf         ,bool const fAdd);
/*...................................................................*/

/*...*/
  void cellTransientSimple(DOUBLE *RESTRICT volume ,INT *RESTRICT id
                        ,DOUBLE *RESTRICT u0     ,DOUBLE *RESTRICT u
                        ,LevelTime density       ,DOUBLE *RESTRICT f
                        ,Temporal const ddt      ,INT const nEq
                        ,INT const numel         ,short const ndf
                        ,bool const fAdd);
/*...................................................................*/

/*...*/
  void cellTransientSimpleInc(DOUBLE *RESTRICT volume ,INT *RESTRICT id
                        ,DOUBLE *RESTRICT prop      ,short  *RESTRICT mat  
                        ,DOUBLE *RESTRICT u0        ,DOUBLE *RESTRICT u
                        ,DOUBLE *RESTRICT f
                        ,Temporal const ddt         ,INT const nEq
                        ,INT const numel            ,short const ndf
                        ,bool const fAdd);
/*...................................................................*/

/*...*/
  void cellTransientPrime(DOUBLE *RESTRICT volume 
                         ,DOUBLE *RESTRICT u0     ,DOUBLE *RESTRICT u
                         ,DOUBLE *RESTRICT density,DOUBLE *RESTRICT f
                         ,Temporal const ddt      
                         ,INT const numel         ,short const ndf
                         ,bool const fAdd); 
/*...................................................................*/

/*...*/
  void cellTransientEnergy(DOUBLE *RESTRICT volume ,INT *RESTRICT id
                          ,DOUBLE *RESTRICT u0     ,DOUBLE *RESTRICT u
                          ,LevelTime density       ,LevelTime sHeat 
                          ,DOUBLE *RESTRICT f
                          ,Temporal const ddt      ,INT const numel
                          ,bool const fAdd);
/*...................................................................*/

/*...*/
  void updateTime(Temporal *ddt, FileOpt *opt
                , Macros *mm   , short const myId );
/*...................................................................*/
#endif /*_TRANSIENT_H_*/
