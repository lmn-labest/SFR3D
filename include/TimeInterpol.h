#ifndef _TI_H_
  #define _TI_H_
/*...*/
  #include<Erro.h>
  #include<HccaStdBool.h>
  #include<Define.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<ParallelMpi.h>

  void initTimeStruct(Memoria *m      ,TimeInterpol *ti
                   ,Mesh *mesh0       ,Mesh *mesh
                   ,Combustion *cModel,FileOpt *opt);
  void interPolTime(DOUBLE *RESTRICT ui ,DOUBLE *RESTRICT u
                 ,DOUBLE *RESTRICT u0   ,DOUBLE const ts
                 ,DOUBLE const t1       ,DOUBLE const t0
                 ,INT const nl          ,short const nc
                 ,bool const fTimePlot  ,bool const fStepPlot);

  void updateTimeStruct(Memoria *m        ,TimeInterpol *ti
                       ,Mesh *mesh
                       ,Combustion *cModel,FileOpt *opt);
#endif /*_FILE_H*/
