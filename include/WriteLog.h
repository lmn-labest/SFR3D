#ifndef _WRITELOG_H_
  #define _WRITELOG_H_
/*...*/  
  #include<Mesh.h>
  #include<Sisteq.h>
  #include<Solv.h>
  #include<HccaTime.h>
  #include<Define.h>
/*...*/
  #include<stdio.h>

  void writeLog(Mesh mesh       ,Scheme sc  
             ,Solv *SolvD1      ,SistEq *sistEqD1
             ,Solv *SolvT1      ,SistEq *sistEqT1
             ,Time t
             ,bool const fSolvD1,bool const fSolvT1
             ,char *nameIn      ,FILE *file);


#endif/*_WRITELOG_H_*/
