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

  void writeLog(Mesh mesh    
             ,Solv SolvD1 ,SistEq sistEqD1
             ,Time t
             ,char *nameIn,FILE *file);


#endif/*_WRITELOG_H_*/
