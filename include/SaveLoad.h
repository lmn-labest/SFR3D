#ifndef _SAVELOAD_H_
  #define _SAVELOAD_H_
/*...*/ 
  #include<stdlib.h>
  #include<stdio.h>
 
/*...*/
  #include<Define.h>
  #include<Mesh.h>
  #include<File.h>
  #include<Properties.h>
  #include<Turbulence.h>
  #include<Residual.h>
/*...................................................................*/

/*...*/
  typedef struct{
    bool fLoad,fSave;
    short step[2];
  }Save;
/*...................................................................*/

/*...*/
  void wSave(PropVarFluid *prop           ,Turbulence *tModel
            ,ThermoDynamic *thDynamic
            ,Mesh *mesh              ,Temporal *ddt
            ,Mean *media
            ,char *preName           ,char *nameOut);
/*...................................................................*/

/*...*/
  void load(PropVarFluid *prop           ,Turbulence *tModel
           ,ThermoDynamic *thDynamic
           ,Mesh *mesh              ,Temporal *ddt
           ,Mean *media
           ,FILE *fileIn);
/*...................................................................*/

#endif/*_SAVELOAD_H_*/