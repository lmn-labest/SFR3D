#ifndef _REACTION_H_
  #define _REACTION_H_
/*...*/
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<HccaStdBool.h>
  #include<Define.h>
  #include<StructDefine.h>
/*...................................................................*/

/*...*/
  DOUBLE arrhenius(DOUBLE const y1   ,DOUBLE const y2
                ,DOUBLE const e1     ,DOUBLE const e2
                ,DOUBLE const mW1    ,DOUBLE const mW2
                ,DOUBLE const t      ,DOUBLE const alpha
                ,DOUBLE const density,DOUBLE const tA    
                ,DOUBLE const coefA  ,bool const fKelvin);
/*...................................................................*/

/*...*/
  DOUBLE massActionMass(Reaction *reac     ,Prop *sHeatPol 
                     ,DOUBLE  *RESTRICT c
                     ,DOUBLE const T     ,unsigned short const nSp);
/*...................................................................*/

/*...*/
  void massRateReaction(Chemical *chem  ,DOUBLE *RESTRICT Q
                     ,DOUBLE *RESTRICT w);
/*...................................................................*/
#endif/*_REACTION_H_*/