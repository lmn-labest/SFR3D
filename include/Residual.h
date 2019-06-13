#ifndef _RESIDUAL_H_
  #define _RESIDUAL_H_
/*...*/
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<Define.h>
  #include<Erro.h>
  #include<HccaBlas.h>
/*...................................................................*/

 void residualCombustion(DOUBLE *RESTRICT vel ,DOUBLE *RESTRICT energy
            ,DOUBLE *RESTRICT zComb         
            ,DOUBLE *RESTRICT rCellVel    ,DOUBLE *RESTRICT rCellMass
            ,DOUBLE *RESTRICT rCellEnergy ,DOUBLE *RESTRICT rCellComb      
            ,DOUBLE *RESTRICT adVel       ,DOUBLE *RESTRICT adEnergy 
            ,DOUBLE *RESTRICT adComb   
            ,DOUBLE *RESTRICT rU          ,DOUBLE *rMass
            ,DOUBLE *rEnergy              ,DOUBLE *rComb    
            ,INT  *RESTRICT idVel         ,INT  *RESTRICT idEnergy
            ,INT  *RESTRICT idComb    
            ,INT const nEl                ,INT const nEqVel
            ,INT const nEqComb
            ,short const ndm              ,short const nComb
            ,short *iCod);

  void residualType(DOUBLE *RESTRICT u      ,DOUBLE *RESTRICT rCellU
                 ,DOUBLE *RESTRICT adU    ,DOUBLE *RESTRICT rU
                 ,INT  *RESTRICT idU      ,INT const nEl   
                 ,INT const nEqU          ,short const ndf 
                 ,short iCod);

  void rScaled(DOUBLE *RESTRICT u    ,DOUBLE *RESTRICT ad
              ,DOUBLE *RESTRICT rCell,DOUBLE *RESTRICT rU
              ,INT  *RESTRICT id
              ,INT const nEl         ,INT const nEq
              ,short const ndf); 

  void rScaledSum(DOUBLE *RESTRICT u    ,DOUBLE *RESTRICT ad
               ,DOUBLE *RESTRICT rCell,DOUBLE *RESTRICT rU
               ,INT  *RESTRICT id
               ,INT const nEl         ,INT const nEq
               ,short const ndf);

  void rScaledSumMax(DOUBLE *RESTRICT u    ,DOUBLE *RESTRICT ad
                  ,DOUBLE *RESTRICT rCell,DOUBLE *RESTRICT rU
                  ,INT  *RESTRICT id
                  ,INT const nEl         ,INT const nEq
                  ,short const ndf); 

  DOUBLE residualMass(DOUBLE *RESTRICT rCellMass,INT const nEl
                  , short const iCod);
#endif/*_RESIDUAL_H_*/

