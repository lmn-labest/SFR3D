#ifndef _LOADS_H_
  #define _LOADS_H_
/*...*/
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<Define.h>
  #include<Mesh.h>
  #include<Properties.h>
  #include<Turbulence.h> 
/*...................................................................*/

/*...*/
  void initLoads();
  void getLoads(DOUBLE *par, Loads *ld, DOUBLE *xm);
/*...................................................................*/

/*...*/
  void loadSenProd(DOUBLE *tA,DOUBLE *par,DOUBLE *xm);
/*...................................................................*/

/*... carga por elmento e condicoes pescritas por celula no 
      metodo simple*/
  void pLoadSimple(DOUBLE *RESTRICT sP, DOUBLE *RESTRICT p
          , DOUBLE *RESTRICT tA     , DOUBLE *RESTRICT xmcc               
          , DOUBLE *RESTRICT velC   , DOUBLE *RESTRICT gradVel
          , DOUBLE const presC      , DOUBLE *RESTRICT gradPresC
          , DOUBLE const viscosityC , DOUBLE const effViscosityC  
          , DOUBLE *RESTRICT xx
          , DOUBLE *RESTRICT sl     , DOUBLE *RESTRICT e        
          , DOUBLE *RESTRICT t      , DOUBLE *RESTRICT n  
          , DOUBLE const densityC   , DOUBLE *RESTRICT wallPar 
          , DOUBLE const fArea      , DOUBLE const dcca
          , Loads *ld               , short  const ndm 
          , INT const nEl
          , bool const fCalVel      , bool const fCalPres
          , bool const fWallModel   , short const wallType);

  void pLoadSimplePres(DOUBLE *RESTRICT sP, DOUBLE *RESTRICT p
          , DOUBLE *RESTRICT tA         , DOUBLE *RESTRICT ksi
          , DOUBLE const presC          , DOUBLE *RESTRICT gradPresC 
          , DOUBLE *RESTRICT sl         , DOUBLE *RESTRICT e        
          , DOUBLE *RESTRICT t          , DOUBLE *RESTRICT n      
          , DOUBLE const densityC       , DOUBLE *RESTRICT velC                                             
          , DOUBLE const fArea          , DOUBLE const dd
          , Loads *ld                   , short  const ndm 
          , bool const fCal);

  void pLoadEnergy(PropVarFluid *vProp
               , DOUBLE *RESTRICT sP   , DOUBLE *RESTRICT p
               , DOUBLE *RESTRICT tA     , DOUBLE *RESTRICT velC
               , DOUBLE const uC         , DOUBLE *RESTRICT n  
               , DOUBLE const thermCoef  , DOUBLE const densityC
               , DOUBLE const viscosityC , DOUBLE const sHeatC
               , DOUBLE const prT        , DOUBLE *RESTRICT xm                   
               , DOUBLE const fArea      , DOUBLE const dcca
               , Loads *ld               , Loads *ldVel 
               , DOUBLE *RESTRICT wallPar, short  const ndm        
               , bool const fCal         , bool const fTemp
               , bool const iKelvin      , bool const fSheat
               , bool const fWallModel   , short const wallType);

  void pLoadOneEqK(DOUBLE *RESTRICT sP   , DOUBLE *RESTRICT p
               , DOUBLE *RESTRICT tA     , DOUBLE *RESTRICT velC
               , DOUBLE const uC         , DOUBLE *RESTRICT n  
               , DOUBLE const densityC   , DOUBLE const viscosityC
               , DOUBLE const prT        , DOUBLE *RESTRICT xm                   
               , DOUBLE const fArea      , DOUBLE const dcca
               , Loads *ld               , Loads *ldVel
               , short const lFaceReK  
               , DOUBLE *RESTRICT wallPar, short  const ndm          
               , bool const fCal         , bool const fWallModel   
               , short const wallType);

  void pLoad(DOUBLE *RESTRICT sP  , DOUBLE *RESTRICT p
           , DOUBLE *RESTRICT tA  , DOUBLE *RESTRICT velC
           , DOUBLE *RESTRICT n
           , DOUBLE const coefDifC, DOUBLE const densityC
           , DOUBLE *RESTRICT xm  , DOUBLE const fArea
           , DOUBLE const dcca    , Loads *ld
           , short const ndm      , bool const fCal);

  void pLoadCombustion(PropVarFluid *vProp
               , DOUBLE *RESTRICT sP      , DOUBLE *RESTRICT p
               , DOUBLE *RESTRICT tA      , DOUBLE *RESTRICT velC
               , DOUBLE *RESTRICT uC      , DOUBLE *RESTRICT n  
               , DOUBLE *RESTRICT diffCoef, DOUBLE const densityC
               , DOUBLE const viscosityC 
               , DOUBLE const prT         , DOUBLE *RESTRICT xx                   
               , DOUBLE const fArea       , DOUBLE const dcca
               , Loads *ld                , Loads *ldVel 
               , DOUBLE *RESTRICT wallPar , short  const ndm          
               , bool const fCal          , bool const fWallModel  
               , short const nComb        , short const wallType);
/*...................................................................*/


#endif/*_LOADS_H_*/
