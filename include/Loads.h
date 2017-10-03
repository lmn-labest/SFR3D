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
  void getLoads(DOUBLE *par, Loads ld);
/*...................................................................*/

/*...*/
  void loadSenProd(DOUBLE *tA,DOUBLE *par,DOUBLE *xm);
/*...................................................................*/

/*... carga por elmento e condicoes pescritas por celula no 
      metodo simple*/
  void pLoadSimple(DOUBLE *RESTRICT sP  , DOUBLE *RESTRICT p
          , DOUBLE *RESTRICT tA          , DOUBLE *RESTRICT velC
          , DOUBLE *RESTRICT n             
          , DOUBLE *RESTRICT gradVel     , DOUBLE *RESTRICT xmcc
          , DOUBLE const viscosityC      , DOUBLE const effViscosityC 
          , DOUBLE const densityC          
          , DOUBLE const fArea           , DOUBLE const dcca
          , Loads ld                     , short const ndm
          , bool const fCal1             , bool const fCal2
          , bool const fWallModel        , short const wallType);

  void pLoadSimplePres(DOUBLE *RESTRICT sP  ,DOUBLE *RESTRICT p
          ,DOUBLE *RESTRICT tA
          ,DOUBLE const viscosityC,DOUBLE const densityC
          ,DOUBLE const wfn
          ,DOUBLE const fArea     ,DOUBLE const dcca 
          ,Loads ld               ,bool const fCal);

  void pLoadEnergy(DOUBLE *RESTRICT sP   , DOUBLE *RESTRICT p
               , DOUBLE *RESTRICT tA     , DOUBLE *RESTRICT velC
               , DOUBLE const uC         , DOUBLE *RESTRICT n  
               , DOUBLE const thermCoef  , DOUBLE const densityC
               , DOUBLE const viscosityC , DOUBLE const sHeatC
               , DOUBLE const prT        , DOUBLE *RESTRICT xm                   
               , DOUBLE const fArea      , DOUBLE const dcca
               , Loads ld                , Loads ldVel 
               , short  const ndm
               , bool const fCal         , bool const fTemp
               , bool const iKelvin      , bool const fSheat
               , bool const fWallModel   , short const wallType);
/*...................................................................*/


#endif/*_LOADS_H_*/
