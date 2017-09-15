#ifndef _TURBULENCE_H_
  #define  _TURBULENCE_H_
/*...*/
  #include<CellLoop.h>
  #include<Define.h>
  #include<Mesh.h>
/*...................................................................*/

/*...*/
    void turbulence(Turbulence tModel     
     , INT    *RESTRICT el                  , INT    *RESTRICT nelcon 
     , short  *RESTRICT nen                 , short  *RESTRICT nFace
     , short  *RESTRICT geomType            , DOUBLE *RESTRICT prop
     , short  *RESTRICT calType             , short  *RESTRICT mat     
     , DOUBLE *RESTRICT cc                      
     , DOUBLE *RESTRICT gKsi                , DOUBLE *RESTRICT gmKsi 
     , DOUBLE *RESTRICT gEta                , DOUBLE *RESTRICT gfArea 
     , DOUBLE *RESTRICT gNormal             , DOUBLE *RESTRICT gVolume
     , DOUBLE *RESTRICT gXm                 , DOUBLE *RESTRICT gXmcc
     , DOUBLE *RESTRICT gvSkew              , DOUBLE *RESTRICT gmvSkew 
     , DOUBLE *RESTRICT gDcca                 
     , DOUBLE *RESTRICT gradVel             , DOUBLE *RESTRICT density
     , DOUBLE *RESTRICT eddyViscosity     
     , short const maxNo                    , short const maxViz
     , short const ndm                      , INT const numel
     , short const ndf);        
/*...................................................................*/

/*...*/
      void cellLes(Turbulence tModel,           
            short *RESTRICT lGeomType  , DOUBLE *RESTRICT prop,
            INT *RESTRICT lViz         , DOUBLE *RESTRICT ksi,
            DOUBLE *RESTRICT mKsi,
            DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT mEta,
            DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT area,
            DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc,
            DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT cc,
            DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew,
            DOUBLE *RESTRICT gradVel   , DOUBLE *RESTRICT lDensity,  
            DOUBLE *viscosity,
            const short nEn            , short const nFace,
            const short ndm            , INT const nel);
/*...................................................................*/

#endif /*_TURBULENCE_H_*/