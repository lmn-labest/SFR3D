#ifndef _CELLLLIB_H_
  #define _CELLLIB_H_
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<Assbly.h>
  #include<Define.h>
  #include<HccaBlas.h>
  #include<HccaStdBool.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<ParallelMpi.h>
  #include<Loads.h>
  #include<Turbulence.h>
  #include<Properties.h>
/*...................................................................*/

/*... chamada da biblioteca de elementos (difusao)*/
  void cellLibDif(Loads *loads               ,Diffusion *diff
                 ,DiffModel *dModel          
                 ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT lprop
                 ,INT   *RESTRICT lViz       ,INT    *RESTRICT lId
                 ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mksi    
                 ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT fArea 
                 ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT volume 
                 ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
                 ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity 
                 ,DOUBLE *RESTRICT lCeofDiffD
                 ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
                 ,DOUBLE *RESTRICT lA        ,DOUBLE *RESTRICT lB
                 ,DOUBLE *RESTRICT lRcell    ,Temporal *ddt
                 ,short  *RESTRICT lFaceR    ,short  *RESTRICT lFaceL
                 ,DOUBLE *RESTRICT u0        ,DOUBLE *RESTRICT lGradU0 
                 ,short const nEn            ,short const nFace
                 ,short const ndm            ,short const lib 
                 ,INT const nel);
/*...................................................................*/

/*... chamada da biblioteca de elementos (transporte)*/
  void cellLibTrans(Loads *loads
                  , Advection *advT          , Diffusion *diffT
                  , TransModel *tModel
                  , short *RESTRICT lGeomType, DOUBLE *RESTRICT lprop
                  , INT   *RESTRICT lViz     , INT *RESTRICT lId
                  , DOUBLE *RESTRICT ksi     , DOUBLE *RESTRICT mKsi
                  , DOUBLE *RESTRICT eta     , DOUBLE *RESTRICT fArea
                  , DOUBLE *RESTRICT normal  , DOUBLE *RESTRICT volume
                  , DOUBLE *RESTRICT xm      , DOUBLE *RESTRICT xmcc
                  , DOUBLE *RESTRICT dcca
                  , DOUBLE *RESTRICT lDensity, DOUBLE *RESTRICT lCoefDiff
                  , DOUBLE *RESTRICT vSkew   , DOUBLE *RESTRICT mvSkew
                  , DOUBLE *RESTRICT lA      , DOUBLE *RESTRICT lB
                  , DOUBLE *RESTRICT lRcell  , Temporal *ddt
                  , short  *RESTRICT lFaceR  , short  *RESTRICT lFaceL
                  , DOUBLE *RESTRICT u0      , DOUBLE *RESTRICT gradU0
                  , DOUBLE *RESTRICT vel     , DOUBLE *RESTRICT cc
                  , short const nEn          , short  const nFace
                  , short const ndm          , short const lib
                  , INT const nel);
/*...................................................................*/

/*... chamada da biblioteca de elementos (Energia)*/
  void cellLibEnergy(Loads *loads                , Loads *loadsVel
                   , Advection  *adv             , Diffusion *diff
                   , Turbulence *tModel          , EnergyModel *model  
                   , Combustion *cModel          , PropVarFluid *vProp
                   , short *RESTRICT lGeomType   
                   , INT   *RESTRICT lViz        , INT *RESTRICT lId
                   , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
                   , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT fArea
                   , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume
                   , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
                   , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
                   , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
                   , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
                   , DOUBLE *RESTRICT lRcell     , Temporal *ddt
                   , short  *RESTRICT lFaceR     , short  *RESTRICT lFaceL
                   , short  *RESTRICT lFaceVelR  , short  *RESTRICT lFaceVelL    
                   , DOUBLE *RESTRICT u          , DOUBLE *RESTRICT gradU
                   , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
                   , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres 
                   , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lSheat
                   , DOUBLE *RESTRICT lDviscosity, DOUBLE *RESTRICT ltConductvity
                   , DOUBLE *RESTRICT lEnthalpyk , DOUBLE *RESTRICT lGradY 
                   , DOUBLE *RESTRICT diffY      , DOUBLE *RESTRICT yFrac 
                   , DOUBLE const lRateHeat  
                   , DOUBLE *RESTRICT dField     , DOUBLE *RESTRICT wallPar
                   , DOUBLE const underU
                   , short const nEn             , short  const nFace
                   , short const ndm             , short const lib
                   , INT const nel);
/*...................................................................*/

/*...*/
  void cellLibCombustion(Loads *lComb      , Loads *lVel
               , Advection  *adv           , Diffusion *diff
               , Turbulence *tModel        , Combustion *cModel
               , PropVarFluid *vProp       
               , short *RESTRICT lGeomType 
               , INT   *RESTRICT lViz      , INT *RESTRICT lId
               , DOUBLE *RESTRICT ksi      , DOUBLE *RESTRICT mKsi
               , DOUBLE *RESTRICT eta      , DOUBLE *RESTRICT fArea
               , DOUBLE *RESTRICT normal   , DOUBLE *RESTRICT volume
               , DOUBLE *RESTRICT xm       , DOUBLE *RESTRICT xmcc
               , DOUBLE *RESTRICT dcca     , DOUBLE *RESTRICT cc
               , DOUBLE *RESTRICT vSkew    , DOUBLE *RESTRICT mvSkew
               , DOUBLE *RESTRICT lA       , DOUBLE *RESTRICT lB
               , DOUBLE *RESTRICT lRcell   , Temporal *ddt
               , short  *RESTRICT lFaceR   , short  *RESTRICT lFaceL
               , short  *RESTRICT lFaceVelR, short  *RESTRICT lFaceVelL
               , DOUBLE *RESTRICT u        , DOUBLE *RESTRICT gradU
               , DOUBLE *RESTRICT Q        , DOUBLE *RESTRICT vel
               , DOUBLE *RESTRICT pres     , DOUBLE *RESTRICT gradPres
               , DOUBLE *RESTRICT lDensity , DOUBLE *RESTRICT lDiff
               , DOUBLE *RESTRICT lEddyVisc
               , DOUBLE *RESTRICT dField   , DOUBLE *RESTRICT wallPar
               , DOUBLE const underU
               , short const nEn           , short  const nFace
               , short const ndm           , short const lib
               , INT const nel);
/*...................................................................*/

/*.......................... PRIME ..................................*/
/*...*/
  void cellLibVelExp(Loads *loadsVel, Loads *loadsPres
        ,Advection advVel           ,Diffusion diffVel
        ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT lprop
        ,INT   *RESTRICT lViz
        ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mKsi
        ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT fArea
        ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT volume
        ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
        ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity
        ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
        ,DOUBLE *RESTRICT lB        ,Temporal const ddt
        ,short  *RESTRICT lFaceVelR ,short  *RESTRICT lFaceVelL
        ,short  *RESTRICT lFacePresR,short  *RESTRICT lFacePresL
        ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres
        ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT gradVel
        ,DOUBLE *RESTRICT dField    ,DOUBLE *RESTRICT cc
        ,DOUBLE *RESTRICT bT        ,DOUBLE const underU
        ,const bool sPressure       ,const bool fResidual
        ,const short nEn            ,short const nFace
        ,const short ndm            ,INT const nel);
/*...................................................................*/

/*.......................... SIMPLE .................................*/
/*... chamada da biblioteca de elementos (escoamento-vel)*/
  void cellLibSimpleVel(Loads *loadsVel     ,Loads *loadsPres 
          ,Advection *advVel          ,Diffusion *diffVel
          ,short const typeSimple 
          ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT lprop
          ,INT   *RESTRICT lViz       ,INT *RESTRICT lId  
          ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mKsi
          ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT fArea
          ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT volume
          ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
          ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity
          ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
          ,DOUBLE *RESTRICT lA        ,DOUBLE *RESTRICT lB
          ,DOUBLE *RESTRICT lRcell    ,Temporal *ddt
          ,short  *RESTRICT lFaceVelR ,short  *RESTRICT lFaceVelL
          ,short  *RESTRICT lFacePresR,short  *RESTRICT lFacePresL
          ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres 
          ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT gradVel 
          ,DOUBLE *RESTRICT dField    ,DOUBLE *RESTRICT cc 
          ,DOUBLE const underU        ,const bool sPressure
          ,short const nEn            ,short  const nFace     
          ,short const ndm            ,short const lib    
          ,INT const nel);
/*... chamada da biblioteca de elementos (escoamento-pres)*/
  void cellLibSimplePres(Loads *loadsVel     ,Loads *loadsPres
	           	 ,Diffusion *diffPres
               ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT lprop
               ,INT   *RESTRICT lViz       ,INT *RESTRICT lId  
               ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mKsi
               ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT fArea
               ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
               ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity
               ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
               ,DOUBLE *RESTRICT lA        ,DOUBLE *RESTRICT lB
               ,DOUBLE *RESTRICT lRcell    ,Temporal *ddt
               ,short  *RESTRICT lFaceVelR ,short  *RESTRICT lFaceVelL
               ,short  *RESTRICT lFacePresR,short  *RESTRICT lFacePresL
               ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres 
               ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT dField 
               ,short const nEn            ,short  const nFace     
               ,short const ndm            ,short const lib    
               ,INT const nel);
/*...................................................................*/

/*.......................... SIMPLE-LM ..............................*/
/*... chamada da biblioteca de elementos (escoamento-vel-Low Mach)*/
  void cellLibSimpleVelLm(Loads *loadsVel, Loads *loadsPres
            , Advection  *advVel         , Diffusion *diffVel 
            , Turbulence *tModel         , MomentumModel *ModelMomentum
            , short const typeSimple     , short *RESTRICT lGeomType  
            , INT   *RESTRICT lViz       , INT *RESTRICT lId
            , DOUBLE *RESTRICT ksi       , DOUBLE *RESTRICT mKsi
            , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT fArea
            , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume
            , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc
            , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT cc
            , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew
            , DOUBLE *RESTRICT lA        , DOUBLE *RESTRICT lB
            , DOUBLE *RESTRICT lRcell    , Temporal *ddt
            , short  *RESTRICT lFaceVelR , short  *RESTRICT lFaceVelL
            , short  *RESTRICT lFacePresR, short  *RESTRICT lFacePresL
            , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres 
            , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT gradVel
            , DOUBLE *RESTRICT lDensity  , DOUBLE *RESTRICT lViscosity
            , DOUBLE *RESTRICT dField    , DOUBLE *RESTRICT stressR
            , DOUBLE *RESTRICT wallPar   , DOUBLE const densityMed
            , DOUBLE const underU        , const bool sPressure
            , short const nEn            , short  const nFace
            , short const ndm            , short const lib
            , INT const nel);
/*... chamada da biblioteca de elementos (escoamento-pres-Low Mach)*/
  void cellLibSimplePresLm(Loads *loadsVel  , Loads *loadsPres
	             , Diffusion *diffPres        , MassEqModel *eMass    
               , short *RESTRICT lGeomType  
               , INT   *RESTRICT lViz       , INT *RESTRICT lId  
               , DOUBLE *RESTRICT ksi       , DOUBLE *RESTRICT mKsi
               , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT fArea
               , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume
               , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc
               , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT lDensity
               , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew
               , DOUBLE *RESTRICT lA        , DOUBLE *RESTRICT lB
               , DOUBLE *RESTRICT lRcell    , Temporal *ddt
               , short  *RESTRICT lFaceVelR , short  *RESTRICT lFaceVelL
               , short  *RESTRICT lFacePresR, short  *RESTRICT lFacePresL
               , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres 
               , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT dField 
               , DOUBLE *RESTRICT temp      , DOUBLE *RESTRICT wallPar
               , short const nEn            , short  const nFace     
               , short const ndm            , short const lib    
               , INT const nel);
/*...................................................................*/

/*... chamada da biblioteca de correcao nao ortogonal para correcao 
      de pressao*/
  void cellLibSimpleNonOrthPres(Diffusion diffPres
               ,short *RESTRICT lGeomType
               ,DOUBLE *RESTRICT lprop   ,INT   *RESTRICT lViz
               ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
               ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea
               ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
               ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
               ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
               ,DOUBLE *RESTRICT lB      
               ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres 
               ,DOUBLE *RESTRICT dField  ,DOUBLE *RESTRICT cc
               ,short const nEn          ,short  const nFace     
               ,short const ndm          ,INT const nel);
/*...................................................................*/

/*... */
  void cellLibTurbulence(Loads *lVel    , Turbulence tModel 
      , short *RESTRICT lGeomType     , DOUBLE *RESTRICT lprop
      , INT   *RESTRICT lViz          , DOUBLE *RESTRICT ksi
      , DOUBLE *RESTRICT mKsi         , DOUBLE *RESTRICT eta
      , DOUBLE *RESTRICT fArea        , DOUBLE *RESTRICT normal    
      , DOUBLE *RESTRICT volume       , DOUBLE *RESTRICT xm        
      , DOUBLE *RESTRICT xmcc         , DOUBLE *RESTRICT dcca    
      , DOUBLE *RESTRICT cc           , DOUBLE *RESTRICT vSkew   
      , DOUBLE *RESTRICT mvSkew       , short  *RESTRICT faceVelR 
      , short *RESTRICT faceVelL      , DOUBLE *RESTRICT vel    
      , DOUBLE *RESTRICT gradVel      , DOUBLE *RESTRICT lDensity  
      , DOUBLE const dViscosity       , DOUBLE *viscosity          
      , DOUBLE *RESTRICT wallPar      , DOUBLE const cDyn
      , short const nEn               , short  const nFace
      , short const ndm               , short const lib
      , INT const nel);
/*...................................................................*/

/*...*/
  void cellLibOneEqK(Loads *ldsK     , Loads *ldsVel
                   , Turbulence *tModel                       
                   , Advection  *adv             , Diffusion *diff              
                   , short *RESTRICT lGeomType   , DOUBLE *RESTRICT lprop
                   , INT   *RESTRICT lViz        , INT *RESTRICT lId
                   , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
                   , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT fArea
                   , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume
                   , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
                   , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
                   , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
                   , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
                   , DOUBLE *RESTRICT lRcell     , Temporal *ddt
                   , short  *RESTRICT lFaceReK   , short  *RESTRICT lFaceLdK 
                   , short  *RESTRICT lFaceReVel , short  *RESTRICT lFaceLdVel 
                   , DOUBLE *RESTRICT u          , DOUBLE *RESTRICT gradU
                   , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
                   , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres  
                   , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lDviscosity
                   , DOUBLE *RESTRICT dField     , DOUBLE *RESTRICT wallPar    
                   , DOUBLE *RESTRICT cDyn
                   , short const nEn             , short  const nFace
                   , short const ndm             , INT const nel);
/*...................................................................*/

/*... funcoes geometricas*/
  void cellGeom2D(DOUBLE *RESTRICT lx       ,short *RESTRICT lnFace
                 ,short  *RESTRICT lGeomType,DOUBLE *RESTRICT xc
                 ,DOUBLE *RESTRICT ksi      ,DOUBLE *RESTRICT mksi
                 ,DOUBLE *RESTRICT eta      ,DOUBLE *RESTRICT mEta
                 ,DOUBLE *RESTRICT normal   ,DOUBLE *RESTRICT volume
                 ,DOUBLE *RESTRICT xm       ,DOUBLE *RESTRICT xmcc
                 ,DOUBLE *RESTRICT dcca
                 ,DOUBLE *RESTRICT vSkew    ,DOUBLE *RESTRICT mvSkew
                 ,short  *RESTRICT sn       ,short const maxNo   
                 ,short const maxViz        ,short const ndm
                 ,INT const nel);

   void cellGeom3D(DOUBLE *RESTRICT lx      ,short *RESTRICT lGeomType
                 ,short *RESTRICT lnFac     ,short *RESTRICT lnEn
                 ,DOUBLE *RESTRICT xc
                 ,DOUBLE *RESTRICT ksi      ,DOUBLE *RESTRICT mksi
                 ,DOUBLE *RESTRICT eta      ,DOUBLE *RESTRICT fArea
                 ,DOUBLE *RESTRICT normal   ,DOUBLE *RESTRICT volume
                 ,DOUBLE *RESTRICT xm       ,DOUBLE *RESTRICT xmcc
                 ,DOUBLE *RESTRICT dcca     
                 ,DOUBLE *RESTRICT vSkew    ,DOUBLE *RESTRICT mvSkew
                 ,short  *RESTRICT sn
                 ,short const maxNo         ,short const maxViz
                 ,short const ndm           ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (difusao)*/
  void cellDif2D(Loads *loads
                ,short *RESTRICT lGeomType,DOUBLE *RESTRICT lprop
                ,INT   *RESTRICT lViz     ,INT *RESTRICT lId
                ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mksi
                ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT mEta
                ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
                ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
                ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
                ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
                ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
                ,DOUBLE *RESTRICT lRcell  ,Temporal const ddt
                ,short  *RESTRICT lFaceR  ,short *RESTRICT lFaceL
                ,DOUBLE *RESTRICT u0      ,DOUBLE *RESTRICT lGradU0
                ,short const nen          ,short const nFace
                ,short const ndm          ,INT const nel);
  
  void cellDif3D(Loads *loads               ,Diffusion *diff 
                ,DiffModel *dModel          
                ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT prop
                ,INT *RESTRICT lViz         ,INT *RESTRICT lId  
                ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mKsi
                ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT fArea
                ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT volume
                ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
                ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity 
                ,DOUBLE *RESTRICT lCeofDiffD
                ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
                ,DOUBLE *RESTRICT lA        ,DOUBLE *RESTRICT lB
                ,DOUBLE *RESTRICT lRcell    ,Temporal *ddt   
                ,short  *RESTRICT lFaceR    ,short  *RESTRICT lFaceL
                ,DOUBLE *RESTRICT u0        ,DOUBLE *RESTRICT gradU0
                ,const short nEn            ,short const nFace    
                ,const short ndm            ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (transporte)*/
  void cellTrans2D(Loads *loads           
                ,Advection advT           ,Diffusion diffT
                ,short *RESTRICT lGeomType,DOUBLE *RESTRICT lprop
                ,INT   *RESTRICT lViz     ,INT *RESTRICT lId
                ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mksi
                ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT mEta
                ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
                ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
                ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity   
                ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
                ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
                ,DOUBLE *RESTRICT lRcell  ,Temporal const ddt
                ,short  *RESTRICT lFaceR  ,short *RESTRICT lFaceL
                ,DOUBLE *RESTRICT u0      ,DOUBLE *RESTRICT lGradU0
                ,DOUBLE *RESTRICT vel     ,DOUBLE *RESTRICT cc
                ,short const nen          ,short const nFace
                ,short const ndm          ,INT const nel);
  
  void cellTrans3D(Loads *loads           
                ,Advection *advT          ,Diffusion *diffT
                ,TransModel *tModel
                ,short *RESTRICT lGeomType,DOUBLE *RESTRICT lprop
                ,INT   *RESTRICT lViz     ,INT *RESTRICT lId
                ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mksi
                ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT mEta
                ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
                ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
                ,DOUBLE *RESTRICT dcca    
                ,DOUBLE *RESTRICT lDensity,DOUBLE *RESTRICT lCoefDiff
                ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
                ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
                ,DOUBLE *RESTRICT lRcell  ,Temporal *ddt
                ,short  *RESTRICT lFaceR  ,short *RESTRICT lFaceL
                ,DOUBLE *RESTRICT u0      ,DOUBLE *RESTRICT lGradU0
                ,DOUBLE *RESTRICT vel     ,DOUBLE *RESTRICT cc
                ,short const nen          ,short const nFace    
                ,short const ndm          ,INT const nel);
/*...................................................................*/

/*.......................... PRIME .................................*/

/*... biblioteca de celulas (simple - vel - Explicito)*/
  void cellVelExp2D(Loads *loadsVel ,Loads *loadsPres
              ,Advection advVel           ,Diffusion diffVel
              ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT prop
              ,INT *RESTRICT lViz
              ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mKsi
              ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT mEta
              ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT area
              ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity
              ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT vC        ,Temporal const ddt
              ,short  *RESTRICT lFaceVelR ,short *RESTRICT lFaceVelL
              ,short  *RESTRICT lFacePresR,short *RESTRICT lFacePresL
              ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres
              ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT gradVel
              ,DOUBLE *RESTRICT dField    ,DOUBLE *RESTRICT cc
              ,DOUBLE *RESTRICT bT        ,DOUBLE const underU
              ,const bool sPressure       ,const bool fResidual
              ,const short nEn            ,short const nFace
              ,const short ndm            ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (simple - vel - Explicito)*/
  void cellVelExp3D(Loads *loadsVel      ,Loads *loadsPres
             ,Advection advVel           ,Diffusion diffVel
             ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT prop
             ,INT *RESTRICT lViz
             ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mKsi
             ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT mEta
             ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT area
             ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
             ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity
             ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
             ,DOUBLE *RESTRICT vC        ,Temporal const ddt
             ,short  *RESTRICT lFaceVelR ,short *RESTRICT lFaceVelL
             ,short  *RESTRICT lFacePresR,short *RESTRICT lFacePresL
             ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres
             ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT gradVel
             ,DOUBLE *RESTRICT dField    ,DOUBLE *RESTRICT cc
             ,DOUBLE *RESTRICT bT        ,DOUBLE const underU
             ,const bool sPressure       ,const bool fResidual
             ,const short nEn            ,short const nFace
             ,const short ndm            ,INT const nel);
/*...................................................................*/

/*.......................... ENERGIA ................................*/
  void cellEnergy2D(Loads *loads           , Loads *loadsVel 
            , Advection *adv               , Diffusion *diff
            , Turbulence *tModel           , EnergyModel *model 
            , PropVarFluid *vProp          , short *RESTRICT lGeomType    
            , INT *RESTRICT lViz           , INT *RESTRICT lId
            , DOUBLE *RESTRICT ksi         , DOUBLE *RESTRICT mKsi
            , DOUBLE *RESTRICT eta         , DOUBLE *RESTRICT mEta
            , DOUBLE *RESTRICT normal      , DOUBLE *RESTRICT volume
            , DOUBLE *RESTRICT xm          , DOUBLE *RESTRICT xmcc
            , DOUBLE *RESTRICT dcca        , DOUBLE *RESTRICT cc
            , DOUBLE *RESTRICT vSkew       , DOUBLE *RESTRICT mvSkew
            , DOUBLE *RESTRICT lA          , DOUBLE *RESTRICT lB
            , DOUBLE *RESTRICT lRcell      , Temporal const ddt
            , short  *RESTRICT lFaceR      , short *RESTRICT lFaceL
            , short  *RESTRICT lFaceVelR   , short *RESTRICT lFacevelL
            , DOUBLE *RESTRICT u0          , DOUBLE *RESTRICT gradU0           
            , DOUBLE *RESTRICT vel         , DOUBLE *RESTRICT gradVel  
            , DOUBLE *RESTRICT pres        , DOUBLE *RESTRICT gradPres         
            , DOUBLE *RESTRICT lDensity    , DOUBLE *RESTRICT lSheat
            , DOUBLE *RESTRICT lDviscosity , DOUBLE *RESTRICT ltConductvity
            , DOUBLE *RESTRICT dField      , DOUBLE *RESTRICT wallPar
            , DOUBLE const underU            
            , const short nEn              , short const nFace
            , const short ndm              , INT const nel);
  
  void cellEnergy3D(Loads *loads          , Loads *lVel
            , Advection *advT             , Diffusion *diffT
            , Turbulence *tModel          , EnergyModel *eModel
            , Combustion *cModel          , PropVarFluid *vProp
            , short *RESTRICT lGeomType   
            , INT *RESTRICT lViz          , INT *RESTRICT lId
            , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
            , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT fArea
            , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume
            , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
            , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
            , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
            , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
            , DOUBLE *RESTRICT lRcell     , Temporal *ddt
            , short  *RESTRICT lFaceR     , short *RESTRICT lFaceL
            , short  *RESTRICT lFaceVelR  , short *RESTRICT lFaceVelL
            , DOUBLE *RESTRICT u0         , DOUBLE *RESTRICT gradU0
            , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
            , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres  
            , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lSheat
            , DOUBLE *RESTRICT lViscosity , DOUBLE *RESTRICT lTconductivity
            , DOUBLE *RESTRICT enthalpyk  , DOUBLE *RESTRICT gradY 
            , DOUBLE *RESTRICT diffY      , DOUBLE *RESTRICT yFrac     
            , DOUBLE const rateHeat    
            , DOUBLE *RESTRICT dField     , DOUBLE *RESTRICT wallPar
            , DOUBLE const underU           
            , const short nEn             , short const nFace
            , const short ndm             , INT const nel);
/*...................................................................*/

/*........................ TURBULENCIA ..............................*/
void cellKinectTurb3D(Loads *ldsK         , Loads *ldsVel        
            , Turbulence *tModel
            , Advection *advT             , Diffusion *diffT
            , short *RESTRICT lGeomType   , DOUBLE *RESTRICT prop
            , INT *RESTRICT lViz          , INT *RESTRICT lId
            , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
            , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT fArea
            , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume
            , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
            , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
            , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
            , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
            , DOUBLE *RESTRICT lRcell     , Temporal *ddt
            , short  *RESTRICT lFaceReK   , short *RESTRICT lFaceLdK
            , short  *RESTRICT lFaceReVel , short *RESTRICT lFaceLdVel
            , DOUBLE *RESTRICT u0         , DOUBLE *RESTRICT gradU0
            , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
            , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres  
            , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lViscosity 
            , DOUBLE *RESTRICT dField     , DOUBLE *RESTRICT wallPar      
            , DOUBLE *RESTRICT cDyn
            , const short nEn             , short const nFace
            , const short ndm             , INT const nel);
/*...................................................................*/

/*.......................... SIMPLE .................................*/

/*... biblioteca de celulas (simple - vel)*/
  void cellSimpleVel2D(Loads *loadsVel    ,Loads *loadsPres 
              ,Advection  *advVel         ,Diffusion *diffVel
              ,short const typeSimple 
              ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT lprop
              ,INT   *RESTRICT lViz       ,INT *RESTRICT lId
              ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mksi
              ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT mEta
              ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT volume
              ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity
              ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT lA        ,DOUBLE *RESTRICT lB
              ,DOUBLE *RESTRICT lRcell    ,Temporal const ddt
              ,short  *RESTRICT lFaceVelR ,short  *RESTRICT lFaceVelL
              ,short  *RESTRICT lFacePresR,short  *RESTRICT lFacePresL
              ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres 
              ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT gradVel
              ,DOUBLE *RESTRICT dField    ,DOUBLE *RESTRICT cc
              ,DOUBLE const underU        ,bool const sPressure
              ,short const nen            ,short const nFace
              ,short const ndm            ,INT const nel);

/*... biblioteca de celulas (simple - vel)*/
  void cellSimpleVel3D(Loads *loadsVel    ,Loads *loadsPres
              ,Advection  *advVel         ,Diffusion  *diffVel					
          		,short const typeSimple 
              ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT lprop
              ,INT   *RESTRICT lViz       ,INT *RESTRICT lId
              ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mksi
              ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT mEta
              ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT volume
              ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity
              ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT lA        ,DOUBLE *RESTRICT lB
              ,DOUBLE *RESTRICT lRcell    ,Temporal *ddt
              ,short  *RESTRICT lFaceVelR ,short  *RESTRICT lFaceVelL
              ,short  *RESTRICT lFacePresR,short  *RESTRICT lFacePresL
              ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres 
              ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT gradVel 
              ,DOUBLE *RESTRICT dField    ,DOUBLE *RESTRICT cc 
              ,DOUBLE const underU        ,bool const sPressure
              ,short const nen            ,short const nFace
              ,short const ndm            ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (simple - pres)*/
  void cellSimplePres2D(Loads *loadsVel ,Loads *loadsPres 						
             	,Diffusion *diffPres        
              ,short *RESTRICT lGeomType,DOUBLE *RESTRICT prop
              ,INT *RESTRICT lViz       ,INT *RESTRICT lId  
              ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
              ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT mEta
              ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT area   
              ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
              ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
              ,DOUBLE *RESTRICT lRcell  
              ,short  *RESTRICT lFaceVelR ,short  *RESTRICT lFaceVelL
              ,short  *RESTRICT lFacePresR,short  *RESTRICT lFacePresL
              ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres 
              ,DOUBLE *RESTRICT vel     ,DOUBLE *RESTRICT dField  
              ,const short nEn          ,short const nFace    
              ,const short ndm          ,INT const nel);
  
/*... biblioteca de celulas (simple - pres)*/
  void cellSimplePres3D(Loads *loadsVel         ,Loads *loadsPres
                      ,Diffusion *diffPres
                      ,short *RESTRICT lGeomType,DOUBLE *RESTRICT prop
                      ,INT *RESTRICT lViz       ,INT *RESTRICT lId  
                      ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
                      ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT mEta
                      ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT area   
                      ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
                      ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
                      ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
                      ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
                      ,DOUBLE *RESTRICT lRcell  
                      ,short  *RESTRICT lFaceVelR ,short  *RESTRICT lFaceVelL
                      ,short  *RESTRICT lFacePresR,short  *RESTRICT lFacePresL
                      ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres 
                      ,DOUBLE *RESTRICT vel     ,DOUBLE *RESTRICT dField  
                      ,const short nEn          ,short const nFace    
                      ,const short ndm          ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (simple - pres-correcao non-ortogonal)*/
  void cellSimpleNonOrthPres2D(Diffusion diffPres
              ,short *RESTRICT lGeomType
              ,DOUBLE *RESTRICT prop    ,INT *RESTRICT lViz
              ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
              ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT mEta
              ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT area   
              ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
              ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT lB      
              ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres 
              ,DOUBLE *RESTRICT dField  ,DOUBLE *RESTRICT cc
              ,const short nEn          ,short const nFace
              ,const short ndm          ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (simple - pres-correcao non-ortogonal)*/
  void cellSimpleNonOrthPres3D(Diffusion diffPres
              ,short *RESTRICT lGeomType
              ,DOUBLE *RESTRICT prop    ,INT *RESTRICT lViz
              ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
              ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT mEta
              ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT area   
              ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
              ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
              ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
              ,DOUBLE *RESTRICT lB      
              ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres 
              ,DOUBLE *RESTRICT dField  ,DOUBLE *RESTRICT cc
              ,const short nEn          ,short const nFace
              ,const short ndm          ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (simple - vel - levemente compressivel)*/

/*... biblioteca de celulas 2D(simple - vel - low mach)*/
  void cellSimpleVel2DLm(Loads *loadsVel , Loads *loadsPres  
           , Advection advVel           , Diffusion diffVel
           , Turbulence tModel          , MomentumModel ModelMomentum
           , short const typeSimple     , short *RESTRICT lGeomType  
           , INT *RESTRICT lViz         , INT *RESTRICT lId 
           , DOUBLE *RESTRICT ksi       , DOUBLE *RESTRICT mKsi 
           , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT mEta 
           , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT area 
           , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc 
           , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT cc 
           , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew 
           , DOUBLE *RESTRICT lA        , DOUBLE *RESTRICT lB 
           , DOUBLE *RESTRICT lRcell    , Temporal const ddt 
           , short  *RESTRICT lFaceVelR , short *RESTRICT lFaceVelL 
           , short  *RESTRICT lFacePresR, short *RESTRICT lFacePresL 
           , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres 
           , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT gradVel 
           , DOUBLE *RESTRICT lDensity  , DOUBLE *RESTRICT lDviscosity 
           , DOUBLE *RESTRICT dField    , DOUBLE *RESTRICT stressR
           , DOUBLE *RESTRICT wallPar   , DOUBLE const densityMed
           , DOUBLE const underU        , const bool sPressure 
           , const short nEn            , short const nFace 
           , const short ndm            , INT const nel);

/*... biblioteca de celulas 3D(simple - pres - low mach)*/
  void cellSimplePres2DLm(Loads *loadsVel, Loads *loadsPres		
    			 , Diffusion *diffPres       	, MassEqModel *eMass
           , short *RESTRICT lGeomType 
           , INT *RESTRICT lViz         , INT *RESTRICT lId  
           , DOUBLE *RESTRICT ksi       , DOUBLE *RESTRICT mKsi
           , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT mEta
           , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT area   
           , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc
           , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT lDensity
           , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew
           , DOUBLE *RESTRICT lA        , DOUBLE *RESTRICT lB
           , DOUBLE *RESTRICT lRcell    , Temporal const ddt
           , short  *RESTRICT lFaceVelR , short *RESTRICT lFaceVelL
           , short  *RESTRICT lFacePresR, short *RESTRICT lFacePresL
           , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres 
           , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT dField  
           , DOUBLE *RESTRICT temp
           , const short nEn            , short const nFace    
           , const short ndm            , INT const nel);

/*... biblioteca de celulas 3D(simple - vel - low mach)*/
  void cellSimpleVel3DLm(Loads *lVel      , Loads *lPres 
            , Advection *advVel           , Diffusion *diffVel
            , Turbulence *tModel          , MomentumModel *ModelMomentum
            , short const typeSimple      , short *RESTRICT lGeomType   
            , INT *RESTRICT lViz          , INT *RESTRICT lId  
            , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
            , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT fArea
            , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume 
            , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
            , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
            , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
            , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
            , DOUBLE *RESTRICT lRcell     , Temporal *ddt
            , short  *RESTRICT lFaceVelR  , short *RESTRICT lFaceVelL
            , short  *RESTRICT lFacePresR , short *RESTRICT lFacePresL
            , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres 
            , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
            , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lViscosity 
            , DOUBLE *RESTRICT dField     , DOUBLE *RESTRICT stressR
            , DOUBLE *RESTRICT wallPar    , DOUBLE const densityMed
            , DOUBLE const underU         , const bool sPressure
            , const short nEn             , short const nFace    
            , const short ndm             , INT const nel);
/*... biblioteca de celulas 3D(simple - pres - low mach)*/
  void cellSimplePres3DLm(Loads *lVel      , Loads *lPres 
							, Diffusion *diffPres        , MassEqModel *eMass
              , short *RESTRICT lGeomType  
              , INT *RESTRICT lViz         , INT *RESTRICT lId  
              , DOUBLE *RESTRICT ksi       , DOUBLE *RESTRICT mKsi
              , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT fArea
              , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume
              , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc
              , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT lDensity
              , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew
              , DOUBLE *RESTRICT lA        , DOUBLE *RESTRICT lB
              , DOUBLE *RESTRICT lRcell    , Temporal *ddt 
              , short  *RESTRICT lFaceVelR , short *RESTRICT lFaceVelL
              , short  *RESTRICT lFacePresR, short *RESTRICT lFacePresL
              , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres
              , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT dField
              , DOUBLE *RESTRICT temp      , DOUBLE *RESTRICT wallPar
              , const short nEn            , short const nFace
              , const short ndm            , INT const nel);
/*...................................................................*/

/*...*/
  void cellCombustion3D(Loads *loads          , Loads *lVel
                  , Advection *advT           , Diffusion *diffT
                  , Turbulence *tModel        , Combustion *cModel
                  , PropVarFluid *vProp       , short *RESTRICT lGeomType 
                  , INT *RESTRICT lViz        , INT *RESTRICT lId
                  , DOUBLE *RESTRICT ksi      , DOUBLE *RESTRICT mKsi
                  , DOUBLE *RESTRICT eta      , DOUBLE *RESTRICT fArea
                  , DOUBLE *RESTRICT normal   , DOUBLE *RESTRICT volume
                  , DOUBLE *RESTRICT xm       , DOUBLE *RESTRICT xmcc
                  , DOUBLE *RESTRICT dcca     , DOUBLE *RESTRICT cc
                  , DOUBLE *RESTRICT vSkew    , DOUBLE *RESTRICT mvSkew
                  , DOUBLE *RESTRICT lA       , DOUBLE *RESTRICT lB
                  , DOUBLE *RESTRICT lRcell   , Temporal *ddt
                  , short  *RESTRICT lFaceR   , short *RESTRICT lFaceL
                  , short  *RESTRICT lFaceVelR, short *RESTRICT lFaceVelL
                  , DOUBLE *RESTRICT u0       , DOUBLE *RESTRICT gradU0
                  , DOUBLE *RESTRICT Q        , DOUBLE *RESTRICT vel
                  , DOUBLE *RESTRICT pres     , DOUBLE *RESTRICT gradPres
                  , DOUBLE *RESTRICT lDensity , DOUBLE *RESTRICT lDiff
                  , DOUBLE *RESTRICT lEddyVisc
                  , DOUBLE *RESTRICT dField   , DOUBLE *RESTRICT wallPar
                  , DOUBLE const underU
                  , const short nEn           , short const nFace
                  , const short ndm           , INT const nel);
/*...................................................................*/

/*... least square*/
  void leastSquareMatrix(DOUBLE *RESTRICT lKsi, DOUBLE *RESTRICT lmKsi
                 , DOUBLE *RESTRICT lLsquare, DOUBLE *RESTRICT lSquareR
                 , short const type         
                 , short const lnFace       , short const ndm);

  void leastSquare(Loads *loads
                 ,DOUBLE *RESTRICT lLsquare,INT *RESTRICT lViz 
                 ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc 
                 ,DOUBLE *RESTRICT lProp   ,DOUBLE *RESTRICT lDcca
                 ,DOUBLE *RESTRICT u       ,DOUBLE *RESTRICT gradU
                 ,short  *RESTRICT lFaceR  ,short *RESTRICT lFaceL
                 ,short const nFace        ,short const ndf
                 ,short const ndm          ,INT const nel);

  void leastSquareQR(Loads *loads
                 ,DOUBLE *RESTRICT lLs     ,DOUBLE *RESTRICT lLsR
                 ,DOUBLE *RESTRICT lProp   ,DOUBLE *RESTRICT lDcca
                 ,INT *RESTRICT lViz       
                 ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc  
                 ,DOUBLE *RESTRICT u       ,DOUBLE *RESTRICT gradU
                 ,short  *RESTRICT lFaceR  ,short *RESTRICT lFaceL
                 ,short const nFace        ,short const ndf
                 ,short const ndm);
/*... biblioteca de reconstrucao de gradiente*/
  void cellLibRcGrad(Loads *loads
           ,INT   *RESTRICT lViz,DOUBLE *RESTRICT lProp
           ,DOUBLE *RESTRICT leastSquare ,DOUBLE *RESTRICT leastSquareR
           ,DOUBLE *RESTRICT ksi         ,DOUBLE *RESTRICT mKsi
           ,DOUBLE *RESTRICT eta         ,DOUBLE *RESTRICT fArea
           ,DOUBLE *RESTRICT normal      ,DOUBLE *RESTRICT volume
           ,DOUBLE *RESTRICT vSkew       
           ,DOUBLE *RESTRICT xm          ,DOUBLE *RESTRICT xmcc 
           ,DOUBLE *RESTRICT lDcca  
           ,short  *RESTRICT lFaceR      ,short *RESTRICT lFaceL
           ,DOUBLE *RESTRICT u           ,DOUBLE *RESTRICT gradU 
           ,DOUBLE *RESTRICT nU          ,short const ty 
           ,short const nFace            ,short const ndm  
           ,short const lib              ,short const ndf  
           ,short *RESTRICT  isNodIN    ,INT const nel);

  void rcLeastSquare(INT *RESTRICT cellFace  , INT *RESTRICT fOwner
                   , DOUBLE *RESTRICT fModKsi, DOUBLE *RESTRICT fKsi
                   , DOUBLE *RESTRICT lSquare, DOUBLE *RESTRICT lSquareR
                   , short *RESTRICT nFace
                   , INT const numel         , short const maxViz
                   , short const type        , short const ndm);
/*...................................................................*/

/*...................................................................*/
  
/*... green-gauss*/
  void greenGaussCell(Loads *loads
               ,INT *RESTRICT lViz ,DOUBLE *RESTRICT mKsi
               ,DOUBLE *RESTRICT lProp   ,DOUBLE *RESTRICT lDcca
               ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT mfArea
               ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT mvSkew
               ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc    
               ,short  *RESTRICT lFaceR  ,short  *RESTRICT lFaceS
               ,DOUBLE *RESTRICT u       ,DOUBLE *RESTRICT gradU 
               ,short const nFace        ,short const ndm   
               ,short const ndf          ,INT const nel);
  void greenGaussNode(INT *RESTRICT lViz   ,DOUBLE *RESTRICT mEta
               ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT u       ,DOUBLE *RESTRICT gradU 
               ,short *RESTRICT isNod       
               ,short const nFace        ,short const ndm   
               ,short const ndf          ,short const ty);
/*...................................................................*/

/*....*/
  void meshQuality(MeshQuality *mq
                 , INT *RESTRICT cellFace    , INT *RESTRICT fOwner
                 , short  *RESTRICT nFace    , DOUBLE *RESTRICT volume
                 , DOUBLE *RESTRICT fKsi     , DOUBLE *RESTRICT fNormal
                 , DOUBLE *RESTRICT fModvSkew, DOUBLE *RESTRICT gDcca
                 , short const maxViz        , short const ndm
                 , INT const numel);
/*...................................................................*/

/*...*/
  bool openDomain(Loads *loadVel
                 , short  *RESTRICT faceL, short  *RESTRICT nFace
                 , INT const numel      , short const maxViz);

  void wallFluid(short *RESTRICT faceR ,INT *RESTRICT nelcon
              ,short *RESTRICT nFace     
              ,INT const nEl         ,short const maxViz);
  
  DOUBLE totalMass(DOUBLE *RESTRICT density  , DOUBLE *RESTRICT volume
                  ,INT const nEl) ;

  void massFluxOpenDomain(Loads *loadVel              , Temporal const ddt
                          , INT *RESTRICT cellFace      , INT *RESTRICT fOwner
                          , short  *RESTRICT faceVelLoad, short  *RESTRICT nFace
                          , DOUBLE *RESTRICT fArea      , DOUBLE *RESTRICT fNormal
                          , DOUBLE *RESTRICT fXm
                          , DOUBLE *RESTRICT density    , DOUBLE *RESTRICT vel
                          , DOUBLE *massInOut           , DOUBLE *deltaMass
                          , INT const numel             , short const ndm
                          , short const maxViz);
 
  void hPres(DOUBLE *RESTRICT pres0, DOUBLE *RESTRICT pres
         , DOUBLE *RESTRICT dFluid, DOUBLE *RESTRICT cc
         , DOUBLE *RESTRICT gravity, DOUBLE *RESTRICT xRef
         , INT const nEl           , short const ndm); 
/*...................................................................*/

/*... funcoes limitadoras de fluxo*/
  DOUBLE limitFaceBase(DOUBLE const rr,short const iCod);
  
  DOUBLE faceBaseTvd(short const nAresta    ,short const idCell
                 ,DOUBLE *RESTRICT u0
                 ,DOUBLE *RESTRICT gradUv,DOUBLE *RESTRICT gradUp
                 ,DOUBLE *RESTRICT lKsi  ,DOUBLE const lModKsi 
                 ,DOUBLE const cv
                 ,short const iCod       ,short const ndm);

  DOUBLE faceBaseTvdV1(DOUBLE const uC        ,DOUBLE const uV
                 ,DOUBLE *RESTRICT gradUc,DOUBLE *RESTRICT gradUv
                 ,DOUBLE *RESTRICT lKsi  ,DOUBLE const lModKsi 
                 ,DOUBLE const cv
                 ,short const iCod       ,short const ndm);
/*...................................................................*/

/*...*/
  DOUBLE interpolFaceVel(DOUBLE *RESTRICT velC,DOUBLE *RESTRICT velV
           ,DOUBLE const presC        ,DOUBLE const presV
           ,DOUBLE *RESTRICT gradPresC,DOUBLE *RESTRICT gradPresV
           ,DOUBLE *RESTRICT lNormal  ,DOUBLE *RESTRICT vKsi
           ,DOUBLE const mKsi         ,DOUBLE *RESTRICT dFieldF
           ,DOUBLE const alphaMenosUm ,DOUBLE const alpha
           ,short const ndm);
  /*...................................................................*/

/*... NVD-HR-CDC*/
  DOUBLE faceBaseNvd(DOUBLE const uP, DOUBLE const uV
                     , DOUBLE *RESTRICT gradUp, DOUBLE *RESTRICT gradUv
                     , DOUBLE *RESTRICT lKsi, DOUBLE const lModKsi
                     , DOUBLE const m
                     , short const iCod, short const ndm);
  void nvd(DOUBLE const phiTil,DOUBLE *p, short const iCod);
/*...................................................................*/

/*... parametro fisicos do escoamento*/
  void parameterCell(DOUBLE *RESTRICT vel, DOUBLE *RESTRICT prop
                    , DOUBLE *RESTRICT density, DOUBLE *RESTRICT volume
                    , short  *RESTRICT mat
                    , DOUBLE *cfl, DOUBLE *reynolds
                    , bool *fParameter, DOUBLE const dt
                    , INT const nEl, short const ndm);
/*...................................................................*/

/*... parametro fisicos do escoamento*/  
  void parameterCellLm(DOUBLE *RESTRICT vel , DOUBLE *RESTRICT prop
                , DOUBLE *RESTRICT density  , DOUBLE *RESTRICT sHeat
                , DOUBLE *RESTRICT tCond    , DOUBLE *RESTRICT viscosity
                , DOUBLE *RESTRICT volume   , short  *RESTRICT mat        
                , DOUBLE *cfl               , DOUBLE *reynolds
                , DOUBLE *peclet            , DOUBLE *mass     
                , bool *fParameter          , DOUBLE const dt
                , INT const nEl             , short const ndm);
/*...................................................................*/

/*... funcoes de correcao atrasada*/
  void difusionScheme(DOUBLE *RESTRICT n ,DOUBLE *RESTRICT ksi 
                    ,DOUBLE const fArea ,DOUBLE const lModKsi  
                    ,DOUBLE *RESTRICT e ,DOUBLE *RESTRICT t   
                    ,short const ndm    ,short const iCod);

  void difusionSchemeNew(DOUBLE *RESTRICT s , DOUBLE *RESTRICT ksi 
                      , DOUBLE *RESTRICT e , DOUBLE *RESTRICT t   
                      , short const ndm    , short const iCod);

  void difusionSchemeAnisotropic(DOUBLE *RESTRICT s,DOUBLE *RESTRICT ksi
                                ,DOUBLE *RESTRICT e,DOUBLE *RESTRICT t
                                ,short const ndm   ,short const iCod);

  void advectiveSchemeNdim(DOUBLE *RESTRICT uC ,DOUBLE *RESTRICT uV
                ,DOUBLE *RESTRICT gradUc     ,DOUBLE *RESTRICT gradUv
                ,DOUBLE *RESTRICT gradUComp  ,DOUBLE *RESTRICT vSkew
                ,DOUBLE *RESTRICT rC         ,DOUBLE *RESTRICT rV
                ,DOUBLE *RESTRICT ksi        ,DOUBLE const modKsi
                ,DOUBLE const wfn            ,DOUBLE *RESTRICT cvc 
                ,DOUBLE const alphaMenosUm   ,DOUBLE const alpha
                ,DOUBLE *RESTRICT parameters ,short const ndm
                ,short const ndf               
                ,short const iCod1           ,short const iCod2);

  void advectiveScheme(DOUBLE *RESTRICT velC      ,DOUBLE *RESTRICT velV
                     ,DOUBLE *RESTRICT gradVelC   ,DOUBLE *RESTRICT gradVelV
                     ,DOUBLE *RESTRICT gradVelComp,DOUBLE *RESTRICT vSkew
                     ,DOUBLE *RESTRICT rC         ,DOUBLE *RESTRICT rV
                     ,DOUBLE *RESTRICT ksi        ,DOUBLE const modKsi
                     ,DOUBLE const wfn            ,DOUBLE *RESTRICT cvc
                     ,DOUBLE const alphaMenosUm   ,DOUBLE const alpha
                     ,DOUBLE *RESTRICT parameters ,short const ndm
                     ,short const iCod1           ,short const iCod2);
   
  void advectiveSchemeScalar(DOUBLE const uC, DOUBLE const uV
              ,DOUBLE *RESTRICT gradUc    ,DOUBLE *RESTRICT gradUv
              ,DOUBLE *RESTRICT gradUcomp ,DOUBLE *RESTRICT vSkew
              ,DOUBLE *RESTRICT rC        ,DOUBLE *RESTRICT rV
              ,DOUBLE *RESTRICT ksi       ,DOUBLE const modKsi
              ,DOUBLE const wfn           ,DOUBLE *cvc
              ,DOUBLE const alphaMenosUm  ,DOUBLE const alpha          
              ,DOUBLE *RESTRICT parameters,short const ndm
              ,short const iCod1          ,short const iCod2);

  DOUBLE deferredCd(DOUBLE const velC,DOUBLE const velV
                   ,DOUBLE const wfn);
  
  DOUBLE deferredLust(DOUBLE const uC         ,DOUBLE const uV
                   ,DOUBLE *RESTRICT gradUc   ,DOUBLE *RESTRICT gradUv
                   ,DOUBLE *RESTRICT rC       ,DOUBLE *RESTRICT rV
                   ,DOUBLE const alphaMenosUm ,DOUBLE const alpha   
                   ,DOUBLE const beta             
                   ,DOUBLE const wfn          ,short const ndm);


  DOUBLE upwindLinearV1(DOUBLE const uC  ,DOUBLE const uV
                 ,DOUBLE *RESTRICT gradUc,DOUBLE *RESTRICT gradUv
                 ,DOUBLE *RESTRICT rC    ,DOUBLE *RESTRICT rV
                 ,DOUBLE const m         ,short const ndm);
/*...................................................................*/

/*...*/
  DOUBLE interpolFace(DOUBLE *RESTRICT lvSkew, DOUBLE *RESTRICT lXmcc
                  , DOUBLE volP            , DOUBLE volV
                  , DOUBLE lModKsi         , short ndm           
                  , short iCod    ) ;
/*...................................................................*/

/*... funcoes de apoio*/ 
  void gradFaceNull(DOUBLE *RESTRICT gradVelFace
                   ,DOUBLE *RESTRICT gradVelCell
                   ,DOUBLE *RESTRICT xmcc       ,short const ndm);

  short sn(short *s,short ty, INT nel);
  DOUBLE areaQuadCell(DOUBLE *RESTRICT eta,short ndm);
  DOUBLE areaTriaCell(DOUBLE *RESTRICT eta,short ndm);
  DOUBLE areaCell(DOUBLE *eta,short ty,short ndm,INT nel);
  DOUBLE volume3DGreenGauss(DOUBLE *RESTRICT xm,DOUBLE *RESTRICT normal
                         ,DOUBLE *RESTRICT fArea
                         ,short const nFace);
  DOUBLE volumeHexa(DOUBLE *RESTRICT x);
  DOUBLE volumeTetra(DOUBLE *RESTRICT x);
  void vectorKm2d(DOUBLE *RESTRICT x      ,DOUBLE *RESTRICT xc
                 ,DOUBLE *RESTRICT xm
                  ,DOUBLE *RESTRICT vSkew ,DOUBLE *RESTRICT mvSkew
                 ,short  *RESTRICT sn     ,short const nFace
                 ,short const maxViz      ,short const maxNo       
                 ,short const ndm         ,INT const nel);

  void vectorKm3d(DOUBLE *RESTRICT xc     ,DOUBLE *RESTRICT xm
                 ,DOUBLE *RESTRICT ksi    ,DOUBLE *RESTRICT normal
                 ,DOUBLE *RESTRICT vSkew ,DOUBLE *RESTRICT mvSkew
                 ,short const nFace      ,short const ndm
                 ,short const maxViz     ,INT nel);
  
  void setNvd(char *word,short *iCod);
  void setTvd(char *word,short *iCod);
  void  setDiffusionScheme(char *word,short *iCod);
  void  setAdvectionScheme(char *word, Advection *adv,FILE *fileIn);
  DOUBLE sizeCar(DOUBLE const volume,short const ndm);
  void vorticity(DOUBLE *RESTRICT w,DOUBLE *RESTRICT gradVel
                ,const short ndm);
  void stress(DOUBLE *RESTRICT s,DOUBLE *RESTRICT gradVel
           ,DOUBLE const nu   ,DOUBLE const lambda
           ,short const ndm);
  void stressEddyViscosity(DOUBLE *RESTRICT s, DOUBLE *RESTRICT gradVel
                         , DOUBLE const nut  , short const ndm);
  DOUBLE qCriterion(DOUBLE *RESTRICT gradVel, short const ndm);
  void velCorrectCombustion(DOUBLE *RESTRICT diff, DOUBLE *RESTRICT gradZ
                  , DOUBLE *RESTRICT velC      , short const ndm           
                  , short const ns    ); 
/*...................................................................*/


#endif/*_CELLLIB_H_*/