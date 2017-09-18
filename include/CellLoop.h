#ifndef _CELLLOOP_H_
  #define _CELLLOOP_H_
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
/*...................................................................*/

/*...*/
  void convTempForKelvin(DOUBLE *RESTRICT uC,INT const n
                        ,bool const fKelvin);
/*...................................................................*/

/* ... calculo das propriedade geometicas*/
  void pGeomForm(DOUBLE *RESTRICT x    ,INT    *RESTRICT el
              ,INT    *RESTRICT nelcon ,short  *RESTRICT nen    
              ,short  *RESTRICT nFace  ,short  *RESTRICT geomType
              ,DOUBLE *RESTRICT gCc    ,DOUBLE *RESTRICT gKsi  
              ,DOUBLE *RESTRICT gmKsi  ,DOUBLE *RESTRICT gEta 
              ,DOUBLE *RESTRICT gmEta  ,DOUBLE *RESTRICT gNormal
              ,DOUBLE *RESTRICT gVolume,DOUBLE *RESTRICT gXm    
              ,DOUBLE *RESTRICT gXmcc  
              ,DOUBLE *RESTRICT gvSkew ,DOUBLE *RESTRICT gmSkew
              ,DOUBLE *RESTRICT gDcca               
              ,short maxNo             ,short maxViz
              ,short ndm               ,INT numel);
/*...................................................................*/

/*.. reconstrucao de gradiente*/
  void rcGradU(Memoria *m             ,Loads *loads
         ,INT    *RESTRICT el         ,INT    *RESTRICT nelcon
         ,DOUBLE *RESTRICT cc         ,DOUBLE *RESTRICT x     
         ,short  *RESTRICT nen        ,short  *RESTRICT nFace
         ,short  *RESTRICT geomType   ,DOUBLE *RESTRICT prop 
         ,short  *RESTRICT mat 
         ,DOUBLE *RESTRICT leastSquare,DOUBLE *RESTRICT leastSquareR   
         ,DOUBLE *RESTRICT gKsi       ,DOUBLE *RESTRICT gmKsi 
         ,DOUBLE *RESTRICT gEta       ,DOUBLE *RESTRICT gmEta 
         ,DOUBLE *RESTRICT gNormal    ,DOUBLE *RESTRICT gVolume
         ,DOUBLE *RESTRICT gvSkew     
         ,DOUBLE *RESTRICT gXm        ,DOUBLE *RESTRICT gXmcc
         ,DOUBLE *RESTRICT gDcca 
         ,short  *RESTRICT faceR      ,short *RESTRICT faceL  
         ,DOUBLE *RESTRICT u          ,DOUBLE *RESTRICT gradU               
         ,DOUBLE *RESTRICT nU         ,short const lib
         ,short const maxNo           ,short const maxViz
         ,short const ndf             ,short const ndm
         ,InterfaceNo *iNo            ,Interface *iCel 
         ,INT const numelNov          ,INT const numel 
         ,INT const nNodeNov          ,INT const nNode);
/*...................................................................*/

/*... */
  void interCellNode(Memoria *m           ,Loads *loads
                   ,DOUBLE *RESTRICT noU  ,DOUBLE *RESTRICT elU
                   ,INT *RESTRICT el      ,short  *RESTRICT geomType 
                   ,DOUBLE *RESTRICT cc   ,DOUBLE *RESTRICT x
                   ,DOUBLE *RESTRICT xm   
                   ,short *RESTRICT nen   ,short *RESTRICT nFace
                   ,short  *RESTRICT faceR,short *RESTRICT faceL  
                   ,InterfaceNo *iNo        
                   ,INT const numelNov    ,INT const numel       
                   ,INT const nNodeNov    ,INT const nNode 
                   ,short const maxNo     ,short const maxViz    
                   ,short const ndf1      ,short const ndf2  
                   ,short const ndm
                   ,bool const fBc        ,short const type);
/*...................................................................*/

/*... */
  void cellLibTurbulence(Turbulence tModel, 
             short *RESTRICT lGeomType  , DOUBLE *RESTRICT lprop,
             INT   *RESTRICT lViz       , DOUBLE *RESTRICT ksi,
             DOUBLE *RESTRICT mKsi      ,
             DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT fArea,
             DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume,
             DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc,
             DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT cc,
             DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew,
             DOUBLE *RESTRICT gradVel   , DOUBLE *RESTRICT lDensity, 
             DOUBLE *lViscosity         ,
             short const nEn            , short  const nFace,
             short const ndm            , short const lib,
             INT const nel);
/*...................................................................*/

/* ... montagem do sistemas de equacoes (difusao)*/
 void systFormDif(Loads *loads
             ,INT    *RESTRICT el  ,INT    *RESTRICT nelcon 
             ,short  *RESTRICT nen     ,short  *RESTRICT nFace
             ,short  *RESTRICT geomType,DOUBLE *RESTRICT prop 
             ,short  *RESTRICT calType ,short  *RESTRICT mat  
             ,DOUBLE *RESTRICT gKsi   ,DOUBLE *RESTRICT gmKsi
             ,DOUBLE *RESTRICT gEta   ,DOUBLE *RESTRICT gfArea
             ,DOUBLE *RESTRICT gNormal,DOUBLE *RESTRICT gVolume
             ,DOUBLE *RESTRICT gXm    ,DOUBLE *RESTRICT gXmcc
             ,DOUBLE *RESTRICT gvSkew ,DOUBLE *RESTRICT gmvSkew
             ,DOUBLE *RESTRICT gDcca  ,DOUBLE *RESTRICT density 
             ,INT    *RESTRICT ia     ,INT    *RESTRICT ja
             ,DOUBLE *RESTRICT al     ,DOUBLE *RESTRICT ad
             ,DOUBLE *RESTRICT b      ,INT    *RESTRICT id
             ,short  *RESTRICT faceR  ,short  *RESTRICT faceL
             ,DOUBLE *RESTRICT u0     ,DOUBLE *RESTRICT gradU0
             ,DOUBLE *RESTRICT rCell  ,Temporal const ddt
             ,INT const nEq           ,INT const nEqNov
             ,INT const nAd           ,INT const nAdR
             ,short const maxNo       ,short const maxViz
             ,short const ndm         ,INT const numel
             ,short const ndf         ,short const storage
             ,bool const forces       ,bool const matrix 
             ,bool const calRcell     ,bool const unsym); 
/*...................................................................*/
 
/* ... montagem do sistemas de equacoes (Transporte)*/
  void systFormTrans(Loads *loads      
             ,Advection advT           ,Diffusion diffT
             ,INT    *RESTRICT el      ,INT    *RESTRICT nelcon 
             ,short  *RESTRICT nen     ,short  *RESTRICT nFace
             ,short  *RESTRICT geomType,DOUBLE *RESTRICT prop 
             ,short  *RESTRICT calType ,short  *RESTRICT mat  
             ,DOUBLE *RESTRICT gCc
             ,DOUBLE *RESTRICT gKsi   ,DOUBLE *RESTRICT gmKsi
             ,DOUBLE *RESTRICT gEta   ,DOUBLE *RESTRICT gfArea
             ,DOUBLE *RESTRICT gNormal,DOUBLE *RESTRICT gVolume
             ,DOUBLE *RESTRICT gXm    ,DOUBLE *RESTRICT gXmcc   
             ,DOUBLE *RESTRICT gvSkew ,DOUBLE *RESTRICT gmvSkew   
             ,DOUBLE *RESTRICT gDcca  ,DOUBLE *RESTRICT density
             ,INT    *RESTRICT ia     ,INT    *RESTRICT ja
             ,DOUBLE *RESTRICT al     ,DOUBLE *RESTRICT ad 
             ,DOUBLE *RESTRICT b      ,INT    *RESTRICT id 
             ,short  *RESTRICT faceR  ,short  *RESTRICT faceL
             ,DOUBLE *RESTRICT u0     ,DOUBLE *RESTRICT gradU0 
             ,DOUBLE *RESTRICT vel 
             ,DOUBLE *RESTRICT rCell  ,Temporal const ddt
             ,INT const nEq           ,INT const nEqNov
             ,INT const nAd           ,INT const nAdR          
             ,short const maxNo       ,short const maxViz
             ,short const ndm         ,INT const numel
             ,short const ndf         ,short const storage
             ,bool const forces       ,bool const matrix 
             ,bool const calRcell     ,bool const unsym); 
/*...................................................................*/

/* ... montagem do sistemas de equacoes (Simple - VEl)*/
  void systFormSimpleVel(Loads *loadsVel,Loads *loadsPres 
             ,Advection advVel          ,Diffusion diffVel      
             ,short const typeSimple 
             ,INT    *RESTRICT el       ,INT    *RESTRICT nelcon 
             ,short  *RESTRICT nen      ,short  *RESTRICT nFace
             ,short  *RESTRICT geomType ,DOUBLE *RESTRICT prop 
             ,short  *RESTRICT calType  ,short  *RESTRICT mat  
             ,DOUBLE *RESTRICT gCc
             ,DOUBLE *RESTRICT gKsi     ,DOUBLE *RESTRICT gmKsi
             ,DOUBLE *RESTRICT gEta     ,DOUBLE *RESTRICT gfArea
             ,DOUBLE *RESTRICT gNormal  ,DOUBLE *RESTRICT gVolume
             ,DOUBLE *RESTRICT gXm      ,DOUBLE *RESTRICT gXmcc
             ,DOUBLE *RESTRICT gvSkew   ,DOUBLE *RESTRICT gmvSkew 
             ,DOUBLE *RESTRICT gDcca    ,DOUBLE *RESTRICT density
             ,INT    *RESTRICT ia       ,INT    *RESTRICT ja
             ,DOUBLE *RESTRICT al       ,DOUBLE *RESTRICT ad
             ,DOUBLE *RESTRICT b        ,INT    *RESTRICT id 
             ,short  *RESTRICT faceVelR ,short  *RESTRICT faceVelL 
             ,short  *RESTRICT facePresR,short  *RESTRICT facePresL 
             ,DOUBLE *RESTRICT pres     ,DOUBLE *RESTRICT gradPres
             ,DOUBLE *RESTRICT vel      ,DOUBLE *RESTRICT gradVel
             ,DOUBLE *RESTRICT dField   ,DOUBLE const underU
             ,DOUBLE *RESTRICT rCell    ,Temporal const ddt
             ,INT const nEq             ,INT const nEqNov
             ,INT const nAd             ,INT const nAdR 
             ,short const maxNo         ,short const maxViz
             ,short const ndm           ,INT const numel
             ,short const ndf           ,short const storage
             ,bool const forces         ,bool const matrix 
             ,bool const calRcell       ,bool const unsym
             ,const bool sPressure);  
/*...................................................................*/

/* ... montagem do sistemas de equacoes (Simple - VEl - low mach)*/
  void systFormSimpleVelLm(Loads *loadsVel,Loads *loadsPres 
      ,Advection advVel               ,Diffusion diffVel  
      ,Turbulence tModel              ,short const typeSimple              
      ,INT    *RESTRICT el            ,INT    *RESTRICT nelcon 
      ,short  *RESTRICT nen           ,short  *RESTRICT nFace
      ,short  *RESTRICT geomType      ,DOUBLE *RESTRICT prop 
      ,short  *RESTRICT calType       ,short  *RESTRICT mat  
      ,DOUBLE *RESTRICT gCc           ,DOUBLE *RESTRICT gKsi     
      ,DOUBLE *RESTRICT gmKsi         ,DOUBLE *RESTRICT gEta     
      ,DOUBLE *RESTRICT gfArea        ,DOUBLE *RESTRICT gNormal 
      ,DOUBLE *RESTRICT gVolume       ,DOUBLE *RESTRICT gXm      
      ,DOUBLE *RESTRICT gXmcc         ,DOUBLE *RESTRICT gvSkew   
      ,DOUBLE *RESTRICT gmvSkew       ,DOUBLE *RESTRICT gDcca    
      ,INT    *RESTRICT ia            ,INT    *RESTRICT ja
      ,DOUBLE *RESTRICT al            ,DOUBLE *RESTRICT ad
      ,DOUBLE *RESTRICT b             ,INT    *RESTRICT id 
      ,short  *RESTRICT faceVelR      ,short  *RESTRICT faceVelL 
      ,short  *RESTRICT facePresR     ,short  *RESTRICT facePresL 
      ,DOUBLE *RESTRICT pres          ,DOUBLE *RESTRICT gradPres
      ,DOUBLE *RESTRICT vel           ,DOUBLE *RESTRICT gradVel
      ,DOUBLE *RESTRICT dField        ,DOUBLE const underU
      ,DOUBLE *RESTRICT rCell  
      ,DOUBLE *RESTRICT density       ,DOUBLE *RESTRICT dViscosity   
       ,DOUBLE *RESTRICT eddyViscosity
      ,Temporal const ddt
      ,INT const nEq                  ,INT const nEqNov
      ,INT const nAd                  ,INT const nAdR 
      ,short const maxNo              ,short const maxViz
      ,short const ndm                ,INT const numel
      ,short const ndf                ,short const storage
      ,bool const forces              ,bool const matrix 
      ,bool const calRcell            ,bool const unsym
      ,const bool sPressure);  

/* ... montagem do sistemas de equacoes (Simple - Pres)*/
  void systFormSimplePresLm(Loads *loadsVel  ,Loads *loadsPres 
							 ,Diffusion diffPres
               ,INT    *RESTRICT el      ,INT    *RESTRICT nelcon 
               ,short  *RESTRICT nen     ,short  *RESTRICT nFace
               ,short  *RESTRICT geomType,DOUBLE *RESTRICT prop 
               ,short  *RESTRICT calType ,short  *RESTRICT mat     
               ,DOUBLE *RESTRICT gKsi    ,DOUBLE *RESTRICT gmKsi 
               ,DOUBLE *RESTRICT gEta    ,DOUBLE *RESTRICT gfArea 
               ,DOUBLE *RESTRICT gNormal ,DOUBLE *RESTRICT gVolume
               ,DOUBLE *RESTRICT gXm     ,DOUBLE *RESTRICT gXmcc 
               ,DOUBLE *RESTRICT gvSkew  ,DOUBLE *RESTRICT gmvSkew 
               ,DOUBLE *RESTRICT gDcca   ,DOUBLE *RESTRICT density
               ,INT    *RESTRICT ia      ,INT    *RESTRICT ja
               ,DOUBLE *RESTRICT a       ,DOUBLE *RESTRICT ad 
               ,DOUBLE *RESTRICT b       ,INT    *RESTRICT id
               ,short  *RESTRICT faceVelR ,short  *RESTRICT faceVelL       
               ,short  *RESTRICT facePresR,short  *RESTRICT facePresL      
               ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres
               ,DOUBLE *RESTRICT vel     ,DOUBLE *RESTRICT dField    
               ,DOUBLE *RESTRICT rCell   ,Temporal const ddt 
               ,INT const nEq            ,INT const nEqNov
               ,INT const nAd            ,INT const nAdR                  
               ,short const maxNo        ,short const maxViz
               ,short const ndm          ,INT const numel
               ,short const ndf          ,short const storage
               ,bool  const forces       ,bool const matrix 
               ,bool const calRcell      ,bool  const  unsym); 
/*...................................................................*/


/*...*/
  void systFormEnergy(Loads *loads        , EnergyModel model
           , Advection advE               , Diffusion diffE
           , Turbulence tModel              
           , INT    *RESTRICT el          , INT    *RESTRICT nelcon
           , short  *RESTRICT nen         , short  *RESTRICT nFace
           , short  *RESTRICT geomType    , DOUBLE *RESTRICT prop
           , short  *RESTRICT calType     , short  *RESTRICT mat
           , DOUBLE *RESTRICT gCc           
           , DOUBLE *RESTRICT gKsi        , DOUBLE *RESTRICT gmKsi
           , DOUBLE *RESTRICT gEta        , DOUBLE *RESTRICT gfArea
           , DOUBLE *RESTRICT gNormal     , DOUBLE *RESTRICT gVolume
           , DOUBLE *RESTRICT gXm         , DOUBLE *RESTRICT gXmcc
           , DOUBLE *RESTRICT gvSkew      , DOUBLE *RESTRICT gmvSkew
           , DOUBLE *RESTRICT gDcca         
           , INT    *RESTRICT ia          , INT    *RESTRICT ja
           , DOUBLE *RESTRICT a           , DOUBLE *RESTRICT ad
           , DOUBLE *RESTRICT b           , INT    *RESTRICT id
           , short  *RESTRICT faceR       , short  *RESTRICT faceL
           , DOUBLE *RESTRICT u0          , DOUBLE *RESTRICT gradU0
           , DOUBLE *RESTRICT vel         , DOUBLE *RESTRICT gradVel
           , DOUBLE *RESTRICT pres0       , DOUBLE *RESTRICT pres   
           , DOUBLE *RESTRICT gradPres    , DOUBLE *RESTRICT rCell                   
           , DOUBLE *RESTRICT density     , DOUBLE *RESTRICT sHeat
           , DOUBLE *RESTRICT dViscosity  , DOUBLE *RESTRICT eddyViscosity
           , DOUBLE *RESTRICT tConductvity, DOUBLE *RESTRICT dField
           , Temporal const ddt           , DOUBLE const underU
           , INT const nEq                , INT const nEqNov
           , INT const nAd                , INT const nAdR
           , short const maxNo            , short const maxViz
           , short const ndm              , INT const numel
           , short const ndf              , short const storage
           , bool  const forces           , bool const matrix
           , bool const calRcell          , bool  const  unsym);
/*...................................................................*/

/*... */
   void velExp(Loads *loadsVel           ,Loads *loadsPres
              ,Advection advVel          ,Diffusion diffVel
              ,INT    *RESTRICT el       ,INT    *RESTRICT nelcon
              ,short  *RESTRICT nen      ,short  *RESTRICT nFace
              ,short  *RESTRICT geomType ,DOUBLE *RESTRICT prop
              ,short  *RESTRICT calType  ,short  *RESTRICT mat
              ,DOUBLE *RESTRICT cc
              ,DOUBLE *RESTRICT gKsi     ,DOUBLE *RESTRICT gmKsi
              ,DOUBLE *RESTRICT gEta     ,DOUBLE *RESTRICT gfArea
              ,DOUBLE *RESTRICT gNormal  ,DOUBLE *RESTRICT gVolume
              ,DOUBLE *RESTRICT gXm      ,DOUBLE *RESTRICT gXmcc
              ,DOUBLE *RESTRICT gvSkew   ,DOUBLE *RESTRICT gmvSkew
              ,DOUBLE *RESTRICT gDcca    ,DOUBLE *RESTRICT density
              ,short  *RESTRICT faceVelR ,short  *RESTRICT faceVelL
              ,short  *RESTRICT facePresR,short  *RESTRICT facePresL
              ,DOUBLE *RESTRICT pres     ,DOUBLE *RESTRICT gradPres
              ,DOUBLE *RESTRICT vel      ,DOUBLE *RESTRICT velUp
              ,DOUBLE *RESTRICT gradVel  ,DOUBLE *RESTRICT bT
              ,DOUBLE *RESTRICT dField
              ,DOUBLE const underU       ,Temporal const ddt
              ,short const maxNo         ,short const maxViz
              ,short const ndm           ,INT const numel
              ,short const ndf           ,bool const sPressure
              ,bool const fResidual);
/*.....................................................................*/

/*...*/
   void velResidual(Loads *loadsVel      ,Loads *loadsPres
              ,Advection advVel          ,Diffusion diffVel
              ,INT    *RESTRICT el       ,INT    *RESTRICT nelcon
              ,short  *RESTRICT nen      ,short  *RESTRICT nFace
              ,short  *RESTRICT geomType ,DOUBLE *RESTRICT prop
              ,short  *RESTRICT calType  ,short  *RESTRICT mat
              ,DOUBLE *RESTRICT cc       ,DOUBLE *RESTRICT ad
              ,DOUBLE *RESTRICT gKsi     ,DOUBLE *RESTRICT gmKsi
              ,DOUBLE *RESTRICT gEta     ,DOUBLE *RESTRICT gfArea
              ,DOUBLE *RESTRICT gNormal  ,DOUBLE *RESTRICT gVolume
              ,DOUBLE *RESTRICT gXm      ,DOUBLE *RESTRICT gXmcc
              ,DOUBLE *RESTRICT gvSkew   ,DOUBLE *RESTRICT gmvSkew
              ,DOUBLE *RESTRICT gDcca    ,DOUBLE *RESTRICT density
              ,short  *RESTRICT faceVelR ,short  *RESTRICT faceVelL
              ,short  *RESTRICT facePresR,short  *RESTRICT facePresL
              ,DOUBLE *RESTRICT pres     ,DOUBLE *RESTRICT gradPres
              ,DOUBLE *RESTRICT vel      ,DOUBLE *RESTRICT res
              ,DOUBLE *RESTRICT gradVel  ,DOUBLE *RESTRICT bT
              ,DOUBLE *RESTRICT dField   ,DOUBLE const underU
              ,Temporal const ddt
              ,short const maxNo         ,short const maxViz
              ,short const ndm           ,INT const numel
              ,short const ndf           ,bool const sPressure
              ,bool const fResidual);
/*.....................................................................*/

/* ... montagem do sistemas de equacoes (Simple - Pres)*/
  void systFormSimplePres(Loads *loadsVel,Loads *loadsPres 
						 ,Diffusion diffPres
             ,INT    *RESTRICT el      ,INT    *RESTRICT nelcon 
             ,short  *RESTRICT nen     ,short  *RESTRICT nFace
             ,short  *RESTRICT geomType,DOUBLE *RESTRICT prop 
             ,short  *RESTRICT calType ,short  *RESTRICT mat  
             ,DOUBLE *RESTRICT gKsi   ,DOUBLE *RESTRICT gmKsi
             ,DOUBLE *RESTRICT gEta   ,DOUBLE *RESTRICT gfArea
             ,DOUBLE *RESTRICT gNormal,DOUBLE *RESTRICT gVolume
             ,DOUBLE *RESTRICT gXm    ,DOUBLE *RESTRICT gXmcc   
             ,DOUBLE *RESTRICT gvSkew ,DOUBLE *RESTRICT gmvSkew   
             ,DOUBLE *RESTRICT gDcca  ,DOUBLE *RESTRICT density
             ,INT    *RESTRICT ia     ,INT    *RESTRICT ja
             ,DOUBLE *RESTRICT al     ,DOUBLE *RESTRICT ad
             ,DOUBLE *RESTRICT b      ,INT    *RESTRICT id
             ,short  *RESTRICT faceVelR ,short  *RESTRICT faceVelL
             ,short  *RESTRICT facePresR,short  *RESTRICT facePresL
             ,DOUBLE *RESTRICT pres   ,DOUBLE *RESTRICT gradPres
             ,DOUBLE *RESTRICT vel    ,DOUBLE *RESTRICT dField 
             ,DOUBLE *RESTRICT rCell  ,Temporal const ddt
             ,INT const nEq           ,INT const nEqNov
             ,INT const nAd           ,INT const nAdR
             ,short const maxNo       ,short const maxViz
             ,short const ndm         ,INT const numel
             ,short const ndf         ,short const storage
             ,bool const forces       ,bool const matrix 
             ,bool const calRcell     ,bool const unsym); 
/*...................................................................*/

/*... correcao nao ortogonal para correcao de pressao*/
  void simpleNonOrthPres(Diffusion diffPres
               ,INT *RESTRICT el           ,INT *RESTRICT nelcon
							 ,short  *RESTRICT nen       ,short  *RESTRICT nFace
               ,short  *RESTRICT geomType  ,DOUBLE *RESTRICT prop 
               ,short  *RESTRICT calType   ,short  *RESTRICT mat
               ,DOUBLE *RESTRICT cc
               ,DOUBLE *RESTRICT gKsi      ,DOUBLE *RESTRICT gmKsi 
               ,DOUBLE *RESTRICT gEta      ,DOUBLE *RESTRICT gfArea 
               ,DOUBLE *RESTRICT gNormal   ,DOUBLE *RESTRICT gVolume
               ,DOUBLE *RESTRICT gXm       ,DOUBLE *RESTRICT gXmcc 
               ,DOUBLE *RESTRICT gvSkew    ,DOUBLE *RESTRICT gmvSkew 
               ,DOUBLE *RESTRICT gDcca     ,DOUBLE *RESTRICT density
               ,DOUBLE *RESTRICT b         ,INT    *RESTRICT id
               ,short  *RESTRICT facePresR ,DOUBLE *RESTRICT pres 
               ,DOUBLE *RESTRICT gradPres  ,DOUBLE *RESTRICT dField 
               ,short const maxNo          ,short const maxViz
               ,short const ndm            ,INT const numel);
/*...................................................................*/

/*... chamada da biblioteca de elementos (difusao)*/
  void cellLibDif(Loads *loads
                 ,short *RESTRICT lGeomType,DOUBLE *RESTRICT lprop
                 ,INT   *RESTRICT lViz     ,INT    *RESTRICT lId
                 ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mksi    
                 ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea 
                 ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume 
                 ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
                 ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity 
                 ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
                 ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
                 ,DOUBLE *RESTRICT lRcell  ,Temporal const ddt
                 ,short  *RESTRICT lFaceR  ,short  *RESTRICT lFaceL
                 ,DOUBLE *RESTRICT u0      ,DOUBLE *RESTRICT lGradU0 
                 ,short const nEn          ,short const nFace
                 ,short const ndm          ,short const lib 
                 ,INT const nel);
/*...................................................................*/

/*... chamada da biblioteca de elementos (transporte)*/
  void cellLibTrans(Loads *loads           
                 ,Advection advT           ,Diffusion diffT
                 ,short *RESTRICT lGeomType,DOUBLE *RESTRICT lprop
                 ,INT   *RESTRICT lViz     ,INT    *RESTRICT lId
                 ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mksi
                 ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea 
                 ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume 
                 ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
                 ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity 
                 ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
                 ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
                 ,DOUBLE *RESTRICT lRcell  ,Temporal const ddt
                 ,short  *RESTRICT lFaceR  ,short  *RESTRICT lFaceL
                 ,DOUBLE *RESTRICT lu0     ,DOUBLE *RESTRICT lGradU0 
                 ,DOUBLE *RESTRICT lVel    ,DOUBLE *RESTRICT cc
                 ,short const nEn          ,short const nFace
                 ,short const ndm          ,short const lib 
                 ,INT const nel);
/*...................................................................*/

/*... chamada da biblioteca de elementos (transporte)*/
  void cellLibEnergy(Loads *loads  , EnergyModel model
     , Advection  adv              , Diffusion diff
     , Turbulence tModel             
     , short *RESTRICT lGeomType   , DOUBLE *RESTRICT lprop
     , INT   *RESTRICT lViz        , INT *RESTRICT lId
     , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
     , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT fArea
     , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume
     , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
     , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
     , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
     , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
     , DOUBLE *RESTRICT lRcell     , Temporal const ddt
     , short  *RESTRICT lFaceR     , short  *RESTRICT lFaceL
     , DOUBLE *RESTRICT u          , DOUBLE *RESTRICT gradU
     , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
     , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres 
     , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lSheat
     , DOUBLE *RESTRICT lDviscosity, DOUBLE *RESTRICT ltConductvity
     , DOUBLE *RESTRICT dField
     , DOUBLE const underU
     , short const nEn             , short  const nFace
     , short const ndm             , short const lib
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
          ,Advection advVel           ,Diffusion diffVel
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
          ,DOUBLE *RESTRICT lRcell    ,Temporal const ddt
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
	           	 ,Diffusion diffPres
               ,short *RESTRICT lGeomType,DOUBLE *RESTRICT lprop
               ,INT   *RESTRICT lViz     ,INT *RESTRICT lId  
               ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
               ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea
               ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
               ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
               ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
               ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
               ,DOUBLE *RESTRICT lRcell  ,Temporal const ddt
               ,short  *RESTRICT lFaceVelR ,short  *RESTRICT lFaceVelL
               ,short  *RESTRICT lFacePresR,short  *RESTRICT lFacePresL
               ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres 
               ,DOUBLE *RESTRICT vel     ,DOUBLE *RESTRICT dField 
               ,short const nEn          ,short  const nFace     
               ,short const ndm          ,short const lib    
               ,INT const nel);
/*...................................................................*/

/*.......................... SIMPLE-LM ..............................*/
/*... chamada da biblioteca de elementos (escoamento-vel-Low Mach)*/
  void cellLibSimpleVelLm(Loads *loadsVel, Loads *loadsPres,
              Advection  advVel          , Diffusion diffVel, 
              Turbulence tModel          , short const typeSimple,       
              short *RESTRICT lGeomType  , DOUBLE *RESTRICT lprop,
              INT   *RESTRICT lViz       , INT *RESTRICT lId,
              DOUBLE *RESTRICT ksi       , DOUBLE *RESTRICT mKsi,
              DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT fArea,
              DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume,
              DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc,
              DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT cc,
              DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew,
              DOUBLE *RESTRICT lA        , DOUBLE *RESTRICT lB,
              DOUBLE *RESTRICT lRcell    , Temporal const ddt,
              short  *RESTRICT lFaceVelR , short  *RESTRICT lFaceVelL,
              short  *RESTRICT lFacePresR, short  *RESTRICT lFacePresL,
              DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres, 
              DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT gradVel,
              DOUBLE *RESTRICT lDensity  , DOUBLE *RESTRICT lViscosity,
              DOUBLE *RESTRICT dField    ,      
              DOUBLE const underU        , const bool sPressure,
              short const nEn            , short  const nFace,
              short const ndm            , short const lib,
              INT const nel);
/*... chamada da biblioteca de elementos (escoamento-pres-Low Mach)*/
  void cellLibSimplePresLm(Loads *loadsVel   , Loads *loadsPres
	             , Diffusion diffPres           
               , short *RESTRICT lGeomType  , DOUBLE *RESTRICT lprop
               , INT   *RESTRICT lViz       , INT *RESTRICT lId  
               , DOUBLE *RESTRICT ksi       , DOUBLE *RESTRICT mKsi
               , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT fArea
               , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume
               , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc
               , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT lDensity
               , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew
               , DOUBLE *RESTRICT lA        , DOUBLE *RESTRICT lB
               , DOUBLE *RESTRICT lRcell    , Temporal const ddt
               , short  *RESTRICT lFaceVelR , short  *RESTRICT lFaceVelL
               , short  *RESTRICT lFacePresR, short  *RESTRICT lFacePresL
               , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres 
               , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT dField 
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

/*...................................................................*/

/*... carga por elmento e condicoes pescritas por celula*/
  void cellPload(Loads *loads           ,DOUBLE *RESTRICT gCc
                ,short  *RESTRICT faceR ,short *RESTRICT faceS
                ,DOUBLE *RESTRICT volume,INT *RESTRICT id 
                ,DOUBLE *RESTRICT u     ,DOUBLE *RESTRICT f
                ,INT const numel        ,short const ndf
                ,short const ndm        ,short const maxViz);

/*... carga por elmento e condicoes pescritas por celula no 
      metodo simple*/
  void pLoadSimple(DOUBLE *RESTRICT sP  ,DOUBLE *RESTRICT p
          ,DOUBLE *RESTRICT tA          ,DOUBLE *RESTRICT velC
          ,DOUBLE *RESTRICT n
          ,DOUBLE *RESTRICT gradVel     ,DOUBLE *RESTRICT xmcc
          ,DOUBLE const viscosityC      ,DOUBLE const densityC
          ,DOUBLE const fArea           ,DOUBLE const dcca
          ,Loads ld                     ,short const ndm
          ,bool const fCal1             ,bool const fCal2);

  void pLoadSimplePres(DOUBLE *RESTRICT sP  ,DOUBLE *RESTRICT p
          ,DOUBLE *RESTRICT tA
          ,DOUBLE const viscosityC,DOUBLE const densityC
          ,DOUBLE const wfn
          ,DOUBLE const fArea     ,DOUBLE const dcca 
          ,Loads ld               ,bool const fCal);
/*...................................................................*/

/*... carga por elmento e condicoes pescritas por celula no metodo 
      simple*/
  void cellPloadSimple(Loads *loads           ,DOUBLE *RESTRICT cc 
                    ,short  *RESTRICT faceR ,short *RESTRICT faceL
                    ,DOUBLE *RESTRICT volume
                    ,INT *RESTRICT idVel    ,INT *RESTRICT idPres
                    ,DOUBLE *RESTRICT vel   ,DOUBLE *RESTRICT pres
                    ,DOUBLE *RESTRICT fVel  ,DOUBLE *RESTRICT fPres
                    ,INT const numel        ,short const ndf
                    ,short const ndm        ,short const maxViz);
/*...................................................................*/

/*... */
  void updateCellValue(DOUBLE *RESTRICT u,DOUBLE *RESTRICT x
                      ,INT *RESTRICT id  ,Interface *iNeq
                      ,INT const numel   ,short const ndf
                      ,bool const fAdd   ,short const fCom);
  
  void updateCellValueSimple(DOUBLE *RESTRICT u,DOUBLE *RESTRICT x
                      ,INT *RESTRICT id  ,Interface *iNeq
                      ,INT const numel   ,INT const nEq 
                      ,short const ndf
                      ,bool const fAdd   ,short const fCom);
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
  
  void cellDif3D(Loads *loads
                ,short *RESTRICT lGeomType,DOUBLE *RESTRICT prop
                ,INT *RESTRICT lViz       ,INT *RESTRICT lId  
                ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
                ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea
                ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
                ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
                ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity 
                ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
                ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
                ,DOUBLE *RESTRICT lRcell  ,Temporal const ddt   
                ,short  *RESTRICT lFaceR  ,short  *RESTRICT lFaceL
                ,DOUBLE *RESTRICT u0      ,DOUBLE *RESTRICT gradU0
                ,const short nEn          ,short const nFace    
                ,const short ndm          ,INT const nel);
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
                ,Advection advT           , Diffusion diffT
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

/*.......................... ENRGIA .................................*/
  void cellEnergy2D(Loads *loads           , EnergyModel model
            , Advection adv                , Diffusion diff
            , Turbulence tModel  
            , short *RESTRICT lGeomType    , DOUBLE *RESTRICT prop
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
            , DOUBLE *RESTRICT u0          , DOUBLE *RESTRICT gradU0           
            , DOUBLE *RESTRICT vel         , DOUBLE *RESTRICT gradVel  
            , DOUBLE *RESTRICT pres        , DOUBLE *RESTRICT gradPres         
            , DOUBLE *RESTRICT lDensity    , DOUBLE *RESTRICT lSheat
            , DOUBLE *RESTRICT lDviscosity , DOUBLE *RESTRICT ltConductvity
            , DOUBLE *RESTRICT dField
            , DOUBLE const underU            
            , const short nEn              , short const nFace
            , const short ndm              , INT const nel);
/*...................................................................*/

/*.......................... SIMPLE .................................*/

/*... biblioteca de celulas (simple - vel)*/
  void cellSimpleVel2D(Loads *loadsVel    ,Loads *loadsPres 
              ,Advection  advVel          ,Diffusion diffVel
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
              ,Advection  advVel          ,Diffusion  diffVel
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
/*...................................................................*/

/*... biblioteca de celulas (simple - pres)*/
  void cellSimplePres2D(Loads *loadsVel    ,Loads *loadsPres   
							,Diffusion diffPres        
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
  void cellSimplePres3D(Loads *loadsVel    ,Loads *loadsPres
							,Diffusion diffPres
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
 /*... biblioteca de celulas (simple - vel - low mach)*/
  void cellSimpleVel2DLm(Loads *loadsVel , Loads *loadsPres  
           , Advection advVel           , Diffusion diffVel 
           , Turbulence tModel          , short const typeSimple      
           , short *RESTRICT lGeomType  , DOUBLE *RESTRICT prop 
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
           , DOUBLE *RESTRICT dField       
           , DOUBLE const underU        , const bool sPressure 
           , const short nEn            , short const nFace 
           , const short ndm            , INT const nel);

  void cellSimplePres2DLm(Loads *loadsVel, Loads *loadsPres 
					 , Diffusion diffPres	
           , short *RESTRICT lGeomType  , DOUBLE *RESTRICT prop
           , INT *RESTRICT lViz         , INT *RESTRICT lId  
           , DOUBLE *RESTRICT ksi       , DOUBLE *RESTRICT mKsi
           , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT mEta
           , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT area   
           , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc
           , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT lDensity
           , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew
           , DOUBLE *RESTRICT lA        , DOUBLE *RESTRICT lB
           , DOUBLE *RESTRICT lRcell    ,  Temporal const ddt
           , short  *RESTRICT lFaceVelR , short *RESTRICT lFaceVelL
           , short  *RESTRICT lFacePresR, short *RESTRICT lFacePresL
           , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres 
           , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT dField  
           , const short nEn            , short const nFace    
           , const short ndm            , INT const nel);
/*...................................................................*/

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

  void rcLeastSquare(DOUBLE *RESTRICT gKsi    ,DOUBLE *RESTRICT gmKsi
                  ,DOUBLE *RESTRICT lSquare   ,DOUBLE *RESTRICT lSquareR
                  ,short *RESTRICT nFace       
                  ,INT const numel          ,short const maxViz
                  ,short const type         ,short const ndm);
/*...................................................................*/

/*... least square*/
  void leastSquareMatrix(DOUBLE *RESTRICT lKsi    ,DOUBLE *RESTRICT lmKsi
                    ,DOUBLE *RESTRICT lLsquare,DOUBLE *RESTRICT lSquareR
                    ,short const type         
                    ,short const lnFace       ,short const ndm);


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
                ,short  *RESTRICT nFace   ,DOUBLE *RESTRICT volume
                ,DOUBLE *RESTRICT gKsi    ,DOUBLE *RESTRICT gNormal
                ,DOUBLE *RESTRICT gmvSkew
                ,short const maxViz      ,short const ndm
                ,INT const numel); 
/*...................................................................*/

/*...*/
  bool openDomain(short  *RESTRICT faceR, short  *RESTRICT nFace
                 , INT const numel      , short const maxViz);

  void wallFluid(short *RESTRICT faceR ,INT *RESTRICT nelcon
              ,short *RESTRICT nFace     
              ,INT const nEl         ,short const maxViz);

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
  DOUBLE nvd(DOUBLE const phiTil, short const iCod);
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

  void difusionSchemeAnisotropic(DOUBLE *RESTRICT s,DOUBLE *RESTRICT ksi
                                ,DOUBLE *RESTRICT e,DOUBLE *RESTRICT t
                                ,short const ndm   ,short const iCod);
  void advectiveScheme(DOUBLE *RESTRICT velC       ,DOUBLE *RESTRICT velV
                     ,DOUBLE *RESTRICT gradVelC   ,DOUBLE *RESTRICT gradVelV
                     ,DOUBLE *RESTRICT gradVelComp,DOUBLE *RESTRICT vSkew
                     ,DOUBLE *RESTRICT rC         ,DOUBLE *RESTRICT rV
                     ,DOUBLE *RESTRICT ksi        ,DOUBLE const modKsi
                     ,DOUBLE const m              ,DOUBLE *RESTRICT cvc
                     ,short const ndm
                     ,short const iCod1, short const iCod2);
   
  void advectiveSchemeScalar(DOUBLE const uC     ,DOUBLE const uV
                     ,DOUBLE *RESTRICT gradUc   ,DOUBLE *RESTRICT gradUv
                     ,DOUBLE *RESTRICT gradUcomp,DOUBLE *RESTRICT vSkew
                     ,DOUBLE *RESTRICT rC       ,DOUBLE *RESTRICT rV
                     ,DOUBLE *RESTRICT ksi      ,DOUBLE const modKsi
                     ,DOUBLE const m            ,DOUBLE *cvc
                     ,short const ndm
                     ,short const iCod1, short const iCod2);

  DOUBLE deferredCd(DOUBLE const velC,DOUBLE const velV
                   ,DOUBLE const wfn);

  DOUBLE upwindLinearV1(DOUBLE const uC  ,DOUBLE const uV
                 ,DOUBLE *RESTRICT gradUc,DOUBLE *RESTRICT gradUv
                 ,DOUBLE *RESTRICT rC    ,DOUBLE *RESTRICT rV
                 ,DOUBLE const m         ,short const ndm);
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
  
  void pLoad(DOUBLE *RESTRICT sP  ,DOUBLE *RESTRICT p
            ,DOUBLE *RESTRICT tA
            ,DOUBLE const coefDifC,DOUBLE const densityC
            ,DOUBLE const wfn     ,DOUBLE *RESTRICT xm                   
            ,DOUBLE const fArea   ,DOUBLE const dcca
            ,Loads ld             ,bool const fCal);
  void setNvd(char *word,short *iCod);
  void setTvd(char *word,short *iCod);
  void  setDiffusionScheme(char *word,short *iCod);
  void  setAdvectionScheme(char *word, Advection *adv,FILE *fileIn);
  DOUBLE sizeCar(DOUBLE const volume,short const ndm);
/*...................................................................*/

#endif/*_CELLLOOP_H_*/
