#ifndef _CELLLOOP_H_
  #define _CELLLOOP_H_
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<Assbly.h>
  #include<CellLib.h>
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
  void rcGradU(Memoria *m             , Loads *loads
         , INT    *RESTRICT el         , INT    *RESTRICT nelcon
         , DOUBLE *RESTRICT cc         , DOUBLE *RESTRICT x     
         , short  *RESTRICT nen        , short  *RESTRICT nFace
         , short  *RESTRICT geomType   , DOUBLE *RESTRICT prop 
         , short  *RESTRICT mat          
         , DOUBLE *RESTRICT leastSquare, DOUBLE *RESTRICT leastSquareR   
         , DOUBLE *RESTRICT gKsi       , DOUBLE *RESTRICT gmKsi 
         , DOUBLE *RESTRICT gEta       , DOUBLE *RESTRICT gmEta 
         , DOUBLE *RESTRICT gNormal    , DOUBLE *RESTRICT gVolume
         , DOUBLE *RESTRICT gvSkew       
         , DOUBLE *RESTRICT gXm        , DOUBLE *RESTRICT gXmcc
         , DOUBLE *RESTRICT gDcca        
         , short  *RESTRICT faceR      , short *RESTRICT faceL  
         , DOUBLE *RESTRICT u          , DOUBLE *RESTRICT gradU               
         , DOUBLE *RESTRICT nU         , short lib
         , short maxNo                 , short maxViz
         , short ndf                   , short ndm
         , InterfaceNo *iNo            , Interface *iCel 
         , INT numelNov                , INT numel 
         , INT nNodeNov                , INT nNode);
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
               ,DOUBLE *RESTRICT gKsi    ,DOUBLE *RESTRICT gmKsi 
               ,DOUBLE *RESTRICT gEta    ,DOUBLE *RESTRICT gfArea 
               ,DOUBLE *RESTRICT gNormal ,DOUBLE *RESTRICT gVolume
               ,DOUBLE *RESTRICT gXm     ,DOUBLE *RESTRICT gXmcc 
               ,DOUBLE *RESTRICT gvSkew  ,DOUBLE *RESTRICT gmvSkew 
               ,DOUBLE *RESTRICT gDcca   ,DOUBLE *RESTRICT density
               ,INT    *RESTRICT ia      ,INT    *RESTRICT ja
               ,DOUBLE *RESTRICT a       ,DOUBLE *RESTRICT ad 
               ,DOUBLE *RESTRICT b       ,INT    *RESTRICT id
               ,short  *RESTRICT faceR   ,short  *RESTRICT faceL        
               ,DOUBLE *RESTRICT u0      ,DOUBLE *RESTRICT gradU0 
               ,DOUBLE *RESTRICT vel                               
               ,DOUBLE *RESTRICT rCell   ,Temporal ddt 
               ,INT nEq                  ,INT nEqNov
               ,INT nAd                  ,INT nAdR                     
               ,short maxNo              ,short maxViz
               ,short ndm                ,INT numel
               ,short ndf                ,short storage
               ,bool forces              ,bool matrix 
               ,bool calRcell            ,bool unsym);
/*...................................................................*/

/* ... montagem do sistemas de equacoes (Simple - VEl)*/
  void systFormSimpleVel(Loads *loadsVel   ,Loads *loadsPres    
             ,Advection advVel           ,Diffusion diffVel
             ,short typeSimple 
             ,INT    *RESTRICT el        ,INT    *RESTRICT nelcon 
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
             ,INT    *RESTRICT ia        ,INT    *RESTRICT ja
             ,DOUBLE *RESTRICT a         ,DOUBLE *RESTRICT ad 
             ,DOUBLE *RESTRICT b         ,INT    *RESTRICT id
             ,short  *RESTRICT faceVelR ,short  *RESTRICT faceVelL       
             ,short  *RESTRICT facePresR,short  *RESTRICT facePresL             
             ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres
             ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT gradVel
             ,DOUBLE *RESTRICT dField    ,DOUBLE underU  
             ,DOUBLE *RESTRICT rCell     ,Temporal ddt 
             ,INT nEq                    ,INT nEqNov
             ,INT nAd                    ,INT nAdR                  
             ,short maxNo                ,short maxViz
             ,short ndm                  ,INT numel
             ,short ndf                  ,short storage
             ,bool forces                ,bool matrix 
             ,bool calRcell              ,bool unsym
             ,bool sPressure); 
/*...................................................................*/

/* ... montagem do sistemas de equacoes (Simple - VEl - low mach)*/
void systFormSimpleVelLm(Loads *loadsVel   , Loads *loadsPres    
    , Advection advVel                     , Diffusion diffVel
    , Turbulence tModel                    , MomentumModel eMomentum
    , short typeSimple     
    , INT    *RESTRICT el                  , INT    *RESTRICT nelcon 
    , short  *RESTRICT nen                 , short  *RESTRICT nFace
    , short  *RESTRICT geomType            , DOUBLE *RESTRICT prop
    , short  *RESTRICT calType             , short  *RESTRICT mat     
    , DOUBLE *RESTRICT cc                  , DOUBLE *RESTRICT gKsi
    , DOUBLE *RESTRICT gmKsi               , DOUBLE *RESTRICT gEta                
    , DOUBLE *RESTRICT gfArea              , DOUBLE *RESTRICT gNormal             
    , DOUBLE *RESTRICT gVolume             , DOUBLE *RESTRICT gXm
    , DOUBLE *RESTRICT gXmcc               , DOUBLE *RESTRICT gvSkew              
    , DOUBLE *RESTRICT gmvSkew             , DOUBLE *RESTRICT gDcca                   
    , INT    *RESTRICT ia                  , INT    *RESTRICT ja
    , DOUBLE *RESTRICT a                   , DOUBLE *RESTRICT ad
    , DOUBLE *RESTRICT b                   , INT    *RESTRICT id
    , short  *RESTRICT faceVelR            , short  *RESTRICT faceVelL       
    , short  *RESTRICT facePresR           , short  *RESTRICT facePresL             
    , DOUBLE *RESTRICT pres                , DOUBLE *RESTRICT gradPres
    , DOUBLE *RESTRICT vel                 , DOUBLE *RESTRICT gradVel
    , DOUBLE *RESTRICT dField              , DOUBLE underU 
    , DOUBLE *RESTRICT rCell               , DOUBLE *RESTRICT stressR  
    , DOUBLE *RESTRICT density             , DOUBLE *RESTRICT dViscosity 
    , DOUBLE *RESTRICT eddyViscosity       , DOUBLE *RESTRICT wallPar
    , Temporal ddt                     
    , INT const nEq                        , INT const nEqNov
    , INT const nAd                        , INT const nAdR                 
    , short const maxNo                    , short const maxViz
    , short const ndm                      , INT const numel
    , short const ndf                      , short const storage
    , short const ntn                      , bool const forces      
    , bool const matrix                    , bool const calRcell
    , const bool unsym                     , bool const sPressure); 

/* ... montagem do sistemas de equacoes (Simple - Pres)*/
  void systFormSimplePresLm(Loads *loadsVel, Loads *loadsPres 
							 , Diffusion diffPres        , MassEqModel eMass  
               , INT    *RESTRICT el       , INT    *RESTRICT nelcon 
               , short  *RESTRICT nen      , short  *RESTRICT nFace
               , short  *RESTRICT geomType , DOUBLE *RESTRICT prop 
               , short  *RESTRICT calType  , short  *RESTRICT mat     
               , DOUBLE *RESTRICT gKsi     , DOUBLE *RESTRICT gmKsi 
               , DOUBLE *RESTRICT gEta     , DOUBLE *RESTRICT gfArea 
               , DOUBLE *RESTRICT gNormal  , DOUBLE *RESTRICT gVolume
               , DOUBLE *RESTRICT gXm      , DOUBLE *RESTRICT gXmcc 
               , DOUBLE *RESTRICT gvSkew   , DOUBLE *RESTRICT gmvSkew 
               , DOUBLE *RESTRICT gDcca    , DOUBLE *RESTRICT density
               , INT    *RESTRICT ia       , INT    *RESTRICT ja
               , DOUBLE *RESTRICT a        , DOUBLE *RESTRICT ad 
               , DOUBLE *RESTRICT b        , INT    *RESTRICT id
               , short  *RESTRICT faceVelR , short  *RESTRICT faceVelL       
               , short  *RESTRICT facePresR, short  *RESTRICT facePresL      
               , DOUBLE *RESTRICT pres     , DOUBLE *RESTRICT gradPres
               , DOUBLE *RESTRICT vel      , DOUBLE *RESTRICT dField   
               , DOUBLE *RESTRICT temp     , DOUBLE *RESTRICT wallPar  
               , DOUBLE *RESTRICT rCell    , Temporal const ddt 
               , INT const nEq             , INT const nEqNov
               , INT const nAd             , INT const nAdR                  
               , short const maxNo         , short const maxViz
               , short const ndm           , INT const numel
               , short const ndf           , short const storage
               , bool  const forces        , bool const matrix 
               , bool const calRcell       , bool  const  unsym); 
/*...................................................................*/

/*...*/
void systFormEnergy(Loads *loads       , Loads *ldVel  
       , Advection adv                 , Diffusion diff 
       , Turbulence tModel             , EnergyModel eModel 
       , PropVar vProp
       , INT    *RESTRICT el           , INT    *RESTRICT nelcon
       , short  *RESTRICT nen          , short  *RESTRICT nFace
       , short  *RESTRICT geomType     , DOUBLE *RESTRICT prop
       , short  *RESTRICT calType      , short  *RESTRICT mat
       , DOUBLE *RESTRICT gCc          , DOUBLE *RESTRICT gKsi
       , DOUBLE *RESTRICT gmKsi        , DOUBLE *RESTRICT gEta 
       , DOUBLE *RESTRICT gfArea       , DOUBLE *RESTRICT gNormal      
       , DOUBLE *RESTRICT gVolume      , DOUBLE *RESTRICT gXm          
       , DOUBLE *RESTRICT gXmcc        , DOUBLE *RESTRICT gvSkew 
       , DOUBLE *RESTRICT gmvSkew      , DOUBLE *RESTRICT gDcca  
       , INT    *RESTRICT ia           , INT    *RESTRICT ja
       , DOUBLE *RESTRICT a            , DOUBLE *RESTRICT ad
       , DOUBLE *RESTRICT b            , INT    *RESTRICT id
       , short  *RESTRICT faceR        , short  *RESTRICT faceL
       , short  *RESTRICT faceVelR     , short  *RESTRICT faceVelL
       , DOUBLE *RESTRICT u0           , DOUBLE *RESTRICT gradU0
       , DOUBLE *RESTRICT vel          , DOUBLE *RESTRICT gradVel
       , DOUBLE *RESTRICT pres0        , DOUBLE *RESTRICT pres 
       , DOUBLE *RESTRICT gradPres     , DOUBLE *RESTRICT rCell
       , DOUBLE *RESTRICT density      , DOUBLE *RESTRICT sHeat
       , DOUBLE *RESTRICT dViscosity   , DOUBLE *RESTRICT eddyViscosity
       , DOUBLE *RESTRICT tConductivity, DOUBLE *RESTRICT dField
       , DOUBLE *RESTRICT wallPar
       , Temporal ddt                  , DOUBLE underU
       , INT nEq                       , INT nEqNov
       , INT nAd                       , INT nAdR
       , short maxNo                   , short maxViz
       , short ndm                     , INT numel
       , short ndf                     , short storage
       , bool forces                   , bool matrix
       , bool calRcell                 , bool unsym);
/*...................................................................*/

/*...*/
  void systFormOneEqK(Loads *lds     ,Loads *ldsVel
       , Turbulence *tModel
       , Advection adv                 , Diffusion diff 
       , INT    *RESTRICT el           , INT    *RESTRICT nelcon
       , short  *RESTRICT nen          , short  *RESTRICT nFace
       , short  *RESTRICT geomType     , DOUBLE *RESTRICT prop
       , short  *RESTRICT calType      , short  *RESTRICT mat
       , DOUBLE *RESTRICT gCc            
       , DOUBLE *RESTRICT gKsi         , DOUBLE *RESTRICT gmKsi
       , DOUBLE *RESTRICT gEta         , DOUBLE *RESTRICT gfArea
       , DOUBLE *RESTRICT gNormal      , DOUBLE *RESTRICT gVolume
       , DOUBLE *RESTRICT gXm          , DOUBLE *RESTRICT gXmcc
       , DOUBLE *RESTRICT gvSkew       , DOUBLE *RESTRICT gmvSkew
       , DOUBLE *RESTRICT gDcca          
       , INT    *RESTRICT ia           , INT    *RESTRICT ja
       , DOUBLE *RESTRICT a            , DOUBLE *RESTRICT ad
       , DOUBLE *RESTRICT b            , INT    *RESTRICT id
       , short  *RESTRICT faceReK      , short  *RESTRICT faceLdK
       , short  *RESTRICT faceReVel    , short  *RESTRICT faceLdVel
       , DOUBLE *RESTRICT u0           , DOUBLE *RESTRICT gradU0
       , DOUBLE *RESTRICT vel          , DOUBLE *RESTRICT gradVel
       , DOUBLE *RESTRICT pres         , DOUBLE *RESTRICT gradPres    
       , DOUBLE *RESTRICT density      , DOUBLE *RESTRICT dViscosity   
       , DOUBLE *RESTRICT eddyViscosity, DOUBLE *RESTRICT dField
       , DOUBLE *RESTRICT rCell        , DOUBLE *RESTRICT wallPar
       , DOUBLE *RESTRICT cDyn         , Temporal ddt  
       , INT nEq                       , INT nEqNov
       , INT nAd                       , INT nAdR
       , short maxNo                   , short maxViz
       , short ndm                     , INT numel
       , short ndf                     , short storage
       , bool forces                   , bool matrix
       , bool calRcell                 , bool unsym);
/*...................................................................*/

/*... */
 void velExp(Loads *loadsVel        ,Loads *loadsPres
             ,Advection advVel           ,Diffusion diffVel
             ,INT    *RESTRICT el        ,INT    *RESTRICT nelcon
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
             ,short  *RESTRICT faceVelR  ,short  *RESTRICT faceVelL
             ,short  *RESTRICT facePresR ,short  *RESTRICT facePresL
             ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres
             ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT velUp
             ,DOUBLE *RESTRICT gradVel   ,DOUBLE *RESTRICT bT
             ,DOUBLE *RESTRICT dField    ,DOUBLE underU
             ,Temporal ddt
             ,short maxNo                ,short maxViz
             ,short ndm                  ,INT numel
             ,short ndf                  ,bool sPressure
             ,bool fResidual);
/*.....................................................................*/

/*...*/
   void velResidual(Loads *loadsVel      , Loads *loadsPres
              ,Advection advVel          , Diffusion diffVel
              ,INT    *RESTRICT el       , INT    *RESTRICT nelcon
              ,short  *RESTRICT nen      , short  *RESTRICT nFace
              ,short  *RESTRICT geomType , DOUBLE *RESTRICT prop
              ,short  *RESTRICT calType  , short  *RESTRICT mat
              ,DOUBLE *RESTRICT cc       , DOUBLE *RESTRICT ad
              ,DOUBLE *RESTRICT gKsi     , DOUBLE *RESTRICT gmKsi
              ,DOUBLE *RESTRICT gEta     , DOUBLE *RESTRICT gfArea
              ,DOUBLE *RESTRICT gNormal  , DOUBLE *RESTRICT gVolume
              ,DOUBLE *RESTRICT gXm      , DOUBLE *RESTRICT gXmcc
              ,DOUBLE *RESTRICT gvSkew   , DOUBLE *RESTRICT gmvSkew
              ,DOUBLE *RESTRICT gDcca    , DOUBLE *RESTRICT density
              ,short  *RESTRICT faceVelR , short  *RESTRICT faceVelL
              ,short  *RESTRICT facePresR, short  *RESTRICT facePresL
              ,DOUBLE *RESTRICT pres     , DOUBLE *RESTRICT gradPres
              ,DOUBLE *RESTRICT vel      , DOUBLE *RESTRICT res
              ,DOUBLE *RESTRICT gradVel  , DOUBLE *RESTRICT bT
              ,DOUBLE *RESTRICT dField   , DOUBLE underU
              ,Temporal ddt                
              ,short maxNo               , short maxViz
              ,short ndm                 , INT numel
              ,short ndf                 , bool sPressure
              ,bool fResidual);
/*.....................................................................*/

/* ... montagem do sistemas de equacoes (Simple - Pres)*/
  void systFormSimplePres(Loads *lVel    , Loads *lPres 
						 , Diffusion diffPres           
             , INT    *RESTRICT el        , INT    *RESTRICT nelcon 
             , short  *RESTRICT nen       , short  *RESTRICT nFace
             , short  *RESTRICT geomType  , DOUBLE *RESTRICT prop 
             , short  *RESTRICT calType   , short  *RESTRICT mat  
             , DOUBLE *RESTRICT gKsi      , DOUBLE *RESTRICT gmKsi
             , DOUBLE *RESTRICT gEta      , DOUBLE *RESTRICT gfArea
             , DOUBLE *RESTRICT gNormal   , DOUBLE *RESTRICT gVolume
             , DOUBLE *RESTRICT gXm       , DOUBLE *RESTRICT gXmcc   
             , DOUBLE *RESTRICT gvSkew    , DOUBLE *RESTRICT gmvSkew   
             , DOUBLE *RESTRICT gDcca     , DOUBLE *RESTRICT density
             , INT    *RESTRICT ia        , INT    *RESTRICT ja
             , DOUBLE *RESTRICT al        , DOUBLE *RESTRICT ad
             , DOUBLE *RESTRICT b         , INT    *RESTRICT id
             , short  *RESTRICT faceVelR  , short  *RESTRICT faceVelL
             , short  *RESTRICT facePresR , short  *RESTRICT facePresL
             , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres
             , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT dField 
             , DOUBLE *RESTRICT rCell     , Temporal ddt
             , INT nEq                    , INT nEqNov
             , INT nAd                    , INT nAdR
             , short maxNo                , short maxViz
             , short ndm                  , INT cnumel
             , short ndf                  , short storage
             , bool forces                , bool matrix 
             , bool calRcell              , bool unsym); 
/*...................................................................*/

/*... correcao nao ortogonal para correcao de pressao*/
  void simpleNonOrthPres(Diffusion diffPres
               , INT *RESTRICT el           , INT *RESTRICT nelcon
							 , short  *RESTRICT nen       , short  *RESTRICT nFace
               , short  *RESTRICT geomType  , DOUBLE *RESTRICT prop 
               , short  *RESTRICT calType   , short  *RESTRICT mat
               , DOUBLE *RESTRICT cc          
               , DOUBLE *RESTRICT gKsi      , DOUBLE *RESTRICT gmKsi 
               , DOUBLE *RESTRICT gEta      , DOUBLE *RESTRICT gfArea 
               , DOUBLE *RESTRICT gNormal   , DOUBLE *RESTRICT gVolume
               , DOUBLE *RESTRICT gXm       , DOUBLE *RESTRICT gXmcc 
               , DOUBLE *RESTRICT gvSkew    , DOUBLE *RESTRICT gmvSkew 
               , DOUBLE *RESTRICT gDcca     , DOUBLE *RESTRICT density
               , DOUBLE *RESTRICT b         , INT    *RESTRICT id
               , short  *RESTRICT facePresR , DOUBLE *RESTRICT pres 
               , DOUBLE *RESTRICT gradPres  , DOUBLE *RESTRICT dField 
               , short maxNo                , short maxViz
               , short ndm                  , INT numel);
/*...................................................................*/

/*... carga por elmento e condicoes pescritas por celula*/
  void cellPload(Loads *loads           ,DOUBLE *RESTRICT gCc
                ,short  *RESTRICT faceR ,short *RESTRICT faceS
                ,DOUBLE *RESTRICT volume,INT *RESTRICT id 
                ,DOUBLE *RESTRICT u     ,DOUBLE *RESTRICT f
                ,INT const numel        ,short const ndf
                ,short const ndm        ,short const maxViz);
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

#endif/*_CELLLOOP_H_*/
