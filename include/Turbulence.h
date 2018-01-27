#ifndef _TURBULENCE_H_
  #define  _TURBULENCE_H_
/*...*/
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<CellLoop.h>
  #include<Define.h>
  #include<Erro.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<NumInt.h>
  #include<Jacobi.h>
  #include<Sisteq.h>
  #include<Solv.h>
/*...................................................................*/

/*...*/
  #define E_WALLMODEL   9.793e0
  #define VANDRIEST    26.e0
  #define VONKARMAN     0.4187e0
  #define AKADER        0.01e0
  #define BKADER        5.e0
/*...................................................................*/

/*...*/
  #define SMAGORINSKY 1
  #define WALEMODEL   2
  #define VREMAN      3
  #define DYNAMIC     4
  #define SIGMAMODEL  5
  #define MIXED       6
  #define BARDINA     7
  #define CLARK       8 
  #define BARDINAMOD  9
  #define ONEEQK     10
/*...................................................................*/
  
/*...*/
  #define LESFUNCMODEL       0
  #define LESSTRUMODEL       1
  #define LESMIXEDMODEL      2
  #define LESMIXEDTWOMODEL   3
  #define LESFUNCMODELONEEQK 4
/*...................................................................*/

/*...*/
  #define LES 0
/*...*/
  #define ESTMODEL 0
  #define FUNMODEL 1
/*...................................................................*/

/*... NEAR-WALL_MODEL*/
  #define STANDARDWALL 1
  #define ENHANCEDWALL 2
/*...................................................................*/

/*...*/
  #define NWALLPAR     4
/*...................................................................*/

/*...*/
  #define LDYNAMIC      1
  #define GDYNAMIC      2
  #define GDYNAMICMOD   3
  #define TWOPARDYNAMIC 4 
/*...................................................................*/

/*...*/
  void turbulence(Memoria *m            
                , Loads *lKturb  , Loads *lVel
                , PartMesh *pMesh, Turbulence *tModel
                , Mesh *mesh     , Scheme *sc
                , Simple *sp     , SistEq *sistEq    
                , Solv *solv     , short const ndf);   
/*...................................................................*/

/*...*/
  void lesDynamicMean(Memoria *m              , Turbulence *tModel
                   , INT    *RESTRICT nelcon  , short  *RESTRICT nen    
                   , short  *RESTRICT nFace   , DOUBLE *RESTRICT gVolume 
                   , DOUBLE *RESTRICT vel     , DOUBLE *RESTRICT gradVel 
                   , DOUBLE *RESTRICT density , DOUBLE *RESTRICT dViscosity
                   , DOUBLE *RESTRICT cDyn                                   
                   , short const maxNo        , short const maxViz
                   , short const ndm          , INT const numel     
                   , short const ndf          , short const iCod) ; 
/*...................................................................*/

/*...*/
  void turbulenceCellLoop(Loads *lVel       , Turbulence *tModel             
      , INT    *RESTRICT el                 , INT    *RESTRICT nelcon 
      , short  *RESTRICT nen                , short  *RESTRICT nFace 
      , short  *RESTRICT geomType           , DOUBLE *RESTRICT prop  
      , short  *RESTRICT calType            , short  *RESTRICT mat        
      , DOUBLE *RESTRICT cc                 , DOUBLE *RESTRICT gKsi                      
      , DOUBLE *RESTRICT gmKsi              , DOUBLE *RESTRICT gEta  
      , DOUBLE *RESTRICT gfArea             , DOUBLE *RESTRICT gNormal  
      , DOUBLE *RESTRICT gVolume            , DOUBLE *RESTRICT gXm  
      , DOUBLE *RESTRICT gXmcc              , DOUBLE *RESTRICT gvSkew    
      , DOUBLE *RESTRICT gmvSkew            , DOUBLE *RESTRICT gDcca  
      , short  *RESTRICT faceVelR           , short *RESTRICT faceVelL               
      , DOUBLE *RESTRICT vel                , DOUBLE *RESTRICT gradVel
      , DOUBLE *RESTRICT density            , DOUBLE *RESTRICT dViscosity        
      , DOUBLE *RESTRICT eddyViscosity      , DOUBLE *RESTRICT wallPar  
      , DOUBLE *RESTRICT dynamic                 
      , short const maxNo                   , short const maxViz
      , short const ndm                     , INT const numel     
      , short const ndf);  
/*...................................................................*/


/*...*/
  void sLesModCellLoop(Turbulence *tModel      
                   , DOUBLE *RESTRICT x       , INT *RESTRICT el       
                   , INT *RESTRICT nelcon     , short  *RESTRICT nen    
                   , short *RESTRICT nFace    , DOUBLE *RESTRICT gVolume 
                   , DOUBLE *RESTRICT nVel    , DOUBLE *RESTRICT eVel 
                   , DOUBLE *RESTRICT nDensity, DOUBLE *RESTRICT eDensity  
                   , DOUBLE *RESTRICT gradVel , DOUBLE *RESTRICT stressR                                  
                   , DOUBLE *RESTRICT cd
                   , short const maxNo        , short const maxViz
                   , short const ndm          , INT const numel     
                   , short const ndf          , bool const fNode);  
/*...................................................................*/

/*...*/
  void cellLes(Loads *loadsVel         , Turbulence tModel           
          , short *RESTRICT lGeomType  , DOUBLE *RESTRICT prop
          , INT *RESTRICT lViz         , DOUBLE *RESTRICT ksi
          , DOUBLE *RESTRICT mKsi
          , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT mEta
          , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT area
          , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc
          , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT cc
          , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew
          , short  *RESTRICT faceVelR  , short *RESTRICT faceVelL      
          , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT gradVel   
          , DOUBLE *RESTRICT lDensity  , DOUBLE const dViscosity
          , DOUBLE *viscosity          
          , const short nEn            , short const nFace
          , const short ndm            , INT const nel);

  void eddyViscosity3D(Loads *lVel             , Turbulence tModel           
          , short *RESTRICT lGeomType  , DOUBLE *RESTRICT prop 
          , INT *RESTRICT lViz         , DOUBLE *RESTRICT fArea
          , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume 
          , DOUBLE *RESTRICT dcca      , short *RESTRICT lFaceVelR  
          , short *RESTRICT lFaceVelL  , DOUBLE *RESTRICT vel       
          , DOUBLE *RESTRICT gradVel   , DOUBLE *RESTRICT lDensity  
          , DOUBLE const dViscosity    , DOUBLE *viscosity          
          , DOUBLE *RESTRICT wallPar   , DOUBLE const nDyn
          , const short nEn            , short const nFace 
          , const short ndm            , INT const nel); 
/*...................................................................*/

/*...*/
  void lesDynOnePar(INT *RESTRICT lViz       , DOUBLE *RESTRICT volume
                , DOUBLE *RESTRICT lDensity  , DOUBLE *RESTRICT vel
                , DOUBLE *RESTRICT gradVel   , DOUBLE *RESTRICT lDynamic
                , short const typeLesFunc    , short const nFace  );
  
  void lesDynOneParG(INT *RESTRICT lViz      , DOUBLE *RESTRICT volume
                , DOUBLE *RESTRICT gradVel , DOUBLE *RESTRICT lDynamic
                , short const typeLesFunc  , short const nFace  );
  
  void lesDynTwoPar(INT *RESTRICT lViz       , DOUBLE *RESTRICT volume
                , DOUBLE *RESTRICT lDensity, DOUBLE *RESTRICT vel
                , DOUBLE *RESTRICT gradVel , DOUBLE *RESTRICT lDynamic
                , short const nFace  );
/*...................................................................*/

/*...*/
  void wallModel(DOUBLE const vt     , DOUBLE const viscosity
               , DOUBLE const density, DOUBLE const dWall 
               , DOUBLE *yP          , DOUBLE *uP
               , short const iCod);
  DOUBLE wallModelHeat(DOUBLE const yPlus,DOUBLE const prM
                     , DOUBLE const prT);
/*...................................................................*/

/*...*/
  void wallParametersCellLoop(Loads *lVel , Turbulence *tModel             
      , INT    *RESTRICT nelcon         , short  *RESTRICT nFace 
      , DOUBLE *RESTRICT gNormal        , DOUBLE *RESTRICT gVolume
      , DOUBLE *RESTRICT gDcca  
      , short  *RESTRICT faceVelR       , short *RESTRICT faceVelL               
      , DOUBLE *RESTRICT vel            , DOUBLE *RESTRICT density
      , DOUBLE *RESTRICT dViscosity     , DOUBLE *RESTRICT wallPar  
      , short const maxNo               , short const maxViz
      , short const ndm                 , INT const numel     
      , short const ndf);                      
/*...................................................................*/

/*...*/
  void oneEquationK(Memoria *m        
                , Loads *lKturb     , Loads *lVel
                , PartMesh *pMesh   , Turbulence *tModel
                , Mesh *mesh        , Scheme *sc
                , Simple *sp  
                , SistEq *sistEq    , Solv *solv);
  
  void lesOneKeqLoc(INT *RESTRICT lViz     , DOUBLE *RESTRICT volume
                , DOUBLE *RESTRICT lDensity, DOUBLE *RESTRICT vel
                , DOUBLE *RESTRICT gradVel , DOUBLE *RESTRICT lDynamic
                , DOUBLE const dViscosity  , DOUBLE const eddyViscosity
                , INT const nel            , short const nFace  );

  void lesOneKeqPar(Turbulence *tModel, INT *RESTRICT lViz
         , DOUBLE *RESTRICT dynamic , DOUBLE *RESTRICT cDyn
         , DOUBLE *RESTRICT volume  , short const nFace
         , INT const nel            , short const nPar 
         , bool const fCap); 

  void lesDynamicOneKeq(Memoria *m             , Turbulence *tModel
                  , INT    *RESTRICT nelcon  , short  *RESTRICT nen    
                  , short  *RESTRICT nFace   , DOUBLE *RESTRICT gVolume 
                  , DOUBLE *RESTRICT vel     , DOUBLE *RESTRICT gradVel 
                  , DOUBLE *RESTRICT density , DOUBLE *RESTRICT dViscosity
                  , DOUBLE *RESTRICT eddyVis , DOUBLE *RESTRICT cDyn                                      
                  , short const maxNo        , short const maxViz
                  , short const ndm          , INT const numel     
                  , short const ndf          , short const iCod) ;
/*...................................................................*/

/*...*/
  bool wallParameters(Loads *lVel               
             , INT *RESTRICT lViz       , DOUBLE *RESTRICT v 
             , DOUBLE *RESTRICT normal  , DOUBLE *RESTRICT dcca
             , short *RESTRICT lFaceVelR, short *RESTRICT lFaceVelL
             , DOUBLE const viscosity   , DOUBLE const density
             , DOUBLE *RESTRICT wallPar , DOUBLE *dWall     
             , short const wallType     , short const nFace );
/*...................................................................*/

/*...*/
  DOUBLE sigmaModel( DOUBLE *RESTRICT gradVel 
                   , short const nFace , short const ndm);
  DOUBLE waleModel(DOUBLE *RESTRICT gradVel);
  DOUBLE svs(DOUBLE *RESTRICT gradVel);
/*...................................................................*/

/*...*/
  void bardinaModel(DOUBLE *RESTRICT xl   , DOUBLE *RESTRICT stressR  
                , DOUBLE *RESTRICT nVel , DOUBLE *RESTRICT nDensity
                , DOUBLE const vol
                , DOUBLE const cs       , INT const nEl
                , short const ndm       , short const nFace);

  void bardinaModelMod(INT *RESTRICT lViz  , DOUBLE *RESTRICT stressR  
                      , DOUBLE *RESTRICT vel, DOUBLE *RESTRICT gradVel
                      , DOUBLE *RESTRICT vol, DOUBLE *RESTRICT lDensity
                      , DOUBLE const cs 
                      , short const ndm     , short const nFace);

  void structuralStress(Turbulence tModel , DOUBLE *RESTRICT xl
                , INT *RESTRICT lViz      , DOUBLE *RESTRICT stressR  
                , DOUBLE *RESTRICT nVel   , DOUBLE *RESTRICT eVel
                , DOUBLE *RESTRICT nDen   , DOUBLE *RESTRICT eDen
                , DOUBLE *RESTRICT gradVel, DOUBLE *RESTRICT vol 
                , DOUBLE const cd                   
                , short const ndm         , short const nFace                
                , INT const nEl);
/*...................................................................*/

/*...*/
  DOUBLE doubleDot(DOUBLE *t);
  DOUBLE doubleDotSym(DOUBLE *t);
  DOUBLE doubleDotSym2(DOUBLE *RESTRICT t,DOUBLE *RESTRICT q);
/*...................................................................*/

/*...*/
  DOUBLE oneParLes(Turbulence *tModel
                 , INT *RESTRICT lViz, DOUBLE *RESTRICT dynamic
                 , short const nFace , bool const cap);
  DOUBLE oneGoParLes(Turbulence *tModel , DOUBLE *gmm
                    , bool const cap);
  DOUBLE oneGoParLesMod(Turbulence *tModel  , DOUBLE *gmm
                    , DOUBLE const density, DOUBLE const viscosity
                    , bool const cap);
  void twoParLes(Turbulence *tModel      , INT *RESTRICT lViz
               , DOUBLE *RESTRICT dynamic, DOUBLE *RESTRICT cDyn             
               , short const nFace       , bool const cap);
  void tensorS( DOUBLE *RESTRICT s,  DOUBLE *RESTRICT gradVel
              , bool const flag );
/*...................................................................*/

/*...*/
  DOUBLE vremanModel(DOUBLE *RESTRICT gradVel);
  DOUBLE waleModel(DOUBLE *RESTRICT gradVel);
/*...................................................................*/

/*...*/
  void sfHexa8(DOUBLE const eps ,DOUBLE const nn
           , DOUBLE const ze  
           , DOUBLE *RESTRICT N        ,DOUBLE *RESTRICT de
           , DOUBLE *RESTRICT dn       ,DOUBLE *RESTRICT dz
           , bool const ninter         ,bool const dev);


  DOUBLE jacob3d(DOUBLE *RESTRICT xl, DOUBLE *RESTRICT de
           , DOUBLE *RESTRICT dn  , DOUBLE *RESTRICT dz 
           , DOUBLE *RESTRICT hx  , DOUBLE *RESTRICT hy
           , DOUBLE *RESTRICT hz   
 	         , short const nen      , bool const afl
           , bool const dev       , INT const nel );
/*...................................................................*/

#endif /*_TURBULENCE_H_*/