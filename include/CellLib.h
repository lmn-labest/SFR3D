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
      , DOUBLE *RESTRICT wallPar      , DOUBLE *RESTRICT dynamic
      , short const nEn               , short  const nFace
      , short const ndm               , short const lib
      , INT const nel);
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
  bool openDomain(Loads *loadVel
                 , short  *RESTRICT faceL, short  *RESTRICT nFace
                 , INT const numel      , short const maxViz);

  void wallFluid(short *RESTRICT faceR ,INT *RESTRICT nelcon
              ,short *RESTRICT nFace     
              ,INT const nEl         ,short const maxViz);
  
  DOUBLE totalMass(DOUBLE *RESTRICT density  , DOUBLE *RESTRICT volume
                  ,INT const nEl) ;

  DOUBLE massFluxOpenDomain(Loads *loadVel    , Temporal const ddt
              , short  *RESTRICT faceVelLoad, short  *RESTRICT nFace
              , DOUBLE *RESTRICT gfArea     , DOUBLE *RESTRICT gNormal
              , DOUBLE *RESTRICT density    , DOUBLE *RESTRICT vel 
              , INT const numel             , short const ndm  
              , short const maxViz  ) ;
 
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

  void difusionSchemeAnisotropic(DOUBLE *RESTRICT s,DOUBLE *RESTRICT ksi
                                ,DOUBLE *RESTRICT e,DOUBLE *RESTRICT t
                                ,short const ndm   ,short const iCod);
  void advectiveScheme(DOUBLE *RESTRICT velC       ,DOUBLE *RESTRICT velV
                     ,DOUBLE *RESTRICT gradVelC   ,DOUBLE *RESTRICT gradVelV
                     ,DOUBLE *RESTRICT gradVelComp,DOUBLE *RESTRICT vSkew
                     ,DOUBLE *RESTRICT rC         ,DOUBLE *RESTRICT rV
                     ,DOUBLE *RESTRICT ksi        ,DOUBLE const modKsi
                     ,DOUBLE const wfn            ,DOUBLE *RESTRICT cvc
                     ,DOUBLE const alphaMenosUm   ,DOUBLE const alpha
                     ,short const ndm
                     ,short const iCod1           ,short const iCod2);
   
  void advectiveSchemeScalar(DOUBLE const uC, DOUBLE const uV
              ,DOUBLE *RESTRICT gradUc    ,DOUBLE *RESTRICT gradUv
              ,DOUBLE *RESTRICT gradUcomp ,DOUBLE *RESTRICT vSkew
              ,DOUBLE *RESTRICT rC        ,DOUBLE *RESTRICT rV
              ,DOUBLE *RESTRICT ksi       ,DOUBLE const modKsi
              ,DOUBLE const wfn           ,DOUBLE *cvc
              ,DOUBLE const alphaMenosUm  ,DOUBLE const alpha          
              ,short const ndm
              ,short const iCod1          ,short const iCod2);

  DOUBLE deferredCd(DOUBLE const velC,DOUBLE const velV
                   ,DOUBLE const wfn);
  
  DOUBLE deferredLust(DOUBLE const uC         ,DOUBLE const uV
                   ,DOUBLE *RESTRICT gradUc   ,DOUBLE *RESTRICT gradUv
                   ,DOUBLE *RESTRICT rC       ,DOUBLE *RESTRICT rV
                   ,DOUBLE const alphaMenosUm ,DOUBLE const alpha               
                   ,DOUBLE const wfn          ,short const ndm);


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
  void vorticity(DOUBLE *RESTRICT w,DOUBLE *RESTRICT gradVel
                ,const short ndm);
  void stress(DOUBLE *RESTRICT s,DOUBLE *RESTRICT gradVel
           ,DOUBLE const nu   ,DOUBLE const lambda
           ,short const ndm);
/*...................................................................*/


#endif/*_CELLLIB_H_*/