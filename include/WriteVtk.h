#ifndef _WRITE_VTK_H_
  #define _WRITE_VTK_H_
  

/*...*/
  #include<stdio.h>
  #include<stdlib.h>
/*...*/
  #include<HccaStdBool.h>
  #include<CellLoop.h>
  #include<Combustion.h>
  #include<Vtk.h>
  #include<Mesh.h>
  #include<Memoria.h>
  #include<File.h>
  #include<Properties.h>
  #include<Erro.h>
/*... particionamento*/
  void wPartVtk(Memoria *m   ,double *x      
               ,INT *el            
               ,short *nen   ,short *typeGeom
               ,INT   *np    ,INT   *ep      
               ,INT nnode    ,INT numel    
               ,short ndm      
               ,short maxNo  ,short maxViz  
               ,char *nameOut,bool iws
               ,FILE *f);

  void wMeshPartVtk(Memoria *m     
            ,double *x      ,INT *el            
            ,short *nen     ,short *typeGeom
            ,INT nnode      ,INT numel    
            ,short ndm      
            ,short maxNo    ,short maxViz  
            ,char *nameOut  ,bool iws
            ,FILE *f);
/*...................................................................*/

/*... geom*/  
  void wGeoVtk(Memoria *m        ,double *x  
              ,INT *el           ,short *mat    
              ,short *nen        ,short *typeGeom
              ,double *prop      ,short *typeCal
              ,short *faceRd1    ,short *faceSd1
              ,short *faceRt1    ,short *faceSt1
              ,short *faceRfluid ,short *faceSfluid
              ,short *faceRenergy,short *faceLenergy 
              ,INT nnode         ,INT numel    
              ,short ndm      
              ,short maxno       ,short maxIt 
              ,short numat    
              ,short *ndfD       ,short *ndfT   
              ,short const ndfF  ,short const ndfFt
              ,char *nameOut     ,bool iws      
              ,FILE *f);

  void wGeoVtk2(Memoria *m   , DOUBLE *x
              , DOUBLE *cc   , INT *el
              , short *nen   , short *typeGeom
              , INT nnode    , INT numel
              , short ndm    , short maxNo
              , short maxViz
              , char *nameOut, bool iws
              , FILE *f);

  void wGeoFaceVtk(Memoria *m     ,DOUBLE *x      
            ,INT *el              ,short *nen     
            ,short *typeGeom
            ,short *faceRd1       ,short *faceSd1
            ,short *faceRt1       ,short *faceSt1
            ,short *faceRfluid    ,short *faceSfluid
            ,short *faceRtemp     ,short *faceLtemp
            ,INT const nnode      ,INT const numel    
            ,short const ndm      
            ,short const ndfD     ,short const ndfT 
            ,short const ndfF     ,short const ndfFt
            ,short const maxViz   ,short const maxNo
            ,char *nameOut        ,bool iws
            ,FILE *f);

  void wGeoFaceVtk2(Memoria *m          , DOUBLE *x
                  , INT *el             , short *nen
                  , short *typeGeom     
                  , short *faceRd       , short *faceLd
                  , INT const nnode     , INT const numel
                  , short const ndm     , short const maxViz
                  , short const ndf     , short const maxNo
                  , char *nameOut       , bool iws
                  , bool const fWallVel , FILE *f);
/*...................................................................*/

/*... resultados*/  
  void wResVtk(Memoria *m     ,double *x      
            ,INT *el        ,short *mat    
            ,short *nen     ,short *typeGeom
            ,DOUBLE *elU    ,DOUBLE *nU
            ,INT nnode      ,INT numel    
            ,short ndm      ,short maxNo 
            ,short numat    ,short *ndfD   
            ,char *nameOut  ,bool iws
            ,FILE *f);

  void wResVtkDif(Memoria *m        , double *x
                , INT *el           , short *mat
                , short *nen        , short *typeGeom
                , DOUBLE *elU       , DOUBLE *nU
                , DOUBLE *elGradU   , DOUBLE *nGradU
                , DOUBLE *elDensity , DOUBLE *nDensity
                , DOUBLE *elCoefDiff, DOUBLE *nCoefDiff
                , INT nnode         , INT numel
                , short ndm         , short maxNo
                , short numat       , short ndf
                , char *uResEl      , char *uResNo
                , char *gradResEl   , char *gradResNo
                , char *nameOut     , FileOpt *opt
                , Temporal *ddt     , FILE *f);
  
  void wResVtkTrans(Memoria *m        , double *x
                  , INT *el           , short *mat
                  , short *nen        , short *typeGeom
                  , DOUBLE *elU       , DOUBLE *nU
                  , DOUBLE *elGradU   , DOUBLE *nGradU
                  , DOUBLE *elVel     , DOUBLE *nVel
                  , DOUBLE *elDensity , DOUBLE *nDensity
                  , DOUBLE *elCoefDiff, DOUBLE *nCoefDiff
                  , INT nnode         , INT numel
                  , short ndm         , short maxNo
                  , short numat       , short ndf
                  , char **ps    
                  , char *nameOut     , FileOpt *opt
                  , Temporal *ddt     , FILE *f);

  void wResVtkFluid(Memoria *m   , DOUBLE *x
          , DOUBLE *cc      
          , INT *el              , short *mat    
          , short *nen           , short *typeGeom
          , DOUBLE *elPres       , DOUBLE *nPres
          , DOUBLE *elGradPres   , DOUBLE *nGradPres
          , DOUBLE *elVel        , DOUBLE *nVel      
          , DOUBLE *elGradVel    , DOUBLE *nGradVel 
          , DOUBLE *elEnergy     , DOUBLE *nEnergy
          , DOUBLE *elGradEnergy , DOUBLE *nGradEnergy
          , DOUBLE *elEddyVis    , DOUBLE *nEddyVis
          , DOUBLE *eDensityFluid, DOUBLE *nDensityFluid
          , DOUBLE *eDyViscosity , DOUBLE *nDyViscosity
          , DOUBLE *eStressR     , DOUBLE *nStressR
          , DOUBLE *eCd          , DOUBLE *nCd 
          , DOUBLE *eWallPar     , DOUBLE *nWallPar 
          , DOUBLE *eKturb       , DOUBLE *nKturb
          , DOUBLE *eMedVel      , DOUBLE *nMedVel
          , DOUBLE *specificHeat , DOUBLE *tConductivity               
          , INT nnode            , INT numel    
          , short const ndm      , short const maxNo 
          , short const numat    , short const ndf  
          , short const ntn   
          , char *nameOut        , FileOpt opt
          , bool fKelvin         , Mean *media  
          , Temporal ddt         , FILE *f);
/*...................................................................*/

/*...*/
  void wResVtkCombustion(Memoria *m,Combustion *cModel     
          , DOUBLE *x            , DOUBLE *cc     
          , INT *el              , short *mat    
          , short *nen           , short *typeGeom
          , DOUBLE *elPres       , DOUBLE *nPres
          , DOUBLE *elGradPres   , DOUBLE *nGradPres
          , DOUBLE *elVel        , DOUBLE *nVel      
          , DOUBLE *elGradVel    , DOUBLE *nGradVel 
          , DOUBLE *elTemp       , DOUBLE *nTemp
          , DOUBLE *elGradEnergy , DOUBLE *nGradEnergy
          , DOUBLE *elZcomb      , DOUBLE *nZcomb  
          , DOUBLE *elGradZcomb  , DOUBLE *nGradZcomb 
          , DOUBLE *elEddyVis    , DOUBLE *nEddyVis
          , DOUBLE *eDensityFluid, DOUBLE *nDensityFluid
          , DOUBLE *eDyViscosity , DOUBLE *nDyViscosity
          , DOUBLE *eStressR     , DOUBLE *nStressR
          , DOUBLE *eCd          , DOUBLE *nCd
          , DOUBLE *eWallPar     , DOUBLE *nWallPar
          , DOUBLE *eKturb       , DOUBLE *nKturb
          , DOUBLE *eRateFuel    , DOUBLE *nRateFuel
          , DOUBLE *eYfrac       , DOUBLE *nYfrac 
          , DOUBLE *eHeatRe      , DOUBLE *nHeatRe     
          , DOUBLE *eMedVel      , DOUBLE *nMedVel
          , DOUBLE *specificHeat , DOUBLE *tConductivity
          , DOUBLE *cDiffSp
          , INT nnode            , INT numel    
          , short const ndm      , short const maxNo 
          , short const numat    , short const ndf
          , short const ntn      
          , char *nameOut        , FileOpt *opt
          , bool fKelvin         , Mean *media  
          , Temporal ddt         , FILE *f);
/*...................................................................*/

/*...*/
  void makeFace(INT *el            ,short *faceR       ,short *faceL 
             ,short *typeGeom
             ,INT *face          ,int    *lFaceL     ,INT *idFace
             ,short *typeGeomFace,short *nenFace
             ,short const maxViz ,short const maxNo
             ,short const ndf    ,INT const numel   
             ,INT *nFace         ,bool const fWallVel);

  void makeVorticity(DOUBLE *RESTRICT w, DOUBLE *RESTRICT gradVel
                    ,INT const n       , const short ndm);
  void makeStress(DOUBLE *RESTRICT str      , DOUBLE *RESTRICT gradVel
                 , DOUBLE *RESTRICT viscosity,INT const n          
                 , short const ndm           , short const ntn           
                 , bool const flag );
  void makeKineticEnergy(DOUBLE *RESTRICT e    , DOUBLE *RESTRICT vel
                    ,DOUBLE *RESTRICT density 
                    ,INT const n               , short const ndm);

  void makeQcriterion(DOUBLE *RESTRICT gradVel, DOUBLE *RESTRICT q
                    , INT const n             , short const ndm);
  void makePresTotal(DOUBLE *RESTRICT presT, DOUBLE *RESTRICT pres
                 , DOUBLE *RESTRICT vel  , DOUBLE *RESTRICT density 
                 , INT const n           , short const ndm);
  void makeModuleVel(DOUBLE *RESTRICT p,DOUBLE *RESTRICT elVel
                 , INT const n           , short const ndm);
/*...................................................................*/

#endif/*_WRITE_VTK_H_*/
