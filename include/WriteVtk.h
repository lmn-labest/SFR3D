#ifndef _WRITE_VTK_H_
  #define _WRITE_VTK_H_
  

/*...*/
  #include<stdio.h>
  #include<stdlib.h>
/*...*/
  #include<HccaStdBool.h>
  #include<CellLoop.h>
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

  void wResVtkDif(Memoria *m        ,double *x      
               ,INT *el           ,short *mat    
               ,short *nen        ,short *typeGeom
               ,DOUBLE *elU       ,DOUBLE *nU
               ,DOUBLE *elGradU   ,DOUBLE *nGradU
               ,INT nnode         ,INT numel    
               ,short ndm         ,short maxNo 
               ,short numat       ,short ndf   
               ,char *uResEl      ,char *uResNo 
               ,char *gradResEl   ,char *gradResNo 
               ,char *nameOut     ,bool iws
               ,Temporal ddt      ,FILE *f);
  
  void wResVtkTrans(Memoria *m        ,double *x      
               ,INT *el           ,short *mat    
               ,short *nen        ,short *typeGeom
               ,DOUBLE *elU       ,DOUBLE *nU
               ,DOUBLE *elGradU   ,DOUBLE *nGradU
               ,DOUBLE *elVel     ,DOUBLE *nVel  
               ,INT nnode         ,INT numel    
               ,short ndm         ,short maxNo 
               ,short numat       ,short ndf   
               ,char *uResEl      ,char *uResNo 
               ,char *gradResEl   ,char *gradResNo 
               ,char *velEl       ,char *velNo     
               ,char *nameOut     ,bool iws
               ,Temporal ddt      ,FILE *f);

  void wResVtkFluid(Memoria *m    ,DOUBLE *x      
          ,INT *el              ,short *mat    
          ,short *nen           ,short *typeGeom
          ,DOUBLE *elPres       ,DOUBLE *nPres
          ,DOUBLE *elGradPres   ,DOUBLE *nGradPres
          ,DOUBLE *elVel        ,DOUBLE *nVel      
          ,DOUBLE *elGradVel    ,DOUBLE *nGradVel 
          ,DOUBLE *elEnergy     ,DOUBLE *nEnergy
          ,DOUBLE *elGradEnergy ,DOUBLE *nGradEnergy
          ,DOUBLE *elEddyVis    ,DOUBLE *nEddyVis
          ,DOUBLE *eDensityFluid,DOUBLE *nDensityFluid
          ,DOUBLE *eDyViscosity ,DOUBLE *nDyViscosity
          ,DOUBLE *specificHeat ,DOUBLE *tConductivity
          ,DOUBLE *wallPar
          ,INT nnode            ,INT numel    
          ,short const ndm      ,short const maxNo 
          ,short const numat    ,short const ndf   
          ,char *nameOut        ,FileOpt opt
          ,bool fKelvin
          ,Temporal ddt         ,FILE *f);
/*...................................................................*/

/*...*/
  void makeFace(INT *el            ,short *faceR       ,short *faceL     
             ,short *typeGeom
             ,INT *face          ,INT *lFaceS    ,INT *idFace
             ,short *typeGeomFace,short *nenFace
             ,short const maxViz ,short const maxNo
             ,short const ndf 
             ,int   const numel  ,int *nFace);
  void makeVorticity(DOUBLE *RESTRICT w, DOUBLE *RESTRICT gradVel
                    ,INT const n       , const short ndm);
  void makeStress(DOUBLE *RESTRICT stress , DOUBLE *RESTRICT gradVel
               ,DOUBLE *RESTRICT viscosity 
               ,INT const n               , short const ndm);
  void makeKineticEnergy(DOUBLE *RESTRICT e    , DOUBLE *RESTRICT vel
                    ,DOUBLE *RESTRICT density 
                    ,INT const n               , short const ndm);
/*...................................................................*/

#endif/*_WRITE_VTK_H_*/
