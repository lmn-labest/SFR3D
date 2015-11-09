#ifndef _WRITE_VTK_H
  #define _WRITE_VTK_H
  #define INTEGER 1
  #define DOUBLEV 2
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
  void wGeoVtk(Memoria *m     ,double *x      
              ,INT *el        ,short *mat    
              ,short *nen     ,short *typeGeom
              ,double *prop   ,short *typeCal
              ,short *faceRd1 ,short *faceSd1
              ,short *faceRt1 ,short *faceSt1
              ,INT nnode      ,INT numel    
              ,short ndm      
              ,short maxno    ,short maxIt 
              ,short numat    
              ,short *ndfD    ,short *ndfT   
              ,char *nameOut  ,bool iws      
              ,FILE *f);

  void wGeoFaceVtk(Memoria *m     ,DOUBLE *x      
            ,INT *el              ,short *nen     
            ,short *typeGeom
            ,short *faceRd1       ,short *faceSd1
            ,short *faceRt1       ,short *faceSt1
            ,INT const nnode      ,INT const numel    
            ,short const ndm      
            ,short const ndfD     ,short const ndfT 
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
/*...................................................................*/

/*...*/
void makeFace(INT *el            ,short *faceR       ,short *faceL     
             ,short *typeGeom
             ,INT *face          ,INT *lFaceS    ,INT *idFace
             ,short *typeGeomFace,short *nenFace
             ,short const maxViz ,short const maxNo
             ,short const ndf 
             ,int   const numel  ,int *nFace);
/*...................................................................*/

#endif/*_WRITE_VTK_H*/
