#ifndef _CELLLOOP_H
  #define _CELLLOOP_H
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<Assbly.h>
  #include<HccaBlas.h>
  #include<HccaStdBool.h>
  #include<Define.h>
/*...................................................................*/

/* ... calculo das propriedade geometicas*/
  void pGeomForm(double *restrict x    ,INT    *restrict el
              ,INT    *restrict nelcon ,short  *restrict nen    
              ,short  *restrict nFace  ,short  *restrict geomType
              ,double *restrict gCc    ,double *restrict gKsi  
              ,double *restrict gmKsi  ,double *restrict gEta 
              ,double *restrict gmEta  ,double *restrict gNormal
              ,double *restrict gVolume,double *restrict gXm    
              ,double *restrict gXmcc  ,double *restrict gmKm   
              ,double *restrict gDcca               
              ,short maxNo             ,short maxViz
              ,short ndm               ,INT numel);
/*...................................................................*/

/* ... montagem do sistemas de equacoes*/
 void systForm(double *restrict x      ,INT    *restrict el
             ,INT    *restrict nelcon ,short  *restrict nen    
             ,short  *restrict nFace  ,short  *restrict geomType
             ,double *restrict prop   ,short  *restrict calType
             ,short  *restrict mat    ,double *restrict gCc
             ,double *restrict gKsi   ,double *restrict gmKsi
             ,double *restrict gEta   ,double *restrict gmEta
             ,double *restrict gNormal,double *restrict gVolume
             ,double *restrict gXm    ,double *restrict gXmcc   
             ,double *restrict gmKm   ,double *restrict gDcca                 
             ,INT    *restrict ia     ,INT    *restrict ja                    
             ,double *restrict ad     ,double *restrict al                    
             ,double *restrict b      ,INT    *restrict id 
             ,short  *restrict faceR  ,double *restrict faceS 
             ,double *restrict u0                                  
             ,short const maxNo       ,short const maxViz
             ,short const ndm         ,INT const numel
             ,short const ndf         ,short const storage
             ,bool const forces       ,bool const matrix 
             ,bool const unsym); 
/*...................................................................*/

/*... chamada da biblioteca de elementos*/
  void cellLib(double *restrict lx                     
              ,short *restrict lGeomType,double *restrict lprop
              ,INT   *restrict lViz     ,double *restrict xc                             
              ,double *restrict ksi     ,double *restrict mksi
              ,double *restrict eta     ,double *restrict meta
              ,double *restrict normal  ,double *restrict volume
              ,double *restrict xm      ,double *restrict xmcc
              ,double *restrict dcca    ,double *restrict mkm
              ,double *restrict lA      ,double *restrict lB
              ,short  *restrict lFaceR  ,double *restrict lFaceS
              ,double *restrict u0      ,short const nen     
              ,short const nFace        ,short const ndm
              ,short const lib          ,INT const nel);
/*...................................................................*/

/*... carga por elmento e condicoes pescritas por celula*/
  void cellPload(short  *restrict faceR ,double *restrict faceS
                ,double *restrict volume
                ,double *restrict u     ,double *restrict f
                ,INT const numel        ,short const ndf
                ,short const maxViz);
/*...................................................................*/

/*... funcoes geometricas*/
  void cellGeom2D(double *restrict lx       ,short *restrict lnFace
                 ,short  *restrict lGeomType,double *restrict xc
                 ,double *restrict ksi      ,double *restrict mksi
                 ,double *restrict eta      ,double *restrict meta
                 ,double *restrict normal   ,double *restrict volume
                 ,double *restrict xm       ,double *restrict xmcc
                 ,double *restrict dcca     ,double *restrict mkm
                 ,short  *restrict sn       ,short maxNo   
                 ,short maxViz              ,short ndm
                 ,INT nel);
/*...................................................................*/

/*biblioteca de celulas*/
  void cellDif2D(double *restrict lx      
                ,short *restrict lGeomType,double *restrict lprop
                ,INT   *restrict lViz     
                ,double *restrict xc                             
                ,double *restrict ksi     ,double *restrict mksi
                ,double *restrict eta     ,double *restrict meta
                ,double *restrict normal  ,double *restrict volume
                ,double *restrict xm      ,double *restrict xmcc
                ,double *restrict dcca    ,double *restrict mkm
                ,double *restrict lA      ,double *restrict lB
                ,short  *restrict lFaceR  ,double *restrict lFaceS
                ,double *restrict u0      ,short const nen     
                ,short const nFace        ,short const ndm
                ,INT const nel);
/*...................................................................*/


/*... funcoes de apoio*/  
  void sn(short *s,short ty, INT nel);
  double areaQuadCell(double *restrict eta,short ndm);
  double areaTriaCell(double *restrict eta,short ndm);
  double volumeCell(double *eta,short ty,short ndm,INT nel);
  void vectorKm2d(double *restrict x ,double *restrict xc
                 ,double *restrict xm,double *restrict mkm
                 ,short  *restrict sn,short nFace
                 ,short maxViz       ,short  maxNo       
                 ,short ndm          ,INT nel);
/*...................................................................*/

#endif
