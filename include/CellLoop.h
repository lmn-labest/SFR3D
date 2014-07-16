#ifndef _CELLLOOP_
  #define _CELLLOOP_
  #include<Define.h>
  #include<stdio.h>
  #include<stdlib.h>
  #include<math.h>
  #include<HccaBlas.h>
  typedef enum Lib{geometria  = 1
                  ,gradiente
                  ,matriz
                  }Lib; 

  void pGeomForm(double *restrict x      ,INT    *restrict el
              ,INT    *restrict nelcon ,short  *restrict nen    
              ,short  *restrict nFace  ,short  *restrict geomType
              ,double *restrict gCc
              ,double *restrict gKsi   ,double *restrict gmKsi
              ,double *restrict gEta   ,double *restrict gmEta
              ,double *restrict gNormal,double *restrict gVolume
              ,double *restrict gXm    ,double *restrict gXmcc   
              ,double *restrict gmKm   ,double *restrict gDcca               
              ,short maxNo             ,short maxViz
              ,short ndm               ,INT numel);

/*chamada da biblioteca de elementos*/
  void cellLib(double *lx,short *lnFace, short *lGeomType
              ,double *xc  
              ,short maxNo, short maxViz,short ndm
              ,short lib,INT nel);

/*funcoes geometricas*/
  void cellGeom2D(double *restrict lx       ,short *restrict lnFace
                 ,short  *restrict lGeomType,double *restrict xc
                 ,double *restrict ksi      ,double *restrict mksi
                 ,double *restrict eta      ,double *restrict meta
                 ,double *restrict normal   ,double *restrict volume
                 ,double *restrict xm       ,double *restrict xmcc
                 ,double *restrict dcca     ,double *restrict mkm
                 ,short  *restrict sn
                 ,short maxNo               ,short maxViz
                 ,short ndm                 ,INT nel);

/*funcoes de apoio*/  
  void sn(short *s,short ty, INT nel);
  double areaQuadCell(double *restrict eta,short ndm);
  double areaTriaCell(double *restrict eta,short ndm);
  double volumeCell(double *eta,short ty,short ndm,INT nel);
  void vectorKm2d(double *restrict x ,double *restrict xc
                 ,double *restrict xm,double *restrict mkm
                 ,short  *restrict sn,short nFace
                 ,short maxViz       ,short  maxNo       
                 ,short ndm          ,INT nel);

#endif
