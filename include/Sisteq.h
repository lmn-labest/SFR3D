#ifndef _SISTEQ_H
  #define _SISTEQ_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
/*...*/
  #include<Csr.h>
  #include<EllPack.h>
  #include<HccaStdBool.h>
  #include<Memoria.h>
  #include<Mesh.h>

  
  typedef struct{
    INT *ja,*ia;
    double *al,*ad,*au;
    double  *b,*b0,*x;
    INT *id;  /*numeracao da esquacoes por celula*/
    bool unsym;/*matriz nao simetrica*/
    INT neq; /*numero de eq*/
    INT nad;
    short storage; /*tecnica de armazenamenro*/
                   /*1 - csr*/
                   /*2 - csrd*/
                   /*3 - csrc*/
                   /*4 - ellpack*/
  }SistEq;
  
  INT numeq(INT  *restrict id  ,INT *restrict num
          ,short *restrict rt  ,short *restrict nen
          ,INT const numel     ,short const nViz
          ,short const ndf);

  void dataStruct(Memoria *m      ,INT *id  
                 ,INT *num        ,INT *nelcon
                 ,short *nViz 
                 ,INT const numel ,short const maxViz
                 ,short const ndf
                 ,char  *strIa    ,char *strJa
                 ,char *strAd     ,char *strA  
                 ,SistEq *SistEqX);
  
  void setDataStruct(char *word,short *data);

/*  void datastruct(Memoria *,int *,INT *,INT,INT,short,short
                 ,Sisteq*);*/
/*SKYLINE*/
/*  void profil(int*,INT*,int*,INT,INT,short,short,INT
             ,INT*);*/

#endif/*_SISTEQ_H*/
