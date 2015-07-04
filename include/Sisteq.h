#ifndef _SISTEQ_H
  #define _SISTEQ_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
  #include<stdio.h>
  #include<stdlib.h>
  #include<HccaStdBool.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<Csr.h>

  
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
  }SistEq;
  
  INT numeq(Memoria *m ,INT *id  ,INT *num, short *rt
            , short *nen,INT numel,short nViz,short ndf);

  void dataStruct(Memoria *m  ,INT *id  ,INT *num   ,INT *nelcon
                 ,short *nViz ,INT numel,short maxViz,short ndf
                 ,char  *strIa,char *strJa,char *strAd ,char *strA  
                 ,SistEq *SistEqX);

/*  void datastruct(Memoria *,int *,INT *,INT,INT,short,short
                 ,Sisteq*);*/
/*SKYLINE*/
/*  void profil(int*,INT*,int*,INT,INT,short,short,INT
             ,INT*);*/

#endif/*_SISTEQ_H*/
