#ifndef _SISTEQ_H
  #define _SISTEQ_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
  #include<stdio.h>
  #include<stdlib.h>
  #include<Mystdbool.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<Csr.h>

  #define min(a, b)  (((a) < (b)) ? (a) : (b))
  #define max(a, b)  (((a) > (b)) ? (a) : (b))

  typedef struct SistEq{
    long *ja,*ia;
    double *a,*ad,*b,*x,*y;
    long *id;  /*numeracao da esquacoes por celula*/
    bool unsym;/*matriz nao simetrica*/
    long neq; /*numero de eq*/
    long nad;
    short storage; /*tecnica de armazenamenro*/
                  /*1 - csr*/
  }SistEq;
  
  long numeq(Memoria *m ,long *id  ,long *num, short *rt
            , short *nen,long numel,short nViz, char *str);

  void dataStruct(Memoria *m  ,long *id  ,long *num   ,long *nelcon
                 ,short *nViz ,long numel,short maxViz,short ndf
                 ,char  *strIa,char *strJa,char *strAd ,char *strA  
                 ,SistEq *SistEqX);

/*  void datastruct(Memoria *,int *,long *,long,long,short,short
                 ,Sisteq*);*/
/*SKYLINE*/
/*  void profil(int*,long*,int*,long,long,short,short,long
             ,long*);*/

#endif/*_SISTEQ_H*/
