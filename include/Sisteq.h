#ifndef _SISTEQ_H_
  #define _SISTEQ_H_
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
/*...*/
  #include<Csr.h>
  #include<Coo.h>
  #include<EllPack.h>
  #include<Erro.h>
  #include<HccaStdBool.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<ParallelMpi.h>

/*... bufferEq*/
  typedef struct {
    INT *thBegin, *thEnd, *thHeight, *thSize;
    DOUBLE *thY;
  }BufferOmp;
/*...................................................................*/
   
  typedef struct{
    INT *ja,*ia;
    double *al,*ad,*au;
    double  *b,*b0,*x;
    INT *id;       /*numeracao da esquacoes por celula*/
    bool unsym;    /*matriz nao simetrica*/
    INT neq;       /*numero de eq*/
    INT neqNov;    /*numero de eq em elementos sem sobreposicao*/
    INT nad;       /*numero de nao zeros fora da diagonal principal*/
    INT nadr;      /*numero de não zeros na parte retangular da matriz*/
    INT bandCsr[3]; /*numero de não zeros na parte retangular da matriz*/
    short storage; /*tecnica de armazenamenro*/
                   /*1 - csr*/
                   /*2 - csrd*/
                   /*3 - csrc*/
                   /*4 - ellpack*/
    Interface iNeq;/*comunicao MPI*/
    BufferOmp omp; /*OPenMp*/
  }SistEq;
  
  INT numeq(INT  *RESTRICT id  ,INT *RESTRICT num
          ,short *RESTRICT rt  ,short *RESTRICT nen
          ,INT const numel     ,short const nViz
          ,short const ndf);
  
  INT numEqV1(INT  *RESTRICT id  ,INT *RESTRICT num
             ,INT const numel);                       
  
  INT numEqV2(INT  *RESTRICT id  ,INT *RESTRICT num
             ,short *RESTRICT rt  ,short *RESTRICT nen
             ,INT const numel     ,short const nViz);
  
  INT countEq(INT *RESTRICT num
           ,short *RESTRICT rt     ,short *RESTRICT nFace 
           ,INT const numel        ,short const maxViz  
           ,short const ndf);

  void dataStruct(Memoria *m      ,INT *id  
                 ,INT *num        ,INT *nelcon
                 ,short *nViz 
                 ,INT const numel ,short const maxViz
                 ,short const ndf
                 ,char  *strIa    ,char *strJa
                 ,char *strAd     ,char *strA  
                 ,SistEq *SistEqX);

  void dataStructSimple(Memoria *m     ,INT *id
                       ,INT *num       ,INT *nelcon
                       ,short *nViz
                       ,INT const numel,short const maxViz
                       ,short const ndf
                       ,char  *strIa   ,char *strJa
                       ,char *strAd    ,char *strA
                       ,SistEq *sistEqX);
  
  void setDataStruct(char *word,short *data);

  void front(Memoria *m
          ,PartMesh *pMesh, SistEq *sistEq, short const ndf);


/*  void datastruct(Memoria *,int *,INT *,INT,INT,short,short
                 ,Sisteq*);*/
/*SKYLINE*/
/*  void profil(int*,INT*,int*,INT,INT,short,short,INT
             ,INT*);*/

#endif/*_SISTEQ_H*/
