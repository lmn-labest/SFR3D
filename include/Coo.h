#ifndef _COO_H  
  #define _COO_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
/*...................................................................*/ 

/*...*/
  #include<Memoria.h>
  #include<Define.h>
  #include<Graph.h>
  #include<HccaStdBool.h>
  #include<File.h>
  #ifdef _MMIO_
    #include<mmio/mmio.h>
  #endif
/*...................................................................*/ 
  void sortCoo(INT *RESTRICT ia0
              ,INT *RESTRICT ia ,INT *RESTRICT ja
              ,INT const nnz    ,INT const nEq);

  void cooIaJaR(INT *RESTRICT ia0
           ,INT *RESTRICT ia     ,INT *RESTRICT ja
           ,INT *RESTRICT id     ,INT *RESTRICT num   
           ,INT *RESTRICT adj    ,short *RESTRICT nViz
           ,INT const numel      ,INT const nEq
           ,short const maxViz   ,short  const ndf);

  INT cooNnzR(INT *RESTRICT id   
            ,INT *RESTRICT num    ,INT *RESTRICT adj
            ,short *RESTRICT nViz
            ,INT const numel      ,INT const nEq
            ,short const maxViz   ,short  const ndf);


/*...*/
  void writeCoo(Memoria *m,INT *ia   ,INT *ja,INT neq
               ,double *au,double*ad ,double *al 
               ,INT nad   ,short type
               ,bool unsym,bool bin
               ,char* name);
  void writeCooB(double *b,INT const neq,char *name);
/*...................................................................*/ 

/*...*/
  void csrToCoo(int  *linha,int  *col ,double *val
               ,INT *ia    ,INT *ja   ,double *au
               ,double *ad ,double *al,INT neq  
               ,INT nad    ,bool unsym,bool bin);
  void cooToFull(int *lin, int *col,double *val,double *a,INT neq
                ,INT nad);
  void cooToCsr(int *lin  , int *col ,double *val
               ,INT *ia   ,INT *ja
               ,double *au,double *ad,double *al
               ,INT neq   ,INT nad   ,short type  
               ,int *aux  ,bool upper,bool diag ,bool lower);
/*...................................................................*/ 

/*...*/
  void matrixCheck(double *val,int *lin,int *col, int nl,int nnz);
/*...................................................................*/ 


#endif/*_COO_H*/
