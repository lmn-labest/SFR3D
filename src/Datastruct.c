#include<Sisteq.h>
/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 03/12/2017                                   *
 * ------------------------------------------------------------------*
 * DATASTRUC: estrutura de dados para a matriz de coeficientes       * 
 *-------------------------------------------------------------------* 
 * Parametros iniciais:                                              * 
 *-------------------------------------------------------------------* 
 * m       -> vetor de memoria                                       * 
 * id      -> numeracao das equacoes por celula                      * 
 * num     -> renumeracao dos elementos                              * 
 * nelcon  -> adjacencia dos elementos                               * 
 * nViz    -> numero de vizinhos por elemento                        * 
 * numel   -> numero de elementos                                    * 
 * maxViz  -> numero maximo de vizinhos                              * 
 * ndf     -> grau de liberdade                                      * 
 *-------------------------------------------------------------------* 
 * Parametros saida:                                                 * 
 *-------------------------------------------------------------------* 
 * sistEqX -> estrutura de matriz atualizado ia e ja.                * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 * ad -> alocados                                                    * 
 * al -> alocados                                                    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void dataStruct(Memoria *m      ,INT *id   
               ,INT *num        ,INT *nelcon
               ,short *nViz
               ,INT const numel ,short const maxViz
               ,short const ndf 
               ,char  *strIa    ,char *strJa
               ,char *strAd     ,char *strA  
               ,SistEq *sistEqX)
{
   short type = sistEqX->storage;
   bool lower,diag,upper,fCoo=false;
   INT n0,n1,n2,nad,nadr,nadT=0;
   
   switch(type)
   {
/*... armazenamento CSR(a)*/
     case CSR:
       diag = true;
/*...*/
       if(sistEqX->unsym){
         upper = true;
         lower = true;
       }
/*...*/
       else{
         upper = true;
         lower = false;
       }
/*...................................................................*/
       sistEqX->nadr = 0;
       n1 = sistEqX->neq + 1; 
       HccaAlloc(INT,m,sistEqX->ia,n1,strIa,_AD_);
/*... */
       nad = sistEqX->nad = csrIa(sistEqX->ia ,id 
                             ,num         ,nelcon 
                             ,nViz
                             ,numel       ,sistEqX->neqNov
                             ,maxViz      ,ndf
                             ,upper       ,diag         
                             ,lower       );
/*...................................................................*/

/*...*/
       HccaAlloc(INT,m,sistEqX->ja   ,nad ,strJa   ,_AD_);
       csrJa(sistEqX->ia ,sistEqX->ja,id ,num    
           ,nelcon,nViz ,numel      ,sistEqX->neqNov, maxViz,ndf
           ,upper,diag ,lower);
/*...................................................................*/

/*... reordenando o grafo*/
       sortGraphCsr(sistEqX->ia,sistEqX->ja,sistEqX->neqNov);
/*...................................................................*/

/*... alocacao da matriz*/
       HccaAlloc(DOUBLE         ,m    ,sistEqX->al
               ,nad            ,strA ,false);
       zero(sistEqX->al,nad,DOUBLEC);
       sistEqX->ad = sistEqX->al;
/*...................................................................*/

/*... banda da matriz*/
       sistEqX->bandCsr[BANDCSRMAX] 
                = bandCsr(sistEqX->ia
                         ,sistEqX->ja
                         ,sistEqX->neqNov
                         ,BANDCSRMAX);
       sistEqX->bandCsr[BANDCSRMED] 
               = bandCsr(sistEqX->ia
                        ,sistEqX->ja
                        ,sistEqX->neqNov
                        ,BANDCSRMED);
       
       sistEqX->bandCsr[BANDCSRMIN] 
               = bandCsr(sistEqX->ia
                        ,sistEqX->ja
                        ,sistEqX->neqNov
                        ,BANDCSRMIN);

       if(!mpiVar.myId  ) {
         fprintf(fileLogExc,"band Maxima: %d\n",sistEqX->bandCsr[BANDCSRMAX]);
         fprintf(fileLogExc,"band Media : %d\n",sistEqX->bandCsr[BANDCSRMED]);
         fprintf(fileLogExc,"band Minima: %d\n",sistEqX->bandCsr[BANDCSRMIN]);
       }
/*...................................................................*/
     break;
/*...................................................................*/

/*... armazenamento CSRD(a,ad)*/
/*...CSRD+COO(simetrico)*/
     case CSRDCOO:
       fCoo = true; 
     case CSRD:
       lower = true;
/*...*/
       if(sistEqX->unsym) 
         upper = true;
/*...*/
       else 
         upper = false;
/*...*/
       diag    = false;
/*...................................................................*/

/*... paralelo (CSRD + CSR - simetrico)*/
       if( mpiVar.nPrcs > 1 && !sistEqX->unsym && !fCoo) 
       { 
         n1 = sistEqX->neqNov + 1; 
         n2 = 2*n1; 
         HccaAlloc(INT,m,sistEqX->ia      ,n2  ,strIa   ,_AD_);
         zero(sistEqX->ia,n2    ,INTC);
/*... */
         nad = sistEqX->nad = csrIa(sistEqX->ia ,id 
                                   ,num         ,nelcon 
                                   ,nViz
                                   ,numel       ,sistEqX->neqNov
                                   ,maxViz      ,ndf
                                   ,upper       ,diag         
                                   ,lower     );
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
         nadr = sistEqX->nadr = csrIaR(&sistEqX->ia[n1],id 
                                      ,num             ,nelcon 
                                      ,nViz
                                      ,numel           ,sistEqX->neqNov
                                      ,maxViz          ,ndf);
/*...................................................................*/

/*...*/
         nadT = nad + nadr; 
         HccaAlloc(INT,m,sistEqX->ja     ,nadT ,strJa   ,_AD_);
         zero(sistEqX->ja,nadT  ,INTC);
         csrJa(sistEqX->ia ,sistEqX->ja 
              ,id          ,num     
              ,nelcon      ,nViz 
              ,numel       ,sistEqX->neqNov
              ,maxViz      ,ndf
              ,upper       ,diag 
              ,lower);
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
         csrJaR(&sistEqX->ia[n1] ,&sistEqX->ja[nad] 
               ,id               ,num     
               ,nelcon           ,nViz 
               ,numel            ,sistEqX->neqNov
               ,maxViz           ,ndf);
/*...................................................................*/

/*... reordenando o grafo*/
         sortGraphCsr(sistEqX->ia,sistEqX->ja,sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo*/
         sortGraphCsr(&sistEqX->ia[n1],&sistEqX->ja[nad]
                     ,sistEqX->neqNov);
/*...................................................................*/

/*... alocacao da matriz*/
         HccaAlloc(DOUBLE         ,m    ,sistEqX->ad
                  ,sistEqX->neqNov,strAd,false);
         HccaAlloc(DOUBLE         ,m    ,sistEqX->al
                  ,nadT           ,strA ,false);
         zero(sistEqX->ad,sistEqX->neqNov,DOUBLEC);
         zero(sistEqX->al,nadT,DOUBLEC);
/*...................................................................*/

       }
/*...................................................................*/

/*... paralelo (CSRD + COO - simetrico)*/
      else if(mpiVar.nPrcs > 1 && !sistEqX->unsym && fCoo) 
      {
/*... obtendo o numero de termos nao nulos na parte retangular*/ 
        nadr = sistEqX->nadr = cooNnzR(id     ,num   
                                      ,nelcon ,nViz 
                                      ,numel  ,sistEqX->neqNov
                                      ,maxViz ,ndf);
/*...................................................................*/

/*... alocando o vetor ia(2*neq+1+nadr)*/
        n0 = sistEqX->neqNov + 1; 
        n1 = n0 + sistEqX->neqNov + 1; 
        n2 = n1 + nadr; 
        HccaAlloc(INT,m,sistEqX->ia      ,n2  ,strIa   ,_AD_);
        zero(sistEqX->ia,n2    ,INTC);
/*... vetor ia do CSRD*/
        nad = sistEqX->nad = csrIa(sistEqX->ia ,id 
                                  ,num         ,nelcon 
                                  ,nViz
                                  ,numel       ,sistEqX->neqNov
                                  ,maxViz      ,ndf
                                  ,upper       ,diag         
                                  ,lower);
/*...................................................................*/

/*... vetor ja do CSRD*/
         nadT = nad + nadr; 
         HccaAlloc(INT,m,sistEqX->ja     ,nadT ,strJa   ,_AD_);
         zero(sistEqX->ja,nadT  ,INTC);
         csrJa(sistEqX->ia ,sistEqX->ja 
              ,id          ,num     
              ,nelcon      ,nViz 
              ,numel       ,sistEqX->neqNov
              ,maxViz      ,ndf
              ,upper       ,diag 
              ,lower);
/*...................................................................*/

/*... vetor ia[n0] parte retangular no formato CSR 
      (para a montagem da matriz global) 
      ia[n1] e ja da parte retangular no formato COO*/
         cooIaJaR(&sistEqX->ia[n0]
                 ,&sistEqX->ia[n1]     ,&sistEqX->ja[nad]          
                 ,id                   ,num
                 ,nelcon               ,nViz
                 ,numel                ,sistEqX->neqNov
                 ,maxViz               ,ndf);
/*...................................................................*/

/*... reordenando o grafo csrd ja(nad)*/
         sortGraphCsr(sistEqX->ia,sistEqX->ja,sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo Coo ja(nadr)*/
         sortGraphCsr(&sistEqX->ia[n0],&sistEqX->ja[nad]
                     ,sistEqX->neqNov);
/*...................................................................*/

/*... alocacao da matriz*/
         HccaAlloc(DOUBLE         ,m    ,sistEqX->ad
                  ,sistEqX->neqNov,strAd,false);
         HccaAlloc(DOUBLE         ,m    ,sistEqX->al
                  ,nadT           ,strA ,false);
         zero(sistEqX->ad,sistEqX->neqNov,DOUBLEC);
         zero(sistEqX->al,nadT,DOUBLEC);
/*...................................................................*/

      }
/*...................................................................*/
  
/*... sequencial (CSRD - simetrico e nao simetrico)*/
/*... paralelo   (CSRD - nao simetrico)*/
       else
       { 
         sistEqX->nadr = 0;
         n1 = sistEqX->neq + 1; 
         HccaAlloc(INT,m,sistEqX->ia,n1,strIa,_AD_);
/*... */
         nad = sistEqX->nad = csrIa(sistEqX->ia ,id 
                             ,num         ,nelcon 
                             ,nViz
                             ,numel       ,sistEqX->neqNov
                             ,maxViz      ,ndf
                             ,upper       ,diag         
                             ,lower       );
/*...................................................................*/

/*...*/
         HccaAlloc(INT,m,sistEqX->ja   ,nad ,strJa   ,_AD_);
         csrJa(sistEqX->ia ,sistEqX->ja,id ,num    
              ,nelcon,nViz ,numel      ,sistEqX->neqNov, maxViz,ndf
              ,upper,diag ,lower);
/*...................................................................*/

/*... reordenando o grafo*/
         sortGraphCsr(sistEqX->ia,sistEqX->ja,sistEqX->neqNov);
/*...................................................................*/

/*... alocacao da matriz*/
         HccaAlloc(DOUBLE         ,m    ,sistEqX->ad
                  ,sistEqX->neqNov,strAd,false);
         HccaAlloc(DOUBLE         ,m    ,sistEqX->al
                  ,nad            ,strA ,false);
         zero(sistEqX->ad,sistEqX->neqNov,DOUBLEC);
         zero(sistEqX->al,nad,DOUBLEC);
/*...................................................................*/

       }
/*...................................................................*/

/*... banda da matriz*/
       sistEqX->bandCsr[BANDCSRMAX] 
       =  bandCsr(sistEqX->ia
                        ,sistEqX->ja
                        ,sistEqX->neqNov
                        ,BANDCSRMAX);
       sistEqX->bandCsr[BANDCSRMED] 
       =  bandCsr(sistEqX->ia
                        ,sistEqX->ja
                        ,sistEqX->neqNov
                        ,BANDCSRMED);
       
      sistEqX->bandCsr[BANDCSRMIN] 
      =  bandCsr(sistEqX->ia
                        ,sistEqX->ja
                        ,sistEqX->neqNov
                        ,BANDCSRMIN);

       if(!mpiVar.myId  ) 
       {
         fprintf(fileLogExc,"band Maxima: %d\n",sistEqX->bandCsr[BANDCSRMAX]);
         fprintf(fileLogExc,"band Media : %d\n",sistEqX->bandCsr[BANDCSRMED]);
         fprintf(fileLogExc,"band Minima: %d\n",sistEqX->bandCsr[BANDCSRMIN]);
       }
/*...................................................................*/
     break;
/*...................................................................*/

/*... armazenamento CSRC(ad,au,al)*/
/*...CSRC+COO*/
     case CSRCCOO:
       fCoo = true; 
     case CSRC:
/*...*/
       lower = true;
/*...*/
       upper = false;
/*...*/
       diag  = false;
/*...................................................................*/

/*... paralelo (CSRC + CSR)*/
       if(mpiVar.nPrcs > 1 && !fCoo)
       {
         n1 = sistEqX->neqNov + 1; 
         n2 = 2*n1; 
         HccaAlloc(INT,m,sistEqX->ia      ,n2  ,strIa   ,_AD_);
         zero(sistEqX->ia,n2    ,INTC);

/*... */
         nad = sistEqX->nad = csrIa(sistEqX->ia ,id 
                                   ,num         ,nelcon 
                                   ,nViz
                                   ,numel       ,sistEqX->neqNov
                                   ,maxViz      ,ndf
                                   ,upper       ,diag         
                                   ,lower       );
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
         nadr = sistEqX->nadr = csrIaR(&sistEqX->ia[n1],id 
                                      ,num             ,nelcon 
                                      ,nViz
                                      ,numel           ,sistEqX->neqNov
                                      ,maxViz          ,ndf);
/*...................................................................*/

/*...*/
         nadT = nad + nadr;
         HccaAlloc(INT,m,sistEqX->ja     ,nadT ,strJa   ,_AD_);
         zero(sistEqX->ja,nadT  ,INTC);
         csrJa(sistEqX->ia ,sistEqX->ja 
              ,id          ,num     
              ,nelcon      ,nViz 
              ,numel       ,sistEqX->neqNov
              ,maxViz      ,ndf
              ,upper       ,diag 
              ,lower);
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
         csrJaR(&sistEqX->ia[n1] ,&sistEqX->ja[nad] 
               ,id               ,num     
               ,nelcon           ,nViz 
               ,numel            ,sistEqX->neqNov
               ,maxViz           ,ndf);
/*...................................................................*/

/*... reordenando o grafo*/
         sortGraphCsr(sistEqX->ia,sistEqX->ja,sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo*/
         sortGraphCsr(&sistEqX->ia[n1],&sistEqX->ja[nad]
                     ,sistEqX->neqNov);
/*...................................................................*/

/*...*/
         if(sistEqX->unsym) 
           nadT = 2*nad + nadr;
         else
           nadT = nad + nadr;
/*...................................................................*/

/*... alocacao da matriz*/
         HccaAlloc(DOUBLE         ,m    ,sistEqX->ad
                  ,sistEqX->neqNov,strAd,false);
         HccaAlloc(DOUBLE         ,m    ,sistEqX->al
                  ,nadT           ,strA ,false);
         zero(sistEqX->ad,sistEqX->neqNov,DOUBLEC);
         zero(sistEqX->al,nadT           ,DOUBLEC);
/*...................................................................*/
       }
/*...................................................................*/

/*... paralelo (CSRC + COO)*/
       else if (mpiVar.nPrcs > 1 && fCoo)
       {
/*... obtendo o numero de termos nao nulos na parte retangular*/ 
         nadr = sistEqX->nadr = cooNnzR(id     ,num   
                                      ,nelcon ,nViz 
                                      ,numel  ,sistEqX->neqNov
                                      ,maxViz ,ndf);
/*...................................................................*/

/*... alocando o vetor ia(2*neq+1+nadr)*/
         n0 = sistEqX->neqNov + 1; 
         n1 = n0 + sistEqX->neqNov + 1; 
         n2 = n1 + nadr; 
         HccaAlloc(INT,m,sistEqX->ia      ,n2  ,strIa   ,_AD_);
         zero(sistEqX->ia,n2    ,INTC);
/*... vetor ia do CSRC*/
         nad = sistEqX->nad = csrIa(sistEqX->ia ,id 
                                  ,num         ,nelcon 
                                  ,nViz
                                  ,numel       ,sistEqX->neqNov
                                  ,maxViz      ,ndf
                                  ,upper       ,diag         
                                  ,lower);
/*...................................................................*/

/*... vetor ja do CSRC*/
         nadT = nad + nadr; 
         HccaAlloc(INT,m,sistEqX->ja     ,nadT ,strJa   ,_AD_);
         zero(sistEqX->ja,nadT  ,INTC);
         csrJa(sistEqX->ia ,sistEqX->ja 
              ,id          ,num     
              ,nelcon      ,nViz 
              ,numel       ,sistEqX->neqNov
              ,maxViz      ,ndf
              ,upper       ,diag 
              ,lower);
/*...................................................................*/

/*... vetor ia[n0] parte retangular no formato CSR 
      (para a montagem da matriz global) 
      ia[n1] e ja da parte retangular no formato COO*/
         cooIaJaR(&sistEqX->ia[n0]
                 ,&sistEqX->ia[n1]     ,&sistEqX->ja[nad]          
                 ,id                   ,num
                 ,nelcon               ,nViz
                 ,numel                ,sistEqX->neqNov
                 ,maxViz               ,ndf);
/*...................................................................*/

/*... reordenando o grafo csrd ja(nad)*/
         sortGraphCsr(sistEqX->ia,sistEqX->ja,sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo Coo ja(nadr)*/
         sortGraphCsr(&sistEqX->ia[n0],&sistEqX->ja[nad]
                     ,sistEqX->neqNov);
/*...................................................................*/

/*...*/
         if(sistEqX->unsym) 
           nadT = 2*nad + nadr;
         else
           nadT = nad;
/*...................................................................*/

/*... alocacao da matriz*/
         HccaAlloc(DOUBLE         ,m    ,sistEqX->ad
                  ,sistEqX->neqNov,strAd,false);
         HccaAlloc(DOUBLE         ,m    ,sistEqX->al
                  ,nadT           ,strA ,false);
         zero(sistEqX->ad,sistEqX->neqNov,DOUBLEC);
         zero(sistEqX->al,nadT,DOUBLEC);
/*...................................................................*/

       }
/*...................................................................*/

/*... sequencial (CSRC - simetrico e nao simetrico)*/
       else
       {
/*...*/
         sistEqX->nadr = 0;
         n1 = sistEqX->neq + 1; 
         HccaAlloc(INT,m,sistEqX->ia,n1,strIa,_AD_);
/*... */
         nad = sistEqX->nad = csrIa(sistEqX->ia ,id 
                                   ,num         ,nelcon 
                                   ,nViz
                                   ,numel       ,sistEqX->neqNov
                                   ,maxViz      ,ndf
                                   ,upper       ,diag         
                                   ,lower       );
/*...................................................................*/

/*...*/
         HccaAlloc(INT,m,sistEqX->ja   ,nad ,strJa   ,_AD_);
         zero(sistEqX->ja,nadT  ,INTC);
         csrJa(sistEqX->ia ,sistEqX->ja,id ,num    
              ,nelcon,nViz ,numel      ,sistEqX->neqNov, maxViz,ndf
              ,upper,diag ,lower);
/*...................................................................*/

/*... reordenando o grafo*/
         sortGraphCsr(sistEqX->ia,sistEqX->ja,sistEqX->neqNov);
/*...................................................................*/

/*...*/
         if(sistEqX->unsym) 
           nadT = 2*nad;
         else
           nadT = nad;
/*...................................................................*/
 
/*... alocacao da matriz*/
         HccaAlloc(DOUBLE         ,m    ,sistEqX->ad
                  ,sistEqX->neqNov,strAd,false);
         HccaAlloc(DOUBLE         ,m    ,sistEqX->al
                  ,nadT           ,strA ,false);
         zero(sistEqX->ad,sistEqX->neqNov,DOUBLEC);
         zero(sistEqX->al,nadT           ,DOUBLEC);
/*...................................................................*/
       }
/*...................................................................*/

/*... banda da matriz*/
       sistEqX->bandCsr[BANDCSRMAX] 
       =  bandCsrC(sistEqX->ia
                        ,sistEqX->ja
                        ,sistEqX->neqNov
                        ,BANDCSRMAX);
       sistEqX->bandCsr[BANDCSRMED] 
       =  bandCsrC(sistEqX->ia
                        ,sistEqX->ja
                        ,sistEqX->neqNov
                        ,BANDCSRMED);
       
      sistEqX->bandCsr[BANDCSRMIN] 
      =  bandCsrC(sistEqX->ia
                        ,sistEqX->ja
                        ,sistEqX->neqNov
                        ,BANDCSRMIN);

       if(!mpiVar.myId  ) 
       {
         fprintf(fileLogExc,"band Maxima: %d\n",sistEqX->bandCsr[BANDCSRMAX]);
         fprintf(fileLogExc,"band Media : %d\n",sistEqX->bandCsr[BANDCSRMED]);
         fprintf(fileLogExc,"band Minima: %d\n",sistEqX->bandCsr[BANDCSRMIN]);
       }  
/*...................................................................*/
     break;
/*...................................................................*/

/*... armazenamento ELLPACK(ad,a)*/
     case ELLPACK:
/*...*/
       n1 = sistEqX->neqNov + 1; 
       HccaAlloc(INT,m,sistEqX->ia,2         ,strIa,_AD_);
       HccaAlloc(INT,m,sistEqX->ja,n1*maxViz ,strJa,_AD_);
       zero(sistEqX->ja,n1*maxViz,INTC);
/*...................................................................*/

/*...*/
       sistEqX->nad = ellPackJa(sistEqX->ia     ,sistEqX->ja 
                               ,id              ,num
                               ,nelcon          ,nViz
                               ,numel           ,n1
                               ,maxViz          ,ndf);
/*...................................................................*/

/*... alocacao da matriz*/
       HccaAlloc(DOUBLE,m,sistEqX->ad,n1        ,strAd,_AD_);
       HccaAlloc(DOUBLE,m,sistEqX->al,n1*maxViz ,strA ,_AD_);
       zero(sistEqX->ad,n1,DOUBLEC);
       zero(sistEqX->al,n1*maxViz,DOUBLEC);
/*...................................................................*/

/*... banda da matriz*/
       if(!mpiVar.myId) 
       {
         fprintf(fileLogExc,"band Maxima: %ld\n"
         ,(long) bandEllPack(sistEqX->ia,sistEqX->ja,n1,1));
         fprintf(fileLogExc,"band Minima: %ld\n"
         ,(long) bandEllPack(sistEqX->ia,sistEqX->ja,n1,3));
         fprintf(fileLogExc,"band Media : %ld\n"
         ,(long) bandEllPack(sistEqX->ia,sistEqX->ja,n1,2));
       }
/*...................................................................*/
     break;
/*...................................................................*/

     default:
       ERRO_OP(__FILE__,__func__,type);
     break;
/*...................................................................*/
  }
}
/*********************************************************************/ 

/********************************************************************
 * Data de criacao    : 15/08/2016                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * DATASTRUCSIMPLE: estrutura de dados para a matriz de coeficientes *
 *-------------------------------------------------------------------*
 * Parametros iniciais:                                              *
 *-------------------------------------------------------------------*
 * m       -> vetor de memoria                                       *
 * id      -> numeracao das equacoes por celula                      *
 * num     -> renumeracao dos elementos                              *
 * nelcon  -> adjacencia dos elementos                               *
 * nViz    -> numero de vizinhos por elemento                        *
 * numel   -> numero de elementos                                    *
 * maxViz  -> numero maximo de vizinhos                              *
 * ndf     -> grau de liberdade                                      *
 *-------------------------------------------------------------------*
 * Parametros saida:                                                 *
 *-------------------------------------------------------------------*
 * sistEqX -> estrutura de matriz atualizado ia e ja.                *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * ad -> alocados                                                    *
 * al -> alocados                                                    *
 *-------------------------------------------------------------------*
 *********************************************************************/
void dataStructSimple(Memoria *m      ,INT *id
                     ,INT *num       ,INT *nelcon
                     ,short *nViz
                     ,INT const numel, short const maxViz
                     ,short const ndf
                     ,char  *strIa   ,char *strJa
                     ,char *strAd    ,char *strA
                     ,SistEq *sistEqX)
{
  short type = sistEqX->storage;
  bool lower, diag, upper, fCoo = false;
  INT n0, n1, n2, nad, nadr, nadT = 0;

  switch (type) {
/*... armazenamento CSR(a)*/
  case CSR:
    break;
/*...................................................................*/

/*... armazenamento CSRD(a,ad)*/
/*...CSRD+COO(simetrico)*/
  case CSRDCOO:
    fCoo = true;
  case CSRD:
    lower = true;
/*...*/
    if (sistEqX->unsym)
      upper = true;
/*...*/
    else
      upper = false;
/*...*/
    diag = false;
/*...................................................................*/

/*... paralelo (CSRD + CSR - simetrico)*/
    if (mpiVar.nPrcs > 1 && !sistEqX->unsym && !fCoo) {
      n1 = sistEqX->neqNov + 1;
      n2 = 2 * n1;
      HccaAlloc(INT, m, sistEqX->ia, n2, strIa, _AD_);
      zero(sistEqX->ia, n2, INTC);
/*... */
      nad = sistEqX->nad = csrIa(sistEqX->ia,id
                                ,num        ,nelcon
                                ,nViz
                                ,numel      ,sistEqX->neqNov
                                ,maxViz     ,1
                                ,upper      ,diag
                                ,lower);
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
      nadr = sistEqX->nadr = csrIaR(&sistEqX->ia[n1],id
                                   ,num             ,nelcon
                                   ,nViz
                                   ,numel           ,sistEqX->neqNov
                                   ,maxViz          ,1);
/*...................................................................*/

/*...*/
      nadT = nad + nadr;
      HccaAlloc(INT, m, sistEqX->ja, nadT, strJa, _AD_);
      zero(sistEqX->ja, nadT, INTC);
      csrJa(sistEqX->ia,sistEqX->ja
            ,id        ,num
            ,nelcon    ,nViz
            ,numel     ,sistEqX->neqNov
            ,maxViz    ,1
            ,upper     ,diag
            ,lower);
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
      csrJaR(&sistEqX->ia[n1],&sistEqX->ja[nad]
            ,id              ,num
            ,nelcon          ,nViz
            ,numel           ,sistEqX->neqNov
            ,maxViz          ,1);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(&sistEqX->ia[n1], &sistEqX->ja[nad]
                  ,sistEqX->neqNov);
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m, sistEqX->ad
               ,sistEqX->neqNov*ndf, strAd, false);
      HccaAlloc(DOUBLE, m, sistEqX->al
               ,nadT, strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*ndf, DOUBLEC);
      zero(sistEqX->al, nadT, DOUBLEC);
/*...................................................................*/

    }
/*...................................................................*/

/*... paralelo (CSRD + COO - simetrico)*/
    else if (mpiVar.nPrcs > 1 && !sistEqX->unsym && fCoo) {
/*... obtendo o numero de termos nao nulos na parte retangular*/
      nadr = sistEqX->nadr = cooNnzR(id    ,num
                                    ,nelcon,nViz
                                    ,numel ,sistEqX->neqNov
                                    ,maxViz,1);
/*...................................................................*/

/*... alocando o vetor ia(2*neq+1+nadr)*/
      n0 = sistEqX->neqNov + 1;
      n1 = n0 + sistEqX->neqNov + 1;
      n2 = n1 + nadr;
      HccaAlloc(INT, m, sistEqX->ia, n2, strIa, _AD_);
      zero(sistEqX->ia, n2, INTC);
/*... vetor ia do CSRD*/
      nad = sistEqX->nad = csrIa(sistEqX->ia, id
                                ,num        ,nelcon
                                ,nViz
                                ,numel      ,sistEqX->neqNov
                                ,maxViz     ,1
                                ,upper      ,diag
                                ,lower);
/*...................................................................*/

/*... vetor ja do CSRD*/
      nadT = nad + nadr;
      HccaAlloc(INT, m, sistEqX->ja, nadT, strJa, _AD_);
      zero(sistEqX->ja, nadT, INTC);
      csrJa(sistEqX->ia,sistEqX->ja
           ,id         ,num
           ,nelcon     ,nViz
           ,numel      ,sistEqX->neqNov
           ,maxViz     ,1
           ,upper      ,diag
           ,lower);
/*...................................................................*/

/*... vetor ia[n0] parte retangular no formato CSR
      (para a montagem da matriz global)
      ia[n1] e ja da parte retangular no formato COO*/
      cooIaJaR(&sistEqX->ia[n0]
              ,&sistEqX->ia[n1],&sistEqX->ja[nad]
              ,id              ,num
              ,nelcon          ,nViz
              ,numel           ,sistEqX->neqNov
              ,maxViz          ,1);
/*...................................................................*/

/*... reordenando o grafo csrd ja(nad)*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo Coo ja(nadr)*/
      sortGraphCsr(&sistEqX->ia[n0], &sistEqX->ja[nad]
                  ,sistEqX->neqNov);
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m, sistEqX->ad
               ,sistEqX->neqNov*ndf, strAd, false);
      HccaAlloc(DOUBLE, m, sistEqX->al
               ,nadT, strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*ndf, DOUBLEC);
      zero(sistEqX->al, nadT, DOUBLEC);
/*...................................................................*/

    }
/*...................................................................*/

/*... sequencial (CSRD - simetrico e nao simetrico)*/
/*... paralelo   (CSRD - nao simetrico)*/
    else {
      sistEqX->nadr = 0;
      n1 = sistEqX->neq + 1;
      HccaAlloc(INT, m, sistEqX->ia, n1, strIa, _AD_);
/*... */
      nad = sistEqX->nad = csrIa(sistEqX->ia, id
                                , num, nelcon
                                , nViz
                                , numel, sistEqX->neqNov
                                , maxViz, 1
                                , upper, diag
                                , lower);
/*...................................................................*/

/*...*/
      HccaAlloc(INT, m, sistEqX->ja, nad, strJa, _AD_);
      csrJa(sistEqX->ia, sistEqX->ja, id, num
        , nelcon, nViz, numel, sistEqX->neqNov, maxViz, 1
        , upper, diag, lower);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m, sistEqX->ad
        , sistEqX->neqNov*ndf, strAd, false);
      HccaAlloc(DOUBLE, m, sistEqX->al
        , nad, strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*ndf, DOUBLEC);
      zero(sistEqX->al, nad, DOUBLEC);
/*...................................................................*/

    }
/*...................................................................*/

/*... banda da matriz*/
    sistEqX->bandCsr[BANDCSRMAX]
      = bandCsr(sistEqX->ia
        , sistEqX->ja
        , sistEqX->neqNov
        , BANDCSRMAX);
    sistEqX->bandCsr[BANDCSRMED]
      = bandCsr(sistEqX->ia
        , sistEqX->ja
        , sistEqX->neqNov
        , BANDCSRMED);

    sistEqX->bandCsr[BANDCSRMIN]
      = bandCsr(sistEqX->ia
        , sistEqX->ja
        , sistEqX->neqNov
        , BANDCSRMIN);

    if (!mpiVar.myId) {
      fprintf(fileLogExc,"band Maxima: %d\n", sistEqX->bandCsr[BANDCSRMAX]);
      fprintf(fileLogExc,"band Media : %d\n", sistEqX->bandCsr[BANDCSRMED]);
      fprintf(fileLogExc,"band Minima: %d\n", sistEqX->bandCsr[BANDCSRMIN]);
    }
/*...................................................................*/
    break;
/*...................................................................*/

/*... armazenamento CSRC(ad,au,al)*/
/*...CSRC+COO*/
  case CSRCCOO:
    fCoo = true;
  case CSRC:
/*...*/
    lower = true;
/*...*/
    upper = false;
/*...*/
    diag = false;
/*...................................................................*/

/*... paralelo (CSRC + CSR)*/
    if (mpiVar.nPrcs > 1 && !fCoo) {
      n1 = sistEqX->neqNov + 1;
      n2 = 2 * n1;
      HccaAlloc(INT, m, sistEqX->ia, n2, strIa, _AD_);
      zero(sistEqX->ia, n2, INTC);

/*... */
      nad = sistEqX->nad = csrIa(sistEqX->ia, id
                                ,num        ,nelcon
                                ,nViz
                                ,numel      ,sistEqX->neqNov
                                ,maxViz     ,1
                                ,upper      ,diag
                                ,lower);
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
      nadr = sistEqX->nadr = csrIaR(&sistEqX->ia[n1], id
                                   ,num             ,nelcon
                                   ,nViz
                                   ,numel           ,sistEqX->neqNov
                                   ,maxViz          ,1);
/*...................................................................*/

/*...*/
      nadT = nad + nadr;
      HccaAlloc(INT, m, sistEqX->ja, nadT, strJa, _AD_);
      zero(sistEqX->ja, nadT, INTC);
      csrJa(sistEqX->ia,sistEqX->ja
            ,id        ,num
            ,nelcon    ,nViz
            ,numel     ,sistEqX->neqNov
            ,maxViz    ,1
            ,upper     ,diag
            ,lower);
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
      csrJaR(&sistEqX->ia[n1], &sistEqX->ja[nad]
            , id, num
            , nelcon, nViz
            , numel, sistEqX->neqNov
            , maxViz, 1);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(&sistEqX->ia[n1],&sistEqX->ja[nad]
                  ,sistEqX->neqNov);
/*...................................................................*/

/*...*/
      if (sistEqX->unsym)
        nadT = 2 * nad + nadr;
      else
        nadT = nad + nadr;
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m,sistEqX->ad
               ,sistEqX->neqNov*ndf, strAd, false);
      HccaAlloc(DOUBLE,m   , sistEqX->al
               ,nadT  ,strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*ndf, DOUBLEC);
      zero(sistEqX->al, nadT, DOUBLEC);
/*...................................................................*/
    }
/*...................................................................*/

/*... paralelo (CSRC + COO)*/
    else if (mpiVar.nPrcs > 1 && fCoo) {
/*... obtendo o numero de termos nao nulos na parte retangular*/
      nadr = sistEqX->nadr = cooNnzR(id    ,num
                                    ,nelcon,nViz
                                    ,numel ,sistEqX->neqNov
                                    ,maxViz,1);
/*...................................................................*/

/*... alocando o vetor ia(2*neq+1+nadr)*/
      n0 = sistEqX->neqNov + 1;
      n1 = n0 + sistEqX->neqNov + 1;
      n2 = n1 + nadr;
      HccaAlloc(INT, m, sistEqX->ia, n2, strIa, _AD_);
      zero(sistEqX->ia, n2, INTC);
/*... vetor ia do CSRC*/
      nad = sistEqX->nad = csrIa(sistEqX->ia,id
                                ,num        ,nelcon
                                ,nViz
                                ,numel      ,sistEqX->neqNov
                                ,maxViz     ,1
                                ,upper      ,diag
                                ,lower);
/*...................................................................*/

/*... vetor ja do CSRC*/
      nadT = nad + nadr;
      HccaAlloc(INT, m, sistEqX->ja, nadT, strJa, _AD_);
      zero(sistEqX->ja, nadT, INTC);
      csrJa(sistEqX->ia,sistEqX->ja
           ,id         ,num
           ,nelcon     ,nViz
           ,numel      ,sistEqX->neqNov
           ,maxViz     ,1
           ,upper      ,diag
           ,lower);
/*...................................................................*/

/*... vetor ia[n0] parte retangular no formato CSR
      (para a montagem da matriz global)
      ia[n1] e ja da parte retangular no formato COO*/
      cooIaJaR(&sistEqX->ia[n0]
              ,&sistEqX->ia[n1],&sistEqX->ja[nad]
              ,id              ,num
              ,nelcon          ,nViz
              ,numel           ,sistEqX->neqNov
              ,maxViz          ,1);
/*...................................................................*/

/*... reordenando o grafo csrd ja(nad)*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo Coo ja(nadr)*/
      sortGraphCsr(&sistEqX->ia[n0], &sistEqX->ja[nad]
                  ,sistEqX->neqNov);
/*...................................................................*/

/*...*/
      if (sistEqX->unsym)
        nadT = 2 * nad + nadr;
      else
        nadT = nad;
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m, sistEqX->ad
               ,sistEqX->neqNov*ndf, strAd, false);
      HccaAlloc(DOUBLE, m, sistEqX->al
               ,nadT, strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*ndf, DOUBLEC);
      zero(sistEqX->al, nadT, DOUBLEC);
/*...................................................................*/

    }
/*...................................................................*/

/*... sequencial (CSRC - simetrico e nao simetrico)*/
    else {
/*...*/
      sistEqX->nadr = 0;
      n1 = sistEqX->neq + 1;
      HccaAlloc(INT, m, sistEqX->ia, n1, strIa, _AD_);
/*... */
      nad = sistEqX->nad = csrIa(sistEqX->ia,id
                                ,num        ,nelcon
                                ,nViz
                                ,numel      ,sistEqX->neqNov
                                ,maxViz     ,1
                                ,upper      ,diag
                                ,lower);
/*...................................................................*/

/*...*/
      HccaAlloc(INT, m, sistEqX->ja, nad, strJa, _AD_);
      zero(sistEqX->ja, nadT, INTC);
      csrJa(sistEqX->ia, sistEqX->ja, id, num
           ,nelcon,nViz, numel, sistEqX->neqNov, maxViz, 1
           ,upper ,diag, lower);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

 /*...*/
      if (sistEqX->unsym)
        nadT = 2 * nad;
      else
        nadT = nad;
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m, sistEqX->ad
        , sistEqX->neqNov*ndf, strAd, false);
      HccaAlloc(DOUBLE, m, sistEqX->al
        , nadT, strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*ndf, DOUBLEC);
      zero(sistEqX->al, nadT, DOUBLEC);
/*...................................................................*/
    }
/*...................................................................*/

/*... banda da matriz*/
    sistEqX->bandCsr[BANDCSRMAX]
      = bandCsrC(sistEqX->ia
        , sistEqX->ja
        , sistEqX->neqNov
        , BANDCSRMAX);
    sistEqX->bandCsr[BANDCSRMED]
      = bandCsrC(sistEqX->ia
        , sistEqX->ja
        , sistEqX->neqNov
        , BANDCSRMED);

    sistEqX->bandCsr[BANDCSRMIN]
      = bandCsrC(sistEqX->ia
        , sistEqX->ja
        , sistEqX->neqNov
        , BANDCSRMIN);

    if (!mpiVar.myId) {
      fprintf(fileLogExc,"band Maxima: %d\n", sistEqX->bandCsr[BANDCSRMAX]);
      fprintf(fileLogExc,"band Media : %d\n", sistEqX->bandCsr[BANDCSRMED]);
      fprintf(fileLogExc,"band Minima: %d\n", sistEqX->bandCsr[BANDCSRMIN]);
    }
/*...................................................................*/
    break;
/*...................................................................*/

/*... armazenamento ELLPACK(ad,a)*/
  case ELLPACK:
/*...*/
    n1 = sistEqX->neqNov + 1;
    HccaAlloc(INT, m, sistEqX->ia, 2, strIa, _AD_);
    HccaAlloc(INT, m, sistEqX->ja, n1*maxViz, strJa, _AD_);
    zero(sistEqX->ja, n1*maxViz, INTC);
/*...................................................................*/

/*...*/
    sistEqX->nad = ellPackJa(sistEqX->ia,sistEqX->ja
                            ,id         ,num
                            ,nelcon     ,nViz
                            ,numel      ,n1
                            ,maxViz     ,1);
/*...................................................................*/

/*... alocacao da matriz*/
    HccaAlloc(DOUBLE, m, sistEqX->ad, n1*ndf, strAd, _AD_);
    HccaAlloc(DOUBLE, m, sistEqX->al, n1*maxViz, strA, _AD_);
    zero(sistEqX->ad, n1, DOUBLEC);
    zero(sistEqX->al, n1*maxViz, DOUBLEC);
/*...................................................................*/

/*... banda da matriz*/
    if (!mpiVar.myId) {
      fprintf(fileLogExc,"band Maxima: %ld\n"
        , (long)bandEllPack(sistEqX->ia, sistEqX->ja, n1, 1));
      fprintf(fileLogExc,"band Minima: %ld\n"
        , (long)bandEllPack(sistEqX->ia, sistEqX->ja, n1, 3));
      fprintf(fileLogExc,"band Media : %ld\n"
        , (long)bandEllPack(sistEqX->ia, sistEqX->ja, n1, 2));
    }
 /*...................................................................*/
    break;
/*...................................................................*/

  default:
    ERRO_OP(__FILE__, __func__, type);
    break;
/*...................................................................*/
  }
}
/*********************************************************************/

/********************************************************************
* Data de criacao    : 05/08/2016                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* DataStructBlockE: estrutura de dados para a matriz de coeficientes*
*-------------------------------------------------------------------*
* Parametros iniciais:                                              *
*-------------------------------------------------------------------*
* m       -> vetor de memoria                                       *
* id      -> numeracao das equacoes por celula                      *
* num     -> renumeracao dos elementos                              *
* nelcon  -> adjacencia dos elementos                               *
* nViz    -> numero de vizinhos por elemento                        *
* numel   -> numero de elementos                                    *
* nAd     -> multiplica o numero de termos no vetor ad              *
* nAl     -> multiplica o numero de termos no vetor al              *
* strIa   -> nome do arranjo ia                                     *
* strJa   -> nome do arranjo ja                                     *
* strAd   -> nome do arranjo Ad                                     *
* strAl   -> nome do arranjo Al                                     *
* maxViz  -> numero maximo de vizinhos                              *
* ndf     -> grau de liberdade                                      *
*-------------------------------------------------------------------*
* Parametros saida:                                                 *
*-------------------------------------------------------------------*
* sistEqX -> estrutura de matriz atualizado ia e ja.                *
*-------------------------------------------------------------------*
* OBS:                                                              *
* ad -> alocados[nAd*neq]                                           *
* al -> alocados[nAl*nad]                                           *
*-------------------------------------------------------------------*
*********************************************************************/
void dataStructBlock(Memoria *m             , INT *id
                   , INT *num               , INT *nelcon
                   , short *nViz
                   , INT const numel        , short const maxViz
                   , short const nAd        , short const nAl
                   , const char *const strIa, const char *const strJa
                   , const char *const strAd, const char *const strA
                   , SistEq *sistEqX)
{
  short type = sistEqX->storage;
  bool lower, diag, upper, fCoo = false;
  INT n0, n1, n2, nad, nadr, nadT = 0;

  switch (type) {
/*... armazenamento CSR(a)*/
  case CSR:
    break;
/*...................................................................*/

/*... armazenamento CSRD(a,ad)*/
/*...CSRD+COO(simetrico)*/
  case CSRDCOO:
    fCoo = true;
  case CSRD:
    lower = true;
    /*...*/
    if (sistEqX->unsym)
      upper = true;
    /*...*/
    else
      upper = false;
    /*...*/
    diag = false;
    /*...................................................................*/

/*... paralelo (CSRD + CSR - simetrico)*/
    if (mpiVar.nPrcs > 1 && !sistEqX->unsym && !fCoo) {
      n1 = sistEqX->neqNov + 1;
      n2 = 2 * n1;
      HccaAlloc(INT, m, sistEqX->ia, n2, strIa, _AD_);
      zero(sistEqX->ia, n2, INTC);
/*... */
      nad = sistEqX->nad = csrIa(sistEqX->ia, id
        , num, nelcon
        , nViz
        , numel, sistEqX->neqNov
        , maxViz, 1
        , upper, diag
        , lower);
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
      nadr = sistEqX->nadr = csrIaR(&sistEqX->ia[n1], id
                                  , num             , nelcon
                                  , nViz
                                  , numel           , sistEqX->neqNov
                                  , maxViz          , 1);
/*...................................................................*/

/*...*/
      nadT = nad + nadr;
      HccaAlloc(INT, m, sistEqX->ja, nadT, strJa, _AD_);
      zero(sistEqX->ja, nadT, INTC);
      csrJa(sistEqX->ia, sistEqX->ja
          , id         , num
          , nelcon     , nViz
          , numel      , sistEqX->neqNov
          , maxViz     , 1
          , upper      , diag
          , lower);
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
      csrJaR(&sistEqX->ia[n1], &sistEqX->ja[nad]
        , id, num
        , nelcon, nViz
        , numel, sistEqX->neqNov
        , maxViz, 1);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(&sistEqX->ia[n1], &sistEqX->ja[nad]
                  , sistEqX->neqNov);
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m          , sistEqX->ad
              , sistEqX->neqNov*nAd, strAd, false);
      HccaAlloc(DOUBLE, m   , sistEqX->al
               , nadT*nAl   , strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*nAd, DOUBLEC);
      zero(sistEqX->al, nadT*nAl, DOUBLEC);
/*...................................................................*/

    }
/*...................................................................*/

/*... paralelo (CSRD + COO - simetrico)*/
    else if (mpiVar.nPrcs > 1 && !sistEqX->unsym && fCoo) 
    {
/*... obtendo o numero de termos nao nulos na parte retangular*/
      nadr = sistEqX->nadr = cooNnzR(id   , num
                                  , nelcon, nViz
                                  , numel , sistEqX->neqNov
                                  , maxViz, 1);
/*...................................................................*/

/*... alocando o vetor ia(2*neq+1+nadr)*/
      n0 = sistEqX->neqNov + 1;
      n1 = n0 + sistEqX->neqNov + 1;
      n2 = n1 + nadr;
      HccaAlloc(INT, m, sistEqX->ia, n2, strIa, _AD_);
      zero(sistEqX->ia, n2, INTC);
      /*... vetor ia do CSRD*/
      nad = sistEqX->nad = csrIa(sistEqX->ia, id
                               , num        , nelcon
                               , nViz
                               , numel      , sistEqX->neqNov
                               , maxViz     , 1
                               , upper      , diag
                               , lower);
/*...................................................................*/

/*... vetor ja do CSRD*/
      nadT = nad + nadr;
      HccaAlloc(INT, m, sistEqX->ja, nadT, strJa, _AD_);
      zero(sistEqX->ja, nadT, INTC);
      csrJa(sistEqX->ia, sistEqX->ja
          , id         , num
          , nelcon     , nViz
          , numel      , sistEqX->neqNov
          , maxViz     , 1
          , upper      , diag
          , lower);
/*...................................................................*/

/*... vetor ia[n0] parte retangular no formato CSR
      (para a montagem da matriz global)
      ia[n1] e ja da parte retangular no formato COO*/
      cooIaJaR(&sistEqX->ia[n0]
             , &sistEqX->ia[n1], &sistEqX->ja[nad]
             , id              , num
             , nelcon          , nViz
             , numel           , sistEqX->neqNov
             , maxViz          , 1);
/*...................................................................*/

/*... reordenando o grafo csrd ja(nad)*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo Coo ja(nadr)*/
      sortGraphCsr(&sistEqX->ia[n0], &sistEqX->ja[nad]
                 , sistEqX->neqNov);
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m, sistEqX->ad
              , sistEqX->neqNov*nAd, strAd, false);
      HccaAlloc(DOUBLE, m, sistEqX->al
              , nadT*nAl, strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*nAd, DOUBLEC);
      zero(sistEqX->al, nadT*nAl, DOUBLEC);
/*...................................................................*/

    }
/*...................................................................*/

/*... sequencial (CSRD - simetrico e nao simetrico)*/
/*... paralelo   (CSRD - nao simetrico)*/
    else {
      sistEqX->nadr = 0;
      n1 = sistEqX->neq + 1;
      HccaAlloc(INT, m, sistEqX->ia, n1, strIa, _AD_);
/*... */
      nad = sistEqX->nad = csrIa(sistEqX->ia, id
        , num, nelcon
        , nViz
        , numel, sistEqX->neqNov
        , maxViz, 1
        , upper, diag
        , lower);
/*...................................................................*/

/*...*/
      HccaAlloc(INT, m, sistEqX->ja, nad, strJa, _AD_);
      csrJa(sistEqX->ia, sistEqX->ja, id, num
        , nelcon, nViz, numel, sistEqX->neqNov, maxViz, 1
        , upper, diag, lower);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m, sistEqX->ad
              , sistEqX->neqNov*nAd, strAd, false);
      HccaAlloc(DOUBLE, m, sistEqX->al
              , nad*nAl, strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*nAd, DOUBLEC);
      zero(sistEqX->al, nad*nAl, DOUBLEC);
/*...................................................................*/

    }
/*...................................................................*/

/*... banda da matriz*/
    sistEqX->bandCsr[BANDCSRMAX] = bandCsr(sistEqX->ia      , sistEqX->ja
                                         , sistEqX->neqNov  , BANDCSRMAX);
    sistEqX->bandCsr[BANDCSRMED] = bandCsr(sistEqX->ia      , sistEqX->ja
                                         , sistEqX->neqNov  , BANDCSRMED);

    sistEqX->bandCsr[BANDCSRMIN] = bandCsr(sistEqX->ia      , sistEqX->ja
                                         , sistEqX->neqNov  , BANDCSRMIN);

    if (!mpiVar.myId) 
    {
      fprintf(fileLogExc, "band Maxima: %d\n", sistEqX->bandCsr[BANDCSRMAX]);
      fprintf(fileLogExc, "band Media : %d\n", sistEqX->bandCsr[BANDCSRMED]);
      fprintf(fileLogExc, "band Minima: %d\n", sistEqX->bandCsr[BANDCSRMIN]);
    }
/*...................................................................*/
    break;
/*...................................................................*/

/*... armazenamento CSRC(ad,au,al)*/
/*...CSRC+COO*/
  case CSRCCOO:
    fCoo = true;
  case CSRC:
/*...*/
    lower = true;
/*...*/
    upper = false;
/*...*/
    diag = false;
/*...................................................................*/

/*... paralelo (CSRC + CSR)*/
    if (mpiVar.nPrcs > 1 && !fCoo) 
    {
      n1 = sistEqX->neqNov + 1;
      n2 = 2 * n1;
      HccaAlloc(INT, m, sistEqX->ia, n2, strIa, _AD_);
      zero(sistEqX->ia, n2, INTC);

/*... */
      nad = sistEqX->nad = csrIa(sistEqX->ia, id
        , num, nelcon
        , nViz
        , numel, sistEqX->neqNov
        , maxViz, 1
        , upper, diag
        , lower);
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
      nadr = sistEqX->nadr = csrIaR(&sistEqX->ia[n1], id
        , num, nelcon
        , nViz
        , numel, sistEqX->neqNov
        , maxViz, 1);
/*...................................................................*/

/*...*/
      nadT = nad + nadr;
      HccaAlloc(INT, m, sistEqX->ja, nadT, strJa, _AD_);
      zero(sistEqX->ja, nadT, INTC);
      csrJa(sistEqX->ia, sistEqX->ja
        , id, num
        , nelcon, nViz
        , numel, sistEqX->neqNov
        , maxViz, 1
        , upper, diag
        , lower);
/*...................................................................*/

/*... parte retangular*/
/* ... CSR*/
      csrJaR(&sistEqX->ia[n1], &sistEqX->ja[nad]
        , id, num
        , nelcon, nViz
        , numel, sistEqX->neqNov
        , maxViz, 1);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(&sistEqX->ia[n1], &sistEqX->ja[nad]
        , sistEqX->neqNov);
/*...................................................................*/

/*...*/
      if (sistEqX->unsym)
        nadT = 2 * nad + nadr;
      else
        nadT = nad + nadr;
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m, sistEqX->ad
              , sistEqX->neqNov*nAd, strAd, false);
      HccaAlloc(DOUBLE, m, sistEqX->al
              , nadT*nAl, strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*nAd, DOUBLEC);
      zero(sistEqX->al, nadT*nAl, DOUBLEC);
/*...................................................................*/
    }
/*...................................................................*/

/*... paralelo (CSRC + COO)*/
    else if (mpiVar.nPrcs > 1 && fCoo) 
    {
/*... obtendo o numero de termos nao nulos na parte retangular*/
      nadr = sistEqX->nadr = cooNnzR(id    , num
                                   , nelcon, nViz
                                   , numel , sistEqX->neqNov
                                   , maxViz, 1);
/*...................................................................*/

/*... alocando o vetor ia(2*neq+1+nadr)*/
      n0 = sistEqX->neqNov + 1;
      n1 = n0 + sistEqX->neqNov + 1;
      n2 = n1 + nadr;
      HccaAlloc(INT, m, sistEqX->ia, n2, strIa, _AD_);
      zero(sistEqX->ia, n2, INTC);
/*... vetor ia do CSRC*/
      nad = sistEqX->nad = csrIa(sistEqX->ia, id
                               , num        , nelcon
                               , nViz
                               , numel      , sistEqX->neqNov
                               , maxViz     , 1
                               , upper      , diag
                               , lower);
/*...................................................................*/

/*... vetor ja do CSRC*/
      nadT = nad + nadr;
      HccaAlloc(INT, m, sistEqX->ja, nadT, strJa, _AD_);
      zero(sistEqX->ja, nadT, INTC);
      csrJa(sistEqX->ia, sistEqX->ja
        , id, num
        , nelcon, nViz
        , numel, sistEqX->neqNov
        , maxViz, 1
        , upper, diag
        , lower);
/*...................................................................*/

/*... vetor ia[n0] parte retangular no formato CSR
      (para a montagem da matriz global)
      ia[n1] e ja da parte retangular no formato COO*/
      cooIaJaR(&sistEqX->ia[n0]
        , &sistEqX->ia[n1], &sistEqX->ja[nad]
        , id, num
        , nelcon, nViz
        , numel, sistEqX->neqNov
        , maxViz, 1);
/*...................................................................*/

/*... reordenando o grafo csrd ja(nad)*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

/*... reordenando o grafo Coo ja(nadr)*/
      sortGraphCsr(&sistEqX->ia[n0], &sistEqX->ja[nad]
        , sistEqX->neqNov);
/*...................................................................*/

/*...*/
      if (sistEqX->unsym)
        nadT = 2 * nad + nadr;
      else
        nadT = nad;
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m, sistEqX->ad
        , sistEqX->neqNov*nAd, strAd, false);
      HccaAlloc(DOUBLE, m, sistEqX->al
        , nadT*nAl, strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*nAd, DOUBLEC);
      zero(sistEqX->al, nadT*nAl, DOUBLEC);
/*...................................................................*/

    }
/*...................................................................*/

/*... sequencial (CSRC - simetrico e nao simetrico)*/
    else 
    {
/*...*/
      sistEqX->nadr = 0;
      n1 = sistEqX->neq + 1;
      HccaAlloc(INT, m, sistEqX->ia, n1, strIa, _AD_);
/*... */
      nad = sistEqX->nad = csrIa(sistEqX->ia, id
        , num, nelcon
        , nViz
        , numel, sistEqX->neqNov
        , maxViz, 1
        , upper, diag
        , lower);
/*...................................................................*/

/*...*/
      HccaAlloc(INT, m, sistEqX->ja, nad, strJa, _AD_);
      zero(sistEqX->ja, nadT, INTC);
      csrJa(sistEqX->ia, sistEqX->ja, id, num
        , nelcon, nViz, numel, sistEqX->neqNov, maxViz, 1
        , upper, diag, lower);
/*...................................................................*/

/*... reordenando o grafo*/
      sortGraphCsr(sistEqX->ia, sistEqX->ja, sistEqX->neqNov);
/*...................................................................*/

/*...*/
      if (sistEqX->unsym)
        nadT = 2 * nad;
      else
        nadT = nad;
/*...................................................................*/

/*... alocacao da matriz*/
      HccaAlloc(DOUBLE, m, sistEqX->ad
        , sistEqX->neqNov*nAd, strAd, false);
      HccaAlloc(DOUBLE, m, sistEqX->al
        , nadT*nAl, strA, false);
      zero(sistEqX->ad, sistEqX->neqNov*nAd, DOUBLEC);
      zero(sistEqX->al, nadT*nAl, DOUBLEC);
/*...................................................................*/
    }
/*...................................................................*/

/*... banda da matriz*/
    sistEqX->bandCsr[BANDCSRMAX] = bandCsrC(sistEqX->ia    , sistEqX->ja
                                          , sistEqX->neqNov, BANDCSRMAX);
    sistEqX->bandCsr[BANDCSRMED] = bandCsrC(sistEqX->ia    , sistEqX->ja
                                          , sistEqX->neqNov, BANDCSRMED);

    sistEqX->bandCsr[BANDCSRMIN] = bandCsrC(sistEqX->ia    , sistEqX->ja
                                          , sistEqX->neqNov, BANDCSRMIN);

    if (!mpiVar.myId) 
    {
      fprintf(fileLogExc, "band Maxima: %d\n", sistEqX->bandCsr[BANDCSRMAX]);
      fprintf(fileLogExc, "band Media : %d\n", sistEqX->bandCsr[BANDCSRMED]);
      fprintf(fileLogExc, "band Minima: %d\n", sistEqX->bandCsr[BANDCSRMIN]);
    }
/*...................................................................*/
    break;
/*...................................................................*/

/*... armazenamento ELLPACK(ad,a)*/
  case ELLPACK:
/*...*/
    n1 = sistEqX->neqNov + 1;
    HccaAlloc(INT, m, sistEqX->ia, 2, strIa, _AD_);
    HccaAlloc(INT, m, sistEqX->ja, n1*maxViz, strJa, _AD_);
    zero(sistEqX->ja, n1*maxViz, INTC);
/*...................................................................*/

/*...*/
    sistEqX->nad = ellPackJa(sistEqX->ia, sistEqX->ja
                           , id, num
                           , nelcon, nViz
                           , numel, n1
                           , maxViz, 1);
/*...................................................................*/

/*... alocacao da matriz*/
    HccaAlloc(DOUBLE, m, sistEqX->ad, n1*nAd, strAd, _AD_);
    HccaAlloc(DOUBLE, m, sistEqX->al, n1*maxViz*nAl, strA, _AD_);
    zero(sistEqX->ad, n1*nAd, DOUBLEC);
    zero(sistEqX->al, n1*maxViz*nAl, DOUBLEC);
/*...................................................................*/

/*... banda da matriz*/
    if (!mpiVar.myId) 
    {
      fprintf(fileLogExc, "band Maxima: %ld\n"
        , (long)bandEllPack(sistEqX->ia, sistEqX->ja, n1, 1));
      fprintf(fileLogExc, "band Minima: %ld\n"
        , (long)bandEllPack(sistEqX->ia, sistEqX->ja, n1, 3));
      fprintf(fileLogExc, "band Media : %ld\n"
        , (long)bandEllPack(sistEqX->ia, sistEqX->ja, n1, 2));
    }
/*...................................................................*/
    break;
/*...................................................................*/

  default:
    ERRO_OP(__FILE__, __func__, type);
    break;
/*...................................................................*/
  }
}
/*********************************************************************/

/**********************************************************************
 *SETDATASTRUCT : escolhe a estruturade dados                         *
 **********************************************************************/
void setDataStruct(char *word,short *data)
{

/*... CSR*/
  if(!strcmp(word,"csr")){
    fprintf(fileLogExc,"%-20s : %s","DataStruct","CSR\n");
    *data = CSR;
  }
/*...................................................................*/

/*... CSRD*/
  else if(!strcmp(word,"card")){
    fprintf(fileLogExc,"%-20s : %s","DataStruct","CSRD\n");
    *data = CSRD;
  }
/*...................................................................*/

/*... ELLPACK*/
  else if(!strcmp(word,"ellpack")){
    fprintf(fileLogExc,"%-20s : %s","DataStruct","ELLPACK\n");
    *data = ELLPACK;
  }
/*...................................................................*/

/*... CSRDCOO*/
  else if(!strcmp(word,"csrdcoo")){
    fprintf(fileLogExc,"%-20s : %s","DataStruct","CSRDCOO\n");
    *data = CSRDCOO;
  }
/*...................................................................*/

/*... CSRC*/
  else if(!strcmp(word,"csrc")){
    fprintf(fileLogExc,"%-20s : %s","DataStruct","CSRC\n");
    *data = CSRC;
  }
/*...................................................................*/

/*... CSRCCOO*/
  else if(!strcmp(word,"csrccoo")){
    fprintf(fileLogExc,"%-20s : %s","DataStruct","CSRCCOO\n");
    *data = CSRCCOO;
  }
/*...................................................................*/

} 
/*********************************************************************/      
