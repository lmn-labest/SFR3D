#include<Sisteq.h>
/********************************************************************* 
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
 * au -> alocados                                                    * 
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
   INT n1,n2,nad,nadr,nadT=0;
   
   switch(type){
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
       if(sistEqX->unsym) 
         upper = true;
/*...*/
       else 
         upper = false;
/*...*/
       diag    = false;
/*...................................................................*/

/*... paralelo (CSRD + CSR - simetrico)*/
       if( mpiVar.nPrcs > 1 && !sistEqX->unsym && !fCoo) { 
         n1 = sistEqX->neqNov + 1; 
         n2 = 2*n1; 
         HccaAlloc(INT,m,sistEqX->ia      ,n2  ,strIa   ,_AD_);
/*... */
         nad = sistEqX->nad = csrIa(sistEqX->ia ,id 
                                   ,num       ,nelcon 
                                   ,nViz
                                   ,numel     ,sistEqX->neqNov
                                   ,maxViz    ,ndf
                                   ,upper     ,diag         
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
      else if(mpiVar.nPrcs > 1 && !sistEqX->unsym && fCoo) {
/*... obtendo o numero de termos nao nulos na parte retangular*/ 
        nadr = sistEqX->nadr = cooNnzR(id     ,num   
                                      ,nelcon ,nViz 
                                      ,numel  ,sistEqX->neqNov
                                      ,maxViz ,ndf);
/*...................................................................*/

/*... alocando o vetor ia(neq+1+nadr)*/
        n1 = sistEqX->neqNov + 1; 
        n2 = n1+nadr; 
        HccaAlloc(INT,m,sistEqX->ia      ,n2  ,strIa   ,_AD_);
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
         csrJa(sistEqX->ia ,sistEqX->ja 
              ,id          ,num     
              ,nelcon      ,nViz 
              ,numel       ,sistEqX->neqNov
              ,maxViz      ,ndf
              ,upper       ,diag 
              ,lower);
/*...................................................................*/

/*... vetor ia e ja do COO*/
         cooIaJaR(&sistEqX->ia[n1]     ,&sistEqX->ja[nad]          
                 ,id                   ,num
                 ,nelcon               ,nViz
                 ,numel                ,sistEqX->neqNov
                 ,maxViz               ,ndf);
/*...................................................................*/

/*... reordenando o grafo ja(nad)*/
         sortGraphCsr(sistEqX->ia,sistEqX->ja,sistEqX->neqNov);
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
       else{ 
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

       if(!mpiVar.myId  ) {
         printf("band Maxima: %d\n",sistEqX->bandCsr[BANDCSRMAX]);
         printf("band Media : %d\n",sistEqX->bandCsr[BANDCSRMED]);
         printf("band Minima: %d\n",sistEqX->bandCsr[BANDCSRMIN]);
       }
/*...................................................................*/
     break;
/*...................................................................*/

/*... armazenamento CSRC(ad,au,al)*/
     case CSRC:
     break;
/*...................................................................*/

/*... armazenamento ELLPACK(ad,d)*/
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
       if(!mpiVar.myId) {
         printf("band Maxima: %ld\n"
         ,(long) bandEllPack(sistEqX->ia,sistEqX->ja,n1,1));
         printf("band Minima: %ld\n"
         ,(long) bandEllPack(sistEqX->ia,sistEqX->ja,n1,3));
         printf("band Media : %ld\n"
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

/**********************************************************************
 *SETDATASTRUCT : escolhe a estruturade dados                         *
 **********************************************************************/
void setDataStruct(char *word,short *data)
{

/*... CSRD*/
  if(!strcmp(word,"CSRD"))
   *data = CSRD;
/*... ELLPACK*/
  else if(!strcmp(word,"ELLPACK"))
   *data = ELLPACK;
/*... CSRDCOO*/
  else if(!strcmp(word,"CSRDCOO"))
   *data = CSRDCOO;

} 
/*********************************************************************/      
