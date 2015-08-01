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
 * sistEqX -> estrutura de matriz atualizado ia e ja. 
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
   bool lower,diag,upper;
   
   switch(type){
/*... armazenamento CSR(a)*/
     case CSR:
     break;
/*...................................................................*/

/*... armazenamento CSRD(a,ad)*/
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

/*...*/ 
       HccaAlloc(INT,m,sistEqX->ia     ,sistEqX->neq+1 ,strIa   ,_AD_);
       sistEqX->nad = csrIa(sistEqX->ia ,id ,num    ,nelcon,nViz
                           ,numel       ,sistEqX->neq, maxViz,ndf
                           ,upper       ,diag        ,lower);
/*...................................................................*/

/*...*/
       HccaAlloc(INT,m,sistEqX->ja     ,sistEqX->nad ,strJa   ,_AD_);
       csrJa(sistEqX->ia ,sistEqX->ja ,id ,num    
            ,nelcon,nViz ,numel       ,sistEqX->neq, maxViz,ndf
            ,upper,diag ,lower);
/*...................................................................*/

/*... reordenando o grafo*/
       sortGraphCsr(sistEqX->ia,sistEqX->ja,sistEqX->neq);
/*...................................................................*/

/*... alocacao da matriz*/
       HccaAlloc(DOUBLE,m,sistEqX->ad     ,sistEqX->neq ,strAd   ,_AD_);
       HccaAlloc(DOUBLE,m,sistEqX->al     ,sistEqX->nad ,strA    ,_AD_);
       zero(sistEqX->ad,sistEqX->neq,DOUBLEC);
       zero(sistEqX->al,sistEqX->nad,DOUBLEC);
/*...................................................................*/
     
/*... banda da matriz*/
       printf("band Maxima: %ld\n"
            ,(long) bandCsr(sistEqX->ia,sistEqX->ja,sistEqX->neq,1));
       printf("band Minima: %ld\n"
            ,(long) bandCsr(sistEqX->ia,sistEqX->ja,sistEqX->neq,3));
       printf("band Media : %ld\n"
            ,(long) bandCsr(sistEqX->ia,sistEqX->ja,sistEqX->neq,2));
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
       HccaAlloc(INT,m,sistEqX->ia,2             ,strIa   ,_AD_);
       HccaAlloc(INT,m,sistEqX->ja,sistEqX->neq*maxViz ,strJa   ,_AD_);
       zero(sistEqX->ja,sistEqX->neq*maxViz,INTC);
/*...................................................................*/

/*...*/
       sistEqX->nad = ellPackJa(sistEqX->ia     ,sistEqX->ja 
                              ,id               ,num
                              ,nelcon           ,nViz
                              ,numel            ,sistEqX->neq
                              ,maxViz           ,ndf);
/*...................................................................*/

/*... alocacao da matriz*/
       HccaAlloc(DOUBLE,m,sistEqX->ad,sistEqX->neq ,strAd   ,_AD_);
       HccaAlloc(DOUBLE,m,sistEqX->al
                ,sistEqX->neq*maxViz ,strA    ,_AD_);
       zero(sistEqX->ad,sistEqX->neq,DOUBLEC);
       zero(sistEqX->al,sistEqX->neq*maxViz,DOUBLEC);
/*...................................................................*/

/*... banda da matriz*/
       printf("band Maxima: %ld\n"
            ,(long) bandEllPack(sistEqX->ia,sistEqX->ja,sistEqX->neq,1));
       printf("band Minima: %ld\n"
            ,(long) bandEllPack(sistEqX->ia,sistEqX->ja,sistEqX->neq,3));
       printf("band Media : %ld\n"
            ,(long) bandEllPack(sistEqX->ia,sistEqX->ja,sistEqX->neq,2));
/*...................................................................*/
     break;
/*...................................................................*/
     
     default:
       ERRO_OP(__FILE__,__func__,type);
     break;
/*...................................................................*/
  }
}
/**********************************************************************
 *SETDATASTRUCT : escolhe a estruturade dados                         *
 **********************************************************************/
void setDataStruct(char *word,short *data)
{

  if(!strcmp(word,"CSRD"))
   *data = CSRD;
  else if(!strcmp(word,"ELLPACK"))
   *data = ELLPACK;

} 
/*********************************************************************/      
