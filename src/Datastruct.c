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
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void dataStruct(Memoria *m ,long *id   ,long *num   ,long *nelcon
               ,short *nViz,long numel ,short maxViz
               ,short ndf 
               ,char  *strIa,char *strJa,char *strAd,char *strA  
               ,SistEq *sistEqX)
{
   short type = sistEqX->storage;
   bool lower,diag,upper;
   
   switch(type){
/*... armazenamento CSR/CSRC(a,ad)*/
     case 1:
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
      Myalloc(long,m,sistEqX->ia     ,sistEqX->neq+1 ,strIa   ,_AD_);
      sistEqX->nad = csrIa(sistEqX->ia ,id ,num    ,nelcon,nViz
                          ,numel       ,sistEqX->neq, maxViz, upper
                          ,diag       ,lower);
/*...................................................................*/

/*...*/
      Myalloc(long,m,sistEqX->ja     ,sistEqX->nad ,strJa   ,_AD_);
      csrJa(sistEqX->ia ,sistEqX->ja ,id ,num    
           ,nelcon,nViz ,numel       ,sistEqX->neq, maxViz
           , upper,diag ,lower);
/*...................................................................*/

/*... reordenando o grafo*/
     sortGraphCsr(sistEqX->ia,sistEqX->ja,sistEqX->neq);
/*...................................................................*/

/*...*/
     Myalloc(double,m,sistEqX->ad     ,sistEqX->neq ,strAd   ,_AD_);
     Myalloc(double,m,sistEqX->a      ,sistEqX->nad ,strA    ,_AD_);
     zero(sistEqX->ad,sistEqX->neq,"double");
     zero(sistEqX->a ,sistEqX->nad,"double");

/*...................................................................*/
     break;
/*...................................................................*/

/*...*/
     case 2:
     break;
/*...................................................................*/

/*...*/
     default:
       printf("\n opcao invalida\n"
           "funcao fname(*,*,*)\narquivo = %s\n",__FILE__);
       exit(EXIT_FAILURE);
     break;
/*...................................................................*/
  }
}
