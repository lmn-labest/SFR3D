#include<Coo.h>
/********************************************************************* 
 * WRITECOO : escreve o grafo da matrix no formato MM (MATRIX MARKET)* 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m     -> memoria principal                                        * 
 * ia    -> vetor da matrix                                          * 
 * ja    -> vetor da matrix                                          * 
 * neq   -> numero de equacoes                                       * 
 * nad   -> numero de termos nao nulos fora da diagonal principal    * 
 * type  -> tipo de armazenamento ( 1 - CSR/CSRC)                    * 
 * unsym -> simetria da matrix                                       * 
 * name  -> nume do arquivo de saida                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void writeCoo(Memoria *m,long *ia   ,long *ja,long neq
             ,long nad  ,short type
             ,bool unsym,char *name){
  
  MM_typecode matcode;
  int  *lin,*col;
  double *val;
  long nTotal; 
  int i;

/*...*/    
  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_coordinate(&matcode);
  mm_set_real(&matcode);
  if(unsym) 
    mm_set_general(&matcode);
  else
    mm_set_symmetric(&matcode);
  printf("%s %ld %ld\n",name,neq,nad);
/*...................................................................*/

/*...*/
  nTotal = neq + nad;
  Myalloc(int   ,m,lin,nTotal,"lin",_AD_);
  Myalloc(int   ,m,col,nTotal,"col",_AD_);
  Myalloc(double,m,val,nTotal,"val",_AD_);
  zero(lin,nTotal,"int");
  zero(col,nTotal,"int");
  zero(val,nTotal,"double");
/*...................................................................*/
  
  switch(type){
/*... CSR/CSRC*/
    case 1:
      csrToCoo(lin,col,val,ia,ja,neq,nad);
      for(i=0;i<nTotal;i++)
        printf("%d %d %d %lf\n",i,lin[i],col[i],val[i]);
      mm_write_mtx_crd(name,neq,neq,nTotal,lin,col,val,matcode);
    break;
/*...................................................................*/

/*...*/
    default:
      printf("\n opcao invalida\n"
           "funcao fname(*,*,*)\narquivo = %s\n",__FILE__);
      exit(EXIT_FAILURE);
    break;
/*...................................................................*/
    Mydealloc(m,val,"val",false);
    Mydealloc(m,col,"col",false);
    Mydealloc(m,lin,"lin",false);
  }

}
/*********************************************************************/ 

/********************************************************************* 
 * CSRTOCOO : conveter do formato CSR para COO                       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lin -> indefinido                                                 * 
 * col -> indefinido                                                 * 
 * val -> indefinido                                                 * 
 * ia  -> vetor CSR                                                  * 
 * ja  -> vetor CSR                                                  * 
 * neq -> numero de equacoes                                         * 
 * nad -> numero de termos nao nulos                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lin -> numero da linha                                            * 
 * col -> numero da coluna                                           * 
 * val -> valor                                                      * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void csrToCoo(int *lin, int *col,double *val,long *ia,long *ja
             ,long neq,long nad){
  
  long i,j,kk,nl1,n,ipoint;

  kk = 0;
  for(i=0;i<neq;i++){
/*... diagonal principal*/
    nl1     = i + 1; 
    lin[kk] = col[kk] = nl1;
    val[kk] = 1.0;
    kk++;
/*...................................................................*/
    n  = ia[i+1] - ia[i];
    ipoint = ia[i];
    for(j=0;j<n ;j++){
      lin[kk] = nl1;
      col[kk]  = ja[ipoint+j]+1;
      val[kk]  = 1.0;
      kk++;
    }
  }
}
/*********************************************************************/ 

