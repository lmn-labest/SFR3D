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
 * au    -> matrix de coeficientes                                   * 
 * ad    -> matrix de coeficientes                                   * 
 * al    -> matrix de coeficientes                                   * 
 * nad   -> numero de termos nao nulos fora da diagonal principal    * 
 * type  -> tipo de armazenamento ( 1 - CSR/CSRC)                    * 
 * unsym -> simetria da matrix                                       * 
 * bin   -> matriz binaria                                           * 
 * name  -> nume do arquivo de saida                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void writeCoo(Memoria *m,INT *ia   ,INT *ja,INT neq
             ,double *au,double *ad,double *al
             ,INT nad  ,short type
             ,bool unsym,bool bin, char *name){
#ifdef _MMIO_  
  MM_typecode matcode;
  int  *lin,*col;
  double *val;
  INT nTotal; 

/*...*/    
  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_coordinate(&matcode);
  mm_set_real(&matcode);
  if(unsym) 
    mm_set_general(&matcode);
  else
    mm_set_symmetric(&matcode);
/*...................................................................*/

/*...*/
  nTotal = neq + nad;
  HccaAlloc(int   ,m,lin,nTotal,"lin",_AD_);
  HccaAlloc(int   ,m,col,nTotal,"col",_AD_);
  HccaAlloc(double,m,val,nTotal,"val",_AD_);
  zero(lin,nTotal,"int");
  zero(col,nTotal,"int");
  zero(val,nTotal,"double");
/*...................................................................*/
  
  switch(type){
/*... CSR/CSRC*/
    case CSRD:
      csrToCoo(lin,col  ,val
              ,ia ,ja   ,au
              ,ad ,al   ,neq
              ,nad,unsym,bin);
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
  }
    HccaDealloc(m,val,"val",false);
    HccaDealloc(m,col,"col",false);
    HccaDealloc(m,lin,"lin",false);
#else
      printf("\nEscrita da matriz no formato MM não disponivel.\n"
           "funcao %s\narquivo = %s\n",__func__,__FILE__);
      exit(EXIT_FAILURE);
#endif
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
 * au  -> matrix de coeficientes                                     * 
 * ad  -> matrix de coeficientes                                     * 
 * al  -> matrix de coeficientes                                     * 
 * neq -> numero de equacoes                                         * 
 * bin -> matriz binaria                                             * 
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
void csrToCoo(int *lin   ,int *col  ,double *val
             ,INT *ia    ,INT *ja   ,double *au
             ,double *ad ,double *al,INT neq   
             ,INT nad    ,bool unsym,bool bin){
  
  INT i,j,kk,nl1,n,ipoint;

  kk = 0;
  for(i=0;i<neq;i++){
/*... diagonal principal*/
    nl1     = i + 1; 
    lin[kk] = col[kk] = nl1;
    if(bin) 
      val[kk] = 1.0;
    else
      val[kk] = ad[i];      
    kk++;
/*...................................................................*/
    n  = ia[i+1] - ia[i];
    ipoint = ia[i];
    for(j=0;j<n ;j++){
      lin[kk]  = nl1;
      col[kk]  = ja[ipoint+j]+1;
      if(bin)
        val[kk]  = 1.0;
      else{
        if( col[kk] < nl1 ) 
          val[kk] = al[ipoint+j];      
        if( col[kk] > nl1 && unsym )
          val[kk]  = au[ipoint+j];
      }
      kk++;
    }
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * WRITECOOB: escreve o vetor de forcao no formato COO               * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * b   -> vetor de forcas                                            * 
 * neq -> numero de equacoes                                         * 
 * name-> nome do arquivo de saida                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
  void writeCooB(double *b,INT const neq,char *name){
 
  FILE *fileOut;
  int i;
   
  fileOut = openFile(name,"w");

  fprintf(fileOut,"%%MatrixMarket matrix array real general\n");
  fprintf(fileOut,"%d %d\n",neq,1);
  for(i=0;i<neq;i++)
    fprintf(fileOut,"%.14e\n",b[i]);

}
/*********************************************************************/ 

