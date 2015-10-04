#include<Csr.h>
/********************************************************************* 
 * CSRIA: Geracao do vetor ia do CSR/CSRC                            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia     -> indefinido                                              * 
 * id     -> numeracao das equacoes por elemento                     * 
 * num    -> renumeracao dos elementos                               * 
 * adj    -> adjacencia dos elementos                                * 
 * nViz   -> numero de vizinhos por elemento                         * 
 * numel  -> numero de elementos                                     * 
 * neq    -> numero de equacoes                                      * 
 * ndf    -> numero de graus de liberade                             * 
 * maxViz -> numero maximo de vizinho da malha                       * 
 * upper  -> armazenamento da parte superior da matriz (CSR/CSRC)    * 
 * diag   -> armazenamento da diagonal (CSR/CSRC)                    * 
 * lower  -> armazenamenro da parte inferior (CSR/CSRC)              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ia     -> ponteiro do CSR                                         * 
 *-------------------------------------------------------------------* 
 * OBS: a funcao retorna o numero do termos nao nulor no CSR/CSRC    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
INT csrIa(INT *restrict ia     ,INT *restrict id   
         ,INT *restrict num    ,INT *restrict adj
         ,short *restrict nViz
         ,INT const numel      ,INT const neq
         ,short const maxViz   ,short  const ndf
         ,bool const upper     ,bool const diag   
         ,bool const lower     ){
  
  INT  i,nel1,neq1,neq2,viz1,col,aux;
  short jNdf,kNdf,j;
/*... gerando arranjo ia*/  
  ia[0] = 0;
  for(i=0;i<numel;i++){
    nel1= num[i]-1;
    for(jNdf=0;jNdf<ndf;jNdf++){
      aux = 0;
      neq1 = MAT2D(nel1,jNdf,id,ndf)-1;
      if(neq1 != -2){
/*... equacoes ligadas ao proprio elemento
        (mais de um grau de liberdade por celula)*/
        for(kNdf=0;kNdf<ndf;kNdf++){
          neq2 = MAT2D(nel1,kNdf,id,ndf)-1;
          if(neq2 != -2){
/*... parte superior*/
            if(lower && neq1 > neq2) 
              aux++;
/*... parte inferior*/            
            else if(upper && neq1 < neq2)
              aux++;
/*... diagonal princial*/      
            else if(diag && neq1 == neq2 ) 
              aux++;
          }
        }
/*...................................................................*/
  
/*... equacoes ligadas as celulas vizinhas*/
        for(j=0;j<nViz[nel1];j++){
          viz1 = MAT2D(nel1,j,adj,maxViz) - 1;
          if( viz1 != -2) {
            for(kNdf=0;kNdf<ndf;kNdf++){
              col   = MAT2D(viz1,kNdf,id,ndf)-1;
              if( col != -2){
/*... parte inferior*/
                if(lower && col < neq1) 
                  aux++;
/*...................................................................*/

/*... parte superior*/            
                else if(upper && col > neq1)
                  aux++;
              }
/*...................................................................*/
            }    
/*...................................................................*/
          }
/*...................................................................*/
        }
/*...................................................................*/
        ia[neq1+1] = ia[neq1] + aux;
      }
/*...................................................................*/
    }
  }
/*...................................................................*/
  return ia[neq] - ia[0];
}
/*********************************************************************/ 

/********************************************************************* 
 * CSRIAR: geracao do vetor ia da parte retangular CSRC-MPI          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia     -> indefinido                                              * 
 * id     -> numeracao das equacoes por elemento                     * 
 * num    -> renumeracao dos elementos                               * 
 * adj    -> adjacencia dos elementos                                * 
 * nViz   -> numero de vizinhos por elemento                         * 
 * numel  -> numero de elementos                                     * 
 * neq    -> numero de equacoes                                      * 
 * ndf    -> numero de graus de liberade                             * 
 * maxViz -> numero maximo de vizinho da malha                       * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ia     -> ponteiro do CSRC-CSR                                    * 
 *-------------------------------------------------------------------* 
 * OBS: a funcao retorna o numero do termos nao nulos no CSR         * 
 * retangular                                                        * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
INT csrIaR(INT *restrict ia     ,INT *restrict id   
          ,INT *restrict num    ,INT *restrict adj
          ,short *restrict nViz
          ,INT const numel      ,INT const nEq
          ,short const maxViz   ,short  const ndf){
  
  INT  i,nel1,neq,neq1,viz1,col,aux,nad=0;
  short jNdf,kNdf,j;
/*... gerando arranjo ia*/  
  ia[0] = 0;
  neq = nEq-1;
  for(i=0;i<numel;i++){
    nel1= num[i]-1;
    for(jNdf=0;jNdf<ndf;jNdf++){
      aux = 0;
      neq1 = MAT2D(nel1,jNdf,id,ndf)-1;
      if(neq1 != -2){
/*... equacoes ligadas as celulas vizinhas*/
        for(j=0;j<nViz[nel1];j++){
          viz1 = MAT2D(nel1,j,adj,maxViz) - 1;
          if( viz1 != -2) {
            for(kNdf=0;kNdf<ndf;kNdf++){
              col   = MAT2D(viz1,kNdf,id,ndf)-1;
              if( col != -2){
/*... parte superior*/
                if(col > neq){ 
                  aux++;
                  nad++;
                }
/*...................................................................*/
              }
/*...................................................................*/
            }    
/*...................................................................*/
          }
/*...................................................................*/
        }
/*...................................................................*/
        ia[neq1+1] = ia[neq1] + aux;
      }
/*...................................................................*/
    }
  }
/*...................................................................*/
  return nad;
}
/*********************************************************************/ 

/********************************************************************* 
 * CSRJA:                                                            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia     -> arranjo CSR/CSRC                                        * 
 * ja     -> indefinido                                              * 
 * id     -> numeracao das equacoes por elemento                     * 
 * num    -> renumeracao dos elementos                               * 
 * adj    -> adjacencia dos elementos                                * 
 * nViz   -> numero de vizinhos por elemento                         * 
 * numel  -> numero de elementos                                     * 
 * neq    -> numero de equacoes                                      * 
 * maxViz -> numero maximo de vizinho da malha                       * 
 * ndf    -> numero de graus de liberdade                            * 
 * upper  -> armazenamento da parte superior da matriz (CSR/CSRC)    * 
 * diag   -> armazenamento da diagonal (CSR/CSRC)                    * 
 * lower  -> armazenamenro da parte inferior (CSR/CSRC)              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ja     -> ponteiro do CSR                                         * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void csrJa(INT *restrict ia    ,INT *restrict ja 
          ,INT *restrict id    ,INT *restrict num 
          ,INT *restrict adj   ,short *restrict nViz
          ,INT const numel     ,INT const neq 
          ,short const maxViz  ,short const ndf
          ,bool const upper    ,bool const diag 
          ,bool const lower){
  
  INT  i,nel1,neq1,neq2,viz1,col,aux,ipont;
  short j,jNdf,kNdf;

/*... gerando arranjo ja*/  
  for(i=0;i<numel;i++){
    nel1= num[i]-1;
    for(jNdf=0;jNdf<ndf;jNdf++){
      aux = 0;
      neq1= MAT2D(nel1,jNdf,id,ndf)-1;
      ipont = ia[neq1];
      if(neq1 != -2){
/*... equacoes ligadas ao proprio elemento
        (mais de um grau de liberdade por celula)*/
        for(kNdf=0;kNdf<ndf;kNdf++){
          neq2 = MAT2D(nel1,kNdf,id,ndf)-1;
          if(neq2 != -2){
/*... parte inferior*/
            if(lower && neq1 > neq2){ 
              ja[ipont+aux] = neq2; 
              aux++;
            }
/*... parte superior*/            
            else if(upper && neq1 < neq2){
              ja[ipont+aux] = neq2; 
              aux++;
            }
/*... diagonal princial*/      
            else if(diag && neq1 == neq2 ){ 
              ja[ipont+aux] = neq1;  
              aux++;
            }
/*...................................................................*/
          }
/*...................................................................*/
        }
/*...................................................................*/

/*... equacoes ligadas as celulas vizinhas*/
        for(j=0;j<nViz[nel1];j++){
          viz1 = MAT2D(nel1,j,adj,maxViz) - 1;
          if( viz1 != -2) {
            for(kNdf=0;kNdf<ndf;kNdf++){
              col= MAT2D(viz1,kNdf,id,ndf)-1;
              if( col != -2){
/*... parte inferior*/
                if(lower && col < neq1){
                  ja[ipont+aux] = col; 
                  aux++;
                }
/*...................................................................*/

/*... parte superior*/            
                else if(upper && col > neq1){
                  ja[ipont+aux] = col; 
                  aux++;
                }
/*...................................................................*/
              }
/*...................................................................*/
            }
/*...................................................................*/
          }    
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * CSRJAR: geracao do vetor ja da parte retangular CSRC-MPI          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia     -> arranjo CSR/CSRC                                        * 
 * ja     -> indefinido                                              * 
 * id     -> numeracao das equacoes por elemento                     * 
 * num    -> renumeracao dos elementos                               * 
 * adj    -> adjacencia dos elementos                                * 
 * nViz   -> numero de vizinhos por elemento                         * 
 * numel  -> numero de elementos                                     * 
 * neq    -> numero de equacoes                                      * 
 * maxViz -> numero maximo de vizinho da malha                       * 
 * ndf    -> numero de graus de liberdade                            * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ja     -> ponteiro do CSR                                         * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void csrJaR(INT *restrict ia    ,INT *restrict ja 
           ,INT *restrict id    ,INT *restrict num 
           ,INT *restrict adj   ,short *restrict nViz
           ,INT const numel     ,INT const nEq 
           ,short const maxViz  ,short const ndf){
  
  INT  i,nel1,neq,neq1,viz1,col,aux,ipont;
  short j,jNdf,kNdf;

  neq = nEq-1;
/*... gerando arranjo ja*/  
  for(i=0;i<numel;i++){
    nel1= num[i]-1;
    for(jNdf=0;jNdf<ndf;jNdf++){
      aux = 0;
      neq1= MAT2D(nel1,jNdf,id,ndf)-1;
      if(neq1 != -2){
        ipont = ia[neq1];
/*... equacoes ligadas as celulas vizinhas*/
        for(j=0;j<nViz[nel1];j++){
          viz1 = MAT2D(nel1,j,adj,maxViz) - 1;
          if( viz1 != -2) {
            for(kNdf=0;kNdf<ndf;kNdf++){
              col= MAT2D(viz1,kNdf,id,ndf)-1;
              if( col != -2){
/*... parte inferior*/
                if(col > neq){
                  ja[ipont+aux] = col; 
                  aux++;
                }
/*...................................................................*/

              }
/*...................................................................*/
            }
/*...................................................................*/
          }    
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/ 


/********************************************************************* 
 * BANDCSR: banda da matriz no formato CSR                           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ia  - ponteiro CSR                                                * 
 * ja  - ponteiro CSR                                                * 
 * neq - numero de equacoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS: retorna a banda da matrix                                    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
INT bandCsr(INT *ia,INT *ja,INT  neq,short type){

  INT i,j,aux;
  unsigned long bandL=0;  

  switch(type){
/*... banda maxima da matriz*/
    case BANDCSRMAX:
      for(i=0;i<neq;i++){
        for(j=ia[i];j<ia[i+1];j++){
          bandL = max(bandL,abs(i-ja[j]));
        }
      }
    break;
/*...................................................................*/ 

/*... banda media da matriz*/
    case BANDCSRMED:
      for(i=0;i<neq;i++){
        aux = 0;
        for(j=ia[i];j<ia[i+1];j++){
          aux = max(aux,abs(i-ja[j]));
        }
        bandL += aux;
      }
      bandL = bandL/neq;
    break;
/*...................................................................*/ 

/*... banda minima da matriz*/
    case BANDCSRMIN:
      bandL = neq;
      for(i=0;i<neq;i++){
        for(j=ia[i];j<ia[i+1];j++){
          bandL = min(bandL,abs(i-ja[j]));
        }
      }
    break;
/*...................................................................*/ 
  }
  return bandL;

}
/*********************************************************************/ 

/********************************************************************* 
 * CSR : Montagem do sistema global do CSR                           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * ia      -> ponteiro para as linhas da matriz esparsa              * 
 * ja      -> ponteiro para as colunas da matriz esparsa             * 
 * au      -> matriz de coeficientes esparsa                         * 
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- triangular superior               )       *
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 * al      -> matriz de coeficientes esparsa                         * 
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- triangular inferior               )       *
 * b       -> vetor de forcas                                        *
 * lId     -> numeracao das equacoes dessa celula                    *
 * lA      -> coficientes da celula                                  *
 * lB      -> vetor de forca da celula                               *
 * nEq     -> numero de equacoes                                     *
 * nAd     -> numero de termos nao nulos                             *
 * nAdR    -> numero de termos nao nulos na parte retangular         *
 * nFace   -> numero de faces da celula                              *
 * ndf     -> graus de liberdade                                     *
 * storage -> tecnica de armazenamento da matriz esparsa             * 
 * forces  -> mantagem no vetor de forcas                            * 
 * matrix  -> mantagem da matriz de coeficientes                     * 
 * unsym   -> matiz nao simetrica                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * au,a,al   -> coeficiente da linha i     (matriz = true)           *
 * b         -> vetor de forca da linha i  (forces = true)           *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 * lA | AP1X AP1Y AP1XY AP1YX|                                       * 
 *    | AP2X AP2Y AP2XY AP2YX|                                       * 
 *    | AP3X AP3Y AP3XY AP3YX|                                       * 
 *    |          ...         |                                       * 
 *    | APX  APY  APXY  APYX |                                       * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void csr(INT    *restrict  ia,INT *restrict ja 
        ,DOUBLE *restrict au ,DOUBLE *restrict ad
        ,DOUBLE *restrict al ,DOUBLE *restrict b
        ,INT *restrict lId                       
        ,DOUBLE *restrict lA ,DOUBLE *restrict lB
        ,INT const nEq       ,INT const nAd 
        ,INT const nAdR                        
        ,short const nFace   ,short const ndf  
        ,short const storage ,bool  const forces
        ,bool const matrix   ,bool  const  unsym)
{
  INT lNeq,lCol=0,iak,jak,iPont,iaKneq,neqS;
  INT *iar,*jar;
  DOUBLE *ar;
  unsigned short i,j,k,jLa,nst;
  bool fCoo = false;

  neqS = nEq - 1;

  switch (storage){
/*... estrutura CSR(ia,ja,a,b)*/
    case CSR:
    break;
/*...................................................................*/

/*... estrutura CSRD(ia,ja,a,al,b)*/
/*...CSRD+COO(symetric)*/
    case CSRDCOO:
      fCoo = true;
    case CSRD:
      if(ndf == 1) {
        lNeq = lId[nFace];
/*... vetor de forcas*/
        if(forces) 
          b[lNeq] = lB[0];
/*...................................................................*/

/*... matriz de coefiencientes*/
        if(matrix){ 
          lNeq     = lId[nFace];
          ad[lNeq] = lA[nFace];
          iPont    = ia[lNeq];
          iaKneq   = ia[lNeq+1];
          for(k=0;k<nFace;k++){
            lCol     = lId[k];
            if(lCol > -1){
              for(iak = iPont; iak < iaKneq;iak++){
                jak = ja[iak];
                if( lCol == jak ) 
                  al[iak] = lA[k];  
              }
            }      
          }
/*...................................................................*/

/*... parte retagunlar para matrizes simetricas (MPI)*/
          if(unsym == false && mpiVar.nPrcs > 1){
/*... parte retangular em COO*/
            if(fCoo){
              iar = &ia[nEq+1];
              jar = &ja[nAd];
              ar  = &al[nAd];
              for(k=0;k<nFace;k++){
                lCol     = lId[k];
                if(lCol > -1 && (lCol > neqS))
                  for(iak = 0; iak < nAdR;iak++)
                    if(iar[iak] == lNeq && jar[iak] == lCol)
                      ar[iak] = lA[k];     
              }
            }   
/*...................................................................*/

/*... parte retangular em CSR*/
            else{
              iar = &ia[nEq+1];
              jar = &ja[nAd];
              ar  = &al[nAd];
              iPont    = iar[lNeq];
              iaKneq   = iar[lNeq+1];
              for(k=0;k<nFace;k++){
                lCol     = lId[k];
                if(lCol > -1 && (lCol > neqS)){
                  for(iak = iPont; iak < iaKneq;iak++){
                    jak = jar[iak];
                    if( lCol == jak ){ 
                      ar[iak] = lA[k];
                    }  
                  }
                }       
              }
            }
/*...................................................................*/
          }  
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/

/*...*/
      else{
        nst  = (nFace+1)*ndf;
/*... vetor de forcas*/
        if(forces) 
          for(j=0;j<ndf;j++){
            lNeq = MAT2D(nFace,j,lId,ndf);
            b[lNeq] = lB[j];
          }
/*...................................................................*/

/*... matriz de coefiencientes*/
        if(matrix) 
          for(i=0;i<ndf;i++){
            jLa      = nFace*ndf + i;
            lNeq     = MAT2D(nFace,i,lId,ndf);
            ad[lNeq] = MAT2D(i,jLa,lA,nst);
            iPont    = ia[lNeq];
            iaKneq   = ia[lNeq+1];
            for(k=0;k<nFace;k++){
              lCol     = MAT2D(k,i,lId,ndf);
              if(lCol > -1)
                for(iak = iPont; iak < iaKneq;iak++){
                  jak = ja[iak];
                  if( lCol == jak ) 
                    for(j=0;j<ndf;j++){
                      jLa = k*ndf+j;
                      al[iak+j] = MAT2D(j,jLa,lA,nst);  
                    }      
                }
            }   
          }
      }
/*...................................................................*/
    break;
/*...................................................................*/

/*... estrutura CSR(ia,ja,au,a,al,b)*/
    case CSRC:
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
/*********************************************************************/ 
