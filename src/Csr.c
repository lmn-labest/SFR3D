#include<Csr.h>
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
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
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
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
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
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
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
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
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
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
  INT bandL=0;  

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
 * Data de criacao    : 16/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * BANDCSRC: banda da matriz no formato CSRC                         * 
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
INT bandCsrC(INT *ia,INT *ja,INT  neq,short type){

  INT i,j,aux;
  INT bandL=0,col=0;
  LONG_INT bandLong=0;

  switch(type){
/*... banda maxima da matriz*/
    case BANDCSRMAX:
/*... loop nas linhas(parte inferior)*/
      for(i=0;i<neq;i++){
        for(j=ia[i];j<ia[i+1];j++){
          bandL = max(bandL,abs(i-ja[j]));
        }
      }
/*... loop nas colunas(parte superior)*/
      for(i=0;i<neq;i++){
        for(j=ia[i];j<ia[i+1];j++){
          if( ja[j] == i)
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
        bandLong += aux;
      }
      bandL = bandLong/neq;
    break;
/*...................................................................*/ 

/*... banda minima da matriz*/
    case BANDCSRMIN:
      bandL = neq;
      col   = neq;
/*... loop nas linhas(parte inferior)*/
      for(i=0;i<neq;i++){
        for(j=ia[i];j<ia[i+1];j++){
          col = min(col,ja[j]);
        }
/*... loop nas colunas(parte superior)*/
        for(j=ia[i];j<ia[i+1];j++){
          if( ja[j] == i)
            col = min(col,ja[j]);
        }
        bandL = min(bandL,abs(i-col));
      }
    break;
/*...................................................................*/ 
  }
  return bandL;

}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
*-------------------------------------------------------------------* 
 * CSR : Montagem do sistema global do CSR                           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * ia      -> ponteiro para as linhas da matriz esparsa              * 
 * ja      -> ponteiro para as colunas da matriz esparsa             * 
 * a       -> matriz de coeficientes esparsa                         * 
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- triangular inferior/triangular superior ) *
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 * b       -> vetor de forcas                                        *
 * lId     -> numeracao das equacoes dessa celula                    *
 * lA      -> coficientes da celula                                  *
 * lB      -> vetor de forca da celula                               *
 * nEq     -> numero de equacoes                                     *
 * neqNov  -> numero de equacoes nao sobrepostas                     *
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
        ,DOUBLE *restrict a  ,DOUBLE *restrict ad
        ,DOUBLE *restrict b
        ,INT *restrict lId                       
        ,DOUBLE *restrict lA ,DOUBLE *restrict lB
        ,INT const nEq       ,INT const nEqNov 
        ,INT const nAd       ,INT const nAdR                        
        ,short const nFace   ,short const ndf  
        ,short const storage ,bool  const forces
        ,bool const matrix   ,bool  const  unsym)
{
  INT lNeq,lCol=0,iak,jak,iPoint,iaKneq,jPoint,iaJneq,neqS,n0;
  INT *iar,*jar;
  DOUBLE *restrict ar=NULL;
  DOUBLE *restrict au=NULL;
  DOUBLE *restrict al=NULL;
  unsigned short i,j,k,jLa,nst;
  bool fCoo = false;

  neqS = nEqNov - 1;

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
          al       = a;
          lNeq     = lId[nFace];
          ad[lNeq] = lA[nFace];
          iPoint   = ia[lNeq];
          iaKneq   = ia[lNeq+1];
          for(k=0;k<nFace;k++){
            lCol     = lId[k];
            if(lCol > -1){
              for(iak = iPoint; iak < iaKneq;iak++){
                jak = ja[iak];
                if( lCol == jak ) 
                  al[iak] = lA[k];  
              }
/*...................................................................*/
            }      
/*...................................................................*/
          }
/*...................................................................*/

/*... parte retagunlar para matrizes simetricas (MPI)*/
          if(unsym == false && mpiVar.nPrcs > 1){
/*... parte retangular em COO*/
            if(fCoo){
              n0     = nEqNov+1;
              iar    = &ia[n0];
              jar    = &ja[nAd];
              ar     = &a[nAd];
              iPoint = iar[lNeq];
              iaKneq = iar[lNeq+1];
              for(k=0;k<nFace;k++){
                lCol     = lId[k];
                if(lCol > -1 && lCol > neqS)
                  for(iak = iPoint; iak < iaKneq;iak++)
                    if(jar[iak] == lCol)
                      ar[iak] = lA[k];     
              }
/*...................................................................*/
            }   
/*...................................................................*/

/*... parte retangular em CSR*/
            else{
              iar    = &ia[nEqNov+1];
              jar    = &ja[nAd];
              ar     = &a[nAd];
              iPoint = iar[lNeq];
              iaKneq = iar[lNeq+1];
              for(k=0;k<nFace;k++){
                lCol     = lId[k];
                if(lCol > -1 && lCol > neqS){
                  for(iak = iPoint; iak < iaKneq;iak++){
                    jak = jar[iak];
                    if( lCol == jak ){ 
                      ar[iak] = lA[k];
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
          al  = a;
          for(i=0;i<ndf;i++){
            jLa      = nFace*ndf + i;
            lNeq     = MAT2D(nFace,i,lId,ndf);
            ad[lNeq] = MAT2D(i,jLa,lA,nst);
            iPoint   = ia[lNeq];
            iaKneq   = ia[lNeq+1];
            for(k=0;k<nFace;k++){
              lCol     = MAT2D(k,i,lId,ndf);
              if(lCol > -1)
                for(iak = iPoint; iak < iaKneq;iak++){
                  jak = ja[iak];
                  if( lCol == jak ) 
                    for(j=0;j<ndf;j++){
                      jLa = k*ndf+j;
                      al[iak+j] = MAT2D(j,jLa,lA,nst);  
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
    break;
/*...................................................................*/

/*... estrutura CSRC(ia,ja,au,a,al,b)*/
/*...CSRC+COO*/
    case CSRCCOO:
      fCoo = true;
    case CSRC:
/*...*/
      if(ndf == 1) {
        lNeq = lId[nFace];
/*... vetor de forcas*/
        if(forces) 
          b[lNeq] = lB[0];
/*...................................................................*/

/*... matriz de coefiencientes*/
        if(matrix){ 
          al       = a;
          lNeq     = lId[nFace];
		      ad[lNeq] = lA[nFace];
          iPoint   = ia[lNeq];
          iaKneq   = ia[lNeq+1];
/*... paralelo MPI*/                                                        
          if(mpiVar.nPrcs > 1){
/*... nao simetrico*/
            if(unsym){
              au = &a[nAd];
              for(k=0;k<nFace;k++){
                lCol     = lId[k];
                if(lCol > -1){
/*... parte inferior*/
                  if( lCol < lNeq)
                    for(iak = iPoint; iak < iaKneq;iak++){
                      jak = ja[iak];
                      if( lCol == jak ) 
                        al[iak] = lA[k];
                    }
/*... parte superior*/
                  else if( lCol > lNeq && lCol <= neqS){
                    jPoint   = ia[lCol];
                    iaJneq   = ia[lCol+1];
                    for(iak = jPoint; iak < iaJneq;iak++){
                      jak = ja[iak];
                      if( lNeq == jak )
                        au[iak] = lA[k];
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

/*... simetrico*/
            else{
              for(k=0;k<nFace;k++){
                lCol     = lId[k];
                if(lCol > -1){
                  for(iak = iPoint; iak < iaKneq;iak++){
                    jak = ja[iak];
                    if( lCol == jak ) 
                      al[iak] = lA[k];  
                  }
/*...................................................................*/
                }      
/*...................................................................*/
              }
/*...................................................................*/
            }
/*...................................................................*/

/*... parte retangular em COO*/
            if(fCoo){
              n0     = nEqNov+1;
              iar    = &ia[n0];
              jar    = &ja[nAd];
/*...*/
              if(unsym)
                ar     = &a[2*nAd];
              else
                ar     = &a[nAd];
/*...................................................................*/
              iPoint = iar[lNeq];
              iaKneq = iar[lNeq+1];
              for(k=0;k<nFace;k++){
                lCol     = lId[k];
                if(lCol > -1 && lCol > neqS)
                  for(iak = iPoint; iak < iaKneq;iak++)
                    if(jar[iak] == lCol) 
                      ar[iak] = lA[k]; 
              }
/*...................................................................*/
            }   
/*...................................................................*/

/*... parte retangular em CSR*/
            else{
              n0     = nEqNov+1;
              iar    = &ia[n0];
              jar    = &ja[nAd];
/*...*/
              if(unsym)
                ar     = &a[2*nAd];
              else
                ar     = &a[nAd];
/*...................................................................*/
              iPoint = iar[lNeq];
              iaKneq = iar[lNeq+1];
              for(k=0;k<nFace;k++){
                lCol     = lId[k];
                if(lCol > -1 && lCol > neqS){
                  for(iak = iPoint; iak < iaKneq;iak++){
                    jak = jar[iak];
                    if( lCol == jak ) 
                      ar[iak] = lA[k];
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

/*... sequencial */                                                        
          else{
/*... nao simetrico*/
            if(unsym){
              au = &a[nAd];
              for(k=0;k<nFace;k++){
                lCol     = lId[k];
                if(lCol > -1){
/*... parte inferior*/
                  if( lCol < lNeq)
                    for(iak = iPoint; iak < iaKneq;iak++){
                      jak = ja[iak];
                      if( lCol == jak ) 
                        al[iak] = lA[k];
                    }
/*... parte superior*/
                  else if( lCol > lNeq){
                    jPoint   = ia[lCol];
                    iaJneq   = ia[lCol+1];
                    for(iak = jPoint; iak < iaJneq;iak++){
                      jak = ja[iak];
                      if( lNeq == jak ) {
                        au[iak] = lA[k];
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

/*... simetrico*/
            else{
              for(k=0;k<nFace;k++){
                lCol     = lId[k];
                if(lCol > -1){
                  for(iak = iPoint; iak < iaKneq;iak++){
                    jak = ja[iak];
                    if( lCol == jak ) 
                      al[iak] = lA[k];  
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
        if(matrix){
          for(i=0;i<ndf;i++){
            jLa      = nFace*ndf + i;
            lNeq     = MAT2D(nFace,i,lId,ndf);
            ad[lNeq] = MAT2D(i,jLa,lA,nst);
            iPoint   = ia[lNeq];
            iaKneq   = ia[lNeq+1];
            for(k=0;k<nFace;k++){
              lCol     = MAT2D(k,i,lId,ndf);
              if(lCol > -1)
                for(iak = iPoint; iak < iaKneq;iak++){
                  jak = ja[iak];
                  if( lCol == jak ) 
                    for(j=0;j<ndf;j++){
                      jLa = k*ndf+j;
                      a[iak+j] = MAT2D(j,jLa,lA,nst);  
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

/********************************************************************* 
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * CSRSIMPLE : Montagem do sistema global do CSR   para sistems      *
 * de velocidades do metodo simple                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * ia      -> ponteiro para as linhas da matriz esparsa              * 
 * ja      -> ponteiro para as colunas da matriz esparsa             * 
 * a       -> matriz de coeficientes esparsa                         * 
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- triangular inferior/triangular superior ) *
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 * b       -> vetor de forcas                                        *
 * lId     -> numeracao das equacoes dessa celula                    *
 * lA      -> coficientes da celula                                  *
 * lB      -> vetor de forca da celula                               *
 * nEq     -> numero de equacoes                                     *
 * neqNov  -> numero de equacoes nao sobrepostas                     *
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
 c b = | b1x b2x ... b(neq)x b1y b2y ... b(neq)y b1z b2z ... b(neq)z * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void csrSimple(INT    *restrict  ia,INT *restrict ja 
        ,DOUBLE *restrict a  ,DOUBLE *restrict ad
        ,DOUBLE *restrict b
        ,INT *restrict lId                       
        ,DOUBLE *restrict lA ,DOUBLE *restrict lB
        ,INT const nEq       ,INT const nEqNov 
        ,INT const nAd       ,INT const nAdR                        
        ,short const nFace   ,short const ndf  
        ,short const storage ,bool  const forces
        ,bool const matrix   ,bool  const  unsym)
{
  INT lNeq,lCol=0,iak,jak,iPoint,iaKneq,jPoint,iaJneq,neqS,n0;
  INT *iar,*jar;
  DOUBLE *restrict ar=NULL;
  DOUBLE *restrict au=NULL;
  DOUBLE *restrict al=NULL;
  unsigned short k;
  bool fCoo = false;

  neqS = nEqNov - 1;

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
      lNeq = lId[nFace];
/*... vetor de forcas*/
      if(forces){
        b[lNeq]     = lB[0];
        b[nEq+lNeq] = lB[1];
        if( ndf == 3) b[2*nEq+lNeq] = lB[2];
      }
/*...................................................................*/

/*... matriz de coefiencientes*/
      if(matrix){ 
        al       = a;
        lNeq     = lId[nFace];
/*...*/
        ad[lNeq]       = lA[nFace];
     		ad[nEq+lNeq]   = lA[nFace+1];
    	  if( ndf == 3) ad[2*nEq+lNeq] = lA[nFace+2];
/*...................................................................*/
        iPoint   = ia[lNeq];
        iaKneq   = ia[lNeq+1];
        for(k=0;k<nFace;k++){
          lCol     = lId[k];
          if(lCol > -1){
            for(iak = iPoint; iak < iaKneq;iak++){
              jak = ja[iak];
              if( lCol == jak ) 
                al[iak] = lA[k];  
            }
/*...................................................................*/
          }      
/*...................................................................*/
        }
/*...................................................................*/

/*... parte retagunlar para matrizes simetricas (MPI)*/
        if(unsym == false && mpiVar.nPrcs > 1){
/*... parte retangular em COO*/
          if(fCoo){
            n0     = nEqNov+1;
            iar    = &ia[n0];
            jar    = &ja[nAd];
            ar     = &a[nAd];
            iPoint = iar[lNeq];
            iaKneq = iar[lNeq+1];
            for(k=0;k<nFace;k++){
              lCol     = lId[k];
              if(lCol > -1 && lCol > neqS)
                for(iak = iPoint; iak < iaKneq;iak++)
                  if(jar[iak] == lCol)
                    ar[iak] = lA[k];     
            }
/*...................................................................*/
          }   
/*...................................................................*/

/*... parte retangular em CSR*/
          else{
            iar    = &ia[nEqNov+1];
            jar    = &ja[nAd];
            ar     = &a[nAd];
            iPoint = iar[lNeq];
            iaKneq = iar[lNeq+1];
            for(k=0;k<nFace;k++){
              lCol     = lId[k];
              if(lCol > -1 && lCol > neqS){
                for(iak = iPoint; iak < iaKneq;iak++){
                  jak = jar[iak];
                  if( lCol == jak ) 
                    ar[iak] = lA[k];
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
    break;
/*...................................................................*/

/*... estrutura CSRC(ia,ja,au,a,al,b)*/
/*...CSRC+COO*/
    case CSRCCOO:
      fCoo = true;
    case CSRC:
/*...*/
      lNeq = lId[nFace];
/*... vetor de forcas*/
      if(forces){ 
        b[lNeq]     = lB[0];
        b[nEq+lNeq] = lB[1];
        if( ndf == 3){ 
          b[2*nEq+lNeq] = lB[2];
        }
      }
/*...................................................................*/

/*... matriz de coefiencientes*/
      if(matrix){ 
        al       = a;
        lNeq     = lId[nFace];
/*...*/
				ad[lNeq]     = lA[nFace];
  		  ad[nEq+lNeq] = lA[nFace+1];
  		  if (ndf == 3) ad[2*nEq+lNeq] = lA[nFace+2];
/*...................................................................*/
        iPoint   = ia[lNeq];
        iaKneq   = ia[lNeq+1];
/*... paralelo MPI*/                                                        
        if(mpiVar.nPrcs > 1){
/*... nao simetrico*/
          if(unsym){
            au = &a[nAd];
            for(k=0;k<nFace;k++){
              lCol     = lId[k];
              if(lCol > -1){
/*... parte inferior*/
                if( lCol < lNeq)
                  for(iak = iPoint; iak < iaKneq;iak++){
                    jak = ja[iak];
                    if( lCol == jak ) 
                      al[iak] = lA[k];
                  }
/*... parte superior*/
                else if( lCol > lNeq && lCol <= neqS){
                  jPoint   = ia[lCol];
                  iaJneq   = ia[lCol+1];
                  for(iak = jPoint; iak < iaJneq;iak++){
                    jak = ja[iak];
                    if( lNeq == jak )
                      au[iak] = lA[k];
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

/*... simetrico*/
          else{
            for(k=0;k<nFace;k++){
              lCol     = lId[k];
              if(lCol > -1){
                for(iak = iPoint; iak < iaKneq;iak++){
                  jak = ja[iak];
                  if( lCol == jak ) 
                    al[iak] = lA[k];  
                }
/*...................................................................*/
              }      
/*...................................................................*/
            }
/*...................................................................*/
          }
/*...................................................................*/

/*... parte retangular em COO*/
          if(fCoo){
            n0     = nEqNov+1;
            iar    = &ia[n0];
            jar    = &ja[nAd];
/*...*/
            if(unsym)
              ar     = &a[2*nAd];
            else
              ar     = &a[nAd];
/*...................................................................*/
            iPoint = iar[lNeq];
            iaKneq = iar[lNeq+1];
            for(k=0;k<nFace;k++){
              lCol     = lId[k];
              if(lCol > -1 && lCol > neqS)
                for(iak = iPoint; iak < iaKneq;iak++)
                  if(jar[iak] == lCol) 
                    ar[iak] = lA[k]; 
            }
/*...................................................................*/
          }   
/*...................................................................*/

/*... parte retangular em CSR*/
          else{
            n0     = nEqNov+1;
            iar    = &ia[n0];
            jar    = &ja[nAd];
/*...*/
            if(unsym)
              ar     = &a[2*nAd];
            else
              ar     = &a[nAd];
/*...................................................................*/
            iPoint = iar[lNeq];
            iaKneq = iar[lNeq+1];
            for(k=0;k<nFace;k++){
              lCol     = lId[k];
              if(lCol > -1 && lCol > neqS){
                for(iak = iPoint; iak < iaKneq;iak++){
                  jak = jar[iak];
                  if( lCol == jak ) 
                    ar[iak] = lA[k];
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

/*... sequencial */                                                        
        else{
/*... nao simetrico*/
          if(unsym){
            au = &a[nAd];
            for(k=0;k<nFace;k++){
              lCol     = lId[k];
              if(lCol > -1){
/*... parte inferior*/
                if( lCol < lNeq)
                  for(iak = iPoint; iak < iaKneq;iak++){
                    jak = ja[iak];
                    if( lCol == jak ) 
                      al[iak] = lA[k];
                  }
/*... parte superior*/
                else if( lCol > lNeq){
                  jPoint   = ia[lCol];
                  iaJneq   = ia[lCol+1];
                  for(iak = jPoint; iak < iaJneq;iak++){
                    jak = ja[iak];
                    if( lNeq == jak ) {
                      au[iak] = lA[k];
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

/*... simetrico*/
          else{
            for(k=0;k<nFace;k++){
              lCol     = lId[k];
              if(lCol > -1){
                for(iak = iPoint; iak < iaKneq;iak++){
                  jak = ja[iak];
                  if( lCol == jak ) 
                    al[iak] = lA[k];  
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

//#ifdef _OPENMP
/*********************************************************************
* Data de criacao    : 27/07/2016                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* PARTITIOBCSRBYNOZEROS: divisa do trabalho por threads para matriz *
* no formato CSR                                                    *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* ia       -> ponteiro CSR                                          *
* ja       -> ponteiro CSR                                          *
* neq      -> numero de equacoes                                    *
* thBegin  -> nao definido                                          *
* thEnd    -> nao definido                                          *
* thSize   -> nao definido                                          *
* thHeigth -> nao definido                                          *
* type     -> CSR,CSRD,CSRC                                         *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* thBegin  -> primeira linha do sistema do thread i                 *
* thEnd    -> ultima linha do sistema do thread i                   *
* thSize   -> numero de termo nao nulos no thread i                 *
* thHeight -> altura efetiva do thread                              *
*-------------------------------------------------------------------*
* OBS: retorna o numero de elmentos nao nulos                       *
*-------------------------------------------------------------------*
*********************************************************************/
void partitionCsrByNonzeros(INT *restrict ia      ,INT *restrict ja
                           ,INT const neq
                           ,int nThreads          ,INT *restrict thBegin
                           ,INT *restrict thEnd   ,INT *restrict thSize
                           ,INT *restrict thHeight,short type) {

  INT nad, meanVariables, line;
  INT tam;
  int i, nTh;

  for (i = 0; i<nThreads; i++) {
    thBegin[i] = 0;
    thEnd[i] = 0;
    thSize[i] = 0;
    thHeight[i] = 0;
  }
/*se o numero de threads for maior que o numero de equacoes*/
  if (neq < nThreads)
    nThreads = neq;
  
  nad = ia[neq];
  switch (type) {
/*... CSR padrao*/
    case CSR:
      meanVariables = nad / nThreads;
      line = 1;
      thBegin[0] = 0;
      for(i = 0; i<nThreads - 1; i++) {
        thSize[i] = 0;
        tam = 0;
        for (;;) {
          tam = ia[line] - ia[line - 1];
          thSize[i] += tam;
          thEnd[i] = line - 1;
          thBegin[i + 1] = line;
          line++;
          if ((thSize[i] + tam) > meanVariables) break;
          if (line > neq) {
            ERRO_GERAL(__FILE__, __func__,__LINE__
                      ,"numero de linhas excedido");
          }
        }
      }
      nTh = nThreads - 1;
      thSize[nTh] = ia[neq] - ia[thBegin[nTh]];
      thEnd[nTh] = neq - 1;
    break;
/*...................................................................*/

/*... CSRD - csr sem a diagonal principal*/
    case CSRD:
      meanVariables = (nad + neq) / nThreads;
      line = 1;
      thBegin[0] = 0;
      for (i = 0; i<nThreads - 1; i++) {
        thSize[i] = 0;
        tam = 0;
        for (;;) {
          tam = ia[line] - ia[line - 1] + 1;
          thSize[i] += tam;
          thEnd[i] = line;
          thBegin[i + 1] = line;
          line++;
          if ((thSize[i] + tam) > meanVariables) break;
          if (line > neq) {
            ERRO_GERAL(__FILE__, __func__, __LINE__
                       , "numero de linhas excedido");
          }
        }
      }
      nTh = nThreads - 1;
      thEnd[nTh] = neq;
      thSize[nTh] = ia[neq] - ia[thBegin[nTh]] +
                   (thEnd[nTh] - thBegin[nTh] + 1);
    break;
/*...................................................................*/

/*... CSRC*/
    case CSRC:
      meanVariables = (2 * nad + neq) / nThreads;
      line = 1;
      thBegin[0] = 0;
      for (i = 0; i<nThreads - 1; i++) {
        thSize[i] = 0;
        tam = 0;
        for (;;) {
          tam = 2 * (ia[line] - ia[line - 1]) + 1;
          thSize[i] += tam;
          thEnd[i] = line;
          thBegin[i + 1] = line;
          line++;
          if ((thSize[i] + tam) > meanVariables) break;
          if (line > neq) {
            ERRO_GERAL(__FILE__, __func__, __LINE__
                      , "numero de linhas excedido");
          }
        }
      }
      nTh = nThreads - 1;
      thEnd[nTh] = neq;
      thSize[nTh] = 2 * (ia[neq] - ia[thBegin[nTh]]) + 1;

/*... calcula o tamanho efetivo do buffer*/
      computeEffectiveWork(ia     ,ja
                          ,neq
                          ,thBegin,thEnd
                          ,thSize ,thHeight);

    break;
/*...................................................................*/

/*...*/
    default:
      ERRO_OP(__FILE__, __func__, type);
    break;
  }
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 27/07/2016                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* COMPUTEEFFECTIVEWORK: calcula o trabalho efetivo por thread       *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* ia      -> ponteiro CSR                                           *
* ja      -> ponteiro CSR                                           *
* neq     -> numero de equacoes                                     *
* thBegin -> nao definido                                           *
* thEnd   -> nao definido                                           *
* thSize  -> nao definido                                           *
* thHeight-> nao definido                                           *
* type    -> CSR,CSRD,CSRC                                          *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* thBegin  -> primeira linha do sistema do thread i                 *
* thEnd    -> ultima linha do sistema do thread i                   *
* thSize   -> numero de termo nao nulos no thread i                 *
* thHeight -> altura efetiva do thread                              *
*-------------------------------------------------------------------*
* OBS: retorna o numero de elmentos nao nulos                       *
*-------------------------------------------------------------------*
*********************************************************************/
void computeEffectiveWork(INT *restrict ia, INT *restrict ja
                    ,INT const nEq
                    ,INT *restrict thBegin,INT *restrict thEnd
                    ,INT *restrict thSize ,INT *restrict thHeight) {
  int i, id = 0, h = 0;

#pragma omp parallel private(id,h) 
  {
    id = omp_get_thread_num();
    h = thBegin[id];
    for (i = thBegin[id]; i<thEnd[id]; i++)
      h = min(h, ja[ia[i]]);
    thHeight[id] = h;
  }

}
/*********************************************************************/
//#endif