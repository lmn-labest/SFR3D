#include<EllPack.h>
/********************************************************************* 
 * ELLPACKJA:                                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ifEllPack  -> mao definido                                        * 
 * ja         -> indefinido                                          * 
 * id         -> numeracao das equacoes por elemento                 * 
 * num        -> renumeracao dos elementos                           * 
 * adj        -> adjacencia dos elementos                            * 
 * nViz       -> numero de vizinhos por elemento                     * 
 * numel      -> numero de elementos                                 * 
 * neq        -> numero de equacoes                                  * 
 * maxLineNzr -> numero maximo de nao zeros por linha                * 
 * ndf        -> numero de graus de liberdade                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ifEllPack -> informacao do ellPack                                * 
 * ja        -> ponteiro do ELLPACK                                  * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
INT ellPackJa(INT *RESTRICT ifEllPack,INT *RESTRICT ja
             ,INT *RESTRICT id       ,INT *RESTRICT num
             ,INT *RESTRICT adj      ,short *RESTRICT nViz
             ,INT const numel        ,INT const nEq    
             ,short const maxLineNzr ,const short ndf)
{

  INT i,nel1,neq1,viz1,col,nad;
  short jNdf,kNdf,j;

  ifEllPack[0] = maxLineNzr;
  ifEllPack[1] = nEq*maxLineNzr;

  nad = 0;
/*... loop na celulas*/
  for(i=0;i<numel;i++){
/*... numeracao da celula*/
    nel1 = num[i] - 1;
    for(jNdf=0;jNdf<ndf;jNdf++){
      neq1 = MAT2D(nel1,jNdf,id,ndf) - 1;
/*... celulas com equacoes associada*/
      if( neq1 != -2){
/*...loop nas faces da celula*/
        for(j=0;j<nViz[nel1];j++){
          viz1 = MAT2D(nel1,j,adj,maxLineNzr) - 1;
          if( viz1 != -2){
/*...equacao associada a face/celula vizinha*/
            for(kNdf=0;kNdf<ndf;kNdf++){
              col = MAT2D(viz1,kNdf,id,ndf) - 1;
/*... termos nao nulos*/
              if( col != -2){
                MAT2D(neq1,j,ja,maxLineNzr) = col; 
                nad++; 
              }
/*...................................................................*/    
            } 
/*...................................................................*/    
          }
/*... termos nulos do contorno*/
          else
            MAT2D(neq1,j,ja,maxLineNzr) = neq1;   
/*...................................................................*/    
        }
/*...................................................................*/    
     }
/*...................................................................*/    
    }
/*...................................................................*/    
  }
/*...................................................................*/

  return nad;

}
/*********************************************************************/ 

/********************************************************************* 
 * ELLPACK : Montagem do sistema global do ELLPACK                   * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * ifEllPack -> informacao do ellPack                                * 
 * ad      -> matrix de coeficientes esparsa (diagonal principal)    *
 * a       -> matriz de coeficientes esparsa (off diagonal)          *
 * b       -> vetor de forcas                                        *
 * lA      -> coficientes da celula                                  *
 * lB      -> vetor de forca da celula                               *
 * lId     -> numeracao das equacoes dessa celula                    *
 * nFace   -> numero de faces da celula                              *
 * ndf     -> graus de liberdade                                     *
 * storage -> tecnica de armazenamento da matriz esparsa             * 
 * forces  -> mantagem no vetor de forcas                            * 
 * matrix  -> mantagem da matriz de coeficientes                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ad,a      -> coeficiente da linha i     (matriz = true)           *
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
void ellPack(INT *RESTRICT  ifEllPack
        ,DOUBLE *RESTRICT ad            ,DOUBLE *RESTRICT a
        ,DOUBLE *RESTRICT b
        ,INT *RESTRICT lId                       
        ,DOUBLE *RESTRICT lA            ,DOUBLE *RESTRICT lB 
        ,short const nFace              ,short const ndf  
        ,short const storage            ,bool  const forces
        ,bool const matrix              )
{
  INT lNeq;
  unsigned short k,maxLineNzr;
  char str[300];

  maxLineNzr = ifEllPack[0];
     
/*... */
  if(ndf == 1){
    lNeq    = lId[nFace];
/*... vetor de forcas*/
    if(forces) 
      b[lNeq] = lB[0];
/*...................................................................*/
    if(matrix){
      ad[lNeq] = lA[nFace];
      for(k=0;k<nFace;k++)
        MAT2D(lNeq,k,a,maxLineNzr) = lA[k];
    } 
  }
/*...................................................................*/

/*... */
  else{
    strcpy(str,"Nao implentados ndf > 1");
    ERRO_GERAL(__FILE__,__func__,__LINE__,str);
  }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * BANDELLPACK: banda da matriz no formato ELLPACK                   * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ifEllPack -> informacao do ellPack                                * 
 * ja  - ponteiro ellPack                                            * 
 * neq - numero de equacoes                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS: retorna a banda da matrix                                    * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
INT bandEllPack(INT *ifEllPack,INT *ja,INT  neq,short type){

  INT i,j,aux,col;
  long bandL=0;
  INT maxLineNzr = ifEllPack[0];

  switch(type){
/*... banda maxima da matriz*/
    case 1:
      for(i=0;i<neq;i++)
        for(j=0;j<maxLineNzr;j++){
          col = MAT2D(i,j,ja,maxLineNzr);
          if( col != i) 
            bandL = max(bandL,abs(i-col));
        }
    break;
/*...................................................................*/ 

/*... banda media da matriz*/
    case 2: 
      for(i=0;i<neq;i++){
        aux = 0;
        for(j=0;j<maxLineNzr;j++){
          col = MAT2D(i,j,ja,maxLineNzr);
          if( col != i) 
            aux = max(aux,abs(i-col));
        }
        bandL += aux;
      }
      bandL = bandL/neq;
    break;
/*...................................................................*/ 

/*... banda minima da matriz*/
    case 3: 
      bandL = neq;
      for(i=0;i<neq;i++)
        for(j=0;j<maxLineNzr;j++){
          col = MAT2D(i,j,ja,maxLineNzr);
          if( col != i) 
            bandL = min(bandL,abs(i-col));
       }
    break;
/*...................................................................*/ 
  }
  return bandL;

}
/*********************************************************************/ 
