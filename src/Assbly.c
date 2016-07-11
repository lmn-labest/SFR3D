#include<Assbly.h>
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * ASSBLY: Montagem do sistema global                                * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * ia      -> ponteiro para as linhas da matriz esparsa              * 
 * ja      -> ponteiro para as colunas da matriz esparsa             * 
 * a       -> matriz de coeficientes esparsa                         * 
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- fora da diagonal principal        )       *
 *            ( ELLPACK  - fora da diagonal principal + zeros)       *
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 *            ( ELLPACK  - diagonal principal                )       *
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
 * au,a,al   -> coeficiente da celula i    (matriz = true)           *
 * b         -> vetor de forca da celula i (forces = true)           *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void assbly(INT    *restrict  ia,INT *restrict ja 
           ,double *restrict a  ,double *restrict ad
           ,double *restrict b
           ,INT *restrict lId
           ,double *restrict lA ,double *restrict lB
           ,INT const nEq       ,INT const nEqNov
           ,INT const nAd       ,INT const nAdR  
           ,short const nFace   ,short const ndf 
           ,short const storage ,bool  const forces  
           ,bool const matrix   ,bool  const  unsym)
{
  

  switch (storage){
/*... estrutura CSR/CSRD/CSRDCOO/CSRC*/
    case CSR:
    case CSRD:
    case CSRC:
    case CSRDCOO:
    case CSRCCOO:
      csr(ia     ,ja
         ,a      ,ad
         ,b      ,lId
         ,lA     ,lB
         ,nEq    ,nEqNov   
         ,nAd    ,nAdR  
         ,nFace  ,ndf
         ,storage,forces
         ,matrix ,unsym);
    break;
/*...................................................................*/

/*... ellPack*/
    case ELLPACK:
      ellPack(ia     
             ,ad     ,a    
             ,b
             ,lId
             ,lA     ,lB
             ,nFace  ,ndf
             ,storage,forces
             ,matrix);
    break;
/*...................................................................*/

/*...*/
    default:
      ERRO_OP(__FILE__,__func__,storage);
    break;
/*...................................................................*/
  }
/*...................................................................*/
} 
/*********************************************************************/


/********************************************************************* 
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * ASSBLY: Montagem do sistema global para sistems de velocidades    *
 * do metodo simple                                                  *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * ia      -> ponteiro para as linhas da matriz esparsa              * 
 * ja      -> ponteiro para as colunas da matriz esparsa             * 
 * a       -> matriz de coeficientes esparsa                         * 
 *            ( CSR - nao utiliza                            )       *
 *            ( CSRD/CSRC- fora da diagonal principal        )       *
 *            ( ELLPACK  - fora da diagonal principal + zeros)       *
 * ad      -> matrix de coeficientes esparsa                         *
 *            ( CSR - matriz completa                        )       *
 *            ( CSRD/CSRC- diagonal principal                )       *
 *            ( ELLPACK  - diagonal principal                )       *
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
 * au,a,al   -> coeficiente da celula i    (matriz = true)           *
 * b         -> vetor de forca da celula i (forces = true)           *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------*
 * Os coeficientes das matrizes são indentencios para as tres        *
 * componentes de velocidades (u1,u2,u3), porém o vetor de forcas e  *
 * diferente                                                         * 
 *********************************************************************/
void assblySimple(INT    *restrict  ia,INT *restrict ja 
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
  

  switch (storage){
/*... estrutura CSR/CSRD/CSRDCOO/CSRC*/
    case CSR:
    case CSRD:
    case CSRC:
    case CSRDCOO:
    case CSRCCOO:
      csrSimple(ia     ,ja
         ,a      ,ad
         ,b      ,lId
         ,lA     ,lB
         ,nEq    ,nEqNov   
         ,nAd    ,nAdR  
         ,nFace  ,ndf
         ,storage,forces
         ,matrix ,unsym);
    break;
/*...................................................................*/

/*... ellPack*/
    case ELLPACK:
      ellPack(ia     
             ,ad     ,a    
             ,b
             ,lId
             ,lA     ,lB
             ,nFace  ,ndf
             ,storage,forces
             ,matrix);
    break;
/*...................................................................*/

/*...*/
    default:
      ERRO_OP(__FILE__,__func__,storage);
    break;
/*...................................................................*/
  }
/*...................................................................*/
} 
    
