#include<Assbly.h>
/********************************************************************* 
 * ASSBLY: Montagem do sistema global                                * 
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
           ,double *restrict au ,double *restrict ad
           ,double *restrict al ,double *restrict b
           ,INT *restrict lId
           ,double *restrict lA ,double *restrict lB
           ,short const nFace   ,short const ndf 
           ,short const storage ,bool  const forces  
           ,bool const matrix   ,bool  const  unsym)
{
  switch (storage){
/*... estrutura CSR/CSRD/CSRC*/
    case CSR:
    case CSRD:
    case CSRC:
      csr(ia     ,ja
         ,au     ,ad
         ,al     ,b
         ,lId
         ,lA     ,lB
         ,nFace  ,ndf
         ,storage,forces
         ,matrix ,unsym);
    break;
/*...................................................................*/

/*... ellPack*/
    case ELLPACK:
      ellPack(ia     
             ,ad     ,al   
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
    
