#include<WriteCsv.h>
/********************************************************************* 
 * WRITECSVCELL: escreve os resultados por celula no formato csv     *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * u         -> solucao conhecida                                    * 
 * gradU     -> gradiente rescontruido da solucao conhecida          * 
 * cc        -> centroide da celula                                  * 
 * nCell     -> numero de celulas                                    * 
 * ndm       -> numero de dimensoes                                  * 
 * ndf       -> grauss de liberdade                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lSquare   -> matriz para a reconstrucao least Square              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void writeCsvCell(DOUBLE *u      ,DOUBLE *gradU
                 ,DOUBLE *cc 
                 ,INT const nCell,short const ndf
                 ,short const ndm,FILE *file){
  INT nEl;
  short j,k;


  fprintf(file,"#nCell");
  for(j=0;j<ndm;j++)
    fprintf(file,",xcc(x%d)",j+1);
  
  for(j=0;j<ndf;j++)
    fprintf(file,",u%d",j+1);
  
  for(j=0;j<ndf;j++)
    for(k=0;k<ndm;k++)
      fprintf(file,",gradU(u%d)(x%d)",j+1,k+1);
  
  fprintf(file,"\n");
  for(nEl=0;nEl<nCell;nEl++){
    fprintf(file,"%9d",nEl+1);

/*... centroide*/    
    for(j=0;j<ndm;j++)
      fprintf(file,",%e",MAT2D(nEl,j,cc,ndm));

/*... valores de u1*/    
    for(j=0;j<ndf;j++)
      fprintf(file,",%e",MAT2D(nEl,j,u,ndf));

/*... valores de gradU1*/    
    for(j=0;j<ndf;j++)
      for(k=0;k<ndm;k++)
      fprintf(file,",%e",MAT3D(nEl,j,k,gradU,ndf,ndm));

    fprintf(file,"\n");
    
  }

}
/*********************************************************************/
