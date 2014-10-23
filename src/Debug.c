#include<Debug.h>
/*********************************************************************/
void testeGeom(double *cc
              ,double *ksi   ,double *mksi 
              ,double *eta   ,double *meta 
              ,double *normal,double *volume
              ,double *xm    ,double *xmcc   
              ,double *mkm   ,double *dcca                 
              ,INT numel     ,short ndm
              ,short maxViz)
{
#ifdef _DEBUG_GEOM_
  INT i;
  short j,k;
  
  fprintf(stderr,"Centroide.\n");
  for(i=0;i<numel;i++){
    fprintf(stderr,"nel = %ld ",i+1);
    for(j=0;j<ndm;j++)
      fprintf(stderr,"%8.4lf ",MAT2D(i,j,cc,ndm));
    fprintf(stderr,"\n");
  }

  fprintf(stderr,"\nVetor entre os centroides.\n");
  for(i=0;i<numel;i++){
    fprintf(stderr,"nel = %ld ",i+1);
    for(j=0;j<maxViz;j++)
      for(k=0;k<ndm;k++)
        fprintf(stderr,"%8.4lf ",MAT3D(i,j,k,ksi,maxViz,ndm));
    fprintf(stderr,"\n");
  }

  fprintf(stderr,"\nDistancia entre os centroides.\n");
  for(i=0;i<numel;i++){
    fprintf(stderr,"nel = %ld ",i+1);
    for(j=0;j<maxViz;j++)
      fprintf(stderr,"%8.4lf ",MAT2D(i,j,mksi,maxViz));
    fprintf(stderr,"\n");
  }
  
  fprintf(stderr,"\nVetor paralelo a face.\n");
  for(i=0;i<numel;i++){
    fprintf(stderr,"nel = %ld ",i+1);
    for(j=0;j<maxViz;j++)
      for(k=0;k<ndm;k++)
        fprintf(stderr,"%8.4lf ",MAT3D(i,j,k,eta,maxViz,ndm));
    fprintf(stderr,"\n");
  }
  
  fprintf(stderr,"\narea da face.\n");
  for(i=0;i<numel;i++){
    fprintf(stderr,"nel = %ld ",i+1);
    for(j=0;j<maxViz;j++)
      fprintf(stderr,"%8.4lf ",MAT2D(i,j,meta,maxViz));
    fprintf(stderr,"\n");
  }
  
  
  fprintf(stderr,"\nVetor normal a face.\n");
  for(i=0;i<numel;i++){
    fprintf(stderr,"nel = %ld ",i+1);
    for(j=0;j<maxViz;j++)
      for(k=0;k<ndm;k++)
        fprintf(stderr,"%8.4lf ",MAT3D(i,j,k,normal,maxViz,ndm));
    fprintf(stderr,"\n");
  }
  
  fprintf(stderr,"\nVolume da celula.\n");
  for(i=0;i<numel;i++){
    fprintf(stderr,"nel = %ld ",i+1);
    fprintf(stderr,"%8.4lf ",volume[i]);
    fprintf(stderr,"\n");
  }
  
  fprintf(stderr,"\nPonto medio\n");
  for(i=0;i<numel;i++){
    fprintf(stderr,"nel = %ld ",i+1);
    for(j=0;j<maxViz;j++)
      for(k=0;k<ndm;k++)
        fprintf(stderr,"%8.4lf ",MAT3D(i,j,k,xm,maxViz,ndm));
    fprintf(stderr,"\n");
  }
  
  fprintf(stderr,"\nVetor ponto medio centroide\n");
  for(i=0;i<numel;i++){
    fprintf(stderr,"nel = %ld ",i+1);
    for(j=0;j<maxViz;j++)
      for(k=0;k<ndm;k++)
        fprintf(stderr,"%8.4lf ",MAT3D(i,j,k,xmcc,maxViz,ndm));
    fprintf(stderr,"\n");
  }
  
  fprintf(stderr,"\nMenor distancia do centroide a Face.\n");
  for(i=0;i<numel;i++){
    fprintf(stderr,"nel = %ld ",i+1);
    for(j=0;j<maxViz;j++)
      fprintf(stderr,"%8.4lf ",MAT2D(i,j,dcca,maxViz));
    fprintf(stderr,"\n");
  }
  
  fprintf(stderr,"\nDistancia entre o ponto medio da face e ponto.\n");
  for(i=0;i<numel;i++){
    fprintf(stderr,"nel = %ld ",i+1);
    for(j=0;j<maxViz;j++)
      fprintf(stderr,"%8.4lf ",MAT2D(i,j,mkm,maxViz));
    fprintf(stderr,"\n");
  }
#endif
}
/*********************************************************************/

