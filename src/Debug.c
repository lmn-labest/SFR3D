#include<Debug.h>
/*********************************************************************/
void testeGeom(double *cc
              ,double *ksi   ,double *mksi 
              ,double *eta   ,double *meta 
              ,double *normal,double *volume
              ,double *xm    ,double *xmcc   
              ,double *vSkew ,double *mvSkew   
              ,double *dcca                 
              ,INT numel     ,short ndm
              ,short maxViz)
{
#ifdef _DEBUG_
  INT i;
  short j,k;
  
  printf("Centroide.\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    for(j=0;j<ndm;j++)
      printf("%8.4lf ",MAT2D(i,j,cc,ndm));
    printf("\n");
  }

  printf("\nVetor entre os centroides.\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    for(j=0;j<maxViz;j++)
      for(k=0;k<ndm;k++)
        printf("%8.4lf ",MAT3D(i,j,k,ksi,maxViz,ndm));
    printf("\n");
  }

  printf("\nDistancia entre os centroides.\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    for(j=0;j<maxViz;j++)
      printf("%8.4lf ",MAT2D(i,j,mksi,maxViz));
    printf("\n");
  }
  
  printf("\nVetor paralelo a face.\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    for(j=0;j<maxViz;j++)
      for(k=0;k<ndm;k++)
        printf("%8.4lf ",MAT3D(i,j,k,eta,maxViz,ndm));
    printf("\n");
  }
  
  printf("\narea da face.\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    for(j=0;j<maxViz;j++)
      printf("%8.4lf ",MAT2D(i,j,meta,maxViz));
    printf("\n");
  }
  
  
  printf("\nVetor normal a face.\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    for(j=0;j<maxViz;j++)
      for(k=0;k<ndm;k++)
        printf("%8.4lf ",MAT3D(i,j,k,normal,maxViz,ndm));
    printf("\n");
  }
  
  printf("\nVolume da celula.\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    printf("%8.4lf ",volume[i]);
    printf("\n");
  }
  
  printf("\nPonto medio\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    for(j=0;j<maxViz;j++)
      for(k=0;k<ndm;k++)
        printf("%8.4lf ",MAT3D(i,j,k,xm,maxViz,ndm));
    printf("\n");
  }
  
  
  printf("\nVetor que une o centroide ao ponto medio\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    for(j=0;j<maxViz;j++)
      for(k=0;k<ndm;k++)
        printf("%8.4lf ",MAT3D(i,j,k,xmcc,maxViz,ndm));
    printf("\n");
  }

  
  printf("\nMenor distancia do centroide a Face.\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    for(j=0;j<maxViz;j++)
      printf("%8.4lf ",MAT2D(i,j,dcca,maxViz));
    printf("\n");
  }
  
  printf("\nVetor entre o ponto medio da face e ponto de intersecao.\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    for(j=0;j<maxViz;j++)
      for(k=0;k<ndm;k++)
        printf("%8.4lf ",MAT3D(i,j,k,vSkew,maxViz,ndm));
    printf("\n");
  }
  
  printf("\nDistancia entre o ponto medio da face e ponto de intersecao.\n");
  for(i=0;i<numel;i++){
    printf("nel = %d ",i+1);
    for(j=0;j<maxViz;j++)
      printf("%8.4lf ",MAT2D(i,j,mvSkew,maxViz));
    printf("\n");
  }
#endif
}
/*********************************************************************/

void testeSist(INT *ia      ,INT *ja
              ,double *au   ,double *ad
              ,double *al   ,double *b
              ,INT const neq,bool const unsym)
{
#ifdef _DEBUG_
  INT i,j;
  
  printf("ia.\n");
  for(i=0;i<neq;i++){
    printf("i   = %9d ia = ",i+1);
    printf("%9d\n",ia[i]);
  }
  
  printf("ja.\n");
  for(i=0;i<neq;i++){
    printf("i   = %9d ja = ",i+1);
    for(j=ia[i];j<ia[i+1];j++)
      printf("%9d ",ja[j]);
    printf("\n");
  }
  
  printf("ad.\n");
  for(i=0;i<neq;i++){
    printf("i   = %9d ",i+1);
    printf("%lf\n",ad[i]);
  }

  printf("al.\n");
  for(i=0;i<neq;i++){
    printf("i   = %d al = ",i+1);
    for(j=ia[i];j<ia[i+1];j++)
      printf("%lf ",al[j]);
    printf("\n");
  }

  if(unsym) {
    printf("au.\n");
    for(i=0;i<neq;i++){
      printf("i   = %9d au= ",i+1);
      for(j=ia[i];j<ia[i+1];j++)
        printf("%lf ",al[j]);
      printf("\n");
    }
  }
  
  printf("b.\n");
  for(i=0;i<neq;i++){
    printf("i   = %9d b =",i+1);
    printf("%lf\n",b[i]);
  }

#endif
}

