#include<WriteVtk.h>
/********************************************************************** 
 * WGEOVTK : escreve a malha com os resultados e condicao             *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * m       -> arranjo da menoria principal                            *  
 * x       -> coordenadas                                             * 
 * el      -> conectividade                                           * 
 * mat     -> materias                                                *  
 * nen     -> conectividades por elemento                             *  
 * mat     -> material por elemento                                   *
 * typeGeom-> tipo geometrico do elemento                             *
 * typeCal -> tipo calculo do elemento                                *
 * nel     -> numeracao do elemento                                   *
 * faceRt1 -> tipo de condicao de contorno T1                         * 
 * faceSt1 -> valor de condicao de contorno T1                        * 
 * nnode   -> numero de nos                                           *  
 * numel   -> numero de elementos                                     *
 * ndm     -> numero de dimensao                                      *
 * maxNo   -> numero maximo de nos por elemento                       *
 * numat   -> numero maximo de nos por elemento                       *
 * ndfT    -> graus de liberdade das equacoes de tranporte T1         *
 * nameOut -> nome de arquivo de saida                                *
 * iws     -> vtk binario                                             *
 * f       -> arquivlo                                                *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 **********************************************************************/
void wGeoVtk(Memoria *m     ,double *x      
            ,INT *el        ,short *mat    
            ,short *nen     ,short *typeGeom
            ,double *prop   ,short *typeCal
            ,short *faceRt1 ,double *faceSt1
            ,INT nnode      ,INT numel    
            ,short ndm      ,short maxNo 
            ,short numat    ,short *ndfT   
            ,char *nameOut  ,bool iws
            ,FILE *f)
{
  int    *lel=NULL;
  double *aux=NULL;
  INT i;
  short j,k;
  char head[]={"GEOM_VOLUME_FINITO"};
  double ddum;
  int    idum;

  if(iws)
    f = openFile(nameOut,"wb");
  else
    f = openFile(nameOut,"w");

/* ...*/
  headVtk(head,iws,f);
/* ..................................................................*/

/*... coordenadas*/
  writeVtkCoor(x,nnode,ndm,iws,f);
/*...................................................................*/
  
/*... conectividades*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++){
    for(j=0;j<maxNo;j++){
     k = i*maxNo+j;
     lel[k]=el[k]-1;
    }
  }  
  writeVtkCell(lel,nen,typeGeom,numel,maxNo,iws,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) numel);
/*...................................................................*/

/*... material*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++)
    lel[i]=(int) mat[i];
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... numero do elemento*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++)
    lel[i]= i+1;
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... tipo celula para o calculo*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  
  for(i=0;i<numel;i++){
    idum = mat[i]-1;
    lel[i] = typeCal[idum];
  }
   
  writeVtkProp(lel,&ddum,numel,1,"elTyCal",iws,INTEGER,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... propriedades dos elementos*/
  HccaAlloc(double,m,aux,numel*numat,"el",_AD_);
  if( aux == NULL){
    fprintf(stderr,"Erro na alocação de aux.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++){
    for(j=0;j<MAXPROP;j++){
      idum = mat[i]-1;
      MAT2D(i,j,aux,MAXPROP) = MAT2D(idum,j,prop,MAXPROP);
    }
  } 
  writeVtkProp(&idum,aux,numel,MAXPROP,"elProp",iws,DOUBLEV,f);
  HccaDealloc(m,aux,"el",_AD_);
/*...................................................................*/

/*...*/
  if(ndfT[0] > 0 ){
/*... faceRt1*/
    HccaAlloc(int,m,lel,numel*(maxNo+1),"el",_AD_);
    if( lel == NULL){
      fprintf(stderr,"Erro na alocação de lel.\n"
                     "Nome do arquivo: %s.\n"
                    ,__FILE__);
      exit(EXIT_FAILURE);
    }
    for(i=0;i<numel*(maxNo+1);i++)
      lel[i]=(int) faceRt1[i];
   
    writeVtkProp(lel,&ddum,numel,maxNo+1,"faceRt1",iws,INTEGER,f);
    HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... faceSt1*/
    writeVtkProp(&idum,faceSt1,numel,maxNo+1,"faceSt1",iws
                    ,DOUBLEV,f);
/*...................................................................*/
  }
/*...................................................................*/

/*.... campo por no*/
  fprintf(f,"POINT_DATA %ld\n",(long) nnode);
/*...................................................................*/

/*... numero do no*/
  HccaAlloc(int,m,lel,nnode,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<nnode;i++)
    lel[i]=i+1;
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  fclose(f);
}
/*********************************************************************/ 

/********************************************************************** 
 * WRESVTK : escreve a malha com os resultados e condicao             *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * m       -> arranjo da menoria principal                            *  
 * x       -> coordenadas                                             * 
 * el      -> conectividade                                           * 
 * mat     -> materias                                                *  
 * nen     -> conectividades por elemento                             *  
 * mat     -> material por elemento                                   *
 * typeGeom-> tipo geometrico do elemento                             *
 * typeCal -> tipo calculo do elemento                                *
 * nel     -> numeracao do elemento                                   *
 * faceRt1 -> tipo de condicao de contorno T1                         * 
 * faceSt1 -> valor de condicao de contorno T1                        * 
 * nnode   -> numero de nos                                           *  
 * numel   -> numero de elementos                                     *
 * ndm     -> numero de dimensao                                      *
 * maxNo   -> numero maximo de nos por elemento                       *
 * numat   -> numero maximo de nos por elemento                       *
 * ndfT    -> graus de liberdade das equacoes de tranporte T1         *
 * nameOut -> nome de arquivo de saida                                *
 * iws     -> vtk binario                                             *
 * f       -> arquivlo                                                *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 **********************************************************************/
void wResVtk(Memoria *m     ,double *x      
            ,INT *el        ,short *mat    
            ,short *nen     ,short *typeGeom
            ,DOUBLE *elU    ,DOUBLE *nU
            ,INT nnode      ,INT numel    
            ,short ndm      ,short maxNo 
            ,short numat    ,short *ndfD   
            ,char *nameOut  ,bool iws
            ,FILE *f)
{
  int    *lel=NULL;
  INT i;
  short j,k;
  char head[]={"GEOM_VOLUME_FINITO"};
  double ddum;
  int    idum;

  if(iws)
    f = openFile(nameOut,"wb");
  else
    f = openFile(nameOut,"w");

/* ...*/
  headVtk(head,iws,f);
/* ..................................................................*/

/*... coordenadas*/
  writeVtkCoor(x,nnode,ndm,iws,f);
/*...................................................................*/
  
/*... conectividades*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++){
    for(j=0;j<maxNo;j++){
     k = i*maxNo+j;
     lel[k]=el[k]-1;
    }
  }  
  writeVtkCell(lel,nen,typeGeom,numel,maxNo,iws,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) numel);
/*...................................................................*/

/*... material*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++)
    lel[i]=(int) mat[i];
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... numero do elemento*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++)
    lel[i]= i+1;
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... escrever resuldos por celula*/  
  if(ndfD[0] > 0){
    writeVtkProp(&idum,elU,numel,ndfD[0],"elD1",iws,DOUBLEV,f);
  }
/*...................................................................*/

/*.... campo por no*/
  fprintf(f,"POINT_DATA %ld\n",(long) nnode);
/*...................................................................*/

/*... numero do no*/
  HccaAlloc(int,m,lel,nnode,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<nnode;i++)
    lel[i]=i+1;
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... escrever resuldos por nos*/  
  if(ndfD[0] > 0){
    writeVtkProp(&idum,nU,nnode,ndfD[0],"nD1",iws,DOUBLEV,f);
  }
/*...................................................................*/

  fclose(f);
}
/*********************************************************************/
