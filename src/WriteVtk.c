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
 * maxno   -> numero maximo de nos por elemento                       *
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
            ,short ndm      ,short maxno 
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
  Myalloc(int,m,lel,numel*maxno,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++){
    for(j=0;j<maxno;j++){
     k = i*maxno+j;
     lel[k]=el[k]-1;
    }
  }  
  writeVtkCell(lel,nen,typeGeom,numel,maxno,iws,f);
  Mydealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) numel);
/*...................................................................*/

/*... material*/
  Myalloc(int,m,lel,numel*maxno,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++)
    lel[i]=(int) mat[i];
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER,f);
  Mydealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... numero do elemento*/
  Myalloc(int,m,lel,numel*maxno,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++)
    lel[i]= i+1;
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER,f);
  Mydealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... tipo celula para o calculo*/
  Myalloc(int,m,lel,numel*maxno,"el",_AD_);
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
  Mydealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... propriedades dos elementos*/
  Myalloc(double,m,aux,numel*numat,"el",_AD_);
  if( aux == NULL){
    fprintf(stderr,"Erro na alocação de aux.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++){
    for(j=0;j<numat;j++){
      idum = mat[i]-1;
      MAT2D(i,j,aux,numat) = MAT2D(idum,j,prop,MAXPROP);
    }
  } 
  writeVtkProp(&idum,aux,numel,numat,"elProp",iws,DOUBLE,f);
  Mydealloc(m,aux,"el",_AD_);
/*...................................................................*/

/*...*/
  if(ndfT[0] > 0 ){
/*... faceRt1*/
    Myalloc(int,m,lel,numel*(maxno+1),"el",_AD_);
    if( lel == NULL){
      fprintf(stderr,"Erro na alocação de lel.\n"
                     "Nome do arquivo: %s.\n"
                    ,__FILE__);
      exit(EXIT_FAILURE);
    }
    for(i=0;i<numel*(maxno+1);i++)
      lel[i]=(int) faceRt1[i];
   
    writeVtkProp(lel,&ddum,numel,maxno+1,"faceRt1",iws,INTEGER,f);
    Mydealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... faceSt1*/
    writeVtkProp(&idum,faceSt1,numel,maxno+1,"faceSt1",iws
                    ,DOUBLE,f);
/*...................................................................*/
  }
/*...................................................................*/

/*.... campo por no*/
  fprintf(f,"POINT_DATA %ld\n",(long) nnode);
/*...................................................................*/

/*... numero do elemento*/
  Myalloc(int,m,lel,nnode,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<nnode;i++)
    lel[i]=i+1;
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER,f);
  Mydealloc(m,lel,"el",_AD_);
/*...................................................................*/
  fclose(f);
}
