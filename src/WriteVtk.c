#include<WriteVtk.h>
/********************************************************************** 
 * WMESHPARTVTK: escreve a malha divida em particoes                  *  
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
 * nel     -> numeracao do elemento                                   *
 * nnode   -> numero de nos                                           *  
 * numel   -> numero de elementos                                     *
 * ndm     -> numero de dimensao                                      *
 * maxNo   -> numero maximo de nos por elemento                       *
 * maxNo   -> numero maximo de vizinhos por elemento                  *
 * numat   -> numero de materias                                      *
 * nameOut -> nome de arquivo de saida                                *
 * iws     -> vtk binario                                             *
 * f       -> arquivlo                                                *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 **********************************************************************/
void wMeshPartVtk(Memoria *m     
                  ,double *x      ,INT *el            
                  ,short *nen     ,short *typeGeom
                  ,INT nnode      ,INT numel    
                  ,short ndm      
                  ,short maxNo    ,short maxViz  
                  ,char *nameOut  ,bool iws
                  ,FILE *f)
{
  int    *lel=NULL;
  INT i;
  short j;
  char head[]={"GEOM_VOLUME_FINITO"};
  double ddum;

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
      MAT2D(i,j,lel,maxNo) = MAT2D(i,j,el,maxNo)-1;
    }
  }  
  writeVtkCell(lel,nen,typeGeom,numel,maxNo,iws,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) numel);
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
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
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
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  fclose(f);
}
/*********************************************************************/ 

/********************************************************************** 
 * WPARTVTK: escreve a malha particionda                              *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * m       -> arranjo da menoria principal                            *  
 * x       -> coordenadas                                             * 
 * el      -> conectividade                                           * 
 * mat     -> materias                                                *  
 * nen     -> conectividades por elemento                             *  
 * typeGeom-> tipo geometrico do elemento                             *
 * nel     -> numeracao do elemento                                   *
 * np      -> no por particao                                         * 
 * ep      -> elementos por particao                                  * 
 * nnode   -> numero de nos                                           *  
 * numel   -> numero de elementos                                     *
 * ndm     -> numero de dimensao                                      *
 * maxNo   -> numero maximo de nos por elemento                       *
 * maxNo   -> numero maximo de vizinhos por elemento                  *
 * numat   -> numero de materias                                      *
 * nameOut -> nome de arquivo de saida                                *
 * iws     -> vtk binario                                             *
 * f       -> arquivlo                                                *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 **********************************************************************/
void wPartVtk(Memoria *m     
            ,double *x      ,INT *el            
            ,short *nen     ,short *typeGeom
            ,INT *np        ,INT *ep
            ,INT nnode      ,INT numel    
            ,short ndm      
            ,short maxNo    ,short maxViz  
            ,char *nameOut  ,bool iws
            ,FILE *f)
{
  int    *lel=NULL;
  INT i;
  short j;
  char head[]={"GEOM_VOLUME_FINITO"};
  double ddum;

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
      MAT2D(i,j,lel,maxNo) = MAT2D(i,j,el,maxNo)-1;
    }
  }  
  writeVtkCell(lel,nen,typeGeom,numel,maxNo,iws,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) numel);
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
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

  writeVtkProp(ep,&ddum,numel,1,"elPart"  ,iws,INTEGER,1,f);

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
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... numero do no*/
  writeVtkProp(np,&ddum,nnode,1,"noPart",iws,INTEGER,1,f);
/*...................................................................*/
  fclose(f);
}
/*********************************************************************/ 

/********************************************************************** 
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 22/08/2017                                    *
 *------------------------------------------------------------------- *
 * WGEOVTK : escreve a malha com os resultados e condicao             *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * m          -> arranjo da menoria principal                         *  
 * x          -> coordenadas                                          * 
 * el         -> conectividade                                        * 
 * mat        -> materias                                             *  
 * nen        -> conectividades por elemento                          *  
 * mat        -> material por elemento                                *
 * typeGeom   -> tipo geometrico do elemento                          *
 * typeCal    -> tipo calculo do elemento                             *
 * nel        -> numeracao do elemento                                *
 * faceRd1    -> condicao de contorno D1                              * 
 * faceLd1    -> tipo da condicao de contorno D1                      * 
 * faceRt1    -> condicao de contorno T1                              * 
 * faceLt1    -> tipo da condicao de contorno T1                      * 
 * faceRfluid -> condicao de contorno fluido                          * 
 * faceLfluid -> tipo da condicao de contorno fluido                  * 
 * nnode      -> numero de nos                                        *  
 * numel      -> numero de elementos                                  *
 * ndm        -> numero de dimensao                                   *
 * maxNo      -> numero maximo de nos por elemento                    *
 * maxNo      -> numero maximo de vizinhos por elemento               *
 * numat      -> numero de materias                                   *
 * ndfD       -> graus de liberdade das equacoes de difusao D1        *
 * ndfT       -> graus de liberdade das equacoes de tranporte T1      *
 * ndfF       -> graus de liberdade das equacoes de fluidos           *
 * ndfFt      -> graus de liberdade das equacoes de fluidos           *
 * nameOut    -> nome de arquivo de saida                             *
 * iws        -> vtk binario                                          *
 * f          -> arquivlo                                             *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 **********************************************************************/
void wGeoVtk(Memoria *m        ,double *x      
            ,INT *el           ,short *mat    
            ,short *nen        ,short *typeGeom
            ,double *prop      ,short *typeCal
            ,short *faceRd1    ,short *faceLd1
            ,short *faceRt1    ,short *faceLt1
            ,short *faceRfluid ,short *faceLfluid
            ,INT nnode         ,INT numel    
            ,short ndm      
            ,short maxNo       ,short maxViz  
            ,short numat    
            ,short *ndfD       ,short *ndfT
            ,short const ndfF  ,short const ndfFt
            ,char *nameOut     ,bool iws
            ,FILE *f)
{
  int    *lel=NULL;
  double *aux=NULL;
  INT i;
  short j;
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
      MAT2D(i,j,lel,maxNo) = MAT2D(i,j,el,maxNo)-1;
    }
  }  
  writeVtkCell(lel,nen,typeGeom,numel,maxNo,iws,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) numel);
/*...................................................................*/

/*... material*/
  HccaAlloc(int,m,lel,numel,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++)
    lel[i]=(int) mat[i];
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER,1,f);
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
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER,1,f);
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
   
  writeVtkProp(lel,&ddum,numel,1,"elTyCal",iws,INTEGER,1,f);
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
  writeVtkProp(&idum,aux,numel,MAXPROP,"elProp",iws,DOUBLEV,1,f);
  HccaDealloc(m,aux,"el",_AD_);
/*...................................................................*/

/*...*/
  if(ndfD[0] > 0 ){
/*... faceRd1*/
    HccaAlloc(int,m,lel,numel*(maxNo+1),"el",_AD_);
    if( lel == NULL){
      fprintf(stderr,"Erro na alocação de lel.\n"
                     "Nome do arquivo: %s.\n"
                    ,__FILE__);
      exit(EXIT_FAILURE);
    }
    for(i=0;i<numel*(maxNo+1);i++)
      lel[i]=(int) faceRd1[i];
   
    writeVtkProp(lel,&ddum,numel,maxViz+1,"faceRd1",iws,INTEGER,1,f);
    HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... faceLoadT1*/
    HccaAlloc(int,m,lel,numel*(maxNo+1),"el",_AD_);
    if( lel == NULL){
      fprintf(stderr,"Erro na alocação de lel.\n"
                     "Nome do arquivo: %s.\n"
                    ,__FILE__);
      exit(EXIT_FAILURE);
    }
    for(i=0;i<numel*(maxNo+1);i++)
      lel[i]=(int) faceLd1[i];
    writeVtkProp(lel,&ddum,numel,maxViz+1,"faceLd1",iws
                    ,INTEGER,1,f);
    HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  }
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
   
    writeVtkProp(lel,&ddum,numel,maxViz+1,"faceRt1",iws,INTEGER,1,f);
    HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... faceLoadT1*/
    HccaAlloc(int,m,lel,numel*(maxNo+1),"el",_AD_);
    if( lel == NULL){
      fprintf(stderr,"Erro na alocação de lel.\n"
                     "Nome do arquivo: %s.\n"
                    ,__FILE__);
      exit(EXIT_FAILURE);
    }
    for(i=0;i<numel*(maxNo+1);i++)
      lel[i]=(int) faceLt1[i];
    writeVtkProp(lel,&ddum,numel,maxViz+1,"faceLt1",iws
                    ,INTEGER,1,f);
    HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if(ndfF > 0 || ndfFt > 0){
/*... faceRfluid*/
    HccaAlloc(int,m,lel,numel*(maxNo+1),"el",_AD_);
    if( lel == NULL){
      fprintf(stderr,"Erro na alocação de lel.\n"
                     "Nome do arquivo: %s.\n"
                    ,__FILE__);
      exit(EXIT_FAILURE);
    }
    for(i=0;i<numel*(maxNo+1);i++)
      lel[i]=(int) faceRfluid[i];
   
    writeVtkProp(lel,&ddum,numel,maxViz+1,"faceRfuild",iws,INTEGER,1,f);
    HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... faceLoadFluid*/
    HccaAlloc(int,m,lel,numel*(maxNo+1),"el",_AD_);
    if( lel == NULL){
      fprintf(stderr,"Erro na alocação de lel.\n"
                     "Nome do arquivo: %s.\n"
                    ,__FILE__);
      exit(EXIT_FAILURE);
    }
    for(i=0;i<numel*(maxNo+1);i++)
      lel[i]=(int) faceLfluid[i];
    writeVtkProp(lel,&ddum,numel,maxViz+1,"faceLfluid",iws
                    ,INTEGER,1,f);
    HccaDealloc(m,lel,"el",_AD_);
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
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  fclose(f);
}
/*********************************************************************/ 

/**********************************************************************
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 22/08/2017                                    *
 *------------------------------------------------------------------- *
 * WGEOFACEVTK : escreve a malha apenas com os faces com condicao     *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * m          -> arranjo da menoria principal                         *  
 * x          -> coordenadas                                          * 
 * el         -> conectividade                                        * 
 * mat        -> materias                                             *  
 * nen        -> conectividades por elemento                          *  
 * mat        -> material por elemento                                *
 * typeGeom   -> tipo geometrico do elemento                          *
 * typeCal    -> tipo calculo do elemento                             *
 * faceRd1    -> condicao de contorno D1                              * 
 * faceLd1    -> tipo da condicao de contorno D1                      * 
 * faceRt1    -> condicao de contorno T1                              * 
 * faceLt1    -> tipo da condicao de contorno T1                      * 
 * faceRfluid -> condicao de contorno fluido                          * 
 * faceLfluid -> tipo da condicao de contorno fluido                  * 
 * nel        -> numeracao do elemento                                *
 * nnode      -> numero de nos                                        *  
 * numel      -> numero de elementos                                  *
 * ndm        -> numero de dimensao                                   *
 * maxNo      -> numero maximo de nos por elemento                    *
 * numat      -> numero maximo de nos por elemento                    *
 * ndfD1      -> graus de liberdade das equacoes difusao pura D1      *
 * ndfT1      -> graus de liberdade das equacoes de transporte T1     *
 * ndfF       -> graus de liberdade das equacoes de fluidos           *
 * ndfFt      -> graus de liberdade das equacoes de fluidos           *
 * nameOut    -> nome de arquivo de saida                             *
 * iws        -> vtk binario                                          *
 * f          -> arquivlo                                             *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 **********************************************************************/
void wGeoFaceVtk(Memoria *m       ,DOUBLE *x      
            ,INT *el              ,short *nen     
            ,short *typeGeom
            ,short *faceRd1       ,short *faceLd1
            ,short *faceRt1       ,short *faceLt1
            ,short *faceRfluid    ,short *faceLfluid
            ,INT const nnode      ,INT const numel    
            ,short const ndm      
            ,short const ndfD1    ,short const ndfT1
            ,short const ndfF     ,short const ndfFt
            ,short const maxViz   ,short const maxNo
            ,char *nameOut        ,bool iws
            ,FILE *f)
{
  char head[]={"FACE_VOLUME_FINITO"};
  int nFace=0;
  int *face = NULL,*idFace=NULL;
  int *lfaceLd1 = NULL,*lfaceLt1=NULL,*lfaceLfluid=NULL;
  short *typeGeomFace = NULL,*nenFace=NULL;
  int i;  
  int *aux=NULL;
  double ddum;

  HccaAlloc(INT   ,m,face ,numel*MAX_NUM_FACE*MAX_NUM_NODE_FACE
           ,"lFace"   ,_AD_);
  HccaAlloc(INT   ,m,idFace      ,numel*MAX_NUM_FACE,"iDFace"  ,_AD_);
  HccaAlloc(short ,m,typeGeomFace,numel*MAX_NUM_FACE,"ltGface" ,_AD_);
  HccaAlloc(short ,m,nenFace     ,numel*MAX_NUM_FACE,"lnenFace",_AD_);

/*... ndfD1*/
  if(ndfD1 > 0){
    HccaAlloc(int   ,m,lfaceLd1    ,numel*MAX_NUM_FACE,"lfaceSd1",_AD_);
    makeFace(el          ,faceRd1   ,faceLd1  
            ,typeGeom    
            ,face        ,lfaceLd1  ,idFace
            ,typeGeomFace,nenFace
            ,maxViz      ,maxNo
            ,ndfD1       
            ,numel       ,&nFace);
  }
/*...................................................................*/

/*... ndfT1*/
  if(ndfT1 > 0){
    HccaAlloc(int   ,m,lfaceLt1    ,numel*MAX_NUM_FACE,"lfaceSt1",_AD_);
    makeFace(el          ,faceRt1   ,faceLt1  
            ,typeGeom    
            ,face        ,lfaceLt1  ,idFace
            ,typeGeomFace,nenFace
            ,maxViz      ,maxNo
            ,ndfT1       
            ,numel       ,&nFace);
  }
/*...................................................................*/

/*... ndfF*/
  if(ndfF > 0 || ndfFt > 0){
    HccaAlloc(int,m,lfaceLfluid,numel*MAX_NUM_FACE,"lfaceSfluid",_AD_);
/*  makeFace(el          ,faceRfluid   ,faceLfluid
            ,typeGeom    
            ,face        ,lfaceLfluid  ,idFace
            ,typeGeomFace,nenFace
            ,maxViz      ,maxNo
            ,ndfT1       
            ,numel       ,&nFace);*/
  }
/*...................................................................*/

/*...*/
  if(nFace)
    if(iws)
      f = openFile(nameOut,"wb");
    else
      f = openFile(nameOut,"w");
/*...................................................................*/

/*... malha sem condicao de contorno na face*/
  else
    return;
/*...................................................................*/

/* ...*/
  headVtk(head,iws,f);
/* ..................................................................*/

/*... coordenadas*/
  writeVtkCoor(x,nnode,ndm,iws,f);
/*...................................................................*/
  
/*... faces*/
  writeVtkCell(face,nenFace,typeGeomFace,nFace,MAX_NUM_NODE_FACE,iws,f);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) nFace);
/*...................................................................*/
  
/*... relacao face celula*/
  writeVtkProp(idFace,&ddum,nFace,1,"idCellFace",iws,INTEGER,1,f);
/*...................................................................*/
  
/*... valores das cargas por celula*/
  if(ndfD1 > 0) 
    writeVtkProp(lfaceLd1,&ddum,nFace,1    ,"lFaceLd1",iws,INTEGER,1,f);
  if(ndfT1 > 0) 
    writeVtkProp(lfaceLt1,&ddum,nFace,1    ,"lFaceLt1",iws,INTEGER,1,f);
  if(ndfF  > 0 || ndfFt > 0) 
    writeVtkProp(lfaceLfluid  ,&ddum,nFace  ,1   
                ,"lFaceLfluid",iws  ,INTEGER,1   ,f);
/*...................................................................*/

/*.... campo por no*/
  fprintf(f,"POINT_DATA %ld\n",(long) nnode);
/*...................................................................*/

/*... numero do no*/
  HccaAlloc(int,m,aux,nnode,"el",_AD_);
  if( aux == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<nnode;i++)
    aux[i]=i+1;
   
  writeVtkProp(aux,&ddum,nnode,1,"pNode",iws,INTEGER,1,f);
  HccaDealloc(m,aux,"el",_AD_);
/*...................................................................*/

/*... dealloc*/
  if(ndfT1 > 0) 
    HccaDealloc(m,lfaceLt1    ,"lfaceSt1",_AD_);
  if(ndfD1 > 0) 
    HccaDealloc(m,lfaceLd1    ,"lfaceSd1",_AD_);
  if(ndfF  > 0) 
    HccaDealloc(m,lfaceLfluid ,"lfaceSfluid",_AD_);

  HccaDealloc(m,nenFace     ,"lnenFace",_AD_);
  HccaDealloc(m,typeGeomFace,"ltGface" ,_AD_);
  HccaDealloc(m,idFace      ,"iDFace"  ,_AD_);
  HccaDealloc(m,face        ,"lFace"   ,_AD_);
/*...................................................................*/
  
  fclose(f);
}
/*********************************************************************/ 

/********************************************************************** 
 * WRESVTK : escreve a malha com os resultados                        *  
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
  char head[]={"RES_VOLUME_FINITO"};
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
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER,1,f);
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
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... escrever resuldos por celula*/  
  if(ndfD[0] > 0){
    writeVtkProp(&idum,elU,numel,ndfD[0],"elD1",iws,DOUBLEV,1,f);
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
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... escrever resuldos por nos*/  
  if(ndfD[0] > 0){
    writeVtkProp(&idum,nU,nnode,ndfD[0],"nD1",iws,DOUBLEV,1,f);
  }
/*...................................................................*/

  fclose(f);
}
/*********************************************************************/

/********************************************************************** 
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 15/02/2017                                    * 
 * -------------------------------------------------------------------* 
 * WRESVTKDIF : escreve a malha com os resultados para problemas de   *  
 * difusao pura                                                       *  
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
 * elU     -> resultados por elementos                                *
 * nU      -> resultados por nos                                      *
 * elGradU -> gradientes dos resultados por elementos                 *
 * nGradU  -> gradientes dos resultados por nos                       *
 * nel     -> numeracao do elemento                                   *
 * nnode   -> numero de nos                                           *  
 * numel   -> numero de elementos                                     *
 * ndm     -> numero de dimensao                                      *
 * maxNo   -> numero maximo de nos por elemento                       *
 * numat   -> numero maximo de nos por elemento                       *
 * ndf     -> graus de liberdade das equacoes                         *
 * nameOut -> nome de arquivo de saida                                *
 * uResEl  -> nome dos resultados                                     *
 * uResNo  -> nome dos resultados                                     *
 * gradResEl -> nome dos resultados                                   *
 * gradResNo -> nome dos resultados                                   *
 * nameOut -> nome de arquivo de saida                                *
 * iws     -> vtk binario                                             *
 * f       -> arquivlo                                                *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * -------------------------------------------------------------------*
 * OBS:                                                               *
 * -------------------------------------------------------------------*
 **********************************************************************/
void wResVtkDif(Memoria *m        ,double *x      
               ,INT *el           ,short *mat    
               ,short *nen        ,short *typeGeom
               ,DOUBLE *elU       ,DOUBLE *nU
               ,DOUBLE *elGradU   ,DOUBLE *nGradU
               ,INT nnode         ,INT numel    
               ,short ndm         ,short maxNo 
               ,short numat       ,short ndf   
               ,char *uResEl      ,char *uResNo 
               ,char *gradResEl   ,char *gradResNo 
               ,char *nameOut     ,bool iws
               ,Temporal ddt      ,FILE *f)
{
  int    *lel=NULL;
  INT i;
  short j;
  char head[]={"DIF_VOLUME_FINITO"};
  double ddum;
  int    idum;

  if(iws)
    f = openFile(nameOut,"wb");
  else
    f = openFile(nameOut,"w");

/* ...*/
  headVtk(head,iws,f);
/* ..................................................................*/

/* ...*/
  if(ddt.flag)
    timeVtk(ddt.t,ddt.timeStep,iws,f);
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
      MAT2D(i,j,lel,maxNo) = MAT2D(i,j,el,maxNo)-1;
    }
  }  
  writeVtkCell(lel,nen,typeGeom,numel,maxNo,iws,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) numel);
/*...................................................................*/

/*... material*/
  HccaAlloc(int,m,lel,numel,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++)
    lel[i]=(int) mat[i];
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER,1,f);
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
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  
/*... escrever resultados por celula*/  
  writeVtkProp(&idum,elU,numel,ndf,uResEl,iws,DOUBLEV,1,f);
/*...................................................................*/

/*... escrever gradiente por celula*/  
  writeVtkProp(&idum,elGradU,numel,ndm,gradResEl,iws,DOUBLEV,2,f);
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
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... escrever resuldos por nos*/  
  writeVtkProp(&idum,nU    ,nnode,ndf,uResNo,iws,DOUBLEV,1,f);
/*...................................................................*/
  
/*... escrever os gradiente por nos*/  
  writeVtkProp(&idum,nGradU,nnode,ndm,gradResNo,iws,DOUBLEV,2,f);
/*...................................................................*/

  fclose(f);
}
/*********************************************************************/

/********************************************************************** 
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 12/09/2017                                   * 
 *-------------------------------------------------------------------* 
 * WRESVTKFLUID:escreve a malha com os resultados para problemas de   *  
 * de escomentos de fluidos imcompressivel                            *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * m            -> arranjo da menoria principal                       *  
 * x            -> coordenadas                                        * 
 * el           -> conectividade                                      * 
 * mat          -> materias                                           *  
 * nen          -> conectividades por elemento                        *  
 * mat          -> material por elemento                              *
 * typeGeom     -> tipo geometrico do elemento                        *
 * elPres       -> resultados por elementos                           *
 * nPres        -> resultados por nos                                 *
 * elGradPres   -> gradientes dos resultados por elementos            *
 * nGradPres    -> gradientes dos resultados por nos                  *
 * elVel        -> campo de velocidade por elementos                  *
 * nVel         -> campo de velocidade por nos                        *
 * elGradVel    -> gradientes dos resultados por elementos            *
 * nGradVel     -> gradientes dos resultados por nos                  *
 * nel          -> numeracao do elemento                              *
 * nnode        -> numero de nos                                      *  
 * numel        -> numero de elementos                                *
 * ndm          -> numero de dimensao                                 *
 * maxNo        -> numero maximo de nos por elemento                  *
 * numat        -> numero maximo de nos por elemento                  *
 * ndf          -> graus de liberdade das equacoes                    *
 * nameOut      -> nome de arquivo de saida                           *
 * presResEl    -> nome dos resultados                                *
 * presResNo    -> nome dos resultados                                *
 * presgradResEl-> nome dos resultados                                *
 * presgradResNo-> nome dos resultados                                *
 * nameOut      -> nome de arquivo de saida                           *
 * iws          -> vtk binario                                        *
 * fVel         -> escreve as velocidades                             *
 * fGradVel     -> escreve os gradientes da velocidades               *
 * fPres        -> escreve as pressoes                                *
 * fGradPes     -> escrever os gradientes da pressoes                 *
 * feddyViscosity> escrever a viscosidade turbulenta                  *
 * f            -> arquivlo                                           *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * OBS:                                                               * 
 *------------------------------------------------------------------- * 
 *                                                                    *
 *           | du1dx1 du1dx2 du1dx3 |                                 *
 * gradVel = | du2dx1 du2dx2 du2dx3 |                                 * 
 *           | du3dx1 du3dx2 du3dx3 |                                 *
 *                                                                    *
 **********************************************************************/
void wResVtkFluid(Memoria *m    ,DOUBLE *x      
          ,INT *el              ,short *mat    
          ,short *nen           ,short *typeGeom
          ,DOUBLE *elPres       ,DOUBLE *nPres
          ,DOUBLE *elGradPres   ,DOUBLE *nGradPres
          ,DOUBLE *elVel        ,DOUBLE *nVel      
          ,DOUBLE *elGradVel    ,DOUBLE *nGradVel 
          ,DOUBLE *elEnergy     ,DOUBLE *nEnergy
          ,DOUBLE *elGradEnergy ,DOUBLE *nGradEnergy
          ,DOUBLE *eddyVis
          ,INT nnode            ,INT numel    
          ,short const ndm      ,short const maxNo 
          ,short const numat    ,short const ndf   
          ,char *presResEl      ,char *presResNo 
          ,char *gradPresResEl  ,char *gradPresResNo 
          ,char *velEl          ,char *velNo       
          ,char *gradVelResEl   ,char *gradVelResNo 
          ,char *energyResEl    ,char *energyResNo
          ,char *gradEnergyResEl,char *gradEnergyResNo
          ,char *eddyVisRes
          ,char *nameOut        ,bool iws
          ,bool fVel            ,bool fGradVel 
          ,bool fPres           ,bool fGradPres
          ,bool fEnergy         ,bool fGradEnergy
          ,bool fEddyViscosity  ,bool fKelvin
          ,Temporal ddt         ,FILE *f)
{
  int    *lel=NULL;
  DOUBLE *p=NULL;
  INT i;
  short j;
  char head[]={"FLUID_VOLUME_FINITO"};
  double ddum;
  int    idum;

  if(iws)
    f = openFile(nameOut,"wb");
  else
    f = openFile(nameOut,"w");

/* ...*/
  headVtk(head,iws,f);
/* ..................................................................*/

/* ...*/
  if(ddt.flag)
    timeVtk(ddt.t,ddt.timeStep,iws,f);
/* ..................................................................*/

/*... coordenadas*/
  writeVtkCoor(x,nnode,ndm,iws,f);
/*...................................................................*/
  
/*... conectividades*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  ERRO_MALLOC(lel,"el",__LINE__,__FILE__,__func__)
  for(i=0;i<numel;i++){
    for(j=0;j<maxNo;j++){
      MAT2D(i,j,lel,maxNo) = MAT2D(i,j,el,maxNo)-1;
    }
  }  
  writeVtkCell(lel,nen,typeGeom,numel,maxNo,iws,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) numel);
/*...................................................................*/

/*... material*/
  HccaAlloc(int,m,lel,numel,"el",_AD_);
  ERRO_MALLOC(lel,"el",__LINE__,__FILE__,__func__)
  for(i=0;i<numel;i++)
    lel[i]=(int) mat[i];
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... numero do elemento*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  ERRO_MALLOC(lel,"el",__LINE__,__FILE__,__func__)
  for(i=0;i<numel;i++)
    lel[i]= i+1;
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  
/*... escrever resultados de pressao por celula*/ 
  if(fPres)
    writeVtkProp(&idum,elPres,numel,1,presResEl,iws,DOUBLEV,1,f);
/*...................................................................*/

/*... escrever gradiente da pressao por celula*/
  if(fGradPres)  
    writeVtkProp(&idum,elGradPres,numel,ndm,gradPresResEl
                ,iws,DOUBLEV,2,f);
/*...................................................................*/


/*... escrever campo de velociade por celula*/  
  if(fVel)
    writeVtkProp(&idum,elVel,numel,ndm,velEl,iws,DOUBLEV,2,f);
/*...................................................................*/

/*... escrever gradiente de velocidade por celula*/  
  if(fGradVel){  
    if( ndm == 2) 
      writeVtkProp(&idum,elGradVel,numel,2*ndm,gradVelResEl
                  ,iws,DOUBLEV,1,f);
    else if( ndm ==3 )
      writeVtkProp(&idum,elGradVel,numel,3*ndm,gradVelResEl
                  ,iws,DOUBLEV,1,f);
  }
/*...................................................................*/

/*... escrever campo de energia por celula*/
  if (fEnergy){
    if(fKelvin ){
      HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
      ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__)
      alphaProdVector(1.e0,elEnergy,numel,p);
      convTempForKelvin(p, numel,false);
      writeVtkProp(&idum, p, numel,   1, energyResEl, iws
                   , DOUBLEV, 1, f);
      HccaDealloc(m,p,"p",_AD_);
    }
    else
      writeVtkProp(&idum, elEnergy, numel,   1, energyResEl, iws
                   , DOUBLEV, 1, f);
   }
/*...................................................................*/

/*... escrever gradiente de velocidade por celula*/
  if (fGradEnergy) {
    if (ndm == 2)
      writeVtkProp(&idum, elGradEnergy, numel,ndm, gradEnergyResEl
                   ,iws, DOUBLEV, 2, f);
  }
/*...................................................................*/

/*... escrever gradiente de velocidade por celula*/  
  if(fEddyViscosity)  
    writeVtkProp(&idum,eddyVis,numel,1,eddyVisRes
                ,iws,DOUBLEV,1,f);
/*...................................................................*/


/*.... campo por no*/
  fprintf(f,"POINT_DATA %ld\n",(long) nnode);
/*...................................................................*/

/*... numero do no*/
  HccaAlloc(int,m,lel,nnode,"el",_AD_);
  ERRO_MALLOC(lel,"el",__LINE__,__FILE__,__func__)
  for(i=0;i<nnode;i++)
    lel[i]=i+1;
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... escrever resultados de pressao por nos*/  
  if(fPres)
    writeVtkProp(&idum,nPres ,nnode,1,presResNo,iws,DOUBLEV,1,f);
/*...................................................................*/
  
/*... escrever os gradiente de pressao por nos*/  
  if(fGradPres)  
    writeVtkProp(&idum,nGradPres,nnode,ndm,gradPresResNo
                ,iws,DOUBLEV,2,f);
/*...................................................................*/

/*... escrever as velocidade por nos*/  
  if(fVel)
    writeVtkProp(&idum,nVel,nnode,ndm,velNo,iws,DOUBLEV,2,f);
/*...................................................................*/

/*... escrever as gradiente de velocidade por nos*/ 
  if(fGradVel){  
    if( ndm == 2) 
      writeVtkProp(&idum,nGradVel,nnode,2*ndm,gradVelResNo
                  ,iws,DOUBLEV,1,f);
    else if( ndm == 3) 
      writeVtkProp(&idum,nGradVel,nnode,3*ndm,gradVelResNo
                  ,iws,DOUBLEV,1,f);
  }
/*...................................................................*/

/*... escrever resultados de energia por nos*/
  if (fEnergy){
    if(fKelvin ){
      HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
      ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__)
      alphaProdVector(1.e0,nEnergy,nnode,p);
      convTempForKelvin(p, nnode,false);
      writeVtkProp(&idum, p, nnode, 1, energyResNo, iws
                , DOUBLEV, 1, f);
      HccaDealloc(m,p,"p",_AD_);
    }
    else
      writeVtkProp(&idum, nEnergy, nnode, 1, energyResNo, iws
                  , DOUBLEV, 1, f);     
  }
/*...................................................................*/

/*... escrever os gradiente de pressao por nos*/
  if (fGradEnergy)
    writeVtkProp(&idum, nGradEnergy, nnode, ndm, gradEnergyResNo
                 ,iws, DOUBLEV, 2, f);
/*...................................................................*/

  fclose(f);
}
/*********************************************************************/

/********************************************************************** 
 * WRESVTKTRANS:escreve a malha com os resultados para problemas de   *  
 * transporte                                                         *  
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
 * elU     -> resultados por elementos                                *
 * nU      -> resultados por nos                                      *
 * elGradU -> gradientes dos resultados por elementos                 *
 * nGradU  -> gradientes dos resultados por nos                       *
 * elVel   -> campo de velocidade por elementos                       *
 * nVel    -> campo de velocidade por nos                             *
 * nel     -> numeracao do elemento                                   *
 * nnode   -> numero de nos                                           *  
 * numel   -> numero de elementos                                     *
 * ndm     -> numero de dimensao                                      *
 * maxNo   -> numero maximo de nos por elemento                       *
 * numat   -> numero maximo de nos por elemento                       *
 * ndf     -> graus de liberdade das equacoes                         *
 * nameOut -> nome de arquivo de saida                                *
 * uResEl  -> nome dos resultados                                     *
 * uResNo  -> nome dos resultados                                     *
 * gradResEl -> nome dos resultados                                   *
 * gradResNo -> nome dos resultados                                   *
 * nameOut -> nome de arquivo de saida                                *
 * iws     -> vtk binario                                             *
 * f       -> arquivlo                                                *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 **********************************************************************/
void wResVtkTrans(Memoria *m        ,double *x      
                 ,INT *el           ,short *mat    
                 ,short *nen        ,short *typeGeom
                 ,DOUBLE *elU       ,DOUBLE *nU
                 ,DOUBLE *elGradU   ,DOUBLE *nGradU
                 ,DOUBLE *elVel     ,DOUBLE *nVel      
                 ,INT nnode         ,INT numel    
                 ,short ndm         ,short maxNo 
                 ,short numat       ,short ndf   
                 ,char *uResEl      ,char *uResNo 
                 ,char *gradResEl   ,char *gradResNo 
                 ,char *velEl       ,char *velNo       
                 ,char *nameOut     ,bool iws
                 ,Temporal ddt      ,FILE *f)
{
  int    *lel=NULL;
  INT i;
  short j;
  char head[]={"DIF_VOLUME_FINITO"};
  double ddum;
  int    idum;

  if(iws)
    f = openFile(nameOut,"wb");
  else
    f = openFile(nameOut,"w");

/* ...*/
  headVtk(head,iws,f);
/* ..................................................................*/

/* ...*/
  if(ddt.flag)
    timeVtk(ddt.t,ddt.timeStep,iws,f);
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
      MAT2D(i,j,lel,maxNo) = MAT2D(i,j,el,maxNo)-1;
    }
  }  
  writeVtkCell(lel,nen,typeGeom,numel,maxNo,iws,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) numel);
/*...................................................................*/

/*... material*/
  HccaAlloc(int,m,lel,numel,"el",_AD_);
  if( lel == NULL){
    fprintf(stderr,"Erro na alocação de lel.\n"
                   "Nome do arquivo: %s.\n"
                  ,__FILE__);
    exit(EXIT_FAILURE);
  }
  for(i=0;i<numel;i++)
    lel[i]=(int) mat[i];
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER,1,f);
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
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  
/*... escrever resultados por celula*/  
  writeVtkProp(&idum,elU,numel,ndf,uResEl,iws,DOUBLEV,1,f);
/*...................................................................*/

/*... escrever gradiente por celula*/  
  writeVtkProp(&idum,elGradU,numel,ndm,gradResEl,iws,DOUBLEV,2,f);
/*...................................................................*/

/*... escrever campo de velociade  por celula*/  
  writeVtkProp(&idum,elVel,numel,ndm,velEl,iws,DOUBLEV,2,f);
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
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... escrever resuldos por nos*/  
  writeVtkProp(&idum,nU    ,nnode,ndf,uResNo,iws,DOUBLEV,1,f);
/*...................................................................*/
  
/*... escrever os gradiente por nos*/  
  writeVtkProp(&idum,nGradU,nnode,ndm,gradResNo,iws,DOUBLEV,2,f);
/*...................................................................*/

/*... escrever as velocidade por nos*/  
  writeVtkProp(&idum,nVel,nnode,ndm,velNo,iws,DOUBLEV,2,f);
/*...................................................................*/

  fclose(f);
}
/*********************************************************************/

/********************************************************************** 
 * MAKEFACE : gera as faces/aresta onde ha carregamento               *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * el      -> conectividade                                           * 
 * faceR   -> tipo de condicao de contorno                            * 
 * faceS   -> valor de condicao de contorno                           * 
 * typeGeom-> tipo geometrico do elemento                             *
 * face    -> indefinido                                              *
 * lfaceL  -> indefinido                                              *
 * idFace  -> indefinido                                              *
 * tyGeomF -> indefinido                                              *
 * nenFace -> indefinido                                              *
 * maxViz  -> numero maximo de vizinho por celula                     *
 * maxNo   -> numero maximo de nos por celula                         *
 * numel   -> numero de elementos                                     *
 * nFace   -> indefinido                                              *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * face    -> faces onde ha condicao de contorno                      *
 * lFaceL  -> tipo de condico de contorno na face                     *
 * tyGeomF -> tipo geometrico da face                                 *
 * nenFace -> numero de nos por face                                  *
 * nFace   -> numero total de faces                                   *
 * ------------------------------------------------------------------ *
 **********************************************************************/
void makeFace(INT *el            ,short *faceR       ,short *faceL 
             ,short *typeGeom
             ,INT *face          ,int    *lFaceL     ,INT *idFace
             ,short *typeGeomFace,short *nenFace
             ,short const maxViz ,short const maxNo
             ,short const ndf     
             ,INT const numel    ,INT *nFace){

  short  isnod[MAX_SN],nenFaceMax=MAX_NUM_NODE_FACE;
  short j,k,ty,tmp;
  unsigned int nf = 0;
  short cCell = maxViz + 1;
  int no,nel;
  
  for(nel=0;nel<numel;nel++){
    ty = typeGeom[nel];

/*... triangulos*/
    if( ty == TRIACELL){
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for(j=0;j<maxViz;j++)
        tmp += abs(MAT2D(nel,j,faceR,cCell));
      if(tmp){
/*... loop nas arestas*/
        for(j=0;j<3;j++){
          tmp = MAT2D(nel,j,faceR,cCell);
/*...*/
          if(tmp){
            idFace[nf]       = nel + 1;
            typeGeomFace[nf] = LINECELL;
            nenFace[nf]      = sn(isnod,ty,nel);
            lFaceL[nf]       = MAT2D(nel,j,faceL,cCell); 
/*... loop nos nos da face*/
            for(k=0;k<nenFace[nf];k++){
              no  = MAT2D(j,k,isnod,nenFace[nf]);
              MAT2D(nf,k,face,nenFaceMax) = MAT2D(nel,no,el,maxNo)-1;
            }    
            nf++;
          } 
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... quadrilatero*/
    else if( ty == QUADCELL){
/*...*/
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for(j=0;j<maxViz;j++)
        tmp += abs(MAT2D(nel,j,faceR,cCell));
/*...*/
      if(tmp){
/*... loop nas arestas*/
        for(j=0;j<4;j++){
          tmp = MAT2D(nel,j,faceR,cCell);
/*...*/
          if(tmp){
            idFace[nf]       = nel + 1;
            typeGeomFace[nf] = LINECELL;
            nenFace[nf]      = sn(isnod,ty,nel);
            lFaceL[nf]       = MAT2D(nel,j,faceL,cCell); 
/*... loop nos nos da face*/
            for(k=0;k<nenFace[nf];k++){
              no  = MAT2D(j,k,isnod,nenFace[nf]);
              MAT2D(nf,k,face,nenFaceMax) = MAT2D(nel,no,el,maxNo)-1;
            }    
            nf++;
          } 
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... tetraedro*/
    else if( ty == TETRCELL){
/*...*/
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for(j=0;j<maxViz;j++)
        tmp += abs(MAT2D(nel,j,faceR,cCell));
/*...*/
      if(tmp){
/*... loop nas faces*/
        for(j=0;j<4;j++){
          tmp = MAT2D(nel,j,faceR,cCell);
/*...*/
          if(tmp){
            idFace[nf]       = nel + 1;
            typeGeomFace[nf] = TRIACELL;
            nenFace[nf]      = sn(isnod,ty,nel);
            lFaceL[nf]       = MAT2D(nel,j,faceL,cCell); 
/*... loop nos nos da face*/
            for(k=0;k<nenFace[nf];k++){
              no  = MAT2D(j,k,isnod,nenFace[nf]);
              MAT2D(nf,k,face,nenFaceMax) = MAT2D(nel,no,el,maxNo)-1;
            }    
            nf++;
          } 
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... hexaedro*/
    else if( ty == HEXACELL){
/*...*/
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for(j=0;j<maxViz;j++)
        tmp += abs(MAT2D(nel,j,faceR,cCell));
/*...*/
      if(tmp){
/*... loop nas faces*/
        for(j=0;j<6;j++){
          tmp = MAT2D(nel,j,faceR,cCell);
/*...*/
          if(tmp){
            idFace[nf]       = nel + 1;
            typeGeomFace[nf] = QUADCELL;
            nenFace[nf]      = sn(isnod,ty,nel); 
            lFaceL[nf]       = MAT2D(nel,j,faceL,cCell); 
/*... loop nos nos da face*/
            for(k=0;k<nenFace[nf];k++){
              no  = MAT2D(j,k,isnod,nenFace[nf]);
              MAT2D(nf,k,face,nenFaceMax) = MAT2D(nel,no,el,maxNo)-1;
            }    
            nf++;
          } 
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
  
  *nFace = nf;

}
/*********************************************************************/ 
