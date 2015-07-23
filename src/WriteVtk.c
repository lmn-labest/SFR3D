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
 * faceRd1 -> tipo de condicao de contorno D1                         * 
 * faceSd1 -> valor de condicao de contorno D1                        * 
 * nnode   -> numero de nos                                           *  
 * numel   -> numero de elementos                                     *
 * ndm     -> numero de dimensao                                      *
 * maxNo   -> numero maximo de nos por elemento                       *
 * maxNo   -> numero maximo de vizinhos por elemento                  *
 * numat   -> numero de materias                                      *
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
            ,short *faceRd1 ,double *faceSd1
            ,INT nnode      ,INT numel    
            ,short ndm      
            ,short maxNo    ,short maxViz  
            ,short numat    ,short *ndfT   
            ,char *nameOut  ,bool iws
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
  if(ndfT[0] > 0 ){
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

/*... faceSt1*/
    writeVtkProp(&idum,faceSd1,numel,maxViz+1,"faceSd1",iws
                    ,DOUBLEV,1,f);
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
 * WGEOFACEVTK : escreve a malha apenas com os faces com condicao     *  
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
 * faceRd1 -> tipo de condicao de contorno T1                         * 
 * faceSd1 -> valor de condicao de contorno T1                        * 
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
void wGeoFaceVtk(Memoria *m       ,DOUBLE *x      
            ,INT *el              ,short *nen     
            ,short *typeGeom
            ,short *faceRd1       ,DOUBLE *faceSd1
            ,INT const nnode      ,INT const numel    
            ,short const ndm      ,short const ndf
            ,short const maxViz   ,short const maxNo
            ,char *nameOut        ,bool iws
            ,FILE *f)
{
  char head[]={"FACE_VOLUME_FINITO"};
  int nFace=0;
  int *face = NULL,*idFace=NULL;
  double *lfaceSd1 = NULL;
  short *typeGeomFace = NULL,*nenFace=NULL;
  int i,idum;  
  int *aux=NULL;
  double ddum;

  HccaAlloc(INT   ,m,face ,numel*MAX_NUM_FACE*MAX_NUM_NODE_FACE
           ,"lFace"   ,_AD_);
  HccaAlloc(INT   ,m,idFace      ,numel*MAX_NUM_FACE,"iDFace"  ,_AD_);
  HccaAlloc(short ,m,typeGeomFace,numel*MAX_NUM_FACE,"ltGface" ,_AD_);
  HccaAlloc(short ,m,nenFace     ,numel*MAX_NUM_FACE,"lnenFace",_AD_);
  HccaAlloc(DOUBLE,m,lfaceSd1    ,numel*MAX_NUM_FACE*ndf 
                                 ,"lfaceSd1",_AD_);
  
  makeFace(el          ,faceRd1   ,faceSd1  
          ,typeGeom    
          ,face        ,lfaceSd1  ,idFace
          ,typeGeomFace,nenFace
          ,maxViz      ,maxNo
          ,ndf         
          ,numel       ,&nFace);
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
  writeVtkProp(&idum,lfaceSd1,nFace,ndf,"lFaceS1",iws,DOUBLEV,1,f);
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
  HccaDealloc(m,lfaceSd1    ,"lfaceSd1",_AD_);
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
               ,FILE *f)
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
 * MAKEFACE : gera as faces/aresta onde ha carregamento               *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * el      -> conectividade                                           * 
 * faceR   -> tipo de condicao de contorno                            * 
 * faceS   -> valor de condicao de contorno                           * 
 * typeGeom-> tipo geometrico do elemento                             *
 * face    -> indefinido                                              *
 * lfaceS  -> indefinido                                              *
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
 * lFaceS  -> valor de contorno por face                              *
 * tyGeomF -> tipo geometrico da face                                 *
 * nenFace -> numero de nos por face                                  *
 * nFace   -> numero total de faces                                   *
 * ------------------------------------------------------------------ *
 **********************************************************************/
void makeFace(INT *el            ,short *faceR       ,DOUBLE *faceS 
             ,short *typeGeom
             ,INT *face          ,DOUBLE *lFaceS     ,INT *idFace
             ,short *typeGeomFace,short *nenFace
             ,short const maxViz ,short const maxNo
             ,short const ndf     
             ,INT const numel    ,INT *nFace){

  short  isnod[MAX_SN],nenFaceMax=MAX_NUM_NODE_FACE;
  short i,j,k,ty,tmp;
  unsigned INT nf = 0;
  short cCell = maxViz + 1;
  int no,nel;
  
  for(nel=0;nel<numel;nel++){
    ty = typeGeom[nel];

/*... triangulos*/
    if( ty == TRIACELL){
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for(j=0;j<maxViz;j++)
        tmp += MAT2D(nel,j,faceR,cCell);
      if(tmp){
/*... loop nas arestas*/
        for(j=0;j<3;j++){
          tmp = MAT2D(nel,j,faceR,cCell);
/*...*/
          if(tmp){
            idFace[nf]       = nel + 1;
            typeGeomFace[nf] = LINECELL;
            nenFace[nf]      = sn(isnod,ty,nel);
            for(i=0;i<ndf;i++)
              MAT2D(nf,i,lFaceS,ndf) = MAT3D(nel,j,i,faceS,cCell,ndf); 
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
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for(j=0;j<maxViz;j++)
        tmp += MAT2D(nel,j,faceR,cCell);
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
            for(i=0;i<ndf;i++)
              MAT2D(nf,i,lFaceS,ndf) = MAT3D(nel,j,i,faceS,cCell,ndf); 
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
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for(j=0;j<maxViz;j++)
        tmp += MAT2D(nel,j,faceR,cCell);
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
            for(i=0;i<ndf;i++)
              MAT2D(nf,i,lFaceS,ndf) = MAT3D(nel,j,i,faceS,cCell,ndf); 
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
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for(j=0;j<maxViz;j++)
        tmp += MAT2D(nel,j,faceR,cCell);
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
            for(i=0;i<ndf;i++)
              MAT2D(nf,i,lFaceS,ndf) = MAT3D(nel,j,i,faceS,cCell,ndf); 
/*... loop nos nos da face*/
            for(k=0;k<1;k++){
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
