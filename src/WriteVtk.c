#include<WriteVtk.h>

/********************************************************************** 
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 00/00/0000                                    * 
 * -------------------------------------------------------------------* 
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
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER_VTK,1,f);
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
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER_VTK,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  fclose(f);
}
/*********************************************************************/ 

/********************************************************************** 
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 00/00/0000                                    *
 * ------------------------------------------------------------------ *
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
  int *lel=NULL;
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
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER_VTK,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

  writeVtkProp(ep,&ddum,numel,1,"elPart"  ,iws,INTEGER_VTK,1,f);

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
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER_VTK,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... numero do no*/
  writeVtkProp(np,&ddum,nnode,1,"noPart",iws,INTEGER_VTK,1,f);
/*...................................................................*/
  fclose(f);
}
/*********************************************************************/ 

/********************************************************************** 
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 06/10/2019                                    *
 *------------------------------------------------------------------- *
 * WGEOVTK : escreve a malha com os resultados e condicao             *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * m          -> arranjo da menoria principal                         *  
 * x          -> coordenadas                                          *
 * cc           -> centro geomentrico das celulas                     *
 * el         -> conectividade                                        * 
 * mat        -> materias                                             *  
 * nen        -> conectividades por elemento                          *  
 * mat        -> material por elemento                                *
 * typeGeom   -> tipo geometrico do elemento                          *
 * typeCal    -> tipo calculo do elemento                             *
 * nel        -> numeracao do elemento                                *
 * faceRd1    -> condicao de contorno D1                              * 
 * faceRt1    -> condicao de contorno T1                              * 
 * faceRfluid -> condicao de contorno fluido                          * 
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
            ,short *faceRd1    ,short *faceRt1
            ,short *faceRfluid ,short *faceRenergy
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
  double *daux=NULL;
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
  writeVtkCoor(x, nnode, ndm, iws, f);
/* ..................................................................*/
 
/*... conectividades*/
  HccaAlloc(int,m,lel,numel*maxNo,"lel",_AD_);
  ERRO_MALLOC(lel,"lel",__LINE__,__FILE__,__func__);
  for(i=0;i<numel;i++)
  {
    for(j=0;j<maxNo;j++)
    {
      MAT2D(i,j,lel,maxNo) = MAT2D(i,j,el,maxNo)-1;
    }
  }
  writeVtkCell(lel,nen,typeGeom,numel,maxNo,iws,f);
  HccaDealloc(m, lel ,"lel"   ,_AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) numel );
/*...................................................................*/

/*... material*/
  HccaAlloc(int,m,lel,numel,"el",_AD_);
  ERRO_MALLOC(lel,"lel",__LINE__,__FILE__,__func__);
  for(i=0;i<numel;i++)
    lel[i]=(int) mat[i];
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER_VTK,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... numero do elemento*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  ERRO_MALLOC(lel,"lel",__LINE__,__FILE__,__func__);
  
  for(i=0;i<numel;i++)
    lel[i]= i+1;
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER_VTK,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... tipo celula para o calculo*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  ERRO_MALLOC(lel,"lel",__LINE__,__FILE__,__func__);
  
  for(i=0;i<numel;i++){
    idum = mat[i]-1;
    lel[i] = typeCal[idum];
  }
   
  writeVtkProp(lel,&ddum,numel,1,"elTyCal",iws,INTEGER_VTK,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... propriedades dos elementos*/
  HccaAlloc(double,m,daux,numel*numat,"el",_AD_);
  ERRO_MALLOC(daux,"daux",__LINE__,__FILE__,__func__);
  for(i=0;i<numel;i++){
    for(j=0;j<MAXPROP;j++){
      idum = mat[i]-1;
      MAT2D(i,j,daux,MAXPROP) = MAT2D(idum,j,prop,MAXPROP);
    }
  } 
  writeVtkProp(&idum,daux,numel,MAXPROP,"elProp",iws,DOUBLE_VTK,1,f);
  HccaDealloc(m,daux,"el",_AD_);
/*...................................................................*/

/*...*/
  if(ndfD[0] > 0 ){
/*... faceRd1*/
    HccaAlloc(int,m,lel,numel*(maxNo+1),"el",_AD_);
    ERRO_MALLOC(lel,"lel",__LINE__,__FILE__,__func__);
    for(i=0;i<numel*(maxNo+1);i++)
      lel[i]=(int) faceRd1[i];
   
    writeVtkProp(lel,&ddum,numel,maxViz+1,"faceRd1",iws,INTEGER_VTK,1,f);
    HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

  }
/*...................................................................*/

/*...*/
  if(ndfT[0] > 0 ){
/*... faceRt1*/
    HccaAlloc(int,m,lel,numel*(maxNo+1),"el",_AD_);
    ERRO_MALLOC(lel,"lel",__LINE__,__FILE__,__func__);
    for(i=0;i<numel*(maxNo+1);i++)
      lel[i]=(int) faceRt1[i];
   
    writeVtkProp(lel,&ddum,numel,maxViz+1,"faceRt1",iws,INTEGER_VTK,1,f);
    HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

  }
/*...................................................................*/

/*...*/
  if(ndfF > 0){
/*... faceRfluid*/
    HccaAlloc(int,m,lel,numel*(maxNo+1),"el",_AD_);
    ERRO_MALLOC(lel,"lel",__LINE__,__FILE__,__func__);
    for(i=0;i<numel*(maxNo+1);i++)
      lel[i]=(int) faceRfluid[i];
   
    writeVtkProp(lel,&ddum,numel,maxViz+1,"faceRfuild",iws,INTEGER_VTK,1,f);
    HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

  }
/*...................................................................*/

/*...*/
  if(ndfFt > 0){
/*... faceRfluid*/
    HccaAlloc(int,m,lel,numel*(maxNo+1),"el",_AD_);
    ERRO_MALLOC(lel,"lel",__LINE__,__FILE__,__func__);
    for(i=0;i<numel*(maxNo+1);i++)
      lel[i]=(int) faceRfluid[i];
   
    writeVtkProp(lel,&ddum,numel,maxViz+1,"faceRfuild",iws,INTEGER_VTK,1,f);
    HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... faceRenergy*/
    HccaAlloc(int,m,lel,numel*(maxNo+1),"el",_AD_);
    ERRO_MALLOC(lel,"lel",__LINE__,__FILE__,__func__);
    for(i=0;i<numel*(maxNo+1);i++)
      lel[i]=(int) faceRenergy[i];
   
    writeVtkProp(lel,&ddum,numel,maxViz+1,"faceRtemp",iws,INTEGER_VTK,1,f);
    HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

  }
/*...................................................................*/

/*.... campo por no*/
  fprintf(f,"POINT_DATA %ld\n",(long) nnode);
/*...................................................................*/

/*... numero do no*/
  HccaAlloc(int,m,lel,nnode,"el",_AD_);
  ERRO_MALLOC(lel,"lel",__LINE__,__FILE__,__func__);
  for(i=0;i<nnode;i++)
    lel[i]=i+1;
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER_VTK,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  fclose(f);
}
/*********************************************************************/ 

/**********************************************************************
* Data de criacao    : 00/00/0000                                    *
* Data de modificaco : 30/04/2018                                    *
*------------------------------------------------------------------- *
* WGEOVTK : escreve a malha com os resultados e condicao             *
* ------------------------------------------------------------------ *
* parametros de entrada:                                             *
* ------------------------------------------------------------------ *
* m          -> arranjo da menoria principal                         *
* x          -> coordenadas                                          *
* cc           -> centro geomentrico das celulas                     *
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
void wGeoVtk2(Memoria *m        , DOUBLE *x
            , DOUBLE *cc        , INT *el         
            , short *nen        , short *typeGeom
            , INT nnode         , INT numel
            , short ndm         , short maxNo       
            , short maxViz      
            , char *nameOut     , bool iws
            , FILE *f)
{
  short  *iaux1 = NULL,*iaux2 = NULL;
  int    *lel = NULL;
  double *daux = NULL;
  INT i;
  short j;
  char head[] = { "GEOM_VOLUME_FINITO" };

  if (iws)
    f = openFile(nameOut, "wb");
  else
    f = openFile(nameOut, "w");

/* ...*/
  headVtk(head, iws, f);
/* ..................................................................*/

/*... coordenadas*/
  HccaAlloc(double, m, daux, (nnode+numel)*ndm, "daux", _AD_);
  ERRO_MALLOC(daux, "daux", __LINE__, __FILE__, __func__);
  for (i = 0; i< nnode; i++)
  {
    for (j = 0; j<ndm; j++)
    {
      MAT2D(i + nnode, j, daux, ndm) = MAT2D(i, j, cc, ndm);
      MAT2D(i, j, daux, ndm) = MAT2D(i, j, x, ndm);
    }
  }

  writeVtkCoor(daux, nnode + numel, ndm, iws, f);
  HccaDealloc(m, daux, "daux", _AD_);
/* ..................................................................*/

/*... conectividades*/
  HccaAlloc(int, m, lel, (numel + numel) * maxNo, "lel", _AD_);
  ERRO_MALLOC(lel, "lel", __LINE__, __FILE__, __func__);
  HccaAlloc(short, m, iaux1, numel + numel, "laux1", _AD_);
  ERRO_MALLOC(iaux1, "laux1", __LINE__, __FILE__, __func__);
  HccaAlloc(short, m, iaux2, numel + numel, "laux2", _AD_);
  ERRO_MALLOC(iaux2, "laux2", __LINE__, __FILE__, __func__);
  for (i = 0; i<numel; i++) {
    iaux1[i] = nen[i];
    iaux2[i] = typeGeom[i];
    for (j = 0; j<maxNo; j++) {
      MAT2D(i, j, lel, maxNo) = MAT2D(i, j, el, maxNo) - 1;
    }
  }

  for (i = 0; i<numel; i++) {
    iaux1[numel+i] = 1;
    iaux2[numel+i] = DOTCELL;
    MAT2D(numel + i, 0, lel, maxNo) = nnode + i;
  }
  writeVtkCell(lel, iaux1, iaux2, 2*numel, maxNo, iws, f);
  HccaDealloc(m, iaux2, "laux2", _AD_);
  HccaDealloc(m, iaux1, "laux1", _AD_);
  HccaDealloc(m, lel  , "lel"   , _AD_);
/*...................................................................*/
  fclose(f);
}
/*********************************************************************/

/**********************************************************************
* Data de criacao    : 20/07/2018                                    *
* Data de modificaco : 07/11/2019                                    *
*------------------------------------------------------------------- *
* WGEOFACEVTK2: escreve a malha apenas com os faces com condicao     *
* ------------------------------------------------------------------ *
* parametros de entrada:                                             *
* ------------------------------------------------------------------ *
* m          -> arranjo da menoria principal                         *
* ld         -> cargas                                               * 
* x          -> coordenadas                                          *
* el         -> conectividade                                        *
* mat        -> materias                                             *
* nen        -> conectividades por elemento                          *
* mat        -> material por elemento                                *
* typeGeom   -> tipo geometrico do elemento                          *
* typeCal    -> tipo calculo do elemento                             *
* faceRd1    -> condicao de contorno D1                              *
* faceRt1    -> condicao de contorno T1                              *
* faceRfluid -> condicao de contorno fluido                          *
* faceRtemp  -> condicao de contorno temp                            *
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
void wGeoFaceVtk2(Memoria *m         , Loads *ld
                , DOUBLE *x
                , INT *el            , short *nen
                , short *typeGeom    , short *faceRd      
                , INT const nnode    , INT const numel
                , short const ndm    , short const maxViz 
                , short const ndf    , short const maxNo
                , char *nameOut      , bool iws
                , FILE *f)
{
  char head[] = { "FACE_VOLUME_FINITO" };
  int nFace = 0;
  int *face = NULL, *idFace = NULL;
  int *lfaceL = NULL, *lfaceLty = NULL; 
  short *typeGeomFace = NULL, *nenFace = NULL;
  int i;
  int *aux = NULL;
  double ddum;

  HccaAlloc(INT, m, face, numel*MAX_NUM_FACE*MAX_NUM_NODE_FACE
           , "lFace", _AD_);
  HccaAlloc(INT, m, idFace, numel*MAX_NUM_FACE, "iDFace", _AD_);
  HccaAlloc(short, m, typeGeomFace, numel*MAX_NUM_FACE, "ltGface", _AD_);
  HccaAlloc(short, m, nenFace, numel*MAX_NUM_FACE, "lnenFace", _AD_);

/*... ndf*/
  HccaAlloc(int, m, lfaceL, numel*MAX_NUM_FACE, "lfaceS", _AD_);
/*... type*/
  HccaAlloc(int, m,lfaceLty, numel*MAX_NUM_FACE, "lfaceLty", _AD_);

/*...*/
  makeFace( el          , ld
          , faceRd      , typeGeom
          , face        , lfaceL
          , lfaceLty    , idFace
          , typeGeomFace, nenFace
          , maxViz      , maxNo
          , ndf         , numel      
          , &nFace      );
/*...................................................................*/

/*...*/
  if (nFace)
  {
    if (iws)
      f = openFile(nameOut, "wb");
    else
      f = openFile(nameOut, "w");
/*...................................................................*/

/* ...*/
    headVtk(head, iws, f);
/* ..................................................................*/

/*... coordenadas*/
    writeVtkCoor(x, nnode, ndm, iws, f);
/*...................................................................*/

/*... faces*/
    writeVtkCell(face, nenFace, typeGeomFace, nFace
                , MAX_NUM_NODE_FACE, iws, f);
/*...................................................................*/

/*... campo por elemento*/
    fprintf(f, "CELL_DATA %ld\n", (long)nFace);
/*...................................................................*/

/*... relacao face celula*/
    writeVtkProp(idFace, &ddum, nFace, 1, "idCellFace", iws
               , INTEGER_VTK, 1, f);
/*...................................................................*/

/*... valores das cargas por celula*/
    writeVtkProp(lfaceL  , &ddum, nFace, 1
              , "lFaceL", iws, INTEGER_VTK, 1, f);
/*...................................................................*/

/*... valores do tipo de cargas por celula*/
    writeVtkProp(lfaceLty , &ddum, nFace, 1
              , "lfaceLty", iws, INTEGER_VTK, 1, f);
/*...................................................................*/

/*.... campo por no*/
    fprintf(f, "POINT_DATA %ld\n", (long)nnode);
/*...................................................................*/

/*... numero do no*/
    HccaAlloc(int, m, aux, nnode, "el", _AD_);
    if (aux == NULL) {
      fprintf(stderr, "Erro na alocação de lel.\n"
        "Nome do arquivo: %s.\n"
        , __FILE__);
      exit(EXIT_FAILURE);
    }
    for (i = 0; i<nnode; i++)
      aux[i] = i + 1;

    writeVtkProp(aux, &ddum, nnode, 1, "pNode", iws, INTEGER_VTK, 1, f);
    HccaDealloc(m, aux, "el", _AD_);
/*...................................................................*/

    fclose(f);
  }
/*...................................................................*/

/*... dealloc*/
  HccaDealloc(m, lfaceL, "lfaceLty", _AD_);
  HccaDealloc(m, lfaceL, "lfaceS", _AD_);
  HccaDealloc(m, nenFace, "lnenFace", _AD_);
  HccaDealloc(m, typeGeomFace, "ltGface", _AD_);
  HccaDealloc(m, idFace, "iDFace", _AD_);
  HccaDealloc(m, face, "lFace", _AD_);
/*...................................................................*/

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
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER_VTK,1,f);
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
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER_VTK,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... escrever resuldos por celula*/  
  if(ndfD[0] > 0){
    writeVtkProp(&idum,elU,numel,ndfD[0],"elD1",iws,DOUBLE_VTK,1,f);
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
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER_VTK,1,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... escrever resuldos por nos*/  
  if(ndfD[0] > 0){
    writeVtkProp(&idum,nU,nnode,ndfD[0],"nD1",iws,DOUBLE_VTK,1,f);
  }
/*...................................................................*/

  fclose(f);
}
/*********************************************************************/

/********************************************************************** 
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 13/05/2018                                    * 
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
 * opt     -> opcoes de arquivo                                       *
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
               ,DOUBLE *elDensity ,DOUBLE *nDensity
               ,DOUBLE *elCoefDiff,DOUBLE *nCoefDiff
               ,INT nnode         ,INT numel    
               ,short ndm         ,short maxNo 
               ,short numat       ,short ndf   
               ,char *uResEl      ,char *uResNo 
               ,char *gradResEl   ,char *gradResNo 
               ,char *nameOut     ,FileOpt *opt
               ,Temporal *ddt     ,FILE *f)
{
  bool iws = opt->bVtk;
  char str[50];
  short j;
  int    *lel=NULL;
  INT i;
  DOUBLE *p=NULL;

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
  if(ddt->flag)
    timeVtk(ddt->t,ddt->timeStep,iws,f);
/* ..................................................................*/

/*... coordenadas*/
  writeVtkCoor(x,nnode,ndm,iws,f);
/*...................................................................*/
  
/*... conectividades*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  ERRO_MALLOC(lel, "el", __LINE__, __FILE__, __func__);

  for(i=0;i<numel;i++)
    for(j=0;j<maxNo;j++)
      MAT2D(i,j,lel,maxNo) = MAT2D(i,j,el,maxNo)-1;
  
  writeVtkCell(lel,nen,typeGeom,numel,maxNo,iws,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f,"CELL_DATA %ld\n",(long) numel);
/*...................................................................*/

/*... material*/
  HccaAlloc(int,m,lel,numel,"el",_AD_);
  ERRO_MALLOC(lel, "el", __LINE__, __FILE__, __func__);
  for(i=0;i<numel;i++)
    lel[i]=(int) mat[i];
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER_VTK, SCALARS_VTK,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... numero do elemento*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  ERRO_MALLOC(lel, "el", __LINE__, __FILE__, __func__);
  for(i=0;i<numel;i++)
    lel[i]= i+1;
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws,INTEGER_VTK
             , SCALARS_VTK,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  
/*... escrever resultados por celula*/  
  if(opt->fCell && opt->uD1)
    writeVtkProp(&idum,elU,numel,ndf,uResEl,iws,DOUBLE_VTK
                , SCALARS_VTK,f);
/*...................................................................*/

/*... escrever os gradiente por nos*/
  if (opt->fCell && opt->graduD1)
    writeVtkProp(&idum, elGradU, numel, ndm, gradResNo, iws, DOUBLE_VTK
               , VECTORS_VTK, f);
/*...................................................................*/


/*... escrever gradiente por celula*/  
  if (opt->fCell && opt->densityD1)
  {
    strcpy(str, "eDensity");
    writeVtkProp(&idum, elDensity, numel, 1, str, iws
               , DOUBLE_VTK, SCALARS_VTK, f);
  }  
/*...................................................................*/

/*... escrever gradiente por celula*/
  if (opt->fCell && opt->coefDiffD1)
  {
    strcpy(str, "eCeofDiff");
    writeVtkProp(&idum, elCoefDiff, numel, 1, str, iws
               , DOUBLE_VTK, SCALARS_VTK, f);
  }
/*...................................................................*/

/*.... campo por no*/
  fprintf(f,"POINT_DATA %ld\n",(long) nnode);
/*...................................................................*/

/*... numero do no*/
  HccaAlloc(int,m,lel,nnode,"el",_AD_);
  ERRO_MALLOC(lel, "el", __LINE__, __FILE__, __func__);
  for(i=0;i<nnode;i++)
    lel[i]=i+1;
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws,INTEGER_VTK
              , SCALARS_VTK,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... escrever resuldos por nos*/  
  if(opt->fNode && opt->uD1)
    writeVtkProp(&idum,nU    ,nnode,ndf,uResNo,iws,DOUBLE_VTK
               , SCALARS_VTK,f);
/*...................................................................*/
  
/*... escrever os gradiente por nos*/
  if(opt->fNode && opt->graduD1)
    writeVtkProp(&idum,nGradU,nnode,ndm,gradResNo,iws,DOUBLE_VTK
                ,VECTORS_VTK,f);
/*...................................................................*/

/*... escrever gradiente por celula*/
  if (opt->fNode && opt->densityD1)
  {
    strcpy(str, "nDensity");
    writeVtkProp(&idum, nDensity, nnode, 1, str, iws
                , DOUBLE_VTK, SCALARS_VTK, f);
  }
/*...................................................................*/

/*... escrever gradiente por celula*/
  if (opt->fNode && opt->coefDiffD1)
  {
    strcpy(str, "nCeofDiff");
    writeVtkProp(&idum, nCoefDiff, nnode, 1, str, iws
      , DOUBLE_VTK, SCALARS_VTK, f);
  }
/*...................................................................*/
  fclose(f);
}
/*********************************************************************/

/**********************************************************************
* Data de criacao    : 00/00/0000                                    *
* Data de modificaco : 13/05/2018                                    *
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
* elVel   -> campo de velocidades por elementos                      *
* nVel    -> campo de velocidades por nos                            *
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
* opt     -> opcoes de arquivo                                       *
* f       -> arquivlo                                                *
* ------------------------------------------------------------------ *
* parametros de saida  :                                             *
* ------------------------------------------------------------------ *
* -------------------------------------------------------------------*
* OBS:                                                               *
* -------------------------------------------------------------------*
**********************************************************************/
void wResVtkTrans(Memoria *m        , double *x
                , INT *el           , short *mat
                , short *nen        , short *typeGeom
                , DOUBLE *elU       , DOUBLE *noU
                , DOUBLE *elGradU   , DOUBLE *noGradU
                , DOUBLE *elVel     , DOUBLE *noVel
                , DOUBLE *elDensity , DOUBLE *noDensity
                , DOUBLE *elCoefDiff, DOUBLE *noCoefDiff
                , INT nnode         , INT numel
                , short ndm         , short maxNo
                , short numat       , short ndf
                , char **ps
                , char *nameOut     , FileOpt *opt
                , Temporal *ddt     , FILE *f)
{
  bool iws = opt->bVtk;
  short j;
  int    *lel = NULL;
  INT i;
  DOUBLE *p = NULL;
  enum { eT1     = 0, nT1      = 1 
       , eGradT1 = 2, nGradT1  = 3   
       , eVel    = 4, nVel     = 5
       , eDen    = 6, nDen     = 7
       , eCoefDif= 8, nCoefDiff= 9};    

  char head[] = { "DIF_VOLUME_FINITO" };
  double ddum;
  int    idum;

  if (iws)
    f = openFile(nameOut, "wb");
  else
    f = openFile(nameOut, "w");

/* ...*/
  headVtk(head, iws, f);
/* ..................................................................*/

/* ...*/
  if (ddt->flag)
    timeVtk(ddt->t, ddt->timeStep, iws, f);
/* ..................................................................*/

/*... coordenadas*/
  writeVtkCoor(x, nnode, ndm, iws, f);
/*...................................................................*/

/*... conectividades*/
  HccaAlloc(int, m, lel, numel*maxNo, "el", _AD_);
  ERRO_MALLOC(lel, "el", __LINE__, __FILE__, __func__);

  for (i = 0; i<numel; i++)
    for (j = 0; j<maxNo; j++)
      MAT2D(i, j, lel, maxNo) = MAT2D(i, j, el, maxNo) - 1;

  writeVtkCell(lel, nen, typeGeom, numel, maxNo, iws, f);
  HccaDealloc(m, lel, "el", _AD_);
/*...................................................................*/

/*... campo por elemento*/
  fprintf(f, "CELL_DATA %ld\n", (long)numel);
/*...................................................................*/

/*... material*/
  HccaAlloc(int, m, lel, numel, "el", _AD_);
  ERRO_MALLOC(lel, "el", __LINE__, __FILE__, __func__);
  for (i = 0; i<numel; i++)
    lel[i] = (int)mat[i];

  writeVtkProp(lel, &ddum, numel, 1, "mat", iws, INTEGER_VTK, SCALARS_VTK, f);
  HccaDealloc(m, lel, "el", _AD_);
/*...................................................................*/

/*... numero do elemento*/
  HccaAlloc(int, m, lel, numel*maxNo, "el", _AD_);
  ERRO_MALLOC(lel, "el", __LINE__, __FILE__, __func__);
  for (i = 0; i<numel; i++)
    lel[i] = i + 1;

  writeVtkProp(lel, &ddum, numel, 1, "elGlobal", iws, INTEGER_VTK
    , SCALARS_VTK, f);
  HccaDealloc(m, lel, "el", _AD_);
/*...................................................................*/

/*... escrever resultados por celula*/
  if (opt->fCell && opt->uT1)
    writeVtkProp(&idum, elU, numel, ndf, ps[eT1], iws, DOUBLE_VTK
               , SCALARS_VTK, f);
/*...................................................................*/

/*... escrever os gradiente por nos*/
  if (opt->fCell && opt->graduT1)
    writeVtkProp(&idum, elGradU   , numel      , ndm, ps[eGradT1]
                , iws , DOUBLE_VTK, VECTORS_VTK, f);
/*...................................................................*/


/*... escrever gradiente por celula*/
  if (opt->fCell && opt->densityT1)
  {
    writeVtkProp(&idum      , elDensity, numel, 1, ps[eDen], iws
                , DOUBLE_VTK, SCALARS_VTK, f);
  }
/*...................................................................*/

/*... escrever gradiente por celula*/
  if (opt->fCell && opt->coefDiffT1)
  {
    writeVtkProp(&idum, elCoefDiff, numel, 1, ps[eCoefDif], iws
      , DOUBLE_VTK, SCALARS_VTK, f);
  }
/*...................................................................*/

/*... escrever gradiente por celula*/
  if (opt->fCell && opt->vel)
  {
    writeVtkProp(&idum, elVel, numel, ndm, ps[eVel], iws
      , DOUBLE_VTK, VECTORS_VTK, f);
  }
/*...................................................................*/

/*.... campo por no*/
  fprintf(f, "POINT_DATA %ld\n", (long)nnode);
/*...................................................................*/

/*... numero do no*/
  HccaAlloc(int, m, lel, nnode, "el", _AD_);
  ERRO_MALLOC(lel, "el", __LINE__, __FILE__, __func__);
  for (i = 0; i<nnode; i++)
    lel[i] = i + 1;

  writeVtkProp(lel, &ddum, nnode, 1, "pNode", iws, INTEGER_VTK
    , SCALARS_VTK, f);
  HccaDealloc(m, lel, "el", _AD_);
/*...................................................................*/

/*... escrever resuldos por nos*/
  if (opt->fNode && opt->uT1)
    writeVtkProp(&idum, noU, nnode, ndf, ps[nT1], iws, DOUBLE_VTK
      , SCALARS_VTK, f);
/*...................................................................*/

/*... escrever os gradiente por nos*/
  if (opt->fNode && opt->graduT1)
    writeVtkProp(&idum, noGradU, nnode, ndm, ps[nGradT1], iws, DOUBLE_VTK
               , VECTORS_VTK, f);
/*...................................................................*/

/*... escrever gradiente por celula*/
  if (opt->fNode && opt->densityT1)
  {
    writeVtkProp(&idum, noDensity, nnode, 1, ps[nDen], iws
      , DOUBLE_VTK, SCALARS_VTK, f);
  }
/*...................................................................*/

/*... escrever gradiente por celula*/
  if (opt->fNode && opt->coefDiffT1)
  {
    writeVtkProp(&idum, noCoefDiff, nnode, 1, ps[nCoefDiff], iws
      , DOUBLE_VTK, SCALARS_VTK, f);
  }
/*...................................................................*/

/*... escrever gradiente por celula*/
  if (opt->fNode && opt->vel)
  {
    writeVtkProp(&idum, noVel, nnode, ndm, ps[nVel], iws
      , DOUBLE_VTK, VECTORS_VTK, f);
  }
/*...................................................................*/

  fclose(f);
}
/*********************************************************************/

/********************************************************************** 
 * Data de criacao    : 30/06/2016                                    *
 * Data de modificaco : 20/10/2019                                    * 
 *------------------------------------------------------------------- * 
 * WRESVTKFLUID:escreve a malha com os resultados para problemas de   *  
 * de escomentos de fluidos imcompressivel                            *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * m            -> arranjo da menoria principal                       *  
 * x            -> coordenadas                                        * 
 * cc           -> centro geomentrico das celulas                     * 
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
 * elEddyVis    -> viscosidade turbulenta                             *
 * eDensityFluid-> densidade do fluido (cell)                         *
 * nDensityFluid-> densidade do fluido (node)                         *
 * eDyViscosity -> viscosidade molecular (cell)                       *
 * nDyViscosity   > viscosidade molecular (cell)                      *
 * eCd          -> coeficientes dincamicamente calculados (cell)      *
 * nCd          -> coeficientes dincamicamente calculados (node)      *
 * eWallPar   -> parametros de parede  ( yPlus, uPlus, uFri) (cell)   *
 * nWallPar   -> parametros de parede  ( yPlus, uPlus, uFri) (node)   *
 * tConductivity-> condutividade termica                              *
 * cDiffSp      -> coeficiente de difusao das especies
 * nel          -> numeracao do elemento                              *
 * nnode        -> numero de nos                                      *  
 * numel        -> numero de elementos                                *
 * ndm          -> numero de dimensao                                 *
 * maxNo        -> numero maximo de nos por elemento                  *
 * numat        -> numero de materias                                 *
 * ndf          -> graus de liberdade das equacoes                    *
 * ntn          -> numero de termos no tensor ( 4 ; 6)                *
 * nameOut      -> nome de arquivo de saida                           *
 * opt          -> opcoes do arquivo                                  *
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
void wResVtkFluid(Memoria *m     , DOUBLE *x 
          , DOUBLE *cc     
          , INT *el              , short *mat    
          , short *nen           , short *typeGeom
          , DOUBLE *elPres       , DOUBLE *nPres
          , DOUBLE *elGradPres   , DOUBLE *nGradPres
          , DOUBLE *elVel        , DOUBLE *nVel      
          , DOUBLE *elGradVel    , DOUBLE *nGradVel 
          , DOUBLE *elTemp       , DOUBLE *nTemp
          , DOUBLE *elGradTemp   , DOUBLE *nGradTemp  
          , DOUBLE *elEddyVis    , DOUBLE *nEddyVis
          , DOUBLE *eDensityFluid, DOUBLE *nDensityFluid
          , DOUBLE *eDyViscosity , DOUBLE *nDyViscosity
          , DOUBLE *eStressR     , DOUBLE *nStressR
          , DOUBLE *eCd          , DOUBLE *nCd
          , DOUBLE *eWallPar     , DOUBLE *nWallPar
          , DOUBLE *eKturb       , DOUBLE *nKturb
          , DOUBLE *eMedVel      , DOUBLE *nMedVel
          , DOUBLE *eSheat       , DOUBLE *nSheat
          , DOUBLE *eTCond       , DOUBLE *nTCond
          , DOUBLE *eGradRho     , DOUBLE *nGradRho
          , INT nnode            , INT numel    
          , short const ndm      , short const maxNo 
          , short const numat    , short const ndf
          , short const ntn        
          , char *nameOut        , FileOpt *opt
          , bool fKelvin         , Mean *media  
          , Temporal *ddt         , FILE *f)
{
  bool iws = opt->bVtk;
  char str[50];
  int    *lel=NULL;
  DOUBLE *p=NULL,*w=NULL;
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
  if(ddt->flag)
    timeVtk(ddt->t,ddt->timeStep,iws,f);
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
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER_VTK,SCALARS_VTK,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... numero do elemento*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  ERRO_MALLOC(lel,"el",__LINE__,__FILE__,__func__)
  for(i=0;i<numel;i++)
    lel[i]= i+1;
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws
              ,INTEGER_VTK,SCALARS_VTK,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  
/*...*/
  if(opt->cc){
    strcpy(str,"cc");
    writeVtkProp(&idum,cc    ,numel,ndm,str,iws
                 ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*...*/
  if(media->fVel && opt->fCell){
    strcpy(str,"<eVel>");
    writeVtkProp(&idum,eMedVel,numel,ndm,str,iws
                 ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*... escrever resultados de pressao por celula*/   
  if(opt->pres && opt->fCell){
    strcpy(str,"ePres");
    writeVtkProp(&idum,elPres,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escrever gradiente da pressao por celula*/  
  if(opt->gradPres && opt->fCell){
    strcpy(str,"eGradPres");
    writeVtkProp(&idum,elGradPres,numel,ndm,str ,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/


/*... escrever campo de velociade por celula*/  
  if(opt->vel && opt->fCell ){
    strcpy(str,"eVel");
    writeVtkProp(&idum,elVel,numel,ndm,str,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
    
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__)
    strcpy(str,"modCellVel");
    makeModuleVel(p,elVel,numel,ndm);
    writeVtkProp(&idum,p,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escrever gradiente de velocidade por celula*/  
  if(opt->gradVel && opt->fCell){  
    strcpy(str,"eGradVel");
    if( ndm == 2) 
      writeVtkProp(&idum,elGradVel,numel,2*ndm,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
    else if( ndm == 3 )
      writeVtkProp(&idum,elGradVel,numel,3*ndm,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve o campo de temperatura  por celula*/
  if (opt->temp && opt->fCell){
    strcpy(str,"eTemp");
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    alphaProdVector(1.e0,elTemp,numel,p);
    if(fKelvin){
      if(!opt->pKelvin) convTempForKelvin(p, numel,false);
      writeVtkProp(&idum,p,numel,1,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
    }
    else{
      if(opt->pKelvin) convTempForKelvin(p, numel,true);
      writeVtkProp(&idum,p,numel,1,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
    }
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a gradiente de temperatura por celula*/
  if (opt->gradTemp && opt->fCell) {
    strcpy(str,"eGradTemp");
    writeVtkProp(&idum,elGradTemp,numel,ndm,str,iws
               ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*... escreve a viscosidade turbulenta por celula*/  
  if(opt->eddyViscosity && opt->fCell){
    strcpy(str,"eEddyViscosity");
    writeVtkProp(&idum,elEddyVis,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a viscosidade dinamica por celula*/  
  if(opt->dViscosity && opt->fCell){
    strcpy(str,"eDinamicyViscosity");
    writeVtkProp(&idum,eDyViscosity,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a condutividade termica por celula*/  
  if(opt->tConductivity && opt->fCell ){
    strcpy(str,"eThermoCondutivity");
    writeVtkProp(&idum,eTCond,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve o calor especifico por celula*/  
  if(opt->specificHeat && opt->fCell ){
    strcpy(str,"eSpecificHeat");
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    writeVtkProp(&idum,eSheat,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a massa especifica por celula*/  
  if(opt->densityFluid && opt->fCell ){
    strcpy(str,"eRho");
    writeVtkProp(&idum,eDensityFluid,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a vorticidade */  
  if(opt->vorticity && opt->fCell ){
    strcpy(str,"eVorticity");
    HccaAlloc(DOUBLE,m,p,numel*3,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeVorticity(p,elGradVel,numel,ndm);
    writeVtkProp(&idum,p,numel,ndm,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a yPlus */  
  if(opt->wallParameters && opt->fCell){
    strcpy(str,"eWallParameters(y+|u+|uf|sW)");
    writeVtkProp(&idum,eWallPar,numel,4,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a energia cinetica */  
  if(opt->kinetic && opt->fCell ){
    strcpy(str,"eKinetic");
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeKineticEnergy(p,elVel,eDensityFluid,numel,ndm);
    writeVtkProp(&idum,p,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve o tensor residual */  
  if(opt->stressR && opt->fCell ){
    HccaAlloc(DOUBLE,m,p,numel*ntn,"p",_AD_);
    HccaAlloc(DOUBLE,m,w,numel*ntn,"w",_AD_);
/*... estrutural*/
    strcpy(str,"eStressRs");
    writeVtkProp(&idum,eStressR,numel,ntn,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    makeStress(p,elGradVel,eDyViscosity,numel,ndm,ntn,false);
/*... funcional*/
    strcpy(str,"eStressRf");
    writeVtkProp(&idum,p,numel,ntn,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
/*... total*/
    strcpy(str,"eStressRtotal");
    addVector(1.e0     , eStressR
             ,1.e0     , p
             ,numel*ntn, w);
    writeVtkProp(&idum,w,numel,ntn,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,w,"w",_AD_);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... coeficiente dinamicamente calculados */  
  if(opt->cDynamic && opt->fCell){
    strcpy(str,"eCdyn");
    writeVtkProp(&idum,eCd,numel,2,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a energia cinetica */  
  if(opt->Qcriterion && opt->fCell ){
    strcpy(str,"eQCriterion");
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeQcriterion(p,elGradVel,numel,ndm);
    writeVtkProp(&idum,p,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a pressao total */  
  if(opt->presTotal &&  opt->fCell ){
    strcpy(str,"ePresTotal");
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makePresTotal(p,elPres,elVel,eDensityFluid,numel,ndm);
    writeVtkProp(&idum,p,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a energia cinetica turbulenta */  
  if(opt->kTurb &&  opt->fCell ){
    strcpy(str,"eKTurbl");
    writeVtkProp(&idum,eKturb,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escrever gradiente da pressao por celula*/  
  if(opt->gradRho && opt->fCell){
    strcpy(str,"eGradRho");
    writeVtkProp(&idum,eGradRho,numel,ndm,str ,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*.... campo por no*/
  fprintf(f,"POINT_DATA %ld\n",(long) nnode);
/*...................................................................*/

/*... numero do no*/
  HccaAlloc(int,m,lel,nnode,"el",_AD_);
  ERRO_MALLOC(lel,"el",__LINE__,__FILE__,__func__)
  for(i=0;i<nnode;i++)
    lel[i]=i+1;
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws
              ,INTEGER_VTK,SCALARS_VTK,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*...*/
  if(media->fVel && opt->fNode){
    strcpy(str,"<nVel>");
    writeVtkProp(&idum,nMedVel,nnode,ndm,str,iws
                 ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*... escrever resultados de pressao por nos*/  
  if(opt->pres && opt->fNode){
    strcpy(str,"NodePres");
    writeVtkProp(&idum,nPres ,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/
  
/*... escrever os gradiente de pressao por nos*/  
  if(opt->gradPres && opt->fNode){
    strcpy(str,"NodeGradPres");
    writeVtkProp(&idum,nGradPres,nnode,ndm,str,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*... escrever as velocidade por nos*/  
  if(opt->vel && opt->fNode){
    strcpy(str,"NodeVel");
    writeVtkProp(&idum,nVel,nnode,ndm,str,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__)
    strcpy(str,"modNodeVel");
    makeModuleVel(p,nVel,nnode,ndm);
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escrever as gradiente de velocidade por nos*/ 
  if(opt->gradVel && opt->fNode){  
    strcpy(str,"NodeGradVel");
    if( ndm == 2) 
      writeVtkProp(&idum,nGradVel,nnode,2*ndm,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
    else if( ndm == 3) 
      writeVtkProp(&idum,nGradVel,nnode,3*ndm,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escrever resultados de energia por nos*/
  if (opt->temp && opt->fNode){
    strcpy(str,"NodeTemp");
    HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    alphaProdVector(1.e0,nTemp,nnode,p);
    if(fKelvin){
      if(!opt->pKelvin) convTempForKelvin(p, nnode,false);
      writeVtkProp(&idum, p, nnode, 1, str, iws
                  , DOUBLE_VTK, SCALARS_VTK, f);
    }
    else{
      if(opt->pKelvin) convTempForKelvin(p, nnode,true);      
      writeVtkProp(&idum, p, nnode, 1, str, iws
                  , DOUBLE_VTK, SCALARS_VTK, f);     
    }
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escrever a viscosidade turbulenta por celula*/  
  if(opt->eddyViscosity && opt->fNode){
    strcpy(str,"NodeEddyViscosity");
    writeVtkProp(&idum, nEddyVis, nnode, 1, str, iws
                , DOUBLE_VTK, SCALARS_VTK, f);
  }
/*...................................................................*/

/*... escrever os gradiente de pressao por nos*/
  if (opt->gradTemp && opt->fNode){ 
    strcpy(str,"NodeGradTemp");
    writeVtkProp(&idum, nGradTemp, nnode, ndm, str,iws
               , DOUBLE_VTK, VECTORS_VTK, f);
  }
/*...................................................................*/

/*... escreve a vorticidade por no*/  
  if(opt->vorticity && opt->fNode ){
    strcpy(str,"NodeVorticity");
    HccaAlloc(DOUBLE,m,p,nnode*3,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeVorticity(p,nGradVel,nnode,ndm);
    writeVtkProp(&idum,p,nnode,ndm,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a tensor desviador por no*/ 
  if(opt->stress && opt->fNode ){
/*  strcpy(str,"NodeStress");
    HccaAlloc(DOUBLE,m,p,nnode*ntn,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeStress(p,nGradVel,nDyViscosity,nnode,ndm,true);
    writeVtkProp(&idum,p,nnode,ntn,str,iws
                ,DOUBLE_VTK,TENSORS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);*/
  }
/*...................................................................*/

/*... escreve a energia cinetica por no*/  
  if(opt->kinetic && opt->fNode ){
    strcpy(str,"NodeKinetic");
    HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeKineticEnergy(p,nVel,nDensityFluid,nnode,ndm);
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve o tensor residual por no*/  
  if(opt->stressR && opt->fNode ){
    strcpy(str,"nStressRs");
    writeVtkProp(&idum,nStressR,nnode,6,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaAlloc(DOUBLE,m,p,numel*ntn,"p",_AD_);
    makeStress(p,nGradVel,nDyViscosity,nnode,ndm,ntn,false);
    strcpy(str,"nStressRf");
    writeVtkProp(&idum,p,nnode,ntn,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... coeficiente dinmicamente calculados por no*/  
  if(opt->cDynamic && opt->fNode){
    strcpy(str,"nCdyn");
    writeVtkProp(&idum,nCd,nnode,2,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a yPlus por no*/   
  if(opt->wallParameters && opt->fNode){
    strcpy(str,"nWallParameters(y+|u+|uf|sW)");
    writeVtkProp(&idum,nWallPar,nnode,4,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a energia cinetica por no*/   
  if(opt->Qcriterion &&  opt->fNode ){
    strcpy(str,"nQCriterion");
    HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeQcriterion(p,nGradVel,nnode,ndm);
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a pressao total por no*/  
  if(opt->presTotal &&  opt->fNode ){
    strcpy(str,"nPresTotal");
    HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makePresTotal(p,nPres,nVel,nDensityFluid,nnode,ndm);
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a energia turbuleta por no*/  
  if(opt->kTurb &&  opt->fNode ){
    strcpy(str,"nKturbl");
    writeVtkProp(&idum,nKturb,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a condutividae termica por no*/  
  if(opt->tConductivity && opt->fNode ){
    strcpy(str,"nThermoCondutivity");
    writeVtkProp(&idum,nTCond,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve o calor especifico por no*/  
  if(opt->specificHeat && opt->fNode ){
    strcpy(str,"nSpecificHeat");
    writeVtkProp(&idum,nSheat,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*...  escreve a massa especifica por no*/  
  if(opt->densityFluid && opt->fNode ){
    strcpy(str,"nDensityFluid");
    writeVtkProp(&idum,nDensityFluid,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escrever gradiente da densidade por celula*/  
  if(opt->gradRho && opt->fNode){
    strcpy(str,"nGradRho");
    writeVtkProp(&idum,nGradRho,nnode,ndm,str ,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

  fclose(f);
}
/*********************************************************************/

/********************************************************************** 
 * Data de criacao    : 05/08/2018                                    *
 * Data de modificaco : 09/11/2019                                    * 
 *------------------------------------------------------------------- * 
 * wResVtkCombustion :                                                * 
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * m            -> arranjo da menoria principal                       *  
 * x            -> coordenadas                                        * 
 * cc           -> centro geomentrico das celulas                     * 
 * el           -> conectividade                                      * 
 * mat          -> materias                                           *  
 * nen          -> conectividades por elemento                        *  
 * mat          -> material por elemento                              *
 * typeGeom     -> tipo geometrico do elemento                        *
 * elPres       -> pressao (cell)                                     *
 * nPres        -> pressao (node)                                     *
 * elGradPres   -> gradientes da pressao (cell)                       *
 * nGradPres    -> gradientes da pressao (node)                       *
 * elVel        -> velocidade (cell)                                  *
 * nVel         -> velocidade (node)                                  *
 * elGradVel    -> gradiente das velocidades (cell)                   *
 * nGradVel     -> gradiente das velocidades (node)                   *
 * elTemp       -> temperatura (cell)                                 *
 * nTemp        -> temperatura (node)                                 *
 * elGradEnergy -> gradientes da energia(Temp) (cell)                 *
 * nGradEnergy  -> gradientes da energia(Temp) (node)                 *
 * elZomb       -> fracao massica das especies agrupadas (cell)       *
 * nZomb        -> fracao massica das especies agrupadas (node)       *
 * elGradZomb   -> campo da grad fra massica das especies agraupadas  *
 * nGradZomb    -> campo da grad fra massica das especies agraupadas  *
 * elEddyVis    -> viscosidade turbulenta (cell)                      *
 * nEddyVis     -> viscosidade turbulenta (node)                      *
 * eDensityFluid-> densidade do fluido (cell)                         *
 * nDensityFluid-> densidade do fluido (node)                         *
 * eDyViscosity -> viscosidade molecular (cell)                       *
 * nDyViscosity   > viscosidade molecular (cell)                      *
 * eCd          -> coeficientes dincamicamente calculados (cell)      *
 * nCd          -> coeficientes dincamicamente calculados (node)      *
 * eWallPar   -> parametros de parede  ( yPlus, uPlus, uFri) (cell)   *
 * nWallPar   -> parametros de parede  ( yPlus, uPlus, uFri) (node)   *
 * eKturb       -> energia cinetrica  turbulenta (Cell)               *
 * nKturb       -> energia cinetrica  turbulenta (node)               *
 * eWk          -> taxa de reacao (Cell)                              *
 * nWk          -> taxa de reacao (Node)                              *
 * eYfrac       -> fracao massica das especies primitivas (cell)      *
 * nYfrac       -> fracao massica das especies primitivas (node)      *
 * elGradY      -> campo da grad fra massica das especies primitivas  *
 * nGradY       -> campo da grad fra massica das especies primitivas  *
 * eMedVel      -> media das velocidades (cell)                       *
 * nMedVel      -> media das velocidades (node)                       *
 * eEntkalpyK   -> entalpia por especies (cell)                       *
 * nEntkalpyK   -> entalpia por especies (node)                       *
 * specificHeat -> massa especifica (cell)                            *
 * tConductivity-> condutividade termica (cell)                       *   
 * nel          -> numeracao do elemento                              *
 * nnode        -> numero de nos                                      *  
 * numel        -> numero de elementos                                *
 * ndm          -> numero de dimensao                                 *
 * maxNo        -> numero maximo de nos por elemento                  *
 * numat        -> numero de materias                                 *
 * ndf          -> graus de liberdade das equacoes                    *
 * ntn          -> numero de termos no tensor ( 4 ; 6)                *
 * nOfPrSp      -> numero de especies primitivas                      *
 * nComb        -> numero de especies resolvidas                      *
 * posN2        -> posisao do N2 nos vetores de especies              *
 * fSpecies     -> especies existentes                                * 
 * nameOut      -> nome de arquivo de saida                           *
 * opt          -> opcoes do arquivo                                  *
 * f            -> arquivo                                            *
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
void wResVtkCombustion(Memoria *m , Combustion *cModel     
          , DOUBLE *x             , DOUBLE *cc     
          , INT *el               , short *mat    
          , short *nen            , short *typeGeom
          , DOUBLE *elPres        , DOUBLE *nPres
          , DOUBLE *elGradPres    , DOUBLE *nGradPres
          , DOUBLE *elVel         , DOUBLE *nVel      
          , DOUBLE *elGradVel     , DOUBLE *nGradVel 
          , DOUBLE *elTemp        , DOUBLE *nTemp
          , DOUBLE *elGradTemp    , DOUBLE *nGradTemp  
          , DOUBLE *elZcomb       , DOUBLE *nZcomb  
          , DOUBLE *elGradZcomb   , DOUBLE *nGradZcomb 
          , DOUBLE *elEddyVis     , DOUBLE *nEddyVis
          , DOUBLE *eDensityFluid , DOUBLE *nDensityFluid
          , DOUBLE *eDyViscosity  , DOUBLE *nDyViscosity
          , DOUBLE *eStressR      , DOUBLE *nStressR
          , DOUBLE *eCd           , DOUBLE *nCd
          , DOUBLE *eWallPar      , DOUBLE *nWallPar
          , DOUBLE *eKturb        , DOUBLE *nKturb
          , DOUBLE *eWk           , DOUBLE *nWk      
          , DOUBLE *eYfrac        , DOUBLE *nYfrac
          , DOUBLE *eGradY        , DOUBLE *nGradY 
          , DOUBLE *eHeatRe       , DOUBLE *nHeatRe     
          , DOUBLE *eMedVel       , DOUBLE *nMedVel
          , DOUBLE *eEnthalpyK    , DOUBLE *nEnthalpyK   
          , DOUBLE *eSheat        , DOUBLE *nSheat
          , DOUBLE *eTCond        , DOUBLE *nTCond
          , DOUBLE *eDiffSp       , DOUBLE *nDiffSp  
          , DOUBLE *eGradRho      , DOUBLE *nGradRho
          , DOUBLE *eMolar        , DOUBLE *nMolar
          , DOUBLE *eTreactor  
          , INT nnode             , INT numel    
          , short const ndm       , short const maxNo 
          , short const numat     , short const ndf
          , short const ntn       
          , char *nameOut         , FileOpt *opt
          , bool fKelvin          , Mean *media  
          , DOUBLE const ts       , FILE *f)
{
  bool iws = opt->bVtk;
  char str[50], st[MAX_STR_NUMBER];;
  int    *lel=NULL;
  DOUBLE *p=NULL,*w=NULL;
  INT i;
  short j,
        nSp   = cModel->nOfSpecies;
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
  timeVtk2(ts,iws,f);
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
   
  writeVtkProp(lel,&ddum,numel,1,"mat",iws,INTEGER_VTK,SCALARS_VTK,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*... numero do elemento*/
  HccaAlloc(int,m,lel,numel*maxNo,"el",_AD_);
  ERRO_MALLOC(lel,"el",__LINE__,__FILE__,__func__)
  for(i=0;i<numel;i++)
    lel[i]= i+1;
   
  writeVtkProp(lel,&ddum,numel,1,"elGlobal",iws
              ,INTEGER_VTK,SCALARS_VTK,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/
  
/*...*/
  if(opt->cc)
  {
    strcpy(str,"cc");
    writeVtkProp(&idum,cc    ,numel,ndm,str,iws
                 ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*...*/
  if(media->fVel && opt->fCell)
  {
    strcpy(str,"<eVel>");
    writeVtkProp(&idum,eMedVel,numel,ndm,str,iws
                 ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*... escrever resultados de pressao por celula*/   
  if(opt->pres && opt->fCell)
  {
    strcpy(str,"ePres");
    writeVtkProp(&idum,elPres,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escrever gradiente da pressao por celula*/  
  if(opt->gradPres && opt->fCell)
  {
    strcpy(str,"eGradPres");
    writeVtkProp(&idum,elGradPres,numel,ndm,str ,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*... escrever o modulo do campo de velociade por celula*/  
  if(opt->vel && opt->fCell )
  {
    strcpy(str,"eVel");
    writeVtkProp(&idum,elVel,numel,ndm,str,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
    
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__)
    strcpy(str,"modCellVel");
    makeModuleVel(p,elVel,numel,ndm);
    writeVtkProp(&idum,p,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escrever gradiente de velocidade por celula*/  
  if(opt->gradVel && opt->fCell)
  {  
    strcpy(str,"eGradVel");
    if( ndm == 2) 
      writeVtkProp(&idum,elGradVel,numel,2*ndm,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
    else if( ndm == 3 )
      writeVtkProp(&idum,elGradVel,numel,3*ndm,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escrever campo de energia por celula*/
  if (opt->temp && opt->fCell)
  {
    strcpy(str,"eTemp");
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    alphaProdVector(1.e0,elTemp,numel,p);
    if(fKelvin)
    {
      if(!opt->pKelvin) convTempForKelvin(p, numel,false);
      writeVtkProp(&idum,p,numel,1,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
    }
    else
    {
      if(opt->pKelvin) convTempForKelvin(p, numel,true);
      writeVtkProp(&idum,p,numel,1,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
    }
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escrever gradiente de velocidade por celula*/
  if (opt->gradTemp && opt->fCell) 
  {
    strcpy(str,"eGradTemp");
    writeVtkProp(&idum,elGradTemp,numel,ndm,str,iws
               ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*... escreve a viscosidade turbulenta por celula*/  
  if(opt->eddyViscosity && opt->fCell){
    strcpy(str,"eEddyViscosity");
    writeVtkProp(&idum,elEddyVis,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a viscosidade dinamica por celula*/  
  if(opt->dViscosity && opt->fCell){
    strcpy(str,"eDviscosity");
    writeVtkProp(&idum,eDyViscosity,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a condutividade termica por celula*/  
  if(opt->tConductivity && opt->fCell ){
    strcpy(str,"eThermoCondutivity");
    writeVtkProp(&idum,eTCond,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve o calor especifico por celula*/  
  if(opt->specificHeat && opt->fCell ){
    strcpy(str,"eSpecificHeat");
    writeVtkProp(&idum,eSheat,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a massa especifica por celula*/  
  if(opt->densityFluid && opt->fCell ){
    strcpy(str,"eRho");
    writeVtkProp(&idum,eDensityFluid,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a vorticidade */  
  if(opt->vorticity && opt->fCell )
  {
    strcpy(str,"eVorticity");
    HccaAlloc(DOUBLE,m,p,numel*3,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeVorticity(p,elGradVel,numel,ndm);
    writeVtkProp(&idum,p,numel,ndm,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a yPlus */  
  if(opt->wallParameters && opt->fCell)
  {
    strcpy(str,"eWallParameters(y+|u+|uf|sW)");
    writeVtkProp(&idum,eWallPar,numel,4,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a energia cinetica */  
  if(opt->kinetic && opt->fCell )
  {
    strcpy(str,"eKinetic");
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeKineticEnergy(p,elVel,eDensityFluid,numel,ndm);
    writeVtkProp(&idum,p,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve o tensor residual */  
  if(opt->stressR && opt->fCell )
  {
    HccaAlloc(DOUBLE,m,p,numel*ntn,"p",_AD_);
    HccaAlloc(DOUBLE,m,w,numel*ntn,"w",_AD_);
/*... estrutural*/
    strcpy(str,"eStressRs");
    writeVtkProp(&idum,eStressR,numel,ntn,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    makeStress(p,elGradVel,eDyViscosity,numel,ndm,ntn,false);
/*... funcional*/
    strcpy(str,"eStressRf");
    writeVtkProp(&idum,p,numel,ntn,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
/*... total*/
    strcpy(str,"eStressRtotal");
    addVector(1.e0     , eStressR
             ,1.e0     , p
             ,numel*ntn, w);
    writeVtkProp(&idum,w,numel,ntn,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,w,"w",_AD_);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... coeficiente dinmicamente calculados */  
  if(opt->cDynamic && opt->fCell)
  {
    strcpy(str,"eCdyn");
    writeVtkProp(&idum,eCd,numel,2,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a energia cinetica */  
  if(opt->Qcriterion && opt->fCell )
  {
    strcpy(str,"eQCriterion");
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeQcriterion(p,elGradVel,numel,ndm);
    writeVtkProp(&idum,p,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a pressao total */  
  if(opt->presTotal &&  opt->fCell )
  {
    strcpy(str,"ePresTotal");
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makePresTotal(p,elPres,elVel,eDensityFluid,numel,ndm);
    writeVtkProp(&idum,p,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a zComb*/  
  if(opt->zComb &&  opt->fCell )
  {
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);

    for(i=0;i<cModel->nComb;i++)
    {
      getColFromMatrix(p,elZcomb,numel,cModel->nComb,i); 
      iota(i,st);
      strcpy(str,"eZ");
      strcat(str,st);
      writeVtkProp(&idum,p,numel,1,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
    }

    sumFracZ(p,elZcomb,numel,cModel->nComb);
    strcpy(str,"eZTotal");
    writeVtkProp(&idum,p,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);

    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a zComb*/  
  if(opt->gradZcomb &&  opt->fCell )
  {
    strcpy(str,"eGradZcomb");
    writeVtkProp(&idum,elGradZcomb,numel,cModel->nComb*ndm,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a energia cinetica turbulenta */  
  if(opt->kTurb &&  opt->fCell )
  {
    strcpy(str,"eKTurbl");
    writeVtkProp(&idum,eKturb,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a taxa de consumo da especies*/  
  if(opt->wk &&  opt->fCell )
  {
    strcpy(str,"eWk");
    writeVtkProp(&idum,eWk,numel,nSp,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a taxa de liberacao de calor pela reacao quimica */  
  if(opt->rateHeatComb &&  opt->fCell )
  {
    strcpy(str,"eRateHeatComb");
    writeVtkProp(&idum,eHeatRe,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve as fracoes massicas por elmento */  
  if(opt->yFrac &&  opt->fCell )
  {
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
/*... Fuel*/
    for(i=0;i<cModel->chem.nSp;i++)
    {
      strcpy(str,"eY");
      strcat(str,cModel->chem.sp[i].name);
      getColFromMatrix(p,eYfrac,numel,cModel->chem.nSp,i); 
      writeVtkProp(&idum,p,numel,1,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
    }

    sumFracZ(p,eYfrac,numel, cModel->chem.nSp);
    strcpy(str,"eYTotal");
    writeVtkProp(&idum,p,numel,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);

    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... coeficiente de difusao das especies */  
  if(opt->coefDiffSp &&  opt->fCell )
  {
    strcpy(str,"eCoefDiffSp");
    writeVtkProp(&idum,eDiffSp,numel,cModel->nOfSpecies,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... entalpiada individuais por especies */  
  if(opt->enthalpyk &&  opt->fCell )
  {
    strcpy(str,"eEnthalpyK");
    writeVtkProp(&idum,eEnthalpyK,numel,cModel->nOfSpecies,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a gradY*/  
  if(opt->gradY &&  opt->fCell )
  {
    strcpy(str,"eGradY");
    writeVtkProp(&idum,eGradY,numel,cModel->nOfSpecies*ndm,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a gradY*/  
  if(opt->tReactor &&  opt->fCell )
  {
    strcpy(str,"tReactor");
    writeVtkProp(&idum,eTreactor,numel,N_TERMS_REACTOR,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escrever gradiente da pressao por celula*/  
  if(opt->gradRho && opt->fCell){
    strcpy(str,"eGradRho");
    writeVtkProp(&idum,eGradRho,numel,ndm,str ,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*... escrever a massa molar da mistura por celula*/
  if (opt->mMolar && opt->fCell) 
  {
    strcpy(str,"eMolar");
    writeVtkProp(&idum,eMolar,numel,1,str,iws
               ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*.... campo por no*/
  fprintf(f,"POINT_DATA %ld\n",(long) nnode);
/*...................................................................*/

/*... numero do no*/
  HccaAlloc(int,m,lel,nnode,"el",_AD_);
  ERRO_MALLOC(lel,"el",__LINE__,__FILE__,__func__)
  for(i=0;i<nnode;i++)
    lel[i]=i+1;
   
  writeVtkProp(lel,&ddum,nnode,1,"pNode",iws
              ,INTEGER_VTK,SCALARS_VTK,f);
  HccaDealloc(m,lel,"el",_AD_);
/*...................................................................*/

/*...*/
  if(media->fVel && opt->fNode)
  {
    strcpy(str,"<nVel>");
    writeVtkProp(&idum,nMedVel,nnode,ndm,str,iws
                 ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*... escrever resultados de pressao por nos*/  
  if(opt->pres && opt->fNode)
  {
    strcpy(str,"NodePres");
    writeVtkProp(&idum,nPres ,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/
  
/*... escrever os gradiente de pressao por nos*/  
  if(opt->gradPres && opt->fNode)
  {
    strcpy(str,"NodeGradPres");
    writeVtkProp(&idum,nGradPres,nnode,ndm,str,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*... escrever as velocidade por nos*/  
  if(opt->vel && opt->fNode)
  {
    strcpy(str,"NodeVel");
    writeVtkProp(&idum,nVel,nnode,ndm,str,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
    HccaAlloc(DOUBLE,m,p,numel,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__)
    strcpy(str,"modNodeVel");
    makeModuleVel(p,nVel,nnode,ndm);
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escrever as gradiente de velocidade por nos*/ 
  if(opt->gradVel && opt->fNode)
  {  
    strcpy(str,"NodeGradVel");
    if( ndm == 2) 
      writeVtkProp(&idum,nGradVel,nnode,2*ndm,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
    else if( ndm == 3) 
      writeVtkProp(&idum,nGradVel,nnode,3*ndm,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escrever resultados de energia por nos*/
  if (opt->temp && opt->fNode)
  {
    strcpy(str,"NodeTemp");
    HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    alphaProdVector(1.e0,nTemp,nnode,p);
    if(fKelvin)
    {
      if(!opt->pKelvin) convTempForKelvin(p, nnode,false);
      writeVtkProp(&idum, p, nnode, 1, str, iws
                  , DOUBLE_VTK, SCALARS_VTK, f);
    }
    else
    {
      if(opt->pKelvin) convTempForKelvin(p, nnode,true);      
      writeVtkProp(&idum, p, nnode, 1, str, iws
                  , DOUBLE_VTK, SCALARS_VTK, f);     
    }
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escrever a viscosidade turbulenta por celula*/  
  if(opt->eddyViscosity && opt->fNode)
  {
    strcpy(str,"NodeEddyViscosity");
    writeVtkProp(&idum, nEddyVis, nnode, 1, str, iws
                , DOUBLE_VTK, SCALARS_VTK, f);
  }
/*...................................................................*/

/*... escrever os gradiente de pressao por nos*/
  if (opt->gradTemp && opt->fNode)
  { 
    strcpy(str,"NodeGradTemp");
    writeVtkProp(&idum, nGradTemp, nnode, ndm, str,iws
               , DOUBLE_VTK, VECTORS_VTK, f);
  }
/*...................................................................*/

/*... escreve a vorticidade por no*/  
  if(opt->vorticity && opt->fNode)
  {
    strcpy(str,"NodeVorticity");
    HccaAlloc(DOUBLE,m,p,nnode*3,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeVorticity(p,nGradVel,nnode,ndm);
    writeVtkProp(&idum,p,nnode,ndm,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a tensor desviador por no*/ 
  if(opt->stress && opt->fNode ){
/*  strcpy(str,"NodeStress");
    HccaAlloc(DOUBLE,m,p,nnode*ntn,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeStress(p,nGradVel,nDyViscosity,nnode,ndm,true);
    writeVtkProp(&idum,p,nnode,ntn,str,iws
                ,DOUBLE_VTK,TENSORS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);*/
  }
/*...................................................................*/

/*... escreve a energia cinetica por no*/  
  if(opt->kinetic && opt->fNode )
  {
    strcpy(str,"NodeKinetic");
    HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeKineticEnergy(p,nVel,nDensityFluid,nnode,ndm);
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve o tensor residual por no*/  
  if(opt->stressR && opt->fNode )
  {
    strcpy(str,"nStressRs");
    writeVtkProp(&idum,nStressR,nnode,6,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaAlloc(DOUBLE,m,p,numel*ntn,"p",_AD_);
    makeStress(p,nGradVel,nDyViscosity,nnode,ndm,ntn,false);
    strcpy(str,"nStressRf");
    writeVtkProp(&idum,p,nnode,ntn,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... coeficiente dinmicamente calculados por no*/  
  if(opt->cDynamic && opt->fNode)
  {
    strcpy(str,"nCdyn");
    writeVtkProp(&idum,nCd,nnode,2,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a yPlus por no*/   
  if(opt->wallParameters && opt->fNode)
  {
    strcpy(str,"nWallParameters(y+|u+|uf|sW)");
    writeVtkProp(&idum,nWallPar,nnode,4,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a energia cinetica por no*/   
  if(opt->Qcriterion &&  opt->fNode)
  {
    strcpy(str,"nQCriterion");
    HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makeQcriterion(p,nGradVel,nnode,ndm);
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a pressao total por no*/  
  if(opt->presTotal &&  opt->fNode )
  {
    strcpy(str,"nPresTotal");
    HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
    makePresTotal(p,nPres,nVel,nDensityFluid,nnode,ndm);
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a pressao total por no*/  
  if(opt->kTurb &&  opt->fNode )
  {
    strcpy(str,"nKturbl");
    writeVtkProp(&idum,nKturb,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a vicosidade  por no*/  
  if(opt->dViscosity && opt->fNode ){
    strcpy(str,"nDviscosity");
    writeVtkProp(&idum,nDyViscosity,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a condutividae termica por no*/  
  if(opt->tConductivity && opt->fNode ){
    strcpy(str,"nThermoCondutivity");
    writeVtkProp(&idum,nTCond,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve o calor especifico por no*/  
  if(opt->specificHeat && opt->fNode ){
    strcpy(str,"nSpecificHeat");
    writeVtkProp(&idum,nSheat,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*...  escreve a massa especifica por no*/  
  if(opt->densityFluid && opt->fNode ){
    strcpy(str,"nRho");
    writeVtkProp(&idum,nDensityFluid,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);

  }
/*...................................................................*/

/*... escreve a zComb*/  
  if(opt->zComb &&  opt->fNode )
  {
    HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);

    getColFromMatrix(p,nZcomb,nnode,3,SL_FUEL); 
    strcpy(str,"nZair");
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);

    getColFromMatrix(p,nZcomb,nnode,3,SL_AIR); 
    strcpy(str,"nZfuel");
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);

    getColFromMatrix(p,nZcomb,nnode,3,SL_PROD); 
    strcpy(str,"nZprod");
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);

    sumFracZ(p,nZcomb,nnode,3);
    strcpy(str,"nZTotal");
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);

    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... escreve a zComb*/  
  if(opt->gradZcomb &&  opt->fNode)
  {
    strcpy(str,"nGradZcomb");
    writeVtkProp(&idum,nGradZcomb,nnode,3*ndm,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a energia cinetica turbulenta */  
  if(opt->wk &&  opt->fNode )
  {
    strcpy(str,"nWk");
    writeVtkProp(&idum,nWk,nnode,nSp,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a energia cinetica turbulenta */  
  if(opt->rateHeatComb &&  opt->fNode )
  {
    strcpy(str,"nRateHeatComb");
    writeVtkProp(&idum,nHeatRe,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escreve a fracao massica por elemento*/  
  if(opt->yFrac && opt->fNode)
  {
    HccaAlloc(DOUBLE,m,p,nnode,"p",_AD_);
    ERRO_MALLOC(p,"p",__LINE__,__FILE__,__func__);
/*... Fuel*/
    for(i=0;i<cModel->chem.nSp;i++)
    {
      strcpy(str,"nY");
      strcat(str,cModel->chem.sp[i].name);
      getColFromMatrix(p,nYfrac,nnode,cModel->chem.nSp,i); 
      writeVtkProp(&idum,p,nnode,1,str,iws
                  ,DOUBLE_VTK,SCALARS_VTK,f);
    }

    sumFracZ(p,nYfrac,nnode, cModel->chem.nSp);
    strcpy(str,"nYTotal");
    writeVtkProp(&idum,p,nnode,1,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);

    HccaDealloc(m,p,"p",_AD_);
  }
/*...................................................................*/

/*... coeficiente de difusao das especies */  
  if(opt->coefDiffSp &&  opt->fNode )
  {
    strcpy(str,"nCoefDiffSp");
    writeVtkProp(&idum,nDiffSp,nnode,cModel->nOfSpecies,str,iws
                ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*... escrever gradiente da densidade por celula*/  
  if(opt->gradRho && opt->fNode){
    strcpy(str,"nGradRho");
    writeVtkProp(&idum,nGradRho,nnode,ndm,str ,iws
                ,DOUBLE_VTK,VECTORS_VTK,f);
  }
/*...................................................................*/

/*... escrever a massa molar da mistura por celula*/
  if (opt->mMolar && opt->fNode) 
  {
    strcpy(str,"nMolar");
    writeVtkProp(&idum,nMolar,nnode,1,str,iws
               ,DOUBLE_VTK,SCALARS_VTK,f);
  }
/*...................................................................*/

/*...*/
  fclose(f);
/*...................................................................*/
}
/*********************************************************************/

/**********************************************************************
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 06/10/2019                                    *
 * ------------------------------------------------------------------ *
 * MAKEFACE : gera as faces/aresta onde ha carregamento               *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * el      -> conectividade                                           * 
 * faceR   -> tipo de condicao de contorno                            * 
 * faceS   -> valor de condicao de contorno                           * 
 * typeGeom-> tipo geometrico do elemento                             *
 * face    -> indefinido                                              *
 * idFace  -> indefinido                                              *
 * tyGeomF -> indefinido                                              *
 * nenFace -> indefinido                                              *
 * maxViz  -> numero maximo de vizinho por celula                     *
 * maxNo   -> numero maximo de nos por celula                         *
 * numel   -> numero de elementos                                     *
 * nFace   -> indefinido                                              *
 * fWallVel-> true considera a parede impermeavel                     *
 *            false desconsidera                                      *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * face    -> faces onde ha condicao de contorno                      *
 * lFaceL  -> tipo de condico de contorno na face                     *
 * tyGeomF -> tipo geometrico da face                                 *
 * nenFace -> numero de nos por face                                  *
 * nFace   -> numero total de faces                                   *
 * ------------------------------------------------------------------ *
 **********************************************************************/
void makeFace(INT *el            ,Loads *ld
             ,short *faceR       ,short *typeGeom
             ,INT *face          ,int    *lFaceL   
             ,int    *lFaceTy    ,INT *idFace
             ,short *typeGeomFace,short *nenFace
             ,short const maxViz ,short const maxNo
             ,short const ndf    ,INT const numel  
             ,INT *nFace         )
{

  short  isnod[MAX_SN],nenFaceMax=MAX_NUM_NODE_FACE;
  short j,k,ty;
  unsigned int nf = 0;
  short cCell = maxViz + 1;
  int no,nel,tmp;
  
  for(nel=0;nel<numel;nel++)
  {
    ty = typeGeom[nel];

/*... triangulos*/
    if( ty == TRIACELL)
    {
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for(j=0;j<maxViz;j++)
        tmp += abs(MAT2D(nel,j,faceR,cCell));
      if(tmp){
/*... loop nas arestas*/
        for(j=0;j<3;j++)
        {
          tmp = MAT2D(nel,j,faceR,cCell);
/*...*/
          if(tmp){
            idFace[nf]       = nel + 1;
            typeGeomFace[nf] = LINECELL;
            nenFace[nf]      = sn(isnod,ty,nel);
            lFaceL[nf]       = MAT2D(nel,j,faceR,cCell); 
            lFaceTy[nf]      = ld[lFaceL[nf]-1].type;
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
    else if( ty == QUADCELL)
    {
/*...*/
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for(j=0;j<maxViz;j++)
        tmp += abs(MAT2D(nel,j,faceR,cCell));
/*...*/
      if(tmp)
      {
/*... loop nas arestas*/
        for(j=0;j<4;j++)
        {
          tmp = MAT2D(nel,j,faceR,cCell);
/*...*/
          if(tmp)
          {
            idFace[nf]       = nel + 1;
            typeGeomFace[nf] = LINECELL;
            nenFace[nf]      = sn(isnod,ty,nel);
            lFaceL[nf]       = MAT2D(nel,j,faceR,cCell); 
            lFaceTy[nf]      = ld[lFaceL[nf]-1].type;
/*... loop nos nos da face*/
            for(k=0;k<nenFace[nf];k++)
            {
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
    else if( ty == TETRCELL)
    {
/*...*/
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for (j = 0; j<maxViz; j++)
        tmp += abs(MAT2D(nel, j, faceR, cCell));
/*...................................................................*/

/*...*/
      if(tmp)
      {
/*... loop nas faces*/
        for(j=0;j<4;j++)
        {
          tmp = MAT2D(nel, j, faceR, cCell);
/*...*/
          if(tmp)
          {
            idFace[nf]       = nel + 1;
            typeGeomFace[nf] = TRIACELL;
            nenFace[nf]      = sn(isnod,ty,nel);
            lFaceL[nf]       = MAT2D(nel,j,faceR,cCell); 
            lFaceTy[nf]      = ld[lFaceL[nf]-1].type;
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
    else if( ty == HEXACELL)
    {
/*...*/
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for(j=0;j<maxViz;j++)
        tmp += abs(MAT2D(nel,j,faceR,cCell));
/*...................................................................*/
  
/*...*/
      if(tmp)
      {
/*... loop nas faces*/
        for(j=0;j<6;j++)
        {
          tmp = MAT2D(nel,j,faceR,cCell);
/*...*/
          if(tmp)
          {
            idFace[nf]       = nel + 1;
            typeGeomFace[nf] = QUADCELL;
            nenFace[nf]      = sn(isnod,ty,nel); 
            lFaceL[nf]       = MAT2D(nel,j,faceR,cCell); 
            lFaceTy[nf]      = ld[lFaceL[nf]-1].type;
/*... loop nos nos da face*/
            for(k=0;k<nenFace[nf];k++)
            {
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
    else if (ty == PIRACELL) 
    {
/*...*/
      tmp = 0;
/*... checa se ha carga nas faces das celula*/
      for (j = 0; j<maxViz; j++)
        tmp += abs(MAT2D(nel, j, faceR, cCell));
/*...*/
      if (tmp)
      {
/*... loop nas faces*/
        for (j = 0; j<5; j++)
        {
          tmp = MAT2D(nel, j, faceR, cCell);
/*...*/
          if (tmp)
          {
            idFace[nf] = nel + 1;
            if(j)
              typeGeomFace[nf] = TRIACELL;
            else
              typeGeomFace[nf] = QUADCELL;
            nenFace[nf] = sn(isnod, ty, nel);
            lFaceL[nf] = MAT2D(nel, j, faceR, cCell);
            lFaceTy[nf]      = ld[lFaceL[nf]-1].type;
/*... loop nos nos da face*/
            for (k = 0; k<nenFace[nf]; k++)
            {
              no = MAT2D(j, k, isnod, nenFace[nf]);
              MAT2D(nf, k, face, nenFaceMax) = MAT2D(nel, no, el, maxNo) - 1;
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

/**********************************************************************
 * Data de criacao    : 11/11/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 *------------------------------------------------------------------- * 
 * makeVorticity : calculo do campo de vorticidade                    *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * w       -> nao definido                                            * 
 * gradVel -> gradienta das velocidades                               * 
 * n       -> numero de elementos/nos                                 *
 * ndm     -> dimensao                                                * 
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * w       -> vorticidade                                             *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *------------------------------------------------------------------- *
 *                                                                    *
 * gradVel(*,ndf,ndm)                                                 *
 *                                                                    *
 *               | du1dx1 du1dx2 du1dx3 |                             * 
 * grad(*,*,*) = | du2dx1 du2dx2 du2dx3 |                             *
 *               | du3dx1 du3dx2 du3dx3 |                             *
 *                                                                    *
 **********************************************************************/
void makeVorticity(DOUBLE *RESTRICT w, DOUBLE *RESTRICT gradVel
                  ,INT const n       , const short ndm) {
  short j;
  INT i;
  DOUBLE v[3],*p;

  for (i = 0; i < n; i++) {
//      vorticity(v,&MAT3D(i,0,0,gradVel,ndm,ndm),ndm);
      p = gradVel + i*ndm*ndm;
      vorticity(v,p,ndm);
      for (j=0;j<ndm;j++)
        MAT2D(i,j,w,ndm) =  v[j];
  }

}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 11/11/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 *------------------------------------------------------------------- * 
 * makeStress : campo de tensoes viscosas desviadoras                 *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * str       -> nao definido                                          * 
 * gradVel   -> gradienta das velocidades                             * 
 * viscosity -> viscosidade molecular                                 * 
 * n         -> numero de pontos                                      * 
 * ndm       -> dimensao                                              * 
 * ntn       -> numetro de termos no tensor simetrico (4;6)           *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * str     -> tensor de forcas viscosas                               *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *------------------------------------------------------------------- *
 *                                                                    *
 * gradVel(*,ndf,ndm)                                                 *
 *                                                                    *
 *               | du1dx1 du1dx2 du1dx3 |                             * 
 * grad(*,*,*) = | du2dx1 du2dx2 du2dx3 |                             *
 *               | du3dx1 du3dx2 du3dx3 |                             *
 *                                                                    *
 **********************************************************************/
void makeStress(DOUBLE *RESTRICT str      , DOUBLE *RESTRICT gradVel
              , DOUBLE *RESTRICT viscosity,INT const n          
              , short const ndm           , short const ntn           
              , bool const flag ) {
  INT i;
  DOUBLE *p,*s,tmp;

/*... tensao desviador laminar*/
  if(flag)
    for (i = 0; i < n; i++) {
      p = gradVel + i*ndm*ndm;
      s = str     + i*ntn;
      tmp = -D2DIV3*viscosity[i];
      stress(s           ,p
            ,viscosity[i],tmp
            ,ndm);
    }
/*.....................................................................*/

/*... tensao desviador turbulenta*/
  else
    for (i = 0; i < n; i++) {
      p = gradVel + i*ndm*ndm;
      s = str     + i*ntn;
      stressEddyViscosity(s           ,p
                         ,viscosity[i],ndm);
    }
/*.....................................................................*/
}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 30 11/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 *------------------------------------------------------------------- * 
 * makeKineticEnergy : campo de energia cinetica                      *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * e         -> energia cinetica                                      * 
 * vel       -> campo de velocidade                                   * 
 * density   -> densidade do fluido                                   * 
 * n         -> numero de pontos                                      * 
 * ndm       -> dimensao                                              * 
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * e       -> energia cinetica espexifica                             *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *------------------------------------------------------------------- *
 **********************************************************************/
void makeKineticEnergy(DOUBLE *RESTRICT e      , DOUBLE *RESTRICT vel
                     , DOUBLE *RESTRICT density 
                     , INT const n             , short const ndm) {
  INT i;
  DOUBLE vv,den,v[3];

  for (i = 0; i < n; i++) {
    den  =  density[i];
    v[0] =  MAT2D(i, 0, vel, ndm);
    v[1] =  MAT2D(i, 1, vel, ndm);
    vv = v[0]*v[0] + v[1]*v[1];
    if(ndm == 3){
      v[2] =  MAT2D(i, 2, vel, ndm);
      vv += v[2]*v[2];
    }
    e[i] = 0.5e0*den*vv;
  }

}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 13/12/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 *------------------------------------------------------------------- * 
 * makeQcriterion : criterio de vorticiade Q                          *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * gradVel -> gradienbte de velocidades                               *
 * n         -> numero de pontos                                      * 
 * ndm       -> dimensao                                              * 
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * q -> retorna o campo Q                                             *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *------------------------------------------------------------------- *
 **********************************************************************/
void makeQcriterion( DOUBLE *RESTRICT q, DOUBLE *RESTRICT gradVel
                  , INT const n        , short const ndm)
{
  INT i;
  DOUBLE *p;

  for (i = 0; i < n; i++) {
    p = gradVel + i*ndm*ndm;
    q[i] = qCriterion(p,ndm);
  }

}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 08/01/2018                                    *
 * Data de modificaco : 00/00/0000                                    *
 *------------------------------------------------------------------- * 
 * makeKineticEnergy : campo de energia cinetica                      *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * presT     -> nao definido                                          * 
 * vel       -> campo de velocidade                                   * 
 * pres      -> pressao estatica                                      * 
 * density   -> densidade do fluido                                   * 
 * n         -> numero de pontos                                      * 
 * ndm       -> dimensao                                              * 
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * presT   -> pressao total                                           *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *------------------------------------------------------------------- *
 **********************************************************************/
void makePresTotal(DOUBLE *RESTRICT presT, DOUBLE *RESTRICT pres
                 , DOUBLE *RESTRICT vel  , DOUBLE *RESTRICT density 
                 , INT const n           , short const ndm) {
  INT i;
  DOUBLE vv,den,v[3];

  for (i = 0; i < n; i++) {
    den  =  density[i];
    v[0] =  MAT2D(i, 0, vel, ndm);
    v[1] =  MAT2D(i, 1, vel, ndm);
    vv = v[0]*v[0] + v[1]*v[1];
    if(ndm == 3){
      v[2] =  MAT2D(i, 2, vel, ndm);
      vv += v[2]*v[2];
    }
    presT[i] = pres[i] +0.5e0*den*vv;
  }

}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 29/01/2018                                    *
 * Data de modificaco : 00/00/0000                                    *
 *------------------------------------------------------------------- * 
 * makeModuleVel : campo de energia cinetica                          *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * p         -> nao definido                                          * 
 * vel       -> campo de velocidade                                   * 
 * n         -> numero de pontos                                      * 
 * ndm       -> dimensao                                              * 
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * p -> sqrt(v1*v1+v1*v1+v1*v1)                                       *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *------------------------------------------------------------------- *
 **********************************************************************/
void makeModuleVel(DOUBLE *RESTRICT p,DOUBLE *RESTRICT vel
                 , INT const n           , short const ndm) {
  INT i;
  DOUBLE vv,v[3];

  for (i = 0; i < n; i++) {
    v[0] =  MAT2D(i, 0, vel, ndm);
    v[1] =  MAT2D(i, 1, vel, ndm);
    vv = v[0]*v[0] + v[1]*v[1];
    if(ndm == 3){
      v[2] =  MAT2D(i, 2, vel, ndm);
      vv += v[2]*v[2];
    }
    p[i] = sqrt(vv);
  }

}
/**********************************************************************/
