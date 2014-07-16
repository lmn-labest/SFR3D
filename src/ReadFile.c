#include<ReadFile.h>
/*********************************************************************
 * readFileFc : leitura de arquivo de dados em volume finitos        *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * m     -> memoria principal                                        * 
 * mesh  -> malha                                                    *
 * FILE  -> arquivo de entrada                                       * 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * *******************************************************************/
void readFileFvMesh(Memoria *m,Mesh *mesh, FILE* file)
{
  char word[WORD_SIZE],str[WORD_SIZE];
  char macro[NMACROS][WORD_SIZE]={
            "coordinates","endMesh"  ,"insert"            /*0,1,2*/
           ,"return"     ,"cells"    ,"faceRt1"           /*3,4,5*/
           ,"faceSt1"    
	   };                                             
  bool rflag[NMACROS],macroFlag;
  INT nn,nel;
  short maxno,ndm,numat,maxViz,ndf;
  int i;

/* leitura dos parametros principais da malha*/
  parametros(&nn,&nel,&maxno,&ndm,&numat,&ndf,file);
  mesh->nnode  = nn;
  mesh->numel  = nel;
  mesh->maxNo  = maxno;
  mesh->maxViz = maxViz = maxno;
  mesh->ndm    = ndm;
  mesh->numat  = numat;
  mesh->ndfT[0] = ndf;

/*... alocando variavies de elementos*/
/*... conectividade*/ 
  Myalloc(INT,m,mesh->elm.node       ,nel*maxno        ,"elnode",_AD_);
/*... materiais*/ 
  Myalloc(short,m,mesh->elm.mat      ,nel              ,"elmat" ,_AD_);
/*... nos por elementos*/
  Myalloc(short,m,mesh->elm.nen      ,nel              ,"elnen" ,_AD_);
/*... tipo geometrico */
  Myalloc(short,m,mesh->elm.geomType ,nel              ,"elgT"  ,_AD_);
/*... centroide */
  Myalloc(double,m,mesh->elm.geom.cc ,nel*ndm          ,"elCc"  ,_AD_);
/*... vetor que une os centroides dos elementos */
  Myalloc(double               ,m       ,mesh->elm.geom.ksi
         ,nel*ndm*maxViz,"elksi"  ,_AD_);
/*... modulo do vetor que une os centroides dos elementos */
  Myalloc(double               ,m       ,mesh->elm.geom.mksi
        ,nel*maxViz     ,"elmksi",_AD_);
/*... vetor paralelo a face da celula */
  Myalloc(double               ,m       ,mesh->elm.geom.eta
         ,nel*ndm*maxViz,"eleta"  ,_AD_);
/*... modulo do vetor paralelo a face da celula */
  Myalloc(double               ,m       ,mesh->elm.geom.meta
        ,nel*maxViz     ,"elmeta",_AD_);
/*... volume da celula*/                           
  Myalloc(double               ,m       ,mesh->elm.geom.volume
        ,nel            ,"elVol",_AD_);
/*... vetor normal a face da celula*/                           
  Myalloc(double               ,m       ,mesh->elm.geom.normal
         ,nel*maxViz*ndm       ,"elnorm",_AD_);
/*... ponto medio da face*/                           
  Myalloc(double               ,m       ,mesh->elm.geom.xm
         ,nel*maxViz*ndm       ,"elxm",_AD_);
/*... vetor que une o centroide ao ponto medio*/                           
  Myalloc(double               ,m       ,mesh->elm.geom.xmcc
         ,nel*maxViz*ndm       ,"elxmcc",_AD_);
/*... distancia entro o ponto medio da face e ponto de intercao 
      entre a linha entre os centroides e a face*/         
  Myalloc(double               ,m       ,mesh->elm.geom.mkm  
         ,nel*maxViz           ,"elmkm" ,_AD_);
/*... distancia entro o ponto medio da face e ponto de intercao 
      entre a linha entre os centroides e a face*/         
  Myalloc(double               ,m       ,mesh->elm.geom.dcca 
         ,nel*maxViz           ,"eldcca",_AD_);

/*... zerando os variavies*/
  zero(mesh->elm.node         ,nel*maxno     ,INTC);
  zero(mesh->elm.mat          ,nel           ,"short"  );
  zero(mesh->elm.nen          ,nel           ,"short"  );
  zero(mesh->elm.geomType     ,nel           ,"short"  );
  zero(mesh->elm.geom.cc      ,nel*ndm       ,"double" );
  zero(mesh->elm.geom.ksi     ,nel*ndm*maxViz,"double" );
  zero(mesh->elm.geom.mksi    ,nel*maxViz    ,"double" );
  zero(mesh->elm.geom.eta     ,nel*ndm*maxViz,"double" );
  zero(mesh->elm.geom.meta    ,nel*maxViz    ,"double" );
  zero(mesh->elm.geom.volume  ,nel           ,"double" );
  zero(mesh->elm.geom.normal  ,nel*ndm*maxViz,"double" );
  zero(mesh->elm.geom.xm      ,nel*ndm*maxViz,"double" );
  zero(mesh->elm.geom.xmcc    ,nel*ndm*maxViz,"double" );
  zero(mesh->elm.geom.mkm     ,nel*maxViz    ,"double" );
  zero(mesh->elm.geom.dcca    ,nel*maxViz    ,"double" );
/*...................................................................*/


/*... alocando materiais*/
/*... Prop*/ 
  Myalloc(double,m,mesh->material.prop,MAXPROP*numat     
         ,"prop" ,_AD_);
/*... type*/ 
  Myalloc(short,m,mesh->material.type,numat     
         ,"type" ,_AD_);
/*... zerando os variavies*/
  zero(mesh->material.prop,MAXPROP*numat,"double");
  zero(mesh->material.type,numat,"short");
/*...................................................................*/

/*... alocando estruturas para vizinhos*/
/*... nelcon*/ 
  Myalloc(INT,m,mesh->adj.nelcon,nel*maxno,"adj" ,_AD_);
/*... type*/ 
  Myalloc(short,m,mesh->adj.nViz,nel       ,"nViz" ,_AD_);
/*... zerando os variavies*/
  zero(mesh->adj.nelcon,nel*maxno,INTC);
  zero(mesh->adj.nViz  ,nel      ,"short");
/*...................................................................*/

/*---alocando variaveis nodais */      
/*---alocando coordenadas      */      
  Myalloc(double,m,mesh->node.x,ndm*nn ,"xnode",_AD_);   
     
/*... zerando os variavies*/
  zero(mesh->node.x,ndm*nn,"double" );
/*...................................................................*/
   if(mesh->ndfT[0] > 0) {   
/*... alocando memoria*/
     Myalloc(short,m,mesh->elm.faceRt1
            ,nel*(maxno+1)*mesh->ndfT[0],"faceRt1"  ,_AD_);
     Myalloc(double,m,mesh->elm.faceSt1
            ,nel*(maxno+1)*ndf,"faceSt1"  ,_AD_);
     zero(mesh->elm.faceRt1  ,nel*(maxno+1)*mesh->ndfT[0],"short"  );
     zero(mesh->elm.faceSt1  ,nel*(maxno+1)*mesh->ndfT[0],"double" );
/*...................................................................*/
   }
/*...................................................................*/

/*... leitura das macros*/
  for(i=0;i<NMACROS;i++)
    rflag[i] = false;
  
  nmacro = 0;
  for(i=0;i<MAX_LINE;i++)
    strcpy(macros[i],"");
  
  macroFlag = true;
  do{
/*...*/
    readMacro(file,word,false);
/*...................................................................*/

/*... coordinates*/
    if((!strcmp(word,macro[0])) && (!rflag[0])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[0] = true;
      printf("loading coordinates...\n");
      readVfCoor(mesh->node.x,mesh->nnode,mesh->ndm,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... endMesh*/
    else if((!strcmp(word,macro[1])) && (!rflag[1])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      printf("%s\n\n",DIF);
      strcpy(macros[nmacro++],word);
      rflag[1] = true;
      macroFlag = false;    
    }
/*...................................................................*/


/*... cells  */
    else if((!strcmp(word,macro[4])) && (!rflag[4])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[4] = true;
      printf("loading cells ...\n");
      readVfElmt(mesh->elm.node    ,mesh->elm.mat
                ,mesh->elm.nen     ,mesh->adj.nViz
                ,mesh->elm.geomType,mesh->numel
                ,mesh->maxNo       ,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/


/*... faceRt1 */
    else if((!strcmp(word,macro[5])) && (!rflag[5])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[5] = true;
      strcpy(str,"endFaceRt1");
      printf("loading faceRt1 ...\n");
      readVfRes(mesh->elm.faceRt1,mesh->numel,mesh->maxViz,str,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceSt1 */
    else if((!strcmp(word,macro[6])) && (!rflag[6])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[6] = true;
      strcpy(str,"endFaceSt1");
      printf("loading faceSt1 ...\n");
      readVfSource(mesh->elm.faceSt1,mesh->numel,mesh->maxViz
                  ,str,file);  
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/
  }while(macroFlag && (!feof(file)));

}

/*********************************************************************
 * PARAMETROS: leitura dos parametros do problema                    *
 *********************************************************************
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * nn    - numero de nos                                             *
 * nel   - numero de elementos                                       *
 * nen   - numero maximo de nos por elemento                         *
 * numat - numero de materias                                        *
 * ndf   - graus de liberdade (onda)                                 *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 *********************************************************************/
void parametros(INT  *nn   ,INT *nel  ,short *maxno
               ,short *ndm  ,short *numat,short *ndf
               ,FILE* file)
{
  char parameter[][WORD_SIZE]={"nnode","numel","numat"
                              ,"maxno","ndf"  ,"ndm"};
  
  char word[WORD_SIZE];
  bool flag[NPARAMETROS];
  int i=0,j;
  long aux;


  *nn    = 0;
  *nel   = 0;
  *numat = 0;
  *maxno = 0;
  *ndf   = 0;
  *ndm   = 0;

  for(j=0;j<NPARAMETROS;j++)
    flag[j] = false;

  while(i < NPARAMETROS && i < 20){
    readMacro(file,word,false);
/*... macro nnode*/   
    if(!strcmp(word,parameter[0])){
      fscanf(file,"%ld",&aux);
      *nn = (INT) aux;
#ifdef _DEBUG_MESH_ 
      printf("nnode %ld\n",*nn);
#endif      
      flag[0] = true;
      i++;
    }
/*...................................................................*/

/*... macro numel*/   
    else if(!strcmp(word,parameter[1])){
      fscanf(file,"%ld",&aux);
      *nel = (INT) aux;
#ifdef _DEBUG_MESH_ 
      printf("numel %ld\n",*nel);
#endif      
      flag[1] = true;
      i++;
    }
/*... macro numat*/   
    else if(!strcmp(word,parameter[2])){
      fscanf(file,"%hd",numat);
#ifdef _DEBUG_MESH_ 
      printf("numat %hd\n",*numat);
#endif      
      flag[2] = true;
      i++;
    }
/*...................................................................*/

/*... macro maxno*/   
    else if(!strcmp(word,parameter[3])){
      fscanf(file,"%hd",maxno);
#ifdef _DEBUG_MESH_ 
      printf("maxno %hd\n",*maxno);
#endif      
      flag[3] = true;
      i++;
    }
/*...................................................................*/

/*... macro ndf*/   
    else if(!strcmp(word,parameter[4])){
      fscanf(file,"%hd",ndf);
#ifdef _DEBUG_MESH_ 
      printf("ndf %hd\n",*ndf);
#endif      
      flag[4] = true;
      i++;
    }
/*...................................................................*/

/*... ndm*/
    else if(!strcmp(word,parameter[5])){
      fscanf(file,"%hd",ndm);
#ifdef _DEBUG_MESH_ 
      printf("ndm %hd\n",*ndm);
#endif      
      flag[5] = true;
      i++;
    }
/*...................................................................*/

    else
      i++;
  }
  
  for(j=0;j<NPARAMETROS;j++){
    if(!flag[j]){
      fprintf(stderr,"parametro: %s faltando.\n"
             "fonte: %s \n",parameter[j],__FILE__);
      exit(EXIT_FAILURE);
    }
  }
}
/*********************************************************************/

/*********************************************************************
 * READVFCOOR: leitura das coordenadas                               *
 *********************************************************************
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*
 * x    -> nao definido                                              * 
 * nn   -> numero de nos                                             *
 * ndm  -> dimensao do problema                                      *
 * file -> arquivo de entrada                                        *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * x    -> coordenada                                                *
 * ------------------------------------------------------------------*
 *********************************************************************/
void readVfCoor(double *x,INT nn, short ndm,FILE *file){
  INT i,k;
  long idum;
  int j;
 
  for(i=0;i<nn;i++){
    fscanf(file,"%ld",&idum);
    k = (INT) idum -1;
    for(j=0;j<ndm;j++){
      fscanf(file,"%lf",&MAT2D(k,j,x,ndm));
    }
  }
#ifdef _DEBUG_MESH_ 
  for(i=0;i<nn;i++){
    fprintf(stderr,"%ld",i+1);
    for(j=0;j<ndm;j++){
      k = i*ndm + j;
      fprintf(stderr," %lf ",x[k]);
    }
    printf("\n");
  }
#endif
}
/*********************************************************************/

/*********************************************************************
 * READVFMAT2DT: leitura das elementos                                 *
 *********************************************************************
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*
 * el   -> nao definido                                              * 
 * mat  -> nao definido                                              *
 * nen  -> nao definido                                              *
 * nFace-> nao definido                                              *
 * ty   -> nao definido                                              *
 * nel  -> numero de elementos                                       *
 * maxno-> numero maximo de nos por elemento                         *
 * file -> arquivo de entrada                                        *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * el   -> conectividade                                             *
 * mat  -> material                                                  *
 * nen  -> numero de conectividade por elemento                      *
 * nFace-> numero de face de um elemento                             *
 * ty   -> numero do tipo geometrico do elemento                     *
 * ------------------------------------------------------------------*
 *********************************************************************/
void readVfElmt(INT *el   ,short *mat ,short *nen,short *nFace 
               ,short *ty ,INT nel    ,short maxno
               ,FILE *file){
  INT i;
  long idum,aux;
  short nenl,face;
  int j;
 
  for(i=0;i<nel;i++){
    fscanf(file,"%ld",&idum);
    idum--;
    fscanf(file,"%hd",&mat[idum]);
    fscanf(file,"%hd",&ty[idum]);
    fscanf(file,"%hd",&nenl);
    nen[idum] = nenl;
    fscanf(file,"%hd",&face);
    nFace[idum] = face;
    for(j=0;j<nenl;j++){
      fscanf(file,"%ld",&aux);
      MAT2D(idum,j,el,maxno) = (INT) aux;
    }
  }
}
/*********************************************************************/

/*********************************************************************
 * READVFRES : leitura das restricoes                                *
 *********************************************************************
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * id    - indefinido                                                *
 * numel - numero de elementos                                       *
 * maxno - numero maximo de no por elementos                         *
 * str   - macro de terminada o fim da secao                         *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * id    - tipos de restricoes                                       *
 *********************************************************************/
void readVfRes(short *id,INT numel,short maxno
              ,char *str    ,FILE* file){
  
  char word[WORD_SIZE];
  int  j,k,kk,nTerm;
  INT nel;
  long aux;  

  readMacro(file,word,false);
  do{
    nel = atol(word);  
    fscanf(file,"%d",&nTerm);
    kk = (nel-1)*nTerm;
    for(j = 0;j < nTerm;j++){
      k = kk + j;
      if ( nel > 0 && nel <= numel){ 
        fscanf(file,"%hd",&id[k]);
      } 
      else{
        aux = (long) nel;
        printf("Erro: numero do elemento nao exitentes. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
               ,aux,__FILE__,__func__);
        exit(EXIT_FAILURE);
      }
    }
  readMacro(file,word,false);
  }while(strcmp(word,str));
}
/*********************************************************************/

/*********************************************************************
 * READVFSOURCE : leitura dos valores das restricoes                 *
 *********************************************************************
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * f     - indefinido                                                *
 * numel - numero de elementos                                       *
 * maxno - numero maximo de no por elementos                         *
 * str   - macro de terminada o fim da secao                         *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * f     - valores das restricoes                                    *
 *********************************************************************/
void readVfSource(double *f,INT numel, short int maxno,char *str
                 ,FILE* file){
  
  char word[WORD_SIZE];
  int  j,k,kk,nTerm;
  INT nel;
  long aux;  

  readMacro(file,word,false);
  do{
    nel = atol(word);  
    fscanf(file,"%d",&nTerm);
    kk = (nel-1)*nTerm;
    for(j = 0;j < nTerm;j++){
      k = kk + j;
      if ( nel > 0 && nel <= numel){ 
        fscanf(file,"%lf",&f[k]);
      }  
      else{
        aux = (long) nel;
        printf("Erro: numero do elemento nao exitentes. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
              ,aux,__FILE__,__func__);
        exit(EXIT_FAILURE);
      }
    }
  readMacro(file,word,false);
  }while(strcmp(word,str));
}


/*********************************************************************
 * CONFIG : configuraceos gerais                                     *
 *********************************************************************
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * bvtk      - arquivo binario para o vtk                            *
 * reordMesh - reordenacao do malha                                  *
 * m         - memoria principal                                     *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 *********************************************************************/
void config(bool *bvtk,Reord *reordMesh,FILE* file)
{
  char config[][WORD_SIZE]={"bvtk","reord","mem"};
  
  char word[WORD_SIZE];
  bool flag[NCONFIG];
  int i=0,j,temp;
  double conv;


  *bvtk            = false;
  reordMesh -> flag= false;

  for(j=0;j<NCONFIG;j++)
    flag[j] = false;

  while(i < NCONFIG && i < 20){
    readMacro(file,word,false);
/*... bvtk*/   
    if(!strcmp(word,config[0])){
      fscanf(file,"%d",&temp);
      *bvtk = (bool) temp;
      flag[0] = true;
      i++;
      if(*bvtk)
        printf("bvtk: true\n");
      else
        printf("bvtk: false\n");
    }
/*... reord*/   
    else if(!strcmp(word,config[1])){
      fscanf(file,"%d",&temp);
      reordMesh -> flag= (bool) temp;
      flag[1] = true;
      i++;
      if(reordMesh->flag)
        printf("Reord: true\n");
      else
        printf("Reod: false\n");
    }
/*... mem*/   
    else if(!strcmp(word,config[2])){
      fscanf(file,"%d",&temp);
      conv    = CONV_BYTES*CONV_BYTES;
      nmax    = temp*conv;                  
      flag[2] = true;
      i++;
      printf("Memoria principal: %d MBytes\n"
            ,(int)(nmax/conv));
    }
    else
      i++;
/*...................................................................*/
  }
  
  for(j=0;j<NCONFIG;j++){
    if(!flag[j]){
      fprintf(stderr,"%s: %s faltando.\n"
             "fonte: %s \n",__func__,config[j],__FILE__);
      exit(EXIT_FAILURE);
    }
  }
}
/*********************************************************************/
