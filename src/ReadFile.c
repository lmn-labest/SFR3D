#include<ReadFile.h>

/*...funcao de apoio*/
  static void getword(char *line, char*word);
  static int getnumprop2(char *line);
/*..................................................................*/

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
            "coordinates","endMesh"  ,"insert"            /* 0, 1, 2*/
           ,"return"     ,"cells"    ,"faceRt1"           /* 3, 4, 5*/
           ,"faceSt1"    ,""         ,""                  /* 6, 7, 8*/ 
           ,""           ,""         ,""                  /* 9,10,11*/ 
           ,"faceRd1"    ,"faceSd1"   ,""                 /*12,13,14*/ 
           ,""           ,""         ,""                  /*15,16,17*/ 
           ,""           ,""         ,""                  /*18,19,20*/ 
           ,"materials"  ,""         ,""                  /*21,22,23*/ 
	   };                                             
  bool rflag[NMACROS],macroFlag;
  INT nn,nel;
  short maxno,ndm,numat,maxViz,ndf;
  int i;

/* leitura dos parametros principais da malha*/
  parametros(&nn,&nel,&maxno,&maxViz,&ndm,&numat,&ndf,file);
  mesh->nnode  = nn;
  mesh->numel  = nel;
  mesh->maxNo  = maxno;
  mesh->maxViz = maxViz; 
  mesh->ndm    = ndm;
  mesh->numat  = numat;
/*mesh->ndfT[0] = ndf;*/
  mesh->ndfD[0] = ndf;

/*... alocando variavies de elementos*/
/*... conectividade*/ 
  HccaAlloc(INT,m,mesh->elm.node       ,nel*maxno        ,"elnode",_AD_);
/*... materiais*/ 
  HccaAlloc(short,m,mesh->elm.mat      ,nel              ,"elmat" ,_AD_);
/*... nos por elementos*/
  HccaAlloc(short,m,mesh->elm.nen      ,nel              ,"elnen" ,_AD_);
/*... tipo geometrico */
  HccaAlloc(short,m,mesh->elm.geomType ,nel              ,"elgT"  ,_AD_);
/*... centroide */
  HccaAlloc(double,m,mesh->elm.geom.cc ,nel*ndm          ,"elCc"  ,_AD_);
/*... vetor que une os centroides dos elementos */
  HccaAlloc(double               ,m       ,mesh->elm.geom.ksi
         ,nel*ndm*maxViz,"elksi"  ,_AD_);
/*... modulo do vetor que une os centroides dos elementos */
  HccaAlloc(double               ,m       ,mesh->elm.geom.mksi
        ,nel*maxViz     ,"elmksi",_AD_);
/*... vetor paralelo a face da celula */
  HccaAlloc(double               ,m       ,mesh->elm.geom.eta
         ,nel*ndm*maxViz,"eleta"  ,_AD_);
/*... modulo do vetor paralelo a face da celula */
  HccaAlloc(double               ,m       ,mesh->elm.geom.meta
        ,nel*maxViz     ,"elmeta",_AD_);
/*... volume da celula*/                           
  HccaAlloc(double               ,m       ,mesh->elm.geom.volume
        ,nel            ,"elVol",_AD_);
/*... vetor normal a face da celula*/                           
  HccaAlloc(double               ,m       ,mesh->elm.geom.normal
         ,nel*maxViz*ndm       ,"elnorm",_AD_);
/*... ponto medio da face*/                           
  HccaAlloc(double               ,m       ,mesh->elm.geom.xm
         ,nel*maxViz*ndm       ,"elxm",_AD_);
/*... vetor que une o centroide ao ponto medio*/                           
  HccaAlloc(double               ,m       ,mesh->elm.geom.xmcc
         ,nel*maxViz*ndm       ,"elxmcc",_AD_);
/*... distancia entro o ponto medio da face e ponto de intercao 
      entre a linha entre os centroides e a face*/         
  HccaAlloc(double               ,m       ,mesh->elm.geom.mkm  
         ,nel*maxViz           ,"elmkm" ,_AD_);
/*... distancia entro o ponto medio da face e ponto de intercao 
      entre a linha entre os centroides e a face*/         
  HccaAlloc(double               ,m       ,mesh->elm.geom.dcca 
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
  HccaAlloc(double,m,mesh->elm.material.prop,MAXPROP*numat     
         ,"prop" ,_AD_);
/*... type*/ 
  HccaAlloc(short,m,mesh->elm.material.type,numat     
         ,"type" ,_AD_);
/*... zerando os variavies*/
  zero(mesh->elm.material.prop,MAXPROP*numat,"double");
  zero(mesh->elm.material.type,numat,"short");
/*...................................................................*/

/*... alocando estruturas para vizinhos*/
/*... nelcon*/ 
  HccaAlloc(INT,m,mesh->elm.adj.nelcon,nel*maxno  ,"adj" ,_AD_);
/*... type*/ 
  HccaAlloc(short,m,mesh->elm.adj.nViz,nel       ,"nViz" ,_AD_);
/*... zerando os variavies*/
  zero(mesh->elm.adj.nelcon,nel*maxno,INTC);
  zero(mesh->elm.adj.nViz  ,nel      ,"short");
/*...................................................................*/

/*---alocando variaveis nodais */      
/*---alocando coordenadas      */      
  HccaAlloc(double,m,mesh->node.x,ndm*nn ,"xnode",_AD_);   
     
/*... zerando os variavies*/
  zero(mesh->node.x,ndm*nn,"double" );
/*...................................................................*/

/*... problema de transporte*/   
  if(mesh->ndfT[0] > 0) {   
/*... alocando memoria*/
     HccaAlloc(short,m,mesh->elm.faceRt1
            ,nel*(maxViz+1)*mesh->ndfT[0],"faceRt1"  ,_AD_);
     HccaAlloc(double,m,mesh->elm.faceSt1
            ,nel*(maxViz+1)*mesh->ndfT[0],"faceSt1"  ,_AD_);
     zero(mesh->elm.faceRt1  ,nel*(maxViz+1)*mesh->ndfT[0],"short"  );
     zero(mesh->elm.faceSt1  ,nel*(maxViz+1)*mesh->ndfT[0],"double" );
/*...................................................................*/
   }
/*...................................................................*/
   
/*... problema de difusa pura*/   
  if(mesh->ndfD[0] > 0) {   
/*... alocando memoria*/
     HccaAlloc(short,m,mesh->elm.faceRd1
            ,nel*(maxViz+1)*mesh->ndfD[0],"faceRd1"  ,_AD_);
     HccaAlloc(double,m,mesh->elm.faceSd1
            ,nel*(maxViz+1)*mesh->ndfD[0],"faceSd1"  ,_AD_);
     HccaAlloc(double,m,mesh->node.temp
            ,nel*mesh->ndfD[0],"temp"     ,_AD_);
     zero(mesh->elm.faceRd1  ,nel*(maxViz+1)*mesh->ndfD[0],"short"  );
     zero(mesh->elm.faceSd1  ,nel*(maxViz+1)*mesh->ndfD[0],"double" );
     zero(mesh->node.temp    ,nel*mesh->ndfD[0]           ,"double" );
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
                ,mesh->elm.nen     ,mesh->elm.adj.nViz
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

/*... faceRd1- condicao de contorno para problemas de difusa pura */
    else if((!strcmp(word,macro[12])) && (!rflag[12])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[12] = true;
      strcpy(str,"endFaceRd1");
      printf("loading faceRd1 ...\n");
      readVfRes(mesh->elm.faceRd1,mesh->numel,mesh->maxViz+1,str,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceSd - condicao de contorno para problemas de difusa pura */
    else if((!strcmp(word,macro[13])) && (!rflag[13])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[13] = true;
      strcpy(str,"endFaceSd1");
      printf("loading faceSd1 ...\n");
      readVfSource(mesh->elm.faceSd1,mesh->numel
                  ,(mesh->maxViz+1)*(mesh->ndfD[0])
                  ,str,file);  
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... materiais */
    else if((!strcmp(word,macro[21])) && (!rflag[21])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[21] = true;
      strcpy(str,"endMaterials");
      printf("loading materials ...\n");
      readVfMat(mesh->elm.material.prop,mesh->elm.material.type
               ,numat,file);
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
 * maxNo   -> numero de nos por celula maximo da malha               * 
 * maxViz  -> numero vizinhos por celula maximo da malha             * 
 * ndf   - graus de liberdade (onda)                                 *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 *********************************************************************/
void parametros(INT  *nn   ,INT *nel  ,short *maxNo, short *maxViz
               ,short *ndm  ,short *numat,short *ndf
               ,FILE* file)
{
  char parameter[][WORD_SIZE]={"nnode","numel","numat"
                              ,"maxno","ndf"  ,"ndm"
                              ,"maxviz"
                              };
  
  char word[WORD_SIZE];
  bool flag[NPARAMETROS];
  int i=0,j;
  long aux;


  *nn    = 0;
  *nel   = 0;
  *numat = 0;
  *maxNo = 0;
  *maxViz= 0;
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
      fscanf(file,"%hd",maxNo);
#ifdef _DEBUG_MESH_ 
      printf("maxno %hd\n",*maxNo);
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

/*... macro maxViz*/  
    else if(!strcmp(word,parameter[6])){
      fscanf(file,"%hd",maxViz);
#ifdef _DEBUG_MESH_ 
      printf("maxno %hd\n",*maxViz);
#endif      
      flag[6] = true;
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
 * READVFELMT: leitura das elementos                                 *
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
 * maxRes- numero maximo de restricoes por elemento                  *
 * str   - macro de terminada o fim da secao                         *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * id    - tipos de restricoes                                       *
 *********************************************************************/
void readVfRes(short *id,INT numel,short maxRes
              ,char *str    ,FILE* file){
  
  char word[WORD_SIZE];
  int   j,kk;
  int   nTerm;
  short res;
  INT   nel;
  long  aux; 
  int error=0; 

  readMacro(file,word,false);
  do{
    nel = atol(word);  
    aux = (long) nel;
    error = fscanf(file,"%d",&nTerm);
    if( error != 1) {
      printf("erro: leitura do numero de termos de restricao. "
             "nel = %ld.\n"
             "arquivo fonte:  \"%s\".\n"
             "nome da funcao: \"%s\".\n"
             ,aux,__FILE__,__func__);
      exit(EXIT_FAILURE);

    }
    for(j = 0;j < nTerm;j++){
      if ( nel > 0 && nel <= numel){ 
        error = fscanf(file,"%hd",&res);
        if( error != 1) {
          printf("erro: leitura da restricao. "
                 "nel = %ld.\n"
                 "res = %hd.\n"
                 "arquivo fonte:  \"%s\".\n"
                 "nome da funcao: \"%s\".\n"
                 ,aux,j,__FILE__,__func__);
          exit(EXIT_FAILURE);
        }
        kk = nel-1;
        MAT2D(kk,j,id,maxRes) = res;   
      } 
      else{
        printf("erro: numero do elemento nao exitentes. nel = %ld.\n"
               "arquivo fonte:  \"%s\".\n"
               "nome da funcao: \"%s\".\n"
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
 * maxRes- numero maximo de carga por elemento                       *
 * str   - macro de terminada o fim da secao                         *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * f     - valores das restricoes                                    *
 *********************************************************************/
void readVfSource(double *f          ,INT numel
                 ,short int maxCarga,char *str
                 ,FILE* file){
  
  char word[WORD_SIZE];
  int  j,kk,nTerm;
  double carga;
  INT nel;
  long aux;  
  int error=0; 

  readMacro(file,word,false);
  do{
    nel = atol(word);  
    aux = (long) nel;
    error = fscanf(file,"%d",&nTerm);
    if( error != 1) {
      printf("erro: leitura do numero de termos de cargas. "
             "nel = %ld.\n"
             "arquivo fonte:  \"%s\".\n"
             "nome da funcao: \"%s\".\n"
             ,aux,__FILE__,__func__);
      exit(EXIT_FAILURE);

    }
    for(j = 0;j < nTerm;j++){
      if ( nel > 0 && nel <= numel){ 
        error = fscanf(file,"%lf",&carga);
        if( error != 1) {
          printf("erro: leitura da carg. "
                 "nel   = %ld.\n"
                 "carga = %hd.\n"
                 "arquivo fonte:  \"%s\".\n"
                 "nome da funcao: \"%s\".\n"
                 ,aux,j,__FILE__,__func__);
          exit(EXIT_FAILURE);
        }
        kk = nel-1;
        MAT2D(kk,j,f,maxCarga) = carga;   
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

/*********************************************************************/
/* Leitura dos materiais                                             */
/*********************************************************************/
void readVfMat(double *prop,short *type,short numat,FILE* file)
{
  
    short i;
    int k,dum;
    int nprop;
    char word[WORD_SIZE];
    int j;  
    char line[WORD_SIZE];
    
    
    
    for(i=0;i<numat;i++){
      fscanf(file,"%d",&k);
      fscanf(file,"%d",&dum);
      --k ;
      type[k] = dum;
      readMacro(file,line,true);
      nprop = getnumprop2(line);
      if( nprop > MAXPROP){
        printf("%s\n"
	       "*** Numero maximo de prorpiedades excedidos\n"
	       "MAXPROPMAX:                 %d\n"
	       "numero de propiedades lidas:%d\n"
	       "Nome do arquivo fonte %s.\n"
	       "Funcao %s.\n"
	       "%s\n"
	       ,DIF,MAXPROP,nprop,__FILE__,__func__,DIF);
	       exit(EXIT_FAILURE);       
      }
  
     for(j=0;j<nprop;j++){
       getword(line,word);
/*       printf("mat = %d prop = %d word = %s\n",k,j,word);*/
       prop[MAXPROP*k+j]=atof(word);
     }
      
   }
}
/*********************************************************************/

/********************************************************************* 
 * GETNUMPROP2: contar o numero de propriedades em um linha          * 
 * ----------------------------------------------------------------- * 
 * Parametros de entrada:                                            * 
 * line - string de caracter                                         * 
 * ----------------------------------------------------------------- * 
 * Parametros de saida:                                              * 
 *      - numero de valores na linha                                 * 
 * ----------------------------------------------------------------- * 
 *********************************************************************/
static int getnumprop2(char *line){

  int n,i;
  bool cont;
  char c;   
  
  n = 0;
  i = 0;
  cont = true;
  while( (c=line[i++]) != '\0'){
    
    if( c != ' ' && cont ){
      cont =false;
      n++;
    }  
    else if( c == ' ' )
      cont = true;
  }
  
  return n;

}
/*********************************************************************/

/********************************************************************* 
 * GETWORD: obtem o numero de um linha separado por espaco           * 
 * ----------------------------------------------------------------- * 
 * Parametros de entrada:                                            * 
 * line - string de caracter                                         * 
 * ----------------------------------------------------------------- * 
 * Parametros de saida:                                              * 
 * word - com o valor numerico                                       * 
 * ----------------------------------------------------------------- * 
 *********************************************************************/
static void getword(char *line, char*word){

    int n,i;
    char c;   
    bool flag;
    
    n = 0;
    i = 0;
    flag =false;
    while( (c=line[i]) != '\0'){
/*primeiro diferente de espaco*/      
      if(c != ' ')
        flag = true;
      if(c == ' ' && flag)
        break;
      word[n++]=c;
      line[i++]  =' ';  
    }
    word[n]='\0';  
    
}
/*********************************************************************/


