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
void readFileFv(Memoria *m,Mesh *mesh, FILE* file)
{
  char word[WORD_SIZE],str[WORD_SIZE];
  char macro[NMACROS][WORD_SIZE]={
            "coordinates","endMesh"  ,"insert"            /*0,1,2*/
           ,"return"     ,"cells"    ,"faceRt1"           /*3,4,5*/
           ,"faceSt1"    
	   };                                             
  bool rflag[NMACROS],macroFlag;
  long int nn,nel;
  short maxno,ndm,numat,ndf;
  int i;

/* leitura dos parametros principais da malha*/
  parametros(&nn,&nel,&maxno,&ndm,&numat,&ndf,file);
  mesh->nnode  = nn;
  mesh->numel  = nel;
  mesh->maxNo  = maxno;
  mesh->maxViz = maxno;
  mesh->ndm    = ndm;
  mesh->numat  = numat;
  mesh->ndfT[0] = ndf;

/*... alocando variavies de elementos*/
/*... conectividade*/ 
  Myalloc(long,m,mesh->elm.node     ,nel*maxno ,"elnode",_AD_);
/*... materiais*/ 
  Myalloc(short,m,mesh->elm.mat     ,nel       ,"elmat" ,_AD_);
/*... nos por elementos*/
  Myalloc(short,m,mesh->elm.nen     ,nel       ,"elnen" ,_AD_);
/*... tipo geometrico */
  Myalloc(short,m,mesh->elm.geomType,nel       ,"elgT"  ,_AD_);

/*... zerando os variavies*/
  zero(mesh->elm.node     ,nel*maxno    ,"long"  );
  zero(mesh->elm.mat      ,nel          ,"short" );
  zero(mesh->elm.nen      ,nel          ,"short" );
  zero(mesh->elm.geomType ,nel          ,"short" );
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
  Myalloc(long,m,mesh->adj.nelcon,nel*maxno,"adj" ,_AD_);
/*... type*/ 
  Myalloc(short,m,mesh->adj.nViz,nel       ,"nViz" ,_AD_);
/*... zerando os variavies*/
  zero(mesh->adj.nelcon,nel*maxno,"long");
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
     Myalloc(short,m,mesh->elm.faceRt1,nel*(maxno+1),"faceRt1"  ,_AD_);
     Myalloc(double,m,mesh->elm.faceSt1
            ,nel*(maxno+1)*ndf,"faceSt1"  ,_AD_);
     zero(mesh->elm.faceRt1  ,nel*(maxno+1)              ,"short"  );
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
      readVfElmt(mesh->elm.node    ,mesh->elm.mat,mesh->elm.nen 
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
      readVfRes(mesh->elm.faceRt1,mesh->numel,mesh->maxNo,str,file);
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
      readVfSource(mesh->elm.faceSt1,mesh->numel,mesh->maxNo
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
void parametros(long  int *nn   ,long  int *nel  ,short int *maxno
               ,short int *ndm  ,short int *numat,short int *ndf
               ,FILE* file)
{
  char parameter[][WORD_SIZE]={"nnode","numel","numat"
                              ,"maxno","ndf"  ,"ndm"};
  
  char word[WORD_SIZE];
  bool flag[NPARAMETROS];
  int i=0,j;


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
      fscanf(file,"%ld",nn);
#ifdef _DEBUG_MESH_ 
      printf("nnode %ld\n",*nn);
#endif      
      flag[0] = true;
      i+=1;
    }
/*...................................................................*/

/*... macro numel*/   
    else if(!strcmp(word,parameter[1])){
      fscanf(file,"%ld",nel);
#ifdef _DEBUG_MESH_ 
      printf("numel %ld\n",*nel);
#endif      
      flag[1] = true;
      i+=1;
    }
/*... macro numat*/   
    else if(!strcmp(word,parameter[2])){
      fscanf(file,"%hd",numat);
#ifdef _DEBUG_MESH_ 
      printf("numat %hd\n",*numat);
#endif      
      flag[2] = true;
      i+=1;
    }
/*...................................................................*/

/*... macro maxno*/   
    else if(!strcmp(word,parameter[3])){
      fscanf(file,"%hd",maxno);
#ifdef _DEBUG_MESH_ 
      printf("maxno %hd\n",*maxno);
#endif      
      flag[3] = true;
      i+=1;
    }
/*...................................................................*/

/*... macro ndf*/   
    else if(!strcmp(word,parameter[4])){
      fscanf(file,"%hd",ndf);
#ifdef _DEBUG_MESH_ 
      printf("ndf %hd\n",*ndf);
#endif      
      flag[4] = true;
      i+=1;
    }
/*...................................................................*/

/*... ndm*/
    else if(!strcmp(word,parameter[5])){
      fscanf(file,"%hd",ndm);
#ifdef _DEBUG_MESH_ 
      printf("ndm %hd\n",*ndm);
#endif      
      flag[5] = true;
      i+=1;
    }
/*...................................................................*/

    else
      i+=1;
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
void readVfCoor(double *x,long int nn, short int ndm,FILE *file){
  long int i,idum,k;
  int j;
 
  for(i=0;i<nn;i++){
    fscanf(file,"%ld",&idum);
    k = idum -1;
    for(j=0;j<ndm;j++){
      fscanf(file,"%lf",&COOR(k,j,x,ndm));
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
 * ty   -> numero do tipo geometrico do elemento                     *
 * ------------------------------------------------------------------*
 *********************************************************************/
void readVfElmt(long int *el   ,short int *mat  ,short int *nen 
               ,short int *ty  ,long int nel    ,short int maxno
               ,FILE *file){
  long int i,idum;
  short int nenl;
  int j;
 
  for(i=0;i<nel;i++){
    fscanf(file,"%ld",&idum);
    idum--;
    fscanf(file,"%hd",&mat[idum]);
    fscanf(file,"%hd",&ty[idum]);
    fscanf(file,"%hd",&nenl);
    nen[idum] = nenl;
    for(j=0;j<nenl;j++){
      fscanf(file,"%ld",&ELM(idum,j,el,maxno));
    }
  }
#ifdef _DEBUG_MESH_ 
  for(i=0;i<nel;i++){
    fprintf(stderr,"%ld",i+1);
    for(j=0;j<ndm;j++){
      k = i*ndm + j;
      fprintf(stderr," %lf ",el[k]);
    }
    printf("\n");
  }
#endif
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
void readVfRes(short int *id,long int numel,short int maxno
              ,char *str    ,FILE* file){
  
  char word[WORD_SIZE];
  int  j,k,kk,nTerm;
  long nel;  

  readMacro(file,word,false);
  do{
    nel = atol(word);  
    fscanf(file,"%d",&nTerm);
    kk = (nel-1)*(maxno+1);
    for(j = 0;j < nTerm;j++){
      k = kk + j;
      if ( nel > 0 && nel <= numel){ 
        fscanf(file,"%hd",&id[k]);
      } 
      else{
        printf("Erro: numero do elemento nao exitentes. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
              ,nel,__FILE__,"readVfSource");
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
void readVfSource(double *f,long numel, short int maxno,char *str
                 ,FILE* file){
  
  char word[WORD_SIZE];
  int  j,k,kk,nTerm;
  long nel;  

  readMacro(file,word,false);
  do{
    nel = atol(word);  
    fscanf(file,"%d",&nTerm);
    kk = (nel-1)*(maxno+1);
    for(j = 0;j < nTerm;j++){
      k = kk + j;
      if ( nel > 0 && nel <= numel){ 
        fscanf(file,"%lf",&f[k]);
      }  
      else{
        printf("Erro: numero do elemento nao exitentes. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
              ,nel,__FILE__,"readVfSource");
        exit(EXIT_FAILURE);
      }
    }
  readMacro(file,word,false);
  }while(strcmp(word,str));
}


/*********************************************************************
 * CONFIG : cinfiguracao gerais                                      *
 *********************************************************************
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * bvtk  - arquivo binario para o vtk                                *
 * nel   - numero de elementos                                       *
 * nen   - numero maximo de nos por elemento                         *
 * numat - numero de materias                                        *
 * ndf   - graus de liberdade (onda)                                 *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 *********************************************************************/
void config(bool *bvtk,FILE* file)
{
  char config[][WORD_SIZE]={"bvtk"};
  
  char word[WORD_SIZE];
  bool flag[NCONFIG];
  int i=0,j,temp;


  *bvtk  = 0;

  for(j=0;j<NCONFIG;j++)
    flag[j] = false;

  while(i < NCONFIG && i < 20){
    readMacro(file,word,false);
/*... macro nnode*/   
    if(!strcmp(word,config[0])){
      fscanf(file,"%d",&temp);
      *bvtk = (bool) temp;
      flag[0] = true;
      i+=1;
    }
/*...................................................................*/

    else
      i+=1;
  }
  
  for(j=0;j<NCONFIG;j++){
    if(!flag[j]){
      fprintf(stderr,"config: %s faltando.\n"
             "fonte: %s \n",config[j],__FILE__);
      exit(EXIT_FAILURE);
    }
  }
}
/*********************************************************************/
