#include<Read_mef.h>
/*********************************************************************/
void read_mef(Memoria *m,Mesh *mesh,FILE* file)
{
  char word[WORD_SIZE];
  char macro[NMACROS][WORD_SIZE]={
            "mesh"           ,"coordinates"    ,"tria3"       
            ,"quad4"         ,"tetra4"         ,"hexa8"
            ,"tria6"         ,"quad8"          ,""      
            ,""              ,"barra2"         ,""      
	    ,"materials"     ,"constraindisp"  ,"nodalforces"
            ,""              ,"nodalloads"     ,"loads"
	    ,"constraintemp" ,"nodalsources"   ,"nodalthermloads"
	    ,"elmtthermloads","constrainconc"  ,"nodalconc"    
	    ,"nodalconcloads","elmtconcloadas" ,"velocityfield" 
	    ,"initialdisp"  ,"initialvel"      ,"initialconc"              
	    ,"initialtemp"  ,"insert"          ,"return"              
	    ,"end"
	   };
  bool rflag[NMACROS];			  
  long nn,nel;
  int nen,dim,numat,ndf,ndfd,ndfcd;
  bool macroflag;
  int i;
  FILE *faux=NULL;

/*...*/  
  for(i=0;i<NMACROS;i++){
    rflag[i] = false;
  }
  nmacro = 0;
  for(i=0;i<MAX_LINE;i++){
    strcpy(macros[i],"");
  }

  macroflag = true;
  while(macroflag && (!feof(file))){
    readmacro(file,word,false);
/*    printf("%d %s\n",nmacro,word);*/
/*    fprintf(stderr,"!%s!\n",word);*/
/*... macro mesh*/    
    if((!strcmp(word,macro[0])) && (!rflag[0])){
      printf("%s\n",DIF);
      printf("Mesh\n");
      strcpy(macros[nmacro++],word);
      rflag[0] = true;
      macroflag = true;
      parametros(&nn,&nel,&nen,&dim,&numat,&ndf,&ndfd,&ndfcd,file);
/*---alocando variaveis de elementos*/
/*... conectividade*/ 
      Myalloc(int,m,mesh->elm.node,nel*nen ,"elnode",_AD_);
/*... materiais*/ 
      Myalloc(int,m,mesh->elm.mat ,nel     ,"elmat" ,_AD_);
/*... nos por elementos*/
      Myalloc(int,m,mesh->elm.nen ,nel     ,"elnen" ,_AD_);
/*... tipo geometrico*/      
      Myalloc(int,m,mesh->elm.type ,nel     ,"geotype" ,_AD_);
/*... Prop*/ 
      Myalloc(double,m,mesh->material.prop,MAXPROP*numat     
             ,"prop" ,_AD_);
/*... tipo*/      
      Myalloc(short,m,mesh->material.type,numat     
             ,"type" ,_AD_);
/*...................................................................*/

/*...*/
      zero(mesh->elm.node     ,nel*nen      ,"int"   );
      zero(mesh->elm.mat      ,nel          ,"int"   );
      zero(mesh->elm.nen      ,nel          ,"int"   );
      zero(mesh->elm.type     ,nel          ,"int"   );
      zero(mesh->material.prop,MAXPROP*numat,"double");
      zero(mesh->material.type,numat        ,"short" );
/*...................................................................*/
     
/*---alocando variaveis nodais */      
/*---alocando coordenadas      */      
      Myalloc(double,m,mesh->node.x,dim*nn ,"xnode",_AD_);   
      
      zero(mesh->node.x,dim*nn,"double" );
/*...................................................................*/

/*... onda*/
      if(ndf){
/*---alocando deslocamentos    */      
        Myalloc(double,m,mesh->node.u,nn*ndf,"unode",_AD_);
/*---alocando velocidades      */      
        Myalloc(double,m,mesh->node.v,nn*ndf,"vnode",_AD_);   
/*---alocando aceleracoes     */      
        Myalloc(double,m,mesh->node.a,nn*ndf,"anode",_AD_);   
/*---alocando forcas           */      
        Myalloc(double,m,mesh->node.f,nn*ndf,"fnode",_AD_);   
/*---alocando restricoes       */      
        Myalloc(long,m,mesh->node.eqn,nn*ndf,"eqnode",_AD_);   
/*---alocando tipo cargas variaveis no tempo por no*/
        Myalloc(int,m,mesh->node.nloads,nn*ndf,"nloads",_AD_);   
/*---alocando tipo cargas variaveis no tempo por elementos*/
        Myalloc(int,m,mesh->elm.eloads,nel*ndf,"eloads",_AD_);   
/*---*/
        zero(mesh->node.u          ,nn*ndf       ,"double");
        zero(mesh->node.v          ,nn*ndf       ,"double");
        zero(mesh->node.a          ,nn*ndf       ,"double");
        zero(mesh->node.f          ,nn*ndf       ,"double");
        zero(mesh->node.eqn        ,nn*ndf       ,"long"  );
        zero(mesh->node.nloads     ,nn*ndf       ,"int"   );
        zero(mesh->elm.eloads     ,nel*ndf      ,"int"   );
/*...................................................................*/
      }
/*...................................................................*/

/*... conveccao-difusao*/
      else if(ndfcd){
/*---alocando concentracoes    */      
        Myalloc(double,m,mesh->node.u,nn*ndfcd,"utnode",_AD_);
/*---alocando velocidades      */      
        Myalloc(double,m,mesh->node.v,nn*ndfcd,"vtnode",_AD_);   
/*---alocando forcas           */      
        Myalloc(double,m,mesh->node.f,nn*ndfcd,"ftnode",_AD_);   
/*---alocando restricoes       */      
        Myalloc(long,m,mesh->node.eqn,nn*ndfcd,"eqtnode",_AD_);   
/*---alocando campo de velocidades*/	
        Myalloc(double,m,mesh->node.w,dim*nn ,"wnode",_AD_);   
/*---alocando tipo cargas variaveis no tempo*/
        Myalloc(int,m,mesh->node.nloads,nn*ndfcd,"ntloads",_AD_);   
/*---alocando tipo cargas variaveis no tempo por elementos*/
        Myalloc(int,m,mesh->elm.eloads,nel*MAXELOADS,"etloads",_AD_);   
/*---*/
        zero(mesh->node.u          ,nn*ndfcd      ,"double");
        zero(mesh->node.v          ,nn*ndfcd      ,"double");
        zero(mesh->node.f          ,nn*ndfcd      ,"double");
        zero(mesh->node.eqn        ,nn*ndfcd      ,"long"  );
        zero(mesh->node.w          ,dim*nn        ,"double");
        zero(mesh->node.nloads     ,nn*ndfcd      ,"int"   );
        zero(mesh->elm.eloads      ,nel*MAXELOADS ,"int"   );
/*...................................................................*/
      }
/*...................................................................*/

/*... difusao*/
      else if(ndfd){
/*---alocando concentracoes    */      
        Myalloc(double,m,mesh->node.u,nn*ndfd,"utnode",_AD_);
/*---alocando velocidades      */      
        Myalloc(double,m,mesh->node.v,nn*ndfd,"vtnode",_AD_);   
/*---alocando forcas           */      
        Myalloc(double,m,mesh->node.f,nn*ndfd,"ftnode",_AD_);   
/*---alocando restricoes       */      
        Myalloc(long,m,mesh->node.eqn,nn*ndfd,"eqtnode",_AD_);   
/*---alocando tipo cargas variaveis no tempo*/
        Myalloc(int,m,mesh->node.nloads,nn*ndfd,"ntloads",_AD_);   
/*---alocando tipo cargas variaveis no tempo por elementos*/
        Myalloc(int,m,mesh->elm.eloads,nel*MAXELOADS,"etloads",_AD_);   
/*---*/
        zero(mesh->node.u          ,nn*ndfd      ,"double");
        zero(mesh->node.v          ,nn*ndfd      ,"double");
        zero(mesh->node.f          ,nn*ndfd      ,"double");
        zero(mesh->node.eqn        ,nn*ndfd      ,"long"  );
        zero(mesh->node.nloads     ,nn*ndfd      ,"int"   );
        zero(mesh->elm.eloads      ,nel*MAXELOADS ,"int"   );
/*...................................................................*/
      }
/*...................................................................*/
      mesh->nnode     =  nn;
      mesh->numel     = nel;
      mesh->dim       = dim;
      mesh->ndf       = ndf;
      mesh->numat     = numat;
      mesh->maxno     = nen;
      mesh->ndfd      = ndfd;
      mesh->ndfcd     = ndfcd;
    }
/* ... coordenadas*/
    else if((!strcmp(word,macro[1]))&&(!rflag[1])){
      printf("loading coordinates...\n");
      strcpy(macros[nmacro++],word);
      rflag[1] = true;
      read_mef_coor(mesh,file);
      printf("load.\n");
    }  
/* ... tria3*/
    else if((!strcmp(word,macro[2]))&&(!rflag[2])){
      printf("loading tria3 ...\n");
      strcpy(macros[nmacro++],word);
      rflag[2] = true;
      read_mef_elm(mesh,nen,TTRIA3,file);
      printf("load.\n");
    }
/* ... quad4*/
    else if((!strcmp(word,macro[3]))&&(!rflag[3])){
      printf("loading quad4 ...\n");
      strcpy(macros[nmacro++],word);
      rflag[3] = true;
      read_mef_elm(mesh,nen,TQUAD4,file);
      printf("load.\n");
    }
/* ... tetra4*/
    else if((!strcmp(word,macro[4]))&&(!rflag[4])){
      printf("loading tetra4 ...\n");
      strcpy(macros[nmacro++],word);
      rflag[4] = true;
      read_mef_elm(mesh,nen,TTETRA4,file);
      printf("load.\n");
    }
/* ... hexa8 - Hexaedro de 8 nos*/
    else if((!strcmp(word,macro[5]))&&(!rflag[5])){
      printf("loading hexa8 ...\n");
      strcpy(macros[nmacro++],word);
      rflag[5] = true;
      read_mef_elm(mesh,nen,THEXA8,file);
      printf("load.\n");
    }
/* ... tria6 - Triangulo de 6 nos */
    else if((!strcmp(word,macro[6]))&&(!rflag[6])){
      printf("loading tria6 ...\n");
      strcpy(macros[nmacro++],word);
      rflag[6] = true;
      read_mef_elm(mesh,nen,TTRIA6,file);
      printf("load.\n");
    }
/* ... quad8 - Quadrilatero de 8 nos*/
    else if((!strcmp(word,macro[7]))&&(!rflag[7])){
      printf("loading quad8 ...\n");
      strcpy(macros[nmacro++],word);
      rflag[7] = true;
      read_mef_elm(mesh,nen,TQUAD8,file);
      printf("load.\n");
    }
/* ... barra2*/
    else if((!strcmp(word,macro[10]))&&(!rflag[10])){
      printf("loading barra2 ...\n");
      strcpy(macros[nmacro++],word);
      rflag[10] = true;
      read_mef_elm(mesh,nen,TBARRA2,file);
      printf("load.\n");
    }
/* ... materials*/
    else if((!strcmp(word,macro[12]))&&(!rflag[12])){
      printf("loading materials ...\n");
      strcpy(macros[nmacro++],word);
      rflag[12] = true;
      read_mef_mat(mesh->material.prop,mesh->material.type,numat,file);
      printf("load.\n");
    }
/* ... constraindisp*/
    else if((!strcmp(word,macro[13]))&&(!rflag[13])){
      printf("loading constraindisp ...\n");
      strcpy(macros[nmacro++],word);
      rflag[13] = true;
      read_mef_bound(mesh->node.eqn,ndf,file);
      printf("load.\n");
    }
/*.....................................................................*/

/* ...              */
    else if((!strcmp(word,macro[14]))&&(!rflag[14])){
      printf("loading nodalforces ...\n");
      strcpy(macros[nmacro++],word);
      rflag[14] = true;
      read_mef_forces(mesh->node.f,ndf,file);
      printf("load.\n");
    }
/*.....................................................................*/

/* ...              */
    else if((!strcmp(word,macro[16]))&&(!rflag[16])){
      printf("loading nodalloads ...\n");
      strcpy(macros[nmacro++],word);
      rflag[16] = true;
      read_mef_nloads(mesh->node.nloads,ndf,file);
      printf("load.\n");
    }
/*.....................................................................*/

/* ...              */
    else if((!strcmp(word,macro[17]))&&(!rflag[17])){
      printf("loading loads ...\n");
      strcpy(macros[nmacro++],word);
      rflag[17] = true;
      read_mef_loads(mesh->loads,file);
      printf("load.\n");
    }
/*.....................................................................*/

/* ...              */
    else if((!strcmp(word,macro[18]))&&(!rflag[18])){
      printf("loading constraintemp ... \n");
      strcpy(macros[nmacro++],word);
      rflag[18] = true;
      read_mef_bound(mesh->node.eqn,ndfd,file);
      printf("load.\n");
    }
/*.....................................................................*/

/* ...              */
    else if((!strcmp(word,macro[19]))&&(!rflag[19])){
      printf("loading nodalsources ...\n");
      strcpy(macros[nmacro++],word);
      rflag[19] = true;
      read_mef_forces(mesh->node.f,ndfd,file);
      printf("load. \n");
    }
/*.....................................................................*/

/* ... nodalthermloads*/
    else if((!strcmp(word,macro[20]))&&(!rflag[20])){
      printf("loading nodalthermloads ...\n");
      strcpy(macros[nmacro++],word);
      rflag[20] = true;
      read_mef_nloads(mesh->node.nloads,ndfd,file);
      printf("load. \n");
    }
/*.....................................................................*/

/* ... elmtthermloads*/
    else if((!strcmp(word,macro[21]))&&(!rflag[21])){
      printf("loading elmthermloads ...\n");
      strcpy(macros[nmacro++],word);
      rflag[21] = true;
      read_mef_eloads(mesh->elm.eloads,MAXELOADS,file);
      printf("load. \n");
    }
/*.....................................................................*/

/* ... concentracoes prescritas*/
    else if((!strcmp(word,macro[22]))&&(!rflag[22])){
      printf("loading constrainconc ... \n");
      strcpy(macros[nmacro++],word);
      rflag[22] = true;
      read_mef_bound(mesh->node.eqn,ndfcd,file);
      printf("load.\n");
    }
/*.....................................................................*/

/* ... valores de concentracoes prescritas ou fontes nodais*/
    else if((!strcmp(word,macro[23]))&&(!rflag[23])){
      printf("loading nodalconc ...\n");
      strcpy(macros[nmacro++],word);
      rflag[23] = true;
      read_mef_forces(mesh->node.f,ndfcd,file);
      printf("load. \n");
    }
/*.....................................................................*/

/* ... nodalconcloads*/
    else if((!strcmp(word,macro[24]))&&(!rflag[24])){
      printf("loading nodalconcloads ...\n");
      strcpy(macros[nmacro++],word);
      rflag[24] = true;
      read_mef_nloads(mesh->node.nloads,ndfcd,file);
      printf("load. \n");
    }
/*.....................................................................*/

/* ... elmtconcloads*/
    else if((!strcmp(word,macro[25]))&&(!rflag[25])){
      printf("loading elmtconcloads ...\n");
      strcpy(macros[nmacro++],word);
      rflag[25] = true;
      read_mef_eloads(mesh->elm.eloads,MAXELOADS,file);
      printf("load. \n");
    }
/*.....................................................................*/

/* ... campo de velocidades*/
    else if((!strcmp(word,macro[26]))&&(!rflag[26])){
      printf("loading velocityfield ...\n");
      strcpy(macros[nmacro++],word);
      rflag[26] = true;
      read_mef_forces(mesh->node.w,dim,file);
      printf("load. \n");
    }
/*.....................................................................*/

/* ...              */
    else if((!strcmp(word,macro[27]))&&(!rflag[27])){
      printf("loading initialdisp ...\n");
      strcpy(macros[nmacro++],word);
      rflag[27] = true;
      read_mef_forces(mesh->node.u,ndf,file);
      printf("load. \n");
    }
/*.....................................................................*/

/* ...  velacidades iniciais (v)             */
    else if((!strcmp(word,macro[28]))&&(!rflag[28])){
      printf("loading initialvel ...\n");
      strcpy(macros[nmacro++],word);
      rflag[28] = true;
      read_mef_forces(mesh->node.v,ndf,file);
      printf("load. \n");
    }
/*.....................................................................*/

/* ...  velacidades iniciais (v)             */
    else if((!strcmp(word,macro[29]))&&(!rflag[29])){
      printf("loading initialconc ...\n");
      strcpy(macros[nmacro++],word);
      rflag[29] = true;
      read_mef_forces(mesh->node.u,ndfcd,file);
      printf("load. \n");
    }
/*.....................................................................*/
/*... temperaturas iniciais (u)*/    
    else if((!strcmp(word,macro[30]))&&(!rflag[30])){
      printf("loading initialtemp ...\n");
      strcpy(macros[nmacro++],word);
      rflag[30] = true;
      read_mef_forces(mesh->node.u,ndfd,file);
      printf("load. \n");
    }
/*.....................................................................*/

/*...arquivo secundario*/    
    else if((!strcmp(word,macro[31]))&&(!rflag[31])){
      strcpy(macros[nmacro++],word);
      rflag[31] = true;
      rflag[32] = false;
      readmacro(file,word,false);
      printf("insert file: %s ...\n",word);
      faux = file; 
      file = open_file(word,"r");
    }
/*.....................................................................*/

/*...retornando ao arquivo principal*/    
    else if((!strcmp(word,macro[32]))&&(!rflag[32])){
      printf("return .\n");
      strcpy(macros[nmacro++],word);
      rflag[31] = false;
      rflag[32] = true;
      file = faux;
    }
/*.....................................................................*/

/*... end mesh*/
    else if(!strcmp(word,macro[33])){
/*      printf("end\n");*/
      readmacro(file,word,false);
      if(!strcmp(word,"mesh")){
/*        printf("mesh\n");*/
        strcpy(macros[nmacro++],"end mesh");
        macroflag = false;
      }
    }
  }
/*... Macros comandos apos o end mesh*/
  i = 0;
  readmacro(file,word,true);
  while(!feof(file) && i < MAX_LINE && strcmp(word,"stop")){
    readmacro(file,word,true);
/*    printf("%d!%s!\n",nmacro,word);*/
    strcpy(macros[nmacro++],word);
    i++;
  }

/*...*/
  if(mesh->ndf) 
    boundc(mesh->nnode,mesh->ndf,mesh->node.eqn,mesh->node.f
          ,mesh->node.u);
  if(mesh->ndfd) 
    boundc(mesh->nnode,mesh->ndfd,mesh->node.eqn,mesh->node.f
          ,mesh->node.u);
  if(mesh->ndfcd) 
    boundc(mesh->nnode,mesh->ndfcd,mesh->node.eqn,mesh->node.f
          ,mesh->node.u);
  
/*...................................................................*/  

  printf("%s\n",DIF);
/*...................................................................*/  
}
/*********************************************************************/

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
 * ndfd  - graus de liberdade problemas difusivos                    *
 * ndfcd - graus de liberdade problemas difusivos-convectivo         *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 *********************************************************************/
void parametros(long *nn,long  *nel,int *nen,int *dim
               ,int *numat,int *ndf,int *therm,int *ndfcd,FILE* file)
{
  char parameter[][WORD_SIZE]={"nnode","numel","numat"
                              ,"maxno","ndf"  ,"therm"
			      ,"ndfcd","dim"};
  char word[WORD_SIZE];
  bool flag[NPARAMETROS];
  int i=0,j;


  *nn    = 0;
  *nel   = 0;
  *numat = 0;
  *nen   = 0;
  *ndf   = 0;
  *therm = 0;
  *ndfcd = 0;
  *dim   = 0;

  for(j=0;j<NPARAMETROS;j++)
    flag[j] = false;
  while(i < NPARAMETROS && i < 20){
    readmacro(file,word,false);
/*... macro mesh*/    
    if(!strcmp(word,parameter[0])){
      fscanf(file,"%ld",nn);
#ifdef _DEBUG_MESH_ 
      printf("nnode %ld\n",*nn);
#endif      
      flag[0] = true;
      i+=1;
    }
    else if(!strcmp(word,parameter[1])){
      fscanf(file,"%ld",nel);
#ifdef _DEBUG_MESH_ 
      printf("numel %ld\n",*nel);
#endif      
      flag[1] = true;
      i+=1;
    }
    else if(!strcmp(word,parameter[2])){
      fscanf(file,"%d",numat);
#ifdef _DEBUG_MESH_ 
      printf("numat %d\n",*numat);
#endif      
      flag[2] = true;
      i+=1;
    }
    else if(!strcmp(word,parameter[3])){
      fscanf(file,"%d",nen);
#ifdef _DEBUG_MESH_ 
      printf("maxno %d\n",*nen);
#endif      
      flag[3] = true;
      i+=1;
    }
    else if(!strcmp(word,parameter[4])){
      fscanf(file,"%d",ndf);
#ifdef _DEBUG_MESH_ 
      printf("ndf %d\n",*ndf);
#endif      
      flag[4] = true;
      i+=1;
    }
    else if(!strcmp(word,parameter[5])){
      fscanf(file,"%d",therm);
#ifdef _DEBUG_MESH_ 
      printf("therm %d\n",*therm);
#endif      
      flag[5] = true;
      i+=1;
    }
    else if(!strcmp(word,parameter[6])){
      fscanf(file,"%d",ndfcd);
#ifdef _DEBUG_MESH_ 
      printf("ndfcd %d\n",*ndfcd);
#endif      
      flag[6] = true;
      i+=1;
    }
    else if(!strcmp(word,parameter[7])){
      fscanf(file,"%d",dim);
#ifdef _DEBUG_MESH_ 
      printf("dim %d\n",*dim);
#endif      
      flag[7] = true;
      i+=1;
    }
    else
      i+=1;
  }
  
  for(j=0;j<NPARAMETROS;j++){
    if(!flag[j]){
      fprintf(stderr,"parametro: %s faltando.\n"
             "fonte: %s \n",parameter[j],SRC_NAME);
      exit(EXIT_FAILURE);
    }
  }
}
/*********************************************************************/

/*********************************************************************
 * BOUNDC: valores prescritos                                        *
 *********************************************************************
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * nnode - numero de nos                                             *
 * ndf   - graus de liberades                                        *
 * eqn   - numero maximo de nos por elemento                         *
 * f     - forcas ou u precritos                                     *
 * u     - vetor solucao                                             *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 *********************************************************************/
  void boundc(int nnode,int ndf,long*eqn,double*f,double *u){

   int i;

    for(i=0;i<nnode*ndf;i++)
/* se for diferente de zero, u prescrito*/    
      if(eqn[i])
        u[i]=f[i];
  }

/*********************************************************************
 * READ_MEF_COOR:                                                    *
 *********************************************************************/
void read_mef_coor(Mesh *mesh,FILE *file){
  long int i,idum,k;
  long int nn;
  int dim;
  int j;
 
   nn = mesh->nnode;
  dim = mesh->dim;
  
  for(i=0;i<nn;i++){
    fscanf(file,"%ld",&idum);
    for(j=0;j<dim;j++){
      k = (idum-1)*dim + j;
      fscanf(file,"%lf",&(mesh->node.x[k]));
    }
  }
#ifdef _DEBUG_MESH_ 
  for(i=0;i<nn;i++){
    fprintf(stderr,"%ld",i+1);
    for(j=0;j<dim;j++){
      k = i*dim + j;
      fprintf(stderr," %lf ",mesh->node.x[k]);
    }
    printf("\n");
  }
#endif
}
/*********************************************************************/

/*********************************************************************/
void read_mef_elm(Mesh *mesh,int nen,int type,FILE *file){
  
  long int i,k,idum;
  int j;
  long int ne;
   
  ne = mesh->numel;
  
  for(i=0;i<ne;i++){
    fscanf(file,"%ld",&idum);
    for(j=0;j<nen;j++){
      k = nen*(idum-1)+j;
      fscanf(file,"%d",&(mesh->elm.node[k]));
    }
    fscanf(file,"%d",&(mesh->elm.mat[i]));
    mesh->elm.nen[i]   = nen;
    mesh->elm.type[i]  = type;
  }
  
#ifdef _DEBUG_MESH_ 
  for(i=0;i<mesh->numel;i++){
    fprintf(stderr,"%ld",i+1);
    for(j=0;j<nen;j++){
      k = nen*i+j;
      fprintf(stderr," %d ",mesh->elm.node[k]);
    }  
    printf("\n");
  }
#endif
}
/*********************************************************************/
/* Leitura das restricoes:                                           */
/*********************************************************************/
void read_mef_bound(long *eqn, int ndf,FILE* file)
{
    char word[WORD_SIZE]; 
    int j,k;
    int no;
    
    readmacro(file,word,false);
    while(strcmp(word,"end")){
      no = atol(word);
      for(j = 0;j < ndf ;j ++){
        k = (no-1)*ndf + j;
        fscanf(file,"%ld",&eqn[k]);
      }
      readmacro(file,word,false);
    }
  
}
/*********************************************************************/

/*********************************************************************/
/* Leitura das forcas por no:                                        */
/*********************************************************************/
void read_mef_forces(double *f, int ndf,FILE* file)
{
  
    char word[WORD_SIZE];
    int j,k;
    int no;
    
    
    readmacro(file,word,false);
    while(strcmp(word,"end")){
      no = atol(word);
      for(j = 0;j < ndf ;j ++){
        k = (no-1)*ndf + j;
        fscanf(file,"%lf",&f[k]);
      }
      readmacro(file,word,false);
    }
    
}
/*********************************************************************/

/*********************************************************************/
/* Leitura dos materiais                                             */
/*********************************************************************/
void read_mef_mat(double *prop,short *type,int numat,FILE* file)
{
  
    int i,k,dum;
    int nprop;
    char word[WORD_SIZE];
    int j;  
    char line[WORD_SIZE];
    
    
    
    for(i=0;i<numat;i++){
      fscanf(file,"%d",&k);
      fscanf(file,"%d",&dum);
      --k ;
      type[k]=dum;
      readmacro(file,line,true);
      nprop = getnumprop2(line);
      if( nprop > MAXPROP){
        printf("%s\n"
	       "*** Numero maximo de prorpiedades excedidos\n"
	       "MAXPROPMAX:                 %d\n"
	       "numero de propiedades lidas:%d\n"
	       "Name do arquivo fonte %s.\n"
	       "Funcao void read_mef_mat(double*,int,FILE*)\n"
	       "%s\n"
	       ,DIF,MAXPROP,nprop,SRC_NAME,DIF);
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

/*********************************************************************/
/* Leitura das defincoes de cargas por no                            */
/*********************************************************************/
void read_mef_nloads(int *nloads,int ndf,FILE* file)
{
  
    char word[WORD_SIZE];
    int j,k;
    int no;
    
    
    readmacro(file,word,false);
    while(strcmp(word,"end")){
      no = atol(word);
      for(j = 0;j < ndf ;j ++){
        k = (no-1)*ndf + j;
        fscanf(file,"%d",&nloads[k]);
      }
      readmacro(file,word,false);
    }


}
/*********************************************************************/

/*********************************************************************/
/* Leitura das defincoes de cargas por elemento                      */
/*********************************************************************/
void read_mef_eloads(int *eloads,int maxeloads,FILE* file)
{
  
    char line[WORD_SIZE];
    char word[WORD_SIZE];
    int j;
    int el;
    int ncarg;
    
    readmacro(file,word,false);
    while(strcmp(word,"end")){
      el= atol(word);
      el--;
/*le o restante da linha*/      
      readmacro(file,line,true);
/*conta o numero de cargas nessa linha*/
      ncarg= getnumprop2(line);
      if( ncarg > maxeloads){
        printf("%s\n"
        "*** Numero maximo de cargas excedidos\n"
	       "MAXELOADS:                  %d\n"
	       "numero de caragas:          %d\n"
	       "Name do arquivo fonte %s.\n"
	       "Funcao void read_mef_mat(double*,int,FILE*)\n"
	       "%s\n"
	      ,DIF,maxeloads,ncarg,SRC_NAME,DIF);
	       exit(EXIT_FAILURE);       
      }
  
     for(j=0;j<ncarg;j++){
/*obtem os valores das cargas na linha lida*/     
       getword(line,word);
/*       printf("el = %d carg = %d word = %s\n",el,j,word);*/
       eloads[maxeloads*el+j]=atof(word);
     }
    readmacro(file,word,false);
   }  

}
/*********************************************************************/
/*********************************************************************/
/* Leitura das defincoes das forcas                                  */
/*********************************************************************/
void read_mef_loads(Loads *loads,FILE* file)
{
  
    char word[WORD_SIZE];
    int i;
    int nload,nparc,type;
    
    readmacro(file,word,false);
    while(strcmp(word,"end")){
      nload = atol(word);
      nload --;
      fscanf(file,"%d",&type);
      fscanf(file,"%d",&nparc);
      loads[nload].type=type;
      loads[nload].nparc=nparc;
      if(nparc >  MAXPARCLOAD){
        printf("***Numero maximo de parcelas no load excedido\n");
        exit(EXIT_FAILURE);
      }	
      for(i=0;i<nparc;i++)
        fscanf(file,"%lf",&loads[nload].parc[i]);
      readmacro(file,word,false);
    }
}
/*********************************************************************/

/********************************************************************* 
 * GETNUMPROP : contar o numero de propriedades em um linha          * 
 * ----------------------------------------------------------------- * 
 * Parametros de entrada:                                            * 
 * file - arquivo de dados                                           * 
 * ----------------------------------------------------------------- * 
 * Parametros de saida:                                              * 
 *      - numero de valores na linha                                 * 
 * ----------------------------------------------------------------- * 
 * OBS:retorna a posicao da arquivo original                         *
 * (problemas entre diferentes OS                                    *  
 *********************************************************************/
int getnumprop(FILE *file){

  int n,i;
  char word[WORD_SIZE],c;
  bool cont;
  long int posicao;

/*... obtendo a posicao do arquivo*/
  posicao=ftell(file);
/*...................................................................*/  
  readmacro(file,word,true);
  
  n = 0;
  i = 0;
  cont = true;
  while( (c=word[i++]) != '\0'){
    
    if( c != ' ' && cont ){
      cont =false;
      n++;
    }  
    else if( c == ' ' )
      cont = true;
  }
  
/*... retornando a posicao*/
  fseek(file,posicao,0);
/*...................................................................*/  
  
  return n;

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
int getnumprop2(char *line){

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
void getword(char *line, char*word){

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
