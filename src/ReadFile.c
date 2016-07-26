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
            "coordinates","endMesh"   ,"insert"       /* 0, 1, 2*/
           ,"return"     ,"cells"     ,"faceRt1"      /* 3, 4, 5*/
           ,"faceLoadT1" ,"loadsT1"   ,""             /* 6, 7, 8*/ 
           ,""           ,""          ,""             /* 9,10,11*/ 
           ,"faceRd1"    ,""          ,"loadsD1"      /*12,13,14*/ 
           ,"faceLoadD1" ,""          ,""             /*15,16,17*/ 
           ,"faceRvel"   ,"loadsVel"  ,"faceLoadVel"  /*18,19,20*/ 
           ,"faceRpres"  ,"loadsPres" ,"faceLoadPres" /*21,22,23*/ 
           ,"materials"  ,""          ,"initialVel"   /*24,25,26*/ 
	   };                                             
  bool rflag[NMACROS],macroFlag;
  INT nn,nel;
  short maxno,ndm,numat,maxViz,ndfVel;
  char nameAux[MAX_STR_LEN_IN];
  FILE *fileAux=NULL;
  int i;

/*... leitura dos parametros principais da malha*/
  parametros(&nn   ,&nel
            ,&maxno,&maxViz
            ,&ndm  ,&numat
            ,file);
/*...................................................................*/

/*...*/
  mesh->nnode    = nn;
  mesh->nnodeNov = nn;
  mesh->nnodeOv  =  0;
  mesh->numel    = nel;
  mesh->numelNov = nel;
  mesh->maxNo    = maxno;
  mesh->maxViz   = maxViz; 
  mesh->ndm      = ndm;
  mesh->numat    = numat;

/*... alocando variavies de elementos*/
/*... conectividade*/ 
  HccaAlloc(INT,m,mesh->elm.node       ,nel*maxno,"elnode"  ,_AD_);
/*... materiais*/ 
  HccaAlloc(short,m,mesh->elm.mat      ,nel      ,"elmat"   ,_AD_);
/*... nos por elementos*/
  HccaAlloc(short,m,mesh->elm.nen      ,nel      ,"elnen"   ,_AD_);
/*... tipo geometrico */
  HccaAlloc(short,m,mesh->elm.geomType ,nel      ,"elgT"    ,_AD_);

/*... centroide */
  HccaAlloc(DOUBLE,m,mesh->elm.geom.cc ,nel*ndm          ,"elCc"    ,_AD_);

  if( mpiVar.nPrcs < 2){

/*... vetor que une os centroides dos elementos */
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.ksi
             ,nel*ndm*maxViz,"elksi"  ,_AD_);
/*... modulo do vetor que une os centroides dos elementos */
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.mksi
              ,nel*maxViz     ,"elmksi",_AD_);
/*... vetor paralelo a face da celula */
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.eta
             ,nel*ndm*maxViz,"eleta"  ,_AD_);
/*... modulo do vetor paralelo a face da celula */
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.fArea
            ,nel*maxViz     ,"elfArea",_AD_);
/*... volume da celula*/                           
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.volume
            ,nel            ,"elVol",_AD_);
/*... vetor normal a face da celula*/                           
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.normal
             ,nel*maxViz*ndm       ,"elnorm",_AD_);
/*... ponto medio da face*/                           
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.xm
             ,nel*maxViz*ndm       ,"elxm",_AD_);
/*... vetor que une o centroide ao ponto medio*/                           
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.xmcc
             ,nel*maxViz*ndm       ,"elxmcc",_AD_);
/*... vetor entre o ponto medio da face e ponto de intercao 
      entre a linha entre os centroides e a face*/         
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.vSkew  
             ,nel*ndm*maxViz         ,"elvSkew" ,_AD_);
/*... distancia entro o ponto medio da face e ponto de intercao 
      entre a linha entre os centroides e a face*/         
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.mvSkew  
             ,nel*maxViz           ,"elmvSkew" ,_AD_);
/*... distancia entro o ponto medio da face e ponto de intercao 
      entre a linha entre os centroides e a face*/         
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.dcca 
             ,nel*maxViz           ,"eldcca",_AD_);

/*... zerando os variavies*/
    zero(mesh->elm.node         ,nel*maxno     ,INTC);
    zero(mesh->elm.mat          ,nel           ,"short"  );
    zero(mesh->elm.nen          ,nel           ,"short"  );
    zero(mesh->elm.geomType     ,nel           ,"short"  );
    zero(mesh->elm.geom.cc      ,nel*ndm       ,DOUBLEC);
    zero(mesh->elm.geom.ksi     ,nel*ndm*maxViz,DOUBLEC);
    zero(mesh->elm.geom.mksi    ,nel*maxViz    ,DOUBLEC);
    zero(mesh->elm.geom.eta     ,nel*ndm*maxViz,DOUBLEC);
    zero(mesh->elm.geom.fArea   ,nel*maxViz    ,DOUBLEC);
    zero(mesh->elm.geom.volume  ,nel           ,DOUBLEC);
    zero(mesh->elm.geom.normal  ,nel*ndm*maxViz,DOUBLEC);
    zero(mesh->elm.geom.xm      ,nel*ndm*maxViz,DOUBLEC);
    zero(mesh->elm.geom.xmcc    ,nel*ndm*maxViz,DOUBLEC);
    zero(mesh->elm.geom.mvSkew  ,nel*maxViz    ,DOUBLEC);
    zero(mesh->elm.geom.mvSkew  ,nel*ndm*maxViz,DOUBLEC);
    zero(mesh->elm.geom.dcca    ,nel*maxViz    ,DOUBLEC);
/*...................................................................*/
  }
/*...................................................................*/

/*... alocando materiais*/
/*... Prop*/ 
  HccaAlloc(DOUBLE,m,mesh->elm.material.prop,MAXPROP*numat     
         ,"prop" ,_AD_);
/*... type*/ 
  HccaAlloc(short,m,mesh->elm.material.type,numat     
         ,"type" ,_AD_);
/*... zerando os variavies*/
  zero(mesh->elm.material.prop,MAXPROP*numat,DOUBLEC);
  zero(mesh->elm.material.type,numat,"short");
/*...................................................................*/

/*... alocando estruturas para vizinhos*/
/*... nelcon*/ 
  HccaAlloc(INT,m,mesh->elm.adj.nelcon,nel*maxViz ,"adj" ,_AD_);
/*... type*/ 
  HccaAlloc(short,m,mesh->elm.adj.nViz,nel       ,"nViz" ,_AD_);
/*... zerando os variavies*/
  zero(mesh->elm.adj.nelcon,nel*maxViz,INTC);
  zero(mesh->elm.adj.nViz  ,nel       ,"short");
/*...................................................................*/

/*---alocando variaveis nodais */      
/*---alocando coordenadas      */      
  HccaAlloc(DOUBLE,m,mesh->node.x,ndm*nn ,"xnode",_AD_);   
     
/*... zerando os variavies*/
  zero(mesh->node.x,ndm*nn,DOUBLEC);
/*...................................................................*/

/*... transporte e fluido*/
  if(mesh->ndfT[0] > 0 || mesh->ndfF > 0) {     
/*... eVel*/
    HccaAlloc(DOUBLE,m,mesh->elm.vel 
             ,nel*mesh->ndm    ,"eVel"              ,_AD_);
    zero(mesh->elm.vel       ,nel*mesh->ndm         ,DOUBLEC);
/*...................................................................*/

/*... eVel0*/
    HccaAlloc(DOUBLE,m,mesh->elm.vel0 
             ,nel*mesh->ndm    ,"eVel0"             ,_AD_);
    zero(mesh->elm.vel0      ,nel*mesh->ndm         ,DOUBLEC);
/*...................................................................*/

/*... nVel*/
     HccaAlloc(DOUBLE,m,mesh->node.vel 
              ,nn*mesh->ndm     ,"nVel"              ,_AD_);
     zero(mesh->node.vel      ,nn*mesh->ndm          ,DOUBLEC);
/*...................................................................*/

  }
/*...................................................................*/

/*... fluido*/
  if(mesh->ndfF > 0) {     
/*... alocando memoria*/
     HccaAlloc(short,m,mesh->elm.faceRvel  
            ,nel*(maxViz+1),"faceRvel"    ,_AD_);
     zero(mesh->elm.faceRvel  ,nel*(maxViz+1),"short"  );
     
     HccaAlloc(short,m,mesh->elm.faceLoadVel  
            ,nel*(maxViz+1),"faceLVel"    ,_AD_);
     zero(mesh->elm.faceLoadVel  ,nel*(maxViz+1),"short"  );
     
     HccaAlloc(short,m,mesh->elm.faceRpres 
            ,nel*(maxViz+1),"faceRpres"   ,_AD_);
     zero(mesh->elm.faceRpres ,nel*(maxViz+1),"short"  );
     
     HccaAlloc(short,m,mesh->elm.faceLoadPres 
            ,nel*(maxViz+1),"faceLPres"   ,_AD_);
     zero(mesh->elm.faceLoadPres ,nel*(maxViz+1),"short"  );

/*... ePres*/
     HccaAlloc(DOUBLE,m,mesh->elm.pressure 
            ,nel              ,"pressure"          ,_AD_);
     zero(mesh->elm.pressure  ,nel                         ,DOUBLEC);
     
/*... densityFluid*/ 
     HccaAlloc(DOUBLE,m,mesh->elm.densityFluid
              ,nel*2            ,"denFluid" ,_AD_);
     zero(mesh->elm.densityFluid ,nel*2                       ,DOUBLEC);

/*... pres*/
     HccaAlloc(DOUBLE,m,mesh->node.pressure
              ,nn    ,"npressure",_AD_);
     zero(mesh->node.pressure     ,nn                         ,DOUBLEC);

/*... nGradVel*/
     ndfVel = mesh->ndfF-1;
     HccaAlloc(DOUBLE,m,mesh->node.gradVel  
              ,nn*ndm*ndfVel ,"nGradVel"     ,_AD_);
     zero(mesh->node.gradVel  ,nn*ndm*ndfVel        ,DOUBLEC);
     
/*... eGradVel*/
     HccaAlloc(DOUBLE,m,mesh->elm.gradVel 
              ,nel*ndm*ndfVel ,"eGradVel"     ,_AD_);
     zero(mesh->elm.gradVel   ,nel*ndm*ndfVel       ,DOUBLEC);

/*... nGradPres*/
     HccaAlloc(DOUBLE,m,mesh->node.gradPres 
              ,nn*ndm        ,"nGradPres"    ,_AD_);
     zero(mesh->node.gradPres ,nn*ndm               ,DOUBLEC);
     
/*... eGradPresl*/
     HccaAlloc(DOUBLE,m,mesh->elm.gradPres
              ,nel*ndm        ,"eGradPres"     ,_AD_);
     zero(mesh->elm.gradPres  ,nel*ndm              ,DOUBLEC);

     if( mpiVar.nPrcs < 2){
       ndfVel = mesh->ndfF-1;
/*... rCellVel*/
       HccaAlloc(DOUBLE,m,mesh->elm.rCellVel  
                ,nel*ndm*ndfVel       ,"rCellVel"     ,_AD_);
       zero(mesh->elm.rCellVel  ,nel*ndm*ndfVel   ,DOUBLEC);
/*... rCellPres*/
       HccaAlloc(DOUBLE,m,mesh->elm.rCellPres 
                ,nel                 ,"rCellPres"    ,_AD_);
       zero(mesh->elm.rCellPres ,nel              ,DOUBLEC);
     }
/*...................................................................*/
  }
/*...................................................................*/

/*... problema de transporte*/   
  if(mesh->ndfT[0] > 0) {     
/*... alocando memoria*/
     HccaAlloc(short,m,mesh->elm.faceRt1
            ,nel*(maxViz+1),"faceRt1"  ,_AD_);
     zero(mesh->elm.faceRt1   ,nel*(maxViz+1),"short"  );
     
     HccaAlloc(short,m,mesh->elm.faceLoadT1
            ,nel*(maxViz+1),"faceLt1"  ,_AD_);
     zero(mesh->elm.faceLoadT1,nel*(maxViz+1),"short"  );
/*...................................................................*/

/*... eUt1*/
     HccaAlloc(DOUBLE,m,mesh->elm.uT1 
            ,nel*mesh->ndfT[0],"eUt1"              ,_AD_);
     zero(mesh->elm.uT1       ,nel*mesh->ndfT[0]           ,DOUBLEC);
     
/*... eU0t1*/
     HccaAlloc(DOUBLE,m,mesh->elm.u0T1
              ,nel*mesh->ndfT[0],"eU0t1"             ,_AD_);
     zero(mesh->elm.u0T1        ,nel*mesh->ndfT[0]           ,DOUBLEC);
/*... densityUt1*/ 
     HccaAlloc(DOUBLE,m,mesh->elm.densityUt1
              ,nel*2            ,"densityUt1" ,_AD_);
     zero(mesh->elm.densityUt1  ,nel*2                       ,DOUBLEC);

/*... uT1*/
     HccaAlloc(DOUBLE,m,mesh->node.uT1 
              ,nn*mesh->ndfT[0] ,"nUt1"              ,_AD_);
     zero(mesh->node.uT1      ,nn*mesh->ndfT[0]            ,DOUBLEC);

/*... nGradUt1*/
     HccaAlloc(DOUBLE,m,mesh->node.gradUt1  
              ,nn*ndm*mesh->ndfT[0] ,"nGradUt1"     ,_AD_);
     zero(mesh->node.gradUt1  ,nn*ndm*mesh->ndfT[0]        ,DOUBLEC);
     
/*... eGradU1*/
     HccaAlloc(DOUBLE,m,mesh->elm.gradUt1 
              ,nel*ndm*mesh->ndfT[0],"eGradUt1"     ,_AD_);
     zero(mesh->elm.gradUt1   ,nel*ndm*mesh->ndfT[0]       ,DOUBLEC);
      
     if( mpiVar.nPrcs < 2){
/*... rCell*/
       HccaAlloc(DOUBLE,m,mesh->elm.rCellUt1  
                ,nel*ndm*mesh->ndfT[0],"rCellUt1"     ,_AD_);
       zero(mesh->elm.rCellUt1  ,nel*mesh->ndfT[0]       ,DOUBLEC);
     }
/*...................................................................*/
   }
/*...................................................................*/
   
/*... problema de difusao pura*/   
  if(mesh->ndfD[0] > 0) {   
/*... alocando memoria*/
     HccaAlloc(short,m,mesh->elm.faceRd1
            ,nel*(maxViz+1),"faceRd1"  ,_AD_);
     zero(mesh->elm.faceRd1   ,nel*(maxViz+1),"short"  );
     
     HccaAlloc(short,m,mesh->elm.faceLoadD1
            ,nel*(maxViz+1),"faceLd1"  ,_AD_);
     zero(mesh->elm.faceLoadD1,nel*(maxViz+1),"short"  );
     
/*... eUd1*/
     HccaAlloc(DOUBLE,m,mesh->elm.uD1 
            ,nel*mesh->ndfD[0],"eUd1"              ,_AD_);
     zero(mesh->elm.uD1       ,nel*mesh->ndfD[0]           ,DOUBLEC);
     
/*... eU0d1*/
     HccaAlloc(DOUBLE,m,mesh->elm.u0D1
              ,nel*mesh->ndfD[0],"eU0d1"             ,_AD_);
     zero(mesh->elm.u0D1        ,nel*mesh->ndfD[0]           ,DOUBLEC);

/*... densityUd1*/ 
     HccaAlloc(DOUBLE,m,mesh->elm.densityUd1
              ,nel*2            ,"densityUd1" ,_AD_);
     zero(mesh->elm.densityUd1  ,nel*2                       ,DOUBLEC);

/*... uD1*/
     HccaAlloc(DOUBLE,m,mesh->node.uD1 
              ,nn*mesh->ndfD[0] ,"nUd1"              ,_AD_);
     zero(mesh->node.uD1      ,nn*mesh->ndfD[0]            ,DOUBLEC);

/*... nGradU1*/
     HccaAlloc(DOUBLE,m,mesh->node.gradUd1  
              ,nn*ndm*mesh->ndfD[0] ,"nGradUd1"     ,_AD_);
     zero(mesh->node.gradUd1  ,nn*ndm*mesh->ndfD[0]        ,DOUBLEC);
     
/*... eGradU1*/
     HccaAlloc(DOUBLE,m,mesh->elm.gradUd1 
              ,nel*ndm*mesh->ndfD[0],"eGradUd1"     ,_AD_);
     zero(mesh->elm.gradUd1   ,nel*ndm*mesh->ndfD[0]       ,DOUBLEC);
      
     if( mpiVar.nPrcs < 2){
/*... rCell*/
       HccaAlloc(DOUBLE,m,mesh->elm.rCellUd1  
                ,nel*ndm*mesh->ndfD[0],"rCellUd1"     ,_AD_);
       zero(mesh->elm.rCellUd1  ,nel*mesh->ndfD[0]       ,DOUBLEC);
     }
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
/*  printf("%s\n",word);*/
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

/*... insert*/
    else if((!strcmp(word,macro[2])) && (!rflag[2])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      printf("%s\n\n",DIF);
      fscanf(file,"%s",nameAux);
      fileAux = file;
      file    = openFile(nameAux,"r");
      rflag[2] = true;
    }
/*...................................................................*/

/*... return*/
    else if((!strcmp(word,macro[3]))){
      printf("%s\n",DIF);
      printf("%s\n",word);
      printf("%s\n\n",DIF);
      if(!rflag[2]){
        ERRO_GERAL(__FILE__,__func__,__LINE__
                  ,"Erro: macro return sem um insert associado!!");
      }
      fclose(file);
      file = fileAux;
      rflag[2] = false;
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
      readVfRes(mesh->elm.faceRt1,mesh->numel,mesh->maxViz+1,str,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceLoadT1 - cargas nas faces transporte */
    else if((!strcmp(word,macro[6])) && (!rflag[6])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[6] = true;
      strcpy(str,"endFaceLoadT1");
      printf("loading faceLoadT1 ...\n");
      readVfRes(mesh->elm.faceLoadT1,mesh->numel
               ,mesh->maxViz+1       ,str       ,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... loadT1 - definicao de cargar transporte */
    else if((!strcmp(word,macro[7])) && (!rflag[7])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[7] = true;
      strcpy(str,"endLoadsT1");
      printf("loading loadsT1 ...\n");
      readVfLoads(loadsT1,str       ,file);
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

/*... loadD1 - definicao de cargar difusao pura */
    else if((!strcmp(word,macro[14])) && (!rflag[14])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[14] = true;
      strcpy(str,"endLoadsD1");
      printf("loading loadsD1 ...\n");
      readVfLoads(loadsD1,str       ,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceLoadD1 - cargas nas faces difusao pura */
    else if((!strcmp(word,macro[15])) && (!rflag[15])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[15] = true;
      strcpy(str,"endFaceLoadD1");
      printf("loading faceLoadD1 ...\n");
      readVfRes(mesh->elm.faceLoadD1,mesh->numel
               ,mesh->maxViz+1       ,str       ,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceRvel- condicao de contorno para problemas fluidos (Vel) */
    else if((!strcmp(word,macro[18])) && (!rflag[18])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[18] = true;
      strcpy(str,"endFaceRvel");
      printf("loading faceRvel ...\n");
      readVfRes(mesh->elm.faceRvel,mesh->numel,mesh->maxViz+1,str,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... loadVel - definicao de cargar fluidos (Vel)*/
    else if((!strcmp(word,macro[19])) && (!rflag[19])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[19] = true;
      strcpy(str,"endLoadsVel");
      printf("loading loadsVel ...\n");
      readVfLoads(loadsVel,str       ,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceLoadVel - cargas nas faces fluido (Vel)*/
    else if((!strcmp(word,macro[20])) && (!rflag[20])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[20] = true;
      strcpy(str,"endFaceLoadVel");
      printf("loading faceLoadVel ...\n");
      readVfRes(mesh->elm.faceLoadVel  ,mesh->numel
               ,mesh->maxViz+1         ,str        ,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceRpres - condicao de contorno para problemas fluidos (Pres)*/
    else if((!strcmp(word,macro[21])) && (!rflag[21])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[21] = true;
      strcpy(str,"endFaceRpres");
      printf("loading faceRpres ...\n");
      readVfRes(mesh->elm.faceRpres,mesh->numel
               ,mesh->maxViz+1     ,str        ,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... loadPres - definicao de cargar fluidos (Pres)*/
    else if((!strcmp(word,macro[22])) && (!rflag[22])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[22] = true;
      strcpy(str,"endLoadsPres");
      printf("loading loadsPres ...\n");
      readVfLoads(loadsPres,str       ,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceLoadPres - cargas nas faces fluido (Pres)*/
    else if((!strcmp(word,macro[23])) && (!rflag[23])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[23] = true;
      strcpy(str,"endFaceLoadPres");
      printf("loading faceLoadPres ...\n");
      readVfRes(mesh->elm.faceLoadPres ,mesh->numel
               ,mesh->maxViz+1         ,str       ,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... materiais */
    else if((!strcmp(word,macro[24])) && (!rflag[24])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[24] = true;
      strcpy(str,"endMaterials");
      printf("loading materials ...\n");
      readVfMat(mesh->elm.material.prop,mesh->elm.material.type
               ,numat,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/

/*... initialVel */
    else if((!strcmp(word,macro[26])) && (!rflag[26])){
      printf("%s\n",DIF);
      printf("%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[25] = true;
      strcpy(str,"endInitialVel");
      printf("loading initialVel ...\n");
      readVfInitial(mesh->elm.vel,mesh->numel,mesh->ndm,str,file);
      printf("load.\n");
      printf("%s\n\n",DIF);
    }
/*...................................................................*/
  }while(macroFlag && (!feof(file)));
  
/*... iniciacao de propriedades do material varaivel com o tempo*/
  if(mesh->ndfD[0] > 0) {   
    initProp(mesh->elm.densityUd1
            ,mesh->elm.material.prop,mesh->elm.mat
            ,2                      ,mesh->numel
            ,DENSITY);         
  }
/*...................................................................*/

/*...*/
  if(mesh->ndfT[0] > 0) {   
    initProp(mesh->elm.densityUt1
            ,mesh->elm.material.prop,mesh->elm.mat
            ,2                      ,mesh->numel
            ,DENSITY);         
  }
/*...................................................................*/

/*...*/
  if(mesh->ndfF > 0) {   
    initProp(mesh->elm.densityFluid
            ,mesh->elm.material.prop,mesh->elm.mat
            ,2                      ,mesh->numel
            ,DENSITY);         
  }
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
 * maxNo   -> numero de nos por celula maximo da malha               * 
 * maxViz  -> numero vizinhos por celula maximo da malha             * 
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 *********************************************************************/
void parametros(INT  *nn    ,INT *nel    
               ,short *maxNo,short *maxViz
               ,short *ndm  ,short *numat
               ,FILE* file)
{
  char parameter[][WORD_SIZE]={"nnode","numel","numat"
                              ,"maxno","ndm"  ,"maxviz"
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

/*... ndm*/
    else if(!strcmp(word,parameter[4])){
      fscanf(file,"%hd",ndm);
#ifdef _DEBUG_MESH_ 
      printf("ndm %hd\n",*ndm);
#endif      
      flag[4] = true;
      i++;
    }
/*...................................................................*/

/*... macro maxViz*/  
    else if(!strcmp(word,parameter[5])){
      fscanf(file,"%hd",maxViz);
#ifdef _DEBUG_MESH_ 
      printf("maxno %hd\n",*maxViz);
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
void readVfCoor(DOUBLE *x,INT nn, short ndm,FILE *file){
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
void readVfSource(DOUBLE *f          ,INT numel
                 ,short const maxCarga,char *str
                 ,FILE* file){
  
  char word[WORD_SIZE];
  int  j,kk,nTerm;
  DOUBLE carga;
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
          printf("erro: leitura da carg.\n"
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
/*********************************************************************/

/*********************************************************************
 * READVFINITIAL: leitura dos valores iniciais                       *
 *********************************************************************
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * f     - indefinido                                                *
 * numel - numero de elementos                                       *
 * gdl   - graus de liberdade                                        *
 * str   - macro de terminada o fim da secao                         *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * f     - valores iniciais                                          *
 *********************************************************************/
void readVfInitial(DOUBLE *f          ,INT numel
                  ,short const ndf     ,char *str
                  ,FILE* file){
  
  char word[WORD_SIZE];
  int  j,kk;
  DOUBLE carga;
  INT nel;
  long aux;  
  int error=0; 

  readMacro(file,word,false);
  do{
    nel = atol(word);  
    aux = (long) nel;
    for(j = 0;j < ndf;j++){
      if ( nel > 0 && nel <= numel){ 
        error = fscanf(file,"%lf",&carga);
        if( error != 1) {
          printf("erro: leitura da carg.\n"
                 "nel   = %ld.\n"
                 "ndf   = %hd.\n"
                 "arquivo fonte:  \"%s\".\n"
                 "nome da funcao: \"%s\".\n"
                 ,aux,j,__FILE__,__func__);
          exit(EXIT_FAILURE);
        }
        kk = nel-1;
        MAT2D(kk,j,f,ndf) = carga;   
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
 * READVFLOADS : leitura da definicoes das cargas                    *
 *********************************************************************
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * load  - indefinido                                                *
 * str   - macro de terminada o fim da secao                         *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * f     - valores das restricoes                                    *
 *********************************************************************/
void readVfLoads(Loads *loads,char *str,FILE* file){
  
  char word[WORD_SIZE];
  int  j,nTerm,nLoad,type;
  DOUBLE par;
  int error=0; 

  readMacro(file,word,false);
  do{
    nLoad = atol(word);  
/*...*/
    error = fscanf(file,"%d",&type);
    if( error != 1) {
      printf("erro: leitura do tipo da carga. "
             "nLoad = %d.\n"
             "arquivo fonte:  \"%s\".\n"
             "nome da funcao: \"%s\".\n"
             ,nLoad,__FILE__,__func__);
      exit(EXIT_FAILURE);

    }
    loads[nLoad-1].type = type;

/*...*/
    error = fscanf(file,"%d",&nTerm);
    if( error != 1) {
      printf("erro: leitura do numero de termos de cargas. "
             "nLoad = %d.\n"
             "arquivo fonte:  \"%s\".\n"
             "nome da funcao: \"%s\".\n"
             ,nLoad,__FILE__,__func__);
      exit(EXIT_FAILURE);

    }
    loads[nLoad-1].np = nTerm;

    if(MAXLOADPARAMETER < nTerm){
      printf("erro: Numeroa de parametros excediso. "
             "MAX  = %d.\n"
             "nPar = %d.\n"
             "arquivo fonte:  \"%s\".\n"
             "nome da funcao: \"%s\".\n"
             ,MAXLOADPARAMETER,nTerm
             ,__FILE__        ,__func__);
      exit(EXIT_FAILURE);
    }

/*...*/
    for(j = 0;j < nTerm;j++){
      error = fscanf(file,"%lf",&par);
      if( error != 1) {
        printf("erro: leitura da carga. "
               "nLoad = %d.\n"
               "nTerm = %d.\n"
               "arquivo fonte:  \"%s\".\n"
               "nome da funcao: \"%s\".\n"
               ,nLoad,j,__FILE__,__func__);
          exit(EXIT_FAILURE);
      }
      loads[nLoad-1].par[j] = par;
    } 
    readMacro(file,word,false);
  }while(strcmp(word,str));
}
/*********************************************************************/

/*********************************************************************
 * CONFIG : configuraceos gerais                                     *
 *********************************************************************
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * fileOpt   - opcoes de arquivo                                     *
 * reordMesh - reordenacao do malha                                  *
 * rcGrad    - tipo de tecnica de rescontrucao de gradiente          *
 * m         - memoria principal                                     *
 * file  - ponteiro para o arquivo de dados                          *
 * ----------------------------------------------------------------- *
 *********************************************************************/
void config(FileOpt *opt,Reord *reordMesh
           ,short *rcGrad
           ,FILE* file)
{
  char config[][WORD_SIZE]={"bvtk"  ,"reord"     ,"mem" 
                           ,"rcGrad","fItPlotRes","fItPlot"};
  
  char word[WORD_SIZE];
  char s[WORD_SIZE];
  bool flag[NCONFIG];
  int i=0,j,temp;
  int conv;


  for(j=0;j<NCONFIG;j++)
    flag[j] = false;

  while(i < NCONFIG && i < 20){
    readMacro(file,word,false);
/*... bvtk*/   
    if(!strcmp(word,config[0])){
      fscanf(file,"%s",s);
      if(!strcmp(s,"true"))
        opt->bVtk = true;
      else
        opt->bVtk = false;
      flag[0] = true;
      i++;
      if(!mpiVar.myId){
        if(opt->bVtk)
          printf("bVtk: true\n");
        else
          printf("bVtk: false\n");
      }
    }
/*... reord*/   
    else if(!strcmp(word,config[1])){
      fscanf(file,"%s",s);
      if(!strcmp(s,"true"))
        reordMesh -> flag= true;
      else
        reordMesh -> flag= false;
          
      flag[1] = true;
      i++;
      if(!mpiVar.myId){
        if(reordMesh->flag)
          printf("Reord: true\n");
        else
          printf("Reod: false\n");
      }
    }
/*... mem*/   
    else if(!strcmp(word,config[2])){
      fscanf(file,"%d",&temp);
//    conv    = CONV_BYTES*CONV_BYTES;
      conv    = 1024*1024;
 			nmax    = (iptx) temp;
			nmax    =  nmax*conv;
      flag[2] = true;
      i++;
      
      if(!mpiVar.myId)
        printf("Memoria principal: %d MBytes\n"
              ,(int)(nmax/conv));
    }
/*... rcGrad*/   
    else if(!strcmp(word,config[3])){
      fscanf(file,"%s",s);
      if(!strcmp(s,"gglc")){
        *rcGrad = RCGRADGAUSSC;
        if(!mpiVar.myId)
          printf("rcGrad: GreenGaussCell\n");
      }
      else if(!strcmp(s,"ggln")){
        *rcGrad = RCGRADGAUSSN;
        if(!mpiVar.myId)
          printf("rcGrad: GreenGaussNode\n");
      }
      else if(!strcmp(s,"lSquare")){
        *rcGrad = RCLSQUARE;
        if(!mpiVar.myId)
          printf("rcGrad: LeastSquare\n");
      }
      
      else if(!strcmp(s,"lSquareQR")){
        *rcGrad = RCLSQUAREQR;
        if(!mpiVar.myId)
          printf("rcGrad: LeastSquareQR\n");
      }

      flag[3] = true;
      i++;
    }
/*... fItPlotRes*/   
    else if(!strcmp(word,config[4])){
      fscanf(file,"%s",s);
      if(!strcmp(s,"true")){
        opt->fItPlotRes = true;
        if(!mpiVar.myId)
          printf("fItPlotRes: true\n");
      }
      else{
        opt->fItPlotRes = false;
        if(!mpiVar.myId)
          printf("fItPlotRes: false\n");
      }
      flag[4] = true;
      i++;
    }
/*... fItPlot*/   
    else if(!strcmp(word,config[5])){
      fscanf(file,"%s",s);
      if(!strcmp(s,"true")){
        opt->fItPlot = true;
        if(!mpiVar.myId)
          printf("fItPlot: true\n");
      }
      else{
        opt->fItPlot = false;
        if(!mpiVar.myId)
          printf("fItPlot: false\n");
      }
      flag[5] = true;
      i++;
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
void readVfMat(DOUBLE *prop,short *type,short numat,FILE* file)
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

/********************************************************************* 
 * INITPROP: inicializao de propriedades com variacao temporal       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * prop    -> nao definido                                           * 
 * propMat -> propriedade de referencia por material                 * 
 * mat     -> material por celula                                    * 
 * np      -> numero de propriedades                                 * 
 * nCell   -> numero de celulas                                      * 
 * iProp   -> numero da propriedade                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * prop    -> propriedade iniciacializada                            * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void initProp(DOUBLE *restrict prop 
             ,DOUBLE *restrict propMat,short *restrict mat
             ,short const np          ,INT    const nCell 
             ,short const iProp)
{    
  INT i;
  unsigned short j,lMat;         
  for(i=0;i<nCell;i++){    
    lMat               = mat[i]-1;
    for(j=0;j<np;j++){
      MAT2D(i,j,prop,np) = MAT2D(lMat,iProp,propMat,MAXPROP);
    }
  }
}
/*********************************************************************/
      
/********************************************************************* 
 * READEDP : graus de liberdade das equacoes diferencias             * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * mesh    ->                                                        * 
 * file    -> arquivo de arquivo                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * mesh    ->                                                        * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void readEdo(Mesh *mesh,FILE *file){

  char *str={"endEdp"};
  char word[WORD_SIZE];
	short i;

/*...*/
	for (i = 0; i < MAX_TRANS_EQ; i++)
		mesh->ndfT[i] = 0;
	for (i = 0; i < MAX_DIF_EQ; i++)
		mesh->ndfD[i] = 0;
	mesh->ndfF = 0;
	mesh->ndfFt = 0;
/*...................................................................*/

  readMacro(file,word,false);
  do{
/*... transport1*/
    if(!strcmp(word,"transport1")){
      readMacro(file,word,false);
      mesh->ndfT[0]=(short) atol(word);
      if(!mpiVar.myId ) printf("transport1 ndf %d\n"
                              ,mesh->ndfT[0]);
    }
/*...................................................................*/

/*... transport2*/
    else if(!strcmp(word,"transport2")){
      readMacro(file,word,false);
      mesh->ndfT[1]= (short) atol(word);
      if(!mpiVar.myId ) printf("transport2 ndf %d\n"
                               ,mesh->ndfT[1]);
    }
/*...................................................................*/

/*... transport3*/
    else if(!strcmp(word,"transport3")){
      readMacro(file,word,false);
      mesh->ndfT[2]= (short) atol(word);
      if(!mpiVar.myId ) printf("transport3 ndf %d\n"
                               ,mesh->ndfT[2]);
    }
/*...................................................................*/

/*... diffusion1*/
    else if(!strcmp(word,"diffusion1")){
      readMacro(file,word,false);
      mesh->ndfD[0]= (short) atol(word);
      if(!mpiVar.myId ) printf("diffusion1 ndf %d\n"
                              ,mesh->ndfD[0]);
    }
/*...................................................................*/

/*... diffusion2*/
    else if(!strcmp(word,"diffusion2")){
      readMacro(file,word,false);
      mesh->ndfD[1]= (short) atol(word);
      if(!mpiVar.myId ) printf("diffusion2 ndf %d\n"
                              ,mesh->ndfD[1]);
    }
/*...................................................................*/

/*... diffusion3*/
    else if(!strcmp(word,"diffusion3")){
      readMacro(file,word,false);
      mesh->ndfD[2]= (short) atol(word);
      if(!mpiVar.myId ) printf("diffusion3 ndf %d\n"
                              ,mesh->ndfD[2]);
    }
/*...................................................................*/

/*... fluido*/
    else if(!strcmp(word,"fluid")){
      readMacro(file,word,false);
      mesh->ndfF= (short) atol(word);
      if(!mpiVar.myId ) printf("fluid ndf %d\n"
                              ,mesh->ndfF);
    }
/*...................................................................*/

/*... fluido-termo ativado*/
    else if(!strcmp(word,"fluidt")){
      readMacro(file,word,false);
      mesh->ndfFt= (short) atol(word);
      if(!mpiVar.myId ) printf("fluidt ndf %d\n"
                              ,mesh->ndfFt);
    }
/*...................................................................*/
    
    readMacro(file,word,false);
  }while(strcmp(word,str));


}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 17/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * SETPPRINTFLUID : Seleciona as veriaves que serao impressas na     *
 * macro pFluid                                                      *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * opt     -> opcoes de arquivos de saida                            * 
 * file    -> arquivo de entrada                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * opt     -> opcoes de arquivos de saida atualizados                * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void setPrintFluid(FileOpt *opt,FILE *file){

  char str[]={"end"};
  char word[WORD_SIZE];
  
  readMacro(file,word,false);
  while(strcmp(word,str)){
/*... fPrintMesh*/        
    if(!strcmp(word,"vel")){ 
      opt->vel = true;
      if(!mpiVar.myId ) printf("print : vel\n");
    }
/*.....................................................................*/

/*...*/
    else if(!strcmp(word,"pres")){ 
      opt->pres = true;
      if(!mpiVar.myId ) printf("print : pres\n");
    }
/*.....................................................................*/

/*...*/
    else if(!strcmp(word,"gradVel")){ 
      opt->gradVel = true;
      if(!mpiVar.myId ) printf("print : gradVel\n");
    }
/*.....................................................................*/

/*...*/
    else if(!strcmp(word,"gradPres")){ 
      opt->gradPres = true;
      if(!mpiVar.myId ) printf("print : gradPres\n");
    }
/*.....................................................................*/
    readMacro(file,word,false);
  }
/*.....................................................................*/

} 
/*********************************************************************/ 
