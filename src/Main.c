/********************* includes padroes ******************************/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
/*********************************************************************/

/*********************************************************************/
#include<Coo.h>
#include<File.h>
#include<Memoria.h>
#include<Mesh.h>
#include<ReadFile.h>
#include<WriteVtk.h>
#include<Sisteq.h>
#include<Reord.h>
/*********************************************************************/


int main(int argc,char**argv){

/*...Memoria principal*/  
  Memoria m;
/*... estrutura da dados para a malha*/
  Mesh *mesh=NULL;
/*... Sistema de equacao*/
  SistEq *sistEqT1=NULL;
/*... reordenacao da malha*/
  Reord  *reordMesh=NULL;
/*...*/
  long i;
/*Estrutura de dados*/
  char strIa[MNOMEPONTEIRO],strJa[MNOMEPONTEIRO];
  char strA[MNOMEPONTEIRO],strAd[MNOMEPONTEIRO];

/*... arquivo*/
  char *nameIn=NULL,*nameOut=NULL,*preName=NULL;
  FILE *fileIn=NULL,*fileOut=NULL;
  bool bvtk=true;

/* ... macro camandos de leitura*/
  bool macroFlag; 
  char word[WORD_SIZE],str[WORD_SIZE];
  char macro[][WORD_SIZE] = {"mesh","stop","config"
                            ,"pgeo","pcoo"
                            }; 
/* ..................................................................*/

/*... abrindo ar quivo de entrada*/ 
  reordMesh = (Reord*) malloc(sizeof(Reord));
  if(reordMesh == NULL){
    printf("Erro ponteiro reord\n");
    exit(EXIT_FAILURE);
  }
  reordMesh->flag = false;  
/* ..................................................................*/

/*... abrindo ar quivo de entrada*/ 
  mesh = (Mesh*) malloc(sizeof(Mesh));
  if(mesh == NULL){
    printf("Erro ponteiro mesh\n");
    exit(EXIT_FAILURE);
  }  
    
  nameIn = (char *) malloc(sizeof(char)*MAX_STR_LEN_IN);
    
  if( argc > 1)
    strcpy(nameIn,argv[1]);
  else{
    printf("Arquivo de dados:\n");
    scanf("%s",nameIn);
  }

  fileIn = openFile(nameIn,"r");
 
/*...................................................................*/

/*... arquivos de saida*/
  preName = (char *) malloc(sizeof(char)*MAX_STR_LEN_SUFIXO);
  if(preName == NULL){
    printf("Erro ponteiro prename\n");
    exit(EXIT_FAILURE);
  }
    
  nameOut = (char *) malloc(sizeof(char)*(SIZEMAX));
  if(nameOut == NULL){
    printf("Erro ponteiro nameout\n");
    exit(EXIT_FAILURE);
  }
    
  if( argc > 2)
    strcpy(preName,argv[2]);
  else{
    printf("Prefixo do arquivo de saida:\n");
    scanf("%s",preName);
  }
/*...................................................................*/
  
/*loop de execucao*/
  macroFlag = true;
  do{
/*... leitura da macro*/
  readMacro(fileIn,word,false);
/*...................................................................*/

/*===================================================================*
 * macro: mesh
 *===================================================================*/
  if((!strcmp(word,macro[0]))){
    printf("%s\n",DIF);
    printf("%s\n",word); 
    printf("%s\n",DIF);
/*... sistema de memoria*/
    initMem(&m,true);
/*... leitura da malha*/
    readFileFv(&m,mesh,fileIn);
/*...................................................................*/

/*... calcula a vizinhaca do elementos*/
    viz(&m             ,mesh->elm.node ,mesh->adj.nelcon
       ,mesh->adj.nViz ,mesh->elm.nen  ,mesh->nnode
       ,mesh->numel    ,mesh->maxNo    ,mesh->maxViz);
/*...................................................................*/

/*...*/
    sistEqT1 = (SistEq*) malloc(sizeof(SistEq));
    if(sistEqT1 == NULL){
      printf("Erro ponteiro sistEqT1\n");
      exit(EXIT_FAILURE);
    }
    sistEqT1->storage = 1;
    sistEqT1->unsym   = false;    
/*...................................................................*/

/*... reodenando as celulas para dimuincao da banda*/
    Myalloc(long,&m,reordMesh->num,mesh->numel,"rNum" ,_AD_);
    reord(reordMesh->num,mesh->adj.nelcon,mesh->numel,reordMesh->flag);
/*...................................................................*/

/*... numeracao das equacoes*/
    Myalloc(long,&m,sistEqT1->id,mesh->numel,"sistT1id",_AD_);
    sistEqT1->neq = numeq(&m,sistEqT1->id  ,reordMesh->num
                         ,mesh->elm.faceRt1,mesh->elm.nen
                         ,mesh->numel      ,mesh->maxViz
                         ,str);
    for(i=0;i<mesh->numel;i++)
      printf("%3ld %3ld %3ld\n",i+1,sistEqT1->id[i],reordMesh->num[i]); 
/*...................................................................*/

/*...*/
    strcpy(strIa,"iaT1");
    strcpy(strJa,"JaT1");
    strcpy(strAd,"aDT1");
    strcpy(strA ,"aT1");
    dataStruct(&m,sistEqT1->id,reordMesh->num,mesh->adj.nelcon
              ,mesh->adj.nViz,mesh->numel,mesh->maxViz,mesh->ndfT[0]
              ,strIa,strJa,strAd,strA
              ,sistEqT1);
    for(i=0;i<=sistEqT1->neq;i++)
      printf("%3ld %3ld\n",i+1,sistEqT1->ia[i]); 
    for(i=0;i<sistEqT1->nad;i++)
      printf("%3ld %3ld\n",i+1,sistEqT1->ja[i]); 
/*...................................................................*/  
    strcpy(str,"B");
    memoriaTotal(str);
    usoMemoria(&m,str);
    mapVector(&m);
  }   
/*===================================================================*/

/*===================================================================*
 * macro: stop
 *===================================================================*/
  else if((!strcmp(word,macro[1]))){
    printf("%s\n",DIF);
    printf("%s\n",word); 
    printf("%s\n\n",DIF);
    finalizeMem(&m,true);
    macroFlag = false;
  }   
/*===================================================================*/

/*===================================================================*
 * macro: config
 *===================================================================*/
  else if((!strcmp(word,macro[2]))){
    printf("%s\n",DIF);
    printf("%s\n",word);
    config(&bvtk,fileIn);
    printf("%s\n\n",DIF);
  }   
/*===================================================================*/

/*===================================================================*
 * macro: pgeo
 *===================================================================*/
  else if((!strcmp(word,macro[3]))){
    printf("%s\n",DIF);
    printf("%s\n",word);
    fName(preName,0,6,&nameOut);
    wGeoVtk(&m               ,mesh->node.x     ,mesh->elm.node 
           ,mesh->elm.mat    ,mesh->elm.nen    ,mesh->elm.geomType
           ,mesh->elm.faceRt1,mesh->elm.faceSt1,mesh->nnode   
           ,mesh->numel      ,mesh->ndm        ,mesh->maxNo  
           ,mesh->ndfT       ,nameOut          ,bvtk             
           ,fileOut);
    printf("%s\n\n",DIF);
  }   
/*===================================================================*/
   
/*===================================================================*
 * macro: pcoo
 *===================================================================*/
  else if((!strcmp(word,macro[4]))){
    printf("%s\n",DIF);
    printf("%s\n",word);
    fName(preName,0,12,&nameOut);
/*...*/
    writeCoo(&m,sistEqT1->ia,sistEqT1->ja,sistEqT1->neq
            ,sistEqT1->nad  ,sistEqT1->storage
            ,sistEqT1->unsym,nameOut);
/*...................................................................*/
    printf("%s\n\n",DIF);
  }   
/*===================================================================*/
  }while(macroFlag && (!feof(fileIn)));

 
  return EXIT_SUCCESS;
}    
/*********************************************************************/
 
 
