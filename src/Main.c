/********************* includes padroes ******************************/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
/*********************************************************************/

/*********************************************************************/
#include<File.h>
#include<Memoria.h>
#include<Mesh.h>
#include<ReadFile.h>
#include<WriteVtk.h>
/*********************************************************************/


int main(int argc,char**argv){

/*...Memoria principal*/  
  Memoria m;
/*... estrutura da dados para a malha*/
  Mesh *mesh=NULL;
/*... arquivo*/
  char *nameIn=NULL,*nameOut=NULL,*preName=NULL;
  FILE *fileIn=NULL,*fileOut=NULL;
  bool bvtk=true;

/* ... macro camandos de leitura*/
  bool macroFlag; 
  char word[WORD_SIZE];
  char macro[][WORD_SIZE] = {"mesh","stop","config"
                            ,"pgeo"
                            }; 
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
           ,mesh->numel      ,mesh->ndm        ,mesh->maxno  
           ,mesh->ndfT       ,nameOut          ,bvtk             
           ,fileOut);
    printf("%s\n\n",DIF);
  }   
/*===================================================================*/
   
  }while(macroFlag && (!feof(fileIn)));

 
  return EXIT_SUCCESS;
}    
/*********************************************************************/
 
 
