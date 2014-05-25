/*==================== includes padroes =============================*/
#include<stdio.h>
#include<stdlib.h>
/*===================================================================*/

/*===================================================================*/
#include<File.h>
#include<Memoria.h>
#include<Mesh.h>
/*===================================================================*/


int main(int argc,char**argv){

/*...Memoria principal*/  
  Memoria m;
/*...*/
  Mesh* mesh=NULL;
/*... arquivo*/
  char *nameIn=NULL,*nameOut=NULL,*preName=NULL;
  FILE *fileIn=NULL,*fileOut=NULL;

/*=== sistema de memoria*/
  init_mem(&m,true);
/*===================================================================*/

/*=== abrindo arquivo de entrada*/ 
    
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

/*  fileIn = open_file(nameIn,"r");*/
 
/*===================================================================*/

/*=== arquivos de saida*/
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

/*===================================================================*/

   finalize_mem(&m,true);
   return EXIT_SUCCESS;
}    
/*********************************************************************/
 
 
