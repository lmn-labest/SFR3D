/********************* includes padroes ******************************/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
/*********************************************************************/

/*********************************************************************/
#include<CellLoop.h>
#include<Coo.h>
#include<File.h>
#include<Memoria.h>
#include<Mesh.h>
#include<WriteVtk.h>
#include<Sisteq.h>
#include<Solv.h>
#include<ReadFile.h>
#include<Reord.h>
/*********************************************************************/

/*********************************************************************/
#ifdef _DEBUG_ 
  #include<Debug.h>
#endif
/*********************************************************************/

int main(int argc,char**argv){

/*...Memoria principal*/  
  Memoria m;
/*... estrutura da dados para a malha*/
  Mesh *mesh=NULL;
/*... Sistema de equacao*/
  SistEq *sistEqT1=NULL;
  SistEq *sistEqD1=NULL;
/*... solver*/
  Solv *solvD1=NULL;
/*... reordenacao da malha*/
  Reord  *reordMesh=NULL;

/*Estrutura de dados*/
  char strIa[MNOMEPONTEIRO],strJa[MNOMEPONTEIRO];
  char strA[MNOMEPONTEIRO],strAd[MNOMEPONTEIRO];

/*... arquivo*/
  char *nameIn=NULL,*nameOut=NULL,*preName=NULL;
  FILE *fileIn=NULL,*fileOut=NULL;
  bool bvtk=true;

/*... loop nas celulas*/
/*Lib lib;*/


/* ... macro camandos de leitura*/
  bool macroFlag; 
  char word[WORD_SIZE],str[WORD_SIZE];
  char macro[][WORD_SIZE] = {"mesh" ,"stop","config"
                            ,"pgeo" ,"pcoob","presolvd"
                            ,"solvd","pcoo" ,"ptemp"}; 
/* ..................................................................*/

/*... Memoria principal(valor padrao - bytes)*/
  nmax = 10000;
/* ..................................................................*/


/*... estrutura de dados para malha*/
  mesh = (Mesh*) malloc(sizeof(Mesh));
  if(mesh == NULL){
    printf("Erro ponteiro mesh\n");
    exit(EXIT_FAILURE);
  }
/* ..................................................................*/

/*...*/  
  reordMesh = (Reord*) malloc(sizeof(Reord));
  if(reordMesh == NULL){
    printf("Erro ponteiro reord\n");
    exit(EXIT_FAILURE);
  }
  reordMesh->flag = false; 
/* ..................................................................*/
    
/*... abrindo ar quivo de entrada*/ 
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
 * macro: mesh - leitura da malha e inicializa das estruturas        *
 * de resolucao do problema                                          * 
 *===================================================================*/
    if((!strcmp(word,macro[0]))){
      printf("%s\n",DIF);
      printf("%s\n",word); 
      printf("%s\n",DIF);
/*... sistema de memoria*/
      initMem(&m,nmax,true);
/*... leitura da malha*/
      readFileFvMesh(&m,mesh,fileIn);
/*...................................................................*/

/*... calcula a vizinhaca do elementos*/
      viz(&m                 ,mesh->elm.node ,mesh->elm.adj.nelcon
         ,mesh->elm.nen      ,mesh->nnode
         ,mesh->numel        ,mesh->maxNo    ,mesh->maxViz);
/*...................................................................*/


/*... reodenando as celulas para dimuincao da banda*/
/* ..................................................................*/

/*...*/
      HccaAlloc(INT,&m,reordMesh->num,mesh->numel,"rNum" ,_AD_);
      printf("%s\n",DIF);
      printf("Reordenando a malha.\n");
      reord(&m                ,reordMesh->num,mesh->elm.adj.nelcon
           ,mesh->elm.adj.nViz,mesh->maxViz  ,mesh->numel   
           ,reordMesh->flag);
      printf("Malha reordenada.\n");
      printf("%s\n",DIF);
/*...................................................................*/


/*... calculo de propriedades geometricas recorrentes*/
      pGeomForm(mesh->node.x         ,mesh->elm.node
               ,mesh->elm.adj.nelcon ,mesh->elm.nen 
               ,mesh->elm.adj.nViz   ,mesh->elm.geomType
               ,mesh->elm.geom.cc    ,mesh->elm.geom.ksi 
               ,mesh->elm.geom.mksi  ,mesh->elm.geom.eta   
               ,mesh->elm.geom.meta  ,mesh->elm.geom.normal
               ,mesh->elm.geom.volume,mesh->elm.geom.xm   
               ,mesh->elm.geom.xmcc  ,mesh->elm.geom.mkm 
               ,mesh->elm.geom.dcca
               ,mesh->maxNo          ,mesh->maxViz
               ,mesh->ndm            ,mesh->numel);
#ifdef _DEBUG_GEOM
      testeGeom(mesh->elm.geom.cc  
               ,mesh->elm.geom.ksi   ,mesh->elm.geom.mksi
               ,mesh->elm.geom.eta   ,mesh->elm.geom.meta
               ,mesh->elm.geom.normal,mesh->elm.geom.volume
               ,mesh->elm.geom.xm    ,mesh->elm.geom.xmcc  
               ,mesh->elm.geom.mkm   ,mesh->elm.geom.dcca
               ,mesh->numel          ,mesh->ndm
               ,mesh->maxViz);
#endif
/*...................................................................*/
      strcpy(str,"KB");
      memoriaTotal(str);
      usoMemoria(&m,str);
//    mapVector(&m);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: stop - finalizacao do programa
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
 * macro: config - configuracao basica de excucao
 *===================================================================*/
    else if((!strcmp(word,macro[2]))){
      printf("%s\n",DIF);
      printf("%s\n",word);
      config(&bvtk,reordMesh,fileIn);
      printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: pgeo - escreve a geometria com os carregamentos
 *===================================================================*/
    else if((!strcmp(word,macro[3]))){
      printf("%s\n",DIF);
      printf("%s\n",word);
      fName(preName,0,6,&nameOut);
      wGeoVtk(&m                      ,mesh->node.x   
             ,mesh->elm.node          ,mesh->elm.mat    
             ,mesh->elm.nen           ,mesh->elm.geomType
             ,mesh->elm.material.prop ,mesh->elm.material.type 
             ,mesh->elm.faceRd1       ,mesh->elm.faceSd1
             ,mesh->nnode             ,mesh->numel    
             ,mesh->ndm               ,mesh->maxNo
             ,mesh->numat             ,mesh->ndfD    
             ,nameOut                 ,bvtk             
             ,fileOut);
      printf("%s\n\n",DIF);
    }   
/*===================================================================*/
   
/*===================================================================*
 * macro: pcoob- escreve a matriz de coeficientes no formato COO
 *===================================================================*/
    else if((!strcmp(word,macro[4]))){
      printf("%s\n",DIF);
      printf("%s\n",word);
      fName(preName,0,13,&nameOut);
/*...*/
      writeCoo(&m,sistEqD1->ia,sistEqD1->ja,sistEqD1->neq
              ,sistEqD1->au   ,sistEqD1->ad,sistEqD1->al        
              ,sistEqD1->nad  ,sistEqD1->storage
              ,sistEqD1->unsym,true
              ,nameOut);
/*...................................................................*/
      printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: presolvd - problema de difusao pura
 *===================================================================*/
    else if((!strcmp(word,macro[5]))){
      printf("%s\n",DIF);
      printf("%s\n",word);
/*... inicializando a estrutura de equacoes do problema*/
      solvD1 = (Solv*) malloc(sizeof(Solv));
      if(solvD1 == NULL){
        printf("Erro ponteiro solvD1\n");
        exit(EXIT_FAILURE);
      }
      solvD1->solver   = PCG;
      solvD1->tol      = 1.0e-14;
      solvD1->maxIt    = 5000;    
      solvD1->fileSolv = NULL;
      solvD1->log      = true;
/*..*/
      if(solvD1->log){  
        strcpy(nameOut,preName);
        strcat(nameOut,"_pcg");
        fName(preName,0,11,&nameOut);
        solvD1->fileSolv = openFile(nameOut,"w");
      }
/*...................................................................*/

/*...................................................................*/

/*... inicializa a estrutura do solver*/
      sistEqD1 = (SistEq*) malloc(sizeof(SistEq));
      if(sistEqD1 == NULL){
        printf("Erro ponteiro sistEqD1\n");
        exit(EXIT_FAILURE);
      }
      sistEqD1->storage = CSRD;
      sistEqD1->unsym   = false;    
/*...................................................................*/

/*... numeracao das equacoes*/
      HccaAlloc(INT,&m,sistEqD1->id
               ,mesh->numel*mesh->ndfD[0]
               ,"sistD1id",_AD_);
      printf("%s\n",DIF);
      printf("Numerando as equacoes.\n");
      sistEqD1->neq = numeq(&m,sistEqD1->id  ,reordMesh->num
                           ,mesh->elm.faceRd1  ,mesh->elm.nen
                           ,mesh->numel        ,mesh->maxViz
                           ,mesh->ndfD[0]);
      printf("Equacoes numeradas.\n");
      printf("%s\n",DIF);
/*...................................................................*/

/*...*/
      HccaAlloc(double                   ,&m        ,sistEqD1->b0
               ,sistEqD1->neq            ,"sistD1b0",_AD_);
      HccaAlloc(double                   ,&m        ,sistEqD1->b 
               ,sistEqD1->neq            ,"sistD1b ",_AD_);
      HccaAlloc(double                   ,&m        ,sistEqD1->x 
               ,sistEqD1->neq            ,"sistD1x ",_AD_);
      zero(sistEqD1->b0,sistEqD1->neq,"double");
      zero(sistEqD1->b ,sistEqD1->neq,"double");
      zero(sistEqD1->x ,sistEqD1->neq,"double");
/*...................................................................*/

/*... Estrutura de Dados*/
      strcpy(strIa,"iaD1");
      strcpy(strJa,"JaD1");
      strcpy(strAd,"aDD1");
      strcpy(strA ,"aD1");
      printf("Montagem da estrura de dados esparsa.\n");
      dataStruct(&m,sistEqD1->id   ,reordMesh->num,mesh->elm.adj.nelcon
                ,mesh->elm.adj.nViz,mesh->numel   ,mesh->maxViz
                ,mesh->ndfD[0]     ,strIa         ,strJa
                ,strAd             ,strA          ,sistEqD1);
      printf("Estrutuda montada.\n");
/*...................................................................*/
      printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: solvd - problema de difusao pura
 *===================================================================*/
    else if((!strcmp(word,macro[6]))){
      printf("%s\n",DIF);
      printf("%s\n",word);
      printf("%s\n",DIF);
/*... restricoes por centro de celula u0 e cargas por volume b0*/
      cellPload(mesh->elm.faceRd1    ,mesh->elm.faceSd1
               ,mesh->elm.geom.volume
               ,mesh->elm.temp       ,sistEqD1->b0
               ,mesh->numel          ,mesh->ndfD[0]
               ,mesh->maxViz);
/*...................................................................*/

/*... calculo de: A(i),b(i)
                  R(i) = b(i) - A(i)x(i)*/
      printf("%s\n",DIF);
      printf("Montagem do sistema de equacoes.\n");
      systForm(mesh->node.x           ,mesh->elm.node       
             ,mesh->elm.adj.nelcon    ,mesh->elm.nen           
             ,mesh->elm.adj.nViz      ,mesh->elm.geomType          
             ,mesh->elm.material.prop ,mesh->elm.material.type 
             ,mesh->elm.mat           ,mesh->elm.geom.cc 
             ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
             ,mesh->elm.geom.eta      ,mesh->elm.geom.meta    
             ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
             ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc    
             ,mesh->elm.geom.mkm      ,mesh->elm.geom.dcca
             ,sistEqD1->ia            ,sistEqD1->ja      
             ,sistEqD1->ad            ,sistEqD1->al       
             ,sistEqD1->b             ,sistEqD1->id       
             ,mesh->elm.faceRd1       ,mesh->elm.faceSd1       
             ,mesh->elm.temp                                               
             ,mesh->maxNo             ,mesh->maxViz
             ,mesh->ndm               ,mesh->numel
             ,mesh->ndfD[0]           ,sistEqD1->storage
             ,true                    ,true   
             ,sistEqD1->unsym);   
      printf("Sistema montado.\n");
      printf("%s\n",DIF);
/*...................................................................*/

/*... soma o vetor b = b + b0*/
      addVector(1.0e0        ,sistEqD1->b
               ,1.0e0        ,sistEqD1->b0
               ,sistEqD1->neq,sistEqD1->b);
/*...................................................................*/

/*...*/
      printf("%s\n",DIF);
      printf("Resolucao do sistema de equacoes.\n");
//    solverC(&m               ,sistEqD1->neq ,sistEqD1->nad
//           ,sistEqD1->ia     ,sistEqD1->ja  
//           ,sistEqD1->al     ,sistEqD1->ad,sistEqD1->au
//           ,sistEqD1->b      ,sistEqD1->x
//           ,solvD1->tol      ,solvD1->maxIt     
//           ,sistEqD1->storage,solvD1->solver
//           ,solvD1->fileSolv ,solvD1->log  
//           ,false            ,false
//           ,sistEqD1->unsym  ,false);
      printf("Sistema resolvido.\n");
      printf("%s\n",DIF);
/*...................................................................*/

/*...*/
#ifdef _DEBUG_
      testeSist(sistEqD1->ia       ,sistEqD1->ja
               ,sistEqD1->au       ,sistEqD1->ad
               ,sistEqD1->al       ,sistEqD1->b
               ,sistEqD1->neq      ,sistEqD1->unsym);
#endif
/*...................................................................*/

/*... x -> temp*/
      updateCellU(mesh->elm.temp,sistEqD1->x
            ,sistEqD1->id  
            ,mesh->numel   , mesh->ndfD[0]);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos*/
     interCellNode(&m
                  ,mesh->node.temp,mesh->elm.temp
                  ,mesh->elm.node                 
                  ,mesh->elm.nen
                  ,mesh->numel    ,mesh->nnode    
                  ,mesh->maxNo    ,mesh->ndfD[0]   
                  ,1);
/*...................................................................*/
      printf("%s\n\n",DIF);
    }
/*===================================================================*/

/*===================================================================*
 * macro: pcoo - escreve a matriz de coeficientes no formato COO
 *===================================================================*/
    else if((!strcmp(word,macro[7]))){
      printf("%s\n",DIF);
      printf("%s\n",word);
      fName(preName,0,12,&nameOut);
/*... matriz A*/
      writeCoo(&m,sistEqD1->ia,sistEqD1->ja,sistEqD1->neq
              ,sistEqD1->au   ,sistEqD1->ad,sistEqD1->al        
              ,sistEqD1->nad  ,sistEqD1->storage
              ,sistEqD1->unsym,false 
              ,nameOut);
/*...................................................................*/

/*... vetor de forcas b*/      
      fName(preName,0,14,&nameOut);
      for(int i=0;i<sistEqD1->neq;i++)
        sistEqD1->b[i] /= sistEqD1->ad[i];

      writeCooB(sistEqD1->b,sistEqD1->neq,nameOut);
/*...................................................................*/
      printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: ptemp - escreve os arquivos dos resultados da temperatura
 *===================================================================*/
    else if((!strcmp(word,macro[8]))){
      printf("%s\n",DIF);
      printf("%s\n",word);
      fName(preName,0,8,&nameOut);
/*...*/
      wResVtk(&m             ,mesh->node.x      
             ,mesh->elm.node ,mesh->elm.mat    
             ,mesh->elm.nen  ,mesh->elm.geomType
             ,mesh->elm.temp ,mesh->node.temp
             ,mesh->nnode    ,mesh->numel  
             ,mesh->ndm      ,mesh->maxNo 
             ,mesh->numat    ,mesh->ndfD    
             ,nameOut        ,bvtk    
             ,fileOut);
/*...................................................................*/
      printf("%s\n\n",DIF);
    }   
/*===================================================================*/
  }while(macroFlag && (!feof(fileIn)));

 
  return EXIT_SUCCESS;
}    
/*********************************************************************/
 
 
