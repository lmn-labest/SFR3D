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
#include<WriteCsv.h>
#include<WriteLog.h>
#include<PartMesh.h>
#include<ParallelMpi.h>
#include<Sisteq.h>
#include<Solv.h>
#include<ReadFile.h>
#include<Reord.h>
#include<Transient.h>
/*********************************************************************/

/*********************************************************************/
#ifdef _DEBUG_ 
  #include<Debug.h>
#endif
/*********************************************************************/

/*********************************************************************/


int main(int argc,char**argv){

/*...Memoria principal*/  
  Memoria m;
/*... estrutura da dados para a malha*/
  Mesh *mesh=NULL,*mesh0=NULL;
/*... Sistema de equacao*/
//SistEq *sistEqT1=NULL;
  SistEq *sistEqD1=NULL;
/*... solver*/
  Solv *solvD1=NULL;
  bool fSolvD1 = false;
/*... reordenacao da malha*/
  Reord  *reordMesh=NULL;

/*... particionamento*/
  PartMesh *pMesh = NULL;

/*...*/
  char loopWord[100][MAX_LINE];
  unsigned short kLoop,jLoop;
  bool flWord=false;

/*... Estrutura de dados*/
  char strIa[MNOMEPONTEIRO],strJa[MNOMEPONTEIRO];
  char strA[MNOMEPONTEIRO],strAd[MNOMEPONTEIRO];

/*... arquivo*/
  char *nameIn=NULL,*nameOut=NULL,*preName=NULL,*auxName=NULL;
  FILE *fileIn=NULL,*fileOut=NULL,*fileLog=NULL;
  char str1[100],str2[100],str3[100],str4[100];
  FileOpt opt;

/*...*/
  Scheme sc;
/*...................................................................*/
 
/*... loop nas celulas*/
/*Lib lib;*/
  
/*...*/
  DOUBLE rCell,rCell0,conv;
/*...*/
  int i; 

/* ... macro camandos de leitura*/
  bool macroFlag; 
  char word[WORD_SIZE],str[WORD_SIZE];
  char macro[][WORD_SIZE] = {"mesh"      ,"stop"      ,"config"
                            ,"pgeo"      ,"pcoob"     ,"presolvD1"
                            ,"solvD1"    ,"pcoo"      ,"pD1"
                            ,"nlItD1"    ,"pD1CsvCell","pD1CsvNode"
                            ,"transient" ,"timeUpdate","partd"}; 
/* ..................................................................*/

/*... Memoria principal(valor padrao - bytes)*/
  nmax = 200000;
/* ..................................................................*/

/*... definicao de variaveis globais*/
  oneDivTree = 1.e0/3.e0;
/* ..................................................................*/

/* ... opcoes de arquivos */                                           
  opt.bVtk       = false;
  opt.fItPlotRes = false;
  opt.fItPlot    = false;
/* ..................................................................*/

/*... Mpi*/
  mpiStart(&argc,argv);
/*...................................................................*/

/*... Time*/
  tm.adjcency          = 0.e0;
  tm.geom              = 0.e0;
  tm.leastSquareMatrix = 0.e0;
  tm.reord             = 0.e0;
/*...*/
  tm.solvD1            = 0.e0;
  tm.numeqD1           = 0.e0;
  tm.dataStructD1      = 0.e0;
  tm.CellPloadD1       = 0.e0;
  tm.CellTransientD1   = 0.e0;
  tm.systFormD1        = 0.e0;
  tm.rcGradD1          = 0.e0;
  tm.solvEdoD1         = 0.e0;
/*... Blas*/
  tm.matVecOverHeadMpi = 0.e0;
  tm.matVecSparse      = 0.e0;
  tm.dot               = 0.e0;
  tm.dotOverHeadMpi    = 0.e0;
/*... Iterativos  */
  tm.pcg               = 0.e0;
/*... particionamento*/
  tm.partdMesh         = 0.e0;
  tm.partdMeshCom      = 0.e0;
/*...*/
  tm.overHeadCelMpi    = 0.e0;
  tm.overHeadNodMpi    = 0.e0;
  tm.overHeadNeqMpi    = 0.e0;
  tm.overHeadGCelMpi   = 0.e0;
  tm.overHeadGNodMpi   = 0.e0;
  tm.overHeadTotalMpi  = 0.e0;

/*... precondicionador*/
  tm.precondDiag       = 0.e0;
  tm.total             = getTimeC();
/*...................................................................*/

/*... particionamento da malha*/
  pMesh = (PartMesh*) malloc(sizeof(PartMesh));
  ERRO_MALLOC(pMesh,"pMesh",__LINE__,__FILE__,__func__);
/*...*/
  pMesh->fPartMesh      = false;
  pMesh->fPrintMesh     = false;
  pMesh->fPrintMeshPart = false;
/*...................................................................*/

/*... estrutura de dados para malha*/
  mesh0 = (Mesh*) malloc(sizeof(Mesh));
  ERRO_MALLOC(mesh0,"mesh0",__LINE__,__FILE__,__func__);

/*... tecnica padrao de resconstrucao de gradiente*/
  sc.rcGrad = RCGRADGAUSSN; 
  sc.nlD1.maxIt   = 100; 
  sc.nlD1.tol     = 1.e-6; 
/* ..................................................................*/

/*...*/  
  sc.ddt.flag     = false;
  sc.ddt.dt       = 1.e0;
  sc.ddt.total    = 1.e0;
  sc.ddt.timeStep = 1;
/* ..................................................................*/

/*...*/  
  reordMesh = (Reord*) malloc(sizeof(Reord));
  ERRO_MALLOC(reordMesh,"reordMesh",__LINE__,__FILE__,__func__);
  reordMesh->flag = false; 
/* ..................................................................*/
    
/*... abrindo ar quivo de entrada*/ 
  nameIn = (char *) malloc(sizeof(char)*MAX_STR_LEN_IN);
    
  if( argc > 1)
    strcpy(nameIn,argv[1]);
  else{
    if(!mpiVar.myId ) printf("Arquivo de dados:\n");
    if(!mpiVar.myId ) scanf("%s",nameIn);
  }

  fileIn = openFile(nameIn,"r");
/*...................................................................*/

/*... arquivos de saida*/
  preName = (char *) malloc(sizeof(char)*MAX_STR_LEN_SUFIXO);
  if(preName == NULL){
    printf("Erro ponteiro prename\n");
    exit(EXIT_FAILURE);
  }
  
  auxName = (char *) malloc(sizeof(char)*MAX_STR_LEN_SUFIXO);
  if(preName == NULL){
    printf("Erro ponteiro auxName\n");
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
    if(!mpiVar.myId ) printf("Prefixo do arquivo de saida:\n");
    if(!mpiVar.myId ) scanf("%s",preName);
  }
/*...................................................................*/
  
/*loop de execucao*/
  macroFlag = true;
  do{
/*... macros na marcro trasient*/
    if(flWord){
      if( jLoop > kLoop) { 
        ERRO_GERAL(__FILE__,__func__,__LINE__,
                   "Numero de comandos na string trasient execedido"); 
      }
      strcpy(word,loopWord[jLoop]);
      jLoop++;
    }
/*... leitura da macro*/
    else 
      readMacro(fileIn,word,false);
/*...................................................................*/

/*===================================================================*
 * macro: mesh - leitura da malha e inicializa das estruturas        *
 * de resolucao do problema                                          * 
 *===================================================================*/
    if((!strcmp(word,macro[0]))){
      if(!mpiVar.myId) printf("%s\n",DIF);
      if(!mpiVar.myId) printf("%s\n",word); 
      if(!mpiVar.myId) printf("%s\n",DIF);
/*... sistema de memoria*/
      initMem(&m,nmax,false);
/*... leitura da malha*/
      if(!mpiVar.myId)
        readFileFvMesh(&m,mesh0,fileIn);
      mpiWait();
/*...................................................................*/

/*... calcula a vizinhaca do elementos*/
      if(!mpiVar.myId){
        tm.adjcency = getTimeC();
        viz(&m                 ,mesh0->elm.node ,mesh0->elm.adj.nelcon
           ,mesh0->elm.nen     ,mesh0->nnode    ,mesh0->numel     
           ,mesh0->maxNo       ,mesh0->maxViz   ,mesh0->ndm);
        tm.adjcency = getTimeC() - tm.adjcency;
      }
/*...................................................................*/
/*... particionamento da malha*/
      if(pMesh->fPartMesh && mpiVar.nPrcs > 1){
        tm.partdMesh = getTimeC() - tm.partdMesh;
        if(!mpiVar.myId){
          partMesh(&m         
                  ,mesh0->node.x  ,mesh0->elm.node
                  ,mesh0->elm.nen
                  ,mesh0->nnode   ,mesh0->numel
                  ,pMesh   
                  ,mesh0->ndm     ,mesh0->maxNo 
                  ,mpiVar.nPrcs);
/*... */
          if(pMesh->fPrintMesh){
            fName(preName,mpiVar.nPrcs,0,1,&nameOut);
            wPartVtk(&m            
                    ,mesh0->node.x      ,mesh0->elm.node              
                    ,mesh0->elm.nen     ,mesh0->elm.geomType
                    ,pMesh->np          ,pMesh->ep          
                    ,mesh0->nnode       ,mesh0->numel    
                    ,mesh0->ndm               
                    ,mesh0->maxNo       ,mesh0->maxViz
                    ,nameOut            ,opt.bVtk             
                    ,fileOut);
          }
/*...................................................................*/
        }
/*...*/
        mesh = (Mesh*) malloc(sizeof(Mesh));
        ERRO_MALLOC(mesh,"mesh",__LINE__,__FILE__,__func__);
/*...*/
        tm.partdMeshCom = getTimeC() - tm.partdMeshCom;
        comunicateMesh(&m
                      ,mesh0     ,mesh
                      ,pMesh
                      ,loadsD1);
        tm.partdMeshCom = getTimeC() - tm.partdMeshCom;
/*...................................................................*/

        if(pMesh->fPrintMeshPart){
          fName(preName,mpiVar.nPrcs,mpiVar.myId,2,&nameOut);
          wMeshPartVtk(&m            
                      ,mesh->node.x      ,mesh->elm.node              
                      ,mesh->elm.nen     ,mesh->elm.geomType
                      ,mesh->nnode       ,mesh->numel    
                      ,mesh->ndm               
                      ,mesh->maxNo       ,mesh->maxViz
                      ,nameOut           ,opt.bVtk             
                      ,fileOut);
          fName(preName,mpiVar.nPrcs,mpiVar.myId,18,&nameOut);
          printMap(*pMesh
                  ,mesh->nnode ,mesh->numel
                  ,mpiVar.myId ,nameOut
                  ,fileOut);
        }
/*...................................................................*/
        tm.partdMesh = getTimeC() - tm.partdMesh;
      }
/*...................................................................*/

/*... sequencial*/
      else{
        mesh = mesh0;
      }
/*...................................................................*/

/*... calculo de propriedades geometricas recorrentes*/
      tm.geom = getTimeC() - tm.geom;
      pGeomForm(mesh->node.x         ,mesh->elm.node
               ,mesh->elm.adj.nelcon ,mesh->elm.nen 
               ,mesh->elm.adj.nViz   ,mesh->elm.geomType
               ,mesh->elm.geom.cc    ,mesh->elm.geom.ksi 
               ,mesh->elm.geom.mksi  ,mesh->elm.geom.eta   
               ,mesh->elm.geom.fArea ,mesh->elm.geom.normal
               ,mesh->elm.geom.volume,mesh->elm.geom.xm   
               ,mesh->elm.geom.xmcc                 
               ,mesh->elm.geom.vSkew ,mesh->elm.geom.mvSkew  
               ,mesh->elm.geom.dcca
               ,mesh->maxNo          ,mesh->maxViz
               ,mesh->ndm            ,mesh->numelNov);
      tm.geom = getTimeC() - tm.geom;
/*...................................................................*/

/*... geom(Cel)*/
      dGlobalCel(&m                  ,pMesh
                ,mesh0->elm.geom.cc  ,mesh->elm.geom.cc
                ,mesh->numelNov 
                ,mesh->ndm           ,1);
/*...................................................................*/
      
/*... reconstrucao de gradiente least square*/
      if(sc.rcGrad ==  RCLSQUARE || sc.rcGrad ==  RCLSQUAREQR){
        if(!mpiVar.myId ){
          printf("%s\n",DIF);
          printf("Least Square ...\n");
        }
/*... wleastSqaure*/
        HccaAlloc(DOUBLE,&m,mesh->elm.leastSquare
                 ,mesh->numelNov*mesh->maxViz*mesh->ndm
                 ,"leastSquare",_AD_);
/*... leastSqaure QR*/
        if(sc.rcGrad ==  RCLSQUAREQR){
          if(mesh->ndm == 2){ 
            HccaAlloc(DOUBLE,&m,mesh->elm.leastSquareR            
                     ,mesh->numelNov*3
                     ,"lSquareR   ",_AD_);
          }
          else if(mesh->ndm == 3){ 
            HccaAlloc(DOUBLE,&m,mesh->elm.leastSquareR            
                     ,mesh->numelNov*6
                     ,"lSquareR   ",_AD_);
          }
        }
/*...*/
        tm.leastSquareMatrix = getTimeC() - tm.leastSquareMatrix;
        rcLeastSquare(mesh->elm.geom.ksi   ,mesh->elm.geom.mksi
                     ,mesh->elm.leastSquare,mesh->elm.leastSquareR
                     ,mesh->elm.adj.nViz       
                     ,mesh->numelNov       ,mesh->maxViz
                     ,sc.rcGrad            ,mesh->ndm);
        tm.leastSquareMatrix = getTimeC() - tm.leastSquareMatrix;
/*...................................................................*/
        if(!mpiVar.myId ){
          printf("Least Square.\n");
          printf("%s\n",DIF);
        }
      }
/*...................................................................*/

/*... qualidade da malha*/
      meshQuality(&mesh->mQuality
                 ,mesh->elm.adj.nViz     ,mesh->elm.geom.volume
                 ,mesh->elm.geom.ksi     ,mesh->elm.geom.normal
                 ,mesh->elm.geom.mvSkew 
                 ,mesh->maxViz           ,mesh->ndm
                 ,mesh->numelNov);
/*... qualidade da malha global*/
      if(mpiVar.nPrcs > 1){
        globalMeshQuality(&mesh->mQuality);
      }  
/*...................................................................*/
         
/*... reodenando as celulas para dimuincao da banda*/
      HccaAlloc(INT,&m,reordMesh->num,mesh->numel,"rNum" ,_AD_);
      printf("%s\n",DIF);
      if(!mpiVar.myId ) printf("Reordenando a malha ...\n");
        tm.reord = getTimeC() - tm.reord;
        reord(&m                ,reordMesh->num,mesh->elm.adj.nelcon
             ,mesh->elm.adj.nViz,mesh->maxViz  
             ,mesh->numel       ,mesh->numelNov
             ,reordMesh->flag   ,mpiVar.nPrcs);
        tm.reord = getTimeC() - tm.reord;
      if(!mpiVar.myId ){
        printf("Malha reordenada.\n");
        printf("%s\n",DIF);
      }
/*...................................................................*/
      
/*...*/
      if(!mpiVar.myId ){
        strcpy(str,"MB");
        memoriaTotal(str);
        usoMemoria(&m,str);
      }
/*...................................................................*/

    }   
/*===================================================================*/

/*===================================================================*
 * macro: stop : finalizacao do programa
 *===================================================================*/
    else if((!strcmp(word,macro[1]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word); 
        printf("%s\n\n",DIF);
      }
      tm.total = getTimeC() - tm.total;
/*... */
      fName(preName,mpiVar.nPrcs,mpiVar.myId,7,&nameOut);
      fileLog = openFile(nameOut,"w");
      writeLog(*mesh     ,sc
              ,solvD1    ,sistEqD1
              ,tm       
              ,solvD1         
              ,nameIn   ,fileLog);
      fclose(fileLog);
/*...................................................................*/

/*... fechando o arquivo log pcg D1*/
      if(fSolvD1 && solvD1->log && !mpiVar.myId)  
        fclose(solvD1->fileSolv);
/*... fechando o arquivo log nao linear D1*/      
      if(fSolvD1 && opt.fItPlot && !mpiVar.myId)  
        fclose(opt.fileItPlot[FITPLOTD1]);
      finalizeMem(&m,false);
      macroFlag = false;
    }    
/*===================================================================*/

/*===================================================================*
 * macro: config : configuracao basica de excucao
 *===================================================================*/
    else if((!strcmp(word,macro[2]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
      
      config(&opt          ,reordMesh
            ,&sc.rcGrad ,fileIn);
      
      if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: pgeo : escreve a geometria com os carregamentos
 *===================================================================*/
    else if((!strcmp(word,macro[3]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
/*... geometrica completa*/
        fName(preName,0,0,6,&nameOut);
        wGeoVtk(&m                       ,mesh0->node.x   
               ,mesh0->elm.node          ,mesh0->elm.mat    
               ,mesh0->elm.nen           ,mesh0->elm.geomType
               ,mesh0->elm.material.prop ,mesh0->elm.material.type 
               ,mesh0->elm.faceRd1       ,mesh0->elm.faceLoadD1
               ,mesh0->nnode             ,mesh0->numel    
               ,mesh0->ndm               
               ,mesh0->maxNo             ,mesh0->maxViz
               ,mesh0->numat             ,mesh0->ndfD    
               ,nameOut                  ,opt.bVtk             
               ,fileOut);  
/*... face com cargas*/
       fName(preName,0,0,17,&nameOut);
       wGeoFaceVtk(&m                    ,mesh0->node.x        
             ,mesh0->elm.node          ,mesh0->elm.nen      
             ,mesh0->elm.geomType
             ,mesh0->elm.faceRd1       ,mesh0->elm.faceLoadD1
             ,mesh0->nnode             ,mesh0->numel    
             ,mesh0->ndm               ,mesh0->ndfD[0]
             ,mesh0->maxViz            ,mesh0->maxNo
             ,nameOut                  ,opt.bVtk             
             ,fileOut);  
      }
    }   
/*===================================================================*/
   
/*===================================================================*
 * macro: pcoob : escreve a matriz de coeficientes no formato COO
 *===================================================================*/
    else if((!strcmp(word,macro[4]))){
      if(mpiVar.nPrcs == 1){
        printf("%s\n",DIF);
        printf("%s\n",word);
       
        fName(preName,0,0,13,&nameOut);
/*...*/
        writeCoo(&m,sistEqD1->ia,sistEqD1->ja,sistEqD1->neq
                ,sistEqD1->au   ,sistEqD1->ad,sistEqD1->al        
                ,sistEqD1->nad  ,sistEqD1->storage
                ,sistEqD1->unsym,true
                ,nameOut);
/*...................................................................*/
        printf("%s\n\n",DIF);
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: presolvd : problema de difusao pura
 *===================================================================*/
    else if((!strcmp(word,macro[5]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
      

/*... inicializando a estrutura de equacoes do problema*/
      solvD1 = (Solv*) malloc(sizeof(Solv));
      if(solvD1 == NULL){
        printf("Erro ponteiro solvD1\n");
        exit(EXIT_FAILURE);
      }
      fSolvD1          = true;
      solvD1->solver   = PCG;
      solvD1->tol      = smachn();
      solvD1->maxIt    = 50000;    
      solvD1->fileSolv = NULL;
      solvD1->log      = true;
      solvD1->flag     = true;
/*...................................................................*/

/*...*/
      if(solvD1->log && !mpiVar.myId){  
        strcpy(auxName,preName);
        strcat(auxName,"_pcg_D1");
        fName(auxName,mpiVar.nPrcs,0,11,&nameOut);
        solvD1->fileSolv = openFile(nameOut,"w");
      }
/*...................................................................*/

/*...*/
      if(opt.fItPlot && !mpiVar.myId){  
        strcpy(auxName,preName);
        strcat(auxName,"_D1");
        fName(auxName,mpiVar.nPrcs,0,10,&nameOut);
        opt.fileItPlot[FITPLOTD1] = openFile(nameOut,"w");
        fprintf(opt.fileItPlot[FITPLOTD1]
               ,"#D1\n#it ||b||/||b0|| ||b||\n");
      }
/*...................................................................*/

/*... inicializa a estrutura do solver*/
      sistEqD1 = (SistEq*) malloc(sizeof(SistEq));
      if(sistEqD1 == NULL){
        printf("Erro ponteiro sistEqD1\n");
        exit(EXIT_FAILURE);
      }
      sistEqD1->storage = CSRD;
//    sistEqD1->storage = ELLPACK;
      sistEqD1->unsym   = false; 
/*...................................................................*/

/*... config*/
      readMacro(fileIn,word,false);
      if(!strcmp(word,"config:")){
/*... solver*/        
        readMacro(fileIn,word,false);
        setSolver(word,&solvD1->solver); 
        
/*... DataStruct*/    
        readMacro(fileIn,word,false);
        setDataStruct(word,&sistEqD1->storage); 

/*... simetria*/        
        readMacro(fileIn,word,false);
        if(!strcmp(word,"sym"))
          sistEqD1->unsym   = false;    
        else if(!strcmp("unSym",word))
          sistEqD1->unsym   = true;    
        
/*... */        
        fscanf(fileIn,"%u" ,&solvD1->maxIt);
        fscanf(fileIn,"%lf",&solvD1->tol);

        if( solvD1->tol == 0.e0) 
          solvD1->tol = smachn();

        if(!mpiVar.myId ) printf("MaxIt     : %d\n",solvD1->maxIt);
        if(!mpiVar.myId ) printf("Tol       : %e\n",solvD1->tol);
      
        if(solvD1->solver == PCG && !mpiVar.myId ){     
          printf("Solver    : PCG\n");
        }

        if(sistEqD1->storage == CSRD && !mpiVar.myId ){     
          printf("DataStruct: CSRD\n");
        }
        else if(sistEqD1->storage == ELLPACK && !mpiVar.myId ){     
          printf("DataStruct: ELLPACK\n");
        }
        
        else if(sistEqD1->storage == CSRDCOO && !mpiVar.myId ){     
          printf("DataStruct: CSRDCOO\n");
        }
        
        if(sistEqD1->unsym && !mpiVar.myId ){
          printf("Matrix    : unsymetric\n");
        } 
        else if (!mpiVar.myId){
          printf("Matrix    : symetric\n");
        } 
      }
/*...................................................................*/

/*... numeracao das equacoes*/
      HccaAlloc(INT,&m,sistEqD1->id
               ,mesh->numel*mesh->ndfD[0]
               ,"sistD1id",_AD_);
      if(!mpiVar.myId){
        printf("%s\n",DIF);
        printf("Numerando as equacoes.\n");
      }
      tm.numeqD1 = getTimeC() - tm.numeqD1;
      sistEqD1->neq = numeq(sistEqD1->id       ,reordMesh->num
                           ,mesh->elm.faceRd1  ,mesh->elm.adj.nViz
                           ,mesh->numel        ,mesh->maxViz
                           ,mesh->ndfD[0]);
      tm.numeqD1 = getTimeC() - tm.numeqD1;
      if(!mpiVar.myId){
        printf("Equacoes numeradas.\n");
        printf("%s\n",DIF);
      }
/*...................................................................*/

/*...*/
      if( mpiVar.nPrcs > 1) {      
        tm.numeqD1 = getTimeC() - tm.numeqD1;
        sistEqD1->neqNov = countEq(reordMesh->num
                            ,mesh->elm.faceRd1  ,mesh->elm.adj.nViz
                            ,mesh->numelNov     ,mesh->maxViz
                            ,mesh->ndfD[0]);
        tm.numeqD1 = getTimeC() - tm.numeqD1;
      }
      else{
        sistEqD1->neqNov = sistEqD1->neq;
      }
/*...................................................................*/

/*...*/
      HccaAlloc(DOUBLE                   ,&m        ,sistEqD1->b0
               ,sistEqD1->neq                  ,"sistD1b0",_AD_);
      HccaAlloc(DOUBLE                   ,&m        ,sistEqD1->b 
               ,sistEqD1->neq                 ,"sistD1b ",_AD_);
      HccaAlloc(DOUBLE                   ,&m        ,sistEqD1->x 
               ,sistEqD1->neq            ,"sistD1x ",_AD_);
      zero(sistEqD1->b0,sistEqD1->neq    ,DOUBLEC);
      zero(sistEqD1->b ,sistEqD1->neq    ,DOUBLEC);
      zero(sistEqD1->x ,sistEqD1->neq    ,DOUBLEC);
/*...................................................................*/

/*... Estrutura de Dados*/
      strcpy(strIa,"iaD1");
      strcpy(strJa,"JaD1");
      strcpy(strAd,"aDD1");
      strcpy(strA ,"aD1");
      if(!mpiVar.myId) printf("Montagem da estrura de dados esparsa.\n");
      tm.dataStructD1 = getTimeC() - tm.dataStructD1 ;
      dataStruct(&m,sistEqD1->id   ,reordMesh->num,mesh->elm.adj.nelcon
                ,mesh->elm.adj.nViz,mesh->numelNov,mesh->maxViz
                ,mesh->ndfD[0]     ,strIa         ,strJa
                ,strAd             ,strA          ,sistEqD1);
      tm.dataStructD1 = getTimeC() - tm.dataStructD1 ;
      if(!mpiVar.myId) printf("Estrutuda montada.\n");
/*...................................................................*/


/*... mapa de equacoes para comunicacao*/
      if( mpiVar.nPrcs > 1) {    
        front(&m,pMesh,sistEqD1,mesh->ndfD[0]);  
      } 
/*...................................................................*/

/*... informacao da memoria total usada*/
      if(!mpiVar.myId  ) {
        strcpy(str,"MB");
        memoriaTotal(str);
        usoMemoria(&m,str);
      }
/*...................................................................*/
      if(!mpiVar.myId  ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: solvd : problema de difusao pura
 *===================================================================*/
    else if((!strcmp(word,macro[6]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
        printf("%s\n",DIF);
      }
      mpiWait();
      tm.solvEdoD1    = getTimeC() - tm.solvEdoD1;
/*...*/
      if(solvD1 == NULL){
        printf("Estrutara de dados nao montada para o solvD1!!!\n");
        exit(EXIT_FAILURE);
      }
/*...................................................................*/
      
      zero(sistEqD1->b0,sistEqD1->neqNov,DOUBLEC);
/*... restricoes por centro de celula u0 e cargas por volume b0*/
      tm.CellPloadD1 = getTimeC() - tm.CellPloadD1;
      cellPload(loadsD1
               ,mesh->elm.faceRd1    ,mesh->elm.faceLoadD1
               ,mesh->elm.geom.volume,sistEqD1->id 
               ,mesh->elm.uD1        ,sistEqD1->b0
               ,mesh->numelNov       ,mesh->ndfD[0]
               ,mesh->maxViz);
      tm.CellPloadD1 = getTimeC() - tm.CellPloadD1;
/*...................................................................*/

/*... discretizacao temporal*/
      if(sc.ddt.flag){
        tm.CellTransientD1 = getTimeC() - tm.CellTransientD1;
        cellTransient(mesh->elm.geom.volume   ,sistEqD1->id     
                     ,mesh->elm.u0D1          ,mesh->elm.uD1
                     ,mesh->elm.densityUd1        
                     ,sistEqD1->b0
                     ,mesh->numelNov           ,mesh->ndfD[0]
                     ,sc.ddt.type              ,true);
/*... u(n-1) = u(n)*/
        alphaProdVector(1.e0,mesh->elm.uD1
                       ,mesh->numel       ,mesh->elm.u0D1); 
/*...................................................................*/
        tm.CellTransientD1 = getTimeC() - tm.CellTransientD1;
      }
/*...................................................................*/

/*... correcao nao ortoganal*/ 
      for(i=0;i<sc.nlD1.maxIt;i++){

/*... calculo de: A(i),b(i)*/
        tm.systFormD1 = getTimeC() - tm.systFormD1;
        systFormDif(mesh->elm.node          ,mesh->elm.adj.nelcon  
                   ,mesh->elm.nen           ,mesh->elm.adj.nViz   
                   ,mesh->elm.geomType      ,mesh->elm.material.prop 
                   ,mesh->elm.material.type ,mesh->elm.mat   
                   ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
                   ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea    
                   ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
                   ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc    
                   ,mesh->elm.geom.vSkew    ,mesh->elm.geom.mvSkew   
                   ,mesh->elm.geom.dcca     ,mesh->elm.densityUd1
                   ,sistEqD1->ia            ,sistEqD1->ja      
                   ,sistEqD1->ad            ,sistEqD1->al       
                   ,sistEqD1->b             ,sistEqD1->id       
                   ,mesh->elm.faceRd1       ,mesh->elm.faceLoadD1  
                   ,mesh->elm.uD1           ,mesh->elm.gradUd1           
                   ,mesh->elm.rCellUd1      ,sc.ddt
                   ,sistEqD1->neqNov        ,sistEqD1->nad
                   ,sistEqD1->nadr      
                   ,mesh->maxNo             ,mesh->maxViz
                   ,mesh->ndm               ,mesh->numelNov
                   ,mesh->ndfD[0]           ,sistEqD1->storage
                   ,true                    ,true   
                   ,true                    ,sistEqD1->unsym);   
        tm.systFormD1 = getTimeC() - tm.systFormD1;
/*...................................................................*/

/*... soma o vetor b(i) = b(i) + b0(i)*/
        addVector(1.0e0           ,sistEqD1->b
                 ,1.0e0           ,sistEqD1->b0
                 ,sistEqD1->neqNov,sistEqD1->b);
/*...................................................................*/

/*... soma o vetor R(i) = R(i) + b0(i)*/
        updateCellValue(mesh->elm.rCellUd1 ,sistEqD1->b0
                       ,sistEqD1->id       ,&sistEqD1->iNeq
                       ,mesh->numelNov     ,mesh->ndfD[0]
                       ,true               ,false);
/*...................................................................*/

/*...*/ 
        if( i == 0 ){
          rCell  = rCell0 = sqrt(dot(mesh->elm.rCellUd1
                                ,mesh->elm.rCellUd1 
                                ,mesh->numelNov));
          conv   = rCell0*sc.nlD1.tol;
        }
        else
          rCell  = sqrt(dot(mesh->elm.rCellUd1 
                       ,mesh->elm.rCellUd1 
                       ,mesh->numelNov));
        
        if(rCell < conv) break;
        if(!mpiVar.myId ){
          printf("it: %8d %.6e\n",i,rCell/rCell0);  
          if(opt.fItPlot)  
            fprintf(opt.fileItPlot[FITPLOTD1]
                   ,"%9d %.6e %0.6e\n",i,rCell/rCell0,rCell);
        }
/*...................................................................*/

/*...*/
#ifdef _DEBUG_
        testeSist(sistEqD1->ia       ,sistEqD1->ja
                 ,sistEqD1->au       ,sistEqD1->ad
                 ,sistEqD1->al       ,sistEqD1->b
                 ,sistEqD1->neq      ,sistEqD1->unsym);
#endif
/*...................................................................*/

/*...*/
        tm.solvD1 = getTimeC() - tm.solvD1;
        solverC(&m               
               ,sistEqD1->neq       ,sistEqD1->neqNov  
               ,sistEqD1->nad       ,sistEqD1->nadr
               ,sistEqD1->ia        ,sistEqD1->ja  
               ,sistEqD1->al        ,sistEqD1->ad,sistEqD1->au
               ,sistEqD1->b         ,sistEqD1->x
               ,&sistEqD1->iNeq
               ,solvD1->tol         ,solvD1->maxIt     
               ,sistEqD1->storage   ,solvD1->solver
               ,solvD1->fileSolv    ,solvD1->log  
               ,false               ,false
               ,sistEqD1->unsym     ,false);  
        tm.solvD1 = getTimeC() - tm.solvD1;
/*...................................................................*/


/*... x -> uD1*/
        updateCellValue(mesh->elm.uD1 ,sistEqD1->x
                       ,sistEqD1->id  ,&sistEqD1->iNeq
                       ,mesh->numel   ,mesh->ndfD[0]
                       ,false         ,true);
/*...................................................................*/

/*... reconstruindo do gradiente*/
       tm.rcGradD1 = getTimeC() - tm.rcGradD1;
       rcGradU(&m                       ,loadsD1
               ,mesh->elm.node          ,mesh->elm.adj.nelcon
               ,mesh->elm.geom.cc       ,mesh->node.x   
               ,mesh->elm.nen           ,mesh->elm.adj.nViz 
               ,mesh->elm.geomType      ,mesh->elm.material.prop 
               ,mesh->elm.mat 
               ,mesh->elm.leastSquare   ,mesh->elm.leastSquareR
               ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
               ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea    
               ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
               ,mesh->elm.geom.vSkew    ,mesh->elm.geom.xmcc  
               ,mesh->elm.geom.dcca
               ,mesh->elm.faceRd1       ,mesh->elm.faceLoadD1    
               ,mesh->elm.uD1           ,mesh->elm.gradUd1                 
               ,mesh->node.uD1          ,sc.rcGrad
               ,mesh->maxNo             ,mesh->maxViz
               ,mesh->ndfD[0]           ,mesh->ndm
               ,&pMesh->iNo             ,&pMesh->iEl  
               ,mesh->numelNov          ,mesh->numel        
               ,mesh->nnodeNov          ,mesh->nnode); 
       tm.rcGradD1 = getTimeC() - tm.rcGradD1;
/*...................................................................*/

       if(opt.fItPlotRes){  

/*... interpolacao das variaveis da celulas para pos nos (Grad)*/
         interCellNode(&m                 ,loadsD1
                      ,mesh->node.gradUd1 ,mesh->elm.gradUd1 
                      ,mesh->elm.node     ,mesh->elm.geomType            
                      ,mesh->elm.geom.cc  ,mesh->node.x                  
                      ,mesh->elm.nen      ,mesh->elm.adj.nViz
                      ,mesh->elm.faceRd1  ,mesh->elm.faceLoadD1    
                      ,&pMesh->iNo          
                      ,mesh->numelNov     ,mesh->numel 
                      ,mesh->nnodeNov     ,mesh->nnode 
                      ,mesh->maxNo        ,mesh->maxViz   
                      ,mesh->ndm          ,1 
                      ,mesh->ndm
                      ,false              ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (uD1)*/
         interCellNode(&m               ,loadsD1
                      ,mesh->node.uD1   ,mesh->elm.uD1 
                      ,mesh->elm.node   ,mesh->elm.geomType
                      ,mesh->elm.geom.cc,mesh->node.x                  
                      ,mesh->elm.nen    ,mesh->elm.adj.nViz
                      ,mesh->elm.faceRd1,mesh->elm.faceLoadD1    
                      ,&pMesh->iNo             
                      ,mesh->numelNov   ,mesh->numel
                      ,mesh->nnodeNov   ,mesh->nnode 
                      ,mesh->maxNo      ,mesh->maxViz 
                      ,mesh->ndfD[0]    ,1
                      ,mesh->ndm
                      ,true             ,2);
/*...................................................................*/

/*... globalizacao das variaveis*/
/*... uD1(Node)*/
          dGlobalNode(&m                 ,pMesh
                     ,mesh0->node.uD1    ,mesh->node.uD1     
                     ,mesh->ndfD[0]      ,1               );
          
/*... gradUd1(Node)*/
          dGlobalNode(&m                 ,pMesh
                     ,mesh0->node.gradUd1,mesh->node.gradUd1     
                     ,mesh->ndm          ,1               );
/*... uD1(Cel)*/
          dGlobalCel(&m                  ,pMesh
                    ,mesh0->elm.uD1      ,mesh->elm.uD1
                    ,mesh->numelNov 
                    ,mesh->ndfD[0]      ,1);
/*... gradUd1(Cel)*/
          dGlobalCel(&m                  ,pMesh
                    ,mesh0->elm.gradUd1  ,mesh->elm.gradUd1
                    ,mesh->numelNov 
                    ,mesh->ndm           ,1);
/*...................................................................*/

/*...*/
         if(!mpiVar.myId){
           fName(preName,sc.ddt.timeStep,i,15,&nameOut);
           strcpy(str1,"elD1");
           strcpy(str2,"noD1");
           strcpy(str3,"elGradD1");
           strcpy(str4,"noGradD1");
/*...*/
           wResVtkDif(&m                 ,mesh0->node.x      
                    ,mesh0->elm.node     ,mesh0->elm.mat    
                    ,mesh0->elm.nen      ,mesh0->elm.geomType
                    ,mesh0->elm.uD1      ,mesh0->node.uD1 
                    ,mesh0->elm.gradUd1  ,mesh0->node.gradUd1 
                    ,mesh0->nnode        ,mesh0->numel  
                    ,mesh0->ndm          ,mesh0->maxNo 
                    ,mesh0->numat        ,mesh0->ndfD[0]
                    ,str1                ,str2         
                    ,str3                ,str4         
                    ,nameOut             ,opt.bVtk    
                    ,sc.ddt              ,fileOut);
         }
/*...................................................................*/
       }
/*...................................................................*/
     }
/*...*/
     tm.solvEdoD1    = getTimeC() - tm.solvEdoD1;
/*...................................................................*/
     if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }
/*===================================================================*/

/*===================================================================*
 * macro: pcoo : escreve a matriz de coeficientes no formato COO
 *===================================================================*/
    else if((!strcmp(word,macro[7]))){
      if( mpiVar.nPrcs == 1 ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      
        fName(preName,0,0,12,&nameOut);
/*... matriz A*/
        writeCoo(&m,sistEqD1->ia,sistEqD1->ja,sistEqD1->neq
                ,sistEqD1->au   ,sistEqD1->ad,sistEqD1->al        
                ,sistEqD1->nad  ,sistEqD1->storage
                ,sistEqD1->unsym,false 
                ,nameOut);
/*...................................................................*/

/*... vetor de forcas b*/      
        fName(preName,0,0,14,&nameOut);
        for(int i=0;i<sistEqD1->neq;i++)
          sistEqD1->b[i] /= sistEqD1->ad[i];

        writeCooB(sistEqD1->b,sistEqD1->neq,nameOut);
/*...................................................................*/
        printf("%s\n\n",DIF);
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: puD1 : escreve os arquivos dos resultados da uD1          
 *===================================================================*/
    else if((!strcmp(word,macro[8]))){
/*... reconstruindo do gradiente*/
      tm.rcGradD1 = getTimeC() - tm.rcGradD1;
      rcGradU(&m                      ,loadsD1
             ,mesh->elm.node          ,mesh->elm.adj.nelcon
             ,mesh->elm.geom.cc       ,mesh->node.x   
             ,mesh->elm.nen           ,mesh->elm.adj.nViz 
             ,mesh->elm.geomType      ,mesh->elm.material.prop 
             ,mesh->elm.mat 
             ,mesh->elm.leastSquare   ,mesh->elm.leastSquareR
             ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
             ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea    
             ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
             ,mesh->elm.geom.vSkew    ,mesh->elm.geom.xmcc 
             ,mesh->elm.geom.dcca
             ,mesh->elm.faceRd1       ,mesh->elm.faceLoadD1    
             ,mesh->elm.uD1           ,mesh->elm.gradUd1                 
             ,mesh->node.uD1          ,sc.rcGrad
             ,mesh->maxNo             ,mesh->maxViz
             ,mesh->ndfD[0]           ,mesh->ndm       
             ,&pMesh->iNo             ,&pMesh->iEl  
             ,mesh->numelNov          ,mesh->numel        
             ,mesh->nnodeNov          ,mesh->nnode);  
      tm.rcGradD1 = getTimeC() - tm.rcGradD1;
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (Grad)*/
      interCellNode(&m                 ,loadsD1
                     ,mesh->node.gradUd1 ,mesh->elm.gradUd1 
                     ,mesh->elm.node     ,mesh->elm.geomType            
                     ,mesh->elm.geom.cc  ,mesh->node.x                  
                     ,mesh->elm.nen      ,mesh->elm.adj.nViz
                     ,mesh->elm.faceRd1  ,mesh->elm.faceLoadD1  
                     ,&pMesh->iNo          
                     ,mesh->numelNov     ,mesh->numel        
                     ,mesh->nnodeNov     ,mesh->nnode 
                     ,mesh->maxNo        ,mesh->maxViz   
                     ,mesh->ndm          ,1
                     ,mesh->ndm      
                     ,false              ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (uD1)*/
      interCellNode(&m               ,loadsD1
                    ,mesh->node.uD1    ,mesh->elm.uD1 
                    ,mesh->elm.node    ,mesh->elm.geomType
                    ,mesh->elm.geom.cc ,mesh->node.x                  
                    ,mesh->elm.nen     ,mesh->elm.adj.nViz
                    ,mesh->elm.faceRd1 ,mesh->elm.faceLoadD1   
                    ,&pMesh->iNo          
                    ,mesh->numelNov    ,mesh->numel        
                    ,mesh->nnodeNov    ,mesh->nnode 
                    ,mesh->maxNo       ,mesh->maxViz 
                    ,mesh->ndfD[0]     ,1 
                    ,mesh->ndm
                    ,true              ,2);
/*...................................................................*/

/*... globalizacao das variaveis*/
/*... uD1(Node)*/
      dGlobalNode(&m                 ,pMesh
                 ,mesh0->node.uD1    ,mesh->node.uD1     
                 ,mesh->ndfD[0]      ,1               );
          
/*... gradUd1(Node)*/
      dGlobalNode(&m                 ,pMesh
                 ,mesh0->node.gradUd1,mesh->node.gradUd1     
                 ,mesh->ndm          ,1               );
/*... uD1(Cel)*/
      dGlobalCel(&m                  ,pMesh
                ,mesh0->elm.uD1      ,mesh->elm.uD1
                ,mesh->numelNov 
                ,mesh->ndfD[0]      ,1);
/*... gradUd1(Cel)*/
      dGlobalCel(&m                  ,pMesh
                ,mesh0->elm.gradUd1  ,mesh->elm.gradUd1
                ,mesh->numelNov 
                ,mesh->ndm           ,1);
/*...................................................................*/

/*...*/
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
        fName(preName,sc.ddt.timeStep,0,8,&nameOut);

        strcpy(str1,"elD1");
        strcpy(str2,"noD1");
        strcpy(str3,"elGradD1");
        strcpy(str4,"noGradD1");
/*...*/
        wResVtkDif(&m                 ,mesh0->node.x      
                  ,mesh0->elm.node    ,mesh0->elm.mat    
                  ,mesh0->elm.nen     ,mesh0->elm.geomType
                  ,mesh0->elm.uD1     ,mesh0->node.uD1 
                  ,mesh0->elm.gradUd1 ,mesh0->node.gradUd1  
                  ,mesh0->nnode       ,mesh0->numel  
                  ,mesh0->ndm         ,mesh0->maxNo 
                  ,mesh0->numat       ,mesh0->ndfD[0]
                  ,str1               ,str2         
                  ,str3               ,str4         
                  ,nameOut            ,opt.bVtk    
                  ,sc.ddt             ,fileOut);
/*...................................................................*/
        printf("%s\n\n",DIF);
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: nlItD1: configura das iteracoes nao lineares             
 *===================================================================*/
    else if((!strcmp(word,macro[9]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
      fscanf(fileIn,"%d",&sc.nlD1.maxIt);
      fscanf(fileIn,"%lf",&sc.nlD1.tol);
      if(!mpiVar.myId ){
        printf("MaxIt: %d\n",sc.nlD1.maxIt);
        printf("Tol  : %e\n",sc.nlD1.tol);
      }
      readMacro(fileIn,word,false);
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: pD1CellCsv:imprime os resultados no formato csv                  
 *===================================================================*/
    else if((!strcmp(word,macro[10]))){
/*... globalizacao das variaveis*/
/*... uD1(Cel)*/
      dGlobalCel(&m                  ,pMesh
                ,mesh0->elm.uD1      ,mesh->elm.uD1
                ,mesh->numelNov 
                ,mesh->ndfD[0]      ,1);
/*... gradUd1(Cel)*/
      dGlobalCel(&m                  ,pMesh
                ,mesh0->elm.gradUd1  ,mesh->elm.gradUd1
                ,mesh->numelNov 
                ,mesh->ndm           ,1);
        
/*...*/
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      
/*...*/
        strcpy(auxName,preName);
        strcat(auxName,"_D1_cell_");
        fName(auxName,sc.ddt.timeStep,0,16,&nameOut);
        fileOut = openFile(nameOut,"w");
/*...*/
        writeCsvCell(mesh0->elm.uD1    ,mesh0->elm.gradUd1
                    ,mesh0->elm.geom.cc                  
                    ,mesh0->numel      ,mesh0->ndfD[0]
                    ,mesh0->ndm        ,fileOut);
/*...*/
        fclose(fileOut);
/*...................................................................*/
        printf("%s\n\n",DIF);
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: pD1CsvNode:imprime os resultados no formato csv                  
 *===================================================================*/
    else if((!strcmp(word,macro[11]))){
/*... globalizacao das variaveis*/
/*... uD1(Node)*/
      dGlobalNode(&m                 ,pMesh
                 ,mesh0->node.uD1    ,mesh->node.uD1     
                 ,mesh->ndfD[0]      ,1               );
          
/*... gradUd1(Node)*/
      dGlobalNode(&m                 ,pMesh
                 ,mesh0->node.gradUd1,mesh->node.gradUd1     
                 ,mesh->ndm          ,1               );

/*...................................................................*/

/*...*/
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
/*...*/
        strcpy(auxName,preName);
        strcat(auxName,"_D1_node_");
        fName(auxName,sc.ddt.timeStep,0,16,&nameOut);
        fileOut = openFile(nameOut,"w");
/*...*/
        writeCsvNode(mesh0->node.uD1    ,mesh0->node.gradUd1
                  ,mesh0->node.x                  
                  ,mesh0->nnode         ,mesh0->ndfD[0]
                  ,mesh0->ndm           ,fileOut);
/*...*/
        fclose(fileOut);
/*...................................................................*/
        printf("%s\n\n",DIF);
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: transient:configuracao da discretizacao temporal                 
 *===================================================================*/
    else if((!strcmp(word,macro[12]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
/*...*/
      sc.ddt.flag = true;
      readMacro(fileIn,word,false);
      if(!strcmp(word,"config:")){
/*... timer*/        
        readMacro(fileIn,word,false);
        setTransientScheme(word,&sc.ddt.type);
/*...*/ 
        fscanf(fileIn,"%lf",&sc.ddt.dt);
        fscanf(fileIn,"%lf",&sc.ddt.total);
/*...*/        
        if(!mpiVar.myId ) printf("dt(s)     : %lf\n",sc.ddt.dt);
        if(!mpiVar.myId ) printf("Total(s)  : %lf\n",sc.ddt.total);
      
        if(sc.ddt.type == EULER && !mpiVar.myId)     
          printf("ddtScheme : EULER\n");
        else if(sc.ddt.type == BACKWARD && !mpiVar.myId )     
          printf("ddtScheme : BACKWARD\n");

        sc.ddt.t        = 0.e0;
        sc.ddt.timeStep =    0;
      }
/*...................................................................*/

/*...*/
      flWord = true;
      kLoop  = 0;
      jLoop  = 0;
      do{
        readMacro(fileIn,word,false);
        strcpy(loopWord[kLoop],word);
        kLoop++;
        if(kLoop > 100){
          ERRO_GERAL(__FILE__,__func__,__LINE__,
                   "Numero de comandos na macro trasient execedido"); 
        }
      }while(strcmp(word,"endTransient"));
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: timeUpdate : macro de atualizaco do tempo                       
 *===================================================================*/
    else if((!strcmp(word,macro[13]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
/*...*/
      jLoop            = 0;
      sc.ddt.t        += sc.ddt.dt;
      sc.ddt.timeStep ++;
/*...................................................................*/
      
/*...*/
      if(sc.ddt.t > sc.ddt.total)
        flWord = false;
/*...................................................................*/
      
/*...*/
      else{
        if(!mpiVar.myId ) printf("t(s) = %lf\n",sc.ddt.t);
        if(!mpiVar.myId ) printf("step = %d\n" ,sc.ddt.timeStep);
      }
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: partd particionamento da malha                                   
 *===================================================================*/
    else if((!strcmp(word,macro[14]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
/*...*/
      pMesh->fPartMesh = true;
      readMacro(fileIn,word,false);
      if(!strcmp(word,"config:")){
/*... fPrintMesh*/        
        readMacro(fileIn,word,false);
        if(!strcmp(word,"true")){ 
          pMesh->fPrintMesh = true;
          if(!mpiVar.myId ) printf("fPrintMesh    : true\n");
        }

/*... fPrintMeshPart*/        
        readMacro(fileIn,word,false);
        if(!strcmp(word,"true")){ 
          pMesh->fPrintMeshPart = true;
          if(!mpiVar.myId ) printf("fPrintMeshPart: true\n");
        }
      }
/*...................................................................*/
/*...................................................................*/
       if(!mpiVar.myId ) printf("%s\n",DIF);
    }   
/*===================================================================*/
  }while(macroFlag && (!feof(fileIn)));

  mpiStop();
  return EXIT_SUCCESS;
}    
/*********************************************************************/
 
 
