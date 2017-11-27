/********************* includes padroes ******************************/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
/*********************************************************************/

/*********************************************************************/
#include<Adjcency.h>
#include<CellLoop.h>
#include<Coo.h>
#include<Diffusion.h>
#include<File.h>
#include<Memoria.h>
#include<Mesh.h>
#include<WriteVtk.h>
#include<WriteCsv.h>
#include<WriteLog.h>
#include<PartMesh.h>
#include<ParallelMpi.h>
#include<Properties.h>
#include<Prime.h>
#include<Sisteq.h>
#include<Solv.h>
#include<Simple.h>
#include<ReadFile.h>
#include<Reord.h>
#include<Transient.h>
#include<Transport.h>
#include<OpenMp.h>
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
  SistEq *sistEqD1 =NULL, *sistEqT1=NULL;
  SistEq *sistEqVel=NULL, *sistEqPres=NULL, *sistEqEnergy = NULL;
/*... metodo de acoplamento pressao-velocidade*/
  Simple *simple = NULL;
  Prime  *prime  = NULL;
/*... tubulence*/
  Turbulence turbModel;
/*...*/
  EnergyModel eModel;
  MassEqModel eMass;
  MomentumModel eMomentum;
/*... propriedade variaveis*/
  PropVar propVarFluid;

/*... solver*/
  INT nEqMax;
  Solv *solvD1=NULL,*solvT1=NULL,*solvVel=NULL,*solvPres=NULL,*solvEnergy=NULL;
  bool fSolvD1 = false, fSolvT1 = false;
  bool fSolvVel = false,fSolvPres = false, fSolvEnergy = false;
  bool fSolvSimple = false,fSolvPrime = false;
/*... reordenacao da malha*/
  Reord  *reordMesh=NULL;

/*... particionamento*/
  PartMesh *pMesh = NULL;

/*...*/
  char loopWord[100][MAX_LINE];
  unsigned short kLoop = 0 ,jLoop = 0,ndfVel;
  bool flWord=false;
  unsigned short nScheme,nOmp,nSistEq; 

/*... Estrutura de dados*/
  char strIa[MNOMEPONTEIRO],strJa[MNOMEPONTEIRO];
  char strA[MNOMEPONTEIRO],strAd[MNOMEPONTEIRO];

/*... arquivo*/
  char *nameIn=NULL,*nameOut=NULL,*preName=NULL,*auxName=NULL;
  FILE *fileIn=NULL,*fileOut=NULL,*fileLog=NULL;
  char str1[100],str2[100],str3[100],str4[100],str5[100],str6[100];
  FileOpt opt;

/*...*/
  Scheme sc;
/*...................................................................*/
 
/*... loop nas celulas*/
/*Lib lib;*/
  
/* ... macro camandos de leitura*/
  bool macroFlag; 
  char word[WORD_SIZE],str[WORD_SIZE];
  char macro[][WORD_SIZE] =
  {"help"        ,"mesh"         ,"stop"          /* 0, 1, 2*/
  ,"config"      ,""             ,""              /* 3, 4, 5*/
  ,"pgeo"        ,"pcoob"        ,"pcoo"          /* 6, 7, 8*/ 
  ,"presolvD1"   ,"presolvT1"    ,"openmp"        /* 9,10,11*/
  ,"solvD1"      ,""             ,"pD1"           /*12,13,14*/
  ,"nlItD1"      ,"pD1CsvCell"   ,"pD1CsvNode"    /*15,16,17*/
  ,"solvT1"      ,""             ,"pT1"           /*18,19,20*/
  ,"nlItT1"      ,"pT1CsvCell"   ,"pT1CsvNode"    /*21,22,23*/
  ,"setSolv"     ,"simple"       ,"setSimple"     /*24,25,26*/
  ,"transient"   ,"timeUpdate"   ,"partd"         /*27,28,29*/
  ,"advection"   ,"edp"          ,"diffusion"     /*30,31,32*/
  ,"pFluid"      ,"setPrintFluid" ,""             /*33,34,35*/
  ,"setPrime"    ,"prime"         ,"propVar"      /*36,37,38*/
  ,"gravity"     ,"model"         ,""        };   /*39,40,41*/
/* ..................................................................*/

/*... Memoria principal(valor padrao - bytes)*/
  nmax = 200000;
/* ..................................................................*/

/*...*/  
  thDynamic.pTh[0]      = PREREF;
  thDynamic.pTh[1]      = PREREF;
  thDynamic.pTh[2]      = PREREF;
/*....................................................................*/

/*...*/
  turbModel.fWall    = false;
  turbModel.wallType = STANDARDWALL;
  turbModel.fTurb    = false;
  turbModel.type     = SMAGORINSKY;
  turbModel.cs       = 0.2e0;
  turbModel.PrandltT = 0.5e0;
/*....................................................................*/

/*...*/
  eModel.fPresWork     = false;
  eModel.fDissipation  = false;
  eModel.fRes          = true;
  eModel.fTemperature  = false;
  eModel.fKelvin       = false;
/*...................................................................*/

/*...*/
  eMomentum.fRes             = true;
  eMomentum.fRhieChowInt     = false;
  eMomentum.fAbsultePressure = false;
/*...................................................................*/

/*...*/
  eMass.RhsDensity = false;
  eMass.LhsDensity = false;
/*...................................................................*/

/*... OpenMP*/
  ompVar.flag          = false;
  ompVar.nThreadsSolver = 1;
  ompVar.fSolver        = false;

  ompVar.nThreadsCell   = 1;
  ompVar.fCell          = false;

  ompVar.nThreadsUpdate = 1;
  ompVar.fUpdate        = false;
/* ..................................................................*/

/* ... opcoes de arquivos */                                           
  opt.bVtk          = false;
  opt.fCell         = false;
  opt.fNode         = true;
  opt.fItPlotRes    = false;
  opt.fItPlot       = false;
  opt.vel           = true;
  opt.pres          = true;
  opt.energy        = true;
  opt.gradVel       = false;
  opt.gradPres      = false;
  opt.gradEnergy    = false;
  opt.eddyViscosity = false;
  opt.densityFluid  = false;
  opt.specificHeat  = false;
  opt.dViscosity    = false;
  opt.tConductivity = false;
  opt.vorticity     = false;
  opt.bconditions   = true;
  opt.wallParameters= false;
  opt.stepPlotFluid[0] =  5;
  opt.stepPlotFluid[1] = opt.stepPlotFluid[0];
/* ..................................................................*/

/*... propriedades variaveis*/
  propVarFluid.fDensity           = false;
  propVarFluid.fSpecificHeat      = false;
  propVarFluid.fDynamicViscosity  = false;
  propVarFluid.fThermalCondutivty = false;
/*...................................................................*/

/*... Mpi*/
  mpiStart(&argc,argv);
/*...................................................................*/

/*... Time*/
  tm.adjcency          = 0.e0;
  tm.geom              = 0.e0;
  tm.leastSquareMatrix = 0.e0;
  tm.reord             = 0.e0;
/*... D1*/
  tm.solvD1            = 0.e0;
  tm.numeqD1           = 0.e0;
  tm.dataStructD1      = 0.e0;
  tm.CellPloadD1       = 0.e0;
  tm.CellTransientD1   = 0.e0;
  tm.systFormD1        = 0.e0;
  tm.rcGradD1          = 0.e0;
  tm.solvEdpD1         = 0.e0;
/*... T1*/
  tm.solvT1            = 0.e0;
  tm.numeqT1           = 0.e0;
  tm.dataStructT1      = 0.e0;
  tm.CellPloadT1       = 0.e0;
  tm.CellTransientT1   = 0.e0;
  tm.systFormT1        = 0.e0;
  tm.rcGradT1          = 0.e0;
  tm.solvEdpT1         = 0.e0;
/*... fluid*/
  tm.solvPres           = 0.e0;
  tm.solvVel            = 0.e0;
  tm.numeqPres          = 0.e0;
  tm.numeqVel           = 0.e0;
  tm.dataStructVel      = 0.e0;
  tm.dataStructPres     = 0.e0;
  tm.solvEdpFluid       = 0.e0;
  tm.cellPloadSimple    = 0.e0;
  tm.cellTransientSimple= 0.e0;
  tm.systFormPres       = 0.e0;
  tm.systFormVel        = 0.e0;
  tm.velExp             = 0.e0;
  tm.rcGradPres         = 0.e0;
  tm.rcGradVel          = 0.e0;

/*... Blas*/
  tm.matVecOverHeadMpi = 0.e0;
  tm.matVecSparse      = 0.e0;
  tm.dot               = 0.e0;
  tm.dotOverHeadMpi    = 0.e0;
/*... Iterativos  */
  tm.pcg               = 0.e0;
  tm.pbicgstab         = 0.e0;
  tm.gmres             = 0.e0;
  tm.minres            = 0.e0;
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
/*...*/
  tm.tempForEnergy     = 0.e0;
/*...*/
  tm.turbulence = 0.e0;

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
  sc.rcGrad       = RCGRADGAUSSN;
/*... D1,T1 */ 
  sc.nlD1.maxIt   = sc.nlT1.maxIt =   100; 
  sc.nlD1.tol     = sc.nlT1.tol   = 1.e-6; 
/*...................................................................*/

/*...*/
  sc.advT1.iCod1 = TVD;
  sc.advT1.iCod2 = VANLEERFACE;
  sc.diffT1.iCod = OVERRELAXED;
/*...................................................................*/

/*...*/
  sc.advVel.iCod1 = TVD;
  sc.advVel.iCod2 = VANLEERFACE;
  sc.advEnergy.iCod1 = TVD;
  sc.advEnergy.iCod2 = VANLEERFACE;
/*...................................................................*/

/*...*/  
  sc.diffVel.iCod  = OVERRELAXED;
  sc.diffPres.iCod = OVERRELAXED;
  sc.diffEnergy.iCod = OVERRELAXED;
/*...................................................................*/

/*...*/  
  sc.ddt.flag     = false;
  sc.ddt.fDynamic = false;
  sc.ddt.dt[0]    = 1.e0;
  sc.ddt.dt[1]    = 1.e0;
  sc.ddt.dt[2]    = 1.e0;
  sc.ddt.dtInicial= 1.e0;
  sc.ddt.total    = 1.e0;
  sc.ddt.timeStep = 0;
  sc.ddt.type     = BACKWARD;
  sc.ddt.iCod     = TDF;
/*...................................................................*/

/*...*/
  simple = NULL;
/*...................................................................*/

/*...*/  
  reordMesh = (Reord*) malloc(sizeof(Reord));
  ERRO_MALLOC(reordMesh,"reordMesh",__LINE__,__FILE__,__func__);
  reordMesh->flag = false; 
/* ..................................................................*/
    
/*... abrindo ar quivo de entrada*/ 
  nameIn = (char *) malloc(sizeof(char)*MAX_STR_LEN_IN);
    
  if( argc > 1){    
    strcpy(nameIn,argv[1]);  
  }
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
 * macro: help - leitura da malha e inicializa das estruturas        *
 * de resolucao do problema                                          * 
 *===================================================================*/
    if(!strcmp(word,macro[0])){
      help(fileIn);
    }
/*===================================================================*/

/*===================================================================*
 * macro: mesh - leitura da malha e inicializa das estruturas        *
 * de resolucao do problema                                          * 
 *===================================================================*/
    else if((!strcmp(word,macro[1]))){
      if(!mpiVar.myId) printf("%s\n",DIF);
      if(!mpiVar.myId) printf("%s\n",word); 
      if(!mpiVar.myId) printf("%s\n",DIF);
/*...*/
      if(mpiVar.nPrcs>1 && !pMesh->fPartMesh){
        printf("Erro: Falta a macro partd!!\n");
        mpiStop();
        exit(EXIT_FAILURE);
      }
/*...................................................................*/

/*... sistema de memoria*/
      initMem(&m,nmax,false);
/*... leitura da malha*/
      if(!mpiVar.myId)
        readFileFvMesh(&m          ,mesh0
                      ,propVarFluid,eModel
                      ,fileIn);
      mpiWait();
/*...................................................................*/
 
/*...*/
      if(!mpiVar.myId){
        printf("%s\n",DIF);
        usoMemoria(&m,"GB");
        printf("%s\n",DIF);
      }
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

/*... identifica parede impermevais*/
      if(mesh0->ndfF > 0 || mesh0->ndfFt > 0){
        wallFluid(mesh0->elm.faceRvel,mesh0->elm.adj.nelcon
                 ,mesh0->elm.adj.nViz   
                 ,mesh0->numel       ,mesh0->maxViz);  
/*... verifica se o dominio e aberto ou nao*/
        mesh0->fOpen = openDomain(loadsVel
                         , mesh0->elm.faceLoadVel , mesh0->elm.adj.nViz
                         , mesh0->numelNov ,mesh0->maxViz);
/*...*/
        thDynamic.fDensityRef = !mesh0->fOpen;
        thDynamic.fPresTh     = !mesh0->fOpen;
/*...................................................................*/   
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
                      ,loadsD1   ,loadsT1);
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

/*...*/
     if(mesh->ndfFt > 0){
/*...*/
        if(thDynamic.fDensityRef)
          specificMassRef(mesh->elm.densityFluid, mesh->elm.geom.volume                  
                       , mesh->elm.material.prop, mesh->elm.mat
                       , mesh->numel);
/*...................................................................*/

/*...*/
        if(thDynamic.fPresTh)
          initPresRef(mesh->elm.temp 
                    , mesh->elm.geom.volume  , thDynamic.pTh                  
                    , mesh->elm.material.prop, mesh->elm.mat
                    , mesh->numel            , eModel.fKelvin);
/*...................................................................*/

/*...*/
        mesh->mass[0] = totalMass(mesh->elm.densityFluid
                                , mesh->elm.geom.volume
                                , mesh->numel); 
        mesh->mass[1] =  mesh->mass[0];
/*....................................................................*/

/*... gera a pressao inicial hidrostatica*/
/*      hPres(mesh->elm.pressure0   , mesh->elm.pressure
            , mesh->elm.densityFluid, mesh->elm.geom.cc
            , gravity               , mesh->xRef
            , mesh->numel           , mesh->ndm );    */
/*...................................................................*/
      }
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
        globalMeshQuality(&mesh->mQuality,&mesh0->mQuality);
      }  
/*...................................................................*/
         
/*... reodenando as celulas para dimuincao da banda*/
      HccaAlloc(INT,&m,reordMesh->num,mesh->numel,"rNum" ,_AD_);
      if(!mpiVar.myId ) printf("%s\n",DIF);
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
        strcpy(str,"GB");
        memoriaTotal(str);
        usoMemoria(&m,str);
      }
/*...................................................................*/

    }   
/*===================================================================*/

/*===================================================================*
 * macro: stop : finalizacao do programa
 *===================================================================*/
    else if((!strcmp(word,macro[2]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word); 
        printf("%s\n\n",DIF);
      }
      tm.total = getTimeC() - tm.total;
/*... */
      fName(preName,mpiVar.nPrcs,mpiVar.myId,7,&nameOut);
      fileLog = openFile(nameOut,"w");
      writeLog(*mesh      ,sc
              ,solvD1     ,sistEqD1
              ,solvT1     ,sistEqT1
              ,solvVel    ,sistEqVel
              ,solvPres   ,sistEqPres
              ,tm       
              ,fSolvD1    ,fSolvT1     
              ,fSolvVel   ,fSolvPres 
              ,fSolvEnergy,turbModel.fTurb  
              ,ompVar
              ,nameIn     ,fileLog);
      fclose(fileLog);
/*...................................................................*/

/*... medias do tempo dos processos Mpi*/
      if(mpiVar.nPrcs > 1) {
        if(!mpiVar.myId){
          fName(preName,mpiVar.nPrcs,mpiVar.myId,60,&nameOut);
          fileLog = openFile(nameOut,"w");
        }
        writeLogMeanTime(*mesh0    ,sc
                        ,solvD1    ,sistEqD1
                        ,solvT1    ,sistEqT1
                        ,tm       
                        ,fSolvD1   ,fSolvT1     
                        ,nameIn    ,fileLog);
        
        if(!mpiVar.myId) fclose(fileLog);
      } 
/*...................................................................*/

/*... fechando o arquivo do log do solver linear Pres*/
      if(fSolvPres && solvPres->log && !mpiVar.myId)  
        fclose(solvPres->fileSolv);
/*... fechando o arquivo do log do solver linear Vel*/
      if(fSolvVel && solvVel->log && !mpiVar.myId)  
        fclose(solvVel->fileSolv);
/*... fechando o arquivo do log do solver linear T1*/
      if(fSolvD1 && solvD1->log && !mpiVar.myId)  
        fclose(solvD1->fileSolv);
/*... fechando o arquivo do log do solver linear T1*/
      if(fSolvT1 && solvT1->log && !mpiVar.myId)  
        fclose(solvT1->fileSolv);
/*...................................................................*/

/*... fechando o arquivo do log nao linear D1*/      
      if(fSolvD1 && opt.fItPlot && !mpiVar.myId)  
        fclose(opt.fileItPlot[FITPLOTD1]);
/*... fechando o arquivo do log nao linear T1*/      
      if(fSolvT1 && opt.fItPlot && !mpiVar.myId)  
        fclose(opt.fileItPlot[FITPLOTT1]);
/*... fechando o arquivo do log nao linear do simple */      
      if(fSolvSimple && opt.fItPlot && !mpiVar.myId)  
        fclose(opt.fileItPlot[FITPLOTSIMPLE]);
/*...................................................................*/
      finalizeMem(&m,false);
      macroFlag = false;
    }    
/*===================================================================*/

/*===================================================================*
 * macro: config : configuracao basica de excucao
 *===================================================================*/
    else if((!strcmp(word,macro[3]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
      
      config(&opt          ,reordMesh
            ,&sc.rcGrad    ,fileIn);
      
      if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: pgeo : escreve a geometria com os carregamentos
 *===================================================================*/
    else if((!strcmp(word,macro[6]))){
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
               ,mesh0->elm.faceRt1       ,mesh0->elm.faceLoadT1
               ,mesh0->elm.faceRvel      ,mesh0->elm.faceLoadVel  
               ,mesh0->elm.faceRenergy   ,mesh0->elm.faceLoadEnergy
               ,mesh0->nnode             ,mesh0->numel    
               ,mesh0->ndm               
               ,mesh0->maxNo             ,mesh0->maxViz
               ,mesh0->numat             
               ,mesh0->ndfD              ,mesh0->ndfT 
               ,mesh0->ndfF              ,mesh0->ndfFt
               ,nameOut                  ,opt.bVtk             
               ,fileOut);  
/*... face com cargas*/
       fName(preName,0,0,17,&nameOut);
       wGeoFaceVtk(&m                  ,mesh0->node.x        
             ,mesh0->elm.node          ,mesh0->elm.nen      
             ,mesh0->elm.geomType
             ,mesh0->elm.faceRd1       ,mesh0->elm.faceLoadD1
             ,mesh0->elm.faceRt1       ,mesh0->elm.faceLoadT1
             ,mesh0->elm.faceRvel      ,mesh0->elm.faceLoadVel  
             ,mesh0->elm.faceRenergy   ,mesh0->elm.faceLoadEnergy 
             ,mesh0->nnode             ,mesh0->numel    
             ,mesh0->ndm               
             ,mesh0->ndfD[0]           ,mesh0->ndfT[0]
             ,mesh0->ndfF              ,mesh0->ndfFt
             ,mesh0->maxViz            ,mesh0->maxNo
             ,nameOut                  ,opt.bVtk             
             ,fileOut);  
      }
    }   
/*===================================================================*/
   
/*===================================================================*
 * macro: pcoob : escreve a matriz de coeficientes no formato COO
 *===================================================================*/
    else if((!strcmp(word,macro[7]))){
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
 * macro: pcoo : escreve a matriz de coeficientes no formato COO
 *===================================================================*/
    else if((!strcmp(word,macro[8]))){
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
 * macro: presolvd1 : problema de difusao pura
 *===================================================================*/
    else if((!strcmp(word,macro[9]))){
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
        strcat(auxName,"_D1");
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
               ,sistEqD1->neq            ,"sistD1b0",_AD_);
      HccaAlloc(DOUBLE                   ,&m        ,sistEqD1->b 
               ,sistEqD1->neq            ,"sistD1b ",_AD_);
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
        strcpy(str,"GB");
        memoriaTotal(str);
        usoMemoria(&m,str);
      }
/*...................................................................*/
      if(!mpiVar.myId  ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: presolvt1 : problema de transporte  
 *===================================================================*/
    else if((!strcmp(word,macro[10]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
      

/*... inicializando a estrutura de equacoes do problema*/
      solvT1 = (Solv*) malloc(sizeof(Solv));
      if(solvT1 == NULL){
        printf("Erro ponteiro solvT1\n");
        exit(EXIT_FAILURE);
      }
      fSolvT1          = true;
      solvT1->solver   = PBICGSTAB;
      solvT1->tol      = smachn();
      solvT1->maxIt    = 50000;    
      solvT1->fileSolv = NULL;
      solvT1->log      = true;
      solvT1->flag     = true;
/*...................................................................*/

/*...*/
      if(solvT1->log && !mpiVar.myId){  
        strcpy(auxName,preName);
        strcat(auxName,"_T1");
        fName(auxName,mpiVar.nPrcs,0,11,&nameOut);
        solvT1->fileSolv = openFile(nameOut,"w");
      }
/*...................................................................*/

/*...*/
      if(opt.fItPlot && !mpiVar.myId){  
        strcpy(auxName,preName);
        strcat(auxName,"_T1");
        fName(auxName,mpiVar.nPrcs,0,10,&nameOut);
        opt.fileItPlot[FITPLOTT1] = openFile(nameOut,"w");
        fprintf(opt.fileItPlot[FITPLOTT1]
               ,"#T1\n#it ||b||/||b0|| ||b||\n");
      }
/*...................................................................*/

/*... inicializa a estrutura do solver*/
      sistEqT1 = (SistEq*) malloc(sizeof(SistEq));
      if(sistEqT1 == NULL){
        printf("Erro ponteiro sistEqT1\n");
        exit(EXIT_FAILURE);
      }
      sistEqT1->storage = CSRD;
      sistEqT1->unsym   = true; 
/*...................................................................*/

/*... config*/
      readMacro(fileIn,word,false);
      if(!strcmp(word,"config:")){
/*... solver*/        
        readMacro(fileIn,word,false);
        setSolver(word,&solvT1->solver); 
        
/*... DataStruct*/    
        readMacro(fileIn,word,false);
        setDataStruct(word,&sistEqT1->storage); 

/*... */        
        fscanf(fileIn,"%u" ,&solvT1->maxIt);
        fscanf(fileIn,"%lf",&solvT1->tol);

        if( solvT1->tol == 0.e0) 
          solvT1->tol = smachn();

        if(!mpiVar.myId ) printf("MaxIt     : %d\n",solvT1->maxIt);
        if(!mpiVar.myId ) printf("Tol       : %e\n",solvT1->tol);
      }
/*...................................................................*/

/*... numeracao das equacoes*/
      HccaAlloc(INT,&m,sistEqT1->id
               ,mesh->numel*mesh->ndfT[0]
               ,"sistT1id",_AD_);
      if(!mpiVar.myId){
        printf("%s\n",DIF);
        printf("Numerando as equacoes.\n");
      }
      tm.numeqT1 = getTimeC() - tm.numeqT1;
      sistEqT1->neq = numeq(sistEqT1->id       ,reordMesh->num
                           ,mesh->elm.faceRt1  ,mesh->elm.adj.nViz
                           ,mesh->numel        ,mesh->maxViz
                           ,mesh->ndfT[0]);
      tm.numeqT1 = getTimeC() - tm.numeqT1;
      if(!mpiVar.myId){
        printf("Equacoes numeradas.\n");
        printf("%s\n",DIF);
      }
/*...................................................................*/

/*...*/
      if( mpiVar.nPrcs > 1) {      
        tm.numeqT1 = getTimeC() - tm.numeqT1;
        sistEqT1->neqNov = countEq(reordMesh->num
                            ,mesh->elm.faceRt1  ,mesh->elm.adj.nViz
                            ,mesh->numelNov     ,mesh->maxViz
                            ,mesh->ndfT[0]);
        tm.numeqT1 = getTimeC() - tm.numeqT1;
      }
      else{
        sistEqT1->neqNov = sistEqT1->neq;
      }
/*...................................................................*/

/*...*/
      HccaAlloc(DOUBLE                   ,&m        ,sistEqT1->b0
               ,sistEqT1->neq            ,"sistT1b0",_AD_);
      HccaAlloc(DOUBLE                   ,&m        ,sistEqT1->b 
               ,sistEqT1->neq            ,"sistT1b ",_AD_);
      HccaAlloc(DOUBLE                   ,&m        ,sistEqT1->x 
               ,sistEqT1->neq            ,"sistT1x ",_AD_);
      zero(sistEqT1->b0,sistEqT1->neq    ,DOUBLEC);
      zero(sistEqT1->b ,sistEqT1->neq    ,DOUBLEC);
      zero(sistEqT1->x ,sistEqT1->neq    ,DOUBLEC);
/*...................................................................*/

/*... Estrutura de Dados*/
      strcpy(strIa,"iaT1");
      strcpy(strJa,"jaT1");
      strcpy(strAd,"adT1");
      strcpy(strA ,"aT1");
      if(!mpiVar.myId) printf("Montagem da estrura de dados esparsa.\n");
      tm.dataStructT1 = getTimeC() - tm.dataStructT1 ;
      dataStruct(&m,sistEqT1->id   ,reordMesh->num,mesh->elm.adj.nelcon
                ,mesh->elm.adj.nViz,mesh->numelNov,mesh->maxViz
                ,mesh->ndfT[0]     ,strIa         ,strJa
                ,strAd             ,strA          ,sistEqT1);
      tm.dataStructT1 = getTimeC() - tm.dataStructT1 ;
      if(!mpiVar.myId) printf("Estrutuda montada.\n");
/*...................................................................*/

/*... Openmp(T1)*/
      if (ompVar.fSolver) {
/*... dividindo a matriz*/
        strcpy(str1,"thBeginT1");
        strcpy(str2,"thEndT1");
        strcpy(str3,"thSizeT1");
        strcpy(str4,"thHeightT1");
        pMatrixSolverOmp(&m,sistEqT1,str1,str2,str3,str4);
/*...................................................................*/

/*alocando o buffer*/
        nEqMax = sistEqT1->neq;
        HccaAlloc(DOUBLE, &m, ompVar.buffer
                 ,nEqMax*ompVar.nThreadsSolver,"bufferOmp", false);
        zero(ompVar.buffer, nEqMax*ompVar.nThreadsSolver, DOUBLEC);
        sistEqT1->omp.thY = ompVar.buffer;
      }
/*...................................................................*/

/*... mapa de equacoes para comunicacao*/
//    if( mpiVar.nPrcs > 1) {    
//      front(&m,pMesh,sistEqT1,mesh->ndfT[0]);  
//    } 
/*...................................................................*/

/*... informacao da memoria total usada*/
      if(!mpiVar.myId  ) {
        strcpy(str,"GB");
        memoriaTotal(str);
        usoMemoria(&m,str);
      }
/*...................................................................*/
      if(!mpiVar.myId  ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
* macro: openmp: configuracao do openmp  
*===================================================================*/
    else if ((!strcmp(word, macro[11]))) {
/*... tecnica de adveccao*/
      readMacro(fileIn, word, false);
      printf("OpenMp:\n");
      nOmp = (short)atol(word);
      ompVar.flag = true;
      do {
        readMacro(fileIn, word, false);
/*... solver*/
        if (!strcmp(word, "solver") || !strcmp(word, "Solver")) {
          readMacro(fileIn, word, false);
/*... codigo da da funcao limitadora de fluxo*/
          ompVar.nThreadsSolver = (short)atol(word);
          ompVar.fSolver        = true;
/*...................................................................*/

/*...*/    
          printf("Solver nThreads: %d\n", ompVar.nThreadsSolver);
/*...................................................................*/
          nOmp--;
        }
/*...................................................................*/

/*... cell*/
        else if (!strcmp(word, "Cell") || !strcmp(word, "cell")) {
          readMacro(fileIn, word, false);
/*...*/
          ompVar.nThreadsCell = (short)atol(word);
          ompVar.fCell = true;
/*...................................................................*/

/*...*/       
          printf("Cell nThreads: %d\n", ompVar.nThreadsCell);
/*...................................................................*/
          nOmp--;
        }
/*...................................................................*/

/*... update*/
        else if (!strcmp(word, "Update") || !strcmp(word, "update")) {
          readMacro(fileIn, word, false);
/*...*/
          ompVar.nThreadsUpdate = (short)atol(word);
          ompVar.fUpdate = true;
/*...................................................................*/

/*...*/       
          printf("Update nThreads: %d\n", ompVar.nThreadsUpdate);
/*...................................................................*/
          nOmp--;
        }
/*...................................................................*/
      } while (nOmp);

      openMpCheck(ompVar.flag);
    }
/*===================================================================*/

/*===================================================================*
 * macro: solvd1: problema de difusao pura
 *===================================================================*/
    else if((!strcmp(word,macro[12]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
        printf("%s\n",DIF);
      }
      mpiWait();
      tm.solvEdpD1    = getTimeC() - tm.solvEdpD1;
/*...*/
      if(solvD1 == NULL){
        printf("Estrutara de dados nao montada para o solvD1!!!\n");
        exit(EXIT_FAILURE);
      }
/*...................................................................*/
     
/*...*/
      diffusion(&m         ,loadsD1
               ,mesh0      ,mesh           ,sistEqD1
               ,solvD1     ,sc             ,pMesh
               ,opt        ,preName        ,nameOut
               ,fileOut);
/*...................................................................*/

/*...*/
     tm.solvEdpD1    = getTimeC() - tm.solvEdpD1;
/*...................................................................*/
     if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }
/*===================================================================*/

/*===================================================================*
 * macro: puD1 : escreve os arquivos dos resultados da uD1          
 *===================================================================*/
    else if((!strcmp(word,macro[14]))){
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
             ,mesh->elm.geom.vSkew     
             ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc    
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
                     ,mesh->elm.geom.xm
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
                    ,mesh->elm.geom.xm                  
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
    else if((!strcmp(word,macro[15]))){
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
    else if((!strcmp(word,macro[16]))){
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
    else if((!strcmp(word,macro[17]))){
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
        writeCsvNode(mesh0->node.uD1,mesh0->node.gradUd1
                    ,mesh0->node.x                  
                    ,mesh0->nnode   ,mesh0->ndfD[0]
                    ,mesh0->ndm     ,fileOut);
/*...*/
        fclose(fileOut);
/*...................................................................*/
        printf("%s\n\n",DIF);
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: solvt1: problema de transporte  
 *===================================================================*/
    else if((!strcmp(word,macro[18]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
        printf("%s\n",DIF);
      }
      mpiWait();
      tm.solvEdpT1    = getTimeC() - tm.solvEdpT1;
/*...*/
      if(solvT1 == NULL){
        printf("Estrutara de dados nao montada para o solvT1!!!\n");
        exit(EXIT_FAILURE);
      }
/*...................................................................*/
     
/*...*/
      transport(&m         ,loadsT1
               ,mesh0      ,mesh           ,sistEqT1
               ,solvT1     ,sc             ,pMesh
               ,opt        ,preName        ,nameOut
               ,fileOut);
/*...................................................................*/

/*...*/
     tm.solvEdpT1    = getTimeC() - tm.solvEdpT1;
/*...................................................................*/
     if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }
/*===================================================================*/

/*===================================================================*
 * macro: puT1 : escreve os arquivos dos resultados da uT1          
 *===================================================================*/
    else if((!strcmp(word,macro[20]))){
/*... reconstruindo do gradiente*/
      tm.rcGradT1 = getTimeC() - tm.rcGradT1;
      rcGradU(&m                      ,loadsT1
             ,mesh->elm.node          ,mesh->elm.adj.nelcon
             ,mesh->elm.geom.cc       ,mesh->node.x   
             ,mesh->elm.nen           ,mesh->elm.adj.nViz 
             ,mesh->elm.geomType      ,mesh->elm.material.prop 
             ,mesh->elm.mat 
             ,mesh->elm.leastSquare   ,mesh->elm.leastSquareR
             ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
             ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea    
             ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
             ,mesh->elm.geom.vSkew     
             ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc    
             ,mesh->elm.geom.dcca
             ,mesh->elm.faceRt1       ,mesh->elm.faceLoadT1    
             ,mesh->elm.uT1           ,mesh->elm.gradUt1                 
             ,mesh->node.uT1          ,sc.rcGrad
             ,mesh->maxNo             ,mesh->maxViz
             ,mesh->ndfT[0]           ,mesh->ndm       
             ,&pMesh->iNo             ,&pMesh->iEl  
             ,mesh->numelNov          ,mesh->numel        
             ,mesh->nnodeNov          ,mesh->nnode);  
      tm.rcGradT1 = getTimeC() - tm.rcGradT1;
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (Grad)*/
      interCellNode(&m                   ,loadsT1
                     ,mesh->node.gradUt1 ,mesh->elm.gradUt1 
                     ,mesh->elm.node     ,mesh->elm.geomType            
                     ,mesh->elm.geom.cc  ,mesh->node.x  
                     ,mesh->elm.geom.xm
                     ,mesh->elm.nen      ,mesh->elm.adj.nViz
                     ,mesh->elm.faceRt1  ,mesh->elm.faceLoadT1  
                     ,&pMesh->iNo          
                     ,mesh->numelNov     ,mesh->numel        
                     ,mesh->nnodeNov     ,mesh->nnode 
                     ,mesh->maxNo        ,mesh->maxViz   
                     ,mesh->ndm          ,1
                     ,mesh->ndm      
                     ,false              ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (vel)*/
      interCellNode(&m                 ,loadsT1
                   ,mesh->node.vel     ,mesh->elm.vel        
                   ,mesh->elm.node     ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc  ,mesh->node.x  
                   ,mesh->elm.geom.xm
                   ,mesh->elm.nen      ,mesh->elm.adj.nViz
                   ,mesh->elm.faceRvel ,mesh->elm.faceLoadVel  
                   ,&pMesh->iNo          
                   ,mesh->numelNov     ,mesh->numel        
                   ,mesh->nnodeNov     ,mesh->nnode 
                   ,mesh->maxNo        ,mesh->maxViz   
                   ,mesh->ndm          ,1
                   ,mesh->ndm      
                   ,false              ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (uT1)*/
      interCellNode(&m                 ,loadsT1
                    ,mesh->node.uT1    ,mesh->elm.uT1 
                    ,mesh->elm.node    ,mesh->elm.geomType
                    ,mesh->elm.geom.cc ,mesh->node.x
                    ,mesh->elm.geom.xm                  
                    ,mesh->elm.nen     ,mesh->elm.adj.nViz
                    ,mesh->elm.faceRt1 ,mesh->elm.faceLoadT1   
                    ,&pMesh->iNo          
                    ,mesh->numelNov    ,mesh->numel        
                    ,mesh->nnodeNov    ,mesh->nnode 
                    ,mesh->maxNo       ,mesh->maxViz 
                    ,mesh->ndfT[0]     ,1 
                    ,mesh->ndm
                    ,true              ,2);
/*...................................................................*/

/*... globalizacao das variaveis*/
/*... uT1(Node)*/
      dGlobalNode(&m                 ,pMesh
                 ,mesh0->node.uT1    ,mesh->node.uT1     
                 ,mesh->ndfT[0]      ,1               );
          
/*... gradUt1(Node)*/
      dGlobalNode(&m                 ,pMesh
                 ,mesh0->node.gradUt1,mesh->node.gradUt1     
                 ,mesh->ndm          ,1               );
/*... vel(Node)*/
      dGlobalNode(&m                 ,pMesh
                 ,mesh0->node.vel    ,mesh->node.vel         
                 ,mesh->ndm          ,1               );
/*... uT1(Cel)*/
      dGlobalCel(&m                  ,pMesh
                ,mesh0->elm.uT1      ,mesh->elm.uT1
                ,mesh->numelNov 
                ,mesh->ndfT[0]      ,1);
/*... gradUt1(Cel)*/
      dGlobalCel(&m                  ,pMesh
                ,mesh0->elm.gradUt1  ,mesh->elm.gradUt1
                ,mesh->numelNov 
                ,mesh->ndm           ,1);
/*... vel(Cel)*/
      dGlobalCel(&m                  ,pMesh
                ,mesh0->elm.vel      ,mesh->elm.vel       
                ,mesh->numelNov 
                ,mesh->ndm           ,1);
/*...................................................................*/

/*...*/
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
        fName(preName,sc.ddt.timeStep,0,20,&nameOut);

        strcpy(str1,"elT1");
        strcpy(str2,"noT1");
        strcpy(str3,"elGradT1");
        strcpy(str4,"noGradT1");
        strcpy(str5,"elVel");
        strcpy(str6,"noVel");
/*...*/
        wResVtkTrans(&m                 ,mesh0->node.x      
                    ,mesh0->elm.node    ,mesh0->elm.mat    
                    ,mesh0->elm.nen     ,mesh0->elm.geomType
                    ,mesh0->elm.uT1     ,mesh0->node.uT1 
                    ,mesh0->elm.gradUt1 ,mesh0->node.gradUt1  
                    ,mesh0->elm.vel     ,mesh0->node.vel      
                    ,mesh0->nnode       ,mesh0->numel  
                    ,mesh0->ndm         ,mesh0->maxNo 
                    ,mesh0->numat       ,mesh0->ndfT[0]
                    ,str1               ,str2         
                    ,str3               ,str4         
                    ,str5               ,str6         
                    ,nameOut            ,opt.bVtk    
                    ,sc.ddt             ,fileOut);
/*...................................................................*/
        printf("%s\n\n",DIF);
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: nlItT1: configura das iteracoes nao lineares             
 *===================================================================*/
    else if((!strcmp(word,macro[21]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
      fscanf(fileIn,"%d",&sc.nlT1.maxIt);
      fscanf(fileIn,"%lf",&sc.nlT1.tol);
      if(!mpiVar.myId ){
        printf("MaxIt: %d\n",sc.nlT1.maxIt);
        printf("Tol  : %e\n",sc.nlT1.tol);
      }
      readMacro(fileIn,word,false);
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: pT1CellCsv:imprime os resultados no formato csv                  
 *===================================================================*/
    else if((!strcmp(word,macro[22]))){
/*... globalizacao das variaveis*/
/*... uT1(Cel)*/
      dGlobalCel(&m                  ,pMesh
                ,mesh0->elm.uT1      ,mesh->elm.uT1
                ,mesh->numelNov 
                ,mesh->ndfT[0]      ,1);
/*... gradUt1(Cel)*/
      dGlobalCel(&m                  ,pMesh
                ,mesh0->elm.gradUt1  ,mesh->elm.gradUt1
                ,mesh->numelNov 
                ,mesh->ndm           ,1);
        
/*...*/
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      
/*...*/
        strcpy(auxName,preName);
        strcat(auxName,"_T1_cell_");
        fName(auxName,sc.ddt.timeStep,0,16,&nameOut);
        fileOut = openFile(nameOut,"w");
/*...*/
        writeCsvCell(mesh0->elm.uT1    ,mesh0->elm.gradUt1
                    ,mesh0->elm.geom.cc                  
                    ,mesh0->numel      ,mesh0->ndfT[0]
                    ,mesh0->ndm        ,fileOut);
/*...*/
        fclose(fileOut);
/*...................................................................*/
        printf("%s\n\n",DIF);
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: pT1CsvNode:imprime os resultados no formato csv                  
 *===================================================================*/
    else if((!strcmp(word,macro[23]))){
/*... globalizacao das variaveis*/
/*... uT1(Node)*/
      dGlobalNode(&m                 ,pMesh
                 ,mesh0->node.uT1    ,mesh->node.uT1     
                 ,mesh->ndfT[0]      ,1               );
          
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
        strcat(auxName,"_T1_node_");
        fName(auxName,sc.ddt.timeStep,0,16,&nameOut);
        fileOut = openFile(nameOut,"w");
/*...*/
        writeCsvNode(mesh0->node.uT1    ,mesh0->node.gradUt1
                  ,mesh0->node.x                  
                  ,mesh0->nnode         ,mesh0->ndfT[0]
                  ,mesh0->ndm           ,fileOut);
/*...*/
        fclose(fileOut);
/*...................................................................*/
        printf("%s\n\n",DIF);
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: setSolv : escoamento de fluidos  
 *===================================================================*/
    else if((!strcmp(word,macro[24]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
/*... tecnica de adveccao*/
      readMacro(fileIn,word,false);
      nSistEq = (short) atol(word);
      do{   
        readMacro(fileIn,word,false);
/*... velocidade*/
        if(!strcmp(word,"Vel") || !strcmp(word, "vel")){
          nSistEq--;
/*... inicializando a estrutura de equacoes do problema (VELOCIDADE)*/
          solvVel = (Solv*) malloc(sizeof(Solv));
          if(solvVel == NULL){
            printf("Erro ponteiro solvVel\n");
            exit(EXIT_FAILURE);
          }
          fSolvVel          = true;
          solvVel->solver   = PBICGSTAB;
          solvVel->tol      = smachn();
          solvVel->maxIt    = 50000;    
          solvVel->fileSolv = NULL;
          solvVel->log      = true;
          solvVel->flag     = true;  
/*...................................................................*/

/*...*/
          if(solvVel->log && !mpiVar.myId){  
            strcpy(auxName,preName);
            strcat(auxName,"_fluid_vel");
            fName(auxName,mpiVar.nPrcs,0,11,&nameOut);
            solvVel->fileSolv = openFile(nameOut,"w");
          }  
/*...................................................................*/

/*... inicializa a estrutura do solver(VELOCIDADES)*/
          sistEqVel = (SistEq*) malloc(sizeof(SistEq));
          if(sistEqVel == NULL){
            printf("Erro ponteiro sistEqVel\n");
            exit(EXIT_FAILURE);
          }
          sistEqVel->unsym   = true;   
/*...................................................................*/

/*... solver*/        
          readMacro(fileIn,word,false);
          setSolverConfig(word,solvVel,fileIn);
/*...................................................................*/ 
        
/*... DataStruct*/    
          readMacro(fileIn,word,false);
          setDataStruct(word,&sistEqVel->storage);    
/*...................................................................*/

/*... numeracao das equacoes das velocidades*/
          HccaAlloc(INT,&m,sistEqVel->id
                   ,mesh->numel              
                   ,"sistVelid",_AD_);
          if(!mpiVar.myId){
            printf("%s\n",DIF);
            printf("Numerando as equacoes.\n");
          }
          tm.numeqVel = getTimeC() - tm.numeqVel;
          sistEqVel->neq = numEqV1(sistEqVel->id,reordMesh->num
                                  ,mesh->numel);
          tm.numeqVel = getTimeC() - tm.numeqVel;
          
          if(!mpiVar.myId){
            printf("Equacoes numeradas.\n");
            printf("%s\n",DIF);  
          }
/*...................................................................*/

/*...*/
          if( mpiVar.nPrcs > 1) {      
//          tm.numeqPres = getTimeC() - tm.numeqPres;
//      sistEqT1->neqNov = countEq(reordMesh->num
//                          ,mesh->elm.faceRt1  ,mesh->elm.adj.nViz
//                          ,mesh->numelNov     ,mesh->maxViz
//                          ,mesh->ndfT[0]);
//          tm.numeqPres = getTimeC() - tm.numeqPres;
          }
          else
            sistEqVel->neqNov = sistEqVel->neq;          
/*...................................................................*/

/*... velovidades*/
          ndfVel = max(mesh->ndfF - 1, mesh->ndfFt - 2);
/*...................................................................*/
          HccaAlloc(DOUBLE        ,&m      ,sistEqVel->b0
               ,sistEqVel->neq*ndfVel,"sistVelb0",_AD_);
          HccaAlloc(DOUBLE        ,&m     ,sistEqVel->b 
               ,sistEqVel->neq*ndfVel,"sistVelb ",_AD_);
          HccaAlloc(DOUBLE        ,&m     ,sistEqVel->x 
               ,sistEqVel->neq*ndfVel,"sistVelx ",_AD_);
          zero(sistEqVel->b0,sistEqVel->neq*ndfVel,DOUBLEC);
          zero(sistEqVel->b ,sistEqVel->neq*ndfVel,DOUBLEC);
          zero(sistEqVel->x ,sistEqVel->neq*ndfVel,DOUBLEC);    
/*...................................................................*/

/*... Estrutura de dados velocidades*/
          strcpy(strIa,"iaVel");
          strcpy(strJa,"jaVel");
          strcpy(strAd,"adVel");
          strcpy(strA ,"aVel");
          
          if(!mpiVar.myId ) printf("Vel:\n");
          if(!mpiVar.myId)
            printf("Montagem da estrura de dados esparsa.\n");
          
          tm.dataStructVel = getTimeC() - tm.dataStructVel;
          dataStructSimple(&m,sistEqVel->id,reordMesh->num
                ,mesh->elm.adj.nelcon
                ,mesh->elm.adj.nViz,mesh->numelNov,mesh->maxViz
                ,ndfVel            ,strIa         ,strJa
                ,strAd             ,strA          ,sistEqVel);
          tm.dataStructVel = getTimeC() - tm.dataStructVel;

          if(!mpiVar.myId) printf("Estrutuda montada.\n");  
/*...................................................................*/

/*... dividindo a matriz*/
/*... Openmp(Vel,Pres)*/
          if(ompVar.fSolver){
            strcpy(str1,"thBeginVel");
            strcpy(str2,"thEndVel");
            strcpy(str3,"thSizeVel");
            strcpy(str4,"thHeightVel");
            pMatrixSolverOmp(&m,sistEqVel,str1,str2,str3,str4);
          }
/*...................................................................*/
        }
/*...................................................................*/

/*...*/
        else if(!strcmp(word,"Pres") || !strcmp(word, "pres")){ 
          nSistEq--;
/*... inicializando a estrutura de equacoes do problema (PRESSAO)*/
          solvPres = (Solv*) malloc(sizeof(Solv));
          if(solvPres == NULL){
            printf("Erro ponteiro solvPres\n");
            exit(EXIT_FAILURE);
          }
          fSolvPres          = true;
          solvPres->solver   = PCG;
          solvPres->tol      = smachn();
          solvPres->maxIt    = 50000;    
          solvPres->fileSolv = NULL;
          solvPres->log      = true;
          solvPres->flag     = true;  
/*...................................................................*/

/*...*/
          if(solvPres->log && !mpiVar.myId){  
            strcpy(auxName,preName);
            strcat(auxName,"_fluid_pres");
            fName(auxName,mpiVar.nPrcs,0,11,&nameOut);
            solvPres->fileSolv = openFile(nameOut,"w");
          }  
/*...................................................................*/

/*... inicializa a estrutura do solver(PRESSAO)*/
          sistEqPres = (SistEq*) malloc(sizeof(SistEq));
          if(sistEqPres == NULL){
            printf("Erro ponteiro sistEqPres\n");
            exit(EXIT_FAILURE);
          }
          sistEqPres->unsym   = false;   
/*...................................................................*/

/*... solver*/        
          readMacro(fileIn,word,false);
          setSolverConfig(word,solvPres,fileIn); 
/*...................................................................*/
        
/*... DataStruct*/    
          readMacro(fileIn,word,false);
          setDataStruct(word,&sistEqPres->storage);    
/*...................................................................*/

/*... numeracao das equacoes das pressoes*/
          HccaAlloc(INT,&m,sistEqPres->id
                   ,mesh->numel              
                   ,"sistPresid",_AD_);
          if(!mpiVar.myId){
            printf("%s\n",DIF);
            printf("Numerando as equacoes.\n");
          }
          tm.numeqPres = getTimeC() - tm.numeqPres;
          sistEqPres->neq = numEqV2(sistEqPres->id ,reordMesh->num
                               ,mesh->elm.faceRpres,mesh->elm.adj.nViz
                               ,mesh->numel        ,mesh->maxViz);
          tm.numeqPres = getTimeC() - tm.numeqPres;
          if(!mpiVar.myId){
            printf("Equacoes numeradas.\n");
            printf("%s\n",DIF);
          }  
/*...................................................................*/

/*...*/
          if( mpiVar.nPrcs > 1) {      
//          tm.numeqPres = getTimeC() - tm.numeqPres;
//      sistEqT1->neqNov = countEq(reordMesh->num
//                          ,mesh->elm.faceRt1  ,mesh->elm.adj.nViz
//                          ,mesh->numelNov     ,mesh->maxViz
//                          ,mesh->ndfT[0]);
//          tm.numeqPres = getTimeC() - tm.numeqPres;
          }
          else
            sistEqPres->neqNov = sistEqPres->neq;          
/*...................................................................*/

/*... pressoes*/
          HccaAlloc(DOUBLE     ,&m        ,sistEqPres->b0
               ,sistEqPres->neq,"sistPresb0",_AD_);
          HccaAlloc(DOUBLE     ,&m        ,sistEqPres->b 
               ,sistEqPres->neq,"sistPresb ",_AD_);
          HccaAlloc(DOUBLE     ,&m        ,sistEqPres->x 
               ,sistEqPres->neq,"sistPresx ",_AD_);
          zero(sistEqPres->b0,sistEqPres->neq    ,DOUBLEC);
          zero(sistEqPres->b ,sistEqPres->neq    ,DOUBLEC);
          zero(sistEqPres->x ,sistEqPres->neq    ,DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados pressoes*/
          strcpy(strIa,"iaPres");
          strcpy(strJa,"japres");
          strcpy(strAd,"adPres");
          strcpy(strA ,"aPres");
          
          if(!mpiVar.myId ) printf("Pres:\n");
          if(!mpiVar.myId) 
            printf("Montagem da estrura de dados esparsa.\n");
          
          tm.dataStructPres = getTimeC() - tm.dataStructPres;
          dataStruct(&m,sistEqPres->id ,reordMesh->num,mesh->elm.adj.nelcon
                ,mesh->elm.adj.nViz,mesh->numelNov,mesh->maxViz
                ,1                 ,strIa         ,strJa
                ,strAd             ,strA          ,sistEqPres);
          tm.dataStructPres = getTimeC() - tm.dataStructPres;
          if(!mpiVar.myId) printf("Estrutuda montada.\n");
/*...................................................................*/

/*... Openmp(Vel,Pres)*/
          if(ompVar.fSolver){
/*... dividindo a matriz*/
            strcpy(str1,"thBeginPres");
            strcpy(str2,"thEndPres");
            strcpy(str3,"thSizePres");
            strcpy(str4,"thHeightPres");
            pMatrixSolverOmp(&m,sistEqPres,str1,str2,str3,str4);
          }
/*...................................................................*/
        }
/*...................................................................*/

/*... energy*/
        else if (!strcmp(word, "Energy") || !strcmp(word, "energy")) {
          nSistEq--;
/*... inicializando a estrutura de equacoes do problema (VELOCIDADE)*/
          solvEnergy = (Solv*)malloc(sizeof(Solv));
          if (solvEnergy == NULL) {
            printf("Erro ponteiro solvEnergy\n");
            exit(EXIT_FAILURE);
          }
          fSolvEnergy = true;
          solvEnergy->solver = PBICGSTAB;
          solvEnergy->tol = smachn();
          solvEnergy->maxIt = 50000;
          solvEnergy->fileSolv = NULL;
          solvEnergy->log = true;
          solvEnergy->flag = true;
/*...................................................................*/

/*...*/
          if (solvEnergy->log && !mpiVar.myId) {
            strcpy(auxName, preName);
            strcat(auxName, "_fluid_energy");
            fName(auxName, mpiVar.nPrcs, 0, 11, &nameOut);
            solvEnergy->fileSolv = openFile(nameOut, "w");
          }
/*...................................................................*/

/*... inicializa a estrutura do solver(Energia)*/
          sistEqEnergy = (SistEq*)malloc(sizeof(SistEq));
          if (sistEqEnergy == NULL) {
            printf("Erro ponteiro sistEqEnergia\n");
            exit(EXIT_FAILURE);
          }
          sistEqEnergy->unsym = true;
/*...................................................................*/

/*... solver*/
          readMacro(fileIn, word, false);
          setSolverConfig(word, solvEnergy, fileIn);
/*...................................................................*/

/*... DataStruct*/
          readMacro(fileIn, word, false);
          setDataStruct(word, &sistEqEnergy->storage);
/*...................................................................*/

/*... numeracao das equacoes das Energy*/
          HccaAlloc(INT           ,&m  ,sistEqEnergy->id
                   ,mesh->numel
                   ,"sistEnergyId",_AD_);
          if (!mpiVar.myId) { 
            printf("%s\n", DIF);
            printf("Numerando as equacoes.\n");
          }
          tm.numeqEnergy = getTimeC() - tm.numeqEnergy;
          sistEqEnergy->neq = numEqV1(sistEqEnergy->id, reordMesh->num
                                     ,mesh->numel);
          tm.numeqEnergy = getTimeC() - tm.numeqEnergy;

          if (!mpiVar.myId) {
            printf("Equacoes numeradas.\n");
            printf("%s\n", DIF);
          }
/*...................................................................*/

/*...*/
          if (mpiVar.nPrcs > 1) {
            //          tm.numeqPres = getTimeC() - tm.numeqPres;
            //      sistEqT1->neqNov = countEq(reordMesh->num
            //                          ,mesh->elm.faceRt1  ,mesh->elm.adj.nViz
            //                          ,mesh->numelNov     ,mesh->maxViz
            //                          ,mesh->ndfT[0]);
            //          tm.numeqPres = getTimeC() - tm.numeqPres;
          }
          else
            sistEqEnergy->neqNov = sistEqEnergy->neq;
/*...................................................................*/

/*... energia*/
          HccaAlloc(DOUBLE           ,&m           ,sistEqEnergy->b0
                   ,sistEqEnergy->neq,"sistEnergy0",_AD_);
          HccaAlloc(DOUBLE             , &m        ,sistEqEnergy->b
                   ,sistEqEnergy->neq  , "sistEnergyb ", _AD_);
          HccaAlloc(DOUBLE           ,&m            ,sistEqEnergy->x
                   ,sistEqEnergy->neq,"sistEnergyx ",_AD_);  
          zero(sistEqEnergy->b0,sistEqEnergy->neq, DOUBLEC);
          zero(sistEqEnergy->b ,sistEqEnergy->neq, DOUBLEC);
          zero(sistEqEnergy->x ,sistEqEnergy->neq, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados energia*/
          strcpy(strIa, "iaEnergy");
          strcpy(strJa, "jaEnergy");
          strcpy(strAd, "adEnergy");
          strcpy(strA , "aEnergy");

          if (!mpiVar.myId) { 
            printf("Energy:\n");
            printf("Montagem da estrura de dados esparsa.\n");
          }

          tm.dataStructEnergy = getTimeC() - tm.dataStructEnergy;
          dataStructSimple(&m,sistEqEnergy->id ,reordMesh->num
                          ,mesh->elm.adj.nelcon
                          ,mesh->elm.adj.nViz  ,mesh->numelNov,mesh->maxViz  
                          ,1                   ,strIa         ,strJa
                          ,strAd               ,strA          ,sistEqEnergy);
          tm.dataStructEnergy = getTimeC() - tm.dataStructEnergy;

          if (!mpiVar.myId) printf("Estrutuda montada.\n");
/*...................................................................*/

/*... dividindo a matriz*/
/*... Openmp(energy)*/
          if (ompVar.fSolver) {
            strcpy(str1, "thBeginEnergy");
            strcpy(str2, "thEndEnergy");
            strcpy(str3, "thSizeEnergy");
            strcpy(str4, "thHeightEnergy");
            pMatrixSolverOmp(&m,sistEqEnergy,str1,str2,str3,str4);
          }
/*...................................................................*/
        }
/*...................................................................*/

      }while(nSistEq);    

/*...*/
      if(opt.fItPlot && !mpiVar.myId){  
        strcpy(auxName,preName);
        fName(auxName,mpiVar.nPrcs,0,22,&nameOut);
        opt.fileItPlot[FITPLOTSIMPLE] = openFile(nameOut,"w");
        if(mesh->ndfF == 3)
          fprintf(opt.fileItPlot[FITPLOTSIMPLE]
          ,"#VelPres\n#it ||rU1||| ||rU2|| ||rMass||\n");
        else if(mesh->ndfF == 4)
          fprintf(opt.fileItPlot[FITPLOTSIMPLE]
          ,"#VelPres\n#it ||rU1|| ||rU2|| ||rU3|| ||rMass||\n");
        else if (mesh->ndfFt == 4)
          fprintf(opt.fileItPlot[FITPLOTSIMPLE]
            , "#VelPres\n#it ||rU1|| ||rU2|| ||rEnergy|| ||rMass||\n");
        else if (mesh->ndfFt == 5)
          fprintf(opt.fileItPlot[FITPLOTSIMPLE]
            , "#VelPres\n#it ||rU1|| ||rU2|| ||rU3|| ||rEnergy|| ||rMass||\n");
      }
/*...................................................................*/

/*... Openmp(Vel,Pres)*/
      if(ompVar.fSolver){

/*... alocando o buffer*/
        if(fSolvPres && fSolvVel && fSolvEnergy){
          nEqMax = max(sistEqPres->neqNov,sistEqVel->neqNov);
          nEqMax = max(sistEqEnergy->neqNov,nEqMax);
          HccaAlloc(DOUBLE, &m, ompVar.buffer
                   ,nEqMax*ompVar.nThreadsSolver,"bufferOmp",false);
          zero(ompVar.buffer,nEqMax*ompVar.nThreadsSolver,DOUBLEC);
          sistEqPres->omp.thY   = ompVar.buffer;
          sistEqVel->omp.thY    = ompVar.buffer;
          sistEqEnergy->omp.thY = ompVar.buffer;
        }
/*... alocando o buffer*/
        else if(fSolvPres && fSolvVel){
          nEqMax = max(sistEqPres->neqNov,sistEqVel->neqNov);
          HccaAlloc(DOUBLE, &m, ompVar.buffer
                   ,nEqMax*ompVar.nThreadsSolver,"bufferOmp",false);
          zero(ompVar.buffer,nEqMax*ompVar.nThreadsSolver,DOUBLEC);
          sistEqPres->omp.thY = ompVar.buffer;
          sistEqVel ->omp.thY = ompVar.buffer;
        }
/*... alocando o buffer*/
        else if(fSolvPres){
          nEqMax = sistEqPres->neqNov;
          HccaAlloc(DOUBLE, &m, ompVar.buffer
                   ,nEqMax*ompVar.nThreadsSolver,"bufferOmp",false);
          zero(ompVar.buffer,nEqMax*ompVar.nThreadsSolver,DOUBLEC);
          sistEqPres->omp.thY = ompVar.buffer;
        }
      }
/*...................................................................*/

/*... mapa de equacoes para comunicacao*/
//    if( mpiVar.nPrcs > 1) {    
//      front(&m,pMesh,sistEqT1,mesh->ndfT[0]);  
//    } 
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
 * macro: simple: escoamento de fluidos (SIMPLE)
 *===================================================================*/
    else if((!strcmp(word,macro[25]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
        printf("%s\n",DIF);
      }
      mpiWait();
      tm.solvEdpFluid = getTimeC() - tm.solvEdpFluid;
/*...*/
      if(solvVel == NULL || solvPres == NULL){
        printf("Estrutara de dados nao montada para o solvFluid!!!\n");
        exit(EXIT_FAILURE);
      }
/*...................................................................*/

/*...*/
     if(!fSolvSimple) {
        printf("Simple nao configurado ainda!!!\n");
        exit(EXIT_FAILURE);
     }  
/*...................................................................*/
     
/*...*/
      if(mesh->ndfF)
        simpleSolver3D(&m         
                      ,loadsVel   ,loadsPres 
                      ,mesh0      ,mesh           
                      ,sistEqVel  ,sistEqPres
                      ,solvVel    ,solvPres
                      ,simple
                      ,sc         ,pMesh
                      ,opt        ,preName        
                      ,nameOut    ,fileOut);
/*...................................................................*/

/*...*/
      else if(mesh->ndfFt)
        simpleSolverLm(&m          , propVarFluid
                     , loadsVel    , loadsPres 
                     , loadsEnergy , eModel
                     , eMass       , eMomentum
                     , turbModel   , &thDynamic   
                     , mesh0       , mesh
                     , sistEqVel   , sistEqPres
                     , sistEqEnergy  
                     , solvVel     , solvPres
                     , solvEnergy    
                     , simple        
                     , &sc         , pMesh
                     , opt         , preName
                     , nameOut     , fileOut);  
/*...................................................................*/

/*...*/
     tm.solvEdpFluid    = getTimeC() - tm.solvEdpFluid;
/*...................................................................*/
     if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }
/*===================================================================*/

/*===================================================================*
 * macro: preSimple: configuracoe do metodo simple
 *===================================================================*/
    else if((!strcmp(word,macro[26]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
        printf("%s\n",DIF);
      }
/*...*/
      simple = (Simple*) malloc(sizeof(Simple));
      if(simple == NULL){
        printf("Erro ponteiro simple\n");
        exit(EXIT_FAILURE);
      }
      fSolvSimple             = true;  
      simple->maxIt           = 1000;
      simple->alphaPres       = 0.3e0; 
      simple->alphaVel        = 0.7e0; 
      simple->type            = SIMPLE;
      simple->kZeroVel        = 4;
      simple->kZeroPres       = 0;
      simple->sPressure       = true;
      simple->faceInterpolVel = 1;
      simple->nNonOrth        = 0;
      simple->tolPres         = 1.e-06;
      simple->tolVel          = 1.e-06;
      if (mesh->ndfFt){
        simple->kZeroEnergy  = 0;
        simple->tolEnergy    = 1.e-06;
        simple->alphaEnergy  = 1.e0;
        simple->alphaDensity = 1.0e0; 
      }
      simple->pSimple         = 500;
/*...................................................................*/
      
/*...*/
      readMacro(fileIn,word,false);
      if(!strcmp(word,"config:")){
/*... timer*/        
        readMacro(fileIn,word,false);
/*... levemente compressivel*/       
        if (mesh->ndfFt)
          setSimpleLmScheme(word,simple,fileIn);
/*... imcompressivel*/
        else
          setSimpleScheme(word, simple, fileIn);
/*...*/        
        if(simple->type == SIMPLE && !mpiVar.myId)     
          printf("PRES-VEL  : SIMPLE\n");
        else if(simple->type == SIMPLEC && !mpiVar.myId )     
          printf("PRES-VEL  : SIMPLEC\n");

/*...*/        
        if(!mpiVar.myId ){ 
          printf("Maxit     : %d\n",simple->maxIt);
          printf("alphaPres : %lf\n",simple->alphaPres);
          printf("alphaVel  : %lf\n",simple->alphaVel);
          printf("tolPres   : %e\n",simple->tolPres);
          printf("tolVel    : %e\n",simple->tolVel);
          if(mesh->ndfFt)
            printf("tolEnergy : %e\n", simple->tolEnergy);
          printf("nNonOrth  : %d\n",simple->nNonOrth);
          printf("pSimple   : %d\n",simple->pSimple);
        }
      }
/*...................................................................*/

/*...*/
      HccaAlloc(DOUBLE     ,&m       ,simple->d
               ,mesh->numel*mesh->ndm,"dField" ,false);
      zero(simple->d    ,mesh->numel*mesh->ndm,DOUBLEC);

      HccaAlloc(DOUBLE     ,&m       ,simple->ePresC
               ,mesh->numel,"ePresC" ,false);
      zero(simple->ePresC,mesh->numel  ,DOUBLEC);

      HccaAlloc(DOUBLE     ,&m       ,simple->nPresC
               ,mesh->nnode,"nPresC" ,false);
      zero(simple->nPresC    ,mesh->numel  ,DOUBLEC);

      HccaAlloc(DOUBLE     ,&m      ,simple->eGradPresC
               ,mesh->numel*mesh->ndm,"eGradPresC",false);
      zero(simple->eGradPresC,mesh->numel*mesh->ndm  ,DOUBLEC);
      
      HccaAlloc(DOUBLE     ,&m       ,simple->ePresC1
               ,mesh->numel,"ePresC1",false);
      zero(simple->ePresC,mesh->numel  ,DOUBLEC);
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }
/*===================================================================*/

/*===================================================================*
 * macro: transient:configuracao da discretizacao temporal                 
 *===================================================================*/
    else if((!strcmp(word,macro[27]))){
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
        fscanf(fileIn,"%lf",&sc.ddt.dt[0]);
        fscanf(fileIn,"%lf",&sc.ddt.total);
/*...*/
        readMacro(fileIn,word,false);
        if(!strcmp(word,"dynamic"))     
          sc.ddt.fDynamic = true;
        if(sc.ddt.fDynamic){
          if(!mpiVar.myId) printf("dynamic : True\n");
        }
        else {
          if(!mpiVar.myId) printf("dynamic : False\n");
        }          
/*...*/        
        if(!mpiVar.myId ) printf("dt(s)     : %lf\n",sc.ddt.dt[0]);
        if(!mpiVar.myId ) printf("Total(s)  : %lf\n",sc.ddt.total);
      
        if(sc.ddt.type == EULER && !mpiVar.myId)     
          printf("ddtScheme : EULER\n");
        else if(sc.ddt.type == BACKWARD && !mpiVar.myId )     
          printf("ddtScheme : BACKWARD\n");

        sc.ddt.t         = 0.e0;
        sc.ddt.dtInicial = sc.ddt.dt[0];
        sc.ddt.dt[1]     = sc.ddt.dt[0];
        sc.ddt.dt[2]     = sc.ddt.dt[0];
        sc.ddt.timeStep  = 0;
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
    else if((!strcmp(word,macro[28]))){
      if(!mpiVar.myId ){
        printf("\n%s\n",DIF);
        printf("%s\n",word);
      }
/*...*/
      jLoop            = 0;
      sc.ddt.t        += sc.ddt.dt[0];
      sc.ddt.timeStep ++; 
/*...................................................................*/

/*...*/
      if(sc.ddt.t > sc.ddt.total + 0.1e0*sc.ddt.dt[0])
        flWord = false;  
/*    if(sc.ddt.t > sc.ddt.total)
        flWord = false;  */
/*...................................................................*/
      
/*...*/
      else{
        if(!mpiVar.myId ){
          printf("dt(n-2) = %lf\n",sc.ddt.dt[2]);
          printf("dt(n-1) = %lf\n",sc.ddt.dt[1]);
          printf("dt(n)   = %lf\n",sc.ddt.dt[0]);
          printf("t(s)    = %lf\n",sc.ddt.t);
          printf("step    = %d\n" ,sc.ddt.timeStep);
        } 
      }
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: partd particionamento da malha                                   
 *===================================================================*/
    else if((!strcmp(word,macro[29]))){
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
       if(!mpiVar.myId ) printf("%s\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: advection tecnica aplicada no termo advectivo          
 *===================================================================*/
    else if((!strcmp(word,macro[30]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
 /*... tecnica de adveccao*/
      readMacro(fileIn, word, false);
      nScheme = (short) atol(word);
      do {
        readMacro(fileIn, word, false);
 /*... velocidade*/
        if (!strcmp(word, "Vel") || !strcmp(word, "vel")) {
          printf("%s:\n", word);
          readMacro(fileIn, word, false);
 /*... codigo da da funcao limitadora de fluxo*/
          setAdvectionScheme(word, &sc.advVel,fileIn);
          nScheme--;
        }
 /*... T1*/
        else if (!strcmp(word, "T1") || !strcmp(word, "t1")) {
          printf("%s:\n", word);
          readMacro(fileIn, word, false);
 /*... codigo da da funcao limitadora de fluxo*/
          setAdvectionScheme(word, &sc.advT1,fileIn);
          nScheme--;
        }
 /*... Energy*/
        else if (!strcmp(word, "Energy") || !strcmp(word, "energy")) {
          printf("%s:\n", word);
          readMacro(fileIn, word, false);
/*... codigo da da funcao limitadora de fluxo*/
          setAdvectionScheme(word, &sc.advEnergy, fileIn);
          nScheme--;
        }
      } while (nScheme);
/*...................................................................*/
      if (!mpiVar.myId) printf("%s\n", DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: edo equacoes diferencias resolvidas                    
 *===================================================================*/
    else if((!strcmp(word,macro[31]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
      readEdo(mesh0,fileIn);
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: diffusion tecnica aplicada no termo difusivo           
 *===================================================================*/
    else if((!strcmp(word,macro[32]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
/*... tecnica de adveccao*/
      readMacro(fileIn,word,false);
      nScheme = (short) atol(word);
      do{ 
        readMacro(fileIn,word,false);
/*... velocidade*/
        if(!strcmp(word,"Vel") || !strcmp(word, "vel")){
          printf("%s:\n",word);
          readMacro(fileIn,word,false);
/*... codigo da da funcao limitadora de fluxo*/        
          setDiffusionScheme(word,&sc.diffVel.iCod);
          nScheme--;
        }
/*... Pressao*/
        else if(!strcmp(word,"Pres") || !strcmp(word, "pres")){
          printf("%s:\n",word);
          readMacro(fileIn,word,false);
/*... codigo da da funcao limitadora de fluxo*/        
          setDiffusionScheme(word,&sc.diffPres.iCod);
          nScheme--;
        }
 /*... transporte T1*/
        else if (!strcmp(word,"T1") || !strcmp(word, "t1")){
          printf("%s:\n", word);
          readMacro(fileIn, word, false);
 /*... codigo da da funcao limitadora de fluxo*/
          setDiffusionScheme(word, &sc.diffT1.iCod);
          nScheme--;
        }
 /*... Energy*/
        else if (!strcmp(word, "Energy") || !strcmp(word, "energy")) {
          printf("%s:\n", word);
          readMacro(fileIn, word, false);
 /*... codigo da da funcao limitadora de fluxo*/
          setDiffusionScheme(word, &sc.diffEnergy.iCod);
          nScheme--;
        }
      }while(nScheme);
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: pFluid : escreve os arquivos dos resultados do escoamente
 * imcompressivel                    
 *===================================================================*/
    else if( (!strcmp(word,macro[33])) && 
             (opt.stepPlotFluid[1]++) == opt.stepPlotFluid[0]){      
/*...*/
      opt.stepPlotFluid[1] = 1;
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
/*...................................................................*/

/*...*/
      ndfVel = max(mesh->ndfF - 1, mesh->ndfFt - 2);
/*...................................................................*/

/*... reconstruindo do gradiente (Pres)*/
      tm.rcGradPres = getTimeC() - tm.rcGradPres;
      rcGradU(&m                      ,loadsPres
             ,mesh->elm.node          ,mesh->elm.adj.nelcon
             ,mesh->elm.geom.cc       ,mesh->node.x   
             ,mesh->elm.nen           ,mesh->elm.adj.nViz 
             ,mesh->elm.geomType      ,mesh->elm.material.prop 
             ,mesh->elm.mat 
             ,mesh->elm.leastSquare   ,mesh->elm.leastSquareR
             ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
             ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea    
             ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
             ,mesh->elm.geom.vSkew     
             ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc    
             ,mesh->elm.geom.dcca
             ,mesh->elm.faceRpres     ,mesh->elm.faceLoadPres  
             ,mesh->elm.pressure      ,mesh->elm.gradPres                
             ,mesh->node.pressure     ,sc.rcGrad
             ,mesh->maxNo             ,mesh->maxViz
             ,1                       ,mesh->ndm       
             ,&pMesh->iNo             ,&pMesh->iEl  
             ,mesh->numelNov          ,mesh->numel        
             ,mesh->nnodeNov          ,mesh->nnode);  
      tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (GradPres)*/
      interCellNode(&m                   ,loadsPres
                   ,mesh->node.gradPres,mesh->elm.gradPres
                   ,mesh->elm.node     ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc  ,mesh->node.x  
                   ,mesh->elm.geom.xm
                   ,mesh->elm.nen      ,mesh->elm.adj.nViz
                   ,mesh->elm.faceRpres,mesh->elm.faceLoadPres
                   ,&pMesh->iNo          
                   ,mesh->numelNov     ,mesh->numel        
                   ,mesh->nnodeNov     ,mesh->nnode 
                   ,mesh->maxNo        ,mesh->maxViz   
                   ,mesh->ndm          ,1
                   ,mesh->ndm      
                   ,false              ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (pres)*/
      interCellNode(&m                 ,loadsPres
                   ,mesh->node.pressure,mesh->elm.pressure   
                   ,mesh->elm.node     ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc  ,mesh->node.x  
                   ,mesh->elm.geom.xm
                   ,mesh->elm.nen      ,mesh->elm.adj.nViz
                   ,mesh->elm.faceRpres,mesh->elm.faceLoadPres  
                   ,&pMesh->iNo          
                   ,mesh->numelNov     ,mesh->numel        
                   ,mesh->nnodeNov     ,mesh->nnode 
                   ,mesh->maxNo        ,mesh->maxViz   
                   ,1                  ,1
                   ,mesh->ndm      
                   ,opt.bconditions    ,2);
/*...................................................................*/

/*... calculo da matrix jacobiana das velocidades
                            | du1dx1 du1dx2 du1dx3 |   
                            | du2dx1 du2dx2 du2dx3 |   
                            | du3dx1 du3dx2 du3dx3 |
*/        
     tm.rcGradVel  = getTimeC() - tm.rcGradVel;
     rcGradU(&m                     ,loadsVel
           ,mesh->elm.node          ,mesh->elm.adj.nelcon
           ,mesh->elm.geom.cc       ,mesh->node.x   
           ,mesh->elm.nen           ,mesh->elm.adj.nViz 
           ,mesh->elm.geomType      ,mesh->elm.material.prop 
           ,mesh->elm.mat 
           ,mesh->elm.leastSquare   ,mesh->elm.leastSquareR
           ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
           ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea    
           ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
           ,mesh->elm.geom.vSkew      
           ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc    
           ,mesh->elm.geom.dcca
           ,mesh->elm.faceRvel      ,mesh->elm.faceLoadVel   
           ,mesh->elm.vel           ,mesh->elm.gradVel                           
           ,mesh->node.vel          ,sc.rcGrad
           ,mesh->maxNo             ,mesh->maxViz
           ,ndfVel                  ,mesh->ndm
           ,&pMesh->iNo             ,&pMesh->iEl  
           ,mesh->numelNov          ,mesh->numel        
           ,mesh->nnodeNov          ,mesh->nnode); 
     tm.rcGradVel = getTimeC() - tm.rcGradVel;  
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (gradVel)*/
      interCellNode(&m                 ,loadsVel
                   ,mesh->node.gradVel ,mesh->elm.gradVel    
                   ,mesh->elm.node     ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc  ,mesh->node.x  
                   ,mesh->elm.geom.xm
                   ,mesh->elm.nen      ,mesh->elm.adj.nViz
                   ,mesh->elm.faceRvel ,mesh->elm.faceLoadVel 
                   ,&pMesh->iNo          
                   ,mesh->numelNov     ,mesh->numel        
                   ,mesh->nnodeNov     ,mesh->nnode 
                   ,mesh->maxNo        ,mesh->maxViz   
                   ,ndfVel             ,mesh->ndm
                   ,mesh->ndm      
                   ,false              ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (vel)*/
      interCellNode(&m                 ,loadsVel
                   ,mesh->node.vel     ,mesh->elm.vel        
                   ,mesh->elm.node     ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc  ,mesh->node.x  
                   ,mesh->elm.geom.xm
                   ,mesh->elm.nen      ,mesh->elm.adj.nViz
                   ,mesh->elm.faceRvel ,mesh->elm.faceLoadVel 
                   ,&pMesh->iNo          
                   ,mesh->numelNov     ,mesh->numel        
                   ,mesh->nnodeNov     ,mesh->nnode 
                   ,mesh->maxNo        ,mesh->maxViz   
                   ,mesh->ndm          ,1
                   ,mesh->ndm      
                   ,opt.bconditions    ,2);
/*...................................................................*/

/*... reconstruindo do gradiente (Energia)*/
      if(mesh->ndfFt){
         rcGradU(&m                 ,loadsTemp
              ,mesh->elm.node       ,mesh->elm.adj.nelcon
              ,mesh->elm.geom.cc    ,mesh->node.x
              ,mesh->elm.nen        ,mesh->elm.adj.nViz
              ,mesh->elm.geomType   ,mesh->elm.material.prop
              ,mesh->elm.mat
              ,mesh->elm.leastSquare,mesh->elm.leastSquareR
              ,mesh->elm.geom.ksi   ,mesh->elm.geom.mksi
              ,mesh->elm.geom.eta   ,mesh->elm.geom.fArea
              ,mesh->elm.geom.normal,mesh->elm.geom.volume
              ,mesh->elm.geom.vSkew
              ,mesh->elm.geom.xm    ,mesh->elm.geom.xmcc
              ,mesh->elm.geom.dcca
              ,mesh->elm.faceRenergy,mesh->elm.faceLoadEnergy
              ,mesh->elm.temp       ,mesh->elm.gradTemp  
              ,mesh->node.temp      ,sc.rcGrad
              ,mesh->maxNo          ,mesh->maxViz
              ,1, mesh->ndm
              ,&pMesh->iNo          ,&pMesh->iEl
              ,mesh->numelNov       ,mesh->numel
              ,mesh->nnodeNov       ,mesh->nnode); 
/*.................................................................. */


/*... interpolacao das variaveis da celulas para pos nos (GradEnergy)*/
        interCellNode(&m            ,loadsTemp  
              ,mesh->node.gradTemp  ,mesh->elm.gradTemp  
              ,mesh->elm.node       ,mesh->elm.geomType
              ,mesh->elm.geom.cc    ,mesh->node.x
              ,mesh->elm.geom.xm
              ,mesh->elm.nen        ,mesh->elm.adj.nViz
              ,mesh->elm.faceRenergy,mesh->elm.faceLoadEnergy
              ,&pMesh->iNo
              ,mesh->numelNov       ,mesh->numel
              ,mesh->nnodeNov       ,mesh->nnode
              ,mesh->maxNo          ,mesh->maxViz
              ,mesh->ndm            ,1
              ,mesh->ndm
              ,false                ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (energy)*/
        interCellNode(&m           ,loadsTemp
             ,mesh->node.temp      ,mesh->elm.temp
             ,mesh->elm.node       ,mesh->elm.geomType
             ,mesh->elm.geom.cc    ,mesh->node.x
             ,mesh->elm.geom.xm
             ,mesh->elm.nen        ,mesh->elm.adj.nViz
             ,mesh->elm.faceRenergy,mesh->elm.faceLoadEnergy
             ,&pMesh->iNo
             ,mesh->numelNov       ,mesh->numel
             ,mesh->nnodeNov       ,mesh->nnode
             ,mesh->maxNo          ,mesh->maxViz
             ,1                    ,1
             ,mesh->ndm
             ,opt.bconditions      ,2);
      }
/*...................................................................*/

/*...*/
      if (turbModel.fTurb) {
        interCellNode(&m                     ,loadsVel
                   ,mesh->node.eddyViscosity ,mesh->elm.eddyViscosity        
                   ,mesh->elm.node           ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc        ,mesh->node.x  
                   ,mesh->elm.geom.xm       
                   ,mesh->elm.nen            ,mesh->elm.adj.nViz
                   ,mesh->elm.faceRvel       ,mesh->elm.faceLoadVel 
                   ,&pMesh->iNo                
                   ,mesh->numelNov           ,mesh->numel        
                   ,mesh->nnodeNov           ,mesh->nnode 
                   ,mesh->maxNo              ,mesh->maxViz   
                   ,1                        ,1
                   ,mesh->ndm               
                   ,false                    ,2);
      }                                 
/*...................................................................*/

/*...*/
      if(!mpiVar.myId ){
        fName(preName,sc.ddt.timeStep,0,21,&nameOut);

/*...*/
        wResVtkFluid(&m                  , mesh0->node.x      
               , mesh0->elm.node         , mesh0->elm.mat    
               , mesh0->elm.nen          , mesh0->elm.geomType
               , mesh0->elm.pressure     , mesh0->node.pressure
               , mesh0->elm.gradPres     , mesh0->node.gradPres  
               , mesh0->elm.vel          , mesh0->node.vel      
               , mesh0->elm.gradVel      , mesh0->node.gradVel 
               , mesh0->elm.temp         , mesh0->node.temp   
               , mesh0->elm.gradTemp     , mesh0->node.gradTemp
               , mesh0->elm.eddyViscosity, mesh0->node.eddyViscosity
               , mesh0->elm.densityFluid , mesh0->elm.specificHeat
               , mesh0->elm.dViscosity   , mesh0->elm.tConductivity
               , mesh0->elm.wallParameters                  
               , mesh0->nnode            , mesh0->numel  
               , mesh0->ndm              , mesh0->maxNo 
               , mesh0->numat            , ndfVel
               , nameOut                 , opt
               , eModel.fKelvin
               , sc.ddt                  , fileOut);  
/*...................................................................*/
      }
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: setPrintFluid                                    
 *===================================================================*/
    else if((!strcmp(word,macro[34]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
/*...*/
      setPrintFluid(&opt,fileIn);
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
* macro: setPrime
*===================================================================*/
    else if ((!strcmp(word, macro[36]))) {
      if (!mpiVar.myId) {
        printf("%s\n", DIF);
        printf("%s\n", word);
        printf("%s\n", DIF);
      }
/*...*/
      prime = (Prime*)malloc(sizeof(Prime));
      if (prime == NULL) {
        printf("Erro ponteiro Prime\n");
        exit(EXIT_FAILURE);
      }
      fSolvPrime       = true;
      prime->maxIt     = 2000;
      prime->alphaPres = 1.0e0;
      prime->alphaVel  = 0.9e0;
      prime->kZeroVel  = 1;
      prime->kZeroPres = 0;
      prime->sPressure = true;
      prime->nNonOrth  = 0;
      prime->tolPres   = 1.e-04;
      prime->tolVel    = 1.e-04;
      prime->pPrime    = 10;
/*...................................................................*/

/*...*/
      readMacro(fileIn, word, false);
      if (!strcmp(word, "config:")) 
        setPrimeScheme(word, prime, fileIn);
/*...................................................................*/

/*...*/
      HccaAlloc(DOUBLE, &m, prime->d
               ,mesh->numel*mesh->ndm, "dField", false);
      zero(prime->d, mesh->numel*mesh->ndm, DOUBLEC);

      HccaAlloc(DOUBLE, &m, prime->velUp
               ,mesh->numel*mesh->ndm, "velUp", false);
      zero(prime->d, mesh->numel*mesh->ndm, DOUBLEC);

      HccaAlloc(DOUBLE, &m, prime->ePresC
               ,mesh->numel, "ePresC", false);
      zero(prime->ePresC, mesh->numel, DOUBLEC);

      HccaAlloc(DOUBLE, &m, prime->nPresC
               ,mesh->nnode, "nPresC", false);
      zero(prime->nPresC, mesh->numel, DOUBLEC);

      HccaAlloc(DOUBLE, &m, prime->eGradPresC
                , mesh->numel*mesh->ndm, "eGradPresC", false);
      zero(prime->eGradPresC, mesh->numel*mesh->ndm, DOUBLEC);

      HccaAlloc(DOUBLE, &m, prime->ePresC1
                , mesh->numel, "ePresC1", false);
      zero(prime->ePresC, mesh->numel, DOUBLEC);

      HccaAlloc(DOUBLE, &m, prime->aD
               ,mesh->numel*mesh->ndm, "aD", false);
      zero(prime->aD, mesh->numel*mesh->ndm, DOUBLEC);

      HccaAlloc(DOUBLE, &m, prime->bTemporal
               ,mesh->numel*mesh->ndm, "bT", false);
      zero(prime->bTemporal, mesh->numel*mesh->ndm, DOUBLEC);
/*...................................................................*/

      if (!mpiVar.myId) printf("%s\n", DIF);
    }
/*===================================================================*/

/*===================================================================*
* macro: prime
*===================================================================*/
    else if ((!strcmp(word, macro[37]))) {
      if (!mpiVar.myId) {
        printf("%s\n", DIF);
        printf("%s\n", word);
        printf("%s\n", DIF);
      }
      mpiWait();
      tm.solvEdpFluid = getTimeC() - tm.solvEdpFluid;
/*...*/
      if (solvPres == NULL) {
        printf("Estrutara de dados nao montada para o solvFluid!!!\n");
        exit(EXIT_FAILURE);
      }
/*...................................................................*/

/*...*/
     if(!fSolvPrime) {
        printf("Prime nao configurado ainda!!!\n");
        exit(EXIT_FAILURE);
     }  
/*...................................................................*/

/*...*/
      primeSolver(&m
                 ,loadsVel  ,loadsPres
                 ,mesh0     ,mesh
                 ,sistEqPres,solvPres
                 ,prime     ,sc 
                 ,pMesh     ,opt
                 ,preName   ,nameOut
                 ,fileOut   );
/*...................................................................*/

/*...*/
      tm.solvEdpFluid = getTimeC() - tm.solvEdpFluid;
/*...................................................................*/
      if (!mpiVar.myId) printf("%s\n\n", DIF);
    }
/*===================================================================*/

/*===================================================================*
 * macro: vprop propriedades variaveis                             
 *===================================================================*/
    else if((!strcmp(word,macro[38]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
      readPropVar(&propVarFluid,fileIn);
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: gravity                                                  
 *===================================================================*/
    else if((!strcmp(word,macro[39]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }   
      readGravity(gravity,fileIn);
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: model modelos utilizados nas equacao diferencias         
 *===================================================================*/
    else if((!strcmp(word,macro[40]))){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
      readModel(&eModel,&turbModel
              , &eMass ,&eMomentum
              , fileIn);
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*/
  }while(macroFlag && (!feof(fileIn)));

  mpiStop();
  return EXIT_SUCCESS;
}    
/*********************************************************************/
 
 
