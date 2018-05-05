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
#include<Print.h>
#include<OpenMp.h>
#include<SaveLoad.h>
#include<FaceStruct.h>
/*********************************************************************/

/*********************************************************************/
#ifdef _DEBUG_ 
  #include<Debug.h>
#endif
/*********************************************************************/

/*********************************************************************/
void testeFace(Geom *geom, Face *face
             , INT *cellFace, short *nFace
             , INT numel, short maxViz
             , short ndm, short maxNo);

/*********************************************************************/

int main(int argc,char**argv){

/*...Memoria principal*/  
  Memoria m;
/*... estrutura da dados para a malha*/
  Mesh *mesh=NULL,*mesh0=NULL;
/*... Sistema de equacao*/
  SistEq sistEqD1, *sistEqT1=NULL;
  SistEq sistEqVel, sistEqPres, sistEqEnergy, sistEqKturb;
/*... metodo de acoplamento pressao-velocidade*/
  Simple simple;
  Prime  prime;
/*... tubulence*/
  Turbulence turbModel;
/*...*/
  Mean media;
/*...*/
  EnergyModel eModel;
  MassEqModel eMass;
  MomentumModel ModelMomentum;
/*... propriedade variaveis*/
  PropVar propVarFluid;

/*... solver*/
  INT nEqMax;
  Solv solvD1, *solvT1 = NULL;
  Solv solvVel,solvPres,solvEnergy,solvKturb;
  bool fSolvD1 = false, fSolvT1 = false;
  bool fSolvVel = false,fSolvPres = false, fSolvEnergy = false;
  bool fSolvKturb = false;
  bool fSolvSimple = false,fSolvPrime = false;
/*... reordenacao da malha*/
  Reord  *reordMesh=NULL;

/*... particionamento*/
  PartMesh *pMesh = NULL;
/*...*/
  Save save;

/*...*/
  char loopWord[100][MAX_LINE];
  unsigned short kLoop = 0 ,jLoop = 0;
  bool flWord=false;

/*... Estrutura de dados*/
  char strIa[MNOMEPONTEIRO],strJa[MNOMEPONTEIRO];
  char strA[MNOMEPONTEIRO],strAd[MNOMEPONTEIRO];

/*... arquivo*/
  char nameIn[MAX_STR_LEN_IN], nameOut[SIZEMAX];
  char auxName[MAX_STR_LEN_SUFIXO], preName[MAX_STR_LEN_SUFIXO];
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
  {"help"        ,"mesh"         ,"stop"         /* 0, 1, 2*/
  ,"config"      ,"nextLoop"     ,"rcGrad"       /* 3, 4, 5*/
  ,"pgeo"        ,"pcoob"        ,"pcoo"         /* 6, 7, 8*/ 
  ,"setSolvDiff" ,"presolvT1"    ,"openmp"       /* 9,10,11*/
  ,"solvD1"      ,""             ,"pD1"          /*12,13,14*/
  ,"nlIt"        ,"pD1CsvCell"   ,"pD1CsvNode"   /*15,16,17*/
  ,"solvT1"      ,""             ,"pT1"          /*18,19,20*/
  ,""            ,"pT1CsvCell"   ,"pT1CsvNode"   /*21,22,23*/
  ,"setSolv"     ,"simple"       ,"setSimple"    /*24,25,26*/
  ,"transient"   ,"timeUpdate"   ,"partd"        /*27,28,29*/
  ,"advection"   ,"edp"          ,"diffusion"    /*30,31,32*/
  ,"pFluid"      ,"setPrint"     ,""             /*33,34,35*/
  ,"setPrime"    ,"prime"        ,"propVar"      /*36,37,38*/
  ,"gravity"     ,"model"        ,"mean"         /*39,40,41*/
  ,"setMean"     ,"save"         ,"load"};       /*42,43,44*/
/* ..................................................................*/

/*... Memoria principal(valor padrao - bytes)*/
  nmax = 200000;
/* ..................................................................*/

/*...*/
  media.fMedia = false;
  media.fVel   = false;
/*..................................................................*/

/*...*/
  save.fLoad = false;
  save.fSave = true;
  save.step[0] =  1;
  save.step[1] = save.step[0];
/*..................................................................*/

/*...*/  
  thDynamic.pTh[0]      = PREREF;
  thDynamic.pTh[1]      = PREREF;
  thDynamic.pTh[2]      = PREREF;
/*....................................................................*/

/*...*/
  turbModel.fWall               = false;
  turbModel.fTurb               = false;
  turbModel.fDynamic            = false;
  turbModel.wallType            = STANDARDWALL;
  turbModel.type                = LES;
  turbModel.typeLes             = LESFUNCMODEL;
  turbModel.typeMixed[FUNMODEL] = SMAGORINSKY;
  turbModel.typeMixed[ESTMODEL] = CLARK;
  turbModel.cs                  = 1.0e0;
  turbModel.cf                  = 0.2e0;
  turbModel.c                   = 1.0e0;
  turbModel.PrandltTwall        = 0.85e0;
  turbModel.PrandltTsgs         = 0.85e0;
/*....................................................................*/

/*...*/
  eModel.fPresWork     = false;
  eModel.fDissipation  = false;
  eModel.fRes          = true;
  eModel.fTemperature  = false;
  eModel.fKelvin       = false;
/*...................................................................*/

/*...*/
  ModelMomentum.fRes             = true;
  ModelMomentum.fRhieChowInt     = false;
  ModelMomentum.iCodBuoyant      = BUOYANT_RHOREF;
  ModelMomentum.fViscosity       = true;
  ModelMomentum.fDiv             = true;
/*...................................................................*/

/*...*/
  eMass.RhsDensity = false;
  eMass.LhsDensity = false;
/*...................................................................*/

/*... OpenMP*/
  ompVar.flag           = false;
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
  opt.presTotal     = true;
  opt.uD1           = true;
  opt.uT1           = true;
  opt.energy        = true;
  opt.gradVel       = false;
  opt.gradPres      = false;
  opt.gradEnergy    = false;
  opt.graduD1       = true;
  opt.graduT1       = true;
  opt.eddyViscosity = false;
  opt.densityFluid  = false;
  opt.specificHeat  = false;
  opt.dViscosity    = false;
  opt.tConductivity = false;
  opt.vorticity     = false;
  opt.wallParameters= false;
  opt.kinetic       = false;
  opt.stressR       = false;
  opt.cDynamic      = false;
  opt.bconditions   = true;
  opt.cc            = false;
  opt.pKelvin       = false;
 
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

/*...*/
  initTime(&tm);
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
  sc.nlD1.pPlot   = sc.nlT1.pPlot = 1;
/*...................................................................*/

/*...*/
  sc.advT1.iCod1 = TVD;
  sc.advT1.iCod2 = VANLEERFACE;
  sc.diffT1.iCod = OVERRELAXED;

/*...*/
  sc.diffD1.iCod = OVERRELAXED;
/*...................................................................*/

/*...*/
  sc.advVel.iCod1    = TVD;
  sc.advVel.iCod2    = VANLEERFACE;
  sc.advVel.par[0]   = 0.e0;
  sc.advEnergy.iCod1 = TVD;
  sc.advEnergy.iCod2 = VANLEERFACE;
  sc.advVel.par[0]   = 0.e0;
  sc.advKturb.iCod1  = TVD;
  sc.advKturb.iCod2  = MIDMODFACE;
  sc.advVel.par[0]   = 0.e0;
/*...................................................................*/

/*...*/  
  sc.diffVel.iCod    = OVERRELAXED;
  sc.diffPres.iCod   = OVERRELAXED;
//  sc.diffPres.iCod   = ORTHOGONAL;
  sc.diffEnergy.iCod = OVERRELAXED;
  sc.diffKturb.iCod  = OVERRELAXED;
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
  reordMesh = (Reord*) malloc(sizeof(Reord));
  ERRO_MALLOC(reordMesh,"reordMesh",__LINE__,__FILE__,__func__);
  reordMesh->flag = false; 
/* ..................................................................*/
    
/*... abrindo ar quivo de entrada*/ 

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
  if( argc > 2)
    strcpy(preName,argv[2]);
  else{
    if(!mpiVar.myId ) printf("Prefixo do arquivo de saida:\n");
    if(!mpiVar.myId ) scanf("%s",preName);
  }
/*...................................................................*/
 
/*...*/
  fName(preName,0,0, 23 ,nameOut);
  fileLogExc = openFileBuffer(nameOut,"w",false);
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
      if(!mpiVar.myId){ 
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word); 
        fprintf(fileLogExc,"%s\n\n",DIF);
      }
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
                      ,&turbModel  ,&media       
                      ,fileIn);
      mpiWait();
/*...................................................................*/
 
/*...*/
      if(!mpiVar.myId){
        tm.adjcency = getTimeC();
/*... calcula a vizinhaca do elementos*/
        neighbors(&m                   ,mesh0->elm.node 
                 ,mesh0->elm.adj.nelcon,mesh0->elm.adj.nViz
                 ,mesh0->elm.nen       ,&mesh0->nFaces
                 ,mesh0->nnode         ,mesh0->numel     
                 ,mesh0->maxNo         ,mesh0->maxViz   
                 ,mesh0->ndm);        
/*...................................................................*/
        
/*...*/        
        faceStruct(&m, mesh0);
/*...................................................................*/
        tm.adjcency = getTimeC() - tm.adjcency;
      }
/*...................................................................*/

/*... identifica parede impermevais*/
      if(mesh0->ndfF > 0 || mesh0->ndfFt > 0){
        wallFluid(mesh0->elm.faceRvel,mesh0->elm.adj.nelcon
                 ,mesh0->elm.adj.nViz   
                 ,mesh0->numel       ,mesh0->maxViz); 
        if(turbModel.fTurb)
          wallFluid(mesh0->elm.faceReKturb,mesh0->elm.adj.nelcon
                   ,mesh0->elm.adj.nViz   
                   ,mesh0->numel          ,mesh0->maxViz);  
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
            fName(preName,mpiVar.nPrcs,0,1,nameOut);
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
          fName(preName,mpiVar.nPrcs,mpiVar.myId,2,nameOut);
          wMeshPartVtk(&m            
                      ,mesh->node.x      ,mesh->elm.node              
                      ,mesh->elm.nen     ,mesh->elm.geomType
                      ,mesh->nnode       ,mesh->numel    
                      ,mesh->ndm               
                      ,mesh->maxNo       ,mesh->maxViz
                      ,nameOut           ,opt.bVtk             
                      ,fileOut);
          fName(preName,mpiVar.nPrcs,mpiVar.myId,18,nameOut);
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
               ,mesh->elm.adj.nelcon ,mesh->elm.adj.nViz
               ,mesh->elm.geomType   ,mesh->elm.nen        
               ,mesh->elm.cellFace
               ,&mesh->elm.geom      ,&mesh->face   
               ,mesh->maxNo          ,mesh->maxViz
               ,mesh->ndm            ,mesh->numelNov);
      tm.geom = getTimeC() - tm.geom;
/*    testeFace(&mesh->elm.geom    , &mesh->face
               , mesh->elm.cellFace, mesh->elm.adj.nViz
               , mesh->numel       , mesh->maxViz
               , mesh->ndm         , mesh->maxNo); */
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
        if(ModelMomentum.iCodBuoyant == BUOYANT_HYDROSTATIC)
          hPres(mesh->elm.pressure0   , mesh->elm.pressure
              , mesh->elm.densityFluid, mesh->elm.geom.cc
              , gravity               , mesh->xRef
              , mesh->numel           , mesh->ndm );      
/*...................................................................*/
      }
/*...................................................................*/  
 
/*... reconstrucao de gradiente least square*/
      if(sc.rcGrad ==  RCLSQUARE || sc.rcGrad ==  RCLSQUAREQR){
        if(!mpiVar.myId ){
          fprintf(fileLogExc,"%s\n",DIF);
          fprintf(fileLogExc,"Least Square ...\n");
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
          fprintf(fileLogExc,"Least Square.\n");
          fprintf(fileLogExc,"%s\n",DIF);
        }
      }
/*...................................................................*/

/*... qualidade da malha*/
      meshQuality(&mesh->mQuality
                 ,mesh->elm.adj.nViz     ,mesh->elm.geom.volume
                 ,mesh->elm.geom.ksi     ,mesh->elm.geom.normal
                 ,mesh->elm.geom.mvSkew  ,mesh->elm.geom.dcca
                 ,mesh->maxViz           ,mesh->ndm
                 ,mesh->numelNov);
/*... qualidade da malha global*/
      if(mpiVar.nPrcs > 1){
        globalMeshQuality(&mesh->mQuality,&mesh0->mQuality);
      }  
/*...................................................................*/
         
/*... reodenando as celulas para dimuincao da banda*/
      HccaAlloc(INT,&m,reordMesh->num,mesh->numel,"rNum" ,_AD_);
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n",DIF);
      if(!mpiVar.myId )fprintf(fileLogExc,"Reordenando a malha ...\n");
      tm.reord = getTimeC() - tm.reord;
      reord(&m                ,reordMesh->num,mesh->elm.adj.nelcon
           ,mesh->elm.adj.nViz,mesh->maxViz  
           ,mesh->numel       ,mesh->numelNov
           ,reordMesh->flag   ,mpiVar.nPrcs);
      tm.reord = getTimeC() - tm.reord;
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"Malha reordenada.\n");
        fprintf(fileLogExc,"%s\n\n",DIF);
      }
/*...................................................................*/

    }   
/*===================================================================*/

/*===================================================================*
 * macro: stop : finalizacao do programa
 *===================================================================*/
    else if((!strcmp(word,macro[2]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word); 
        fprintf(fileLogExc,"%s\n\n",DIF);
      }
      tm.total = getTimeC() - tm.total;
/*... */
      fName(preName,mpiVar.nPrcs,mpiVar.myId,7,nameOut);
      fileLog = openFile(nameOut,"w");
      writeLog(*mesh      ,sc
              ,&solvD1    ,&sistEqD1
              ,solvT1     ,sistEqT1
              ,&solvVel   ,&sistEqVel
              ,&solvPres  ,&sistEqPres
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
          fName(preName,mpiVar.nPrcs,mpiVar.myId,60,nameOut);
          fileLog = openFile(nameOut,"w");
        }
        writeLogMeanTime(*mesh0    ,sc
                        ,&solvD1   ,&sistEqD1
                        ,solvT1    ,sistEqT1
                        ,tm       
                        ,fSolvD1   ,fSolvT1     
                        ,nameIn    ,fileLog);
        
        if(!mpiVar.myId) fclose(fileLog);
      } 
/*...................................................................*/

/*... fechando o arquivo do log do solver linear Pres*/
      if(fSolvPres && solvPres.log && !mpiVar.myId)  
        fclose(solvPres.fileSolv);
/*... fechando o arquivo do log do solver linear Vel*/
      if(fSolvVel && solvVel.log && !mpiVar.myId)  
        fclose(solvVel.fileSolv);
/*... fechando o arquivo do log do solver linear T1*/
      if(fSolvD1 && solvD1.log && !mpiVar.myId)  
        fclose(solvD1.fileSolv);
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
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
      
      config(&opt          ,reordMesh ,fileIn);
      
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: nextLoop : proximo loop tempo            
 *===================================================================*/
    else if((!strcmp(word,macro[4]))){
      jLoop = 0;
    }   
/*===================================================================*/

/*===================================================================*
* macro: rcGrad: tipo de reconstrucao de gradeiente
*===================================================================*/
    else if ((!strcmp(word, macro[5]))) {
      if (!mpiVar.myId) {
        fprintf(fileLogExc, "%s\n", DIF);
        fprintf(fileLogExc, "%s\n", word);
      }
      setReGrad(&sc.rcGrad, fileIn);
      if (!mpiVar.myId) fprintf(fileLogExc, "%s\n\n", DIF);
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
        fName(preName,0,0,6,nameOut);
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
/*...................................................................*/

/*... geometrica completa*/
        fName(preName, 0, 0, 25, nameOut);
        wGeoVtk2(&m                   , mesh0->node.x
               , mesh0->elm.geom.cc   , mesh0->elm.node
               , mesh0->elm.nen       , mesh0->elm.geomType
               , mesh0->nnode         , mesh0->numel
               , mesh0->ndm           , mesh0->maxNo 
               , mesh0->maxViz
               , nameOut              , opt.bVtk
               , fileOut);
/*...................................................................*/

/*... face com cargas*/
       fName(preName,0,0,17,nameOut);
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
/*...................................................................*/
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
       
        fName(preName,0,0,13,nameOut);
/*...*/
        writeCoo(&m,sistEqD1.ia,sistEqD1.ja,sistEqD1.neq
                ,sistEqD1.au   ,sistEqD1.ad,sistEqD1.al        
                ,sistEqD1.nad  ,sistEqD1.storage
                ,sistEqD1.unsym,true
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
      
        fName(preName,0,0,12,nameOut);
/*... matriz A*/
        writeCoo(&m,sistEqD1.ia,sistEqD1.ja,sistEqD1.neq
                ,sistEqD1.au   ,sistEqD1.ad,sistEqD1.al        
                ,sistEqD1.nad  ,sistEqD1.storage
                ,sistEqD1.unsym,false 
                ,nameOut);
/*...................................................................*/

/*... vetor de forcas b*/      
        fName(preName,0,0,14,nameOut);
        for(int i=0;i<sistEqD1.neq;i++)
          sistEqD1.b[i] /= sistEqD1.ad[i];

        writeCooB(sistEqD1.b,sistEqD1.neq,nameOut);
/*...................................................................*/
        printf("%s\n\n",DIF);
      }
    }   
/*===================================================================*/


/*===================================================================*
 * macro: setSolvDiff : problema de difusao pura
 *===================================================================*/
    else if((!strcmp(word,macro[9])))
    {
      if(!mpiVar.myId )
      {
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
      readSolvDiff(&m     , mesh     , reordMesh
                  ,&solvD1, &sistEqD1, &fSolvD1
                  ,auxName, preName  , nameOut
                  ,fileIn , &opt);

/*...................................................................*/
      if(!mpiVar.myId  ) fprintf(fileLogExc,"%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: presolvt1 : problema de transporte  
 *===================================================================*/
    else if((!strcmp(word,macro[10]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
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
        fName(auxName,mpiVar.nPrcs,0,11,nameOut);
        solvT1->fileSolv = openFile(nameOut,"w");
      }
/*...................................................................*/

/*...*/
      if(opt.fItPlot && !mpiVar.myId){  
        strcpy(auxName,preName);
        strcat(auxName,"_T1");
        fName(auxName,mpiVar.nPrcs,0,10,nameOut);
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

        if(!mpiVar.myId ) 
          fprintf(fileLogExc,"MaxIt     : %d\n",solvT1->maxIt);
        if(!mpiVar.myId ) 
          fprintf(fileLogExc,"Tol       : %e\n",solvT1->tol);
      }
/*...................................................................*/

/*... numeracao das equacoes*/
      HccaAlloc(INT,&m,sistEqT1->id
               ,mesh->numel*mesh->ndfT[0]
               ,"sistT1id",_AD_);
      if(!mpiVar.myId){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"Numerando as equacoes.\n");
      }
      tm.numeqT1 = getTimeC() - tm.numeqT1;
      sistEqT1->neq = numeq(sistEqT1->id       ,reordMesh->num
                           ,mesh->elm.faceRt1  ,mesh->elm.adj.nViz
                           ,mesh->numel        ,mesh->maxViz
                           ,mesh->ndfT[0]);
      tm.numeqT1 = getTimeC() - tm.numeqT1;
      if(!mpiVar.myId){
        fprintf(fileLogExc,"Equacoes numeradas.\n");
        fprintf(fileLogExc,"%s\n",DIF);
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
      if(!mpiVar.myId) 
        fprintf(fileLogExc,"Montagem da estrura de dados esparsa.\n");
      tm.dataStructT1 = getTimeC() - tm.dataStructT1 ;
      dataStruct(&m,sistEqT1->id   ,reordMesh->num,mesh->elm.adj.nelcon
                ,mesh->elm.adj.nViz,mesh->numelNov,mesh->maxViz
                ,mesh->ndfT[0]     ,strIa         ,strJa
                ,strAd             ,strA          ,sistEqT1);
      tm.dataStructT1 = getTimeC() - tm.dataStructT1 ;
      if(!mpiVar.myId) fprintf(fileLogExc,"Estrutuda montada.\n");
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
      if(!mpiVar.myId  ) fprintf(fileLogExc,"%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
* macro: openmp: configuracao do openmp  
*===================================================================*/
    else if ((!strcmp(word, macro[11]))) {
      if(!mpiVar.myId  ) fprintf(fileLogExc,"%s\n\n",DIF);
/*... */
      openMpSet(fileIn,&ompVar);
/*..................................................................*/
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*===================================================================*/

/*===================================================================*
 * macro: solvd1: problema de difusao pura
 *===================================================================*/
    else if((!strcmp(word,macro[12])))
    {
      if(!mpiVar.myId )
      {
        printf("%s\n",DIF);
        printf("%s\n",word);
        printf("%s\n",DIF);
      }
      mpiWait();
      tm.solvEdpD1    = getTimeC() - tm.solvEdpD1;
/*...*/
      if(!fSolvD1)
      {
        fprintf(fileLogExc
                ,"Estrutara de dados nao montada para o solvD1!!!\n");
        exit(EXIT_FAILURE);
      }
/*...................................................................*/
     
/*...*/
      diffusion(&m         ,loadsD1        
               ,mesh0      ,mesh           ,&sistEqD1
               ,&solvD1    ,&sc            ,pMesh
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
      if (!mpiVar.myId) {
        printf("%s\n", DIF);
        printf("%s\n", word);
      }
      printDiff(&m
              , pMesh   , &sc
              , loadsD1 , &opt
              , mesh0   , mesh
              , preName , nameOut);
      if (!mpiVar.myId) printf("%s\n\n", DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: nlIt: configura das iteracoes nao lineares             
 *===================================================================*/
    else if((!strcmp(word,macro[15]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }

      readNlIt(&sc, fileIn);
/*...................................................................*/
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
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
        fName(auxName,sc.ddt.timeStep,0,16,nameOut);
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
        fName(auxName,sc.ddt.timeStep,0,16,nameOut);
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
        fName(preName,sc.ddt.timeStep,0,20,nameOut);

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
 * macro:             
 *===================================================================*/
    else if((!strcmp(word,macro[21]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
/*...................................................................*/
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
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
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      
/*...*/
        strcpy(auxName,preName);
        strcat(auxName,"_T1_cell_");
        fName(auxName,sc.ddt.timeStep,0,16,nameOut);
        fileOut = openFile(nameOut,"w");
/*...*/
        writeCsvCell(mesh0->elm.uT1    ,mesh0->elm.gradUt1
                    ,mesh0->elm.geom.cc                  
                    ,mesh0->numel      ,mesh0->ndfT[0]
                    ,mesh0->ndm        ,fileOut);
/*...*/
        fclose(fileOut);
/*...................................................................*/
        fprintf(fileLogExc,"%s\n\n",DIF);
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
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
/*...*/
        strcpy(auxName,preName);
        strcat(auxName,"_T1_node_");
        fName(auxName,sc.ddt.timeStep,0,16,nameOut);
        fileOut = openFile(nameOut,"w");
/*...*/
        writeCsvNode(mesh0->node.uT1    ,mesh0->node.gradUt1
                  ,mesh0->node.x                  
                  ,mesh0->nnode         ,mesh0->ndfT[0]
                  ,mesh0->ndm           ,fileOut);
/*...*/
        fclose(fileOut);
/*...................................................................*/
        fprintf(fileLogExc,"%s\n\n",DIF);
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: setSolv : escoamento de fluidos  
 *===================================================================*/
    else if((!strcmp(word,macro[24]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
/*...................................................................*/

/*...*/
      readSolvFluid(&m          , mesh          , reordMesh
                   , &solvVel   , &sistEqVel    , &fSolvVel
                   , &solvPres  , &sistEqPres   , &fSolvPres
                   , &solvEnergy, &sistEqEnergy , &fSolvEnergy
                   , &solvKturb , &sistEqKturb  , &fSolvKturb
                   , auxName    , preName       , nameOut
                   , fileIn                     , &opt);
/*...................................................................*/
      if(!mpiVar.myId  ) fprintf(fileLogExc,"%s\n\n",DIF);
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
      if(!fSolvVel || !fSolvPres){
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
                      ,&sistEqVel ,&sistEqPres
                      ,&solvVel   ,&solvPres
                      ,&simple
                      ,sc         ,pMesh
                      ,opt        ,preName        
                      ,nameOut    ,fileOut);
/*...................................................................*/

/*...*/
      else if(mesh->ndfFt)
        simpleSolverLm(&m           , propVarFluid
                     , loadsVel     , loadsPres 
                     , loadsEnergy  , loadsKturb
                     , eModel       
                     , eMass        , ModelMomentum
                     , &turbModel   , &thDynamic   
                     , mesh0        , mesh
                     , &sistEqVel   , &sistEqPres
                     , &sistEqEnergy, &sistEqKturb   
                     , &solvVel     , &solvPres
                     , &solvEnergy  , &solvKturb  
                     , &simple      , &sc     
                     , pMesh        , &media 
                     , opt          , preName
                     , nameOut      , fileOut);  
/*...................................................................*/

/*...*/
     tm.solvEdpFluid    = getTimeC() - tm.solvEdpFluid;
/*...................................................................*/
     if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }
/*===================================================================*/

/*===================================================================*
 * macro: setSimple: configuracoe do metodo simple
 *===================================================================*/
    else if((!strcmp(word,macro[26]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
        fprintf(fileLogExc,"%s\n",DIF);
      }
/*...*/
      readSetSimple(&m     , fileIn
                  , mesh0  , mesh
                  , &simple, &fSolvSimple);
/*...................................................................*/

      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*===================================================================*/

/*===================================================================*
 * macro: transient:configuracao da discretizacao temporal                 
 *===================================================================*/
    else if((!strcmp(word,macro[27]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
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
          if(!mpiVar.myId) fprintf(fileLogExc,"dynamic : True\n");
        }
        else {
          if(!mpiVar.myId) fprintf(fileLogExc,"dynamic : False\n");
        }          
/*...*/        
        if(!mpiVar.myId ) 
          fprintf(fileLogExc,"dt(s)     : %lf\n",sc.ddt.dt[0]);
        if(!mpiVar.myId ) 
          fprintf(fileLogExc,"Total(s)  : %lf\n",sc.ddt.total);
      
        if(sc.ddt.type == EULER && !mpiVar.myId)     
          fprintf(fileLogExc,"ddtScheme : EULER\n");
        else if(sc.ddt.type == BACKWARD && !mpiVar.myId )     
          fprintf(fileLogExc,"ddtScheme : BACKWARD\n");

        if(!save.fLoad){
          sc.ddt.t         = 0.e0;
          sc.ddt.dtInicial = sc.ddt.dt[0];
          sc.ddt.dt[1]     = sc.ddt.dt[0];
          sc.ddt.dt[2]     = sc.ddt.dt[0];
          sc.ddt.timeStep  = 0;
        }
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
      strcpy(loopWord[kLoop-1],"nextLoop");
/*...................................................................*/
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
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
//    jLoop            = 0;
      sc.ddt.t        += sc.ddt.dt[0];
      sc.ddt.timeStep ++; 
/*...................................................................*/

/*...*/
      if(sc.ddt.t > sc.ddt.total + 0.1e0*sc.ddt.dt[TIME_N])
        flWord = false;  
/*    if(sc.ddt.t > sc.ddt.total)
        flWord = false;  */
/*...................................................................*/
      
/*...*/
      else{
        if(!mpiVar.myId ){
          printf("dt(n-2) = %lf\n",sc.ddt.dt[TIME_N_MINUS_2]);
          printf("dt(n-1) = %lf\n",sc.ddt.dt[TIME_N_MINUS_1]);
          printf("dt(n)   = %lf\n",sc.ddt.dt[TIME_N ]);
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
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
/*...*/
      pMesh->fPartMesh = true;
      readMacro(fileIn,word,false);
      if(!strcmp(word,"config:")){
/*... fPrintMesh*/        
        readMacro(fileIn,word,false);
        if(!strcmp(word,"true")){ 
          pMesh->fPrintMesh = true;
          if(!mpiVar.myId ) fprintf(fileLogExc,"fPrintMesh    : true\n");
        }

/*... fPrintMeshPart*/        
        readMacro(fileIn,word,false);
        if(!strcmp(word,"true")){ 
          pMesh->fPrintMeshPart = true;
          if(!mpiVar.myId ) fprintf(fileLogExc,"fPrintMeshPart: true\n");
        }
      }
/*...................................................................*/
       if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: advection tecnica aplicada no termo advectivo          
 *===================================================================*/
    else if((!strcmp(word,macro[30]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
/*...*/
      readAdvectionScheme(fileIn, &sc);
 /*...................................................................*/
      if (!mpiVar.myId) fprintf(fileLogExc,"%s\n\n", DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: edo equacoes diferencias resolvidas                    
 *===================================================================*/
    else if((!strcmp(word,macro[31]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
      readEdo(mesh0,fileIn);
/*...................................................................*/
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: diffusion tecnica aplicada no termo difusivo           
 *===================================================================*/
    else if((!strcmp(word,macro[32]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
/*... */
      readDiffusionScheme(fileIn, &sc);
/*...................................................................*/
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
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
      printFluid(&m         
               , &turbModel, &eModel
               , pMesh     , sc
               , loadsVel  , loadsPres 
               , loadsTemp , opt
               , mesh0     , mesh 
               , &media      
               , preName   , nameOut);
/*...................................................................*/

      if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: setPrint                                   
 *===================================================================*/
    else if((!strcmp(word,macro[34]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
/*...*/
      setPrint(&opt,fileIn);
/*...................................................................*/
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
* macro: setPrime
*===================================================================*/
    else if ((!strcmp(word, macro[36]))) {
      if (!mpiVar.myId) {
        fprintf(fileLogExc,"%s\n", DIF);
        fprintf(fileLogExc,"%s\n", word);
        fprintf(fileLogExc,"%s\n", DIF);
      }
/*...*/
      readSetPrime(&m    , fileIn
                  , mesh0, mesh
                  , &prime, &fSolvPrime);

      if (!mpiVar.myId) printf("%s\n\n", DIF);
    }
/*===================================================================*/

/*===================================================================*
* macro: prime
*===================================================================*/
    else if ((!strcmp(word, macro[37]))) {
      if (!mpiVar.myId) {
        fprintf(fileLogExc,"%s\n", DIF);
        fprintf(fileLogExc,"%s\n", word);
        fprintf(fileLogExc,"%s\n", DIF);
      }
      mpiWait();
      tm.solvEdpFluid = getTimeC() - tm.solvEdpFluid;
/*...*/
      if (!fSolvPres) {
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
                 ,loadsVel   ,loadsPres
                 ,mesh0      ,mesh
                 ,&sistEqPres,&solvPres
                 ,&prime     ,sc 
                 ,pMesh      ,opt
                 ,preName    ,nameOut
                 ,fileOut   );
/*...................................................................*/

/*...*/
      tm.solvEdpFluid = getTimeC() - tm.solvEdpFluid;
/*...................................................................*/
      if (!mpiVar.myId) fprintf(fileLogExc,"%s\n\n", DIF);
    }
/*===================================================================*/

/*===================================================================*
 * macro: vprop: propriedades variaveis                             
 *===================================================================*/
    else if((!strcmp(word,macro[38]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
      readPropVar(&propVarFluid,fileIn);
/*...................................................................*/
      if(!mpiVar.myId )fprintf(fileLogExc,"%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: gravity : gravidade                                                 
 *===================================================================*/
    else if((!strcmp(word,macro[39]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }   
      readGravity(gravity,fileIn);
/*...................................................................*/
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: model: modelos utilizados nas equacao diferencias         
 *===================================================================*/
    else if((!strcmp(word,macro[40]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
      readModel(&eModel,&turbModel
              , &eMass ,&ModelMomentum
              , fileIn);
/*...................................................................*/
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: mean : calcula a media                                                        
 *===================================================================*/
    else if(!strcmp(word,macro[41])){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }    
      calMean(&media,mesh,sc.ddt.dt[TIME_N_MINUS_1]
            , sc.ddt.t,sc.ddt.timeStep);
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: setMean : configuracao das medias temporais                                   
 *===================================================================*/
    else if(!strcmp(word,macro[42])){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
      readMean(&m,fileIn,mesh,&media);
/*...................................................................*/
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: save : salva o estado do programa                                          
 *===================================================================*/
    else if( (!strcmp(word,macro[43])) &&
             (save.step[1]++ == save.step[0]) ){
      if(!mpiVar.myId ){
        printf("%s\n",DIF);
        printf("%s\n",word);
      }
      save.step[1] = 1;
      wSave(&propVarFluid,&turbModel
          ,&thDynamic
          ,mesh         ,&sc.ddt
          ,&media
          ,preName      ,nameOut);
/*...................................................................*/
      if(!mpiVar.myId ) printf("%s\n\n",DIF);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: load : carrega um estodo do sistema especifico                             
 *===================================================================*/
    else if((!strcmp(word,macro[44]))){
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"%s\n",DIF);
        fprintf(fileLogExc,"%s\n",word);
      }
      save.fLoad = true;
      load(&propVarFluid,&turbModel
          ,&thDynamic
          ,mesh         ,&sc.ddt
          ,&media
          ,fileIn);
/*...................................................................*/
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n\n",DIF);
    }   
/*===================================================================*/


/*===================================================================*/
  }while(macroFlag && (!feof(fileIn)));

  mpiStop();
  return EXIT_SUCCESS;
}    
/*********************************************************************/


/*********************************************************************/
void testeFace(Geom *geom, Face *face
              , INT *cellFace, short *nFace
              , INT numel, short maxViz
              , short ndm, short maxNo) {
  short i, j;
  INT nel, idFace, cellOwner;
  DOUBLE lG[6], lF[6];
  DOUBLE lG1[6][3], lF1[6][3];



  for (i = 0; i < 6; i++) {
    lG[i] = lF[i] = 0.0;
    for (j = 0; j < 3; j++)
      lG1[i][j] = lF1[i][j] = 0.0;
  }

  for (nel = 0; nel < numel; nel++) {

    printf("volume\nnel %d\n", nel + 1);
    printf("%lf\n",geom->volume[nel]);
 
    printf("cc\nnel %d\n", nel + 1);
    printf("%lf %lf %lf\n", MAT2D(nel, 0,geom->cc,ndm)
                          , MAT2D(nel, 1, geom->cc, ndm)
                          , MAT2D(nel, 2, geom->cc, ndm));

    for (i = 0; i < nFace[nel]; i++) {
      lG[i]  = MAT2D(nel, i, geom->mksi, maxViz);
      idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
      lF[i] = face->mksi[idFace];
    }
    printf("mksi\nnel %d\n", nel + 1);
    printf("%lf %lf %lf %lf %lf %lf\n"
          , lG[0], lG[1], lG[2]
          , lG[3], lG[4], lG[5]);
    printf("%lf %lf %lf %lf %lf %lf\n"
          , lF[0], lF[1], lF[2]
          , lF[3], lF[4], lF[5]); 

    for (i = 0; i < nFace[nel]; i++) {
      lG[i] = MAT2D(nel, i, geom->fArea, maxViz);
      idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
      lF[i] = face->area[idFace];
    }
    printf("area\nnel %d\n", nel + 1);
    printf("%lf %lf %lf %lf %lf %lf\n"
      , lG[0], lG[1], lG[2]
      , lG[3], lG[4], lG[5]);
    printf("%lf %lf %lf %lf %lf %lf\n"
      , lF[0], lF[1], lF[2]
      , lF[3], lF[4], lF[5]);

    for (i = 0; i < nFace[nel]; i++) {
      lG[i] = MAT2D(nel, i, geom->mvSkew, maxViz);
      idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
      lF[i] = face->mvSkew[idFace];
    }
    printf("mvSkew\nnel %d\n", nel + 1);
    printf("%lf %lf %lf %lf %lf %lf\n"
      , lG[0], lG[1], lG[2]
      , lG[3], lG[4], lG[5]);
    printf("%lf %lf %lf %lf %lf %lf\n"
      , lF[0], lF[1], lF[2]
      , lF[3], lF[4], lF[5]);

    for (i = 0; i < nFace[nel]; i++) {
      idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
      cellOwner = MAT2D(idFace, 0, face->owner, 2) - 1;
      for (j = 0; j < ndm; j++) {
        lG1[i][j] = MAT3D(nel, i, j, geom->ksi, maxViz, ndm);
        lF1[i][j] = OWNER(cellOwner, nel)*MAT2D(idFace, j, face->ksi, ndm);
      }
    }
    printf("ksi\nnel %d\n", nel + 1);
    printf("(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n\n"
      , lG1[0][0], lG1[0][1], lG1[0][2]
      , lG1[1][0], lG1[1][1], lG1[1][2]
      , lG1[2][0], lG1[2][1], lG1[2][2]
      , lG1[3][0], lG1[3][1], lG1[3][2]
      , lG1[4][0], lG1[4][1], lG1[4][2]
      , lG1[5][0], lG1[5][1], lG1[5][2]);

    printf("(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      , lF1[0][0], lF1[0][1], lF1[0][2]
      , lF1[1][0], lF1[1][1], lF1[1][2]
      , lF1[2][0], lF1[2][1], lF1[2][2]
      , lF1[3][0], lF1[3][1], lF1[3][2]
      , lF1[4][0], lF1[4][1], lF1[4][2]
      , lF1[5][0], lF1[5][1], lF1[5][2]);

    for (i = 0; i < nFace[nel]; i++) {
      idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
      cellOwner = MAT2D(idFace, 0, face->owner, 2) - 1;
      for (j = 0; j < ndm; j++) {
        lG1[i][j] = MAT3D(nel, i, j, geom->eta, maxViz, ndm);
        lF1[i][j] = OWNER(cellOwner, nel)*MAT2D(idFace, j, face->eta, ndm);
      }
    }
    printf("eta\nnel %d\n", nel + 1);
    printf("(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n\n"
      , lG1[0][0], lG1[0][1], lG1[0][2]
      , lG1[1][0], lG1[1][1], lG1[1][2]
      , lG1[2][0], lG1[2][1], lG1[2][2]
      , lG1[3][0], lG1[3][1], lG1[3][2]
      , lG1[4][0], lG1[4][1], lG1[4][2]
      , lG1[5][0], lG1[5][1], lG1[5][2]);

    printf("(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      , lF1[0][0], lF1[0][1], lF1[0][2]
      , lF1[1][0], lF1[1][1], lF1[1][2]
      , lF1[2][0], lF1[2][1], lF1[2][2]
      , lF1[3][0], lF1[3][1], lF1[3][2]
      , lF1[4][0], lF1[4][1], lF1[4][2]
      , lF1[5][0], lF1[5][1], lF1[5][2]);
    for (i = 0; i < nFace[nel]; i++) {
      idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
      cellOwner = MAT2D(idFace, 0, face->owner, 2) - 1;
      for (j = 0; j < ndm; j++) {
        lG1[i][j] = MAT3D(nel, i, j, geom->normal, maxViz, ndm);
        lF1[i][j] = OWNER(cellOwner, nel)*MAT2D(idFace, j, face->normal, ndm);
      }
    }
    printf("normal\nnel %d\n", nel + 1);
    printf("(%9lf %9lf %9lf)\n"
           "(%9lf %9lf %9lf)\n"
           "(%9lf %9lf %9lf)\n"
           "(%9lf %9lf %9lf)\n"
           "(%9lf %9lf %9lf)\n"
           "(%9lf %9lf %9lf)\n\n"
           , lG1[0][0], lG1[0][1], lG1[0][2]
           , lG1[1][0], lG1[1][1], lG1[1][2]
           , lG1[2][0], lG1[2][1], lG1[2][2]
           , lG1[3][0], lG1[3][1], lG1[3][2]
           , lG1[4][0], lG1[4][1], lG1[4][2]
           , lG1[5][0], lG1[5][1], lG1[5][2]);

    printf("(%9lf %9lf %9lf)\n"
           "(%9lf %9lf %9lf)\n"
           "(%9lf %9lf %9lf)\n"
           "(%9lf %9lf %9lf)\n"
           "(%9lf %9lf %9lf)\n"
           "(%9lf %9lf %9lf)\n"
           , lF1[0][0], lF1[0][1], lF1[0][2]
           , lF1[1][0], lF1[1][1], lF1[1][2]
           , lF1[2][0], lF1[2][1], lF1[2][2]
           , lF1[3][0], lF1[3][1], lF1[3][2]
           , lF1[4][0], lF1[4][1], lF1[4][2]
           , lF1[5][0], lF1[5][1], lF1[5][2]);

    for (i = 0; i < nFace[nel]; i++) {
      idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
      cellOwner = MAT2D(idFace, 0, face->owner, 2) - 1;
      for (j = 0; j < ndm; j++) {
        lG1[i][j] = MAT3D(nel, i, j, geom->xm, maxViz, ndm);
        lF1[i][j] = MAT2D(idFace, j, face->xm, ndm);
      }
    }
    printf("xm\nnel %d\n", nel + 1);
    printf("(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n\n"
      , lG1[0][0], lG1[0][1], lG1[0][2]
      , lG1[1][0], lG1[1][1], lG1[1][2]
      , lG1[2][0], lG1[2][1], lG1[2][2]
      , lG1[3][0], lG1[3][1], lG1[3][2]
      , lG1[4][0], lG1[4][1], lG1[4][2]
      , lG1[5][0], lG1[5][1], lG1[5][2]);

    printf("(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      , lF1[0][0], lF1[0][1], lF1[0][2]
      , lF1[1][0], lF1[1][1], lF1[1][2]
      , lF1[2][0], lF1[2][1], lF1[2][2]
      , lF1[3][0], lF1[3][1], lF1[3][2]
      , lF1[4][0], lF1[4][1], lF1[4][2]
      , lF1[5][0], lF1[5][1], lF1[5][2]);

    for (i = 0; i < nFace[nel]; i++) {
      idFace = MAT2D(nel, i, cellFace, maxViz) - 1;
      cellOwner = MAT2D(idFace, 0, face->owner, 2) - 1;
      for (j = 0; j < ndm; j++) {
        lG1[i][j] = MAT3D(nel, i, j, geom->vSkew, maxViz, ndm);
        lF1[i][j] = MAT2D(idFace, j, face->vSkew, ndm);
      }
    }

    printf("vSkew\nnel %d\n", nel + 1);
    printf("(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n\n"
      , lG1[0][0], lG1[0][1], lG1[0][2]
      , lG1[1][0], lG1[1][1], lG1[1][2]
      , lG1[2][0], lG1[2][1], lG1[2][2]
      , lG1[3][0], lG1[3][1], lG1[3][2]
      , lG1[4][0], lG1[4][1], lG1[4][2]
      , lG1[5][0], lG1[5][1], lG1[5][2]);

    printf("(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      "(%9lf %9lf %9lf)\n"
      , lF1[0][0], lF1[0][1], lF1[0][2]
      , lF1[1][0], lF1[1][1], lF1[1][2]
      , lF1[2][0], lF1[2][1], lF1[2][2]
      , lF1[3][0], lF1[3][1], lF1[3][2]
      , lF1[4][0], lF1[4][1], lF1[4][2]
      , lF1[5][0], lF1[5][1], lF1[5][2]);
  }

}
/***************************************************************************/
 
