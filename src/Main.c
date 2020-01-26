/********************* includes padroes ******************************/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
/*********************************************************************/

/*********************************************************************/
#include<Adjcency.h>
#include<CellLoop.h>
#include<Coo.h>
#include<Combustion.h>
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
#include<TimeInterpol.h>
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

static void initSec(char *word,short icod) 
{

  switch(icod)
  {
    case OUTPUT_FOR_FILE:
      fprintf(fileLogExc, "%s\n", DIF);
      fprintf(fileLogExc, "%s\n", word);
    break;
    case OUTPUT_FOR_SCREEN:
      if (!mpiVar.myId)
      {
        printf("%s\n", DIF);
        printf("%s\n", word);
      }
    break;
  }
}

static void endSec(short icod)
{
  switch (icod) 
  {
    case OUTPUT_FOR_FILE:
      fprintf(fileLogExc, "%s\n\n", DIF);
    break;
    case OUTPUT_FOR_SCREEN:
      if (!mpiVar.myId) printf("%s\n\n", DIF);
    break;
  }
}

/*********************************************************************/
void grad(Mesh *mesh);
void gradErro(Mesh *mesh);
void vel(Mesh *mesh);
void initialZ(Mesh *mesh);
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
  SistEq sistEqD1, sistEqT1;
  SistEq sistEqVel, sistEqPres, sistEqEnergy, sistEqKturb;
  SistEq sistEqComb;
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
  MomentumModel momentumModel;
  DiffModel diffModel[3];
  TransModel transModel[3];
  Combustion combModel;
/*... propriedade variaveis*/
  PropVarFluid propVarFluid;
  PropVarCD propVarD[3],propVarT[3];
/*...*/
  TimeInterpol tInterpol;

/*... solver*/
  Solv solvD1, solvT1;
  Solv solvVel,solvPres,solvEnergy,solvKturb, solvComb;
  bool fSolvD1 = false, fSolvT1 = false;
  bool fSolvVel = false,fSolvPres = false, fSolvEnergy = false;
  bool fSolvComb = false;
  bool fSolvKturb = false;
  bool fSolvSimple = false,fSolvPrime = false;  
/*... reordenacao da malha*/
  Reord  reordMesh;

/*... particionamento*/
  PartMesh *pMesh = NULL;
/*...*/
  Save save;

/*...*/
  Macros mm;
//bool fReadModels=false;

/*... arquivo*/
  char nameIn[MAX_STR_LEN_IN], nameOut[SIZEMAX];
  char auxName[MAX_STR_LEN_SUFIXO], preName[MAX_STR_LEN_SUFIXO];
  FILE *fileIn=NULL,*fileOut=NULL,*fileLog=NULL;
  FileOpt opt;

/*...*/
  Scheme sc;
/*...................................................................*/

/*... loop nas celulas*/
/*Lib lib;*/
  
/* ... macro camandos de leitura*/
  bool macroFlag; 
  char word[WORD_SIZE];
  char macro[][WORD_SIZE] =
  {"help"         ,"mesh"         ,"stop"           /* 0, 1, 2*/
  ,"config"       ,"nextLoop"     ,"rcGrad"         /* 3, 4, 5*/
  ,"pgeo"         ,"pcoob"        ,"pcoo"           /* 6, 7, 8*/ 
  ,"setSolvDiff"  ,"setSolvTrans" ,"openmp"         /* 9,10,11*/
  ,"solvD1"       ,"propVar"      ,"pD1"            /*12,13,14*/
  ,"nlIt"         ,"pD1CsvCell"   ,"pD1CsvNode"     /*15,16,17*/
  ,"solvT1"       ,"edo"          ,"pT1"            /*18,19,20*/
  ,"simpleLm"     ,"pT1CsvCell"   ,"pT1CsvNode"     /*21,22,23*/
  ,"setSolvFluid" ,"simple"       ,"setSimple"      /*24,25,26*/
  ,"transient"    ,"timeUpdate"   ,"partd"          /*27,28,29*/
  ,"advection"    ,"edp"          ,"diffusion"      /*30,31,32*/
  ,"pFluid"       ,"setPrint"     ,"reScaleMesh"    /*33,34,35*/
  ,"setPrime"     ,"prime"        ,""               /*36,37,38*/
  ,"setSolvComb"  ,"pCombustion"  ,"simpleComb"     /*39,40,41*/
  ,"pRes"         ,"chemical"     ,"residual"       /*42,43,44*/ 
  ,"gravity"      ,"model"        ,"mean"           /*45,46,47*/
  ,"setMean"      ,"save"         ,"load"};         /*48,49,50*/
/* ..................................................................*/

/*... Memoria principal(valor padrao - bytes)*/
  nmax = 200000;
/* ..................................................................*/

  xRef[0] =  0.e0;
  xRef[1] =  0.e0;
  xRef[2] =  0.e0;

/*...*/
  mm.flWord = false;
/*...*/
  diffModel[0].fRes = true;
/*...*/
  media.fMedia = false;
  media.fVel   = false;
/*..................................................................*/

/*...*/
  combModel.fCombustion     = false;
  combModel.nComb           = 0;
  combModel.fRes            = true;
  combModel.fCorrectVel     = false;
  combModel.typeHeatRealese = HFORMATION; 
  combModel.edc.tMix        = 0.125;
  combModel.edc.cGamma      = 1.01;
  combModel.edc.cTau        = 2.1377;
  combModel.edc.type        = PANJWANI_CONST_TMIX_EDC;
//combModel.edc.type        = FLUENT_CONST_TMIX_EDC;
  combModel.totalHeat       = 0.e0;
  combModel.edc.edo.maxIt   = 1000;
  combModel.edc.edo.tol     = 1.0e-05;
  combModel.edc.edo.type    = EDO_SIE;
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
  turbModel.fOneEq              = false;
  turbModel.fTurbStruct         = false;
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
  turbModel.SchmidtTsgs         = 0.5e0;
/*....................................................................*/

/*...*/
  eModel.fPresWork     = false;
  eModel.fDissipation  = false;
  eModel.fRes          = true;
  eModel.fTemperature  = false;
  eModel.fKelvin       = false;
/*...................................................................*/

/*...*/
  momentumModel.fRes             = true;
  momentumModel.fRhieChowInt     = false;
  momentumModel.iCodBuoyant      = BUOYANT_RHOREF;
  momentumModel.fViscosity       = true;
  momentumModel.fDiv             = true;
  momentumModel.fSoPressure     = true;
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

  ompVar.nThreadsGrad = 1;
  ompVar.fGrad        = false;

  ompVar.nThreadsReaction = 1;
  ompVar.fReaction        = false;
/* ..................................................................*/

/* ... opcoes de arquivos */
  initPrintVtk(&opt);  
/* ..................................................................*/

/*... propriedades variaveis*/
  propVarFluid.fDensity             = false;
  propVarFluid.fSpecificHeat        = false;
  propVarFluid.fDynamicViscosity    = false;
  propVarFluid.fThermalConductivity = false;
  propVarFluid.fDiffusion           = false;
/*...................................................................*/

/*... propriedades variaveis*/
  initPropStructCD(propVarD,3);
  initPropStructCD(propVarT,3);
/*...................................................................*/

/*...*/
  initLoads();
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
//mesh0 = (Mesh*) malloc(sizeof(Mesh));
  mesh0 = (Mesh*) calloc(1,sizeof(Mesh));
  ERRO_MALLOC(mesh0,"mesh0",__LINE__,__FILE__,__func__);
  mesh0->mass[2] =  mesh0->mass[1] = mesh0->mass[0] = 0.e0;
  mesh0->scaleX[0] = mesh0->scaleX[1] = mesh0->scaleX[2] = 1.0; 
/*... tecnica padrao de resconstrucao de gradiente*/
  sc.rcGrad.type      = RCGRADGAUSSN;
  sc.rcGrad.func      = GL_BARTH_MOD;
  sc.rcGrad.beta      = 1.0;
  sc.rcGrad.fLimiter  = false;
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
  sc.advVel.iCod2    = MIDMODFACE;
  sc.advVel.par[0]   = 0.e0;
  sc.advEnergy.iCod1 = TVD;
  sc.advEnergy.iCod2 = MIDMODFACE;
  sc.advVel.par[0]   = 0.e0;
  sc.advKturb.iCod1  = TVD;
  sc.advKturb.iCod2  = MIDMODFACE;
  sc.advVel.par[0]   = 0.e0;
/*...................................................................*/

/*...*/  
  sc.diffVel.iCod    = OVERRELAXED;
  sc.diffPres.iCod   = OVERRELAXED;
  sc.diffEnergy.iCod = OVERRELAXED;
  sc.diffKturb.iCod  = OVERRELAXED;
/*...................................................................*/

/*...*/
  sc.diffComb.iCod = OVERRELAXED;
  sc.advComb.iCod1 = TVD;
//sc.advComb.iCod1 = FOUP;
  sc.advComb.iCod2 = SUPERBEEFACE;
  sc.advComb.par[0] = 0.e0;
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
  reordMesh.flag = false; 
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
  fName(preName,mpiVar.myId,0, 23 ,nameOut);
  fileLogExc = openFileBuffer(nameOut,"w",false);
/*...................................................................*/

/*...*/
  fName(preName,mpiVar.myId,0, 31 ,nameOut);
  fileLogDebug = openFileBuffer(nameOut,"w",false);
/*...................................................................*/

/*loop de execucao*/
  macroFlag = true;
  do{
/*... macros na marcro trasient*/
    if(mm.flWord){
      if( mm.jLoop > mm.kLoop)
      { 
        ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__,
                   "Numero de comandos na string trasient execedido"
                  ,EXIT_PROG); 
      }
      strcpy(word,mm.loopWord[mm.jLoop]);
      mm.jLoop++;
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
      initSec(word, OUTPUT_FOR_SCREEN);
      help(fileIn);
      endSec(OUTPUT_FOR_SCREEN);
    }
/*===================================================================*/

/*===================================================================*
 * macro: mesh - leitura da malha e inicializa das estruturas        *
 * de resolucao do problema                                          * 
 *===================================================================*/
    else if((!strcmp(word,macro[1]))){
      initSec(word, OUTPUT_FOR_FILE);
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
        readFileFvMesh(&m           ,mesh0
                      ,&propVarFluid
                      ,propVarD     ,propVarT
                      ,&eModel      ,&turbModel 
                      ,&combModel   ,&momentumModel
                      ,&media       
                      ,fileIn);
      mpiWait(); 
/*...................................................................*/

/*... calcula a vizinhaca do elementos da malha completa*/
      if(!mpiVar.myId)
        neighbors(&m                 ,mesh0->elm.node 
                 ,mesh0->elm.adj.nelcon,mesh0->elm.adj.nViz
                 ,mesh0->elm.nen       ,&mesh0->nFaces
                 ,mesh0->nnode         ,mesh0->numel     
                 ,mesh0->maxNo         ,mesh0->maxViz   
                 ,mesh0->ndm);        
/*...................................................................*/

/*... identifica parede impermevais (malha completa)*/
      if(!mpiVar.myId)
      {
        if(mesh0->ndfF > 0 || mesh0->ndfFt > 0)
        {
          wallFluidVelPres(loadsVel
                   ,loadsPresC           ,loadsPres           
                   ,mesh0->elm.faceRvel  ,mesh0->elm.faceRpres
                   ,mesh0->elm.adj.nelcon,mesh0->elm.adj.nViz   
                   ,mesh0->numel         ,mesh0->maxViz
                   ,momentumModel.iCodBuoyant); 
          if(turbModel.fOneEq)
            wallFluid(loadsKturb                  
                     ,mesh0->elm.faceReKturb,mesh0->elm.adj.nelcon
                     ,mesh0->elm.adj.nViz   
                     ,mesh0->numel          ,mesh0->maxViz);   
/*... verifica se o dominio e aberto ou nao*/
          mesh0->fOpen = openDomain(loadsVel
                         , mesh0->elm.faceRvel , mesh0->elm.adj.nViz
                         , mesh0->numelNov ,mesh0->maxViz);
/*...*/
          thDynamic.fDensityRef = !mesh0->fOpen;          
//        thDynamic.fDensityRef = true;
          thDynamic.fPresTh     = !mesh0->fOpen;
/*...................................................................*/   
        }
      }
#ifdef _MPI_
      MPI_Bcast(&thDynamic.fDensityRef,1,MPI_C_BOOL,0,mpiVar.comm);
      MPI_Bcast(&thDynamic.fPresTh    ,1,MPI_C_BOOL,0,mpiVar.comm);
#endif
/*...................................................................*/

/*... particionamento da malha*/
      if(pMesh->fPartMesh && mpiVar.nPrcs > 1)
      {
        tm.partdMesh = getTimeC() - tm.partdMesh;
        if(!mpiVar.myId)
        {
/*...*/
          fprintf(fileLogExc, "%s\n", "PartMesh");
          fPartMesh(&m            ,mesh0->node.x  
                  ,mesh0->elm.node,mesh0->elm.nen
                  ,mesh0->nnode   ,mesh0->numel
                  ,pMesh   
                  ,mesh0->ndm     ,mesh0->maxNo 
                  ,mpiVar.nPrcs);
/*...................................................................*/

/*... */
          if(pMesh->fPrintMesh)
          {
            fprintf(fileLogExc, "%s\n", "WritePartMesh");
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
/*...................................................................*/

/*...*/
//      mesh = (Mesh*) malloc(sizeof(Mesh));
        mesh = (Mesh*) calloc(1,sizeof(Mesh));
        ERRO_MALLOC(mesh,"mesh",__LINE__,__FILE__,__func__);
/*...*/
        tm.partdMeshCom = getTimeC() - tm.partdMeshCom;
        comunicateMesh(&m        ,&combModel
                      ,&turbModel
                      ,mesh0     ,mesh
                      ,pMesh
                      ,loadsD1   ,loadsT1
                      ,loadsVel  ,loadsPres 
                      ,loadsPresC,loadsEnergy
                      ,loadsTemp ,loadsKturb
                      ,loadsZcomb);
        tm.partdMeshCom = getTimeC() - tm.partdMeshCom;
/*...................................................................*/

/*... individuais*/
        if(pMesh->fPrintMeshPart)
        {
          fprintf(fileLogExc, "%s\n", "WritePartMeshes");
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

/*...*/
      tm.adjcency = getTimeC();
/*... calcula a vizinhaca do elementos*/
      neighbors(&m                  ,mesh->elm.node 
               ,mesh->elm.adj.nelcon,mesh->elm.adj.nViz
               ,mesh->elm.nen       ,&mesh->nFaces
               ,mesh->nnode         ,mesh->numel     
               ,mesh->maxNo         ,mesh->maxViz   
               ,mesh->ndm);        
/*...................................................................*/
      
/*...*/        
      faceStruct(&m, mesh);
/*...................................................................*/
      tm.adjcency = getTimeC() - tm.adjcency;
/*...................................................................*/
     
/*... calculo de propriedades geometricas recorrentes*/
      tm.geom = getTimeC() - tm.geom;
      pGeomForm(mesh->node.x         ,mesh->elm.node
               ,mesh->elm.adj.nelcon ,mesh->elm.adj.nViz
               ,mesh->elm.geomType   ,mesh->elm.nen        
               ,mesh->elm.cellFace
               ,&mesh->elm.geom      ,&mesh->face 
               ,&pMesh->iEl   
               ,mesh->maxNo          ,mesh->maxViz  
               ,mesh->ndm            ,mesh->numelNov);
      tm.geom = getTimeC() - tm.geom;
/*    testeFace(&mesh->elm.geom    , &mesh->face
               , mesh->elm.cellFace, mesh->elm.adj.nViz
               , mesh->numel       , mesh->maxViz
               , mesh->ndm         , mesh->maxNo);*/
/*...................................................................*/

/*... geom(Cel)*/
      dGlobalCel(&m                  ,pMesh
                ,mesh0->elm.geom.cc  ,mesh->elm.geom.cc
                ,mesh->numelNov 
                ,mesh->ndm           ,1);
/*...................................................................*/

/*...*/
      if(mesh->ndfFt > 0)
      {
        if(combModel.fCombustion)
          propVarFluid.molarMass = mixtureMolarMassMedMpi(&combModel
                              ,mesh->elm.yFrac
                              ,mesh->elm.geom.volume,mesh->numelNov);
/*...*/
//      if(thDynamic.fDensityRef)
          propVarFluid.densityRef = specificMassRef(mesh->elm.densityFluid.t
                                                , mesh->elm.geom.volume                  
                                                , mesh->numelNov);
/*...................................................................*/

/*...*/
        if(thDynamic.fPresTh)
          initPresRef(mesh->elm.temp        , mesh->elm.geom.volume  
                    , thDynamic.pTh         , propVarFluid.densityRef
                    , propVarFluid.molarMass                   
                    , mesh->numelNov        , eModel.fKelvin);
/*...................................................................*/

/*...*/
        mesh->mass[0] = totalMass(mesh->elm.densityFluid.t
                              , mesh->elm.geom.volume
                              , mesh->numelNov); 
        mesh->mass[1] =  mesh->mass[0];
/*....................................................................*/

/*... gera a pressao inicial hidrostatica*/
        if(momentumModel.iCodBuoyant == BUOYANT_HYDROSTATIC)
          hPres(mesh->elm.pressure0     , mesh->elm.pressure
              , mesh->elm.densityFluid.t, mesh->elm.geom.cc
              , gravity                 , mesh->xRef
              , mesh->numel             , mesh->ndm );      
/*...................................................................*/
      }
/*...................................................................*/  
 
/*... reconstrucao de gradiente least square*/
      if(sc.rcGrad.type ==  RCLSQUARE || sc.rcGrad.type  ==  RCLSQUAREQR){
        if(!mpiVar.myId ){
          fprintf(fileLogExc,"%s\n",DIF);
          fprintf(fileLogExc,"Least Square ...\n");
        }
/*... wleastSqaure*/
        HccaAlloc(DOUBLE,&m,mesh->elm.leastSquare
                 ,mesh->numelNov*mesh->maxViz*mesh->ndm
                 ,"leastSquare",_AD_);
/*... leastSqaure QR*/
        if(sc.rcGrad.type  ==  RCLSQUAREQR){
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
        rcLeastSquare(mesh->elm.cellFace   , mesh->face.owner
                    , mesh->face.mksi      , mesh->face.ksi
                    , mesh->elm.leastSquare, mesh->elm.leastSquareR
                    , mesh->elm.adj.nViz       
                    , mesh->numelNov       , mesh->maxViz
                    , sc.rcGrad.type       , mesh->ndm);
        tm.leastSquareMatrix = getTimeC() - tm.leastSquareMatrix;
/*...................................................................*/
        if(!mpiVar.myId ){
          fprintf(fileLogExc,"Least Square.\n");
          fprintf(fileLogExc,"%s\n",DIF);
        }
      }
/*...................................................................*/

/*... qualidade da malha*/
      meshQuality(&mesh->mQuality        , mesh->elm.cellFace     
                , mesh->face.owner       , mesh->elm.adj.nViz  
                , mesh->elm.geom.volume  , mesh->face.area
                , mesh->face.ksi         , mesh->face.normal
                , mesh->face.mvSkew      , mesh->elm.geom.dcca
                , mesh->maxViz           , mesh->ndm
                , mesh->numelNov);
/*... qualidade da malha global*/
      if(mpiVar.nPrcs > 1)
        globalMeshQuality(&mesh->mQuality,&mesh0->mQuality);
/*...................................................................*/
         
/*... reodenando as celulas para dimuincao da banda*/
      HccaAlloc(INT,&m,reordMesh.num,mesh->numel,"rNum" ,_AD_);
      if(!mpiVar.myId ) fprintf(fileLogExc,"%s\n",DIF);
      if(!mpiVar.myId ) fprintf(fileLogExc,"Reordenando a malha ...\n");
      tm.reord = getTimeC() - tm.reord;
      reord(&m                ,reordMesh.num,mesh->elm.adj.nelcon
           ,mesh->elm.adj.nViz,mesh->maxViz  
           ,mesh->numel       ,mesh->numelNov
           ,reordMesh.flag   ,mpiVar.nPrcs);
      tm.reord = getTimeC() - tm.reord;
      if(!mpiVar.myId ){
        fprintf(fileLogExc,"Malha reordenada.\n");
        fprintf(fileLogExc,"%s\n\n",DIF);
      }
/*...................................................................*/

/*...*/
//    initialZ(mesh0);
//    vel(mesh0);
//    grad(mesh0);
//    writeMeshPart(mesh,&combModel);
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: stop : finalizacao do programa
 *===================================================================*/
    else if((!strcmp(word,macro[2]))){
      initSec(word, OUTPUT_FOR_FILE);
      tm.total = getTimeC() - tm.total;

/*... */
      fName(preName,mpiVar.nPrcs,mpiVar.myId,7,nameOut);
      fileLog = openFile(nameOut,"w");
      writeLog(mesh           ,&sc
              ,&solvD1        ,&sistEqD1
              ,&solvT1        ,&sistEqT1
              ,&solvVel       ,&sistEqVel
              ,&solvPres      ,&sistEqPres
              ,&solvEnergy    ,&sistEqEnergy 
              ,&solvComb      ,&sistEqComb 
              ,&tm             
              ,fSolvD1        ,fSolvT1     
              ,fSolvVel       ,fSolvPres 
              ,fSolvEnergy    ,turbModel.fTurb  
              ,fSolvComb      ,ompVar
              ,nameIn         ,fileLog);
      fclose(fileLog);
/*...................................................................*/

/*... medias do tempo dos processos Mpi*/
      if(mpiVar.nPrcs > 1) 
      {
        if(!mpiVar.myId){
          fName(preName,mpiVar.nPrcs,mpiVar.myId,60,nameOut);
          fileLog = openFile(nameOut,"w");
        }
        writeLogMeanTime(mesh0          ,&sc
              ,&solvD1        ,&sistEqD1
              ,&solvT1        ,&sistEqT1
              ,&solvVel       ,&sistEqVel
              ,&solvPres      ,&sistEqPres
              ,&solvEnergy    ,&sistEqEnergy 
              ,&solvComb      ,&sistEqComb 
              ,&tm            ,&ompVar 
              ,fSolvD1        ,fSolvT1     
              ,fSolvVel       ,fSolvPres 
              ,fSolvEnergy    ,turbModel.fTurb  
              ,fSolvComb 
              ,nameIn         ,fileLog);        
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
      if(fSolvT1 && solvT1.log && !mpiVar.myId)  
        fclose(solvT1.fileSolv);
/*...................................................................*/

/*... fechando o arquivo do log nao linear D1*/      
      if(fSolvD1 && opt.fItPlot && !mpiVar.myId)  
        fclose(opt.fileItPlot[FITPLOTD1]);
/*... fechando o arquivo do log nao linear T1*/      
      if(fSolvT1 && opt.fItPlot && !mpiVar.myId)  
        fclose(opt.fileItPlot[FITPLOTT1]);
/*... fechando o arquivo do log nao linear do simple */      
      if(fSolvSimple && opt.fItPlot && !mpiVar.myId)
      {
        fclose(opt.fileItPlot[FITPLOTSIMPLE]);
        fclose(opt.fileParameters);
      }

      fclose(fileLogDebug);
      finalizeMem(&m,false);
      macroFlag = false;
      endSec(OUTPUT_FOR_FILE);
      fclose(fileLogExc);
    }    
/*===================================================================*/

/*===================================================================*
 * macro: config : configuracao basica de excucao
 *===================================================================*/
    else if((!strcmp(word,macro[3]))){
      initSec(word, OUTPUT_FOR_FILE);
      config(&opt,&reordMesh ,fileIn);      
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: nextLoop : proximo loop tempo            
 *===================================================================*/
    else if((!strcmp(word,macro[4]))){
      mm.jLoop = 0;
    }   
/*===================================================================*/

/*===================================================================*
* macro: rcGrad: tipo de reconstrucao de gradeiente
*===================================================================*/
    else if ((!strcmp(word, macro[5]))){
      initSec(word, OUTPUT_FOR_FILE);
      setReGrad(&sc.rcGrad, fileIn);
      endSec(OUTPUT_FOR_FILE);
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
               ,mesh0->elm.faceRd1       ,mesh0->elm.faceRt1       
               ,mesh0->elm.faceRvel      ,mesh0->elm.faceRenergy  
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
        printFace(&m      , mesh0
              , preName , combModel.fCombustion
              , opt.bVtk, fileOut);
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
//      fName(preName,0,0,14,nameOut);
//      for(i=0;i<sistEqD1.neq;i++)
//        sistEqD1.b[i] /= sistEqD1.ad[i];
//
//      writeCooB(sistEqD1.b,sistEqD1.neq,nameOut);
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
      initSec(word, OUTPUT_FOR_FILE);
      readSolvDiff(&m     , mesh     , &reordMesh
                  ,&solvD1, &sistEqD1, &fSolvD1
                  ,auxName, preName  , nameOut
                  ,fileIn , &opt);
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: setSolvTrans : problema de transporte  
 *===================================================================*/
    else if((!strcmp(word,macro[10]))){
      initSec(word, OUTPUT_FOR_FILE);
/*... inicializando a estrutura de equacoes do problema*/
      readSolvTrans(&m     , mesh     , &reordMesh
                  , &solvT1, &sistEqT1, &fSolvT1
                  , auxName, preName  , nameOut
                  , fileIn , &opt);
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
* macro: openmp: configuracao do openmp  
*===================================================================*/
    else if ((!strcmp(word, macro[11]))) {
      initSec(word, OUTPUT_FOR_FILE);
/*... */
      openMpSet(fileIn,&ompVar);
/*..................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }
/*===================================================================*/

/*===================================================================*
 * macro: solvd1: problema de difusao pura
 *===================================================================*/
    else if((!strcmp(word,macro[12])))
    {
      initSec(word, OUTPUT_FOR_SCREEN);
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
      diffusion(&m          ,loadsD1        ,&diffModel[0]
               ,mesh0       ,mesh           ,&sistEqD1
               ,&solvD1     ,&sc            ,pMesh
               ,&propVarD[0],&opt           ,preName    
               ,nameOut     ,fileOut);
/*...................................................................*/

/*...*/
     tm.solvEdpD1    = getTimeC() - tm.solvEdpD1;
/*...................................................................*/
     endSec(OUTPUT_FOR_SCREEN);
    }
/*===================================================================*/

/*===================================================================*
  * macro: vprop: propriedades variaveis
*===================================================================*/
    else if ((!strcmp(word, macro[13])))
    {
      initSec(word, OUTPUT_FOR_FILE);
      readPropVar(&propVarFluid,propVarD,propVarT,&combModel,fileIn);
      endSec(OUTPUT_FOR_FILE);
    }
/*===================================================================*/

/*===================================================================*
 * macro: puD1 : escreve os arquivos dos resultados da uD1          
 *===================================================================*/
    else if((!strcmp(word,macro[14]) && 
            (opt.stepPlot[1]++) == opt.stepPlot[0]))
    {
/*...*/
      opt.stepPlot[1] = 1;
      initSec(word, OUTPUT_FOR_SCREEN);
      printDiff(&m
              , pMesh   , &sc
              , loadsD1 , &opt
              , mesh0   , mesh
              , preName , nameOut);
      endSec(OUTPUT_FOR_SCREEN);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: nlIt: configura das iteracoes nao lineares             
 *===================================================================*/
    else if((!strcmp(word,macro[15])))
    {
      initSec(word, OUTPUT_FOR_FILE);
      readNlIt(&sc, fileIn);
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: pD1CellCsv:imprime os resultados no formato csv                  
 *===================================================================*/
    else if((!strcmp(word,macro[16])))
    {
      initSec(word, OUTPUT_FOR_SCREEN);
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
      if(!mpiVar.myId )
      {    
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
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: pD1CsvNode:imprime os resultados no formato csv                  
 *===================================================================*/
    else if((!strcmp(word,macro[17])))
    {
      initSec(word, OUTPUT_FOR_SCREEN);
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
      if(!mpiVar.myId )
      {
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
      }
    }   
/*===================================================================*/

/*===================================================================*
 * macro: solvt1: problema de transporte  
 *===================================================================*/
    else if((!strcmp(word,macro[18])))
    {
      initSec(word, OUTPUT_FOR_SCREEN);
      mpiWait();
      tm.solvEdpT1    = getTimeC() - tm.solvEdpT1;
     
/*...*/
      transport(&m          ,loadsT1        ,&transModel[0]
               ,mesh0       ,mesh           ,&sistEqT1
               ,&solvT1     ,&sc            ,pMesh
               ,&propVarT[0],&opt           ,preName
               ,nameOut     ,fileOut);
/*...................................................................*/

/*...*/
     tm.solvEdpT1    = getTimeC() - tm.solvEdpT1;
/*...................................................................*/
     endSec(OUTPUT_FOR_SCREEN);
    }
/*===================================================================*/

/*===================================================================*
 * macro: edo: configura do solver de edo   
 *===================================================================*/
    else if((!strcmp(word,macro[19])))
    {
      initSec(word, OUTPUT_FOR_SCREEN);
/*...*/
      configEdo(&combModel.edc.edo,fileIn);
/*...................................................................*/
      endSec(OUTPUT_FOR_SCREEN);
    }
/*===================================================================*/

/*===================================================================*
 * macro: puT1 : escreve os arquivos dos resultados da uT1          
 *===================================================================*/
    else if((!strcmp(word,macro[20]) && 
            (opt.stepPlot[1]++) == opt.stepPlot[0]) )
    {
      opt.stepPlot[1] = 1;
      initSec(word, OUTPUT_FOR_SCREEN);
      printTrans(&m
               , pMesh   ,&sc
               , loadsT1 ,&opt
               , mesh0   ,mesh
               , preName ,nameOut);
      endSec(OUTPUT_FOR_SCREEN);
    }   
/*===================================================================*/

/*===================================================================*
* macro: simpleLm: escoamento de fluidos (SIMPLELM)
*===================================================================*/
    else if ((!strcmp(word, macro[21])))
    {
      initSec(word, OUTPUT_FOR_SCREEN);
      mpiWait();
      tm.solvEdpFluid = getTimeC() - tm.solvEdpFluid;
/*...*/
      if (!fSolvVel || !fSolvPres) {
        printf("Estrutara de dados nao montada para o solvFluid!!!\n");
        exit(EXIT_FAILURE);
      }
/*...................................................................*/

/*...*/
      if (!fSolvSimple) {
        printf("Simple nao configurado ainda!!!\n");
        exit(EXIT_FAILURE);
      }
/*...................................................................*/

/*...*/
      simpleSolverLm(&m           , &propVarFluid
                   , loadsVel     , loadsPres
                   , loadsEnergy  , loadsKturb
                   , &eModel      , &combModel
                   , &eMass       , &momentumModel
                   , &turbModel   , &thDynamic
                   , mesh0        , mesh
                   , &sistEqVel   , &sistEqPres
                   , &sistEqEnergy, &sistEqKturb
                   , &solvVel     , &solvPres
                   , &solvEnergy  , &solvKturb
                   , &simple      , &sc
                   , pMesh        , &media
                   , &opt         , preName
                   , nameOut      , fileOut);
/*...................................................................*/

/*...*/
      tm.solvEdpFluid = getTimeC() - tm.solvEdpFluid;
/*...................................................................*/
      endSec(OUTPUT_FOR_SCREEN);
    }
/*===================================================================*/

/*===================================================================*
 * macro: pT1CellCsv:imprime os resultados no formato csv                  
 *===================================================================*/
    else if((!strcmp(word,macro[22])))
    {
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
    else if((!strcmp(word,macro[23])))
    {
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
    else if((!strcmp(word,macro[24])))
    {
      initSec(word, OUTPUT_FOR_FILE);
/*...................................................................*/

/*...*/
      readSolvFluid(&m          , mesh           , &reordMesh
                   , &combModel 
                   , &solvVel   , &sistEqVel    , &fSolvVel
                   , &solvPres  , &sistEqPres   , &fSolvPres
                   , &solvEnergy, &sistEqEnergy , &fSolvEnergy
                   , &solvKturb , &sistEqKturb  , &fSolvKturb
                   , &solvComb  , &sistEqComb   , &fSolvComb 
                   , pMesh
                   , auxName    , preName       , nameOut
                   , fileIn                     , &opt);
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: simple: escoamento de fluidos (SIMPLE)
 *===================================================================*/
    else if((!strcmp(word,macro[25])))
    {
      initSec(word, OUTPUT_FOR_SCREEN);
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
     simpleSolver(&m         
                 ,loadsVel   ,loadsPres 
                 ,&eMass     ,&momentumModel
                 ,&turbModel 
                 ,mesh0      ,mesh           
                 ,&sistEqVel ,&sistEqPres
                 ,&solvVel   ,&solvPres
                 ,&simple
                 ,&sc        ,pMesh
                 ,&opt       ,preName        
                 ,nameOut    ,fileOut);
/*...................................................................*/

/*...*/
      tm.solvEdpFluid    = getTimeC() - tm.solvEdpFluid;
/*...................................................................*/
      endSec(OUTPUT_FOR_SCREEN);
    }
/*===================================================================*/

/*===================================================================*
 * macro: setSimple: configuracoe do metodo simple
 *===================================================================*/
    else if((!strcmp(word,macro[26])))
    {
    initSec(word, OUTPUT_FOR_FILE);
/*...*/
    if(!mpiVar.myId )
    {
      fName(preName,mpiVar.nPrcs,0,27 ,nameOut);
      opt.fileParameters = openFileBuffer(nameOut,"w",true);
      fprintf(opt.fileParameters,"%s %s %s\n"
                            ,"#step t cfl reynolds peclet P0 "
                            ,"mass(Inc) mass(Avg) massIn massOut "
                            ,"totalHeat Heat temMax temMed");
    }
/*...................................................................*/

/*...*/
      fSolvSimple = readSetSimple(&m     , fileIn
                                , mesh0  , mesh
                                , &simple);
/*...................................................................*/

      endSec(OUTPUT_FOR_FILE);
    }
/*===================================================================*/

/*===================================================================*
 * macro: transient:configuracao da discretizacao temporal                 
 *===================================================================*/
    else if((!strcmp(word,macro[27])))
    {
      initSec(word, OUTPUT_FOR_FILE);
/*...*/
      readTransConfig(&sc.ddt,&save,&mm,fileIn);
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
/*...*/
      if(opt.fPolt)
        printCall(&m        , &propVarFluid
              , &turbModel, &eModel 
              , &combModel, &tInterpol
              , pMesh     , &sc
              , loadsVel  , loadsPres 
              , loadsTemp , loadsZcomb
              , &opt       
              , mesh0     , mesh 
              , &media    , PINITIAL_TIME  
              , preName   , nameOut); 
/*...................................................................*/

    }   
/*===================================================================*/

/*===================================================================*
 * macro: timeUpdate : macro de atualizaco do tempo                       
 *===================================================================*/
    else if((!strcmp(word,macro[28])))
    {
      initSec(word, OUTPUT_FOR_SCREEN);
/*...*/
      updateTime(&sc.ddt,&opt, &mm, mpiVar.myId );
/*...................................................................*/

/*...*/
      if(opt.fPolt)
      {
        printCall(&m        , &propVarFluid
                , &turbModel, &eModel 
                , &combModel, &tInterpol
                , pMesh     , &sc
                , loadsVel  , loadsPres 
                , loadsTemp , loadsZcomb
                , &opt       
                , mesh0     , mesh 
                , &media    , PINT_TIME  
                , preName   , nameOut); 
      }
      updateTimeStruct(&m       ,&tInterpol
                     ,mesh
                     ,&combModel,&opt); 
/*...................................................................*/
      endSec(OUTPUT_FOR_SCREEN);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: partd particionamento da malha                                   
 *===================================================================*/
    else if((!strcmp(word,macro[29])))
    {
      initSec(word, OUTPUT_FOR_FILE);
/*...*/
      pMesh->fPartMesh = true;
      readMacro(fileIn,word,false);
      if(!strcmp(word,"config:"))
      {
/*... fPrintMesh*/        
        readMacro(fileIn,word,false);
        if(!strcmp(word,"true"))
        { 
          pMesh->fPrintMesh = true;
          if(!mpiVar.myId ) fprintf(fileLogExc,"fPrintMesh    : true\n");
        }

/*... fPrintMeshPart*/        
        readMacro(fileIn,word,false);
        if(!strcmp(word,"true"))
        { 
          pMesh->fPrintMeshPart = true;
          if(!mpiVar.myId ) fprintf(fileLogExc,"fPrintMeshPart: true\n");
        }
      }
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: advection tecnica aplicada no termo advectivo          
 *===================================================================*/
    else if((!strcmp(word,macro[30])))
    {
      initSec(word, OUTPUT_FOR_FILE);
/*...*/
      readAdvectionScheme(fileIn, &sc);
 /*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: edp equacoes diferencias resolvidas                    
 *===================================================================*/
    else if((!strcmp(word,macro[31])))
    {
      initSec(word, OUTPUT_FOR_FILE);
      readEdo(mesh0,fileIn);
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: diffusion tecnica aplicada no termo difusivo           
 *===================================================================*/
    else if((!strcmp(word,macro[32])))
    {
      initSec(word, OUTPUT_FOR_FILE);
/*... */
      readDiffusionScheme(fileIn, &sc);
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: pFluid : escreve os arquivos dos resultados do escoamente
   no regime permanente               
 *===================================================================*/
    else if( (!strcmp(word,macro[33])) && 
             (opt.stepPlot[1]++) == opt.stepPlot[0])
    {      
/*...*/
      opt.stepPlot[1] = 1;
      initSec(word, OUTPUT_FOR_SCREEN);
/*...................................................................*/

/*...*/                     
      printCall(&m        , &propVarFluid
              , &turbModel, &eModel 
              , &combModel, &tInterpol
              , pMesh     , &sc
              , loadsVel  , loadsPres 
              , loadsTemp , loadsZcomb
              , &opt       
              , mesh0     , mesh 
              , &media    , PLAST_TIME  
              , preName   , nameOut); 

/*...................................................................*/
      endSec(OUTPUT_FOR_SCREEN);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: setPrint                                   
 *===================================================================*/
    else if((!strcmp(word,macro[34])))
    {
      initSec(word, OUTPUT_FOR_FILE);
/*...*/
      setPrint(&opt,fileIn);
      initTimeStruct(&m        ,&tInterpol
                    ,mesh0     ,mesh 
                    ,&combModel,&opt);
      updateTimeStruct(&m       ,&tInterpol
                     ,mesh
                     ,&combModel,&opt); 
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
* macro: muda a escala da malha
*===================================================================*/
    else if ((!strcmp(word, macro[35])))
    {
      initSec(word, OUTPUT_FOR_FILE);
/*...*/
      reScaleMesh(mesh0,fileIn);
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }
/*===================================================================*/

/*===================================================================*
* macro: setPrime
*===================================================================*/
    else if ((!strcmp(word, macro[36])))
    {
      initSec(word, OUTPUT_FOR_FILE);
/*...*/
      readSetPrime(&m    , fileIn
                  , mesh0, mesh
                  , &prime, &fSolvPrime);
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }
/*===================================================================*/

/*===================================================================*
* macro: prime
*===================================================================*/
    else if ((!strcmp(word, macro[37])))
    {
      initSec(word, OUTPUT_FOR_FILE);
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
      endSec(OUTPUT_FOR_FILE);
    }
/*===================================================================*/

/*===================================================================*
* macro:                  
*===================================================================*/
    else if ((!strcmp(word, macro[38])))
    {
      initSec(word, OUTPUT_FOR_FILE);
      endSec(OUTPUT_FOR_FILE);
    }
/*===================================================================*/

/*===================================================================*
* setSolvComb : configuracao da combustao
*===================================================================*/
    else if ((!strcmp(word, macro[39])))
    {
      initSec(word, OUTPUT_FOR_FILE);
/*...*/
      readSolvFluid(&m         , mesh         , &reordMesh
                 , &combModel
                 , &solvVel   , &sistEqVel   , &fSolvVel
                 , &solvPres  , &sistEqPres  , &fSolvPres
                 , &solvEnergy, &sistEqEnergy, &fSolvEnergy
                 , &solvKturb , &sistEqKturb , &fSolvKturb
                 , &solvComb  , &sistEqComb  , &fSolvComb  
                 , pMesh       
                 , auxName    , preName      , nameOut
                 , fileIn                    , &opt);
/*...................................................................*/

/*...*/
      initEntalpyOfFormation(&combModel,&propVarFluid.sHeat);
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }
/*===================================================================*/

/*===================================================================*
 * macro: pCombustion : escreve os arquivos dos resultados do escoamente
 * reativos (combustao)              
 *===================================================================*/
    else if( (!strcmp(word,macro[40])) && 
             (opt.stepPlot[1]++) == opt.stepPlot[0])
    {      
/*...*/
      opt.stepPlot[1] = 1;
      initSec(word, OUTPUT_FOR_SCREEN);
/*...................................................................*/

/*...*/                     
      printCombustion(&m      , &propVarFluid
                  , &turbModel
                  , &eModel   , &combModel
                  , pMesh     , &sc
                  , loadsVel  , loadsPres 
                  , loadsTemp , loadsZcomb
                  , &opt       
                  , mesh0     , mesh 
                  , &media      
                  , preName   , nameOut);
/*...................................................................*/
      endSec(OUTPUT_FOR_SCREEN);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: simpleComb : escreve os arquivos dos resultados do escoamente
 * reativos (combustao)              
 *===================================================================*/
    else if ((!strcmp(word, macro[41])))
    {      
      initSec(word, OUTPUT_FOR_SCREEN);
      mpiWait();
      tm.solvEdpFluid = getTimeC() - tm.solvEdpFluid;
/*...*/
      if(!fSolvVel || !fSolvPres || !fSolvEnergy || !fSolvComb)
      {
        printf("Estrutara de dados nao montada para o setSolvComb!!!\n");
        exit(EXIT_FAILURE);
      }
/*...................................................................*/

/*...*/
     if(!fSolvSimple) 
     {
        printf("Simple nao configurado ainda!!!\n");
        exit(EXIT_FAILURE);
     }  
/*...................................................................*/

/*...*/
     if(!combModel.fCombustion) 
     {
        printf("Modelo de combustao nao configurado ainda!!!\n");
        exit(EXIT_FAILURE);
     }  
/*...................................................................*/

/*...*/
     combustionSolver(&m             , &propVarFluid
                    , loadsVel       , loadsPres
                    , loadsEnergy    , loadsKturb
                    , loadsZcomb
                    , &eModel        , &combModel
                    , &eMass         , &momentumModel
                    , &turbModel     , &thDynamic
                    , mesh0          , mesh
                    , &sistEqVel      , &sistEqPres
                    , &sistEqEnergy   , &sistEqKturb
                    , &sistEqComb
                    , &solvVel       , &solvPres
                    , &solvEnergy    , &solvKturb
                    , &solvComb
                    , &simple        , &sc
                    , pMesh          , &media
                    , &opt           , preName
                    , nameOut        , fileOut);
/*...................................................................*/

/*...*/
      tm.solvEdpFluid    = getTimeC() - tm.solvEdpFluid;
/*...................................................................*/
      endSec(OUTPUT_FOR_SCREEN);
   }   
/*===================================================================*/

/*===================================================================*
 *             
 *===================================================================*/
    else if((!strcmp(word,macro[42])))
    {
/*...*/
      initSec(word, OUTPUT_FOR_SCREEN);
      endSec(OUTPUT_FOR_SCREEN);
/*...................................................................*/
    }   
/*===================================================================*/

/*===================================================================*
 * macro:                       
 *===================================================================*/
    else if((!strcmp(word,macro[43])))
    {
      initSec(word, OUTPUT_FOR_FILE);
/*...*/
      readChemical(&combModel,fileIn);
/*...................................................................*/

      endSec(OUTPUT_FOR_FILE);
    }
/*===================================================================*/

/*===================================================================*
 * macro: residual : leitura dos residuos                      
 *===================================================================*/
    else if((!strcmp(word,macro[44])))
    {
      initSec(word, OUTPUT_FOR_FILE);
/*...*/
      readResidual(&simple, fileIn);
/*...................................................................*/

      endSec(OUTPUT_FOR_FILE);
    }
/*===================================================================*/


/*===================================================================*
 * macro: gravity : gravidade                                                 
 *===================================================================*/
    else if((!strcmp(word,macro[45])))
    {
      initSec(word, OUTPUT_FOR_FILE);
      readGravity(gravity,fileIn);
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: model: modelos utilizados nas equacao diferencias         
 *===================================================================*/
    else if((!strcmp(word,macro[46])))
    {
      initSec(word, OUTPUT_FOR_FILE);
      readModel(&eModel   ,&turbModel
              , &eMass    ,&momentumModel
              , diffModel ,transModel
              , &combModel  
              , fileIn);
//    fReadModels = true;
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: mean : calcula a media                                                        
 *===================================================================*/
    else if(!strcmp(word,macro[47]))
    {
      initSec(word, OUTPUT_FOR_SCREEN);
      calMean(&media,mesh,sc.ddt.dt[TIME_N_MINUS_1]
            , sc.ddt.t,sc.ddt.timeStep);
/*...................................................................*/
      endSec(OUTPUT_FOR_SCREEN);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: setMean : configuracao das medias temporais                                   
 *===================================================================*/
    else if(!strcmp(word,macro[48])){
      initSec(word, OUTPUT_FOR_FILE);
      readMean(&m,fileIn,mesh,&media);
/*...................................................................*/
      endSec(OUTPUT_FOR_FILE);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: save : salva o estado do programa                                          
 *===================================================================*/
    else if( (!strcmp(word,macro[49])) &&
             (save.step[1]++ == save.step[0]) ){
      initSec(word, OUTPUT_FOR_SCREEN);
      save.step[1] = 1;
      wSave(&propVarFluid,&turbModel
          ,&thDynamic
          ,mesh         ,&sc.ddt
          ,&media
          ,preName      ,nameOut);
/*...................................................................*/
      endSec(OUTPUT_FOR_SCREEN);
    }   
/*===================================================================*/

/*===================================================================*
 * macro: load : carrega um estado do sistema especifico                             
 *===================================================================*/
    else if((!strcmp(word,macro[50]))){
      initSec(word, OUTPUT_FOR_SCREEN);
      save.fLoad = true;
      load(&propVarFluid,&turbModel
          ,&thDynamic
          ,mesh         ,&sc.ddt
          ,&media
          ,fileIn);
/*...................................................................*/
      endSec(OUTPUT_FOR_SCREEN);
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
//      lG[i]  = MAT2D(nel, i, geom->mksi, maxViz);
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
//      lG[i] =   MAT2D(nel, i, geom->fArea, maxViz);
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
//      lG[i] = MAT2D(nel, i, geom->mvSkew, maxViz);
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
//        lG1[i][j] = MAT3D(nel, i, j, geom->ksi, maxViz, ndm);
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
//        lG1[i][j] = MAT3D(nel, i, j, geom->eta, maxViz, ndm);
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
//        lG1[i][j] = MAT3D(nel, i, j, geom->normal, maxViz, ndm);
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
      for (j = 0; j < ndm; j++) {
//        lG1[i][j] = MAT3D(nel, i, j, geom->xm, maxViz, ndm);
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
//        lG1[i][j] = MAT3D(nel, i, j, geom->vSkew, maxViz, ndm);
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
 
