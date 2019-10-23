#include<ReadFile.h>

/*...funcao de apoio*/
  static void getword(char *line, char*word);
  static int getnumprop2(char *line);
  static void setLoadsEnergy(PropVarFluid *pF
                           ,Loads *loadsEnergy   ,Loads *loadsTemp
                           ,Loads *loadsVel      ,DOUBLE *RESTRICT prop
                           ,bool const fTemp     ,bool const iKelvin);
 static void setLoadsEnergyMix(Combustion *cModel   ,Prop *pDen
                              ,Prop *sHeatProp   ,Loads *loadsEnergy
                              ,Loads *loadsTemp     ,Loads *loadsZ  
                              ,Loads *loadsVel      ,DOUBLE *RESTRICT prop
                              ,bool const fTemp     ,bool const fSheat 
                              ,bool const iKelvin   ,bool const fDensity
                              ,bool const fGrouped);  
 static void convLoadsZcombMix(Combustion *cModel  ,Prop *pDen
                             ,Prop *sHeatProp   ,Loads *loadsTemp 
                             ,Loads *loadsZ        ,Loads *loadsVel
                             ,DOUBLE *RESTRICT prop
                             ,bool const fTemp     ,bool const fSheat 
                             ,bool const iKelvin   ,bool const fDensity
                             ,bool const fGrouped);  

 static void convLoadsVelMix(Combustion *cModel   ,Prop *pDen
                            ,Prop *sHeatProp   ,Loads *loadsVel   
                            ,Loads *loadsTemp     ,Loads *loadsZ     
                            ,DOUBLE *RESTRICT prop
                            ,bool const fTemp     ,bool const fSheat 
                            ,bool const iKelvin   ,bool const fDensity
                            ,bool const fGrouped); 
 static void setLoadsVel(PropVarFluid *propFluid ,Loads *loadsVel
                        , Loads *loadsTemp       , DOUBLE *RESTRICT prop
                        , bool const fTemp       , bool const iKelvin);

 static void setLoadsPresC(Loads *loadsPres,Loads *loadsPresC
                           ,short const iCod);

 static void setLoadsRho(PropVarFluid *pf       ,Loads *loadsRho      
                          ,Loads *loadsVel     ,Loads *loadsTemp
                          ,bool const iKelvin );
 static void initFaceRrho(Loads *ldVel
                , short  *RESTRICT faceRvel, short  *RESTRICT faceRrho
                , short  *RESTRICT nFace   , short const maxViz
                , INT numel);
/*..................................................................*/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 16/20/2019                                   *
 *-------------------------------------------------------------------*
 * readFileFvMesh : leitura de arquivo de dados em volume finitos    *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * m     -> memoria principal                                        * 
 * mesh  -> malha                                                    *
 * propF -> propriedades variavei do fluido                          *
 * propD -> propriedades variavei do problema de difusao             *
 * propT -> propriedades variavei do problema de transporte          *
 * eModel-> modelo da eq energia                                     *
 * tModel-> modelo de turbulencia                                    *
 * cModel-> modelo de combustao                                      *
 * mModel-> modelo da equacao de momento                             *
 * file  -> arquivo de entrada                                       * 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * *******************************************************************/
void readFileFvMesh( Memoria *m              , Mesh *mesh
                   , PropVarFluid *propF           
                   , PropVarCD *propD        , PropVarCD *propT
                   , EnergyModel *energyModel, Turbulence *tModel     
                   , Combustion *cModel      , MomentumModel *mModel     
                   , Mean *media             , FILE* file)
{
  char word[WORD_SIZE],str[WORD_SIZE];
  char macro[NMACROS][WORD_SIZE]={
        "coordinates"  ,"endMesh"    ,"insert"         /* 0, 1, 2*/
       ,"return"       ,"cells"      ,"faceResT1"      /* 3, 4, 5*/
       ,""             ,"loadsT1"    ,"uniformT1"      /* 6, 7, 8*/ 
       ,"faceResZ"     ,"faceLoadZ"  ,"loadsZ"         /* 9,10,11*/ 
       ,"faceResD1"    ,"uniformD1"  ,"loadsD1"        /*12,13,14*/ 
       ,""             ,"initialD1"  ,""               /*15,16,17*/ 
       ,"faceResVel"   ,"loadsVel"   ,""               /*18,19,20*/ 
       ,"faceResPres"  ,"loadsPres"  ,"faceRpres"      /*21,22,23*/
       ,"faceResTemp"  ,"loadsTemp"  ,""               /*24,25,26*/
       ,"materials"    ,"uniformPres","initialVel"     /*27,28,29*/
       ,"uniformTemp"  ,"uniformVel" ,"uniformZ"       /*30,31,32*/
       ,"faceReKturb"  ,"loadsKturb" ,"faceLoadKturb"  /*33,34,35*/
       ,"fileMaterials","initialTemp",""               /*36,37,38*/
	   };                                             
  bool rflag[NMACROS],macroFlag;
  bool fOneEqK,fComb,fWallModel,fTurbStruct,fDynamic;
  INT nn,nel;
  short maxno,ndm,numat,maxViz,ndfVel,nComb,nSpPri,nReac;
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
  ndfVel = max(mesh->ndfF - 1,mesh->ndfFt - 2);
/*...................................................................*/

/*...*/   
  fComb       = cModel->fCombustion;
  nComb       = cModel->nComb;
  nSpPri      = cModel->nOfSpecies;
  nReac       = cModel->chem.nReac;
  fOneEqK     = tModel->fOneEq;
  fWallModel  = tModel->fWall;
  fTurbStruct = tModel->fTurbStruct;
  fDynamic    = tModel->fDynamic;
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

/*...*/
  if(mesh->ndm == 3)
    mesh->ntn = 6;
  else
    mesh->ntn = 4;
/*...................................................................*/

/*... alocando variavies de elementos*/
/*... conectividade*/ 
  HccaAlloc(INT,m,mesh->elm.node       ,nel*maxno,"elnode"  ,_AD_);
/*... materiais*/ 
  HccaAlloc(short,m,mesh->elm.mat      ,nel       ,"elmat"   ,_AD_);
/*... nos por elementos*/
  HccaAlloc(short,m,mesh->elm.nen      ,nel      ,"elnen"   ,_AD_);
/*... tipo geometrico */
  HccaAlloc(short,m,mesh->elm.geomType ,nel      ,"elgT"    ,_AD_);
/*... centroide */
  HccaAlloc(DOUBLE,m,mesh->elm.geom.cc ,nel*ndm,"elCc"    ,_AD_);
/*... face por elemento*/
  HccaAlloc(INT, m, mesh->elm.cellFace, nel*maxViz,"cellface", _AD_);

  if( mpiVar.nPrcs < 2){

/*... volume da celula*/                           
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.volume
            ,nel            ,"elVol",_AD_);
/*... vetor que une o centroide ao ponto medio*/                           
    HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.xmcc
             ,nel*maxViz*ndm       ,"elxmcc",_AD_);
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
    zero(mesh->elm.geom.volume  ,nel           ,DOUBLEC);
    zero(mesh->elm.geom.xmcc    ,nel*ndm*maxViz,DOUBLEC);
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
  if(mesh->ndfT[0] > 0 || mesh->ndfF > 0 || mesh->ndfFt > 0) {
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
  if(mesh->ndfF > 0 || mesh->ndfFt > 0) 
  {
/*... alocando memoria*/
/*... cc da equacao de velocidades*/
     HccaAlloc(short,m,mesh->elm.faceRvel  
            ,nel*(maxViz+1),"faceRvel"    ,_AD_);
     zero(mesh->elm.faceRvel  ,nel*(maxViz+1),"short"  );
     
     HccaAlloc(short,m,mesh->elm.faceRvel  
            ,nel*(maxViz+1),"faceLVel"    ,_AD_);
     zero(mesh->elm.faceRvel  ,nel*(maxViz+1),"short"  );

/*... cc da equacao de pressao*/
     HccaAlloc(short,m,mesh->elm.faceRpres 
            ,nel*(maxViz+1),"faceRpres"   ,_AD_);
     zero(mesh->elm.faceRpres ,nel*(maxViz+1),"short"  );
     
     HccaAlloc(short,m,mesh->elm.faceRpres 
            ,nel*(maxViz+1),"faceLPres"   ,_AD_);
     zero(mesh->elm.faceRpres ,nel*(maxViz+1),"short"  );

/*... viscosidade turbulenta*/
     HccaAlloc(DOUBLE, m, mesh->elm.eddyViscosity
              , nel  , "eddyVis", _AD_);
     zero(mesh->elm.eddyViscosity, nel, DOUBLEC);

/*... tensor residual turbulenta*/
     if(fTurbStruct)
     {
       HccaAlloc(DOUBLE, m, mesh->elm.stressR        
              , nel*mesh->ntn , "stressR", _AD_);
       zero(mesh->elm.stressR, nel*mesh->ntn, DOUBLEC);
     } 

/*... paramentros de parede*/
     if(fWallModel)
     {
        HccaAlloc(DOUBLE, m, mesh->elm.wallParameters
                 , nel*NWALLPAR   , "wallParm", _AD_);
        zero(mesh->elm.wallParameters, nel*NWALLPAR, DOUBLEC);
     }

/*... coeficientes dinamicos locais*/
     if(fDynamic)
     {
       HccaAlloc(DOUBLE, m, mesh->elm.cd
              , nel*2   , "cDynamic"     , _AD_);
       zero(mesh->elm.cd, nel*2, DOUBLEC);
     } 

/*... energia cinetica turbulenta*/
     if(fOneEqK){
       HccaAlloc(DOUBLE, m           , mesh->elm.kTurb
              , nel  , "kTurb"     , _AD_);
       zero(mesh->elm.kTurb, nel, DOUBLEC);
       
       HccaAlloc(DOUBLE, m        , mesh->elm.kTurb0
                , nel  , "kTurb0" , _AD_);
       zero(mesh->elm.kTurb0, nel, DOUBLEC);
       
       HccaAlloc(DOUBLE  , m           , mesh->elm.gradKturb
                , nel*ndm, "gradKturb" , _AD_);
       zero(mesh->elm.gradKturb, nel*ndm, DOUBLEC);
       
       HccaAlloc(DOUBLE, m        , mesh->node.kTurb
                , nn   , "nkTurb" , _AD_);
       zero(mesh->node.kTurb, nn  , DOUBLEC);

/*... eGradTurb*/
       HccaAlloc(DOUBLE,m,mesh->elm.gradKturb
                ,nel*ndm        ,"eGradKturb",_AD_);
       zero(mesh->elm.gradKturb ,nel*ndm,DOUBLEC);
       
/*... cc da equacao energia cinetica turbulenta*/
       HccaAlloc(short,m,mesh->elm.faceReKturb
              ,nel*(maxViz+1),"faceReKturb",_AD_);
       zero(mesh->elm.faceReKturb ,nel*(maxViz+1),"short"  );
       
       HccaAlloc(short,m,mesh->elm.faceLoadKturb 
              ,nel*(maxViz+1),"faceLKturb",_AD_);
       zero(mesh->elm.faceLoadKturb ,nel*(maxViz+1),"short"  );
     }
 
/*... densityFluid*/
     HccaAlloc(DOUBLE , m         , mesh->elm.densityFluid
              ,nel * DENSITY_LEVEL, "eDenFluid", _AD_);
     zero(mesh->elm.densityFluid, nel * DENSITY_LEVEL, DOUBLEC);

/*... ePres(n+1)*/
     HccaAlloc(DOUBLE,m,mesh->elm.pressure 
            ,nel              ,"pressure"          ,_AD_);
     zero(mesh->elm.pressure  ,nel ,DOUBLEC);

/*... ePres(n)*/
     HccaAlloc(DOUBLE,m,mesh->elm.pressure0
            ,nel              ,"pressure0"         ,_AD_);
     zero(mesh->elm.pressure0 ,nel,DOUBLEC);

/*... pres*/
     HccaAlloc(DOUBLE,m,mesh->node.pressure
              ,nn    ,"npressure",_AD_);
     zero(mesh->node.pressure     ,nn,DOUBLEC);

/*... nGradVel*/
     HccaAlloc(DOUBLE,m,mesh->node.gradVel  
              ,nn*ndm*ndfVel ,"nGradVel"     ,_AD_);
     zero(mesh->node.gradVel  ,nn*ndm*ndfVel,DOUBLEC);
     
/*... eGradVel*/
     HccaAlloc(DOUBLE,m,mesh->elm.gradVel 
              ,nel*ndm*ndfVel ,"eGradVel"     ,_AD_);
     zero(mesh->elm.gradVel   ,nel*ndm*ndfVel,DOUBLEC);

/*... nGradPres*/
     HccaAlloc(DOUBLE,m,mesh->node.gradPres 
              ,nn*ndm        ,"nGradPres"    ,_AD_);
     zero(mesh->node.gradPres ,nn*ndm,DOUBLEC);
     
/*... eGradPres*/
     HccaAlloc(DOUBLE,m,mesh->elm.gradPres
              ,nel*ndm        ,"eGradPres"     ,_AD_);
     zero(mesh->elm.gradPres  ,nel*ndm,DOUBLEC);

/*... */
     if (mesh->ndfFt > 0) 
     {
/*...*/
       HccaAlloc(short, m, mesh->elm.faceRenergy
                , nel*(maxViz + 1), "faceRenergy", _AD_);
       zero(mesh->elm.faceRenergy, nel*(maxViz + 1), "short");
/*...*/
       HccaAlloc(short, m, mesh->elm.faceRrho    
                , nel*(maxViz + 1), "faceRrho", _AD_);
       zero(mesh->elm.faceRrho, nel*(maxViz + 1), "short");

/*... eEnergy0*/
       HccaAlloc(DOUBLE,m         ,mesh->elm.energy0
                ,nel   ,"eEnergy0",_AD_);
       zero(mesh->elm.energy0,nel,DOUBLEC);

/*... eEnergy*/
       HccaAlloc(DOUBLE,m        ,mesh->elm.energy
                ,nel   ,"eEnergy",_AD_);
       zero(mesh->elm.energy, nel, DOUBLEC);

/*... nEnergy*/
       HccaAlloc(DOUBLE,m        ,mesh->node.energy
                ,nn    ,"nEnergy",_AD_);
       zero(mesh->node.energy, nn, DOUBLEC);

/*... nGradEnergy*/
       HccaAlloc(DOUBLE, m, mesh->node.gradEnergy
                ,nn*ndm, "nGradEnergy", _AD_);
       zero(mesh->node.gradEnergy, nn*ndm, DOUBLEC);

/*... eGradEnergy*/
       HccaAlloc(DOUBLE ,m            ,mesh->elm.gradEnergy
                ,nel*ndm,"eGradEnergy",_AD_);
       zero(mesh->elm.gradEnergy, nel*ndm, DOUBLEC);

/*... eTemp*/
       HccaAlloc(DOUBLE,m        ,mesh->elm.temp
                ,nel   ,"eTemp",_AD_);
       zero(mesh->elm.temp, nel, DOUBLEC);

/*... eTemp*/
       HccaAlloc(DOUBLE,m        ,mesh->elm.temp0
                ,nel   ,"eTemp0",_AD_);
       zero(mesh->elm.temp0, nel, DOUBLEC);

/*... nTemp*/
       HccaAlloc(DOUBLE,m        ,mesh->node.temp
                ,nn    ,"nTemp",_AD_);
       zero(mesh->node.energy, nn, DOUBLEC);

/*... nGradTemp*/
       HccaAlloc(DOUBLE, m, mesh->node.gradTemp
                ,nn*ndm, "nGradTemp", _AD_);
       zero(mesh->node.gradTemp, nn*ndm, DOUBLEC);

/*... eGradTemp*/
       HccaAlloc(DOUBLE ,m            ,mesh->elm.gradTemp
                ,nel*ndm,"eGradTemp",_AD_);
       zero(mesh->elm.gradTemp, nel*ndm, DOUBLEC);

/*... calor especifico*/
        HccaAlloc(DOUBLE, m, mesh->elm.specificHeat
                 , nel * SHEAT_LEVEL, "sHeat", _AD_);
        zero(mesh->elm.specificHeat, nel * SHEAT_LEVEL, DOUBLEC);

/*... viscosidade dinamica*/
        HccaAlloc(DOUBLE, m, mesh->elm.dViscosity
                 , nel  , "dVis", _AD_);
        zero(mesh->elm.dViscosity, nel, DOUBLEC);
        
/*... condutividade termica*/
        HccaAlloc(DOUBLE, m, mesh->elm.tConductivity
                 , nel  , "tCon", _AD_);
        zero(mesh->elm.tConductivity, nel, DOUBLEC);

/*... nRhoFluid*/
       HccaAlloc(DOUBLE, m, mesh->node.rhoFluid 
                ,nn*DENSITY_LEVEL, "nRhoFluid", _AD_);
       zero(mesh->node.rhoFluid, nn*DENSITY_LEVEL, DOUBLEC);

/*... eGradRho*/
       HccaAlloc(DOUBLE ,m            ,mesh->elm.gradRhoFluid
                ,nel*ndm,"eGradFluid",_AD_);
       zero(mesh->elm.gradRhoFluid, nel*ndm, DOUBLEC);

     }
/*...................................................................*/

/*...*/
     if( mpiVar.nPrcs < 2)
     {
/*... rCellVel*/
       HccaAlloc(DOUBLE,m,mesh->elm.rCellVel  
                ,nel*ndm*ndfVel       ,"rCellVel"     ,_AD_);
       zero(mesh->elm.rCellVel  ,nel*ndm*ndfVel   ,DOUBLEC);
/*... rCellPres*/
       HccaAlloc(DOUBLE,m,mesh->elm.rCellPres 
                ,nel                 ,"rCellPres"    ,_AD_);
       zero(mesh->elm.rCellPres ,nel              ,DOUBLEC);
/*... rCellEnergy*/
       if (mesh->ndfFt > 0)
       {
          HccaAlloc(DOUBLE,m            ,mesh->elm.rCellEnergy
                   ,nel   ,"rCellEnergy",_AD_);
          zero(mesh->elm.rCellEnergy, nel, DOUBLEC);
       }
/*... rCellKturb*/
       if(fOneEqK)
       {
         HccaAlloc(DOUBLE,m,mesh->elm.rCellKturb
                ,nel  ,"rCellKturb",_AD_);
         zero(mesh->elm.rCellKturb,nel,DOUBLEC);
       }
     }
/*...................................................................*/
  }
/*...................................................................*/

/*... modelo de combustao*/
  if (fComb)
  {
    HccaAlloc(short, m, mesh->elm.faceResZcomb
            , nel*(maxViz + 1), "faceResComb", _AD_);
    zero(mesh->elm.faceResZcomb, nel*(maxViz + 1), "short");

    HccaAlloc(short, m, mesh->elm.faceLoadZcomb
      , nel*(maxViz + 1), "faceLoadComb", _AD_);
    zero(mesh->elm.faceLoadZcomb, nel*(maxViz + 1), "short");

/*... eZcomb*/
    HccaAlloc(DOUBLE, m, mesh->elm.zComb
      , nel*nComb, "zComb", _AD_);
    zero(mesh->elm.zComb, nel*nComb, DOUBLEC);

/*... eZcomb0*/
    HccaAlloc(DOUBLE, m, mesh->elm.zComb0
      , nel*nComb, "zComb0", _AD_);
    zero(mesh->elm.zComb0, nel*nComb, DOUBLEC);

/*... eGradZcomb*/
    HccaAlloc(DOUBLE, m, mesh->elm.gradZcomb
      , nel*ndm*nComb, "gradComb", _AD_);
    zero(mesh->elm.gradZcomb, nel*ndm*nComb, DOUBLEC);

/*... nZcomb*/
    HccaAlloc(DOUBLE  ,m       ,mesh->node.zComb
             ,nn*nComb,"nZcomb",_AD_);
    zero(mesh->node.zComb, nn*nComb, DOUBLEC);

/*... nGradZcomb*/
    HccaAlloc(DOUBLE      , m, mesh->node.gradZcomb
             ,nn*ndm*nComb, "nGradZcomb", _AD_);
    zero(mesh->node.gradZcomb, nn*ndm*nComb, DOUBLEC);

/*... */
    HccaAlloc(DOUBLE, m, mesh->elm.cDiffComb
      , nel*nSpPri, "cDiffComb", _AD_);
    zero(mesh->elm.cDiffComb, nel*nSpPri, DOUBLEC);

/*... */
    HccaAlloc(DOUBLE, m, mesh->elm.wk          
            , nel* nSpPri,"wk", _AD_);
    zero(mesh->elm.wk, nel* nSpPri, DOUBLEC);

/*... */
    HccaAlloc(DOUBLE, m, mesh->elm.rateHeatReComb
            , nel, "rateHeatCom", _AD_);
    zero(mesh->elm.rateHeatReComb, nel, DOUBLEC);

/*... yFrac*/
    HccaAlloc(DOUBLE, m, mesh->elm.yFrac
            , nel*nSpPri, "yFrac", _AD_);
    zero(mesh->elm.yFrac, nel*nSpPri, DOUBLEC);
/*... yFrac0*/
    HccaAlloc(DOUBLE, m, mesh->elm.yFrac0
            , nel*nSpPri, "yFrac0", _AD_);
    zero(mesh->elm.yFrac0, nel*nSpPri, DOUBLEC);

/*... enthalpyK*/
    HccaAlloc(DOUBLE, m, mesh->elm.enthalpyk
            , nel*nSpPri, "enthalpyk", _AD_);
    zero(mesh->elm.enthalpyk, nel*nSpPri, DOUBLEC);

/*... eGradY*/
    HccaAlloc(DOUBLE, m, mesh->elm.gradY
      , nel*ndm*nSpPri, "eGradY", _AD_);
    zero(mesh->elm.gradY, nel*ndm*nSpPri, DOUBLEC);

/*... timeReactor*/
    HccaAlloc(DOUBLE, m, mesh->elm.tReactor
      , nel*N_TERMS_REACTOR, "tReactor", _AD_);
    zero(mesh->elm.tReactor, nel*N_TERMS_REACTOR, DOUBLEC);

    if (mpiVar.nPrcs < 2)
    {
      if (fComb)
      {
        HccaAlloc(DOUBLE, m, mesh->elm.rCellComb
          , nel*nComb, "rCellComb", _AD_);
        zero(mesh->elm.rCellComb, nel*nComb, DOUBLEC);
      }
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
     HccaAlloc(DOUBLE, m, mesh->elm.densityUt1
             , nel*DENSITY_LEVEL, "densityUt1", _AD_);
     zero(mesh->elm.densityUt1, nel*DENSITY_LEVEL, DOUBLEC);

/*... ceoficiente de diffusividae T1*/
     HccaAlloc(DOUBLE, m, mesh->elm.cDiffT1, nel, "cDiffT1", _AD_);
     zero(mesh->elm.cDiffT1, nel, DOUBLEC);

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
              ,nel*DENSITY_LEVEL ,"densityUd1" ,_AD_);
     zero(mesh->elm.densityUd1  ,nel*DENSITY_LEVEL,DOUBLEC);

/*... ceoficiente de diffusividae D1*/
     HccaAlloc(DOUBLE, m, mesh->elm.cDiffD1, nel, "cDiffD1", _AD_);
     zero(mesh->elm.cDiffD1, nel, DOUBLEC);

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
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[0] = true;
      fprintf(fileLogExc,"loading coordinates...\n");
      readVfCoor(mesh->node.x,mesh->nnode,mesh->ndm,file);      
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... endMesh*/
    else if((!strcmp(word,macro[1])) && (!rflag[1])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      fprintf(fileLogExc,"%s\n\n",DIF);
      strcpy(macros[nmacro++],word);
      rflag[1] = true;
      macroFlag = false;  
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... insert*/
    else if((!strcmp(word,macro[2])) && (!rflag[2])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      fprintf(fileLogExc,"%s\n\n",DIF);
      fscanf(file,"%s",nameAux);
      fileAux = file;
      file    = openFile(nameAux,"r");
      rflag[2] = true;
    }
/*...................................................................*/

/*... return*/
    else if((!strcmp(word,macro[3]))){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      fprintf(fileLogExc,"%s\n\n",DIF);
      if(!rflag[2]){
        ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__
                  ,"Erro: macro return sem um insert associado!!"
                  ,EXIT_PROG);
      }
      fclose(file);
      file = fileAux;
      rflag[2] = false;
    }
/*...................................................................*/

/*... cells  */
    else if((!strcmp(word,macro[4])) && (!rflag[4])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[4] = true;
      fprintf(fileLogExc,"loading cells ...\n");
      readVfElmt(mesh->elm.node    ,mesh->elm.mat
                ,mesh->elm.nen     ,mesh->elm.adj.nViz
                ,mesh->elm.geomType,mesh->numel
                ,mesh->maxNo       ,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... faceResT1 */
    else if((!strcmp(word,macro[5])) && (!rflag[5])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[5] = true;
      strcpy(str,"endFaceResT1");
      fprintf(fileLogExc,"loading faceResT1 ...\n");
      readVfRes(mesh->elm.faceRt1,mesh->numel,mesh->maxViz+1,str,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... loadT1 - definicao de cargar transporte */
    else if((!strcmp(word,macro[7])) && (!rflag[7])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[7] = true;
      strcpy(str,"endLoadsT1");
      fprintf(fileLogExc,"loading loadsT1 ...\n");
      readVfLoads(loadsT1,str       ,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... uniformT1 */
    else if ((!strcmp(word, macro[8])) && (!rflag[8])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++], word);
      rflag[13] = true;
      uniformField(mesh->elm.u0T1, mesh->numel, 1, file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/


/*... faceResZ */
    else if((!strcmp(word,macro[9])) && (!rflag[9])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[9] = true;
      strcpy(str,"endFaceResZ");
      fprintf(fileLogExc,"loading faceResZ ...\n");
      readVfRes(mesh->elm.faceResZcomb,mesh->numel,mesh->maxViz+1,str,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... faceLoadZ*/
    else if((!strcmp(word,macro[10])) && (!rflag[10])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[10] = true;
      strcpy(str,"endFaceLoadZ");
      fprintf(fileLogExc,"loading faceLoadZ ...\n");
      readVfRes(mesh->elm.faceLoadZcomb,mesh->numel
               ,mesh->maxViz+1       ,str       ,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... loadZ*/
    else if((!strcmp(word,macro[11])) && (!rflag[11])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[11] = true;
      strcpy(str,"endLoadsZ");
      fprintf(fileLogExc,"loading loadsZ ...\n");
      readVfLoads(loadsZcomb,str       ,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... faceResD1- condicao de contorno para problemas de difusa pura */
    else if((!strcmp(word,macro[12])) && (!rflag[12])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[12] = true;
      strcpy(str,"endFaceResD1");
      fprintf(fileLogExc,"loading faceResD1 ...\n");
      readVfRes(mesh->elm.faceRd1,mesh->numel,mesh->maxViz+1,str,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... uniformD1 */
    else if ((!strcmp(word, macro[13])) && (!rflag[13])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++], word);
      rflag[13] = true;
      uniformField(mesh->elm.u0D1, mesh->numel, 1, file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... loadD1 - definicao de cargar difusao pura */
    else if((!strcmp(word,macro[14])) && (!rflag[14])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[14] = true;
      strcpy(str,"endLoadsD1");
      fprintf(fileLogExc,"loading loadsD1 ...\n");
      readVfLoads(loadsD1,str       ,file);
//    loadsD1[0].nTypeVar = LFUNC;
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... initialD1 */
    else if ((!strcmp(word, macro[16])) && (!rflag[16])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++], word);
      rflag[29] = true;
      strcpy(str, "endInitialD1");
      fprintf(fileLogExc, "loading initialD1 ...\n");
      readVfInitial(mesh->elm.u0D1, mesh->numel, 1, str, file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... faceResvel- condicao de contorno para problemas fluidos (Vel) */
    else if((!strcmp(word,macro[18])) && (!rflag[18])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[18] = true;
      strcpy(str,"endFaceResVel");
      fprintf(fileLogExc,"loading faceRvel ...\n");
      readVfRes(mesh->elm.faceRvel,mesh->numel,mesh->maxViz+1,str,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... loadVel - definicao de cargar fluidos (Vel)*/
    else if((!strcmp(word,macro[19])) && (!rflag[19])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[19] = true;
      strcpy(str,"endLoadsVel");
      fprintf(fileLogExc,"loading loadsVel ...\n");
      readVfLoads(loadsVel,str       ,file);      
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... faceRespres - condicao de contorno para problemas fluidos (Pres)*/
    else if((!strcmp(word,macro[21])) && (!rflag[21])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[21] = true;
      strcpy(str,"endFaceResPres");
      fprintf(fileLogExc,"loading faceResPres ...\n");
      readVfRes(mesh->elm.faceRpres,mesh->numel
               ,mesh->maxViz+1     ,str        ,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... loadPres - definicao de cargar fluidos (Pres)*/
    else if((!strcmp(word,macro[22])) && (!rflag[22])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[22] = true;
      strcpy(str,"endLoadsPres");
      fprintf(fileLogExc,"loading loadsPres ...\n");
      readVfLoads(loadsPres,str       ,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... faceRenergy - condicao de contorno para problemas fluidos (Energy)*/
    else if ((!strcmp(word, macro[24])) && (!rflag[24])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++], word);
      rflag[24] = true;
      strcpy(str, "endFaceResTemp");
      fprintf(fileLogExc,"loading faceRtemp ...\n");
      readVfRes(mesh->elm.faceRenergy,mesh->numel
               ,mesh->maxViz + 1     ,str        ,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... loadEnergy - definicao de cargar fluidos (Energy)*/
    else if ((!strcmp(word, macro[25])) && (!rflag[25])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++], word);
      rflag[25] = true;
      strcpy(str, "endLoadsTemp");
      fprintf(fileLogExc,"loading loadsTemp ...\n");
      readVfLoads(loadsTemp, str, file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... materiais */
    else if((!strcmp(word,macro[27])) && (!rflag[27])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[27] = true;
      strcpy(str,"endMaterials");
      fprintf(fileLogExc,"loading materials ...\n");
      readVfMat(mesh->elm.material.prop,mesh->elm.material.type
               ,numat,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... uniformPres */
    else if ((!strcmp(word, macro[28])) && (!rflag[28])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++], word);
      rflag[28] = true;
      uniformField(mesh->elm.pressure0, mesh->numel,1, file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... initialVel */
    else if((!strcmp(word,macro[29])) && (!rflag[29])){
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++],word);
      rflag[29] = true;
      strcpy(str,"endInitialVel");
      fprintf(fileLogExc,"loading initialVel ...\n");
      readVfInitial(mesh->elm.vel0,mesh->numel,mesh->ndm,str,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... uniformTemp */
    else if ((!strcmp(word, macro[30])) && (!rflag[30])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++], word);
      rflag[30] = true;
      fprintf(fileLogExc,"loading uniformTemp ...\n");
      uniformField(mesh->elm.temp0, mesh->numel, 1, file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... uniformVel */
    else if ((!strcmp(word, macro[31])) && (!rflag[31])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++], word);
      rflag[31] = true;
      fprintf(fileLogExc,"loading uniformVel ...\n");
      uniformField(mesh->elm.vel0, mesh->numel, ndm, file);        
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... uniformZ */
    else if ((!strcmp(word, macro[32])) && (!rflag[32])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      if(!fComb)
      {
        printf("Combustion model disable!\n");
        exit(EXIT_FAILURE);
      }
      strcpy(macros[nmacro++], word);
      rflag[32] = true;
      fprintf(fileLogExc,"loading uniformZ ...\n");
      uniformField(mesh->elm.zComb0, mesh->numel, nComb, file);        
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... faceReKturb - condicao de contorno para problemas fluidos (Kturb)*/
    else if ((!strcmp(word, macro[33])) && (!rflag[33])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++], word);
      rflag[33] = true;
      strcpy(str, "endFaceReKturb");
      fprintf(fileLogExc,"loading faceReKturb ...\n");
      readVfRes(mesh->elm.faceReKturb,mesh->numel
               ,mesh->maxViz + 1     ,str        ,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... loadKturb - definicao de cargar fluidos (Kturb)*/
    else if ((!strcmp(word, macro[34])) && (!rflag[34])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++], word);
      rflag[34] = true;
      strcpy(str, "endLoadsKturb");
      fprintf(fileLogExc,"loading loadsKturb ...\n");
      readVfLoads(loadsKturb, str, file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... faceLoadKturb - cargas nas faces fluido (Kturb)*/
    else if ((!strcmp(word, macro[35])) && (!rflag[35])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF, word);
      strcpy(macros[nmacro++], word);
      rflag[35] = true;
      strcpy(str, "endFaceLoadKturb");
      fprintf(fileLogExc,"loading faceLoadKturb ...\n");
      readVfRes(mesh->elm.faceLoadKturb,mesh->numel
               ,mesh->maxViz + 1        ,str        ,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... fileMaterials - materias dados vias arquivos*/
    else if ((!strcmp(word, macro[36])) && (!rflag[36])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF,word);
      strcpy(macros[nmacro++], word);
      rflag[36] = true;
      strcpy(str, "endFileMaterials");
      fprintf(fileLogExc, "loading endFileMaterials ...\n");
      readFileMat(mesh->elm.material.prop, mesh->elm.material.type
                 , numat, file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/

/*... initTemp - valore de temperatura inicias*/
    else if ((!strcmp(word, macro[37])) && (!rflag[37])) {
      fprintf(fileLogExc, "%s\n%s\n", DIF,word);
      strcpy(macros[nmacro++], word);
      rflag[36] = true;
      strcpy(str, "endInitialTemp");
      fprintf(fileLogExc, "loading endInitTemp ...\n");
      readVfInitial(mesh->elm.temp0,mesh->numel,1,str,file);
      fprintf(fileLogExc, "done.\n%s\n\n", DIF);
    }
/*...................................................................*/



  }while(macroFlag && (!feof(file)));
  
/*... propriedades de referencia*/
  initPropRef(propF,mesh->elm.material.prop,0);
/*...................................................................*/

/*... incializando as condicoes de contorno da PresVel*/
  if(rflag[22])
    setLoadsPresC(loadsPres, loadsPresC,mModel->iCodBuoyant);
/*...................................................................*/

/*... incializando as condicoes de contorno da Vel*/
  if(rflag[19])
  {
    if(fComb)
      convLoadsVelMix(cModel                   ,&propF->den 
                    ,&propF->sHeat            ,loadsVel              
                    ,loadsTemp                ,loadsZcomb
                    ,mesh->elm.material.prop
                    ,energyModel->fTemperature,propF->fSpecificHeat
                    ,energyModel->fKelvin     ,propF->fDensity
                    ,cModel->fLump);   
    else
      setLoadsVel(propF                     ,loadsVel              
                 ,loadsTemp                 ,mesh->elm.material.prop
                 ,energyModel->fTemperature ,energyModel->fKelvin);   
  }
/*...................................................................*/

/*... incializando as condicoes de contorno da energia*/
  if(rflag[25])
  {
    if(fComb)
      setLoadsEnergyMix(cModel                 ,&propF->den 
                      ,&propF->sHeat            ,loadsEnergy              
                      ,loadsTemp                ,loadsZcomb
                      ,loadsVel                 ,mesh->elm.material.prop
                      ,energyModel->fTemperature,propF->fSpecificHeat
                      ,energyModel->fKelvin     ,propF->fDensity
                      ,cModel->fLump);  
    else
      setLoadsEnergy(propF
                     ,loadsEnergy              ,loadsTemp
                     ,loadsVel                 ,mesh->elm.material.prop
                     ,energyModel->fTemperature,energyModel->fKelvin);  
  }
/*...................................................................*/

/*... incializando as condicoes de contorno da Zcomb*/
  if(rflag[11])
  {
    convLoadsZcombMix(cModel               ,&propF->den 
                 ,&propF->sHeat            ,loadsTemp   
                 ,loadsZcomb               ,loadsVel    
                 ,mesh->elm.material.prop
                 ,energyModel->fTemperature,propF->fSpecificHeat
                 ,energyModel->fKelvin     ,propF->fDensity
                 ,cModel->fLump);  
  }
/*...................................................................*/

/*... incializando as condicoes de contorno o rho*/
  setLoadsRho(propF               , loadsRhoFluid
            , loadsVel            , loadsTemp
            , energyModel->fKelvin );
  if(rflag[18])
    initFaceRrho( loadsVel
                , mesh->elm.faceRvel       , mesh->elm.faceRrho  
                , mesh->elm.adj.nViz       , mesh->maxViz
                , mesh->numelNov);
/*...................................................................*/

/*... combustao*/
 if(fComb)
  {
/*...*/
    alphaProdVector(1.e0             , mesh->elm.zComb0
                 , mesh->numel*nComb , mesh->elm.zComb);
/*...................................................................*/

/*...*/
    getSpeciesPrimitives(cModel
                      ,mesh->elm.yFrac,mesh->elm.zComb
                      ,mesh->numelNov);
/*...................................................................*/  

/*...*/
    alphaProdVector(1.e0              , mesh->elm.yFrac
                 , mesh->numel*nSpPri , mesh->elm.yFrac0);
/*...................................................................*/          
  }
/*...................................................................*/

/*... iniciacao de propriedades do material varaivel com o tempo*/
  if(mesh->ndfD[0] > 0) 
  {
/*...*/
    alphaProdVector(1.e0, mesh->elm.u0D1
                  , mesh->numel*mesh->ndfD[0], mesh->elm.uD1);
    if(propD[0].fDensity)
/*... inicia a massa especifica com o campo inicial*/
      initPropCD(&propD[0].den          , mesh->elm.densityUd1
               , mesh->elm.u0D1         , mesh->elm.material.prop
               , mesh->elm.mat
               , DENSITY_LEVEL          , mesh->numel
               , DENSITY);
    else
      initProp(mesh->elm.densityUd1
            ,mesh->elm.material.prop,mesh->elm.mat
            ,DENSITY_LEVEL          ,mesh->numel
            ,DENSITY);    


/*... inicializando o coeficiente de difusao*/
    if (propD[0].fCeofDiff)
      initPropCD(&propD[0].ceofDiff, mesh->elm.cDiffD1
        , mesh->elm.u0D1, mesh->elm.material.prop
        , mesh->elm.mat
        , COEFDIFF_LEVEL, mesh->numel
        , COEFDIF);
    else
      initProp(mesh->elm.cDiffD1
             , mesh->elm.material.prop, mesh->elm.mat
             , COEFDIFF_LEVEL         , mesh->numel
             , COEFDIF);
/*...................................................................*/

  }
/*...................................................................*/

/*...*/
  if(mesh->ndfT[0] > 0) 
  {
/*...*/
    alphaProdVector(1.e0             , mesh->elm.vel0
                 , mesh->numel*ndm   , mesh->elm.vel);
/*...*/
    alphaProdVector(1.e0                     , mesh->elm.u0T1
                  , mesh->numel*mesh->ndfT[0], mesh->elm.uT1);
    if (propT[0].fDensity)
/*... inicia a massa especifica com o campo inicial*/
      initPropCD(&propT[0].den , mesh->elm.densityUt1
                , mesh->elm.u0T1, mesh->elm.material.prop
                , mesh->elm.mat
                , DENSITY_LEVEL , mesh->numel
                , DENSITY);
    else
      initProp(mesh->elm.densityUt1
             , mesh->elm.material.prop, mesh->elm.mat
             , DENSITY_LEVEL, mesh->numel
             , DENSITY);


/*... inicializando o coeficiente de difusao*/
    if (propT[0].fCeofDiff)
      initPropCD(&propT[0].ceofDiff, mesh->elm.cDiffT1
               , mesh->elm.u0T1    , mesh->elm.material.prop
               , mesh->elm.mat
               , COEFDIFF_LEVEL    , mesh->numel
               , COEFDIF);
  else
    initProp(mesh->elm.cDiffT1
           , mesh->elm.material.prop, mesh->elm.mat
           , COEFDIFF_LEVEL         , mesh->numel
           , COEFDIF);
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if(mesh->ndfF > 0) 
  {   
    initProp(mesh->elm.densityFluid
            ,mesh->elm.material.prop,mesh->elm.mat
            ,DENSITY_LEVEL          ,mesh->numel
            ,DENSITY); 

/*...*/
    alphaProdVector(1.e0, mesh->elm.vel0
                  , mesh->numel*ndfVel, mesh->elm.vel);
/*...*/
    alphaProdVector(1.e0, mesh->elm.pressure0
                  , mesh->numel, mesh->elm.pressure);
        
  }
/*...................................................................*/

/*...*/
  if (mesh->ndfFt > 0) 
  {
/*...*/
    alphaProdVector(1.e0          ,mesh->elm.temp0
                     ,mesh->numel ,mesh->elm.temp);
/*...................................................................*/

    if(energyModel->fTemperature)
    {
/*... convertendo temperatura para kelvin*/
      if(energyModel->fKelvin)
        convTempForKelvin(mesh->elm.temp0,mesh->numel,true); 

      alphaProdVector(1.e0        ,mesh->elm.temp
                     ,mesh->numel ,mesh->elm.energy0);

      alphaProdVector(1.e0        ,mesh->elm.energy0
                     ,mesh->numel ,mesh->elm.energy);
    }
    else
    {
      if(fComb)
        getEnergyFromTheTempMix(&propF->sHeat          ,mesh->elm.yFrac
                           ,mesh->elm.temp         ,mesh->elm.energy0
                           ,mesh->elm.material.prop,mesh->elm.mat                        
                           ,mesh->numel            ,cModel->nOfSpecies
                           ,propF->fSpecificHeat   ,energyModel->fKelvin
                           ,ompVar.fUpdate         ,ompVar.nThreadsUpdate);
      else
        getEnergyForTemp(propF
                        ,mesh->elm.temp         ,mesh->elm.energy0                     
                        ,mesh->numel            
                        ,propF->fSpecificHeat   ,energyModel->fKelvin
                        ,ompVar.fUpdate         ,ompVar.nThreadsUpdate);
/*...*/
      alphaProdVector(1.e0        ,mesh->elm.energy0
                     ,mesh->numel ,mesh->elm.energy);
    }   
/*...................................................................*/


/*...*/
    alphaProdVector(1.e0              ,mesh->elm.vel0
                   ,mesh->numel*ndfVel,mesh->elm.vel);
/*...*/
    alphaProdVector(1.e0        ,mesh->elm.pressure0
                   ,mesh->numel ,mesh->elm.pressure);
  
/*... inicializando a densidade*/
/*... mixtura*/
    if(fComb)
      initPropTempMix(propF              ,cModel
                     ,mesh->elm.densityFluid,mesh->elm.temp
                     ,mesh->elm.pressure0   ,mesh->elm.yFrac
                     ,mesh->elm.material.prop ,mesh->elm.mat
                     ,cModel->nOfSpecies    ,DENSITY_LEVEL   
                     ,mesh->numel           ,energyModel->fKelvin 
                     ,DENSITY);
/*...*/
    else
    {
      if(propF->fDensity)
        initPropTemp(propF
                  ,mesh->elm.densityFluid ,mesh->elm.temp 
                  ,mesh->elm.pressure0    ,mesh->elm.material.prop
                  ,mesh->elm.mat
                  ,DENSITY_LEVEL          ,mesh->numel
                  ,energyModel->fKelvin   ,DENSITY);
      else
        initProp(mesh->elm.densityFluid 
                ,mesh->elm.material.prop,mesh->elm.mat
                ,DENSITY_LEVEL          ,mesh->numel
                ,DENSITY);
    }
/*...................................................................*/

/*... inicializando o calor especifico*/
/*... mixtura*/
    if(fComb)
      initPropTempMix(propF              ,cModel
                    ,mesh->elm.specificHeat,mesh->elm.temp   
                    ,mesh->elm.pressure0   ,mesh->elm.yFrac  
                    ,mesh->elm.material.prop ,mesh->elm.mat  
                    ,cModel->nOfSpecies    ,SHEAT_LEVEL 
                    ,mesh->numel           ,energyModel->fKelvin 
                    ,SPECIFICHEATCAPACITYFLUID); 
    else
    {
      if(propF->fSpecificHeat)
      {
        initPropTemp(propF
                    ,mesh->elm.specificHeat   ,mesh->elm.temp 
                    ,mesh->elm.pressure0      ,mesh->elm.material.prop  
                    ,mesh->elm.mat
                    ,SHEAT_LEVEL              ,mesh->numel
                    ,energyModel->fKelvin     ,SPECIFICHEATCAPACITYFLUID);
      }
      else
        initProp(mesh->elm.specificHeat   
             ,mesh->elm.material.prop  ,mesh->elm.mat
             ,SHEAT_LEVEL              ,mesh->numel
             ,SPECIFICHEATCAPACITYFLUID);
    }
/*...................................................................*/

/*... inicializando a viscosidade dinamica*/   
/*... mixtura*/
    if(fComb)
      initPropTempMix(propF                ,cModel
                    ,mesh->elm.dViscosity    ,mesh->elm.temp   
                    ,mesh->elm.pressure0     ,mesh->elm.yFrac  
                    ,mesh->elm.material.prop ,mesh->elm.mat
                    ,cModel->nOfSpecies      ,DVISCOSITY_LEVEL
                    ,mesh->numel             ,energyModel->fKelvin  
                    ,DYNAMICVISCOSITY);
    else
    {
      if(propF->fDynamicViscosity)
      {
        initPropTemp(propF
                ,mesh->elm.dViscosity     ,mesh->elm.temp 
                ,mesh->elm.pressure0      ,mesh->elm.material.prop  
                ,mesh->elm.mat
                ,DVISCOSITY_LEVEL         ,mesh->numel
                ,energyModel->fKelvin     ,DYNAMICVISCOSITY);
      }
      else
        initProp(mesh->elm.dViscosity 
              ,mesh->elm.material.prop  ,mesh->elm.mat
              ,DVISCOSITY_LEVEL         ,mesh->numel
              ,DYNAMICVISCOSITY);
    }
/*...................................................................*/

/*... inicializando a condutividade termica*/
/*... mixtura*/
    if(fComb)
        initPropTempMix(propF               ,cModel
                    ,mesh->elm.tConductivity,mesh->elm.temp      
                    ,mesh->elm.pressure0    ,mesh->elm.yFrac 
                    ,mesh->elm.material.prop ,mesh->elm.mat  
                    ,cModel->nOfSpecies     ,TCONDUCTIVITY_LEVEL
                    ,mesh->numel            ,energyModel->fKelvin
                    ,THERMALCONDUCTIVITY);
    else
    {
      if(propF->fThermalConductivity)
        initPropTemp(propF
                    ,mesh->elm.tConductivity ,mesh->elm.temp 
                    ,mesh->elm.pressure0 
                    ,mesh->elm.material.prop   ,mesh->elm.mat
                    ,TCONDUCTIVITY_LEVEL       ,mesh->numel
                    ,energyModel->fKelvin      ,THERMALCONDUCTIVITY);
      else
        initProp(mesh->elm.tConductivity 
              ,mesh->elm.material.prop     ,mesh->elm.mat
              ,TCONDUCTIVITY_LEVEL         ,mesh->numel
              ,THERMALCONDUCTIVITY);
    }
/*...................................................................*/
 
/*...*/
    if(fComb)
    {
/*... inicializando a difusividade das especies*/
      initDiffMix(propF               , cModel
             ,mesh->elm.cDiffComb     , mesh->elm.temp 
             ,mesh->elm.pressure0     , mesh->elm.yFrac 
             ,mesh->elm.material.prop ,mesh->elm.mat   
             ,cModel->nOfSpecies      ,cModel->nComb   
             ,mesh->numel             ,energyModel->fKelvin);
/*... inicializando as entalpia das especies*/
      getEnthalpySpecies(cModel         , propF
                   , mesh->elm.enthalpyk, mesh->elm.temp 
                   , mesh->numel        , energyModel->fKelvin
                   , ompVar.fUpdate     , ompVar.nThreadsUpdate);
/*...................................................................*/

    }
/*...................................................................*/
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
      fprintf(fileLogExc,"parametro: %s faltando.\n"
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
 
  for(i=0;i<nn;i++)
  {
    fscanf(file,"%ld",&idum);
    k = (INT) idum -1;
    for(j=0;j<ndm;j++)
    {
      fscanf(file,"%lf",&MAT2D(k,j,x,ndm));
    }
  }

//alphaProdVector(0.1,x,nn*ndm,x);
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
/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 15/09/2019                                   *
 * ------------------------------------------------------------------*
 * READVFLOADS : leitura da definicoes das cargas                    *
 * ------------------------------------------------------------------*
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
    loads[nLoad-1].type        = type;
    loads[nLoad-1].fUse        = true;
    error = fscanf(file,"%d",&type);
    loads[nLoad-1].nTypeVar = type;
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
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 04/08/2019                                   *
 * ------------------------------------------------------------------*
 * CONFIG : configuraceos gerais                                     *
 * ------------------------------------------------------------------*
 * Parametro de entrada:                                             *
 * ----------------------------------------------------------------- *
 * fileOpt   - opcoes de arquivo                                     *
 * reordMesh - reordenacao do malha                                  *
 * m         - memoria principal                                     *
 * file  - ponteiro para o arquivo de dados                          *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 * lRcell    -> residuo por celula                                   *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void config(FileOpt *opt,Reord *reordMesh,FILE* file)
{
  char config[][WORD_SIZE]={"reord"     ,"memory" 
                           ,"fItPlotRes","fItPlot"};
  
  char word[WORD_SIZE];
  char s[WORD_SIZE];
  bool flag[NCONFIG];
  int i=0,j,temp;
  int conv;


  for(j=0;j<NCONFIG;j++)
    flag[j] = false;

  while(i < NCONFIG && i < 20){
    readMacro(file,word,false);
/*... reord*/   
    if(!strcmp(word,config[0])){
      fscanf(file,"%s",s);
      if(!strcmp(s,"true"))
        reordMesh -> flag= true;
      else
        reordMesh -> flag= false;
          
      flag[0] = true;
      i++;
      if(reordMesh->flag)
        fprintf(fileLogExc,"%-20s: %s\n","Reord","true");
      else
        fprintf(fileLogExc,"%-20s: %s\n","Reod","false");
    }
/*... mem*/   
    else if(!strcmp(word,config[1])){
      fscanf(file,"%d",&temp);
//    conv    = CONV_BYTES*CONV_BYTES;
      conv    = 1024*1024;
 			nmax    = (iptx) temp;
			nmax    =  nmax*conv;
      flag[1] = true;
      i++;
      
      fprintf(fileLogExc,"%-20s: %d\n","Memory(MBytes)"
              ,(int)(nmax/conv));
    }

/*... fItPlotRes*/   
    else if(!strcmp(word,config[2])){
      fscanf(file,"%s",s);
      if(!strcmp(s,"true")){
        opt->fItPlotRes = true;
        fprintf(fileLogExc,"%-20s: %s\n","fItPlotRes","true");
      }
      else{
        opt->fItPlotRes = false;
        fprintf(fileLogExc,"%-20s: %s\n","fItPlotRes","false");
      }
      flag[2] = true;
      i++;
    }
/*... fItPlot*/   
    else if(!strcmp(word,config[3])){
      fscanf(file,"%s",s);
      if(!strcmp(s,"true")){
        opt->fItPlot = true;
        fprintf(fileLogExc,"%-20s: %s\n","fItPlot","true");
      }
      else{
        opt->fItPlot = false;
        fprintf(fileLogExc,"%-20s: %s\n","fItPlot","false");
      }
      flag[3] = true;
      i++;
    }
    else
      i++;
/*...................................................................*/
  }
  readMacro(file,word,true);
  readMacro(file,word,true);
  for(j=0;j<NCONFIG;j++){
    if(!flag[j]){
      fprintf(fileLogExc,"%s: %s faltando.\n"
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
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------*
 * INITPROP: inicializao de propriedades com variacao temporal       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * prop    -> nao definido                                           * 
 * propMat -> propriedade de referencia por material                 * 
 * mat     -> material por celula                                    * 
 * np      -> numero niveis de tempos                                * 
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
void initProp(DOUBLE *RESTRICT prop 
             ,DOUBLE *RESTRICT propMat,short *RESTRICT mat
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
 * Data de criacao    : 11/11/2017                                   *
 * Data de modificaco : 05/05/2018                                   *
 *-------------------------------------------------------------------*
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

  short i;
  char *str={"endEdp"};
  char word[WORD_SIZE];;
  char macros[][WORD_SIZE] =
       { "transport1" ,"transport2","transport3"    /* 0, 1, 2*/
         ,"diffusion1","diffusion2","diffusion3"    /* 3, 4, 5*/
         ,"fluid"     ,"fluidt"    ,""};            /* 6, 7, 8*/
 
/*...*/
	for (i = 0; i < MAX_TRANS_EQ; i++)
		mesh->ndfT[i] = 0;
	for (i = 0; i < MAX_DIF_EQ; i++)
		mesh->ndfD[i] = 0;
	mesh->ndfF  = 0;
	mesh->ndfFt = 0;
/*...................................................................*/

  readMacro(file,word,false);
  do{
/*... transport1*/
    if(!strcmp(word, macros[0])){
      readMacro(file,word,false);
      mesh->ndfT[0]=(short) atol(word);
      fprintf(fileLogExc,"transport1 ndf %d\n",mesh->ndfT[0]);
    }
/*...................................................................*/

/*... transport2*/
    else if(!strcmp(word, macros[1])){
      readMacro(file,word,false);
      mesh->ndfT[1]= (short) atol(word);
      fprintf(fileLogExc,"transport2 ndf %d\n",mesh->ndfT[1]);
    }
/*...................................................................*/

/*... transport3*/
    else if(!strcmp(word, macros[2])){
      readMacro(file,word,false);
      mesh->ndfT[2]= (short) atol(word);
      fprintf(fileLogExc,"transport3 ndf %d\n",mesh->ndfT[2]);
    }
/*...................................................................*/

/*... diffusion1*/
    else if(!strcmp(word, macros[3])){
      readMacro(file,word,false);
      mesh->ndfD[0]= (short) atol(word);
      fprintf(fileLogExc,"diffusion1 ndf %d\n",mesh->ndfD[0]);
    }
/*...................................................................*/

/*... diffusion2*/
    else if(!strcmp(word, macros[4])){
      readMacro(file,word,false);
      mesh->ndfD[1]= (short) atol(word);
      fprintf(fileLogExc,"diffusion2 ndf %d\n",mesh->ndfD[1]);
    }
/*...................................................................*/

/*... diffusion3*/
    else if(!strcmp(word, macros[5])){
      readMacro(file,word,false);
      mesh->ndfD[2]= (short) atol(word);
      fprintf(fileLogExc,"diffusion3 ndf %d\n",mesh->ndfD[2]);
    }
/*...................................................................*/

/*... fluido*/
    else if(!strcmp(word,macros[6])){
      readMacro(file,word,false);
      mesh->ndfF= (short) atol(word);
      fprintf(fileLogExc,"fluid ndf %d\n",mesh->ndfF);
    }
/*...................................................................*/

/*... fluido-termo ativado*/
    else if(!strcmp(word, macros[7])){
      readMacro(file,word,false);
      mesh->ndfFt= (short) atol(word);
      fprintf(fileLogExc,"fluidt ndf %d\n",mesh->ndfFt);
    }
/*...................................................................*/
    
    readMacro(file,word,false);
  }while(strcmp(word,str));


}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 29/08/2017                                   *
 * Data de modificaco : 10/07/2018                                   *
 *-------------------------------------------------------------------* 
 * READPROPVAR : propriedades variaveis                              * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * pFluid  ->                                                        *
 * file    -> arquivo de arquivo                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * p       ->                                                        * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void readPropVarFluid(PropVarFluid *p,FILE *file){

  short kk = 0; 
  char *str={"endfluid"};
  char macros[][WORD_SIZE] =
       { "sheat"       ,"density"   ,"dviscosity"   /* 0, 1, 2*/
        ,"tcondutivity",""          ,""         };  /* 3, 4, 5*/
  char word[WORD_SIZE];

/*...*/
  p->fDensity           = false;
  p->fSpecificHeat      = false;
  p->fDynamicViscosity  = false;
  p->fThermalConductivity = false;
/*...................................................................*/

  readMacroV2(file,word,false,true);
  do{    
/*... specific heat*/
    if(!strcmp(word,macros[0])){
      readMacroV2(file, word, false, true);
      p->fSpecificHeat = true;
      if(p->fSpecificHeat) 
        fprintf(fileLogExc,"%-25s: %s\n","sHeat variation","Enable\n");
      initSheatPol(&p->sHeat, word, file);
    }
/*...................................................................*/

/*... densidade*/
    else if(!strcmp(word,macros[1])){
      readMacroV2(file, word, false, true);
      p->fDensity = true;
      if(p->fDensity) 
        fprintf(fileLogExc,"%-25s: %s\n","Density variation","Enable");
      initDensityPol(&p->den,word,file);
    }
/*...................................................................*/

/*... viscosidade dinamica*/
    else if(!strcmp(word, macros[2])){
      readMacroV2(file, word, false, true);
      p->fDynamicViscosity = true;  
      if(p->fDynamicViscosity)
        fprintf(fileLogExc,"%-25s: %s\n","dViscosity variation"
                          ,"Enable");
      initDviscosityPol(&p->dVisc, word, file);
    }
/*...................................................................*/

/*... condutiviade termica*/
    else if(!strcmp(word, macros[3])){
      readMacroV2(file, word, false, true);
      p->fThermalConductivity = true;
      if(p->fThermalConductivity)
        fprintf(fileLogExc,"%-25s: %s\n","tCondutivity variation"
                                        ,"Enable"); 
      initThCondPol(&p->thCond, word, file);
    }
/*...................................................................*/
    
    readMacroV2(file, word, false, true);
    kk++;
  }while(strcmp(word,str) && kk < 1000);

/*...*/
  if(kk == 1000)
    ERRO_READ_LOOP(__LINE__,__FILE__,__func__,EXIT_READ_LOOP);
  /*...................................................................*/
}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 19/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------* 
 * readPropVarMixture : propriedades variaveis de uma mistura gasosa * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * p       ->                                                        *
 * file    -> arquivo de arquivo                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * p       ->                                                        * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void readPropVarMixture(PropVarFluid *p,Combustion *cModel,FILE *file)
{

  char *str={"endmixture"};
  char macros[][WORD_SIZE] =
       { "sheat"       ,"density"   ,"dviscosity"   /* 0, 1, 2*/
        ,"tcondutivity","diffusion" ,""         };  /* 3, 4, 5*/
  char word[WORD_SIZE];

/*...*/
  p->fDensity             = false;
  p->fSpecificHeat        = false;
  p->fDynamicViscosity    = false;
  p->fThermalConductivity = false;
  p->fDiffusion           = false;
/*...................................................................*/

  readMacroV2(file,word,false,true);
  do
  {
/*... specific heat*/
    if(!strcmp(word,macros[0]))
    {
      readMacroV2(file, word, false, true);
      p->fSpecificHeat = true;
      if(p->fSpecificHeat) 
        fprintf(fileLogExc,"%-25s: %s\n","sHeat variation","Enable\n");
      initMixtureSpeciesfiHeat(&p->sHeat, word,cModel, file);
      
    }
/*...................................................................*/

/*... densidade*/
    else if(!strcmp(word,macros[1]))
    {
      readMacroV2(file, word, false, true);
      p->fDensity = true;
      if(p->fDensity) 
        fprintf(fileLogExc,"%-25s: %s\n","Density variation","Enable");
      initDensityPol(&p->den,word,file);
    }
/*...................................................................*/

/*... viscosiade dinamica*/
    else if(!strcmp(word, macros[2]))
    {
      readMacroV2(file, word, false, true);
      p->fDynamicViscosity = true;  
      if(p->fDynamicViscosity)
        fprintf(fileLogExc,"%-25s: %s\n","dViscosity variation"
                          ,word);
      initDviscosityPol(&p->dVisc, word, file);
    }
/*...................................................................*/

/*... condutividade termica*/
    else if(!strcmp(word, macros[3])){
      readMacroV2(file, word, false, true);
      p->fThermalConductivity = true;
      if(p->fThermalConductivity)
        fprintf(fileLogExc,"%-25s: %s\n","tCondutivity variation"
                                        ,word); 
      initThCondPol(&p->thCond, word, file);
    }
/*...................................................................*/
    
/*... coeficient de difusao*/
    else if(!strcmp(word, macros[4])){
      readMacroV2(file, word, false, true);
      p->fDiffusion = true;
      if(p->fDiffusion)
        fprintf(fileLogExc,"%-25s: %s\n","Diffusion variation"
                                        ,word); 
      initDiffSp(&p->diff, word, file);
    }
/*...................................................................*/

    readMacroV2(file, word, false, true);
  }while(strcmp(word,str));

}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 12/05/2018                                   *
 * Data de modificaco : 10/07/2018                                   *
 *-------------------------------------------------------------------*
 * READPROPVAR : propriedades variaveis                              *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * p       ->                                                        *
 * file    -> arquivo de arquivo                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * p       ->                                                        *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void readPropVarDiff(PropVarCD *p, FILE *file)
{

  char *str = { "enddiff" };
  char macros[][WORD_SIZE] =
             { "densityd1"   ,"cdiffd1"};  /* 0, 1, 2*/
  char word[WORD_SIZE];

/*...*/
  p[0].fDensity = false;
  p[0].fCeofDiff = false;
/*...................................................................*/

  readMacroV2(file, word, false, true);
  do 
  {
/*... density D1*/
    if (!strcmp(word, macros[0]))
    {
      readMacroV2(file, word, false, true);
      p[0].fDensity = true;
      initCdPol(&p[0].den, word,file);
      if (p[0].fDensity)
        fprintf(fileLogExc, "%-25s: %s\n", "DensityD1 variation"
                                         , "Enable");
    }
/*...................................................................*/

/*... condutividade termica D1*/
    else if (!strcmp(word, macros[1]))
    {
      readMacroV2(file, word, false, true);
      p[0].fCeofDiff = true;
      initCdPol(&p[0].ceofDiff, word, file);
      if (p[0].fCeofDiff)
        fprintf(fileLogExc, "%-25s: %s\n", "CeofDiff D1 variation"
                                         , "Enable");
    }
/*...................................................................*/
    readMacroV2(file, word, false, true);
  } while (strcmp(word, str));

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/06/2018                                   *
 * Data de modificaco : 10/07/2018                                   *
 *-------------------------------------------------------------------*
 * READPROPVAR : propriedades variaveis                              *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * p       ->                                                        *
 * file    -> arquivo de arquivo                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * p       ->                                                        *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void readPropVarTrans(PropVarCD *p, FILE *file)
{

  char *str = { "endtrans" };
  char macros[][WORD_SIZE] =  { "densityt1","cdifft1" };  /* 0, 1, 2*/
  char word[WORD_SIZE];

/*...*/
  p[0].fDensity = false;
  p[0].fCeofDiff = false;
/*...................................................................*/

  readMacroV2(file, word, false, true);
  do
  {
 /*... density T1*/
    if (!strcmp(word, macros[0]))
    {
      readMacroV2(file, word, false, true);
      p[0].fDensity = true;
      initCdPol(&p[0].den, word, file);
      if (p[0].fDensity)
        fprintf(fileLogExc, "%-25s: %s\n", "DensityD1 variation"
          , "Enable");
    }
/*...................................................................*/

/*... condutividade termica T1*/
    else if (!strcmp(word, macros[1]))
    {
      readMacroV2(file, word, false, true);
      p[0].fCeofDiff = true;
      initCdPol(&p[0].ceofDiff, word, file);
      if (p[0].fCeofDiff)
        fprintf(fileLogExc, "%-25s: %s\n", "CeofDiff D1 variation"
          , "Enable");
    }
/*...................................................................*/
    readMacroV2(file, word, false, true);
  } while (strcmp(word, str));
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/06/2018                                   *
 * Data de modificaco : 26/05/2019                                   *
 *-------------------------------------------------------------------*
 * READPROPVAR : propriedades variaveis                              *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * pf      -> fluido                                                 *
 * pd      -> difusao                                                *
 * pt      -> transporte                                             *
 * file    -> arquivo de arquivo                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * pd      -> difusao                                                *
 * pt      -> transporte                                             *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void readPropVar(PropVarFluid *pf  , PropVarCD *pd, PropVarCD *pt
               , Combustion *cModel,FILE *file)
{

  char *str = { "endpropvar" };
  char macros[][WORD_SIZE] = { "diff"   ,"trans","fluid"  /* 0, 1, 2*/
                              ,"mixture"};                /* 3      */
  char word[WORD_SIZE];

  readMacroV2(file, word, false, true);
  do
  {
/*... diff*/
    if (!strcmp(word, macros[0]))
      readPropVarDiff(pd,file);
/*...................................................................*/

/*... Trans*/
    else if (!strcmp(word, macros[1]))
      readPropVarTrans(pt, file);
 /*...................................................................*/

/*... Fluid*/
    else if (!strcmp(word, macros[2]))
      readPropVarFluid(pf, file);
/*...................................................................*/

/*... mixute*/
    else if (!strcmp(word, macros[3]))
      readPropVarMixture(pf,cModel, file);
/*...................................................................*/

    readMacroV2(file, word, false, true);
  } while (strcmp(word, str));

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 04/09/2017                                   *
 * Data de modificaco : 07/06/2019                                   *
 *-------------------------------------------------------------------* 
 * readMode : le as configuraoes dos modelos                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * e       -> modelos/termos usandos na eq energia                   * 
 * t       -> modelo de turbilencia                                  *
 * eMass   -> modelos/termos usados na equacao da conv de mass       * 
 * momentumModel  -> modelos/termos usados na equacao da conv de mass*
 * dModel  -> modelos/termos usados nas equacoes de diffusao         *
 * tModel  -> modelos/termos usados nas equacoes de transporte       *
 * cModel  -> modelos/termos usados na equacoes de combustao         *
 * file    -> arquivo de arquivo                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * e       ->                                                        * 
 * t       ->                                                        * 
 * eModel  ->                                                        * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void readModel(EnergyModel *e         , Turbulence *t
             , MassEqModel *eMass     , MomentumModel *momentumModel
             , DiffModel   *dModel    , TransModel *tModel
             , Combustion *cModel
             , FILE *file){

  char *str={"endmodel"};
  char format[1024];
  char word[WORD_SIZE];
  char macros[][WORD_SIZE] = {"energy"  ,"turbulence","mass"
                             ,"momentum","diffusion" ,"transport"
                             ,"combustion"};

  char energy[][WORD_SIZE] = { "preswork"  , "dissipation" /*0,1*/
                             , "residual"  , "absolute"    /*2,3*/
                             ,"temperature", "entalphy"    /*4,5*/
                             , "diffenergy" };             /*6, */  

  char turb[][WORD_SIZE] = { "smagorinsky","wallmodel" , "wale"     
                            ,"vreman"     ,"ldynamic"  , "sigmamodel"
                            ,"mixed"      ,"bardina"   , "clark"
                            ,"bardinaMod" ,"towdynamic", "gdynamic"
                            ,"gdynamicmod","oneeqk"}; 

  char mass[][WORD_SIZE] = { "lhsdensity","rhsdensity"}; 

  char momentum[][WORD_SIZE] = {"residual"   ,"absolute"       /*0,1*/                    
                               ,"rhiechow"   ,"viscosity"      /*2,3*/
                               ,"div"        ,"buoyanthy"      /*4,5*/  
                               ,"buoyantprgh","buoyantrhoref"  /*6,7*/
                               ,"sopressure" };                /*8*/

  char combustion[][WORD_SIZE] = {"residual"   ,"absolute"      /*0,1*/                    
                                 ,"grouped"    ,"ungrouped"     /*2,3*/
                                 ,"edc"        ,"arrhenius"     /*4,5*/ 
                                 ,"hcombustion","hformation"    /*6,7*/   
                                 ,"correctvel" ,"edm"};         /*8*/
  
  char diff[][WORD_SIZE] = { "residual","absolute"};        /*0,1*/
  char tran[][WORD_SIZE] = { "residual","absolute" };       /*0,1*/

  char typeWallModel[][WORD_SIZE] ={"standard","enhanced"};
  
  int i,nPar,id;

  readMacroV2(file, word, false, true);
  do{
/*... equacao da energia*/
    if(!strcmp(word,macros[0])){
/*...*/
      e->fPresWork    = false;
      e->fDissipation = false;
      e->fRes         = false;
      e->fTemperature = false;
      e->fDiffEnergy  = false;
      fprintf(fileLogExc,"\n%-20s: \n","EnergyModel");  
/*...................................................................*/      
      fscanf(file,"%d",&nPar);
      for(i=0;i<nPar;i++){
        strcpy(format,"%-20s: %s\n");
        readMacroV2(file,word,false,true);
/*... PresWork*/
        if(!strcmp(word,energy[0])){    
          e->fPresWork = true;
          if(e->fPresWork) 
            fprintf(fileLogExc,format,"PresWork","Enable");
        }
/*...................................................................*/

/*... Dissipation*/
        else if(!strcmp(word,energy[1])){        
          e->fDissipation = true;
          if(e->fDissipation) 
            fprintf(fileLogExc,format,"Dissipation","Enable");                           
        }
/*...................................................................*/

/*... Residual*/
        else if(!strcmp(word,energy[2])){
          e->fRes = true;
          if(e->fRes)
            fprintf(fileLogExc,format,"Residual","Enable");
        }
/*...................................................................*/

/*... Absolute*/
        else if(!strcmp(word,energy[3])){
          e->fRes = false;
          if(!e->fRes)
            fprintf(fileLogExc,format,"Absolute","Enable");
        }
/*...................................................................*/

/*... Temperatura*/
        else if(!strcmp(word,energy[4])){
          e->fTemperature = true;
          if(e->fTemperature)
            fprintf(fileLogExc,format,"Temperatura","Enable");
        }
/*...................................................................*/

/*... Entalphy*/
        else if(!strcmp(word,energy[5])){
          e->fTemperature = false;
          if(!e->fTemperature)
            fprintf(fileLogExc,format,"Entalphy","Enable");
        }
/*...................................................................*/

/*... Diffusion energy source*/
        else if(!strcmp(word,energy[6])){
          e->fDiffEnergy = true;
          if(!e->fTemperature)
            fprintf(fileLogExc,format,"Diffusion Energy","Enable");
        }
/*...................................................................*/

      }
/*...................................................................*/
    }
/*...................................................................*/

/*... turbulencia*/
    else if(!strcmp(word,macros[1])){   
      fprintf(fileLogExc,"\n%-20s: \n","TurbulenceModel");   
      fscanf(file,"%d",&nPar);
      for(i=0;i<nPar;i++){
        readMacroV2(file, word, false, true);
/*... Smagorinsky*/
        if(!strcmp(word,turb[0])){
          t->fDynamic             = false;
          t->fTurb               = true;      
          t->type                = LES;
          t->typeMixed[FUNMODEL] = SMAGORINSKY;
          t->typeLes             = LESFUNCMODEL; 
          fscanf(file,"%lf",&t->cf);    
          fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[0],t->cf); 
        }
/*...................................................................*/    

/*... wallModel*/
        else if(!strcmp(word,turb[1])){
          t->fWall = true;         
          readMacroV2(file, word, false, true);
/*...*/
          if(!strcmp(word,typeWallModel[0])){
            t->wallType  = STANDARDWALL; 
          }
/*...................................................................*/

/*...*/
          else if(!strcmp(word,typeWallModel[1])){
            t->wallType  = ENHANCEDWALL; 
          }
/*...................................................................*/
          fprintf(fileLogExc,"%-20s: %s\n",turb[1]
                            ,typeWallModel[t->wallType-1]);
        }
/*...................................................................*/

/*... Wale*/
        else if(!strcmp(word,turb[2])){
          t->fDynamic                = false;
          t->fTurb                  = true;     
          t->type                   = LES;
          t->typeMixed[FUNMODEL]    = WALEMODEL; 
          t->typeLes = LESFUNCMODEL;
          fscanf(file,"%lf",&t->cf);    
          fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[2],t->cf); 
        }
/*...................................................................*/ 

/*... Vreman*/
        else if(!strcmp(word,turb[3])){
          t->fDynamic               = false;
          t->fTurb                  = true;      
          t->type                   = LES;
          t->typeMixed[FUNMODEL]    = VREMAN;
          t->typeLes = LESFUNCMODEL;
          fscanf(file,"%lf",&t->cf);    
          fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[3],t->cf); 
        }
/*...................................................................*/

/*... Dynamic*/
        else if(!strcmp(word,turb[4])){
          t->fDynamic    = true;
          t->fTurb       = true;      
          t->type        = LES;
          t->typeDynamic = LDYNAMIC;
          fscanf(file,"%lf",&t->c);  
          fprintf(fileLogExc,"%-20s:\n", "lDynamic"); 
          setDynamicModelLes(t,file);
        }
/*...................................................................*/ 

/*... sigmaModel*/
        else if(!strcmp(word,turb[5])){
          t->fDynamic             = false;
          t->fTurb               = true;      
          t->type                = LES;
          t->typeMixed[FUNMODEL] = SIGMAMODEL;
          t->typeLes = LESFUNCMODEL;
          fscanf(file,"%lf",&t->cf);   
          fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[5],t->cf); 
        }
/*...................................................................*/ 

/*... mixed*/
        else if(!strcmp(word,turb[6])){
          t->fDynamic    = false;
          t->fTurb       = true;      
          t->type        = LES;
          t->typeLes     = LESMIXEDMODEL; 
          t->fTurbStruct = true;
          fscanf(file,"%lf",&t->c);  
          fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[6],t->c); 
          setMixedModelLes(t,file);
        }
/*...................................................................*/

/*... bardina*/
        else if(!strcmp(word,turb[7])){
          t->fDynamic            = false;
          t->fTurb               = true;   
          t->type                = LES;
          t->typeMixed[ESTMODEL] = BARDINA;
          t->typeLes             = LESSTRUMODEL;
          t->fTurbStruct         = true;
          fscanf(file,"%lf",&t->cs);  
          fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[7],t->cs); 
        }
/*...................................................................*/  

/*... clark*/
        else if(!strcmp(word,turb[8])){
          t->fDynamic            = false;
          t->fTurb               = true;      
          t->type                = LES;
          t->typeMixed[ESTMODEL] = CLARK;
          t->typeLes             = LESSTRUMODEL;
          t->fTurbStruct         = true;
          fscanf(file,"%lf",&t->cs);  
          fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[8],t->cs); 
        }
/*...................................................................*/ 

/*... bardinaMod*/
        else if(!strcmp(word,turb[9])){
          t->fDynamic            = false;
          t->fTurb               = true;      
          t->type                = LES;
          t->typeMixed[ESTMODEL] = BARDINAMOD;
          t->typeLes             = LESSTRUMODEL;
          t->fTurbStruct         = true;
          fscanf(file,"%lf",&t->cs);  
          fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[9],t->cs); 
        }
/*...................................................................*/   

/*... mixed 2 paramentros*/
        else if(!strcmp(word,turb[10])){
          t->fDynamic    = true;
          t->fTurb       = true;      
          t->type        = LES;
          t->typeLes     = LESMIXEDTWOMODEL; 
          t->typeDynamic = TWOPARDYNAMIC;
          t->fTurbStruct = true;
          fscanf(file,"%lf",&t->c);  
          fprintf(fileLogExc,"%-20s:\n", "TowDynamic"); 
          setMixedModelLes(t,file);
        }
/*...................................................................*/

/*... Gdynamic*/
        else if(!strcmp(word,turb[11])){
          t->fDynamic    = true;
          t->fTurb       = true;      
          t->type        = LES;
          t->typeDynamic = GDYNAMIC;
          fscanf(file,"%lf",&t->c);  
          fprintf(fileLogExc,"%-20s:\n","gDynamic"); 
          setDynamicModelLes(t,file);
        }
/*...................................................................*/ 

/*... Gdynamic*/
        else if(!strcmp(word,turb[12])){
          t->fDynamic    = true;
          t->fTurb       = true;      
          t->type        = LES;
          t->typeDynamic = GDYNAMICMOD;
          fscanf(file,"%lf",&t->c);  
          fprintf(fileLogExc,"%-20s:\n","gDynamicMod"); 
          setDynamicModelLes(t,file);
        }
/*...................................................................*/

/*... oneEqk*/
        else if(!strcmp(word,turb[13])){
          t->fDynamic    = true;
          t->fTurb       = true;      
          t->type        = LES;
          t->fOneEq      = true;
          t->typeMixed[FUNMODEL] = ONEEQK;
          t->typeLes = LESFUNCMODELONEEQK;
          t->eK.fRes =  true;
          fscanf(file,"%lf %lf %lf %d %lf",&(t->eK.ck)
                                          ,&(t->eK.ce)
                                          ,&(t->eK.sk)
                                          ,&(t->eK.maxIt)
                                          ,&(t->eK.tol) );  
          strcpy(format,"%-20s: Ck = %.3lf Ce = %.3lf Sk = %.3lf\n");
          fprintf(fileLogExc,format,"oneEqK",t->eK.ck,t->eK.ce,t->eK.sk); 
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... mass*/
    else if(!strcmp(word,macros[2])){ 
      strcpy(format,"%-20s: %s\n");
      fprintf(fileLogExc,"\n%-20s: \n","MassEqModel");    
      eMass->LhsDensity = false;
      eMass->RhsDensity = false;
      fscanf(file,"%d",&nPar);
      for(i=0;i<nPar;i++){
        readMacroV2(file, word, false, true);
/*... LhsDensity*/
        if(!strcmp(word,mass[0])){
          eMass->LhsDensity = true;          
          if(eMass->LhsDensity){ 
            fprintf(fileLogExc,format,"LhsDensity","Enable");
          }
        }
/*...................................................................*/

/*... RhsDensity*/
        else if(!strcmp(word,mass[1])){
          eMass->RhsDensity = true;          
          if(eMass->RhsDensity){ 
            fprintf(fileLogExc,format,"RhsDensity","Enable");
          }
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... Momentum*/
    else if(!strcmp(word,macros[3]))
    { 
      strcpy(format,"%-20s: %s\n");
      fprintf(fileLogExc,"\n%-20s: \n","MomentumEqModel");    
      momentumModel->fRes             = false;
      momentumModel->fRhieChowInt     = false;
      momentumModel->fViscosity       = false;
      momentumModel->fDiv             = false;
      momentumModel->fSoPressure      = false;
      fscanf(file,"%d",&nPar);
      for(i=0;i<nPar;i++)
      {
        readMacroV2(file, word, false, true);
/*... residual*/
        if(!strcmp(word,momentum[0]))
        {
          momentumModel->fRes = true;          
          if(momentumModel->fRes) 
            fprintf(fileLogExc,format,"Residual","Enable");
        }
/*...................................................................*/

/*... Absolute*/
        else if(!strcmp(word,momentum[1]))
        {
          momentumModel->fRes = false;            
          if(!momentumModel->fRes) 
            fprintf(fileLogExc,format,"Absolute","Enable");
        }
/*...................................................................*/

/*... RhieChow*/
        else if(!strcmp(word,momentum[2]))
        {
          momentumModel->fRhieChowInt = true;            
          if(momentumModel->fRhieChowInt) 
            fprintf(fileLogExc,format,"RhieChowInt","Enable");
        }
/*...................................................................*/

/*... Viscosity*/
        else if(!strcmp(word,momentum[3]))
        {
          momentumModel->fViscosity = true;            
          if(momentumModel->fViscosity) 
            fprintf(fileLogExc,format,"Viscosity","Enable");
        }
/*...................................................................*/

/*... div*/
        else if(!strcmp(word,momentum[4]))
        {
          momentumModel->fDiv = true;            
          if(momentumModel->fDiv) 
            fprintf(fileLogExc,format,"Divergente","Enable");
        }
/*...................................................................*/

/*... bouyant_hydrostatic*/
        else if(!strcmp(word,momentum[5]))
        {
          momentumModel->iCodBuoyant = BUOYANT_HYDROSTATIC;            
          fprintf(fileLogExc,format,"bouyant_hydrostatic","Enable");
        }
/*...................................................................*/

/*... bouyant_prgh*/
        else if(!strcmp(word,momentum[6]))
        {
          momentumModel->iCodBuoyant = BUOYANT_PRGH;            
          fprintf(fileLogExc,format,"bouyant_prgh","Enable");
        }
/*...................................................................*/

/*... bouyant_prgh*/
        else if(!strcmp(word,momentum[7]))
        {
          momentumModel->iCodBuoyant = BUOYANT_RHOREF;            
          fprintf(fileLogExc,format,"bouyant_rofref","Enable");
        }
/*...................................................................*/

/*... soPressure*/
        else if(!strcmp(word,momentum[8]))
        {
          momentumModel->fSoPressure = true;            
          fprintf(fileLogExc,format,"SoPressure","Enable");
        }
/*...................................................................*/
      }
    }
/*...................................................................*/

/*... Diffusion*/
    else if (!strcmp(word, macros[4])) 
    {
      strcpy(format, "%-20s: %s\n");
      readMacroV2(file, word, false, true);
/*...*/
      if(!strcmp(word,"d1"))
        id = 0;
      else if (!strcmp(word,"d2"))
        id = 1;
      else if (!strcmp(word, "d3"))
        id = 2;  
/*...................................................................*/

      fprintf(fileLogExc, "\n%-20s: D%1d\n", "DiffusionMode",id+1);

      dModel[id].fRes = false;
      fscanf(file, "%d", &nPar);
      for (i = 0; i<nPar; i++) 
      {
        readMacroV2(file, word, false, true);
/*... residual*/
        if (!strcmp(word, diff[0])) 
        {
          dModel[id].fRes = true;
          if (dModel[id].fRes)
            fprintf(fileLogExc, format, "Residual", "Enable");
        }
/*...................................................................*/

/*... Absolute*/
        else if (!strcmp(word, diff[1]))
        {
          dModel[id].fRes = false;
          if (!dModel[id].fRes)
            fprintf(fileLogExc, format, "Absolute", "Enable");
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... Transport*/
    else if (!strcmp(word, macros[5]))
    {
      strcpy(format, "%-20s: %s\n");
      readMacroV2(file, word, false, true);
/*...*/
      if (!strcmp(word, "t1"))
        id = 0;
      else if (!strcmp(word, "t2"))
        id = 1;
      else if (!strcmp(word, "t3"))
        id = 2;
/*...................................................................*/

      fprintf(fileLogExc, "\n%-20s: D%1d\n", "TransportMode", id + 1);

      tModel[id].fRes = false;
      fscanf(file, "%d", &nPar);
      for (i = 0; i<nPar; i++)
      {
        readMacroV2(file, word, false, true);
/*... residual*/
        if (!strcmp(word, tran[0]))
        {
          tModel[id].fRes = true;
          if (tModel[id].fRes)
            fprintf(fileLogExc, format, "Residual", "Enable");
        }
/*...................................................................*/

/*... Absolute*/
        else if (!strcmp(word, tran[1]))
        {
          tModel[id].fRes = false;
          if (!tModel[id].fRes)
            fprintf(fileLogExc, format, "Absolute", "Enable");
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... Combustion*/
    else if (!strcmp(word, macros[6]))
    {
      strcpy(format,"%-20s: %s\n");
      fprintf(fileLogExc, "\n%-20s:\n", "CombustionModel");

      cModel->fCombustion  = true;
      fscanf(file, "%d", &nPar);
      for (i = 0; i<nPar; i++)
      {
        strcpy(format,"%-20s: %s\n");
        readMacroV2(file, word, false, true);
/*... residual*/
        if (!strcmp(word, combustion[0]))
        {
          cModel->fRes = true;
          if (cModel->fRes)
            fprintf(fileLogExc, format, "Residual", "Enable");
        }
/*...................................................................*/

/*... Absolute*/
        else if (!strcmp(word, combustion[1]))
        {
          cModel->fRes = false;
          if (!cModel->fRes)
            fprintf(fileLogExc, format, "Absolute", "Enable");
        }
/*...................................................................*/

/*... grouped (lumped)*/
        else if (!strcmp(word, combustion[2]))
        {
          cModel->fLump                = true;
          fscanf(file,"%hd %hd %hd",&(cModel->nComb)
                                   ,&(cModel->nOfSpecies)
                                   ,&(cModel->nOfSpeciesLump));
          strcpy(format,"%-20s: nComb= %d nOfSp = %d nOfSpLp = %d\n");
          fprintf(fileLogExc,format,"grouped species"
                              ,cModel->nComb
                              ,cModel->nOfSpecies 
                              ,cModel->nOfSpeciesLump); 
        }
/*...................................................................*/

/*... ungrouped*/
        else if (!strcmp(word, combustion[3]))
        {
          cModel->fLump                = false;
          fscanf(file,"%hd %hd %hd",&(cModel->nComb)
                                   ,&(cModel->nOfSpecies)
                                   ,&(cModel->nOfSpeciesLump)); 

          strcpy(format,"%-20s: nComb = %d nOfSp = %d nOfSpLp = %d\n");
          fprintf(fileLogExc,format,"ungrouped species"
                              ,cModel->nComb
                              ,cModel->nOfSpecies 
                              ,cModel->nOfSpeciesLump);     
        }
/*...................................................................*/

/*... edu*/
        else if (!strcmp(word, combustion[4]))
        {          
          cModel->reactionKinetic = EDC;
          setEdc(&cModel->edc,file);
          fprintf(fileLogExc, format, "EDC", "Enable");          
        }
/*...................................................................*/

/*... arrhenius*/
        else if (!strcmp(word, combustion[5]))
        {
          cModel->reactionKinetic           = ARRHENIUS;
          fprintf(fileLogExc, format, "ARRHENIUS", "Enable");          
        }
/*...................................................................*/

/*... hcombustion*/
        else if (!strcmp(word, combustion[6]))
        {
          cModel->typeHeatRealese  = HCOMBUSTION; 
          fprintf(fileLogExc, format, "Hcombustion", "Enable");          
        }
/*...................................................................*/

/*... hformation*/
        else if (!strcmp(word, combustion[7]))
        {
          cModel->typeHeatRealese  = HFORMATION; 
          fprintf(fileLogExc, format, "HFormation", "Enable");          
        }
/*...................................................................*/

/*... correctVel*/
        else if (!strcmp(word, combustion[8]))
        {
          cModel->fCorrectVel  = true; 
          if (cModel->fCorrectVel)
            fprintf(fileLogExc, format, "coorectVel", "Enable");
        }
/*...................................................................*/

/*... edm*/
        else if (!strcmp(word, combustion[9]))
        {
          cModel->reactionKinetic = EDM;
          setEdm(&cModel->edm,file);
          fprintf(fileLogExc, format, "EDM", "Enable");
        }
/*...................................................................*/
      }
/*...................................................................*/
      if (cModel->fLump) initLumpedMatrix(cModel);
      initEntalpyOfCombustion(cModel); 
    }
/*...................................................................*/
    readMacroV2(file, word, false, true);
  }while(strcmp(word,str));

}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 30/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------* 
 * readGravity : lendo o campo gravitacional                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * p       ->                                                        * 
 * file    -> arquivo de arquivo                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * p       ->                                                        * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void readGravity(DOUBLE *gravity,FILE *file){

  short i,n;
  char word[WORD_SIZE];

  readMacro(file,word,false);
  n = (short)atol(word);
  for(i=0;i<n;i++)
    fscanf(file,"%lf",gravity+i);

  fprintf(fileLogExc,"g = (%lf,%lf,%lf)\n"
         ,gravity[0],gravity[1],gravity[2]);

}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 17/07/2016                                   *
 * Data de modificaco : 11/09/2019                                   * 
 *-------------------------------------------------------------------* 
 * SETPPRINT : Seleciona as veriaves que serao impressas na          *
 * macro pFluid, puD1, puT1                                          *
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
void setPrint(FileOpt *opt,FILE *file){

  char str[]={"end"};
  char format[1024];
  char word[WORD_SIZE];
  char macro[][WORD_SIZE] = 
               {"cell"         ,"node"        ,"vel"             /* 0, 1, 2*/
               ,"pres"         ,"gradvel"     ,"gradpres"        /* 3, 4, 5*/
               ,"temp"         ,"gradtemp"    ,"eddyviscosity"   /* 6, 7, 8*/
               ,"densityfluid" ,"specificheat","dviscosity"      /* 9,10,11*/
               ,"tconductivity","vorticity"   ,"wallparameters"  /*12,13,14*/
               ,"stress"       ,"kinecit"     ,"stressr"         /*15,16,17*/
               ,"cdynamic"     ,"qcriterion"  ,"prestotal"       /*18,19,20*/
               ,"kturb"        ,"pkelvin"     ,"ud1"             /*21,22,23*/
               ,"gradud1"      ,"ut1"         ,"gradut1"         /*24,25,26*/
               ,"densityd1"    ,"coefdiffd1"  ,"densityt1"       /*27,28,29*/
               ,"coefdifft1"   ,"zcomb"       ,"gradzcomb"       /*30,31,32*/
               ,"qchemical"    ,"yfrac"       ,"rateheatcomb"    /*33,34,35*/                  
               ,"coefdiffsp"   ,"enthalpyk"   ,"grady"           /*36,37,38*/
               ,"treactor"     ,"binary"      ,"gradrho"         /*39,40,41*/
               ,""             ,""            ,""};              /*42,43,44*/
  int tmp,i=0,maxWord=100;

  strcpy(format,"%-20s: %s\n");

  initPrintVtk(opt);
  opt->bconditions   = true;

  fscanf(file,"%d",&tmp);
  opt->stepPlot[0] = opt->stepPlot[1] = (short) tmp;
  readMacroV2(file, word, false, true);
  while(strcmp(word,str) && i < maxWord)
  {
/*... cell*/        
    if(!strcmp(word,macro[0]))
    { 
      opt->fCell = true;
      fprintf(fileLogExc,format,"print","cell");
    }
/*.....................................................................*/

/*... node*/        
    else if(!strcmp(word,macro[1]))
    { 
      opt->fNode = true;
      fprintf(fileLogExc,format,"print","node");
    }
/*.....................................................................*/

/*... vel*/        
    else if(!strcmp(word,macro[2]))
    { 
      opt->vel = true;
      fprintf(fileLogExc,format,"print","vel");
    }
/*.....................................................................*/

/*...*/
    else if(!strcmp(word,macro[3]))
    { 
      opt->pres = true;
      fprintf(fileLogExc,format,"print","pres");
    }
/*.....................................................................*/

/*...*/
    else if(!strcmp(word,macro[4]))
    { 
      opt->gradVel = true;
      fprintf(fileLogExc,format,"print","gradVel");
    }
/*.....................................................................*/

/*...*/
    else if(!strcmp(word,macro[5]))
    { 
      opt->gradPres = true;
      fprintf(fileLogExc,format,"print","gradPres");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[6])) 
    {
      opt->temp = true;
      fprintf(fileLogExc,format,"print","temp");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[7])) 
    {
      opt->gradTemp = true;
      fprintf(fileLogExc,format,"print","gradTemp");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[8])) 
    {
      opt->eddyViscosity = true;
      fprintf(fileLogExc,format,"print","eddyViscosity");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[9])) 
    {
      opt->densityFluid = true;
      fprintf(fileLogExc,format,"print","densityFluid");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[10])) 
    {
      opt->specificHeat = true;
      fprintf(fileLogExc,format,"print","specificHeat");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[11])) 
    {
      opt->dViscosity = true;
      fprintf(fileLogExc,format,"print","dViscosity");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[12])) 
    {
      opt->tConductivity = true;
      fprintf(fileLogExc,format,"print","tConductivity");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[13])) 
    {
      opt->vorticity = true;
      fprintf(fileLogExc,format,"print","vorticity");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[14])) 
    {
      opt->wallParameters = true;
      fprintf(fileLogExc,format,"print","wallParameters");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[15])) 
    {
      opt->stress = true;
      fprintf(fileLogExc,format,"print","stress");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[16])) 
    {
      opt->kinetic = true;
      fprintf(fileLogExc,format,"print","kinecit");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[17])) 
    {
      opt->stressR = true;
      fprintf(fileLogExc,format,"print","stressR");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[18])) 
    {
      opt->cDynamic = true;
      fprintf(fileLogExc,format,"print","cDynamic");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[19])) 
    {
      opt->Qcriterion = true;
      fprintf(fileLogExc,format,"print","Qcriterion");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[20])) 
    {
      opt->presTotal = true;
      fprintf(fileLogExc,format,"print","presTotal");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[21])) 
    {
      opt->kTurb = true;
      fprintf(fileLogExc,format,"print","kTurb");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[22])) 
    {
      opt->pKelvin = true;
      fprintf(fileLogExc,format,"print","pKelvin");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[23])) 
    {
      opt->uD1 = true;
      fprintf(fileLogExc, format, "print", "uD1");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[24])) 
    {
      opt->graduD1 = true;
      fprintf(fileLogExc, format, "print", "graduD1");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[25]))  
    {
      opt->uT1 = true;
      fprintf(fileLogExc, format, "print", "uT1");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[26])) 
    {
      opt->graduT1 = true;
      fprintf(fileLogExc, format, "print", "graduT1");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[27]))
    {
      opt->densityD1 = true;
      fprintf(fileLogExc, format, "print", "densityD1");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[28]))
    {
      opt->coefDiffD1 = true;
      fprintf(fileLogExc, format, "print", "coefDiffD1");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[29]))
    {
      opt->densityT1 = true;
      fprintf(fileLogExc, format, "print", "densityT1");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[30]))
    {
      opt->coefDiffT1 = true;
      fprintf(fileLogExc, format, "print", "coefDiffT1");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[31]))
    {
      opt->zComb = true;
      fprintf(fileLogExc, format, "print", "zComb");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[32]))
    {
      opt->gradZcomb = true;
      fprintf(fileLogExc, format, "print", "gradZcomb");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[33]))
    {
      opt->wk = true;
      fprintf(fileLogExc, format, "print", "QChemical");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[34]))
    {
      opt->yFrac = true;
      fprintf(fileLogExc, format, "print", "yFrac");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[35]))
    {
      opt->rateHeatComb = true;
      fprintf(fileLogExc, format, "print"
                               , "rateHeatCombustion");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[36]))
    {
      opt->coefDiffSp = true;
      fprintf(fileLogExc, format, "print", "coefDiffSp");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[37]))
    {
      opt->enthalpyk = true;
      fprintf(fileLogExc, format, "print", macro[37]);
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[38]))
    {
      opt->gradY = true;
      fprintf(fileLogExc, format, "print", macro[38]);
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[39]))
    {
      opt->tReactor = true;
      fprintf(fileLogExc, format, "print", macro[39]);
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[40]))
    {
      opt->bVtk = true;
      fprintf(fileLogExc, format, "print", macro[40]);
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word, macro[41]))
    {
      opt->gradRho = true;
      fprintf(fileLogExc, format, "print", "gradRho");
    }
/*.....................................................................*/

    readMacroV2(file, word, false, true);
    i += 1;
  }
/*.....................................................................*/

/*...*/
  if (maxWord == i)
  {
    fprintf(fileLogDebug
          ,"Numero maximo de palavras lida na macro setPrint!!\n");
    fprintf(fileLogDebug
          ,"Possivel falta na palavra end!!\n");
    mpiStop();
    exit(EXIT_FAILURE);
  }
/*.....................................................................*/

} 
/*********************************************************************/ 

/*********************************************************************/
/* Leitura dos materiais                                             */
/*********************************************************************/
void uniformField(DOUBLE *field, INT const n, short const ndf
                ,FILE* file)
{

  DOUBLE value[MAXSPECIES];
  INT i;
  short j;
  
  for (j=0;j<ndf;j++)
    fscanf(file, "%lf", value+j);

  for (i = 0; i<n; i++) 
    for (j = 0; j<ndf; j++)
      MAT2D(i,j,field,ndf) = value[j];   
  
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 11/11/2017                                   *
 * Data de modificaco : 25/09/2019                                   *
 *-------------------------------------------------------------------*
 * help : Ajuda em relação a algumas macros                          *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * f - arquivo                                                       *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void help(FILE *f){

  char word[WORD_SIZE];
  short iHelp = 12;
  char help [][WORD_SIZE] = 
               {"macros"   ,"setprint" ,"advection"    /* 0, 1, 2*/
               ,"model"    ,"diffusion","nlit"         /* 3, 4, 5*/
               ,"transient","rcgrad"   ,"config"       /* 6, 7, 8*/
               ,"edp"      ,"openmp"   ,"propvar"      /* 9,10,11*/
               ,"mesh"     ,""         ,""    };       /*12,13,14*/

  short iMesh = 39;
  char mesh[][WORD_SIZE] ={
       "coordinates"  ,"endMesh"    ,"insert"         /* 0, 1, 2*/
       ,"return"       ,"cells"      ,"faceResT1"      /* 3, 4, 5*/
       ,"faceLoadT1"   ,"loadsT1"    ,"uniformT1"      /* 6, 7, 8*/ 
       ,"faceResZ"     ,"faceLoadZ"  ,"loadsZ"         /* 9,10,11*/ 
       ,"faceResD1"    ,"uniformD1"  ,"loadsD1"        /*12,13,14*/ 
       ,"faceLoadD1"   ,"initialD1"  ,""               /*15,16,17*/ 
       ,"faceResVel"   ,"loadsVel"   ,"faceRvel"    /*18,19,20*/ 
       ,"faceResPres"  ,"loadsPres"  ,"faceRpres"   /*21,22,23*/
       ,"faceResTemp"  ,"loadsTemp"  ,"faceLoadTemp"   /*24,25,26*/
       ,"materials"    ,"uniformPres","initialVel"     /*27,28,29*/
       ,"uniformTemp"  ,"uniformVel" ,"uniformZ"       /*30,31,32*/
       ,"faceReKturb"  ,"loadsKturb" ,"faceLoadKturb"  /*33,34,35*/
       ,"fileMaterials","initialTemp",""               /*36,37,38*/
       };

  short iMacros = 48;
  char macros[][WORD_SIZE] = 
       { ""            ,"mesh"         ,"stop"          /* 0, 1, 2*/
       ,"config"      ,"nextLoop"     ,"rcGrad"         /* 3, 4, 5*/
       ,"pgeo"        ,"pcoob"        ,"pcoo"           /* 6, 7, 8*/
       ,"setSolvDiff" ,"presolvT1"    ,"openmp"         /* 9,10,11*/
       ,"solvD1"      ,""             ,"pD1"            /*12,13,14*/
       ,"nlIt"        ,"pD1CsvCell"   ,"pD1CsvNode"     /*15,16,17*/
       ,"solvT1"      ,""             ,"pT1"            /*18,19,20*/
       ,""            ,"pT1CsvCell"   ,"pT1CsvNode"     /*21,22,23*/
       ,"setSolv"     ,"simple"       ,"setSimple"      /*24,25,26*/
       ,"transient"   ,"timeUpdate"   ,"partd"          /*27,28,29*/
       ,"advection"   ,"edp"          ,"diffusion"      /*30,31,32*/
       ,"pFluid"      ,"setPrint"     ,"reScaleMesh"    /*33,34,35*/
       ,"setPrime"    ,"prime"        ,"propVar"        /*36,37,38*/
       ,"setSolvComb" ,"pCombustion"  ,"simpleComb"     /*39,40,41*/
       ,"gravity"     ,"model"        ,"mean"           /*42,43,44*/
       ,"setMean"     ,"save"         ,"load"};         /*45,46,47*/
 
  short iPrint = 45;
  char print[][WORD_SIZE] = 
               {"cell"         ,"node"        ,"vel"             /* 0, 1, 2*/
               ,"pres"         ,"gradvel"     ,"gradpres"        /* 3, 4, 5*/
               ,"temp"         ,"gradtemp"    ,"eddyviscosity"   /* 6, 7, 8*/
               ,"densityfluid" ,"specificheat","dviscosity"      /* 9,10,11*/
               ,"tconductivity","vorticity"   ,"wallparameters"  /*12,13,14*/
               ,"stress"       ,"kinecit"     ,"stressr"         /*15,16,17*/
               ,"cdynamic"     ,"qcriterion"  ,"prestotal"       /*18,19,20*/
               ,"kturb"        ,"pkelvin"     ,"ud1"             /*21,22,23*/
               ,"gradud1"      ,"ut1"         ,"gradut1"         /*24,25,26*/
               ,"densityd1"    ,"coefdiffd1"  ,"densityt1"       /*27,28,29*/
               ,"coefdifft1"   ,"zcomb"       ,"gradzcomb"       /*30,31,32*/
               ,"ratefuel"     ,"yfrac"       ,"rateheatcomb"    /*33,34,35*/                  
               ,""             ,""            ,""                /*36,37,38*/
               ,""             ,""            ,""                /*39,40,41*/
               ,""             ,""            ,""};              /*42,43,44*/
;
/*... adveccao*/
  char fAdv[][WORD_SIZE] =                                   
                         { "FoUp","CD" ,"SoUp"                /* 0, 1, 2*/
                          ,"TVD" ,"NVD","LUST 0.25"};         /* 3, 4, 5*/  
  char tvd[][WORD_SIZE]=
                         {"VanLeer" ,"VanAlbada","MidMod "   /* 0, 1, 2*/
                         ,"Osher"   ,"SuperBee"};            /* 3, 4*/

  char scheme_nvd[][20] ={"BCD"        ,"MUSCL"   ,"Smart"       /* 0, 1, 2*/
                         ,"ModSmart"   ,"SuperBee","ModSuperBee" /* 3, 4, 5*/
                         ,"Stoic"      ,"MinMod"  ,"ModBCD"};    /* 5, 6, 7*/
/*....................................................................*/

/*... diffusion*/
  short iDiff = 4;
  char fDif[][WORD_SIZE]={"Orthogonal"  ,"Minimal","OrthogonalC" /* 0, 1, 2*/
                         ,"OverRelaxed"};                        /* 3*/
/*....................................................................*/

/*... model*/
  short iModels = 4;
  char models[][WORD_SIZE] = {"energy"  ,"turbulence","mass"     /* 0, 1, 2*/
                             ,"momentum","combustion"};          /* 3, 4*/
  short iEnergy = 7;  
  char energy[][WORD_SIZE] = { "preswork"   , "dissipation" /*0,1*/
                             , "residual"   , "absolute"    /*2,3*/
                             , "temperature", "entalphy"    /*4,5*/
                             , "diffenergy" };              /*6, */  

  short iTurb = 14;
  char turbulence[][WORD_SIZE] = {"wallmodel type"                       
                                 ,"smagorinsky 0.2"
                                 ,"wale 0.325"      
                                 ,"vreman 0.2" 
                                 ,"sigmamodel 0.136"     
                                 ,"ldynamic smagorinsky 0.2" 
                                 ,"gdynamic smagorinsky 0.2" 
                                 ,"gdynamicmod wale 0.325" 
                                 ,"bardina 1.0"
                                 ,"bardinaMod 1.0"
                                 ,"clark 1.0"
                                 ,"mixed 0.5 clark 1.0 smagorinsky 0.2"
                                 ,"towdynamic clark 1.0 smagorinsky 0.2"
                                 ,"oneeqk 0.094 1.048 1.0"};   
  short iMass = 2;
  char mass[][WORD_SIZE] = { "lhsDensity","rhsDensity"}; 

  short iMom = 8;  
  char momentum[][WORD_SIZE] =  {"residual"   ,"absolute"       /*0,1*/                    
                               ,"rhiechow"    ,"viscosity"      /*2,3*/
                               ,"div"         ,"buoyanthy"      /*4,5*/  
                               ,"buoyantprgh" ,"buoyantrhoref"  /*6,7*/
                               ,"sopressure"  };                /*6*/ 
  short iComb = 9;  
  char combustion[][WORD_SIZE] = {"residual"        ,"absolute"          /*0,1*/                    
                                 ,"grouped ns np nl","ungrouped ns np nl"/*2,3*/
                                 ,"edc options"     ,"arrhenius "        /*4,5*/ 
                                 ,"hcombustion"     ,"hformation"        /*6,7*/   
                                 ,"correctVel"};                         /*8*/
  short iEdc = 4;  
  char edc[][WORD_SIZE] = { "pan 1.01 1.0","panct 1.01 0.125"
                           ,"fluent 1.0"  ,"fluentct 0.125"}; 

  short iWall = 2;
  char typeWallModel[][WORD_SIZE] ={"standard","enhanced"};

  short iNlIt = 2;
  char sNlTy[][WORD_SIZE] = 
    { "nlIt 1 d1 500 1.e-06 50"
     ,"nlIt 2 d1 500 1.e-06 50 t1 500 1.e-06 50"};

  short iTrans = 2;
  char sTrans[][WORD_SIZE] = 
    { "transient config: BACKWARD 1.0e-02 1.e+01 dynamic"
     ,"transient config: EULER    1.0e-02 1.e+01 static" };


  short iRcGrad = 4;
  char sRcGrad[][WORD_SIZE] = 
    { "greenGaussCell"  ,"greenGaussNode"
     ,"leastSquare"     ,"leastSquareQR" };

  short iEdp= 8;
  char sEdp[][WORD_SIZE] =
    { "transport1" ,"transport2","transport3"   
     ,"diffusion1","diffusion2","diffusion3"    
     ,"fluid"     ,"fluidt"    ,"" };
/*....................................................................*/

  int i;                                                     

  printf("Help options:\n");
  readMacro(f,word,false);
  convStringLower(word);
/*... macros*/
  if(!strcmp(word,help[0])){
    printf("Macros:\n");
    for(i=0;i<iMacros;i++)
      printf("%3d - %s\n",i+1,macros[i]);
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... setPrintFluid*/        
  else if(!strcmp(word,help[1])){     
    printf("setPrint 10 options end\n");
    printf("10 time steps to print results\n");
    printf("options:\n");
    for(i=0;i<iPrint;i++)
      printf("%3d - %s\n",i+1,print[i]);
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... advection*/        
  else if(!strcmp(word,help[2])){   
    printf("Ex1:\n"); 
    printf("advection 3 Vel FoUp Energy CD Temp SoUp \n"); 
    printf("Ex2:\n"); 
    printf("advection 2 Vel NVD [options] Energy TVD [opitions]\n");
    printf("Options:\n");
    for(i=0;i<6;i++)
      printf("%3d - %s\n",i+1,fAdv[i]);
    printf("TVD options:\n");
    for(i=0;i<NFUNCLIMTFACE;i++)
      printf("%3d - %s\n",i+1,tvd[i]);
    printf("NVD options:\n");
    for(i=0;i<NFUNCNVD;i++)
      printf("%3d - %s\n",i+1,scheme_nvd[i]);
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... model*/        
  else if(!strcmp(word,help[3])){   
    printf("Ex1:\n"); 
    printf("model\nenergy 2 presWork dissipation\n"
           "turbulence 1 dynamic\nendModel\n\n"); 

    printf("Models options:\n");
    for(i=0;i<iModels;i++)
      printf("%3d - %s\n",i+1,models[i]);

    printf("Energy options:\n");
    for(i=0;i<iEnergy;i++)
      printf("%3d - %s\n",i+1,energy[i]);

    printf("Turbulence options:\n");
    for(i=0;i<iTurb;i++)
      printf("%3d - %s\n",i+1,turbulence[i]);

    printf("WallModel options:\n");
    for(i=0;i<iWall;i++)
      printf("%3d - %s\n",i+1,typeWallModel[i]);

    printf("Mass options:\n");
    for(i=0;i<iMass;i++)
      printf("%3d - %s\n",i+1,mass[i]);

    printf("Momentum options:\n");
    for(i=0;i<iMom;i++)
      printf("%3d - %s\n",i+1,momentum[i]);

    printf("Combustion options:\n");
    for(i=0;i<iComb;i++)
      printf("%3d - %s\n",i+1,combustion[i]);
    printf("EDC options:\n");  
    for (i = 0; i < iEdc; i++)
      printf("%3d - %s\n",i+1,edc[i]);
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... advection*/        
  else if(!strcmp(word,help[4])){   
    printf("Ex1:\n"); 
    printf("diffusion 2 Vel Orthogonal Pres Orthogonal \n"); 
    printf("Options:\n");
    for(i=0;i<iDiff;i++)
      printf("%3d - %s\n",i+1,fDif[i]);
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... nlIt*/
  else if (!strcmp(word, help[5])) {
    printf("Exemplos:\n");
    for (i = 0; i<iNlIt; i++)
      printf("%3d - %s\n", i + 1, sNlTy[i]);
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... transient*/
  else if (!strcmp(word, help[6])) {
    printf("Exemplos:\n");
    for (i = 0; i<iTrans; i++)
      printf("%3d - %s\n", i + 1, sTrans[i]);
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... rcGrad*/
  else if (!strcmp(word, help[7])) {
    printf("rcGrad |options|\n");
    printf("Exemplos:\n");
    printf("rcGrad greenGaussCell\n");
    printf("Options:\n");
    for (i = 0; i<iRcGrad; i++)
      printf("%3d - %s\n", i + 1, sRcGrad[i]);
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... config*/
  else if (!strcmp(word, help[8])) {
    printf("Ex:\n");
    printf("config\n");
    printf("reord false bvtk false memory 100 "
           "fItPlotRes false fItPlot true\n");
    printf("endConfig\n");
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... edp*/
  else if (!strcmp(word, help[9])) {
    printf("Ex:\n");
    printf("edp\nfluid 4\ndiffusion1 1\n");
    printf("endEdp\n");
    for (i = 0; i<iEdp; i++)
      printf("%3d - %s\n", i + 1, sEdp[i]);
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... openmp*/
  else if (!strcmp(word, help[10])) {
    printf("Ex:\n");
    printf("openmp 5 solver 2 update 2 cell 4 grad 2 reaction 2\n");
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... propVar*/
  else if (!strcmp(word, help[11])) {
    printf("Ex:\n");
    printf("PropVar\n");
    printf("Fluid\n");
    printf("density      polinomial mat/denPol.dat\n");
    printf("tcondutivity polinomial mat/ThCondPol.dat\n");
    printf("dviscosity   polinomial mat/dViscPol.dat\n");
    printf("sHeat        polinomial mat/sHeatPol.dat\n");
    printf("endFluid\n");
    printf("endPropVar\n");    
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... mesh*/
  else if (!strcmp(word, help[12])) {
    printf("Ex:\n");
    printf("Mesh options:\n");
    for (i = 0; i<iMesh; i++)
      printf("%3d - %s\n", i + 1, mesh[i]);
    exit(EXIT_HELP);
  }
/*.....................................................................*/

/*... opcoes*/        
  else{    
    printf("Options:\n");
    for(i=0;i<iHelp;i++)
      printf("%3d - %s\n",i+1,help[i]);
    exit(EXIT_HELP);
  }
/*.....................................................................*/

  exit(EXIT_SUCCESS);
}
/***********************************************************************/

/**********************************************************************
 * Data de criacao    : 09/12/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 *--------------------------------------------------------------------* 
 * setMixedModelLes : modelos mistos riedades variaveis               *                * 
 *--------------------------------------------------------------------* 
 * Parametros de entrada:                                             * 
 *--------------------------------------------------------------------* 
 * t       -> modelo de turbilencia                                   *
 * file    -> arquivo de arquivo                                      * 
 *--------------------------------------------------------------------* 
 * Parametros de saida:                                               * 
 *--------------------------------------------------------------------* 
 * t       -> atualizado                                              * 
 *--------------------------------------------------------------------* 
 * OBS:                                                               * 
 *--------------------------------------------------------------------*
 * Ex: mixed 1.0 bardina 1.0 wale 0.325                               *
 **********************************************************************/
void setMixedModelLes(Turbulence *t       , FILE *file) {

  char word[WORD_SIZE];
  char turb[][WORD_SIZE] = { "smagorinsky","wale"
                            ,"sigmamodel" ,"bardina" 
                            ,"clark"      ,"bardinamod" }; 

  short k = 0, ii=0;

  do{
    ii++;
    readMacro(file,word,false);
    convStringLower(word);
/*... Smagorinsky*/
    if(!strcmp(word,turb[0])){ 
      k++;
      t->typeMixed[FUNMODEL] = SMAGORINSKY;
      fscanf(file,"%lf",&t->cf);  
      fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[0],t->cf);    
    }
/*...................................................................*/ 

/*... Wale*/
    else if(!strcmp(word,turb[1])){ 
      k++;
      t->typeMixed[FUNMODEL] = WALEMODEL;
      fscanf(file,"%lf",&t->cf);        
      fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[1],t->cf);    
    }
/*...................................................................*/ 

/*... sigmaModel*/
    else if(!strcmp(word,turb[2])){
      k++; 
      t->typeMixed[FUNMODEL]  = SIGMAMODEL;
      fscanf(file,"%lf",&t->cf);   
      fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[2],t->cf);
    }
/*...................................................................*/ 

/*... bardina*/
    else if(!strcmp(word,turb[3])){
      k++;   
      t->typeMixed[ESTMODEL] = BARDINA;
      fscanf(file,"%lf",&t->cs);  
      fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[3],t->cs);
    }
/*...................................................................*/  

/*... clark*/
    else if(!strcmp(word,turb[4])){
      k++;  
      t->typeMixed[ESTMODEL] = CLARK;
      fscanf(file,"%lf",&t->cs);  
      fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[4],t->cs);
    } 
/*...................................................................*/ 

/*... bardinaMod*/
    else if(!strcmp(word,turb[5])){
      k++;   
      t->typeMixed[ESTMODEL] = BARDINAMOD;
      fscanf(file,"%lf",&t->cs);  
      fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[5],t->cs);
    }
/*...................................................................*/ 
  
  }while( k!=2 && ii != 10);


  if (ii == 10)
    ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__
              ,"Erro na leitura",EXIT_PROG);

}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 12/12/2017                                    *
 * Data de modificaco : 13/01/2018                                    *
 *--------------------------------------------------------------------* 
 * setDynamicModelLes : modelos dinamicos com um parementro           *
 *--------------------------------------------------------------------* 
 * Parametros de entrada:                                             * 
 *--------------------------------------------------------------------* 
 * t       -> modelo de turbilencia                                   *
 * file    -> arquivo de arquivo                                      * 
 *--------------------------------------------------------------------* 
 * Parametros de saida:                                               * 
 *--------------------------------------------------------------------* 
 * t       -> atualizado                                              * 
 *--------------------------------------------------------------------* 
 * OBS:                                                               * 
 *--------------------------------------------------------------------*
 **********************************************************************/
void setDynamicModelLes(Turbulence *t       , FILE *file) {

  char word[WORD_SIZE];
  char turb[][WORD_SIZE] = { "smagorinsky","sigmamodel"
                            ,"wale"       ,"vreman"}; 

  short k = 0;

  readMacro(file,word,false);
  convStringLower(word);
/*... Smagorinsky*/
  if(!strcmp(word,turb[0])){ 
    k++;
    t->typeLes             = LESFUNCMODEL;
    t->typeMixed[FUNMODEL] = SMAGORINSKY;
    fscanf(file,"%lf",&t->cf);  
    fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[0],t->cf);    
  }
/*...................................................................*/ 

/*... SigmaModel*/
  else if(!strcmp(word,turb[1])){ 
    k++;
    t->typeLes             = LESFUNCMODEL;
    t->typeMixed[FUNMODEL] = SIGMAMODEL;
    fscanf(file,"%lf",&t->cf);  
    fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[1],t->cf);    
  }
/*...................................................................*/

/*... WaleModel*/
  else if(!strcmp(word,turb[2])){ 
    k++;
    t->typeLes             = LESFUNCMODEL;
    t->typeMixed[FUNMODEL] = WALEMODEL;
    fscanf(file,"%lf",&t->cf);  
    fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[2],t->cf);    
  }
/*...................................................................*/ 

/*... VREMAN*/
  else if(!strcmp(word,turb[3])){ 
    k++;
    t->typeLes             = LESFUNCMODEL;
    t->typeMixed[FUNMODEL] = VREMAN;
    fscanf(file,"%lf",&t->cf);  
    fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[3],t->cf);    
  }
/*...................................................................*/ 
}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 26/05/2019                                    *
 * Data de modificaco : 00/00/0000                                    *
 *--------------------------------------------------------------------* 
 * setEdc : Eddy Disspantion concept                                  *                * 
 *--------------------------------------------------------------------* 
 * Parametros de entrada:                                             * 
 *--------------------------------------------------------------------* 
 * e       -> modelo de EDC      ia                                   *
 * file    -> arquivo de arquivo                                      * 
 *--------------------------------------------------------------------* 
 * Parametros de saida:                                               * 
 *--------------------------------------------------------------------* 
 * e       -> atualizado                                              * 
 *--------------------------------------------------------------------* 
 * OBS:                                                               * 
 *--------------------------------------------------------------------*
 **********************************************************************/
void setEdc(Edc *e       , FILE *file) 
{

  char word[WORD_SIZE];
  char edc[][WORD_SIZE] = { "pan"        ,"panct"
                           ,"fluent"     ,"fluentct"
                           ,"fds"                   }; 

  readMacro(file,word,false);
  convStringLower(word);
/*... PANJWANIN (2010)*/
  if(!strcmp(word,edc[0]))
  { 
    e->type = PANJWANI_EDC;
    fscanf(file,"%lf %lf",&e->cGamma,&e->cTau);  
    fprintf(fileLogExc,"%-20s: cGamma = %lf  cTau = %lf\n"
                        , edc[0],e->cGamma,e->cTau);    
  }
/*...................................................................*/ 

/*... PANJWANIN (2010) com tempo de mistura constante*/
  else if(!strcmp(word,edc[1]))
  { 
    e->type = PANJWANI_CONST_TMIX_EDC;
    fscanf(file,"%lf %lf",&e->cGamma,&e->tMix);  
    fprintf(fileLogExc,"%-20s: cGamma = %lf  tMix = %lf\n"
                        , edc[1],e->cGamma,e->tMix);      
  }
/*...................................................................*/ 

/*... Fluent*/
  else if(!strcmp(word,edc[2]))
  {
    e->type = FLUENT_EDC;
    fscanf(file,"%lf",&e->cTau);  
    fprintf(fileLogExc,"%-20s:  cTau = %lf\n"
                        , edc[2], e->cTau);
  }   
/*...................................................................*/ 

/*... Fluent com tempo de mistura constante*/
  else if(!strcmp(word,edc[3])) 
  {
    e->type = FLUENT_CONST_TMIX_EDC;
    fscanf(file,"%lf",&e->tMix);  
    fprintf(fileLogExc,"%-20s: tMix = %lf\n"
                        , edc[3],e->tMix);   
  }
/*...................................................................*/ 

/*... Fluent com tempo de mistura constante*/
  else if(!strcmp(word,edc[4])) 
  {
    e->type = FDS_EDC;
    fscanf(file,"%lf %lf",&e->cGamma,&e->cTau);   
    fprintf(fileLogExc,"%-20s: tc = %lf tf = %lf\n"
                        , edc[3],e->cGamma,e->cTau);   
  }
/*...................................................................*/ 

}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 26/05/2019                                    *
 * Data de modificaco : 26/07/2019                                    *
 *--------------------------------------------------------------------* 
 * setEdm : Eddy Disspantion model                                    *                * 
 *--------------------------------------------------------------------* 
 * Parametros de entrada:                                             * 
 *--------------------------------------------------------------------* 
 * e       -> modelo de EDm                                           *
 * file    -> arquivo de arquivo                                      * 
 *--------------------------------------------------------------------* 
 * Parametros de saida:                                               * 
 *--------------------------------------------------------------------* 
 * e       -> atualizado                                              * 
 *--------------------------------------------------------------------* 
 * OBS:                                                               * 
 *--------------------------------------------------------------------*
 **********************************************************************/
void setEdm(Edm *e       , FILE *file) 
{

  char word[WORD_SIZE];
  char edm[][WORD_SIZE] = { "a"          ,"b"
                           ,"tmix"};        
  short i;
  int nTerms;
  e->tMixConst = false;
  e->fProd     = false;
  fscanf(file,"%d",&nTerms);
  for(i=0;i<nTerms;i++)
  {
    readMacro(file,word,false);
    convStringLower(word);
/*...*/
    if(!strcmp(word,edm[0]))
      fscanf(file,"%lf",&e->coef[0]);
/*...................................................................*/ 

/*...*/
    else if(!strcmp(word,edm[1]))
      fscanf(file,"%lf",&e->coef[1]);     
/*...................................................................*/ 

/*...*/
    else if(!strcmp(word,edm[2]))
    {
      e->tMixConst = true;
      fscanf(file,"%lf",&e->coef[2]);
    }
/*...................................................................*/ 
  }
/*...................................................................*/

/*...*/
   if(e->tMixConst)
    fprintf(fileLogExc,"%-20s: A = %lf  B = %lf tMix = %lf\n"
                        , "edm",e->coef[0],e->coef[1],e->coef[2]);
   else
    fprintf(fileLogExc,"%-20s: A = %lf  B = %lf\n"
                        , "edm",e->coef[0],e->coef[1]);
/*...................................................................*/
}
/**********************************************************************/


/**********************************************************************
 * Data de criacao    : 28/01/2018                                    *
 * Data de modificaco : 09/05/2019                                    *
 *--------------------------------------------------------------------* 
 * readAdvectionScheme:                                               *
 *--------------------------------------------------------------------* 
 * Parametros de entrada:                                             * 
 *--------------------------------------------------------------------* 
 * file    -> arquivo de arquivo                                      * 
 *--------------------------------------------------------------------* 
 * Parametros de saida:                                               * 
 *--------------------------------------------------------------------* 
 *--------------------------------------------------------------------* 
 * OBS:                                                               * 
 *--------------------------------------------------------------------*
 **********************************************************************/
void readAdvectionScheme(FILE *fileIn, Scheme *sc) {

  char word[WORD_SIZE];
  unsigned short nScheme;
/*... tecnica de adveccao*/
  readMacro(fileIn, word, false);
  nScheme = (short) atol(word);
  do {
    readMacro(fileIn, word, false);
/*... velocidade*/
    if (!strcmp(word, "Vel") || !strcmp(word, "vel")) {
      fprintf(fileLogExc,"%s:\n", word);
      readMacro(fileIn, word, false);
 /*... codigo da da funcao limitadora de fluxo*/
      setAdvectionScheme(word, &sc->advVel,fileIn);
      nScheme--;
    }
 /*... T1*/
    else if (!strcmp(word, "T1") || !strcmp(word, "t1")) {
      fprintf(fileLogExc,"%s:\n", word);
      readMacro(fileIn, word, false);
 /*... codigo da da funcao limitadora de fluxo*/
      setAdvectionScheme(word, &sc->advT1,fileIn);
      nScheme--;
    }
 /*... Energy*/
    else if (!strcmp(word, "Energy") || !strcmp(word, "energy")) {
      fprintf(fileLogExc,"%s:\n", word);
      readMacro(fileIn, word, false);
/*... codigo da da funcao limitadora de fluxo*/
      setAdvectionScheme(word, &sc->advEnergy, fileIn);
      nScheme--;
    }
/*... kTutb*/
    else if (!strcmp(word, "kTurb") || !strcmp(word, "kturb")) {
      fprintf(fileLogExc,"%s:\n", word);
      readMacro(fileIn, word, false);
/*... codigo da da funcao limitadora de fluxo*/
      setAdvectionScheme(word, &sc->advKturb, fileIn);
      nScheme--;
    }
/*... Zcomb*/
    else if (!strcmp(word, "Zcomb") || !strcmp(word, "zcomb")) {
      fprintf(fileLogExc,"%s:\n", word);
      readMacro(fileIn, word, false);
 /*... codigo da da funcao limitadora de fluxo*/
      setAdvectionScheme(word, &sc->advComb, fileIn);
      nScheme--;
    }
  } while (nScheme);
/*...................................................................*/
}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 28/01/2018                                    *
 * Data de modificaco : 09/05/2019                                    *
 *--------------------------------------------------------------------* 
 * readDiffusionScheme:                                               *
 *--------------------------------------------------------------------* 
 * Parametros de entrada:                                             * 
 *--------------------------------------------------------------------* 
 * file    -> arquivo de arquivo                                      * 
 *--------------------------------------------------------------------* 
 * Parametros de saida:                                               * 
 *--------------------------------------------------------------------* 
 *--------------------------------------------------------------------* 
 * OBS:                                                               * 
 *--------------------------------------------------------------------*
 **********************************************************************/
void readDiffusionScheme(FILE *fileIn, Scheme *sc) {

  char word[WORD_SIZE];
  unsigned short nScheme;

/*... tecnica de adveccao*/
  readMacro(fileIn, word, false);
  nScheme = (short) atol(word);
  do {
    readMacro(fileIn,word,false);
/*... velocidade*/
    if(!strcmp(word,"Vel") || !strcmp(word, "vel")){
      fprintf(fileLogExc,"%s:\n",word);
      readMacro(fileIn,word,false);
/*... codigo da da funcao limitadora de fluxo*/        
      setDiffusionScheme(word,&sc->diffVel.iCod);
      nScheme--;
    }
/*... Pressao*/
    else if(!strcmp(word,"Pres") || !strcmp(word, "pres")){
      fprintf(fileLogExc,"%s:\n",word);
      readMacro(fileIn,word,false);
/*... codigo da da funcao limitadora de fluxo*/        
      setDiffusionScheme(word,&sc->diffPres.iCod);
      nScheme--;
    }
 /*... transporte T1*/
    else if (!strcmp(word,"T1") || !strcmp(word, "t1")){
      fprintf(fileLogExc,"%s:\n", word);
      readMacro(fileIn, word, false);
 /*... codigo da da funcao limitadora de fluxo*/
      setDiffusionScheme(word, &sc->diffT1.iCod);
      nScheme--;
    }
/*... difusao D1*/
    else if (!strcmp(word, "D1") || !strcmp(word, "D1")) {
      fprintf(fileLogExc, "%s:\n", word);
      readMacro(fileIn, word, false);
/*... codigo da da funcao limitadora de fluxo*/
      setDiffusionScheme(word, &sc->diffD1.iCod);
      nScheme--;
    }
 /*... Energy*/
    else if (!strcmp(word, "Energy") || !strcmp(word, "energy")) {
      fprintf(fileLogExc,"%s:\n", word);
      readMacro(fileIn, word, false);
 /*... codigo da da funcao limitadora de fluxo*/
      setDiffusionScheme(word, &sc->diffEnergy.iCod);
      nScheme--;
    }
/*... Zcomb*/
    else if (!strcmp(word, "Zcomb") || !strcmp(word, "zcomb")) {
      fprintf(fileLogExc,"%s:\n", word);
      readMacro(fileIn, word, false);
 /*... codigo da da funcao limitadora de fluxo*/
      setDiffusionScheme(word, &sc->diffComb.iCod);
      nScheme--;
    }
  } while (nScheme);
/*...................................................................*/

}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 28/01/2018                                    *
 * Data de modificaco : 31/03/2018                                    *
 *--------------------------------------------------------------------* 
 * readSetSimple:                                                     *
 *--------------------------------------------------------------------* 
 * Parametros de entrada:                                             * 
 *--------------------------------------------------------------------* 
 * file    -> arquivo de arquivo                                      * 
 *--------------------------------------------------------------------* 
 * Parametros de saida:                                               * 
 *--------------------------------------------------------------------* 
 *--------------------------------------------------------------------* 
 * OBS:                                                               * 
 *--------------------------------------------------------------------*
 **********************************************************************/
void readSetSimple(Memoria *m    , FILE *fileIn
                 , Mesh *mesh0   , Mesh *mesh
                 , Simple *simple, bool *fSolvSimple) {

  char word[WORD_SIZE];

/*...*/
  *fSolvSimple            = true;  
  simple->maxIt           = 1000;
  simple->alphaPres       = 0.3e0; 
  simple->alphaVel        = 0.7e0; 
  simple->type            = SIMPLE;
  simple->sPressure       = true;
  simple->faceInterpolVel = 1;
  simple->nNonOrth        = 0;
  if (mesh->ndfFt){
    simple->alphaEnergy  = 1.e0;
    simple->alphaDensity = 1.0e0; 
    simple->alphaComb    = 1.0e0;
  }
  simple->pSimple         = 500;
/*...................................................................*/
      
/*...*/
  readMacro(fileIn,word,false);
  if(!strcmp(word,"config:")){
/*...*/        
    readMacro(fileIn,word,false);
    setSimpleScheme(word, simple);
/*...................................................................*/    
  }
/*...................................................................*/

/*...*/
  HccaAlloc(DOUBLE     ,m       ,simple->d
           ,mesh->numel*mesh->ndm,"dField" ,false);
  zero(simple->d    ,mesh->numel*mesh->ndm,DOUBLEC);

  HccaAlloc(DOUBLE     ,m       ,simple->ePresC
           ,mesh->numel,"ePresC" ,false);
  zero(simple->ePresC,mesh->numel  ,DOUBLEC);

  HccaAlloc(DOUBLE     ,m       ,simple->nPresC
           ,mesh->nnode,"nPresC" ,false);
  zero(simple->nPresC    ,mesh->numel  ,DOUBLEC);

  HccaAlloc(DOUBLE     ,m      ,simple->eGradPresC
           ,mesh->numel*mesh->ndm,"eGradPresC",false);
  zero(simple->eGradPresC,mesh->numel*mesh->ndm  ,DOUBLEC);
  
  HccaAlloc(DOUBLE     ,m       ,simple->ePresC1
           ,mesh->numel,"ePresC1",false);
  zero(simple->ePresC,mesh->numel  ,DOUBLEC);
/*...................................................................*/
}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 08/05/2019                                    *
 * Data de modificaco : 00/00/0000                                    *
 *--------------------------------------------------------------------* 
 * readSetSimpleComb:                                                 *
 *--------------------------------------------------------------------* 
 * Parametros de entrada:                                             * 
 *--------------------------------------------------------------------* 
 * file    -> arquivo de arquivo                                      * 
 *--------------------------------------------------------------------* 
 * Parametros de saida:                                               * 
 *--------------------------------------------------------------------* 
 *--------------------------------------------------------------------* 
 * OBS:                                                               * 
 *--------------------------------------------------------------------*
 **********************************************************************/
void readSetSimpleComb(Memoria *m    , FILE *fileIn
                     , Mesh *mesh0   , Mesh *mesh
                    , Simple *simple, bool *fSolvComb) {

  char word[WORD_SIZE];

/*...*/
  *fSolvComb              = true;  
  simple->maxIt           = 1000;
  simple->alphaPres       = 0.3e0; 
  simple->alphaVel        = 0.7e0; 
  simple->type            = SIMPLE;
  simple->sPressure       = true;
  simple->faceInterpolVel = 1;
  simple->nNonOrth        = 0;
  simple->alphaEnergy     = 1.e0;
  simple->alphaDensity    = 1.0e0; 
  simple->alphaComb       = 1.0e0;

  simple->pSimple         = 500;
/*...................................................................*/
      
/*...*/
  readMacro(fileIn,word,false);
  if(!strcmp(word,"config:")){
/*... timer*/        
    readMacro(fileIn,word,false);
/*... levemente compressivel*/       
    setSimpleCombustionScheme(word,mesh0->ndm,simple,fileIn);
/*...*/        
    if(simple->type == SIMPLE && !mpiVar.myId)     
      fprintf(fileLogExc,"PRES-VEL  : SIMPLE\n");
    else if(simple->type == SIMPLEC && !mpiVar.myId )     
      fprintf(fileLogExc,"PRES-VEL  : SIMPLEC\n");

/*...*/        
    fprintf(fileLogExc,"%-15s : %d\n","Maxit"       ,simple->maxIt);
    fprintf(fileLogExc,"%-15s : %lf\n","alphaPres"  ,simple->alphaPres);
    fprintf(fileLogExc,"%-15s : %lf\n","alphaVel"   ,simple->alphaVel);
    fprintf(fileLogExc,"%-15s : %lf\n","alphaEnergy",simple->alphaEnergy);
    fprintf(fileLogExc,"%-15s : %lf\n","alphaComb"  ,simple->alphaComb);
    fprintf(fileLogExc,"%-15s : %d\n","nNonOrth"    ,simple->nNonOrth);
    fprintf(fileLogExc,"%-15s : %d\n","pSimple"     ,simple->pSimple);
  }
/*...................................................................*/

/*...*/
  HccaAlloc(DOUBLE     ,m       ,simple->d
           ,mesh->numel*mesh->ndm,"dField" ,false);
  zero(simple->d    ,mesh->numel*mesh->ndm,DOUBLEC);

  HccaAlloc(DOUBLE     ,m       ,simple->ePresC
           ,mesh->numel,"ePresC" ,false);
  zero(simple->ePresC,mesh->numel  ,DOUBLEC);

  HccaAlloc(DOUBLE     ,m       ,simple->nPresC
           ,mesh->nnode,"nPresC" ,false);
  zero(simple->nPresC    ,mesh->numel  ,DOUBLEC);

  HccaAlloc(DOUBLE     ,m      ,simple->eGradPresC
           ,mesh->numel*mesh->ndm,"eGradPresC",false);
  zero(simple->eGradPresC,mesh->numel*mesh->ndm  ,DOUBLEC);
  
  HccaAlloc(DOUBLE     ,m       ,simple->ePresC1
           ,mesh->numel,"ePresC1",false);
  zero(simple->ePresC,mesh->numel  ,DOUBLEC);
/*...................................................................*/
}
/**********************************************************************/

/**********************************************************************
* Data de criacao    : 31/03/2018                                    *
* Data de modificaco : 00/00/0000                                    *
*--------------------------------------------------------------------*
* readSetSimple:                                                     *
*--------------------------------------------------------------------*
* Parametros de entrada:                                             *
*--------------------------------------------------------------------*
*--------------------------------------------------------------------*
* Parametros de saida:                                               *
*--------------------------------------------------------------------*
*--------------------------------------------------------------------*
* OBS:                                                               *
*--------------------------------------------------------------------*
**********************************************************************/
void readSetPrime(Memoria *m    , FILE *fileIn
                 , Mesh *mesh0  , Mesh *mesh
                 , Prime  *prime, bool *fSolvPrime) {

  char word[WORD_SIZE];

  *fSolvPrime = true;
  prime->maxIt = 2000;
  prime->alphaPres = 1.0e0;
  prime->alphaVel = 0.9e0;
  prime->kZeroVel = 1;
  prime->kZeroPres = 0;
  prime->sPressure = true;
  prime->nNonOrth = 0;
  prime->tolPres = 1.e-04;
  prime->tolVel = 1.e-04;
  prime->pPrime = 10;
/*...................................................................*/

/*...*/
  readMacro(fileIn, word, false);
  if (!strcmp(word, "config:"))
    setPrimeScheme(word, prime, fileIn);
/*...................................................................*/

/*...*/
  HccaAlloc(DOUBLE, m, prime->d
    , mesh->numel*mesh->ndm, "dField", false);
  zero(prime->d, mesh->numel*mesh->ndm, DOUBLEC);

  HccaAlloc(DOUBLE, m, prime->velUp
    , mesh->numel*mesh->ndm, "velUp", false);
  zero(prime->d, mesh->numel*mesh->ndm, DOUBLEC);

  HccaAlloc(DOUBLE, m, prime->ePresC
    , mesh->numel, "ePresC", false);
  zero(prime->ePresC, mesh->numel, DOUBLEC);

  HccaAlloc(DOUBLE, m, prime->nPresC
    , mesh->nnode, "nPresC", false);
  zero(prime->nPresC, mesh->numel, DOUBLEC);

  HccaAlloc(DOUBLE, m, prime->eGradPresC
    , mesh->numel*mesh->ndm, "eGradPresC", false);
  zero(prime->eGradPresC, mesh->numel*mesh->ndm, DOUBLEC);

  HccaAlloc(DOUBLE, m, prime->ePresC1
    , mesh->numel, "ePresC1", false);
  zero(prime->ePresC, mesh->numel, DOUBLEC);

  HccaAlloc(DOUBLE, m, prime->aD
    , mesh->numel*mesh->ndm, "aD", false);
  zero(prime->aD, mesh->numel*mesh->ndm, DOUBLEC);

  HccaAlloc(DOUBLE, m, prime->bTemporal
    , mesh->numel*mesh->ndm, "bT", false);
  zero(prime->bTemporal, mesh->numel*mesh->ndm, DOUBLEC);
  /*...................................................................*/
}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 31/03/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * readSolvFluid: leitura das configuracoes do solvers utilizados    *
 * no escoamento de fluidos                                          *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void readSolvFluid(Memoria *m      , Mesh *mesh          , Reord *reordMesh                 
                 , Solv *solvVel   , SistEq* sistEqVel   , bool *fSolvVel
                 , Solv *solvPres  , SistEq* sistEqPres  , bool *fSolvPres
                 , Solv *solvEnergy, SistEq* sistEqEnergy, bool *fSolvEnergy
                 , Solv *solvKturb , SistEq* sistEqKturb , bool *fSolvKturb                 
                 , char* auxName   , char* preName       , char* nameOut
                 , FILE *fileIn    , FileOpt *opt)
{

  unsigned short ndfVel, nSistEq;

  INT nEqMax;
/*... Estrutura de dados*/
  char strIa[MNOMEPONTEIRO], strJa[MNOMEPONTEIRO];
  char strA[MNOMEPONTEIRO], strAd[MNOMEPONTEIRO];

  char str1[100], str2[100], str3[100], str4[100];

  char word[WORD_SIZE];

/*... solver*/
  readMacro(fileIn, word, false);
  nSistEq = (short)atol(word);
  do 
  {
    readMacro(fileIn, word, false);
/*... velocidade*/
    if (!strcmp(word, "Vel") || !strcmp(word, "vel")) 
    {
      fprintf(fileLogExc, "%s.\n",word);
      nSistEq--;
      *fSolvVel = true;
      solvVel->solver = PBICGSTAB;
      solvVel->tol = smachn();
      solvVel->maxIt = 50000;
      solvVel->fileSolv = NULL;
      solvVel->log = true;
      solvVel->flag = true;
/*...................................................................*/

/*...*/
      if (solvVel->log && !mpiVar.myId) 
      {
        strcpy(auxName, preName);
        strcat(auxName, "_fluid_vel");
        fName(auxName, mpiVar.nPrcs, 0, 11, nameOut);
        solvVel->fileSolv = openFile(nameOut, "w");
      }
/*...................................................................*/

      sistEqVel->unsym = true;
/*...................................................................*/

/*... solver*/
      readMacro(fileIn, word, false);
      setSolverConfig(word, solvVel, &sistEqVel->storage);
/*...................................................................*/

/*... numeracao das equacoes das velocidades*/
      HccaAlloc(INT, m, sistEqVel->id
        , mesh->numel
        , "sistVelid", _AD_);

      fprintf(fileLogExc, "%s\n", DIF);
      fprintf(fileLogExc, "Numerando as equacoes.\n");

      tm.numeqVel = getTimeC() - tm.numeqVel;
      sistEqVel->neq = numEqV1(sistEqVel->id, reordMesh->num
                             , mesh->numel);
      tm.numeqVel = getTimeC() - tm.numeqVel;

      fprintf(fileLogExc, "Equacoes numeradas.\n");
      fprintf(fileLogExc, "%s\n", DIF);
/*...................................................................*/

/*...*/
      if (mpiVar.nPrcs > 1) 
      {
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
      HccaAlloc(DOUBLE, m, sistEqVel->b0
        , sistEqVel->neq*ndfVel, "sistVelb0", _AD_);
      HccaAlloc(DOUBLE, m, sistEqVel->b
        , sistEqVel->neq*ndfVel, "sistVelb ", _AD_);
      HccaAlloc(DOUBLE, m, sistEqVel->x
        , sistEqVel->neq*ndfVel, "sistVelx ", _AD_);
      zero(sistEqVel->b0, sistEqVel->neq*ndfVel, DOUBLEC);
      zero(sistEqVel->b, sistEqVel->neq*ndfVel, DOUBLEC);
      zero(sistEqVel->x, sistEqVel->neq*ndfVel, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados velocidades*/
      strcpy(strIa, "iaVel");
      strcpy(strJa, "jaVel");
      strcpy(strAd, "adVel");
      strcpy(strA, "aVel");

      fprintf(fileLogExc, "Vel:\n");
      fprintf(fileLogExc, "Montagem da estrura de dados esparsa.\n");

      tm.dataStructVel = getTimeC() - tm.dataStructVel;
      dataStructSimple(m, sistEqVel->id, reordMesh->num
        , mesh->elm.adj.nelcon
        , mesh->elm.adj.nViz, mesh->numelNov, mesh->maxViz
        , ndfVel, strIa, strJa
        , strAd, strA, sistEqVel);
      tm.dataStructVel = getTimeC() - tm.dataStructVel;

      fprintf(fileLogExc, "Estrutuda montada.\n");
/*...................................................................*/

/*... dividindo a matriz*/
/*... Openmp(Vel,Pres)*/
      if (ompVar.fSolver) 
      {
        strcpy(str1, "thBeginVel");
        strcpy(str2, "thEndVel");
        strcpy(str3, "thSizeVel");
        strcpy(str4, "thHeightVel");
        pMatrixSolverOmp(m, sistEqVel, str1, str2, str3, str4);
      }
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    else if (!strcmp(word, "Pres") || !strcmp(word, "pres")) 
    {
      fprintf(fileLogExc, "%s.\n",word);
      nSistEq--;
      *fSolvPres = true;
      solvPres->solver = PCG;
      solvPres->tol = smachn();
      solvPres->maxIt = 50000;
      solvPres->fileSolv = NULL;
      solvPres->log = true;
      solvPres->flag = true;
/*...................................................................*/

/*...*/
      if (solvPres->log && !mpiVar.myId) 
      {
        strcpy(auxName, preName);
        strcat(auxName, "_fluid_pres");
        fName(auxName, mpiVar.nPrcs, 0, 11, nameOut);
        solvPres->fileSolv = openFile(nameOut, "w");
      }
/*...................................................................*/

/*... inicializa a estrutura do solver(PRESSAO)*/
      sistEqPres->unsym = false;
/*...................................................................*/

/*... solver*/
      readMacro(fileIn, word, false);
      setSolverConfig(word, solvPres, &sistEqPres->storage);
/*...................................................................*/

/*... numeracao das equacoes das pressoes*/
      HccaAlloc(INT, m, sistEqPres->id
        , mesh->numel
        , "sistPresid", _AD_);

      fprintf(fileLogExc, "%s\n", DIF);
      fprintf(fileLogExc, "Numerando as equacoes.\n");

      tm.numeqPres = getTimeC() - tm.numeqPres;
      sistEqPres->neq = numEqV2(sistEqPres->id, reordMesh->num
        , mesh->elm.faceRpres, mesh->elm.adj.nViz
        , mesh->numel, mesh->maxViz);
      tm.numeqPres = getTimeC() - tm.numeqPres;
      fprintf(fileLogExc, "Equacoes numeradas.\n");
      fprintf(fileLogExc, "%s\n", DIF);
/*...................................................................*/

/*...*/
      if (mpiVar.nPrcs > 1) 
      {
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
      HccaAlloc(DOUBLE, m, sistEqPres->b0
        , sistEqPres->neq, "sistPresb0", _AD_);
      HccaAlloc(DOUBLE, m, sistEqPres->b
        , sistEqPres->neq, "sistPresb ", _AD_);
      HccaAlloc(DOUBLE, m, sistEqPres->x
        , sistEqPres->neq, "sistPresx ", _AD_);
      zero(sistEqPres->b0, sistEqPres->neq, DOUBLEC);
      zero(sistEqPres->b, sistEqPres->neq, DOUBLEC);
      zero(sistEqPres->x, sistEqPres->neq, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados pressoes*/
      strcpy(strIa, "iaPres");
      strcpy(strJa, "japres");
      strcpy(strAd, "adPres");
      strcpy(strA, "aPres");

      fprintf(fileLogExc, "Pres:\n");
      fprintf(fileLogExc, "Montagem da estrura de dados esparsa.\n");

      tm.dataStructPres = getTimeC() - tm.dataStructPres;
      dataStruct(m, sistEqPres->id, reordMesh->num, mesh->elm.adj.nelcon
        , mesh->elm.adj.nViz, mesh->numelNov, mesh->maxViz
        , 1, strIa, strJa
        , strAd, strA, sistEqPres);
      tm.dataStructPres = getTimeC() - tm.dataStructPres;
      fprintf(fileLogExc, "Estrutuda montada.\n");
/*...................................................................*/

/*... Openmp(Vel,Pres)*/
      if (ompVar.fSolver) 
      {
/*... dividindo a matriz*/
        strcpy(str1, "thBeginPres");
        strcpy(str2, "thEndPres");
        strcpy(str3, "thSizePres");
        strcpy(str4, "thHeightPres");
        pMatrixSolverOmp(m, sistEqPres, str1, str2, str3, str4);
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... energy*/
    else if (!strcmp(word, "Energy") || !strcmp(word, "energy")) 
    {
      fprintf(fileLogExc, "%s.\n",word);
      nSistEq--;
/*... inicializando a estrutura de equacoes do problema (ENERGY)*/
      *fSolvEnergy = true;
      solvEnergy->solver = PBICGSTAB;
      solvEnergy->tol = smachn();
      solvEnergy->maxIt = 50000;
      solvEnergy->fileSolv = NULL;
      solvEnergy->log = true;
      solvEnergy->flag = true;
/*...................................................................*/

/*...*/
      if (solvEnergy->log && !mpiVar.myId) 
      {
        strcpy(auxName, preName);
        strcat(auxName, "_fluid_energy");
        fName(auxName, mpiVar.nPrcs, 0, 11, nameOut);
        solvEnergy->fileSolv = openFile(nameOut, "w");
      }
/*...................................................................*/

/*... inicializa a estrutura do solver(Energia)*/
      sistEqEnergy->unsym = true;
/*...................................................................*/

/*... solver*/
      readMacro(fileIn, word, false);
      setSolverConfig(word, solvEnergy, &sistEqEnergy->storage);
/*...................................................................*/

/*... numeracao das equacoes das Energy*/
      HccaAlloc(INT, m, sistEqEnergy->id
        , mesh->numel
        , "sistEnergyId", _AD_);

      fprintf(fileLogExc, "%s\n", DIF);
      fprintf(fileLogExc, "Numerando as equacoes.\n");
      tm.numeqEnergy = getTimeC() - tm.numeqEnergy;
      sistEqEnergy->neq = numEqV1(sistEqEnergy->id, reordMesh->num
        , mesh->numel);
      tm.numeqEnergy = getTimeC() - tm.numeqEnergy;

      fprintf(fileLogExc, "Equacoes numeradas.\n");
      fprintf(fileLogExc, "%s\n", DIF);
/*...................................................................*/

/*...*/
      if (mpiVar.nPrcs > 1) 
      {
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
      HccaAlloc(DOUBLE, m, sistEqEnergy->b0
        , sistEqEnergy->neq, "sistEnergy0", _AD_);
      HccaAlloc(DOUBLE, m, sistEqEnergy->b
        , sistEqEnergy->neq, "sistEnergyb ", _AD_);
      HccaAlloc(DOUBLE, m, sistEqEnergy->x
        , sistEqEnergy->neq, "sistEnergyx ", _AD_);
      zero(sistEqEnergy->b0, sistEqEnergy->neq, DOUBLEC);
      zero(sistEqEnergy->b, sistEqEnergy->neq, DOUBLEC);
      zero(sistEqEnergy->x, sistEqEnergy->neq, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados energia*/
      strcpy(strIa, "iaEnergy");
      strcpy(strJa, "jaEnergy");
      strcpy(strAd, "adEnergy");
      strcpy(strA, "aEnergy");

      fprintf(fileLogExc, "Energy:\n");
      fprintf(fileLogExc, "Montagem da estrura de dados esparsa.\n");

      tm.dataStructEnergy = getTimeC() - tm.dataStructEnergy;
      //        dataStructSimple(&m,sistEqEnergy->id ,reordMesh->num
      //                        ,mesh->elm.adj.nelcon
      //                        ,mesh->elm.adj.nViz  ,mesh->numelNov,mesh->maxViz  
      //                       ,1                   ,strIa         ,strJa
      //                       ,strAd               ,strA          ,sistEqEnergy);
      dataStruct(m, sistEqEnergy->id, reordMesh->num
        , mesh->elm.adj.nelcon
        , mesh->elm.adj.nViz, mesh->numelNov, mesh->maxViz
        , 1, strIa, strJa
        , strAd, strA, sistEqEnergy);
      tm.dataStructEnergy = getTimeC() - tm.dataStructEnergy;

      fprintf(fileLogExc, "Estrutuda montada.\n");
/*...................................................................*/

/*... dividindo a matriz*/
/*... Openmp(energy)*/
      if (ompVar.fSolver) 
      {
        strcpy(str1, "thBeginEnergy");
        strcpy(str2, "thEndEnergy");
        strcpy(str3, "thSizeEnergy");
        strcpy(str4, "thHeightEnergy");
        pMatrixSolverOmp(m, sistEqEnergy, str1, str2, str3, str4);
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... kTurb*/
    else if (!strcmp(word, "kTurb") || !strcmp(word, "kturb")) 
    {
      fprintf(fileLogExc, "%s.\n",word);
      nSistEq--;
/*... inicializando a estrutura de equacoes do problema
      (energica cinetica turbulenta)*/
      *fSolvKturb = true;
      solvKturb->solver = PBICGSTAB;
      solvKturb->tol = smachn();
      solvKturb->maxIt = 50000;
      solvKturb->fileSolv = NULL;
      solvKturb->log = true;
      solvKturb->flag = true;
/*...................................................................*/

/*...*/
      if (solvKturb->log && !mpiVar.myId) 
      {
        strcpy(auxName, preName);
        strcat(auxName, "_fluid_kturb");
        fName(auxName, mpiVar.nPrcs, 0, 11, nameOut);
        solvKturb->fileSolv = openFile(nameOut, "w");
      }
/*...................................................................*/

/*... inicializa a estrutura do solver(Energia cinetica turbulenta)*/
      sistEqKturb->unsym = true;
/*...................................................................*/

/*... solver*/
      readMacro(fileIn, word, false);
      setSolverConfig(word, solvKturb, &sistEqKturb->storage);
/*...................................................................*/

/*... numeracao das equacoes (Energia cinetica turbulenta)*/
      HccaAlloc(INT, m, sistEqKturb->id
        , mesh->numel
        , "sistKturbId", _AD_);

      fprintf(fileLogExc, "%s\n", DIF);
      fprintf(fileLogExc, "Numerando as equacoes.\n");

      tm.turbulence = getTimeC() - tm.turbulence;
      sistEqKturb->neq = numEqV1(sistEqKturb->id, reordMesh->num
        , mesh->numel);
      tm.turbulence = getTimeC() - tm.turbulence;

      fprintf(fileLogExc, "Equacoes numeradas.\n");
      fprintf(fileLogExc, "%s\n", DIF);
/*...................................................................*/

/*...*/
      if (mpiVar.nPrcs > 1)
      {
        //          tm.numeqPres = getTimeC() - tm.numeqPres;
        //      sistEqT1->neqNov = countEq(reordMesh->num
        //                          ,mesh->elm.faceRt1  ,mesh->elm.adj.nViz
        //                          ,mesh->numelNov     ,mesh->maxViz
        //                          ,mesh->ndfT[0]);
        //          tm.numeqPres = getTimeC() - tm.numeqPres;
      }
      else
        sistEqKturb->neqNov = sistEqKturb->neq;
/*...................................................................*/

/*... energia*/
      HccaAlloc(DOUBLE, m, sistEqKturb->b0
        , sistEqKturb->neq, "sistKturb0", _AD_);
      HccaAlloc(DOUBLE, m, sistEqKturb->b
        , sistEqKturb->neq, "sistKturb", _AD_);
      HccaAlloc(DOUBLE, m, sistEqKturb->x
        , sistEqKturb->neq, "sistKturbx ", _AD_);
      zero(sistEqKturb->b0, sistEqKturb->neq, DOUBLEC);
      zero(sistEqKturb->b, sistEqKturb->neq, DOUBLEC);
      zero(sistEqKturb->x, sistEqKturb->neq, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados energia*/
      strcpy(strIa, "iaKturb");
      strcpy(strJa, "jaKturb");
      strcpy(strAd, "adKturb");
      strcpy(strA, "akTurb");

      fprintf(fileLogExc, "Kturb:\n");
      fprintf(fileLogExc, "Montagem da estrura de dados esparsa.\n");

      tm.turbulence = getTimeC() - tm.turbulence;
//        dataStructSimple(&m,sistEqEnergy->id ,reordMesh->num
//                        ,mesh->elm.adj.nelcon
//                        ,mesh->elm.adj.nViz  ,mesh->numelNov,mesh->maxViz  
//                       ,1                   ,strIa         ,strJa
//                       ,strAd               ,strA          ,sistEqEnergy);
      dataStruct(m, sistEqKturb->id, reordMesh->num
        , mesh->elm.adj.nelcon
        , mesh->elm.adj.nViz, mesh->numelNov, mesh->maxViz
        , 1, strIa, strJa
        , strAd, strA, sistEqKturb);
      tm.turbulence = getTimeC() - tm.turbulence;

      fprintf(fileLogExc, "Estrutuda montada.\n");
/*...................................................................*/

/*... dividindo a matriz*/
/*... Openmp(energy)*/
      if (ompVar.fSolver)
      {
        strcpy(str1, "thBeginKturb");
        strcpy(str2, "thEndKturb");
        strcpy(str3, "thSizeKturb");
        strcpy(str4, "thHeightKturb");
        pMatrixSolverOmp(m, sistEqKturb, str1, str2, str3, str4);
      }
/*...................................................................*/
    }
/*...................................................................*/
  } while (nSistEq);
/*...................................................................*/

/*...*/
  if (opt->fItPlot && !mpiVar.myId) {
    strcpy(auxName, preName);
    fName(auxName, mpiVar.nPrcs, 0, 22, nameOut);
    opt->fileItPlot[FITPLOTSIMPLE] = openFile(nameOut, "w");
    if (mesh->ndfF == 3)
      fprintf(opt->fileItPlot[FITPLOTSIMPLE]
        , "#VelPres\n#it ||rU1||| ||rU2|| ||rMass||\n");
    else if (mesh->ndfF == 4)
      fprintf(opt->fileItPlot[FITPLOTSIMPLE]
        , "#VelPres\n#it ||rU1|| ||rU2|| ||rU3|| ||rMass||\n");
    else if (mesh->ndfFt == 4)
      fprintf(opt->fileItPlot[FITPLOTSIMPLE]
        , "#VelPres\n#it ||rU1|| ||rU2|| ||rMass|| ||rEnergy||\n");
    else if (mesh->ndfFt == 5)
      fprintf(opt->fileItPlot[FITPLOTSIMPLE]
        , "#VelPres\n#it ||rU1|| ||rU2|| ||rU3|| ||rMass|| ||rEnergy||\n");
  }
/*...................................................................*/

/*... Openmp(Vel,Pres)*/
  if (ompVar.fSolver) 
  {
/*... alocando o buffer*/
    if (*fSolvPres && *fSolvVel && *fSolvEnergy && *fSolvKturb) 
    {   
      nEqMax = max(sistEqPres->neqNov, sistEqVel->neqNov);
      nEqMax = max(sistEqEnergy->neqNov, nEqMax);
      nEqMax = max(sistEqKturb->neqNov, nEqMax);
      nEqMax = max(nEqMax, ompVar.nEqMax);
      ompVar.nEqMax = nEqMax;
      HccaAlloc(DOUBLE, m, ompVar.buffer
        , nEqMax*ompVar.nThreadsSolver, "bufferOmp", false);
      zero(ompVar.buffer, nEqMax*ompVar.nThreadsSolver, DOUBLEC);
      sistEqPres->omp.thY   = ompVar.buffer;
      sistEqVel->omp.thY    = ompVar.buffer;
      sistEqEnergy->omp.thY = ompVar.buffer;
      sistEqKturb->omp.thY  = ompVar.buffer;
    }
/*... alocando o buffer*/
    else if (*fSolvPres && *fSolvVel && *fSolvEnergy)
    {
      nEqMax = max(sistEqPres->neqNov, sistEqVel->neqNov);
      nEqMax = max(sistEqEnergy->neqNov, nEqMax);
      nEqMax = max(nEqMax, ompVar.nEqMax);
      ompVar.nEqMax = nEqMax;
      HccaAlloc(DOUBLE, m, ompVar.buffer
        , nEqMax*ompVar.nThreadsSolver, "bufferOmp", false);
      zero(ompVar.buffer, nEqMax*ompVar.nThreadsSolver, DOUBLEC);
      sistEqPres->omp.thY   = ompVar.buffer;
      sistEqVel->omp.thY    = ompVar.buffer;
      sistEqEnergy->omp.thY = ompVar.buffer;
    }

/*... alocando o buffer*/
    else if (*fSolvPres && *fSolvVel) 
    {
      nEqMax = max(sistEqPres->neqNov, sistEqVel->neqNov);
      nEqMax = max(nEqMax, ompVar.nEqMax);
      ompVar.nEqMax = nEqMax;
      HccaAlloc(DOUBLE, m, ompVar.buffer
        , nEqMax*ompVar.nThreadsSolver, "bufferOmp", false);
      zero(ompVar.buffer, nEqMax*ompVar.nThreadsSolver, DOUBLEC);
      ompVar.nEqMax = nEqMax;
      sistEqPres->omp.thY = ompVar.buffer;
      sistEqVel->omp.thY  = ompVar.buffer;
    }
/*... alocando o buffer*/
    else if (*fSolvPres)
    {
      nEqMax = sistEqPres->neqNov;
      nEqMax = max(nEqMax, ompVar.nEqMax);
      ompVar.nEqMax = nEqMax;
      HccaAlloc(DOUBLE, m, ompVar.buffer
        , nEqMax*ompVar.nThreadsSolver, "bufferOmp", false);
      zero(ompVar.buffer, nEqMax*ompVar.nThreadsSolver, DOUBLEC);
      ompVar.nEqMax = nEqMax;
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
  if (!mpiVar.myId) 
  {
    strcpy(str1, "MB");
    memoriaTotal(str1);
    usoMemoria(m, str1);
  }
/*...................................................................*/

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 30/07/2018                                   *
 * Data de modificaco : 20/08/2019                                   *
 *-------------------------------------------------------------------*
 * readSolvComb   leitura das configuracoes do solvers utilizados    *
 * no modelo de combustao                                            *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void readSolvComb(Memoria *m      , Mesh *mesh          , Reord *reordMesh
                , Combustion *cModel 
                , Solv *solvVel   , SistEq* sistEqVel   , bool *fSolvVel
                , Solv *solvPres  , SistEq* sistEqPres  , bool *fSolvPres
                , Solv *solvEnergy, SistEq* sistEqEnergy, bool *fSolvEnergy
                , Solv *solvKturb , SistEq* sistEqKturb , bool *fSolvKturb
                , Solv *solvComb  , SistEq* sistEqComb  , bool *fSolvComb
                , PartMesh *pMesh 
                , char* auxName   , char* preName       , char* nameOut
                , FILE *fileIn    , FileOpt *opt)
{

  unsigned short i, ndfVel, nSistEq, ndfComb = cModel->nComb;

  INT nEqMax;
/*... Estrutura de dados*/
  char strIa[MNOMEPONTEIRO], strJa[MNOMEPONTEIRO];
  char strA[MNOMEPONTEIRO], strAd[MNOMEPONTEIRO];

  char str1[100], str2[100], str3[100], str4[100];

  char word[WORD_SIZE];

/*... solver*/
  readMacro(fileIn, word, false);
  nSistEq = (short)atol(word);
  do
  {
    readMacro(fileIn, word, false);
/*... velocidade*/
    if (!strcmp(word, "Vel") || !strcmp(word, "vel"))
    {
      nSistEq--;
      *fSolvVel = true;
      solvVel->solver = PBICGSTAB;
      solvVel->tol = smachn();
      solvVel->maxIt = 50000;
      solvVel->fileSolv = NULL;
      solvVel->log = true;
      solvVel->flag = true;
/*...................................................................*/

/*...*/
      if (solvVel->log && !mpiVar.myId)
      {
        strcpy(auxName, preName);
        strcat(auxName, "_fluid_vel");
        fName(auxName, mpiVar.nPrcs, 0, 11, nameOut);
        solvVel->fileSolv = openFile(nameOut, "w");
      }
/*...................................................................*/

      sistEqVel->unsym = true;
/*...................................................................*/

/*... solver*/
      readMacro(fileIn, word, false);
      setSolverConfig(word, solvVel, &sistEqVel->storage);
/*...................................................................*/

/*... DataStruct*/
      readMacro(fileIn, word, false);
      setDataStruct(word, &sistEqVel->storage);
/*...................................................................*/

/*... numeracao das equacoes das velocidades*/
      HccaAlloc(INT, m, sistEqVel->id
        , mesh->numel
        , "sistVelid", _AD_);

      fprintf(fileLogExc, "%s\n", DIF);
      fprintf(fileLogExc, "Numerando as equacoes.\n");

      tm.numeqVel = getTimeC() - tm.numeqVel;
      sistEqVel->neq = numEqV1(sistEqVel->id, reordMesh->num
        , mesh->numel);
      tm.numeqVel = getTimeC() - tm.numeqVel;

      fprintf(fileLogExc, "Equacoes numeradas.\n");
      fprintf(fileLogExc, "%s\n", DIF);
/*...................................................................*/

/*...*/
      if (mpiVar.nPrcs > 1)
      {
        tm.numeqVel = getTimeC() - tm.numeqVel;
        sistEqVel->neqNov = countEq(reordMesh->num
                                  ,mesh->elm.faceRvel ,mesh->elm.adj.nViz
                                  ,mesh->numelNov     ,mesh->maxViz
                                  ,1);
        tm.numeqVel = getTimeC() - tm.numeqVel;
      }
      else
        sistEqVel->neqNov = sistEqVel->neq;
/*...................................................................*/

/*... velovidades*/
      ndfVel = max(mesh->ndfF - 1, mesh->ndfFt - 2);
/*...................................................................*/
      HccaAlloc(DOUBLE, m, sistEqVel->b0
        , sistEqVel->neqNov*ndfVel, "sistVelb0", _AD_);
      HccaAlloc(DOUBLE, m, sistEqVel->b
        , sistEqVel->neqNov*ndfVel, "sistVelb ", _AD_);
      HccaAlloc(DOUBLE, m, sistEqVel->x
        , sistEqVel->neq*ndfVel, "sistVelx ", _AD_);
      zero(sistEqVel->b0, sistEqVel->neqNov*ndfVel, DOUBLEC);
      zero(sistEqVel->b, sistEqVel->neqNov*ndfVel, DOUBLEC);
      zero(sistEqVel->x, sistEqVel->neq*ndfVel, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados velocidades*/
      strcpy(strIa, "iaVel");
      strcpy(strJa, "jaVel");
      strcpy(strAd, "adVel");
      strcpy(strA, "aVel");

      fprintf(fileLogExc, "Vel:\n");
      fprintf(fileLogExc, "Montagem da estrura de dados esparsa.\n");

      tm.dataStructVel = getTimeC() - tm.dataStructVel;
      dataStructSimple(m, sistEqVel->id, reordMesh->num
        , mesh->elm.adj.nelcon
        , mesh->elm.adj.nViz, mesh->numelNov, mesh->maxViz
        , ndfVel, strIa, strJa
        , strAd, strA, sistEqVel);
      tm.dataStructVel = getTimeC() - tm.dataStructVel;

      fprintf(fileLogExc, "Estrutuda montada.\n");
/*...................................................................*/

/*... dividindo a matriz*/
/*... Openmp(Vel,Pres)*/
      if (ompVar.fSolver)
      {
        strcpy(str1, "thBeginVel");
        strcpy(str2, "thEndVel");
        strcpy(str3, "thSizeVel");
        strcpy(str4, "thHeightVel");
        pMatrixSolverOmp(m, sistEqVel, str1, str2, str3, str4);
      }
/*...................................................................*/

/*... mapa de equacoes para comunicacao*/
      if( mpiVar.nPrcs > 1)
        front(m           ,pMesh
             ,sistEqVel   
             ,"fMapNeqVel"      ,"iaSendsNeqVel"
             ,"iaRcvsNeqVel"    ,"xBufferMpiNeqVel"
             ,"nVizPartNeqVel"  ,1);      
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    else if (!strcmp(word, "Pres") || !strcmp(word, "pres"))
    {
      nSistEq--;
      *fSolvPres = true;
      solvPres->solver = PCG;
      solvPres->tol = smachn();
      solvPres->maxIt = 50000;
      solvPres->fileSolv = NULL;
      solvPres->log = true;
      solvPres->flag = true;
/*...................................................................*/

/*...*/
      if (solvPres->log && !mpiVar.myId)
      {
        strcpy(auxName, preName);
        strcat(auxName, "_fluid_pres");
        fName(auxName, mpiVar.nPrcs, 0, 11, nameOut);
        solvPres->fileSolv = openFile(nameOut, "w");
      }
/*...................................................................*/

/*... inicializa a estrutura do solver(PRESSAO)*/
      sistEqPres->unsym = false;
/*...................................................................*/

/*... solver*/
      readMacro(fileIn, word, false);
      setSolverConfig(word, solvPres, &sistEqPres->storage);
/*...................................................................*/

/*... DataStruct*/
      readMacro(fileIn, word, false);
      setDataStruct(word, &sistEqPres->storage);
/*...................................................................*/

/*... numeracao das equacoes das pressoes*/
      HccaAlloc(INT, m, sistEqPres->id
        , mesh->numel
        , "sistPresid", _AD_);

      fprintf(fileLogExc, "%s\n", DIF);
      fprintf(fileLogExc, "Numerando as equacoes.\n");

      tm.numeqPres = getTimeC() - tm.numeqPres;
      sistEqPres->neq = numEqV2(sistEqPres->id, reordMesh->num
        , mesh->elm.faceRpres, mesh->elm.adj.nViz
        , mesh->numel, mesh->maxViz);
      tm.numeqPres = getTimeC() - tm.numeqPres;

      fprintf(fileLogExc, "Equacoes numeradas.\n");
      fprintf(fileLogExc, "%s\n", DIF);
/*...................................................................*/

/*...*/
      if (mpiVar.nPrcs > 1)
      {
        tm.numeqPres = getTimeC() - tm.numeqPres;
        sistEqPres->neqNov = countEq(reordMesh->num
                                  ,mesh->elm.faceRpres,mesh->elm.adj.nViz
                                  ,mesh->numelNov     ,mesh->maxViz
                                  ,1);
        tm.numeqPres = getTimeC() - tm.numeqPres;
      }
      else
        sistEqPres->neqNov = sistEqPres->neq;
/*...................................................................*/

/*... pressoes*/
      HccaAlloc(DOUBLE, m, sistEqPres->b0
        , sistEqPres->neqNov, "sistPresb0", _AD_);
      HccaAlloc(DOUBLE, m, sistEqPres->b
        , sistEqPres->neqNov, "sistPresb", _AD_);
      HccaAlloc(DOUBLE, m, sistEqPres->x
        , sistEqPres->neq, "sistPresx", _AD_);
      zero(sistEqPres->b0, sistEqPres->neqNov, DOUBLEC);
      zero(sistEqPres->b, sistEqPres->neqNov, DOUBLEC);
      zero(sistEqPres->x, sistEqPres->neq, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados pressoes*/
      strcpy(strIa, "iaPres");
      strcpy(strJa, "japres");
      strcpy(strAd, "adPres");
      strcpy(strA, "aPres");

      fprintf(fileLogExc, "Pres:\n");
      fprintf(fileLogExc, "Montagem da estrura de dados esparsa.\n");

      tm.dataStructPres = getTimeC() - tm.dataStructPres;
      dataStruct(m, sistEqPres->id, reordMesh->num, mesh->elm.adj.nelcon
        , mesh->elm.adj.nViz, mesh->numelNov, mesh->maxViz
        , 1, strIa, strJa
        , strAd, strA, sistEqPres);
      tm.dataStructPres = getTimeC() - tm.dataStructPres;
      fprintf(fileLogExc, "Estrutuda montada.\n");
/*...................................................................*/

/*... Openmp(Vel,Pres)*/
      if (ompVar.fSolver)
      {
/*... dividindo a matriz*/
        strcpy(str1, "thBeginPres");
        strcpy(str2, "thEndPres");
        strcpy(str3, "thSizePres");
        strcpy(str4, "thHeightPres");
        pMatrixSolverOmp(m, sistEqPres, str1, str2, str3, str4);
      }
/*...................................................................*/

/*... mapa de equacoes para comunicacao*/
      if( mpiVar.nPrcs > 1)
        front(m           ,pMesh
             ,sistEqPres  
             ,"fMapNeqP"    ,"iaSendsNeqP"
             ,"iaRcvsNeqP"  ,"xBufferMpiNeqP"
             ,"nVizPartNeqP"  ,1);    
/*...................................................................*/
    }
/*...................................................................*/

/*... energy*/
    else if (!strcmp(word, "Energy") || !strcmp(word, "energy"))
    {
      nSistEq--;
/*... inicializando a estrutura de equacoes do problema (ENERGY)*/
      *fSolvEnergy = true;
      solvEnergy->solver = PBICGSTAB;
      solvEnergy->tol = smachn();
      solvEnergy->maxIt = 50000;
      solvEnergy->fileSolv = NULL;
      solvEnergy->log = true;
      solvEnergy->flag = true;
/*...................................................................*/

/*...*/
      if (solvEnergy->log && !mpiVar.myId)
      {
        strcpy(auxName, preName);
        strcat(auxName, "_fluid_energy");
        fName(auxName, mpiVar.nPrcs, 0, 11, nameOut);
        solvEnergy->fileSolv = openFile(nameOut, "w");
      }
/*...................................................................*/

/*... inicializa a estrutura do solver(Energia)*/
      sistEqEnergy->unsym = true;
/*...................................................................*/

/*... solver*/
      readMacro(fileIn, word, false);
      setSolverConfig(word, solvEnergy, &sistEqEnergy->storage);
/*...................................................................*/

/*... DataStruct*/
      readMacro(fileIn, word, false);
      setDataStruct(word, &sistEqEnergy->storage);
/*...................................................................*/

/*... numeracao das equacoes das Energy*/
      HccaAlloc(INT, m, sistEqEnergy->id
        , mesh->numel
        , "sistEnergyId", _AD_);

      fprintf(fileLogExc, "%s\n", DIF);
      fprintf(fileLogExc, "Numerando as equacoes.\n");

      tm.numeqEnergy = getTimeC() - tm.numeqEnergy;
      sistEqEnergy->neq = numEqV1(sistEqEnergy->id, reordMesh->num
        , mesh->numel);
      tm.numeqEnergy = getTimeC() - tm.numeqEnergy;

      fprintf(fileLogExc, "Equacoes numeradas.\n");
      fprintf(fileLogExc, "%s\n", DIF);
/*...................................................................*/

/*...*/
      if (mpiVar.nPrcs > 1)
      {
        tm.numeqEnergy = getTimeC() - tm.numeqEnergy;
        sistEqEnergy->neqNov = countEq(reordMesh->num
                                  ,mesh->elm.faceRenergy,mesh->elm.adj.nViz
                                  ,mesh->numelNov       ,mesh->maxViz
                                  ,1);
        tm.numeqEnergy = getTimeC() - tm.numeqEnergy;
      }
      else
        sistEqEnergy->neqNov = sistEqEnergy->neq;
/*...................................................................*/

/*... energia*/
      HccaAlloc(DOUBLE, m, sistEqEnergy->b0
        , sistEqEnergy->neqNov, "sistEnergy0", _AD_);
      HccaAlloc(DOUBLE, m, sistEqEnergy->b
        , sistEqEnergy->neqNov, "sistEnergyb", _AD_);
      HccaAlloc(DOUBLE, m, sistEqEnergy->x
        , sistEqEnergy->neq, "sistEnergyx ", _AD_);
      zero(sistEqEnergy->b0, sistEqEnergy->neq, DOUBLEC);
      zero(sistEqEnergy->b, sistEqEnergy->neq, DOUBLEC);
      zero(sistEqEnergy->x, sistEqEnergy->neq, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados energia*/
      strcpy(strIa, "iaEnergy");
      strcpy(strJa, "jaEnergy");
      strcpy(strAd, "adEnergy");
      strcpy(strA, "aEnergy");

      fprintf(fileLogExc, "Energy:\n");
      fprintf(fileLogExc, "Montagem da estrura de dados esparsa.\n");

      tm.dataStructEnergy = getTimeC() - tm.dataStructEnergy;
      dataStruct(m, sistEqEnergy->id, reordMesh->num
                , mesh->elm.adj.nelcon
                , mesh->elm.adj.nViz, mesh->numelNov, mesh->maxViz
                , 1, strIa, strJa
                , strAd, strA, sistEqEnergy);
      tm.dataStructEnergy = getTimeC() - tm.dataStructEnergy;

      fprintf(fileLogExc, "Estrutuda montada.\n");
/*...................................................................*/

/*... dividindo a matriz*/
/*... Openmp(energy)*/
      if (ompVar.fSolver)
      {
        strcpy(str1, "thBeginEnergy");
        strcpy(str2, "thEndEnergy");
        strcpy(str3, "thSizeEnergy");
        strcpy(str4, "thHeightEnergy");
        pMatrixSolverOmp(m, sistEqEnergy, str1, str2, str3, str4);
      }
/*...................................................................*/

/*... mapa de equacoes para comunicacao*/
      if( mpiVar.nPrcs > 1)
        front(m           ,pMesh
             ,sistEqEnergy
             ,"fMapNeqE"    ,"iaSendsNeqE"
             ,"iaRcvsNeqE"  ,"xBufferMpiNeqE"
             ,"nVizPartNeqE"  ,1);   
/*...................................................................*/

    }
/*...................................................................*/

/*... kTurb*/
    else if (!strcmp(word, "kTurb") || !strcmp(word, "kturb"))
    {
      nSistEq--;
/*... inicializando a estrutura de equacoes do problema
      (energica cinetica turbulenta)*/
      *fSolvKturb = true;
      solvKturb->solver = PBICGSTAB;
      solvKturb->tol = smachn();
      solvKturb->maxIt = 50000;
      solvKturb->fileSolv = NULL;
      solvKturb->log = true;
      solvKturb->flag = true;
/*...................................................................*/

/*...*/
      if (solvKturb->log && !mpiVar.myId)
      {
        strcpy(auxName, preName);
        strcat(auxName, "_fluid_kturb");
        fName(auxName, mpiVar.nPrcs, 0, 11, nameOut);
        solvKturb->fileSolv = openFile(nameOut, "w");
      }
/*...................................................................*/

/*... inicializa a estrutura do solver(Energia cinetica turbulenta)*/
      sistEqKturb->unsym = true;
/*...................................................................*/

/*... solver*/
      readMacro(fileIn, word, false);
      setSolverConfig(word, solvKturb, &sistEqKturb->storage);
/*...................................................................*/

/*... numeracao das equacoes (Energia cinetica turbulenta)*/
      HccaAlloc(INT, m, sistEqKturb->id
        , mesh->numel
        , "sistKturbId", _AD_);

      fprintf(fileLogExc, "%s\n", DIF);
      fprintf(fileLogExc, "Numerando as equacoes.\n");

      tm.turbulence = getTimeC() - tm.turbulence;
      sistEqKturb->neq = numEqV1(sistEqKturb->id, reordMesh->num
        , mesh->numel);
      tm.turbulence = getTimeC() - tm.turbulence;

/*...................................................................*/

/*...*/
      if (mpiVar.nPrcs > 1)
      {
        tm.turbulence = getTimeC() - tm.turbulence;
        sistEqKturb->neqNov = countEq(reordMesh->num
                                  ,mesh->elm.faceReKturb  ,mesh->elm.adj.nViz
                                  ,mesh->numelNov         ,mesh->maxViz
                                  ,1);
        tm.turbulence = getTimeC() - tm.turbulence;
      }
      else
        sistEqKturb->neqNov = sistEqKturb->neq;
/*...................................................................*/

/*... energia*/
      HccaAlloc(DOUBLE, m, sistEqKturb->b0
        , sistEqKturb->neqNov, "sistKturb0", _AD_);
      HccaAlloc(DOUBLE, m, sistEqKturb->b
        , sistEqKturb->neqNov, "sistKturb", _AD_);
      HccaAlloc(DOUBLE, m, sistEqKturb->x
        , sistEqKturb->neq, "sistKturbx ", _AD_);
      zero(sistEqKturb->b0, sistEqKturb->neqNov, DOUBLEC);
      zero(sistEqKturb->b, sistEqKturb->neqNov, DOUBLEC);
      zero(sistEqKturb->x, sistEqKturb->neq, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados energia*/
      strcpy(strIa, "iaKturb");
      strcpy(strJa, "jaKturb");
      strcpy(strAd, "adKturb");
      strcpy(strA, "akTurb");

      tm.turbulence = getTimeC() - tm.turbulence;
         dataStruct(m, sistEqKturb->id, reordMesh->num
        , mesh->elm.adj.nelcon
        , mesh->elm.adj.nViz, mesh->numelNov, mesh->maxViz
        , 1, strIa, strJa
        , strAd, strA, sistEqKturb);
      tm.turbulence = getTimeC() - tm.turbulence;

      fprintf(fileLogExc, "Estrutuda montada.\n");
/*...................................................................*/

/*... dividindo a matriz*/
/*... Openmp(energy)*/
      if (ompVar.fSolver)
      {
        strcpy(str1, "thBeginKturb");
        strcpy(str2, "thEndKturb");
        strcpy(str3, "thSizeKturb");
        strcpy(str4, "thHeightKturb");
        pMatrixSolverOmp(m, sistEqKturb, str1, str2, str3, str4);
      }
/*...................................................................*/

/*... mapa de equacoes para comunicacao*/
      if( mpiVar.nPrcs > 1)
        front(m                 ,pMesh
             ,sistEqKturb   
             ,"fMapNeqKt"      ,"iaSendsNeqKt"
             ,"iaRcvsNeqKt"    ,"xBufferMpiNeqKt"
             ,"nVizPartNeqKt"  ,1);      
/*...................................................................*/

    }
/*...................................................................*/

/*... fracao massiva*/
    else if (!strcmp(word, "z") || !strcmp(word, "Z")) 
    {
      nSistEq--;
      *fSolvComb = true;
      solvComb->solver = PBICGSTAB;
      solvComb->tol = smachn();
      solvComb->maxIt = 50000;
      solvComb->fileSolv = NULL;
      solvComb->log = true;
      solvComb->flag = true;
/*...................................................................*/

/*...*/
      if (solvComb->log && !mpiVar.myId) 
      {
        strcpy(auxName, preName);
        strcat(auxName, "_fluid_comb");
        fName(auxName, mpiVar.nPrcs, 0, 11, nameOut);
        solvComb->fileSolv = openFile(nameOut, "w");
      }
/*...................................................................*/

      sistEqComb->unsym = true;
/*...................................................................*/

/*... solver*/
      readMacro(fileIn, word, false);
      setSolverConfig(word, solvComb, &sistEqComb->storage);
/*...................................................................*/

/*... numeracao das equacoes das velocidades*/
      HccaAlloc(INT, m, sistEqComb->id
              , mesh->numel
              , "sistCombId", _AD_);

      tm.numeqComb = getTimeC() - tm.numeqComb;

      sistEqComb->neq = numEqV1(sistEqComb->id, reordMesh->num
                              , mesh->numel);    

      tm.numeqComb = getTimeC() - tm.numeqComb;
/*...................................................................*/

/*...*/
      if (mpiVar.nPrcs > 1) 
      {
        tm.numeqComb = getTimeC() - tm.numeqComb;
        sistEqComb->neqNov = countEq(reordMesh->num
                            ,mesh->elm.faceResZcomb  ,mesh->elm.adj.nViz
                            ,mesh->numelNov          ,mesh->maxViz
                            ,1            );
        tm.numeqComb = getTimeC() - tm.numeqComb;
      }
      else
        sistEqComb->neqNov = sistEqComb->neq;
/*...................................................................*/

/*... */
      HccaAlloc(DOUBLE, m, sistEqComb->b0
              , sistEqComb->neqNov*ndfComb, "sistCombb0", _AD_);
      HccaAlloc(DOUBLE, m, sistEqComb->b
              , sistEqComb->neqNov*ndfComb, "sistCombb ", _AD_);
      HccaAlloc(DOUBLE, m, sistEqComb->x
              , sistEqComb->neq*ndfComb, "sistCombx ", _AD_);
      zero(sistEqComb->b0, sistEqComb->neqNov*ndfComb, DOUBLEC);
      zero(sistEqComb->b , sistEqComb->neqNov*ndfComb, DOUBLEC);
      zero(sistEqComb->x , sistEqComb->neq*ndfComb, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados velocidades*/
      strcpy(strIa, "iaComb");
      strcpy(strJa, "jaComb");
      strcpy(strAd, "adComb");
      strcpy(strA , "aComb");

      fprintf(fileLogExc, "Combustion:\n");
      fprintf(fileLogExc, "Montagem da estrura de dados esparsa.\n");

      tm.dataStructComb = getTimeC() - tm.dataStructComb;
      dataStructBlock(m                , sistEqComb->id   
                    , reordMesh->num   , mesh->elm.adj.nelcon
                    , mesh->elm.adj.nViz
                    , mesh->numelNov   , mesh->maxViz
                    , ndfComb          , ndfComb
                    , strIa            , strJa
                    , strAd            , strA
                    , sistEqComb);
      tm.dataStructComb = getTimeC() - tm.dataStructComb;

      fprintf(fileLogExc, "Estrutuda montada.\n");
/*...................................................................*/

/*... dividindo a matriz*/
/*... Openmp(Vel,Pres)*/
      if (ompVar.fSolver) 
      {
        strcpy(str1, "thBeginComb");
        strcpy(str2, "thEndComb");
        strcpy(str3, "thSizeComb");
        strcpy(str4, "thHeightComb");
        pMatrixSolverOmp(m, sistEqComb, str1, str2, str3, str4);
      }
/*...................................................................*/

/*... mapa de equacoes para comunicacao*/
      if( mpiVar.nPrcs > 1)
        front(m                  ,pMesh
             ,sistEqComb   
             ,"fMapNeqComb"      ,"iaSendsNeqComb"
             ,"iaRcvsNeqComb"    ,"xBufferMpiNeqComb"
             ,"nVizPartNeqComb"  ,1);        
/*...................................................................*/

    }
/*...................................................................*/

  } while (nSistEq);
/*...................................................................*/

/*...*/
  if (opt->fItPlot && !mpiVar.myId)
  {
    strcpy(auxName, preName);
    fName(auxName, mpiVar.nPrcs, 0, 22, nameOut);
    opt->fileItPlot[FITPLOTSIMPLE] = openFile(nameOut, "w");
/*...*/
    if (mesh->ndfFt == 4)
      fprintf(opt->fileItPlot[FITPLOTSIMPLE]
        , "#VelPres\n#it ||rU1|| ||rU2|| ||rMass|| ||rEnergy||\n");
    else if (mesh->ndfFt == 5)
      fprintf(opt->fileItPlot[FITPLOTSIMPLE]
        , "#VelPres\n#it ||rU1|| ||rU2|| ||rU3|| ||rMass|| ||rEnergy|| ");
/*...................................................................*/

/*...nComb*/
    if(ndfComb > 0)
    {
      fprintf(opt->fileItPlot[FITPLOTSIMPLE],"#Combustion\n#it");
      for(i=0;i<ndfComb;i++)
        fprintf(opt->fileItPlot[FITPLOTSIMPLE]," ||%s|| "
                                              ,cModel->chem.sp[i].name);  
      fprintf(opt->fileItPlot[FITPLOTSIMPLE],"\n");
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... Openmp(Vel,Pres)*/
  if (ompVar.fSolver)
  {
/*... alocando o buffer*/
    if (*fSolvKturb)
    {
      nEqMax = max(sistEqPres->neqNov  , sistEqVel->neqNov);
      nEqMax = max(sistEqEnergy->neqNov, nEqMax);
      nEqMax = max(sistEqKturb->neqNov , nEqMax);
      nEqMax = max(sistEqComb->neqNov  , nEqMax);
      nEqMax = max(nEqMax, ompVar.nEqMax);
      ompVar.nEqMax = nEqMax;
      HccaAlloc(DOUBLE, m, ompVar.buffer
        , nEqMax*ompVar.nThreadsSolver, "bufferOmp", false);
      zero(ompVar.buffer, nEqMax*ompVar.nThreadsSolver, DOUBLEC);
      sistEqPres->omp.thY   = ompVar.buffer;
      sistEqVel->omp.thY    = ompVar.buffer;
      sistEqEnergy->omp.thY = ompVar.buffer;
      sistEqKturb->omp.thY  = ompVar.buffer;
      sistEqComb->omp.thY   = ompVar.buffer;
    }
/*...................................................................*/

/*...*/
    else
    {
      nEqMax = max(sistEqPres->neqNov  , sistEqVel->neqNov);
      nEqMax = max(sistEqEnergy->neqNov, nEqMax);
      nEqMax = max(sistEqComb->neqNov  , nEqMax);
      nEqMax = max(nEqMax, ompVar.nEqMax);
      ompVar.nEqMax = nEqMax;
      HccaAlloc(DOUBLE, m, ompVar.buffer
              , nEqMax*ompVar.nThreadsSolver, "bufferOmp", false);
      zero(ompVar.buffer, nEqMax*ompVar.nThreadsSolver, DOUBLEC);
      sistEqPres->omp.thY   = ompVar.buffer;
      sistEqVel->omp.thY    = ompVar.buffer;
      sistEqEnergy->omp.thY = ompVar.buffer;
      sistEqComb->omp.thY   = ompVar.buffer;  
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... informacao da memoria total usada*/
  if (!mpiVar.myId)
  {
    strcpy(str1, "MB");
    memoriaTotal(str1);
    usoMemoria(m, str1);
  }
/*...................................................................*/

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 01/05/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * readNlIt  leitura das configuracoes do solvers utilizados         *
 * no escoamento de fluidos                                          *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void readNlIt(Scheme *sc, FILE *fileIn)
{
  unsigned short nCount;
  char word[WORD_SIZE];

/*... solver*/
  readMacro(fileIn, word, false);
  nCount = (short)atol(word);
  do 
  {
    readMacro(fileIn, word, false);
/*... D1*/
    if (!strcmp(word, "D1") || !strcmp(word, "d1")) 
    {
      nCount--;
      fscanf(fileIn, "%d" ,&sc->nlD1.maxIt);
      fscanf(fileIn, "%lf",&sc->nlD1.tol);
      fscanf(fileIn, "%d" ,&sc->nlD1.pPlot);
      fprintf(fileLogExc, "MaxIt: %d\n", sc->nlD1.maxIt);
      fprintf(fileLogExc, "Tol  : %e\n", sc->nlD1.tol);
      fprintf(fileLogExc, "pPlot: %d\n", sc->nlD1.pPlot);
    }
/*...................................................................*/

/*... T1*/
    if (!strcmp(word, "T1") || !strcmp(word, "t1"))
    {
      nCount--;
      fscanf(fileIn, "%d" , &sc->nlT1.maxIt);
      fscanf(fileIn, "%lf", &sc->nlT1.tol);
      fscanf(fileIn, "%d" , &sc->nlT1.pPlot);
      fprintf(fileLogExc, "MaxIt: %d\n", sc->nlT1.maxIt);
      fprintf(fileLogExc, "Tol  : %e\n", sc->nlT1.tol);
      fprintf(fileLogExc, "pPlot: %d\n", sc->nlT1.pPlot);
    }
/*...................................................................*/

  } while (nCount);
/*...................................................................*/

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 30/04/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * readSolvDifff: leitura das configuracoes do solvers utilizados    *
 * nos problemas de difusoes                                         *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void readSolvDiff(Memoria *m   , Mesh *mesh, Reord *reordMesh
                 , Solv *solvD1, SistEq* sistEqD1, bool *fSolvD1
                 , char* auxName, char* preName, char* nameOut
                 , FILE *fileIn, FileOpt *opt)
{

  unsigned short nSistEq;

  INT nEqMax;
/*... Estrutura de dados*/
  char strIa[MNOMEPONTEIRO], strJa[MNOMEPONTEIRO];
  char strA[MNOMEPONTEIRO], strAd[MNOMEPONTEIRO];

  char str1[100], str2[100], str3[100], str4[100];

  char word[WORD_SIZE];

/*... solver*/
  readMacro(fileIn, word, false);
  nSistEq = (short)atol(word);
  do 
  {
    readMacro(fileIn, word, false);
/*... equaocao de difusao d1*/
    if (!strcmp(word, "D1") || !strcmp(word, "d1")) 
    {
      nSistEq--;
      *fSolvD1         = true;
      solvD1->solver   = PCG;
      solvD1->tol      = smachn();
      solvD1->maxIt    = 50000;
      solvD1->fileSolv = NULL;
      solvD1->log      = true;
      solvD1->flag     = true;
/*...................................................................*/

/*...*/
      if (solvD1->log && !mpiVar.myId) 
      {
        strcpy(auxName, preName);
        strcat(auxName, "_D1");
        fName(auxName, mpiVar.nPrcs, 0, 11, nameOut);
        solvD1->fileSolv = openFile(nameOut, "w");
      }
/*...................................................................*/
      
/*...*/
      sistEqD1->storage = CSRD;
      sistEqD1->unsym = false;
/*...................................................................*/

/*... solver*/
      readMacro(fileIn, word, false);
      setSolverConfig(word, solvD1, &sistEqD1->storage);
/*...................................................................*/

/*... numeracao das equacoes das velocidades*/
      HccaAlloc(INT, m, sistEqD1->id, mesh->numel*mesh->ndfD[0]
               ,"sistD1id", _AD_);
      fprintf(fileLogExc, "%s\n", DIF);
      fprintf(fileLogExc, "Numerando as equacoes.\n");
      tm.numeqD1 = getTimeC() - tm.numeqD1;
      sistEqD1->neq = numeq(sistEqD1->id     , reordMesh->num
                          , mesh->elm.faceRd1, mesh->elm.adj.nViz
                          , mesh->numel      , mesh->maxViz
                          , mesh->ndfD[0]);
      tm.numeqD1 = getTimeC() - tm.numeqD1;

      fprintf(fileLogExc, "Equacoes numeradas.\n");
      fprintf(fileLogExc, "%s\n", DIF);
/*...................................................................*/

/*...*/
      if (mpiVar.nPrcs > 1) 
      {
        //          tm.numeqPres = getTimeC() - tm.numeqPres;
        //      sistEqT1->neqNov = countEq(reordMesh->num
        //                          ,mesh->elm.faceRt1  ,mesh->elm.adj.nViz
        //                          ,mesh->numelNov     ,mesh->maxViz
        //                          ,mesh->ndfT[0]);
        //          tm.numeqPres = getTimeC() - tm.numeqPres;
      }
      else
        sistEqD1->neqNov = sistEqD1->neq;
/*...................................................................*/

/*... */
      HccaAlloc(DOUBLE        , m         , sistEqD1->b0
               , sistEqD1->neq, "sistD1b0", _AD_);
      HccaAlloc(DOUBLE        , m         , sistEqD1->b
               , sistEqD1->neq, "sistD1b ", _AD_);
      HccaAlloc(DOUBLE       , m         , sistEqD1->x
              , sistEqD1->neq, "sistD1x ", _AD_);
      zero(sistEqD1->b0, sistEqD1->neq, DOUBLEC);
      zero(sistEqD1->b , sistEqD1->neq, DOUBLEC);
      zero(sistEqD1->x , sistEqD1->neq, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados velocidades*/
      strcpy(strIa, "iaD1");
      strcpy(strJa, "jaD1");
      strcpy(strAd, "adD1");
      strcpy(strA , "aD1");

      fprintf(fileLogExc, "D1:\n");
      fprintf(fileLogExc, "Montagem da estrura de dados esparsa.\n");

      tm.dataStructD1 = getTimeC() - tm.dataStructD1;
      dataStruct(m, sistEqD1->id, reordMesh->num, mesh->elm.adj.nelcon
               , mesh->elm.adj.nViz, mesh->numelNov, mesh->maxViz
               , mesh->ndfD[0], strIa, strJa
               , strAd, strA, sistEqD1);
      tm.dataStructD1 = getTimeC() - tm.dataStructD1;

      fprintf(fileLogExc, "Estrutuda montada.\n");
/*...................................................................*/

/*... dividindo a matriz*/
/*... Openmp(Vel,Pres)*/
      if (ompVar.fSolver) 
      {
        strcpy(str1, "thBeginD1");
        strcpy(str2, "thEndD1");
        strcpy(str3, "thSizeD1");
        strcpy(str4, "thHeightD1");
        pMatrixSolverOmp(m, sistEqD1, str1, str2, str3, str4);
      }
/*...................................................................*/
    }
/*...................................................................*/
  } while (nSistEq);
/*...................................................................*/

/*...*/
  if (opt->fItPlot && !mpiVar.myId) 
  {
    strcpy(auxName, preName);
    strcat(auxName, "_D1");
    fName(auxName, mpiVar.nPrcs, 0, 10, nameOut);
    opt->fileItPlot[FITPLOTD1] = openFile(nameOut, "w");
    fprintf(opt->fileItPlot[FITPLOTD1], "#D1\n#it ||b||/||b0|| ||b||\n");  
  }
/*...................................................................*/

/*... Openmp(Vel,Pres)*/
  if (ompVar.fSolver) 
  {
    if (*fSolvD1) 
    {
      nEqMax = sistEqD1->neqNov;
      HccaAlloc(DOUBLE, m, ompVar.buffer
               , nEqMax*ompVar.nThreadsSolver, "bufferOmp", false);
      ompVar.nEqMax = nEqMax;
      zero(ompVar.buffer, nEqMax*ompVar.nThreadsSolver, DOUBLEC);
      sistEqD1->omp.thY = ompVar.buffer;
    }
  }
/*...................................................................*/

/*... mapa de equacoes para comunicacao*/
//    if( mpiVar.nPrcs > 1) {    
//      front(&m,pMesh,sistEqT1,mesh->ndfT[0]);  
//    } 
/*...................................................................*/

/*... informacao da memoria total usada*/
  if (!mpiVar.myId) 
  {
    strcpy(str1, "MB");
    memoriaTotal(str1);
    usoMemoria(m, str1);
  }
/*...................................................................*/

}
/**********************************************************************/

/*********************************************************************
* Data de criacao    : 02/06/2018                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* readSolvTrans: leitura das configuracoes do solvers utilizados    *
* nos problemas de transporte                                       *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
void readSolvTrans(Memoria *m  , Mesh *mesh      , Reord *reordMesh
                , Solv *solvT1 , SistEq* sistEqT1, bool *fSolvT1
                , char* auxName, char* preName   , char* nameOut
                , FILE *fileIn , FileOpt *opt)
{

  unsigned short nSistEq;

  INT nEqMax;
  /*... Estrutura de dados*/
  char strIa[MNOMEPONTEIRO], strJa[MNOMEPONTEIRO];
  char strA[MNOMEPONTEIRO], strAd[MNOMEPONTEIRO];

  char str1[100], str2[100], str3[100], str4[100];

  char word[WORD_SIZE];

/*... solver*/
  readMacro(fileIn, word, false);
  nSistEq = (short)atol(word);
  do
  {
    readMacro(fileIn, word, false);
/*... equaocao de difusao d1*/
    if (!strcmp(word, "T1") || !strcmp(word, "t1"))
    {
      nSistEq--;
      *fSolvT1         = true;
      solvT1->solver   = PBICGSTAB;
      solvT1->tol      = smachn();
      solvT1->maxIt    = 50000;
      solvT1->fileSolv = NULL;
      solvT1->log      = true;
      solvT1->flag     = true;
/*...................................................................*/

/*...*/
      if (solvT1->log && !mpiVar.myId)
      {
        strcpy(auxName, preName);
        strcat(auxName, "_T1");
        fName(auxName, mpiVar.nPrcs, 0, 11, nameOut);
        solvT1->fileSolv = openFile(nameOut, "w");
      }
/*...................................................................*/

/*...*/
      sistEqT1->storage = CSRD;
      sistEqT1->unsym   = true;
/*...................................................................*/

/*... solver*/
      readMacro(fileIn, word, false);
      setSolverConfig(word, solvT1, &sistEqT1->storage);
/*...................................................................*/

/*... numeracao das equacoes das velocidades*/
      HccaAlloc(INT, m, sistEqT1->id, mesh->numel*mesh->ndfT[0]
              , "sistT1id", _AD_);

      fprintf(fileLogExc, "%s\n", DIF);
      fprintf(fileLogExc, "Numerando as equacoes.\n");

      tm.numeqT1 = getTimeC() - tm.numeqT1;
      sistEqT1->neq = numeq(sistEqT1->id     , reordMesh->num
                          , mesh->elm.faceRt1, mesh->elm.adj.nViz
                          , mesh->numel      , mesh->maxViz
                          , mesh->ndfT[0]);
      tm.numeqT1 = getTimeC() - tm.numeqT1;

      fprintf(fileLogExc, "Equacoes numeradas.\n");
      fprintf(fileLogExc, "%s\n", DIF);
/*...................................................................*/

/*...*/
      if (mpiVar.nPrcs > 1)
      {
        //          tm.numeqPres = getTimeC() - tm.numeqPres;
        //      sistEqT1->neqNov = countEq(reordMesh->num
        //                          ,mesh->elm.faceRt1  ,mesh->elm.adj.nViz
        //                          ,mesh->numelNov     ,mesh->maxViz
        //                          ,mesh->ndfT[0]);
        //          tm.numeqPres = getTimeC() - tm.numeqPres;
      }
      else
        sistEqT1->neqNov = sistEqT1->neq;
/*...................................................................*/

/*...*/
      HccaAlloc(DOUBLE        , m        , sistEqT1->b0
              , sistEqT1->neq, "sistT1b0", _AD_);
      HccaAlloc(DOUBLE      , m         , sistEqT1->b
             , sistEqT1->neq, "sistT1b ", _AD_);
      HccaAlloc(DOUBLE       , m         , sistEqT1->x
              , sistEqT1->neq, "sistT1x ", _AD_);
      zero(sistEqT1->b0, sistEqT1->neq, DOUBLEC);
      zero(sistEqT1->b , sistEqT1->neq, DOUBLEC);
      zero(sistEqT1->x , sistEqT1->neq, DOUBLEC);
/*...................................................................*/

/*... Estrutura de dados velocidades*/
      strcpy(strIa, "iaT1");
      strcpy(strJa, "jaT1");
      strcpy(strAd, "adT1");
      strcpy(strA , "aT1");

      fprintf(fileLogExc, "T1:\n");
      fprintf(fileLogExc, "Montagem da estrura de dados esparsa.\n");

      tm.dataStructT1 = getTimeC() - tm.dataStructT1;
      dataStruct(m               
               , sistEqT1->id      , reordMesh->num, mesh->elm.adj.nelcon
               , mesh->elm.adj.nViz, mesh->numelNov, mesh->maxViz
               , mesh->ndfT[0]     , strIa         , strJa
               , strAd             , strA          , sistEqT1);
      tm.dataStructT1 = getTimeC() - tm.dataStructT1;

      fprintf(fileLogExc, "Estrutuda montada.\n");
/*...................................................................*/

/*... dividindo a matriz*/
/*... Openmp(T1)*/
      if (ompVar.fSolver)
      {
        strcpy(str1, "thBeginT1");
        strcpy(str2, "thEndT1");
        strcpy(str3, "thSizeT1");
        strcpy(str4, "thHeightT1");
        pMatrixSolverOmp(m, sistEqT1, str1, str2, str3, str4);
      }
/*...................................................................*/
    }
/*...................................................................*/
  } while (nSistEq);
/*...................................................................*/

  /*...*/
  if (opt->fItPlot && !mpiVar.myId)
  {
    strcpy(auxName, preName);
    strcat(auxName, "_T1");
    fName(auxName, mpiVar.nPrcs, 0, 10, nameOut);
    opt->fileItPlot[FITPLOTT1] = openFile(nameOut, "w");
    fprintf(opt->fileItPlot[FITPLOTT1], "#T1\n#it ||b||/||b0|| ||b||\n");
  }
  /*...................................................................*/

/*... Openmp(Vel,Pres)*/
  if (ompVar.fSolver) {
    if (*fSolvT1) {
      nEqMax = sistEqT1->neqNov;
      HccaAlloc(DOUBLE, m, ompVar.buffer
              , nEqMax*ompVar.nThreadsSolver, "bufferOmp", false);
      ompVar.nEqMax = nEqMax;
      zero(ompVar.buffer, nEqMax*ompVar.nThreadsSolver, DOUBLEC);
      sistEqT1->omp.thY = ompVar.buffer;
    }
  }
/*...................................................................*/

/*... mapa de equacoes para comunicacao*/
  //    if( mpiVar.nPrcs > 1) {    
  //      front(&m,pMesh,sistEqT1,mesh->ndfT[0]);  
  //    } 
/*...................................................................*/

/*... informacao da memoria total usada*/
  if (!mpiVar.myId) 
  {
    strcpy(str1, "MB");
    memoriaTotal(str1);
    usoMemoria(m, str1);
  }
/*...................................................................*/

}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 29/01/2018                                    *
 * Data de modificaco : 02/02/2018                                    *
 *--------------------------------------------------------------------* 
 * readMean:                                                          *
 *--------------------------------------------------------------------* 
 * Parametros de entrada:                                             * 
 *--------------------------------------------------------------------* 
 * file    -> arquivo de arquivo                                      * 
 *--------------------------------------------------------------------* 
 * Parametros de saida:                                               * 
 *--------------------------------------------------------------------* 
 *--------------------------------------------------------------------* 
 * OBS:                                                               * 
 *--------------------------------------------------------------------*
 **********************************************************************/
void readMean(Memoria *m, FILE *fileIn
            , Mesh *mesh, Mean *media) {

  char word[WORD_SIZE];
  unsigned short nTerm,ndfVel;
  INT nel=mesh->numel;
/*...*/
  ndfVel = max(mesh->ndfF - 1,mesh->ndfFt - 2);
/*...................................................................*/

/*...*/
  media->t0 = 0.e0;
  fscanf(fileIn,"%d %d",&media->startSample,&media->endSample);
  fprintf(fileLogExc,"%-20s: (%d-%d)\n", "Samples"
                                       , media->startSample 
                                       , media->endSample);    
/*...................................................................*/

/*...*/
  media->fMedia = true;
  media->fInit  = false;
/*... tecnica de adveccao*/
  readMacro(fileIn, word, false);
  nTerm = (short) atol(word);
  do {
    readMacro(fileIn, word, false);
/*... velocidade*/
    if (!strcmp(word, "Vel") || !strcmp(word, "vel")) {
      fprintf(fileLogExc,"%-20s: enable\n", word);
      media->fVel = true;
      nTerm--;
/*... media das velcidades*/
      HccaAlloc(DOUBLE,m,media->mVel,nel*ndfVel,"medVel",_AD_);
      zero(media->mVel, nel*ndfVel, DOUBLEC);
/*... integral acumulada da media*/
      HccaAlloc(DOUBLE,m,media->sVel,nel*ndfVel,"medsVel",_AD_);
      zero(media->sVel, nel*ndfVel, DOUBLEC);      
    }
/*...................................................................*/

  } while (nTerm);
/*....................................................................*/

}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 05/05/2018                                    *
 * Data de modificaco : 19/09/2019                                    *
 *--------------------------------------------------------------------*
 * setReGrad:                                                         *
 *--------------------------------------------------------------------*
 * Parametros de entrada:                                             *
 *--------------------------------------------------------------------*
 * rcGrad  -> nao definido                                            *
 * file    -> arquivo de arquivo                                      *
 *--------------------------------------------------------------------*
 * Parametros de saida:                                               *
 *--------------------------------------------------------------------*
 * rcGrad  -> tipo de tecnica de reconstrucao de gradiente            *
 *--------------------------------------------------------------------*
 * OBS:                                                               *
 *--------------------------------------------------------------------*
 **********************************************************************/
void setReGrad(RcGrad *rcGrad, FILE *file)
{
  char word[WORD_SIZE];
  char macro[][WORD_SIZE] = { "greengausscell"  ,"greengaussnode"
                             ,"leastsquare"     ,"leastsquareqr"
                             ,"limiter"};
  char name[][WORD_SIZE] = {"barth","barthmod"};
  short i;
/*...*/
  readMacroV2(file, word, false, true);
/*.....................................................................*/

/*...*/
  if(!strcmp(word, macro[0])) 
  {
    rcGrad->type = RCGRADGAUSSC;
    fprintf(fileLogExc, "%-20s: %s\n", "Gradient", "GreenGaussCell");
  }
/*.....................................................................*/

/*...*/
  else if(!strcmp(word, macro[1]))
  {
    rcGrad->type = RCGRADGAUSSN;
    fprintf(fileLogExc, "%-20s: %s\n", "Gradient", "GreenGaussNode");
  }
/*.....................................................................*/

/*...*/
  else if(!strcmp(word, macro[2]))
  {
    rcGrad->type = RCLSQUARE;
    fprintf(fileLogExc, "%-20s: %s\n", "Gradient", "LeatSquare");
  }
/*.....................................................................*/

/*...*/
  else if(!strcmp(word, macro[3]))
  {
    rcGrad->type = RCLSQUAREQR;
    fprintf(fileLogExc, "%-20s: %s\n", "Gradient", "LeatSquareQR");
  }
/*.....................................................................*/
  
  rcGrad->fLimiter = false;
/*...*/
  readMacroV2(file, word, false, true);
/*.....................................................................*/
  
/*...*/
  if(!strcmp(word, macro[4])) 
  {
    rcGrad->fLimiter = true;
    fscanf(file,"%lf",&rcGrad->beta);
    readMacroV2(file, word, false, true);
    for(i=0;i<2;i++)
      if(!strcmp(word, name[i]))
        rcGrad->func = i;
      
    fprintf(fileLogExc, "%-20s: %s beta %lf\n", "Limiter"
                                         , name[rcGrad->func]
                                         , rcGrad->beta);
  }
/*.....................................................................*/


}
/**********************************************************************/

/**********************************************************************
* Data de criacao    : 19/05/2018                                    *
* Data de modificaco : 27/09/2019                                    *
*--------------------------------------------------------------------*
* setReGrad:                                                         *
*--------------------------------------------------------------------*
* Parametros de entrada:                                             *
*--------------------------------------------------------------------*
* prop    -> nao definido                                            *
* type    -> tipo vo de arquivo                                      *
*--------------------------------------------------------------------*
* Parametros de saida:                                               *
*--------------------------------------------------------------------*
* rcGrad  -> tipo de tecnica de reconstrucao de gradiente            *
*--------------------------------------------------------------------*
* OBS:                                                               *
*--------------------------------------------------------------------*
**********************************************************************/
void readFileMat(DOUBLE *prop, short *type, short numat,FILE *file)
{
  FILE *fileIn = NULL;
  char str[] = {"end"};
  char word[WORD_SIZE];
  char macro[][WORD_SIZE] = { "celltype"                       /* 0  */
                            ,"densityd1"   ,"coefdiffd1"       /* 1, 2*/   
                            ,"densityt1"   ,"coefdifft1"       /* 3, 4*/
                            ,"densityfluid","dviscosity"       /* 5, 6*/
                            ,"ctherm"      ,"cp"               /* 7, 8*/
                            ,"mmolar"        };                /* 9,*/
  short i, j,iMat;
  int aux; 
  DOUBLE v;

  for(i=0;i<numat;i++)
  {
    fscanf(file, "%hd", &iMat);
    iMat--;    
    readMacro(file, word, false);
    fileIn = openFile(word,"r");
/*...*/
    readMacro(fileIn, word, false);
    do
    {
      convStringLower(word);
      j = 0;
/*... cellType*/
      if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%d", &aux);
        type[iMat] = (short) aux;
        fprintf(fileLogExc, "%-20s: %d\n",macro[j-1], type[iMat]);
      }
/*.....................................................................*/

/*... densityd1*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn,"%lf",&v);
        prop[DENSITY] = v;
        fprintf(fileLogExc, "%-20s: %e\n", macro[j-1], prop[DENSITY]);
      }
/*.....................................................................*/

/*... coefdiffd1*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%lf", &v);
        prop[COEFDIF] = v;
        fprintf(fileLogExc, "%-20s: %e\n", macro[j-1], prop[COEFDIF]);
      }
/*.....................................................................*/

/*... densityt1*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%lf", &v);
        prop[DENSITY] = v;
        fprintf(fileLogExc, "%-20s: %e\n", macro[j - 1], prop[DENSITY]);
      }
/*.....................................................................*/

/*... coefdifft1*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%lf", &v);
        prop[COEFDIF] = v;
        fprintf(fileLogExc, "%-20s: %e\n", macro[j - 1], prop[COEFDIF]);
      }
/*.....................................................................*/

/*... densitytFluid*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%lf", &v);
        prop[DENSITY] = v;
        fprintf(fileLogExc, "%-20s: %e\n", macro[j - 1], prop[DENSITY]);
      }
/*.....................................................................*/

/*... viscosidade dinamica*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%lf", &v);
        prop[DYNAMICVISCOSITY ] = v;
        fprintf(fileLogExc, "%-20s: %e\n", macro[j - 1], prop[DYNAMICVISCOSITY]);
      }
/*.....................................................................*/

/*... condutividade termica*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%lf", &v);
        prop[THERMALCONDUCTIVITY] = v;
        fprintf(fileLogExc, "%-20s: %e\n", macro[j - 1], prop[THERMALCONDUCTIVITY]);
      }
/*.....................................................................*/

/*... condutividade termica*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%lf", &v);
        prop[SPECIFICHEATCAPACITYFLUID] = v;
        fprintf(fileLogExc, "%-20s: %e\n", macro[j - 1], prop[SPECIFICHEATCAPACITYFLUID]);
      }
/*.....................................................................*/

/*... massa molar*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%lf", &v);
        prop[MMOLARMASS] = v;
        fprintf(fileLogExc, "%-20s: %e\n", macro[j - 1], prop[MMOLARMASS]);
      }
/*.....................................................................*/
      readMacro(fileIn, word, false);
    }while(strcmp(word,str));
/*.....................................................................*/
    fclose(fileIn);
  }
/*.....................................................................*/
}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 20/07/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * readChemical : leitura do sistema quimico                         *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * c       ->                                                        *
 * file    -> arquivo de arquivo                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void readChemical(Combustion *c, FILE *file)
{
  FILE *fileAux;
  char word[WORD_SIZE];
  char macros[][WORD_SIZE] =
       { "elements"       ,"species"   ,"reactions"  /* 0, 1, 2*/
       , "fuel"           ,"ox"        ,"product"};  /* 3, 4, 5*/

  short id;
  int i,j,n,nn;
  DOUBLE value,ru;

  readMacroV2(file, word, false, true);
  fileAux = openFile(word, "r");

 /*... incializando com zero*/
  for(i=0;i<MAXREAC;i++)
  {
    c->chem.reac[i].reverse = false;
    for(j=0;j<MAXSPECIES;j++)
    {
      c->chem.reac[i].stch[0][j] = 0.e0;
      c->chem.reac[i].stch[1][j] = 0.e0;
      c->chem.reac[i].exp[0][j] = 0.e0;
      c->chem.reac[i].exp[1][j] = 0.e0;
    }
  }
  c->chem.sH2O  = -1;
  c->chem.sO2   = -1;
  c->chem.sN2   = -1;
  c->chem.sCO2  = -1;
  c->chem.sCO   = -1;
  c->chem.sCH4  = -1;
/*.....................................................................*/

/*...*/
  do
  {
/*... lendo*/
    readMacroV2(fileAux, word, false, true);
/*... elements*/
    if(!strcmp(word,macros[0]))    
    {
      fscanf(fileAux, "%d", &n);
      c->chem.nEp = n;
      for(i=0;i<n;i++)
      {
        readMacroV2(fileAux, word, false, true);  
/*... C*/
        if(!strcmp(word,"c"))
        {
          c->chem.eC    = i; 
          fscanf(fileAux,"%lf",&value);
          c->chem.mE[i] = value;
        } 
/*... H*/
        else if(!strcmp(word,"h"))
        {
          c->chem.eH    = i; 
          fscanf(fileAux,"%lf",&value);
          c->chem.mE[i] = value;
        }   
/*... O*/
        else if(!strcmp(word,"o"))
        {
          c->chem.eO    = i; 
          fscanf(fileAux,"%lf",&value);
          c->chem.mE[i] = value;
        }  
/*... N*/
        else if(!strcmp(word,"n"))
        {
          c->chem.eN    = i; 
          fscanf(fileAux,"%lf",&value);
          c->chem.mE[i] = value;
        }  
      }
    }
/*...................................................................*/

/*... species*/ 
    else if(!strcmp(word,macros[1]))
    {
      fscanf(fileAux, "%d", &n);
      c->chem.nSp = n;
      for(i=0;i<n;i++)
      {
        readMacroV2(fileAux, word, false, true);  
/*... O2*/
        if(!strcmp(word,"o2"))
          c->chem.sO2 = i;           
/*... H2O*/
        else if(!strcmp(word,"h2o"))
          c->chem.sH2O = i;           
/*... CO2*/
        else if(!strcmp(word,"co2"))
          c->chem.sCO2 = i; 
/*... CO*/
        else if(!strcmp(word,"co"))
          c->chem.sCO = i; 
/*... N2*/
        else if(!strcmp(word,"n2"))
          c->chem.sN2 = i;
/*... CH4*/
        else if(!strcmp(word,"ch4"))
          c->chem.sCH4 = i;
/*... C3H8*/
        else if(!strcmp(word,"c3h8"))
          c->chem.sC3H8 = i;
/*...*/
        strcpy(c->chem.sp[i].name,word);
        fscanf(fileAux,"%hd",&c->chem.sp[i].nC);
        fscanf(fileAux,"%hd",&c->chem.sp[i].nH);
        fscanf(fileAux,"%hd",&c->chem.sp[i].nO);
        fscanf(fileAux,"%hd",&c->chem.sp[i].nN);
        fscanf(fileAux,"%lf",&c->chem.sp[i].leornadJones[0]);
        fscanf(fileAux,"%lf",&c->chem.sp[i].leornadJones[1]);
/*...................................................................*/
      }
    }
/*...................................................................*/

/*... reaction*/
    else if(!strcmp(word,macros[2]))
    {
      fscanf(fileAux, "%d", &n);
      c->chem.nReac = n;
      for(i=0;i<n;i++)
      {
        do
        {
          readMacroV2(fileAux, word, false, true);  
/*... */
          if (!strcmp(word, "r"))
          {
            fscanf(fileAux,"%d",&nn);
            c->chem.reac[i].nPartSp[0] = nn;
            for(j=0;j<nn;j++)
            {
              fscanf(fileAux,"%lf %s",&value,word);
              id = searchSpeciesId(&c->chem,word);
              c->chem.reac[i].stch[0][id]  = value;
              c->chem.reac[i].partSp[0][j] = id;
            }
          }
/*...................................................................*/

/*... */
          else if (!strcmp(word, "p"))
          {
            fscanf(fileAux,"%d",&nn);
            c->chem.reac[i].nPartSp[1] = nn;
            for(j=0;j<nn;j++)
            {
              fscanf(fileAux,"%lf %s",&value,word);
              id = searchSpeciesId(&c->chem,word);
              c->chem.reac[i].stch[1][id] = value;
              c->chem.reac[i].partSp[1][j] = id;
            }
          }
/*...................................................................*/

/*... */
          else if (!strcmp(word,"arrd"))
          {
            fscanf(fileAux,"%lf %lf %lf",&c->chem.reac[i].ArrF.A
                                        ,&c->chem.reac[i].ArrF.beta
                                        ,&c->chem.reac[i].ArrF.E);
            
/*... cal/(mol*kelvin)*/
            ru = IDEALGASRC;
/*J/(mol.kelvin) */
//          ru = IDEALGASR*1.e-03;
            c->chem.reac[i].ArrF.Ta = c->chem.reac[i].ArrF.E/ru;

            fscanf(fileAux,"%d",&nn);
            for(j=0;j<nn;j++)
            {
              fscanf(fileAux,"%s %lf",word,&value);
              id = searchSpeciesId(&c->chem,word);
              c->chem.reac[i].exp[0][id] = value;
            }
          }
/*...................................................................*/

/*... */
          else if (!strcmp(word,"arrr"))
          {
            c->chem.reac[i].reverse = true;
            fscanf(fileAux,"%lf %lf %lf",&c->chem.reac[i].ArrR.A
                                        ,&c->chem.reac[i].ArrR.beta
                                        ,&c->chem.reac[i].ArrR.E);
/*... cal/(mol*kelvin)*/
            ru = IDEALGASRC;
/*J/(mol.kelvin) */
//          ru = IDEALGASR*1.e-03;
            c->chem.reac[i].ArrR.Ta = c->chem.reac[i].ArrR.E/ru;

            fscanf(fileAux,"%d",&nn);
            for(j=0;j<nn;j++)
            {
              fscanf(fileAux,"%s %lf",word,&value);
              id = searchSpeciesId(&c->chem,word);
              c->chem.reac[i].exp[1][id] = value;
            }
/*...................................................................*/
          }
/*...................................................................*/
        }while(strcmp(word,"end"));
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... fuel*/ 
    else if(!strcmp(word,macros[3]))
    {
      fscanf(fileAux, "%d", &n);
      c->chem.nFuel = n;
      for(i=0;i<n;i++)
      {
        readMacroV2(fileAux, word, false, true);  
        c->chem.fuel[i] = searchSpeciesId(&c->chem,word);
      }
    }
/*...................................................................*/

/*... ox*/ 
    else if(!strcmp(word,macros[4]))
    {
      fscanf(fileAux, "%d", &n);
      c->chem.nOx = n;
      for(i=0;i<n;i++)
      {
        readMacroV2(fileAux, word, false, true);  
        c->chem.ox[i] = searchSpeciesId(&c->chem,word);
      }
    }
/*...................................................................*/

/*... product*/ 
    else if(!strcmp(word,macros[5]))
    {
      fscanf(fileAux, "%d", &n);
      c->chem.nProd = n;
      for(i=0;i<n;i++)
      {
        readMacroV2(fileAux, word, false, true);  
        c->chem.prod[i] = searchSpeciesId(&c->chem,word);
      }
    }
/*...................................................................*/

/*...................................................................*/
  }while(strcmp(word,"end"));
/*...................................................................*/

/*...*/
  initMolarMass(c);
/*...................................................................*/

/*...*/
  for(i=0;i<MAXREAC;i++)
    for(j=0;j<MAXSPECIES;j++)
      c->chem.reac[i].stch[2][j] = c->chem.reac[i].stch[1][j]
                                 - c->chem.reac[i].stch[0][j];
/*.....................................................................*/

/*...*/
  fprintf(fileLogExc, "%-20s:\n","Elements");
  fprintf(fileLogExc, "%-20s: %d %lf\n","C",c->chem.eC
                                           ,c->chem.mE[c->chem.eC]); 
  fprintf(fileLogExc, "%-20s: %d %lf\n","O",c->chem.eO
                                           ,c->chem.mE[c->chem.eO]);
  fprintf(fileLogExc, "%-20s: %d %lf\n","H",c->chem.eH
                                           ,c->chem.mE[c->chem.eH]);
  fprintf(fileLogExc, "%-20s: %d %lf\n","N",c->chem.eN
                                           ,c->chem.mE[c->chem.eN]);
/*...................................................................*/

/*...*/
  fprintf(fileLogExc, "\n");
  fprintf(fileLogExc, "%-20s:\n","Species");
  for(i=0;i<c->chem.nSp;i++)
    fprintf(fileLogExc, "%-20s: %d C%dH%dO%dN%d %lf\n",c->chem.sp[i].name
                                     ,i
                                     ,c->chem.sp[i].nC,c->chem.sp[i].nH
                                     ,c->chem.sp[i].nO,c->chem.sp[i].nN
                                     ,c->chem.sp[i].mW); 
/*...................................................................*/

/*...*/
  fprintf(fileLogExc, "\n");
  fprintf(fileLogExc, "%-20s:\n","Reactions");
/*...*/
  for(i=0;i<c->chem.nReac;i++)
  {
    fprintf(fileLogExc, "%d ",i);
    for(j=0;j<c->chem.nSp;j++)
    {
      value = c->chem.reac[i].stch[0][j];
      if(value != 0.e0) 
        fprintf(fileLogExc, " %lf %s ",value,c->chem.sp[j].name);
    }
    if(c->chem.reac[i].reverse) fprintf(fileLogExc, " <=> ");
    else fprintf(fileLogExc, " => ");
    for(j=0;j<c->chem.nSp;j++)
    {
      value = c->chem.reac[i].stch[1][j];
      if(value != 0.e0)
        fprintf(fileLogExc, " %lf %s ",value,c->chem.sp[j].name);
      
    }
/*...................................................................*/

/*... Lei de arrhenius da reacao direta*/
    fprintf(fileLogExc, "arrD %e %e %e %e",c->chem.reac[i].ArrF.A
                                           ,c->chem.reac[i].ArrF.beta
                                           ,c->chem.reac[i].ArrF.E
                                           ,c->chem.reac[i].ArrF.Ta); 
    for(j=0;j<c->chem.nSp;j++)
    {
      value = c->chem.reac[i].exp[0][j];
      if(value != 0.e0)
        fprintf(fileLogExc, " %s  %lf ",c->chem.sp[j].name,value);      
    }
/*...................................................................*/

/*... Lei de arrhenius da reacao inversa*/
    if(c->chem.reac[i].reverse)
    {
      fprintf(fileLogExc, "arrR %e %e %e %e",c->chem.reac[i].ArrR.A
                                       ,c->chem.reac[i].ArrR.beta
                                       ,c->chem.reac[i].ArrR.E
                                       ,c->chem.reac[i].ArrR.Ta);  
      for(j=0;j<c->chem.nSp;j++)
      {
        value = c->chem.reac[i].exp[1][j];
        if(value != 0.e0)
          fprintf(fileLogExc, " %s  %lf ",c->chem.sp[j].name,value);
      }
    }
/*...................................................................*/
    fprintf(fileLogExc,"\n");
  }                                 
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 11/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * readResidual : leitura dos paramentros do modelo de combustao     *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * c       ->                                                        *
 * file    -> arquivo de arquivo                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void readResidual(Simple *sc, FILE *file)
{

  char *str = { "endresidual" };
  char word[WORD_SIZE];

  readMacroV2(file, word, false, true);
  do 
  {

/*... vel*/
      if (!strcmp(word, "vel"))
      {
        fscanf(file, "%lf %lf %lf %hd", &sc->vel.tol[0]
                                      , &sc->vel.tol[1]
                                      , &sc->vel.tol[2]
                                      , &sc->vel.k);
        readMacroV2(file, word, false, true);
        sc->vel.fRel = readBool(word);
        readMacroV2(file, word, false, true);
        typeResidual(word,&sc->vel);

        fprintf(fileLogExc, "%-20s:\n"   ,"Vel");
        fprintf(fileLogExc, "%-20s: %e\n","tolU1",sc->vel.tol[0]);
        fprintf(fileLogExc, "%-20s: %e\n","tolU2",sc->vel.tol[1]);
        fprintf(fileLogExc, "%-20s: %e\n","tolU3",sc->vel.tol[2]);
        fprintf(fileLogExc, "%-20s: %hd\n","k"    ,sc->vel.k);
        if(sc->vel.fRel)
          fprintf(fileLogExc, "%-20s: %s\n","fRelative","true"); 
        else
          fprintf(fileLogExc, "%-20s: %s\n","fRelative","false");
        fprintf(fileLogExc, "%-20s: %s\n\n","type"       ,sc->vel.name);
      }
/*.....................................................................*/

/*... mass*/
      else if (!strcmp(word, "mass"))
      {
        fscanf(file, "%lf %hd", &sc->mass.tol[0]
                              , &sc->mass.k);
        readMacroV2(file, word, false, true);
        sc->mass.fRel = readBool(word);
        readMacroV2(file, word, false, true);
        typeResidual(word,&sc->mass);

        fprintf(fileLogExc, "%-20s:\n"   ,"Mass");
        fprintf(fileLogExc, "%-20s: %e\n","tolMass",sc->mass.tol[0]);
        fprintf(fileLogExc, "%-20s: %hd\n","k"      ,sc->mass.k);
        if(sc->vel.fRel)
          fprintf(fileLogExc, "%-20s: %s\n","fRelative","true"); 
        else
          fprintf(fileLogExc, "%-20s: %s\n","fRelative","false");
        fprintf(fileLogExc, "%-20s: %s\n\n","type"       ,sc->mass.name);
      }
/*.....................................................................*/

/*... energy*/
      else if (!strcmp(word, "energy"))
      {
        fscanf(file, "%lf %hd", &sc->energy.tol[0]
                              , &sc->energy.k);
        readMacroV2(file, word, false, true);
        sc->energy.fRel = readBool(word);
        readMacroV2(file, word, false, true);
        typeResidual(word,&sc->energy);

        fprintf(fileLogExc, "%-20s:\n"   ,"Energy");
        fprintf(fileLogExc, "%-20s: %e\n","tolMass",sc->energy.tol[0]);
        fprintf(fileLogExc, "%-20s: %hd\n","k"      ,sc->energy.k);
        if(sc->vel.fRel)
          fprintf(fileLogExc, "%-20s: %s\n","fRelative","true"); 
        else
          fprintf(fileLogExc, "%-20s: %s\n","fRelative","false");
        fprintf(fileLogExc, "%-20s: %s\n\n","type"       ,sc->energy.name);
      }
/*.....................................................................*/

/*... energy*/
      else if (!strcmp(word, "z"))
      {
        fscanf(file, "%lf %hd", &sc->z.tol[0]
                              , &sc->z.k);
        readMacroV2(file, word, false, true);
        sc->z.fRel = readBool(word);
        readMacroV2(file, word, false, true);
        typeResidual(word,&sc->z);

        fprintf(fileLogExc, "%-20s:\n"   ,"Z");
        fprintf(fileLogExc, "%-20s: %e\n","tolZ",sc->z.tol[0]);
        fprintf(fileLogExc, "%-20s: %hd\n","k"      ,sc->z.k);
        if(sc->vel.fRel)
          fprintf(fileLogExc, "%-20s: %s\n","fRelative","true"); 
        else
          fprintf(fileLogExc, "%-20s: %s\n","fRelative","false");
        fprintf(fileLogExc, "%-20s: %s\n\n","type"       ,sc->z.name);
      }
/*.....................................................................*/

    readMacroV2(file, word, false, true);
  } while (strcmp(word, str));

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 13/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * readBool:                                                         *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
bool readBool(char *word)
{

  if( word[0] == 'F' || word[0] == 'f')
    return false;
  else if( word[0] == 'T' || word[0] == 'f')
    return true;
  
  return -1;
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 13/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * typeResidual:                                                     *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * word  ->                                                          *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void typeResidual(char *word,Residual *re)
{
  char macros[][WORD_SIZE] =
             {"rscaledsummax","rscaledsum"   /* 0, 1*/
             ,"rscaled"      ,"rsqrt"     }; /* 2, 3*/  

  if(!strcmp(word, macros[0]))
  {
    strcpy(re->name,word);
    re->type = RSCALEDSUMMAX;
  } 
  else if(!strcmp(word, macros[1]))
  {
    strcpy(re->name,word);
    re->type = RSCALEDSUM;
  }
  else if(!strcmp(word, macros[2]))
  {
   strcpy(re->name,word);
    re->type = RSCALED;
  }  
  else if(!strcmp(word, macros[3]))
  {
    strcpy(re->name,word);
    re->type = RSQRT;
  }
  else
  {
    ERRO_OP_WORD(fileLogDebug,__FILE__,__func__,__LINE__,"Residual"
                 ,word,EXIT_FAILURE);
  }
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 26/08/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * typeResidual:                                                     *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void configEdo(Edo *edo, FILE *file)
{
  char *str = { "endedo" };
  char word[WORD_SIZE];

  readMacroV2(file, word, false, true);
  do 
  {
/*... vel*/
    if (!strcmp(word, "sie"))
    {
      edo->type = EDO_SIE;
      fscanf(file, "%d %lf", &edo->maxIt,&edo->tol);
      fprintf(fileLogExc, "%-20s:\n"   ,"Sie");
      fprintf(fileLogExc, "%-20s: %d\n","Maxit",edo->maxIt);
      fprintf(fileLogExc, "%-20s: %e\n","tol"  ,edo->tol);
    }
/*.....................................................................*/
    readMacroV2(file, word, false, true);
  } while (strcmp(word, str));

}
/*********************************************************************/

/********************************************************************
* Data de criacao    : 00/00/0000                                   *
* Data de modificaco : 16/10/2019                                   *
*-------------------------------------------------------------------*
* CONVLOADSPRESC : Converte condicoes de contorno da Temperatura    *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
static void setLoadsPresC(Loads *loadsPres,Loads *loadsPresC
                         ,short const iCod)
{

  short i,j,type;

  for(i=0;i<MAXLOAD;i++){
    loadsPresC[i].type = loadsPres[i].type;
    loadsPresC[i].np   = loadsPres[i].np;
    for(j=0;j<MAXLOADPARAMETER;j++){
      loadsPresC[i].par[j] = loadsPres[i].par[j];
    }

    type = loadsPresC[i].type;
    if(type == DIRICHLETBC)
      loadsPresC[i].par[0] = 0.e0;
    else if(type == INLETTOTALPRES)
      loadsPresC[i].par[1] = 0.e0;
    else if (type == FLUXPRES)
    {
      loadsPresC[i].type    = FLUXPRESC;
      loadsPres[i].iCod[0]  = iCod;
      loadsPresC[i].iCod[0] = iCod;
    }
  }
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 09/05/2019                                   *
 *-------------------------------------------------------------------*
 * setLoadsEnergy:                                                   *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
static void setLoadsEnergy(PropVarFluid *pF
                           ,Loads *loadsEnergy   ,Loads *loadsTemp
                           ,Loads *loadsVel      ,DOUBLE *RESTRICT prop
                           ,bool const fTemp     ,bool const iKelvin   )
{
  bool fDensity = pF->fDensity, fSheat = pF->fSpecificHeat;
  short i,j,type;
  DOUBLE t,sHeat,tmp,tmp1;
  
/*... cc da equacao da energia e em temperatura*/
  if(fTemp){

    for(i=0;i<MAXLOAD;i++){
      loadsEnergy[i].fUse    = loadsTemp[i].fUse;
      loadsEnergy[i].type    = loadsTemp[i].type;
      loadsEnergy[i].np      = loadsTemp[i].np;
      for(j=0;j<MAXLOADPARAMETER;j++)
        loadsEnergy[i].par[j] = loadsTemp[i].par[j];
    }
/*... converte c para kelvin*/    
//  if(iKelvin)
//    for(i=0;i<MAXLOADFLUID;i++){
//      type = loadsEnergy[i].type;
//      if( type == DIRICHLETBC ||  type == INLET 
//      ||  type == CONVECTIONHEAT)
//        loadsEnergy[i].par[0] 
//                       = CELSIUS_FOR_KELVIN(loadsEnergy[i].par[0]);
//    }
      
  }
/*....................................................................*/

  else{
    sHeat = MAT2D(0,SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
    for(i=0;i<MAXLOAD;i++){
      loadsEnergy[i].fUse    = loadsTemp[i].fUse;
      loadsEnergy[i].type    = loadsTemp[i].type;
      loadsEnergy[i].np      = loadsTemp[i].np;
      type = loadsTemp[i].type;
      for(j=0;j<MAXLOADPARAMETER;j++)
        loadsEnergy[i].par[j] = loadsTemp[i].par[j];

/*...*/
      if( type == DIRICHLETBC ){
        t        = loadsTemp[i].par[0];
        tmp = tempForSpecificEnthalpy(&pF->sHeat,t 
                                     , pF->sHeatRef
                                     , fSheat   , iKelvin);
        loadsEnergy[i].par[0] = tmp;               
      }
/*....................................................................*/

/*...*/
      else if ( type == INLET ||  type == OPEN) {
        t = loadsTemp[i].par[0];
        tmp = tempForSpecificEnthalpy(&pF->sHeat  ,t
                                    , pF->sHeatRef
                                    , fSheat    , iKelvin);       
        loadsEnergy[i].par[0] = tmp;
/*... densidade*/
        if(fDensity)
        {         
          tmp1 = airDensity(&pF->den          
                           ,t               ,thDynamic.pTh[2]
                           ,thDynamic.pTh[2],pF->molarMass
                           ,iKelvin);
          loadsEnergy[i].density = tmp1;
          loadsTemp[i].density   = tmp1;
        }       
        else
        {
          loadsEnergy[i].density = loadsVel[i].par[loadsVel[i].np-1];
          loadsTemp[i].density = loadsVel[i].par[loadsVel[i].np-1];
        }
/*....................................................................*/

/*... velocidades*/
        loadsEnergy[i].vel[0] = loadsVel[i].par[0];
        loadsEnergy[i].vel[1] = loadsVel[i].par[1]; 
        loadsEnergy[i].vel[2] = loadsVel[i].par[2];
   
        loadsTemp[i].vel[0] = loadsVel[i].par[0];
        loadsTemp[i].vel[1] = loadsVel[i].par[1]; 
        loadsTemp[i].vel[2] = loadsVel[i].par[2];
/*....................................................................*/
      }
/*....................................................................*/

/*...*/
      else if (type == NEUMANNBC  ||  type == CONVECTIONHEAT
           ||  type == OUTLET) {
        loadsEnergy[i].par[0] = loadsTemp[i].par[0];
      }
/*....................................................................*/
    }     
  }
/*....................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 09/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * setLoadsEnergy:                                                   *
  *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
static void setLoadsEnergyMix(Combustion *cModel   ,Prop *pDen
                              ,Prop *sHeatProp      ,Loads *loadsEnergy
                              ,Loads *loadsTemp     ,Loads *loadsZ     
                              ,Loads *loadsVel      ,DOUBLE *RESTRICT prop
                              ,bool const fTemp     ,bool const fSheat 
                              ,bool const iKelvin   ,bool const fDensity
                              ,bool const fGrouped)   
{
  short i,j,n,type;
  short ns = cModel->nOfSpecies, nl =cModel->nOfSpeciesLump;
  DOUBLE t,yFrac[MAXSPECIES],sHeat,tmp,tmp1,molarMassMix;
  
/*... cc da equacao da energia e em temperatura*/
  if(fTemp){

    for(i=0;i<MAXLOAD;i++){
      loadsEnergy[i].fUse    = loadsTemp[i].fUse;
      loadsEnergy[i].type    = loadsTemp[i].type;
      loadsEnergy[i].np      = loadsTemp[i].np;
      for(j=0;j<MAXLOADPARAMETER;j++)
        loadsEnergy[i].par[j] = loadsTemp[i].par[j];
    }
/*... converte c para kelvin*/    
//  if(iKelvin)
//    for(i=0;i<MAXLOADFLUID;i++){
//      type = loadsEnergy[i].type;
//      if( type == DIRICHLETBC ||  type == INLET 
//      ||  type == CONVECTIONHEAT)
//        loadsEnergy[i].par[0] 
//                       = CELSIUS_FOR_KELVIN(loadsEnergy[i].par[0]);
//    }
      
  }
/*....................................................................*/

  else{
    sHeat = MAT2D(0,SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
    for(i=0;i<MAXLOAD;i++){
      loadsEnergy[i].fUse      = loadsTemp[i].fUse;
      loadsEnergy[i].type    = loadsTemp[i].type;
      loadsEnergy[i].np      = loadsTemp[i].np;
      type = loadsTemp[i].type;
      for(j=0;j<MAXLOADPARAMETER;j++)
        loadsEnergy[i].par[j] = loadsTemp[i].par[j];

/*...*/
      if( type == DIRICHLETBC ){
        t        = loadsTemp[i].par[0];
        
        n        = loadsZ[i].np;
        getSpeciesPrimitivesCc(cModel,yFrac,loadsZ[i].par);

        tmp = tempToSpecificEnthalpyMix( sHeatProp, yFrac 
                                         , t        , sHeat
                                         , n
                                         , fSheat  , iKelvin);
        loadsEnergy[i].par[0] = tmp;               
      }
/*....................................................................*/

/*...*/
      else if ( type == INLET ||  type == OPEN) {
        t = loadsTemp[i].par[0];

        getSpeciesPrimitivesCc(cModel,yFrac,loadsZ[i].par);
        
        tmp = tempToSpecificEnthalpyMix( sHeatProp  , yFrac 
                                         , t        , sHeat
                                         , ns
                                         , fSheat  , iKelvin);
        loadsEnergy[i].par[0] = tmp;
/*... densidade*/
        if(fDensity)
        {
          molarMassMix =  mixtureMolarMass(cModel,yFrac); 
          tmp1 = mixtureSpeciesDensity(pDen           ,molarMassMix
                                     ,t               ,thDynamic.pTh[2]
                                     ,thDynamic.pTh[2],iKelvin);
          loadsEnergy[i].density = tmp1;
          loadsTemp[i].density   = tmp1;
        }       
        else
        {
          loadsEnergy[i].density = loadsVel[i].par[loadsVel[i].np-1];
          loadsTemp[i].density = loadsVel[i].par[loadsVel[i].np-1];
        }
/*....................................................................*/

/*... velocidades*/
        loadsEnergy[i].vel[0] = loadsVel[i].par[0];
        loadsEnergy[i].vel[1] = loadsVel[i].par[1]; 
        loadsEnergy[i].vel[2] = loadsVel[i].par[2];
   
        loadsTemp[i].vel[0] = loadsVel[i].par[0];
        loadsTemp[i].vel[1] = loadsVel[i].par[1]; 
        loadsTemp[i].vel[2] = loadsVel[i].par[2];
/*....................................................................*/
      }
/*....................................................................*/

/*...*/
      else if (type == NEUMANNBC  ||  type == CONVECTIONHEAT
           ||  type == OUTLET) {
        loadsEnergy[i].par[0] = loadsTemp[i].par[0];
      }
/*....................................................................*/
    }     
  }
/*....................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 08/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * setLoadsEnergy:                                                  *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
static void convLoadsZcombMix(Combustion *cModel   ,Prop *pDen
                             ,Prop *sHeatProp   ,Loads *loadsTemp    
                             ,Loads *loadsZ        ,Loads *loadsVel
                             ,DOUBLE *RESTRICT prop
                             ,bool const fTemp     ,bool const fSheat 
                             ,bool const iKelvin   ,bool const fDensity
                             ,bool const fGrouped)   
{
  short i,type;
  short ns = cModel->nOfSpecies, nl =cModel->nOfSpeciesLump;
  DOUBLE t,yFrac[MAXSPECIES],sHeat,tmp1,molarMassMix;
  
  sHeat = MAT2D(0,SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
  for(i=0;i<MAXLOAD;i++)
  {
    type = loadsZ[i].type;
/*...*/
    if ( type == INLET ||  type == OPEN) 
    {
      t = loadsTemp[i].par[0];

      getSpeciesPrimitivesCc(cModel,yFrac,loadsZ[i].par);

/*... densidade*/
      if(fDensity)
      {
        molarMassMix =  mixtureMolarMass(cModel,yFrac); 
        tmp1 = mixtureSpeciesDensity(pDen           ,molarMassMix
                                    ,t               ,thDynamic.pTh[2]
                                   ,thDynamic.pTh[2],iKelvin);
        loadsZ[i].density = tmp1;
      }
      else
        loadsZ[i].density = loadsVel[i].par[loadsVel[i].np-1];
/*....................................................................*/

/*... velocidade*/
      loadsZ[i].vel[0] = loadsVel[i].par[0];
      loadsZ[i].vel[1] = loadsVel[i].par[1]; 
      loadsZ[i].vel[2] = loadsVel[i].par[2];
/*....................................................................*/
    }     
  }
/*....................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 08/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * setLoadsEnergy: Converte condicoes de contorno da Temperatura    *
 * de C para kelvin                                                  *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
static void convLoadsVelMix(Combustion *cModel   ,Prop *pDen
                           ,Prop *sHeatProp   ,Loads *loadsVel
                           ,Loads *loadsTemp     ,Loads *loadsZ     
                           ,DOUBLE *RESTRICT prop
                           ,bool const fTemp     ,bool const fSheat 
                           ,bool const iKelvin   ,bool const fDensity
                           ,bool const fGrouped)   
{
  short i,type;
  short ns = cModel->nOfSpecies, nl =cModel->nOfSpeciesLump;
  DOUBLE t,yFrac[MAXSPECIES],sHeat,tmp1,molarMassMix;
  
  sHeat = MAT2D(0,SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
  for(i=0;i<MAXLOAD;i++)
  {
    type = loadsVel[i].type;
/*...*/
    if ( type == INLET ||  type == OPEN) 
    {
      t = loadsTemp[i].par[0];

      getSpeciesPrimitivesCc(cModel,yFrac,loadsZ[i].par);
      
      if(fDensity)
      {
        molarMassMix =  mixtureMolarMass(cModel,yFrac); 
        tmp1 = mixtureSpeciesDensity(pDen           ,molarMassMix
                                    ,t               ,thDynamic.pTh[2]
                                   ,thDynamic.pTh[2],iKelvin);
        loadsVel[i].par[loadsVel[i].np-1] = tmp1;  
      }
/*....................................................................*/
    }     
  }
/*....................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 10/05/2019                                   *
 * Data de modificaco : 16/10/2019                                   *
 *-------------------------------------------------------------------*
 * setLoadsVel:                                                      *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
static void setLoadsVel(PropVarFluid *propFluid
                        ,Loads *loadsVel      ,Loads *loadsTemp
                        ,DOUBLE *RESTRICT prop
                        ,bool const fTemp     ,bool const iKelvin)

{

  bool fDensity = propFluid->fDensity, fSheat = propFluid->fSpecificHeat;
  short i,type;
  DOUBLE t,sHeat,tmp1;
  
  sHeat   = MAT2D(0,SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
  for(i=0;i<MAXLOAD;i++)
  {
    type               = loadsVel[i].type;
/*...*/
    if ( type == INLET ||  type == OPEN) 
    {
      t = loadsTemp[i].par[0];
      if(fDensity)
      {
        tmp1 = airDensity(&propFluid->den
                   ,t               , thDynamic.pTh[2]
                   ,thDynamic.pTh[2], propFluid->molarMass
                   , iKelvin);
        loadsVel[i].par[loadsVel[i].np-1] = tmp1; 
        loadsVel[i].density               = tmp1;
      }
/*....................................................................*/
    }     
  }
/*....................................................................*/
}
/*********************************************************************/


/********************************************************************* 
 * Data de criacao    : 10/05/2019                                   *
 * Data de modificaco : 16/10/2019                                   *
 *-------------------------------------------------------------------*
 * setLoadsVel:                                                      *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
static void setLoadsRho(PropVarFluid *pf       ,Loads *loadsRho      
                          ,Loads *loadsVel     ,Loads *loadsTemp
                          ,bool const iKelvin )

{
  short i,type;
  DOUBLE tmp1,t;

  for(i=0;i<MAXLOAD;i++)
  {
    type               = loadsVel[i].type;
/*...*/
    if ( type == INLET ||  type == OPEN) 
    {
      t = loadsTemp[i].par[0];
      loadsRho[i].type  = type == INLET ? INLET:OPEN; 
      loadsRho[i].nTypeVar = loadsVel[i].nTypeVar;
      loadsRho[i].np       = loadsVel[i].np;
      loadsRho[i].fUse     = loadsVel[i].fUse;      
      if(pf->fDensity)
      {
        tmp1 = airDensity(&pf->den
                   ,t               , thDynamic.pTh[2]
                   ,thDynamic.pTh[2], pf->molarMass
                   ,iKelvin);
        loadsRho[i].par[0] = tmp1; 
        
      }
      else
        loadsRho[i].par[0] = loadsVel[i].density; 
/*....................................................................*/
    }     
  }
/*....................................................................*/
}
/*********************************************************************/


/********************************************************************* 
 * Data de criacao    : 27/09/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * readTransConfig: lendo a configura do termpo transisnte           *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void readTransConfig(Temporal *ddt, Save *save
                   , Macros *mm   , FILE *fileIn)
{
      char word[WORD_SIZE];
/*...*/
      ddt->flag = true;
      readMacro(fileIn,word,false);
      if(!strcmp(word,"config:")){
/*... timer*/        
        readMacro(fileIn,word,false);
        setTransientScheme(word,ddt);
/*...*/ 
        fscanf(fileIn,"%lf",&ddt->dt[0]);
        fscanf(fileIn,"%lf",&ddt->total);
/*...*/
        readMacro(fileIn,word,false);
        if(!strcmp(word,"dynamic"))     
          ddt->fDynamic = true;
        if(ddt->fDynamic){
          fprintf(fileLogExc,"dynamic : True\n");
        }
        else {
          fprintf(fileLogExc,"dynamic : False\n");
        }          
/*...*/        
        fprintf(fileLogExc,"dt(s)     : %.10lf\n",ddt->dt[0]);
        fprintf(fileLogExc,"Total(s)  : %.10lf\n",ddt->total);
      
        if(ddt->typeReal == EULER)     
          fprintf(fileLogExc,"ddtScheme : EULER\n");
        else if(ddt->typeReal == BACKWARD)     
          fprintf(fileLogExc,"ddtScheme : BACKWARD\n");

        if(!save->fLoad)
        {
          ddt->t         = 0.e0;
          ddt->dtInicial = ddt->dt[0];
          ddt->dt[1]     = ddt->dt[0];
          ddt->dt[2]     = ddt->dt[0];
          ddt->timeStep  = 0;
        }
      }
/*...................................................................*/

/*...*/
      mm->flWord = true;
      mm->kLoop  = 0;
      mm->jLoop  = 0;
      do{
        readMacro(fileIn,word,false);
        strcpy(mm->loopWord[mm->kLoop],word);
        mm->kLoop++;
        if(mm->kLoop > 100)
        {
          ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__,
                   "Numero de comandos na macro trasient execedido"
                   ,EXIT_PROG); 
        }
      }while(strcmp(word,"endTransient"));
      strcpy(mm->loopWord[mm->kLoop-1],"nextLoop");
/*...................................................................*/


}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 09/07/2016                                   *
 * Data de modificaco : 08/10/2019                                   * 
 *-------------------------------------------------------------------* 
 * SETSIMPLESCHEME : set o metodo simple incompressivel              *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void setSimpleScheme(char *fileName, Simple *sp)
{

  FILE *fileIn=NULL;
  short j=0;
  char str[] = {"end"};
  char word[WORD_SIZE];
  char macro[][WORD_SIZE] = { "name"       ,"alphapres"    /*0,1*/
                             ,"alphavel"   ,"north"        /*2,3*/   
                             ,"presidual"  ,"maxit"        /*4,5*/
                             ,"alphaenergy","alphaz"       /*6,7*/
                             ,"alpharho"   };              /*8*/ 
  
  fileIn = openFile(fileName,"r");

/*...*/
  readMacroV2(fileIn, word, false, true);
  do
  {
      j = 0;
/*... name*/
      if (!strcmp(word, macro[j++]))
      {
        readMacroV2(fileIn, word, false, true);
        if(!strcmp(word,"SIMPLE"))
        {
          sp->type = SIMPLE;
          fprintf(fileLogExc, "%-20s : %s\n", "PRES-VEL", "SIMPLE");
        }
        else if(!strcmp(word,"SIMPLEC"))
        {
          sp->type =  SIMPLEC;
          fprintf(fileLogExc,"%-20s : %s\n","PRES-VEL","SIMPLEC");
        }
      }
/*.....................................................................*/

/*... alphaPres*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn,"%lf",&sp->alphaPres); 
        fprintf(fileLogExc,"%-20s : %lf\n","alphaPres", sp->alphaPres);
      }
/*.....................................................................*/

/*... alphaVel*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn,"%lf" ,&sp->alphaVel);
        fprintf(fileLogExc, "%-20s : %lf\n","alphaVel", sp->alphaVel);
      }
/*.....................................................................*/

/*... nOrth*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn,"%u" ,&sp->nNonOrth);
        fprintf(fileLogExc, "%-20s : %u\n","nNonOrth", sp->nNonOrth);
      }
/*.....................................................................*/

/*... pResidual*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn,"%d" ,&sp->pSimple);
        fprintf(fileLogExc,"%-20s : %d\n","pSimple", sp->pSimple);
      }
/*.....................................................................*/

/*... maxit*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn,"%d" ,&sp->maxIt);
        fprintf(fileLogExc,"%-20s : %d\n","Maxit", sp->maxIt);
      }
/*.....................................................................*/

/*... alphaEnergy*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%lf", &sp->alphaEnergy);
        fprintf(fileLogExc, "%-20s : %lf\n", "alphaEnergy", sp->alphaEnergy);
      }
/*.....................................................................*/

/*... alphaZ*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%lf", &sp->alphaComb);
        fprintf(fileLogExc, "%-20s : %lf\n","alphaZ", sp->alphaComb);
      }
/*.....................................................................*/

/*... rho*/
      else if (!strcmp(word, macro[j++]))
      {
        fscanf(fileIn, "%lf", &sp->alphaDensity);
        fprintf(fileLogExc,"%-20s : %lf\n","alphaEnergy", sp->alphaEnergy);
      }
/*.....................................................................*/

      readMacroV2(fileIn, word, false, true);
    }while(strcmp(word,str));
/*.....................................................................*/
    fclose(fileIn);
}
/*********************************************************************/


/********************************************************************* 
 * Data de criacao    : 21/10/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * initFaceRrho : set o metodo simple incompressivel                 *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void initFaceRrho(Loads *ldVel
                , short  *RESTRICT faceRvel, short  *RESTRICT faceRrho
                , short  *RESTRICT nFace   , short const maxViz
                , INT numel)
{
  short j,ty;
  INT nel,aux1=maxViz+1;

  for(nel=0;nel<numel;nel++)
  {
    for(j=0;j<nFace[nel];j++)
    {
      ty = MAT2D(nel, j, faceRvel, aux1);
      if(ty>0)
      {
        if(ldVel[ty-1].type == INLET)
          MAT2D(nel, j, faceRrho, aux1) = ty;      
      }  
    }
  }
    

}
