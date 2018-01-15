#include<ReadFile.h>

/*...funcao de apoio*/
  static void getword(char *line, char*word);
  static int getnumprop2(char *line);
  static void convLoadsEnergy(Loads *loadsEnergy  ,Loads *loadsTemp
                             ,DOUBLE *RESTRICT prop
                             ,bool const fTemp    ,bool const fSheat                             
                             ,bool const fKelvin );
  static void convLoadsPresC(Loads *loadsPres,Loads *loadsPresC);
/*..................................................................*/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 05/09/2017                                   *
 *-------------------------------------------------------------------*
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
void readFileFvMesh( Memoria *m        , Mesh *mesh
                   , PropVar prop      , EnergyModel energyModel
                   , FILE* file)
{
  char word[WORD_SIZE],str[WORD_SIZE];
  char macro[NMACROS][WORD_SIZE]={
          "coordinates","endMesh"    ,"insert"         /* 0, 1, 2*/
         ,"return"     ,"cells"      ,"faceRt1"        /* 3, 4, 5*/
         ,"faceLoadT1" ,"loadsT1"    ,""               /* 6, 7, 8*/ 
         ,""           ,""           ,""               /* 9,10,11*/ 
         ,"faceRd1"    ,""           ,"loadsD1"        /*12,13,14*/ 
         ,"faceLoadD1" ,""           ,""               /*15,16,17*/ 
         ,"faceRvel"   ,"loadsVel"   ,"faceLoadVel"    /*18,19,20*/ 
         ,"faceRpres"  ,"loadsPres"  ,"faceLoadPres"   /*21,22,23*/
         ,"faceRtemp"  ,"loadsTemp"  ,"faceLoadTemp"   /*24,25,26*/
         ,"materials"  ,"uniformPres","initialVel"     /*27,28,29*/
         ,"uniformTemp","uniformVel" ,""               /*30,31,32*/
         ,""           ,""           ,""               /*33,34,35*/
         ,""           ,""           ,""               /*36,37,38*/
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
  ndfVel = max(mesh->ndfF - 1,mesh->ndfFt - 2);
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
  if(mesh->ndfF > 0 || mesh->ndfFt > 0) {
/*... alocando memoria*/
/*... cc da equacao de velocidades*/
     HccaAlloc(short,m,mesh->elm.faceRvel  
            ,nel*(maxViz+1),"faceRvel"    ,_AD_);
     zero(mesh->elm.faceRvel  ,nel*(maxViz+1),"short"  );
     
     HccaAlloc(short,m,mesh->elm.faceLoadVel  
            ,nel*(maxViz+1),"faceLVel"    ,_AD_);
     zero(mesh->elm.faceLoadVel  ,nel*(maxViz+1),"short"  );

/*... cc da equacao de pressao*/
     HccaAlloc(short,m,mesh->elm.faceRpres 
            ,nel*(maxViz+1),"faceRpres"   ,_AD_);
     zero(mesh->elm.faceRpres ,nel*(maxViz+1),"short"  );
     
     HccaAlloc(short,m,mesh->elm.faceLoadPres 
            ,nel*(maxViz+1),"faceLPres"   ,_AD_);
     zero(mesh->elm.faceLoadPres ,nel*(maxViz+1),"short"  );

/*... cc da equacao de energia*/
     if (mesh->ndfFt > 0 ){
        HccaAlloc(short, m, mesh->elm.faceRenergy
                  , nel*(maxViz + 1), "faceRenergy", _AD_);
        zero(mesh->elm.faceRenergy, nel*(maxViz + 1), "short");

        HccaAlloc(short, m, mesh->elm.faceLoadEnergy
                 , nel*(maxViz + 1), "faceLEnergy", _AD_);
        zero(mesh->elm.faceLoadEnergy, nel*(maxViz + 1), "short");

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

     }
/*... viscosidade turbulenta*/
     HccaAlloc(DOUBLE, m, mesh->elm.eddyViscosity
              , nel  , "eddyVis", _AD_);
     zero(mesh->elm.eddyViscosity, nel, DOUBLEC);

/*... viscosidade turbulenta*/
     HccaAlloc(DOUBLE, m, mesh->elm.stressR        
              , nel*mesh->ntn , "stressR", _AD_);
     zero(mesh->elm.stressR, nel*mesh->ntn, DOUBLEC);

/*... viscosidade turbulenta*/
     HccaAlloc(DOUBLE, m, mesh->elm.wallParameters
              , nel*4   , "wallParm", _AD_);
     zero(mesh->elm.wallParameters, nel*4, DOUBLEC);

/*... viscosidade turbulenta*/
     HccaAlloc(DOUBLE, m, mesh->elm.cd
              , nel*2   , "cDynamic"     , _AD_);
     zero(mesh->elm.cd, nel*2, DOUBLEC);

/*... densityFluid*/
     HccaAlloc(DOUBLE , m         , mesh->elm.densityFluid
              ,nel * DENSITY_LEVEL, "eDenFluid", _AD_);
     zero(mesh->elm.densityFluid, nel * DENSITY_LEVEL, DOUBLEC);


/*... ePres(n+1)*/
     HccaAlloc(DOUBLE,m,mesh->elm.pressure 
            ,nel              ,"pressure"          ,_AD_);
     zero(mesh->elm.pressure  ,nel                         ,DOUBLEC);

/*... ePres(n)*/
     HccaAlloc(DOUBLE,m,mesh->elm.pressure0
            ,nel              ,"pressure0"         ,_AD_);
     zero(mesh->elm.pressure0 ,nel                         ,DOUBLEC);

/*... pres*/
     HccaAlloc(DOUBLE,m,mesh->node.pressure
              ,nn    ,"npressure",_AD_);
     zero(mesh->node.pressure     ,nn                         ,DOUBLEC);

/*... nGradVel*/
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
     
/*... eGradPres*/
     HccaAlloc(DOUBLE,m,mesh->elm.gradPres
              ,nel*ndm        ,"eGradPres"     ,_AD_);
     zero(mesh->elm.gradPres  ,nel*ndm              ,DOUBLEC);

/*... */
     if (mesh->ndfFt > 0) {
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

     }

     if( mpiVar.nPrcs < 2){
/*... rCellVel*/
       HccaAlloc(DOUBLE,m,mesh->elm.rCellVel  
                ,nel*ndm*ndfVel       ,"rCellVel"     ,_AD_);
       zero(mesh->elm.rCellVel  ,nel*ndm*ndfVel   ,DOUBLEC);
/*... rCellPres*/
       HccaAlloc(DOUBLE,m,mesh->elm.rCellPres 
                ,nel                 ,"rCellPres"    ,_AD_);
       zero(mesh->elm.rCellPres ,nel              ,DOUBLEC);
/*... rCellEnergy*/
       if (mesh->ndfFt > 0){
          HccaAlloc(DOUBLE,m            ,mesh->elm.rCellEnergy
                   ,nel   ,"rCellEnergy",_AD_);
          zero(mesh->elm.rCellEnergy, nel, DOUBLEC);
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
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[0] = true;
      fprintf(fileLogExc,"loading coordinates...\n");
      readVfCoor(mesh->node.x,mesh->nnode,mesh->ndm,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
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
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[4] = true;
      fprintf(fileLogExc,"loading cells ...\n");
      readVfElmt(mesh->elm.node    ,mesh->elm.mat
                ,mesh->elm.nen     ,mesh->elm.adj.nViz
                ,mesh->elm.geomType,mesh->numel
                ,mesh->maxNo       ,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceRt1 */
    else if((!strcmp(word,macro[5])) && (!rflag[5])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[5] = true;
      strcpy(str,"endFaceRt1");
      fprintf(fileLogExc,"loading faceRt1 ...\n");
      readVfRes(mesh->elm.faceRt1,mesh->numel,mesh->maxViz+1,str,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceLoadT1 - cargas nas faces transporte */
    else if((!strcmp(word,macro[6])) && (!rflag[6])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[6] = true;
      strcpy(str,"endFaceLoadT1");
      fprintf(fileLogExc,"loading faceLoadT1 ...\n");
      readVfRes(mesh->elm.faceLoadT1,mesh->numel
               ,mesh->maxViz+1       ,str       ,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... loadT1 - definicao de cargar transporte */
    else if((!strcmp(word,macro[7])) && (!rflag[7])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[7] = true;
      strcpy(str,"endLoadsT1");
      fprintf(fileLogExc,"loading loadsT1 ...\n");
      readVfLoads(loadsT1,str       ,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceRd1- condicao de contorno para problemas de difusa pura */
    else if((!strcmp(word,macro[12])) && (!rflag[12])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[12] = true;
      strcpy(str,"endFaceRd1");
      fprintf(fileLogExc,"loading faceRd1 ...\n");
      readVfRes(mesh->elm.faceRd1,mesh->numel,mesh->maxViz+1,str,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... loadD1 - definicao de cargar difusao pura */
    else if((!strcmp(word,macro[14])) && (!rflag[14])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[14] = true;
      strcpy(str,"endLoadsD1");
      fprintf(fileLogExc,"loading loadsD1 ...\n");
      readVfLoads(loadsD1,str       ,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceLoadD1 - cargas nas faces difusao pura */
    else if((!strcmp(word,macro[15])) && (!rflag[15])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[15] = true;
      strcpy(str,"endFaceLoadD1");
      fprintf(fileLogExc,"loading faceLoadD1 ...\n");
      readVfRes(mesh->elm.faceLoadD1,mesh->numel
               ,mesh->maxViz+1       ,str       ,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceRvel- condicao de contorno para problemas fluidos (Vel) */
    else if((!strcmp(word,macro[18])) && (!rflag[18])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[18] = true;
      strcpy(str,"endFaceRvel");
      fprintf(fileLogExc,"loading faceRvel ...\n");
      readVfRes(mesh->elm.faceRvel,mesh->numel,mesh->maxViz+1,str,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... loadVel - definicao de cargar fluidos (Vel)*/
    else if((!strcmp(word,macro[19])) && (!rflag[19])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[19] = true;
      strcpy(str,"endLoadsVel");
      fprintf(fileLogExc,"loading loadsVel ...\n");
      readVfLoads(loadsVel,str       ,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceLoadVel - cargas nas faces fluido (Vel)*/
    else if((!strcmp(word,macro[20])) && (!rflag[20])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[20] = true;
      strcpy(str,"endFaceLoadVel");
      fprintf(fileLogExc,"loading faceLoadVel ...\n");
      readVfRes(mesh->elm.faceLoadVel  ,mesh->numel
               ,mesh->maxViz+1         ,str        ,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceRpres - condicao de contorno para problemas fluidos (Pres)*/
    else if((!strcmp(word,macro[21])) && (!rflag[21])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[21] = true;
      strcpy(str,"endFaceRpres");
      fprintf(fileLogExc,"loading faceRpres ...\n");
      readVfRes(mesh->elm.faceRpres,mesh->numel
               ,mesh->maxViz+1     ,str        ,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... loadPres - definicao de cargar fluidos (Pres)*/
    else if((!strcmp(word,macro[22])) && (!rflag[22])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[22] = true;
      strcpy(str,"endLoadsPres");
      fprintf(fileLogExc,"loading loadsPres ...\n");
      readVfLoads(loadsPres,str       ,file);
      fprintf(fileLogExc,"done.\n");
      convLoadsPresC(loadsPres, loadsPresC);
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceLoadPres - cargas nas faces fluido (Pres)*/
    else if((!strcmp(word,macro[23])) && (!rflag[23])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[23] = true;
      strcpy(str,"endFaceLoadPres");
      fprintf(fileLogExc,"loading faceLoadPres ...\n");
      readVfRes(mesh->elm.faceLoadPres ,mesh->numel
               ,mesh->maxViz+1         ,str        ,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... faceRenergy - condicao de contorno para problemas fluidos (Energy)*/
    else if ((!strcmp(word, macro[24])) && (!rflag[24])) {
      fprintf(fileLogExc,"%s\n", DIF);
      fprintf(fileLogExc,"%s\n", word);
      strcpy(macros[nmacro++], word);
      rflag[24] = true;
      strcpy(str, "endFaceRtemp");
      fprintf(fileLogExc,"loading faceRtemp ...\n");
      readVfRes(mesh->elm.faceRenergy,mesh->numel
               ,mesh->maxViz + 1     ,str        ,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n", DIF);
    }
/*...................................................................*/

/*... loadEnergy - definicao de cargar fluidos (Energy)*/
    else if ((!strcmp(word, macro[25])) && (!rflag[25])) {
      fprintf(fileLogExc,"%s\n", DIF);
      fprintf(fileLogExc,"%s\n", word);
      strcpy(macros[nmacro++], word);
      rflag[25] = true;
      strcpy(str, "endLoadsTemp");
      fprintf(fileLogExc,"loading loadsTemp ...\n");
      readVfLoads(loadsTemp, str, file);
      fprintf(fileLogExc,"done.\n");
      convLoadsEnergy(loadsEnergy             ,loadsTemp
                     ,mesh->elm.material.prop
                     ,energyModel.fTemperature,prop.fSpecificHeat
                     ,energyModel.fKelvin);  
      fprintf(fileLogExc,"%s\n\n", DIF);
    }
/*...................................................................*/

/*... faceLoadEnergy - cargas nas faces fluido (Energy)*/
    else if ((!strcmp(word, macro[26])) && (!rflag[26])) {
      fprintf(fileLogExc,"%s\n", DIF);
      fprintf(fileLogExc,"%s\n", word);
      strcpy(macros[nmacro++], word);
      rflag[26] = true;
      strcpy(str, "endFaceLoadTemp");
      fprintf(fileLogExc,"loading faceLoadTemp ...\n");
      readVfRes(mesh->elm.faceLoadEnergy,mesh->numel
               ,mesh->maxViz + 1        ,str        ,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n", DIF);
    }
/*...................................................................*/

/*... materiais */
    else if((!strcmp(word,macro[27])) && (!rflag[27])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[27] = true;
      strcpy(str,"endMaterials");
      fprintf(fileLogExc,"loading materials ...\n");
      readVfMat(mesh->elm.material.prop,mesh->elm.material.type
               ,numat,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... uniformPres */
    else if ((!strcmp(word, macro[28])) && (!rflag[28])) {
      fprintf(fileLogExc,"%s\n", DIF);
      fprintf(fileLogExc,"%s\n", word);
      strcpy(macros[nmacro++], word);
      rflag[28] = true;
      uniformField(mesh->elm.pressure, mesh->numel,1, file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n", DIF);
    }
/*...................................................................*/

/*... initialVel */
    else if((!strcmp(word,macro[29])) && (!rflag[29])){
      fprintf(fileLogExc,"%s\n",DIF);
      fprintf(fileLogExc,"%s\n",word);
      strcpy(macros[nmacro++],word);
      rflag[29] = true;
      strcpy(str,"endInitialVel");
      fprintf(fileLogExc,"loading initialVel ...\n");
      readVfInitial(mesh->elm.vel,mesh->numel,mesh->ndm,str,file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n",DIF);
    }
/*...................................................................*/

/*... uniformTemp */
    else if ((!strcmp(word, macro[30])) && (!rflag[30])) {
      fprintf(fileLogExc,"%s\n", DIF);
      fprintf(fileLogExc,"%s\n", word);
      strcpy(macros[nmacro++], word);
      rflag[30] = true;
      fprintf(fileLogExc,"loading uniformTemp ...\n");
      uniformField(mesh->elm.temp0, mesh->numel, 1, file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n", DIF);
    }
/*...................................................................*/

/*... uniformVel */
    else if ((!strcmp(word, macro[31])) && (!rflag[31])) {
      fprintf(fileLogExc,"%s\n", DIF);
      fprintf(fileLogExc,"%s\n", word);
      strcpy(macros[nmacro++], word);
      rflag[31] = true;
      fprintf(fileLogExc,"loading uniformVel ...\n");
      uniformField(mesh->elm.vel, mesh->numel, ndm, file);
      fprintf(fileLogExc,"done.\n");
      fprintf(fileLogExc,"%s\n\n", DIF);
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

/*...*/
  if (mesh->ndfFt > 0) {
/*...*/
    alphaProdVector(1.e0        ,mesh->elm.temp0
                     ,mesh->numel ,mesh->elm.temp);
/*...................................................................*/

    if(energyModel.fTemperature){
/*... convertendo temperatura para kelvin*/
      if(energyModel.fKelvin)
        convTempForKelvin(mesh->elm.temp0,mesh->numel,true); 

      alphaProdVector(1.e0        ,mesh->elm.temp
                     ,mesh->numel ,mesh->elm.energy0);

      alphaProdVector(1.e0        ,mesh->elm.energy0
                     ,mesh->numel ,mesh->elm.energy);
    }
    else{
      getEnergyForTemp(mesh->elm.temp         ,mesh->elm.energy0
                      ,mesh->elm.material.prop,mesh->elm.mat                        
                      ,mesh->numel            
                      ,prop.fSpecificHeat     ,energyModel.fKelvin
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
    if(prop.fDensity){
/*... inicia a massa especifica com o campo de temperatura inicial*/
      initPropTemp(mesh->elm.densityFluid ,mesh->elm.temp 
                  ,mesh->elm.material.prop,mesh->elm.mat
                  ,DENSITY_LEVEL          ,mesh->numel
                  ,energyModel.fKelvin    ,DENSITY);
    }
    else
      initProp(mesh->elm.densityFluid 
              ,mesh->elm.material.prop,mesh->elm.mat
              ,DENSITY_LEVEL          ,mesh->numel
              ,DENSITY);
/*...................................................................*/

/*... inicializando o calor especifico*/
    if(prop.fSpecificHeat)
      initPropTemp(mesh->elm.specificHeat   ,mesh->elm.temp 
                  ,mesh->elm.material.prop  ,mesh->elm.mat
                  ,SHEAT_LEVEL              ,mesh->numel
                  ,energyModel.fKelvin      ,SPECIFICHEATCAPACITYFLUID);
    else
      initProp(mesh->elm.specificHeat  
             ,mesh->elm.material.prop,mesh->elm.mat
             ,SHEAT_LEVEL              ,mesh->numel
             ,SPECIFICHEATCAPACITYFLUID);
/*...................................................................*/

/*... inicializando a viscosidade dinamica*/
    if(prop.fDynamicViscosity)
      initPropTemp(mesh->elm.dViscosity     ,mesh->elm.temp 
                ,mesh->elm.material.prop  ,mesh->elm.mat
                ,DVISCOSITY_LEVEL         ,mesh->numel
                ,energyModel.fKelvin      ,DYNAMICVISCOSITY);
   else
      initProp(mesh->elm.dViscosity 
              ,mesh->elm.material.prop  ,mesh->elm.mat
              ,DVISCOSITY_LEVEL         ,mesh->numel
              ,DYNAMICVISCOSITY);
/*...................................................................*/

/*... inicializando a condutividade termica*/
    if(prop.fThermalCondutivty)
      initPropTemp(mesh->elm.tConductivity ,mesh->elm.temp 
                ,mesh->elm.material.prop   ,mesh->elm.mat
                ,TCONDUCTIVITY_LEVEL       ,mesh->numel
                ,energyModel.fKelvin       ,THERMALCONDUCTIVITY);
   else
      initProp(mesh->elm.tConductivity 
              ,mesh->elm.material.prop  ,mesh->elm.mat
              ,TCONDUCTIVITY_LEVEL         ,mesh->numel
              ,THERMALCONDUCTIVITY);
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
 
  for(i=0;i<nn;i++){
    fscanf(file,"%ld",&idum);
    k = (INT) idum -1;
    for(j=0;j<ndm;j++){
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
          fprintf(fileLogExc,"%-20s: %s\n","bVtk","true");
        else
          fprintf(fileLogExc,"%-20s: %s\n","bVtk","false");
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
          fprintf(fileLogExc,"%-20s: %s\n","Reord","true");
        else
          fprintf(fileLogExc,"%-20s: %s\n","Reod","false");
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
        fprintf(fileLogExc,"%-20s: %d\n","Memory(MBytes)"
              ,(int)(nmax/conv));
    }
/*... rcGrad*/   
    else if(!strcmp(word,config[3])){
      fscanf(file,"%s",s);
      if(!strcmp(s,"gglc")){
        *rcGrad = RCGRADGAUSSC;
        if(!mpiVar.myId)
          fprintf(fileLogExc,"%-20s: %s\n","rcGrad","GreenGaussCell");
      }
      else if(!strcmp(s,"ggln")){
        *rcGrad = RCGRADGAUSSN;
        if(!mpiVar.myId)
          fprintf(fileLogExc,"%-20s: %s\n","rcGrad","GreenGaussNode");
      }
      else if(!strcmp(s,"lSquare")){
        *rcGrad = RCLSQUARE;
        if(!mpiVar.myId)
          fprintf(fileLogExc,"%-20s: %s\n","rcGrad","LeastSquare");
      }
      
      else if(!strcmp(s,"lSquareQR")){
        *rcGrad = RCLSQUAREQR;
        if(!mpiVar.myId)
          fprintf(fileLogExc,"%-20s: %s\n","rcGrad","LeastSquareQR");
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
          fprintf(fileLogExc,"%-20s: %s\n","fItPlotRes","true");
      }
      else{
        opt->fItPlotRes = false;
        if(!mpiVar.myId)
          fprintf(fileLogExc,"%-20s: %s\n","fItPlotRes","false");
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
          fprintf(fileLogExc,"%-20s: %s\n","fItPlot","true");
      }
      else{
        opt->fItPlot = false;
        if(!mpiVar.myId)
          fprintf(fileLogExc,"%-20s: %s\n","fItPlot","false");
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
	mesh->ndfF  = 0;
	mesh->ndfFt = 0;
/*...................................................................*/

  readMacro(file,word,false);
  do{
/*... transport1*/
    if(!strcmp(word,"transport1")){
      readMacro(file,word,false);
      mesh->ndfT[0]=(short) atol(word);
      if(!mpiVar.myId ) fprintf(fileLogExc,"transport1 ndf %d\n"
                              ,mesh->ndfT[0]);
    }
/*...................................................................*/

/*... transport2*/
    else if(!strcmp(word,"transport2")){
      readMacro(file,word,false);
      mesh->ndfT[1]= (short) atol(word);
      if(!mpiVar.myId ) fprintf(fileLogExc,"transport2 ndf %d\n"
                               ,mesh->ndfT[1]);
    }
/*...................................................................*/

/*... transport3*/
    else if(!strcmp(word,"transport3")){
      readMacro(file,word,false);
      mesh->ndfT[2]= (short) atol(word);
      if(!mpiVar.myId ) fprintf(fileLogExc,"transport3 ndf %d\n"
                               ,mesh->ndfT[2]);
    }
/*...................................................................*/

/*... diffusion1*/
    else if(!strcmp(word,"diffusion1")){
      readMacro(file,word,false);
      mesh->ndfD[0]= (short) atol(word);
      if(!mpiVar.myId ) fprintf(fileLogExc,"diffusion1 ndf %d\n"
                              ,mesh->ndfD[0]);
    }
/*...................................................................*/

/*... diffusion2*/
    else if(!strcmp(word,"diffusion2")){ 
     readMacro(file,word,false);
      mesh->ndfD[1]= (short) atol(word);
      if(!mpiVar.myId ) fprintf(fileLogExc,"diffusion2 ndf %d\n"
                              ,mesh->ndfD[1]);
    }
/*...................................................................*/

/*... diffusion3*/
    else if(!strcmp(word,"diffusion3")){
      readMacro(file,word,false);
      mesh->ndfD[2]= (short) atol(word);
      if(!mpiVar.myId ) fprintf(fileLogExc,"diffusion3 ndf %d\n"
                              ,mesh->ndfD[2]);
    }
/*...................................................................*/

/*... fluido*/
    else if(!strcmp(word,"fluid")){
      readMacro(file,word,false);
      mesh->ndfF= (short) atol(word);
      if(!mpiVar.myId ) fprintf(fileLogExc,"fluid ndf %d\n"
                              ,mesh->ndfF);
    }
/*...................................................................*/

/*... fluido-termo ativado*/
    else if(!strcmp(word,"fluidt")){
      readMacro(file,word,false);
      mesh->ndfFt= (short) atol(word);
      if(!mpiVar.myId ) fprintf(fileLogExc,"fluidt ndf %d\n"
                              ,mesh->ndfFt);
    }
/*...................................................................*/
    
    readMacro(file,word,false);
  }while(strcmp(word,str));


}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 29/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
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
void readPropVar(PropVar *p,FILE *file){

  char *str={"endPropVar"};
  char word[WORD_SIZE];

/*...*/
  p->fDensity           = false;
  p->fSpecificHeat      = false;
  p->fDynamicViscosity  = false;
  p->fThermalCondutivty = false;
/*...................................................................*/

  readMacro(file,word,false);
  do{
/*... specific heat*/
    if(!strcmp(word,"sHeat")){
      readMacro(file,word,false);
      p->fSpecificHeat = true;
      initSheatPol(); 
      if(!mpiVar.myId && p->fSpecificHeat) 
        fprintf(fileLogExc,"sHeat variation        : Enable\n");
    }
/*...................................................................*/

/*... densidade*/
    else if(!strcmp(word,"density")){
      readMacro(file,word,false);
      p->fDensity = true;
      convStringLower(word);
      initDensityPol(word);
      if(!mpiVar.myId && p->fDensity) 
        fprintf(fileLogExc,"Density variation      : Enable\n");                           
    }
/*...................................................................*/

/*... specific heat*/
    else if(!strcmp(word,"dViscosity")){
      readMacro(file,word,false);
      p->fDynamicViscosity = true;
      initDviscosityPol(word); 
      if(!mpiVar.myId && p->fDynamicViscosity)
        fprintf(fileLogExc,"dViscosity variation   : Enable\n");;
    }
/*...................................................................*/

/*... specific heat*/
    else if(!strcmp(word,"tCondutivity")){
      readMacro(file,word,false);
      p->fThermalCondutivty = true;
      initThCondPol(word);
      if(!mpiVar.myId && p->fThermalCondutivty)
        fprintf(fileLogExc,"tCondutivity variation : Enable\n");                          
    }
/*...................................................................*/
    
    readMacro(file,word,false);
  }while(strcmp(word,str));


}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 04/09/2017                                   *
 * Data de modificaco : 13/01/2018                                   *
 *-------------------------------------------------------------------* 
 * READPROPVAR : propriedades variaveis                              * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * e       -> modelos/termos usandos na eq energia                   * 
 * t       -> modelo de turbilencia                                  *
 * eMass   -> modelos/termos usados na equacao da conv de mass       * 
 * eMomentum  -> modelos/termos usados na equacao da conv de mass    *
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
void readModel(EnergyModel *e     , Turbulence *t
             , MassEqModel *eMass, MomentumModel *eMomentum
             , FILE *file){

  char *str={"endModel"};
  char format[1024];
  char word[WORD_SIZE];
  char macros[][WORD_SIZE] = {"energy"  ,"turbulence","mass"
                             ,"momentum"};

  char energy[][WORD_SIZE] = { "preswork", "dissipation", "residual"  
                             , "absolute", "temperature", "entalphy"}; 

  char turb[][WORD_SIZE] = { "smagorinsky","wallmodel" , "wale"     
                            ,"vreman"     ,"ldynamic"  , "sigmamodel"
                            ,"mixed"      ,"bardina"   , "clark"
                            ,"bardinaMod" ,"towdynamic", "gdynamic"
                            ,"gdynamicmod"}; 

  char mass[][WORD_SIZE] = { "lhsdensity","rhsdensity"}; 

  char momentum[][WORD_SIZE] = {"residual","absolute"
                               ,"presa"   ,"presref" 
                               ,"rhiechow","viscosity"
                               ,"div"};

  char typeWallModel[][WORD_SIZE] ={"standard","enhanced"};
  


  int i,nPar;

  readMacro(file,word,false);
  do{
/*... equacao da energia*/
    if(!strcmp(word,macros[0])){
/*...*/
      e->fPresWork    = false;
      e->fDissipation = false;
      e->fRes         = false;
      e->fTemperature = false;
      if(!mpiVar.myId)
        fprintf(fileLogExc,"\n%-20s: \n","EnergyModel");  
/*...................................................................*/      
      fscanf(file,"%d",&nPar);
      for(i=0;i<nPar;i++){
        strcpy(format,"%-20s: %s\n");
        readMacro(file,word,false);
        convStringLower(word);
/*... PresWork*/
        if(!strcmp(word,energy[0])){    
          e->fPresWork = true;
          if(!mpiVar.myId && e->fPresWork) 
            fprintf(fileLogExc,format,"PresWork","Enable");
        }
/*...................................................................*/

/*... Dissipation*/
        else if(!strcmp(word,energy[1])){        
          e->fDissipation = true;
          if(!mpiVar.myId && e->fDissipation) 
            fprintf(fileLogExc,format,"Dissipation","Enable");                           
        }
/*...................................................................*/

/*... Residual*/
        else if(!strcmp(word,energy[2])){
          e->fRes = true;
          if(!mpiVar.myId && e->fRes)
            fprintf(fileLogExc,format,"Residual","Enable");
        }
/*...................................................................*/

/*... Absolute*/
        else if(!strcmp(word,energy[3])){
          e->fRes = false;
          if(!mpiVar.myId && e->fRes)
            fprintf(fileLogExc,format,"Absolute","Enable");
        }
/*...................................................................*/

/*... Temperatura*/
        else if(!strcmp(word,energy[4])){
          e->fTemperature = true;
          if(!mpiVar.myId && e->fTemperature)
            fprintf(fileLogExc,format,"Temperatura","Enable");
        }
/*...................................................................*/

/*... Entalphy*/
        else if(!strcmp(word,energy[5])){
          e->fTemperature = false;
          if(!mpiVar.myId && e->fTemperature)
            fprintf(fileLogExc,format,"Entalphy","Enable");
        }
/*...................................................................*/

      }
/*...................................................................*/
    }
/*...................................................................*/

/*... turbulencia*/
    else if(!strcmp(word,macros[1])){   
      if(!mpiVar.myId)
        fprintf(fileLogExc,"\n%-20s: \n","TurbulenceModel");   
      fscanf(file,"%d",&nPar);
      for(i=0;i<nPar;i++){
        readMacro(file,word,false);
        convStringLower(word);
/*... Smagorinsky*/
        if(!strcmp(word,turb[0])){
          t->dynamic             = false;
          t->fTurb               = true;      
          t->type                = LES;
          t->typeMixed[FUNMODEL] = SMAGORINSKY;
          t->typeLes             = LESFUNCMODEL; 
          fscanf(file,"%lf",&t->cf);    
          if(!mpiVar.myId){ 
            fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[0],t->cf); 
          }
        }
/*...................................................................*/    

/*... wallModel*/
        else if(!strcmp(word,turb[1])){
          t->fWall = true;         
          readMacro(file,word,false); 
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
          if(!mpiVar.myId)
            fprintf(fileLogExc,"%-20s: %s\n",turb[1]
                              ,typeWallModel[t->wallType-1]);
        }
/*...................................................................*/

/*... Wale*/
        else if(!strcmp(word,turb[2])){
          t->dynamic                = false;
          t->fTurb                  = true;     
          t->type                   = LES;
          t->typeMixed[FUNMODEL]    = WALEMODEL; 
          t->typeLes = LESFUNCMODEL;
          fscanf(file,"%lf",&t->cf);    
          if(!mpiVar.myId){ 
            fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[2],t->cf); 
          }
        }
/*...................................................................*/ 

/*... Vreman*/
        else if(!strcmp(word,turb[3])){
          t->dynamic                = false;
          t->fTurb                  = true;      
          t->type                   = LES;
          t->typeMixed[FUNMODEL]    = VREMAN;
          t->typeLes = LESFUNCMODEL;
          fscanf(file,"%lf",&t->cf);    
          if(!mpiVar.myId){ 
             fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[3],t->cf); 
          }
        }
/*...................................................................*/

/*... Dynamic*/
        else if(!strcmp(word,turb[4])){
          t->dynamic     = true;
          t->fTurb       = true;      
          t->type        = LES;
          t->typeDynamic = LDYNAMIC;
          fscanf(file,"%lf",&t->c);  
          if(!mpiVar.myId){ 
            fprintf(fileLogExc,"%-20s:\n", "lDynamic"); 
          }
          setDynamicModelLes(t,file);
        }
/*...................................................................*/ 

/*... sigmaModel*/
        else if(!strcmp(word,turb[5])){
          t->dynamic             = false;
          t->fTurb               = true;      
          t->type                = LES;
          t->typeMixed[FUNMODEL] = SIGMAMODEL;
          t->typeLes = LESFUNCMODEL;
          fscanf(file,"%lf",&t->cf);   
          if(!mpiVar.myId){ 
             fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[5],t->cf); 
          }
        }
/*...................................................................*/ 

/*... mixed*/
        else if(!strcmp(word,turb[6])){
          t->dynamic = false;
          t->fTurb   = true;      
          t->type    = LES;
          t->typeLes = LESMIXEDMODEL; 
          fscanf(file,"%lf",&t->c);  
          if(!mpiVar.myId){ 
            fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[6],t->c); 
          }
          setMixedModelLes(t,file);
        }
/*...................................................................*/

/*... bardina*/
        else if(!strcmp(word,turb[7])){
          t->dynamic             = false;
          t->fTurb               = true;   
          t->type                = LES;
          t->typeMixed[ESTMODEL] = BARDINA;
          t->typeLes = LESSTRUMODEL;
          fscanf(file,"%lf",&t->cs);  
          if(!mpiVar.myId){ 
             fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[7],t->cs); 
          }
        }
/*...................................................................*/  

/*... clark*/
        else if(!strcmp(word,turb[8])){
          t->dynamic             = false;
          t->fTurb               = true;      
          t->type                = LES;
          t->typeMixed[ESTMODEL] = CLARK;
          t->typeLes             = LESSTRUMODEL;
          fscanf(file,"%lf",&t->cs);  
          if(!mpiVar.myId){ 
             fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[8],t->cs); 
          }
        }
/*...................................................................*/ 

/*... bardinaMod*/
        else if(!strcmp(word,turb[9])){
          t->dynamic             = false;
          t->fTurb               = true;      
          t->type                = LES;
          t->typeMixed[ESTMODEL] = BARDINAMOD;
          t->typeLes             = LESSTRUMODEL;
          fscanf(file,"%lf",&t->cs);  
          if(!mpiVar.myId){ 
             fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[9],t->cs); 
          }
        }
/*...................................................................*/   

/*... mixed 2 paramentros*/
        else if(!strcmp(word,turb[10])){
          t->dynamic     = true;
          t->fTurb       = true;      
          t->type        = LES;
          t->typeLes     = LESMIXEDTWOMODEL; 
          t->typeDynamic = TWOPARDYNAMIC;
          fscanf(file,"%lf",&t->c);  
          if(!mpiVar.myId){ 
            fprintf(fileLogExc,"%-20s:\n", "TowDynamic"); 
          }
          setMixedModelLes(t,file);
        }
/*...................................................................*/

/*... Gdynamic*/
        else if(!strcmp(word,turb[11])){
          t->dynamic     = true;
          t->fTurb       = true;      
          t->type        = LES;
          t->typeDynamic = GDYNAMIC;
          fscanf(file,"%lf",&t->c);  
          if(!mpiVar.myId){ 
            fprintf(fileLogExc,"%-20s:\n","gDynamic"); 
          }
          setDynamicModelLes(t,file);
        }
/*...................................................................*/ 

/*... Gdynamic*/
        else if(!strcmp(word,turb[12])){
          t->dynamic     = true;
          t->fTurb       = true;      
          t->type        = LES;
          t->typeDynamic = GDYNAMICMOD;
          fscanf(file,"%lf",&t->c);  
          if(!mpiVar.myId){ 
            fprintf(fileLogExc,"%-20s:\n","gDynamicMod"); 
          }
          setDynamicModelLes(t,file);
        }
/*...................................................................*/ 

      }
/*...................................................................*/
    }
/*...................................................................*/

/*... mass*/
    else if(!strcmp(word,macros[2])){ 
      strcpy(format,"%-20s: %s\n");
      if(!mpiVar.myId)
        fprintf(fileLogExc,"\n%-20s: \n","MassEqModel");    
      eMass->LhsDensity = false;
      eMass->RhsDensity = false;
      fscanf(file,"%d",&nPar);
      for(i=0;i<nPar;i++){
        readMacro(file,word,false);
        convStringLower(word);
/*... LhsDensity*/
        if(!strcmp(word,mass[0])){
          eMass->LhsDensity = true;          
          if(!mpiVar.myId && eMass->LhsDensity){ 
            fprintf(fileLogExc,format,"LhsDensity","Enable");
          }
        }
/*...................................................................*/

/*... RhsDensity*/
        else if(!strcmp(word,mass[1])){
          eMass->RhsDensity = true;          
          if(!mpiVar.myId && eMass->RhsDensity){ 
            fprintf(fileLogExc,format,"RhsDensity","Enable");
          }
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... Momentum*/
    else if(!strcmp(word,macros[3])){ 
      strcpy(format,"%-20s: %s\n");
      if(!mpiVar.myId)
        fprintf(fileLogExc,"\n%-20s: \n","MomentumEqModel");    
      eMomentum->fRes             = false;
      eMomentum->fAbsultePressure = false;
      eMomentum->fRhieChowInt     = false;
      eMomentum->fViscosity       = false;
      eMomentum->fDiv             = false;
      fscanf(file,"%d",&nPar);
      for(i=0;i<nPar;i++){
        readMacro(file,word,false);
        convStringLower(word);
/*... residual*/
        if(!strcmp(word,momentum[0])){
          eMomentum->fRes = true;          
          if(!mpiVar.myId && eMomentum->fRes){ 
            fprintf(fileLogExc,format,"Residual","Enable");
          }
        }
/*...................................................................*/

/*... Absolute*/
        else if(!strcmp(word,momentum[1])){
          eMomentum->fRes = false;            
          if(!mpiVar.myId && !eMomentum->fRes){ 
            fprintf(fileLogExc,format,"Absolute","Enable");
          }
        }
/*...................................................................*/

/*... presAbs*/
        else if(!strcmp(word,momentum[2])){
          eMomentum->fAbsultePressure = true;            
          if(!mpiVar.myId && eMomentum->fAbsultePressure ){ 
            fprintf(fileLogExc,format,"presAbs","Enable");
          }
        }
/*...................................................................*/

/*... presRef*/
        else if(!strcmp(word,momentum[3])){
          eMomentum->fAbsultePressure = false;            
          if(!mpiVar.myId && !eMomentum->fAbsultePressure ){ 
            fprintf(fileLogExc,format,"presRef","Enable");
          }
        }
/*...................................................................*/

/*... RhieChow*/
        else if(!strcmp(word,momentum[4])){
          eMomentum->fRhieChowInt = true;            
          if(!mpiVar.myId && eMomentum->fRhieChowInt){ 
            fprintf(fileLogExc,format,"RhieChowInt","Enable");
          }
        }
/*...................................................................*/

/*... Viscosity*/
        else if(!strcmp(word,momentum[5])){
          eMomentum->fViscosity = true;            
          if(!mpiVar.myId && eMomentum->fViscosity){ 
            fprintf(fileLogExc,format,"Viscosity","Enable");
          }
        }
/*...................................................................*/

/*... Viscosity*/
        else if(!strcmp(word,momentum[6])){
          eMomentum->fDiv = true;            
          if(!mpiVar.myId && eMomentum->fDiv){ 
            fprintf(fileLogExc,format,"Divergente","Enable");
          }
        }
/*...................................................................*/

      }
/*...................................................................*/
    }
/*...................................................................*/
    readMacro(file,word,false);
  }while(strcmp(word,str));

}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 30/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
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
void readGravity(DOUBLE *gravity,FILE *file){

  short i,n;
  char word[WORD_SIZE];

  readMacro(file,word,false);
  n = (short)atol(word);
  for(i=0;i<n;i++)
    fscanf(file,"%lf",gravity+i);

  if(!mpiVar.myId ) 
    fprintf(fileLogExc,"g = (%lf,%lf,%lf)\n"
           ,gravity[0],gravity[1],gravity[2]);

}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 17/07/2016                                   *
 * Data de modificaco : 13/12/2017                                   * 
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
  char format[1024];
  char word[WORD_SIZE];
  char macro[][WORD_SIZE] = 
               {"cell"         ,"node"        ,"vel"             /* 0, 1, 2*/
               ,"pres"         ,"gradvel"     ,"gradpres"        /* 3, 4, 5*/
               ,"temp"         ,"gradtemp"    ,"eddyviscosity"   /* 6, 7, 8*/
               ,"densityfluid" ,"specificheat","dviscosity"      /* 9,10,11*/
               ,"tconductivity","vorticity"   ,"wallparameters"  /*12,13,14*/
               ,"stress"       ,"kinecit"     ,"stressr"         /*15,16,17*/
               ,"cdynamic"     ,"qcriterion"  ,"prestotal"};     /*18,19,20*/  
  int tmp;

  strcpy(format,"%-20s: %s\n");

  opt->fCell          = false;
  opt->fNode          = false;
  opt->vel            = false;
  opt->pres           = false;
  opt->presTotal      = false;
  opt->energy         = false;
  opt->gradVel        = false;
  opt->gradPres       = false;
  opt->gradEnergy     = false;
  opt->eddyViscosity  = false;
  opt->densityFluid   = false;
  opt->specificHeat   = false;
  opt->dViscosity     = false;
  opt->tConductivity  = false;
  opt->vorticity      = false;
  opt->wallParameters = false;
  opt->stress         = false;
  opt->kinetic        = false;
  opt->stressR        = false;
  opt->cDynamic        = false;  
  opt->Qcriterion     = false;

  fscanf(file,"%d",&tmp);
  opt->stepPlotFluid[0] = opt->stepPlotFluid[1] = (short) tmp;
  readMacro(file,word,false);
  convStringLower(word);
  while(strcmp(word,str)){
/*... cell*/        
    if(!strcmp(word,macro[0])){ 
      opt->fCell = true;
      if(!mpiVar.myId ) fprintf(fileLogExc,format,"print","cell");
    }
/*.....................................................................*/

/*... node*/        
    else if(!strcmp(word,macro[1])){ 
      opt->fNode = true;
      if(!mpiVar.myId ) fprintf(fileLogExc,format,"print","node");
    }
/*.....................................................................*/

/*... vel*/        
    else if(!strcmp(word,macro[2])){ 
      opt->vel = true;
      if(!mpiVar.myId ) fprintf(fileLogExc,format,"print","vel");
    }
/*.....................................................................*/

/*...*/
    else if(!strcmp(word,macro[3])){ 
      opt->pres = true;
      if(!mpiVar.myId ) fprintf(fileLogExc,format,"print","pres");
    }
/*.....................................................................*/

/*...*/
    else if(!strcmp(word,macro[4])){ 
      opt->gradVel = true;
      if(!mpiVar.myId ) fprintf(fileLogExc,format,"print","gradVel");
    }
/*.....................................................................*/

/*...*/
    else if(!strcmp(word,macro[5])){ 
      opt->gradPres = true;
      if(!mpiVar.myId ) fprintf(fileLogExc,format,"print","gradPres");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[6])) {
      opt->energy = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","temp");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[7])) {
      opt->gradEnergy = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","gradTemp");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[8])) {
      opt->eddyViscosity = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","eddyViscosity");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[9])) {
      opt->densityFluid = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","densityFluid");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[10])) {
      opt->specificHeat = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","specificHeat");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[11])) {
      opt->dViscosity = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","dViscosity");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[12])) {
      opt->tConductivity = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","tConductivity");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[13])) {
      opt->vorticity = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","vorticity");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[14])) {
      opt->wallParameters = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","wallParameters");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[15])) {
      opt->stress = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","stress");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[16])) {
      opt->kinetic = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","kinecit");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[17])) {
      opt->stressR = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","stressR");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[18])) {
      opt->cDynamic = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","cDynamic");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[19])) {
      opt->Qcriterion = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","Qcriterion");
    }
/*.....................................................................*/

/*...*/
    else if (!strcmp(word,macro[20])) {
      opt->presTotal = true;
      if (!mpiVar.myId) fprintf(fileLogExc,format,"print","presTotal");
    }
/*.....................................................................*/

    readMacro(file,word,false);
    convStringLower(word);
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

  DOUBLE value[4];
  INT i;
  short j;
  
  for (j=0;j<ndf;j++)
    fscanf(file, "%lf", value+j);

  for (i = 0; i<n; i++) 
    for (j = 0; j<ndf; j++)
      MAT2D(i,j,field,ndf) = value[j];   
  
}
/*********************************************************************/

/*********************************************************************/
/*Converte condicoes de contorno da pressao para presssao de correcao*/
/*********************************************************************/
static void convLoadsPresC(Loads *loadsPres,Loads *loadsPresC){

  short i,j,type;

  for(i=0;i<MAXLOADFLUID;i++){
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
  }
}

/********************************************************************* 
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 05/09/2017                                   *
 *-------------------------------------------------------------------*
 * CONVLOADSENERGY: Converte condicoes de contorno da Temperatura    *
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
static void convLoadsEnergy(Loads *loadsEnergy   ,Loads *loadsTemp
                           ,DOUBLE *RESTRICT prop
                           ,bool const fTemp     ,bool const fSheat 
                           ,bool const iKelvin){

  short i,j,type;
  DOUBLE t,sHeat,tmp;
  
/*... cc da equacao da energia e em temperatura*/
  if(fTemp){

    for(i=0;i<MAXLOADFLUID;i++){
      loadsEnergy[i].type = loadsTemp[i].type;
      loadsEnergy[i].np   = loadsTemp[i].np;
      for(j=0;j<MAXLOADPARAMETER;j++)
        loadsEnergy[i].par[j] = loadsTemp[i].par[j];
    }
/*... converte c para kelvin*/    
    if(iKelvin)
      for(i=0;i<MAXLOADFLUID;i++){
        type = loadsEnergy[i].type;
        if( type == DIRICHLETBC ||  type == INLET 
        ||  type == CONVECTIONHEAT)
          loadsEnergy[i].par[0] 
                         = CELSIUS_FOR_KELVIN(loadsEnergy[i].par[0]);
      }
      
  }
/*....................................................................*/

  else{
    sHeat = MAT2D(0,SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
    for(i=0;i<MAXLOADFLUID;i++){
      loadsEnergy[i].type = loadsTemp[i].type;
      loadsEnergy[i].np   = loadsTemp[i].np;
      type = loadsTemp[i].type;
      for(j=0;j<MAXLOADPARAMETER;j++)
        loadsEnergy[i].par[j] = loadsTemp[i].par[j];

/*...*/
      if( type == DIRICHLETBC ){
        t = loadsTemp[i].par[0];
        tmp = tempForSpecificEnthalpy(t, sHeat, fSheat, iKelvin);
        loadsEnergy[i].par[0] = tmp;               
      }
/*....................................................................*/

/*...*/
      else if ( type == INLET ||  type == OPEN) {
        t = loadsTemp[i].par[1];
        tmp = tempForSpecificEnthalpy(t, sHeat, fSheat, iKelvin);
        loadsEnergy[i].par[1] = tmp;  
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
 * Data de criacao    : 05/11/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * convStringLower: Converte um string toda para caracteres minusculo*
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * s - recebe uma string                                             *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * s - converte os carateres maiuscolos para minuscolos              *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void convStringLower(char *s){
  char *p = s;
  for ( ; *p; ++p) *p = tolower(*p);
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 11/11/2017                                   *
 * Data de modificaco : 08/01/2018                                   *
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
  short iHelp = 5;
  char help [][WORD_SIZE] = 
               {"macros","setprintfluid","advection"    /* 0, 1, 2*/
               ,"model" ,"diffusion"    ,""             /* 3, 4, 5*/
               ,""      ,""             ,""             /* 6, 7, 8*/
               ,""      ,""             ,""             /* 9,10,11*/
               ,""      ,""             ,""};           /*12,13,15*/

  short iMacros = 18;
  char macros[][WORD_SIZE] = 
               {"setPrintFluid","setSolv"     ,"setSimple"    /* 0, 1, 2*/
               ,"pgeo"         ,"transient"   ,"simple"       /* 3, 4, 5*/
               ,"pFluid"       ,"endTransient","stop"         /* 6, 7, 8*/
               ,"model"        ,"propVar"     ,"gravity"      /* 9,10,11*/
               ,"edp"          ,"vorticity"   ,"openmp"       /*12,13,14*/
               ,"advection"    ,"diffusion"   ,"config" };    /*15,16,17*/

  short iPrint = 21;
  char print[][WORD_SIZE] = 
               {"cell"         ,"node"        ,"vel"             /* 0, 1, 2*/
               ,"pres"         ,"gradvel"     ,"gradpres"        /* 3, 4, 5*/
               ,"temp"         ,"gradtemp"    ,"eddyviscosity"   /* 6, 7, 8*/
               ,"densityfluid" ,"specificheat","dviscosity"      /* 9,10,11*/
               ,"tconductivity","vorticity"   ,"wallparameters"  /*12,13,14*/
               ,"stress"       ,"kinecit"     ,"stressr"         /*15,16,17*/
               ,"cdynamic"     ,"qcriterion"  ,"prestotal"};     /*18,19,20*/  
  
/*... adveccao*/
  char fAdv[][WORD_SIZE] =                                   
                         { "FoUp","CD" ,"SoUp"                /* 0, 1, 2*/
                          ,"TVD" ,"NVD","LUST"};              /* 3, 4, 5*/  
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
  short iModels = 11;
  char models[][WORD_SIZE] = {"energy"  ,"turbulence","mass"     /* 0, 1, 2*/
                             ,"momentum"};                       /* 3*/
  short iEnergy = 6;
  char energy[][WORD_SIZE] = { "preswork", "dissipation", "residual"  
                             , "absolute", "temperature", "entalphy"}; 
  short iTurb = 13;
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
                                 ,"towdynamic clark 1.0 smagorinsky 0.2"};   
  short iMass = 2;
  char mass[][WORD_SIZE] = { "lhsDensity","rhsDensity"}; 

  short iMom = 7;  
  char momentum[][WORD_SIZE] = {"residual","absolute"    /* 0, 1*/
                               ,"presa"   ,"presref"     /* 2, 3*/
                               ,"rhiechow","viscosity"   /* 4, 5*/
                               ,"div"};                  /* 6*/

  short iWall = 2;
  char typeWallModel[][WORD_SIZE] ={"standard","enhanced"};
/*....................................................................*/

  int i;                                                     

  printf("Help:\n");
  readMacro(f,word,false);
  convStringLower(word);
/*... macros*/
  if(!strcmp(word,help[0])){
    printf("Macros:\n");
    for(i=0;i<iMacros;i++)
      printf("%3d - %s\n",i+1,macros[i]);
    exit(EXIT_FAILURE);
  }
/*.....................................................................*/

/*... setPrintFluid*/        
  else if(!strcmp(word,help[1])){     
    printf("setPrintFluid 10 options end\n");
    printf("options:\n");
    for(i=0;i<iPrint;i++)
      printf("%3d - %s\n",i+1,print[i]);
    exit(EXIT_FAILURE);
  }
/*.....................................................................*/

/*... advection*/        
  else if(!strcmp(word,help[2])){   
    printf("Ex1:\n"); 
    printf("advection 3 Vel FoUp Energy CD Temp SoUp \n"); 
    printf("Ex2:\n"); 
    printf("advection 2 Vel NVD [options] Energy TVD [opitions]\n");
    printf("TVD options:\n");
    for(i=0;i<NFUNCLIMTFACE;i++)
      printf("%3d - %s\n",i+1,tvd[i]);
    printf("NVD options:\n");
    for(i=0;i<NFUNCNVD;i++)
      printf("%3d - %s\n",i+1,scheme_nvd[i]);
    exit(EXIT_FAILURE);
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

    exit(EXIT_FAILURE);
  }
/*.....................................................................*/

/*... advection*/        
  else if(!strcmp(word,help[4])){   
    printf("Ex1:\n"); 
    printf("diffusion 2 Vel Orthogonal Pres Orthogonal \n"); 
    printf("Options:\n");
    for(i=0;i<iDiff;i++)
      printf("%3d - %s\n",i+1,fDif[i]);
    exit(EXIT_FAILURE);
  }
/*.....................................................................*/

/*... advection*/        
  else{   
    printf("Options:\n");
    for(i=0;i<iHelp;i++)
      printf("%3d - %s\n",i+1,help[i]);
    exit(EXIT_FAILURE);
  }
/*.....................................................................*/

  exit(EXIT_FAILURE);
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

  short k = 0, ii=0,jj;

  do{
    ii++;
    readMacro(file,word,false);
    convStringLower(word);
    jj = 0;
/*... Smagorinsky*/
    if(!strcmp(word,turb[0])){ 
      k++;
      t->typeMixed[FUNMODEL] = SMAGORINSKY;
      fscanf(file,"%lf",&t->cf);  
      if(!mpiVar.myId) 
        fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[0],t->cf);    
    }
/*...................................................................*/ 

/*... Wale*/
    else if(!strcmp(word,turb[1])){ 
      k++;
      t->typeMixed[FUNMODEL] = WALEMODEL;
      fscanf(file,"%lf",&t->cf);    
      if(!mpiVar.myId)         
        fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[1],t->cf);    
    }
/*...................................................................*/ 

/*... sigmaModel*/
    else if(!strcmp(word,turb[2])){
      k++; 
      t->typeMixed[FUNMODEL]  = SIGMAMODEL;
      fscanf(file,"%lf",&t->cf);   
      if(!mpiVar.myId){ 
        fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[2],t->cf);
      }
    }
/*...................................................................*/ 

/*... bardina*/
    else if(!strcmp(word,turb[3])){
      k++;   
      t->typeMixed[ESTMODEL] = BARDINA;
      fscanf(file,"%lf",&t->cs);  
      if(!mpiVar.myId){ 
         fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[3],t->cs);
      }
    }
/*...................................................................*/  

/*... clark*/
    else if(!strcmp(word,turb[4])){
      k++;  
      t->typeMixed[ESTMODEL] = CLARK;
      fscanf(file,"%lf",&t->cs);  
      if(!mpiVar.myId){
        fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[4],t->cs);
      }
    } 
/*...................................................................*/ 

/*... bardinaMod*/
    else if(!strcmp(word,turb[5])){
      k++;   
      t->typeMixed[ESTMODEL] = BARDINAMOD;
      fscanf(file,"%lf",&t->cs);  
      if(!mpiVar.myId){ 
        fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[5],t->cs);
      }
    }
/*...................................................................*/ 
  
  }while( k!=2 && ii != 10);


  if (ii == 10)
    ERRO_GERAL(__FILE__,__func__,__LINE__,"Erro na leitura");

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

  short k = 0, ii=0,jj;

  readMacro(file,word,false);
  convStringLower(word);
/*... Smagorinsky*/
  if(!strcmp(word,turb[0])){ 
    k++;
    t->typeLes             = LESFUNCMODEL;
    t->typeMixed[FUNMODEL] = SMAGORINSKY;
    fscanf(file,"%lf",&t->cf);  
    if(!mpiVar.myId) 
      fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[0],t->cf);    
  }
/*...................................................................*/ 

/*... SigmaModel*/
  else if(!strcmp(word,turb[1])){ 
    k++;
    t->typeLes             = LESFUNCMODEL;
    t->typeMixed[FUNMODEL] = SIGMAMODEL;
    fscanf(file,"%lf",&t->cf);  
    if(!mpiVar.myId) 
      fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[1],t->cf);    
  }
/*...................................................................*/

/*... WaleModel*/
  else if(!strcmp(word,turb[2])){ 
    k++;
    t->typeLes             = LESFUNCMODEL;
    t->typeMixed[FUNMODEL] = WALEMODEL;
    fscanf(file,"%lf",&t->cf);  
    if(!mpiVar.myId) 
      fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[2],t->cf);    
  }
/*...................................................................*/ 

/*... VREMAN*/
  else if(!strcmp(word,turb[3])){ 
    k++;
    t->typeLes             = LESFUNCMODEL;
    t->typeMixed[FUNMODEL] = VREMAN;
    fscanf(file,"%lf",&t->cf);  
    if(!mpiVar.myId) 
      fprintf(fileLogExc,"%-20s: Cf = %lf\n", turb[3],t->cf);    
  }
/*...................................................................*/ 
}
/**********************************************************************/