#include<Simple.h>

static DOUBLE trunkNumber(DOUBLE const a);
void testeSist(INT *ia      ,INT *ja
              ,double *au   ,double *ad
              ,double *al   ,double *b
              ,INT const neq,bool const unsym);

/********************************************************************* 
 * Data de criacao    : 17/07/2016                                   *
 * Data de modificaco : 27/09/2019                                   * 
 *-------------------------------------------------------------------* 
 * SIMPLESOLVER3D: metodo simple e simpleC para escoamentos 3D       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void simpleSolver(Memoria *m        
                   ,Loads *loadsVel   ,Loads *loadsPres 
                   ,MassEqModel *eMass, MomentumModel *momentumModel
                   ,Turbulence *tModel 
                   ,Mesh *mesh0       ,Mesh *mesh       
                   ,SistEq *sistEqVel ,SistEq *sistEqPres
                   ,Solv *solvVel     ,Solv *solvPres 
                   ,Simple *sp
                   ,Scheme *sc        ,PartMesh *pMesh 
                   ,FileOpt *opt       ,char *preName 
                   ,char *nameOut     ,FILE *fileOut)
{
  FILE *fStop=NULL;
	short unsigned ndfVel = mesh->ndfF-1;
  short unsigned conv;
  int itSimple;
  short unsigned kZeroVel  = sp->vel.k;
  short unsigned kZeroPres = sp->mass.k;
  short typeResidual[2];
  INT jj = 1;
  DOUBLE time,timei;
/*...*/
  DOUBLE rU[3],rU0[3],tmp,tb[3],rMass0,rMass;
  DOUBLE *rCellPc;
/*...*/
  DOUBLE tolSimpleU1,tolSimpleU2,tolSimpleU3,tolSimpleMass;
/*...*/
  bool xMomentum ,yMomentum ,zMomentum ,pCor;
  bool relResVel, relResMass;
  bool fPrint  = false;
  DOUBLE cfl,reynolds;
  bool fParameter[2];

/*...*/
  time = getTimeC();

/*... vel*/
  typeResidual[0] = sp->vel.type;
  relResVel       = sp->vel.fRel;
/*... conservacao de massa*/
  typeResidual[1] = sp->mass.type; 
  relResMass      = sp->mass.fRel;
/*...*/
  tolSimpleU1   = sp->vel.tol[0];
  tolSimpleU2   = sp->vel.tol[1];
  tolSimpleU3   = sp->vel.tol[2];
  tolSimpleMass = sp->mass.tol[0];
/*...................................................................*/

/*...*/
  rMass0 = 1.e0;
  rMass  = 0.e0;
  rU[0]  = rU[1]  = rU[2]  = 0.e0;
  rU0[0] = rU0[1] = rU0[2] = 1.e0;
  conv = 0;
  xMomentum = yMomentum = zMomentum = true;
  tmp = tb[0] = tb[1] = tb[2] = 0.e0;
  rCellPc = mesh->elm.rCellPres;
/*...................................................................*/

/*...*/
  zero(sistEqVel->b0 ,sistEqVel->neqNov*ndfVel,DOUBLEC);
  zero(sistEqPres->b0,sistEqPres->neqNov      ,DOUBLEC);
/*...................................................................*/

/*... restricoes por centro de celula u0 e cargas por volume b0*/
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
  cellPloadSimple(loadsPres            ,mesh->elm.geom.cc 
                 ,mesh->elm.faceRpres    
                 ,mesh->elm.geom.volume
                 ,sistEqVel->id        ,sistEqPres->id
                 ,mesh->elm.vel        ,mesh->elm.pressure 
                 ,sistEqVel->b0        ,sistEqPres->b0
                 ,mesh->numelNov       ,ndfVel
                 ,mesh->ndm            ,mesh->maxViz);
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
/*...................................................................*/

/*... discretizacao temporal*/
  if(sc->ddt.flag){
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
/*...*/
    cellTransientSimple(mesh->elm.geom.volume   ,sistEqVel->id     
                       ,mesh->elm.vel0          ,mesh->elm.vel 
                       ,mesh->elm.densityFluid  ,sistEqVel->b0
                       ,sc->ddt                  ,sistEqVel->neqNov
                       ,mesh->numelNov          ,ndfVel        
                       ,true);
/*... vel(n-1) = vel(n)*/
    alphaProdVector(1.e0,mesh->elm.vel
                   ,mesh->numel*ndfVel ,mesh->elm.vel0); 
/*...................................................................*/
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
  }
/*...................................................................*/

/*... arquivo de log*/
  if (opt->fItPlot)
    fprintf(opt->fileItPlot[FITPLOTSIMPLE]
           , "istep = %d time = %lf\n",sc->ddt.timeStep,sc->ddt.t);
/*...................................................................*/

/*...*/
  for(itSimple=0;itSimple<sp->maxIt;itSimple++){
/*...*/
     if ( (fStop=fopen("stopSimple.mvf","r")) !=NULL){
		   fclose(fStop);			
       break;
     } 
/*...................................................................*/

/*... pressao-velociade*/
    velPresCoupling(m           , NULL
                    , loadsVel  , loadsPres
                    , eMass     , momentumModel
                    , tModel
                    , mesh
                    , sistEqVel , sistEqPres
                    , solvVel   , solvPres
                    , sp        , sc
                    , pMesh     , rCellPc
                    , &xMomentum, &yMomentum
                    , &zMomentum, &pCor
                    , fPrint    , itSimple);     
/*...................................................................*/
      
/*... residual*/
     tm.residualSimple = getTimeC() - tm.residualSimple;
     residualSimple(mesh->elm.vel 
                   ,mesh->elm.rCellVel,rCellPc  
                   ,sistEqVel->ad           
                   ,rU                ,&rMass
                   ,sistEqVel->id     
                   ,mesh->numelNov    ,sistEqVel->neq  
                   ,mesh->ndm         ,typeResidual    );  
     tm.residualSimple = getTimeC() - tm.residualSimple;
/*...................................................................*/

/*...*/
     if (itSimple == kZeroPres && relResMass &&  pCor) rMass0 = rMass;
     if (itSimple == kZeroVel && relResVel)
     {
       if (xMomentum) rU0[0] = rU[0];
       if (yMomentum) rU0[1] = rU[1];
       if (zMomentum && ndfVel == 3) rU0[2] = rU[2];
     }
/*...................................................................*/

/*...*/
     timei = getTimeC() -time;
/*... arquivo de log*/
     if(opt->fItPlot)
     {
       if (ndfVel == 3)
         fprintf(opt->fileItPlot[FITPLOTSIMPLE]
                , "%d %20.8e %20.8e %20.8e %20.8e\n"
                ,itSimple + 1, rU[0], rU[1], rU[2], rMass);
       else
         fprintf(opt->fileItPlot[FITPLOTSIMPLE]
                ,"%d %20.8e %20.8e %20.8e\n"
                ,itSimple+1,rU[0],rU[1],rMass);
     }
/*...................................................................*/

/*...*/
     if( jj == sp->pSimple) {
       jj = 0; 
       printf("It simple: %d \n",itSimple+1);
       printf("CPU Time(s)  : %lf \n", timei);
       printf("Residuo:\n");
       printf("%-25s : %20.8e\n","conservacao da massa",rMass/rMass0);
       printf("%-25s : %20.8e\n","momentum x1", rU[0] / rU0[0]);
       printf("%-25s : %20.8e\n","momentum x2", rU[1] / rU0[1]);
       if (ndfVel == 3)
          printf("%-25s : %20.8e\n","momentum x3", rU[2] / rU0[2]);
     } 
     jj++; 
/*...................................................................*/

/*... Massa*/
     conv = 0;
     if(rMass/rMass0<tolSimpleMass || rMass<SZERO) conv++;
/*... Vel*/
     if(rU[0]/rU0[0]<tolSimpleU1 || rU[0]<SZERO) conv++;
     if(rU[1]/rU0[1]<tolSimpleU2 || rU[1]<SZERO) conv++;
     if(rU[2]/rU0[2]<tolSimpleU3 || rU[2]<SZERO) conv++;
/*..*/
     if(conv == 4 ) break;
/*...................................................................*/

  }
/*...................................................................*/
  timei = getTimeC() -time;

/*...*/  
  fParameter[0] = true;
  fParameter[1] = true;
  parameterCell(mesh->elm.vel           ,mesh->elm.material.prop
               ,mesh->elm.densityFluid  ,mesh->elm.geom.volume 
               ,mesh->elm.mat  
               ,&cfl                    ,&reynolds
               ,fParameter              ,sc->ddt.dt[0]
               ,mesh->numelNov          ,mesh->ndm);
/*...................................................................*/

/*...*/
  printf("It simple: %d \n",itSimple+1);
  printf("Time(s)  : %lf \n",timei);
  printf("Reynolds: %lf\n",reynolds);
  if(sc->ddt.flag)
    printf("%-20s : %13.6lf\n","CFL", cfl);
  printf("%-25s : %15s %20s\n","Residuo:","init","final");
  printf("%-25s : %20.8e %20.8e\n","conservacao da massa", rMass0, rMass);
  printf("%-25s : %20.8e %20.8e\n","momentum x1",rU0[0], rU[0]);
  printf("%-25s : %20.8e %20.8e\n","momentum x2 ",rU0[1], rU[1]);
  if (ndfVel == 3)
    printf("%-25s : %20.8e %20.8e\n","momentum x3", rU0[2], rU[2]);
/*...................................................................*/

} 
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 30/07/2018                                   *
 * Data de modificaco : 23/07/2019                                   *
 *-------------------------------------------------------------------*
 * combustionSolver: metodo simple e simpleC para escoamentos        * 
 * 2D/3D termo ativados                                              *
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
void combustionSolver(Memoria *m        , PropVarFluid *propF     
                  , Loads *loadsVel     , Loads *loadsPres
                  , Loads *loadsEnergy  , Loads *loadsKturb
                  , Loads *loadsComb    
                  , EnergyModel *eModel , Combustion *cModel
                  , MassEqModel *eMass  , MomentumModel *momentumModel
                  , Turbulence *tModel  , ThermoDynamic *thDynamic
                  , Mesh *mesh0         , Mesh *mesh
                  , SistEq *sistEqVel   , SistEq *sistEqPres
                  , SistEq *sistEqEnergy, SistEq *sistEqKturb
                  , SistEq *sistEqComb
                  , Solv *solvVel       , Solv *solvPres
                  , Solv *solvEnergy    , Solv *solvKturb
                  , Solv *solvComb     
                  , Simple *sp          , Scheme *sc          
                  , PartMesh *pMesh     , Mean *media
                  , FileOpt *opt        , char *preName
                  , char *nameOut       , FILE *fileOut) {
  FILE *fStop = NULL;
  short unsigned ndfVel = mesh->ndfFt - 2;
  short unsigned nComb = cModel->nComb,
                 nSp   = cModel->nOfSpecies; 
  short unsigned conv, i;
  short unsigned itSimple;
  short unsigned kZeroVel    = sp->vel.k,
                 kZeroPres   = sp->mass.k,
                 kZeroEnergy = sp->energy.k,
                 kZeroComb   = sp->z.k;
  short typeResidual[4];
  INT jj = 1;
  DOUBLE time, timei;
  DOUBLE *rCellPc;

/*...*/
  DOUBLE rU[3], rU0[3], rMass0, rMass, rComb[MAXSPECIES],rComb0[MAXSPECIES]
       , rEnergy0, rEnergy;
/*...*/
  DOUBLE tolSimpleU1, tolSimpleU2, tolSimpleU3, tolComb
       , tolSimpleMass, tolSimpleEnergy;
/*...*/
  bool xMomentum, yMomentum, zMomentum, pCor, fEnergy, fComb[MAX_COMB];
  bool relRes = false;
  bool fPrint = false;
  bool fDensity       = propF->fDensity,
    fSheat            = propF->fSpecificHeat,
    fDvisc            = propF->fDynamicViscosity,
    fTcond            = propF->fThermalConductivity,
    fDiff             = propF->fDiffusion,
    fDensityRef       = thDynamic->fDensityRef,
    fPresRef          = thDynamic->fPresTh,
    fDeltaTimeDynamic = sc->ddt.fDynamic;
  DOUBLE cfl, reynolds, peclet,  deltaMass, prMax;
  bool fParameter[10];

/*...*/
  time = getTimeC();

/*...*/
//relRes          = true;
/*... vel*/
  typeResidual[0] = sp->vel.type;
/*... energia*/  
  typeResidual[1] = sp->energy.type;
/*... conservao de especie*/ 
  typeResidual[2] = sp->z.type;
/*... conservacao de massa*/
  typeResidual[3] = sp->mass.type;      

/*...*/
  tolSimpleU1     = sp->vel.tol[0];
  tolSimpleU2     = sp->vel.tol[1];
  tolSimpleU3     = sp->vel.tol[2];
  tolSimpleMass   = sp->mass.tol[0];
  tolSimpleEnergy = sp->energy.tol[0];
  tolComb         = sp->z.tol[0];
/*...................................................................*/

/*...*/
  rMass0  = 1.e0;
  rMass   = 0.e0;
  rU[0]   = rU[1] = rU[2] = 0.e0;
  rU0[0]  = rU0[1] = rU0[2] = 1.e0;
  for (i = 0; i < nComb; i++)
  {
    rComb[i]  = 0.e0; 
    rComb0[i] = 1.e0;
  }
   
  rEnergy0 = 1.e0;
  rEnergy  = 0.e0;
  conv = 0;
  xMomentum = yMomentum = zMomentum = true;
  rCellPc = mesh->elm.rCellPres;
/*...................................................................*/

/*...*/
  zero(sistEqVel->b0   ,sistEqVel->neqNov*ndfVel,DOUBLEC);
  zero(sistEqPres->b0  ,sistEqPres->neqNov      ,DOUBLEC);
  zero(sistEqEnergy->b0,sistEqEnergy->neqNov    ,DOUBLEC);
  zero(sistEqComb->b0  ,sistEqComb->neqNov*nComb,DOUBLEC)
/*...................................................................*/

/*... restricoes por centro de celula vel e cargas por volume b0*/
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
  cellPloadSimple(loadsPres            ,mesh->elm.geom.cc
                 ,mesh->elm.faceRpres  ,mesh->elm.geom.volume
                 ,sistEqVel->id        ,sistEqPres->id
                 ,mesh->elm.vel        ,mesh->elm.pressure
                 ,sistEqVel->b0        ,sistEqPres->b0
                 ,mesh->numelNov       ,ndfVel
                 ,mesh->ndm            ,mesh->maxViz);
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
/*...................................................................*/

/*... restricoes por centro de celula energy e cargas por volume b0*/
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
  cellPload(loadsEnergy            ,mesh->elm.geom.cc
           ,mesh->elm.faceRenergy  
           ,mesh->elm.geom.volume  ,sistEqEnergy->id
           ,mesh->elm.energy       ,sistEqEnergy->b0
           ,mesh->numelNov         ,1        
           ,mesh->ndm              ,mesh->maxViz);
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
/*...................................................................*/

/* ... restricoes por centro de celula kTurb e cargas por volume b0*/
  if (tModel->typeLes == LESFUNCMODELONEEQK) 
  {
    zero(sistEqKturb->b0 ,sistEqKturb->neqNov     ,DOUBLEC);
    tm.turbulence = getTimeC() - tm.turbulence;
    cellPload(loadsKturb         ,mesh->elm.geom.cc
         ,mesh->elm.faceReKturb  
         ,mesh->elm.geom.volume  ,sistEqKturb->id
         ,mesh->elm.energy       ,sistEqKturb->b0
         ,mesh->numelNov         ,1         
         ,mesh->ndm              ,mesh->maxViz);
    tm.turbulence = getTimeC() - tm.turbulence;
  }
/*...................................................................*/

/*... discretizacao temporal*/
  if (sc->ddt.flag) 
  {
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
/*... velocidade*/
    cellTransientSimple(mesh->elm.geom.volume ,sistEqVel->id
                       ,mesh->elm.vel0        ,mesh->elm.vel
                       ,mesh->elm.densityFluid,sistEqVel->b0
                       ,sc->ddt               ,sistEqVel->neqNov
                       ,mesh->numelNov        ,ndfVel
                       ,true);  
/*...................................................................*/

/*... Energia*/
    if(eModel->fTemperature)
      cellTransientEnergy(mesh->elm.geom.volume  ,sistEqEnergy->id
                         ,mesh->elm.energy0      ,mesh->elm.energy
                         ,mesh->elm.densityFluid ,mesh->elm.specificHeat
                         ,sistEqEnergy->b0
                         ,sc->ddt                 ,mesh->numelNov
                         ,true);
    else
      cellTransient(mesh->elm.geom.volume  ,sistEqEnergy->id
                   ,mesh->elm.energy0      ,mesh->elm.energy
                   ,mesh->elm.densityFluid ,sistEqEnergy->b0
                   ,sc->ddt                ,mesh->numelNov
                   ,1                      ,true);  
/*...................................................................*/

/*... energia cinetica turbulenta*/
    if (tModel->typeLes == LESFUNCMODELONEEQK) 
    {
      tm.turbulence = getTimeC() - tm.turbulence;
      cellTransient(mesh->elm.geom.volume  ,sistEqKturb->id
                   ,mesh->elm.kTurb0       ,mesh->elm.kTurb  
                   ,mesh->elm.densityFluid ,sistEqKturb->b0
                   ,sc->ddt                ,mesh->numelNov
                   ,1                      ,true);
/*... kTurb(n-1) = kTurb(n)*/
      alphaProdVector(1.e0              ,mesh->elm.kTurb
                     ,mesh->numel       ,mesh->elm.kTurb0);
      tm.turbulence = getTimeC() - tm.turbulence;
    }
/*...................................................................*/

/*... vel(n-1) = vel(n)*/
    alphaProdVector(1.e0              ,mesh->elm.vel 
                   ,mesh->numel*ndfVel,mesh->elm.vel0);
/*... energy(n-1) = energy(n)*/
    alphaProdVector(1.e0       ,mesh->elm.energy
                   ,mesh->numel,mesh->elm.energy0);    
/*... pres(n-1) = pres(n)*/
    alphaProdVector(1.e0       ,mesh->elm.pressure
                   ,mesh->numel,mesh->elm.pressure0);
/*...................................................................*/
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;

/*...*/
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
/*... zComb*/
    cellTransientSimple(mesh->elm.geom.volume , sistEqComb->id
                      , mesh->elm.zComb0      , mesh->elm.zComb
                      , mesh->elm.densityFluid, sistEqComb->b0
                      , sc->ddt               , sistEqComb->neqNov
                      , mesh->numelNov        , nComb
                      , true);    
/*...................................................................*/

/*... z(n-1) = z(n)*/
    alphaProdVector(1.e0             , mesh->elm.zComb
                  , mesh->numel*nComb, mesh->elm.zComb0);
/*...................................................................*/

/*... y(n-1) = y(n)*/
    alphaProdVector(1.e0             , mesh->elm.yFrac
                  , mesh->numel*nSp  , mesh->elm.yFrac0);
/*...................................................................*/
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
  }
/*...................................................................*/

/*... arquivo de log*/
  if (opt->fItPlot)
    fprintf(opt->fileItPlot[FITPLOTSIMPLE]
           , "istep = %d time = %lf\n",sc->ddt.timeStep,sc->ddt.t);
/*...................................................................*/

/*...*/
  for (itSimple = 0; itSimple<sp->maxIt; itSimple++) 
  {

/*...*/
    if ((fStop = fopen("stopSimple.mvf", "r")) != NULL) 
    {
      fclose(fStop);      
      break;
    }
/*...................................................................*/

/*... pressao-velociade*/
    velPresCoupling(m           , propF
                    , loadsVel  , loadsPres
                    , eMass     , momentumModel
                    , tModel    , mesh
                    , sistEqVel , sistEqPres
                    , solvVel   , solvPres
                    , sp        , sc
                    , pMesh     , rCellPc
                    , &xMomentum, &yMomentum
                    , &zMomentum, &pCor
                    , fPrint    , itSimple);     
/*...................................................................*/

/*... modelo de turbulencia*/
    if(tModel->fTurb){
      tm.turbulence = getTimeC() - tm.turbulence;
      turbulence(m    
               , loadsKturb, loadsVel  
               , pMesh     , tModel
               , mesh      , sc
               , sp        , sistEqKturb
               , solvKturb , ndfVel);  
      tm.turbulence = getTimeC() - tm.turbulence;
    }
/*...................................................................*/

/*... modelo de combustao*/
   if (fPrint) printf("Modelo de combustao:\n");
   combustionModel(m          , propF
                   , loadsComb , loadsVel
                   , tModel    , cModel
                   , eModel    , mesh   
                   , sistEqComb, solvComb 
                   , sp        , sc       
                   , pMesh     , opt
                   , fComb     , itSimple  );              
/*...................................................................*/

/*... equacao de energia*/
    if (fPrint) printf("Conservacao de Energia:\n");
    fEnergy = energyEquation(m            , propF 
                           , loadsVel     , loadsEnergy  
                           , eModel       , tModel      
                           , cModel       , thDynamic
                           , mesh          
                           , sistEqEnergy, solvEnergy
                           , sp          , sc
                           , pMesh);    
/*...................................................................*/

/*... residual*/
    tm.residualSimple = getTimeC() - tm.residualSimple;
    residualCombustion(mesh->elm.vel      ,mesh->elm.energy
            ,mesh->elm.zComb
            ,mesh->elm.rCellVel           ,rCellPc
            ,mesh->elm.rCellEnergy        ,mesh->elm.rCellComb      
            ,sistEqVel->ad                ,sistEqEnergy->ad 
            ,sistEqComb->ad   
            ,rU                           ,&rMass
            ,&rEnergy                     ,rComb    
            ,sistEqVel->id                ,sistEqEnergy->id
            ,sistEqComb->id   
            ,mesh->numelNov               ,sistEqVel->neqNov
            ,sistEqComb->neqNov
            ,mesh->ndm                    ,nComb       
            ,typeResidual ); 
    tm.residualSimple = getTimeC() - tm.residualSimple;
/*...................................................................*/

/*...*/
    tm.tempFromTheEnergy = getTimeC() - tm.tempFromTheEnergy;  
    getTempFromTheEnergyMix(propF               ,mesh->elm.yFrac
                        ,mesh->elm.temp         ,mesh->elm.energy
                        ,mesh->numel            ,cModel->nOfSpecies
                        ,eModel->fTemperature   ,eModel->fKelvin
                        ,ompVar.fUpdate         ,ompVar.nThreadsUpdate);  
    tm.tempFromTheEnergy = getTimeC() - tm.tempFromTheEnergy;
/*...................................................................*/

/*... pressa de referencia*/
    if(fPresRef)
      presRefMix(cModel
            , mesh->elm.temp0      , mesh->elm.temp
            , mesh->elm.yFrac0      , mesh->elm.yFrac
            , mesh->elm.geom.volume, thDynamic->pTh  
            , mesh->numelNov       , eModel->fKelvin);    
/*...................................................................*/

/*...*/
    tm.updateProp = getTimeC() - tm.updateProp;
    if(fDensity)
      updateMixDensity(&propF->den          , cModel
                 , mesh->elm.temp           , mesh->elm.pressure
                 , mesh->elm.densityFluid   , mesh->elm.yFrac         
                 , sp->alphaDensity         , eModel->fKelvin      
                 , mesh->numel              , PROP_UPDATE_NL_LOOP
                 , ompVar.fUpdate           , ompVar.nThreadsUpdate);    
    if(fSheat)
      updateMixSpecificHeat(&propF->sHeat
                         , mesh->elm.temp        , mesh->elm.yFrac
                         , mesh->elm.specificHeat, cModel->nOfSpecies 
                         , eModel->fKelvin     
                         , mesh->numel           , PROP_UPDATE_NL_LOOP
                         , ompVar.fUpdate        , ompVar.nThreadsUpdate);  
    if(fDvisc)
      updateMixDynamicViscosity(&propF->dVisc    ,cModel
                          ,mesh->elm.temp        ,mesh->elm.yFrac 
                          ,mesh->elm.dViscosity  ,cModel->nOfSpecies   
                          ,eModel->fKelvin       , mesh->numel
                          ,ompVar.fUpdate        , ompVar.nThreadsUpdate);
    if(fTcond)
      updateMixDynamicThermalCond(propF           ,cModel 
                          ,mesh->elm.temp         , mesh->elm.yFrac
                          ,mesh->elm.tConductivity,cModel->nOfSpecies
                          ,eModel->fKelvin        , mesh->numel
                          ,ompVar.fUpdate         , ompVar.nThreadsUpdate);
    
    if(fDiff)
      updateMixDiffusion(propF             , cModel 
                       ,mesh->elm.temp     , mesh->elm.yFrac
                       ,mesh->elm.cDiffComb, cModel->nOfSpecies 
                       ,cModel->nComb       
                       ,eModel->fKelvin    , mesh->numel
                       ,ompVar.fUpdate     , ompVar.nThreadsUpdate); 
    tm.updateProp = getTimeC() - tm.updateProp;   
/*...................................................................*/

/*...*/
    if (itSimple == kZeroPres && relRes &&  pCor) rMass0 = rMass;
    if (itSimple == kZeroEnergy && relRes && fEnergy) rEnergy0 = rEnergy;
    if (itSimple == kZeroVel && relRes)
    {
      if (xMomentum) rU0[0] = rU[0];
      if (yMomentum) rU0[1] = rU[1];
      if (zMomentum && ndfVel == 3) rU0[2] = rU[2];
    }
    if(itSimple == kZeroComb && relRes)
    {
      for(i=0;i<nComb;i++)
        if(fComb[i]) rComb0[i] = rComb[i];
    }
    conv = 0;  
/*...................................................................*/

/*...*/
    timei = getTimeC() - time;
/*... arquivo de log*/
    if (opt->fItPlot)
    {
      if (ndfVel == 3 && !mpiVar.myId)
      {
        fprintf(opt->fileItPlot[FITPLOTSIMPLE]
          , "%10d %20.8e %20.8e %20.8e %20.8e %20.8e "
          , itSimple + 1, rU[0], rU[1], rU[2], rMass, rEnergy);
        for(i=0;i<nComb;i++)
          fprintf(opt->fileItPlot[FITPLOTSIMPLE]
          , " %20.8e ",  rComb[i]);
        fprintf(opt->fileItPlot[FITPLOTSIMPLE],"\n");
      }
    }
/*...................................................................*/

/*...*/
    if (jj == sp->pSimple) 
    {
      jj = 0;
      if(!mpiVar.myId)
      {
        printf("It simple: %d \n", itSimple + 1);
        printf("CPU Time(s)  : %lf \n", timei);
        printf("Residuo:\n");
        printf("%-25s : %20.8e\n","conservacao da massa", rMass / rMass0);
        printf("%-25s : %20.8e\n","momentum x1", rU[0] / rU0[0]);
        printf("%-25s : %20.8e\n","momentum x2", rU[1] / rU0[1]);
        if (ndfVel == 3)
          printf("%-25s : %20.8e\n","momentum x3", rU[2] / rU0[2]);
        printf("%-25s : %20.8e\n","conservacao da energia",rEnergy/rEnergy0);
        for(i=0;i<nComb;i++)
          printf("%-25s : %20.8e\n",cModel->chem.sp[i].name,rComb[i]/rComb0[i]);
      }
    }
    jj++; 
/*...................................................................*/

/*... Energia*/
    if(rEnergy/rEnergy0<tolSimpleEnergy || rEnergy<SZERO) conv++;
/*... Massa*/
    if(rMass/rMass0<tolSimpleMass || rMass<SZERO) conv++;
/*... Vel*/
    if(rU[0]/rU0[0]<tolSimpleU1 || rU[0]<SZERO) conv++;
    if(rU[1]/rU0[1]<tolSimpleU2 || rU[1]<SZERO) conv++;
    if(rU[2]/rU0[2]<tolSimpleU3 || rU[2]<SZERO) conv++;
/*... Combustion*/
    for(i=0;i<nComb;i++)
      if(rComb[i]/rComb0[i]<tolComb || rComb[i]<SZERO) conv++;
/*..*/
    if(conv == (5 + nComb)) break;
/*...................................................................*/

  }
/*...................................................................*/
  timei = getTimeC() - time;

/*...*/
  fParameter[0] = true;
  fParameter[1] = true;
  fParameter[2] = true;
  fParameter[3] = true;
  fParameter[4] = true;
  parameterCellLm(mesh->elm.vel         , mesh->elm.material.prop
               , mesh->elm.densityFluid , mesh->elm.specificHeat
               , mesh->elm.tConductivity, mesh->elm.dViscosity
               , mesh->elm.geom.volume  , mesh->elm.mat            
               , &cfl                   , &reynolds
               , &peclet                , &mesh->mass[2]
               , &prMax
               , fParameter             , sc->ddt.dt[0]
               , mesh->numelNov         , mesh->ndm);  
/*...................................................................*/

/*...*/
  if(fDeltaTimeDynamic)
    dynamicDeltat(mesh->elm.vel          , mesh->elm.geom.volume 
                 , mesh->elm.densityFluid , mesh->elm.specificHeat
                 , mesh->elm.tConductivity, mesh->elm.dViscosity
                 , sc->ddt.dt             , mesh->numelNov
                 , mesh->ndm              , CFL);  
/*...................................................................*/

/*... guardando as propriedades para o proximo passo*/
  tm.updateProp = getTimeC() - tm.updateProp;
  if(fDensity)
    updateMixDensity(&propF->den         , cModel
                 , mesh->elm.temp        , mesh->elm.pressure
                 , mesh->elm.densityFluid, mesh->elm.yFrac         
                 , sp->alphaDensity      , eModel->fKelvin      
                 , mesh->numelNov        , PROP_UPDATE_OLD_TIME
                 , ompVar.fUpdate        , ompVar.nThreadsUpdate); 
  if(fSheat)
    updateMixSpecificHeat(&propF->sHeat
                         , mesh->elm.temp        , mesh->elm.yFrac
                         , mesh->elm.specificHeat, cModel->nOfSpecies
                         , eModel->fKelvin     
                         , mesh->numel           , PROP_UPDATE_OLD_TIME
                         , ompVar.fUpdate        , ompVar.nThreadsUpdate);
  tm.updateProp = getTimeC() - tm.updateProp;
/*...................................................................*/

/*... calcula a massa especifica de referencia*/
  if(fDensityRef)
    propF->densityRef = specificMassRef(mesh->elm.densityFluid
                                      , mesh->elm.geom.volume                  
                                      , mesh->numel);  
/*...................................................................*/

/*... calculo da taxa de massa atravessando o contorno aberto*/
  massFluxOpenDomain(loadsVel         , sc->ddt 
              , mesh->elm.cellFace    , mesh->face.owner
              , mesh->elm.faceRvel , mesh->elm.adj.nViz 
              , mesh->face.area       , mesh->face.normal
              , mesh->face.xm  
              , mesh->elm.densityFluid, mesh->elm.vel 
              , mesh->massInOut       , &deltaMass
              , mesh->numelNov        , mesh->ndm  
              , mesh->maxViz  );
  mesh->mass[1] += deltaMass;  
/*...................................................................*/

/*... temp(n) = temp(n+1)*/
  alphaProdVector(1.e0        ,mesh->elm.temp
                 ,mesh->numel ,mesh->elm.temp0);  
/*...................................................................*/

/*... pth(n-2) = pth(n-1)*/
  thDynamic->pTh[0] = thDynamic->pTh[1];  
/*... pth(n-1) = pth(n)*/
  thDynamic->pTh[1] = thDynamic->pTh[2];  
/*...................................................................*/

/*... Enegia total liberada*/
  cModel->totalHeat += totalHeatRealeseComb(mesh->elm.rateHeatReComb
                      , mesh->elm.geom.volume
                      , sc->ddt.dt[TIME_N], mesh->numelNov); 
/*...................................................................*/

/*...*/
  mesh->tempMax = maxArray(mesh->elm.temp,mesh->numelNov);
  mesh->tempMed = getVolumeMed(mesh->elm.temp,mesh->elm.geom.volume
                              ,mesh->numelNov);
/*...................................................................*/

/*...*/
  if(!mpiVar.myId)
  {
    printf("It simple: %d \n", itSimple + 1);
    printf("Time(s) : %lf \n", timei);
    if (sc->ddt.flag)
      printf("%-20s : %13.6lf\n","CFL", cfl);
    printf("%-20s : %13.6lf\n","Reynolds", reynolds);
    printf("%-20s : %13.6lf\n","Peclet"  , peclet);
    printf("%-20s : %13.6lf\n","Prandtl" , prMax);
    printf("%-20s : %13.6e\n","PresRef"             , thDynamic->pTh[2]);
    printf("%-20s : %13.6e\n","(Inc    ) Mass/Mass0", mesh->mass[1]/mesh->mass[0]);
    printf("%-20s : %13.6e\n","(Average) Mass/Mass0", mesh->mass[2]/mesh->mass[0]);
    printf("%-20s : %13.6e\n","MassIn              ", mesh->massInOut[0]);
    printf("%-20s : %13.6e\n","MassOut             ", mesh->massInOut[1]);
    printf("%-20s : %13.6e\n","Heat Combustion"     ,cModel->totalHeat);
    printf("%-20s : %13.6lf\n","Temperatue Max"     ,mesh->tempMax);
    printf("%-20s : %13.6lf\n","Temperatue Med"     ,mesh->tempMed);
    
    printf("%-25s : %15s %20s\n","Residuo:","init","final");
    printf("%-25s : %20.8e %20.8e\n","conservacao da massa", rMass0, rMass);
    printf("%-25s : %20.8e %20.8e\n","momentum x1",rU0[0], rU[0]);
    printf("%-25s : %20.8e %20.8e\n","momentum x2 ",rU0[1], rU[1]);
    if (ndfVel == 3)
      printf("%-25s : %20.8e %20.8e\n","momentum x3", rU0[2], rU[2]);
    printf("%-25s : %20.8e %20.8e\n","conservacao da energia", rEnergy0, rEnergy);
    for(i=0;i<nComb;i++)
    {
      printf("%-25s : %20.8e %20.8e\n",cModel->chem.sp[i].name
                                      ,rComb0[i],rComb[i]);
    }
/*...................................................................*/

/*...*/
    fprintf(opt->fileParameters
           ,"%9d %e %e %e %e %e %e %e %e %e %e %e %e\n"
                               , sc->ddt.timeStep  , sc->ddt.t
                               , cfl               , reynolds
                               , peclet            , thDynamic->pTh[2]
                               , mesh->mass[1]     , mesh->mass[2]
                               , mesh->massInOut[0], mesh->massInOut[1]
                               , cModel->totalHeat , mesh->tempMax
                               ,mesh->tempMed);
/*...................................................................*/
  }
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 21/08/2017                                   *
* Data de modificaco : 30/07/2018                                   *
*-------------------------------------------------------------------*
* SIMPLESOLVERENERGY: metodo simple e simpleC para escoamentos      *
* 2D/3D termo ativados                                              *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
*********************************************************************/
void simpleSolverLm(Memoria *m         , PropVarFluid *propF
                 , Loads *loadsVel     , Loads *loadsPres
                 , Loads *loadsEnergy  , Loads *loadsKturb
                 , EnergyModel *eModel , Combustion *cModel
                 , MassEqModel *eMass  , MomentumModel *momentumModel
                 , Turbulence *tModel  , ThermoDynamic *thDynamic
                 , Mesh *mesh0         , Mesh *mesh
                 , SistEq *sistEqVel   , SistEq *sistEqPres
                 , SistEq *sistEqEnergy, SistEq *sistEqKturb
                 , Solv *solvVel       , Solv *solvPres
                 , Solv *solvEnergy    , Solv *solvKturb
                 , Simple *sp          , Scheme *sc
                 , PartMesh *pMesh     , Mean *media
                 , FileOpt *opt        , char *preName
                 , char *nameOut       , FILE *fileOut) {
  FILE *fStop = NULL;
  short unsigned ndfVel = mesh->ndfFt - 2;
  short unsigned conv;
  short itSimple;
  short unsigned kZeroVel    = sp->vel.k,
                 kZeroPres   = sp->mass.k,
                 kZeroEnergy = sp->energy.k;
  short typeResidual[3];
  INT jj = 1;
  DOUBLE time, timei;
  DOUBLE *rCellPc;

/*...*/
  DOUBLE rU[3], rU0[3], rMass0, rMass
    , rEnergy0, rEnergy;
/*...*/
  DOUBLE tolSimpleU1, tolSimpleU2, tolSimpleU3
    , tolSimpleMass, tolSimpleEnergy;
/*...*/
  bool xMomentum, yMomentum, zMomentum, pCor, fEnergy;
  bool relResVel, relResMass, relResEnergy;
  bool fPrint = false;
  bool fDensity = propF->fDensity,
    fSheat = propF->fSpecificHeat,
    fDvisc = propF->fDynamicViscosity,
    fTcond = propF->fThermalConductivity,
    fDensityRef = thDynamic->fDensityRef,
    fPresRef = thDynamic->fPresTh,
    fDeltaTimeDynamic = sc->ddt.fDynamic;
  DOUBLE cfl, reynolds, peclet, deltaMass, prMax;
  bool fParameter[10];

  time = getTimeC();

/*...*/
/*... vel*/
  typeResidual[0] = sp->vel.type;
  relResVel       = sp->vel.fRel;
/*... conservacao de massa*/
  typeResidual[1] = sp->mass.type;
  relResMass      = sp->mass.fRel;
  /*... conservacao de energiamassa*/
  typeResidual[2] = sp->energy.type;
  relResEnergy    = sp->energy.fRel;

/*...*/
  tolSimpleU1     = sp->vel.tol[0];
  tolSimpleU2     = sp->vel.tol[1];
  tolSimpleU3     = sp->vel.tol[2];
  tolSimpleMass   = sp->mass.tol[0];
  tolSimpleEnergy = sp->energy.tol[0];
/*...................................................................*/

/*...*/
  rMass0   = 1.e0;
  rMass    = 0.e0;
  rU[0]    = rU[1] = rU[2] = 0.e0;
  rU0[0]   = rU0[1] = rU0[2] = 1.e0;
  rEnergy0 = 1.e0;
  rEnergy  = 0.e0;
  conv     = 0;
  fEnergy  = xMomentum = yMomentum = zMomentum = true;
  rCellPc  = mesh->elm.rCellPres;
/*...................................................................*/

/*...*/
  zero(sistEqVel->b0, sistEqVel->neqNov*ndfVel, DOUBLEC);
  zero(sistEqPres->b0, sistEqPres->neqNov, DOUBLEC);
  zero(sistEqEnergy->b0, sistEqEnergy->neqNov, DOUBLEC);
/*...................................................................*/

/*... restricoes por centro de celula u0 e cargas por volume b0*/
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
  cellPloadSimple(loadsPres            , mesh->elm.geom.cc
                , mesh->elm.faceRpres  , mesh->elm.geom.volume
                , sistEqVel->id        , sistEqPres->id
                , mesh->elm.vel        , mesh->elm.pressure
                , sistEqVel->b0        , sistEqPres->b0
                , mesh->numelNov       , ndfVel
                , mesh->ndm            , mesh->maxViz);
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
/*...................................................................*/

/*... restricoes por centro de celula u0 e cargas por volume b0*/
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
  cellPload(loadsEnergy          , mesh->elm.geom.cc
          , mesh->elm.faceRenergy
          , mesh->elm.geom.volume, sistEqEnergy->id
          , mesh->elm.energy     , sistEqEnergy->b0
          , mesh->numelNov       , 1
          , mesh->ndm            , mesh->maxViz);
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
/*...................................................................*/

/* ... restricoes por centro de celula u0 e cargas por volume b0*/
  if (tModel->typeLes == LESFUNCMODELONEEQK)
  {
    zero(sistEqKturb->b0, sistEqKturb->neqNov, DOUBLEC);
    tm.turbulence = getTimeC() - tm.turbulence;
    cellPload(loadsKturb           , mesh->elm.geom.cc
            , mesh->elm.faceReKturb
            , mesh->elm.geom.volume, sistEqKturb->id
            , mesh->elm.energy     , sistEqKturb->b0
            , mesh->numelNov       , 1
            , mesh->ndm            , mesh->maxViz);
    tm.turbulence = getTimeC() - tm.turbulence;
  }
/*...................................................................*/

/*... discretizacao temporal*/
  if (sc->ddt.flag)
  {
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
/*... velocidade*/
    cellTransientSimple(mesh->elm.geom.volume , sistEqVel->id
                      , mesh->elm.vel0        , mesh->elm.vel
                      , mesh->elm.densityFluid, sistEqVel->b0
                      , sc->ddt               , sistEqVel->neqNov
                      , mesh->numelNov        , ndfVel
                      , true);
/*...................................................................*/

/*... Energia*/
    if (eModel->fTemperature)
      cellTransientEnergy(mesh->elm.geom.volume , sistEqEnergy->id
                        , mesh->elm.energy0     , mesh->elm.energy
                        , mesh->elm.densityFluid, mesh->elm.specificHeat
                        , sistEqEnergy->b0
                        , sc->ddt               , mesh->numelNov
                        , true);
    else
      cellTransient(mesh->elm.geom.volume , sistEqEnergy->id
                  , mesh->elm.energy0     , mesh->elm.energy
                  , mesh->elm.densityFluid, sistEqEnergy->b0
                  , sc->ddt               , mesh->numelNov
                  , 1                     , true);
/*...................................................................*/

/*... energia cinetica turbulenta*/
    if (tModel->typeLes == LESFUNCMODELONEEQK)
    {
      tm.turbulence = getTimeC() - tm.turbulence;
      cellTransient(mesh->elm.geom.volume , sistEqKturb->id
                  , mesh->elm.kTurb0      , mesh->elm.kTurb
                  , mesh->elm.densityFluid, sistEqKturb->b0
                  , sc->ddt               , mesh->numelNov
                  , 1                     , true);
/*... kTurb(n-1) = kTurb(n)*/
      alphaProdVector(1.e0       , mesh->elm.kTurb
                    , mesh->numel, mesh->elm.kTurb0);
      tm.turbulence = getTimeC() - tm.turbulence;
    }
/*...................................................................*/

/*... vel(n-1) = vel(n)*/
    alphaProdVector(1.e0               , mesh->elm.vel
                   , mesh->numel*ndfVel, mesh->elm.vel0);
/*... energy(n-1) = energy(n)*/
    alphaProdVector(1.e0       , mesh->elm.energy
                  , mesh->numel, mesh->elm.energy0);
/*... pres(n-1) = pres(n)*/
    alphaProdVector(1.e0      , mesh->elm.pressure
                 , mesh->numel, mesh->elm.pressure0);
/*...................................................................*/
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
  }
/*...................................................................*/

/*...*/
  for (itSimple = 0; itSimple<sp->maxIt; itSimple++)
  {
/*...*/
    if ((fStop = fopen("stopSimple.mvf", "r")) != NULL)
    {
      fclose(fStop);
      break;
    }
/*...................................................................*/

/*... pressao-velociade*/
    velPresCoupling(m         , propF
                  , loadsVel  , loadsPres
                  , eMass     , momentumModel
                  , tModel    , mesh
                  , sistEqVel , sistEqPres
                  , solvVel   , solvPres
                  , sp        , sc
                  , pMesh     , rCellPc
                  , &xMomentum, &yMomentum
                  , &zMomentum, &pCor
                  , fPrint    , itSimple);
/*...................................................................*/

/*... modelo de turbulencia*/
    if (tModel->fTurb) 
    {
      tm.turbulence = getTimeC() - tm.turbulence;
      turbulence(m
               , loadsKturb, loadsVel
               , pMesh     , tModel
               , mesh      , sc
               , sp        , sistEqKturb
               , solvKturb , ndfVel);

      tm.turbulence = getTimeC() - tm.turbulence;
    }
/*...................................................................*/

/*... equacao de energia*/
    if (fPrint) printf("Consercao de Energia:\n");
    fEnergy = energyEquation(m           , propF
                           , loadsVel    , loadsEnergy
                           , eModel      , tModel  
                           , cModel      , thDynamic
                           , mesh
                           , sistEqEnergy, solvEnergy
                           , sp          , sc
                           , pMesh);
/*...................................................................*/

/*... residual*/
    residualSimpleLm(mesh->elm.vel        , mesh->elm.energy
                   , mesh->elm.rCellVel   , rCellPc
                   , mesh->elm.rCellEnergy
                   , sistEqVel->ad        , sistEqEnergy->ad
                   , rU                   , &rMass
                   , &rEnergy
                   , sistEqVel->id        , sistEqEnergy->id
                   , mesh->numelNov       , sistEqVel->neqNov
                   , mesh->ndm            , typeResidual);
/*...................................................................*/

/*...*/
    tm.tempFromTheEnergy = getTimeC() - tm.tempFromTheEnergy;
    getTempForEnergy( propF
                    , mesh->elm.temp         , mesh->elm.energy
                    , mesh->numel            , eModel->fTemperature
                    , fSheat                 , eModel->fKelvin
                    , ompVar.fUpdate         , ompVar.nThreadsUpdate);
    tm.tempFromTheEnergy = getTimeC() - tm.tempFromTheEnergy;
/*...................................................................*/

/*...*/
    if (fDensity)
      updateDensity(propF
                  , mesh->elm.temp, mesh->elm.pressure0
                  , mesh->elm.densityFluid
                  , sp->alphaDensity, eModel->fKelvin
                  , mesh->numel, PROP_UPDATE_NL_LOOP);
    if (fSheat)
      updateSpecificHeat(&propF->sHeat
                       , mesh->elm.temp, mesh->elm.specificHeat
                       , eModel->fKelvin
                       , mesh->numel, PROP_UPDATE_NL_LOOP);
    if (fDvisc)
      updateDynamicViscosity(&propF->dVisc
                           , mesh->elm.temp, mesh->elm.dViscosity
                           , eModel->fKelvin, mesh->numel);
    if (fTcond)
      updateThermalconductivity(&propF->thCond
                              , mesh->elm.temp, mesh->elm.tConductivity
                              , eModel->fKelvin, mesh->numel);
/*...................................................................*/

/*... pressao de referencia*/
    if (fPresRef)
      presRef(mesh->elm.temp0, mesh->elm.temp
            , mesh->elm.geom.volume, thDynamic->pTh
            , mesh->numel, eModel->fKelvin);
/*...................................................................*/

/*...*/
    if (itSimple == kZeroPres && relResMass &&  pCor) rMass0 = rMass;
    if (itSimple == kZeroEnergy && relResEnergy && fEnergy) rEnergy0 = rEnergy;
    if (itSimple == kZeroVel && relResVel) 
    {
      if (xMomentum) rU0[0] = rU[0];
      if (yMomentum) rU0[1] = rU[1];
      if (zMomentum && ndfVel == 3) rU0[2] = rU[2];
    }
/*...................................................................*/

/*...*/
    timei = getTimeC() - time;
/*... arquivo de log*/
    if (opt->fItPlot) 
    {
      if (ndfVel == 3)
        fprintf(opt->fileItPlot[FITPLOTSIMPLE]
          , "%d %20.8e %20.8e %20.8e %20.8e %20.8e\n"
          , itSimple + 1, rU[0], rU[1], rU[2], rMass, rEnergy);
      else
        fprintf(opt->fileItPlot[FITPLOTSIMPLE]
          , "%d %20.8e %20.8e %20.8e %20.8e\n"
          , itSimple + 1, rU[0], rU[1], rMass, rEnergy);
    }
/*...................................................................*/

/*...*/
    if (jj == sp->pSimple) 
    {
      jj = 0;
      printf("It simple: %d \n", itSimple + 1);
      printf("CPU Time(s)  : %lf \n", timei);
      printf("Residuo:\n");
      printf("conservacao da massa   : %20.8e\n", rMass / rMass0);
      printf("momentum x1            : %20.8e\n", rU[0] / rU0[0]);
      printf("momentum x2            : %20.8e\n", rU[1] / rU0[1]);
      if (ndfVel == 3)
        printf("momentum x3            : %20.8e\n", rU[2] / rU0[2]);
      printf("conservacao da energia : %20.8e\n", rEnergy / rEnergy0);
    }
    jj++;
/*...................................................................*/

/*... 3D*/
    conv = 0;
    if (rEnergy / rEnergy0<tolSimpleEnergy || rEnergy<SZERO) conv++;
    if (rMass / rMass0<tolSimpleMass || rMass<SZERO) conv++;     
    if (rU[0] / rU0[0]<tolSimpleU1 || rU[0]<SZERO) conv++;
    if (rU[1] / rU0[1]<tolSimpleU2 || rU[1]<SZERO) conv++;
    if (rU[2] / rU0[2]<tolSimpleU3 || rU[2]<SZERO) conv++;
    if (conv == 5) break;
/*...................................................................*/

  }
/*...................................................................*/
  timei = getTimeC() - time;

/*...*/
  fParameter[0] = true;
  fParameter[1] = true;
  fParameter[2] = true;
  fParameter[3] = true;
  fParameter[4] = true;
  parameterCellLm(mesh->elm.vel          , mesh->elm.material.prop
                , mesh->elm.densityFluid , mesh->elm.specificHeat
                , mesh->elm.tConductivity, mesh->elm.dViscosity
                , mesh->elm.geom.volume  , mesh->elm.mat
                , &cfl                   , &reynolds
                , &peclet                , &mesh->mass[2]
                , &prMax
                , fParameter             , sc->ddt.dt[0]
                , mesh->numelNov, mesh->ndm);
/*...................................................................*/

/*...*/
  if (fDeltaTimeDynamic)
    dynamicDeltat(mesh->elm.vel, mesh->elm.geom.volume
      , mesh->elm.densityFluid, mesh->elm.specificHeat
      , mesh->elm.tConductivity, mesh->elm.dViscosity
      , sc->ddt.dt, mesh->numelNov
      , mesh->ndm, CFL);
/*...................................................................*/

/*... guardando as propriedades para o proximo passo*/
  if (fDensity) updateDensity(propF
                            , mesh->elm.temp, mesh->elm.pressure0
                            , mesh->elm.densityFluid
                            , sp->alphaDensity, eModel->fKelvin
                            , mesh->numel, PROP_UPDATE_OLD_TIME);
  if (fSheat) updateSpecificHeat(&propF->sHeat
    , mesh->elm.temp, mesh->elm.specificHeat
    , eModel->fKelvin
    , mesh->numel, PROP_UPDATE_OLD_TIME);
/*...................................................................*/

/*... calcula a massa especifica de referencia*/
  if (fDensityRef)
    propF->densityRef = specificMassRef(mesh->elm.densityFluid
                                      , mesh->elm.geom.volume
                                      , mesh->numel);
/*...................................................................*/

/*... calculo da taxa de massa atravessando o contorno aberto*/
  massFluxOpenDomain(loadsVel         , sc->ddt 
              , mesh->elm.cellFace    , mesh->face.owner
              , mesh->elm.faceRvel , mesh->elm.adj.nViz 
              , mesh->face.area       , mesh->face.normal
              , mesh->face.xm  
              , mesh->elm.densityFluid, mesh->elm.vel 
              , mesh->massInOut       , &deltaMass
              , mesh->numelNov        , mesh->ndm  
              , mesh->maxViz  );
  mesh->mass[1] += deltaMass;  
/*...................................................................*/

/*... temp(n) = temp(n+1)*/
  alphaProdVector(1.e0, mesh->elm.temp
    , mesh->numel, mesh->elm.temp0);
/*...................................................................*/

/*... pth(n-2) = pth(n-1)*/
  thDynamic->pTh[0] = thDynamic->pTh[1];
/*... pth(n-1) = pth(n)*/
  thDynamic->pTh[1] = thDynamic->pTh[2];
/*...................................................................*/

/*...*/
  printf("It simple: %d \n", itSimple + 1);
  printf("Time(s) : %lf \n", timei);
  if (sc->ddt.flag)
    printf("%-20s : %13.6lf\n","CFL", cfl);
  printf("%-20s : %13.6lf\n","Reynolds", reynolds);
  printf("%-20s : %13.6lf\n","Peclet"  , peclet);
  printf("%-20s : %13.6lf\n","Prandtl" , prMax);
  printf("%-20s : %13.6e\n","PresRef"  , thDynamic->pTh[2]);
  printf("%-20s : %13.6e\n","(Inc    ) Mass/Mass0", mesh->mass[1]/mesh->mass[0]);
  printf("%-20s : %13.6e\n","(Average) Mass/Mass0", mesh->mass[2]/mesh->mass[0]);
  printf("%-20s : %13.6e\n","MassIn              ", mesh->massInOut[0]);
  printf("%-20s : %13.6e\n","MassOut             ", mesh->massInOut[1]);
  printf("%-25s : %15s %20s\n","Residuo:","init","final");
  printf("%-25s : %20.8e %20.8e\n","conservacao da massa", rMass0, rMass);
  printf("%-25s : %20.8e %20.8e\n","momentum x1",rU0[0], rU[0]);
  printf("%-25s : %20.8e %20.8e\n","momentum x2 ",rU0[1], rU[1]);
  if (ndfVel == 3)
    printf("%-25s : %20.8e %20.8e\n","momentum x3", rU0[2], rU[2]);
  printf("%-25s : %20.8e %20.8e\n","conservacao da energia", rEnergy0, rEnergy);
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 02/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * UPDATECELLSIMPLEVEL : atualizacao dos valores das velocidades     *
 * estimadas com os valores das respectivas equacoes                 *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * w       -> variavel nas celulas                                   * 
 * u1      -> solucao do sistema                                     * 
 * u2      -> solucao do sistema                                     * 
 * id      -> numeracao das equacoes                                 * 
 * numel   -> numero de elementos                                    * 
 * ndm     -> numero de dimensoes                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * w      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void updateCellSimpleVel(DOUBLE  *RESTRICT w
                        ,DOUBLE  *RESTRICT u1 ,DOUBLE  *RESTRICT u2
                        ,INT  *RESTRICT id    ,INT const nEl
                        ,short const ndm)
{
  INT i,lNeq;

  for(i=0;i<nEl;i++){
    lNeq             = id[i] - 1;
    if(lNeq > -1){  
      MAT2D(i,0,w,2) = u1[lNeq];
      MAT2D(i,1,w,2) = u2[lNeq];
    }
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 30/08/2017                                   *
 * Data de modificaco : 22/08/2019                                   * 
 *-------------------------------------------------------------------* 
 * UPDATECELLSIMPLEVELR: atualizacao dos valores das velocidades     *
 * estimadas com os valores das respectivas equacoes                 *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * w       -> variavel nas celulas                                   * 
 * u1      -> solucao do sistema                                     * 
 * u2      -> solucao do sistema                                     *
 * u2      -> solucao do sistema                                     *  
 * id      -> numeracao das equacoes                                 * 
 * iNeq    -> mapa de equacoes de interface                          *
 * numel   -> numero de elementos                                    * 
 * ndm     -> numero de dimensoes                                    * 
 * fCom    -> comunica os valores x entre as particoes               * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * w      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS: Versao residual  v*(n) = v(n) + dv(n+1)                      * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void updateCellSimpleVelR(DOUBLE  *RESTRICT w  ,DOUBLE  *RESTRICT u1
                         ,DOUBLE  *RESTRICT u2 ,DOUBLE  *RESTRICT u3
                         ,INT  *RESTRICT id    ,Interface *iNeq
                         ,INT const nEl        ,short const ndm  
                         ,bool const fRes      ,bool const fCom)
{
  INT i,lNeq;

/*... obtem os valores de x das equacoes em overlaping*/  
  if(fCom)
  {
    comunicateNeq(iNeq,u1);
    comunicateNeq(iNeq,u2);
    if(ndm == 3) 
      comunicateNeq(iNeq,u3);
  }
/*.................................................................*/

/*...*/
  if(fRes)
  {
/*...*/
    if(ndm == 3)
    {    
      for(i=0;i<nEl;i++)
      {
        lNeq             = id[i] - 1;
        if(lNeq > -1)
        {  
          MAT2D(i,0,w,3) += u1[lNeq];
          MAT2D(i,1,w,3) += u2[lNeq];
          MAT2D(i,2,w,3) += u3[lNeq];
        }
      }
    }
/*...................................................................*/

/*...*/
    else
    {
      for(i=0;i<nEl;i++)
      {
        lNeq             = id[i] - 1;
        if(lNeq > -1)
        {  
          MAT2D(i,0,w,2) += u1[lNeq];
          MAT2D(i,1,w,2) += u2[lNeq];
        }
      }
    }
/*...................................................................*/
  }

/*...*/
  else
  {
/*...*/
    if(ndm == 3)
    {    
      for(i=0;i<nEl;i++)
      {
        lNeq             = id[i] - 1;
        if(lNeq > -1){  
          MAT2D(i,0,w,3) = u1[lNeq];
          MAT2D(i,1,w,3) = u2[lNeq];
          MAT2D(i,2,w,3) = u3[lNeq];
        }
      }
    }
/*...................................................................*/

/*...*/
    else
    {
      for(i=0;i<nEl;i++)
      {
        lNeq             = id[i] - 1;
        if(lNeq > -1)
        {  
          MAT2D(i,0,w,2) = u1[lNeq];
          MAT2D(i,1,w,2) = u2[lNeq];
        }
      }
    }
/*...................................................................*/
  }    

}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 11/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * UPDATECELLSIMPLEVEL3D : atualizacao dos valores das velocidades   *
 * estimadas com os valores das respectivas equacoes                 *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * w       -> variavel nas celulas                                   * 
 * u1      -> solucao do sistema                                     * 
 * u2      -> solucao do sistema                                     * 
 * u3      -> solucao do sistema                                     * 
 * id      -> numeracao das equacoes                                 * 
 * numel   -> numero de elementos                                    * 
 * ndm     -> numero de dimensoes                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * w      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void updateCellSimpleVel3D(DOUBLE  *RESTRICT w
                          ,DOUBLE  *RESTRICT u1 
                          ,DOUBLE  *RESTRICT u2
                          ,DOUBLE  *RESTRICT u3
                          ,INT  *RESTRICT id    ,INT const nEl
                          ,short const ndm)
{
  INT i,lNeq;

  for(i=0;i<nEl;i++){
    lNeq             = id[i] - 1;
    if(lNeq > -1){  
      MAT2D(i,0,w,3) = u1[lNeq];
      MAT2D(i,1,w,3) = u2[lNeq];
      MAT2D(i,2,w,3) = u3[lNeq];
    }
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 02/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * UPDATECELLSIMPLEPRES: atualizacao dos valores das pressoes de     *
 * correcao com os valores das respectivas equacoes                  *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * presC   -> variavel nas celulas                                   * 
 * xp      -> solucao do sistema                                     * 
 * id      -> numeracao das equacoes                                 * 
 * iNeq    -> mapa de equacoes de interface                          *
 * numel   -> numero de elementos                                    *
 * fCom    -> comunica os valores x entre as particoes               *  
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * presC  -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void updateCellSimplePres(DOUBLE  *RESTRICT presC,DOUBLE  *RESTRICT xp   
                         ,INT  *RESTRICT id      ,Interface *iNeq
                         ,INT const nEl          ,bool const fCom)
{
  INT i,lNeq;

/*... obtem os valores de x das equacoes em overlaping*/  
  if(fCom)
    comunicateNeq(iNeq,xp);
/*.................................................................*/


  for(i=0;i<nEl;i++){
    lNeq     = id[i] - 1;
    if(lNeq > -1)  
      presC[i] = xp[lNeq];
/*... celulas prescritas*/
    else
      presC[i] = 0.e0;
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 02/07/2016                                   *
 * Data de modificaco : 15/08/2016                                   * 
 *-------------------------------------------------------------------* 
 * SIMPLEUPDATE : atualizacao das variasveis do metodo simple        *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * w         -> velocidades estimadas                                * 
 * pressure  -> pressao desatualizado                                * 
 * presC     -> pressao de correcao                                  * 
 * gradPresC -> gradiente da pressao de correcao                     * 
 * dField    -> matriz D do metodo simple                            * 
 * numel     -> numero de elementos                                  * 
 * ndm       -> numero de dimensoes                                  * 
 * alphaPres -> under-relaxationmensoes                              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * w        -> atualizado                                            * 
 * pressure -> atualizado                                            * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void simpleUpdate(DOUBLE *RESTRICT w     ,DOUBLE *RESTRICT pressure
                 ,DOUBLE *RESTRICT presC ,DOUBLE *RESTRICT gradPresC
                 ,DOUBLE *RESTRICT dField         
                 ,INT const nEl          ,short const ndm
                 ,DOUBLE const alphaPres){
  INT i; 

/*...*/  
  if( ndm == 2)
  {
    for(i=0;i<nEl;i++)
    {
/*... atualizacoes da velocidades*/
      MAT2D(i,0,w,2) -= 
      MAT2D(i,0,dField,2)*MAT2D(i,0,gradPresC,2);
      
      MAT2D(i,1,w,2) -=
      MAT2D(i,1,dField,2)*MAT2D(i,1,gradPresC,2);
/*...................................................................*/

/*... atualizacoes da velocidades*/
      pressure[i] += alphaPres*presC[i];
/*...................................................................*/
    }
  }
/*...................................................................*/

/*...*/  
  else if( ndm == 3)
  {
    for(i=0;i<nEl;i++)
    {
/*... atualizacoes da velocidades*/
      MAT2D(i,0,w,3) -=
      MAT2D(i,0,dField,3)*MAT2D(i,0,gradPresC,3);

      MAT2D(i,1,w,3) -=
      MAT2D(i,1,dField,3)*MAT2D(i,1,gradPresC,3);

      MAT2D(i,2,w,3) -=
      MAT2D(i,2,dField,3)*MAT2D(i,2,gradPresC,3);
  
/*...................................................................*/

/*... atualizacoes da velocidades*/
      pressure[i] += alphaPres*presC[i];
/*...................................................................*/

    }
  }

}               
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 08/05/2019                                   *
* Data de modificaco : 00/00/00000                                  *
*-------------------------------------------------------------------*
* SETSIMPLELMCOMBUSTIONSCHEME : set o metodo simple para baixo      *
* numero mach para modelos de combustao                             *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* word -> str com a discretizacao                                   *
* sp   -> estrutura simple                                          *
* fieIn-> arquivo de entrada                                        *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
void setSimpleCombustionScheme(char *word , short const ndm
                              , Simple *sp, FILE *fileIn) 
{

  if (!strcmp(word, "SIMPLE"))
    sp->type = SIMPLE;
  else if (!strcmp(word, "SIMPLEC"))
    sp->type = SIMPLEC;

  fscanf(fileIn, "%d", &sp->maxIt);

  fscanf(fileIn, "%lf", &sp->alphaPres);
  fscanf(fileIn, "%lf", &sp->alphaVel);
  fscanf(fileIn, "%lf", &sp->alphaEnergy);
  fscanf(fileIn, "%lf", &sp->alphaComb);

  fscanf(fileIn, "%d", &sp->nNonOrth);
  fscanf(fileIn, "%d", &sp->pSimple);


}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 24/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * RESIDUALCOMBUSTION : calculo dos residuos no metodo simple        *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * vel      -  > campo de velocidade                                 *
 * energy     -> campo de energia                                    *  
 * zComb      -> campo com a fracao massica das especies             *
 * rCellVel   -> residuo das equacoes das velocidade por celulas     * 
 * rCellMass  -> residuo de massa da equacoes de pressao por celulas * 
 * rCellEnergy-> residuo de massa da equacoes de energia por celulas * 
 * adVel    -> diagonal da equacoes de velocidades                   *
 * adEnergy -> diagonal da equacoes de energia                       *  
 * rU       -> nao definido                                          * 
 * rMass    -> nao definido                                          * 
 * rEnergy  -> nao definido                                          * 
 * rComb    -> nao definido                                          * 
 * idVel    -> numero da equacao da celula i                         * 
 * idEnergy -> numero da equacao da celula i                         * 
 * idComb   -> numero da equacao da celula i                         * 
 * nEl      -> numero de elementos                                   * 
 * nEqVel   -> numero de equacoes de velocidade                      * 
 * ndm      -> numero de dimensoes                                   * 
 * nComb    -> numero de especies explicitamente resolvidas          * 
 * iCod     -> tipo de residuo                                       * 
 *          RSCALED - residuo com escala de grandeza                 * 
 *          RSQRT   - norma p-2 ( norma euclidiana)                  * 
 *          RSCALEDM- residuo com escala de grandeza                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * rU       -> residuo das velocidades                               * 
 * rMass    -> residuo de mass                                       * 
 * rEnergy  -> residuo de energia                                    * 
 * rComb    -> resido das especies                                   *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *            | rvx1 rvx2 ... rvxnEl |                               *
 * rCellVel = | rvy1 rvy2 ... rvynEl |                               * 
 *            | rvz1 rvz2 ... rvznEl |                               * 
 *                                                                   *
 *         | adx1 adx2 ... adxnEq |                                  *
 * adVel = | ady1 ady2 ... adynEq |                                  *
 *         | adz1 adz2 ... adznEq |                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void residualCombustionOld(DOUBLE *RESTRICT vel ,DOUBLE *RESTRICT energy
            ,DOUBLE *RESTRICT zComb
            ,DOUBLE *RESTRICT rCellVel     ,DOUBLE *RESTRICT rCellMass
            ,DOUBLE *RESTRICT rCellEnergy ,DOUBLE *RESTRICT rCellComb      
            ,DOUBLE *RESTRICT adVel       ,DOUBLE *RESTRICT adEnergy 
            ,DOUBLE *RESTRICT adComb   
            ,DOUBLE *RESTRICT rU          ,DOUBLE *rMass
            ,DOUBLE *rEnergy              ,DOUBLE *rComb    
            ,INT  *RESTRICT idVel         ,INT  *RESTRICT idEnergy
            ,INT  *RESTRICT idComb    
            ,INT const nEl                ,INT const nEqVel
            ,INT const nEqComb
            ,short const ndm              ,short const nComb
            ,short iCod)

{
  DOUBLE maxV[3],sum[MAXSPECIES],maxE,maxComb[MAXSPECIES],mod,tmp,v,rScale,sumMax;
  DOUBLE *p;
  INT i,j,lNeq;
  
/*...*/
  maxV[0] = maxV[1] = maxV[2] = maxE = 0.e0; 
  
  for(j=0;j<nComb;j++)
    maxComb[j] = rComb[j] = sum[j] = 0.e0;

  for(j=0;j<ndm;j++)
    rU[j]   = 0.e0;

  *rEnergy = *rMass = 0.e0;
/*...................................................................*/

/*...*/
  switch(iCod){

/*... scaled*/
    case RSCALED:
/*... max(Ap*velP) */
      for(i=0;i<nEl;i++)
      {
        lNeq = idVel[i] - 1;
        if(lNeq > -1)
        { 
          for(j=0;j<ndm;j++)
          {
            v       = MAT2D(i,j,vel,ndm);
            lNeq   += j*nEqVel;  
            mod     = fabs(adVel[lNeq]*v);
            maxV[j] = max(maxV[j],mod);
          }
        }
      }
/*...................................................................*/
      
/*... max ( | F - Ax |P / max(Ap*velP) )*/
      for(j=0;j<ndm;j++)
      {
        for(i=0;i<nEl;i++)
        {
          mod    = fabs(MAT2D(j,i,rCellVel,nEl));
          if(maxV[j] != 0.e0)
            rScale = mod/maxV[j];
          else
            rScale = 0.e0;
          rU[j]  = max(rU[j],rScale);
        }
      }  
/*...................................................................*/

/*... max(Ap*energyP) */
      for(i=0;i<nEl;i++)
      {
        lNeq = idEnergy[i] - 1;
        if(lNeq > -1)
        { 
          v       = energy[i];
          mod     = fabs(adEnergy[lNeq]*v);
          maxE    = max(maxE,mod);
        }
      }
/*...................................................................*/
      
/*... max ( | F - Ax |P /max(Ap*energyP) )*/
      for(i=0;i<nEl;i++)
      {
        mod    = fabs(MAT2D(j,i,rCellVel,nEl));
        if(maxE != 0.e0)
          rScale = mod/maxE;
        else
          rScale = 0.e0;
        *rEnergy = max(*rEnergy,rScale);
      }  
/*...................................................................*/

/*... max(Ap*zP) */
      for(i=0;i<nEl;i++)
      {
        lNeq = idComb[i] - 1;
        if(lNeq > -1){ 
          for(j=0;j<nComb;j++)
          {
            v       = MAT2D(i,j,zComb,nComb);
            lNeq   += j*nEqComb;  
            mod     = fabs(adComb[lNeq]*v);
            maxComb[j] = max(maxComb[j],mod);
          }
        }
      }
/*...................................................................*/
      
/*... max ( | F - Ax |P / max(Ap*zP) )*/
      for(j=0;j<ndm;j++)
      {
        for(i=0;i<nEl;i++)
        { 
          mod    = fabs(MAT2D(j,i,rCellVel,nEl));
          if(maxComb[j] != 0.e0)
            rScale = mod/maxComb[j];
          else
            rScale = 0.e0;
          rComb[j]  = max(rComb[j],rScale);
        }
      }  
/*...................................................................*/

/*...*/
      tmp = 0.e0;
      for(i=0;i<nEl;i++)
      {
        v    = fabs(rCellMass[i]);
        tmp += v;
      } 
      *rMass = tmp; 
/*...................................................................*/
    break;
/*...................................................................*/

/*... norma euclidiana*/
    case RSQRT:
/*...*/
      for(j=0;j<ndm;j++)
      {
        p     = &rCellVel[j*nEl]; 
        rU[j] = sqrt(dot(p,p,nEl));
      }
/*...................................................................*/

/*...*/
      *rEnergy = sqrt(dot(rCellEnergy,rCellEnergy,nEl));
/*...................................................................*/

/*...*/
      for(j=0;j<nComb;j++)
      {
        p     = &rCellComb[j*nEl]; 
        rComb[j] = sqrt(dot(p,p,nEl));
      }
/*...................................................................*/

/*...*/
      *rMass = sqrt(dot(rCellMass,rCellMass,nEl));
/*...................................................................*/
    break;
/*...................................................................*/

/*... scaled*/
    case RSCALEDSUM:
/*... sum( |Ap*velP |) */
      for(i=0;i<nEl;i++)
      {
        lNeq = idVel[i] - 1;
        if(lNeq > -1)
        { 
          for(j=0;j<ndm;j++)
          {
            v       = MAT2D(i,j,vel,ndm);
            lNeq   += j*nEqVel;
            mod     = fabs(adVel[lNeq]*v);
            sum[j] += mod;
          }
        }
      }
/*...................................................................*/
      
/*... sum ( | F - Ax |P / sum( |Ap*velP| ) )*/
      for(j=0;j<ndm;j++)
      {
        for(i=0;i<nEl;i++)
        {
          mod    = fabs(MAT2D(j,i,rCellVel,nEl));
          rU[j] += mod;
        }
        if( sum[j] > rU[j]*SZERO)
          rU[j]  /=  sum[j];
      }  
/*...................................................................*/

/*... sum( |Ap*energyP| ) */
      sum[0] = 0.e0;
      for(i=0;i<nEl;i++)
      {
        lNeq = idEnergy[i] - 1;
        if(lNeq > -1)
        { 
          v       = energy[i];
          mod     = fabs(adEnergy[lNeq]*v);
          sum[0] += mod;
        }
      }
/*...................................................................*/
      
/*... sum ( | F - Ax |P / sum( |Ap*energyP| ) )*/
      for(i=0;i<nEl;i++)
      {
         mod      = fabs(rCellEnergy[i]);
        *rEnergy += mod;
      }
      if( sum[0] > (*rEnergy)*SZERO)
        *rEnergy /=  sum[0];  
/*...................................................................*/

/*... sum( |Ap*zP |) */
      for(i=0;i<nEl;i++)
      {
        lNeq = idComb[i] - 1;
        if(lNeq > -1)
        { 
          for(j=0;j<nComb;j++)
          {
            v       = MAT2D(i,j,zComb,nComb);
            lNeq   += j*nEqComb;
            mod     = fabs(adComb[lNeq]*v);
            sum[j] += mod;
          }
        }
      }
/*...................................................................*/
      
/*... sum ( | F - Ax |P / sum( |Ap*zP| ) )*/
      for(j=0;j<nComb;j++)
      {
        for(i=0;i<nEl;i++)
        {
          mod    = fabs(MAT2D(j,i,rCellComb,nEl));
          rComb[j] += mod;
        }
        if( sum[j] > rComb[j]*SZERO)
          rComb[j]  /=  sum[j];
      }  
/*...................................................................*/

/*...*/
      tmp = 0.e0;
      for(i=0;i<nEl;i++)
      {
        v    = fabs(rCellMass[i]);
        tmp += v;
      } 
      *rMass = tmp; 
/*...................................................................*/
    break;
/*...................................................................*/

/*... scaled*/
    case RSCALEDSUMMAX:
/*... sum( |Ap*velP |) */
      for(i=0;i<nEl;i++)
      {
        lNeq = idVel[i] - 1;
        if(lNeq > -1)
        { 
          for(j=0;j<ndm;j++)
          {
            v       = MAT2D(i,j,vel,ndm);
            lNeq   += j*nEqVel;
            mod     = fabs(adVel[lNeq]*v);
            sum[j] += mod;
          }
        }
      }
/*...................................................................*/
      
/*...*/
      sumMax = sum[0];
      for (j = 1; j < ndm; j++) 
        sumMax = max(sumMax,sum[j]);
/*...................................................................*/

/*... sum ( | F - Ax |P / sum( |Ap*velP| ) )*/
      for(j=0;j<ndm;j++)
      {
        for(i=0;i<nEl;i++)
        {
          mod    = fabs(MAT2D(j,i,rCellVel,nEl));
          rU[j] += mod;
        }
        rU[j]  /= sumMax;
      }  
/*...................................................................*/

/*... sum( |Ap*energyP| ) */
      sum[0] = 0.e0;
      for(i=0;i<nEl;i++)
      {
        lNeq = idEnergy[i] - 1;
        if(lNeq > -1)
        { 
          v       = energy[i];
          mod     = fabs(adEnergy[lNeq]*v);
          sum[0] += mod;
        }
      }
/*...................................................................*/
      
/*... sum ( | F - Ax |P / sum( |Ap*energyP| ) )*/
      for(i=0;i<nEl;i++)
      {
         mod      = fabs(rCellEnergy[i]);
        *rEnergy += mod;
      }
      if( sum[0] > (*rEnergy)*SZERO)
        *rEnergy /=  sum[0];  
/*...................................................................*/

/*... sum( |Ap*zP |) */
      for(i=0;i<nEl;i++)
      {
        lNeq = idComb[i] - 1;
        if(lNeq > -1)
        { 
          for(j=0;j<nComb;j++)
          {
            v       = MAT2D(i,j,zComb,nComb);
            lNeq   += j*nEqComb;
            mod     = fabs(adComb[lNeq]*v);
            sum[j] += mod;
          }
        }
      }
/*...................................................................*/
      
/*...*/
      sumMax = sum[0];
      for (j = 1; j < nComb; j++) 
        sumMax = max(sumMax,sum[j]);
/*...................................................................*/

/*... sum ( | F - Ax |P / sum( |Ap*zP| ) )*/
      for(j=0;j<nComb;j++)
      {
        for(i=0;i<nEl;i++)
        {
          mod    = fabs(MAT2D(j,i,rCellComb,nEl));
          rComb[j] += mod;
        }
        rComb[j]  /= sumMax;
      }  
/*...................................................................*/

/*...*/
      tmp = 0.e0;
      for(i=0;i<nEl;i++)
      {
        v    = fabs(rCellMass[i]);
        tmp += v;
      } 
      *rMass = tmp; 
/*...................................................................*/
    break;
/*...................................................................*/

/*... */
     default:
       ERRO_OP(__FILE__,__func__,iCod);
     break;
/*...................................................................*/
  }

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 12/11/2017                                   *
 * Data de modificaco : 20/08/2019                                   * 
 *-------------------------------------------------------------------* 
 * dinamicyDeltat : calculo dos residuos no metodo simple            *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * vel      -> campo de velocidade                                   *
 * density   - densidade por celula                                  *
 * sHeat     - calor especifico por celula                           *
 * tCond     - condutividade termica                                 *
 * viscosity -  viscosidade dinamica                                 *
 * volume    - volume                                                *
 * nEl      -> numero de elementos                                   * 
 * dt       -> dt                                                    * 
 * ndm      -> numero de dimensoes                                   * 
 * iCod     -> tipo de residuo                                       * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * ddt.dt   -> dt atualizado                                         *  
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *                                                                   *
 * ddt.dt[0] - delta atual                      (n)                  * 
 * ddt.dt[1] - delta do passo anterior          (n-1)                * 
 * ddt.dt[2] - delta do passo anterior anterior (n-2)                *  
 *-------------------------------------------------------------------* 
 *********************************************************************/
void dynamicDeltat(DOUBLE *RESTRICT vel    , DOUBLE *RESTRICT volume
                 , DOUBLE *RESTRICT density, DOUBLE *RESTRICT sHeat
                 , DOUBLE *RESTRICT tCond  , DOUBLE *RESTRICT dViscosity
                 , DOUBLE *dt              , INT const nEl
                 , short const ndm         , short const iCod)

{
  short nD = DENSITY_LEVEL;
  INT i;
  DOUBLE dtCfl,dtVn,modVel,nCfl,cfl,lc,v[3],deltaT,deltaT0,diff,den;
#ifdef _MPI_
  DOUBLE gg;
#endif  
/*...*/
  switch(iCod)
  {

/*... scaled*/
    case CFL:
      nCfl = 1.5e0;
      cfl  = 0.0e0;
      dtVn = dtCfl = 86400.e0; /*24*3600*/
      for(i=0;i<nEl;i++)
      {
 /*... modulo das velocidades*/
        if (ndm == 2)
        {
          v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
          v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
          modVel = sqrt(v[0] + v[1]);
        }
        else
        {
          v[0] = MAT2D(i, 0, vel, ndm)*MAT2D(i, 0, vel, ndm);
          v[1] = MAT2D(i, 1, vel, ndm)*MAT2D(i, 1, vel, ndm);
          v[2] = MAT2D(i, 2, vel, ndm)*MAT2D(i, 2, vel, ndm);
          modVel = sqrt(v[0] + v[1] + v[2]);
        }
/*...................................................................*/

/*... tamanho caracteristico*/
        lc = sizeCar(volume[i],ndm);
/*...................................................................*/

/*...*/
        diff = tCond[i]/sHeat[i];
        den  = MAT2D(i,2,density,nD);
/*...................................................................*/

/*...*/
        cfl   = max(cfl,lc/modVel);
        dtCfl = min(dtCfl,nCfl*lc/modVel);
        dtVn  = min(dtVn,(den*lc*lc)/(2.0*diff));
/*...................................................................*/
      }
/*...................................................................*/
  
/*....*/
#ifdef _MPI_
  if(mpiVar.nPrcs>1)
  { 
    tm.overHeadMiscMpi = getTimeC() - tm.overHeadMiscMpi;
    MPI_Allreduce(&cfl  ,&gg ,1,MPI_DOUBLE,MPI_MAX,mpiVar.comm);
    cfl = gg;
    MPI_Allreduce(&dtCfl,&gg ,1,MPI_DOUBLE,MPI_MIN,mpiVar.comm);
    dtCfl = gg;
    MPI_Allreduce(&dtVn,&gg ,1,MPI_DOUBLE,MPI_MIN,mpiVar.comm);
    dtVn = gg;
    tm.overHeadMiscMpi = getTimeC() - tm.overHeadMiscMpi;
  }
#endif
/*...................................................................*/

    
/*...*/
      deltaT  = dt[0] ;
      deltaT0 = dt[1];
/*    if (deltaT > dtCfl || deltaT > dtVn ) {
        dt[0] = min(dtCfl,dtVn);
        printf("change dt for: %lf\n",dt[0]);
      }
*/
      if (deltaT > dtCfl) 
      {
        dt[0] = trunkNumber(dtCfl);
        printf("change dt for: %lf\n",dt[0]);
      }  
      dt[1] = deltaT;  
      dt[2] = deltaT0;
/*...................................................................*/
      break;
/*...................................................................*/

/*... */
    default:
       ERRO_OP(__FILE__,__func__,iCod);
    break;
/*...................................................................*/
  }

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 01/04/2018                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * fTrunkNumber: mantem apenas o primeiro digito no numero           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * a        -> numero                                                *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * retorna o numero truncado                                         *  
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *                                                                   *
 * 5.02012  -> 5.0                                                   * 
 * 0.02012  -> 0.02                                                  * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
static DOUBLE trunkNumber(DOUBLE const a){

  double base,aa;
  int firstNumber,j;

  base = 10.0e0;
  for(j=0;j<100;j++){
    aa = a*base;
    firstNumber = (int) aa;
    if(firstNumber) break;
    base*=10.0e0;
  }
  return firstNumber/base;    

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 24/07/2018                                   *
 * Data de modificaco : 21/04/2019                                   *
 *-------------------------------------------------------------------*
 * velPresCouplingLm : Acoplamento pressao-velocidade do metodo      *
 * simple                                                            *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * campos de pressao, velocidade, gradiente de pres, gradiente de    *
 * velocidade e D atualizados                                        * 
 *-------------------------------------------------------------------*
*********************************************************************/
void velPresCoupling(Memoria *m         , PropVarFluid *propF
                    , Loads *loadsVel   , Loads *loadsPres
                    , MassEqModel *eMass, MomentumModel *momentumModel
                    , Turbulence *tModel
                    , Mesh *mesh
                    , SistEq *sistEqVel , SistEq *sistEqPres
                    , Solv *solvVel     , Solv *solvPres
                    , Simple *sp        , Scheme *sc
                    , PartMesh *pMesh   , DOUBLE *rCellPc
                    , bool *xMomentum   , bool *yMomentum
                    , bool *zMomentum   , bool *pCor
                    , bool fPrint       , short itSimple )
{

  short unsigned ndfVel = mesh->ndfFt - 2;
  short nonOrth;
  DOUBLE *rho=NULL;
  DOUBLE *b1, *b2, *b3, *bPc, *xu1, *xu2, *xu3, *xp;
  DOUBLE *adU1, *adU2, *adU3;
/*...*/
  DOUBLE tb[3], tmp;
  tmp = tb[0] = tb[1] = tb[2] = 0.e0;

  ndfVel = max(mesh->ndfF - 1,mesh->ndfFt - 2);

/*...*/
  b1 = b2 = b3 = bPc = NULL;
  xu1 = xu2 = xu3 = NULL;
  adU1 = adU2 = adU3 = NULL;

  b1 = sistEqVel->b;
  b2 = &sistEqVel->b[sistEqVel->neqNov];
  if (ndfVel == 3) b3 = &sistEqVel->b[2 * sistEqVel->neqNov];
  bPc = sistEqPres->b;

  xu1 = sistEqVel->x;
  xu2 = &sistEqVel->x[sistEqVel->neq];
  if (ndfVel == 3) xu3 = &sistEqVel->x[2 * sistEqVel->neq];
  xp = sistEqPres->x;

  adU1 = sistEqVel->ad;
  adU2 = &sistEqVel->ad[sistEqVel->neqNov];
  if (ndfVel == 3) adU3 = &sistEqVel->ad[2 * sistEqVel->neqNov];
/*...................................................................*/

/*... reconstruindo do gradiente da massa especifica*/
  tm.rcGradRho = getTimeC() - tm.rcGradRho;
  HccaAlloc(DOUBLE, m, rho , mesh->numel, "rhot"  , _AD_);
  getColFromMatrix(rho         ,mesh->elm.densityFluid
                  ,mesh->numel ,DENSITY_LEVEL
                  ,TIME_N);
  rcGradU(m                         , loadsRhoFluid
           , mesh->elm.node         , mesh->elm.adj.nelcon
           , mesh->node.x
           , mesh->elm.nen          , mesh->elm.adj.nViz
           , mesh->elm.cellFace     , mesh->face.owner
           , mesh->elm.geom.volume  , mesh->elm.geom.dcca
           , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
           , mesh->face.mksi        , mesh->face.ksi
           , mesh->face.eta         , mesh->face.area
           , mesh->face.normal      , mesh->face.xm
           , mesh->face.mvSkew      , mesh->face.vSkew
           , mesh->elm.geomType     , mesh->elm.material.prop
           , mesh->elm.material.type
           , mesh->elm.mat          , NULL
           , mesh->elm.leastSquare  , mesh->elm.leastSquareR
           , mesh->elm.faceRrho    
           , rho                    , mesh->elm.gradRhoFluid
           , mesh->node.rhoFluid   
           , NULL                   , NULL
           , propF->densityRef
           , &sc->rcGrad
           , mesh->maxNo            , mesh->maxViz
           , 1                      , mesh->ndm
           , &pMesh->iNo            , &pMesh->iEl
           , mesh->numelNov         , mesh->numel
           , mesh->nnodeNov         , mesh->nnode
           , false);  
  HccaDealloc(m, rho,  "rhot", _AD_);
  tm.rcGradRho = getTimeC() - tm.rcGradRho;
/*...................................................................*/

/*...*/
  if (itSimple == 0)
/*... calculo da matrix jacobiana das velocidades
| du1dx1 du1dx2 du1dx3 |
| du2dx1 du2dx2 du2dx3 |
| du3dx1 du3dx2 du3dx3 |
*/
  {
    tm.rcGradVel = getTimeC() - tm.rcGradVel;
    rcGradU(m                       , loadsVel
           , mesh->elm.node         , mesh->elm.adj.nelcon
           , mesh->node.x           
           , mesh->elm.nen          , mesh->elm.adj.nViz
           , mesh->elm.cellFace     , mesh->face.owner
           , mesh->elm.geom.volume  , mesh->elm.geom.dcca
           , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
           , mesh->face.mksi        , mesh->face.ksi
           , mesh->face.eta         , mesh->face.area
           , mesh->face.normal      , mesh->face.xm
           , mesh->face.mvSkew      , mesh->face.vSkew
           , mesh->elm.geomType     , mesh->elm.material.prop
           , mesh->elm.material.type
           , mesh->elm.mat          , NULL
           , mesh->elm.leastSquare  , mesh->elm.leastSquareR
           , mesh->elm.faceRvel     
           , mesh->elm.vel          , mesh->elm.gradVel
           , mesh->node.vel         
           , NULL                   , NULL
           , 0
           , &sc->rcGrad
           , mesh->maxNo            , mesh->maxViz
           , ndfVel                 , mesh->ndm
           , &pMesh->iNo            , &pMesh->iEl
           , mesh->numelNov         , mesh->numel
           , mesh->nnodeNov         , mesh->nnode
           , false);  
    tm.rcGradVel = getTimeC() - tm.rcGradVel;
/*...................................................................*/

/*... reconstruindo do gradiente da pressao*/
    tm.rcGradPres = getTimeC() - tm.rcGradPres;
    rcGradU(m                       , loadsPres
           , mesh->elm.node         , mesh->elm.adj.nelcon
           , mesh->node.x
           , mesh->elm.nen          , mesh->elm.adj.nViz
           , mesh->elm.cellFace     , mesh->face.owner
           , mesh->elm.geom.volume  , mesh->elm.geom.dcca
           , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
           , mesh->face.mksi        , mesh->face.ksi
           , mesh->face.eta         , mesh->face.area
           , mesh->face.normal      , mesh->face.xm
           , mesh->face.mvSkew      , mesh->face.vSkew
           , mesh->elm.geomType     , mesh->elm.material.prop
           , mesh->elm.material.type
           , mesh->elm.mat          , NULL
           , mesh->elm.leastSquare  , mesh->elm.leastSquareR
           , mesh->elm.faceRpres    
           , mesh->elm.pressure     , mesh->elm.gradPres
           , mesh->node.pressure    
           , mesh->elm.densityFluid , mesh->elm.gradRhoFluid 
           , propF->densityRef
           , &sc->rcGrad
           , mesh->maxNo            , mesh->maxViz
           , 1                      , mesh->ndm
           , &pMesh->iNo            , &pMesh->iEl
           , mesh->numelNov         , mesh->numel
           , mesh->nnodeNov         , mesh->nnode
           , false);  
    tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/
  }
/*...................................................................*/

/*... montagem do sistema u, v e w*/
  tm.systFormVel = getTimeC() - tm.systFormVel;
  if(propF == NULL)
    systFormSimpleVel(loadsVel           , loadsPres
               , &sc->advVel             , &sc->diffVel
               , tModel                  , momentumModel
               , &pMesh->iEl             , sp->type          
               , mesh->elm.node          , mesh->elm.adj.nelcon
               , mesh->elm.nen           , mesh->elm.adj.nViz
               , mesh->elm.cellFace      , mesh->face.owner
               , mesh->elm.geom.volume   , mesh->elm.geom.dcca
               , mesh->elm.geom.xmcc     , mesh->elm.geom.cc
               , mesh->face.mksi         , mesh->face.ksi
               , mesh->face.eta          , mesh->face.area
               , mesh->face.normal       , mesh->face.xm
               , mesh->face.mvSkew       , mesh->face.vSkew
               , mesh->elm.geomType      , mesh->elm.material.prop
               , mesh->elm.material.type , mesh->elm.mat
               , sistEqVel->ia           , sistEqVel->ja
               , sistEqVel->al           , sistEqVel->ad
               , sistEqVel->b            , sistEqVel->id
               , mesh->elm.faceRvel      , mesh->elm.faceRpres
               , mesh->elm.pressure      , mesh->elm.gradPres
               , mesh->elm.vel           , mesh->elm.gradVel
               , sp->d                   , sp->alphaVel
               , mesh->elm.rCellVel      , mesh->elm.stressR
               , mesh->elm.eddyViscosity , mesh->elm.wallParameters
               , &sc->ddt
               , sistEqVel->neq          , sistEqVel->neqNov
               , sistEqVel->nad          , sistEqVel->nadr
               , mesh->maxNo             , mesh->maxViz
               , mesh->ndm               , mesh->numelNov
               , ndfVel                  , mesh->ntn
               , sistEqVel->storage
               , true                    , true
               , true                    , sistEqVel->unsym);
  else
    systFormSimpleVelLm(loadsVel, loadsPres
                      , &sc->advVel             , &sc->diffVel
                      , tModel                  , momentumModel
                      , &pMesh->iEl             , sp->type
                      , mesh->elm.node          , mesh->elm.adj.nelcon
                      , mesh->elm.nen           , mesh->elm.adj.nViz
                      , mesh->elm.cellFace      , mesh->face.owner
                      , mesh->elm.geom.volume   , mesh->elm.geom.dcca
                      , mesh->elm.geom.xmcc     , mesh->elm.geom.cc
                      , mesh->face.mksi         , mesh->face.ksi
                      , mesh->face.eta          , mesh->face.area
                      , mesh->face.normal       , mesh->face.xm
                      , mesh->face.mvSkew       , mesh->face.vSkew
                      , mesh->elm.geomType      
                      , mesh->elm.material.type , mesh->elm.mat
                      , sistEqVel->ia           , sistEqVel->ja
                      , sistEqVel->al           , sistEqVel->ad
                      , sistEqVel->b            , sistEqVel->id
                      , mesh->elm.faceRvel      , mesh->elm.faceRpres
                      , mesh->elm.pressure      , mesh->elm.gradPres
                      , mesh->elm.vel           , mesh->elm.gradVel
                      , mesh->elm.gradRhoFluid
                      , sp->d                   , sp->alphaVel
                      , mesh->elm.rCellVel      , mesh->elm.stressR
                      , mesh->elm.densityFluid  , mesh->elm.dViscosity
                      , mesh->elm.eddyViscosity , mesh->elm.wallParameters
                      , propF->densityRef       , &sc->ddt
                      , sistEqVel->neq          , sistEqVel->neqNov
                      , sistEqVel->nad          , sistEqVel->nadr
                      , mesh->maxNo             , mesh->maxViz
                      , mesh->ndm               , mesh->numelNov
                      , ndfVel                  , sistEqVel->storage
                      , mesh->ntn               , true
                      , true                    , true
                      , sistEqVel->unsym        , sp->sPressure);
  tm.systFormVel = getTimeC() - tm.systFormVel;
/*...................................................................*/ 

/*... soma o vetor b(i) = b(i) + b0(i)*/
  addVector(1.0e0                   , sistEqVel->b
          , 1.0e0                   , sistEqVel->b0
          , sistEqVel->neqNov*ndfVel, sistEqVel->b);
/*...................................................................*/

/*... soma o vetor R(i) = R(i) + b0(i)*/
  updateCellValueSimple(mesh->elm.rCellVel, sistEqVel->b0
                      , sistEqVel->id     , &sistEqVel->iNeq
                      , mesh->numelNov    , sistEqVel->neqNov
                      , ndfVel
                      , true              , false);
/*...................................................................*/

/*...*/
  tb[0] = sqrt(dot(b1, b1, sistEqVel->neqNov));
  tb[1] = sqrt(dot(b2, b2, sistEqVel->neqNov));
  if (ndfVel == 3) tb[2] = sqrt(dot(b3, b3, sistEqVel->neqNov));
  if (itSimple == 0)
  {
    tmp = max(tb[0], tb[1]);
    if (ndfVel == 3) tmp = max(tmp, tb[2]);
  }
  /*...*/
  *xMomentum = true;
  if (tb[0] <= tmp * SZERO || tb[0] == 0.e0) *xMomentum = false;
  /*...................................................................*/

/*...*/
  *yMomentum = true;
  if (tb[1] <= tmp * SZERO || tb[1] == 0.e0) *yMomentum = false;
/*...................................................................*/


/*...*/
  if (ndfVel == 3) 
  {
    *zMomentum = true;
    if (tb[2] <= tmp * SZERO || tb[2] == 0.e0) *zMomentum = false;
  }
/*...................................................................*/

/*... solver Au = bu (velocidade estimadas)*/
  if (*xMomentum) {
    if (fPrint) printf("Quantidade de movimento u1:\n");
    tm.solvVel = getTimeC() - tm.solvVel;
    solverC(m
      , sistEqVel->neq    , sistEqVel->neqNov
      , sistEqVel->nad    , sistEqVel->nadr
      , sistEqVel->ia     , sistEqVel->ja
      , sistEqVel->al     , adU1, sistEqVel->au
      , b1                , xu1
      , &sistEqVel->iNeq  , &sistEqVel->omp
      , solvVel->tol      , solvVel->maxIt
      , sistEqVel->storage, solvVel->solver
      , solvVel->fileSolv , solvVel->log
      , true              , sistEqVel->unsym);
    tm.solvVel = getTimeC() - tm.solvVel;
  }
/*...................................................................*/

/*... solver Av = bv (velocidade estimadas)*/
  if (*yMomentum) {
    if (fPrint) printf("Quantidade de movimento u2:\n");
    tm.solvVel = getTimeC() - tm.solvVel;
    solverC(m
      , sistEqVel->neq    , sistEqVel->neqNov
      , sistEqVel->nad    , sistEqVel->nadr
      , sistEqVel->ia     , sistEqVel->ja
      , sistEqVel->al     , adU2, sistEqVel->au
      , b2                , xu2
      , &sistEqVel->iNeq  , &sistEqVel->omp
      , solvVel->tol      , solvVel->maxIt
      , sistEqVel->storage, solvVel->solver
      , solvVel->fileSolv , solvVel->log
      , true              , sistEqVel->unsym);
    tm.solvVel = getTimeC() - tm.solvVel;
  }
/*...................................................................*/

/*... solver Aw = bw (velocidade estimadas)*/
  if (*zMomentum && ndfVel == 3) {
    if (fPrint) printf("Quantidade de movimento u3:\n");
    tm.solvVel = getTimeC() - tm.solvVel;
    solverC(m
      , sistEqVel->neq    , sistEqVel->neqNov
      , sistEqVel->nad    , sistEqVel->nadr
      , sistEqVel->ia     , sistEqVel->ja
      , sistEqVel->al     , adU3, sistEqVel->au
      , b3                , xu3
      , &sistEqVel->iNeq  , &sistEqVel->omp
      , solvVel->tol      , solvVel->maxIt
      , sistEqVel->storage, solvVel->solver
      , solvVel->fileSolv , solvVel->log
      , true              , sistEqVel->unsym);
    tm.solvVel = getTimeC() - tm.solvVel;
  }
/*...................................................................*/

/*... atualizando o campo de velociade estimadas*/
  updateCellSimpleVelR(mesh->elm.vel      , xu1
                     , xu2                , xu3
                     , sistEqVel->id      , &sistEqVel->iNeq 
                     , mesh->numel        , mesh->ndm
                     , momentumModel->fRes, true);                     
/*...................................................................*/

/*...*/
  if (fPrint) printf("Correcao de pressao:\n");

/*... montagem do sistema  da pressao de correca*/
  tm.systFormPres = getTimeC() - tm.systFormPres;
  if(propF == NULL)
    systFormSimplePres(loadsVel            , loadsPresC 
							 , &sc->diffPres             , tModel
               , mesh->elm.node            , mesh->elm.adj.nelcon 
               , mesh->elm.nen             , mesh->elm.adj.nViz
               , mesh->elm.cellFace        , mesh->face.owner
               , mesh->elm.geom.volume     , mesh->elm.geom.dcca
               , mesh->elm.geom.xmcc       , mesh->elm.geom.cc
               , mesh->face.mksi           , mesh->face.ksi
               , mesh->face.eta            , mesh->face.area
               , mesh->face.normal         , mesh->face.xm
               , mesh->face.mvSkew         , mesh->face.vSkew
               , mesh->elm.geomType        , mesh->elm.material.prop
               , mesh->elm.material.type   , mesh->elm.mat
               , sistEqPres->ia            , sistEqPres->ja
               , sistEqPres->al            , sistEqPres->ad
               , bPc                       , sistEqPres->id
               , mesh->elm.faceRvel        , mesh->elm.faceRpres     
               , mesh->elm.pressure        , mesh->elm.gradPres
               , mesh->elm.vel             , sp->d   
               , mesh->elm.wallParameters  , rCellPc
               , &sc->ddt              
               , sistEqPres->neq           , sistEqPres->neqNov
               , sistEqPres->nad           , sistEqPres->nadr                 
               , mesh->maxNo               , mesh->maxViz
               , mesh->ndm                 , mesh->numelNov
               , ndfVel                    , sistEqPres->storage
               , true                      , true
               , true                      , sistEqPres->unsym); 
  else
    systFormSimplePresLm(loadsVel, loadsPresC
            , &sc->diffPres         
            , eMass                  , momentumModel 
            , tModel
            , mesh->elm.node         , mesh->elm.adj.nelcon
            , mesh->elm.nen          , mesh->elm.adj.nViz
            , mesh->elm.cellFace     , mesh->face.owner
            , mesh->elm.geom.volume  , mesh->elm.geom.dcca
            , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
            , mesh->face.mksi        , mesh->face.ksi
            , mesh->face.eta         , mesh->face.area
            , mesh->face.normal      , mesh->face.xm
            , mesh->face.mvSkew      , mesh->face.vSkew
            , mesh->elm.geomType     
            , mesh->elm.material.type, mesh->elm.mat
            , sistEqPres->ia         , sistEqPres->ja
            , sistEqPres->al         , sistEqPres->ad
            , bPc                    , sistEqPres->id
            , mesh->elm.faceRvel     , mesh->elm.faceRpres
            , mesh->elm.pressure     , mesh->elm.gradPres
            , mesh->elm.gradRhoFluid 
            , mesh->elm.vel          , sp->d
            , mesh->elm.temp         , mesh->elm.wallParameters
            , rCellPc                , mesh->elm.densityFluid
            , propF->densityRef      , &sc->ddt
            , sistEqPres->neq        , sistEqPres->neqNov
            , sistEqPres->nad        , sistEqPres->nadr
            , mesh->maxNo            , mesh->maxViz
            , mesh->ndm              , mesh->numelNov
            , ndfVel                 , sistEqPres->storage
            , true                   , true
            , true                   , sistEqPres->unsym);
  tm.systFormPres = getTimeC() - tm.systFormPres;
/*...................................................................*/

/*... residual*/
  tmp = sqrt(dot(bPc, bPc, sistEqPres->neqNov));
/*...................................................................*/

/*...*/
  *pCor = true;
  if (tmp <= SZERO || tmp == 0.e0) *pCor = false;
/*...................................................................*/

/*... solver ApPc = bpC (velocidade estimadas)*/
  if (*pCor) {
/*...*/
    zero(sp->ePresC , mesh->numel, DOUBLEC);
    zero(sp->ePresC1, mesh->numel, DOUBLEC);
/*...................................................................*/

/*...*/
    tm.solvPres = getTimeC() - tm.solvPres;
    solverC(m
      , sistEqPres->neq, sistEqPres->neqNov
      , sistEqPres->nad, sistEqPres->nadr
      , sistEqPres->ia, sistEqPres->ja
      , sistEqPres->al, sistEqPres->ad, sistEqPres->au
      , bPc, xp
      , &sistEqPres->iNeq, &sistEqPres->omp
      , solvPres->tol, solvPres->maxIt
      , sistEqPres->storage, solvPres->solver
      , solvPres->fileSolv, solvPres->log
      , true, sistEqPres->unsym);
    tm.solvPres = getTimeC() - tm.solvPres;
/*...................................................................*/

/*... atualizando da pressao de correcao*/
    updateCellSimplePres(sp->ePresC    , xp
                       , sistEqPres->id, &sistEqPres->iNeq
                       , mesh->numel   ,true);
/*...................................................................*/

/*...*/
    alphaProdVector(1.e0, sp->ePresC, mesh->numel, sp->ePresC1);
/*...................................................................*/

/*... correcao nao ortoganal da pressao de correcao*/
    for (nonOrth = 0; nonOrth < sp->nNonOrth; nonOrth++) {
/*... reconstruindo do gradiente da pressao correcao*/
      tm.rcGradPres = getTimeC() - tm.rcGradPres;
      rcGradU(m                      , loadsPresC
            , mesh->elm.node         , mesh->elm.adj.nelcon
            , mesh->node.x           
            , mesh->elm.nen          , mesh->elm.adj.nViz
            , mesh->elm.cellFace     , mesh->face.owner
            , mesh->elm.geom.volume  , mesh->elm.geom.dcca
            , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
            , mesh->face.mksi        , mesh->face.ksi
            , mesh->face.eta         , mesh->face.area
            , mesh->face.normal      , mesh->face.xm
            , mesh->face.mvSkew      , mesh->face.vSkew
            , mesh->elm.geomType     , mesh->elm.material.prop
            , mesh->elm.material.type
            , mesh->elm.mat          , NULL
            , mesh->elm.leastSquare  , mesh->elm.leastSquareR
            , mesh->elm.faceRpres    
            , sp->ePresC1            , sp->eGradPresC
            , sp->nPresC             
            , mesh->elm.densityFluid , mesh->elm.gradRhoFluid 
            , propF->densityRef
            , &sc->rcGrad
            , mesh->maxNo            , mesh->maxViz
            , 1                      , mesh->ndm
            , &pMesh->iNo            , &pMesh->iEl
            , mesh->numelNov         , mesh->numel
            , mesh->nnodeNov         , mesh->nnode
            , false);
      tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*...*/
      tm.systFormPres = getTimeC() - tm.systFormPres;
      simpleNonOrthPres(&sc->diffPres
        , mesh->elm.node, mesh->elm.adj.nelcon
        , mesh->elm.nen, mesh->elm.adj.nViz
        , mesh->elm.cellFace, mesh->face.owner
        , mesh->elm.geom.volume, mesh->elm.geom.dcca
        , mesh->elm.geom.xmcc, mesh->elm.geom.cc
        , mesh->face.mksi, mesh->face.ksi
        , mesh->face.eta, mesh->face.area
        , mesh->face.normal, mesh->face.xm
        , mesh->face.mvSkew, mesh->face.vSkew
        , mesh->elm.geomType, mesh->elm.material.prop
        , mesh->elm.material.type, mesh->elm.mat
        , mesh->elm.densityFluid
        , bPc, sistEqPres->id
        , mesh->elm.faceRpres, sp->ePresC1
        , sp->eGradPresC, sp->d
        , mesh->maxNo, mesh->maxViz
        , mesh->ndm, mesh->numelNov);
      tm.systFormPres = getTimeC() - tm.systFormPres;
/*...................................................................*/

/*...*/
      tm.solvPres = getTimeC() - tm.solvPres;
      solverC(m
        , sistEqPres->neq, sistEqPres->neqNov
        , sistEqPres->nad, sistEqPres->nadr
        , sistEqPres->ia, sistEqPres->ja
        , sistEqPres->al, sistEqPres->ad, sistEqPres->au
        , bPc, xp
        , &sistEqPres->iNeq, &sistEqPres->omp
        , solvPres->tol, solvPres->maxIt
        , sistEqPres->storage, solvPres->solver
        , solvPres->fileSolv, solvPres->log
        , true, sistEqPres->unsym);
      tm.solvPres = getTimeC() - tm.solvPres;
/*...................................................................*/

/*... atualizando da pressao de correcao*/
      updateCellSimplePres(sp->ePresC1 , xp
                       , sistEqPres->id, &sistEqPres->iNeq
                       , mesh->numel   ,true);
/*...................................................................*/

/*... soma o vetor presC(i) = presC + presC1*/
      addVector(1.0e0       , sp->ePresC
               , 1.0e0      , sp->ePresC1
               , mesh->numel, sp->ePresC);
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... reconstruindo do gradiente da pressao correcao*/
  tm.rcGradPres = getTimeC() - tm.rcGradPres;
  rcGradU(m                      , loadsPresC
        , mesh->elm.node         , mesh->elm.adj.nelcon
        , mesh->node.x           
        , mesh->elm.nen          , mesh->elm.adj.nViz
        , mesh->elm.cellFace     , mesh->face.owner
        , mesh->elm.geom.volume  , mesh->elm.geom.dcca
        , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
        , mesh->face.mksi        , mesh->face.ksi
        , mesh->face.eta         , mesh->face.area
        , mesh->face.normal      , mesh->face.xm
        , mesh->face.mvSkew      , mesh->face.vSkew
        , mesh->elm.geomType     , mesh->elm.material.prop
        , mesh->elm.material.type
        , mesh->elm.mat          , NULL
        , mesh->elm.leastSquare  , mesh->elm.leastSquareR
        , mesh->elm.faceRpres    
        , sp->ePresC             , sp->eGradPresC
        , sp->nPresC             
        , mesh->elm.densityFluid , mesh->elm.gradRhoFluid 
        , propF->densityRef
        , &sc->rcGrad
        , mesh->maxNo            , mesh->maxViz
        , 1                      , mesh->ndm
        , &pMesh->iNo            , &pMesh->iEl
        , mesh->numelNov         , mesh->numel
        , mesh->nnodeNov         , mesh->nnode
        , false);  
  tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*... atualizacao de u, v, w e p*/
  simpleUpdate(mesh->elm.vel, mesh->elm.pressure
             , sp->ePresC   , sp->eGradPresC
             , sp->d
             , mesh->numel  , mesh->ndm
             , sp->alphaPres);
/*...................................................................*/

/*... calculo da matrix jacobiana das velocidades
  | du1dx1 du1dx2 du1dx3 |
  | du2dx1 du2dx2 du2dx3 |
  | du3dx1 du3dx2 du3dx3 |
*/
  tm.rcGradVel = getTimeC() - tm.rcGradVel;
  rcGradU(m                      , loadsVel
        , mesh->elm.node         , mesh->elm.adj.nelcon
        , mesh->node.x           
        , mesh->elm.nen          , mesh->elm.adj.nViz
        , mesh->elm.cellFace     , mesh->face.owner
        , mesh->elm.geom.volume  , mesh->elm.geom.dcca
        , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
        , mesh->face.mksi        , mesh->face.ksi
        , mesh->face.eta         , mesh->face.area
        , mesh->face.normal      , mesh->face.xm
        , mesh->face.mvSkew      , mesh->face.vSkew
        , mesh->elm.geomType     , mesh->elm.material.prop
        , mesh->elm.material.type
        , mesh->elm.mat          , NULL
        , mesh->elm.leastSquare  , mesh->elm.leastSquareR
        , mesh->elm.faceRvel     
        , mesh->elm.vel          , mesh->elm.gradVel
        , mesh->node.vel         
        , NULL                   , NULL
        , 0
        , &sc->rcGrad
        , mesh->maxNo            , mesh->maxViz
        , ndfVel                 , mesh->ndm
        , &pMesh->iNo            , &pMesh->iEl
        , mesh->numelNov         , mesh->numel
        , mesh->nnodeNov         , mesh->nnode
        , false);  
  tm.rcGradVel = getTimeC() - tm.rcGradVel;
/*...................................................................*/

/*... reconstruindo do gradiente da pressao*/
  tm.rcGradPres = getTimeC() - tm.rcGradPres;
  rcGradU(m                      , loadsPres
        , mesh->elm.node         , mesh->elm.adj.nelcon
        , mesh->node.x
        , mesh->elm.nen          , mesh->elm.adj.nViz
        , mesh->elm.cellFace     , mesh->face.owner
        , mesh->elm.geom.volume  , mesh->elm.geom.dcca
        , mesh->elm.geom.xmcc    , mesh->elm.geom.cc
        , mesh->face.mksi        , mesh->face.ksi
        , mesh->face.eta         , mesh->face.area
        , mesh->face.normal      , mesh->face.xm
        , mesh->face.mvSkew      , mesh->face.vSkew
        , mesh->elm.geomType     , mesh->elm.material.prop
        , mesh->elm.material.type
        , mesh->elm.mat          , NULL
        , mesh->elm.leastSquare  , mesh->elm.leastSquareR
        , mesh->elm.faceRpres    
        , mesh->elm.pressure     , mesh->elm.gradPres
        , mesh->node.pressure    
        , mesh->elm.densityFluid , mesh->elm.gradRhoFluid 
        , propF->densityRef
        , &sc->rcGrad
        , mesh->maxNo            , mesh->maxViz
        , 1                      , mesh->ndm
        , &pMesh->iNo            , &pMesh->iEl
        , mesh->numelNov         , mesh->numel
        , mesh->nnodeNov         , mesh->nnode
        , false);  
  tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

}
/*********************************************************************/


