#include<Print.h>
/********************************************************************* 
 * Data de criacao    : 02/12/2017                                   *
 * Data de modificaco : 30/01/2018                                   *
 *-------------------------------------------------------------------*
 * printFluid: impressao do fluido                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*  
 * m       -> vetor de memoria principal                             * 
 * turbModel -> modelo de turbulencia                                *
 * eModel    -> modelo da equacao de energia                         *
 * loadVel   -> deficicao de cargas velocidade                       *
 * loadPres  -> deficicao de cargas pres                             *
 * loadTemp  -> deficicao de cargas temperatura                      *
 * opt       -> opcoes de arquivo                                    *
 * mesh0     -> malha global                                         *
 * mesh      -> malha particionada                                   *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*  
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void printFluid(Memoria *m        
               ,Turbulence *turbModel,EnergyModel *eModel
               ,PartMesh *pMesh      ,Scheme sc
               ,Loads *loadsVel      ,Loads *loadsPres 
               ,Loads *loadsTemp     ,FileOpt opt
               ,Mesh *mesh0          ,Mesh *mesh  
               ,Mean *media      
               ,char *preName        ,char *nameOut){
 
  short ndm = mesh->ndm;
  void *dum=NULL;
  INT ndfVel;
  DOUBLE *nStressR=NULL,*nEddyV=NULL,*nDvisc=NULL,*nDenFluid=NULL;
  DOUBLE *nMedVel=NULL,*nP2Vel=NULL,*nMedP2Vel=NULL;
  DOUBLE *nCdyn=NULL,*nWall=NULL,*nKturb=NULL;
  FILE *fileOut=NULL;

/*...*/
  ndfVel = max(mesh->ndfF - 1, mesh->ndfFt - 2);
/*...................................................................*/

/*...*/
  HccaAlloc(DOUBLE, m, nDenFluid, mesh->nnode*3  , "nDenFluid", _AD_);
  HccaAlloc(DOUBLE, m, nDvisc   , mesh->nnode    , "nVis"     , _AD_); 
  HccaAlloc(DOUBLE, m, nEddyV   , mesh->nnode    , "nEddyV"   , _AD_);
  HccaAlloc(DOUBLE, m, nStressR , mesh->nnode*6  , "nStressR" , _AD_);
  HccaAlloc(DOUBLE, m, nCdyn    , mesh->nnode*2  , "nCdyn"    , _AD_); 
  HccaAlloc(DOUBLE, m, nWall    , mesh->nnode*4  , "nWall"    , _AD_); 
  HccaAlloc(DOUBLE, m, nKturb   , mesh->nnode    , "nKturb"   , _AD_);
/*...................................................................*/

/*...*/
  if (media->fMedia) 
  {
    HccaAlloc(DOUBLE, m, nMedVel  , mesh->nnode*ndm, "nMediaVel", _AD_);
    zero(nMedVel, mesh->nnode*ndm,DOUBLEC);
    HccaAlloc(DOUBLE, m, nP2Vel   , mesh->nnode*ndm, "nP2Vel", _AD_);
    zero(nMedVel, mesh->nnode*ndm,DOUBLEC);
    HccaAlloc(DOUBLE, m,nMedP2Vel, mesh->nnode*ndm, "nMedP2Vel", _AD_);
    zero(nMedVel, mesh->nnode*ndm,DOUBLEC);
  }
/*...................................................................*/

/*... reconstruindo do gradiente (Pres)*/
  tm.rcGradPres = getTimeC() - tm.rcGradPres;
  rcGradU(m                        , loadsPres
         , mesh->elm.node          , mesh->elm.adj.nelcon
         , mesh->node.x
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
         , mesh->elm.leastSquare   , mesh->elm.leastSquareR
         , mesh->elm.faceRpres     , mesh->elm.faceLoadPres
         , mesh->elm.pressure      , mesh->elm.gradPres                
         , mesh->node.pressure     , sc.rcGrad
         , mesh->maxNo             , mesh->maxViz
         , 1                       , mesh->ndm       
         , &pMesh->iNo             , &pMesh->iEl  
         , mesh->numelNov          , mesh->numel        
         , mesh->nnodeNov          , mesh->nnode);  
  tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (GradPres)*/
  interCellNode(m                  , loadsPres
              , mesh->elm.cellFace , mesh->face.owner
              , mesh->node.gradPres, mesh->elm.gradPres
              , mesh->elm.node     , mesh->elm.geomType            
              , mesh->elm.geom.cc  , mesh->node.x  
              , mesh->face.xm
              , mesh->elm.nen      , mesh->elm.adj.nViz
              , mesh->elm.faceRpres, mesh->elm.faceLoadPres
              , &pMesh->iNo          
              , mesh->numelNov     , mesh->numel        
              , mesh->nnodeNov     , mesh->nnode 
              , mesh->maxNo        , mesh->maxViz   
              , mesh->ndm          , 1
              , mesh->ndm      
              , false              , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (pres)*/
  interCellNode(m                   , loadsPres
               , mesh->elm.cellFace , mesh->face.owner
               , mesh->node.pressure, mesh->elm.pressure   
               , mesh->elm.node     , mesh->elm.geomType            
               , mesh->elm.geom.cc  , mesh->node.x  
               , mesh->face.xm        
               , mesh->elm.nen      , mesh->elm.adj.nViz
               , mesh->elm.faceRpres, mesh->elm.faceLoadPres  
               , &pMesh->iNo           
               , mesh->numelNov     , mesh->numel        
               , mesh->nnodeNov     , mesh->nnode 
               , mesh->maxNo        , mesh->maxViz   
               , 1                  , 1
               , mesh->ndm            
               , opt.bconditions    , 2);
/*...................................................................*/

/*... calculo da matrix jacobiana das velocidades
                            | du1dx1 du1dx2 du1dx3 |   
                            | du2dx1 du2dx2 du2dx3 |   
                            | du3dx1 du3dx2 du3dx3 |
*/        
  tm.rcGradVel  = getTimeC() - tm.rcGradVel;
  rcGradU(m                        , loadsVel
         , mesh->elm.node          , mesh->elm.adj.nelcon
         , mesh->node.x
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
         , mesh->elm.leastSquare   , mesh->elm.leastSquareR
         , mesh->elm.faceRvel      , mesh->elm.faceLoadVel   
         , mesh->elm.vel           , mesh->elm.gradVel                           
         , mesh->node.vel          , sc.rcGrad
         , mesh->maxNo             , mesh->maxViz
         , ndfVel                  , mesh->ndm
         , &pMesh->iNo             , &pMesh->iEl  
         , mesh->numelNov          , mesh->numel        
         , mesh->nnodeNov          , mesh->nnode); 
  tm.rcGradVel = getTimeC() - tm.rcGradVel;  
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (gradVel)*/
  interCellNode(m                   , loadsVel
               , mesh->elm.cellFace , mesh->face.owner
               , mesh->node.gradVel , mesh->elm.gradVel    
               , mesh->elm.node     , mesh->elm.geomType            
               , mesh->elm.geom.cc  , mesh->node.x  
               , mesh->face.xm        
               , mesh->elm.nen      , mesh->elm.adj.nViz
               , mesh->elm.faceRvel , mesh->elm.faceLoadVel 
               , &pMesh->iNo           
               , mesh->numelNov     , mesh->numel        
               , mesh->nnodeNov     , mesh->nnode 
               , mesh->maxNo        , mesh->maxViz   
               , ndfVel             , mesh->ndm
               , mesh->ndm            
               , false              , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (vel)*/
  interCellNode(m                  , loadsVel
               , mesh->elm.cellFace, mesh->face.owner
               , mesh->node.vel    , mesh->elm.vel        
               , mesh->elm.node    , mesh->elm.geomType            
               , mesh->elm.geom.cc , mesh->node.x  
               , mesh->face.xm   
               , mesh->elm.nen     , mesh->elm.adj.nViz
               , mesh->elm.faceRvel, mesh->elm.faceLoadVel 
               , &pMesh->iNo          
               , mesh->numelNov    , mesh->numel        
               , mesh->nnodeNov    , mesh->nnode 
               , mesh->maxNo       , mesh->maxViz   
               , mesh->ndm         , 1
               , mesh->ndm      
               , opt.bconditions   , 2);
/*...................................................................*/

/*... medias*/
  if (media->fMedia) {

/*... medias das velociade*/
    if(media->fVel)
      interCellNode(m                   , loadsVel
                   , mesh->elm.cellFace , mesh->face.owner
                   , nMedVel            , media->mVel        
                   , mesh->elm.node     , mesh->elm.geomType            
                   , mesh->elm.geom.cc  , mesh->node.x  
                   , mesh->face.xm        
                   , mesh->elm.nen      , mesh->elm.adj.nViz
                   , mesh->elm.faceRvel , mesh->elm.faceLoadVel 
                   , &pMesh->iNo           
                   , mesh->numelNov     , mesh->numel        
                   , mesh->nnodeNov     , mesh->nnode 
                   , mesh->maxNo        , mesh->maxViz   
                   , mesh->ndm          , 1
                   , mesh->ndm            
                   , false              , 2);
/*...................................................................*/
  }
/*...................................................................*/

/*... reconstruindo do gradiente (Energia)*/
  if(mesh->ndfFt){
    rcGradU(m                       , loadsTemp
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
           , mesh->elm.material.type, mesh->elm.mat
           , mesh->elm.leastSquare  , mesh->elm.leastSquareR
           , mesh->elm.faceRenergy  , mesh->elm.faceLoadEnergy
           , mesh->elm.temp         , mesh->elm.gradTemp  
           , mesh->node.temp        , sc.rcGrad
           , mesh->maxNo            , mesh->maxViz
           , 1                      , mesh->ndm
           , &pMesh->iNo            , &pMesh->iEl
           , mesh->numelNov         , mesh->numel
           , mesh->nnodeNov         , mesh->nnode); 
/*.................................................................. */


/*... interpolacao das variaveis da celulas para pos nos (GradEnergy)*/
    interCellNode(m                  , loadsTemp  
              , mesh->elm.cellFace   , mesh->face.owner
              , mesh->node.gradTemp  , mesh->elm.gradTemp  
              , mesh->elm.node       , mesh->elm.geomType
              , mesh->elm.geom.cc    , mesh->node.x
              , mesh->face.xm          
              , mesh->elm.nen        , mesh->elm.adj.nViz
              , mesh->elm.faceRenergy, mesh->elm.faceLoadEnergy
              , &pMesh->iNo            
              , mesh->numelNov       , mesh->numel
              , mesh->nnodeNov       , mesh->nnode
              , mesh->maxNo          , mesh->maxViz
              , mesh->ndm            , 1
              , mesh->ndm              
              , false                , 2);  
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (energy)*/
    interCellNode(m                 , loadsTemp
             , mesh->elm.cellFace   , mesh->face.owner
             , mesh->node.temp      , mesh->elm.temp
             , mesh->elm.node       , mesh->elm.geomType
             , mesh->elm.geom.cc    , mesh->node.x
             , mesh->face.xm          
             , mesh->elm.nen        , mesh->elm.adj.nViz
             , mesh->elm.faceRenergy, mesh->elm.faceLoadEnergy
             , &pMesh->iNo            
             , mesh->numelNov       , mesh->numel
             , mesh->nnodeNov       , mesh->nnode
             , mesh->maxNo          , mesh->maxViz
             , 1                    , 1
             , mesh->ndm              
             , opt.bconditions      , 2);  
/*...................................................................*/
  }
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (density)*/
  interCellNode(m                , loadsTemp
             , mesh->elm.cellFace, mesh->face.owner
             , nDenFluid         , mesh->elm.densityFluid
             , mesh->elm.node    , mesh->elm.geomType
             , mesh->elm.geom.cc , mesh->node.x
             , mesh->face.xm   
             , mesh->elm.nen     , mesh->elm.adj.nViz
             , dum               , dum                           
             , &pMesh->iNo         
             , mesh->numelNov    , mesh->numel
             , mesh->nnodeNov    , mesh->nnode
             , mesh->maxNo       , mesh->maxViz
             , 3                 , 1
             , mesh->ndm           
             , false             , 2);  
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (dViscosity)*/
  interCellNode(m                    , loadsTemp
              , mesh->elm.cellFace   , mesh->face.owner
              , nDvisc               , mesh->elm.dViscosity
              , mesh->elm.node       , mesh->elm.geomType
              , mesh->elm.geom.cc    , mesh->node.x
              , mesh->face.xm          
              , mesh->elm.nen        , mesh->elm.adj.nViz
              , dum                  , dum                           
              , &pMesh->iNo            
              , mesh->numelNov       , mesh->numel
              , mesh->nnodeNov       , mesh->nnode
              , mesh->maxNo          , mesh->maxViz
              , 1                    , 1
              , mesh->ndm              
              , false                , 2);
/*...................................................................*/

/*...*/
  if (turbModel->fTurb) {
/*... viscisidade turbulenta*/
    if(opt.eddyViscosity)       
      interCellNode(m                  , loadsVel
                   , mesh->elm.cellFace, mesh->face.owner
                   , nEddyV            , mesh->elm.eddyViscosity        
                   , mesh->elm.node    , mesh->elm.geomType            
                   , mesh->elm.geom.cc , mesh->node.x  
                   , mesh->face.xm   
                   , mesh->elm.nen     , mesh->elm.adj.nViz
                   , dum               , dum                   
                   , &pMesh->iNo          
                   , mesh->numelNov    , mesh->numel        
                   , mesh->nnodeNov    , mesh->nnode 
                   , mesh->maxNo       , mesh->maxViz   
                   , 1                 , 1
                   , mesh->ndm           
                   , false             , 2);
/*...................................................................*/

/*... tensor residual (modelos estruturais)*/
    if(opt.stressR)
      interCellNode(m                     , loadsVel
                , mesh->elm.cellFace      , mesh->face.owner
                , nStressR                , mesh->elm.stressR              
                , mesh->elm.node          , mesh->elm.geomType            
                , mesh->elm.geom.cc       , mesh->node.x  
                , mesh->face.xm             
                , mesh->elm.nen           , mesh->elm.adj.nViz
                , dum                     , dum                   
                , &pMesh->iNo                
                , mesh->numelNov          , mesh->numel        
                , mesh->nnodeNov          , mesh->nnode 
                , mesh->maxNo             , mesh->maxViz   
                , 6                       , 1
                , mesh->ndm                 
                , false                   , 2);
/*...................................................................*/

/*... coeficiente dinamico*/
    if(opt.cDynamic)
      interCellNode(m                      , loadsVel
                , mesh->elm.cellFace       , mesh->face.owner
                , nCdyn                    , mesh->elm.cd              
                , mesh->elm.node           , mesh->elm.geomType            
                , mesh->elm.geom.cc        , mesh->node.x  
                , mesh->face.xm              
                , mesh->elm.nen            , mesh->elm.adj.nViz
                , dum                      , dum                   
                , &pMesh->iNo                 
                , mesh->numelNov           , mesh->numel        
                , mesh->nnodeNov           , mesh->nnode 
                , mesh->maxNo              , mesh->maxViz   
                , 2                        , 1
                , mesh->ndm                  
                , false                    , 2);
/*...................................................................*/

/*... parametros de parede*/
    if(opt.wallParameters)
      interCellNode(m                     , loadsVel
                , mesh->elm.cellFace      , mesh->face.owner
                , nWall                   , mesh->elm.wallParameters              
                , mesh->elm.node          , mesh->elm.geomType            
                , mesh->elm.geom.cc       , mesh->node.x  
                , mesh->face.xm             
                , mesh->elm.nen           , mesh->elm.adj.nViz
                , dum                     , dum                   
                , &pMesh->iNo                
                , mesh->numelNov          , mesh->numel        
                , mesh->nnodeNov          , mesh->nnode 
                , mesh->maxNo             , mesh->maxViz   
                , 4                       , 1
                , mesh->ndm                 
                , false                   , 2);
/*...................................................................*/

/*... energia cinetica turbulenta*/
    if(opt.kTurb)       
      interCellNode(m              , loadsVel
            , mesh->elm.cellFace   , mesh->face.owner
            , nKturb               , mesh->elm.kTurb                
            , mesh->elm.node       , mesh->elm.geomType            
            , mesh->elm.geom.cc    , mesh->node.x  
            , mesh->face.xm      
            , mesh->elm.nen        , mesh->elm.adj.nViz
            , mesh->elm.faceReKturb,mesh->elm.faceLoadKturb                   
            , &pMesh->iNo          
            , mesh->numelNov       , mesh->numel        
            , mesh->nnodeNov       , mesh->nnode 
            , mesh->maxNo          , mesh->maxViz   
            , 1                    , 1
            , mesh->ndm              
            , false                , 2);
/*...................................................................*/

  }                                 
/*...................................................................*/

/*...*/
  if(!mpiVar.myId ){
    fName(preName,sc.ddt.timeStep,0,21,nameOut);
/*...*/
    wResVtkFluid(m                        , mesh0->node.x
               , mesh0->elm.geom.cc       
               , mesh0->elm.node          , mesh0->elm.mat    
               , mesh0->elm.nen           , mesh0->elm.geomType
               , mesh0->elm.pressure      , mesh0->node.pressure
               , mesh0->elm.gradPres      , mesh0->node.gradPres  
               , mesh0->elm.vel           , mesh0->node.vel      
               , mesh0->elm.gradVel       , mesh0->node.gradVel 
               , mesh0->elm.temp          , mesh0->node.temp   
               , mesh0->elm.gradTemp      , mesh0->node.gradTemp
               , mesh0->elm.eddyViscosity , nEddyV
               , mesh0->elm.densityFluid  , nDenFluid
               , mesh0->elm.dViscosity    , nDvisc
               , mesh0->elm.stressR       , nStressR
               , mesh0->elm.cd            , nCdyn
               , mesh0->elm.wallParameters, nWall
               , mesh0->elm.kTurb         , nKturb
               , media->mVel              , nMedVel
               , mesh0->elm.specificHeat  , mesh0->elm.tConductivity                                               
               , mesh0->nnode             , mesh0->numel  
               , mesh0->ndm               , mesh0->maxNo 
               , mesh0->numat             , ndfVel
               , mesh0->ntn               
               , nameOut                  , opt
               , eModel->fKelvin          , media
               , sc.ddt                   , fileOut);   
/*...................................................................*/
  }
/*...................................................................*/

/*... desalocando memoria*/
  if (media->fMedia) {
    HccaDealloc(m, nMedP2Vel, "nMedP2Vel", _AD_);
    HccaDealloc(m, nP2Vel   , "nP2Vel"   , _AD_);
    HccaDealloc(m, nMedVel  , "nMediaVel", _AD_);
  }
  HccaDealloc(m, nKturb   , "nKturb"   , _AD_);
  HccaDealloc(m, nWall    , "nWall"    , _AD_);
  HccaDealloc(m, nCdyn    , "nCdyn"    , _AD_); 
  HccaDealloc(m, nStressR , "nStressR" , _AD_); 
  HccaDealloc(m, nEddyV   , "nEddyV"   , _AD_); 
  HccaDealloc(m, nDvisc   , "nVis"     , _AD_);   
  HccaDealloc(m, nDenFluid, "nDenFluid", _AD_);  
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 05/08/2018                                   *
 * Data de modificaco : 20/08/2018                                   *
 *-------------------------------------------------------------------*
 * printCombustion: impressao do fluido                              * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*  
 * m       -> vetor de memoria principal                             * 
 * turbModel -> modelo de turbulencia                                *
 * eModel    -> modelo da equacao de energia                         *
 * cModel    -> modelo de combustao
 * loadVel   -> deficicao de cargas velocidade                       *
 * loadPres  -> deficicao de cargas pres                             *
 * loadTemp  -> deficicao de cargas temperatura                      *
 * opt       -> opcoes de arquivo                                    *
 * mesh0     -> malha global                                         *
 * mesh      -> malha particionada                                   *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*  
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void printCombustion(Memoria *m      ,Turbulence *turbModel
               ,EnergyModel *eModel  ,Combustion *cModel
               ,PartMesh *pMesh      ,Scheme sc
               ,Loads *loadsVel      ,Loads *loadsPres 
               ,Loads *loadsTemp     ,Loads *loadsComb
               ,FileOpt *opt
               ,Mesh *mesh0          ,Mesh *mesh  
               ,Mean *media      
               ,char *preName        ,char *nameOut){
 
  short ndm = mesh->ndm, nOfPrSp = cModel->nOfSpecies;
  void *dum=NULL;
  short ndfVel,nComb;
  DOUBLE *nStressR=NULL,*nEddyV=NULL,*nDvisc=NULL,*nDenFluid=NULL;
  DOUBLE *nMedVel=NULL,*nP2Vel=NULL,*nMedP2Vel=NULL;
  DOUBLE *nCdyn=NULL,*nWall=NULL,*nKturb=NULL,*nRateFuel=NULL;
  DOUBLE *nYfrac=NULL,*nRaHeReComb=NULL;
  FILE *fileOut=NULL;

/*...*/
  nComb   = cModel->nOfSpeciesLump;
  nOfPrSp = cModel->nOfSpecies; 
  ndfVel = max(mesh->ndfF - 1, mesh->ndfFt - 2);
/*...................................................................*/

/*...*/
  HccaAlloc(DOUBLE, m, nDenFluid, mesh->nnode*3      , "nDenFluid", _AD_);
  HccaAlloc(DOUBLE, m, nDvisc   , mesh->nnode        , "nVis"     , _AD_); 
  HccaAlloc(DOUBLE, m, nEddyV   , mesh->nnode        , "nEddyV"   , _AD_);
  HccaAlloc(DOUBLE, m, nStressR , mesh->nnode*6      , "nStressR" , _AD_);
  HccaAlloc(DOUBLE, m, nCdyn    , mesh->nnode*2      , "nCdyn"    , _AD_); 
  HccaAlloc(DOUBLE, m, nWall    , mesh->nnode*4      , "nWall"    , _AD_); 
  HccaAlloc(DOUBLE, m, nKturb   , mesh->nnode        , "nKturb"   , _AD_);
  HccaAlloc(DOUBLE, m, nYfrac   , mesh->nnode*nOfPrSp, "nYfrac"   , _AD_);
/*...................................................................*/

/*...*/
  if(opt->rateFuel)
    HccaAlloc(DOUBLE, m, nRateFuel, mesh->nnode, "nRateFuel"   , _AD_);
/*...................................................................*/

/*...*/
  if(opt->rateHeatComb)
    HccaAlloc(DOUBLE, m, nRaHeReComb, mesh->nnode, "nRaHeComb",_AD_);
/*...................................................................*/

/*...*/
  if (media->fMedia) {
    HccaAlloc(DOUBLE, m, nMedVel  , mesh->nnode*ndm, "nMediaVel", _AD_);
    zero(nMedVel, mesh->nnode*ndm,DOUBLEC);
    HccaAlloc(DOUBLE, m, nP2Vel   , mesh->nnode*ndm, "nP2Vel", _AD_);
    zero(nMedVel, mesh->nnode*ndm,DOUBLEC);
    HccaAlloc(DOUBLE, m,nMedP2Vel, mesh->nnode*ndm, "nMedP2Vel", _AD_);
    zero(nMedVel, mesh->nnode*ndm,DOUBLEC);
  }
/*...................................................................*/

/*... reconstruindo do gradiente (Pres)*/
  tm.rcGradPres = getTimeC() - tm.rcGradPres;
  rcGradU(m                        , loadsPres
         , mesh->elm.node          , mesh->elm.adj.nelcon
         , mesh->node.x
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
         , mesh->elm.leastSquare   , mesh->elm.leastSquareR
         , mesh->elm.faceRpres     , mesh->elm.faceLoadPres
         , mesh->elm.pressure      , mesh->elm.gradPres                
         , mesh->node.pressure     , sc.rcGrad
         , mesh->maxNo             , mesh->maxViz
         , 1                       , mesh->ndm       
         , &pMesh->iNo             , &pMesh->iEl  
         , mesh->numelNov          , mesh->numel        
         , mesh->nnodeNov          , mesh->nnode);  
  tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (GradPres)*/
  interCellNode(m                  , loadsPres
              , mesh->elm.cellFace , mesh->face.owner
              , mesh->node.gradPres, mesh->elm.gradPres
              , mesh->elm.node     , mesh->elm.geomType            
              , mesh->elm.geom.cc  , mesh->node.x  
              , mesh->face.xm
              , mesh->elm.nen      , mesh->elm.adj.nViz
              , mesh->elm.faceRpres, mesh->elm.faceLoadPres
              , &pMesh->iNo          
              , mesh->numelNov     , mesh->numel        
              , mesh->nnodeNov     , mesh->nnode 
              , mesh->maxNo        , mesh->maxViz   
              , mesh->ndm          , 1
              , mesh->ndm      
              , false              , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (pres)*/
  interCellNode(m                   , loadsPres
               , mesh->elm.cellFace , mesh->face.owner
               , mesh->node.pressure, mesh->elm.pressure   
               , mesh->elm.node     , mesh->elm.geomType            
               , mesh->elm.geom.cc  , mesh->node.x  
               , mesh->face.xm        
               , mesh->elm.nen      , mesh->elm.adj.nViz
               , mesh->elm.faceRpres, mesh->elm.faceLoadPres  
               , &pMesh->iNo           
               , mesh->numelNov     , mesh->numel        
               , mesh->nnodeNov     , mesh->nnode 
               , mesh->maxNo        , mesh->maxViz   
               , 1                  , 1
               , mesh->ndm            
               , opt->bconditions   , 2);
/*...................................................................*/

/*... calculo da matrix jacobiana das velocidades
                            | du1dx1 du1dx2 du1dx3 |   
                            | du2dx1 du2dx2 du2dx3 |   
                            | du3dx1 du3dx2 du3dx3 |
*/        
  tm.rcGradVel  = getTimeC() - tm.rcGradVel;
  rcGradU(m                        , loadsVel
         , mesh->elm.node          , mesh->elm.adj.nelcon
         , mesh->node.x
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
         , mesh->elm.leastSquare   , mesh->elm.leastSquareR
         , mesh->elm.faceRvel      , mesh->elm.faceLoadVel   
         , mesh->elm.vel           , mesh->elm.gradVel                           
         , mesh->node.vel          , sc.rcGrad
         , mesh->maxNo             , mesh->maxViz
         , ndfVel                  , mesh->ndm
         , &pMesh->iNo             , &pMesh->iEl  
         , mesh->numelNov          , mesh->numel        
         , mesh->nnodeNov          , mesh->nnode); 
  tm.rcGradVel = getTimeC() - tm.rcGradVel;  
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (gradVel)*/
  interCellNode(m                   , loadsVel
               , mesh->elm.cellFace , mesh->face.owner
               , mesh->node.gradVel , mesh->elm.gradVel    
               , mesh->elm.node     , mesh->elm.geomType            
               , mesh->elm.geom.cc  , mesh->node.x  
               , mesh->face.xm        
               , mesh->elm.nen      , mesh->elm.adj.nViz
               , mesh->elm.faceRvel , mesh->elm.faceLoadVel 
               , &pMesh->iNo           
               , mesh->numelNov     , mesh->numel        
               , mesh->nnodeNov     , mesh->nnode 
               , mesh->maxNo        , mesh->maxViz   
               , ndfVel             , mesh->ndm
               , mesh->ndm            
               , false              , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (vel)*/
  interCellNode(m                  , loadsVel
               , mesh->elm.cellFace, mesh->face.owner
               , mesh->node.vel    , mesh->elm.vel        
               , mesh->elm.node    , mesh->elm.geomType            
               , mesh->elm.geom.cc , mesh->node.x  
               , mesh->face.xm   
               , mesh->elm.nen     , mesh->elm.adj.nViz
               , mesh->elm.faceRvel, mesh->elm.faceLoadVel 
               , &pMesh->iNo          
               , mesh->numelNov    , mesh->numel        
               , mesh->nnodeNov    , mesh->nnode 
               , mesh->maxNo       , mesh->maxViz   
               , mesh->ndm         , 1
               , mesh->ndm      
               , opt->bconditions  , 2);
/*...................................................................*/

/*... calculo da matrix jacobiana das velocidades
                            | dz1dx1 du1dx2 dz1dx3 |   
                            | dz2dx1 du2dx2 dz2dx3 |   
                            | dz3dx1 dz3dx2 dz3dx3 |
*/        
    rcGradU(m                    , loadsComb
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
        , mesh->elm.material.type, mesh->elm.mat
        , mesh->elm.leastSquare  , mesh->elm.leastSquareR
        , mesh->elm.faceResZcomb , mesh->elm.faceLoadZcomb
        , mesh->elm.zComb        , mesh->elm.gradZcomb
        , mesh->node.zComb       , sc.rcGrad
        , mesh->maxNo            , mesh->maxViz
        , nComb                  , mesh->ndm              
        , &pMesh->iNo            , &pMesh->iEl
        , mesh->numelNov         , mesh->numel
        , mesh->nnodeNov         , mesh->nnode);   
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (gradVel)*/
  interCellNode(m                      , loadsComb
               , mesh->elm.cellFace    , mesh->face.owner
               , mesh->node.gradZcomb  , mesh->elm.gradZcomb  
               , mesh->elm.node        , mesh->elm.geomType            
               , mesh->elm.geom.cc     , mesh->node.x  
               , mesh->face.xm           
               , mesh->elm.nen         , mesh->elm.adj.nViz
               , mesh->elm.faceResZcomb, mesh->elm.faceLoadZcomb 
               , &pMesh->iNo           
               , mesh->numelNov        , mesh->numel        
               , mesh->nnodeNov        , mesh->nnode 
               , mesh->maxNo           , mesh->maxViz   
               , nComb                 , mesh->ndm
               , mesh->ndm            
               , false                 , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (zComb)*/
  interCellNode(m                      , loadsComb
               , mesh->elm.cellFace    , mesh->face.owner
               , mesh->node.zComb      , mesh->elm.zComb  
               , mesh->elm.node        , mesh->elm.geomType            
               , mesh->elm.geom.cc     , mesh->node.x  
               , mesh->face.xm           
               , mesh->elm.nen         , mesh->elm.adj.nViz
               , mesh->elm.faceResZcomb, mesh->elm.faceLoadZcomb  
               , &pMesh->iNo           
               , mesh->numelNov        , mesh->numel        
               , mesh->nnodeNov        , mesh->nnode 
               , mesh->maxNo           , mesh->maxViz   
               , nComb                 , 1
               , mesh->ndm               
               , opt->bconditions      , 2);
/*...................................................................*/

/*... medias*/
  if (media->fMedia) {

/*... medias das velociade*/
    if(media->fVel)
      interCellNode(m                   , loadsVel
                   , mesh->elm.cellFace , mesh->face.owner
                   , nMedVel            , media->mVel        
                   , mesh->elm.node     , mesh->elm.geomType            
                   , mesh->elm.geom.cc  , mesh->node.x  
                   , mesh->face.xm        
                   , mesh->elm.nen      , mesh->elm.adj.nViz
                   , mesh->elm.faceRvel , mesh->elm.faceLoadVel 
                   , &pMesh->iNo           
                   , mesh->numelNov     , mesh->numel        
                   , mesh->nnodeNov     , mesh->nnode 
                   , mesh->maxNo        , mesh->maxViz   
                   , mesh->ndm          , 1
                   , mesh->ndm            
                   , false              , 2);
/*...................................................................*/
  }
/*...................................................................*/

/*... reconstruindo do gradiente (Energia)*/
  if(mesh->ndfFt){
    rcGradU(m                       , loadsTemp
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
           , mesh->elm.material.type, mesh->elm.mat
           , mesh->elm.leastSquare  , mesh->elm.leastSquareR
           , mesh->elm.faceRenergy  , mesh->elm.faceLoadEnergy
           , mesh->elm.temp         , mesh->elm.gradTemp  
           , mesh->node.temp        , sc.rcGrad
           , mesh->maxNo            , mesh->maxViz
           , 1                      , mesh->ndm
           , &pMesh->iNo            , &pMesh->iEl
           , mesh->numelNov         , mesh->numel
           , mesh->nnodeNov         , mesh->nnode); 
/*.................................................................. */


/*... interpolacao das variaveis da celulas para pos nos (GradEnergy)*/
    interCellNode(m                  , loadsTemp  
              , mesh->elm.cellFace   , mesh->face.owner
              , mesh->node.gradTemp  , mesh->elm.gradTemp  
              , mesh->elm.node       , mesh->elm.geomType
              , mesh->elm.geom.cc    , mesh->node.x
              , mesh->face.xm          
              , mesh->elm.nen        , mesh->elm.adj.nViz
              , mesh->elm.faceRenergy, mesh->elm.faceLoadEnergy
              , &pMesh->iNo            
              , mesh->numelNov       , mesh->numel
              , mesh->nnodeNov       , mesh->nnode
              , mesh->maxNo          , mesh->maxViz
              , mesh->ndm            , 1
              , mesh->ndm              
              , false                , 2);  
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (energy)*/
    interCellNode(m                 , loadsTemp
             , mesh->elm.cellFace   , mesh->face.owner
             , mesh->node.temp      , mesh->elm.temp
             , mesh->elm.node       , mesh->elm.geomType
             , mesh->elm.geom.cc    , mesh->node.x
             , mesh->face.xm          
             , mesh->elm.nen        , mesh->elm.adj.nViz
             , mesh->elm.faceRenergy, mesh->elm.faceLoadEnergy
             , &pMesh->iNo            
             , mesh->numelNov       , mesh->numel
             , mesh->nnodeNov       , mesh->nnode
             , mesh->maxNo          , mesh->maxViz
             , 1                    , 1
             , mesh->ndm              
             , opt->bconditions     , 2);  
/*...................................................................*/
  }
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (density)*/
  interCellNode(m                , loadsTemp
             , mesh->elm.cellFace, mesh->face.owner
             , nDenFluid         , mesh->elm.densityFluid
             , mesh->elm.node    , mesh->elm.geomType
             , mesh->elm.geom.cc , mesh->node.x
             , mesh->face.xm   
             , mesh->elm.nen     , mesh->elm.adj.nViz
             , dum               , dum                           
             , &pMesh->iNo         
             , mesh->numelNov    , mesh->numel
             , mesh->nnodeNov    , mesh->nnode
             , mesh->maxNo       , mesh->maxViz
             , 3                 , 1
             , mesh->ndm           
             , false             , 2);  
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (dViscosity)*/
  interCellNode(m                    , loadsTemp
              , mesh->elm.cellFace   , mesh->face.owner
              , nDvisc               , mesh->elm.dViscosity
              , mesh->elm.node       , mesh->elm.geomType
              , mesh->elm.geom.cc    , mesh->node.x
              , mesh->face.xm          
              , mesh->elm.nen        , mesh->elm.adj.nViz
              , dum                  , dum                           
              , &pMesh->iNo            
              , mesh->numelNov       , mesh->numel
              , mesh->nnodeNov       , mesh->nnode
              , mesh->maxNo          , mesh->maxViz
              , 1                    , 1
              , mesh->ndm              
              , false                , 2);
/*...................................................................*/

/*...*/
  if (turbModel->fTurb) 
  {
/*... viscisidade turbulenta*/
    if(opt->eddyViscosity)       
      interCellNode(m                  , loadsVel
                   , mesh->elm.cellFace, mesh->face.owner
                   , nEddyV            , mesh->elm.eddyViscosity        
                   , mesh->elm.node    , mesh->elm.geomType            
                   , mesh->elm.geom.cc , mesh->node.x  
                   , mesh->face.xm   
                   , mesh->elm.nen     , mesh->elm.adj.nViz
                   , dum               , dum                   
                   , &pMesh->iNo          
                   , mesh->numelNov    , mesh->numel        
                   , mesh->nnodeNov    , mesh->nnode 
                   , mesh->maxNo       , mesh->maxViz   
                   , 1                 , 1
                   , mesh->ndm           
                   , false             , 2);
/*...................................................................*/

/*... tensor residual (modelos estruturais)*/
    if(opt->stressR)
      interCellNode(m                     , loadsVel
                , mesh->elm.cellFace      , mesh->face.owner
                , nStressR                , mesh->elm.stressR              
                , mesh->elm.node          , mesh->elm.geomType            
                , mesh->elm.geom.cc       , mesh->node.x  
                , mesh->face.xm             
                , mesh->elm.nen           , mesh->elm.adj.nViz
                , dum                     , dum                   
                , &pMesh->iNo                
                , mesh->numelNov          , mesh->numel        
                , mesh->nnodeNov          , mesh->nnode 
                , mesh->maxNo             , mesh->maxViz   
                , 6                       , 1
                , mesh->ndm                 
                , false                   , 2);
/*...................................................................*/

/*... coeficiente dinamico*/
    if(opt->cDynamic)
      interCellNode(m                      , loadsVel
                , mesh->elm.cellFace       , mesh->face.owner
                , nCdyn                    , mesh->elm.cd              
                , mesh->elm.node           , mesh->elm.geomType            
                , mesh->elm.geom.cc        , mesh->node.x  
                , mesh->face.xm              
                , mesh->elm.nen            , mesh->elm.adj.nViz
                , dum                      , dum                   
                , &pMesh->iNo                 
                , mesh->numelNov           , mesh->numel        
                , mesh->nnodeNov           , mesh->nnode 
                , mesh->maxNo              , mesh->maxViz   
                , 2                        , 1
                , mesh->ndm                  
                , false                    , 2);
/*...................................................................*/

/*... parametros de parede*/
    if(opt->wallParameters)
      interCellNode(m                     , loadsVel
                , mesh->elm.cellFace      , mesh->face.owner
                , nWall                   , mesh->elm.wallParameters              
                , mesh->elm.node          , mesh->elm.geomType            
                , mesh->elm.geom.cc       , mesh->node.x  
                , mesh->face.xm             
                , mesh->elm.nen           , mesh->elm.adj.nViz
                , dum                     , dum                   
                , &pMesh->iNo                
                , mesh->numelNov          , mesh->numel        
                , mesh->nnodeNov          , mesh->nnode 
                , mesh->maxNo             , mesh->maxViz   
                , 4                       , 1
                , mesh->ndm                 
                , false                   , 2);
/*...................................................................*/

/*... energia cinetica turbulenta*/
    if(opt->kTurb)       
      interCellNode(m              , loadsVel
            , mesh->elm.cellFace   , mesh->face.owner
            , nKturb               , mesh->elm.kTurb                
            , mesh->elm.node       , mesh->elm.geomType            
            , mesh->elm.geom.cc    , mesh->node.x  
            , mesh->face.xm      
            , mesh->elm.nen        , mesh->elm.adj.nViz
            , mesh->elm.faceReKturb,mesh->elm.faceLoadKturb                   
            , &pMesh->iNo          
            , mesh->numelNov       , mesh->numel        
            , mesh->nnodeNov       , mesh->nnode 
            , mesh->maxNo          , mesh->maxViz   
            , 1                    , 1
            , mesh->ndm              
            , false                , 2);
/*...................................................................*/

  }                                 
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (gradVel)*/
  if(opt->rateFuel)
    interCellNode(m                     , loadsComb
                , mesh->elm.cellFace    , mesh->face.owner
                , nRateFuel             , mesh->elm.rateFuel
                , mesh->elm.node        , mesh->elm.geomType            
                , mesh->elm.geom.cc     , mesh->node.x  
                , mesh->face.xm           
                , mesh->elm.nen         , mesh->elm.adj.nViz
                , mesh->elm.faceResZcomb, mesh->elm.faceLoadZcomb 
                , &pMesh->iNo           
                , mesh->numelNov        , mesh->numel        
                , mesh->nnodeNov        , mesh->nnode 
                , mesh->maxNo           , mesh->maxViz   
                , 1                     , 1            
                , mesh->ndm               
                , false                 , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (gradVel)*/
  if(opt-> rateHeatComb)
    interCellNode(m                     , loadsComb
                , mesh->elm.cellFace    , mesh->face.owner
                , nRaHeReComb           , mesh->elm.rateHeatReComb
                , mesh->elm.node        , mesh->elm.geomType            
                , mesh->elm.geom.cc     , mesh->node.x  
                , mesh->face.xm           
                , mesh->elm.nen         , mesh->elm.adj.nViz
                , mesh->elm.faceResZcomb, mesh->elm.faceLoadZcomb 
                , &pMesh->iNo           
                , mesh->numelNov        , mesh->numel        
                , mesh->nnodeNov        , mesh->nnode 
                , mesh->maxNo           , mesh->maxViz   
                , 1                     , 1            
                , mesh->ndm               
                , false                 , 2);
/*...................................................................*

/*...*/
  if(opt->yFrac)
    interCellNode(m                     , loadsComb
                , mesh->elm.cellFace    , mesh->face.owner
                , nYfrac                , mesh->elm.yFrac
                , mesh->elm.node        , mesh->elm.geomType            
                , mesh->elm.geom.cc     , mesh->node.x  
                , mesh->face.xm           
                , mesh->elm.nen         , mesh->elm.adj.nViz
                , mesh->elm.faceResZcomb, mesh->elm.faceLoadZcomb 
                , &pMesh->iNo           
                , mesh->numelNov        , mesh->numel        
                , mesh->nnodeNov        , mesh->nnode 
                , mesh->maxNo           , mesh->maxViz   
                , nOfPrSp               , 1            
                , mesh->ndm               
                , false                 , 2);
/*...................................................................*/

/*...*/
  if(!mpiVar.myId ){
    fName(preName,sc.ddt.timeStep,0,30,nameOut);
/*...*/
    wResVtkCombustion(m                        , mesh0->node.x
               , mesh0->elm.geom.cc       
               , mesh0->elm.node          , mesh0->elm.mat    
               , mesh0->elm.nen           , mesh0->elm.geomType
               , mesh0->elm.pressure      , mesh0->node.pressure
               , mesh0->elm.gradPres      , mesh0->node.gradPres  
               , mesh0->elm.vel           , mesh0->node.vel      
               , mesh0->elm.gradVel       , mesh0->node.gradVel 
               , mesh0->elm.temp          , mesh0->node.temp   
               , mesh0->elm.gradTemp      , mesh0->node.gradTemp
               , mesh0->elm.zComb         , mesh0->node.zComb
               , mesh0->elm.gradZcomb     , mesh0->node.gradZcomb
               , mesh0->elm.eddyViscosity , nEddyV
               , mesh0->elm.densityFluid  , nDenFluid
               , mesh0->elm.dViscosity    , nDvisc
               , mesh0->elm.stressR       , nStressR
               , mesh0->elm.cd            , nCdyn
               , mesh0->elm.wallParameters, nWall
               , mesh0->elm.kTurb         , nKturb
               , mesh0->elm.rateFuel      , nRateFuel
               , mesh0->elm.yFrac         , nYfrac
               , mesh0->elm.rateHeatReComb, nRaHeReComb
               , media->mVel              , nMedVel 
               , mesh0->elm.specificHeat  , mesh0->elm.tConductivity                                               
               , mesh0->nnode             , mesh0->numel  
               , mesh0->ndm               , mesh0->maxNo 
               , mesh0->numat             , ndfVel
               , mesh0->ntn               , nOfPrSp
               , nameOut                  , opt
               , eModel->fKelvin          , media
               , sc.ddt                   , fileOut);   
/*...................................................................*/
  }
/*...................................................................*/

/*... desalocando memoria*/
  if (media->fMedia) {
    HccaDealloc(m, nMedP2Vel, "nMedP2Vel", _AD_);
    HccaDealloc(m, nP2Vel   , "nP2Vel"   , _AD_);
    HccaDealloc(m, nMedVel  , "nMediaVel", _AD_);
  }

/*...*/
  if(opt->rateHeatComb)
    HccaDealloc( m, nRaHeReComb, "nRaHeComb",_AD_);
/*...................................................................*/

/*...*/
  if(opt->rateFuel)
    HccaDealloc( m, nRateFuel, "nRateFuel"   , _AD_);
/*...................................................................*/

  HccaDealloc(m, nYfrac   , "nYfrac"   , _AD_);
  HccaDealloc(m, nKturb   , "nKturb"   , _AD_);
  HccaDealloc(m, nWall    , "nWall"    , _AD_);
  HccaDealloc(m, nCdyn    , "nCdyn"    , _AD_); 
  HccaDealloc(m, nStressR , "nStressR" , _AD_); 
  HccaDealloc(m, nEddyV   , "nEddyV"   , _AD_); 
  HccaDealloc(m, nDvisc   , "nVis"     , _AD_);   
  HccaDealloc(m, nDenFluid, "nDenFluid", _AD_);  
/*...................................................................*/

}
/*********************************************************************/



/*********************************************************************
 * Data de criacao    : 02/06/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * printTrans:impressao da equacao de de transporte                  *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * m       -> vetor de memoria principal                             *
 * pMesh     -> modelo de turbulencia                                *
 * sc        -> modelo da equacao de energia                         *
 * loadsT1   -> deficicao de cargas                                  *
 * opt       -> opcoes de arquivo                                    *
 * mesh0     -> malha global                                         *
 * mesh      -> malha particionada                                   *
 * preName   -> prefixo do arquivo                                   *
 * nameOut   -> arquivo de saida                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
*********************************************************************/
void printTrans(Memoria *m
              , PartMesh *pMesh, Scheme *sc
              , Loads *loadsT1 , FileOpt *opt
              , Mesh *mesh0    , Mesh *mesh
              , char *preName  , char *nameOut)
{
  void *dum = NULL;
  char i,str[10][100],*ps[10];
  FILE *fileOut = NULL;
  DOUBLE *nDen = NULL, *nCoefDiff = NULL, *nVel = NULL;

/*...*/
  for(i=0;i<10;ps[i] = str[i],i++);
/*....................................................................*/

/*...*/
  HccaAlloc(DOUBLE,m,nDen     ,mesh->nnode*DENSITY_LEVEL,"nDen",_AD_);
  HccaAlloc(DOUBLE,m,nCoefDiff,mesh->nnode          ,"nCeofDiff",_AD_);
/*...................................................................*/

/*... reconstruindo do gradiente*/
  tm.rcGradT1 = getTimeC() - tm.rcGradT1;
  rcGradU(m                     , loadsT1
         , mesh->elm.node       , mesh->elm.adj.nelcon
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
         , mesh->elm.material.type, mesh->elm.mat
         , mesh->elm.leastSquare  , mesh->elm.leastSquareR
         , mesh->elm.faceRt1      , mesh->elm.faceLoadT1
         , mesh->elm.uT1          , mesh->elm.gradUt1
         , mesh->node.uT1         , sc->rcGrad
         , mesh->maxNo            , mesh->maxViz
         , mesh->ndfT[0]          , mesh->ndm
         , &pMesh->iNo            , &pMesh->iEl
         , mesh->numelNov         , mesh->numel
         , mesh->nnodeNov         , mesh->nnode);
  tm.rcGradT1 = getTimeC() - tm.rcGradT1;
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (Grad)*/
  interCellNode(m                 , loadsT1
              , mesh->elm.cellFace, mesh->face.owner
              , mesh->node.gradUt1, mesh->elm.gradUt1
              , mesh->elm.node    , mesh->elm.geomType
              , mesh->elm.geom.cc , mesh->node.x
              , mesh->face.xm
              , mesh->elm.nen     , mesh->elm.adj.nViz
              , mesh->elm.faceRt1 , mesh->elm.faceLoadT1
              , &pMesh->iNo
              , mesh->numelNov    , mesh->numel
              , mesh->nnodeNov    , mesh->nnode
              , mesh->maxNo       , mesh->maxViz
              , mesh->ndm         , 1
              , mesh->ndm
              , false             , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (uT1)*/
  interCellNode(m                 , loadsT1
              , mesh->elm.cellFace, mesh->face.owner
              , mesh->node.uT1    , mesh->elm.uT1
              , mesh->elm.node    , mesh->elm.geomType
              , mesh->elm.geom.cc , mesh->node.x
              , mesh->face.xm
              , mesh->elm.nen     , mesh->elm.adj.nViz
              , mesh->elm.faceRt1 , mesh->elm.faceLoadT1
              , &pMesh->iNo
              , mesh->numelNov    , mesh->numel
              , mesh->nnodeNov    , mesh->nnode
              , mesh->maxNo       , mesh->maxViz
              , mesh->ndfT[0]     , 1
              , mesh->ndm
              , true              , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (density)*/
  interCellNode(m                 , loadsT1
              , mesh->elm.cellFace, mesh->face.owner
              , nDen              , mesh->elm.densityUt1
              , mesh->elm.node    , mesh->elm.geomType
              , mesh->elm.geom.cc , mesh->node.x
              , mesh->face.xm
              , mesh->elm.nen     , mesh->elm.adj.nViz
              , dum               , dum
              , &pMesh->iNo
              , mesh->numelNov    , mesh->numel
              , mesh->nnodeNov    , mesh->nnode
              , mesh->maxNo       , mesh->maxViz
              , 3                 , 1
              , mesh->ndm
              , false             , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (coefDif)*/
  interCellNode(m                 , loadsD1
              , mesh->elm.cellFace, mesh->face.owner
              , nCoefDiff         , mesh->elm.cDiffT1
              , mesh->elm.node    , mesh->elm.geomType
              , mesh->elm.geom.cc , mesh->node.x
              , mesh->face.xm
              , mesh->elm.nen     , mesh->elm.adj.nViz
              , dum               , dum
              , &pMesh->iNo
              , mesh->numelNov    , mesh->numel
              , mesh->nnodeNov    , mesh->nnode
              , mesh->maxNo       , mesh->maxViz
              , 1                 , 1
              , mesh->ndm
              , false             , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (vel)*/
  interCellNode(m                  , loadsT1
               , mesh->elm.cellFace, mesh->face.owner
               , mesh->node.vel    , mesh->elm.vel
               , mesh->elm.node    , mesh->elm.geomType
               , mesh->elm.geom.cc , mesh->node.x
               , mesh->face.xm
               , mesh->elm.nen     , mesh->elm.adj.nViz
               , dum               , dum                     
               , &pMesh->iNo
               , mesh->numelNov    , mesh->numel
               , mesh->nnodeNov    , mesh->nnode
               , mesh->maxNo       , mesh->maxViz
               , mesh->ndm         , 1
               , mesh->ndm
               , false             , 2);
/*...................................................................*/

/*... globalizacao das variaveis*/
/*... uT1(Node)*/
  dGlobalNode(m              , pMesh
            , mesh0->node.uT1, mesh->node.uT1
            , mesh->ndfT[0]  , 1);
  /*... vel(Node)*/
  dGlobalNode(m              , pMesh
            , mesh0->node.vel, mesh->node.vel
            , mesh->ndm, 1);
/*... gradUt1(Node)*/
  dGlobalNode(m                  , pMesh
            , mesh0->node.gradUt1, mesh->node.gradUt1
            , mesh->ndm          , 1);
/*... uT1(Cel)*/
  dGlobalCel(m             , pMesh
           , mesh0->elm.uT1, mesh->elm.uT1
           , mesh->numelNov
           , mesh->ndfT[0] , 1);
/*... gradUt1(Cel)*/
  dGlobalCel(m                 , pMesh
           , mesh0->elm.gradUt1, mesh->elm.gradUt1
           , mesh->numelNov
           , mesh->ndm         , 1);
/*...................................................................*/

/*...*/
  if (!mpiVar.myId)
  {
    fName(preName, sc->ddt.timeStep, 0, 20, nameOut);

    strcpy(str[0], "elT1");
    strcpy(str[1], "noT1");
    strcpy(str[2], "elGradT1");
    strcpy(str[3], "noGradT1");
    strcpy(str[4], "elVel");
    strcpy(str[5], "nVel");
    strcpy(str[6], "eDensityT1");
    strcpy(str[7], "nDensityT1");
    strcpy(str[8], "nCoefDiffT1");
    strcpy(str[9], "nCoefDiffT1");
/*...*/
    wResVtkTrans(m           , mesh0->node.x
      , mesh0->elm.node      , mesh0->elm.mat
      , mesh0->elm.nen       , mesh0->elm.geomType
      , mesh0->elm.uT1       , mesh0->node.uT1
      , mesh0->elm.gradUt1   , mesh0->node.gradUt1
      , mesh0->elm.vel       , nVel
      , mesh0->elm.densityUt1, nDen
      , mesh0->elm.cDiffT1   , nCoefDiff
      , mesh0->nnode         , mesh0->numel
      , mesh0->ndm           , mesh0->maxNo
      , mesh0->numat         , mesh0->ndfT[0]
      , ps    
      , nameOut              , opt
      , &sc->ddt             , fileOut);
/*...................................................................*/
  }
/*...................................................................*/

/*... desalocando memoria*/
  HccaDealloc(m, nCoefDiff, "nCeofDiff", _AD_);
  HccaDealloc(m, nDen     , "nDen"     , _AD_);
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 01/05/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * printDiff: impressao da equacao de diff                           *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * m       -> vetor de memoria principal                             *
 * pMesh     -> modelo de turbulencia                                *
 * sc        -> modelo da equacao de energia                         *
 * loadsD1   -> deficicao de cargas velocidade                       *
 * opt       -> opcoes de arquivo                                    *
 * mesh0     -> malha global                                         *
 * mesh      -> malha particionada                                   *
 * preName   -> prefixo do arquivo                                   *
 * nameOut   -> arquivo de saida                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void printDiff(Memoria *m
             , PartMesh *pMesh      , Scheme *sc
             , Loads *loadsD1       , FileOpt *opt
             , Mesh *mesh0          , Mesh *mesh
             , char *preName        , char *nameOut)
{
  void *dum = NULL;
  char str1[100], str2[100], str3[100], str4[100];
  FILE *fileOut = NULL;
  DOUBLE *nDenD = NULL, *nCoefDiffD = NULL;

/*...*/
  HccaAlloc(DOUBLE,m, nDenD,mesh->nnode * DENSITY_LEVEL,"nDenD", _AD_);
  HccaAlloc(DOUBLE,m,nCoefDiffD,mesh->nnode, "nCeofDiff", _AD_);
  /*...................................................................*/

/*... reconstruindo do gradiente*/
  tm.rcGradD1 = getTimeC() - tm.rcGradD1;
  rcGradU(m                      , loadsD1
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
        , mesh->elm.material.type, mesh->elm.mat
        , mesh->elm.leastSquare  , mesh->elm.leastSquareR
        , mesh->elm.faceRd1      , mesh->elm.faceLoadD1
        , mesh->elm.uD1          , mesh->elm.gradUd1
        , mesh->node.uD1         , sc->rcGrad
        , mesh->maxNo            , mesh->maxViz
        , mesh->ndfD[0]          , mesh->ndm
        , &pMesh->iNo            , &pMesh->iEl
        , mesh->numelNov         , mesh->numel
        , mesh->nnodeNov         , mesh->nnode);
  tm.rcGradD1 = getTimeC() - tm.rcGradD1;
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (Grad)*/
  interCellNode(m                 , loadsD1
               , mesh->elm.cellFace, mesh->face.owner
               , mesh->node.gradUd1, mesh->elm.gradUd1
               , mesh->elm.node    , mesh->elm.geomType
               , mesh->elm.geom.cc , mesh->node.x
               , mesh->face.xm
               , mesh->elm.nen     , mesh->elm.adj.nViz
               , mesh->elm.faceRd1 , mesh->elm.faceLoadD1
               , &pMesh->iNo
               , mesh->numelNov    , mesh->numel
               , mesh->nnodeNov    , mesh->nnode
               , mesh->maxNo       , mesh->maxViz
               , mesh->ndm         , 1
               , mesh->ndm
               , false             , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (uD1)*/
  interCellNode(m                 , loadsD1
              , mesh->elm.cellFace, mesh->face.owner
              , mesh->node.uD1    , mesh->elm.uD1
              , mesh->elm.node    , mesh->elm.geomType
              , mesh->elm.geom.cc , mesh->node.x
              , mesh->face.xm
              , mesh->elm.nen     , mesh->elm.adj.nViz
              , mesh->elm.faceRd1 , mesh->elm.faceLoadD1
              , &pMesh->iNo
              , mesh->numelNov    , mesh->numel
              , mesh->nnodeNov    , mesh->nnode
              , mesh->maxNo       , mesh->maxViz
              , mesh->ndfD[0]     , 1
              , mesh->ndm
              , true              , 2);
/*...................................................................*/


/*... interpolacao das variaveis da celulas para pos nos (density)*/
  interCellNode(m                 , loadsD1
              , mesh->elm.cellFace, mesh->face.owner
              , nDenD             , mesh->elm.densityUd1
              , mesh->elm.node    , mesh->elm.geomType
              , mesh->elm.geom.cc , mesh->node.x
              , mesh->face.xm
              , mesh->elm.nen     , mesh->elm.adj.nViz
              , dum               , dum
              , &pMesh->iNo
              , mesh->numelNov    , mesh->numel
              , mesh->nnodeNov    , mesh->nnode
              , mesh->maxNo       , mesh->maxViz
              , 3                 , 1
              , mesh->ndm
              , false             , 2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (coefDif)*/
  interCellNode(m                 , loadsD1
              , mesh->elm.cellFace, mesh->face.owner
              , nCoefDiffD        , mesh->elm.cDiffD1
              , mesh->elm.node    , mesh->elm.geomType
              , mesh->elm.geom.cc , mesh->node.x
              , mesh->face.xm
              , mesh->elm.nen     , mesh->elm.adj.nViz
              , dum               , dum
              , &pMesh->iNo
              , mesh->numelNov    , mesh->numel
              , mesh->nnodeNov    , mesh->nnode
              , mesh->maxNo       , mesh->maxViz
              , 1                 , 1
              , mesh->ndm
              , false             , 2);
/*...................................................................*/

/*... globalizacao das variaveis*/
/*... uD1(Node)*/
  dGlobalNode(m              , pMesh
            , mesh0->node.uD1, mesh->node.uD1
            , mesh->ndfD[0]  , 1);

/*... gradUd1(Node)*/
  dGlobalNode(m                   , pMesh
             , mesh0->node.gradUd1, mesh->node.gradUd1
             , mesh->ndm          , 1);
/*... uD1(Cel)*/
  dGlobalCel(m             , pMesh
           , mesh0->elm.uD1, mesh->elm.uD1
           , mesh->numelNov
           , mesh->ndfD[0] , 1);
/*... gradUd1(Cel)*/
  dGlobalCel(m                 , pMesh
           , mesh0->elm.gradUd1, mesh->elm.gradUd1
           , mesh->numelNov
           , mesh->ndm         , 1);
/*...................................................................*/

/*...*/
  if (!mpiVar.myId)
  {
    fName(preName, sc->ddt.timeStep, 0, 8, nameOut);

    strcpy(str1, "elD1");
    strcpy(str2, "noD1");
    strcpy(str3, "elGradD1");
    strcpy(str4, "noGradD1");
/*...*/
    wResVtkDif(m                    , mesh0->node.x
             , mesh0->elm.node      , mesh0->elm.mat
             , mesh0->elm.nen       , mesh0->elm.geomType
             , mesh0->elm.uD1       , mesh0->node.uD1
             , mesh0->elm.gradUd1   , mesh0->node.gradUd1
             , mesh0->elm.densityUd1, nDenD
             , mesh0->elm.cDiffD1   , nCoefDiffD
             , mesh0->nnode         , mesh0->numel
             , mesh0->ndm           , mesh0->maxNo
             , mesh0->numat         , mesh0->ndfD[0]
             , str1                 , str2
             , str3                 , str4
             , nameOut              , opt
             , &sc->ddt             , fileOut);
/*...................................................................*/
  }
/*...................................................................*/

/*... desalocando memoria*/
  HccaDealloc(m, nCoefDiffD, "nCeofDiff", _AD_);
  HccaDealloc(m, nDenD     , "nDenD"    , _AD_);
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 01/05/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * reScaleMesh : redimensio as coordenada da matriz                  *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * m       -> vetor de memoria principal                             *
 * pMesh     -> modelo de turbulencia                                *
 * sc        -> modelo da equacao de energia                         *
 * loadsD1   -> deficicao de cargas velocidade                       *
 * opt       -> opcoes de arquivo                                    *
 * mesh0     -> malha global                                         *
 * mesh      -> malha particionada                                   *
 * preName   -> prefixo do arquivo                                   *
 * nameOut   -> arquivo de saida                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void reScaleMesh(DOUBLE *x,INT nnode, short ndm, FILE *fileIn)
{
  char word[WORD_SIZE], fNomeOut[MAX_STR_LEN_SUFIXO];
  short j;
  int i;
  DOUBLE sC[MAX_NDM];

/*... lendo os valores de escala*/  
  if(ndm == 2)
    fscanf(fileIn,"%lf %lf", &sC[0],&sC[1]);
  else if (ndm == 3)
    fscanf(fileIn, "%lf %lf %lf", &sC[0], &sC[1], &sC[2]);
  readMacro(fileIn, word, false);
/*...................................................................*/

/*... lendo o nome do arquivo de saida*/
  strcpy(fNomeOut,word);
/*...................................................................*/

/*...*/
  for(i=0;i<nnode;i++)
  {
    for (j = 0; j < ndm; j++) 
    {
      MAT2D(i,j,x,ndm) *= sC[j];
    }
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 01/05/2018                                   *
 * Data de modificaco : 11/08/2018                                   *
 *-------------------------------------------------------------------*
 * reScaleMesh : redimensio as coordenada da matriz                  *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * m       -> vetor de memoria principal                             *
 * pMesh     -> modelo de turbulencia                                *
 * sc        -> modelo da equacao de energia                         *
 * loadsD1   -> deficicao de cargas velocidade                       *
 * opt       -> opcoes de arquivo                                    *
 * mesh0     -> malha global                                         *
 * mesh      -> malha particionada                                   *
 * preName   -> prefixo do arquivo                                   *
 * nameOut   -> arquivo de saida                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void initPrintVtk(FileOpt *opt)
{
  opt->fCell          = false;
  opt->fNode          = false;
  opt->gradPres       = false;
  opt->gradVel        = false;
  opt->gradEnergy     = false;
  opt->graduD1        = false;
  opt->graduT1        = false;
  opt->gradZcomb      = false;
  opt->uD1            = false;
  opt->uT1            = false;
  opt->vel            = false;
  opt->pres           = false;
  opt->presTotal      = false;
  opt->energy         = false;
  opt->eddyViscosity  = false;
  opt->densityFluid   = false;
  opt->specificHeat   = false;
  opt->dViscosity     = false;
  opt->tConductivity  = false;
  opt->densityD1      = false;
  opt->coefDiffD1     = false;
  opt->densityT1      = false;
  opt->coefDiffT1     = false;
  opt->vorticity      = false;
  opt->wallParameters = false;
  opt->stress         = false;
  opt->kinetic        = false;
  opt->stressR        = false;
  opt->cDynamic       = false;
  opt->Qcriterion     = false;
  opt->kTurb          = false;
  opt->zComb          = false;
  opt->rateFuel       = false;
  opt->bconditions    = false;
  opt->cc             = false;
  opt->pKelvin        = false;

  opt->stepPlot[0] = 5;
  opt->stepPlot[1] = opt->stepPlot[0];

}
/*********************************************************************
* Data de criacao    : 20/07/2018                                   *
* Data de modificaco : 15/08/2018                                   *
*-------------------------------------------------------------------*
* reScaleMesh : redimensio as coordenada da matriz                  *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* m       -> vetor de memoria principal                             *
* mesh      -> malha global                                         *
* preName   -> prefixo do arquivo                                   *
* fComb     -> modelo de combustao                                  *
* bVtk      -> arquivo vtk binario                                  *
* fileOut  -> ponteiro do arquivo de saida                          *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
void printFace(Memoria *m   , Mesh *mesh
             , char* preName, bool fComb 
             , bool bVtk
             , FILE *fileOut)
{

  char nameOut[SIZEMAX],aux[MAX_STR_LEN_SUFIXO];

/*...*/
  if (mesh->ndfF > 0)
  {
    aux[0] = '\0';
    strcpy(aux, preName);
    strcat(aux, "_vel");
    fName(aux, 0, 0, 17, nameOut);
    wGeoFaceVtk2(m        , mesh->node.x
               , mesh->elm.node    , mesh->elm.nen
               , mesh->elm.geomType
               , mesh->elm.faceRvel, mesh->elm.faceLoadVel
               , mesh->nnode       , mesh->numel
               , mesh->ndm         , mesh->maxViz
               , mesh->ndfF - 1    , mesh->maxNo
               , nameOut           , bVtk
               , fileOut);

    aux[0] = '\0';
    strcpy(aux, preName);
    strcat(aux, "_pres");
    fName(aux, 0, 0, 17, nameOut);
    wGeoFaceVtk2(m                  , mesh->node.x
               , mesh->elm.node     , mesh->elm.nen
               , mesh->elm.geomType
               , mesh->elm.faceRpres, mesh->elm.faceLoadPres
               , mesh->nnode        , mesh->numel
               , mesh->ndm          , mesh->maxViz
               , 1                  , mesh->maxNo
               , nameOut            , bVtk
               , fileOut);
  }
/*..................................................................*/

/*...*/
  if(mesh->ndfFt > 0)
  {
    aux[0] = '\0';
    strcpy(aux,preName);
    strcat(aux,"_energy");
    fName(aux, 0, 0, 17, nameOut);
    wGeoFaceVtk2(m                   , mesh->node.x
               , mesh->elm.node       , mesh->elm.nen
              , mesh->elm.geomType
              , mesh->elm.faceRenergy, mesh->elm.faceLoadEnergy
              , mesh->nnode          , mesh->numel
              , mesh->ndm            , mesh->maxViz
              , 1                    , mesh->maxNo          
              , nameOut              , bVtk                 
              , fileOut);
    
    aux[0] = '\0';
    strcpy(aux, preName);
    strcat(aux, "_vel");
    fName(aux, 0, 0, 17, nameOut);
    wGeoFaceVtk2(m                    , mesh->node.x
               , mesh->elm.node       , mesh->elm.nen
               , mesh->elm.geomType
               , mesh->elm.faceRvel   , mesh->elm.faceLoadVel
               , mesh->nnode          , mesh->numel
               , mesh->ndm            , mesh->maxViz
               , mesh->ndfFt-1        , mesh->maxNo
               , nameOut              , bVtk
               , fileOut);

    aux[0] = '\0';
    strcpy(aux, preName);
    strcat(aux, "_pres");
    fName(aux, 0, 0, 17, nameOut);
    wGeoFaceVtk2(m                  , mesh->node.x
               , mesh->elm.node     , mesh->elm.nen
               , mesh->elm.geomType
               , mesh->elm.faceRpres, mesh->elm.faceLoadPres
               , mesh->nnode        , mesh->numel
               , mesh->ndm          , mesh->maxViz
               , 1                  , mesh->maxNo
               , nameOut            , bVtk
               , fileOut);
  }
/*..................................................................*/

/*...*/
  if (mesh->ndfD[0] > 0)
  {
    aux[0] = '\0';
    strcpy(aux, preName);
    strcat(aux, "_D1");
    fName(aux, 0, 0, 17, nameOut);
    wGeoFaceVtk2(m         , mesh->node.x
      , mesh->elm.node     , mesh->elm.nen
      , mesh->elm.geomType
      , mesh->elm.faceRd1  , mesh->elm.faceLoadD1
      , mesh->nnode        , mesh->numel
      , mesh->ndm          , mesh->maxViz
      , 1                  , mesh->maxNo
      , nameOut            , bVtk
      , fileOut);
  }
/*..................................................................*/

/*...*/
  if (mesh->ndfT[0] > 0)
  {
    aux[0] = '\0';
    strcpy(aux, preName);
    strcat(aux, "_T1");
    fName(aux, 0, 0, 17, nameOut);
    wGeoFaceVtk2(m                 , mesh->node.x
               , mesh->elm.node    , mesh->elm.nen
               , mesh->elm.geomType
               , mesh->elm.faceRt1 , mesh->elm.faceLoadT1
               , mesh->nnode       , mesh->numel
               , mesh->ndm         , mesh->maxViz
               , 1                 , mesh->maxNo
               , nameOut           , bVtk
               , fileOut);
  }
/*..................................................................*/

/*...*/
  if (fComb)
  {
    aux[0] = '\0';
    strcpy(aux, preName);
    strcat(aux, "_combustion");
    fName(aux, 0, 0, 17, nameOut);
    wGeoFaceVtk2(m                     , mesh->node.x
               , mesh->elm.node        , mesh->elm.nen
               , mesh->elm.geomType
               , mesh->elm.faceResZcomb, mesh->elm.faceLoadZcomb
               , mesh->nnode           , mesh->numel
               , mesh->ndm             , mesh->maxViz
               , 1                     , mesh->maxNo
               , nameOut               , bVtk
               , fileOut);
  }
/*..................................................................*/

}
/********************************************************************/
