#include<Energy.h>
/********************************************************************* 
 * Data de criacao    : 16/01/2018                                   *
 * Data de modificaco : 15/08/2018                                   * 
 *-------------------------------------------------------------------* 
 * energyEquation: Solucao da equação de energia                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 *********************************************************************/
bool energyEquation(Memoria *m               , PropVarFluid *prop 
                   , Loads *loadsVel         , Loads *loadsEnergy  
                   , EnergyModel *eModel     , Turbulence *tModel  
                   , Combustion *cModel
                   , ThermoDynamic *thDynamic, Mesh *mesh          
                   , SistEq *sistEqEnergy    , Solv *solvEnergy
                   , Simple *sp  
                   , Scheme *sc              , PartMesh *pMesh){

  bool fEnergy;
  DOUBLE tmp;

/*...*/
  tmp = 0.e0;
/*.................................................................. */

/*... reconstruindo do gradiente (Energia)*/
  tm.rcGradEnergy = getTimeC() - tm.rcGradEnergy;
  rcGradU(m                       , loadsEnergy
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
         , mesh->elm.faceRenergy  ,mesh->elm.faceLoadEnergy
         , mesh->elm.energy       ,mesh->elm.gradEnergy
         , mesh->node.energy      ,sc->rcGrad
         , mesh->maxNo            ,mesh->maxViz
         , 1, mesh->ndm
         , &pMesh->iNo            ,&pMesh->iEl
         , mesh->numelNov         ,mesh->numel
         , mesh->nnodeNov         ,mesh->nnode);
  tm.rcGradEnergy = getTimeC() - tm.rcGradEnergy;
/*.................................................................. */

/*... calculo de: A(i),bE(i)*/
  tm.systFormEnergy = getTimeC() - tm.systFormEnergy;
  systFormEnergy(loadsEnergy          , loadsVel
            , &sc->advEnergy          , &sc->diffEnergy
            , tModel                  , eModel  
            , cModel                  , prop                    
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
            , sistEqEnergy->ia        , sistEqEnergy->ja
            , sistEqEnergy->al        , sistEqEnergy->ad
            , sistEqEnergy->b         , sistEqEnergy->id
            , mesh->elm.faceRenergy   , mesh->elm.faceLoadEnergy
            , mesh->elm.faceRvel      , mesh->elm.faceLoadVel   
            , mesh->elm.energy        , mesh->elm.gradEnergy
            , mesh->elm.vel           , mesh->elm.gradVel
            , mesh->elm.pressure0     , mesh->elm.pressure  
            , mesh->elm.gradPres      , mesh->elm.rCellEnergy 
            , mesh->elm.densityFluid  , mesh->elm.specificHeat
            , mesh->elm.dViscosity    , mesh->elm.eddyViscosity 
            , mesh->elm.tConductivity
            , mesh->elm.enthalpyk     , mesh->elm.gradY          
            , mesh->elm.cDiffComb     , mesh->elm.yFrac
            , mesh->elm.rateHeatReComb
            , sp->d                   , mesh->elm.wallParameters
            , sc->ddt                 , sp->alphaEnergy
            , sistEqEnergy->neq       , sistEqEnergy->neqNov
            , sistEqEnergy->nad       , sistEqEnergy->nadr
            , mesh->maxNo             , mesh->maxViz
            , mesh->ndm               , mesh->numelNov
            , 1                       , sistEqEnergy->storage
            , true                    , true
            , true                    , sistEqEnergy->unsym);
  tm.systFormEnergy = getTimeC() - tm.systFormEnergy;
/*...................................................................*/

/*... soma o vetor bE(i) = bE(i) + bE0(i)*/
  addVector(1.0e0               ,sistEqEnergy->b
           ,1.0e0               ,sistEqEnergy->b0
           ,sistEqEnergy->neqNov,sistEqEnergy->b);
/*...................................................................*/

/*... soma o vetor RE(i) = RE(i) + bE0(i)*/
  updateCellValue(mesh->elm.rCellEnergy,sistEqEnergy->b0
                  ,sistEqEnergy->id     ,&sistEqEnergy->iNeq
                  ,mesh->numelNov       ,1
                  ,true                 ,false);
/*...................................................................*/

/*...*/
  tmp = sqrt(dot(sistEqEnergy->b,sistEqEnergy->b
                  ,sistEqEnergy->neqNov));
  fEnergy = true;
  if (tmp == 0.e0) fEnergy = false;  
/*...................................................................*/

/*...Ae=bE*/
  if (fEnergy) {
    tm.solvEnergy = getTimeC() - tm.solvEnergy;
    solverC(m
        ,sistEqEnergy->neq    ,sistEqEnergy->neqNov
        ,sistEqEnergy->nad    ,sistEqEnergy->nadr
        ,sistEqEnergy->ia     ,sistEqEnergy->ja
        ,sistEqEnergy->al     ,sistEqEnergy->ad, sistEqEnergy->au
        ,sistEqEnergy->b      ,sistEqEnergy->x
        ,&sistEqEnergy->iNeq  ,&sistEqEnergy->omp
        ,solvEnergy->tol      ,solvEnergy->maxIt
        ,sistEqEnergy->storage,solvEnergy->solver
        ,solvEnergy->fileSolv ,solvEnergy->log
        ,true                 ,sistEqEnergy->unsym);  
    tm.solvEnergy = getTimeC() - tm.solvEnergy;
  }
/*...................................................................*/

/*... x -> Energy*/
  updateCellValue(mesh->elm.energy,sistEqEnergy->x
                 ,sistEqEnergy->id,&sistEqEnergy->iNeq
                 ,mesh->numel     ,1
                 ,eModel->fRes    ,true);
/*...................................................................*/

  return fEnergy;

}