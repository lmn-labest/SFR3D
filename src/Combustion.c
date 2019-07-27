#include<Combustion.h>

DOUBLE edc2(Combustion *cModel,Prop *sHeatPol 
          ,DOUBLE *y          ,short const iYf
          ,short const iYox      
          ,short *iProd       ,short const nProd  
          ,DOUBLE const s     ,DOUBLE const density
          ,DOUBLE const vol   ,DOUBLE const eddyVisc
          ,DOUBLE *c          ,DOUBLE const modS 
          ,DOUBLE const dVisc ,DOUBLE const df 
          ,DOUBLE const tMix  ,DOUBLE const tChemical
          ,DOUBLE const temp   
          ,short const iCod   ,bool const fKelvin);



/*********************************************************************
 * Data de criacao    : 30/07/2018                                   *
 * Data de modificaco : 24/05/2019                                   *
 *-------------------------------------------------------------------*
 * combustionModel : modelo de combustao                             *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 *********************************************************************/
void combustionModel(Memoria *m         , PropVarFluid *prop
                   , Loads *loadsComb   , Loads *loadsVel
                   , Turbulence *tModel , Combustion *cModel
                   , EnergyModel *eModel, Mesh *mesh
                   , SistEq *sistEqComb , Solv *solvComb  
                   , Simple *sp         , Scheme *sc        
                   , PartMesh *pMesh    , FileOpt *opt
                   , bool *fComb        , short itSimple    ) 
{
  bool fSheat = prop->fSpecificHeat;
  short conv, i, nComb;
  INT desloc;
  DOUBLE tmp,tb[MAX_COMB];
  DOUBLE *b[MAX_COMB], *xu[MAX_COMB];
  DOUBLE *ad[MAX_COMB],*al[MAX_COMB],*au[MAX_COMB],*rCellC[MAX_COMB];
/*... nComb - numero de especies transportadas*/     
  nComb   = cModel->nComb;
/*...................................................................*/

/*...*/
  for(i=0;i<nComb; i++)
  {
    fComb[i] = false; 

    desloc = sistEqComb->neq;
    b[i] = &sistEqComb->b[i*desloc];
  
    xu[i] = &sistEqComb->x[i*desloc]; 

    ad[i] = &sistEqComb->ad[i*desloc];

    desloc = 2*sistEqComb->nad + sistEqComb->nadr;
    al[i] = &sistEqComb->al[i*desloc];
 
    au[i] = &sistEqComb->au[i*desloc];

    desloc = mesh->numel;
    rCellC[i] = &mesh->elm.rCellComb[i*desloc];
  }
/*...................................................................*/

/*... taxa de comsumo do combustivel*/
/*m.fuelConsume  = getTimeC() - tm.fuelConsume;
  timeChemical(cModel                  , tModel
              , mesh->elm.zComb        , mesh->elm.temp      
              , mesh->elm.densityFluid , mesh->elm.gradVel 
              , mesh->elm.eddyViscosity
              , mesh->elm.dViscosity   , mesh->elm.tReactor
              , mesh->ndm 
              , mesh->numel            , eModel->fKelvin ); 
  tm.fuelConsume  = getTimeC() - tm.fuelConsume;  */
/*...................................................................*/


/*... taxa de comsumo do combustivel*/
  tm.fuelConsume  = getTimeC() - tm.fuelConsume;
  rateReaction(cModel                   , tModel
                  , &prop->sHeat
                  , mesh->elm.zComb        , mesh->elm.cDiffComb
                  , mesh->elm.temp         , mesh->elm.wk
                  , mesh->elm.densityFluid , mesh->elm.gradVel
                  , mesh->elm.eddyViscosity,mesh->elm.dViscosity
                  , mesh->elm.geom.volume
                  , mesh->ndm              , mesh->numel
                  , eModel->fKelvin );  
  tm.fuelConsume  = getTimeC() - tm.fuelConsume;  
/*...................................................................*/

/*... reconstruindo do gradiente (gradZ)*/
  tm.rcGradComb   = getTimeC() - tm.rcGradComb;
  rcGradU(m                      , loadsComb
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
        , mesh->node.zComb       , sc->rcGrad
        , mesh->maxNo            , mesh->maxViz
        , nComb                  , mesh->ndm              
        , &pMesh->iNo            , &pMesh->iEl
        , mesh->numelNov         , mesh->numel
        , mesh->nnodeNov         , mesh->nnode);      
  tm.rcGradComb = getTimeC() - tm.rcGradComb;
/*.................................................................. */

/*... calculo de: A(i),b(i)*/
  tm.systFormComb = getTimeC() - tm.systFormComb;
  systFormComb(loadsComb             , loadsVel
             , &sc->advComb            , &sc->diffComb  
             , tModel                  , cModel
             , prop                    
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
             , sistEqComb->ia          , sistEqComb->ja
             , sistEqComb->al          , sistEqComb->ad
             , sistEqComb->b           , sistEqComb->id
             , mesh->elm.faceResZcomb  , mesh->elm.faceLoadZcomb
             , mesh->elm.faceRvel      , mesh->elm.faceLoadVel
             , mesh->elm.zComb         , mesh->elm.gradZcomb  
             , mesh->elm.wk            , mesh->elm.vel           
             , mesh->elm.pressure0     , mesh->elm.pressure
             , mesh->elm.gradPres      , mesh->elm.rCellComb
             , mesh->elm.densityFluid  , mesh->elm.cDiffComb
             , mesh->elm.eddyViscosity , mesh->elm.wallParameters
             , sp->d
             , sc->ddt                , sp->alphaComb
             , sistEqComb->neq        , sistEqComb->neqNov
             , sistEqComb->nad        , sistEqComb->nadr
             , mesh->maxNo            , mesh->maxViz
             , mesh->ndm              , mesh->numelNov
             , cModel->nComb          , sistEqComb->storage
             , true                   , true
             , true                   , sistEqComb->unsym);
  tm.systFormComb = getTimeC() - tm.systFormComb;
/*...................................................................*/

/*... soma o vetor b(i) = b(i) + b0(i)*/
  addVector(1.0e0                   , sistEqComb->b
          , 1.0e0                   , sistEqComb->b0
          , sistEqComb->neqNov*nComb, sistEqComb->b);
/*...................................................................*/

/*... soma o vetor R(i) = R(i) + b0(i)*/
  updateCellValueSimple(mesh->elm.rCellComb  , sistEqComb->b0
                      , sistEqComb->id       , &sistEqComb->iNeq
                      , mesh->numelNov       , sistEqComb->neqNov
                      , nComb
                      , true                 , false);  
/*...................................................................*/
    
/*...*/
  for(i=0,conv = 0;i<nComb;i++)
  {
    tb[i] = sqrt(dot(b[i], b[i], sistEqComb->neqNov));
    if (itSimple == 0) tmp = maxArray(tb,nComb);
/*...................................................................*/

 /*...*/
    fComb[i] = true; 
    if ( tb[i] < SZERO || tb[i] == 0.e0 ) fComb[i] = false;
//  printf("%d %e\n",i,tb[i]);
  }
/*...................................................................*/

/*...*/
  for(i=0,conv = 0;i<nComb;i++)
  {
/*...Ax=b*/
    tm.solvComb = getTimeC() - tm.solvComb;
    if(fComb[i])
      solverC(m
              , sistEqComb->neq    , sistEqComb->neqNov
              , sistEqComb->nad    , sistEqComb->nadr
              , sistEqComb->ia     , sistEqComb->ja
              , al[i]              , ad[i]           , au[i]
              , b[i]               , xu[i]              
              , &sistEqComb->iNeq  , &sistEqComb->omp
              , solvComb->tol      , solvComb->maxIt
              , sistEqComb->storage, solvComb->solver
              , solvComb->fileSolv , solvComb->log
              , true               , sistEqComb->unsym);
    tm.solvComb = getTimeC() - tm.solvComb;    
  }
/*...................................................................*/

/*... x -> zComb*/
  updateCellValueBlock(mesh->elm.zComb , sistEqComb->x
                     , sistEqComb->id  , &sistEqComb->iNeq
                     , mesh->numel     , sistEqComb->neq 
                     , cModel->nComb   
                     , cModel->fRes    , true);  
/*...................................................................*/

/*...*/
  tm.speciesLoop = getTimeC() - tm.speciesLoop;
/*...................................................................*/

/*...*/
  getSpeciesPrimitives(cModel
                      ,mesh->elm.yFrac,mesh->elm.zComb
                      ,mesh->numelNov);
/*...................................................................*/

/*...*/
  regularZ(mesh->elm.yFrac,mesh->numelNov,cModel->nOfSpecies);
/*...................................................................*/

/*...*/
  getGradSpecies(cModel    
                , mesh->elm.gradZcomb, mesh->elm.gradY
                , mesh->numel        , mesh->ndm);
/*...................................................................*/

  tm.speciesLoop = getTimeC() - tm.speciesLoop;
/*...................................................................*/

/*...*/
  tm.enthalpySpecies = getTimeC() - tm.enthalpySpecies;
  getEnthalpySpecies(cModel             , prop
                   , mesh->elm.enthalpyk, mesh->elm.temp 
                   , mesh->numel        , eModel->fKelvin
                   , ompVar.fUpdate     , ompVar.nThreadsUpdate);
  tm.enthalpySpecies = getTimeC() - tm.enthalpySpecies;
/*...................................................................*/

/*...*/
  tm.heatRelease  = getTimeC() - tm.heatRelease; 
  rateHeatReakesedReaction(cModel            , &prop->sHeat                
                    , mesh->elm.rateHeatReComb, mesh->elm.temp     
                    , mesh->elm.zComb0        , mesh->elm.zComb
                    , mesh->elm.densityFluid  , mesh->elm.wk
                    , mesh->elm.material.prop , mesh->elm.mat    
                    , sc->ddt.dt[TIME_N]      , mesh->numelNov
                    , fSheat                  , eModel->fKelvin
                    , ompVar.fUpdate     , ompVar.nThreadsUpdate); 
  tm.heatRelease  = getTimeC() - tm.heatRelease; 
/*...................................................................*/

/*...*/
  prop->molarMass = mixtureMolarMass(cModel,mesh->elm.yFrac);
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 15/08/2018                                   *
 * Data de modificaco : 08/06/2019                                   *
 *-------------------------------------------------------------------*
 * regularZ : mantem o z entre 0 e 1                                 *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * z            -> fracao massica                                    *
 * numel        -> numero de elementos                               *
 * nOfSpecies   -> numero de especie agrupadas                       *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * z -> fracao de massa de regularizadas                             * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * z < 0.0 -> z = 0.e0                                               *
 *********************************************************************/
void regularZ(DOUBLE *RESTRICT y, INT const numel, short const ns)
{
  short j;
  INT nel;

/*... z < 0.0*/
  for(nel = 0 ; nel < numel; nel++)
  {
    for(j=0;j<ns;j++)
      if(MAT2D(nel,j,y,ns) < 0.e0)
        MAT2D(nel,j,y,ns) = 0.e0;
  }

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 04/06/2019                                   *
 * Data de modificaco : 09/06/2019                                   *
 *-------------------------------------------------------------------*
 * getEnthalpySpecies : obetem as entalpias por especies             *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * z            -> fracao massica                                    *
 * enthalpyk    -> nao definido                                      *
 * temp         -> temperatura                                       *
 * numel        -> numero de celulas                                 *
 * fOmp         -> openmp                                            *
 * nThreads     -> numero de threads                                 *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * enthalpyk    -> entalpia por especie                              * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void getEnthalpySpecies(Combustion *cModel        ,  PropVarFluid *propF
                      , DOUBLE *RESTRICT enthalpyk, DOUBLE *RESTRICT temp 
                      , INT const numel           , bool const fKelvin
                      , bool const fOmp           , short const nThreads )
{
  short j,ns=cModel->nOfSpecies;
  INT nel;
  DOUBLE hs,sHeatRef;

  sHeatRef = propF->sHeatRef;
/*... omp*/
  if(fOmp)
  {
#pragma omp parallel  for default(none) num_threads(nThreads)\
        private(nel,j,hs) shared(propF,cModel,enthalpyk,temp,numel,ns,fKelvin,sHeatRef)
    for(nel = 0 ; nel < numel; nel++)
    { 
      for(j=0;j<ns;j++)
      {
        hs = tempForSpecificEnthalpySpecies(&propF->sHeat       , j
                                           , temp[nel]           , sHeatRef
                                           , propF->fSpecificHeat, fKelvin);
        MAT2D(nel,j,enthalpyk,ns) = hs;
      }
    }
  }
/*.....................................................................*/

/*... seq*/
  else
  {
    for(nel = 0 ; nel < numel; nel++)
    { 
      for(j=0;j<ns;j++)
      {
        hs = tempForSpecificEnthalpySpecies(&propF->sHeat       , j
                                           , temp[nel]           , sHeatRef
                                           , propF->fSpecificHeat, fKelvin);
        MAT2D(nel,j,enthalpyk,ns) = hs;
      }
    }
  }
/*.....................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 05/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * getGradSpecies : obetem o gradiente te todas das especies         *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel       -> fracao massica                                    *
 * gradZ        -> gradiente das especie resolvidas explicitamente   *
 * numel        -> numero de celulas                                 *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * gradZ        -> gradiente das especie resolvidas explicitamente   *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void  getGradSpecies(Combustion *cModel   
                   , DOUBLE *RESTRICT gradZ, DOUBLE *RESTRICT gradY
                   , INT const numel       , short const ndm)
{
  short j,k,ns=cModel->nOfSpecies,nc=cModel->nComb;
  INT nel;

  if(ns!=nc)
  {
    for(nel = 0 ; nel < numel; nel++)
    {

      for(k=0;k<ndm;k++)
      {
        MAT3D(nel,ns-1,k,gradY,ns,ndm) = 0.e0;
        for(j=0;j<nc;j++)
        {
          MAT3D(nel,j,k,gradY,ns,ndm) = MAT3D(nel,j,k,gradZ,nc,ndm);
/*... grad(ZN) = - sum(gradZ)*/
          MAT3D(nel,ns-1,k,gradY,ns,ndm) -= MAT3D(nel,j,k,gradZ,nc,ndm);
        }
      }
    }
  }
  else
  {
   for(nel = 0 ; nel < numel; nel++)
     for(j=0;j<ns;j++)
       for(k=0;k<ndm;k++)
      {
        MAT3D(nel,j,k,gradY,ns,ndm) = MAT3D(nel,j,k,gradZ,nc,ndm); 
      }

  }
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 15/08/2018                                   *
 * Data de modificaco : 27/05/2019                                   *
 *-------------------------------------------------------------------*
 * getSpeciesPrimitives : obtem as especies primitivas das especies  *
 * agrupadas                                                         *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> modelo de combustao                                    *
 * y       -> nao definido                                           *
 * zComb   -> fracao massica                                         *
 * numel   -> numero de elementos                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * y -> fracao de massa de especies primitivas                       * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * z(nel,1) -> air (oxidante)                                        *
 * z(nel,2) -> fuel(comburente)                                      *
 * z(nel,3) -> produto                                               *
 *                                                                   * 
 * y[0] - CH4                                                        *
 * y[1] - O2                                                         *
 * y[2] - N2                                                         *
 * y[3] - CO2                                                        *
 * y[4] - H2O                                                        *
 * y[5] - CO                                                         *
 * y[6] - C                                                          *
 *********************************************************************/
void getSpeciesPrimitives(Combustion *cModel 
                        , DOUBLE *RESTRICT y,DOUBLE *RESTRICT z
                        , INT const numel   )
{
  short i,ns,nl,nc;
  INT nel;
  DOUBLE *py,*pz,*a,zLump[3],sum;

  if(cModel->fLump)
  {
    ns = cModel->nOfSpecies; 
    nl = cModel->nOfSpeciesLump;
    nc = cModel->nComb;
    a  = cModel->lumpedMatrix;
/*... resolvendo as N especies*/
    if (nc == nl)
      for(nel = 0 ; nel < numel; nel++)
      {
        py = &MAT2D(nel,0,y,ns);
        pz = &MAT2D(nel,0,z,nc);
        yLumpedMatrixZ(py,a,pz,ns,nl);
      }
    else
      for(nel = 0 ; nel < numel; nel++)
      {
        py = &MAT2D(nel,0,y,ns);
        zLump[SL_FUEL] = MAT2D(nel,SL_FUEL,z,nc);
        zLump[SL_AIR]  = MAT2D(nel,SL_AIR ,z,nc);
        zLump[SL_PROD] = 1.e0 - (zLump[SL_AIR] + zLump[SL_FUEL]);
        yLumpedMatrixZ(py,a,zLump,ns,nl);
      }
  }
  else
  {
    ns = cModel->nOfSpecies; 
    nc = cModel->nComb;
/*... resolvendo as N especies*/
    if (ns == nc)
      for(nel = 0 ; nel < numel; nel++)
        for(i=0;i<nc;i++)
          MAT2D(nel,i,y,ns) = MAT2D(nel,i,z,nc);

/*... resolvendo as N - 1 especies*/
    else
      for(nel = 0 ; nel < numel; nel++)
      {
      
        for(i=0;i<nc;i++)
          MAT2D(nel,i,y,ns) = MAT2D(nel,i,z,nc);

        for(i=0,sum=0.e0;i<nc;i++)
          sum+=MAT2D(nel,i,y,ns);
        
        sum = min(sum,1.e0);

        MAT2D(nel,cModel->chem.sN2,y,ns) = 1.e0 - sum;
      }
  }
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 15/05/2019                                   *
 * Data de modificaco : 27/05/2019                                   *
 *-------------------------------------------------------------------*
 * getSpeciesPrimitives : obtem as especies primitivas das especies  *
 * condicoes de contorno                                             *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> modelo de combustao                                    *
 * y       -> nao definido                                           *
 * zComb   -> fracao massica                                         *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * y -> fracao de massa de especies primitivas                       * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void getSpeciesPrimitivesCc(Combustion *cModel 
                           , DOUBLE *RESTRICT y,DOUBLE *RESTRICT z)
{
  short i,ns,nl,nc;
  DOUBLE zLumped[3],sum;
  if(cModel->fLump)
  {
    ns = cModel->nOfSpecies; 
    nl = cModel->nOfSpeciesLump; 
    nc = cModel->nComb;
    if (nc != nl)
    {
      zLumped[SL_FUEL] = z[SL_FUEL];
      zLumped[SL_AIR]  = z[SL_AIR];
      zLumped[SL_PROD] = 1.e0 - (z[SL_FUEL] + z[SL_AIR]);
    }
    else
    {
      zLumped[SL_FUEL] = z[SL_FUEL];
      zLumped[SL_AIR]  = z[SL_AIR];
      zLumped[SL_PROD] = z[SL_PROD];
    }
    yLumpedMatrixZ(y,cModel->lumpedMatrix,zLumped,ns,nl);
  }
  else
  {
    ns = cModel->nOfSpecies; 
    nc = cModel->nComb;
/*... resolvendo as N especies*/
    if (ns == nc)
    {
      for(i=0;i<nc;i++)
        y[i] = z[i];
    }
/*... resolvendo as N - 1 especies*/
    else
    {
      for(i=0;i<nc;i++)
        y[i] = z[i];

      for(i=0,sum=0.e0;i<nc;i++)
        sum+=y[i];
      y[cModel->chem.sN2] = 1.e0 - sum;

    }
  }
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 25/07/2019                                   *
 *-------------------------------------------------------------------*
 * rateHeatReakesedReaction: calculo da taxa de liberacao de calor  *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> modelo de combustao                                    *
 * sHeatPol - estrutra par o polinoimio do calor especifico          *
 * q       -> nao definido                                           *
 * zComb0  -> fracao de massa agrupada do passo de tempo anterior    *
 * zComb   ->fracao de massa agrupada do passo de tempo atural       *
 * wk        -> taxa de consumo massico das especies kg/(m3 s)       * 
 * prop   - propriedades por material                                *
 * mat    - material da celula                                       * 
 * dt      -> delta dessa passo de tempo                             * 
 * numel   -> numero de elementos                                    *
 * fSheat   - calor especifico com variacao com a Temperatura        *
 * fKelvin  - temperatura dada em kelvin                             *
 * fOmp         -> openmp                                            *
 * nThreads     -> numero de threads                                 *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * q -> taxa de energia liberada KJ/s                                * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * 1) aCH4 + bO2 + cN2 -> dCO + eH2O + fN2                           *
 * nSpeciePart[0][0] = 3                                             *
 * nSpeciePart[0][1] = 3                                             * 
 * Reagente                                                          *
 * speciesPart[0][0][0] = SP_FUEL;                                   *
 * speciesPart[0][0][1] = SP_O2;                                     *
 * speciesPart[0][0][2] = posN2;                                     *
 * Produto                                                           *
 * speciesPart[0][1][0] = SP_CO;                                     *
 * speciesPart[0][1][1] = SP_H2O;                                    *
 * speciesPart[0][1][2] = posN2;                                     *
 * 2) aCO + bO2 -> dCO2                                              *
 * nSpeciePart[1][0] = 2                                             *
 * nSpeciePart[1][1] = 1                                             *                      
 * Reagente                                                          *
 * speciesPart[1][0][0] = SP_CO;                                     *
 * speciesPart[1][0][1] = SP_O2;                                     *
 * Produto                                                           *
 * speciesPart[1][1][0] = SP_CO2;                                    *
 *********************************************************************/
void rateHeatReakesedReaction(Combustion *cModel,Prop *sHeat   
                   , DOUBLE *RESTRICT q      , DOUBLE *RESTRICT temp
                   , DOUBLE *RESTRICT zComb0 , DOUBLE *RESTRICT zComb
                   , DOUBLE *RESTRICT density, DOUBLE *RESTRICT wk
                   , DOUBLE *RESTRICT prop   , short  *RESTRICT mat
                   , DOUBLE const dt         , INT const numel
                   , bool const fsHeat       , bool const fKelvin
                   , bool const fOmp         , short const nThreads)
{

  short i
      , iCod = cModel->typeHeatRealese
      , nReac = cModel->chem.nReac
      , nSp   = cModel->chem.nSp;
        
  INT nel;
  DOUBLE hc,h[MAXSPECIES];
  DOUBLE sum;

  switch(iCod) 
  {
/*... Entalpia de formacao*/
    case  HFORMATION:
    
    for(i=0;i<nSp;i++)
      h[i] = cModel->chem.sp[i].entalphyOfForm;

    if (fOmp)
    {
#pragma omp parallel  for default(none) num_threads(nThreads)\
        private(nel,i,sum,nSp,h)\
        shared(cModel,temp,nReac,wk,q)    
      for(nel = 0; nel < numel; nel++)
      {
/*... KJ/KG*/
         for(i=0,sum=0.e0;i<nSp;i++)
            sum += wk[i]*h[i];
        q[nel] = -sum; 
      }
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    else
    {
      for(nel = 0; nel < numel; nel++)
      {
/*... KJ/KG*/
         for(i=0,sum=0.e0;i<nSp;i++)
            sum += wk[i]*h[i];
        q[nel] = -sum; 
      }
/*...................................................................*/
    }
    break;
/*...................................................................*/

/*... Entalpia de combustao*/
  case HCOMBUSTION:
    printf("Nao implementado!!\n");
    exit(EXIT_FAILURE);
    break;
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * rateReaction: calculo da taxa de consumo do combustivel        *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * q       -> nao definido                                           *
 * rate    -> taxa de consumo de combustivel                         * 
 * numel   -> numero de elementos                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * q -> energia total liberada KJ                                    * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE totalHeatRealeseComb(DOUBLE *RESTRICT q, DOUBLE *RESTRICT vol  
                          , DOUBLE const dt   , INT const numel)
{

  INT nel;
  DOUBLE qTotal =0.e0, vTotal = 0.e0;

  for(nel = 0; nel < numel; nel++)
    qTotal += q[nel]*vol[nel];

  return qTotal*dt;
}
/*********************************************************************/

/**********************************************************************
 * Data de criacao    : 15/08/2018                                    *
 * Data de modificaco : 10/05/2019                                    *
 *------------------------------------------------------------------- * 
 * sumFracZ : Soma todas as fracoes de massa                          *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * z     -> nao definido                                              * 
 * zComb -> fracoes de massa                                          * 
 * n     -> numero de linhas na matriz                                * 
 * nComb -> numero de componentes na mistura                          * 
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * z -> soma ( deve ser igual 1)                                      * 
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *------------------------------------------------------------------- *
 **********************************************************************/
void sumFracZ(DOUBLE *z       ,DOUBLE *zComb 
         ,INT const n     ,short const nComb)
{
  short j;
  INT i;

  for(i=0;i<n;i++)
  {
    z[i] = 0.e0;
    for(j=0;j<nComb;j++)
      z[i] += MAT2D(i,j,zComb,nComb);
  }
}
/**********************************************************************/



/*********************************************************************
 * Data de criacao    : 03/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * maxArray : retorna o valor maximo de vetor                        *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * x -> valores                                                      *
 * y -> dimensao de vetor                                            *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE maxArray(DOUBLE *RESTRICT x,INT const n)
{
  
  INT j;
  DOUBLE tmp;

  tmp = max(x[0], x[1]);
  for(j=2;j<n;j++)
    tmp = max(tmp, x[j]); 

  return tmp;

}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 24/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * edc : Eddy dissipation concept                                    *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * y        -> fracao massica das especies primitivas                *
 * s        -> taxa de consumo do oxidante                           *
 * desnity  -> densidade do fluido dentro do reator/celula           *
 * vol      -> volume do reator/celula                               *
 * eddyVisc -> viscosidae turbulenta                                 *
 * c        -> constantes necessaria                                 *
 * dVics    ->  viscosidae dinamica                                  *
 * df       -> coeficiente de dufusao do fuel                        *
 * tMix     -> tempo de mistura definido pe usuario                  *
 * iCod     -> metodo escolhido                                      *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * y(nel,0) -> Fuel(comburante)                                      * 
 * y(nel,1) -> O2  (oxidante)                                        * 
 * y(nel,2) -> CO2                                                   *
 * y(nel,3) -> H2O                                                   *
 * y(nel,4) -> N2                                                    *
 *********************************************************************/
DOUBLE edcOld(DOUBLE *y          ,short const iYf
          ,short const iYox      
          ,short *iProd       ,short const nProd  
          ,DOUBLE const s     ,DOUBLE const density
          ,DOUBLE const vol   ,DOUBLE const eddyVisc
          ,DOUBLE *c          ,DOUBLE const modS 
          ,DOUBLE const dVisc ,DOUBLE const df 
          ,DOUBLE const tMix  ,DOUBLE const tChemical
          ,short const iCod)
{
  short i;
  DOUBLE omega,r,k,delta,tm,tg,tu,td,tc,tf,itMix;
  DOUBLE x,x1,x2,x3,yF,yOx,yP,yMin,e;
  DOUBLE tmp1,tmp2,tmp3,gEdc,gamma;

/*...*/
  yF    = y[iYf];
  yOx   = y[iYox];
  for(i=0,yP=0.e0;i<nProd;i++)
    yP += y[iProd[i]];
/*..................................................................*/

/*... fds*/
  switch (iCod)
  {
/*...*/
    case FDS_EDC:
      delta = pow(vol,D1DIV3);
/*... estimativa a energia cinetica turbulenta*/
      tmp1 = c[2]*c[2]/0.094;
      k = (tmp1)*(tmp1)*delta*delta*modS;
/*..................................................................*/
      td = delta*delta/df;
      tu = 0.4*delta/sqrt((2.e0/3.e0)*k);
      tg = sqrt(2.e0*delta/9.81);
      tc = tChemical;
      tf = c[1];
      itMix = max(tc,min(tg,min(td,min(tu,tf))));
      itMix = 1.e0/itMix; 
      gEdc  = 1.e0;
      break;
/*..................................................................*/

/*...*/    
    case FLUENT_EDC:
    case FLUENT_CONST_TMIX_EDC:
      gEdc  = 1.e0;
      itMix = c[1]*modS;     
      break;
/*..................................................................*/

/*... Balram Panjwani - 2010*/
    case PANJWANI_EDC:
    case PANJWANI_CONST_TMIX_EDC:
/*...*/
      itMix = c[1]*modS;
/*... calculo */
      gamma = min(c[0]*pow(dVisc/(eddyVisc+dVisc),0.25),0.99);
      yMin = min(yF,yOx/s); 
      tmp1 = s + 1;
      tmp2 = yMin + yP/tmp1;
      tmp3 = yP/tmp1;
      x1 = (tmp2*tmp2)/(yF+tmp3)*(yOx+tmp3);
      x2 = min((tmp2)/(gamma*tmp2),1.e0);
      x3 = (gamma*tmp2/yMin,1.e0);
      x = x1*x2*x3;
      gEdc = x/(1.e0 - x*gamma);
      break;
/*...................................................................*/
    default:
      ERRO_OP_NEW(__FILE__,__func__,__LINE__,"Eddy Dissipation model"
                                          ,iCod)
      break;
  }

/*... tempo de mixutura constante e definida pelo usuario*/
  if(FLUENT_CONST_TMIX_EDC==iCod 
  || PANJWANI_CONST_TMIX_EDC==iCod)  itMix = 1.e0/tMix;
/*....................................................................*/

  r = min(yF,yOx/s);
  omega = density*max(r,0.e0)*itMix*gEdc;
  return omega;
}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 24/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * edc : Eddy dissipation concept                                    *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * y        -> fracao massica das especies primitivas                *
 * s        -> taxa de consumo do oxidante                           *
 * desnity  -> densidade do fluido dentro do reator/celula           *
 * vol      -> volume do reator/celula                               *
 * eddyVisc -> viscosidae turbulenta                                 *
 * c        -> constantes necessaria                                 *
 * dVics    ->  viscosidae dinamica                                  *
 * df       -> coeficiente de dufusao do fuel                        *
 * tMix     -> tempo de mistura definido pe usuario                  *
 * iCod     -> metodo escolhido                                      *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * y(nel,0) -> Fuel(comburante)                                      * 
 * y(nel,1) -> O2  (oxidante)                                        * 
 * y(nel,2) -> CO2                                                   *
 * y(nel,3) -> H2O                                                   *
 * y(nel,4) -> N2                                                    *
 *********************************************************************/
DOUBLE edc(DOUBLE *y          ,short const iYf
          ,short const iYox      
          ,short *iProd       ,short const nProd  
          ,DOUBLE const s     ,DOUBLE const density
          ,DOUBLE const vol   ,DOUBLE const eddyVisc
          ,DOUBLE *c          ,DOUBLE const modS 
          ,DOUBLE const dVisc ,DOUBLE const df 
          ,DOUBLE const tMix  ,DOUBLE const tChemical
          ,short const iCod)
{
  short i;
  DOUBLE omega,r,k,delta,tm,tg,tu,td,tc,tf,itMix;
  DOUBLE x,x1,x2,x3,yF,yOx,yP,yMin,e;
  DOUBLE tmp1,tmp2,tmp3,gEdc,gamma;

/*...*/
  yF    = y[iYf];
  yOx   = y[iYox];
  for(i=0,yP=0.e0;i<nProd;i++)
    yP += y[iProd[i]];
/*..................................................................*/

/*... fds*/
  switch (iCod)
  {
/*...*/
    case FDS_EDC:
      delta = pow(vol,D1DIV3);
/*... estimativa a energia cinetica turbulenta*/
      tmp1 = c[2]*c[2]/0.094;
      k = (tmp1)*(tmp1)*delta*delta*modS;
/*..................................................................*/
      td = delta*delta/df;
      tu = 0.4*delta/sqrt((2.e0/3.e0)*k);
      tg = sqrt(2.e0*delta/9.81);
      tc = tChemical;
      tf = c[1];
      itMix = max(tc,min(tg,min(td,min(tu,tf))));
      itMix = 1.e0/itMix; 
      gEdc  = 1.e0;
      break;
/*..................................................................*/

/*...*/    
    case FLUENT_EDC:
    case FLUENT_CONST_TMIX_EDC:
      gEdc  = 1.e0;
      itMix = c[1]*modS;     
      break;
/*..................................................................*/

/*... Balram Panjwani - 2010*/
    case PANJWANI_EDC:
    case PANJWANI_CONST_TMIX_EDC:
/*...*/
      itMix = c[1]*modS;
/*... calculo */
      gamma = min(c[0]*pow(dVisc/(eddyVisc+dVisc),0.25),0.99);
      yMin = min(yF,yOx/s); 
      tmp1 = s + 1;
      tmp2 = yMin + yP/tmp1;
      tmp3 = yP/tmp1;
      x1 = (tmp2*tmp2)/(yF+tmp3)*(yOx+tmp3);
      x2 = min((tmp2)/(gamma*tmp2),1.e0);
      x3 = (gamma*tmp2/yMin,1.e0);
      x = x1*x2*x3;
      gEdc = x/(1.e0 - x*gamma);
      break;
/*...................................................................*/
    default:
      ERRO_OP_NEW(__FILE__,__func__,__LINE__,"Eddy Dissipation model"
                                          ,iCod)
      break;
  }

/*... tempo de mixutura constante e definida pelo usuario*/
  if(FLUENT_CONST_TMIX_EDC==iCod 
  || PANJWANI_CONST_TMIX_EDC==iCod)  itMix = 1.e0/tMix;
/*....................................................................*/

  r = min(yF,yOx/s);
  omega = density*max(r,0.e0)*itMix*gEdc;
  return omega;
}
/*********************************************************************/ 


/*********************************************************************
 * Data de criacao    : 26/07/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * edm : Eddy dissipation model                                      *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * y        -> fracao massica das especies primitivas                *
 * w        -> nao de definido                                       *
 * desnity  -> densidade do fluido dentro do reator/celula           *
 * tMix     -> tempo de mistura definido pe usuario                  *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * w        -> taxa de racao de massa                                *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void edm(Combustion *c       ,DOUBLE *RESTRICT y  
          ,DOUBLE *RESTRICT w  ,DOUBLE const density
          ,DOUBLE const itMix)
{
  short i,j,id,
        nSp  = c->nOfSpecies;
          
  DOUBLE omega,dy,r1,r2,my,A,B,Q[MAXREAC],itM=itMix;

/*...*/
  A = c->edm.coef[0];
  B = c->edm.coef[1];
  if(c->edm.tMixConst)
    itM = 1.e0/c->edm.coef[2];
/*....................................................................*/

/*... reagentes*/  
  for(i=0;i<c->chem.nReac;i++)
  {
    r1 = 1.e+32;
    for(j=0;j<c->chem.reac[i].nPartSp[0];j++)
    {
      id = c->chem.reac[i].partSp[0][j];
      dy = y[id]/(c->chem.reac[i].stch[0][id]*c->chem.sp[id].mW);
      r1  = min(r1,dy);
    }
/*....................................................................*/

/*... produtos*/  
    r2 = 1.e+32;
    if(c->edm.fProd){
      for(j=0,dy=my=0.e0;j<c->chem.reac[i].nPartSp[1];j++)
      {
        id = c->chem.reac[i].partSp[1][j];
        dy += y[id];
        my += c->chem.reac[i].stch[1][id]*c->chem.sp[id].mW;
      }
      r2 = B*dy/my;
    }
/*....................................................................*/

/*...*/
    Q[i] = density*A*itM*min(r1,r2);
/*....................................................................*/
  }
/*....................................................................*/

  massRateReaction(&c->chem,Q,w);

}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 27/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * globalReac:                                                       *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *                                                                   *
 * 1*(CmHnOl) + nO2*(O2) + (nO2*pN2/pO2)*N2                          * 
 * -> nCO2*(CO2) + nH2O*(H2O) + nN2*(N2)                             *                              *
 *                                                                   *
 * C - m = nCO2   :.                                                 *
 * H - n = 2*nH2O :. nH2O = n/2                                      *
 * O - 2*nO2 +l = 2*nCO2 + nH2O :. nO2 = nCO2 + nH2O/2               *
 * :. nO2 = m + n/4 - l/2                                            *
 * N - 2*nO2*pN2/pO2 = 2*nN2 :. nO2*pN2/pO2 = nN2 :.                 *
 * :.  nN2 = (m + n/4)*pN2/pO2                                       *
 *                                                                   *
 * nAir = nO2/pCO2                                                   *
 *********************************************************************/
void globalReac(Combustion *c, short const iReac)
{

  short m,n,l,posN2=c->chem.sN2;
  DOUBLE nAir,nProd,pO2,pN2,pN2r,nN2p,nCO2p,nH2Op,nCOp,nTotal;
  DOUBLE pCO2r,pCOr,pH2Or,nO2,nN2;

//m = c->fuel[iReac].c;
//n = c->fuel[iReac].h;
//l = c->fuel[iReac].o;

//pO2 = c->O2InAir;
//pN2 = c->N2InAir;


/*... Ar ( pO2*O2 + pN2*N2)*/
//nO2  = m + n/4.0e0 - l/2.0;
//nN2  = nO2*pN2/pO2;
//nAir = nO2*(1.0e0/pO2);

/*... produto  */
//nN2p  = (m + n/4.e0)*pN2/pO2;
//nCO2p = m;
//nH2Op = n/2.e0;  

//nTotal = nCO2p + nH2Op + nN2p;
//pCO2r = nCO2p/nTotal;
//pH2Or = nH2Op/nTotal;
//pN2r = nN2p/nTotal;
//nProd = (nN2p/pN2r+nCO2p/pCO2r+nH2Op/pH2Or)/3.0;

//c->stoich[iReac][0][c->sp_fuel[iReac]]= 1.e0;
//c->stoich[iReac][0][c->sp_O2]  = nO2;
//c->stoich[iReac][0][c->sp_N2]  = nN2;
//c->stoich[iReac][1][c->sp_H2O] = nH2Op;
//c->stoich[iReac][1][c->sp_CO2] = nCO2p;
//c->stoich[iReac][1][c->sp_N2]  = nN2p;
//c->stoichAir                   = nAir;

//c->CO2InProd = pCO2r;
//c->H2OInProd = pH2Or;
//c->N2InProd  = pN2r;

/*.*/
//c->nSpeciesPart[iReac][0] = 2;
//c->nSpeciesPart[iReac][1] = 2;
/*. reagente*/
//c->speciesPart[iReac][0][0] = c->sp_fuel[iReac];
//c->speciesPart[iReac][0][1] = c->sp_O2;
//c->speciesPart[iReac][0][2] = c->sp_N2;
/*. produto*/
//c->speciesPart[iReac][1][0] = c->sp_CO2;
//c->speciesPart[iReac][1][1] = c->sp_H2O;
//c->speciesPart[iReac][1][2] = c->sp_N2;

/*.*/
//c->sMassAir             = nAir*c->mW_Air/c->mW[c->sp_fuel[iReac]];
//c->sMass[iReac][0][c->sp_fuel[iReac]] = 1.0;
//c->sMass[iReac][0][c->sp_O2]   = nO2*c->mW[c->sp_O2]/c->mW[c->sp_fuel[iReac]];
//c->sMass[iReac][0][c->sp_N2]   = nN2*c->mW[c->sp_N2]/c->mW[c->sp_fuel[iReac]];
//c->sMass[iReac][1][c->sp_CO2]  = nCO2p*c->mW[c->sp_CO2]/c->mW[c->sp_fuel[iReac]];
//c->sMass[iReac][1][c->sp_H2O]  = nH2Op*c->mW[c->sp_H2O]/c->mW[c->sp_fuel[iReac]];
//c->sMass[iReac][1][c->sp_N2]   = nN2p*c->mW[c->sp_N2]/c->mW[c->sp_fuel[iReac]];

/* fprintf(fileLogExc,"\nReaction:\n\n"); 
/*...*/
//if(c->fLump)
//{
//  fprintf(fileLogExc,"C%dH%dO%d" 
//         " + %lf (%lf O2 + %lf N2)\n"
//         " -> %lf (%lf CO2 + %lf H2O + %lf N2)\n\n"
//        ,m    ,n     ,l
//        ,nAir ,pO2   ,pN2
//        ,nProd,pCO2r, pH2Or,pN2r);
    
//  fprintf(fileLogExc,"1 kg C%dH%dO%d" 
//         " + %lf (%lf O2 + %lf N2)\n"
//         " -> %lf (%lf CO2 + %lf H2O + %lf N2)\n\n"
//        ,m                ,n     ,l
//        ,c->sMassAir    ,pO2   ,pN2
//        ,c->sMassAir + 1,pCO2r,pH2Or,pN2r);
//} 
/*..................................................................*/

/*...*/ 
//else
//{
//  fprintf(fileLogExc,"C%dH%dO%d" 
//         " + %lf O2 + %lf N2\n"
 //        " -> %lf CO2 + %lf H2O + %lf N2\n\n"
//        ,m    ,n    ,l
//        ,c->stoich[iReac][0][c->sp_O2]
//        ,c->stoich[iReac][0][c->sp_N2]
//        ,c->stoich[iReac][1][c->sp_CO2]
//        ,c->stoich[iReac][1][c->sp_H2O]
//        ,c->stoich[iReac][1][c->sp_N2]);

//  fprintf(fileLogExc,"1 kg C%dH%dO%d" 
//         " + %lf  kg O2 + %lf  kg N2\n"
//         " -> %lf kg CO2 + %lf kg H2O + %lf kg N2\n\n"
//        ,m    ,n    ,l
//        ,c->sMass[iReac][0][c->sp_O2]  
//        ,c->sMass[iReac][0][c->sp_N2]
//        ,c->sMass[iReac][1][c->sp_CO2]
//        ,c->sMass[iReac][1][c->sp_H2O]
//        ,c->sMass[iReac][1][c->sp_N2]);  
//}
/*..................................................................*/

}
/*******************************************************************/

/*********************************************************************
* Data de criacao    : 02/06/2019                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* getVolumeMed : Calula a media volumetria de um gradeza            *
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
DOUBLE getVolumeMed(DOUBLE *RESTRICT x,DOUBLE *RESTRICT vol
                   ,INT const n)
{
  INT i;
  DOUBLE sum,volT;
 
  for (i = 0,sum=volT=0.e0; i < n; i++)
  {
    sum+= x[i]*vol[i];
    volT+=vol[i];
  }
  return sum/volT;
}
/*******************************************************************/



void printt(double *x, int n)
{
  int i;
  for(i=0;i<n;i++)
    if( x[i] != 0.e0) printf("%d %e\n",i,x[i]);
}