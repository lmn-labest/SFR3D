#include<Combustion.h>

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
  rateFuelConsume(cModel                   , tModel
                  , mesh->elm.zComb        , mesh->elm.cDiffComb
                  , mesh->elm.temp         , mesh->elm.rateFuel 
                  , mesh->elm.densityFluid , mesh->elm.gradVel
                  , mesh->elm.eddyViscosity,mesh->elm.dViscosity
                  , mesh->elm.geom.volume
                  , mesh->ndm              , mesh->numel
                  , eModel->fKelvin );    
/*...................................................................*/

/*... reconstruindo do gradiente (Fuel)*/
  tm.rcGradComb   = getTimeC() - tm.rcGradComb;
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
             , mesh->elm.rateFuel      , mesh->elm.vel           
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
//regularZ(mesh->elm.zComb,mesh->numelNov
//        ,cModel->nComb  ,cModel->fLump);
/*...................................................................*/

/*...*/
  getSpeciesPrimitives(cModel
                      ,mesh->elm.yFrac,mesh->elm.zComb
                      ,mesh->numelNov);
/*...................................................................*/

/*...*/
  rateHeatRealeseCombustion(cModel            , &prop->sHeat                
                    , mesh->elm.rateHeatReComb, mesh->elm.temp     
                    , mesh->elm.zComb0        , mesh->elm.zComb
                    , mesh->elm.densityFluid  , mesh->elm.rateFuel 
                    , mesh->elm.material.prop , mesh->elm.mat    
                    , sc->ddt.dt[TIME_N]      , mesh->numelNov
                    , fSheat                  , eModel->fKelvin );  
/*...................................................................*/

/*...*/
  prop->molarMass = mixtureMolarMass(cModel,mesh->elm.yFrac);
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 15/08/2018                                   *
 * Data de modificaco : 15/05/2019                                   *
 *-------------------------------------------------------------------*
 * regularZ : mantem o z entre 0 e 1                                 *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * z            -> fracao massica                                    *
 * numel        -> numero de elementos                               *
 * nOfSpecies   -> numero de especie agrupadas                       *
 * fLump        -> true  - especies agrupadas                        *
 *                 false - especies primitivas                       *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * z -> fracao de massa de regularizadas                             * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * z < 0.0 -> z = 0.e0                                               *
 *********************************************************************/
void regularZ(DOUBLE *RESTRICT z    , INT const numel
            , short const nComb      , bool fLump)
{
  short j;
  INT nel;
/*DOUBLE zSum;*/

/*... z < 0.0*/
  for(nel = 0 ; nel < numel; nel++)
  {
    for(j=0;j<nComb;j++)
      if(MAT2D(nel,j,z,nComb) < 0.e0)
        MAT2D(nel,j,z,nComb) = 0.e0;
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
;
        for(i=0,sum=0.e0;i<ns-1;i++)
          sum+=MAT2D(nel,i,y,ns) ;

        MAT2D(nel,cModel->sp_N2,y,ns) = 1.e0 - sum;

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

      for(i=0,sum=0.e0;i<ns-1;i++)
        sum+=y[i];
      y[cModel->sp_N2] = 1.e0 - sum;

    }
  }
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 24/05/2019                                   *
 *-------------------------------------------------------------------*
 * rateFuelConsume: calculo da taxa de consumo do combustivel        *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> modelo de combustao                                    *
 * zComb   -> fracao massica                                         *
 * temp    -> temperatura                                            *
 * rate    -> nao definido                                           * 
 * density -> densidade                                              *
 * numel   -> numero de elementos                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * rate    -> taxa de consumo do combustivel                         * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * fLump - true                                                      *
 * z(nel,0) -> fuel(comburente)                                      *
 * z(nel,1) -> air (oxidante)                                        *
 * z(nel,2) -> produto                                               *
 * z(nel,0) -> fuel(comburente)                                      *
 * fLump - false                                                     *
 * z(nel,1) -> O2  (oxidante)                                        * 
 * z(nel,2) -> CO2                                                   *
 * z(nel,3) -> H2O                                                   *
 * z(nel,4) -> N2                                                    *
 *********************************************************************/
void rateFuelConsume(Combustion *cModel      , Turbulence *tModel
             , DOUBLE *RESTRICT zComb        , DOUBLE *RESTRICT diffComb
             , DOUBLE *RESTRICT temp         , DOUBLE *RESTRICT rate
             , DOUBLE *RESTRICT density      , DOUBLE *RESTRICT gradVel
             , DOUBLE *RESTRICT eddyViscosity, DOUBLE *RESTRICT dViscosity
             , DOUBLE *RESTRICT volume
             , short const ndm               , INT const numel
             , bool const fKelvin )
{
  short nComb = cModel->nComb
       , iCod = cModel->reactionKinetic
       , nReac=cModel->nReac;
  short iComb,iOx,iProd[3],i,j;
  INT nel;
  DOUBLE s,tMix,eddy,sT[6],*iGradVel,df;
  DOUBLE zFuel, zAir, zO2,zOx, zProp, omega, densityC, alpha, coefA;
  DOUBLE tmp1, tmp2, tmp3, tmp4, modS, c[3], e1, e2;
  DOUBLE mWfuel,mWox, tempA, ru, tc, y[MAXSPECIES]; 

/*...*/
  switch(iCod)
  {
/*...*/
    case ARRHENIUS:
/*... cal/(mol*kelvin)*/
      ru = IDEALGASR*2.39006e-04;
 /*J/(mol.kelvin) */
//    ru = IDEALGASR*1.e-03;
      for(i=0;i<cModel->nReac;i++)
      {
        iComb  = cModel->sp_fuel[i];
        alpha  = cModel->arrhenius[i].alpha;
        mWfuel = cModel->mW[iComb];

        tempA  = cModel->arrhenius[i].energyAtivation/ru;
        coefA  = cModel->arrhenius[i].a;
        e1     = cModel->arrhenius[i].e1; 
        e2     = cModel->arrhenius[i].e2; 
/*...*/
        if(cModel->fLump)
        {
          s      = cModel->stoichAir;
          mWox   = cModel->mW_Air;
        }
        else
        {
          s      = cModel->stoich[i][0][cModel->sp_O2];
          mWox   = cModel->mW[cModel->sp_O2];
        }
/*...................................................................*/

/*...*/
        for(nel = 0; nel < numel; nel++)
        {

          densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
          getSpeciesPrimitivesCc(cModel,y,zComb);
          omega = arrhenius(y[iComb]  ,y[cModel->sp_O2]
                                ,e1       ,e2   
                                ,mWfuel   ,mWox
                               ,temp[nel] ,alpha
                               ,densityC  ,tempA 
                               ,coefA     ,fKelvin);

          MAT2D(nel,i,rate,nReac) = mWfuel*omega;
        }
/*...................................................................*/ 
    } 
/*...................................................................*/    
    break;
/*...................................................................*/

/*...*/
    case EDC:
      tMix = cModel->edc.tMix;
      c[0] = cModel->edc.cGamma;
      c[1] = cModel->edc.cTau;
      c[2] = tModel->cf;
      for(i=0;i<cModel->nReac;i++)
      {
        if(cModel->fLump) 
          s = cModel->sMassAir; 
        else
        {
          iComb   = cModel->sp_fuel[i];
          iOx     = cModel->sp_O2;  
          s       = cModel->sMass[i][0][cModel->sp_O2]; 
          for(j=0;j<cModel->nSpeciesPart[i][1];j++)
            iProd[j] = cModel->speciesPart[i][1][j];
        }
/*...*/
        for(nel = 0; nel < numel; nel++)
        {
          densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
          getSpeciesPrimitivesCc(cModel,y,zComb);
/*.. calculo Sij*/
          iGradVel = &MAT3D(nel,0,0,gradVel,ndm,ndm);
          tensorS(sT,iGradVel,false);
/*... |S| = sqrt(2S:S)*/
          modS = sqrt(2.e0*doubleDotSym(sT));
/*...................................................................*/
          densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);
          df    = MAT2D(nel,cModel->sp_fuel[i],diffComb,nComb);
          omega  = edc(y              ,iComb
                      ,iOx            ,iProd 
                      ,cModel->nSpeciesPart[i][1]
                      ,s              ,densityC
                      ,volume[nel]    ,eddyViscosity[nel]
                      ,c              ,modS
                      ,dViscosity[nel],df
                      ,tMix           ,cModel->edc.type);
          MAT2D(nel,i,rate,nReac) = omega;
        }
      }
/*...................................................................*/
      break;
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 07/05/2019                                   *
 *-------------------------------------------------------------------*
 * rateHeatRealeseCombustion: calculo da taxa de liberacao de calor  *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> modelo de combustao                                    *
 * sHeatPol - estrutra par o polinoimio do calor especifico          *
 * q       -> nao definido                                           *
 * zComb0  -> fracao de massa agrupada do passo de tempo anterior    *
 * zComb   ->fracao de massa agrupada do passo de tempo atural       *
 * rateFuel-> taxa de consumo do combustivel                         *
 * prop   - propriedades por material                                *
 * mat    - material da celula                                       * 
 * dt      -> delta dessa passo de tempo                             * 
 * numel   -> numero de elementos                                    *
 * fSheat   - calor especifico com variacao com a Temperatura        *
 * fKelvin  - temperatura dada em kelvin                             *
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
void rateHeatRealeseCombustion(Combustion *cModel,PropPol *sHeat   
                   , DOUBLE *RESTRICT q      , DOUBLE *RESTRICT temp
                   , DOUBLE *RESTRICT zComb0 , DOUBLE *RESTRICT zComb
                   , DOUBLE *RESTRICT density, DOUBLE *RESTRICT rateFuel 
                   , DOUBLE *RESTRICT prop   , short  *RESTRICT mat
                   , DOUBLE const dt         , INT const numel
                   , bool const fsHeat       , bool const fKelvin)
{

  short lMat, iCod = cModel->typeHeatRealese, nReac = cModel->nReac,iComb,i, j, kSp;
  INT nel;
  DOUBLE *h,hc,H,DH,HP,HR,hs,nSp;
  DOUBLE sHeatRef,sum;

/*... Entalpia de formacao*/
  if(iCod == HFORMATION)
  {
    h     = cModel->entalphyOfForm;
    for(nel = 0; nel < numel; nel++)
    {
      lMat  = mat[nel] - 1;
/*... reacao i*/
      for(i=0,sum=0;i<nReac;i++)
      {   
        iComb = cModel->sp_fuel[i];
/*... reagentes*/
        for(j=0,HR=0.e0;j<cModel->nSpeciesPart[i][0];j++)
        {
          kSp = cModel->speciesPart[i][0][j];
          nSp = cModel->stoich[i][0][kSp];

          sHeatRef = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);

          hs = tempForSpecificEnthalpySpecies(sHeat    , kSp 
                                                     , temp[nel], sHeatRef
                                                     , fsHeat   , fKelvin);
          H = h[kSp]/* + hs*/;
          HR += nSp*H;
         
        } 
/*... reagentes*/
        for(j=0,HP=0.e0;j<cModel->nSpeciesPart[i][1];j++)
        {
          kSp = cModel->speciesPart[i][1][j];
          nSp = cModel->stoich[i][1][kSp];

          sHeatRef = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);

          hs =  tempForSpecificEnthalpySpecies(sHeat    , kSp 
                                             , temp[nel], sHeatRef
                                             , fsHeat   , fKelvin);
          H = h[kSp]/* + hs*/;
          HP += nSp*H;
         
        } 
        sum += (HP-HR)*MAT2D(nel,i,rateFuel,nReac)/cModel->mW[iComb];      
      }
/*...................................................................*/
      q[nel] = -sum; 
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... Entalpia de combustao*/
  else if (iCod == HCOMBUSTION)
  {
    hc = cModel->entalphyOfCombustion;
    for(nel = 0; nel < numel; nel++)
    {
      q[nel] = -rateFuel[nel]*hc;
    }
  }
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * rateFuelConsume: calculo da taxa de consumo do combustivel        *
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

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 12/05/2019                                   *
 *-------------------------------------------------------------------*
 * initLumpedMatrix : inicializa a matriz de especies                *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> modelo de combustao                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void initLumpedMatrix(Combustion *cModel)
{
  
  short nl = cModel->nOfSpeciesLump
      , ns =  cModel->nOfSpecies 
      , posN2 = cModel->sp_N2;
  DOUBLE mO2,mN2,mCO2p,mH2Op,mN2p,mAir,mProd;

  mO2 = mN2 = mCO2p = mH2Op = mN2p = 0.e0; 

  mO2   = cModel->stoich[0][0][cModel->sp_O2]*cModel->mW[cModel->sp_O2];
  mN2   = cModel->stoich[0][0][posN2]*cModel->mW[posN2];
/*...*/
  mCO2p = cModel->stoich[0][1][cModel->sp_CO2]*cModel->mW[cModel->sp_CO2];
  mH2Op = cModel->stoich[0][1][cModel->sp_H2O]*cModel->mW[cModel->sp_H2O];
  mN2p  = cModel->stoich[0][1][posN2]*cModel->mW[posN2];

  mAir  = mO2 + mN2;
  mProd = mCO2p + mH2Op + mN2p;
/*... Fuel*/
  MAT2D(cModel->sp_fuel[0],0,cModel->lumpedMatrix,nl) = 1.e0; 
  MAT2D(cModel->sp_fuel[0],1,cModel->lumpedMatrix,nl) = 0.e0; 
  MAT2D(cModel->sp_fuel[0],2,cModel->lumpedMatrix,nl) = 0.e0; 
/*...................................................................*/

/*... O2*/
  MAT2D(cModel->sp_O2,0,cModel->lumpedMatrix,nl) = 0.e0;
  MAT2D(cModel->sp_O2,1,cModel->lumpedMatrix,nl) = mO2/mAir; 
  MAT2D(cModel->sp_O2,2,cModel->lumpedMatrix,nl) = 0.0e0; 
/*...................................................................*/

/*... CO2*/
  MAT2D(cModel->sp_CO2,0,cModel->lumpedMatrix,nl) = 0.e0; 
  MAT2D(cModel->sp_CO2,1,cModel->lumpedMatrix,nl) = 0.0e0; 
  MAT2D(cModel->sp_CO2,2,cModel->lumpedMatrix,nl) = mCO2p/mProd;  
/*...................................................................*/

/*... H2O*/
  MAT2D(cModel->sp_H2O,0,cModel->lumpedMatrix,nl) = 0.e0; 
  MAT2D(cModel->sp_H2O,1,cModel->lumpedMatrix,nl) = 0.e0; 
  MAT2D(cModel->sp_H2O,2,cModel->lumpedMatrix,nl) =  mH2Op/mProd; 
/*...................................................................*/

/*... N2*/
  MAT2D(cModel->sp_N2,0,cModel->lumpedMatrix,nl) = 0.e0; 
  MAT2D(cModel->sp_N2,1,cModel->lumpedMatrix,nl) = mN2/mAir; 
  MAT2D(cModel->sp_N2,2,cModel->lumpedMatrix,nl) = mN2p/mProd;  
/*...................................................................*/


}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * initLumpedMatrix : inicializa a matriz de especies                *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * y       -> nao definido                                           *
 * a       -> matriz das especies agrupadas                          *
 * z       -> especies agrupadas                                     *
 * ns      -> numero de especies primitivas                          *
 * nl      -> numero de especies agrupadas                           *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * y       -> especies primitivas                                    *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * c[0] - CH4                                                        *
 * c[1] - O2                                                         *
 * c[2] - CO2                                                        *
 * c[3] - H2O                                                        *
 * c[4] - N2                                                         *
 *********************************************************************/
void yLumpedMatrixZ(DOUBLE *RESTRICT y, DOUBLE *RESTRICT a
                  , DOUBLE *RESTRICT z
                  , short const ns    , short const nl)
{
  short i,j;

/*...*/
  for(i=0;i<ns;i++){
    y[i] = 0.e0;
    for(j=0;j<nl;j++)
      y[i] += MAT2D(i,j,a,nl)*z[j];  
  }
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * initMolarMass: inicializa a massa molar das especies              *
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
void initMolarMass(Combustion *cModel)
{
  
  short nO,nC,nH,i;
  DOUBLE pO2,pN2;

  pO2 = cModel->O2InAir;
  pN2 = cModel->N2InAir;

/*... Fuel */
  for(i=0;i<cModel->nReac;i++)
  {
    nC = cModel->fuel[i].c;
    nO = cModel->fuel[i].o;
    nH = cModel->fuel[i].h;
    cModel->mW[cModel->sp_fuel[i]] = nC*MW_C + nH*MW_H + nO*MW_O;
  }
/*...O2*/
  cModel->mW[cModel->sp_O2] =  2.0e0*MW_O;
/*...H2O*/
  cModel->mW[cModel->sp_H2O]  = 2.0e0*MW_H + MW_O;
/*...CO2*/
  cModel->mW[cModel->sp_CO2] = MW_C +  2.0e0*MW_O;
/*...N2*/
  cModel->mW[cModel->sp_N2] =  2.0e0*MW_N;
/*... Air*/
  cModel->mW_Air = pO2* cModel->mW[cModel->sp_O2] 
                  + pN2*cModel->mW[cModel->sp_N2];

}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 23/08/2018                                   *
 * Data de modificaco : 21/05/2019                                   *
 *-------------------------------------------------------------------*
 * stoichiometricCoeff:                                              *
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
void stoichiometricCoeff(Combustion *cModel)
{
  short i;
  for(i=0;i<cModel->nReac;i++)
  {
    fprintf(fileLogExc,"%d)\n",i);
    globalReac(cModel, i);
  }
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/08/2018                                   *
 * Data de modificaco : 27/05/2019                                   *
 *-------------------------------------------------------------------*
 * initMolarMass: inicializa a massa molar das especies              *
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
void initEntalpyOfFormation(Combustion *cModel)
{
  short i;
/*...  Fuel - KJ/kMol*/
  for(i=0;i<cModel->nReac;i++)
    cModel->entalphyOfForm[cModel->sp_fuel[i]] = cModel->fuel[i].hf; 

/*... CO2 - KJ/kMol*/
  cModel-> entalphyOfForm[cModel->sp_CO2] = -393510.e0;
/*... H2O - KJ/kMol*/
  cModel-> entalphyOfForm[cModel->sp_H2O] = -241826.e0;
/*... O2 - KJ/kMol*/
  cModel-> entalphyOfForm[cModel->sp_O2] = 0.e0;
/*... N2 - KJ/kMol*/
  cModel-> entalphyOfForm[cModel->sp_N2] = 0.e0;
    
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/08/2018                                   *
 * Data de modificaco : 05/05/2019                                   *
 *-------------------------------------------------------------------*
 * initMolarMass: inicializa a massa molar das especies              *
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
void initEntalpyOfCombustion(Combustion *cModel)
{
  
  DOUBLE pCO2,pH2O,hFuel,hCO2,hH2O;

  pCO2 = cModel->CO2InProd;
  pH2O = cModel->H2OInProd;

/*...  Fuel - KJ/kMol*/
  hFuel = -74870.e0; 
/*... CO2 - KJ/kMol*/
  hCO2 = -393520.e0;
/*... H2O - KJ/kMol*/
  hH2O = -241830.e0;

/*... KJ/Kmol de fuel*/
  cModel->entalphyOfCombustion = (hCO2 + 2.0*hH2O) - hFuel;

  fprintf(fileLogExc,"Primitives species:\n");

  fprintf(fileLogExc,"Entalphy of combustion (KJ/Kmol) = %lf\n"  
                    ,cModel->entalphyOfCombustion);


/*... KJ/KG de fuel */
  cModel->entalphyOfCombustion  /= cModel->mW[cModel->sp_fuel[0]];

  fprintf(fileLogExc,"Entalphy of combustion (KJ/KG)   = %lf\n\n"  
                    ,cModel->entalphyOfCombustion);
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 02/05/2019                                   *
 *-------------------------------------------------------------------*
 * concetracionOfSpecies:                                            *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> nao definido                                           *
 * z       -> especies agrupadas   (cModel->fLump = true)            *
 *            especies primitivas  (cModel->fLump = false)           *
 * c       -> concetracao das especies primitivas                    *
 * density -> densidade da mistura                                   *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * c       -> concentracao das especies                              *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * c[0]   - CH4                                                      *
 * c[1]   - O2                                                       *
 * c[2]   - CO2                                                      *
 * c[3]   - H2O                                                      *
 * c[4]   - CO                                                       *
 *  ...                                                              *
 * c[N-1] - N2                                                       *
 *********************************************************************/
void concetracionOfSpecies(Combustion *cModel,DOUBLE *RESTRICT z
                          ,DOUBLE *RESTRICT c,DOUBLE const density)
{
  bool fLump = cModel->fLump;
  short i, ns = cModel->nOfSpecies, nl = cModel->nOfSpeciesLump;
  DOUBLE y[MAXSPECIES],*mW;

  mW = cModel->mW;

  if (fLump) yLumpedMatrixZ(y,cModel->lumpedMatrix, z, ns, nl);
  
  for(i=0;i<ns;i++)
    c[i] = y[i]*density/mW[i];

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
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 25/08/2018                                   *
 *-------------------------------------------------------------------*
 * mixtureMolarMass: massa molar da mistura                          *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> nao definido                                           *
 * y       -> especies primitivas                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE mixtureMolarMass(Combustion *cModel,DOUBLE *RESTRICT y) 
{                    

  short i, ns = cModel->nOfSpecies;
  DOUBLE *mW,tmp;

  mW = cModel->mW;

  for(i=0,tmp=0.e0;i<ns;i++)
    tmp += y[i]/mW[i];

  return 1.e0/tmp;

}
/*********************************************************************/

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
DOUBLE edc(DOUBLE *y           ,short const iYf
          ,short const iYox   ,short *iProd
          ,short const nProd
          ,DOUBLE const s     ,DOUBLE const density
          ,DOUBLE const vol   ,DOUBLE const eddyVisc
          ,DOUBLE *c          ,DOUBLE const modS 
          ,DOUBLE const dVisc ,DOUBLE const df 
          ,DOUBLE const tMix  ,short const iCod)
{
  short i;
  DOUBLE omega,r,k,delta,tm,tg,tu,td,tc,tf,itMix;
  DOUBLE x,x1,x2,x3,yF,yOx,yP,yMin,e;
  DOUBLE tmp1,tmp2,tmp3,gEdc,gamma;

/*...*/
  yF    = y[iYf];
  yOx   = y[iYox];
  for(i=0,yP=0;i<nProd;i++)
    yP += y[i];
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
      tc = 1.e-04;
      tf = 0.125;
      itMix =1.e0/max(tc,min(tg,min(td,min(tu,tf))));
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
      gamma = min(c[0]*pow(dVisc/(eddyVisc+dVisc),0.25),0.95);
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
 * Data de criacao    : 26/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * edc : Eddy dissipation concept                                    *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * y1       -> fracao massica                                        *
 * e1       -> expoente                                              *
 * mW1      -> massa molar                                           *
 * y2       -> fracao massica                                        *
 * e2       ->  expoente                                             *
 * mW2      -> massa molar                                           *
 * t        -> temperatura                                           *
 * alpha    -> coeficiente da temperatura                            *
 * desnity  -> densidade do fluido                                   *
 * tA       -> temperatura de ativacao                               *
 * coefA    -> coeficiente                                           *
 * fKelvin  ->                                                       *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * omega :kmol/m^3s                                                  *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * Reacao quimica a1A1 + a2A2 => a3A3 + a4A4                         * 
 * Unidade de Coef                                                   *
 * (N/L^3)^(1-(e1+e2)*(1/theta^alpha)(1/T)                           *
 * Exemplo:                                                          *
 * N     = mol                                                       *
 * L     = cm                                                        *
 * theta = Kelvin                                                    *
 * T     = segundos                                                  *
 * A = ((mol/cm^3)^(1-(e1+e2))/((K^alpha)(1/s))                      *
 *********************************************************************/
DOUBLE arrhenius(DOUBLE const y1     ,DOUBLE const y2
                ,DOUBLE const e1     ,DOUBLE const e2
                ,DOUBLE const mW1    ,DOUBLE const mW2
                ,DOUBLE const t      ,DOUBLE const alpha
                ,DOUBLE const density,DOUBLE const tA    
                ,DOUBLE const coefA  ,bool const fKelvin)
{
  short iCod=3;
  DOUBLE tc;
  DOUBLE omega,al,c1,c2,k,prodC,d;


  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);
/*... A = ((kmol/m^3)^(1-(e1+e2))/((K^alpha)(1/s))  */
  if(iCod == 1)
  {
    c1 = density*y1/mW1;
    c2 = density*y2/mW2;
    d  = 1.e0;
  }
/*..................................................................*/

/*... mol/m3*/
/*... A = ((mol/m^3)^(1-(e1+e2))/((K^alpha)(1/s))  */
  else if(iCod == 2)
  {
/*... kmol/m3 - > mol/m3*/
    c1 = 1.e+03*density*y1/mW1;
    c2 = 1.e+03*density*y2/mW2;
/*... mol/m^3s -> kmol/m^3s*/
    d  = 1.e-03;
  }
/*..................................................................*/

/*... A = ((mol/cm^3)^(1-(e1+e2))/((K^alpha)(1/s))  */
  else if(iCod == 3)
  {
/*... kmol/m3 - > mol/cm3*/
    c1 = 1.e-03*density*y1/mW1;
    c2 = 1.e-03*density*y2/mW2;
/*... mol/cm^3s -> kmol/m^3s*/
    d  = 1.e+03;
  }
/*..................................................................*/

  c1 = max(c1,0.e0);
  c2 = max(c2,0.e0);

/*... (c1^a1)x(c1^a2)*/
  prodC = pow(c1,e1)*pow(c2,e2);
/*... exp(Ea/RT)*/
  k = coefA*pow(tc,alpha)*exp(-tA/tc);

  omega = d*prodC*k;

  return omega; 
}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 30/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * arrheniusA : Conversao de unidades da lei e arrhenius             *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * e1       -> expoente                                              *
 * mW1      -> massa molar                                           *
 * y2       -> fracao massica                                        *
 * e2       ->  expoente                                             *
 * mW2      -> massa molar                                           *
 * desnity  -> densidade do fluido                                   *
 * coefA    -> coeficiente                                           *
 * fKelvin  ->                                                       *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * coefA : modificado para incorporar a mudacao das unidades         *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE arrheniusA(DOUBLE const e1     ,DOUBLE const e2
                 ,DOUBLE const mW1    ,DOUBLE const mW2
                 ,DOUBLE const density,DOUBLE const tA    
                 ,DOUBLE const coefA  ,bool const fKelvin)
{

  short iCod=3;
  DOUBLE d;

/*... A = ((kmol/m^3)^(1-(e1+e2))/((K^alpha)(1/s))  */
  if(iCod == 1)
    d = pow(density,e1+e2)*pow(mW1,-e1)*pow(mW2,-e2);
/*..................................................................*/

/*... A = ((mol/m^3)^(1-(e1+e2))/((K^alpha)(1/s))  */
  else if(iCod == 2)
    d = 1.e-03*pow(1.e+03*density,e1+e2)*pow(mW1,-e1)*pow(mW2,-e2);
/*..................................................................*/

/*... A = ((mol/cm^3)^(1-(e1+e2))/((K^alpha)(1/s))  */
  else if(iCod == 3)
    d = 1.e+03*pow(1.e-03*density,e1+e2)*pow(mW1,-e1)*pow(mW2,-e2);
/*..................................................................*/

  return d*coefA;

}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 30/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * arrheniusO : Lei de arrhenius com otimizacao do numero de calculos*
 * para a conversao de unidades                                      *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * y1       -> fracao massica                                        *
 * e1       -> expoente                                              *
 * y2       -> fracao massica                                        *
 * e2       ->  expoente                                             *
 * t        -> temperatura                                           *
 * alpha    -> coeficiente da temperatura                            *
 * tA       -> temperatura de ativacao                               *
 * coefA    -> coeficiente com conversao de unidades imbutidades     *
 * fKelvin  ->                                                       *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * omega :kmol/m^3s                                                  *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * Uso em conjunto com a funcao arrheniusA                           *
 *********************************************************************/
DOUBLE arrheniusO(DOUBLE const y1     ,DOUBLE const y2
                 ,DOUBLE const e1     ,DOUBLE const e2
                 ,DOUBLE const t      ,DOUBLE const alpha
                 ,DOUBLE const tA     ,DOUBLE const coefAl 
                 ,bool const fKelvin)
{
  DOUBLE tc;
  DOUBLE omega,k,prodC;

  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);

/*... (c1^a1)x(c1^a2)*/
  prodC = pow(y1,e1)*pow(y1,e2);
/*... exp(Ea/RT)*/
  k = coefAl*pow(tc,alpha)*exp(-tA/tc);

  omega = prodC*k;

  return omega; 
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

  short m,n,l,posN2=c->sp_N2;
  DOUBLE nAir,nProd,pO2,pN2,pN2r,nN2p,nCO2p,nH2Op,nCOp,nTotal;
  DOUBLE pCO2r,pCOr,pH2Or,nO2,nN2;

  m = c->fuel[iReac].c;
  n = c->fuel[iReac].h;
  l = c->fuel[iReac].o;

  pO2 = c->O2InAir;
  pN2 = c->N2InAir;


/*... Ar ( pO2*O2 + pN2*N2)*/
  nO2  = m + n/4.0e0 - l/2.0;
  nN2  = nO2*pN2/pO2;
  nAir = nO2*(1.0e0/pO2);

/*... produto  */
  nN2p  = (m + n/4.e0)*pN2/pO2;
  nCO2p = m;
  nH2Op = n/2.e0;  

  nTotal = nCO2p + nH2Op + nN2p;
  pCO2r = nCO2p/nTotal;
  pH2Or = nH2Op/nTotal;
  pN2r = nN2p/nTotal;
  nProd = (nN2p/pN2r+nCO2p/pCO2r+nH2Op/pH2Or)/3.0;

  c->stoich[iReac][0][c->sp_fuel[iReac]]= 1.e0;
  c->stoich[iReac][0][c->sp_O2]  = nO2;
  c->stoich[iReac][0][c->sp_N2]  = nN2;
  c->stoich[iReac][1][c->sp_H2O] = nH2Op;
  c->stoich[iReac][1][c->sp_CO2] = nCO2p;
  c->stoich[iReac][1][c->sp_N2]  = nN2p;
  c->stoichAir                   = nAir;

  c->CO2InProd = pCO2r;
  c->H2OInProd = pH2Or;
  c->N2InProd  = pN2r;

/*.*/
  c->nSpeciesPart[iReac][0] = 3;
  c->nSpeciesPart[iReac][1] = 3;
/*. reagente*/
  c->speciesPart[iReac][0][0] = c->sp_fuel[iReac];
  c->speciesPart[iReac][0][1] = c->sp_O2;
  c->speciesPart[iReac][0][2] = c->sp_N2;
/*. produto*/
  c->speciesPart[iReac][1][0] = c->sp_CO2;
  c->speciesPart[iReac][1][1] = c->sp_H2O;
  c->speciesPart[iReac][1][2] = c->sp_N2;

/*.*/
  c->sMassAir             = nAir*c->mW_Air/c->mW[c->sp_fuel[iReac]];
  c->sMass[iReac][0][c->sp_fuel[iReac]] = 1.0;
  c->sMass[iReac][0][c->sp_O2]   = nO2*c->mW[c->sp_O2]/c->mW[c->sp_fuel[iReac]];
  c->sMass[iReac][0][c->sp_N2]   = nN2*c->mW[c->sp_N2]/c->mW[c->sp_fuel[iReac]];
  c->sMass[iReac][1][c->sp_CO2]  = nCO2p*c->mW[c->sp_CO2]/c->mW[c->sp_fuel[iReac]];
  c->sMass[iReac][1][c->sp_H2O]  = nH2Op*c->mW[c->sp_H2O]/c->mW[c->sp_fuel[iReac]];
  c->sMass[iReac][1][c->sp_N2]   = nN2p*c->mW[c->sp_N2]/c->mW[c->sp_fuel[iReac]];

/* fprintf(fileLogExc,"\nReaction:\n\n"); 
/*...*/
  if(c->fLump)
  {
    fprintf(fileLogExc,"C%dH%dO%d" 
           " + %lf (%lf O2 + %lf N2)\n"
           " -> %lf (%lf CO2 + %lf H2O + %lf N2)\n\n"
          ,m    ,n     ,l
          ,nAir ,pO2   ,pN2
          ,nProd,pCO2r, pH2Or,pN2r);
    
    fprintf(fileLogExc,"1 kg C%dH%dO%d" 
           " + %lf (%lf O2 + %lf N2)\n"
           " -> %lf (%lf CO2 + %lf H2O + %lf N2)\n\n"
          ,m                ,n     ,l
          ,c->sMassAir    ,pO2   ,pN2
          ,c->sMassAir + 1,pCO2r,pH2Or,pN2r);

  } 
/*..................................................................*/

/*...*/
  else
  {
    fprintf(fileLogExc,"C%dH%dO%d" 
           " + %lf O2 + %lf N2\n"
           " -> %lf CO2 + %lf H2O + %lf N2\n\n"
          ,m    ,n    ,l
          ,c->stoich[iReac][0][c->sp_O2]
          ,c->stoich[iReac][0][c->sp_N2]
          ,c->stoich[iReac][1][c->sp_CO2]
          ,c->stoich[iReac][1][c->sp_H2O]
          ,c->stoich[iReac][1][c->sp_N2]);

    fprintf(fileLogExc,"1 kg C%dH%dO%d" 
           " + %lf  kg O2 + %lf  kg N2\n"
           " -> %lf kg CO2 + %lf kg H2O + %lf kg N2\n\n"
          ,m    ,n    ,l
          ,c->sMass[iReac][0][c->sp_O2]  
          ,c->sMass[iReac][0][c->sp_N2]
          ,c->sMass[iReac][1][c->sp_CO2]
          ,c->sMass[iReac][1][c->sp_H2O]
          ,c->sMass[iReac][1][c->sp_N2]);  
  }
/*..................................................................*/

}
/*******************************************************************/
void printt(double *x, int n)
{
  int i;
  for(i=0;i<n;i++)
    if( x[i] != 0.e0) printf("%d %e\n",i,x[i]);
}