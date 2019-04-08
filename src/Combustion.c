#include<Combustion.h>

static void vectorMolarMass(DOUBLE *mW, Combustion *cModel);

/*********************************************************************
 * Data de criacao    : 30/07/2018                                   *
 * Data de modificaco : 26/08/2018                                   *
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
                   , short itSimple    ) 
{
  short itComb, conv, jj = 1;
  short nComp = cModel->nComb;
  INT desloc;
  bool fComb[3];
  DOUBLE tmp,tb[3],rCell[3],rCell0[3];
  DOUBLE *b1, *b2, *b3, *bPc, *xu1, *xu2, *xu3, *rCell1, *rCell2, *rCell3;
  DOUBLE *ad1,*ad2,*ad3,*al1,*al2,*al3,*au1,*au2,*au3;

  fComb[0] = fComb[1] = fComb[2] = false; 

  desloc = sistEqComb->neq;
  b1 = sistEqComb->b;
  b2 = &sistEqComb->b[desloc];
  b3 = &sistEqComb->b[2 * desloc];
  
  xu1 = sistEqComb->x;
  xu2 = &sistEqComb->x[desloc];
  xu3 = &sistEqComb->x[2 * desloc];

  ad1 = sistEqComb->ad;
  ad2 = &sistEqComb->ad[desloc];
  ad3 = &sistEqComb->ad[2 * desloc];

  desloc = 2*sistEqComb->nad + sistEqComb->nadr;
  al1 = sistEqComb->al;
  al2 = &sistEqComb->al[desloc];
  al3 = &sistEqComb->al[2 * desloc];

  au1 = sistEqComb->au;
  au2 = &sistEqComb->au[desloc];
  au3 = &sistEqComb->au[2 * desloc];

  desloc = mesh->numel;
  rCell1 = mesh->elm.rCellComb;
  rCell2 = &mesh->elm.rCellComb[desloc];
  rCell3 = &mesh->elm.rCellComb[ 2 * desloc];

  rCell0[0] = rCell0[1] = rCell0[2] = 1.e0;
/*...*/
  if(opt->fItPlot)  
    fprintf(opt->fileItPlot[FITPLOTCOMB]
           ,"#itSimple = %d t = %lf\n",itSimple+1,sc->ddt.t);
/*...................................................................*/

  printf("1 %lf\n",mixtureMolarMass(cModel,mesh->elm.yFrac));

/*...*/
  for(itComb = 0; itComb < 15;itComb++)
  {
/*... taxa de comsumo do combustivel*/
    rateFuelConsume(cModel                , mesh->elm.zComb
                  , mesh->elm.temp        , mesh->elm.rateFuel 
                  , mesh->elm.densityFluid,eModel->fKelvin
                  , mesh->numel);    
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
        , nComp                  , mesh->ndm              
        , &pMesh->iNo            , &pMesh->iEl
        , mesh->numelNov         , mesh->numel
        , mesh->nnodeNov         , mesh->nnode);    
    tm.rcGradComb = getTimeC() - tm.rcGradComb;
/*.................................................................. */

  /*... calculo de: A(i),bE(i)*/
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

/*... soma o vetor bE(i) = bE(i) + bE0(i)*/
    addVector(1.0e0                   , sistEqComb->b
            , 1.0e0                   , sistEqComb->b0
            , sistEqComb->neqNov*nComp, sistEqComb->b);
/*...................................................................*/

/*... soma o vetor RE(i) = RE(i) + bE0(i)*/
    updateCellValueSimple(mesh->elm.rCellComb, sistEqComb->b0
                      , sistEqComb->id       , &sistEqComb->iNeq
                      , mesh->numelNov       , sistEqComb->neqNov
                      , nComp
                      , true                 , false);  
/*...................................................................*/

/*...*/
    tb[0] = sqrt(dot(b1, b1, sistEqComb->neqNov));
    tb[1] = sqrt(dot(b2, b2, sistEqComb->neqNov));
    tb[2] = sqrt(dot(b3, b3, sistEqComb->neqNov));
    if (itComb == 0)
    {
      tmp = max(tb[0], tb[1]);
      tmp = max(tmp, tb[2]);
    }
/*...................................................................*/

/*...*/
    fComb[0] = fComb[1] = fComb[2] = true; 
    if ( tb[0] < SZERO || tb[0] == 0.e0 ) fComb[0] = false;
    if ( tb[1] < SZERO || tb[1] == 0.e0 ) fComb[1] = false;
    if ( tb[2] < SZERO || tb[2] == 0.e0 ) fComb[2] = false;
/*...................................................................*/

/*...*/ 
    if( itComb == 0 )
    {
      if(fComb[0]) 
        rCell[0]  = rCell0[0] = sqrt(dot(rCell1,rCell1,mesh->numelNov));
      if(fComb[1]) 
        rCell[1]  = rCell0[1] = sqrt(dot(rCell2,rCell2,mesh->numelNov));
      if(fComb[2])
        rCell[2]  = rCell0[2] = sqrt(dot(rCell3,rCell3,mesh->numelNov));
    }
    else
    {
      rCell[0] = sqrt(dot(rCell1,rCell1,mesh->numelNov));
      rCell[1] = sqrt(dot(rCell2,rCell2,mesh->numelNov));
      rCell[2] = sqrt(dot(rCell3,rCell3,mesh->numelNov));
    }   
    if( jj == 20 && !mpiVar.myId )
    {
      jj = 0;
      printf("It Combustion: %d \n", itComb + 1);
      printf("Residuo:\n");
      printf("zFlue : %20.8e\n", rCell[0] / rCell0[0]);
      printf("zAir  : %20.8e\n", rCell[1] / rCell0[1]);
      printf("zPord : %20.8e\n", rCell[2] / rCell0[2]);

    }
    jj++;
    if(opt->fItPlot && !mpiVar.myId)  
        fprintf(opt->fileItPlot[FITPLOTCOMB]
               ,"%9d  %20.8e  %20.8e  %20.8e\n",itComb+1
                                               ,rCell[0]
                                               ,rCell[1]
                                               ,rCell[2]);
/*...................................................................*/

/*...*/
    conv = 0;
    if(rCell[0]/ rCell0[0] < 1.e-6 || rCell[0] < 1.e-16) conv++;
    if(rCell[1]/ rCell0[1] < 1.e-6 || rCell[1] < 1.e-16) conv++;
    if(rCell[2]/ rCell0[2] < 1.e-6 || rCell[2] < 1.e-16) conv++;
    if(conv == 3) break;
/*...................................................................*/

/*...AZ1=bZ1*/
    tm.solvComb = getTimeC() - tm.solvComb;
    if(fComb[0])
      solverC(m
            , sistEqComb->neq    , sistEqComb->neqNov
            , sistEqComb->nad    , sistEqComb->nadr
            , sistEqComb->ia     , sistEqComb->ja
            , al1                , ad1             , au1
            , b1                 , xu1                
            , &sistEqComb->iNeq  , &sistEqComb->omp
            , solvComb->tol      , solvComb->maxIt
            , sistEqComb->storage, solvComb->solver
            , solvComb->fileSolv , solvComb->log
            , true               , sistEqComb->unsym);
/*...AZ2=bZ2*/
    if(fComb[1]) 
      solverC(m
            , sistEqComb->neq, sistEqComb->neqNov
            , sistEqComb->nad, sistEqComb->nadr
            , sistEqComb->ia , sistEqComb->ja
            , al2            , ad2           , au2
            , b2             , xu2
            , &sistEqComb->iNeq, &sistEqComb->omp
            , solvComb->tol, solvComb->maxIt
            , sistEqComb->storage, solvComb->solver
            , solvComb->fileSolv, solvComb->log
            , true, sistEqComb->unsym);        
  
/*...AZ3=bZ3*/
    if(fComb[2]) 
      solverC(m
         , sistEqComb->neq    , sistEqComb->neqNov
         , sistEqComb->nad    , sistEqComb->nadr
         , sistEqComb->ia     , sistEqComb->ja
         , al3                , ad3         , au3
         , b3                 , xu3
         , &sistEqComb->iNeq  , &sistEqComb->omp
         , solvComb->tol      , solvComb->maxIt
         , sistEqComb->storage, solvComb->solver
         , solvComb->fileSolv , solvComb->log
         , true, sistEqComb->unsym);        
    tm.solvComb = getTimeC() - tm.solvComb;    
/*...................................................................*/

/*... x -> zComb*/
    updateCellValueBlock(mesh->elm.zComb , sistEqComb->x
                       , sistEqComb->id  , &sistEqComb->iNeq
                       , mesh->numel     , sistEqComb->neq 
                       , cModel->nComb   
                       , cModel->fRes    , true);  
/*...................................................................*/

/*...*/
//  regularZ(mesh->elm.zComb,mesh->numelNov,cModel->nOfSpeciesLump);
/*...................................................................*/

  }
/*...................................................................*/

/*...*/
  getSpeciesPrimitives(cModel
                      ,mesh->elm.yFrac,mesh->elm.zComb
                      ,mesh->numelNov);
/*...................................................................*/

/*...*/
  rateHeatRealeseCombustion(cModel                
                    , mesh->elm.rateHeatReComb     
                    , mesh->elm.zComb0        , mesh->elm.zComb
                    , mesh->elm.densityFluid  , mesh->elm.rateFuel    
                    , sc->ddt.dt[TIME_N]      , mesh->numelNov);
/*...................................................................*/

  prop->molarMass = mixtureMolarMass(cModel,mesh->elm.yFrac);
//printf("2 %lf\n", mixtureMolarMass(cModel,mesh->elm.yFrac));

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 15/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * regularZ : mantem o z entre 0 e 1                                 *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * zComb   -> fracao massica                                         *
 * numel   -> numero de elementos                                    *
 * nLump   -> numero de especie agrupadas                            *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * y -> fracao de massa de especies primitivas                       * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * z < 0.0 -> z = 0.e0                                               *
 * z > 1.0 -> z = 1.e0 - (zFuel + zAir)                              *
 *********************************************************************/
void regularZ(DOUBLE *RESTRICT z, INT const numel, short const nLump)
{
  short j;
  INT nel;
  DOUBLE zSum;

/*... z < 0.0*/
  for(nel = 0 ; nel < numel; nel++)
  {
    for(j=0;j<nLump;j++)
      if(MAT2D(nel,j,z,nLump) < 0.e0)
        MAT2D(nel,j,z,nLump) = 0.e0;
  }

/*... z <= 1.0*/
  for(nel = 0 ; nel < numel; nel++)
  {
    for(j=0,zSum = 0.e0;j<nLump;j++)
      zSum += MAT2D(nel,j,z,nLump);
    if(zSum > 1.e0)
      MAT2D(nel,2,z,nLump) = 1.e0 
                       - (MAT2D(nel,0,z,nLump) + MAT2D(nel,1,z,nLump));
  }

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 15/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
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
 * y[1] - N2                                                         *
 * y[2] - O2                                                         *
 * y[3] - CO2                                                        *
 * y[4] - H2O                                                        *
 * y[5] - CO                                                         *
 * y[6] - C                                                          *
 *********************************************************************/
void getSpeciesPrimitives(Combustion *cModel 
                        , DOUBLE *RESTRICT y,DOUBLE *RESTRICT z
                        , INT const numel)
{
  short ns,nl;
  INT nel;
  DOUBLE *py,*pz,*a;

  ns = cModel->nOfSpecies; 
  nl = cModel->nOfSpeciesLump;
  a  = cModel->lumpedMatrix;

  for(nel = 0 ; nel < numel; nel++)
  {
    py = &MAT2D(nel,0,y,ns);
    pz = &MAT2D(nel,0,z,nl);
    yLumpedMatrixZ(py,a,pz,ns,nl);
  }

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 26/08/2018                                   *
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
 * z(nel,1) -> air (oxidante)                                        *
 * z(nel,2) -> fuel(comburente)                                      *
 * z(nel,3) -> produto                                               *
 *********************************************************************/
void rateFuelConsume(Combustion *cModel       , DOUBLE *RESTRICT zComb
                   , DOUBLE *RESTRICT temp    , DOUBLE *RESTRICT rate
                   , DOUBLE *RESTRICT density , bool const fKelvin
                   , INT const numel)
{
  short nComb = cModel->nComb,iCod = cModel->reactionKinetic;
  INT nel;
  DOUBLE s,tMix;
  DOUBLE zFuel, zAir, omega, densityC, alpha, coefA;
  DOUBLE tmp1, tmp2, tmp3, tmp4;
  DOUBLE mWfuel,mWox, tempA, eOx, eFuel, ru, tc; 

/*...*/
  switch(iCod)
  {
/*...*/
    case ARRHENIUS:
      s      = cModel->sMolar;
      alpha  = cModel->arrhenius.alpha;
      mWfuel = cModel->mW_Fuel;
      mWox   = cModel->mW_Air;
      coefA  = cModel->arrhenius.a;
/*... KJ/(kmol*kelvin)*/
      ru     = IDEALGASR*1.e-03;
      tempA  = cModel->arrhenius.energyAtivation/ru; 
      coefA *= mWfuel/(mWfuel+pow(mWox,s));
      eOx     = s;
      eFuel   = 1.0;

      for(nel = 0; nel < numel; nel++)
      {

        if(fKelvin)
          tc = temp[nel];  
        else
          tc = CELSIUS_FOR_KELVIN(temp[nel]);

        densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);

        zAir  = MAT2D(nel,0,zComb,nComb);
        zFuel = MAT2D(nel,1,zComb,nComb);
   
        tmp1 = pow(temp[nel],alpha); 
        tmp2 = pow(densityC,s+1.0);
        tmp3 = pow(zFuel,eFuel)*pow(zAir,eOx);
        tmp4 = exp(-tempA/tc);
//      printf("%lf %lf\n",tmp3,tmp4);

        rate[nel] = coefA*tmp1*tmp2*tmp3*tmp4;
      }
    break;
/*...................................................................*/

/*...*/
    case EBU:
      s    = cModel->sMass;
      tMix = cModel->tMix;
      for(nel = 0; nel < numel; nel++)
      {

        densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);

        zAir  = MAT2D(nel,0,zComb,nComb);
        zFuel = MAT2D(nel,1,zComb,nComb);
    
        omega = max(densityC*min(zFuel,zAir/s),0.e0);
        rate[nel] = omega/tMix;
      }
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
 * rateHeatRealeseCombustion: calculo da taxa de liberacao de calor  *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> modelo de combustao                                    *
 * q       -> nao definido                                           *
 * zComb0  -> fracao de massa agrupada do passo de tempo anterior    *
 * zComb   ->fracao de massa agrupada do passo de tempo atural       *
 * rateFuel-> taxa de consumo do combustivel                         * 
 * dt      -> delta dessa passo de tempo                             * 
 * numel   -> numero de elementos                                    *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * q -> taxa de energia liberada KJ/s                                * 
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * entalphyOfFormLumped[0] - Air KJ/Kg                               *
 * entalphyOfFormLumped[1] - Fuel KJ/Kg                              *
 * entalphyOfFormLumped[2] - Froduto KJ/Kg                           *
 *********************************************************************/
void rateHeatRealeseCombustion(Combustion *cModel   
                   , DOUBLE *RESTRICT q
                   , DOUBLE *RESTRICT zComb0 , DOUBLE *RESTRICT zComb
                   , DOUBLE *RESTRICT density, DOUBLE *RESTRICT rateFuel 
                   , DOUBLE const dt         , INT const numel)
{

  short i, nOfSpL = cModel->nOfSpeciesLump;
  short iCod = cModel->typeHeatRealese;
  INT nel;
  DOUBLE s = cModel->sMass;
  DOUBLE *h,hc;
  DOUBLE z0,z1;
  DOUBLE densityC,densityC0,tmp,dRoZ;

  if(iCod == HFORMATION)
  {
    h  = cModel->entalphyOfFormLumped;
    for(nel = 0; nel < numel; nel++)
    {
      tmp = ((1.e0+s)*h[2] - s*h[0] - h[1])*rateFuel[nel];
      q[nel] = -tmp;
    }
  }
  else if (iCod == HCOMBUSTION)
  {
    hc = cModel->entalphyOfCombustion;
    for(nel = 0; nel < numel; nel++)
    {
      tmp = rateFuel[nel]*hc;
      q[nel] = -tmp;
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
 * entalphyOfFormLumped[0] - Air KJ/Kg                               *
 * entalphyOfFormLumped[1] - Fuel KJ/Kg                              *
 * entalphyOfFormLumped[2] - Froduto KJ/Kg                           *
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
 * Data de modificaco : 00/00/0000                                   *
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

  short nComb = cModel->nComb;

/*... Fuel*/
  MAT2D(0,0,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(0,1,cModel->lumpedMatrix,nComb) = 1.e0; 
  MAT2D(0,2,cModel->lumpedMatrix,nComb) = 0.e0; 
/*...................................................................*/

/*... N2*/
  MAT2D(1,0,cModel->lumpedMatrix,nComb) = 0.7670e0; 
  MAT2D(1,1,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(1,2,cModel->lumpedMatrix,nComb) = 0.7248e0; 
/*...................................................................*/

/*... O2*/
  MAT2D(2,0,cModel->lumpedMatrix,nComb) = 0.2330e0; 
  MAT2D(2,1,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(2,2,cModel->lumpedMatrix,nComb) = 0.000000e0; 
/*...................................................................*/

/*... CO2*/
  MAT2D(3,0,cModel->lumpedMatrix,nComb) = 0.000e0; 
  MAT2D(3,1,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(3,2,cModel->lumpedMatrix,nComb) = 0.1514e0; 
/*...................................................................*/

/*... H2O*/
  MAT2D(4,0,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(4,1,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(4,2,cModel->lumpedMatrix,nComb) = 0.1238e0; 
/*...................................................................*/

/*... CO*/
  MAT2D(5,0,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(5,1,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(5,2,cModel->lumpedMatrix,nComb) = 0.e0; 
/*...................................................................*/

/*... C*/
  MAT2D(6,0,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(6,1,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(6,2,cModel->lumpedMatrix,nComb) = 0.e0; 
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
 * y[0] - CH4                                                        *
 * y[1] - N2                                                         *
 * y[2] - O2                                                         *
 * y[3] - CO2                                                        *
 * y[4] - CO                                                         *
 * y[5] - H2O                                                        *
 * y[6] - C                                                          *
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
  
  short nO,nC,nH;
  DOUBLE pO2,pN2;

  nC = cModel->fuel.c;
  nO = cModel->fuel.o;
  nH = cModel->fuel.h;

  pO2 = cModel->O2InAir;
  pN2 = cModel->N2InAir;

/*... Fuel */
  cModel->mW_Fuel = nC*MW_C + nH*MW_H + nO*MW_O;

/*...O2*/
  cModel->mW_O2 = 2.0e0*MW_O;

/*...N2*/
  cModel->mW_N2 = 2.0e0*MW_N;

/*...H2O*/
  cModel->mW_H2O = 2.0e0*MW_H + MW_O;

/*...CO2*/
  cModel->mW_CO2 =MW_C +  2.0e0*MW_O;

/*...CO*/
  cModel->mW_CO = MW_C + MW_O;

/*...C*/
  cModel->mW_C = MW_C;

/*... Air*/
  cModel->mW_Air = pO2*cModel->mW_O2 + pN2*cModel->mW_N2;

}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 23/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
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
 *                                                                   *
 * 1*(CmHnOl) + nO2*(O2) + (nO2*pN2/pO2)*N2 - nCO2*(CO2) + nH2O*(H2O)*
 * nN2*(N2)                                                          *
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
void stoichiometricCoeff(Combustion *cModel)
{
  short m,n,l;
  DOUBLE nAir,nProd,pO2,pN2,pN2r,nN2r,nCO2r,nH2Or,nTotal;
  DOUBLE pCO2r,pH2Or,nO2;

  m = cModel->fuel.c;
  n = cModel->fuel.h;
  l = cModel->fuel.o;

  pO2 = cModel->O2InAir;
  pN2 = cModel->N2InAir;

/* FUEL + nAir*AIR -> nPrpd*PROD*/

/*... Ar ( pO2*O2 + pN2*N2)*/
//nAir = 2.e0*(1.0e0/(pO2));

/*... porduto  */
//nCO2r  = 1.0;
//nH2Or  = 2.0;
//nN2r   = nAir*pN2;
//nTotal = nCO2r + nH2Or + nN2r;
//pCO2 = nCO2r/nTotal;
//pH2O = nH2Or/nTotal;
//pN2  = nN2r/nTotal;
//nProd = nN2r/pN2;

//cModel->CO2InPord = pCO2;
//cModel->H2OInPord = pH2O;

/*... Ar ( pO2*O2 + pN2*N2)*/
  nO2  = m + n/4.0e0 - l/2.0;
  nAir = nO2*(1.0e0/pO2);

/*... produto  */
  nN2r  = (m + n/4.e0)*pN2/pO2;
  nCO2r = m;
  nH2Or = n/2.e0;  

  nTotal = nCO2r + nH2Or + nN2r;
  pCO2r = nCO2r/nTotal;
  pH2Or = nH2Or/nTotal;
  pN2r = nN2r/nTotal;
  nProd = (nN2r/pN2r+nCO2r/pCO2r+nH2Or/pH2Or)/3.0;

  cModel->CO2InPord = pCO2r;
  cModel->H2OInPord = pH2Or;

/*...*/
  cModel->sMolar = nAir;
  cModel->sMass  = nAir*cModel->mW_Air/cModel->mW_Fuel;

}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/08/2018                                   *
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
void initEntalpyOfFormation(Combustion *cModel)
{
  
  DOUBLE hFuelFormation,hCO2Formation;
  DOUBLE pCO2,pH2O;


  pCO2 = cModel->CO2InPord;
  pH2O = cModel->H2OInPord;

/*...  Fuel - KJ/kMol*/
  cModel->H_Fuel = -74870.e0; 
/*... CO2 - KJ/kMol*/
  cModel->H_CO2 = -393510.e0;
/*... H2O - KJ/kMol*/
  cModel->H_H2O = -241826.e0;
/*... CO - KJ/kMol*/
  cModel->H_CO = -110500.e0;

/*... Air - KJ/Kg*/
  cModel->entalphyOfFormLumped[0] = 0.e0;

/*... Fuel */
  cModel->entalphyOfFormLumped[1] = cModel->H_Fuel / cModel->mW_Fuel;
 
/*... Prod */
  cModel->entalphyOfFormLumped[2] =  pCO2*cModel->H_CO2 / cModel->mW_CO2 
                                   + pH2O*cModel->H_H2O / cModel->mW_H2O;
 
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/08/2018                                   *
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
void initEntalpyOfCombustion(Combustion *cModel)
{
  
  DOUBLE hFuelFormation,hCO2Formation;
  DOUBLE pCO2,pH2O;


  pCO2 = cModel->CO2InPord;
  pH2O = cModel->H2OInPord;

/*...  Fuel - KJ/kMol*/
  cModel->H_Fuel = -74870.e0; 
/*... CO2 - KJ/kMol*/
  cModel->H_CO2 = -393520.e0;
/*... H2O - KJ/kMol*/
  cModel->H_H2O = -241830.e0;

/*... KJ/Kmol de fuel*/
  cModel->entalphyOfCombustion = (cModel->H_CO2 + 2.0*cModel->H_H2O) 
                               - (cModel->H_Fuel);
/*... KJ/KG   de fuel */
  cModel->entalphyOfCombustion  /= cModel->mW_Fuel;

}
/********************************************************************/


/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * concetracionOfSpecies:                                            *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * cModel  -> nao definido                                           *
 * z       -> especies agrupadas                                     *
 * c       -> concetracao das especies primitivas                    *
 * density -> densidade da mistura                                   *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * c       -> concentracao das especies                              *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * c[0] - CH4                                                        *
 * c[1] - N2                                                         *
 * c[2] - O2                                                         *
 * c[3] - CO2                                                        *
*  c[4] - H2O                                                        *
 * c[5] - CO                                                         * 
 * c[6] - C                                                          *
 *********************************************************************/
void concetracionOfSpecies(Combustion *cModel,DOUBLE *RESTRICT z
                          ,DOUBLE *RESTRICT c,DOUBLE const density)
{
  
  short i, ns = cModel->nOfSpecies, nl = cModel->nOfSpeciesLump;
  DOUBLE y[7],mW[7];

  vectorMolarMass(mW,cModel);

  yLumpedMatrixZ(y,cModel->lumpedMatrix, z, ns, nl);
  
  for(i=0;i<ns;i++)
    c[i] = y[i]*density/mW[i];

}
/*********************************************************************/

/**********************************************************************
 * Data de criacao    : 15/08/2018                                    *
 * Data de modificaco : 00/00/0000                                    *
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
  short i,j;

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
  DOUBLE mW[7],tmp;

  vectorMolarMass(mW,cModel);

  for(i=0,tmp=0.e0;i<ns;i++)
    tmp += y[i]/mW[i];

  return 1.e0/tmp;

}
/*********************************************************************/


/********************************************************************/
static void vectorMolarMass(DOUBLE *mW, Combustion *cModel)
{

  mW[0] = cModel->mW_Fuel;
  mW[1] = cModel->mW_N2;
  mW[2] = cModel->mW_O2;
  mW[3] = cModel->mW_CO2;
  mW[4] = cModel->mW_H2O;
  mW[5] = cModel->mW_CO;
  mW[6] = cModel->mW_C;

}
/*********************************************************************/