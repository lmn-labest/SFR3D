#include<Combustion.h>

static void vectorMolarMass(DOUBLE *mW, Combustion *cModel);

/*********************************************************************
 * Data de criacao    : 30/07/2018                                   *
 * Data de modificaco : 02/05/2019                                   *
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
  bool fComb[MAX_COMB],rel;
  short itComb, conv, i, j, kZero, nComb, jj = 1, jPrint;
  INT desloc;
  DOUBLE tmp,tb[MAX_COMB],rCell[MAX_COMB],rCell0[MAX_COMB];
  DOUBLE *b[MAX_COMB], *bPc, *xu[MAX_COMB];
  DOUBLE *ad[MAX_COMB],*al[MAX_COMB],*au[MAX_COMB],*rCellC[MAX_COMB];
  char slName[][10] = {"zFuel","zAir","zProd"},
       spName[][10] = {"zFuel","zO2" ,"zN2","zCO2","zH2O"};
/*... rel   - criterio de parada com vlores relativos
      kZero - iteracao na qual o resido e normalizadp  
      nComb - numero de especies transportadas
      jPrint - numero de iteracao que ser impresso o residuo na tela
*/
  rel    = true;          
  kZero  = 0;
  nComb  = cModel->nComb;
  jPrint = 20;
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
  
    rCell0[i] = 1.e0;
  }
/*...................................................................*/

/*...*/
  if(opt->fItPlot && !mpiVar.myId)  
    fprintf(opt->fileItPlot[FITPLOTCOMB]
           ,"#itSimple = %d t = %lf\n",itSimple+1,sc->ddt.t);
/*...................................................................*/

  printf("1 %lf\n",mixtureMolarMass(cModel,mesh->elm.yFrac));

/*...*/
  for(itComb = 0; itComb < 10;itComb++)
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
    updateCellValueSimple(mesh->elm.rCellComb, sistEqComb->b0
                      , sistEqComb->id       , &sistEqComb->iNeq
                      , mesh->numelNov       , sistEqComb->neqNov
                      , nComb
                      , true                 , false);  
/*...................................................................*/

/*...*/
    if(jj == jPrint && !mpiVar.myId)
    {
      printf("\nIt Combustion: %d\n", itComb + 1);
      printf("Residuo:\n");
    }
    if(opt->fItPlot && !mpiVar.myId)
      fprintf(opt->fileItPlot[FITPLOTCOMB],"%9d",itComb+1);
/*...................................................................*/
    
/*...*/
    for(i=0;i<nComb;i++)
      tb[i] = sqrt(dot(b[i], b[i], sistEqComb->neqNov));
    if (itComb == 0) tmp = maxArray(tb,nComb);
/*...................................................................*/

/*...*/
    for(i=0;i<nComb;i++)
    {
 /*...*/
      fComb[i] = true; 
      if ( tb[i] < SZERO || tb[i] == 0.e0 ) fComb[i] = false;
/*...................................................................*/

/*...*/ 
      if( itComb == kZero && rel )
        rCell[i]  = rCell0[i] = sqrt(dot(rCellC[i],rCellC[i]
                                    ,mesh->numelNov));  
      rCell[i] = sqrt(dot(rCellC[i],rCellC[i],mesh->numelNov));
                
      if( jj == jPrint && !mpiVar.myId )
      {
        if(cModel->fLump)
        {
          printf("%5s : %20.8e\n", slName[i],rCell[i] / rCell0[i]);
        }
        else
        {
          printf("%5s : %20.8e\n", spName[i],rCell[i] / rCell0[i]);
        }  
      }
      if(opt->fItPlot && !mpiVar.myId)
        fprintf(opt->fileItPlot[FITPLOTCOMB]," %20.8e ",rCell[i]);
/*...................................................................*/

/*...*/
      conv = 0;
      if(rCell[i]/ rCell0[i] < 1.e-6 || rCell[i] < 1.e-16) conv++;
      if(conv == nComb) break;
/*...................................................................*/

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
      else
        zero(xu[i],sistEqComb->neq,DOUBLEC);
//    if(i==0)printf("%.16e %.16e %.16e %.16e\n",*ad[i],*al[i],*b[i],*xu[i]);
      tm.solvComb = getTimeC() - tm.solvComb;    
    }
/*...................................................................*/

/*...*/
    if( jj == jPrint) jj = 0;
    jj++;
    if(opt->fItPlot && !mpiVar.myId) 
      fprintf(opt->fileItPlot[FITPLOTCOMB],"\n");
/*...................................................................*/

/*... x -> zComb*/
    updateCellValueBlock(mesh->elm.zComb , sistEqComb->x
                       , sistEqComb->id  , &sistEqComb->iNeq
                       , mesh->numel     , sistEqComb->neq 
                       , cModel->nComb   
                       , cModel->fRes    , true);  
/*...................................................................*/

  }
/*...................................................................*/

/*...*/
//  regularZ(mesh->elm.zComb,mesh->numelNov
//          ,cModel->nComb  ,cModel->fLump);
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

/*...*/
  if(opt->fItPlot)  
    fprintf(opt->fileItPlot[FITPLOTCOMB],"\n");
/*...................................................................*/

/*...*/
  prop->molarMass = mixtureMolarMass(cModel,mesh->elm.yFrac);
  printf("2 %lf\n", mixtureMolarMass(cModel,mesh->elm.yFrac));
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 15/08/2018                                   *
 * Data de modificaco : 03/05/2019                                   *
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
 * fLump = true                                                      *
 * sum(z) > 1.0 -> z[zProd] = 1.e0 - (zFuel + zAir)                  *
* fLump = false                                                      *
 * sum(z) > 1.0 -> z[zN2] = 1.e0 - (zFuel + zO2 + zCO2 + zH2O)       * 
 *********************************************************************/
void regularZ(DOUBLE *RESTRICT z    , INT const numel
            , short const nComb      , bool fLump)
{
  short j;
  INT nel;
  DOUBLE zSum;

/*... z < 0.0*/
  for(nel = 0 ; nel < numel; nel++)
  {
    for(j=0;j<nComb;j++)
      if(MAT2D(nel,j,z,nComb) < 0.e0)
        MAT2D(nel,j,z,nComb) = 0.e0;
  }

/*... z <= 1.0*/
  for(nel = 0 ; nel < numel; nel++)
  {
    for(j=0,zSum = 0.e0;j<nComb;j++)
      zSum += MAT2D(nel,j,z,nComb);
    if(zSum > 1.e0)
    {
      if(fLump)
        MAT2D(nel,SL_PROD,z,nComb) = 1.e0 
        - (MAT2D(nel,SL_AIR,z,nComb) + MAT2D(nel,SL_PROD,z,nComb));
      else
        MAT2D(nel,SP_N2,z,nComb) = 1.e0 
        - (MAT2D(nel,SP_FUEL,z,nComb) + MAT2D(nel,SP_O2 ,z,nComb)
        +  MAT2D(nel,SP_CO2 ,z,nComb) + MAT2D(nel,SP_H2O,z,nComb));
    } 
  }

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 15/08/2018                                   *
 * Data de modificaco : 03/05/2019                                   *
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
  short i,ns,nl;
  INT nel;
  DOUBLE *py,*pz,*a;

  if(cModel->fLump)
  {
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
  else
  {
    ns = cModel->nOfSpecies; 
    for(nel = 0 ; nel < numel; nel++)
      for(i=0;i<ns;i++)
        MAT2D(nel,i,y,ns) = MAT2D(nel,i,z,ns);
  }
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 12/08/2018                                   *
 * Data de modificaco : 04/05/2019                                   *
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
  DOUBLE zFuel, zAir, zO2, omega, densityC, alpha, coefA;
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
      if(cModel->fLump)
      {
        s    = cModel->sMassAir;
        tMix = cModel->tMix;
        for(nel = 0; nel < numel; nel++)
        {

          densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);

          zAir  = MAT2D(nel,SL_AIR ,zComb,nComb);
          zFuel = MAT2D(nel,SL_FUEL,zComb,nComb);
    
          omega = max(densityC*min(zFuel,zAir/s),0.e0);
          rate[nel] = omega/tMix;
        }
      }
      else
      {
        s    = cModel->sMassO2;
        tMix = cModel->tMix;
        for(nel = 0; nel < numel; nel++)
        {

          densityC = MAT2D(nel, TIME_N, density, DENSITY_LEVEL);

          zO2   = MAT2D(nel,SP_O2,zComb,nComb);
          zFuel = MAT2D(nel,SP_FUEL,zComb,nComb);
    
          omega = max(densityC*min(zFuel,zO2/s),0.e0);
          rate[nel] = omega/tMix;
        }
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
  DOUBLE s = cModel->sMassAir;
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
 * Data de modificaco : 04/05/2019                                   *
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
  MAT2D(SP_FUEL,0,cModel->lumpedMatrix,nComb) = 1.e0; 
  MAT2D(SP_FUEL,1,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(SP_FUEL,2,cModel->lumpedMatrix,nComb) = 0.e0; 
/*...................................................................*/

/*... O2*/
  MAT2D(SP_O2,0,cModel->lumpedMatrix,nComb) = 0.e0;
  MAT2D(SP_O2,1,cModel->lumpedMatrix,nComb) = 0.2330e0; 
  MAT2D(SP_O2,2,cModel->lumpedMatrix,nComb) = 0.0e0; 
/*...................................................................*/

/*... N2*/
  MAT2D(SP_N2,0,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(SP_N2,1,cModel->lumpedMatrix,nComb) = 0.7670e0; 
  MAT2D(SP_N2,2,cModel->lumpedMatrix,nComb) = 0.7248e0; 
/*...................................................................*/

/*... CO2*/
  MAT2D(SP_CO2,0,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(SP_CO2,1,cModel->lumpedMatrix,nComb) = 0.0e0; 
  MAT2D(SP_CO2,2,cModel->lumpedMatrix,nComb) = 0.1514e0; 
/*...................................................................*/

/*... H2O*/
  MAT2D(SP_H2O,0,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(SP_H2O,1,cModel->lumpedMatrix,nComb) = 0.e0; 
  MAT2D(SP_H2O,2,cModel->lumpedMatrix,nComb) = 0.1238e0; 
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
 * y[1] - O2                                                         *
 * y[2] - N2                                                         *
 * y[3] - CO2                                                        *
 * y[4] - H2O                                                        *
 * y[5] - CO                                                         *
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
 * Data de modificaco : 03/05/2019                                   *
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
void stoichiometricCoeff(Combustion *cModel)
{
  short m,n,l;
  DOUBLE nAir,nProd,pO2,pN2,pN2r,nN2p,nCO2p,nH2Op,nTotal;
  DOUBLE pCO2r,pH2Or,nO2,nN2;

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

  cModel->CO2InProd = pCO2r;
  cModel->H2OInProd = pH2Or;
  cModel->N2InProd  = pN2r;

/*...*/
  cModel->sMolar    = nAir;
  cModel->sMassAir  = nAir*cModel->mW_Air/cModel->mW_Fuel;
  cModel->sMassO2   = nO2*cModel->mW_O2/cModel->mW_Fuel;
  cModel->sMassN2   = nN2*cModel->mW_N2/cModel->mW_Fuel;
  cModel->sMassCO2p = nCO2p*cModel->mW_CO2/cModel->mW_Fuel;
  cModel->sMassH2Op = nH2Op*cModel->mW_H2O/cModel->mW_Fuel;
  cModel->sMassN2p  = nN2p*cModel->mW_N2/cModel->mW_Fuel;
/*..................................................................*/

  printf("Reaction:\n"); 
/*...*/
  if(cModel->fLump)
  {
    printf("C%dH%dO%d" 
           " + %lf (%lf O2 + %lf N2)"
           " -> %lf (%lf CO2 + %lf H2O + %lf N2)\n\n"
          ,m    ,n     ,l
          ,nAir ,pO2   ,pN2
          ,nProd,pCO2r, pH2Or,pN2r);
    
    printf("1 kg C%dH%dO%d" 
           " + %lf (%lf O2 + %lf N2)"
           " -> %lf (%lf CO2 + %lf H2O + %lf N2)\n\n"
          ,m                ,n     ,l
          ,cModel->sMassAir    ,pO2   ,pN2
          ,cModel->sMassAir + 1,pCO2r,pH2Or,pN2r);

  } 
/*..................................................................*/

/*...*/
  else
  {
    printf("C%dH%dO%d" 
           " + %lf O2 + %lf N2"
           " -> %lf CO2 + %lf H2O + %lf N2\n\n"
          ,m    ,n    ,l
          ,nO2  ,nO2*pN2/pO2
          ,nCO2p,nH2Op,nN2p);

     printf("1 kg C%dH%dO%d" 
           " + %lf  kg O2 + %lf  kg N2 "
           " -> %lf kg CO2 + %lf kg H2O + %lf kg N2\n\n"
          ,m    ,n    ,l
          ,cModel->sMassO2  ,cModel->sMassN2
          ,cModel->sMassCO2p,cModel->sMassH2Op,cModel->sMassN2p);
  }
/*..................................................................*/
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


  pCO2 = cModel->CO2InProd;
  pH2O = cModel->H2OInProd;

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


  pCO2 = cModel->CO2InProd;
  pH2O = cModel->H2OInProd;

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
 * c[0] - CH4                                                        *
 * c[1] - O2                                                         *
 * c[2] - N2                                                         *
 * c[3] - CO2                                                        *
*  c[4] - H2O                                                        *
 * c[5] - CO                                                         * 
 * c[6] - C                                                          *
 *********************************************************************/
void concetracionOfSpecies(Combustion *cModel,DOUBLE *RESTRICT z
                          ,DOUBLE *RESTRICT c,DOUBLE const density)
{
  bool fLump = cModel->fLump;
  short i, ns = cModel->nOfSpecies, nl = cModel->nOfSpeciesLump;
  DOUBLE y[MAXSPECIES],mW[MAXSPECIES];

  vectorMolarMass(mW,cModel);

  if (fLump) yLumpedMatrixZ(y,cModel->lumpedMatrix, z, ns, nl);
  
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
  DOUBLE mW[MAXSPECIES],tmp;

  vectorMolarMass(mW,cModel);

  for(i=0,tmp=0.e0;i<ns;i++)
    tmp += y[i]/mW[i];

  return 1.e0/tmp;

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 05/05/2019                                   *
 *-------------------------------------------------------------------*
 * vectorMolarMass : massa molar da mistura                          *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * mW -> nao definido                                                *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * mW -> massa molar                                                 *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
static void vectorMolarMass(DOUBLE *mW, Combustion *cModel)
{

  mW[SP_FUEL] = cModel->mW_Fuel;
  mW[SP_O2  ] = cModel->mW_O2;
  mW[SP_N2  ] = cModel->mW_N2;
  mW[SP_CO2 ] = cModel->mW_CO2;
  mW[SP_H2O ] = cModel->mW_H2O;
/*mW[5] = cModel->mW_CO;
  mW[6] = cModel->mW_C;*/

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

void printt(double *x, short n)
{
  int i;
  for(i=0;i<n;i++)
    printf("%d %e\n",i,x[i]);
}