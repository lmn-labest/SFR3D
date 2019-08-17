#include<HccaTime.h>

double getTimeC(void){
/**********************************************************************
 * GETTIMEC : Obter tempo                                             *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * -------------------------------------------------------------------*
*********************************************************************/

/*funcao MPI */
#if _MPI_
  double t=0.0;
  t = MPI_Wtime();
  return t;
/*funcao do openmp*/
#elif _OPENMP 
  #if _MPI_
    double t=0.0;
    t = MPI_Wtime();
    return t;
  #else
    double t=0.0;
    t = omp_get_wtime();
    return t;
  #endif
/*funcao do nativa do C (baixa precisao)*/
#elif _WIN32
   clock_t start;
   double  t;

   start = clock();
   t = (double)(start) / CLOCKS_PER_SEC;
   return t;
#else
  double t =0.0;
  struct timeval start;

  gettimeofday(&start,NULL);
  t = (double) (start.tv_sec + start.tv_usec/1000000.e0);
  
  return t;
#endif

}

/*********************************************************************/

/**********************************************************************
 * Data de criacao    : 03/12/2017                                    *
 * Data de modificaco : 08/06/2019                                    *
 * -------------------------------------------------------------------*
 * initTime : inicializacao da variaveis de tempo                     *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 * ------------------------------------------------------------------ *
*********************************************************************/
void initTime(Time *tm) {

/*... Time*/
  tm->adjcency          = 0.e0;
  tm->geom              = 0.e0;
  tm->leastSquareMatrix = 0.e0;
  tm->reord             = 0.e0;
/*... D1*/
  tm->solvD1            = 0.e0;
  tm->numeqD1           = 0.e0;
  tm->dataStructD1      = 0.e0;
  tm->CellPloadD1       = 0.e0;
  tm->CellTransientD1   = 0.e0;
  tm->systFormD1        = 0.e0;
  tm->rcGradD1          = 0.e0;
  tm->solvEdpD1         = 0.e0;
/*... T1*/
  tm->solvT1            = 0.e0;
  tm->numeqT1           = 0.e0;
  tm->dataStructT1      = 0.e0;
  tm->CellPloadT1       = 0.e0;
  tm->CellTransientT1   = 0.e0;
  tm->systFormT1        = 0.e0;
  tm->rcGradT1          = 0.e0;
  tm->solvEdpT1         = 0.e0;
/*... fluid*/
  tm->solvPres            = 0.e0;
  tm->solvVel             = 0.e0;
  tm->solvEnergy          = 0.e0;
  tm->numeqPres           = 0.e0;
  tm->numeqVel            = 0.e0;
  tm->numeqEnergy         = 0.e0;
  tm->dataStructVel       = 0.e0;
  tm->dataStructPres      = 0.e0;
  tm->dataStructEnergy    = 0.e0;
  tm->solvEdpFluid        = 0.e0;
  tm->cellPloadSimple     = 0.e0;
  tm->cellTransientSimple = 0.e0; 
  tm->systFormPres        = 0.e0;
  tm->systFormVel         = 0.e0;
  tm->systFormEnergy      = 0.e0;
  tm->velExp              = 0.e0;
  tm->rcGradPres          = 0.e0;
  tm->rcGradVel           = 0.e0;
  tm->rcGradEnergy        = 0.e0;
  tm->updateProp          = 0.e0;
  tm->residualSimple      = 0.e0;
/*... Combustion*/
  tm->tempFromTheEnergy  = 0.e0;
  tm->solvComb           = 0.e0;
  tm->numeqComb          = 0.e0;
  tm->dataStructComb     = 0.e0;
  tm->systFormComb       = 0.e0;
  tm->rcGradComb         = 0.e0;
  tm->rateReaction       = 0.e0;
  tm->timeChemical       = 0.e0;
  tm->heatRelease        = 0.e0;
  tm->speciesLoop        = 0.e0;
  tm->enthalpySpecies    = 0.e0;

/*... Blas*/
  tm->matVecOverHeadMpi = 0.e0;
  tm->matVecSparse      = 0.e0;
  tm->dot               = 0.e0;
  tm->dotOverHeadMpi    = 0.e0;
/*... Iterativos  */
  tm->pcg               = 0.e0;
  tm->pbicgstab         = 0.e0;
  tm->gmres             = 0.e0;
  tm->minres            = 0.e0;
/*... Direto*/
  tm->pardiso           = 0.e0;
/*... particionamento*/
  tm->partdMesh         = 0.e0;
  tm->partdMeshCom      = 0.e0;
/*...*/
  tm->overHeadCelMpi    = 0.e0;
  tm->overHeadNodMpi    = 0.e0;
  tm->overHeadNeqMpi    = 0.e0;
  tm->overHeadGCelMpi   = 0.e0;
  tm->overHeadGNodMpi   = 0.e0;
  tm->overHeadTotalMpi  = 0.e0;
/*...*/

/*...*/
  tm->turbulence        = 0.e0;
/*... precondicionador*/
  tm->precondDiag       = 0.e0;
  tm->total             = getTimeC();
/*...................................................................*/

}