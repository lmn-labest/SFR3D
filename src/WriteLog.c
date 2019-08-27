#include<WriteLog.h>


/********************************************************************* 
 * Data de criacao    : 21/08/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------* 
 * sist :                                                            * 
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
static void sist(Solv *solv,SistEq *sistEq
                ,short maxViz  ,FILE *file)
{
/*... tecnica de armazenamento*/
  if(sistEq->storage == CSR)
    fprintf(file,"Armazenamento      : CSR\n");
  else if(sistEq->storage == CSRD)
    fprintf(file,"Armazenamento      : CSRD\n");
  else if(sistEq->storage == CSRC)
    fprintf(file,"Armazenamento      : CSRC\n");
  else if(sistEq->storage == ELLPACK)
    fprintf(file,"Armazenamento      : ELLPACK\n");
  else if(sistEq->storage == CSRDCOO)
    fprintf(file,"Armazenamento      : CSRDCOO\n");
  else if(sistEq->storage == CSRCCOO)
    fprintf(file,"Armazenamento      : CSRCCOO\n");
  fprintf(file,"nEq                : %d\n",sistEq->neq);
  if(mpiVar.nPrcs > 1)
    fprintf(file,"nEqNov             : %d\n",sistEq->neqNov);
/*... ELLPACK*/
  if(sistEq->storage == ELLPACK){
    fprintf(file,"nAd                : %d\n",sistEq->nad);
    fprintf(file,"nAd Overhead       : %d\n"
           ,sistEq->neq*maxViz);
  }
/*... CSRD*/
  else{
    fprintf(file,"band Min           : %d\n"
           ,sistEq->bandCsr[BANDCSRMIN]);
    fprintf(file,"band Max           : %d\n"
           ,sistEq->bandCsr[BANDCSRMAX]);
    fprintf(file,"band Med           : %d\n"
           ,sistEq->bandCsr[BANDCSRMED]);
    fprintf(file,"nAd                : %d\n"
           ,sistEq->nad);
    if(mpiVar.nPrcs > 1)
      fprintf(file,"nAdR               : %d\n",sistEq->nadr);
  }
/*...................................................................*/
  fprintf(file,"tol                : %e\n",solv->tol);
/*... solver */
  if(solv->solver == PBICGSTAB )
    fprintf(file,"Iterativo          : PBICGSTAB\n");
  else if (solv->solver == PBICGSTABL2)
    fprintf(file, "Iterativo          : PBICGSTAB(2)\n");
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 03/10/2017                                   *
 * Data de modificaco : 21/08/2019                                   *
 *-------------------------------------------------------------------* 
 * WRITELOG : escrita do arquivo log de execuao                      * 
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
void writeLog(Mesh *mesh            ,Scheme *sc
             ,Solv *solvD1          ,SistEq *sistEqD1
             ,Solv *solvT1          ,SistEq *sistEqT1
             ,Solv *solvVel         ,SistEq *sistEqVel
             ,Solv *solvPres        ,SistEq *sistEqPres
             ,Solv *solvEnergy      ,SistEq *sistEqEnergy
             ,Solv *solvComb        ,SistEq *sistEqComb  
             ,Time *t                
             ,bool const fSolvD1    ,bool const fSolvT1
             ,bool const fSolvVel   ,bool const fSolvPres
             ,bool const fEnergy    ,bool const fTurb
             ,bool const fCombustion
             ,Omp omp
             ,char *nameIn          ,FILE *file)
{

  
  fprintf(file,"Log do execucao    : %s\n\n",nameIn); 
  
  if(mpiVar.nPrcs > 1)
  {
    fprintf(file,"myId               :%2d\n",mpiVar.myId); 
    fprintf(file,"nPrcs              :%2d\n\n",mpiVar.nPrcs); 
  } 
  
  fprintf(file,"Time:\n");
  fprintf(file,"%-25s : %13.3lf\n","adjcency",t->adjcency);
  fprintf(file,"%-25s : %13.3lf\n","geom"    ,t->geom);
  if( sc->rcGrad == RCLSQUARE )
    fprintf(file,"%-25s : %13.3lf\n","lSquareMatrix",t->leastSquareMatrix);
  fprintf(file,"%-25s : %13.3lf\n","reord",t->reord);

/*... Blas*/
  fprintf(file,"%-25s : %13.3lf\n","matVecSparse",t->matVecSparse);
  fprintf(file,"%-25s : %13.3lf\n","dot"         ,t->dot); 
/*... Blas overHead do mpi */
  if(mpiVar.nPrcs > 1)
  {
    fprintf(file,"%-25s : %13.3lf\n","matVecOverHead",t->matVecOverHeadMpi);
    fprintf(file,"%-25s : %13.3lf\n","dotOverHead"   ,t->dotOverHeadMpi);
  }

/*... Solver*/
  fprintf(file,"%-25s : %13.3lf\n","Pcg"        ,t->pcg);
  fprintf(file,"%-25s : %13.3lf\n","Pbicgstab"  ,t->pbicgstab);
  fprintf(file,"%-25s : %13.3lf\n","Gmres"      ,t->gmres);
  fprintf(file,"%-25s : %13.3lf\n","Pardiso"    ,t->pardiso);
  fprintf(file,"%-25s : %13.3lf\n","precondDiag",t->precondDiag);

/*... particionamento*/
  if(mpiVar.nPrcs > 1)
  {
    fprintf(file,"%-25s : %13.3lf\n","partdMesh"    ,t->partdMesh);
    fprintf(file,"%-25s : %13.3lf\n","partdMeshComm",t->partdMeshCom);
  }
/*... comunicacao entre as particoes*/
  if(mpiVar.nPrcs > 1)
  {
    fprintf(file,"%-25s : %13.3lf\n","overHeadCelMpi"   ,t->overHeadCelMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadNodMpi"   ,t->overHeadNodMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadNeqMpi"   ,t->overHeadNeqMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadGCelMpi"  ,t->overHeadGCelMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadGNodMpi"  ,t->overHeadGNodMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadMiscMpi"  ,t->overHeadMiscMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadBufferMpi",t->overHeadBufferMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadRecvMpi"  ,t->overHeadRecvMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadSendMpi"  ,t->overHeadSendMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadWaitMpi"  ,t->overHeadWaitMpi);
  }

/*... solvD1*/
  if(fSolvD1)
  {
    fprintf(file,"%-25s : %13.3lf\n","SolvD1"         ,t->solvD1);
    fprintf(file,"%-25s : %13.3lf\n","numeqD1"        ,t->numeqD1);
    fprintf(file,"%-25s : %13.3lf\n","dataStructD1"   ,t->dataStructD1);
    fprintf(file,"%-25s : %13.3lf\n","CellPloadD1"    ,t->CellPloadD1);
    fprintf(file,"%-25s : %13.3lf\n","CellTransientD1",t->CellTransientD1);
    fprintf(file,"%-25s : %13.3lf\n","systFormD1"     ,t->systFormD1);
    fprintf(file,"%-25s : %13.3lf\n","rcGradD1"       ,t->rcGradD1);
    fprintf(file,"%-25s : %13.3lf\n","solvEdpD1"      ,t->solvEdpD1);
  }
/*...................................................................*/

/*... solvT1*/
  if(fSolvT1)
  {
    fprintf(file,"%-25s : %13.3lf\n","SolvT1"         ,t->solvT1);
    fprintf(file,"%-25s : %13.3lf\n","numeqT1"        ,t->numeqT1);
    fprintf(file,"%-25s : %13.3lf\n","dataStructT1"   ,t->dataStructT1);
    fprintf(file,"%-25s : %13.3lf\n","CellPloadT1"    ,t->CellPloadT1);
    fprintf(file,"%-25s : %13.3lf\n","CellTransientT1",t->CellTransientT1);
    fprintf(file,"%-25s : %13.3lf\n","systFormT1"     ,t->systFormT1);
    fprintf(file,"%-25s : %13.3lf\n","rcGradT1"       ,t->rcGradT1);
    fprintf(file,"%-25s : %13.3lf\n","solvEdpT1"      ,t->solvEdpT1);
  }
/*...................................................................*/

/*...*/
  fprintf(file,"%-25s : %13.3lf\n","solvEdpFuild",t->solvEdpFluid);
  fprintf(file,"%-25s : %13.3lf\n","updateProp"         ,t->updateProp);

/*... solvVel*/
  if(fSolvVel)
  {
    fprintf(file,"%-25s : %13.3lf\n","SolvVel"            ,t->solvVel);
    fprintf(file,"%-25s : %13.3lf\n","numeqVel"           ,t->numeqVel);
    fprintf(file,"%-25s : %13.3lf\n","dataStructVel"      ,t->dataStructVel);
    fprintf(file,"%-25s : %13.3lf\n","CellPloadSimple"    ,t->cellPloadSimple);
    fprintf(file,"%-25s : %13.3lf\n","CellTransientSimple",t->cellTransientSimple);
    fprintf(file,"%-25s : %13.3lf\n","systFormVel"        ,t->systFormVel);
    fprintf(file,"%-25s : %13.3lf\n","velExp"             ,t->velExp);
    fprintf(file,"%-25s : %13.3lf\n","rcGradVel"          ,t->rcGradVel);
    fprintf(file,"%-25s : %13.3lf\n","residualSimple"     ,t->residualSimple);
  }
/*...................................................................*/

/*... solvPres*/
  if(fSolvPres)
  {
    fprintf(file,"%-25s : %13.3lf\n","SolvPres"      ,t->solvPres);
    fprintf(file,"%-25s : %13.3lf\n","numeqPres"     ,t->numeqPres);
    fprintf(file,"%-25s : %13.3lf\n","dataStructPres",t->dataStructPres);
    fprintf(file,"%-25s : %13.3lf\n","systFormPres"  ,t->systFormPres);
    fprintf(file,"%-25s : %13.3lf\n","rcGradPres"    ,t->rcGradPres);
  }
/*...................................................................*/

/*... energy*/
  if(fEnergy)
  {
    fprintf(file,"%-25s : %13.3lf\n","SolvEnergy"       ,t->solvEnergy);
    fprintf(file,"%-25s : %13.3lf\n","numeqEnergy"      ,t->numeqEnergy);
    fprintf(file,"%-25s : %13.3lf\n","dataStructEnergy" ,t->dataStructEnergy);
    fprintf(file,"%-25s : %13.3lf\n","systFormEnergy"   ,t->systFormEnergy);
    fprintf(file,"%-25s : %13.3lf\n","rcGradEnergy"     ,t->rcGradEnergy);
    fprintf(file,"%-25s : %13.3lf\n","tempFromTheEnergy",t->tempFromTheEnergy);
  }
/*...................................................................*/

/*... combustion*/
  if(fCombustion)
  {
    fprintf(file,"%-25s : %13.3lf\n","SolvCombustion"      ,t->solvComb);
    fprintf(file,"%-25s : %13.3lf\n","numeqCombustion"     ,t->numeqComb);
    fprintf(file,"%-25s : %13.3lf\n","dataStructCombustion",t->dataStructComb);
    fprintf(file,"%-25s : %13.3lf\n","systFormCombustion"  ,t->systFormComb);
    fprintf(file,"%-25s : %13.3lf\n","rcGradCombustion"    ,t->rcGradComb);
    fprintf(file,"%-25s : %13.3lf\n","rateReaction"        ,t->rateReaction);
    fprintf(file,"%-25s : %13.3lf\n","timeChemical"        ,t->timeChemical);
    fprintf(file,"%-25s : %13.3lf\n","heatRelease"         ,t->heatRelease);
    fprintf(file,"%-25s : %13.3lf\n","speciesLoop"         ,t->speciesLoop);
    fprintf(file,"%-25s : %13.3lf\n","enthalpySpecies "    ,t->enthalpySpecies);
  }
/*...................................................................*/

/*... turbulence*/
  if(fTurb)  
    fprintf(file,"%-25s : %13.3lf\n","turbulence",t->turbulence);
/*...................................................................*/

  fprintf(file,"%-25s : %13.3lf\n","Total",t->total);
/*...*/
  fprintf(file,"\nMesh:\n");
  fprintf(file,"nnode              : %d\n",mesh->nnode);
  fprintf(file,"nCell              : %d\n",mesh->numel);
  fprintf(file,"volume             : %lf m3 \n",mesh->mQuality.volume);
  fprintf(file,"Mass(0)            : %lf kg\n",mesh->mass[0]);
  fprintf(file,"Mass(1)            : %lf kg\n",mesh->mass[1]);
  fprintf(file,"Mass(2)            : %lf kg\n",mesh->mass[2]);
  fprintf(file,"non-OrthMed        : %.1lf°\n",mesh->mQuality.nonOrthMed);
  fprintf(file,"non-OtthMax        : %.1lf°\n",mesh->mQuality.nonOrthMax);
  fprintf(file,"skewMed            : %lf\n",mesh->mQuality.skewMed);
  fprintf(file,"skewMax            : %lf\n",mesh->mQuality.skewMax);
  fprintf(file,"aspect ratio max   : %lf\n",mesh->mQuality.aspectRaMax);
  fprintf(file,"aspect ratio min   : %lf\n",mesh->mQuality.aspectRaMin);

/*...*/
  if(fSolvD1)
  {
    fprintf(file,"\nSistema D1:\n");
    sist(solvD1,sistEqD1,mesh->maxViz,file);
  }
/*...................................................................*/

/*...*/
  if(fSolvT1)
  {
    fprintf(file,"\nSistema T1:\n");
    sist(solvT1,sistEqT1,mesh->maxViz,file);
  }
/*...................................................................*/

/*... solvVel*/
  if(fSolvVel){
    fprintf(file,"\nSistema Vel:\n");
    sist(solvVel,sistEqVel,mesh->maxViz,file);
  }
/*...................................................................*/

/*... solvPres*/
  if(fSolvPres)
  {
    fprintf(file,"\nSistema Pres:\n");
    sist(solvPres,sistEqPres,mesh->maxViz,file);
  }
/*...................................................................*/

/*... solvPres*/
  if(fEnergy)
  {
    fprintf(file,"\nSistema Energy:\n");
    sist(solvEnergy,sistEqEnergy,mesh->maxViz,file);
  }
/*...................................................................*/

/*... solvComb*/
  if(fCombustion)
  {
    fprintf(file,"\nSistema Combustion:\n");
    sist(solvComb,sistEqComb,mesh->maxViz,file);
  }
/*...................................................................*/

/*... OpenMp*/    
  fprintf(file, "\nOpenMp             :\n");
/*...*/
  fprintf(file,   "Solver             :\n");
  if(omp.fSolver)
    fprintf(file, "nThreads           : %d\n",omp.nThreadsSolver);
  else
    fprintf(file, "Disable\n");
/*...*/
  fprintf(file,   "Cell               :\n");
  if (omp.fCell)
    fprintf(file, "nThreads           : %d\n",omp.nThreadsCell);
  else
    fprintf(file, "Disable\n");
/*...*/
  fprintf(file,   "Update             :\n");
  if (omp.fUpdate)
    fprintf(file, "nThreads           : %d\n",omp.nThreadsUpdate);
  else
    fprintf(file, "Disable\n");
/*...*/
  fprintf(file,   "Grad               :\n");
  if (omp.fGrad)
    fprintf(file, "nThreads           : %d\n",omp.nThreadsGrad);
  else
    fprintf(file, "Disable\n");
/*...*/
  fprintf(file,   "Reaction           :\n");
  if (omp.fReaction)
    fprintf(file, "nThreads           : %d\n",omp.nThreadsReaction);
  else
    fprintf(file, "Disable\n");
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 03/10/2017                                   *
 * Data de modificaco : 21/08/2019                                   *
 *-------------------------------------------------------------------*
 * WRITELOGMEANTIME: escrita do arquivo log de execuao com os tempos * 
 * medios dp MPI                                                     * 
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
  void writeLogMeanTime(Mesh *mesh         ,Scheme *sc
             ,Solv *solvD1                 ,SistEq *sistEqD1
             ,Solv *solvT1                 ,SistEq *sistEqT1
             ,Solv *solvVel                ,SistEq *sistEqVel
             ,Solv *solvPres               ,SistEq *sistEqPres
             ,Solv *solvEnergy             ,SistEq *sistEqEnergy
             ,Solv *solvComb               ,SistEq *sistEqComb  
             ,Time *t                      ,Omp *omp 
             ,bool const fSolvD1           ,bool const fSolvT1
             ,bool const fSolvVel          ,bool const fSolvPres
             ,bool const fEnergy           ,bool const fTurb
             ,bool const fCombustion             
             ,char *nameIn                 ,FILE *file)
{

#ifdef _MPI_
  DOUBLE mTime;   
  DOUBLE  nPrcs= (DOUBLE) mpiVar.nPrcs;

  if(!mpiVar.myId){
    fprintf(file,"Log do execucao : %s\n\n",nameIn); 
    fprintf(file,"nPrcs           :%2d\n\n",mpiVar.nPrcs); 
    fprintf(file,"Time:\n");
  }

/*...*/
  MPI_Reduce(&t->leastSquareMatrix, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM              , 0     , mpiVar.comm);
  if( sc->rcGrad == RCLSQUARE ) 
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","lSquareMatrix",mTime/nPrcs);

  MPI_Reduce(&t->reord, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM ,0       , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","reord",mTime/nPrcs);

/*... Blas*/
  MPI_Reduce(&t->matVecSparse, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM       , 0       , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","matVecSparse",mTime/nPrcs);
  
  MPI_Reduce(&t->dot         , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM       , 0       , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","dot",mTime/nPrcs); 

/*... Blas overHead do mpi */
  MPI_Reduce(&t->matVecOverHeadMpi    , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM        , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","matVecOverHead",mTime/nPrcs);
  
  MPI_Reduce(&t->dotOverHeadMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","dotOverHead",mTime/nPrcs);

/*... Solver*/
  MPI_Reduce(&t->pcg           , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","Pcg",mTime/nPrcs);

  MPI_Reduce(&t->pbicgstab   , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","Pbicgstab",mTime/nPrcs);

  MPI_Reduce(&t->gmres         , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","Gmres",mTime/nPrcs);

  MPI_Reduce(&t->precondDiag   , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","precondDiag",mTime/nPrcs);

/*... particionamento*/
  MPI_Reduce(&t->partdMesh    , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","partdMesh",mTime/nPrcs);


  MPI_Reduce(&t->partdMeshCom , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","partdMeshCom",mTime/nPrcs);

/*... comunicacao entre as particoes*/
  MPI_Reduce(&t->overHeadCelMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","overHeadCelMpi",mTime/nPrcs);
  
  MPI_Reduce(&t->overHeadNodMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","overHeadNodMpi",mTime/nPrcs);
  
  MPI_Reduce(&t->overHeadNeqMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","overHeadNeqMpi",mTime/nPrcs);
  
  MPI_Reduce(&t->overHeadGCelMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","overHeadGCelMpi",mTime/nPrcs);

  MPI_Reduce(&t->overHeadGNodMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM           , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","overHeadGNodMpi",mTime/nPrcs);
  
  MPI_Reduce(&t->overHeadMiscMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM            , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","overHeadMisMpi",mTime/nPrcs);

  MPI_Reduce(&t->overHeadBufferMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM              , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","overHeadBufferMpi",mTime/nPrcs);

  MPI_Reduce(&t->overHeadRecvMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM              , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","overHeadRecvMpi",mTime/nPrcs);
 
  MPI_Reduce(&t->overHeadSendMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM              , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","overHeadSendMpi",mTime/nPrcs);

  MPI_Reduce(&t->overHeadWaitMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM            , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","overHeadWaitMpi",mTime/nPrcs);

/*... solvD1*/
  if(fSolvD1){
    MPI_Reduce(&t->solvD1          , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","SolvD1",mTime/nPrcs);
    
    MPI_Reduce(&t->numeqD1          , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","numeqD1",mTime/nPrcs);
  
    MPI_Reduce(&t->dataStructD1    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","dataStructD1",mTime/nPrcs);
    
    MPI_Reduce(&t->CellPloadD1     , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","CellPloadD1",mTime/nPrcs);
    
    MPI_Reduce(&t->CellTransientD1 , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","cellTransientD1",mTime/nPrcs);
    
    MPI_Reduce(&t->systFormD1      , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","systFormD1",mTime/nPrcs);
    
    MPI_Reduce(&t->rcGradD1        , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","rcGradD1",mTime/nPrcs);
    
    MPI_Reduce(&t->solvEdpD1       , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","solvEdpD1",mTime/nPrcs);
  }
/*...................................................................*/

/*... solvT1*/
  if(fSolvT1){
    MPI_Reduce(&t->solvT1          , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","SolvT1",mTime/nPrcs);
    
    MPI_Reduce(&t->numeqT1          , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","numeqT1",mTime/nPrcs);
  
    MPI_Reduce(&t->dataStructT1    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","dataStructT1",mTime/nPrcs);
    
    MPI_Reduce(&t->CellPloadT1     , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","cellPloadT1",mTime/nPrcs);
    
    MPI_Reduce(&t->CellTransientT1 , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","cellTransientT1",mTime/nPrcs);
    
    MPI_Reduce(&t->systFormT1      , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","systFormT1",mTime/nPrcs);
    
    MPI_Reduce(&t->rcGradT1        , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","rcGradT1",mTime/nPrcs);
    
    MPI_Reduce(&t->solvEdpT1       , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","solvEdpT1",mTime/nPrcs);
  }
/*...................................................................*/
  
  MPI_Reduce(&t->solvEdpFluid    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","solvEdpFluid",mTime/nPrcs);

  MPI_Reduce(&t->updateProp       , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","updateProp ",mTime/nPrcs);


/*... solver fluid*/
  if(fSolvVel)
  {
    MPI_Reduce(&t->solvVel         , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","solvVel",mTime/nPrcs);

    MPI_Reduce(&t->numeqVel        , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","numeqVel",mTime/nPrcs);
    
    MPI_Reduce(&t->dataStructVel   , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","dataStructVel",mTime/nPrcs);

    MPI_Reduce(&t->cellPloadSimple, &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","cellStructVel",mTime/nPrcs);

    MPI_Reduce(&t->cellTransientSimple, &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","cellTransientSimple",mTime/nPrcs);

    MPI_Reduce(&t->systFormVel    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","systFormaVeltSimple",mTime/nPrcs);

    MPI_Reduce(&t->rcGradVel      , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","rcGraVel",mTime/nPrcs);

  }
/*...................................................................*/

/*... solver pres*/
  if(fSolvPres)
  {
    MPI_Reduce(&t->solvPres    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","SolvPres",mTime/nPrcs);

    MPI_Reduce(&t->numeqPres         , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","numeqPres",mTime/nPrcs);

    MPI_Reduce(&t->dataStructPres , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","dataStructPres",mTime/nPrcs);
    
    MPI_Reduce(&t->systFormPres   , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","systFormPres",mTime/nPrcs);

    MPI_Reduce(&t->rcGradPres     , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","rcGradPres",mTime/nPrcs);
  }
/*...................................................................*/

/*... solver fluid*/
  if(fEnergy)
  {
    MPI_Reduce(&t->solvEnergy         , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","solvEnergy",mTime/nPrcs);

    MPI_Reduce(&t->numeqEnergy        , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","numeqEnergy",mTime/nPrcs);
    
    MPI_Reduce(&t->dataStructEnergy      , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","dataStructEnergy",mTime/nPrcs);

    MPI_Reduce(&t->systFormEnergy, &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","systFormEnergy",mTime/nPrcs);

    MPI_Reduce(&t->rcGradEnergy , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","rcGradEnergy",mTime/nPrcs);

    MPI_Reduce(&t->tempFromTheEnergy    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","tempFromTheEnergy",mTime/nPrcs);
  }
/*...................................................................*/

/*... solver fluid*/
  if(fCombustion)
  {
    MPI_Reduce(&t->solvComb       , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","solvCombustion",mTime/nPrcs);

    MPI_Reduce(&t->numeqComb      , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","numeqCombustion",mTime/nPrcs);
    
    MPI_Reduce(&t->dataStructComb , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","dataStructCombustion",mTime/nPrcs);

    MPI_Reduce(&t->systFormComb    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","systFormCombustion",mTime/nPrcs);

    MPI_Reduce(&t->rcGradComb , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","rcGradCombustion",mTime/nPrcs);

    MPI_Reduce(&t->rateReaction    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","rateReaction",mTime/nPrcs);

    MPI_Reduce(&t->timeChemical    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","timeChemical",mTime/nPrcs);
    
    MPI_Reduce(&t->heatRelease   , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","heatRelease",mTime/nPrcs);

    MPI_Reduce(&t->speciesLoop    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","speciesLoop",mTime/nPrcs);

    MPI_Reduce(&t->enthalpySpecies   , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"%-25s : %13.3lf\n","enthalpySpecies",mTime/nPrcs);

  }
/*...................................................................*/

  MPI_Reduce(&t->total           , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM            , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"%-25s : %13.3lf\n","Total",mTime/nPrcs);


/*...*/
  if(!mpiVar.myId){
    fprintf(file,"\nMesh:\n");
    fprintf(file,"nnode           : %d\n",mesh->nnode);
    fprintf(file,"nCell           : %d\n",mesh->numel);
    fprintf(file,"volume          : %lf\n",mesh->mQuality.volume);
    fprintf(file,"non-OrthMed     : %.2lfÂ°\n",mesh->mQuality.nonOrthMed);
    fprintf(file,"non-OtthMax     : %.2lfÂ°\n",mesh->mQuality.nonOrthMax);
    fprintf(file,"skewMed         : %lf\n",mesh->mQuality.skewMed);
    fprintf(file,"skewMax         : %lf\n",mesh->mQuality.skewMax);
  }
/*...................................................................*/
#endif

}
/*********************************************************************/ 
