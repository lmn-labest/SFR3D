#include<WriteLog.h>
/********************************************************************* 
 * Data de criacao    : 03/10/2017                                   *
 * Data de modificaco : 08/06/2018                                   *
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
void writeLog(Mesh mesh             ,Scheme sc
             ,Solv *solvD1          ,SistEq *sistEqD1
             ,Solv *solvT1          ,SistEq *sistEqT1
             ,Solv *solvVel         ,SistEq *sistEqVel
             ,Solv *solvPres        ,SistEq *sistEqPres
             ,Time t                
             ,bool const fSolvD1    ,bool const fSolvT1
             ,bool const fSolvVel   ,bool const fSolvPres
             ,bool const fEnergy    ,bool const fTurb
             ,bool const fCombustion
             ,Omp omp
             ,char *nameIn          ,FILE *file){

  
  fprintf(file,"Log do execucao    : %s\n\n",nameIn); 
  
  if(mpiVar.nPrcs > 1){
    fprintf(file,"myId               :%2d\n",mpiVar.myId); 
    fprintf(file,"nPrcs              :%2d\n\n",mpiVar.nPrcs); 
  } 
  
  fprintf(file,"Time:\n");
  fprintf(file,"%-25s : %13.3lf\n","adjcency",t.adjcency);
  fprintf(file,"%-25s : %13.3lf\n","geom"    ,t.geom);
  if( sc.rcGrad == RCLSQUARE )
    fprintf(file,"%-25s : %13.3lf\n","lSquareMatrix",t.leastSquareMatrix);
  fprintf(file,"%-25s : %13.3lf\n","reord",t.reord);

/*... Blas*/
  fprintf(file,"%-25s : %13.3lf\n","matVecSparse",t.matVecSparse);
  fprintf(file,"%-25s : %13.3lf\n","dot"         ,t.dot); 
/*... Blas overHead do mpi */
  if(mpiVar.nPrcs > 1){
    fprintf(file,"%-25s : %13.3lf\n","matVecOverHead",t.matVecOverHeadMpi);
    fprintf(file,"%-25s : %13.3lf\n","dotOverHead"   ,t.dotOverHeadMpi);
  }

/*... Solver*/
  fprintf(file,"%-25s : %13.3lf\n","Pcg"        ,t.pcg);
  fprintf(file,"%-25s : %13.3lf\n","Pbicgstab"  ,t.pbicgstab);
  fprintf(file,"%-25s : %13.3lf\n","Gmres"      ,t.gmres);
  fprintf(file,"%-25s : %13.3lf\n","Pardiso"    ,t.pardiso);
  fprintf(file,"%-25s : %13.3lf\n","precondDiag",t.precondDiag);

/*... particionamento*/
  if(mpiVar.nPrcs > 1){
    fprintf(file,"%-25s : %13.3lf\n","partdMesh"    ,t.partdMesh);
    fprintf(file,"%-25s : %13.3lf\n","partdMeshComm",t.partdMeshCom);
  }
/*... comunicacao entre as particoes*/
  if(mpiVar.nPrcs > 1){
    fprintf(file,"%-25s : %13.3lf\n","overHeadCelMpi"  ,t.overHeadCelMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadNodMpi"  ,t.overHeadNodMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadNeqMpi"  ,t.overHeadNeqMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadGCelMpi" ,t.overHeadGCelMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadGNodMpi" ,t.overHeadGNodMpi);
    fprintf(file,"%-25s : %13.3lf\n","overHeadTotalMpi",t.overHeadTotalMpi);
  }

/*... solvD1*/
  if(fSolvD1){
    fprintf(file,"%-25s : %13.3lf\n","SolvD1"         ,t.solvD1);
    fprintf(file,"%-25s : %13.3lf\n","numeqD1"        ,t.numeqD1);
    fprintf(file,"%-25s : %13.3lf\n","dataStructD1"   ,t.dataStructD1);
    fprintf(file,"%-25s : %13.3lf\n","CellPloadD1"    ,t.CellPloadD1);
    fprintf(file,"%-25s : %13.3lf\n","CellTransientD1",t.CellTransientD1);
    fprintf(file,"%-25s : %13.3lf\n","systFormD1"     ,t.systFormD1);
    fprintf(file,"%-25s : %13.3lf\n","rcGradD1"       ,t.rcGradD1);
    fprintf(file,"%-25s : %13.3lf\n","solvEdpD1"      ,t.solvEdpD1);
  }
/*...................................................................*/

/*... solvT1*/
  if(fSolvT1){
    fprintf(file,"%-25s : %13.3lf\n","SolvT1"         ,t.solvT1);
    fprintf(file,"%-25s : %13.3lf\n","numeqT1"        ,t.numeqT1);
    fprintf(file,"%-25s : %13.3lf\n","dataStructT1"   ,t.dataStructT1);
    fprintf(file,"%-25s : %13.3lf\n","CellPloadT1"    ,t.CellPloadT1);
    fprintf(file,"%-25s : %13.3lf\n","CellTransientT1",t.CellTransientT1);
    fprintf(file,"%-25s : %13.3lf\n","systFormT1"     ,t.systFormT1);
    fprintf(file,"%-25s : %13.3lf\n","rcGradT1"       ,t.rcGradT1);
    fprintf(file,"%-25s : %13.3lf\n","solvEdpT1"      ,t.solvEdpT1);
  }
/*...................................................................*/

/*...*/
  fprintf(file,"%-25s : %13.3lf\n","solvEdpFuild",t.solvEdpFluid);

/*... solvVel*/
  if(fSolvVel){
    fprintf(file,"%-25s : %13.3lf\n","SolvVel"            ,t.solvVel);
    fprintf(file,"%-25s : %13.3lf\n","numeqVel"           ,t.numeqVel);
    fprintf(file,"%-25s : %13.3lf\n","dataStructVel"      ,t.dataStructVel);
    fprintf(file,"%-25s : %13.3lf\n","CellPloadSimple"    ,t.cellPloadSimple);
    fprintf(file,"%-25s : %13.3lf\n","CellTransientSimple",t.cellTransientSimple);
    fprintf(file,"%-25s : %13.3lf\n","systFormVel"        ,t.systFormVel);
    fprintf(file,"%-25s : %13.3lf\n","velExp"             ,t.velExp);
    fprintf(file,"%-25s : %13.3lf\n","rcGradVel"          ,t.rcGradVel);
    fprintf(file,"%-25s : %13.3lf\n","updateProp"         ,t.updateProp);
    fprintf(file,"%-25s : %13.3lf\n","residualSimple"     ,t.residualSimple);
  }
/*...................................................................*/

/*... energy*/
  if(fEnergy){
    fprintf(file,"%-25s : %13.3lf\n","SolvEnergy"       ,t.solvEnergy);
    fprintf(file,"%-25s : %13.3lf\n","numeqEnergy"      ,t.numeqEnergy);
    fprintf(file,"%-25s : %13.3lf\n","dataStructEnergy" ,t.dataStructEnergy);
    fprintf(file,"%-25s : %13.3lf\n","systFormEnergy"   ,t.systFormEnergy);
    fprintf(file,"%-25s : %13.3lf\n","rcGradEnergy"     ,t.rcGradEnergy);
    fprintf(file,"%-25s : %13.3lf\n","tempFromTheEnergy",t.tempFromTheEnergy);
  }
/*...................................................................*/

/*... combustion*/
  if(fCombustion){
    fprintf(file,"%-25s : %13.3lf\n","SolvCombustion"      ,t.solvComb);
    fprintf(file,"%-25s : %13.3lf\n","numeqCombustion"     ,t.numeqComb);
    fprintf(file,"%-25s : %13.3lf\n","dataStructCombustion",t.dataStructComb);
    fprintf(file,"%-25s : %13.3lf\n","systFormCombustion"  ,t.systFormComb);
    fprintf(file,"%-25s : %13.3lf\n","rcGradCombustion"    ,t.rcGradComb);
    fprintf(file,"%-25s : %13.3lf\n","rateReaction"        ,t.rateReaction);
    fprintf(file,"%-25s : %13.3lf\n","timeChemical"        ,t.timeChemical);
    fprintf(file,"%-25s : %13.3lf\n","heatRelease"         ,t.heatRelease);
    fprintf(file,"%-25s : %13.3lf\n","speciesLoop"         ,t.speciesLoop);
    fprintf(file,"%-25s : %13.3lf\n","enthalpySpecies "    ,t.enthalpySpecies);
  }
/*...................................................................*/

/*... solvPres*/
  if(fSolvPres){
    fprintf(file,"%-25s : %13.3lf\n","SolvPres"      ,t.solvPres);
    fprintf(file,"%-25s : %13.3lf\n","numeqPres"     ,t.numeqPres);
    fprintf(file,"%-25s : %13.3lf\n","dataStructPres",t.dataStructPres);
    fprintf(file,"%-25s : %13.3lf\n","systFormPres"  ,t.systFormPres);
    fprintf(file,"%-25s : %13.3lf\n","rcGradPres"    ,t.rcGradPres);
  }
/*...................................................................*/

/*... turbulence*/
  if(fTurb)  
    fprintf(file,"%-25s : %13.3lf\n","turbulence",t.turbulence);
/*...................................................................*/

  fprintf(file,"%-25s : %13.3lf\n","Total",t.total);
/*...*/
  fprintf(file,"\nMesh:\n");
  fprintf(file,"nnode              : %d\n",mesh.nnode);
  fprintf(file,"nCell              : %d\n",mesh.numel);
  fprintf(file,"volume             : %lf m3 \n",mesh.mQuality.volume);
  fprintf(file,"Mass(0)            : %lf kg\n",mesh.mass[0]);
  fprintf(file,"Mass(1)            : %lf kg\n",mesh.mass[1]);
  fprintf(file,"Mass(2)            : %lf kg\n",mesh.mass[2]);
  fprintf(file,"non-OrthMed        : %.1lf°\n",mesh.mQuality.nonOrthMed);
  fprintf(file,"non-OtthMax        : %.1lf°\n",mesh.mQuality.nonOrthMax);
  fprintf(file,"skewMed            : %lf\n",mesh.mQuality.skewMed);
  fprintf(file,"skewMax            : %lf\n",mesh.mQuality.skewMax);
  fprintf(file,"aspect ratio max   : %lf\n",mesh.mQuality.aspectRaMax);
  fprintf(file,"aspect ratio min   : %lf\n",mesh.mQuality.aspectRaMin);

/*...*/
  if(fSolvD1){
    fprintf(file,"\nSistema D1:\n");
/*... tecnica de armazenamento*/
    if(sistEqD1->storage == CSR)
      fprintf(file,"Armazenamento      : CSR\n");
    else if(sistEqD1->storage == CSRD)
      fprintf(file,"Armazenamento      : CSRD\n");
    else if(sistEqD1->storage == CSRC)
      fprintf(file,"Armazenamento      : CSRC\n");
    else if(sistEqD1->storage == ELLPACK)
      fprintf(file,"Armazenamento      : ELLPACK\n");
    else if(sistEqD1->storage == CSRDCOO)
      fprintf(file,"Armazenamento      : CSRDCOO\n");
    else if(sistEqT1->storage == CSRCCOO)
      fprintf(file,"Armazenamento      : CSRCCOO\n");
    fprintf(file,"nEq               : %d\n",sistEqD1->neq);
    if(mpiVar.nPrcs > 1)
      fprintf(file,"nEqNov            : %d\n",sistEqD1->neqNov);
/*... ELLPACK*/
    if(sistEqD1->storage == ELLPACK){
      fprintf(file,"nAd               : %d\n",sistEqD1->nad);
      fprintf(file,"nAd Overhead      : %d\n",sistEqD1->neq*mesh.maxViz);
    }
/*... CSRD*/
    else{
      fprintf(file,"band Min          : %d\n"
             ,sistEqD1->bandCsr[BANDCSRMIN]);
      fprintf(file,"band Max          : %d\n"
             ,sistEqD1->bandCsr[BANDCSRMAX]);
      fprintf(file,"band Med          : %d\n"
             ,sistEqD1->bandCsr[BANDCSRMED]);
      fprintf(file,"nAd               : %d\n",sistEqD1->nad);
      if(mpiVar.nPrcs > 1)
        fprintf(file,"nAdR              : %d\n",sistEqD1->nadr);
    }
/*...................................................................*/
    fprintf(file,"tol               : %e\n",solvD1->tol);
/*... solver */
    if(solvD1->solver == PCG)
      fprintf(file,"Iterativo         : PCG\n");
    else if(solvD1->solver == PBICGSTAB )
      fprintf(file,"Iterativo         : PBICGSTAB\n");
  }
/*...................................................................*/

/*...*/
  if(fSolvT1){
    fprintf(file,"\nSistema T1:\n");
/*... tecnica de armazenamento*/
    if(sistEqT1->storage == CSR)
      fprintf(file,"Armazenamento      : CSR\n");
    else if(sistEqT1->storage == CSRD)
      fprintf(file,"Armazenamento      : CSRD\n");
    else if(sistEqT1->storage == CSRC)
      fprintf(file,"Armazenamento      : CSRC\n");
    else if(sistEqT1->storage == ELLPACK)
      fprintf(file,"Armazenamento      : ELLPACK\n");
    else if(sistEqT1->storage == CSRDCOO)
      fprintf(file,"Armazenamento      : CSRDCOO\n");
    else if(sistEqT1->storage == CSRCCOO)
      fprintf(file,"Armazenamento      : CSRCCOO\n");
    fprintf(file,"nEq                : %d\n",sistEqT1->neq);
    if(mpiVar.nPrcs > 1)
      fprintf(file,"nEqNov             : %d\n",sistEqT1->neqNov);
/*... ELLPACK*/
    if(sistEqT1->storage == ELLPACK){
      fprintf(file,"nAd                : %d\n",sistEqT1->nad);
      fprintf(file,"nAd Overhead       : %d\n",sistEqT1->neq*mesh.maxViz);
    }
/*... CSRD*/
    else{
      fprintf(file,"band Min           : %d\n"
             ,sistEqT1->bandCsr[BANDCSRMIN]);
      fprintf(file,"band Max           : %d\n"
             ,sistEqT1->bandCsr[BANDCSRMAX]);
      fprintf(file,"band Med           : %d\n"
             ,sistEqT1->bandCsr[BANDCSRMED]);
      fprintf(file,"nAd                : %d\n"
             ,sistEqT1->nad);
      if(mpiVar.nPrcs > 1)
        fprintf(file,"nAdR              : %d\n",sistEqT1->nadr);
    }
/*...................................................................*/
    fprintf(file,"tol                : %e\n",solvT1->tol);
/*... solver */
    if(solvT1->solver == PBICGSTAB )
      fprintf(file,"Iterativo          : PBICGSTAB\n");
    else if (solvT1->solver == PBICGSTABL2)
      fprintf(file, "Iterativo          : PBICGSTAB(2)\n");
  }
/*...................................................................*/

/*... solvVel*/
  if(fSolvVel){
    fprintf(file,"\nSistema Vel:\n");
/*... tecnica de armazenamento*/
    if(sistEqVel->storage == CSR)
      fprintf(file,"Armazenamento      : CSR\n");
    else if(sistEqVel->storage == CSRD)
      fprintf(file,"Armazenamento      : CSRD\n");
    else if(sistEqVel->storage == CSRC)
      fprintf(file,"Armazenamento      : CSRC\n");
    else if(sistEqVel->storage == ELLPACK)
      fprintf(file,"Armazenamento      : ELLPACK\n");
    else if(sistEqVel->storage == CSRDCOO)
      fprintf(file,"Armazenamento      : CSRDCOO\n");
    else if(sistEqVel->storage == CSRCCOO)
      fprintf(file,"Armazenamento      : CSRCCOO\n");
    fprintf(file,"nEq                : %d\n",sistEqVel->neq);
    if(mpiVar.nPrcs > 1)
      fprintf(file,"nEqNov             : %d\n",sistEqVel->neqNov);
/*... ELLPACK*/
    if(sistEqVel->storage == ELLPACK){
      fprintf(file,"nAd                : %d\n",sistEqVel->nad);
      fprintf(file,"nAd Overhead       : %d\n"
             ,sistEqVel->neq*mesh.maxViz);
    }
/*... CSRD*/
    else{
      fprintf(file,"band Min           : %d\n"
             ,sistEqVel->bandCsr[BANDCSRMIN]);
      fprintf(file,"band Max           : %d\n"
             ,sistEqVel->bandCsr[BANDCSRMAX]);
      fprintf(file,"band Med           : %d\n"
             ,sistEqVel->bandCsr[BANDCSRMED]);
      fprintf(file,"nAd                : %d\n"
             ,sistEqVel->nad);
      if(mpiVar.nPrcs > 1)
        fprintf(file,"nAdR              : %d\n",sistEqVel->nadr);
    }
/*...................................................................*/
    fprintf(file,"tol                : %e\n",solvVel->tol);
/*... solver */
    if(solvVel->solver == PBICGSTAB )
      fprintf(file,"Iterativo          : PBICGSTAB\n");
    else if (solvVel->solver == PBICGSTABL2)
      fprintf(file, "Iterativo          : PBICGSTAB(2)\n");
  }
/*...................................................................*/

/*... solvPres*/
  if(fSolvPres){
    fprintf(file,"\nSistema Pres:\n");
/*... tecnica de armazenamento*/
    if(sistEqPres->storage == CSR)
      fprintf(file,"Armazenamento      : CSR\n");
    else if(sistEqPres->storage == CSRD)
      fprintf(file,"Armazenamento      : CSRD\n");
    else if(sistEqPres->storage == CSRC)
      fprintf(file,"Armazenamento      : CSRC\n");
    else if(sistEqPres->storage == ELLPACK)
      fprintf(file,"Armazenamento      : ELLPACK\n");
    else if(sistEqPres->storage == CSRDCOO)
      fprintf(file,"Armazenamento      : CSRDCOO\n");
    else if(sistEqPres->storage == CSRCCOO)
      fprintf(file,"Armazenamento      : CSRCCOO\n");
    fprintf(file,"nEq                : %d\n",sistEqPres->neq);
    if(mpiVar.nPrcs > 1)
      fprintf(file,"nEqNov             : %d\n",sistEqPres->neqNov);
/*... ELLPACK*/
    if(sistEqPres->storage == ELLPACK){
      fprintf(file,"nAd                : %d\n",sistEqPres->nad);
      fprintf(file,"nAd Overhead       : %d\n"
             ,sistEqPres->neq*mesh.maxViz);
    }
/*... CSRD*/
    else{
      fprintf(file,"band Min           : %d\n"
             ,sistEqPres->bandCsr[BANDCSRMIN]);
      fprintf(file,"band Max           : %d\n"
             ,sistEqPres->bandCsr[BANDCSRMAX]);
      fprintf(file,"band Med           : %d\n"
             ,sistEqPres->bandCsr[BANDCSRMED]);
      fprintf(file,"nAd                : %d\n"
             ,sistEqPres->nad);
      if(mpiVar.nPrcs > 1)
        fprintf(file,"nAdR              : %d\n",sistEqPres->nadr);
    }
/*...................................................................*/
    fprintf(file,"tol                : %e\n",solvPres->tol);
/*... solver */
    if(solvPres->solver == PCG)
      fprintf(file,"Iterativo          : PCG\n");
    else if(solvPres->solver == PBICGSTAB )
      fprintf(file,"Iterativo          : PBICGSTAB\n");
    else if (solvPres->solver == PBICGSTABL2)
      fprintf(file, "Iterativo          : PBICGSTAB(2)\n");
    else if (solvPres->solver == PARDISO)
      fprintf(file, "Ditero             : PARDISO\n");
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
void writeLogMeanTime(Mesh mesh         ,Scheme sc
             ,Solv *solvD1      ,SistEq *sistEqD1
             ,Solv *solvT1      ,SistEq *sistEqT1
             ,Time t
             ,bool const fSolvD1,bool const fSolvT1
             ,char *nameIn      ,FILE *file){

#ifdef _MPI_
  DOUBLE mTime;   
  DOUBLE  nPrcs= (DOUBLE) mpiVar.nPrcs;

  if(!mpiVar.myId){
    fprintf(file,"Log do execucao : %s\n\n",nameIn); 
    fprintf(file,"nPrcs           :%2d\n\n",mpiVar.nPrcs); 
    fprintf(file,"Time:\n");
  }

/*...*/
  MPI_Reduce(&t.leastSquareMatrix, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM              , 0     , mpiVar.comm);
  if( sc.rcGrad == RCLSQUARE ) 
    if(!mpiVar.myId)
      fprintf(file,"lSquareMatrix   : %lf\n",mTime/nPrcs);

  MPI_Reduce(&t.reord, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM ,0       , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"reord           : %lf\n",mTime/nPrcs);

/*... Blas*/
  MPI_Reduce(&t.matVecSparse, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM       , 0       , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"matVecSparse    : %lf\n",mTime/nPrcs);
  
  MPI_Reduce(&t.dot         , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM       , 0       , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"dot             : %lf\n",mTime/nPrcs); 

/*... Blas overHead do mpi */
  MPI_Reduce(&t.matVecOverHeadMpi    , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM        , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"matVecOverHead  : %lf\n",mTime/nPrcs);
  
  MPI_Reduce(&t.dotOverHeadMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"dotOverHead     : %lf\n",mTime/nPrcs);

/*... Solver*/
  MPI_Reduce(&t.pcg           , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"Pcg             : %lf\n",mTime/nPrcs);

  MPI_Reduce(&t.pbicgstab   , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"Pbicgstab       : %lf\n",mTime/nPrcs);

  MPI_Reduce(&t.precondDiag   , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"precondDiag     : %lf\n",mTime/nPrcs);


/*... comunicacao entre as particoes*/
  MPI_Reduce(&t.overHeadCelMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"overHeadCelMpi  : %lf\n",mTime/nPrcs);
  
  MPI_Reduce(&t.overHeadNodMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"overHeadNodMpi  : %lf\n",mTime/nPrcs);
  
  MPI_Reduce(&t.overHeadNeqMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"overHeadNeqMpi  : %lf\n",mTime/nPrcs);
  
  MPI_Reduce(&t.overHeadGCelMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM          , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"overHeadGCelMpi : %lf\n",mTime/nPrcs);

  MPI_Reduce(&t.overHeadGNodMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM           , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"overHeadGNodMpi : %lf\n",mTime/nPrcs);
  
  MPI_Reduce(&t.overHeadTotalMpi, &mTime , 1, MPI_DOUBLE
            ,MPI_SUM            , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"overHeadTotalMpi: %lf\n",mTime/nPrcs);

/*... solvD1*/
  if(fSolvD1){
    MPI_Reduce(&t.solvD1          , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"SolvD1          : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.numeqD1          , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"numeqD1         : %lf\n",mTime/nPrcs);
  
    MPI_Reduce(&t.dataStructD1    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"dataStructD1    : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.CellPloadD1     , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"CellPloadD1     : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.CellTransientD1 , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"CellTransientD1 : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.systFormD1      , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"systFormD1      : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.rcGradD1        , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"rcGradD1        : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.solvEdpD1       , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"solvEdpD1       : %lf\n",mTime/nPrcs);
  }
/*...................................................................*/

/*... solvT1*/
  if(fSolvT1){
    MPI_Reduce(&t.solvT1          , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"SolvT1          : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.numeqT1          , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"numeqT1         : %lf\n",mTime/nPrcs);
  
    MPI_Reduce(&t.dataStructT1    , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"dataStructT1    : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.CellPloadT1     , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"CellPloadT1     : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.CellTransientT1 , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"CellTransientT1 : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.systFormT1      , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"systFormT1      : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.rcGradT1        , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"rcGradT1        : %lf\n",mTime/nPrcs);
    
    MPI_Reduce(&t.solvEdpT1       , &mTime , 1, MPI_DOUBLE
              ,MPI_SUM            , 0      , mpiVar.comm);
    if(!mpiVar.myId)
      fprintf(file,"solvEdpT1       : %lf\n",mTime/nPrcs);
  }
/*...................................................................*/
  
  MPI_Reduce(&t.total           , &mTime , 1, MPI_DOUBLE
            ,MPI_SUM            , 0      , mpiVar.comm);
  if(!mpiVar.myId)
    fprintf(file,"Total           : %lf\n",mTime/nPrcs);


/*...*/
  if(!mpiVar.myId){
    fprintf(file,"\nMalha:\n");
    fprintf(file,"nnode           : %d\n",mesh.nnode);
    fprintf(file,"nCell           : %d\n",mesh.numel);
    fprintf(file,"volume          : %lf\n",mesh.mQuality.volume);
    fprintf(file,"non-OrthMed     : %.2lfÂ°\n",mesh.mQuality.nonOrthMed);
    fprintf(file,"non-OtthMax     : %.2lfÂ°\n",mesh.mQuality.nonOrthMax);
    fprintf(file,"skewMed         : %lf\n",mesh.mQuality.skewMed);
    fprintf(file,"skewMax         : %lf\n",mesh.mQuality.skewMax);
  }
/*...................................................................*/
#endif

}
/*********************************************************************/ 
