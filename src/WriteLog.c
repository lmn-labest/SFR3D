#include<WriteLog.h>
/********************************************************************* 
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
void writeLog(Mesh mesh          ,Scheme sc
             ,Solv *solvD1       ,SistEq *sistEqD1
             ,Solv *solvT1       ,SistEq *sistEqT1
             ,Solv *solvVel      ,SistEq *sistEqVel
             ,Solv *solvPres     ,SistEq *sistEqPres
             ,Time t
             ,bool const fSolvD1 ,bool const fSolvT1
             ,bool const fSolvVel,bool const fSolvPres
             ,Omp omp
             ,char *nameIn       ,FILE *file){

  
  fprintf(file,"Log do execucao    : %s\n\n",nameIn); 
  
  if(mpiVar.nPrcs > 1){
    fprintf(file,"myId               :%2d\n",mpiVar.myId); 
    fprintf(file,"nPrcs              :%2d\n\n",mpiVar.nPrcs); 
  } 
  
  fprintf(file,"Time:\n");
  fprintf(file,"adjcency           : %lf\n",t.adjcency);
  fprintf(file,"geom               : %lf\n",t.geom);
  if( sc.rcGrad == RCLSQUARE )
    fprintf(file,"lSquareMatrix      : %lf\n",t.leastSquareMatrix);
  fprintf(file,"reord              : %lf\n",t.reord);

/*... Blas*/
  fprintf(file,"matVecSparse       : %lf\n",t.matVecSparse);
  fprintf(file,"dot                : %lf\n",t.dot); 
/*... Blas overHead do mpi */
  if(mpiVar.nPrcs > 1){
    fprintf(file,"matVecOverHead     : %lf\n",t.matVecOverHeadMpi);
    fprintf(file,"dotOverHead        : %lf\n",t.dotOverHeadMpi);
  }

/*... Solver*/
  fprintf(file,"Pcg                : %lf\n",t.pcg);
  fprintf(file,"Pbicgstab          : %lf\n",t.pbicgstab);
  fprintf(file,"precondDiag        : %lf\n",t.precondDiag);

/*... particionamento*/
  if(mpiVar.nPrcs > 1){
    fprintf(file,"partdMesh          : %lf\n",t.partdMesh);
    fprintf(file,"partdMeshComm      : %lf\n",t.partdMeshCom);
  }
/*... comunicacao entre as particoes*/
  if(mpiVar.nPrcs > 1){
    fprintf(file,"overHeadCelMpi     : %lf\n",t.overHeadCelMpi);
    fprintf(file,"overHeadNodMpi     : %lf\n",t.overHeadNodMpi);
    fprintf(file,"overHeadNeqMpi     : %lf\n",t.overHeadNeqMpi);
    fprintf(file,"overHeadGCelMpi    : %lf\n",t.overHeadGCelMpi);
    fprintf(file,"overHeadGNodMpi    : %lf\n",t.overHeadGNodMpi);
    fprintf(file,"overHeadTotalMpi   : %lf\n",t.overHeadTotalMpi);
  }


/*... solvD1*/
  if(fSolvD1){
    fprintf(file,"SolvD1             : %lf\n",t.solvD1);
    fprintf(file,"numeqD1            : %lf\n",t.numeqD1);
    fprintf(file,"dataStructD1       : %lf\n",t.dataStructD1);
    fprintf(file,"CellPloadD1        : %lf\n",t.CellPloadD1);
    fprintf(file,"CellTransientD1    : %lf\n",t.CellTransientD1);
    fprintf(file,"systFormD1         : %lf\n",t.systFormD1);
    fprintf(file,"rcGradD1           : %lf\n",t.rcGradD1);
    fprintf(file,"solvEdpD1          : %lf\n",t.solvEdpD1);
  }
/*...................................................................*/

/*... solvT1*/
  if(fSolvT1){
    fprintf(file,"SolvT1             : %lf\n",t.solvT1);
    fprintf(file,"numeqT1            : %lf\n",t.numeqT1);
    fprintf(file,"dataStructT1       : %lf\n",t.dataStructT1);
    fprintf(file,"CellPloadT1        : %lf\n",t.CellPloadT1);
    fprintf(file,"CellTransientT1    : %lf\n",t.CellTransientT1);
    fprintf(file,"systFormT1         : %lf\n",t.systFormT1);
    fprintf(file,"rcGradT1           : %lf\n",t.rcGradT1);
    fprintf(file,"solvEdpT1          : %lf\n",t.solvEdpT1);
  }
/*...................................................................*/

/*... solvPres*/
  if(fSolvPres){
    fprintf(file,"SolvPres           : %lf\n",t.solvPres);
    fprintf(file,"numeqPres          : %lf\n",t.numeqPres);
    fprintf(file,"dataStructPres     : %lf\n",t.dataStructPres);
    fprintf(file,"systFormPres       : %lf\n",t.systFormPres);
    fprintf(file,"rcGradPres         : %lf\n",t.rcGradPres);
  }
/*...................................................................*/

/*... solvVel*/
  if(fSolvVel){
    fprintf(file,"SolvVel            : %lf\n",t.solvVel);
    fprintf(file,"numeqVel           : %lf\n",t.numeqVel);
    fprintf(file,"dataStructVel      : %lf\n",t.dataStructVel);
    fprintf(file,"CellPloadSimple    : %lf\n",t.cellPloadSimple);
    fprintf(file,"CellTransientSimple: %lf\n",t.cellTransientSimple);
    fprintf(file,"systFormVel        : %lf\n",t.systFormVel);
    fprintf(file,"rcGradVel          : %lf\n",t.rcGradVel);
    fprintf(file,"solvEdpFuild       : %lf\n",t.solvEdpFluid);
  }
/*...................................................................*/
  
  fprintf(file,"Total              : %lf\n",t.total);
/*...*/
  fprintf(file,"\nMalha:\n");
  fprintf(file,"nnode              : %d\n",mesh.nnode);
  fprintf(file,"nCell              : %d\n",mesh.numel);
  fprintf(file,"volume             : %lf\n",mesh.mQuality.volume);
  fprintf(file,"non-OrthMed        : %.1lf°\n",mesh.mQuality.nonOrthMed);
  fprintf(file,"non-OtthMax        : %.1lf°\n",mesh.mQuality.nonOrthMax);
  fprintf(file,"skewMed            : %lf\n",mesh.mQuality.skewMed);
  fprintf(file,"skewMax            : %lf\n",mesh.mQuality.skewMax);

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
  }
/*...................................................................*/

/*... OpenMp*/    
  fprintf(file, "\nOpenMp             :\n");
  fprintf(file,   "Solver             :\n");
  if(omp.fSolver)
    fprintf(file, "nThreads           : %d\n",omp.nThreadsSolver);
  else
    fprintf(file, "Disable\n");
  fprintf(file,   "Cell               :\n");
  if (omp.fCell)
    fprintf(file, "nThreads           : %d\n",omp.nThreadsCell);
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

#ifdef _MPICH_
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
