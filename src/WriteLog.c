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
void writeLog(Mesh mesh    ,Scheme sc
             ,Solv *solvD1 ,SistEq *sistEqD1
             ,Time t
             ,bool const fSolvD1
             ,char *nameIn,FILE *file){

  
  fprintf(file,"Log do execucao : %s\n\n",nameIn); 
  
  if(mpiVar.nPrcs > 1){
    fprintf(file,"myId            :%2d\n",mpiVar.myId); 
    fprintf(file,"nPrcs           :%2d\n\n",mpiVar.nPrcs); 
  } 
  
  fprintf(file,"Time:\n");
  fprintf(file,"adjcency        : %lf\n",t.adjcency);
  fprintf(file,"geom            : %lf\n",t.geom);
  if( sc.rcGrad == RCLSQUARE )
    fprintf(file,"lSquareMatrix   : %lf\n",t.leastSquareMatrix);
  fprintf(file,"reord           : %lf\n",t.reord);

/*... Blas*/
  fprintf(file,"matVecSparse    : %lf\n",t.matVecSparse);
  fprintf(file,"dot             : %lf\n",t.dot); 
/*... Blas overHead do mpi */
  if(mpiVar.nPrcs > 1){
    fprintf(file,"matVecOverHead  : %lf\n",t.matVecOverHeadMpi);
    fprintf(file,"dotOverHead     : %lf\n",t.dotOverHeadMpi);
  }

/*... Solver*/
  fprintf(file,"Pcg             : %lf\n",t.pcg);
  fprintf(file,"precondDiag     : %lf\n",t.precondDiag);

/*... particionamento*/
  if(mpiVar.nPrcs > 1){
    fprintf(file,"partdMesh       : %lf\n",t.partdMesh);
    fprintf(file,"partdMeshComm   : %lf\n",t.partdMeshCom);
  }
/*... comunicacao entre as particoes*/
  if(mpiVar.nPrcs > 1){
    fprintf(file,"overHeadCelMpi  : %lf\n",t.overHeadCelMpi);
    fprintf(file,"overHeadNodMpi  : %lf\n",t.overHeadNodMpi);
    fprintf(file,"overHeadNeqMpi  : %lf\n",t.overHeadNeqMpi);
    fprintf(file,"overHeadGCelMpi : %lf\n",t.overHeadGCelMpi);
    fprintf(file,"overHeadGNodMpi : %lf\n",t.overHeadGNodMpi);
    fprintf(file,"overHeadTotalMpi: %lf\n",t.overHeadTotalMpi);
  }

/*... solvD1*/
  if(fSolvD1){
    fprintf(file,"SolvD1          : %lf\n",t.solvD1);
    fprintf(file,"numeqD1         : %lf\n",t.numeqD1);
    fprintf(file,"dataStructD1    : %lf\n",t.dataStructD1);
    fprintf(file,"CellPloadD1     : %lf\n",t.CellPloadD1);
    fprintf(file,"CellTransientD1 : %lf\n",t.CellTransientD1);
    fprintf(file,"systFormD1      : %lf\n",t.systFormD1);
    fprintf(file,"rcGradD1        : %lf\n",t.rcGradD1);
    fprintf(file,"solvEdoD1       : %lf\n",t.solvEdoD1);
  }
/*...................................................................*/
  
  fprintf(file,"Total           : %lf\n",t.total);
/*...*/
  fprintf(file,"\nMalha:\n");
  fprintf(file,"nnode           : %d\n",mesh.nnode);
  fprintf(file,"nCell           : %d\n",mesh.numel);
  fprintf(file,"volume          : %lf\n",mesh.mQuality.volume);
  fprintf(file,"non-OrthMed     : %.2lf°\n",mesh.mQuality.nonOrthMed);
  fprintf(file,"non-OtthMax     : %.2lf°\n",mesh.mQuality.nonOrthMax);
  fprintf(file,"skewMed         : %lf\n",mesh.mQuality.skewMed);
  fprintf(file,"skewMax         : %lf\n",mesh.mQuality.skewMax);
  
/*...*/
  if(fSolvD1){
    fprintf(file,"\nSistema D1:\n");
/*... tecnica de armazenamento*/
    if(sistEqD1->storage == CSR)
      fprintf(file,"Armazenamento   : CSR\n");
    else if(sistEqD1->storage == CSRD)
      fprintf(file,"Armazenamento   : CSRD\n");
    else if(sistEqD1->storage == CSRC)
      fprintf(file,"Armazenamento   : CSRC\n");
    else if(sistEqD1->storage == ELLPACK)
      fprintf(file,"Armazenamento   : ELLPACK\n");
    fprintf(file,"nEq             : %d\n",sistEqD1->neq);
    if(mpiVar.nPrcs > 1)
      fprintf(file,"nEqNov          : %d\n",sistEqD1->neqNov);
/*... ELLPACK*/
    if(sistEqD1->storage == ELLPACK){
      fprintf(file,"nAd             : %d\n",sistEqD1->nad);
      fprintf(file,"nAd Overhead    : %d\n",sistEqD1->neq*mesh.maxViz);
    }
/*... CSRD*/
    else{
      fprintf(file,"band Min        : %d\n",sistEqD1->bandCsr[BANDCSRMIN]);
      fprintf(file,"band Max        : %d\n",sistEqD1->bandCsr[BANDCSRMAX]);
      fprintf(file,"band Med        : %d\n",sistEqD1->bandCsr[BANDCSRMED]);
      fprintf(file,"nAd             : %d\n",sistEqD1->nad);
      if(mpiVar.nPrcs > 1)
        fprintf(file,"nAdR            : %d\n",sistEqD1->nadr);
    }
/*...................................................................*/
    fprintf(file,"tol             : %e\n",solvD1->tol);
/*... tenica de armazenamento*/
    if(solvD1->solver == PCG)
      fprintf(file,"Iterativo       : PCG\n");
    else if(solvD1->solver == PBICGSTAB )
      fprintf(file,"Iterativo       : PBICGSTAB\n");
  }

}
/*********************************************************************/ 
