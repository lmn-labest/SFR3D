#include<WriteLog.h>

void writeLog(Mesh mesh    
             ,Solv solvD1 ,SistEq sistEqD1
             ,Time t
             ,char *nameIn,FILE *file){

  fprintf(file,"Log do execucao: %s\n\n",nameIn); 

  fprintf(file,"Time:\n");
  fprintf(file,"adjcency     : %lf\n",t.adjcency);
  fprintf(file,"geom         : %lf\n",t.geom);
  if( mesh.rcGrad == RCLSQUARE )
    fprintf(file,"lSquareMatrix: %lf\n",t.leastSquareMatrix);
  fprintf(file,"reord        : %lf\n",t.reord);
  fprintf(file,"matVecSparse : %lf\n",t.matVecSparse);
  fprintf(file,"dot          : %lf\n",t.dot); 
  fprintf(file,"Pcg          : %lf\n",t.pcg);
  fprintf(file,"precondDiag  : %lf\n",t.precondDiag);

/*... solvD1*/
  if(solvD1.flag){
    fprintf(file,"SolvD1       : %lf\n",t.solvD1);
    fprintf(file,"numeqD1      : %lf\n",t.numeqD1);
    fprintf(file,"dataStructD1 : %lf\n",t.dataStructD1);
    fprintf(file,"CellPloadD1  : %lf\n",t.CellPloadD1);
    fprintf(file,"systFormD1   : %lf\n",t.systFormD1);
    fprintf(file,"rcGradD1     : %lf\n",t.rcGradD1);
  }
/*...................................................................*/
  
  fprintf(file,"Total        : %lf\n",t.total);

/*...*/
  fprintf(file,"\nMalha:\n");
  fprintf(file,"nnode        : %d\n",mesh.nnode);
  fprintf(file,"nCell        : %d\n",mesh.numel);
/*...*/
  if(solvD1.flag){
    fprintf(file,"\nSitema D1:\n");
    fprintf(file,"nEq          : %d\n",sistEqD1.neq);
    fprintf(file,"nad          : %d\n",sistEqD1.nad);
    fprintf(file,"tol          : %e\n",solvD1.tol);
/*... tenica de armazenamento*/
    if(sistEqD1.storage == CSR)
      fprintf(file,"Armazenamento: CSR\n");
    else if(sistEqD1.storage == CSRD)
      fprintf(file,"Armazenamento: CSRD\n");
    else if(sistEqD1.storage == CSRC)
      fprintf(file,"Armazenamento: CSRC\n");
/*... tenica de armazenamento*/
    if(solvD1.solver == PCG)
      fprintf(file,"Iterativo    : PCG\n");
    else if(solvD1.solver == PBICGSTAB )
      fprintf(file,"Iterativo    : PBICGSTAB\n");
  }


} 
