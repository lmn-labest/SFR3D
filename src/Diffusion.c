#include<Diffusion.h>
/********************************************************************* 
 * DIFFUSION :Resolucao do problema de difusao pura no passo de tempo* 
 * n+1                                                               * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m        -> nao definido                                          *
 * loadsDif -> definicoes de cargas                                  *
 * mesh0    -> parametros                                            * 
 * mesh     -> parametros                                            * 
 * sistEqD  -> sistema de equacoes                                   * 
 * solvD    -> sistema de equacoes                                   * 
 * sc       -> tecnica de discretizacao temporal                     * 
 * pMesh    -> varaives de particionamento                           * 
 * opt      -> opcoes de arquivos                                    * 
 * preName  -> sufixo dos arquivos de saida                          * 
 * nameOut  -> string para o nome do arquivo de saida                * 
 * fileOut  -> nao definido                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS: variasveis do problema de difusao atualizadas para o passo de* 
 * de tempo n+1                                                      * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void diffusion(Memoria *m   ,Loads *loadsDif
              ,Mesh *mesh0  ,Mesh *mesh     ,SistEq *sistEqD
              ,Solv *solvD  ,Scheme sc      ,PartMesh *pMesh 
              ,FileOpt opt  ,char *preName  ,char *nameOut
              ,FILE *fileOut) 
{   
  
  short unsigned i;
  char str1[100],str2[100],str3[100],str4[100];
  DOUBLE rCell,rCell0,conv;

  zero(sistEqD->b0,sistEqD->neqNov,DOUBLEC);
/*... restricoes por centro de celula u0 e cargas por volume b0*/
  tm.CellPloadD1 = getTimeC() - tm.CellPloadD1;
  cellPload(loadsDif              ,mesh->elm.geom.cc 
           ,mesh->elm.faceRd1    ,mesh->elm.faceLoadD1
           ,mesh->elm.geom.volume,sistEqD->id 
           ,mesh->elm.uD1        ,sistEqD->b0
           ,mesh->numelNov       ,mesh->ndfD[0]
           ,mesh->ndm            ,mesh->maxViz);
  tm.CellPloadD1 = getTimeC() - tm.CellPloadD1;
/*...................................................................*/

/*... discretizacao temporal*/
  if(sc.ddt.flag){
    tm.CellTransientD1 = getTimeC() - tm.CellTransientD1;
    cellTransient(mesh->elm.geom.volume   ,sistEqD->id     
                 ,mesh->elm.u0D1          ,mesh->elm.uD1
                 ,mesh->elm.densityUd1    ,sistEqD->b0
                 ,sc.ddt                  ,mesh->numelNov 
                 ,mesh->ndfD[0]           ,true);
/*... u(n-1) = u(n)*/
    alphaProdVector(1.e0,mesh->elm.uD1
                   ,mesh->numel       ,mesh->elm.u0D1); 
/*...................................................................*/
    tm.CellTransientD1 = getTimeC() - tm.CellTransientD1;
  }
/*...................................................................*/

/*... correcao nao ortoganal*/ 
  for(i=0;i<sc.nlD1.maxIt;i++){

/*... calculo de: A(i),b(i)*/
    tm.systFormD1 = getTimeC() - tm.systFormD1;
    systFormDif(loadsDif
               ,mesh->elm.node          ,mesh->elm.adj.nelcon  
               ,mesh->elm.nen           ,mesh->elm.adj.nViz   
               ,mesh->elm.geomType      ,mesh->elm.material.prop 
               ,mesh->elm.material.type ,mesh->elm.mat   
               ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
               ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea    
               ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
               ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc    
               ,mesh->elm.geom.vSkew    ,mesh->elm.geom.mvSkew   
               ,mesh->elm.geom.dcca     ,mesh->elm.densityUd1
               ,sistEqD->ia            ,sistEqD->ja      
               ,sistEqD->al            ,sistEqD->ad       
               ,sistEqD->b             ,sistEqD->id       
               ,mesh->elm.faceRd1       ,mesh->elm.faceLoadD1  
               ,mesh->elm.uD1           ,mesh->elm.gradUd1           
               ,mesh->elm.rCellUd1      ,sc.ddt
               ,sistEqD->neq            ,sistEqD->neqNov        
               ,sistEqD->nad            ,sistEqD->nadr      
               ,mesh->maxNo             ,mesh->maxViz
               ,mesh->ndm               ,mesh->numelNov
               ,mesh->ndfD[0]           ,sistEqD->storage
               ,true                    ,true   
               ,true                    ,sistEqD->unsym);   
    tm.systFormD1 = getTimeC() - tm.systFormD1;
/*...................................................................*/

/*... soma o vetor b(i) = b(i) + b0(i)*/
    addVector(1.0e0           ,sistEqD->b
             ,1.0e0           ,sistEqD->b0
             ,sistEqD->neqNov,sistEqD->b);
/*...................................................................*/

/*... soma o vetor R(i) = R(i) + b0(i)*/
    updateCellValue(mesh->elm.rCellUd1 ,sistEqD->b0
                   ,sistEqD->id       ,&sistEqD->iNeq
                   ,mesh->numelNov     ,mesh->ndfD[0]
                   ,true               ,false);
/*...................................................................*/

/*...*/ 
    if( i == 0 ){
      rCell  = rCell0 = sqrt(dot(mesh->elm.rCellUd1
                            ,mesh->elm.rCellUd1 
                            ,mesh->numelNov));
      conv   = rCell0*sc.nlD1.tol;
    }
    else
      rCell  = sqrt(dot(mesh->elm.rCellUd1 
                   ,mesh->elm.rCellUd1 
                   ,mesh->numelNov));
    if(!mpiVar.myId ){
      printf("it: %8d %.6e\n",i,rCell/rCell0);  
      if(opt.fItPlot)  
        fprintf(opt.fileItPlot[FITPLOTD1]
               ,"%9d %.6e %0.6e\n",i,rCell/rCell0,rCell);
    }
    if(rCell < conv) break;
/*...................................................................*/

/*...*/
    tm.solvD1 = getTimeC() - tm.solvD1;
    solverC(m               
           ,sistEqD->neq       ,sistEqD->neqNov  
           ,sistEqD->nad       ,sistEqD->nadr
           ,sistEqD->ia        ,sistEqD->ja  
           ,sistEqD->al        ,sistEqD->ad,sistEqD->au
           ,sistEqD->b         ,sistEqD->x
           ,&sistEqD->iNeq
           ,solvD->tol         ,solvD->maxIt     
           ,sistEqD->storage   ,solvD->solver
           ,solvD->fileSolv    ,solvD->log  
           ,false               ,false
           ,sistEqD->unsym     ,false);  
    tm.solvD1 = getTimeC() - tm.solvD1;
/*...................................................................*/


/*... x -> uD1*/
    updateCellValue(mesh->elm.uD1 ,sistEqD->x
                   ,sistEqD->id  ,&sistEqD->iNeq
                   ,mesh->numel   ,mesh->ndfD[0]
                   ,false         ,true);
/*...................................................................*/

/*... reconstruindo do gradiente*/
    tm.rcGradD1 = getTimeC() - tm.rcGradD1;
    rcGradU(m                       ,loadsDif
           ,mesh->elm.node          ,mesh->elm.adj.nelcon
           ,mesh->elm.geom.cc       ,mesh->node.x   
           ,mesh->elm.nen           ,mesh->elm.adj.nViz 
           ,mesh->elm.geomType      ,mesh->elm.material.prop 
           ,mesh->elm.mat 
           ,mesh->elm.leastSquare   ,mesh->elm.leastSquareR
           ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
           ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea    
           ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
           ,mesh->elm.geom.vSkew      
           ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc    
           ,mesh->elm.geom.dcca
           ,mesh->elm.faceRd1       ,mesh->elm.faceLoadD1    
           ,mesh->elm.uD1           ,mesh->elm.gradUd1                 
           ,mesh->node.uD1          ,sc.rcGrad
           ,mesh->maxNo             ,mesh->maxViz
           ,mesh->ndfD[0]           ,mesh->ndm
           ,&pMesh->iNo             ,&pMesh->iEl  
           ,mesh->numelNov          ,mesh->numel        
           ,mesh->nnodeNov          ,mesh->nnode); 
    tm.rcGradD1 = getTimeC() - tm.rcGradD1;
/*...................................................................*/

    if(opt.fItPlotRes){  

/*... interpolacao das variaveis da celulas para pos nos (Grad)*/
      interCellNode(m                  ,loadsDif
                   ,mesh->node.gradUd1 ,mesh->elm.gradUd1 
                   ,mesh->elm.node     ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc  ,mesh->node.x   
                   ,mesh->elm.geom.xm               
                   ,mesh->elm.nen      ,mesh->elm.adj.nViz
                   ,mesh->elm.faceRd1  ,mesh->elm.faceLoadD1    
                   ,&pMesh->iNo          
                   ,mesh->numelNov     ,mesh->numel 
                   ,mesh->nnodeNov     ,mesh->nnode 
                   ,mesh->maxNo        ,mesh->maxViz   
                   ,mesh->ndm          ,1 
                   ,mesh->ndm
                   ,false              ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (uD1)*/
       interCellNode(m                ,loadsDif
                    ,mesh->node.uD1   ,mesh->elm.uD1 
                    ,mesh->elm.node   ,mesh->elm.geomType
                    ,mesh->elm.geom.cc,mesh->node.x    
                    ,mesh->elm.geom.xm              
                    ,mesh->elm.nen    ,mesh->elm.adj.nViz
                    ,mesh->elm.faceRd1,mesh->elm.faceLoadD1    
                    ,&pMesh->iNo             
                    ,mesh->numelNov   ,mesh->numel
                    ,mesh->nnodeNov   ,mesh->nnode 
                    ,mesh->maxNo      ,mesh->maxViz 
                    ,mesh->ndfD[0]    ,1
                    ,mesh->ndm
                    ,true             ,2);
/*...................................................................*/

/*... globalizacao das variaveis*/
/*... uD1(Node)*/
       dGlobalNode(m                  ,pMesh
                  ,mesh0->node.uD1    ,mesh->node.uD1     
                  ,mesh->ndfD[0]      ,1               );
          
/*... gradUd1(Node)*/
       dGlobalNode(m                  ,pMesh
                  ,mesh0->node.gradUd1,mesh->node.gradUd1     
                  ,mesh->ndm          ,1               );
/*... uD1(Cel)*/
       dGlobalCel(m                   ,pMesh
                 ,mesh0->elm.uD1      ,mesh->elm.uD1
                 ,mesh->numelNov 
                 ,mesh->ndfD[0]      ,1);
/*... gradUd1(Cel)*/
       dGlobalCel(m                   ,pMesh
                 ,mesh0->elm.gradUd1  ,mesh->elm.gradUd1
                 ,mesh->numelNov 
                 ,mesh->ndm           ,1);
/*...................................................................*/

/*...*/
       if(!mpiVar.myId){
         fName(preName,sc.ddt.timeStep,i,15,&nameOut);
         strcpy(str1,"elD1");
         strcpy(str2,"noD1");
         strcpy(str3,"elGradD1");
         strcpy(str4,"noGradD1");
/*...*/
         wResVtkDif(m                  ,mesh0->node.x      
                   ,mesh0->elm.node     ,mesh0->elm.mat    
                   ,mesh0->elm.nen      ,mesh0->elm.geomType
                   ,mesh0->elm.uD1      ,mesh0->node.uD1 
                   ,mesh0->elm.gradUd1  ,mesh0->node.gradUd1 
                   ,mesh0->nnode        ,mesh0->numel  
                   ,mesh0->ndm          ,mesh0->maxNo 
                   ,mesh0->numat        ,mesh0->ndfD[0]
                   ,str1                ,str2         
                   ,str3                ,str4         
                   ,nameOut             ,opt.bVtk    
                   ,sc.ddt              ,fileOut);
      }
/*...................................................................*/
    }
/*...................................................................*/
  } 
/*...................................................................*/
}
/*********************************************************************/
