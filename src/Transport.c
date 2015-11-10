#include<Transport.h>
/********************************************************************* 
 * TRANSPORT :Resolucao do problema de transporte no passo de tempo  * 
 * n+1                                                               * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m         -> nao definido                                         *
 * loadsTrans-> definicoes de cargas                                 *
 * mesh0     -> parametros                                           * 
 * mesh      -> parametros                                           * 
 * sistEqD   -> sistema de equacoes                                  * 
 * solvD     -> sistema de equacoes                                  * 
 * sc        -> tecnica de discretizacao temporal                    * 
 * pMesh     -> varaives de particionamento                          * 
 * opt       -> opcoes de arquivos                                   * 
 * preName   -> sufixo dos arquivos de saida                         * 
 * nameOut   -> string para o nome do arquivo de saida               * 
 * fileOut   -> nao definido                                         * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS: variasveis do problema de difusao atualizadas para o passo de* 
 * de tempo n+1                                                      * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void transport(Memoria *m   ,Loads *loadsTrans
              ,Mesh *mesh0  ,Mesh *mesh       ,SistEq *sistEqT
              ,Solv *solvT  ,Scheme sc        ,PartMesh *pMesh 
              ,FileOpt opt  ,char *preName    ,char *nameOut
              ,FILE *fileOut) 
{   
  
  short unsigned i;
  char str1[100],str2[100],str3[100],str4[100],str5[100],str6[100];
  DOUBLE rCell,rCell0,conv;

  zero(sistEqT->b0,sistEqT->neqNov,DOUBLEC);
/*... restricoes por centro de celula u0 e cargas por volume b0*/
  tm.CellPloadT1 = getTimeC() - tm.CellPloadT1;
  cellPload(loadsTrans            ,mesh->elm.geom.cc 
           ,mesh->elm.faceRt1    ,mesh->elm.faceLoadT1
           ,mesh->elm.geom.volume,sistEqT->id 
           ,mesh->elm.uT1        ,sistEqT->b0
           ,mesh->numelNov       ,mesh->ndfT[0]
           ,mesh->ndm            ,mesh->maxViz);
  tm.CellPloadT1 = getTimeC() - tm.CellPloadT1;
/*...................................................................*/

/*... discretizacao temporal*/
  if(sc.ddt.flag){
    tm.CellTransientT1 = getTimeC() - tm.CellTransientT1;
    cellTransient(mesh->elm.geom.volume   ,sistEqT->id     
                 ,mesh->elm.u0T1          ,mesh->elm.uT1
                 ,mesh->elm.densityUt1    ,sistEqT->b0
                 ,sc.ddt.dt     
                 ,mesh->numelNov          ,mesh->ndfT[0]
                 ,sc.ddt.type             ,true);
/*... u(n-1) = u(n)*/
    alphaProdVector(1.e0,mesh->elm.uT1
                   ,mesh->numel       ,mesh->elm.u0T1); 
/*...................................................................*/
    tm.CellTransientT1 = getTimeC() - tm.CellTransientT1;
  }
/*...................................................................*/
  
/*... correcao nao ortoganal*/ 
  for(i=0;i<sc.nlT1.maxIt;i++){

/*... calculo de: A(i),b(i)*/
    tm.systFormT1 = getTimeC() - tm.systFormT1;
    systFormTrans(loadsTrans            ,sc.advT1
               ,mesh->elm.node          ,mesh->elm.adj.nelcon  
               ,mesh->elm.nen           ,mesh->elm.adj.nViz   
               ,mesh->elm.geomType      ,mesh->elm.material.prop 
               ,mesh->elm.material.type ,mesh->elm.mat   
               ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
               ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea    
               ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
               ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc    
               ,mesh->elm.geom.vSkew    ,mesh->elm.geom.mvSkew   
               ,mesh->elm.geom.dcca     ,mesh->elm.densityUt1
               ,sistEqT->ia             ,sistEqT->ja      
               ,sistEqT->ad             ,sistEqT->al       
               ,sistEqT->b              ,sistEqT->id       
               ,mesh->elm.faceRt1       ,mesh->elm.faceLoadT1  
               ,mesh->elm.uT1           ,mesh->elm.gradUt1           
               ,mesh->elm.vel                                        
               ,mesh->elm.rCellUt1      ,sc.ddt
               ,sistEqT->neqNov         ,sistEqT->nad
               ,sistEqT->nadr      
               ,mesh->maxNo             ,mesh->maxViz
               ,mesh->ndm               ,mesh->numelNov
               ,mesh->ndfT[0]           ,sistEqT->storage
               ,true                    ,true   
               ,true                    ,sistEqT->unsym);   
    tm.systFormT1 = getTimeC() - tm.systFormT1;
/*...................................................................*/

/*... soma o vetor b(i) = b(i) + b0(i)*/
    addVector(1.0e0           ,sistEqT->b
             ,1.0e0           ,sistEqT->b0
             ,sistEqT->neqNov ,sistEqT->b);
/*...................................................................*/

/*... soma o vetor R(i) = R(i) + b0(i)*/ 
    updateCellValue(mesh->elm.rCellUt1 ,sistEqT->b0
                   ,sistEqT->id        ,&sistEqT->iNeq
                   ,mesh->numelNov     ,mesh->ndfT[0]
                   ,true               ,false);
/*...................................................................*/

/*...*/ 
    if( i == 0 ){
      rCell  = rCell0 = sqrt(dot(mesh->elm.rCellUt1
                            ,mesh->elm.rCellUt1 
                            ,mesh->numelNov));
      conv   = rCell0*sc.nlT1.tol;
    }
    else
      rCell  = sqrt(dot(mesh->elm.rCellUt1 
                   ,mesh->elm.rCellUt1 
                   ,mesh->numelNov));
        
    if(!mpiVar.myId ){
      printf("it: %8d %.6e\n",i,rCell/rCell0);  
      if(opt.fItPlot)  
        fprintf(opt.fileItPlot[FITPLOTT1]
               ,"%9d %.6e %0.6e\n",i,rCell/rCell0,rCell);
    }
    if(rCell < conv) break;
/*...................................................................*/

/*...*/
    tm.solvT1 = getTimeC() - tm.solvT1;
    solverC(m               
           ,sistEqT->neq       ,sistEqT->neqNov  
           ,sistEqT->nad       ,sistEqT->nadr
           ,sistEqT->ia        ,sistEqT->ja  
           ,sistEqT->al        ,sistEqT->ad,sistEqT->au
           ,sistEqT->b         ,sistEqT->x
           ,&sistEqT->iNeq
           ,solvT->tol         ,solvT->maxIt     
           ,sistEqT->storage   ,solvT->solver
           ,solvT->fileSolv    ,solvT->log  
           ,false              ,false
           ,sistEqT->unsym     ,false);  
    tm.solvT1 = getTimeC() - tm.solvT1;
/*...................................................................*/

/*... x -> uT1*/
    updateCellValue(mesh->elm.uT1 ,sistEqT->x
                   ,sistEqT->id   ,&sistEqT->iNeq
                   ,mesh->numel   ,mesh->ndfT[0]
                   ,false         ,true);
/*...................................................................*/

/*... reconstruindo do gradiente*/
    tm.rcGradT1 = getTimeC() - tm.rcGradT1;
    rcGradU(m                       ,loadsTrans
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
           ,mesh->elm.faceRt1       ,mesh->elm.faceLoadT1    
           ,mesh->elm.uT1           ,mesh->elm.gradUt1                 
           ,mesh->node.uT1          ,sc.rcGrad
           ,mesh->maxNo             ,mesh->maxViz
           ,mesh->ndfT[0]           ,mesh->ndm
           ,&pMesh->iNo             ,&pMesh->iEl  
           ,mesh->numelNov          ,mesh->numel        
           ,mesh->nnodeNov          ,mesh->nnode); 
    tm.rcGradT1 = getTimeC() - tm.rcGradT1;
/*...................................................................*/

    if(opt.fItPlotRes){  

/*... interpolacao das variaveis da celulas para pos nos (Grad)*/
      interCellNode(m                  ,loadsTrans
                   ,mesh->node.gradUt1 ,mesh->elm.gradUt1 
                   ,mesh->elm.node     ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc  ,mesh->node.x   
                   ,mesh->elm.geom.xm               
                   ,mesh->elm.nen      ,mesh->elm.adj.nViz
                   ,mesh->elm.faceRt1  ,mesh->elm.faceLoadT1    
                   ,&pMesh->iNo          
                   ,mesh->numelNov     ,mesh->numel 
                   ,mesh->nnodeNov     ,mesh->nnode 
                   ,mesh->maxNo        ,mesh->maxViz   
                   ,mesh->ndm          ,1 
                   ,mesh->ndm
                   ,false              ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (vel)*/
      interCellNode(m                  ,loadsTrans
                   ,mesh->node.vel     ,mesh->elm.vel        
                   ,mesh->elm.node     ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc  ,mesh->node.x  
                   ,mesh->elm.geom.xm
                   ,mesh->elm.nen      ,mesh->elm.adj.nViz
                   ,mesh->elm.faceRt1  ,mesh->elm.faceLoadT1  
                   ,&pMesh->iNo          
                   ,mesh->numelNov     ,mesh->numel        
                   ,mesh->nnodeNov     ,mesh->nnode 
                   ,mesh->maxNo        ,mesh->maxViz   
                   ,mesh->ndm          ,1
                   ,mesh->ndm      
                   ,false              ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (uT1)*/
      interCellNode(m                ,loadsTrans
                    ,mesh->node.uT1   ,mesh->elm.uT1 
                    ,mesh->elm.node   ,mesh->elm.geomType
                    ,mesh->elm.geom.cc,mesh->node.x    
                    ,mesh->elm.geom.xm              
                    ,mesh->elm.nen    ,mesh->elm.adj.nViz
                    ,mesh->elm.faceRt1,mesh->elm.faceLoadT1    
                    ,&pMesh->iNo             
                    ,mesh->numelNov   ,mesh->numel
                    ,mesh->nnodeNov   ,mesh->nnode 
                    ,mesh->maxNo      ,mesh->maxViz 
                    ,mesh->ndfT[0]    ,1
                    ,mesh->ndm
                    ,true             ,2);
/*...................................................................*/

/*... globalizacao das variaveis*/
/*... uT1(Node)*/
       dGlobalNode(m                  ,pMesh
                  ,mesh0->node.uT1    ,mesh->node.uT1     
                  ,mesh->ndfT[0]      ,1               );
/*... gradUt1(Node)*/
       dGlobalNode(m                  ,pMesh
                  ,mesh0->node.gradUt1,mesh->node.gradUt1     
                  ,mesh->ndm          ,1               );
/*... vel(Node)*/
      dGlobalNode(m                  ,pMesh
                 ,mesh0->node.vel    ,mesh->node.vel         
                 ,mesh->ndm          ,1               );
/*... uT1(Cel)*/
       dGlobalCel(m                   ,pMesh
                 ,mesh0->elm.uT1      ,mesh->elm.uT1
                 ,mesh->numelNov 
                 ,mesh->ndfT[0]      ,1);
/*... gradUt1(Cel)*/
       dGlobalCel(m                   ,pMesh
                 ,mesh0->elm.gradUt1  ,mesh->elm.gradUt1
                 ,mesh->numelNov 
                 ,mesh->ndm           ,1);
/*... vel(Cel)*/
      dGlobalCel(m                   ,pMesh
                ,mesh0->elm.vel      ,mesh->elm.vel       
                ,mesh->numelNov 
                ,mesh->ndm           ,1);
/*...................................................................*/

/*...*/
       if(!mpiVar.myId){
         fName(preName,sc.ddt.timeStep,i,19,&nameOut);
         strcpy(str1,"elT1");
         strcpy(str2,"noT1");
         strcpy(str3,"elGradT1");
         strcpy(str4,"noGradT1");
         strcpy(str5,"elVel");
         strcpy(str6,"noVel");
/*...*/
         wResVtkTrans(m                 ,mesh0->node.x      
                    ,mesh0->elm.node    ,mesh0->elm.mat    
                    ,mesh0->elm.nen     ,mesh0->elm.geomType
                    ,mesh0->elm.uT1     ,mesh0->node.uT1 
                    ,mesh0->elm.gradUt1 ,mesh0->node.gradUt1  
                    ,mesh0->elm.vel     ,mesh0->node.vel      
                    ,mesh0->nnode       ,mesh0->numel  
                    ,mesh0->ndm         ,mesh0->maxNo 
                    ,mesh0->numat       ,mesh0->ndfT[0]
                    ,str1               ,str2         
                    ,str3               ,str4         
                    ,str5               ,str6         
                    ,nameOut            ,opt.bVtk    
                    ,sc.ddt             ,fileOut);
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  } 
/*...................................................................*/
}
/*********************************************************************/
