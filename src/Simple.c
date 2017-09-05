#include<Simple.h>
/********************************************************************* 
 * Data de criacao    : 17/07/2016                                   *
 * Data de modificaco : 20/08/2017                                   * 
 *-------------------------------------------------------------------* 
 * SIMPLESOLVER3D: metodo simple e simpleC para escoamentos 3D       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void simpleSolver3D(Memoria *m        
                   ,Loads *loadsVel   ,Loads *loadsPres 
                   ,Mesh *mesh0       ,Mesh *mesh       
                   ,SistEq *sistEqVel ,SistEq *sistEqPres
                   ,Solv *solvVel     ,Solv *solvPres 
                   ,Simple *sp
                   ,Scheme sc         ,PartMesh *pMesh 
                   ,FileOpt opt       ,char *preName 
                   ,char *nameOut     ,FILE *fileOut){
  FILE *fStop=NULL;
	short unsigned ndfVel = mesh->ndfF-1;
  short unsigned conv;
  int itSimple;
  int nonOrth;
  short unsigned kZeroVel  = sp->kZeroVel;
  short unsigned kZeroPres = sp->kZeroPres;
  INT jj = 1;
  DOUBLE time,timei;
  DOUBLE *b1,*b2,*b3,*bPc,*xu1,*xu2,*xu3,*xp,*adU1,*adU2,*adU3;
  DOUBLE *rCellPc;
/*...*/
  DOUBLE rU[3],rU0[3],tmp,tb[3],rMass0,rMass;
/*...*/
  DOUBLE tolSimpleU1,tolSimpleU2,tolSimpleU3,tolSimpleMass;
/*...*/
  bool xMomentum ,yMomentum ,zMomentum ,pCor;
  bool relRes  = false;
  bool fPrint  = false;
  DOUBLE cfl,reynolds;
  bool fParameter[2];

  time = getTimeC();
/*...*/
  b1 = b2 = b3 = bPc = NULL;
  xu1= xu2 = xu3 = NULL;
  adU1 = adU2 = adU3 = NULL;

  b1       =  sistEqVel->b; 
  b2       = &sistEqVel->b[sistEqVel->neq]; 
  if(ndfVel == 3) b3 = &sistEqVel->b[2*sistEqVel->neq];
  bPc      = sistEqPres->b; 

  xu1      = sistEqVel->x;
  xu2      = &sistEqVel->x[sistEqVel->neq];
  if (ndfVel == 3) xu3 = &sistEqVel->x[2*sistEqVel->neq];
  xp       = sistEqPres->x;

  adU1 = sistEqVel->ad;
  adU2 = &sistEqVel->ad[sistEqVel->neq];
  if (ndfVel == 3) adU3 = &sistEqVel->ad[2*sistEqVel->neq];

  rCellPc  = mesh->elm.rCellPres;
/*...................................................................*/

/*...*/
  tolSimpleU1   = tolSimpleU2 = tolSimpleU3 = sp->tolVel;
  tolSimpleMass = sp->tolPres;
/*...................................................................*/

/*...*/
  rMass0 = 1.e0;
  rMass  = 0.e0;
  rU[0]  = rU[1]  = rU[2]  = 0.e0;
  rU0[0] = rU0[1] = rU0[2] = 1.e0;
  conv = 0;
  xMomentum = yMomentum = zMomentum = true;
/*...................................................................*/

/*...*/
  zero(sistEqVel->b0 ,sistEqVel->neqNov*ndfVel,DOUBLEC);
  zero(sistEqPres->b0,sistEqPres->neqNov      ,DOUBLEC);
/*...................................................................*/

/*... restricoes por centro de celula u0 e cargas por volume b0*/
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
  cellPloadSimple(loadsPres            ,mesh->elm.geom.cc 
                 ,mesh->elm.faceRpres  ,mesh->elm.faceLoadPres  
                 ,mesh->elm.geom.volume
                 ,sistEqVel->id        ,sistEqPres->id
                 ,mesh->elm.vel        ,mesh->elm.pressure 
                 ,sistEqVel->b0        ,sistEqPres->b0
                 ,mesh->numelNov       ,ndfVel
                 ,mesh->ndm            ,mesh->maxViz);
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
/*...................................................................*/

/*... discretizacao temporal*/
  if(sc.ddt.flag){
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
    cellTransientSimple(mesh->elm.geom.volume   ,sistEqVel->id     
                       ,mesh->elm.vel0          ,mesh->elm.vel 
                       ,mesh->elm.densityFluid  ,sistEqVel->b0
                       ,sc.ddt                  ,sistEqVel->neqNov
                       ,mesh->numelNov          ,ndfVel        
                       ,true);
/*... vel(n-1) = vel(n)*/
    alphaProdVector(1.e0,mesh->elm.vel
                   ,mesh->numel*ndfVel ,mesh->elm.vel0); 
/*...................................................................*/
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
  }
/*...................................................................*/

/*...*/
  for(itSimple=0;itSimple<sp->maxIt;itSimple++){
/*...*/
     if ( (fStop=fopen("stopSimple.mvf","r")) !=NULL){
		   fclose(fStop);			
       break;
     } 
/*...................................................................*/

/*... calculo da matrix jacobiana das velocidades
                            | du1dx1 du1dx2 du1dx3 |   
                            | du2dx1 du2dx2 du2dx3 |   
                            | du3dx1 du3dx2 du3dx3 |
*/   
     tm.rcGradVel  = getTimeC() - tm.rcGradVel;
     rcGradU(m                      ,loadsVel
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
           ,mesh->elm.faceRvel      ,mesh->elm.faceLoadVel   
           ,mesh->elm.vel           ,mesh->elm.gradVel
           ,mesh->node.vel          ,sc.rcGrad
           ,mesh->maxNo             ,mesh->maxViz
           ,ndfVel                  ,mesh->ndm
           ,&pMesh->iNo             ,&pMesh->iEl  
           ,mesh->numelNov          ,mesh->numel        
           ,mesh->nnodeNov          ,mesh->nnode); 
     tm.rcGradVel = getTimeC() - tm.rcGradVel;
/*...................................................................*/

/*... reconstruindo do gradiente da pressao*/
     tm.rcGradPres = getTimeC() - tm.rcGradPres;
     rcGradU(m                     ,loadsPres
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
           ,mesh->elm.faceRpres     ,mesh->elm.faceLoadPres  
           ,mesh->elm.pressure      ,mesh->elm.gradPres
           ,mesh->node.pressure     ,sc.rcGrad
           ,mesh->maxNo             ,mesh->maxViz
           ,1                       ,mesh->ndm
           ,&pMesh->iNo             ,&pMesh->iEl  
           ,mesh->numelNov          ,mesh->numel        
           ,mesh->nnodeNov          ,mesh->nnode); 
     tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*... montagem do sistema u, v e w*/
     tm.systFormVel = getTimeC() - tm.systFormVel;
     systFormSimpleVel(loadsVel              ,loadsPres
                    ,sc.advVel               ,sc.diffVel                
                    ,sp->type
                    ,mesh->elm.node          ,mesh->elm.adj.nelcon  
                    ,mesh->elm.nen           ,mesh->elm.adj.nViz   
                    ,mesh->elm.geomType      ,mesh->elm.material.prop 
                    ,mesh->elm.material.type ,mesh->elm.mat  
                    ,mesh->elm.geom.cc 
                    ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
                    ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea    
                    ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
                    ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc    
                    ,mesh->elm.geom.vSkew    ,mesh->elm.geom.mvSkew   
                    ,mesh->elm.geom.dcca     ,mesh->elm.densityFluid
                    ,sistEqVel->ia           ,sistEqVel->ja      
                    ,sistEqVel->al           ,sistEqVel->ad       
                    ,sistEqVel->b            ,sistEqVel->id       
                    ,mesh->elm.faceRvel      ,mesh->elm.faceLoadVel 
                    ,mesh->elm.faceRpres     ,mesh->elm.faceLoadPres
                    ,mesh->elm.pressure      ,mesh->elm.gradPres           
                    ,mesh->elm.vel           ,mesh->elm.gradVel 
                    ,sp->d                   ,sp->alphaVel          
                    ,mesh->elm.rCellVel      ,sc.ddt
                    ,sistEqVel->neq          ,sistEqVel->neqNov      
                    ,sistEqVel->nad          ,sistEqVel->nadr      
                    ,mesh->maxNo             ,mesh->maxViz
                    ,mesh->ndm               ,mesh->numelNov
                    ,ndfVel                  ,sistEqVel->storage
                    ,true                    ,true   
                    ,true                    ,sistEqVel->unsym
                    ,sp->sPressure);   
     tm.systFormVel = getTimeC() - tm.systFormVel;
/*...................................................................*/

/*... soma o vetor b(i) = b(i) + b0(i)*/
     addVector(1.0e0                    ,sistEqVel->b
              ,1.0e0                    ,sistEqVel->b0
              ,sistEqVel->neqNov*ndfVel ,sistEqVel->b);
/*...................................................................*/

/*... soma o vetor R(i) = R(i) + b0(i)*/ 
     updateCellValueSimple(mesh->elm.rCellVel ,sistEqVel->b0
                   ,sistEqVel->id      ,&sistEqVel->iNeq
                   ,mesh->numelNov     ,sistEqVel->neqNov
                   ,ndfVel
                   ,true               ,false);
/*...................................................................*/

/*...*/
     tb[0] = sqrt(dot(b1,b1,sistEqVel->neqNov));
     tb[1] = sqrt(dot(b2,b2,sistEqVel->neqNov));
     if (ndfVel == 3) tb[2] = sqrt(dot(b3,b3,sistEqVel->neqNov));
     if(itSimple == 0){ 
       tmp   = max(tb[0],tb[1]);      
       if (ndfVel == 3) tmp   = max(tmp,tb[2]);
     }       
/*...*/ 
     xMomentum = true;
     if( tb[0] < tmp*SZERO ) xMomentum = false;
/*...................................................................*/

/*...*/
     yMomentum = true;
     if( tb[1] < tmp*SZERO ) yMomentum = false;
/*...................................................................*/

/*...*/
     if (ndfVel == 3) {
       zMomentum = true;
       if (tb[2] < tmp*SZERO) zMomentum = false;
     }
/*...................................................................*/

/*... solver Au = bu (velocidade estimadas)*/
     if(xMomentum){
       if(fPrint) printf("Quantidade de movimento u1:\n");
       tm.solvVel = getTimeC() - tm.solvVel;
       solverC(m               
           ,sistEqVel->neq     ,sistEqVel->neqNov  
           ,sistEqVel->nad     ,sistEqVel->nadr
           ,sistEqVel->ia      ,sistEqVel->ja  
           ,sistEqVel->al      ,adU1          ,sistEqVel->au
           ,b1                 ,xu1          
           ,&sistEqVel->iNeq   ,&sistEqVel->omp
           ,solvVel->tol       ,solvVel->maxIt     
           ,sistEqVel->storage ,solvVel->solver
           ,solvVel->fileSolv  ,solvVel->log  
           ,true               ,sistEqVel->unsym   );   
       tm.solvVel = getTimeC() - tm.solvVel;
     }
/*...................................................................*/

/*... solver Av = bv (velocidade estimadas)*/
     if(yMomentum){
       if(fPrint) printf("Quantidade de movimento u2:\n");
       tm.solvVel = getTimeC() - tm.solvVel;
       solverC(m               
           ,sistEqVel->neq     ,sistEqVel->neqNov  
           ,sistEqVel->nad     ,sistEqVel->nadr
           ,sistEqVel->ia      ,sistEqVel->ja  
           ,sistEqVel->al      ,adU2         ,sistEqVel->au
           ,b2                 ,xu2
           ,&sistEqVel->iNeq   ,&sistEqVel->omp
           ,solvVel->tol       ,solvVel->maxIt     
           ,sistEqVel->storage ,solvVel->solver
           ,solvVel->fileSolv  ,solvVel->log  
           ,true               ,sistEqVel->unsym   );  
       tm.solvVel = getTimeC() - tm.solvVel;
     }
/*...................................................................*/

/*... solver Aw = bw (velocidade estimadas)*/
     if(zMomentum && ndfVel == 3){
       if(fPrint) printf("Quantidade de movimento u3:\n");
       tm.solvVel = getTimeC() - tm.solvVel;
       solverC(m               
           ,sistEqVel->neq     ,sistEqVel->neqNov  
           ,sistEqVel->nad     ,sistEqVel->nadr
           ,sistEqVel->ia      ,sistEqVel->ja  
           ,sistEqVel->al      ,adU3           ,sistEqVel->au
           ,b3                 ,xu3
           ,&sistEqVel->iNeq   ,&sistEqVel->omp
           ,solvVel->tol       ,solvVel->maxIt     
           ,sistEqVel->storage ,solvVel->solver
           ,solvVel->fileSolv  ,solvVel->log  
           ,true               ,sistEqVel->unsym   );  
       tm.solvVel = getTimeC() - tm.solvVel;
     }
/*...................................................................*/

/*... atualizando o campo de velociade estimadas*/
     if (ndfVel == 3) {
       updateCellSimpleVel3D(mesh->elm.vel, xu1, xu2, xu3, sistEqVel->id
         , mesh->numelNov, mesh->ndm);
     }
     else {
       updateCellSimpleVel(mesh->elm.vel, xu1, xu2, sistEqVel->id
         , mesh->numelNov, mesh->ndm);
     }
/*...................................................................*/

/*...*/
     if(fPrint) printf("Correcao de pressao:\n");

/*... montagem do sistema  da pressao de correca*/
     tm.systFormPres = getTimeC() - tm.systFormPres;
     systFormSimplePres(loadsVel             ,loadsPresC
										,sc.diffPres
                    ,mesh->elm.node          ,mesh->elm.adj.nelcon  
                    ,mesh->elm.nen           ,mesh->elm.adj.nViz   
                    ,mesh->elm.geomType      ,mesh->elm.material.prop 
                    ,mesh->elm.material.type ,mesh->elm.mat   
                    ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi  
                    ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea    
                    ,mesh->elm.geom.normal   ,mesh->elm.geom.volume   
                    ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc    
                    ,mesh->elm.geom.vSkew    ,mesh->elm.geom.mvSkew   
                    ,mesh->elm.geom.dcca     ,mesh->elm.densityFluid
                    ,sistEqPres->ia          ,sistEqPres->ja      
                    ,sistEqPres->al          ,sistEqPres->ad       
                    ,bPc                     ,sistEqPres->id       
                    ,mesh->elm.faceRvel      ,mesh->elm.faceLoadVel 
                    ,mesh->elm.faceRpres     ,mesh->elm.faceLoadPres
                    ,mesh->elm.pressure      ,mesh->elm.gradPres           
                    ,mesh->elm.vel           ,sp->d                           
                    ,rCellPc                 ,sc.ddt
                    ,sistEqPres->neq         ,sistEqPres->neqNov      
                    ,sistEqPres->nad         ,sistEqPres->nadr      
                    ,mesh->maxNo             ,mesh->maxViz
                    ,mesh->ndm               ,mesh->numelNov
                    ,ndfVel                  ,sistEqPres->storage
                    ,true                    ,true   
                    ,true                    ,sistEqPres->unsym); 
     tm.systFormPres = getTimeC() - tm.systFormPres;
/*...................................................................*/
      
/*... residual*/
     residualSimple(mesh->elm.vel 
                   ,mesh->elm.rCellVel,rCellPc  
                   ,sistEqVel->ad           
                   ,rU                ,&rMass
                   ,sistEqVel->id     
                   ,mesh->numelNov    ,sistEqVel->neq  
                   ,mesh->ndm         ,3    );  
/*...................................................................*/

/*...*/
     pCor = true;
     if( rMass < tmp*SZERO ) pCor = false;
     if( itSimple == kZeroPres && relRes) rMass0 = rMass;
     if( itSimple == kZeroVel && relRes ){
       rU0[0] = rU[0]; 
       rU0[1] = rU[1]; 
       if (ndfVel == 3) rU0[2] = rU[2];
     } 
     conv = 0;
/*...................................................................*/

/*... solver ApPc = bpC (velocidade estimadas)*/
     if(pCor){
/*...*/
       zero(sp->ePresC ,mesh->numel  ,DOUBLEC);
       zero(sp->ePresC1,mesh->numel  ,DOUBLEC);
/*...................................................................*/

/*...*/
       tm.solvPres = getTimeC() - tm.solvPres;
       solverC(m               
           ,sistEqPres->neq    ,sistEqPres->neqNov  
           ,sistEqPres->nad    ,sistEqPres->nadr
           ,sistEqPres->ia     ,sistEqPres->ja  
           ,sistEqPres->al     ,sistEqPres->ad,sistEqPres->au
           ,bPc                ,xp
           ,&sistEqPres->iNeq  ,&sistEqPres->omp
           ,solvPres->tol      ,solvPres->maxIt     
           ,sistEqPres->storage,solvPres->solver
           ,solvPres->fileSolv ,solvPres->log  
           ,true               ,sistEqPres->unsym  );   
       tm.solvPres = getTimeC() - tm.solvPres;
/*...................................................................*/

/*... atualizando da pressao de correcao*/
       updateCellSimplePres(sp->ePresC,xp,sistEqPres->id,mesh->numelNov); 
/*...................................................................*/

/*...*/
       alphaProdVector(1.e0,sp->ePresC,mesh->numel,sp->ePresC1); 
/*...................................................................*/

/*... correcao nao ortoganal da pressao de correcao*/
       for(nonOrth=0;nonOrth < sp->nNonOrth;nonOrth++){
/*... reconstruindo do gradiente da pressao correcao*/
         tm.rcGradPres = getTimeC() - tm.rcGradPres;
         rcGradU(m                  ,loadsPresC
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
           ,mesh->elm.faceRpres     ,mesh->elm.faceLoadPres
           ,sp->ePresC1             ,sp->eGradPresC
           ,sp->nPresC              ,sc.rcGrad
           ,mesh->maxNo             ,mesh->maxViz
           ,1                       ,mesh->ndm
           ,&pMesh->iNo             ,&pMesh->iEl
           ,mesh->numelNov          ,mesh->numel
           ,mesh->nnodeNov          ,mesh->nnode);
         tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*...*/
         tm.systFormPres = getTimeC() - tm.systFormPres;
         simpleNonOrthPres(sc.diffPres
                    ,mesh->elm.node    ,mesh->elm.adj.nelcon
					          ,mesh->elm.nen           ,mesh->elm.adj.nViz
                    ,mesh->elm.geomType      ,mesh->elm.material.prop
                    ,mesh->elm.material.type ,mesh->elm.mat
                    ,mesh->elm.geom.cc 
                    ,mesh->elm.geom.ksi      ,mesh->elm.geom.mksi
                    ,mesh->elm.geom.eta      ,mesh->elm.geom.fArea
                    ,mesh->elm.geom.normal   ,mesh->elm.geom.volume
                    ,mesh->elm.geom.xm       ,mesh->elm.geom.xmcc
                    ,mesh->elm.geom.vSkew    ,mesh->elm.geom.mvSkew
                    ,mesh->elm.geom.dcca     ,mesh->elm.densityFluid
                    ,bPc                     ,sistEqPres->id
                    ,mesh->elm.faceRpres     ,sp->ePresC1     
                    ,sp->eGradPresC          ,sp->d
                    ,mesh->maxNo             ,mesh->maxViz
                    ,mesh->ndm               ,mesh->numelNov);
         tm.systFormPres = getTimeC() - tm.systFormPres;
/*...................................................................*/

/*...*/
         tm.solvPres = getTimeC() - tm.solvPres;
         solverC(m
            ,sistEqPres->neq    ,sistEqPres->neqNov
            ,sistEqPres->nad    ,sistEqPres->nadr
            ,sistEqPres->ia     ,sistEqPres->ja
            ,sistEqPres->al     ,sistEqPres->ad,sistEqPres->au
            ,bPc                ,xp
            ,&sistEqPres->iNeq  ,&sistEqPres->omp
            ,solvPres->tol      ,solvPres->maxIt
            ,sistEqPres->storage,solvPres->solver
            ,solvPres->fileSolv ,solvPres->log
            ,true               ,sistEqPres->unsym  );
         tm.solvPres = getTimeC() - tm.solvPres;  
/*...................................................................*/

/*... atualizando da pressao de correcao*/
         updateCellSimplePres(sp->ePresC1,xp,sistEqPres->id
                             ,mesh->numelNov);  
/*...................................................................*/

/*... soma o vetor presC(i) = presC + presC1*/
         addVector(1.0e0           ,sp->ePresC
                  ,1.0e0           ,sp->ePresC1
                  ,mesh->numel     ,sp->ePresC);
/*...................................................................*/
       }
/*...................................................................*/
     }
/*...................................................................*/

/*... reconstruindo do gradiente da pressao correcao*/
     tm.rcGradPres = getTimeC() - tm.rcGradPres;
     rcGradU(m                      ,loadsPresC
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
           ,mesh->elm.faceRpres     ,mesh->elm.faceLoadPres
           ,sp->ePresC              ,sp->eGradPresC
           ,sp->nPresC              ,sc.rcGrad
           ,mesh->maxNo             ,mesh->maxViz
           ,1                       ,mesh->ndm
           ,&pMesh->iNo             ,&pMesh->iEl
           ,mesh->numelNov          ,mesh->numel
           ,mesh->nnodeNov          ,mesh->nnode); 
     tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/


/*... atualizacao de u, v, w e p*/
     simpleUpdate(mesh->elm.vel,mesh->elm.pressure
                ,sp->ePresC,sp->eGradPresC
                ,sp->d       
                ,mesh->numel,mesh->ndm
                ,sp->alphaPres);
/*...................................................................*/

/*... 3D*/
     if (ndfVel == 3) {
       if (rMass / rMass0 < tolSimpleMass || rMass < tmp*SZERO) conv++;
/*..*/
       if (rU[0] / rU0[0] < tolSimpleU1 || rU[0] < tmp*SZERO) conv++;
/*...*/
       if (rU[1] / rU0[1] < tolSimpleU2 || rU[1] < tmp*SZERO) conv++;
/*...*/
       if (rU[2] / rU0[2] < tolSimpleU3 || rU[2] < tmp*SZERO) conv++;
/*..*/
       if (conv == 4) break;
     }
/*...................................................................*/

/*... 2D*/
     else {
/*...*/
         if (rMass / rMass0 < tolSimpleMass || rMass < tmp*SZERO) conv++;
/*..*/
       if (rU[0] / rU0[0] < tolSimpleU1 || rU[0] < tmp*SZERO) conv++;
/*...*/
       if (rU[1] / rU0[1] < tolSimpleU2 || rU[1] < tmp*SZERO) conv++;
/*..*/
       if (conv == 3) break;
/*...................................................................*/
     }
/*...................................................................*/

/*...*/
     timei = getTimeC() -time;
/*... arquivo de log*/
     if(opt.fItPlot){
       if (ndfVel == 3)
         fprintf(opt.fileItPlot[FITPLOTSIMPLE]
                , "%d %20.8e %20.8e %20.8e %20.8e\n"
                ,itSimple + 1, rU[0], rU[1], rU[2], rMass);
       else
         fprintf(opt.fileItPlot[FITPLOTSIMPLE]
                ,"%d %20.8e %20.8e %20.8e\n"
                ,itSimple+1,rU[0],rU[1],rMass);
     }
/*...................................................................*/

/*...*/
     if( jj == sp->pSimple) {
       jj = 0; 
       printf("It simple: %d \n",itSimple+1);
       printf("Time(s)  : %lf \n",timei);
       printf("Residuo:\n");
       printf("conservacao da massa: %20.8e\n",rMass/rMass0);
       printf("momentum x1         : %20.8e\n",rU[0]/rU0[0]);
       printf("momentum x2         : %20.8e\n",rU[1]/rU0[1]);
       if (ndfVel == 3)
         printf("momentum x3         : %20.8e\n",rU[2]/rU0[2]);
     } 
     jj++; 
/*...................................................................*/
  }
/*...................................................................*/
  timei = getTimeC() -time;

/*...*/  
  fParameter[0] = true;
  fParameter[1] = true;
  parameterCell(mesh->elm.vel           ,mesh->elm.material.prop
               ,mesh->elm.densityFluid  ,mesh->elm.geom.volume 
               ,mesh->elm.mat  
               ,&cfl                    ,&reynolds
               ,fParameter              ,sc.ddt.dt[0]
               ,mesh->numelNov          ,mesh->ndm);
/*...................................................................*/

/*...*/
  printf("It simple: %d \n",itSimple+1);
  printf("Time(s)  : %lf \n",timei);
  printf("Reynolds: %lf\n",reynolds);
  if(sc.ddt.flag)
    printf("CFL     : %lf\n",cfl);
  printf("Residuo:\n");
  printf("conservacao da massa (init,final): %20.8e %20.8e \n"
        ,rMass0,rMass);
  printf("momentum x1          (init,final): %20.8e %20.8e\n"
        ,rU0[0],rU[0]);
  printf("momentum x2          (init,final): %20.8e %20.8e\n"
        ,rU0[1],rU[1]);
  printf("momentum x3          (init,final): %20.8e %20.8e\n"
        ,rU0[2],rU[2]);
/*...................................................................*/

} 
/*********************************************************************/ 

/*********************************************************************
* Data de criacao    : 21/08/2017                                   *
* Data de modificaco : 03/09/2017                                   *
*-------------------------------------------------------------------*
* SIMPLESOLVERENERGY: metodo simple e simpleC para escoamentos      * 
* 2D/3D termo ativados                                              *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
*********************************************************************/
void simpleSolverLm(Memoria *m          ,PropVar prop     
                   ,Loads *loadsVel     ,Loads *loadsPres
                   ,Loads *loadsEnergy  ,EnergyModel eModel
                   ,Mesh *mesh0         ,Mesh *mesh
                   ,SistEq *sistEqVel   ,SistEq *sistEqPres
                   ,SistEq *sistEqEnergy
                   ,Solv *solvVel       ,Solv *solvPres
                   ,Solv *solvEnergy
                   ,Simple *sp
                   ,Scheme sc           ,PartMesh *pMesh
                   ,FileOpt opt         ,char *preName
                   ,char *nameOut       ,FILE *fileOut) {
  FILE *fStop = NULL;
  short unsigned ndfVel = mesh->ndfFt - 2;
  short unsigned conv;
  int itSimple;
  int nonOrth;
  short unsigned kZeroVel    = sp->kZeroVel,
                 kZeroPres   = sp->kZeroPres,
                 kZeroEnergy = sp->kZeroEnergy;
  INT jj = 1;
  DOUBLE time, timei;
  DOUBLE *b1, *b2, *b3, *bPc, *bE, *xu1, *xu2, *xu3, *xp;
  DOUBLE *adU1, *adU2, *adU3;
  DOUBLE *rCellPc,*rCellE;
/*...*/
  DOUBLE rU[3], rU0[3], tmp, tb[3], rMass0, rMass,rEnergy0, rEnergy;
  /*...*/
  DOUBLE tolSimpleU1, tolSimpleU2, tolSimpleU3, tolSimpleMass, tolSimpleEnergy;
/*...*/
  bool xMomentum, yMomentum, zMomentum, pCor, fEnergy;
  bool relRes = false;
  bool fPrint = false;
  bool fDensity = prop.fDensity,
       fSheat   = prop.fSpecificHeat,
       fDvisc   = prop.fDynamicViscosity,
       fTcond   = prop.fThermalCondutivty;
  DOUBLE cfl, reynolds, peclet;
  bool fParameter[3];

  time = getTimeC();
/*...*/
  b1 = b2 = b3 = bPc = NULL;
  xu1 = xu2 = xu3 = NULL;
  adU1 = adU2 = adU3 = NULL;

  b1 = sistEqVel->b;
  b2 = &sistEqVel->b[sistEqVel->neq];
  if (ndfVel == 3) b3 = &sistEqVel->b[2 * sistEqVel->neq];
  bPc = sistEqPres->b;
  bE  = sistEqEnergy->b;

  xu1 = sistEqVel->x;
  xu2 = &sistEqVel->x[sistEqVel->neq];
  if (ndfVel == 3) xu3 = &sistEqVel->x[2 * sistEqVel->neq];
  xp = sistEqPres->x;

  adU1 = sistEqVel->ad;
  adU2 = &sistEqVel->ad[sistEqVel->neq];
  if (ndfVel == 3) adU3 = &sistEqVel->ad[2 * sistEqVel->neq];

  rCellPc = mesh->elm.rCellPres;
  rCellE  = mesh->elm.rCellEnergy;
/*...................................................................*/

/*...*/
  tolSimpleU1 = tolSimpleU2 = tolSimpleU3 = sp->tolVel;
  tolSimpleMass = sp->tolPres;
  tolSimpleEnergy = sp->tolEnergy;
/*...................................................................*/

/*
*/

/*...*/
  rMass0 = 1.e0;
  rMass = 0.e0;
  rU[0] = rU[1] = rU[2] = 0.e0;
  rU0[0] = rU0[1] = rU0[2] = 1.e0;
  rEnergy0 = 1.e0;
  rEnergy  = 0.e0;
  conv = 0;
  xMomentum = yMomentum = zMomentum = true;
/*...................................................................*/

/*...*/
  zero(sistEqVel->b0   ,sistEqVel->neqNov*ndfVel,DOUBLEC);
  zero(sistEqPres->b0  ,sistEqPres->neqNov      ,DOUBLEC);
  zero(sistEqEnergy->b0,sistEqEnergy->neqNov    ,DOUBLEC);
/*...................................................................*/

/*... restricoes por centro de celula u0 e cargas por volume b0*/
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
  cellPloadSimple(loadsPres            ,mesh->elm.geom.cc
                 ,mesh->elm.faceRpres  ,mesh->elm.faceLoadPres
                 ,mesh->elm.geom.volume
                 ,sistEqVel->id        ,sistEqPres->id
                 ,mesh->elm.vel        ,mesh->elm.pressure
                 ,sistEqVel->b0        ,sistEqPres->b0
                 ,mesh->numelNov       ,ndfVel
                 ,mesh->ndm            ,mesh->maxViz);
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
/*...................................................................*/

/*... restricoes por centro de celula u0 e cargas por volume b0*/
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
  cellPload(loadsEnergy            ,mesh->elm.geom.cc
           ,mesh->elm.faceRenergy  ,mesh->elm.faceLoadEnergy
           ,mesh->elm.geom.volume  ,sistEqEnergy->id
           ,mesh->elm.energy       ,sistEqEnergy->b0
           ,mesh->numelNov         ,ndfVel
           ,mesh->ndm              ,mesh->maxViz);
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
/*...................................................................*/

/*... discretizacao temporal*/
  if (sc.ddt.flag) {
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
/*... velocidade*/
    cellTransientSimple(mesh->elm.geom.volume ,sistEqVel->id
                       ,mesh->elm.vel0        ,mesh->elm.vel
                       ,mesh->elm.densityFluid,sistEqVel->b0
                       ,sc.ddt                ,sistEqVel->neqNov
                       ,mesh->numelNov        ,ndfVel
                       ,true);
/*... Energia*/
    cellTransientEnergy(mesh->elm.geom.volume  ,sistEqEnergy->id
                       ,mesh->elm.energy0      ,mesh->elm.energy
                       ,mesh->elm.densityFluid ,mesh->elm.specificHeat
                       ,sistEqEnergy->b0
                       ,sc.ddt                 ,mesh->numelNov
                       ,true);
/*... vel(n-1) = vel(n)*/
    alphaProdVector(1.e0              ,mesh->elm.vel
                   ,mesh->numel*ndfVel,mesh->elm.vel0);
/*... energy(n-1) = energy(n)*/
    alphaProdVector(1.e0       ,mesh->elm.energy
                   ,mesh->numel,mesh->elm.energy0);
/*... pres(n-1) = pres(n)*/
    alphaProdVector(1.e0       ,mesh->elm.pressure
                   ,mesh->numel,mesh->elm.pressure0);
/*...................................................................*/
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
  }
/*...................................................................*/

/*... calculo da matrix jacobiana das velocidades
    | du1dx1 du1dx2 du1dx3 |
    | du2dx1 du2dx2 du2dx3 |
    | du3dx1 du3dx2 du3dx3 |
*/
  tm.rcGradVel = getTimeC() - tm.rcGradVel;
  rcGradU(m, loadsVel
            , mesh->elm.node, mesh->elm.adj.nelcon
            , mesh->elm.geom.cc, mesh->node.x
            , mesh->elm.nen, mesh->elm.adj.nViz
            , mesh->elm.geomType, mesh->elm.material.prop
            , mesh->elm.mat
            , mesh->elm.leastSquare, mesh->elm.leastSquareR
            , mesh->elm.geom.ksi, mesh->elm.geom.mksi
            , mesh->elm.geom.eta, mesh->elm.geom.fArea
            , mesh->elm.geom.normal, mesh->elm.geom.volume
            , mesh->elm.geom.vSkew
            , mesh->elm.geom.xm, mesh->elm.geom.xmcc
            , mesh->elm.geom.dcca
            , mesh->elm.faceRvel, mesh->elm.faceLoadVel
            , mesh->elm.vel, mesh->elm.gradVel
            , mesh->node.vel, sc.rcGrad
            , mesh->maxNo, mesh->maxViz
            , ndfVel, mesh->ndm
            , &pMesh->iNo, &pMesh->iEl
            , mesh->numelNov, mesh->numel
            , mesh->nnodeNov, mesh->nnode);
  tm.rcGradVel = getTimeC() - tm.rcGradVel;
/*...................................................................*/

/*... reconstruindo do gradiente da pressao*/
  tm.rcGradPres = getTimeC() - tm.rcGradPres;
  rcGradU(m, loadsPres
           ,mesh->elm.node       ,mesh->elm.adj.nelcon
           ,mesh->elm.geom.cc    ,mesh->node.x
           ,mesh->elm.nen        ,mesh->elm.adj.nViz
           ,mesh->elm.geomType   ,mesh->elm.material.prop
           ,mesh->elm.mat
           ,mesh->elm.leastSquare,mesh->elm.leastSquareR
           ,mesh->elm.geom.ksi   ,mesh->elm.geom.mksi
           ,mesh->elm.geom.eta   ,mesh->elm.geom.fArea
           ,mesh->elm.geom.normal,mesh->elm.geom.volume
           ,mesh->elm.geom.vSkew
           ,mesh->elm.geom.xm    ,mesh->elm.geom.xmcc
           ,mesh->elm.geom.dcca
           ,mesh->elm.faceRpres  ,mesh->elm.faceLoadPres
           ,mesh->elm.pressure   ,mesh->elm.gradPres
           ,mesh->node.pressure  ,sc.rcGrad
           ,mesh->maxNo          ,mesh->maxViz
           ,1                    ,mesh->ndm
           ,&pMesh->iNo          ,&pMesh->iEl
           ,mesh->numelNov       ,mesh->numel
           ,mesh->nnodeNov       ,mesh->nnode);
  tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*...*/
  for (itSimple = 0; itSimple<sp->maxIt; itSimple++) {
/*...*/
    if ((fStop = fopen("stopSimple.mvf", "r")) != NULL) {
      fclose(fStop);
      break;
    }
/*...................................................................*/

/*... montagem do sistema u, v e w*/
    tm.systFormVel = getTimeC() - tm.systFormVel;
    systFormSimpleVelLw(loadsVel               ,loadsPres
                     ,sc.advVel              ,sc.diffVel
                     ,sp->type
                     ,mesh->elm.node         ,mesh->elm.adj.nelcon
                     ,mesh->elm.nen          ,mesh->elm.adj.nViz
                     ,mesh->elm.geomType     ,mesh->elm.material.prop
                     ,mesh->elm.material.type,mesh->elm.mat
                     ,mesh->elm.geom.cc      ,mesh->elm.geom.ksi
                     ,mesh->elm.geom.mksi    ,mesh->elm.geom.eta     
                     ,mesh->elm.geom.fArea   ,mesh->elm.geom.normal  
                     ,mesh->elm.geom.volume  ,mesh->elm.geom.xm    
                     ,mesh->elm.geom.xmcc    ,mesh->elm.geom.vSkew   
                     ,mesh->elm.geom.mvSkew  ,mesh->elm.geom.dcca    
                     ,sistEqVel->ia          ,sistEqVel->ja
                     ,sistEqVel->al          ,sistEqVel->ad
                     ,sistEqVel->b           ,sistEqVel->id
                     ,mesh->elm.faceRvel     ,mesh->elm.faceLoadVel
                     ,mesh->elm.faceRpres    ,mesh->elm.faceLoadPres
                     ,mesh->elm.pressure     ,mesh->elm.gradPres
                     ,mesh->elm.vel          ,mesh->elm.gradVel
                     ,sp->d                  ,sp->alphaVel
                     ,mesh->elm.rCellVel
                     ,mesh->elm.densityFluid ,mesh->elm.dViscosity    
                     ,sc.ddt
                     ,sistEqVel->neq         ,sistEqVel->neqNov
                     ,sistEqVel->nad         ,sistEqVel->nadr
                     ,mesh->maxNo            ,mesh->maxViz
                     ,mesh->ndm              ,mesh->numelNov
                     ,ndfVel                 ,sistEqVel->storage
                     ,true                   ,true
                     ,true                   ,sistEqVel->unsym
                     ,sp->sPressure);
    tm.systFormVel = getTimeC() - tm.systFormVel;
/*...................................................................*/

/*... soma o vetor b(i) = b(i) + b0(i)*/
    addVector(1.0e0                   ,sistEqVel->b
             ,1.0e0                   ,sistEqVel->b0
             ,sistEqVel->neqNov*ndfVel,sistEqVel->b);
/*...................................................................*/

/*... soma o vetor R(i) = R(i) + b0(i)*/
   updateCellValueSimple(mesh->elm.rCellVel,sistEqVel->b0
                         ,sistEqVel->id     ,&sistEqVel->iNeq
                         ,mesh->numelNov    ,sistEqVel->neqNov
                         ,ndfVel
                         ,true              ,false);  
/*...................................................................*/

/*...*/
    tb[0] = sqrt(dot(b1, b1, sistEqVel->neqNov));
    tb[1] = sqrt(dot(b2, b2, sistEqVel->neqNov));
    if (ndfVel == 3) tb[2] = sqrt(dot(b3, b3, sistEqVel->neqNov));
    if (itSimple == 0) {
      tmp = max(tb[0], tb[1]);
      if (ndfVel == 3) tmp = max(tmp, tb[2]);
    }
/*...*/
    xMomentum = true;
    if (tb[0] <= tmp*SZERO || tb[0] == 0.e0 ) xMomentum = false;
/*...................................................................*/

/*...*/
    yMomentum = true;
    if (tb[1] <= tmp*SZERO || tb[1] == 0.e0) yMomentum = false;
/*...................................................................*/

/*...*/
    if (ndfVel == 3) {
      zMomentum = true;
      if (tb[2] <= tmp*SZERO || tb[2] == 0.e0 ) zMomentum = false;
    }
/*...................................................................*/

/*... solver Au = bu (velocidade estimadas)*/
    if (xMomentum) {
      if (fPrint) printf("Quantidade de movimento u1:\n");
      tm.solvVel = getTimeC() - tm.solvVel;
      solverC(m
             ,sistEqVel->neq    ,sistEqVel->neqNov
             ,sistEqVel->nad    ,sistEqVel->nadr
             ,sistEqVel->ia     ,sistEqVel->ja
             ,sistEqVel->al     ,adU1             ,sistEqVel->au
             ,b1                ,xu1
             ,&sistEqVel->iNeq  ,&sistEqVel->omp
             ,solvVel->tol      ,solvVel->maxIt
             ,sistEqVel->storage,solvVel->solver
             ,solvVel->fileSolv ,solvVel->log
             ,true              ,sistEqVel->unsym);  
      tm.solvVel = getTimeC() - tm.solvVel;
    }
/*...................................................................*/

/*... solver Av = bv (velocidade estimadas)*/
    if (yMomentum) {
      if (fPrint) printf("Quantidade de movimento u2:\n");
      tm.solvVel = getTimeC() - tm.solvVel;
      solverC(m
             ,sistEqVel->neq    ,sistEqVel->neqNov
             ,sistEqVel->nad    ,sistEqVel->nadr
             ,sistEqVel->ia     ,sistEqVel->ja
             ,sistEqVel->al     ,adU2             ,sistEqVel->au
             ,b2                ,xu2
             ,&sistEqVel->iNeq  ,&sistEqVel->omp
             ,solvVel->tol      ,solvVel->maxIt
             ,sistEqVel->storage,solvVel->solver
             ,solvVel->fileSolv ,solvVel->log
             ,true              ,sistEqVel->unsym);  
      tm.solvVel = getTimeC() - tm.solvVel;
    }
/*...................................................................*/

/*... solver Aw = bw (velocidade estimadas)*/
    if (zMomentum && ndfVel == 3) {
      if (fPrint) printf("Quantidade de movimento u3:\n");
      tm.solvVel = getTimeC() - tm.solvVel;
      solverC(m
             ,sistEqVel->neq    ,sistEqVel->neqNov
             ,sistEqVel->nad    ,sistEqVel->nadr
             ,sistEqVel->ia     ,sistEqVel->ja
             ,sistEqVel->al     ,adU3             ,sistEqVel->au
             ,b3                ,xu3
             ,&sistEqVel->iNeq  ,&sistEqVel->omp
             ,solvVel->tol      ,solvVel->maxIt
             ,sistEqVel->storage,solvVel->solver
             ,solvVel->fileSolv ,solvVel->log
             ,true              ,sistEqVel->unsym);  
      tm.solvVel = getTimeC() - tm.solvVel;
    }
/*...................................................................*/

/*... atualizando o campo de velociade estimadas*/
//  if (ndfVel == 3) {
//    updateCellSimpleVel3D(mesh->elm.vel, xu1, xu2, xu3, sistEqVel->id
//      , mesh->numelNov, mesh->ndm);
//  }
//  else {
//    updateCellSimpleVelR(mesh->elm.vel ,xu1, xu2, sistEqVel->id
//                        ,mesh->numelNov,mesh->ndm);
//  }
    updateCellSimpleVelR(mesh->elm.vel ,xu1
                        ,xu2           ,xu3
                        ,sistEqVel->id ,mesh->numelNov
                        ,mesh->ndm);
/*...................................................................*/

/*...*/
    if (fPrint) printf("Correcao de pressao:\n");

/*... montagem do sistema  da pressao de correca*/
    tm.systFormPres = getTimeC() - tm.systFormPres;
    systFormSimplePres(loadsVel               ,loadsPresC
                      ,sc.diffPres
                      ,mesh->elm.node         ,mesh->elm.adj.nelcon
                      ,mesh->elm.nen          ,mesh->elm.adj.nViz
                      ,mesh->elm.geomType     ,mesh->elm.material.prop
                      ,mesh->elm.material.type,mesh->elm.mat
                      ,mesh->elm.geom.ksi     ,mesh->elm.geom.mksi
                      ,mesh->elm.geom.eta     ,mesh->elm.geom.fArea
                      ,mesh->elm.geom.normal  ,mesh->elm.geom.volume
                      ,mesh->elm.geom.xm      ,mesh->elm.geom.xmcc
                      ,mesh->elm.geom.vSkew   ,mesh->elm.geom.mvSkew
                      ,mesh->elm.geom.dcca    ,mesh->elm.densityFluid
                      ,sistEqPres->ia         ,sistEqPres->ja
                      ,sistEqPres->al         ,sistEqPres->ad
                      ,bPc                    ,sistEqPres->id
                      ,mesh->elm.faceRvel     ,mesh->elm.faceLoadVel
                      ,mesh->elm.faceRpres    ,mesh->elm.faceLoadPres
                      ,mesh->elm.pressure     ,mesh->elm.gradPres
                      ,mesh->elm.vel          ,sp->d
                      ,rCellPc                ,sc.ddt
                      ,sistEqPres->neq        ,sistEqPres->neqNov
                      ,sistEqPres->nad        ,sistEqPres->nadr
                      ,mesh->maxNo            ,mesh->maxViz
                      ,mesh->ndm              ,mesh->numelNov
                      ,ndfVel                 ,sistEqPres->storage
                      ,true                   ,true
                      ,true                   ,sistEqPres->unsym);  
    tm.systFormPres = getTimeC() - tm.systFormPres;
/*...................................................................*/

/*... residual*/
    tmp = sqrt(dot(bPc,bPc,sistEqPres->neqNov));
/*...................................................................*/

/*...*/
    pCor = true;
    if (tmp <= SZERO || tmp == 0.e0) pCor = false;
/*...................................................................*/

/*... solver ApPc = bpC (velocidade estimadas)*/
    if (pCor) {
/*...*/
      zero(sp->ePresC, mesh->numel, DOUBLEC);
      zero(sp->ePresC1, mesh->numel, DOUBLEC);
/*...................................................................*/

/*...*/
      tm.solvPres = getTimeC() - tm.solvPres;
      solverC(m
             ,sistEqPres->neq    ,sistEqPres->neqNov
             ,sistEqPres->nad    ,sistEqPres->nadr
             ,sistEqPres->ia     ,sistEqPres->ja
             ,sistEqPres->al     ,sistEqPres->ad, sistEqPres->au
             ,bPc                ,xp
             ,&sistEqPres->iNeq  ,&sistEqPres->omp
             ,solvPres->tol      ,solvPres->maxIt
             ,sistEqPres->storage,solvPres->solver
             ,solvPres->fileSolv ,solvPres->log
             ,true               ,sistEqPres->unsym);  
      tm.solvPres = getTimeC() - tm.solvPres;
/*...................................................................*/

/*... atualizando da pressao de correcao*/
      updateCellSimplePres(sp->ePresC, xp, sistEqPres->id, mesh->numelNov);
/*...................................................................*/

/*...*/
      alphaProdVector(1.e0, sp->ePresC, mesh->numel, sp->ePresC1);
/*...................................................................*/

/*... correcao nao ortoganal da pressao de correcao*/
      for (nonOrth = 0; nonOrth < sp->nNonOrth; nonOrth++) {
/*... reconstruindo do gradiente da pressao correcao*/
        tm.rcGradPres = getTimeC() - tm.rcGradPres;
        rcGradU(m, loadsPresC
                ,mesh->elm.node       ,mesh->elm.adj.nelcon
                ,mesh->elm.geom.cc    ,mesh->node.x
                ,mesh->elm.nen        ,mesh->elm.adj.nViz
                ,mesh->elm.geomType   ,mesh->elm.material.prop
                ,mesh->elm.mat
                ,mesh->elm.leastSquare,mesh->elm.leastSquareR
                ,mesh->elm.geom.ksi   ,mesh->elm.geom.mksi
                ,mesh->elm.geom.eta   ,mesh->elm.geom.fArea
                ,mesh->elm.geom.normal,mesh->elm.geom.volume
                ,mesh->elm.geom.vSkew
                ,mesh->elm.geom.xm    ,mesh->elm.geom.xmcc
                ,mesh->elm.geom.dcca
                ,mesh->elm.faceRpres  ,mesh->elm.faceLoadPres
                ,sp->ePresC1          ,sp->eGradPresC
                ,sp->nPresC, sc.rcGrad
                ,mesh->maxNo          , mesh->maxViz
                ,1                    ,mesh->ndm
                ,&pMesh->iNo          ,&pMesh->iEl
                ,mesh->numelNov       ,mesh->numel
                ,mesh->nnodeNov       ,mesh->nnode);  
        tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*...*/
        tm.systFormPres = getTimeC() - tm.systFormPres;
        simpleNonOrthPres(sc.diffPres
                         ,mesh->elm.node         ,mesh->elm.adj.nelcon
                         ,mesh->elm.nen          ,mesh->elm.adj.nViz
                         ,mesh->elm.geomType     ,mesh->elm.material.prop
                         ,mesh->elm.material.type,mesh->elm.mat
                         ,mesh->elm.geom.cc
                         ,mesh->elm.geom.ksi     ,mesh->elm.geom.mksi
                         ,mesh->elm.geom.eta     ,mesh->elm.geom.fArea
                         ,mesh->elm.geom.normal  ,mesh->elm.geom.volume
                         ,mesh->elm.geom.xm      ,mesh->elm.geom.xmcc
                         ,mesh->elm.geom.vSkew   ,mesh->elm.geom.mvSkew
                         ,mesh->elm.geom.dcca    ,mesh->elm.densityFluid
                         ,bPc                    ,sistEqPres->id
                         ,mesh->elm.faceRpres    ,sp->ePresC1
                         ,sp->eGradPresC         ,sp->d
                         ,mesh->maxNo            ,mesh->maxViz
                         ,mesh->ndm              ,mesh->numelNov);
        tm.systFormPres = getTimeC() - tm.systFormPres;
/*...................................................................*/

/*...*/
        tm.solvPres = getTimeC() - tm.solvPres;
        solverC(m
               ,sistEqPres->neq    ,sistEqPres->neqNov
               ,sistEqPres->nad    ,sistEqPres->nadr
               ,sistEqPres->ia     ,sistEqPres->ja
               ,sistEqPres->al     ,sistEqPres->ad, sistEqPres->au
               ,bPc                ,xp
               ,&sistEqPres->iNeq  ,&sistEqPres->omp
               ,solvPres->tol      ,solvPres->maxIt
               ,sistEqPres->storage,solvPres->solver
               ,solvPres->fileSolv ,solvPres->log
               ,true               ,sistEqPres->unsym);
        tm.solvPres = getTimeC() - tm.solvPres;
/*...................................................................*/

/*... atualizando da pressao de correcao*/
        updateCellSimplePres(sp->ePresC1, xp, sistEqPres->id
                            ,mesh->numelNov);
/*...................................................................*/

/*... soma o vetor presC(i) = presC + presC1*/
        addVector(1.0e0      ,sp->ePresC
                 ,1.0e0      ,sp->ePresC1
                 ,mesh->numel,sp->ePresC);
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... reconstruindo do gradiente da pressao correcao*/
    tm.rcGradPres = getTimeC() - tm.rcGradPres;
    rcGradU(m, loadsPresC
           ,mesh->elm.node, mesh->elm.adj.nelcon
           ,mesh->elm.geom.cc, mesh->node.x
           ,mesh->elm.nen, mesh->elm.adj.nViz
           ,mesh->elm.geomType, mesh->elm.material.prop
           ,mesh->elm.mat
           ,mesh->elm.leastSquare, mesh->elm.leastSquareR
           ,mesh->elm.geom.ksi, mesh->elm.geom.mksi
           ,mesh->elm.geom.eta, mesh->elm.geom.fArea
           ,mesh->elm.geom.normal, mesh->elm.geom.volume
           ,mesh->elm.geom.vSkew
           ,mesh->elm.geom.xm, mesh->elm.geom.xmcc
           ,mesh->elm.geom.dcca
           ,mesh->elm.faceRpres, mesh->elm.faceLoadPres
           ,sp->ePresC, sp->eGradPresC
           ,sp->nPresC, sc.rcGrad
           ,mesh->maxNo, mesh->maxViz
           ,1, mesh->ndm
           ,&pMesh->iNo, &pMesh->iEl
           ,mesh->numelNov, mesh->numel
           ,mesh->nnodeNov, mesh->nnode);
    tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*... atualizacao de u, v, w e p*/
    simpleUpdate(mesh->elm.vel ,mesh->elm.pressure
                ,sp->ePresC    ,sp->eGradPresC
                ,sp->d
                ,mesh->numel   ,mesh->ndm
                ,sp->alphaPres);
/*...................................................................*/

/*... calculo da matrix jacobiana das velocidades
    | du1dx1 du1dx2 du1dx3 |
    | du2dx1 du2dx2 du2dx3 |
    | du3dx1 du3dx2 du3dx3 |
*/
    tm.rcGradVel = getTimeC() - tm.rcGradVel;
    rcGradU(m, loadsVel
            , mesh->elm.node, mesh->elm.adj.nelcon
            , mesh->elm.geom.cc, mesh->node.x
            , mesh->elm.nen, mesh->elm.adj.nViz
            , mesh->elm.geomType, mesh->elm.material.prop
            , mesh->elm.mat
            , mesh->elm.leastSquare, mesh->elm.leastSquareR
            , mesh->elm.geom.ksi, mesh->elm.geom.mksi
            , mesh->elm.geom.eta, mesh->elm.geom.fArea
            , mesh->elm.geom.normal, mesh->elm.geom.volume
            , mesh->elm.geom.vSkew
            , mesh->elm.geom.xm, mesh->elm.geom.xmcc
            , mesh->elm.geom.dcca
            , mesh->elm.faceRvel, mesh->elm.faceLoadVel
            , mesh->elm.vel, mesh->elm.gradVel
            , mesh->node.vel, sc.rcGrad
            , mesh->maxNo, mesh->maxViz
            , ndfVel, mesh->ndm
            , &pMesh->iNo, &pMesh->iEl
            , mesh->numelNov, mesh->numel
            , mesh->nnodeNov, mesh->nnode);
    tm.rcGradVel = getTimeC() - tm.rcGradVel;
/*...................................................................*/

/*... reconstruindo do gradiente da pressao*/
    tm.rcGradPres = getTimeC() - tm.rcGradPres;
    rcGradU(m, loadsPres
           ,mesh->elm.node       ,mesh->elm.adj.nelcon
           ,mesh->elm.geom.cc    ,mesh->node.x
           ,mesh->elm.nen        ,mesh->elm.adj.nViz
           ,mesh->elm.geomType   ,mesh->elm.material.prop
           ,mesh->elm.mat
           ,mesh->elm.leastSquare,mesh->elm.leastSquareR
           ,mesh->elm.geom.ksi   ,mesh->elm.geom.mksi
           ,mesh->elm.geom.eta   ,mesh->elm.geom.fArea
           ,mesh->elm.geom.normal,mesh->elm.geom.volume
           ,mesh->elm.geom.vSkew
           ,mesh->elm.geom.xm    ,mesh->elm.geom.xmcc
           ,mesh->elm.geom.dcca
           ,mesh->elm.faceRpres  ,mesh->elm.faceLoadPres
           ,mesh->elm.pressure   ,mesh->elm.gradPres
           ,mesh->node.pressure  ,sc.rcGrad
           ,mesh->maxNo          ,mesh->maxViz
           ,1                    ,mesh->ndm
           ,&pMesh->iNo          ,&pMesh->iEl
           ,mesh->numelNov       ,mesh->numel
           ,mesh->nnodeNov       ,mesh->nnode);
    tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*... reconstruindo do gradiente (Energia)*/
    tm.rcGradEnergy = getTimeC() - tm.rcGradEnergy;
    rcGradU(m                    ,loadsEnergy
           ,mesh->elm.node       ,mesh->elm.adj.nelcon
           ,mesh->elm.geom.cc    ,mesh->node.x
           ,mesh->elm.nen        ,mesh->elm.adj.nViz
           ,mesh->elm.geomType   ,mesh->elm.material.prop
           ,mesh->elm.mat
           ,mesh->elm.leastSquare,mesh->elm.leastSquareR
           ,mesh->elm.geom.ksi   ,mesh->elm.geom.mksi
           ,mesh->elm.geom.eta   ,mesh->elm.geom.fArea
           ,mesh->elm.geom.normal,mesh->elm.geom.volume
           ,mesh->elm.geom.vSkew
           ,mesh->elm.geom.xm    ,mesh->elm.geom.xmcc
           ,mesh->elm.geom.dcca
           ,mesh->elm.faceRenergy,mesh->elm.faceLoadEnergy
           ,mesh->elm.energy     ,mesh->elm.gradEnergy
           ,mesh->node.energy    ,sc.rcGrad
           ,mesh->maxNo          ,mesh->maxViz
           ,1, mesh->ndm
           ,&pMesh->iNo          ,&pMesh->iEl
           ,mesh->numelNov       ,mesh->numel
           ,mesh->nnodeNov       ,mesh->nnode);
    tm.rcGradEnergy = getTimeC() - tm.rcGradEnergy;
/*.................................................................. */


/*... equacao de energia*/
/*... calculo de: A(i),bE(i)*/
    tm.systFormEnergy = getTimeC() - tm.systFormEnergy;
    systFormEnergy(loadsEnergy       ,eModel
             ,sc.advEnergy           ,sc.diffEnergy
             ,mesh->elm.node         ,mesh->elm.adj.nelcon
             ,mesh->elm.nen          ,mesh->elm.adj.nViz
             ,mesh->elm.geomType     ,mesh->elm.material.prop
             ,mesh->elm.material.type,mesh->elm.mat
             ,mesh->elm.geom.cc
             ,mesh->elm.geom.ksi     ,mesh->elm.geom.mksi
             ,mesh->elm.geom.eta     ,mesh->elm.geom.fArea
             ,mesh->elm.geom.normal  ,mesh->elm.geom.volume
             ,mesh->elm.geom.xm      ,mesh->elm.geom.xmcc
             ,mesh->elm.geom.vSkew   ,mesh->elm.geom.mvSkew
             ,mesh->elm.geom.dcca    
             ,sistEqEnergy->ia       ,sistEqEnergy->ja
             ,sistEqEnergy->al       ,sistEqEnergy->ad
             ,sistEqEnergy->b        ,sistEqEnergy->id
             ,mesh->elm.faceRenergy  ,mesh->elm.faceLoadEnergy
             ,mesh->elm.energy       ,mesh->elm.gradEnergy
             ,mesh->elm.vel          ,mesh->elm.gradVel
             ,mesh->elm.pressure0    ,mesh->elm.pressure  
             ,mesh->elm.gradPres
             ,mesh->elm.rCellEnergy 
             ,mesh->elm.densityFluid ,mesh->elm.specificHeat
             ,mesh->elm.dViscosity   ,mesh->elm.tConductivity
             ,sc.ddt                 ,sp->alphaEnergy
             ,sistEqEnergy->neq      ,sistEqEnergy->neqNov
             ,sistEqEnergy->nad      ,sistEqEnergy->nadr
             ,mesh->maxNo            ,mesh->maxViz
             ,mesh->ndm              ,mesh->numelNov
             ,1                      ,sistEqEnergy->storage
             ,true                   ,true
             ,true                   ,sistEqEnergy->unsym);
    tm.systFormEnergy = getTimeC() - tm.systFormEnergy;
/*...................................................................*/

/*... soma o vetor bE(i) = bE(i) + bE0(i)*/
    addVector(1.0e0               ,sistEqEnergy->b
             ,1.0e0               ,sistEqEnergy->b0
             ,sistEqEnergy->neqNov,sistEqEnergy->b);
/*...................................................................*/

/*... soma o vetor RE(i) = RE(i) + bE0(i)*/
    updateCellValue(mesh->elm.rCellEnergy,sistEqEnergy->b0
                   ,sistEqEnergy->id     ,&sistEqEnergy->iNeq
                   ,mesh->numelNov       ,1
                   ,true                 ,false);
/*...................................................................*/

/*...*/
    tmp = sqrt(dot(sistEqEnergy->b,sistEqEnergy->b
                  ,sistEqEnergy->neqNov));
    fEnergy = true;
    if (tmp == 0.e0) fEnergy = false;  
/*...................................................................*/

/*...Ae=bE*/
    if (fEnergy) {
      tm.solvEnergy = getTimeC() - tm.solvEnergy;
      solverC(m
        ,sistEqEnergy->neq    ,sistEqEnergy->neqNov
        ,sistEqEnergy->nad    ,sistEqEnergy->nadr
        ,sistEqEnergy->ia     ,sistEqEnergy->ja
        ,sistEqEnergy->al     ,sistEqEnergy->ad, sistEqEnergy->au
        ,sistEqEnergy->b      ,sistEqEnergy->x
        ,&sistEqEnergy->iNeq  ,&sistEqEnergy->omp
        ,solvEnergy->tol      ,solvEnergy->maxIt
        ,sistEqEnergy->storage,solvEnergy->solver
        ,solvEnergy->fileSolv ,solvEnergy->log
        ,true                 ,sistEqEnergy->unsym);
      tm.solvEnergy = getTimeC() - tm.solvEnergy;
    }
/*...................................................................*/

/*... x -> Energy*/
    updateCellValue(mesh->elm.energy,sistEqEnergy->x
                   ,sistEqEnergy->id,&sistEqEnergy->iNeq
                   ,mesh->numel     ,1
                   ,eModel.fRes     ,true);
/*...................................................................*/

/*... residual*/
    residualSimpleLm(mesh->elm.vel     ,mesh->elm.energy
                 ,mesh->elm.rCellVel   ,rCellPc
                 ,mesh->elm.rCellEnergy
                 ,sistEqVel->ad        ,sistEqEnergy->ad
                 ,rU                   ,&rMass
                 ,&rEnergy          
                 ,sistEqVel->id        ,sistEqEnergy->id
                 ,mesh->numelNov       ,sistEqVel->neqNov
                 ,mesh->ndm            ,RSCALEDSUM);  
/*...................................................................*/

/*...*/
    if(fDensity)
      updateDensity(mesh->elm.energy,mesh->elm.densityFluid
                   ,mesh->numel      ,PROP_UPDATE_SIMPLE_LOOP);
    if(fSheat)
      updateSpecificHeat(mesh->elm.energy,mesh->elm.specificHeat
                       ,mesh->numel     ,PROP_UPDATE_SIMPLE_LOOP);
    if(fDvisc)
     updateDynamicViscosity(mesh->elm.energy,mesh->elm.dViscosity
                           ,mesh->numel);
    if(fTcond)
      updateDynamicViscosity(mesh->elm.energy,mesh->elm.tConductivity
                            ,mesh->numel);
/*...................................................................*/

/*...*/
    if (itSimple == kZeroPres && relRes &&  pCor) rMass0 = rMass;
    if (itSimple == kZeroEnergy && relRes && fEnergy) rEnergy0 = rEnergy;
    if (itSimple == kZeroVel && relRes) {
      if (xMomentum) rU0[0] = rU[0];
      if (yMomentum) rU0[1] = rU[1];
      if (zMomentum && ndfVel == 3) rU0[2] = rU[2];
    }
    conv = 0;
/*...................................................................*/

/*...*/
    timei = getTimeC() - time;
/*... arquivo de log*/
    if (opt.fItPlot) {
      if (ndfVel == 3)
        fprintf(opt.fileItPlot[FITPLOTSIMPLE]
          , "%d %20.8e %20.8e %20.8e %20.8e %20.8e\n"
          , itSimple + 1, rU[0], rU[1], rU[2], rMass, rEnergy);
      else
        fprintf(opt.fileItPlot[FITPLOTSIMPLE]
          , "%d %20.8e %20.8e %20.8e %20.8e\n"
          , itSimple + 1, rU[0], rU[1], rMass, rEnergy);
    }
/*...................................................................*/

/*...*/
    if (jj == sp->pSimple) {
      jj = 0;
      printf("It simple: %d \n", itSimple + 1);
      printf("Time(s)  : %lf \n", timei);
      printf("Residuo:\n");
      printf("conservacao da massa   : %20.8e\n", rMass / rMass0);
      printf("momentum x1            : %20.8e\n", rU[0] / rU0[0]);
      printf("momentum x2            : %20.8e\n", rU[1] / rU0[1]);
      if (ndfVel == 3)
        printf("momentum x3           : %20.8e\n", rU[2] / rU0[2]);
      printf("conservacao da energia : %20.8e\n",rEnergy/rEnergy0);
    }
   
/*...................................................................*/

/*... 3D*/
    if (ndfVel == 3) {
      if(rEnergy/rEnergy0<tolSimpleEnergy || rEnergy<SZERO) conv++;
/*...*/
      if(rMass/rMass0<tolSimpleMass || rMass<SZERO) conv++;
/*...*/
      if(rU[0]/rU0[0]<tolSimpleU1 || rU[0]<tmp*SZERO) conv++;
/*...*/
      if(rU[1]/rU0[1]<tolSimpleU2 || rU[1]<tmp*SZERO) conv++;
/*...*/
      if(rU[2]/rU0[2]<tolSimpleU3 || rU[2]<tmp*SZERO) conv++;
/*..*/
      if(conv == 5) break;
    }
/*...................................................................*/

/*... 2D*/
    else {
      if(rEnergy/rEnergy0<tolSimpleEnergy || rEnergy<SZERO) conv++;
/*...*/
      if(rMass/rMass0<tolSimpleMass || rMass<SZERO) conv++;
/*..*/
      if(rU[0]/rU0[0]<tolSimpleU1 || rU[0]<=tmp*SZERO) conv++;
/*...*/
      if(rU[1]/rU0[1]<tolSimpleU2 || rU[1]<=tmp*SZERO) conv++;
/*..*/
      if(conv == 4) break;
/*...................................................................*/
      jj++;
    }
/*...................................................................*/


  }
/*...................................................................*/
  timei = getTimeC() - time;

/*...*/
  fParameter[0] = true;
  fParameter[1] = true;
  fParameter[2] = true;
  parameterCellLm(mesh->elm.vel         ,mesh->elm.material.prop
                 ,mesh->elm.densityFluid,mesh->elm.specificHeat
                 ,mesh->elm.geom.volume ,mesh->elm.mat
                 ,&cfl                  ,&reynolds
                 ,&peclet
                 ,fParameter            ,sc.ddt.dt[0]
                 ,mesh->numelNov        ,mesh->ndm);
/*...................................................................*/

/*... guardando as propriedades para o proximo passo*/
  if(fDensity) updateDensity(mesh->elm.energy,mesh->elm.densityFluid
              ,mesh->numel     ,PROP_UPDATE_OLD_TIME   );
  if(fSheat) updateSpecificHeat(mesh->elm.energy,mesh->elm.specificHeat
                         ,mesh->numel   ,PROP_UPDATE_OLD_TIME);
/*...................................................................*/

/*...*/
  printf("It simple: %d \n", itSimple + 1);
  printf("Time(s)  : %lf \n", timei);
  printf("Reynolds: %lf\n", reynolds);
  printf("Pecelt  : %lf\n", peclet);
  if (sc.ddt.flag)
    printf("CFL     : %lf\n", cfl);
  printf("Residuo:\n");
  printf("conservacao da massa   (init,final): %20.8e %20.8e \n"
        , rMass0, rMass);
  printf("momentum x1            (init,final): %20.8e %20.8e\n"
        , rU0[0], rU[0]);
  printf("momentum x2            (init,final): %20.8e %20.8e\n"
        , rU0[1], rU[1]);
  if (ndfVel == 3)
    printf("momentum x3            (init,final): %20.8e %20.8e\n"
          , rU0[2], rU[2]);
  printf("conservacao da energia (init,final): %20.8e %20.8e\n"
        , rEnergy0, rEnergy);
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 02/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * UPDATECELLSIMPLEVEL : atualizacao dos valores das velocidades     *
 * estimadas com os valores das respectivas equacoes                 *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * w       -> variavel nas celulas                                   * 
 * u1      -> solucao do sistema                                     * 
 * u2      -> solucao do sistema                                     * 
 * id      -> numeracao das equacoes                                 * 
 * numel   -> numero de elementos                                    * 
 * ndm     -> numero de dimensoes                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * w      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void updateCellSimpleVel(DOUBLE  *RESTRICT w
                        ,DOUBLE  *RESTRICT u1 ,DOUBLE  *RESTRICT u2
                        ,INT  *RESTRICT id    ,INT const nEl
                        ,short const ndm)
{
  INT i,lNeq;

  for(i=0;i<nEl;i++){
    lNeq             = id[i] - 1;
    if(lNeq > -1){  
      MAT2D(i,0,w,2) = u1[lNeq];
      MAT2D(i,1,w,2) = u2[lNeq];
    }
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 30/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * UPDATECELLSIMPLEVELR: atualizacao dos valores das velocidades     *
 * estimadas com os valores das respectivas equacoes                 *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * w       -> variavel nas celulas                                   * 
 * u1      -> solucao do sistema                                     * 
 * u2      -> solucao do sistema                                     *
 * u2      -> solucao do sistema                                     *  
 * id      -> numeracao das equacoes                                 * 
 * numel   -> numero de elementos                                    * 
 * ndm     -> numero de dimensoes                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * w      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS: Versao residual  v*(n) = v(n) + dv(n+1)                      * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void updateCellSimpleVelR(DOUBLE  *RESTRICT w  ,DOUBLE  *RESTRICT u1
                         ,DOUBLE  *RESTRICT u2 ,DOUBLE  *RESTRICT u3
                         ,INT  *RESTRICT id    ,INT const nEl
                         ,short const ndm)
{
  INT i,lNeq;

  if(ndm == 3)
    for(i=0;i<nEl;i++){
      lNeq             = id[i] - 1;
      if(lNeq > -1){  
        MAT2D(i,0,w,2) += u1[lNeq];
        MAT2D(i,1,w,2) += u2[lNeq];
        MAT2D(i,1,w,3) += u3[lNeq];
      }
    }
  else
    for(i=0;i<nEl;i++){
      lNeq             = id[i] - 1;
      if(lNeq > -1){  
        MAT2D(i,0,w,2) += u1[lNeq];
        MAT2D(i,1,w,2) += u2[lNeq];
      }
    }
}
/*********************************************************************/ 



/********************************************************************* 
 * Data de criacao    : 11/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * UPDATECELLSIMPLEVEL3D : atualizacao dos valores das velocidades   *
 * estimadas com os valores das respectivas equacoes                 *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * w       -> variavel nas celulas                                   * 
 * u1      -> solucao do sistema                                     * 
 * u2      -> solucao do sistema                                     * 
 * u3      -> solucao do sistema                                     * 
 * id      -> numeracao das equacoes                                 * 
 * numel   -> numero de elementos                                    * 
 * ndm     -> numero de dimensoes                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * w      -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void updateCellSimpleVel3D(DOUBLE  *RESTRICT w
                          ,DOUBLE  *RESTRICT u1 
                          ,DOUBLE  *RESTRICT u2
                          ,DOUBLE  *RESTRICT u3
                          ,INT  *RESTRICT id    ,INT const nEl
                          ,short const ndm)
{
  INT i,lNeq;

  for(i=0;i<nEl;i++){
    lNeq             = id[i] - 1;
    if(lNeq > -1){  
      MAT2D(i,0,w,3) = u1[lNeq];
      MAT2D(i,1,w,3) = u2[lNeq];
      MAT2D(i,2,w,3) = u3[lNeq];
    }
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 02/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * UPDATECELLSIMPLEPRES: atualizacao dos valores das pressoes de     *
 * correcao com os valores das respectivas equacoes                  *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * presC   -> variavel nas celulas                                   * 
 * xp      -> solucao do sistema                                     * 
 * id      -> numeracao das equacoes                                 * 
 * numel   -> numero de elementos                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * presC  -> atualizado                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void updateCellSimplePres(DOUBLE  *RESTRICT presC,DOUBLE  *RESTRICT xp   
                          ,INT  *RESTRICT id      ,INT const nEl)
{
  INT i,lNeq;

  for(i=0;i<nEl;i++){
    lNeq     = id[i] - 1;
    if(lNeq > -1)  
      presC[i] = xp[lNeq];
/*... celulas prescritas*/
    else
      presC[i] = 0.e0;
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 02/07/2016                                   *
 * Data de modificaco : 15/08/2016                                   * 
 *-------------------------------------------------------------------* 
 * SIMPLEUPDATE : atualizacao das variasveis do metodo simple        *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * w         -> velocidades estimadas                                * 
 * pressure  -> pressao desatualizado                                * 
 * presC     -> pressao de correcao                                  * 
 * gradPresC -> gradiente da pressao de correcao                     * 
 * dField    -> matriz D do metodo simple                            * 
 * numel     -> numero de elementos                                  * 
 * ndm       -> numero de dimensoes                                  * 
 * alphaPres -> under-relaxationmensoes                              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * w        -> atualizado                                            * 
 * pressure -> atualizado                                            * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void simpleUpdate(DOUBLE *RESTRICT w     ,DOUBLE *RESTRICT pressure
                 ,DOUBLE *RESTRICT presC ,DOUBLE *RESTRICT gradPresC
                 ,DOUBLE *RESTRICT dField         
                 ,INT const nEl          ,short const ndm
                 ,DOUBLE const alphaPres){
  INT i; 

/*...*/  
  if( ndm == 2){
    for(i=0;i<nEl;i++){
/*... atualizacoes da velocidades*/
      MAT2D(i,0,w,2) -= 
      MAT2D(i,0,dField,2)*MAT2D(i,0,gradPresC,2);
      
      MAT2D(i,1,w,2) -=
      MAT2D(i,1,dField,2)*MAT2D(i,1,gradPresC,2);
/*...................................................................*/

/*... atualizacoes da velocidades*/
      pressure[i] += alphaPres*presC[i];
/*...................................................................*/
    }
  }
/*...................................................................*/

/*...*/  
  else if( ndm == 3){
    for(i=0;i<nEl;i++){
/*... atualizacoes da velocidades*/
      MAT2D(i,0,w,3) -=
      MAT2D(i,0,dField,3)*MAT2D(i,0,gradPresC,3);

      MAT2D(i,1,w,3) -=
      MAT2D(i,1,dField,3)*MAT2D(i,1,gradPresC,3);

      MAT2D(i,2,w,3) -=
      MAT2D(i,2,dField,3)*MAT2D(i,2,gradPresC,3);
/*...................................................................*/

/*... atualizacoes da velocidades*/
      pressure[i] += alphaPres*presC[i];
/*...................................................................*/
    }
  }

}               
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 09/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * SETSIMPLESCHEME : set o metodo simple incompressivel              *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * word -> str com a discretizacao                                   * 
 * sp   -> estrutura simple                                          * 
 * fieIn-> arquivo de entrada                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void setSimpleScheme(char *word,Simple *sp,FILE *fileIn){

  if(!strcmp(word,"SIMPLE"))
   sp->type = SIMPLE;
  else if(!strcmp(word,"SIMPLEC"))
   sp->type =  SIMPLEC;

  fscanf(fileIn,"%d",&sp->maxIt); 
  
  fscanf(fileIn,"%lf",&sp->alphaPres); 
  fscanf(fileIn,"%lf",&sp->alphaVel); 
  fscanf(fileIn,"%lf",&sp->tolPres); 
  fscanf(fileIn,"%lf",&sp->tolVel); 
  
  fscanf(fileIn,"%d",&sp->nNonOrth); 

  fscanf(fileIn, "%d", &sp->pSimple);

}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 24/08/2017                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* SETSIMPLELMSCHEME : set o metodo simple para baixo numero mach    *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* word -> str com a discretizacao                                   *
* sp   -> estrutura simple                                          *
* fieIn-> arquivo de entrada                                        *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
void setSimpleLmScheme(char *word, Simple *sp, FILE *fileIn) {

  if (!strcmp(word, "SIMPLE"))
    sp->type = SIMPLE;
  else if (!strcmp(word, "SIMPLEC"))
    sp->type = SIMPLEC;

  fscanf(fileIn, "%d", &sp->maxIt);

  fscanf(fileIn, "%lf", &sp->alphaPres);
  fscanf(fileIn, "%lf", &sp->alphaVel);
  fscanf(fileIn, "%lf", &sp->alphaEnergy);
  fscanf(fileIn, "%lf", &sp->tolPres);
  fscanf(fileIn, "%lf", &sp->tolVel);
  fscanf(fileIn, "%lf", &sp->tolEnergy);

  fscanf(fileIn, "%d", &sp->nNonOrth);

  fscanf(fileIn, "%d", &sp->pSimple);

}
/*********************************************************************/


/********************************************************************* 
 * Data de criacao    : 18/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * RESIDUALSIMPLE : calculo dos residuos no metodo simple            *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * vel      -> campo de velocidade                                   * 
 * rCellVel -> residuo das equacoes das velocidade por celulas       * 
 * rCellMass-> residuo de massa da equacoes de pressao por celulas   * 
 * adVel    -> diagonal da equacoes de velocidades                   * 
 * rU       -> nao definido                                          * 
 * rMass    -> nao definido                                          * 
 * idVel    -> numero da equacao da celula i                         * 
 * nEl      -> numero de elementos                                   * 
 * ndm      -> numero de dimensoes                                   * 
 * iCod     -> tipo de residuo                                       * 
 *          RSCALED - residuo com escala de grandeza                 * 
 *          RSQRT   - norma p-2 ( norma euclidiana)                  * 
 *          RSCALEDM- residuo com escala de grandeza                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * rU       -> residuo das velocidades                               * 
 * rMass    -> residuo de mass                                       * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 * rCellVel = | rvx1 rvx2 ... rvxn |                                 *
 *            | rvy1 rvy2 ... rvyn |                                 * 
 *            | rvz1 rvz2 ... rvzn |                                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void residualSimple(DOUBLE *RESTRICT vel
                 ,DOUBLE *RESTRICT rCellVel,DOUBLE *RESTRICT rCellMass
                 ,DOUBLE *RESTRICT adVel
                 ,DOUBLE *RESTRICT rU      ,DOUBLE *rMass
                 ,INT  *RESTRICT idVel     
                 ,INT const nEl            ,INT const nEqVel
                 ,short const ndm          ,short iCod)
{

  DOUBLE maxV[3],sum[3],mod,tmp,v,rScale;
  DOUBLE *p;
  INT i,j,lNeq;
  
/*...*/
  maxV[0] = maxV[1] = maxV[2] = 0.e0; 
  sum [0] = sum[1]  = sum[2] = 0.e0; 
  for(j=0;j<ndm;j++){
    rU[j]   = 0.e0;
  }
/*...................................................................*/

/*...*/
  switch(iCod){

/*... scaled*/
    case RSCALED:
/*... max(|Ap*velP|) */
      for(i=0;i<nEl;i++){
        lNeq = idVel[i] - 1;
        if(lNeq > -1){ 
          for(j=0;j<ndm;j++){
            v       = MAT2D(i,j,vel,ndm);
            lNeq   += j*nEqVel;
            mod     = fabs(adVel[lNeq]*v);
            maxV[j] = max(maxV[j],mod);
          }
        }
      }
/*...................................................................*/
      
/*... max ( | F - Ax |P / max( |Ap*velP| ) )*/
      for(j=0;j<ndm;j++){
        for(i=0;i<nEl;i++){
          mod    = fabs(MAT2D(j,i,rCellVel,nEl));
          if(maxV[j] != 0.e0)
            rScale = mod/maxV[j];
          else
            rScale = 0.e0;
          rU[j]  = max(rU[j],rScale);
        }
      }  
/*...................................................................*/

/*...*/
      tmp = 0.e0;
      for(i=0;i<nEl;i++){
        v    = fabs(rCellMass[i]);
        tmp += v;
      } 
      *rMass = tmp; 
/*...................................................................*/
    break;
/*...................................................................*/

/*... norma euclidiana*/
    case RSQRT:
/*...*/
      for(j=0;j<ndm;j++){
        p     = &rCellVel[j*nEl]; 
        rU[j] = sqrt(dot(p,p,nEl));
      }
/*...................................................................*/

/*...*/
      *rMass = sqrt(dot(rCellMass,rCellMass,nEl));
/*...................................................................*/
    break;
/*...................................................................*/

/*... scaled*/
    case RSCALEDSUM:
/*... sum(|Ap*velP|) */
      for(i=0;i<nEl;i++){
        lNeq = idVel[i] - 1;
        if(lNeq > -1){ 
          for(j=0;j<ndm;j++){
            v       = MAT2D(i,j,vel,ndm);
            lNeq   += j*nEqVel;
            mod     = fabs(adVel[lNeq]*v);
            sum[j] += mod;
          }
        }
      }
/*...................................................................*/
      
/*... sum ( | F - Ax |P / sum( |Ap*velP |) )*/
      for(j=0;j<ndm;j++){
        for(i=0;i<nEl;i++){
          mod    = fabs(MAT2D(j,i,rCellVel,nEl));
          rU[j] += mod;
        }
        if( sum[j] > rU[j]*SZERO)
          rU[j]  /=  sum[j];
      }  
/*...................................................................*/

/*...*/
      tmp = 0.e0;
      for(i=0;i<nEl;i++){
        v    = fabs(rCellMass[i]);
        tmp += v;
      } 
      *rMass = tmp; 
/*...................................................................*/
    break;
/*...................................................................*/

/*... */
     default:
       ERRO_OP(__FILE__,__func__,iCod);
     break;
/*...................................................................*/
  }

}

/********************************************************************* 
 * Data de criacao    : 26/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * RESIDUALSIMPLE : calculo dos residuos no metodo simple            *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * vel      -  > campo de velocidade                                 *
 * energy     -> campo de energia                                    *  
 * rCellVel   -> residuo das equacoes das velocidade por celulas     * 
 * rCellMass  -> residuo de massa da equacoes de pressao por celulas * 
 * rCellEnergy-> residuo de massa da equacoes de energia por celulas * 
 * adVel    -> diagonal da equacoes de velocidades                   *
 * adEnergy -> diagonal da equacoes de energia                       *  
 * rU       -> nao definido                                          * 
 * rMass    -> nao definido                                          * 
 * rEnergy  -> nao definido                                          * 
 * idVel    -> numero da equacao da celula i                         * 
 * idEnergy -> numero da equacao da celula i                         * 
 * nEl      -> numero de elementos                                   * 
 * nEqVel   -> numero de equacoes de velocidade                      * 
 * ndm      -> numero de dimensoes                                   * 
 * iCod     -> tipo de residuo                                       * 
 *          RSCALED - residuo com escala de grandeza                 * 
 *          RSQRT   - norma p-2 ( norma euclidiana)                  * 
 *          RSCALEDM- residuo com escala de grandeza                 * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * rU       -> residuo das velocidades                               * 
 * rMass    -> residuo de mass                                       * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *            | rvx1 rvx2 ... rvxnEl |                               *
 * rCellVel = | rvy1 rvy2 ... rvynEl |                               * 
 *            | rvz1 rvz2 ... rvznEl |                               * 
 *                                                                   *
 *         | adx1 adx2 ... adxnEq |                                  *
 * adVel = | ady1 ady2 ... adynEq |                                  *
 *         | adz1 adz2 ... adznEq |                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void residualSimpleLm(DOUBLE *RESTRICT vel ,DOUBLE *RESTRICT energy
            ,DOUBLE *RESTRICT rCellVel   ,DOUBLE *RESTRICT rCellMass
            ,DOUBLE *RESTRICT rCellEnergy
            ,DOUBLE *RESTRICT adVel       ,DOUBLE *RESTRICT adEnergy 
            ,DOUBLE *RESTRICT rU          ,DOUBLE *rMass
            ,DOUBLE *rEnergy 
            ,INT  *RESTRICT idVel         ,INT  *RESTRICT idEnergy 
            ,INT const nEl                ,INT const nEqVel
            ,short const ndm              ,short iCod)

{
  DOUBLE maxV[3],sum[3],maxE,mod,tmp,v,rScale;
  DOUBLE *p;
  INT i,j,lNeq;
  
/*...*/
  maxV[0] = maxV[1] = maxV[2] = maxE = 0.e0; 
  sum [0] = sum[1]  = sum[2] = 0.e0; 
  for(j=0;j<ndm;j++){
    rU[j]   = 0.e0;
  }
  *rEnergy = *rMass = 0.e0;
/*...................................................................*/

/*...*/
  switch(iCod){

/*... scaled*/
    case RSCALED:
/*... max(Ap*velP) */
      for(i=0;i<nEl;i++){
        lNeq = idVel[i] - 1;
        if(lNeq > -1){ 
          for(j=0;j<ndm;j++){
            v       = MAT2D(i,j,vel,ndm);
            lNeq   += j*nEqVel;  
            mod     = fabs(adVel[lNeq]*v);
            maxV[j] = max(maxV[j],mod);
          }
        }
      }
/*...................................................................*/
      
/*... max ( | F - Ax |P / max(Ap*velP) )*/
      for(j=0;j<ndm;j++){
        for(i=0;i<nEl;i++){
          mod    = fabs(MAT2D(j,i,rCellVel,nEl));
          if(maxV[j] != 0.e0)
            rScale = mod/maxV[j];
          else
            rScale = 0.e0;
          rU[j]  = max(rU[j],rScale);
        }
      }  
/*...................................................................*/

/*... max(Ap*energyP) */
      for(i=0;i<nEl;i++){
        lNeq = idEnergy[i] - 1;
        if(lNeq > -1){ 
          v       = energy[i];
          mod     = fabs(adEnergy[lNeq]*v);
          maxE    = max(maxE,mod);
        }
      }
/*...................................................................*/
      
/*... max ( | F - Ax |P /max(Ap*energyP) )*/
      for(i=0;i<nEl;i++){
        mod    = fabs(MAT2D(j,i,rCellVel,nEl));
        if(maxE != 0.e0)
          rScale = mod/maxE;
        else
          rScale = 0.e0;
        *rEnergy = max(*rEnergy,rScale);
      }  
/*...................................................................*/

/*...*/
      tmp = 0.e0;
      for(i=0;i<nEl;i++){
        v    = fabs(rCellMass[i]);
        tmp += v;
      } 
      *rMass = tmp; 
/*...................................................................*/
    break;
/*...................................................................*/

/*... norma euclidiana*/
    case RSQRT:
/*...*/
      for(j=0;j<ndm;j++){
        p     = &rCellVel[j*nEl]; 
        rU[j] = sqrt(dot(p,p,nEl));
      }
/*...................................................................*/

/*...*/
      *rEnergy = sqrt(dot(rCellEnergy,rCellEnergy,nEl));
/*...................................................................*/

/*...*/
      *rMass = sqrt(dot(rCellMass,rCellMass,nEl));
/*...................................................................*/
    break;
/*...................................................................*/

/*... scaled*/
    case RSCALEDSUM:
/*... sum( |Ap*velP |) */
      for(i=0;i<nEl;i++){
        lNeq = idVel[i] - 1;
        if(lNeq > -1){ 
          for(j=0;j<ndm;j++){
            v       = MAT2D(i,j,vel,ndm);
            lNeq   += j*nEqVel;
            mod     = fabs(adVel[lNeq]*v);
            sum[j] += mod;
          }
        }
      }
/*...................................................................*/
      
/*... sum ( | F - Ax |P / sum( |Ap*velP| ) )*/
      for(j=0;j<ndm;j++){
        for(i=0;i<nEl;i++){
          mod    = fabs(MAT2D(j,i,rCellVel,nEl));
          rU[j] += mod;
        }
        if( sum[j] > rU[j]*SZERO)
          rU[j]  /=  sum[j];
      }  
/*...................................................................*/

/*... sum( |Ap*energyP| ) */
      for(i=0;i<nEl;i++){
        lNeq = idEnergy[i] - 1;
        if(lNeq > -1){ 
          v       = energy[i];
          mod     = fabs(adEnergy[lNeq]*v);
          sum[0] += mod;
        }
      }
/*...................................................................*/
      
/*... sum ( | F - Ax |P / sum( |Ap*energyP| ) )*/
      for(i=0;i<nEl;i++){
         mod      = fabs(rCellEnergy[i]);
        *rEnergy += mod;
      }
      if( sum[0] > (*rEnergy)*SZERO)
        *rEnergy /=  sum[0];  
/*...................................................................*/


/*...*/
      tmp = 0.e0;
      for(i=0;i<nEl;i++){
        v    = fabs(rCellMass[i]);
        tmp += v;
      } 
      *rMass = tmp; 
/*...................................................................*/
    break;
/*...................................................................*/

/*... */
     default:
       ERRO_OP(__FILE__,__func__,iCod);
     break;
/*...................................................................*/
  }

}