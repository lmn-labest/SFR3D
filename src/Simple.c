#include<Simple.h>

void simpleSolver(Memoria *m        
                 ,Loads *loadsVel   ,Loads *loadsPres 
                 ,Mesh *mesh0       ,Mesh *mesh       
                 ,SistEq *sistEqVel ,SistEq *sistEqPres
                 ,Solv *solvVel     ,Solv *solvPres 
                 ,Simple *sp
                 ,Scheme sc         ,PartMesh *pMesh 
                 ,FileOpt opt       ,char *preName  
                 ,char *nameOut     ,FILE *fileOut){

  #define SZERO 1.0e-32
  short unsigned ndfVel = mesh->ndfF-1;
  short unsigned itSimple,conv;
  short unsigned kZeroVel  = sp->kZeroVel;
  short unsigned kZeroPres = sp->kZeroPres;
  INT nEl = mesh->numel,i;
  INT jj = 1;
  DOUBLE time;
  DOUBLE *b1,*b2,*bPc,*xu1,*xu2,*xp;
  DOUBLE *rCellU1,*rCellU2,*rCellPc;
/*...*/
  DOUBLE rU1,rU2,rU10,rU20,rPc,rPc0;
/*...*/
  DOUBLE tolSimpleU1,tolSimpleU2,tolSimpleMass;
/*...*/
  bool xMomentum,yMomentum,pCor,fPrint=false;

  time = getTimeC();
/*...*/
  b1       =  sistEqVel->b; 
  b2       = &sistEqVel->b[sistEqVel->neq]; 
  bPc      = sistEqPres->b; 

  xu1      = sistEqVel->x;
  xu2      = &sistEqVel->x[sistEqVel->neq];
  xp       = sistEqVel->x;

  rCellU1  = mesh->elm.rCellVel;
  rCellU2  = &mesh->elm.rCellVel[nEl];
  rCellPc  = mesh->elm.rCellPres;
/*...................................................................*/

/*...*/
  tolSimpleU1   = tolSimpleU2 = sp->tolVel;
  tolSimpleMass = sp->tolPres;
/*...................................................................*/

/*...*/
  rU10 = rU20 = rPc0 = 1.e0;
  rU1  = rU2  = rPc  = 0.e0;
  conv = 0;
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
                 ,mesh->numelNov       ,mesh->ndfF
                 ,mesh->ndm            ,mesh->maxViz);
  tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
/*...................................................................*/

/*... discretizacao temporal*/
  if(sc.ddt.flag){
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
    cellTransientSimple(mesh->elm.geom.volume   ,sistEqVel->id     
                       ,mesh->elm.vel0          ,mesh->elm.vel 
                       ,mesh->elm.densityFluid  ,sistEqVel->b0
                       ,sc.ddt.dt               ,sistEqVel->neqNov
                       ,mesh->numelNov          ,ndfVel     
                       ,sc.ddt.type             ,true);
/*... vel(n-1) = vel(n)*/
    alphaProdVector(1.e0,mesh->elm.vel
                   ,mesh->numel*ndfVel ,mesh->elm.vel0); 
/*...................................................................*/
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
  }
/*...................................................................*/

/*...*/
  for(itSimple=0;itSimple<sp->maxIt;itSimple++){
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
                    ,sc.advVel               ,sp->type
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
     rU1 = sqrt(dot(b1,b1,sistEqVel->neqNov));
     xMomentum = true;
     if( rU1 < SZERO ) xMomentum = false;
     if( itSimple == kZeroVel && xMomentum){
       rU1  = sqrt(dot(rCellU1,rCellU1,sistEqVel->neqNov));
       rU10 = rU1;
     }
/*...................................................................*/

/*...*/
     rU2 = sqrt(dot(b2,b2,sistEqVel->neqNov));
     yMomentum = true;
     if( rU2 < SZERO ) yMomentum = false;
     if( itSimple == kZeroVel && yMomentum){
       rU1  = sqrt(dot(rCellU2,rCellU2,sistEqVel->neqNov));
       rU20 = rU2;
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
           ,sistEqVel->al      ,sistEqVel->ad,sistEqVel->au
           ,b1                 ,xu1          
           ,&sistEqVel->iNeq
           ,solvVel->tol       ,solvVel->maxIt     
           ,sistEqVel->storage ,solvVel->solver
           ,solvVel->fileSolv  ,solvVel->log  
           ,false              ,false
           ,sistEqVel->unsym   ,false);   
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
           ,sistEqVel->al      ,sistEqVel->ad,sistEqVel->au
           ,b2                 ,xu2
           ,&sistEqVel->iNeq
           ,solvVel->tol       ,solvVel->maxIt     
           ,sistEqVel->storage ,solvVel->solver
           ,solvVel->fileSolv  ,solvVel->log  
           ,false              ,false
           ,sistEqVel->unsym   ,false);  
       tm.solvVel = getTimeC() - tm.solvVel;
     }
/*...................................................................*/

/*... atualizando o campo de velociade estimadas*/
     updateCellSimpleVel(mesh->elm.vel,xu1,xu2,sistEqVel->id
                       ,mesh->numelNov,mesh->ndm); 
/*...................................................................*/

/*...*/
     if(fPrint) printf("Correcao de pressao:\n");

/*... montagem do sistema  da pressao de correca*/
     tm.systFormPres = getTimeC() - tm.systFormPres;
     systFormSimplePres(loadsVel              ,loadsPres
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
      
/*...*/
     rPc = sqrt(dot(rCellPc,rCellPc,sistEqPres->neqNov));
     pCor = true;
     if( rPc < SZERO ) pCor = false;
     if( itSimple == kZeroPres) rPc0 = rPc;
     rU1  = sqrt(dot(rCellU1,rCellU1,sistEqVel->neqNov));
     rU2  = sqrt(dot(rCellU2,rCellU2,sistEqVel->neqNov));
     conv = 0;
/*...................................................................*/

/*... solver ApPc = bpC (velocidade estimadas)*/
     if(pCor){
       tm.solvPres = getTimeC() - tm.solvPres;
       solverC(m               
           ,sistEqPres->neq    ,sistEqPres->neqNov  
           ,sistEqPres->nad    ,sistEqPres->nadr
           ,sistEqPres->ia     ,sistEqPres->ja  
           ,sistEqPres->al     ,sistEqPres->ad,sistEqPres->au
           ,bPc                ,xp
           ,&sistEqPres->iNeq
           ,solvPres->tol      ,solvPres->maxIt     
           ,sistEqPres->storage,solvPres->solver
           ,solvPres->fileSolv ,solvPres->log  
           ,false              ,false
           ,sistEqPres->unsym  ,false);   
       tm.solvPres = getTimeC() - tm.solvPres;
/*...................................................................*/

/*... atualizando da pressao de correcao*/
       updateCellSimplePres(sp->ePresC,xp,sistEqPres->id,mesh->numelNov); 
/*...................................................................*/

/*... reconstruindo do gradiente da pressao correcao*/
       tm.rcGradPres = getTimeC() - tm.rcGradPres;
       rcGradU(m                    ,loadsPres
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
     }
/*...................................................................*/

/*... atualizacao de u, v, w e p*/
     simpleUpdate(mesh->elm.vel,mesh->elm.pressure
                ,sp->ePresC,sp->eGradPresC
                ,sp->d       
                ,mesh->numel,mesh->ndm
                ,sp->alphaPres);
/*...................................................................*/

/*...*/
     if( rPc/rPc0 < tolSimpleMass || rPc < SZERO) conv++;
/*..*/
     if( rU1/rU10 < tolSimpleU1   || rU1 < SZERO) conv++;
/*..*/
     if( rU2/rU20 < tolSimpleU2   || rU2 < SZERO) conv++;
/*..*/
     if( conv == 3) break;
/*...................................................................*/

/*...*/
     if( jj == 50) { 
       jj = 0; 
       printf("It simple: %d \n",itSimple+1);
       printf("Residuo:\n");
       printf("conservacao da massa: %20.8e\n",rPc/rPc0);
       printf("momentum x1         : %20.8e\n",rU1/rU10);
       printf("momentum x2         : %20.8e\n",rU2/rU20);
     } 
     jj++; 
/*...................................................................*/
  }
/*...................................................................*/
  time = getTimeC() -time;

/*...*/
  printf("It simple: %d \n",itSimple+1);
  printf("Time(s)  : %lf \n",time);
  printf("Residuo:\n");
  printf("conservacao da massa: %20.8e\n",rPc);
  printf("momentum x1         : %20.8e\n",rU1);
  printf("momentum x2         : %20.8e\n",rU2);
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
void updateCellSimpleVel(DOUBLE  *restrict w
                        ,DOUBLE  *restrict u1 ,DOUBLE  *restrict u2
                        ,INT  *restrict id    ,INT const nEl
                        ,short const ndm)
{
  INT i,lNeq;

  for(i=0;i<nEl;i++){
    lNeq             = id[i] - 1;
    if(lNeq > -1){  
      MAT2D(i,0,w,ndm) = u1[lNeq];
      MAT2D(i,1,w,ndm) = u2[lNeq];
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
void updateCellSimplePres(DOUBLE  *restrict presC,DOUBLE  *restrict xp   
                          ,INT  *restrict id      ,INT const nEl)
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
 * Data de modificaco : 00/00/0000                                   * 
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
void simpleUpdate(DOUBLE *restrict w     ,DOUBLE *restrict pressure
                 ,DOUBLE *restrict presC ,DOUBLE *restrict gradPresC
                 ,DOUBLE *restrict dField         
                 ,INT const nEl          ,short const ndm
                 ,DOUBLE const alphaPres){
  INT i; 

  if( ndm == 2){
    for(i=0;i<nEl;i++){
/*... atualizacoes da velocidades*/
      MAT2D(i,0,w,ndm) -= dField[i]*MAT2D(i,0, gradPresC,ndm);
      MAT2D(i,1,w,ndm) -= dField[i]*MAT2D(i,1, gradPresC,ndm);
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
 * SETSIMPLESCHEME : set o metodo simple                             *
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

}
/*********************************************************************/




 
