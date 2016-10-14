#include<Prime.h>
/*********************************************************************
* Data de criacao    : 05/09/2016                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* PRIME: metodo prime ▓                                             *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
*********************************************************************/
void primeSolver(Memoria *m
                  ,Loads *loadsVel, Loads *loadsPres
                  ,Mesh *mesh0, Mesh *mesh
                  ,SistEq *sistEqPres
                  ,Solv *solvPres
                  ,Prime *pr
                  ,Scheme sc, PartMesh *pMesh
                  ,FileOpt opt, char *preName
                  ,char *nameOut, FILE *fileOut
                  ,short const ndf) {
  FILE *fStop = NULL;
  short unsigned ndfVel = mesh->ndfF - 1;
  short unsigned conv;
  int itPrime;
  short unsigned kZeroVel = pr->kZeroVel;
  short unsigned kZeroPres = pr->kZeroPres;
  INT jj = 1;
  DOUBLE time, timei;
  DOUBLE *bPc,*xp,*v1,*v2,*v3;
  DOUBLE *rCellPc;
/*...*/
  DOUBLE rU[3], rU0[3], tmp, tb[3], rMass0, rMass;
/*...*/
  DOUBLE tolPrimeU1,tolPrimeU2,tolPrimeU3,tolPrimeMass;
/*...*/
  bool xMomentum, yMomentum, pCor;
  bool relRes = false;
  bool fPrint = false;
  DOUBLE cfl, reynolds;
  bool fParameter[2];

  time = getTimeC();
/*...*/
  v1 = pr->velUp;
  v2 = &pr->velUp[mesh->numelNov];
  v3 = &pr->velUp[2*mesh->numelNov];

  bPc = sistEqPres->b;

  xp = sistEqPres->x;

  rCellPc = mesh->elm.rCellPres;
/*...................................................................*/

/*...*/
  tolPrimeU1 = tolPrimeU2 = tolPrimeU3 = pr->tolVel;
  tolPrimeMass = pr->tolPres;
/*...................................................................*/

/*...*/
  rMass0 = 1.e0;
  rMass  = 0.e0;
  rU[0]  = rU[1]  = rU[2]   = 0.e0;
  rU0[0] = rU0[1] = rU0[2] = 1.e0;
  conv = 0;
/*...................................................................*/

/*...*/
  zero(sistEqPres->b0, sistEqPres->neqNov, DOUBLEC);
/*...................................................................*/

/*... restricoes por centro de celula u0 e cargas por volume b0*/
//tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
//cellPloadSimple(loadsPres, mesh->elm.geom.cc
//                , mesh->elm.faceRpres, mesh->elm.faceLoadPres
//                , mesh->elm.geom.volume
//                , sistEqVel->id, sistEqPres->id
//                , mesh->elm.vel, mesh->elm.pressure
//                , sistEqVel->b0, sistEqPres->b0
//                , mesh->numelNov, mesh->ndfF
//                , mesh->ndm, mesh->maxViz);
//tm.cellPloadSimple = getTimeC() - tm.cellPloadSimple;
/*...................................................................*/

/*... discretizacao temporal*/
  if (sc.ddt.flag) {  
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
    cellTransientPrime(mesh->elm.geom.volume
                      ,mesh->elm.vel0        ,mesh->elm.vel
                      ,mesh->elm.densityFluid,pr->bTemporal
                      ,sc.ddt
                      ,mesh->numelNov        ,ndfVel
                      ,false);
/*... vel(n-1) = vel(n)*/
    alphaProdVector(1.e0, mesh->elm.vel
                    , mesh->numel*ndfVel, mesh->elm.vel0);
/*...................................................................*/
    tm.cellTransientSimple = getTimeC() - tm.cellTransientSimple;
  }  
/*...................................................................*/

/*...*/
  for (itPrime = 0; itPrime< pr->maxIt; itPrime++) {
/*...*/
    if ((fStop = fopen("stopPrime.mvf", "r")) != NULL) {
      fclose(fStop);
      break;
    }
/*...................................................................*/

/*... calculo da matrix jacobiana das velocidades
    | du1dx1 du1dx2 du1dx3 |
    | du2dx1 du2dx2 du2dx3 |
    | du3dx1 du3dx2 du3dx3 |
*/
    tm.rcGradVel = getTimeC() - tm.rcGradVel;
    rcGradU(m, loadsVel
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
            ,mesh->elm.faceRvel   ,mesh->elm.faceLoadVel
            ,mesh->elm.vel        ,mesh->elm.gradVel
            ,mesh->node.vel       ,sc.rcGrad
            ,mesh->maxNo          ,mesh->maxViz
            ,ndfVel               ,mesh->ndm
            ,&pMesh->iNo          ,&pMesh->iEl
            ,mesh->numelNov       ,mesh->numel
            ,mesh->nnodeNov       ,mesh->nnode);
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
            ,mesh->elm.faceRpres,mesh->elm.faceLoadPres
            ,mesh->elm.pressure ,mesh->elm.gradPres
            ,mesh->node.pressure,sc.rcGrad
            ,mesh->maxNo        ,mesh->maxViz
            ,1, mesh->ndm
            ,&pMesh->iNo        ,&pMesh->iEl
            ,mesh->numelNov     ,mesh->numel
            ,mesh->nnodeNov     ,mesh->nnode);
    tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*... calculo da velocidade explicita*/
   tm.velExp = getTimeC() - tm.velExp;
   velExp(loadsVel               ,loadsPres
         ,sc.advVel              ,sc.diffVel
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
         ,mesh->elm.faceRvel     ,mesh->elm.faceLoadVel
         ,mesh->elm.faceRpres    ,mesh->elm.faceLoadPres
         ,mesh->elm.pressure     ,mesh->elm.gradPres
         ,mesh->elm.vel          ,pr->velUp
         ,mesh->elm.gradVel      ,pr->bTemporal
         ,pr->d                  ,pr->alphaVel
         ,sc.ddt
         ,mesh->maxNo            ,mesh->maxViz
         ,mesh->ndm              ,mesh->numelNov
         ,ndfVel                 ,pr->sPressure
         ,false);  
     tm.velExp = getTimeC() - tm.velExp;
/*...................................................................*/

/*... calculo do residui R = b - A vExp*/
     tm.velExp = getTimeC() - tm.velExp;
     velResidual(loadsVel             ,loadsPres
              ,sc.advVel              ,sc.diffVel
              ,mesh->elm.node         ,mesh->elm.adj.nelcon
              ,mesh->elm.nen          , mesh->elm.adj.nViz
              ,mesh->elm.geomType     ,mesh->elm.material.prop
              ,mesh->elm.material.type,mesh->elm.mat
              ,mesh->elm.geom.cc      ,pr->aD
              ,mesh->elm.geom.ksi     ,mesh->elm.geom.mksi
              ,mesh->elm.geom.eta     ,mesh->elm.geom.fArea
              ,mesh->elm.geom.normal  ,mesh->elm.geom.volume
              ,mesh->elm.geom.xm      ,mesh->elm.geom.xmcc
              ,mesh->elm.geom.vSkew   ,mesh->elm.geom.mvSkew
              ,mesh->elm.geom.dcca    ,mesh->elm.densityFluid
              ,mesh->elm.faceRvel     ,mesh->elm.faceLoadVel
              ,mesh->elm.faceRpres    ,mesh->elm.faceLoadPres
              ,mesh->elm.pressure     ,mesh->elm.gradPres
              ,pr->velUp              ,mesh->elm.rCellVel
              ,mesh->elm.gradVel      ,pr->bTemporal
              ,pr->d                  ,pr->alphaVel
              ,sc.ddt
              ,mesh->maxNo            ,mesh->maxViz
              ,mesh->ndm              ,mesh->numelNov
              ,ndfVel                 ,pr->sPressure
              ,true);    
     tm.velExp = getTimeC() - tm.velExp;
/*...................................................................*/
 
/*...*/
     tb[0] = sqrt(dot(v1,v1,mesh->numelNov));
     tb[1] = sqrt(dot(v2,v2,mesh->numelNov));
     if(ndf == 3) tb[2] = sqrt(dot(v3,v3,mesh->numelNov));
//     printf("%lf %lf\n",tb[0],tb[1]);
     if (itPrime == 0){
       tmp = max(tb[0], tb[1]);
       if(ndf == 3) tmp = max(tmp, tb[2]);
     }
/*...................................................................*/

/*...*/
    if (fPrint) printf("Correcao de pressao:\n");

/*... montagem do sistema  da pressao de correca*/
    tm.systFormPres = getTimeC() - tm.systFormPres;
    systFormSimplePres(loadsVel    ,loadsPresC
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
           ,pr->velUp              ,pr->d
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
    residualPrime(mesh->elm.vel     ,pr->aD 
                 ,mesh->elm.rCellVel,rCellPc
                 ,rU                ,&rMass
                 ,mesh->numelNov    ,mesh->ndm
                 ,3);
/*...................................................................*/

/*...*/
    pCor = true;
    if (rMass < tmp*SZERO) pCor = false;
    if (itPrime == kZeroPres && relRes) rMass0 = rMass;
    if (itPrime == kZeroVel && relRes) {
      rU0[0] = rU[0];
      rU0[1] = rU[1];
      if (ndf == 3)  rU0[2] = rU[2];
    }
    conv = 0;
/*...................................................................*/

/*... solver ApPc = bpC (velocidade estimadas)*/
    if (pCor) {
      tm.solvPres = getTimeC() - tm.solvPres;
      solverC(m
              , sistEqPres->neq, sistEqPres->neqNov
              , sistEqPres->nad, sistEqPres->nadr
              , sistEqPres->ia, sistEqPres->ja
              , sistEqPres->al, sistEqPres->ad, sistEqPres->au
              , bPc, xp
              , &sistEqPres->iNeq, &sistEqPres->omp
              , solvPres->tol, solvPres->maxIt
              , sistEqPres->storage, solvPres->solver
              , solvPres->fileSolv, solvPres->log
              , false, sistEqPres->unsym);
      tm.solvPres = getTimeC() - tm.solvPres;
    }
/*...................................................................*/

/*... atualizando da pressao de correcao*/
    updateCellPrimePres(pr->ePresC, xp, sistEqPres->id
                       ,mesh->numelNov);
/*...................................................................*/

/*... reconstruindo do gradiente da pressao correcao*/
    tm.rcGradPres = getTimeC() - tm.rcGradPres;
    rcGradU(m                    ,loadsPresC
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
           ,pr->ePresC           ,pr->eGradPresC
           ,pr->nPresC           ,sc.rcGrad
           ,mesh->maxNo          ,mesh->maxViz
           ,1                    ,mesh->ndm
           ,&pMesh->iNo          ,&pMesh->iEl
           ,mesh->numelNov       ,mesh->numel
           ,mesh->nnodeNov       ,mesh->nnode);
    tm.rcGradPres = getTimeC() - tm.rcGradPres;
/*...................................................................*/

/*... atualizacao de u, v, w e p*/
    primeUpdate(mesh->elm.vel      ,pr->velUp
               ,mesh->elm.pressure ,pr->ePresC    
               ,pr->eGradPresC     ,pr->d
               ,mesh->numel        ,mesh->ndm
               ,pr->alphaPres);
/*...................................................................*/

/*...*/
    if( ndf == 2) {
      if (rMass / rMass0 < tolPrimeMass || rMass < tmp*SZERO) conv++;
/*..*/
      if (rU[0] / rU0[0] < tolPrimeU1 || rU[0] < tmp*SZERO) conv++;
/*...*/
      if (rU[1] / rU0[1] < tolPrimeU2 || rU[1] < tmp*SZERO) conv++;
/*..*/
      if (conv == 3) break;
    }
/*...................................................................*/

/*...*/
    else if (ndf == 3) {
      if (rMass / rMass0 < tolPrimeMass || rMass < tmp*SZERO) conv++;
/*...*/
      if (rU[0] / rU0[0] < tolPrimeU1 || rU[0] < tmp*SZERO) conv++;
/*...*/
      if (rU[1] / rU0[1] < tolPrimeU2 || rU[1] < tmp*SZERO) conv++;
/*...*/
      if (rU[2] / rU0[2] < tolPrimeU3 || rU[2] < tmp*SZERO) conv++;
/*..*/
      if (conv == 4) break;
    }
/*...................................................................*/

/*...*/
    timei = getTimeC() - time;
/*... arquivo de log*/
    if (opt.fItPlot){
      if(ndf == 2)
        fprintf(opt.fileItPlot[FITPLOTSIMPLE]
               , "%d %20.8e %20.8e %20.8e\n"
               , itPrime+1, rU[0], rU[1], rMass);
      else if (ndf == 3)
        fprintf(opt.fileItPlot[FITPLOTSIMPLE]
               , "%d %20.8e %20.8e %20.8e %20.8e\n"
               ,itPrime+1,rU[0],rU[1],rU[2],rMass);
    }
/*...................................................................*/

/*...*/
    timei = getTimeC() - time;
    if (jj == pr->pPrime) {
      jj = 0;
      printf("It simple: %d \n", itPrime + 1);
      printf("Time(s)  : %lf \n", timei);
      printf("Residuo:\n");
      printf("conservacao da massa: %20.8e\n", rMass/rMass0);
      printf("momentum x1         : %20.8e\n", rU[0]/rU0[0]);
      printf("momentum x2         : %20.8e\n", rU[1]/rU0[1]);
      if (ndf == 3) 
        printf("momentum x3         : %20.8e\n", rU[2]/rU0[2]);
    }
    jj++;
 /*...................................................................*/
  }
/*...................................................................*/
  time = getTimeC() - time;

/*...*/
  fParameter[0] = true;
  fParameter[1] = true;
  parameterCell(mesh->elm.vel         ,mesh->elm.material.prop
               ,mesh->elm.densityFluid,mesh->elm.geom.volume
               ,mesh->elm.mat
               ,&cfl                  ,&reynolds
               ,fParameter            ,sc.ddt.dt[0]
               ,mesh->numelNov        ,mesh->ndm);
/*...................................................................*/

/*...*/
  printf("It prime: %d \n", itPrime + 1);
  printf("Time(s)  : %lf \n", time);
  printf("Reynolds: %lf\n", reynolds);
  if(sc.ddt.flag)
    printf("CFL     : %lf\n", cfl);
  printf("Residuo:\n");

/*...*/
  if(relRes){
    printf("conservacao da massa (init,final): %20.8e %20.8e \n"
         , rMass0, rMass);
    printf("momentum x1          (init,final): %20.8e %20.8e\n"
         , rU0[0], rU[0]);
    printf("momentum x2          (init,final): %20.8e %20.8e\n"
         , rU0[1], rU[1]);
    if(ndf == 3)
      printf("momentum x3          (init,final): %20.8e %20.8e\n"
             , rU0[2], rU[2]);
  }
/*...................................................................*/

/*...*/
  else{
    printf("conservacao da massa (final): %20.8e \n",rMass);
    printf("momentum x1          (final): %20.8e\n",rU[0]);
    printf("momentum x2          (final): %20.8e\n",rU[1]);
    if(ndf == 3)
      printf("momentum x3          (final): %20.8e\n",rU[2]);
  }
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 05/09/2016                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * UPDATECELLPIMEPRES: atualizacao dos valores das pressoes de       *
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
void updateCellPrimePres(DOUBLE  *restrict presC, DOUBLE  *restrict xp
                          , INT  *restrict id, INT const nEl)
{
  INT i, lNeq;

  for (i = 0; i<nEl; i++) {
    lNeq = id[i] - 1;
    if (lNeq > -1)
      presC[i] = xp[lNeq];
    /*... celulas prescritas*/
    else
      presC[i] = 0.e0;
  }
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 07/09/2016                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * PRIMEUPDATE : atualizacao das variasveis do metodo simple         *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * w         -> velocidades estimadas                                *
 * wUp       -> velocidades obitidas explicitamente                  *
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
void primeUpdate(DOUBLE *restrict w        ,DOUBLE *restrict wUp
                ,DOUBLE *restrict pressure ,DOUBLE *restrict presC 
                ,DOUBLE *restrict gradPresC,DOUBLE *restrict dField
                ,INT const nEl             ,short const ndm
                ,DOUBLE const alphaPres) {
  INT i;

/*...*/
  if (ndm == 2) {
    for (i = 0; i<nEl; i++) {
/*... atualizacoes da velocidades*/
     MAT2D(i,0,w,2) = MAT2D(i,0,wUp,2) -
        MAT2D(i,0,dField,2)*MAT2D(i,0,gradPresC,2);

      MAT2D(i,1,w,2) = MAT2D(i,1,wUp,2) -
        MAT2D(i,1,dField,2)*MAT2D(i,1,gradPresC,2);
/*...................................................................*/

/*... atualizacoes da velocidades*/
      pressure[i] += alphaPres*presC[i];
/*...................................................................*/
    }
  }
/*...................................................................*/

/*...*/
  else if (ndm == 3) {
    for (i = 0; i<nEl; i++) {
/*... atualizacoes da velocidades*/
      MAT2D(i,0,w,3) = MAT2D(i,0,wUp,3) -
        MAT2D(i,0,dField,3)*MAT2D(i,0,gradPresC,3);

      MAT2D(i,1,w,3) = MAT2D(i,1,wUp,3) -
        MAT2D(i,1,dField,3)*MAT2D(i,1,gradPresC,3);

      MAT2D(i,2,w,3) = MAT2D(i,2,wUp,3) -
        MAT2D(i,2,dField,3)*MAT2D(i,2,gradPresC,3);
/*...................................................................*/

/*... atualizacoes da velocidades*/
      pressure[i] += alphaPres*presC[i];
/*...................................................................*/
    }
  }
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 07/09/2016                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* RESIDUALPRIME : calculo dos residuos no metodo prime              *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* vel      -> campo de velocidade                                   *
* advel    -> diagonal principal da equacao de velociades           *
* rCellVel -> residuo das equacoes das velocidade por celulas       *
* rCellMass-> residuo de massa da equacoes de pressao por celulas   *
* rU       -> nao definido                                          *
* rMass    -> nao definido                                          *
* nEl      -> numero de elementos                                   *
* ndm      -> numero de dimensoes                                   *
* iCod     -> tipo de residuo                                       *
*          RSCALED - residuo com escala de grandeza                 *
*          RSQRT   - norma p-2 ( norma euclidiana)                  *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* rU       -> residuo das velocidades                               *
* rMass    -> residuo de mass                                       *
*-------------------------------------------------------------------*
* OBS:                                                              *
*            | rvx1 rvx2 ... rvxn |                                 *
* rCellVel = | rvy1 rvy2 ... rvyn |                                 *
*            | rvz1 rvz2 ... rvzn |                                 *
*                                                                   *
*            | ad1x ad1y ad1z |                                     *    
* adVel    = | ad2x ad2y ad1z |                                     *
*            |   ...          |                                     *
             | adnx adny adnz |                                     *
*-------------------------------------------------------------------*
*********************************************************************/
void residualPrime(DOUBLE *restrict vel     ,DOUBLE *restrict adVel
                  ,DOUBLE *restrict rCellVel,DOUBLE *restrict rCellMass
                  ,DOUBLE *restrict rU      ,DOUBLE *restrict rMass
                  ,INT const nEl            ,short const ndm
                  ,short const iCod)
{

  DOUBLE maxV[3], sum[3], mod, tmp, v, rScale,ad;
  DOUBLE *p;
  INT i, j, lNeq;

/*...*/
  maxV[0] = maxV[1] = maxV[2] = 0.e0;
  sum[0] = sum[1] = sum[2] = 0.e0;
  for (j = 0; j<ndm; j++) {
    rU[j] = 0.e0;
  }
/*...................................................................*/

/*...*/
  switch (iCod) {

/*... norma euclidiana*/
    case RSQRT:
/*...*/
      for (j = 0; j<ndm; j++) {
        p = &rCellVel[j*nEl];
        rU[j] = sqrt(dot(p, p, nEl));
      }
/*...................................................................*/

/*...*/
      *rMass = sqrt(dot(rCellMass, rCellMass, nEl));
/*...................................................................*/
    break;  
/*...................................................................*/

/*... scaled*/
    case RSCALEDSUM:
/*... max(Ap*velP) */
      for (i = 0; i<nEl; i++) {
        for (j = 0; j<ndm; j++) {
          v  = MAT2D(i,j,vel  ,ndm);
          ad = MAT2D(i,j,adVel,ndm);
          mod = fabs(ad*v);
          sum[j] += mod;
        }
      }
/*...................................................................*/

/*... max ( | F - Ax |P / max(Ap*velP) )*/
      for (j = 0; j<ndm; j++) {
        for (i = 0; i<nEl; i++) {
          mod = fabs(MAT2D(j, i, rCellVel, nEl));
          rU[j] += mod;
        }
        if (sum[j] > rU[j] * SZERO)
          rU[j] /= sum[j];
      }
/*...................................................................*/

/*...*/
      tmp = 0.e0;
      for (i = 0; i<nEl; i++) {
        v = fabs(rCellMass[i]);
        tmp += v;
      }
      *rMass = tmp;
/*...................................................................*/
    break;
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 11/08/2016                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* SETPRIMESCHEME : set o metodo prime                               *
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
void setPrimeScheme(char *word, Prime *pr, FILE *fileIn) {

  fscanf(fileIn, "%d", &pr->maxIt);

  fscanf(fileIn, "%lf", &pr->alphaPres);
  fscanf(fileIn, "%lf", &pr->alphaVel);
  fscanf(fileIn, "%lf", &pr->tolPres);
  fscanf(fileIn, "%lf", &pr->tolVel);

  fscanf(fileIn, "%d", &pr->nNonOrth);

  fscanf(fileIn, "%d", &pr->pPrime);

/*...*/
  if (!mpiVar.myId) {
    printf("PRES-VEL  : PRIME\n");
    printf("Maxit     : %d\n" , pr->maxIt);
    printf("alphaPres : %lf\n", pr->alphaPres);
    printf("alphaVel  : %lf\n", pr->alphaVel);
    printf("tolPres   : %e\n" , pr->tolPres);
    printf("tolVel    : %e\n" , pr->tolVel);
    printf("nNonOrth  : %d\n" , pr->nNonOrth);
    printf("pSimple   : %d\n" , pr->pPrime);
  }
/*...................................................................*/

}
/*********************************************************************/