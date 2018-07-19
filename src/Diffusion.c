#include<Diffusion.h>
/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 12/05/2018                                   *
 *-------------------------------------------------------------------*
 * DIFFUSION :Resolucao do problema de difusao pura no passo de tempo* 
 * n+1                                                               * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m        -> nao definido                                          *
 * loadsDif -> definicoes de cargas                                  *
 * dModel   -> configuracoes do modelo difusivo                      *
 * mesh0    -> parametros                                            * 
 * mesh     -> parametros                                            * 
 * sistEqD  -> sistema de equacoes                                   * 
 * solvD    -> sistema de equacoes                                   * 
 * sc       -> tecnica de discretizacao temporal                     * 
 * pMesh    -> varaives de particionamento                           * 
 * prop     -> propeidades variaveis
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
void diffusion(Memoria *m       ,Loads *loadsDif,DiffModel *dModel
              ,Mesh *mesh0      ,Mesh *mesh     ,SistEq *sistEqD
              ,Solv *solvD      ,Scheme *sc     ,PartMesh *pMesh 
              ,PropVarCD *prop  ,FileOpt *opt   ,char *preName
              ,char *nameOut    ,FILE *fileOut) 
{   
  bool fDensity=prop->fDensity,fCoefDiff=prop->fCeofDiff;
  short unsigned i,jj;
  DOUBLE rCell,rCell0,conv;

  rCell = rCell0 = conv = 0.e0;

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
  if(sc->ddt.flag){
    tm.CellTransientD1 = getTimeC() - tm.CellTransientD1;
    cellTransient(mesh->elm.geom.volume   ,sistEqD->id     
                 ,mesh->elm.u0D1          ,mesh->elm.uD1
                 ,mesh->elm.densityUd1    ,sistEqD->b0
                 ,sc->ddt                 ,mesh->numelNov 
                 ,mesh->ndfD[0]           ,true);
/*... u(n-1) = u(n)*/
    alphaProdVector(1.e0              ,mesh->elm.uD1
                   ,mesh->numel       ,mesh->elm.u0D1); 
/*...................................................................*/
    tm.CellTransientD1 = getTimeC() - tm.CellTransientD1;
  }
/*...................................................................*/

/*... correcao nao ortoganal*/ 
  for(i=0,jj=0;i<sc->nlD1.maxIt;i++){

/*... calculo de: A(i),b(i)*/
    tm.systFormD1 = getTimeC() - tm.systFormD1;
    systFormDif(loadsDif                ,&sc->diffD1
               ,dModel 
               ,mesh->elm.node          ,mesh->elm.adj.nelcon  
               ,mesh->elm.nen           ,mesh->elm.adj.nViz   
               ,mesh->elm.cellFace      ,mesh->face.owner
               ,mesh->elm.geom.volume   ,mesh->elm.geom.dcca
               ,mesh->elm.geom.xmcc     
               ,mesh->face.mksi         ,mesh->face.ksi
               ,mesh->face.eta          ,mesh->face.area       
               ,mesh->face.normal       ,mesh->face.xm    
               ,mesh->face.mvSkew       ,mesh->face.vSkew
               ,mesh->elm.geomType      ,mesh->elm.material.prop 
               ,mesh->elm.material.type ,mesh->elm.mat   
               ,mesh->elm.densityUd1    ,mesh->elm.cDiffD1
               ,sistEqD->ia             ,sistEqD->ja      
               ,sistEqD->al             ,sistEqD->ad       
               ,sistEqD->b              ,sistEqD->id       
               ,mesh->elm.faceRd1       ,mesh->elm.faceLoadD1  
               ,mesh->elm.uD1           ,mesh->elm.gradUd1           
               ,mesh->elm.rCellUd1      ,&sc->ddt
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
      conv   = rCell0*sc->nlD1.tol;
    }
    else
      rCell  = sqrt(dot(mesh->elm.rCellUd1 
                   ,mesh->elm.rCellUd1 
                   ,mesh->numelNov));
    if(!mpiVar.myId && jj == sc->nlD1.pPlot)
    {
      jj = 0;
      printf("it: %9d %.6e %0.6e\n",i,rCell/rCell0,rCell);
      if(opt->fItPlot)  
        fprintf(opt->fileItPlot[FITPLOTD1]
               ,"%9d %.6e %0.6e\n",i,rCell/rCell0,rCell);
    }
    jj++;
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
           ,&sistEqD->iNeq     ,&sistEqD->omp
           ,solvD->tol         ,solvD->maxIt     
           ,sistEqD->storage   ,solvD->solver
           ,solvD->fileSolv    ,solvD->log  
           ,false              ,sistEqD->unsym);
    tm.solvD1 = getTimeC() - tm.solvD1;
/*...................................................................*/

/*... x -> uD1*/
    updateCellValue(mesh->elm.uD1 ,sistEqD->x
                   ,sistEqD->id  ,&sistEqD->iNeq
                   ,mesh->numel  ,mesh->ndfD[0]
                   ,dModel->fRes ,true);
/*...................................................................*/

/*... reconstruindo do gradiente*/
    tm.rcGradD1 = getTimeC() - tm.rcGradD1;
    rcGradU( m                       , loadsDif
           , mesh->elm.node          , mesh->elm.adj.nelcon
           , mesh->node.x  
           , mesh->elm.nen           , mesh->elm.adj.nViz
           , mesh->elm.cellFace      , mesh->face.owner
           , mesh->elm.geom.volume   , mesh->elm.geom.dcca
           , mesh->elm.geom.xmcc     , mesh->elm.geom.cc
           , mesh->face.mksi         , mesh->face.ksi
           , mesh->face.eta          , mesh->face.area
           , mesh->face.normal       , mesh->face.xm
           , mesh->face.mvSkew       , mesh->face.vSkew
           , mesh->elm.geomType      , mesh->elm.material.prop
           , mesh->elm.material.type , mesh->elm.mat
           , mesh->elm.leastSquare   , mesh->elm.leastSquareR
           , mesh->elm.faceRd1       , mesh->elm.faceLoadD1    
           , mesh->elm.uD1           , mesh->elm.gradUd1                 
           , mesh->node.uD1          , sc->rcGrad
           , mesh->maxNo             , mesh->maxViz
           , mesh->ndfD[0]           , mesh->ndm
           , &pMesh->iNo             , &pMesh->iEl  
           , mesh->numelNov          , mesh->numel        
           , mesh->nnodeNov          , mesh->nnode);   
    tm.rcGradD1 = getTimeC() - tm.rcGradD1;
/*...................................................................*/

/*...*/
    if(fDensity)
      updateDensityCD( &prop->den
                     , mesh->elm.uD1 , mesh->elm.densityUd1
                     , mesh->numel   , PROP_UPDATE_NL_LOOP);
/*...................................................................*/

/*...*/
    if (fCoefDiff)
      updateProp(&prop->ceofDiff  , mesh->elm.uD1
               , mesh->elm.cDiffD1, mesh->numel);
/*...................................................................*/

  } 
/*...................................................................*/

/*...*/
  if (fDensity)
    updateDensityCD(&prop->den
                  , mesh->elm.uD1, mesh->elm.densityUd1
                  , mesh->numel  , PROP_UPDATE_OLD_TIME);
/*...................................................................*/
}
/*********************************************************************/




/*  systFormDifOld(loadsDif, &sc->diffD1
, mesh->elm.node, mesh->elm.adj.nelcon
, mesh->elm.nen, mesh->elm.adj.nViz
, mesh->elm.geom.volume, mesh->elm.geom.dcca
, mesh->elm.geom.xmcc
, mesh->elm.geom.mksi, mesh->elm.geom.ksi
, mesh->elm.geom.eta, mesh->elm.geom.fArea
, mesh->elm.geom.normal, mesh->elm.geom.xm
, mesh->elm.geom.mvSkew, mesh->elm.geom.vSkew
, mesh->elm.geomType, mesh->elm.material.prop
, mesh->elm.material.type, mesh->elm.mat
, mesh->elm.densityUd1
, sistEqD->ia, sistEqD->ja
, sistEqD->al, sistEqD->ad
, sistEqD->b, sistEqD->id
, mesh->elm.faceRd1, mesh->elm.faceLoadD1
, mesh->elm.uD1, mesh->elm.gradUd1
, mesh->elm.rCellUd1, &sc->ddt
, sistEqD->neq, sistEqD->neqNov
, sistEqD->nad, sistEqD->nadr
, mesh->maxNo, mesh->maxViz
, mesh->ndm, mesh->numelNov
, mesh->ndfD[0], sistEqD->storage
, true, true
, true, sistEqD->unsym);*/