#include<Transport.h>
/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 29/06/2018                                   *
 *-------------------------------------------------------------------*
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
void transport(Memoria *m      ,Loads *loadsTrans,TransModel *tModel
              ,Mesh *mesh0     ,Mesh *mesh       ,SistEq *sistEqT
              ,Solv *solvT     ,Scheme *sc       ,PartMesh *pMesh 
              ,PropVarCD *prop ,FileOpt *opt     ,char *preName  
              ,char *nameOut   ,FILE *fileOut) 
{   
  bool fDensity = prop->fDensity, fCoefDiff = prop->fCeofDiff;
  short unsigned i, jj;
  DOUBLE rCell,rCell0,conv;

  rCell = rCell0 = conv = 0.e0;

  zero(sistEqT->b0,sistEqT->neqNov,DOUBLEC);
  zero(sistEqT->x ,sistEqT->neq   ,DOUBLEC);
/*... restricoes por centro de celula u0 e cargas por volume b0*/
  tm.CellPloadT1 = getTimeC() - tm.CellPloadT1;
  cellPload(loadsTrans            ,mesh->elm.geom.cc 
           ,mesh->elm.faceRt1    
           ,mesh->elm.geom.volume,sistEqT->id 
           ,mesh->elm.uT1        ,sistEqT->b0
           ,mesh->numelNov       ,mesh->ndfT[0]
           ,mesh->ndm            ,mesh->maxViz);
  tm.CellPloadT1 = getTimeC() - tm.CellPloadT1;
/*...................................................................*/

/*... discretizacao temporal*/
  if(sc->ddt.flag){
    tm.CellTransientT1 = getTimeC() - tm.CellTransientT1;
    cellTransient(mesh->elm.geom.volume   ,sistEqT->id     
                 ,mesh->elm.u0T1          ,mesh->elm.uT1
                 ,mesh->elm.densityUt1    ,sistEqT->b0
                 ,sc->ddt                 ,mesh->numelNov   
                 ,mesh->ndfT[0]           ,true);
/*... u(n-1) = u(n)*/
    alphaProdVector(1.e0,mesh->elm.uT1
                   ,mesh->numel       ,mesh->elm.u0T1); 
/*...................................................................*/
    tm.CellTransientT1 = getTimeC() - tm.CellTransientT1;
  }
/*...................................................................*/
  
/*... correcao nao ortoganal*/   
  for(i=0,jj=0;i<sc->nlT1.maxIt;i++){
 
/*... calculo de: A(i),b(i)*/
    tm.systFormT1 = getTimeC() - tm.systFormT1;
    systFormTrans(loadsTrans            
               ,&sc->advT1              ,&sc->diffT1
               ,tModel
               ,mesh->elm.node          ,mesh->elm.adj.nelcon  
               ,mesh->elm.nen           ,mesh->elm.adj.nViz   
               ,mesh->elm.cellFace      ,mesh->face.owner
               ,mesh->elm.geom.volume   ,mesh->elm.geom.dcca
               ,mesh->elm.geom.xmcc     ,mesh->elm.geom.cc
               ,mesh->face.mksi         ,mesh->face.ksi
               ,mesh->face.eta          ,mesh->face.area
               ,mesh->face.normal       ,mesh->face.xm
               ,mesh->face.mvSkew       ,mesh->face.vSkew
               ,mesh->elm.geomType      ,mesh->elm.material.prop
               ,mesh->elm.material.type ,mesh->elm.mat
               ,mesh->elm.densityUt1.t  ,mesh->elm.cDiffT1
               ,sistEqT->ia             ,sistEqT->ja      
               ,sistEqT->al             ,sistEqT->ad       
               ,sistEqT->b              ,sistEqT->id       
               ,mesh->elm.faceRt1       
               ,mesh->elm.uT1           ,mesh->elm.gradUt1           
               ,mesh->elm.vel                                        
               ,mesh->elm.rCellUt1      ,&sc->ddt
               ,sistEqT->neq            ,sistEqT->neqNov      
               ,sistEqT->nad            ,sistEqT->nadr      
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
    if( i == 0 )
    {
      rCell  = rCell0 = sqrt(dot(mesh->elm.rCellUt1
                            ,mesh->elm.rCellUt1 
                            ,mesh->numelNov));
      conv   = rCell0*sc->nlT1.tol;
    }
    else
      rCell  = sqrt(dot(mesh->elm.rCellUt1 
                   ,mesh->elm.rCellUt1 
                   ,mesh->numelNov));
        
    if(!mpiVar.myId )
    {
      jj = 0;
      printf("it: %9d %.6e %0.6e\n", i, rCell / rCell0, rCell);
      if(opt->fItPlot)  
        fprintf(opt->fileItPlot[FITPLOTT1]
               ,"%9d %.6e %0.6e\n",i,rCell/rCell0,rCell);
    }
    jj++;
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
           ,&sistEqT->iNeq     ,&sistEqT->omp 
           ,solvT->tol         ,solvT->maxIt     
           ,sistEqT->storage   ,solvT->solver
           ,solvT->fileSolv    ,solvT->log  
           ,true               ,sistEqT->unsym);  
    tm.solvT1 = getTimeC() - tm.solvT1;
/*...................................................................*/

/*... x -> uT1*/
    updateCellValue(mesh->elm.uT1 ,sistEqT->x
                   ,sistEqT->id   ,&sistEqT->iNeq
                   ,mesh->numel   ,mesh->ndfT[0]
                   ,tModel->fRes  ,true);
/*...................................................................*/

/*... reconstruindo do gradiente*/
    tm.rcGradT1 = getTimeC() - tm.rcGradT1;
    rcGradU(m                        , loadsTrans
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
           , mesh->elm.material.type 
           , mesh->elm.mat           , mesh->elm.cDiffT1
           , mesh->elm.leastSquare   , mesh->elm.leastSquareR
           , mesh->elm.faceRt1         
           , mesh->elm.uT1           , mesh->elm.gradUt1                 
           , mesh->node.uT1          
           , NULL                   , NULL
           , 0
           , &sc->rcGrad
           , mesh->maxNo             , mesh->maxViz
           , mesh->ndfT[0]           , mesh->ndm
           , &pMesh->iNo             , &pMesh->iEl  
           , mesh->numelNov          , mesh->numel        
           , mesh->nnodeNov          , mesh->nnode
           , prop->fCeofDiff); 
    tm.rcGradT1 = getTimeC() - tm.rcGradT1;
/*...................................................................*/

/*...*/
    if (fDensity)
      updateDensityCD(&prop->den
                    , mesh->elm.uT1, mesh->elm.densityUt1
                    , mesh->numel, PROP_UPDATE_NL_LOOP);
/*...................................................................*/

/*...*/
    if (fCoefDiff)
      updateProp(&prop->ceofDiff, mesh->elm.uT1
                , mesh->elm.cDiffT1, mesh->numel);
/*...................................................................*/

  }
/*...................................................................*/

/*...*/
  if (fDensity)
    updateDensityCD(&prop->den
                   , mesh->elm.uT1, mesh->elm.densityUt1
                   , mesh->numel, PROP_UPDATE_OLD_TIME);
/*...................................................................*/
}
/*********************************************************************/
