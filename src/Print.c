#include<Print.h>

/********************************************************************* 
 * Data de criacao    : 02/12/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * printFluid: impressao do fluido                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*  
 * m       -> vetor de memoria principal                             * 
 * turbModel -> modelo de turbulencia                                *
 * eModel    -> modelo da equacao de energia                         *
 * loadVel   -> deficicao de cargas velocidade                       *
 * loadPres  -> deficicao de cargas pres                             *
 * loadTemp  -> deficicao de cargas temperatura                      *
 * opt       -> opcoes de arquivo                                    *
 * mesh0     -> malha global                                         *
 * mesh      -> malha particionada                                   *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*  
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void printFluid(Memoria *m        
               ,Turbulence turbModel,EnergyModel eModel
               ,PartMesh *pMesh     ,Scheme sc
               ,Loads *loadsVel     ,Loads *loadsPres 
               ,Loads *loadsTemp    ,FileOpt opt
               ,Mesh *mesh0         ,Mesh *mesh  
               ,char *preName       ,char *nameOut){
 
  void *dum=NULL;
  DOUBLE ndfVel;
  FILE *fileOut=NULL;

/*...*/
  ndfVel = max(mesh->ndfF - 1, mesh->ndfFt - 2);
/*...................................................................*/

/*... reconstruindo do gradiente (Pres)*/
  tm.rcGradPres = getTimeC() - tm.rcGradPres;
  rcGradU(m                      ,loadsPres
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

/*... interpolacao das variaveis da celulas para pos nos (GradPres)*/
  interCellNode(m                   ,loadsPres
               ,mesh->node.gradPres,mesh->elm.gradPres
               ,mesh->elm.node     ,mesh->elm.geomType            
               ,mesh->elm.geom.cc  ,mesh->node.x  
               ,mesh->elm.geom.xm
               ,mesh->elm.nen      ,mesh->elm.adj.nViz
               ,mesh->elm.faceRpres,mesh->elm.faceLoadPres
               ,&pMesh->iNo          
               ,mesh->numelNov     ,mesh->numel        
               ,mesh->nnodeNov     ,mesh->nnode 
               ,mesh->maxNo        ,mesh->maxViz   
               ,mesh->ndm          ,1
               ,mesh->ndm      
               ,false              ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (pres)*/
  interCellNode(m                 ,loadsPres
               ,mesh->node.pressure,mesh->elm.pressure   
               ,mesh->elm.node     ,mesh->elm.geomType            
               ,mesh->elm.geom.cc  ,mesh->node.x  
               ,mesh->elm.geom.xm
               ,mesh->elm.nen      ,mesh->elm.adj.nViz
               ,mesh->elm.faceRpres,mesh->elm.faceLoadPres  
               ,&pMesh->iNo          
               ,mesh->numelNov     ,mesh->numel        
               ,mesh->nnodeNov     ,mesh->nnode 
               ,mesh->maxNo        ,mesh->maxViz   
               ,1                  ,1
               ,mesh->ndm      
               ,opt.bconditions    ,2);
/*...................................................................*/

/*... calculo da matrix jacobiana das velocidades
                            | du1dx1 du1dx2 du1dx3 |   
                            | du2dx1 du2dx2 du2dx3 |   
                            | du3dx1 du3dx2 du3dx3 |
*/        
  tm.rcGradVel  = getTimeC() - tm.rcGradVel;
  rcGradU(m                     ,loadsVel
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

/*... interpolacao das variaveis da celulas para pos nos (gradVel)*/
  interCellNode(m                 ,loadsVel
                   ,mesh->node.gradVel ,mesh->elm.gradVel    
                   ,mesh->elm.node     ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc  ,mesh->node.x  
                   ,mesh->elm.geom.xm
                   ,mesh->elm.nen      ,mesh->elm.adj.nViz
                   ,mesh->elm.faceRvel ,mesh->elm.faceLoadVel 
                   ,&pMesh->iNo          
                   ,mesh->numelNov     ,mesh->numel        
                   ,mesh->nnodeNov     ,mesh->nnode 
                   ,mesh->maxNo        ,mesh->maxViz   
                   ,ndfVel             ,mesh->ndm
                   ,mesh->ndm      
                   ,false              ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (vel)*/
  interCellNode(m                 ,loadsVel
                   ,mesh->node.vel     ,mesh->elm.vel        
                   ,mesh->elm.node     ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc  ,mesh->node.x  
                   ,mesh->elm.geom.xm
                   ,mesh->elm.nen      ,mesh->elm.adj.nViz
                   ,mesh->elm.faceRvel ,mesh->elm.faceLoadVel 
                   ,&pMesh->iNo          
                   ,mesh->numelNov     ,mesh->numel        
                   ,mesh->nnodeNov     ,mesh->nnode 
                   ,mesh->maxNo        ,mesh->maxViz   
                   ,mesh->ndm          ,1
                   ,mesh->ndm      
                   ,opt.bconditions    ,2);
/*...................................................................*/

/*... reconstruindo do gradiente (Energia)*/
  if(mesh->ndfFt){
    rcGradU(m                 ,loadsTemp
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
              ,mesh->elm.temp       ,mesh->elm.gradTemp  
              ,mesh->node.temp      ,sc.rcGrad
              ,mesh->maxNo          ,mesh->maxViz
              ,1, mesh->ndm
              ,&pMesh->iNo          ,&pMesh->iEl
              ,mesh->numelNov       ,mesh->numel
              ,mesh->nnodeNov       ,mesh->nnode); 
/*.................................................................. */


/*... interpolacao das variaveis da celulas para pos nos (GradEnergy)*/
    interCellNode(m            ,loadsTemp  
              ,mesh->node.gradTemp  ,mesh->elm.gradTemp  
              ,mesh->elm.node       ,mesh->elm.geomType
              ,mesh->elm.geom.cc    ,mesh->node.x
              ,mesh->elm.geom.xm
              ,mesh->elm.nen        ,mesh->elm.adj.nViz
              ,mesh->elm.faceRenergy,mesh->elm.faceLoadEnergy
              ,&pMesh->iNo
              ,mesh->numelNov       ,mesh->numel
              ,mesh->nnodeNov       ,mesh->nnode
              ,mesh->maxNo          ,mesh->maxViz
              ,mesh->ndm            ,1
              ,mesh->ndm
              ,false                ,2);
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (energy)*/
    interCellNode(m           ,loadsTemp
             ,mesh->node.temp      ,mesh->elm.temp
             ,mesh->elm.node       ,mesh->elm.geomType
             ,mesh->elm.geom.cc    ,mesh->node.x
             ,mesh->elm.geom.xm
             ,mesh->elm.nen        ,mesh->elm.adj.nViz
             ,mesh->elm.faceRenergy,mesh->elm.faceLoadEnergy
             ,&pMesh->iNo
             ,mesh->numelNov       ,mesh->numel
             ,mesh->nnodeNov       ,mesh->nnode
             ,mesh->maxNo          ,mesh->maxViz
             ,1                    ,1
             ,mesh->ndm
             ,opt.bconditions      ,2);
  }
/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (density)*/
  interCellNode(m               ,loadsTemp
             ,mesh->node.densityFluid,mesh->elm.densityFluid
             ,mesh->elm.node         ,mesh->elm.geomType
             ,mesh->elm.geom.cc      ,mesh->node.x
             ,mesh->elm.geom.xm      
             ,mesh->elm.nen          ,mesh->elm.adj.nViz
             ,dum                    ,dum                           
             ,&pMesh->iNo            
             ,mesh->numelNov         ,mesh->numel
             ,mesh->nnodeNov         ,mesh->nnode
             ,mesh->maxNo            ,mesh->maxViz
             ,3                      ,1
             ,mesh->ndm              
             ,false                  ,2);

/*...................................................................*/

/*... interpolacao das variaveis da celulas para pos nos (dViscosity)*/
  interCellNode(m               ,loadsTemp
             ,mesh->node.dViscosity  ,mesh->elm.dViscosity
             ,mesh->elm.node         ,mesh->elm.geomType
             ,mesh->elm.geom.cc      ,mesh->node.x
             ,mesh->elm.geom.xm      
             ,mesh->elm.nen          ,mesh->elm.adj.nViz
             ,dum                    ,dum                           
             ,&pMesh->iNo            
             ,mesh->numelNov         ,mesh->numel
             ,mesh->nnodeNov         ,mesh->nnode
             ,mesh->maxNo            ,mesh->maxViz
             ,1                      ,1
             ,mesh->ndm              
             ,false                  ,2);

/*...................................................................*/

/*...*/
  if (turbModel.fTurb) {
    interCellNode(m                     ,loadsVel
                   ,mesh->node.eddyViscosity ,mesh->elm.eddyViscosity        
                   ,mesh->elm.node           ,mesh->elm.geomType            
                   ,mesh->elm.geom.cc        ,mesh->node.x  
                   ,mesh->elm.geom.xm       
                   ,mesh->elm.nen            ,mesh->elm.adj.nViz
                   ,dum                      ,dum                   
                   ,&pMesh->iNo                
                   ,mesh->numelNov           ,mesh->numel        
                   ,mesh->nnodeNov           ,mesh->nnode 
                   ,mesh->maxNo              ,mesh->maxViz   
                   ,1                        ,1
                   ,mesh->ndm               
                   ,false                    ,2);
  }                                 
/*...................................................................*/

/*...*/
  if(!mpiVar.myId ){
    fName(preName,sc.ddt.timeStep,0,21,&nameOut);
/*...*/
    wResVtkFluid(m                  , mesh0->node.x      
               , mesh0->elm.node         , mesh0->elm.mat    
               , mesh0->elm.nen          , mesh0->elm.geomType
               , mesh0->elm.pressure     , mesh0->node.pressure
               , mesh0->elm.gradPres     , mesh0->node.gradPres  
               , mesh0->elm.vel          , mesh0->node.vel      
               , mesh0->elm.gradVel      , mesh0->node.gradVel 
               , mesh0->elm.temp         , mesh0->node.temp   
               , mesh0->elm.gradTemp     , mesh0->node.gradTemp
               , mesh0->elm.eddyViscosity, mesh0->node.eddyViscosity
               , mesh0->elm.densityFluid , mesh0->node.densityFluid
               , mesh0->elm.dViscosity   , mesh0->node.dViscosity
               , mesh0->elm.specificHeat , mesh0->elm.tConductivity
               , mesh0->elm.wallParameters                  
               , mesh0->nnode            , mesh0->numel  
               , mesh0->ndm              , mesh0->maxNo 
               , mesh0->numat            , ndfVel
               , nameOut                 , opt
               , eModel.fKelvin
               , sc.ddt                  , fileOut);  
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/