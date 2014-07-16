#include<Adjcency.h>
/*********************************************************************
 * VIZ : calcula os elementos vinzinhos                              *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * m     -> memoria principal                                        *
 * el    -> conectividade nodal                                      *
 * nelcon-> indefinido                                               *
 * nViz  -> indefinido                                               *
 * nen   -> numero de nos por elemento                               *
 * nnode -> numero de nos da malha                                   *
 * numel -> numero de elmentos                                       *
 * maxNo -> numero de nos por elementos maximo da malha              * 
 * maxViz-> numero vizinhos por elementos maximo da malha            * 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * nelcon-> vizinhas dos elementos                                   *
 * nViz  -> numero de vizinhos por elemento                          *
 * ------------------------------------------------------------------*
 * *******************************************************************/
void viz(Memoria *m ,INT *el    ,INT *nelcon
        ,short *nViz,short *nen ,INT nnode
        ,INT numel ,short maxNo ,short maxViz){
  
  INT *nodcon,nEdge;
  
  Myalloc(INT,m,nodcon      ,nnode ,"nodcon",_AD_);
  
  adj2d(el,nodcon,nelcon,nen,numel,nnode,maxNo,maxViz,&nEdge);
  
//for(i=0;i<numel;i++)
//  nViz[i] = nen[i];
  
  Mydealloc(m,nodcon,"nodcon",_AD_);


}
/*********************************************************************/

/*********************************************************************
 * adj2D : calcula os elementos vinzinhos                            *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * el    -> conectividade nodal                                      *
 * nodcon-> indefinido                                               *
 * nelcon-> indefinido                                               *
 * nen   -> numero de nos por elemento                               *
 * numel -> numero de elmentos                                       *
 * nnode -> numero de nos da malha                                   *
 * maxNo -> numero de nos por elementos maximo da malha              * 
 * maxViz-> numero vizinhos por elementos maximo da malha            * 
 * nEdge -> indefinifo                                               * 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * nelcon-> vizinhas dos elementos                                   *
 * nEdge -> numero de arestas                                        * 
 * ------------------------------------------------------------------*
 * *******************************************************************/
void adj2d(INT *el  ,INT *nodcon,INT *nelcon,short *nen
          ,INT numel, INT nnode ,short maxNo , short maxViz
          ,INT *nEdge){
 
 
  INT  i,nel1,nel2,no1,no2,no21,no22;
  int   j,k;
  short is1,is2; 
  bool imiss;
  short int isnod3[][2]= {{0,1},{1,2},{2,0}};  
  short int isnod4[][2]= {{0,1},{1,2},{2,3},{3,0}};  
  
  for(nel1=0;nel1<numel;nel1++)
    for(j=0;j<maxViz;j++) 
      MAT2D(nel1,j,nelcon,maxViz) = -1;
 
  for(i = 0;i < nnode;i++)
    nodcon[i] = -1; 
  
  *nEdge = 0;
  do{
    imiss = false;
    for(nel1=0;nel1 < numel;nel1++){ 
      is1 = nen[nel1];
      for(j=0;j<is1;j++){
        if(MAT2D(nel1,j,nelcon,maxViz) == -1){
/*... triangulo*/
          if(is1 == 3) {
            no1 = MAT2D(nel1,isnod3[j][0],el,maxNo)-1; 
            no2 = MAT2D(nel1,isnod3[j][1],el,maxNo)-1; 
          }
/*...................................................................*/

/*... quadrilatero*/
          else if(is1 == 4){
            no1 = MAT2D(nel1,isnod4[j][0],el,maxNo)-1; 
            no2 = MAT2D(nel1,isnod4[j][1],el,maxNo)-1; 
          }
/*...................................................................*/
          if( nodcon[no1] == -1 || nodcon[no2] == -1){
            nodcon[no1] = nodcon[no2] = nel1;
            imiss       = true;
          }
        }
      }
    }
    
/*...*/ 
    for(nel1=0;nel1 < numel;nel1++){ 
      is1 = nen[nel1];
      for(j=0;j<is1;j++){
        if(MAT2D(nel1,j,nelcon,maxViz) == -1){
/*... triangulo*/
          if(is1 == 3){
           no1  = MAT2D(nel1,isnod3[j][0],el,maxNo)-1; 
           no2  = MAT2D(nel1,isnod3[j][1],el,maxNo)-1;
          }
/*...................................................................*/ 

/*... quadrilateros*/
          else if(is1 == 4){
           no1  = MAT2D(nel1,isnod4[j][0],el,maxNo)-1; 
           no2  = MAT2D(nel1,isnod4[j][1],el,maxNo)-1;
          }
/*...................................................................*/ 
          nel2 = nodcon[no1];
          if(nel2 > -1){
            if(nel2 == nodcon[no2] && nel2 != nel1){
              is2 = nen[nel2];
/*... triangulos*/
              for(k=0;k<is2;k++){
/*... triangulo*/
                if(is2 == 3){
                  no21  = MAT2D(nel2,isnod3[k][0],el,maxNo)-1; 
                  no22  = MAT2D(nel2,isnod3[k][1],el,maxNo)-1;
                }
/*...................................................................*/

/*... tquadrilateros*/
                else if(is2 == 4){
                  no21  = MAT2D(nel2,isnod4[k][0],el,maxNo)-1; 
                  no22  = MAT2D(nel2,isnod4[k][1],el,maxNo)-1;
                }
/*...................................................................*/
                if( (no21 == no2) && (no22 == no1)){
                      MAT2D(nel1,j,nelcon,maxViz) = nel2;
                      MAT2D(nel2,k,nelcon,maxViz) = nel1;
                      nodcon[no1] = nodcon[no2] = -1;
                      imiss       = true; 
                      (*nEdge)++;
                 } 
                 
              }
/*...................................................................*/
            }
          } 
        }
      }
    }
/*...................................................................*/

    for(nel1=0;nel1 < numel;nel1++){ 
      is1 = nen[nel1];
      for(j=0;j<is1;j++){
        if(MAT2D(nel1,j,nelcon,maxViz) == -1){
/*... triangulo*/
          if( is1 == 3){
            no1  = MAT2D(nel1,isnod3[j][0],el,maxNo)-1; 
            no2  = MAT2D(nel1,isnod3[j][1],el,maxNo)-1;
          }
/*...................................................................*/

/*... quadrilateros*/
          if( is1 == 4){
            no1  = MAT2D(nel1,isnod4[j][0],el,maxNo)-1; 
            no2  = MAT2D(nel1,isnod4[j][1],el,maxNo)-1;
          }
/*...................................................................*/
          if(nodcon[no1] == nodcon[no2] && nodcon[no1] == nel1){
            MAT2D(nel1,j,nelcon,maxViz) = -2;
            nodcon[no1] = nodcon[no2]   =  -1;
            imiss       = true;
            (*nEdge)++;
          }
/*...................................................................*/
        }
      }
    }
/*...................................................................*/
  }while(imiss);  

/*...*/  
  for(nel1=0;nel1<numel;nel1++)
    for(j=0;j<maxViz;j++)
      MAT2D(nel1,j,nelcon,maxViz) += 1;
/*...................................................................*/
 
}
/*********************************************************************/
