#include<Adjcency.h>
/*********************************************************************
 * VIZ : calcula os elementos vinzinhos                              *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * m     -> memoria principal                                        *
 * el    -> conectividade nodal                                      *
 * nelcon-> indefinido                                               *
 * nen   -> numero de nos por elemento                               *
 * nnode -> numero de nos da malha                                   *
 * numel -> numero de elmentos                                       *
 * maxNo -> numero de nos por elementos maximo da malha              * 
 * maxViz-> numero vizinhos por elementos maximo da malha            * 
 * ndm   -> numero deminseoes                                        * 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * nelcon-> vizinhos dos elementos                                   *
 * ------------------------------------------------------------------*
 * *******************************************************************/
void viz(Memoria *m        ,INT *el           ,INT *nelcon
        ,short *nen        ,INT const nnode   ,INT const numel  
        ,short const maxNo ,short const maxViz,short const ndm){
  
  INT *nodcon,nEdge;
  
  HccaAlloc(INT,m,nodcon      ,nnode ,"nodcon",_AD_);
  
  if( ndm == 2)
    adj2d(el,nodcon,nelcon,nen,nnode,numel,maxNo,maxViz,&nEdge);

  else if( ndm == 3)
     adjHex8(el         ,nodcon         ,nelcon
            ,nnode      ,numel 
            ,maxNo      ,maxViz);
  
/*  for(i=0;i<numel;i++)
    printf("nel %d %d %d %d %d %d %d\n",i+1
                  ,MAT2D(i,0,nelcon,maxViz)
                  ,MAT2D(i,1,nelcon,maxViz)
                  ,MAT2D(i,2,nelcon,maxViz)
                  ,MAT2D(i,3,nelcon,maxViz)
                  ,MAT2D(i,4,nelcon,maxViz)
                  ,MAT2D(i,5,nelcon,maxViz));*/
  HccaDealloc(m,nodcon,"nodcon",_AD_);


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
 * nnode -> numero de nos da malha                                   *
 * numel -> numero de elmentos                                       *
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
void adj2d(INT *el         ,INT *nodcon       ,INT *nelcon
          ,short *nen      ,INT const nnode   , INT const numel 
         ,short const maxNo,short const maxViz,INT *nEdge){
 
 
  INT  i,nel1,nel2,no1=0,no2=0,no21=0,no22=0;
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

/*... quadrilateros*/
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

/*********************************************************************
 * adj2D : calcula os elementos vinzinhos                            *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * el    -> conectividade nodal                                      *
 * nodcon-> indefinido                                               *
 * nelcon-> indefinido                                               *
 * nen   -> numero de nos por elemento                               *
 * nnode -> numero de nos da malha                                   *
 * numel -> numero de elmentos                                       *
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
void adjHex8(INT *el         ,INT *nodcon         ,INT *nelcon
            ,INT const nnode    , INT const numel 
            ,short const maxNo,short const maxViz){
  
  INT nel1,nel2,node[4];
  short k,j;
  bool imiss;

  for(nel1=0;nel1<numel;nel1++)
    for(j=0;j<maxViz;j++) 
      MAT2D(nel1,j,nelcon,maxViz) = -1;
  
  for(nel1 = 0;nel1 < nnode;nel1++)
    nodcon[nel1] = -1; 

  do{
    imiss = false;
/*...*/
    for(nel1=0;nel1<numel;nel1++){
/*...*/
      for(j=0;j<6;j++){
        if( MAT2D(nel1,j,nelcon,maxViz) == -1){
          hexa8fNod(nel1        ,j
                   ,el          ,node
                   ,maxNo);                       
          if(  nodcon[node[0]] == -1 || nodcon[node[1]] == -1  
            || nodcon[node[2]] == -1 || nodcon[node[3]] == -1){
            nodcon[node[0]]= nel1;
            nodcon[node[1]]= nel1;
            nodcon[node[2]]= nel1;
            nodcon[node[3]]= nel1;
            imiss = true;
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*...*/    
    for(nel1=0;nel1<numel;nel1++){
      for(j=0;j<6;j++){
        if( MAT2D(nel1,j,nelcon,maxViz) == -2){
          hexa8fNod(nel1        ,j
                   ,el          ,node
                   ,maxNo);                       
/*...................................................................*/
          nel2 = nodcon[node[0]];
          if( nel2 > -1){
              if(nel2 == nodcon[node[1]] && nel2 == nodcon[node[2]] &&
                 nel2 == nodcon[node[3]] && nel2 != nel1){
                 k = hexa8face(nel2,el,node,maxNo);
                 if( k == -1) {
                   printf("adjHex8: Erro na vizinhaca!!!\n");
                   exit(EXIT_FAILURE); 
                 }
                 MAT2D(nel2,k,nelcon,maxViz) = nel1;
                 MAT2D(nel1,j,nelcon,maxViz) = nel2;
                 nodcon[node[0]] = -1;
                 nodcon[node[1]] = -1;
                 nodcon[node[2]] = -1;
                 nodcon[node[3]] = -1;
                 imiss = true;
              }
/*...................................................................*/
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/
    } 
/*...................................................................*/

/*...*/    
    for(nel1=0;nel1<numel;nel1++){
      for(j=0;j<6;j++){
        if( MAT2D(nel1,j,nelcon,maxViz) == -1){
          hexa8fNod(nel1        ,j
                   ,el          ,node
                   ,maxNo);                       
          nel2 = nodcon[node[0]];
          if(nel2 == nodcon[node[1]] && nel2 == nodcon[node[2]] &&
             nel2 == nodcon[node[3]] && nel2 == nel1){
              MAT2D(nel1,j,nelcon,maxViz) = -2;
              nodcon[node[0]] = -1;
              nodcon[node[1]] = -1;
              nodcon[node[2]] = -1;
              nodcon[node[3]] = -1;
              imiss = true;
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/
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

void hexa8fNod(INT const nEl     ,short const face
             ,INT *restrict el   ,INT *restrict nodeFace
             ,short const maxNo)                       
{
  short isNodHex[][4] = {{0,3,2,1}
                        ,{4,5,6,7}
                        ,{0,1,5,4}
                        ,{1,2,6,5}
                        ,{2,3,7,6}
                        ,{3,0,4,7}};
  short i;

  for(i=0;i<4;i++)
    nodeFace[i] = MAT2D(nEl,isNodHex[face][i],el,maxNo) - 1;


}

short hexa8face(INT const nEl         ,INT *restrict el
               ,INT *restrict nodeFace,short const maxNo )
{
  INT no[4];
  short face = -1,nFace,j;
/*possiveis numeracoes de faces*/
  short ind[][3] = {{1,2,3},{2,3,0},{3,0,1},{0,1,2}};

  for(nFace=0;nFace<6;nFace++){
    hexa8fNod(nEl,nFace,el,no,maxNo);
    for(j=0;j<4;j++){
      if(no[0] == nodeFace[j]){
        if( no[1] == nodeFace[ind[j][2]] &&
            no[2] == nodeFace[ind[j][1]] &&
            no[3] == nodeFace[ind[j][0]] ) 
        return nFace; 
      }
    }
  }
  return face; 
}                
