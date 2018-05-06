#include<Adjcency.h>

/*********************************************************************
* Data de criacao    : 14/01/2018                                   *
* Data de modificaco : 08/04/2018                                   *
* ------------------------------------------------------------------*
* face2d : numeracao de faces para malhas 2D                        *
* ------------------------------------------------------------------*
* Parametros de entrada:                                            *
* ------------------------------------------------------------------*
* el    -> conectividade nodal                                      *
* nelcon-> indefinido                                               *
* nViz  -> numero de vizinhos por celulas                           *
* cellFace -> nao definido                                          *
* fNode    -> nao definido                                          *
* fOwner   -> nao definido                                          *
* nnode -> numero de nos da malha                                   *
* numel -> numero de elmentos                                       *
* maxNo -> numero de nos por elementos maximo da malha              *
* maxFace-> numero maximo de faces por celulas                      *
* nFaces-> numero total de faces                                    *
* ------------------------------------------------------------------*
* Paramanetros de saida:                                            *
* ------------------------------------------------------------------*
* nelcon-> vizinhos dos elementos                                   *
* ------------------------------------------------------------------*
* *******************************************************************/
void face2d(INT *RESTRICT el         ,INT *RESTRICT nelcon
          ,short *RESTRICT nViz      ,INT *RESTRICT cellFace
          ,INT *fNode                ,INT *fOwner          
          ,INT nnode                 ,INT numel
          ,short maxNo               ,short maxFace
          ,INT nFace){

  INT  nel1,vizNel,no1,no2,no21,no22,iEdge;
  int   j,k;
  short is1,is2;
  short int isnod3[][2]= {{0,1},{1,2},{2,0}};
  short int isnod4[][2]= {{0,1},{1,2},{2,3},{3,0}};
  

  no1 = no2 = no21 = no22 = vizNel = 0;

  for(j=0;j<nFace;j++){
    MAT2D(j,0,fNode ,2) = -1;
    MAT2D(j,1,fNode ,2) = -1;
    MAT2D(j,0,fOwner,2) = -1;
    MAT2D(j,1,fOwner,2) = -1;
  }

  for(nel1=0;nel1 < numel;nel1++){
    is1 = nViz[nel1];
    for(j=0;j<is1;j++)
      MAT2D(nel1,j,cellFace,maxFace) = -1;      
  }

  iEdge = 0;
/*... aresta no interior do domineo*/
  for(nel1=0;nel1 < numel;nel1++){
    is1 = nViz[nel1];
    for(j=0;j<is1;j++){
/*... triangulo*/
      if(is1 == 3) {
        no1 = MAT2D(nel1,isnod3[j][0],el,maxNo)-1;
        no2 = MAT2D(nel1,isnod3[j][1],el,maxNo)-1;
        vizNel = MAT2D(nel1, j, nelcon, maxFace) - 1;
      }
/*...................................................................*/

/*... quadrilatero*/
      else if(is1 == 4){
        no1    = MAT2D(nel1,isnod4[j][0],el,maxNo)-1;
        no2    = MAT2D(nel1,isnod4[j][1],el,maxNo)-1;
        vizNel = MAT2D(nel1,j,nelcon,maxFace) - 1;
      }
/*...................................................................*/

/*...*/
      if( vizNel != -2 && MAT2D(nel1, j, cellFace, maxFace) == -1){
        MAT2D(iEdge, 0, fNode,  2)         = no1+1;
        MAT2D(iEdge, 1, fNode,  2)         = no2+1;
        MAT2D(iEdge, 0, fOwner, 2)         = nel1+1;
        MAT2D(iEdge, 1, fOwner, 2)         = vizNel+1;
        MAT2D(nel1 , j, cellFace, maxFace) = iEdge+1;
        is2 = nViz[vizNel];
        for (k = 0; k < is2; k++) {
/*... triangulo*/
          if (is2 == 3) {
            no21 = MAT2D(vizNel, isnod3[k][0], el, maxNo) - 1;
            no22 = MAT2D(vizNel, isnod3[k][1], el, maxNo) - 1;
          }
/*...................................................................*/

/*... quadrilatero*/
          else if (is2 == 4) {
            no21 = MAT2D(vizNel, isnod4[k][0], el, maxNo) - 1;
            no22 = MAT2D(vizNel, isnod4[k][1], el, maxNo) - 1;
          }
/*...................................................................*/
          if ((no21 == no2) && (no22 == no1)) {
            MAT2D(vizNel, k, cellFace, maxFace) = iEdge + 1;
          }
/*...................................................................*/
        }
        iEdge++;
/*...................................................................*/
      }
/*...................................................................*/
    }
  }
/*...................................................................*/

/*... aresta no contorno*/
  for (nel1 = 0; nel1 < numel; nel1++) {
    is1 = nViz[nel1];
    for (j = 0; j<is1; j++) {
      vizNel = MAT2D(nel1, j, nelcon, maxFace) - 1;
      if (vizNel == -2) {
/*... triangulo*/
        if (is1 == 3) {
          no1 = MAT2D(nel1, isnod3[j][0], el, maxNo) - 1;
          no2 = MAT2D(nel1, isnod3[j][1], el, maxNo) - 1;
        }
/*...................................................................*/

/*... quadrilatero*/
        else if (is1 == 4) {
          no1 = MAT2D(nel1, isnod4[j][0], el, maxNo) - 1;
          no2 = MAT2D(nel1, isnod4[j][1], el, maxNo) - 1;
        }
/*...................................................................*/

/*...*/
        MAT2D(iEdge, 0, fNode, 2) = no1 + 1;
        MAT2D(iEdge, 1, fNode, 2) = no2 + 1;
        MAT2D(iEdge, 0, fOwner, 2) = nel1 + 1;
        MAT2D(nel1, j, cellFace, maxFace) = iEdge + 1;
        iEdge++;
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 15/04/2018                                   *
* Data de modificaco : 26/04/2018                                   *
* ------------------------------------------------------------------*
* face3d : numeracao de faces para malhas 2D                        *
* ------------------------------------------------------------------*
* Parametros de entrada:                                            *
* ------------------------------------------------------------------*
* el    -> conectividade nodal                                      *
* nelcon-> indefinido                                               *
* nViz  -> numero de vizinhos por celulas                           *
* cellFace -> nao definido                                          *
* fNode    -> nao definido                                          *
* fOwner   -> nao definido                                          *
* nnode -> numero de nos da malha                                   *
* numel -> numero de elmentos                                       *
* maxNo -> numero de nos por elementos maximo da malha              *
* maxFace-> numero maximo de faces por celulas                      *
* nFaces-> numero total de faces                                    *
* ------------------------------------------------------------------*
* Paramanetros de saida:                                            *
* ------------------------------------------------------------------*
* nelcon-> vizinhos dos elementos                                   *
* ------------------------------------------------------------------*
* *******************************************************************/
void face3d(INT *RESTRICT el    , INT *RESTRICT nelcon
          , short *RESTRICT nViz, INT *RESTRICT cellFace
          , INT *fNode          , INT *fOwner
          , INT nnode           , INT numel
          , short maxNo         , short maxNoFace
          , short maxFace       , INT nFace) {
  short is1, is2, noPface;
  INT  nel1, vizNel, iFace,node[4];
  int   i ,j, k, idFace;

  vizNel = 0;
  for (k = 0; k<maxNoFace; k++)
    node[k] = 0;

  for (j = 0; j<nFace; j++) {
    for(k=0;k<maxNoFace;k++)
      MAT2D(j, k, fNode, maxNoFace) = -1;
/*...*/
    MAT2D(j, 0, fOwner, 2) = -1;
    MAT2D(j, 1, fOwner, 2) = -1;
  }

  for (nel1 = 0; nel1 < numel; nel1++) {
    is1 = nViz[nel1];
    for (j = 0; j<is1; j++)
      MAT2D(nel1, j, cellFace, maxFace) = -1;
  }

  iFace = 0;
/*... aresta no interior do domineo*/
  for (nel1 = 0; nel1 < numel; nel1++) {
    is1 = nViz[nel1];
    for (j = 0; j<is1; j++) {
/*... tetraedro*/
      if (is1 == 4) {
        tetra4fNod(nel1, j, el, node, maxNo);
        vizNel = MAT2D(nel1, j, nelcon, maxFace) - 1;
        noPface = 3;
      }
/*...................................................................*/

/*... piramide*/
      else if (is1 == 5) {
        pira5fNod(nel1, j, el, node, maxNo);
        vizNel = MAT2D(nel1, j, nelcon, maxFace) - 1;
        if(j)
          noPface = 3;
        else
          noPface = 4;
      }
/*...................................................................*/

/*... hexaedros*/
      else if (is1 == 6) {
        hexa8fNod(nel1  , j, el   , node, maxNo);
        vizNel = MAT2D(nel1, j, nelcon, maxFace) - 1;
        noPface = 4;
      }
/*...................................................................*/

/*...*/
      if (vizNel != -2 && MAT2D(nel1, j, cellFace, maxFace) == -1) {
        for(i=0;i<noPface;i++)
          MAT2D(iFace, i, fNode, maxNoFace) = node[i] + 1;
        
        MAT2D(iFace, 0, fOwner, 2) = nel1 + 1;
        MAT2D(iFace, 1, fOwner, 2) = vizNel + 1;
        
        MAT2D(nel1, j, cellFace, maxFace) = iFace + 1;
        is2 = nViz[vizNel];
/*... tetraedro*/
        if (is2 == 4) {
          idFace = tetra4face(vizNel, el, node, maxNo);
        }
/*...................................................................*/

/*... piramide*/
        else if (is2 == 5) {
          idFace = pira5face(vizNel, el, node, maxNo);
        }
/*...................................................................*/

/*... hexaedro*/
        else if (is2 == 6) {
          idFace = hexa8face(vizNel, el, node, maxNo);
        }
/*...................................................................*/
        if (idFace!=-1) {
          MAT2D(vizNel, idFace, cellFace, maxFace) = iFace + 1;
        }
/*...................................................................*/
        iFace++;
/*...................................................................*/
      }
/*...................................................................*/
    }
  }
/*...................................................................*/

/*... aresta no contorno*/
  for (nel1 = 0; nel1 < numel; nel1++) {
    is1 = nViz[nel1];
    for (j = 0; j<is1; j++) {
      vizNel = MAT2D(nel1, j, nelcon, maxFace) - 1;
      if (vizNel == -2) {
/*... tetraedro*/
        if (is1 == 4) 
        {
          tetra4fNod(nel1, j, el, node, maxNo);
          noPface = 3;
        }
/*...................................................................*/

/*... piramide*/
        else if (is1 == 5) {
          pira5fNod(nel1, j, el, node, maxNo);
          if (j)
            noPface = 3;
          else
            noPface = 4;
        }
/*...................................................................*/

/*... hexaedro*/
        else if (is1 == 6) {
          hexa8fNod(nel1, j, el, node, maxNo);
          noPface = 4;
        }
/*...................................................................*/

/*...*/
        for (i = 0; i<noPface; i++)
          MAT2D(iFace, i, fNode, maxNoFace) = node[i] + 1;
        MAT2D(iFace, 0, fOwner, 2) = nel1 + 1;
        MAT2D(nel1, j, cellFace, maxFace) = iFace + 1;
        iFace++;
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 26/04/2018                                   *
 * ------------------------------------------------------------------*
 * neighbors : calcula os elementos vinzinhos                        *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*
 * m     -> memoria principal                                        *
 * el    -> conectividade nodal                                      *
 * nelcon-> indefinido                                               *
 * nViz  -> numero de vizinhos por celulas                           *
 * nen   -> numero de nos por celulas                                *
 * nFaces-> numero total de faces                                    *
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
void neighbors(Memoria *m          ,INT *RESTRICT el   
              ,INT *RESTRICT nelcon,short *RESTRICT nViz
              ,short *RESTRICT nen ,INT *nFace
              ,INT nnode           ,INT numel
              ,short maxNo         ,short maxViz
              ,short ndm)
{

  INT *nodcon;

/*...*/
  HccaAlloc(INT,m,nodcon      ,nnode ,"nodcon",_AD_);
/*...................................................................*/

/*...*/
  if( ndm == 2)
  {
    adj2d(el   ,nodcon,nelcon
         ,nen  ,nnode ,numel
         ,maxNo,maxViz,nFace);
  }
/*.....................................................................*/

/*...*/
  else if( ndm == 3)
  {
    adj3d(el   , nodcon, nelcon
        , nViz , nnode , numel
        , maxNo, maxViz, nFace);
  
/*... tetraedro*/
//  if(maxViz == 4)
//    adjTetra4(el         ,nodcon         ,nelcon
//             ,nnode      ,numel
//             ,maxNo      ,maxViz);
/*... hexaedro*/
//  else if(maxViz == 6)
//    adjHex8(el         ,nodcon         ,nelcon
//           ,nnode      ,numel
//           ,maxNo      ,maxViz, nFace);
  }
/*.....................................................................*/

/*for(int i=0;i<numel;i++)
  {
    printf("nel %4d ",i+1);
    for(int j=0;j<maxViz;j++)
      printf("%4d ",MAT2D(i,j,nelcon,maxViz));
    printf("\n");
  }
  printf("nFace = %4d\n", *nFace);*/

  HccaDealloc(m,nodcon,"nodcon",_AD_);
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 00/00/0000                                   *
* Data de modificaco : 07/04/2018                                   *
* ------------------------------------------------------------------*
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
void makeFaces(Memoria *m             , INT *RESTRICT el  
              , INT *RESTRICT nelcon  , short *RESTRICT nViz 
              , INT *RESTRICT cellFace
              , INT *RESTRICT fNode   , INT *RESTRICT owner
              , INT nnode             , INT numel
              , INT nFaces            , short maxNo   
              , short maxViz          , short ndm) 
{

/*...*/
  if (ndm == 2)
  {
    face2d(el    , nelcon
          , nViz , cellFace
          , fNode, owner
          , nnode, numel
          , maxNo, maxViz
          , nFaces);

    for (int i = 0; i<numel; i++)
    {
      printf("nel %4d faces", i + 1);
      for (int j = 0; j<maxViz; j++)
        printf("%4d ", MAT2D(i, j, cellFace, maxViz));
      printf("\n");
    }
    for (int i = 0; i<nFaces; i++)
    {
      printf("%4d ", i + 1);
      for (int j = 0; j<2; j++)
        printf(" %4d ", MAT2D(i, j, fNode, 2));
      for (int j = 0; j<2; j++)
        printf(" %4d ", MAT2D(i, j, owner, 2));
      printf("\n");
    }
  }
/*.....................................................................*/

/*...*/
  if (ndm == 3)
  {
    face3d(el    , nelcon
         , nViz  , cellFace
         , fNode , owner
         , nnode , numel
         , maxNo , 4
         , maxViz, nFaces);
  
/*  for (int i = 0; i<numel; i++)
    {
      printf("nel %4d faces ", i + 1);
      for (int j = 0; j<maxViz; j++)
        printf("%4d ", MAT2D(i, j, cellFace, maxViz));
      printf("\n");
    }
    for (int i = 0; i<nFaces; i++)
    {
      printf("%4d ", i + 1);
      for (int j = 0; j<4; j++)
        printf(" %4d ", MAT2D(i, j, fNode, 4));
      for (int j = 0; j<2; j++)
        printf(" %4d ", MAT2D(i, j, owner, 2));
      printf("\n");
    }
    exit(0);*/
  }
/*.....................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 07/04/2018                                   *
 *-------------------------------------------------------------------*
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
          ,short *nen      ,INT const nnode   ,INT const numel
          ,short const maxNo,short const maxViz,INT *nEdge)
{
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
  do
  {
    imiss = false;
    for(nel1=0;nel1 < numel;nel1++)
      {
      is1 = nen[nel1];
      for(j=0;j<is1;j++){
        if(MAT2D(nel1,j,nelcon,maxViz) == -1)
        {
/*... triangulo*/
          if(is1 == 3)
          {
            no1 = MAT2D(nel1,isnod3[j][0],el,maxNo)-1;
            no2 = MAT2D(nel1,isnod3[j][1],el,maxNo)-1;
          }
/*...................................................................*/

/*... quadrilatero*/
          else if(is1 == 4)
          {
            no1 = MAT2D(nel1,isnod4[j][0],el,maxNo)-1;
            no2 = MAT2D(nel1,isnod4[j][1],el,maxNo)-1;
          }
/*...................................................................*/
          if( nodcon[no1] == -1 || nodcon[no2] == -1)
          {
            nodcon[no1] = nodcon[no2] = nel1;
            imiss       = true;
          }
        }
      }
    }

/*...*/
    for(nel1=0;nel1 < numel;nel1++)
    {
      is1 = nen[nel1];
      for(j=0;j<is1;j++)
      {
        if(MAT2D(nel1,j,nelcon,maxViz) == -1)
        {
/*... triangulo*/
          if(is1 == 3){
           no1  = MAT2D(nel1,isnod3[j][0],el,maxNo)-1;
           no2  = MAT2D(nel1,isnod3[j][1],el,maxNo)-1;
          }
/*...................................................................*/

/*... quadrilateros*/
          else if(is1 == 4)
          {
           no1  = MAT2D(nel1,isnod4[j][0],el,maxNo)-1;
           no2  = MAT2D(nel1,isnod4[j][1],el,maxNo)-1;
          }
/*...................................................................*/
          nel2 = nodcon[no1];
          if(nel2 > -1)
          {
            if(nel2 == nodcon[no2] && nel2 != nel1)
            {
              is2 = nen[nel2];
              for(k=0;k<is2;k++)
              {
/*... triangulo*/
                if(is2 == 3)
                {
                  no21  = MAT2D(nel2,isnod3[k][0],el,maxNo)-1;
                  no22  = MAT2D(nel2,isnod3[k][1],el,maxNo)-1;
                }
/*...................................................................*/

/*... quadrilateros*/
                else if(is2 == 4)
                {
                  no21  = MAT2D(nel2,isnod4[k][0],el,maxNo)-1;
                  no22  = MAT2D(nel2,isnod4[k][1],el,maxNo)-1;
                }
/*...................................................................*/
                if( (no21 == no2) && (no22 == no1))
                {
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

    for(nel1=0;nel1 < numel;nel1++)
    {
      is1 = nen[nel1];
      for(j=0;j<is1;j++)
      {
        if(MAT2D(nel1,j,nelcon,maxViz) == -1)
        {
/*... triangulo*/
          if( is1 == 3)
          {
            no1  = MAT2D(nel1,isnod3[j][0],el,maxNo)-1;
            no2  = MAT2D(nel1,isnod3[j][1],el,maxNo)-1;
          }
/*...................................................................*/

/*... quadrilateros*/
          if( is1 == 4)
          {
            no1  = MAT2D(nel1,isnod4[j][0],el,maxNo)-1;
            no2  = MAT2D(nel1,isnod4[j][1],el,maxNo)-1;
          }
/*...................................................................*/
          if(nodcon[no1] == nodcon[no2] && nodcon[no1] == nel1)
          {
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

/*... numero comecando em 1*/
  for(nel1=0;nel1<numel;nel1++)
    for(j=0;j<maxViz;j++)
      MAT2D(nel1,j,nelcon,maxViz) += 1;
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 26/04/2018                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* adj3D : calcula os elementos vinzinhos                            *
* ------------------------------------------------------------------*
* Parametros de entrada:                                            *
* ------------------------------------------------------------------*
* el    -> conectividade nodal                                      *
* nodcon-> indefinido                                               *
* nelcon-> indefinido                                               *
* nViz  -> numero de nos vizinhos                                   *
* nnode -> numero de nos da malha                                   *
* numel -> numero de elmentos                                       *
* maxNo -> numero de nos por elementos maximo da malha              *
* maxViz-> numero vizinhos por elementos maximo da malha            *
* nFace -> indefinifo                                               *
* ------------------------------------------------------------------*
* Paramanetros de saida:                                            *
* ------------------------------------------------------------------*
* nelcon-> vizinhas dos elementos                                   *
* nFace -> numero de vizinhos                                       *
* ------------------------------------------------------------------*
* *******************************************************************/
void adj3d(INT *el          , INT *nodcon       , INT *nelcon
         , short *nViz      , INT const nnode   , INT const numel
         , short const maxNo, short const maxViz, INT *nFace)
{
  bool imiss;
  short l, k, j, is1, is2;
  INT  i, nel1, nel2, node[4];

  for (nel1 = 0; nel1<numel; nel1++)
    for (j = 0; j<maxViz; j++)
      MAT2D(nel1, j, nelcon, maxViz) = -1;

  for (i = 0; i < nnode; i++)
    nodcon[i] = -1;

  *nFace = 0;
  do
  {
    imiss = false;
    for (nel1 = 0; nel1 < numel; nel1++)
    {
      is1 = nViz[nel1];
      for (j = 0; j<is1; j++) {
        if (MAT2D(nel1, j, nelcon, maxViz) == -1)
        {
/*... tetraedro*/
          if (is1 == 4)
          {
            tetra4fNod(nel1, j
              , el, node
              , maxNo);
            if (nodcon[node[0]] == -1 || nodcon[node[1]] == -1
              || nodcon[node[2]] == -1) {
              nodcon[node[0]] = nel1;
              nodcon[node[1]] = nel1;
              nodcon[node[2]] = nel1;
              imiss = true;
            }
          }
/*...................................................................*/

/*... piramide5*/
          else if (is1 == 5)
          {
            pira5fNod(nel1, j, el , node, maxNo);
            pira5Aux1(j, nodcon, node, nel1, &imiss);
/*...................................................................*/
          }
/*...................................................................*/

/*... hexa*/
          else if (is1 == 6)
          {
            hexa8fNod(nel1   , j
                    , el     , node
                    , maxNo);
            if (nodcon[node[0]] == -1 || nodcon[node[1]] == -1
              || nodcon[node[2]] == -1 || nodcon[node[3]] == -1)
            {
              nodcon[node[0]] = nel1;
              nodcon[node[1]] = nel1;
              nodcon[node[2]] = nel1;
              nodcon[node[3]] = nel1;
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
    for (nel1 = 0; nel1 < numel; nel1++)
    {
      is1 = nViz[nel1];
      for (j = 0; j<is1; j++)
      {
        if (MAT2D(nel1, j, nelcon, maxViz) == -1)
        {
/*... tetra*/
          if (is1 == 4)
            tetra4fNod(nel1, j, el, node, maxNo);
/*...................................................................*/

/*... piramide5*/
          else if (is1 == 5)
            pira5fNod(nel1, j, el, node, maxNo);
/*...................................................................*/

/*... hexa*/
          else if (is1 == 6)
            hexa8fNod(nel1, j, el  , node, maxNo);
/*...................................................................*/
          nel2 = nodcon[node[0]];
          if (nel2 > -1)
          {
            is2 = nViz[nel2];
            for (k = 0; k<is2; k++)
            {
/*... triangulo*/
              if (is2 == 4)
              {
                if (nel2 == nodcon[node[1]] && nel2 == nodcon[node[2]] &&
                  nel2 != nel1) {
                  k = tetra4face(nel2, el, node, maxNo);
                  if (k == -1) {
                    printf("adjTetra4: Erro na vizinhaca!!!\n");
                    exit(EXIT_FAILURE);
                  }
                  MAT2D(nel2, k, nelcon, maxViz) = nel1;
                  MAT2D(nel1, j, nelcon, maxViz) = nel2;
                  nodcon[node[0]] = -1;
                  nodcon[node[1]] = -1;
                  nodcon[node[2]] = -1;
                  imiss = true;
                  (*nFace)++;
                }
              }
/*...................................................................*/

/*... piramide*/
              else if (is2 == 5){
                pira5Aux2(k    , j     
                         ,maxNo, maxViz
                         ,el   , nodcon, nelcon
                         ,node , nel1  , nel2
                         ,nFace, &imiss);
              }
/*...................................................................*/

/*... hexa*/
              else if (is2 == 6)
              {
                if (nel2 == nodcon[node[1]] && nel2 == nodcon[node[2]] &&
                    nel2 == nodcon[node[3]] && nel2 != nel1)
                {
                  l = hexa8face(nel2, el, node, maxNo);
                  if (l == -1)
                  {
                    printf("adjHex8: Erro na vizinhaca!!!\n");
                    exit(EXIT_FAILURE);
                  }
                  MAT2D(nel2, l, nelcon, maxViz) = nel1;
                  MAT2D(nel1, j, nelcon, maxViz) = nel2;
                  nodcon[node[0]] = -1;
                  nodcon[node[1]] = -1;
                  nodcon[node[2]] = -1;
                  nodcon[node[3]] = -1;
                  imiss = true;
                  (*nFace)++;
                }
/*...................................................................*/
              }    
/*.....................................................................*/
            }
          }
        }
      }
    }
/*...................................................................*/

    for (nel1 = 0; nel1 < numel; nel1++)
    {
      is1 = nViz[nel1];
      for (j = 0; j<is1; j++)
      {
        if (MAT2D(nel1, j, nelcon, maxViz) == -1)
        {
/*... tetra*/
          if (is1 == 4)
          {
            tetra4fNod(nel1, j, el, node, maxNo);
            nel2 = nodcon[node[0]];
            if (nel2 == nodcon[node[1]] && nel2 == nodcon[node[2]] &&
              nel2 == nel1) {
              MAT2D(nel1, j, nelcon, maxViz) = -2;
              nodcon[node[0]] = -1;
              nodcon[node[1]] = -1;
              nodcon[node[2]] = -1;
              imiss = true;
              (*nFace)++;
            }
          }
/*...................................................................*/

/*... piramide*/
          else if (is1 == 5)
          {
            pira5fNod(nel1, j, el, node, maxNo);
            nel2 = nodcon[node[0]];
            pira5Aux3(j     , maxViz
                    , nodcon, nelcon, node
                    , nel1  , nel2  , nFace
                    , &imiss);
           }
/*...................................................................*/

/*... hexa*/
          else if (is1 == 6)
          {
            hexa8fNod(nel1, j, el, node, maxNo);
            nel2 = nodcon[node[0]];
            if (nel2 == nodcon[node[1]] && nel2 == nodcon[node[2]] &&
                nel2 == nodcon[node[3]] && nel2 == nel1)
            {
              MAT2D(nel1, j, nelcon, maxViz) = -2;
              nodcon[node[0]] = -1;
              nodcon[node[1]] = -1;
              nodcon[node[2]] = -1;
              nodcon[node[3]] = -1;
              imiss = true;
              (*nFace)++;
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
  } while (imiss);

/*... numero comecando em 1*/
  for (nel1 = 0; nel1<numel; nel1++)
    for (j = 0; j<maxViz; j++)
      MAT2D(nel1, j, nelcon, maxViz) += 1;
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 26/04/2018                                   *
 *-------------------------------------------------------------------*
 * adjHex8 : calcula os elementos vinzinhos de uma malha de hexaedros*
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
 * nFace -> indefinido                                               *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * nelcon-> vizinhas dos elementos                                   *
 * nFace -> numero de faces                                          *  
 * ------------------------------------------------------------------*
 * *******************************************************************/
void adjHex8(INT *el          ,INT *nodcon         ,INT *nelcon
            ,INT const nnode  ,INT const numel
            ,short const maxNo,short const maxViz, INT *nFace)
{
  bool imiss;
  short k, j;
  INT nel1,nel2,node[4];

  for(nel1=0;nel1<numel;nel1++)
    for(j=0;j<maxViz;j++)
      MAT2D(nel1,j,nelcon,maxViz) = -1;

  for(nel1 = 0;nel1 < nnode;nel1++)
    nodcon[nel1] = -1;

  *nFace = 0;
  do{
    imiss = false;
/*...*/
    for(nel1=0;nel1<numel;nel1++)
    {
/*...*/
      for(j=0;j<6;j++)
      {
        if( MAT2D(nel1,j,nelcon,maxViz) == -1)
        {
          hexa8fNod(nel1        ,j
                   ,el          ,node
                   ,maxNo);
          if(  nodcon[node[0]] == -1 || nodcon[node[1]] == -1
            || nodcon[node[2]] == -1 || nodcon[node[3]] == -1)
          {
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
    for(nel1=0;nel1<numel;nel1++)
    {
      for(j=0;j<6;j++)
      {
        if( MAT2D(nel1,j,nelcon,maxViz) == -1)
        {
          hexa8fNod(nel1        ,j
                   ,el          ,node
                   ,maxNo);
/*...................................................................*/
          nel2 = nodcon[node[0]];
          if( nel2 > -1)
          {
             if(nel2 == nodcon[node[1]] && nel2 == nodcon[node[2]] &&
                nel2 == nodcon[node[3]] && nel2 != nel1)
              {
                k = hexa8face(nel2,el,node,maxNo);
                if( k == -1)
                {
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
                (*nFace)++;
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
    for(nel1=0;nel1<numel;nel1++)
    {
      for(j=0;j<6;j++)
      {
        if( MAT2D(nel1,j,nelcon,maxViz) == -1)
        {
          hexa8fNod(nel1        ,j
                   ,el          ,node
                   ,maxNo);
          nel2 = nodcon[node[0]];
          if(nel2 == nodcon[node[1]] && nel2 == nodcon[node[2]] &&
             nel2 == nodcon[node[3]] && nel2 == nel1)
          {
            MAT2D(nel1,j,nelcon,maxViz) = -2;
            nodcon[node[0]] = -1;
            nodcon[node[1]] = -1;
            nodcon[node[2]] = -1;
            nodcon[node[3]] = -1;
            imiss = true;
            (*nFace)++;
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

  }while(imiss);

/*... numero comecando em 1*/
  for(nel1=0;nel1<numel;nel1++)
    for(j=0;j<maxViz;j++)
      MAT2D(nel1,j,nelcon,maxViz) += 1;
/*...................................................................*/

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 07/04/2018                                   *
 *-------------------------------------------------------------------*
 * adjTetra4 : calcula os elementos vinzinhos de uma malha de        *
 * tetraedros                                                        *
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
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * nelcon-> vizinhas dos elementos                                   *
 * ------------------------------------------------------------------*
 * *******************************************************************/
void adjTetra4(INT *RESTRICT el    ,INT *RESTRICT nodcon
              ,INT *RESTRICT nelcon
              ,INT const nnode     ,INT const numel
              ,short const maxNo   ,short const maxViz){

  INT nel1,nel2,node[3];
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
      for(j=0;j<4;j++){
        if( MAT2D(nel1,j,nelcon,maxViz) == -1){
          tetra4fNod(nel1        ,j
                    ,el          ,node
                    ,maxNo);
          if(  nodcon[node[0]] == -1 || nodcon[node[1]] == -1
            || nodcon[node[2]] == -1){
            nodcon[node[0]]= nel1;
            nodcon[node[1]]= nel1;
            nodcon[node[2]]= nel1;
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
      for(j=0;j<4;j++){
        if( MAT2D(nel1,j,nelcon,maxViz) == -1){
          tetra4fNod(nel1        ,j
                    ,el          ,node
                    ,maxNo);
/*...................................................................*/
          nel2 = nodcon[node[0]];
          if( nel2 > -1){
              if(nel2 == nodcon[node[1]] && nel2 == nodcon[node[2]] &&
                 nel2 != nel1){
                 k = tetra4face(nel2,el,node,maxNo);
                 if( k == -1) {
                   printf("adjTetra4: Erro na vizinhaca!!!\n");
                   exit(EXIT_FAILURE);
                 }
                 MAT2D(nel2,k,nelcon,maxViz) = nel1;
                 MAT2D(nel1,j,nelcon,maxViz) = nel2;
                 nodcon[node[0]] = -1;
                 nodcon[node[1]] = -1;
                 nodcon[node[2]] = -1;
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
      for(j=0;j<4;j++){
        if( MAT2D(nel1,j,nelcon,maxViz) == -1){
          tetra4fNod(nel1        ,j
                   ,el          ,node
                   ,maxNo);
          nel2 = nodcon[node[0]];
          if(nel2 == nodcon[node[1]] && nel2 == nodcon[node[2]] &&
             nel2 == nel1){
            MAT2D(nel1,j,nelcon,maxViz) = -2;
            nodcon[node[0]] = -1;
            nodcon[node[1]] = -1;
            nodcon[node[2]] = -1;
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

/***********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * TETRA4FNOD -determina os nos da face j do elemento k                *
 * ------------------------------------------------------------------- *
 * Parametros de entrada:                                              *
 * ------------------------------------------------------------------- *
 * nEl      -> numero do elemento                                      *
 * face     -> numero da face do elemento                              *
 * el       -> conetividades nodal                                     *
 * nodeFace -> nao definido                                            *
 * maxNo    -> numero maximo de nos por elemento na malha              *
 * ------------------------------------------------------------------- *
 * Parametros de saida:                                                *
 * ------------------------------------------------------------------- *
 * nodeFace -> nos da face j (numeracao local)                         *
 **********************************************************************/
void tetra4fNod(INT const nEl     ,short const face
               ,INT *RESTRICT el   ,INT *RESTRICT nodeFace
               ,short const maxNo)
{
  short isNod[][3] = {{1,2,3}
                     ,{0,3,2}
                     ,{0,1,3}
                     ,{0,2,1}};
  short i;

  for(i=0;i<3;i++)
    nodeFace[i] = MAT2D(nEl,isNod[face][i],el,maxNo) - 1;


}
/**********************************************************************/

/***********************************************************************
 * Data de criacao    : 00/00/0000                                     *
 * Data de modificaco : 00/00/0000                                     *
 *-------------------------------------------------------------------- *
 * TETRA4FACE - determina a face do tetraedro k adjacente a face j     *
 * cujos nos estao armazenados em nodeFace                             *
 * ------------------------------------------------------------------- *
 * Parametros de entrada:                                              *
 * ------------------------------------------------------------------- *
 * nEl      -> numero do elemento                                      *
 * face     -> numero da face do elemento                              *
 * el       -> conetividades nodal                                     *
 * nodeFace -> nos da face j (numeracao local)                         *
 * maxNo    -> numero maximo de nos por elemento na malha              *
 * ------------------------------------------------------------------- *
 * Parametros de saida:                                                *
 * ------------------------------------------------------------------- *
 **********************************************************************/
short tetra4face(INT const nEl         ,INT *RESTRICT el
                ,INT *RESTRICT nodeFace,short const maxNo )
{
  INT no[3];
  short nFace,j;
/*possiveis numeracoes de faces*/
  short ind[][2] = {{1,2},{2,0},{0,1}};

  for(nFace=0;nFace<4;nFace++){
    tetra4fNod(nEl,nFace,el,no,maxNo);
    for(j=0;j<3;j++){
      if(no[0] == nodeFace[j]){
        if( no[1] == nodeFace[ind[j][1]] &&
            no[2] == nodeFace[ind[j][0]] )
        return nFace;
      }
    }
  }
  return -1;
}
/*********************************************************************/

/***********************************************************************
 * Data de criacao    : 00/00/0000                                     *
 * Data de modificaco : 00/00/0000                                     *
 * ------------------------------------------------------------------- *
 * HEXA8FNOD - determina os nos da face j do elemento k                *
 * ------------------------------------------------------------------- *
 * Parametros de entrada:                                              *
 * ------------------------------------------------------------------- *
 * nEl      -> numero do elemento                                      *
 * face     -> numero da face do elemento                              *
 * el       -> conetividades nodal                                     *
 * nodeFace -> nao definido                                            *
 * maxNo    -> numero maximo de nos por elemento na malha              *
 * ------------------------------------------------------------------- *
 * Parametros de saida:                                                *
 * ------------------------------------------------------------------- *
 * nodeFace -> nos da face j (numeracao local)                         *
 * normal apontando para fora                                          *
 * face1 1 4 3 2                                                       *
 * face2 5 6 7 8                                                       *
 * face3 1 2 6 5                                                       *
 * face4 2 3 7 6                                                       *
 * face5 3 4 8 7                                                       *
 * face6 4 1 5 8                                                       *
 **********************************************************************/
void hexa8fNod(INT const nEl     ,short const face
             ,INT *RESTRICT el   ,INT *RESTRICT nodeFace
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
/**********************************************************************/

/***********************************************************************
 * Data de criacao    : 00/00/0000                                     *
 * Data de modificaco : 00/00/0000                                     *
 * ------------------------------------------------------------------- *
 * HEXA8FACE - determina a face do hexaedro k adjacente a face j       *
 * cujos nos estao armazenados em nodeFace                             *
 * ------------------------------------------------------------------- *
 * Parametros de entrada:                                              *
 * ------------------------------------------------------------------- *
 * nEl      -> numero do elemento                                      *
 * el       -> conetividades nodal                                     *
 * nodeFace -> nos da face j (numeracao local)                         *
 * maxNo    -> numero maximo de nos por elemento na malha              *
 * ------------------------------------------------------------------- *
 * Parametros de saida:                                                *
 * ------------------------------------------------------------------- *
 **********************************************************************/
short hexa8face(INT const nEl         ,INT *RESTRICT el
               ,INT *RESTRICT nodeFace,short const maxNo )
{
  INT no[4];
  short nFace,j;
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
  return -1;
}
/*********************************************************************/

/***********************************************************************
* Data de criacao    : 26/04/2018                                     *
* Data de modificaco : 00/00/0000                                     *
* ------------------------------------------------------------------- *
* PRISM5FNOD - determina os nos da face j do elemento k               *
* ------------------------------------------------------------------- *
* Parametros de entrada:                                              *
* ------------------------------------------------------------------- *
* nEl      -> numero do elemento                                      *
* face     -> numero da face do elemento                              *
* el       -> conetividades nodal                                     *
* nodeFace -> nao definido                                            *
* maxNo    -> numero maximo de nos por elemento na malha              *
* ------------------------------------------------------------------- *
* Parametros de saida:                                                *
* ------------------------------------------------------------------- *
* nodeFace -> nos da face j (numeracao local)                         *
* normal apontando para fora                                          *
* face1 1 4 3 2                                                       *
* face2 3 6 9                                                         *
* face3 6 7 9                                                         *
* face4 7 8 9                                                         *
* face5 8 3 9                                                         *
**********************************************************************/
void pira5fNod(INT const nEl   , short const face
             , INT *RESTRICT el, INT *RESTRICT nodeFace
             , short const maxNo)
{
  short isNodPira[][4] = { { 0,3,2, 1}
                          ,{ 0,1,4, 0}
                          ,{ 1,2,4, 0}
                          ,{ 2,3,4, 0}
                          ,{ 3,0,4, 0}};

  short i;

  for (i = 0; i<4; i++){
      nodeFace[i] = MAT2D(nEl, isNodPira[face][i], el, maxNo) - 1;
  }

}
/**********************************************************************/

/***********************************************************************
* Data de criacao    : 00/00/0000                                     *
* Data de modificaco : 00/00/0000                                     *
* ------------------------------------------------------------------- *
* prism5face - determina a face do piramide k adjacente a face j      *
* cujos nos estao armazenados em nodeFace                             *
* ------------------------------------------------------------------- *
* Parametros de entrada:                                              *
* ------------------------------------------------------------------- *
* nEl      -> numero do elemento                                      *
* el       -> conetividades nodal                                     *
* nodeFace -> nos da face j (numeracao local)                         *
* maxNo    -> numero maximo de nos por elemento na malha              *
* ------------------------------------------------------------------- *
* Parametros de saida:                                                *
* ------------------------------------------------------------------- *
**********************************************************************/
short pira5face(INT const nEl, INT *RESTRICT el
              , INT *RESTRICT nodeFace, short const maxNo)
{
  INT no[4];
  short nFace, j;
/*possiveis numeracoes de faces*/
  short ind4[][3] = { { 1,2,3 },{ 2,3,0 },{ 3,0,1 },{ 0,1,2 } };
  short ind3[][2] = { { 1,2 }  ,{ 2,0 }  ,{ 0,1 } };

  for (nFace = 0; nFace<5; nFace++) 
  {
    pira5fNod(nEl, nFace, el, no, maxNo);
/*... face 1 - 4 nos*/
    if(!nFace)
    {
      for (j = 0; j<4; j++) 
      {
        if (no[0] == nodeFace[j]) 
        {
          if (no[1] == nodeFace[ind4[j][2]] &&
              no[2] == nodeFace[ind4[j][1]] &&
              no[3] == nodeFace[ind4[j][0]])
            return nFace;
        }
      }
    }
/*...................................................................*/

/*... face 2,3,4,5 - 3 nos*/
    else
    {
      for (j = 0; j<3; j++)
      {
        if (no[0] == nodeFace[j])
        {
          if (no[1] == nodeFace[ind3[j][1]] &&
              no[2] == nodeFace[ind3[j][0]])
            return nFace;
        }
      }
    }
/*...................................................................*/
  }
  return -1;
}
/*********************************************************************/

/*********************************************************************/
void pira5Aux1(short idFace, INT *nodcon, INT *node, INT nel1
             , bool *imiss)
{
/*...*/
  switch (idFace)
  {
/*...*/
    case 0:
      if (nodcon[node[0]] == -1 || nodcon[node[1]] == -1
        || nodcon[node[2]] == -1 || nodcon[node[3]] == -1)
      {
        nodcon[node[0]] = nel1;
        nodcon[node[1]] = nel1;
        nodcon[node[2]] = nel1;
        nodcon[node[3]] = nel1;
        *imiss = true;
      }
      break;
/*...................................................................*/

/*...*/
  case 1:
  case 2:
  case 3:
  case 4:
      if (nodcon[node[0]] == -1 || nodcon[node[1]] == -1
        || nodcon[node[2]] == -1)
      {
        nodcon[node[0]] = nel1;
        nodcon[node[1]] = nel1;
        nodcon[node[2]] = nel1;
        *imiss = true;
      }
    break;
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************/
void pira5Aux2(short k    , short j  
             , short maxNo, short maxViz
             , INT *el    , INT *nodcon, INT *nelcon
             , INT *node  , INT nel1   , INT nel2
             , INT *nFace , bool *imiss)
{

  short l;
/*...*/
  switch (k)
  {
/*...*/
    case 0:
      if (nel2 == nodcon[node[1]]
        && nel2 == nodcon[node[2]]
        && nel2 == nodcon[node[3]]
        && nel2 != nel1)
      {
        l = pira5face(nel2, el, node, maxNo);
        if (l == -1)
        {
          printf("adjPira5: Erro na vizinhaca!!!\n");
          exit(EXIT_FAILURE);
        }
        MAT2D(nel2, l, nelcon, maxViz) = nel1;
        MAT2D(nel1, j, nelcon, maxViz) = nel2;
        nodcon[node[0]] = -1;
        nodcon[node[1]] = -1;
        nodcon[node[2]] = -1;
        nodcon[node[3]] = -1;
        (*nFace)++;
        *imiss = true;
      }
      break;
/*...................................................................*/

/*...*/
    case 1:
    case 2:
    case 3:
    case 4:
      if  (nel2 == nodcon[node[1]]
        && nel2 == nodcon[node[2]]
        && nel2 != nel1)
      {
        l = pira5face(nel2, el, node, maxNo);
        if (k == -1)
        {
          printf("adjPira5: Erro na vizinhaca!!!\n");
          exit(EXIT_FAILURE);
        }
        MAT2D(nel2, l, nelcon, maxViz) = nel1;
        MAT2D(nel1, j, nelcon, maxViz) = nel2;
        nodcon[node[0]] = -1;
        nodcon[node[1]] = -1;
        nodcon[node[2]] = -1;
        (*nFace)++;
        *imiss = true;
    }
    break;
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/*********************************************************************/
void pira5Aux3(short j    , short maxViz
             , INT *nodcon, INT *nelcon, INT *node
             , INT nel1   , INT nel2   , INT *nFace
             , bool *imiss)
{
  switch (j)
  {
/*...*/
  case 0:
    if (nel2 == nodcon[node[1]]
      && nel2 == nodcon[node[2]]
      && nel2 == nodcon[node[3]]
      && nel2 == nel1)
    {
      MAT2D(nel1, j, nelcon, maxViz) = -2;
      nodcon[node[0]] = -1;
      nodcon[node[1]] = -1;
      nodcon[node[2]] = -1;
      nodcon[node[3]] = -1;
      (*nFace)++;
      *imiss = true;
    }
    break;
/*...................................................................*/

/*...*/
  case 1:
  case 2:
  case 3:
  case 4:
    if (nel2 == nodcon[node[1]]
      && nel2 == nodcon[node[2]]
      && nel2 == nel1)
    {
      MAT2D(nel1, j, nelcon, maxViz) = -2;
      nodcon[node[0]] = -1;
      nodcon[node[1]] = -1;
      nodcon[node[2]] = -1;
      (*nFace)++;
      *imiss = true;
    }
    break;
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/