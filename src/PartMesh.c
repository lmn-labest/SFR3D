#include<PartMesh.h>
/********************************************************************* 
 * PARTMESH :                                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m       -> vetor de memoria                                       * 
 * x       -> coordenadas                                            * 
 * el      -> conectividade dos elementos                            * 
 * nen     -> numero de nos por celulas                              * 
 * nNode   -> numero do nos da malha                                 * 
 * nEl     -> numero de elementos da malha                           * 
 * pMesh   -> varaives de particionamento                            * 
 * ndm     -> numero de dimensoes                                    * 
 * maxNo   -> numero maximo de nos por celula                        * 
 * nDiv    -> numero de divisoes                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void partMesh(Memoria *m      
             ,DOUBLE *RESTRICT x   ,INT *RESTRICT el
             ,short  *RESTRICT nen
             ,INT const nNode      ,INT const nEl   
             ,PartMesh *pMesh  
             ,short const ndm      ,short const maxNo     
             ,short const nDiv     ){
#ifdef _METIS_
  INT n,eSize;  
  INT *eptr=NULL,*eind=NULL;
  int *vsize   = NULL; 
  int *vwgt    = NULL; 
  double *tpwgts  = NULL; 
  int *options = NULL; 
  int ne,nn,nCommon,nParts,edgeCut;
#endif


/*... alocando variaveis do particionamento*/
  HccaAlloc(INT,m,pMesh->np     ,nNode ,"nNodeP" ,false);
  HccaAlloc(INT,m,pMesh->ep     ,nEl   ,"elP"    ,false);
/*...................................................................*/

/*... METIS5*/
#ifdef _METIS_
  n = nEl + 1;
/*... grafo CSR da malha*/
  HccaAlloc(INT,m,eptr          ,n     ,"eptr" ,false);
  eSize = meshToCsrIa(nen,eptr,nEl);  

  HccaAlloc(INT,m,eind          ,eSize ,"eind" ,false);
  meshToCsrJa(el  ,nen
             ,eptr,eind
             ,nEl ,maxNo);
/*...................................................................*/ 

/*... metis*/
/*... triangulos*/
  if( maxNo == 3 ) nCommon = 2;
/*... quadrilateros*/
  else if( maxNo == 4 && ndm == 2 ) nCommon = 2;
/*... tetraedros*/
  else if( maxNo == 4 && ndm == 3) nCommon = 3;
/*... hexaedros*/
  else if( maxNo == 8 ) nCommon = 4;

  ne      = nEl;
  nn      = nNode;
  nParts  = nDiv;
  printf("Metis-5 ...\n");
  METIS_PartMeshDual(&ne      ,&nn      ,eptr    ,eind
                    ,vwgt     ,vsize    ,&nCommon,&nParts
                    ,tpwgts   ,options  ,&edgeCut 
                    ,pMesh->ep,pMesh->np);         
  printf("Metis-5.\n");
/*...................................................................*/ 
    
/*
  printf("ep\n");
  for(i=0;i<nEl;i++)
    printf("%3d %3d\n",i+1,ep[i]);
*/
/*... desalocando o grafo CSR da malha*/
  HccaDealloc(m,eind  ,"eind",false);
  HccaDealloc(m,eptr  ,"eptr",false);
/*...................................................................*/ 

#else
/*...*/
  printf("divCoorXY ...\n");
  divCoorXY(x        ,el
           ,nen
           ,nNode    ,nEl
           ,pMesh->np,pMesh->ep
           ,ndm      ,maxNo 
           ,nDiv     ,true);
  printf("divCoorXY.\n");
/*...................................................................*/
#endif

}
/*********************************************************************/

/********************************************************************* 
 * DIVCOORXY : Divisao da malha pelas coordenadas no plano XY        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * coor    -> coordenadas                                            * 
 * el      -> conectividade dos elementos                            * 
 * nen     -> numero de nos por celulas                              * 
 * nNode   -> numero do nos da malha                                 * 
 * nEl     -> numero de elementos da malha                           * 
 * np      -> nao definido                                           * 
 * ep      -> nao definido                                           * 
 * ndm     -> numero de dimensoes                                    * 
 * maxNo   -> numero maximo de nos por celula                        * 
 * nDiv    -> numero de divisoes                                     * 
 * fC      -> true numeracao iniciando em 0                          * 
 *            false numeracao iniciando em 1                         * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * np      -> divisao dos nos                                        * 
 * ep      -> divisao dos elementos                                  * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void divCoorXY(DOUBLE *RESTRICT coor,INT *RESTRICT el
              ,short  *RESTRICT nen                          
              ,INT const nNode      ,INT const nEl   
              ,INT *RESTRICT np     ,INT *RESTRICT ep   
              ,short const ndm      ,short const maxNo
              ,short const nDiv     ,bool const fC){

  DOUBLE yMax,yMin,xMax,xMin;
  DOUBLE lx,ly,dx,dy,x,y;
  int hx,hy;
  short nx=1,ny=1;   
  INT i;
  short j;

  for(i=0;i<nEl;i++)
    ep[i]=nEl;
  
  xMax = xMin = MAT2D(0,0,coor,ndm);
  yMax = yMin = MAT2D(0,1,coor,ndm);
 
  for(i=1;i<nNode;i++){
    xMax = max(MAT2D(i,0,coor,ndm),xMax);
    yMax = max(MAT2D(i,1,coor,ndm),yMax);
    xMin = min(MAT2D(i,0,coor,ndm),xMin);
    yMin = min(MAT2D(i,1,coor,ndm),yMin);
  }

  if( nDiv == 2){
    nx = 1;
    ny = 2;
  }
  
  else if( nDiv == 4){
    nx = 2;
    ny = 2;
  }
  
  else if( nDiv == 6){
    nx = 3;
    ny = 2;
  }
  
  else if( nDiv == 8){
    nx = 4;
    ny = 2;
  }
  
  else if( nDiv == 10){
    nx = 5;
    ny = 2;
  }
  
  else if( nDiv == 12){
    nx = 4;
    ny = 3;
  }
  
  else if( nDiv == 16){
    nx = 4;
    ny = 4;
  }
  
  lx = (xMax - xMin);
  ly = (yMax - yMin);

  for(i=0;i<nNode;i++){
    x  = MAT2D(i,0,coor,ndm);
    y  = MAT2D(i,1,coor,ndm);
    dx = (x - xMin)/lx;
    dy = (y - yMin)/ly;
    hx = (INT) (((DOUBLE) nx)*dx); 
    hy = (INT) (((DOUBLE) ny)*dy);
    hx = min(hx+1,nx); 
    hy = min(hy+1,ny); 
    np[i] = (hy - 1)*nx+hx;
  }

  for(i=0;i<nEl;i++){
    for(j=0;j<nen[i];j++){
      ep[i] = min(ep[i],np[MAT2D(i,j,el,maxNo)-1]);
    }
  } 
  

  if(fC){
    for(i=0;i<nNode;i++)
      np[i]--;
  
    for(i=0;i<nEl;i++)
      ep[i]--;
  }

/*  
  printf("np\n");
  for(i=0;i<nNode;i++)
    printf("%3d %3d\n",i+1,np[i]);
  
  printf("ep\n");
  for(i=0;i<nEl;i++)
    printf("%3d %3d\n",i+1,ep[i]);

  ep[ 0] = 0;
  ep[ 1] = 0;
  ep[ 2] = 0;
  ep[ 3] = 0;
  ep[ 4] = 0;
  ep[ 5] = 0;
  ep[ 6] = 0;
  ep[ 7] = 0;
  ep[ 8] = 0;
  ep[ 9] = 1;
  ep[10] = 1;
  ep[11] = 1;
  ep[12] = 2;
  ep[13] = 2;
  ep[14] = 2;
  ep[15] = 2;
*/
}
/*********************************************************************/
  
/********************************************************************* 
 * GETNUMEBERLOCALMESH : obtem o numero de elementos e nos da malha  * 
 * local                                                             * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ep       -> numero de particao do elemento                        * 
 * fNode    -> vetor axiliar                                         * 
 * fEq      -> vetor axiliar                                         * 
 * el       -> conectividade dos elementos global                    * 
 * nen      -> numero de nos por elementos global                    * 
 * nNode    -> numero de nos global                                  * 
 * nEl      -> numero de elementos gobal                             * 
 * maxNo    -> numero maximo de no por elementos                     * 
 * maxViz   -> numero maximo de vizinhos                             * 
 * rank     -> numero do processo da malha                           * 
 * numelNov -> nao definido                                          * 
 * numelOv  -> nao definido                                          * 
 * nNodeNov -> nao definido                                          * 
 * nNodeOv  -> nao definido                                          * 
 * nno1     -> nao definido                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * numelNov -> elementos em sem sobreposicao                         * 
 * numelOv  -> elementos em sobreposicao                             * 
 * nNodeNov -> numero de nos local                                   * 
 * nNodeOv  -> numero de nos local em elementos em sobre posicao     * 
 * nno1     -> numero de nos pertecentes apenas a minha particao     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void getNumberLocalMesh(INT *RESTRICT ep        ,INT *RESTRICT np                     
                       ,bool *RESTRICT fNode    ,bool *RESTRICT fEp
                       ,INT *RESTRICT el        ,INT *RESTRICT nelcon 
                       ,short *RESTRICT nen     ,short *RESTRICT nFace
                       ,INT const nNode         ,INT const nEl   
                       ,short const maxNo       ,short const maxViz
                       ,INT *numelNov           ,INT *numelOv       
                       ,INT *nNodeNov           ,INT *nNodeOv 
                       ,INT *nno1                         
                       ,short const rank )
{
  INT i,j,k,no,nelViz;

  *numelNov = 0;
  *numelOv  = 0;
  *nNodeNov = 0;
  *nNodeOv  = 0;
  *nno1     = 0;
  
  for(i=0;i<nEl;i++){
    fEp[i] = true;
  }

/*... numero de elementos*/
  for(i=0;i<nEl;i++){
    if( ep[i] == rank){
      (*numelNov)++;
      for(j=0;j<nFace[i];j++){
        nelViz = MAT2D(i,j,nelcon,maxViz) - 1;
/*... elementos vizinhos*/
        if( nelViz != -2 ) 
          if( ep[nelViz] != rank && fEp[nelViz]){
            (*numelOv)++;
            fEp[nelViz] = false;
          }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
  
  for(i=0;i<nNode;i++){
    fNode[i] = true;
    if( np[i] == rank) (*nno1)++;
  }

/*... numero de nos sem sobeprosicao*/
  for(i=0;i<nEl;i++){
    if( ep[i] == rank){
      for(k=0;k<nen[i];k++){
        no = MAT2D(i,k,el,maxNo) - 1;
        if(fNode[no]){
          (*nNodeNov)++;
          fNode[no] = false;
        }
      }
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... numero de nos com sobreposicao*/
  for(i=0;i<nEl;i++){
    if( ep[i] == rank){
      for(j=0;j<nFace[i];j++){
        nelViz = MAT2D(i,j,nelcon,maxViz) - 1;
/*... elementos vizinhos*/
        if( nelViz != -2 ) 
          if( ep[nelViz] != rank)
            for(k=0;k<nen[nelViz];k++){
              no = MAT2D(nelViz,k,el,maxNo) - 1;
              if(fNode[no]){
                (*nNodeOv)++;
                fNode[no] = false;
              }
          }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... contando numero total de nos locais*/
//  for(i=0;i<nNode;i++)
//    if(fNode[i])
//      (*lNnode)++;  
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * GETMAPVIZ  :  gera o mapa de vizinho por particao(face)           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ep         -> numero de particao                                  * 
 * elLG       -> mapa local->global de elementos                     * 
 * vizPart    -> nao definido                                        * 
 * lNel       -> numero de elementos locais                          * 
 * numelNov   -> numero de elementos sem sobreposicao                * 
 * rank       -> numero do processo da malha                         * 
 * nPrcs      -> numero de processos Total                           * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * vizPart    -> particao vizinhas a um particao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 unsigned short getMapViz(INT *RESTRICT ep,INT *RESTRICT  elLG
                ,short *RESTRICT vizPart  
                ,INT const  lNel          ,INT const numelNov
                ,short const rank         ,short const nPrcs)
{
  
  INT i,j,kk;
  unsigned short k = 0;  

//for(i=0;i<nPrcs;i++)
//  vizPart[i] = rank;

  for(j=0;j<nPrcs;j++){
    if( rank != j)
      for(i=numelNov;i<lNel;i++){
        kk = elLG[i];
        if(ep[kk] == j){
          vizPart[k] = ep[kk];
          k++;
          break;
        }
      }
  }
/*... se houver mais de uma particao vizinha ordenda em ordem 
      cresente*/
  if( k > 1)
   sBubblesort(vizPart,k);  

  return k;  

}
/*********************************************************************/ 

/********************************************************************* 
 * GETMAPVIZNO:  gera o mapa de vizinho por particao(no)             * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ep         -> numero de particao                                  * 
 * elLG       -> mapa local->global de elementos                     * 
 * vizPart    -> nao definido                                        * 
 * lNel       -> numero de elementos locais                          * 
 * numelNov   -> numero de elementos sem sobreposicao                * 
 * rank       -> numero do processo da malha                         * 
 * nPrcs      -> numero de processos Total                           * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * vizPart    -> particao vizinhas a um particao                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 unsigned short getMapVizNo(INT *RESTRICT ep,INT *RESTRICT  noLG
                  ,short *RESTRICT vizPart  ,bool *RESTRICT contViz
                  ,INT *RESTRICT nincid     ,INT *RESTRICT incid 
                  ,INT const  nNodeNov      ,INT const maxGrade              
                  ,short const rank         ,short const nPrcs)
{
  
  INT i,ii,el,no,partd;
  unsigned short k = 0;  

  for(i=0;i<nPrcs;i++)
    contViz[i] = true;

  for(i=0;i<nNodeNov;i++){
/*... loop na particoes vizinha*/
    no        = noLG[i];
    for(ii=0;ii<nincid[no];ii++){
      el = MAT2D(no,ii,incid,maxGrade);
      partd = ep[el]; 
      if(partd != rank && contViz[partd]){
        vizPart[k++]  = partd;
        contViz[partd] = false;
      }
    }
  }
/*... se houver mais de uma particao vizinha ordenda em ordem 
      cresente*/
  if( k > 1)
   sBubblesort(vizPart,k);  

  return k;  
}
/*********************************************************************/

/********************************************************************* 
 * GETMAPINTERFACEEL: gera o mapa de interface de elementos entre    *
 * as particoes                                                      * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ep       -> divisao dos elementos                                 * 
 * nelcon   -> vizinhos dos elementos                                * 
 * nFace    -> numero de faces por celulas                           * 
 * elLG     -> mapa local->global de elementos                       * 
 * fEq      -> vetor auxiliar                                        * 
 * nRcvs    -> nao definido                                          * 
 * nSends   -> nao definido                                          * 
 * vizPart  -> particoes vizinhas                                    * 
 * fMap     -> nao definido                                          * 
 * lNel     -> numero total de elementos da particao                 * 
 * numelNov -> numero elementos da particao sem sobreposicoes        * 
 * rank     -> numero da particao                                    * 
 * nPrcs    -> numero de processos total                             * 
 * nVizPart -> numero de particoes vizinhas                          * 
 * maxViz   -> numero de maximo de celulas vizinhas                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * nRcvs    -> numero de elementos no buffer de recebimento          * 
 * nSends   -> numero de elementos no buffer de envio                * 
 * fMap     -> buffer de comunicao ( numeracao por celulas)          * 
 *-------------------------------------------------------------------* 
 * OBS: fMap e gerado considerando a numeracao local das celulas     * 
 * ordenado primeiro por vizinho e depois numeracao global           * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void  getMapInterfaceEl(INT *RESTRICT ep
                       ,INT *RESTRICT nelcon    ,short *RESTRICT nFace
                       ,INT *RESTRICT elLG      ,bool *RESTRICT fEl
                       ,INT *nRcvs              ,INT *nSends
                       ,INT *RESTRICT iaRcvs    ,INT *RESTRICT iaSends
                       ,short *RESTRICT vizPart ,INT *RESTRICT fMap 
                       ,INT const lNel          ,INT const numelNov
                       ,short const rank        ,short const nPrcs
                       ,short const nVizPart    ,short const maxViz)
{

  INT i,el,nelViz,kk=0;
  short partId,k,j;

  for(i=0;i<lNel;i++){
    fMap[i]    = 0;
  }

/*... elementos no buffer de recebimento*/
  for(j=0;j<nVizPart;j++){
    iaRcvs[j] =kk;
    partId    = vizPart[j];
    for(i=numelNov;i<lNel;i++){
/*... loop na particoes vizinha*/
      el        = elLG[i];
      if(ep[el] == partId)
        fMap[kk++] = i;
    }
    v2Bubblesort(&fMap[iaRcvs[j]],elLG,kk-iaRcvs[j]);
/*...................................................................*/
  }
  *nRcvs = kk;  
  iaRcvs[nVizPart] = kk;

/*... elementos no buffer de envio*/
  for(j=0;j<nVizPart;j++){
    for(i=0;i<lNel;i++)
      fEl[i] = true;  
/*... loop na particoes vizinha*/
    iaSends[j] = kk;
    for(i=0;i<numelNov;i++){
/*... loop na elementos vizinho*/
      for(k=0;k<nFace[i];k++){
        el     = elLG[i];
        nelViz = MAT2D(el,k,nelcon,maxViz) - 1;
        if(nelViz != -2){
          partId = vizPart[j];
          if(ep[nelViz] == partId && fEl[i]){
            fMap[kk++]  = i;
            fEl[i]      = false;  
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
 
  *nSends = kk - *nRcvs;  
  iaSends[nVizPart] = kk;
  
}                  
/*********************************************************************/

/********************************************************************* 
 * GETMAPINTERFACENO: gera o mapa de interface de nos entre as       *
 * particoes                                                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ep       -> divisao dos elementos                                 * 
 * noLG     -> mapa local->global de elementos                       * 
 * fNo      -> vetor auxiliar                                        * 
 * nRcvs    -> nao definido                                          * 
 * vizPart  -> particoes vizinhas                                    * 
 * fMap     -> nao definido                                          * 
 * nNodeNov -> numero de nos sem sobreposicao                        * 
 * rank     -> numero da particao                                    * 
 * nPrcs    -> numero de processos total                             * 
 * nVizPart -> numero de particoes vizinhas                          * 
 * maxViz   -> numero de maximo de celulas vizinhas                  * 
 * maxGrade - numero de incidencias maximo na malha                  *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * nRcvs    -> numero de elementos no buffer de comunicacao          * 
 * fMap     -> buffer de comunicao ( numeracao por celulas)          * 
 *-------------------------------------------------------------------* 
 * OBS: fMap e gerado considerando a numeracao local dos nos e que   * 
 * a comunicação nodal e bidirecional e ordenado por vizinho         *
 * e depois por numeracao global                                     * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void  getMapInterfaceNo(INT *RESTRICT ep
                       ,INT *RESTRICT noLG      ,bool *RESTRICT fNod
                       ,INT *nComNo             ,INT *RESTRICT iaComNo
                       ,INT *RESTRICT nincid    ,INT *RESTRICT incid      
                       ,short *RESTRICT vizPart ,INT *RESTRICT fMap 
                       ,INT const nNodeNov       
                       ,short const rank        ,short const nPrcs
                       ,short const nVizPart    ,short const maxGrade)
{

  INT i,ii,no,el,kk=0;
  short partId,j;
  

/*... elementos no buffer de recebimento*/
  for(j=0;j<nVizPart;j++){
    iaComNo[j] =kk;
    partId     = vizPart[j];
    for(i=0;i<nNodeNov;i++)
      fNod[i] = true;  
/*...*/    
    for(i=0;i<nNodeNov;i++){
/*... loop na particoes vizinha*/
      no        = noLG[i];
      for(ii=0;ii<nincid[no];ii++){
        el = MAT2D(no,ii,incid,maxGrade);
        if(ep[el] == partId && fNod[i]){
          fMap[kk++] = i;
          fNod[i]    = false;  
        }
      }
    }
    v2Bubblesort(&fMap[iaComNo[j]],noLG,kk-iaComNo[j]);
/*...................................................................*/
  }
  *nComNo = kk;  
  iaComNo[nVizPart] = kk;

//  for(kk=0;kk<*nComNo;kk++)   
//    printf("%4d %4d %4d\n",kk+1,fMap[kk]+1,noLG[fMap[kk]]);

}                  
/*********************************************************************/

/********************************************************************* 
 * GETMAPELM : gera o mapa Local->global de elementos                * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ep       -> numero de particao                                    * 
 * elLG     -> nao definido                                          * 
 * elGL     -> nao definido                                          * 
 * fEq      -> vetor auxiliar                                        * 
 * el       -> conectividade dos elementos global                    * 
 * nelcon   -> vizinhos dos elementos                                * 
 * nFace    -> numero de faces por celulas                           * 
 * nEl      -> numero de elementos locais                            * 
 * lNel     -> numero de elementos locais                            * 
 * maxViz   -> numero maximo de vizinhos                             * 
 * rank     -> numero do processo da malha                           * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * elLG     -> mapa local->global de elementos                       * 
 * elGL     -> mapa gobal->local  de elementos                       * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void getMapElm(INT *RESTRICT ep
              ,INT *RESTRICT elLG      ,INT *RESTRICT elGL
              ,bool *RESTRICT fEp
              ,INT *RESTRICT el        ,INT *RESTRICT nelcon 
              ,short *RESTRICT nFace
              ,INT const nEl           ,INT const lNel
              ,short const maxViz      ,short const rank )
{
  INT i,j,k=0,nelViz;

/*... elementos da propria particao*/
  for(i=0;i<nEl;i++){
    fEp[i] = true;
    if(ep[i] == rank)
      elLG[k++] = i;
  }
/*...................................................................*/

/*... numero de elementos sem sobreposicao*/
  for(i=0;i<nEl;i++){
    if( ep[i] == rank){
      for(j=0;j<nFace[i];j++){
        nelViz = MAT2D(i,j,nelcon,maxViz) - 1;
/*... elementos vizinhos*/
        if( nelViz != -2 ) 
          if( ep[nelViz] != rank && fEp[nelViz]){
            elLG[k++]   = nelViz;
            fEp[nelViz] = false;
          }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/

/*... mapa global-> local */                
  for(i=0;i<nEl;i++){
    elGL[i] = 0;
  }
 
  for(i=0;i<lNel;i++){
    elGL[elLG[i]] = i;
  }
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * GETMAPNODE: gera o mapa Local->global de nos                      * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * ep       -> numero de particao do elemento                        * 
 * noLG     -> nao definido                                          * 
 * noGL     -> nao definido                                          * 
 * fNode    -> vetor axiliar                                         * 
 * el       -> conectividade dos elementos global                    * 
 * nelcon   -> vizinhos dos elementos                                * 
 * nen     -> numero de nos por celulas                              * 
 * nFace    -> numero de faces por celulas                           * 
 * nEl      -> numero de elementos gobal                             * 
 * nNode    -> numero de nos global                                  * 
 * maxNo   -> numero de nos por celula maximo da malha               * 
 * maxViz  -> numero vizinhos por celula maximo da malha             * 
 * rank     -> numero do processo da malha                           * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * noLG     -> mapa local->global de nos                             * 
 * noGL     -> mapa global->local de nos                             * 
 *-------------------------------------------------------------------* 
 * OBS: primiros os nos pertecente a propria particao, depois os nos * 
 * de outras particoes e nos em elemento em sobreposicao             * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void getMapNode(INT *RESTRICT ep        ,INT *RESTRICT np                     
               ,INT *RESTRICT noLG      ,INT *RESTRICT noGL
               ,bool *RESTRICT fNode
               ,INT *RESTRICT el        ,INT *RESTRICT nelcon 
               ,short *RESTRICT nen     ,short *RESTRICT nFace
               ,INT const nNode         ,INT const nEl           
               ,short const maxNo       ,short const maxViz
               ,INT  const lNode        ,short const rank )
{
  INT i,k,j,no,nelViz,kk=0;
          

  for(i=0;i<nNode;i++)
    fNode[i]  = true;

/*... mapa local-> global */                

/*... loop no elementos sem sobreposicao da particao local*/
/*... nos da propria particapo*/
  for(i=0;i<nEl;i++)
    if( ep[i] == rank)
      for(k=0;k<nen[i];k++){
        no = MAT2D(i,k,el,maxNo) - 1;
        if(fNode[no] && np[no] == rank){
          noLG[kk++] = no;
          fNode[no]  = false;
        }
      }
/*... nos da particapo vizinha*/
  for(i=0;i<nEl;i++)
    if( ep[i] == rank)
      for(k=0;k<nen[i];k++){
        no = MAT2D(i,k,el,maxNo) - 1;
        if(fNode[no]){
          noLG[kk++] = no;
          fNode[no]  = false;
        }
      }
/*...................................................................*/

/*... loop no elementos sobrepostos*/      
  for(i=0;i<nEl;i++)
    if( ep[i] == rank)
      for(j=0;j<nFace[i];j++){
        nelViz = MAT2D(i,j,nelcon,maxViz) - 1;
/*... elementos vizinhos*/
        if( nelViz != -2 && ep[nelViz] != rank) 
          for(k=0;k<nen[nelViz];k++){
            no = MAT2D(nelViz,k,el,maxNo) - 1;
            if(fNode[no]){
              noLG[kk++] = no;
              fNode[no]  = false;
            }
          }
      }  
/*...................................................................*/

/*... mapa global-> local */                
  for(i=0;i<nNode;i++){
    noGL[i] = 0;
  }
 
  for(i=0;i<lNode;i++){
    noGL[noLG[i]] = i;
  }
/*...................................................................*/

}
/*********************************************************************/
        
/********************************************************************* 
 * GETLOCALEL: obtem as conectividades local dos elementos           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * el      -> conectividade dos elementos global                     * 
 * lEl     -> nao definido                                           * 
 * elLG    -> mapa local->global de elementos                        * 
 * noGL    -> mapa global->local de nos                              * 
 * nen     -> numero de nos por celulas                              * 
 * nLel    -> numero de elementos local                              * 
 * maxNo   -> numero de nos maximo por elemento                      * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lEl     -> conectividade local dos elementos                      * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void getLocalEl(INT *RESTRICT el   ,INT *RESTRICT lEl
               ,INT *RESTRICT elLG ,INT *RESTRICT noGL
               ,short *RESTRICT nen   
               ,INT const nLel     ,short const maxNo)
{
  INT i,elG,noG;
  short j;
  
  for(i=0;i<nLel;i++){
    elG     = elLG[i];
    for(j=0;j<nen[elG];j++){ 
      noG   =  MAT2D(elG,j,el,maxNo)-1;
      MAT2D(i,j,lEl,maxNo) = noGL[noG]+1; 
    }
  }

}
/*********************************************************************/ 

/********************************************************************* 
 * GETLOCALADJ: obtem as adjacencia local dos elementos              * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * adj     -> vizinhaca gobal dos elementos                          * 
 * lAdj    -> nao definido                                           * 
 * elLG    -> mapa local->global de elementos                        * 
 * elGL    -> mapa global->local de elementos                        * 
 * nViz    -> numero de viznhos por celula                           * 
 * nLel    -> numero de elementos local                              * 
 * maxViz  -> numero de viznhos maximo por elemento                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lEl     -> conectividade local dos elementos                      * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void getLocalAdj(INT *RESTRICT adj  ,INT *RESTRICT lAdj
               ,INT *RESTRICT elLG  ,INT *RESTRICT elGL
               ,short *RESTRICT nViz  
               ,INT const nLel     ,short const maxViz)
{
  INT i,elG,adjG;
  short j;
  
  for(i=0;i<nLel;i++){
    elG     = elLG[i];
    for(j=0;j<nViz[elG];j++){ 
      adjG   =  MAT2D(elG,j,adj,maxViz)-1;
      if(adjG > -1)
        MAT2D(i,j,lAdj,maxViz) = elGL[adjG]+1;
      else 
        MAT2D(i,j,lAdj,maxViz) = adjG+1;
    }
  }

}
/*********************************************************************/ 

/********************************************************************* 
 * DGETLOCALV: obtem as arranjo local de um variavel global (double) * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * vG      -> arranjo global                                         * 
 * vL      -> nao definido                                           * 
 * mapLG   -> mapa local->global                                     * 
 * noGL    -> mapa global->local de nos                              * 
 * nLin    -> numero de linhas                                       * 
 * nCol    -> numero de colunas                                      * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * vL      -> arranjo local                                          * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void dGetLocalV(DOUBLE *RESTRICT vG  ,DOUBLE *RESTRICT vL
               ,INT *RESTRICT mapLG
               ,INT const nLin     ,INT const nCol)
{
  INT i,j,nG;
  
  if(nCol == 1){
    for(i=0;i<nLin;i++){
      nG      = mapLG[i];
      vL[i]   = vG[nG]; 
    }
  }
  else{
    for(i=0;i<nLin;i++){
      nG      = mapLG[i];
      for(j=0;j<nCol;j++)  
        MAT2D(i,j,vL,nCol) = MAT2D(nG,j,vG,nCol); 
    }
  }

}
/*********************************************************************/ 

/********************************************************************* 
 * SGETLOCALV: obtem as arranjo local de um variavel global (short)  * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * vG      -> arranjo global                                         * 
 * vL      -> nao definido                                           * 
 * mapLG   -> mapa local->global                                     * 
 * noGL    -> mapa global->local de nos                              * 
 * nLin    -> numero de linhas                                       * 
 * nCol    -> numero de colunas                                      * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * vL      -> arranjo local                                          * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void sGetLocalV(short *RESTRICT vG  ,short *RESTRICT vL
               ,INT *RESTRICT mapLG
               ,INT const nLin      ,INT const nCol)
{
  INT i,j,nG;
  
  if(nCol == 1){
    for(i=0;i<nLin;i++){
      nG      = mapLG[i];
      vL[i]   = vG[nG]; 
    }
  }
  else{
    for(i=0;i<nLin;i++){
      nG      = mapLG[i];
      for(j=0;j<nCol;j++)  
        MAT2D(i,j,vL,nCol) = MAT2D(nG,j,vG,nCol); 
    }
  }

}
/*********************************************************************/ 

/*********************************************************************     
 * MESHTOCSRIA : convertendo o grafo da malha para o formato csr     *     
 * ----------------------------------------------------------------- *     
 * Parametros de entrada :                                           *     
 * ----------------------------------------------------------------- *     
 * nen   - no por elementos                                          *     
 * eptr  - nao definido                                              *     
 * numel - numero de elementos                                       *     
 * ----------------------------------------------------------------- *     
 * Parametros de saida :                                             *     
 * ----------------------------------------------------------------- *     
 * eptr  - malha no formato CSR ( ponteiro )                         *     
 * esize - numero total de elementos no CSR                          *     
 *********************************************************************/    
INT meshToCsrIa(short *RESTRICT nen,INT *RESTRICT eptr
                ,INT const numel   ){
  INT i,eSize;
/*...*/
  eptr[0] = 0;
  for(i = 1; i < numel + 1;i++)
    eptr[i] = eptr[i-1] + nen[i-1];
  
  eSize = eptr[numel] - eptr[0];
/*...................................................................*/

  return eSize;
}
/*********************************************************************/    

/*********************************************************************     
 * MESHTOCSRIAJA : convertendo o grafo da malha para o formato csr   *     
 * ----------------------------------------------------------------- *     
 * Parametros de entrada :                                           *     
 * ----------------------------------------------------------------- *     
 * el    - conecitividade dos elementos                              *     
 * nen   - no por elementos                                          *     
 * eptr  - ponteiro CSR                                              *     
 * eind  - nao definido                                              *     
 * numel - numero de elementos                                       *     
 * maxno - numero maximo de nos por elemento                         *     
 * ----------------------------------------------------------------- *     
 * Parametros de saida :                                             *     
 * ----------------------------------------------------------------- *     
 * eind  - malha no formato CSR ( conectividade)                     *     
 *********************************************************************/     
void meshToCsrJa(INT *RESTRICT el  ,short *RESTRICT nen
                ,INT *RESTRICT eptr,INT *RESTRICT eind
                ,INT const numel   ,short const maxNo){
      INT i,iPoint;
      short j;
/*...*/
      for(i=0;i<numel;i++){
        iPoint = eptr[i];
        for(j=0;j<nen[i];j++)
           eind[iPoint+j] = MAT2D(i,j,el,maxNo)-1;
      }    
/*...................................................................*/

}
/*********************************************************************/     


/*********************************************************************/
void printMap(PartMesh pMesh  
             ,INT const nNode ,INT const nEl 
             ,short const myId
             ,char *nameOut   ,FILE *f){

  INT i;

  f = openFile(nameOut,"w");
  
  fprintf(f,"Myid: %d\n",myId);
  fprintf(f,"nPrcViz(El): %d\n",pMesh.iEl.nVizPart);
  
  fprintf(f,"\nmapPrcViz(El)\n");
  for(i=0;i<pMesh.iEl.nVizPart;i++)
    fprintf(f,"%d %hd\n",i,pMesh.iEl.vizPart[i]);
  
  
  fprintf(f,"\nmapInterface: nRcvs %d  nSends %d\n"
         ,pMesh.iEl.nRcvs,pMesh.iEl.nSends);
  
  for(i=0;i<pMesh.iEl.nVizPart+1;i++){
    fprintf(f,"%9d rcvs %9d send %9d\n"
           ,i,pMesh.iEl.iaRcvs[i],pMesh.iEl.iaSends[i]);
  }
  
  fprintf(f,"\nmapInterface: nRcvs %9d  nSends %9d\n"
            ,pMesh.iEl.nRcvs,pMesh.iEl.nSends);
  for(i=0;i<pMesh.iEl.nRcvs+pMesh.iEl.nSends;i++){
    fprintf(f,"fMapEl %9d %9d %9d\n",i,pMesh.iEl.fMap[i]+1
             ,pMesh.elLG[pMesh.iEl.fMap[i]]+1);
  }
  
  fprintf(f,"\nnPrcViz(No): %d\n\n",pMesh.iNo.nVizPart);
  
  fprintf(f,"mapPrcViz(No)\n");
  for(i=0;i<pMesh.iNo.nVizPart;i++)
    fprintf(f,"%d %hd\n",i,pMesh.iNo.vizPart[i]);
  
  fprintf(f,"\nmapInterface: nCom %d\n",pMesh.iNo.nCom);
  
  for(i=0;i<pMesh.iNo.nVizPart+1;i++){
    fprintf(f,"%9d iaComNo %9d\n"
           ,i,pMesh.iNo.iaComNo[i]);
  }
  
  fprintf(f,"\nmapInterface: nCom %9d\n"
            ,pMesh.iNo.nCom);
  for(i=0;i<pMesh.iNo.nCom;i++){
    fprintf(f,"fMapNo %9d %9d %9d\n",i,pMesh.iNo.fMap[i]+1
             ,pMesh.noLG[pMesh.iNo.fMap[i]]+1);
  }
  
  fprintf(f,"\nelLG: %9d\n",nEl);
  for(i=0;i<nEl;i++){
    fprintf(f,"elLG %9d %9d\n",i,pMesh.elLG[i]+1);
  }
  
  fprintf(f,"\nnoLG: %9d\n",nNode);
  for(i=0;i<nNode;i++){
    fprintf(f,"noLG %9d %9d\n",i,pMesh.noLG[i]+1);
  }
  
    

}
/*********************************************************************/

