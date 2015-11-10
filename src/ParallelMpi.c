#include<ParallelMpi.h>
/********************************************************************* 
 * STARTMPI: inicializado do MPI                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * argc, argv -> argumentos de linha de comando                      * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * nPrcs      -> numero de processos                                 * 
 * myId       -> id do processo                                      * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void mpiStart(int *argc,char **argv){

 int np=1,id=0;
 
 mpiVar.comm  = 0;
  
#ifdef _MPICH_
  MPI_Init(argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &mpiVar.comm);
  MPI_Comm_size(mpiVar.comm, &np);
  MPI_Comm_rank(mpiVar.comm, &id);
  
#endif

  mpiVar.nPrcs = np;
  mpiVar.myId  = id;
/*...................................................................*/
  
}
/*********************************************************************/ 

/********************************************************************* 
 * STOPMPI: finaliza os processos MPI                                * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void mpiStop(void){

#ifdef _MPICH_
  MPI_Finalize();
#endif
}
/*********************************************************************/ 

/********************************************************************* 
 * MPIWAIT: sicronza os processos MPI                                * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void mpiWait(void){

#ifdef _MPICH_
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  MPI_Barrier(mpiVar.comm);
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
#endif

}
/*********************************************************************/ 

/********************************************************************* 
 * GOBBALMESHQUALITY : obtem as estatistica globais da malha         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * mQl     -> estatisticas da malha                                  *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * mQl     -> estatisticas da malha atualizada                       *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void globalMeshQuality(MeshQuality *mQl){

#ifdef _MPICH_
  DOUBLE v;
/*... volume */
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  MPI_Allreduce(&mQl->volume,&v ,1           ,MPI_DOUBLE
               ,MPI_SUM     ,mpiVar.comm);
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mQl->volume = v;
/*...................................................................*/

/*... nonOrthMed */
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  MPI_Allreduce(&mQl->nonOrthMed,&v      ,1          , MPI_DOUBLE
               ,MPI_SUM         , mpiVar.comm);
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mQl->nonOrthMed = v/(DOUBLE) mpiVar.nPrcs;  
/*...................................................................*/

/*... nonOrthMAX */
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  MPI_Allreduce(&mQl->nonOrthMax,&v      , 1          , MPI_DOUBLE
            ,MPI_MAX            , mpiVar.comm);
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mQl->nonOrthMax = v;  
/*...................................................................*/

/*... skewMed*/
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  MPI_Allreduce(&mQl->skewMed   ,&v      ,1          , MPI_DOUBLE
            ,MPI_SUM            , mpiVar.comm);
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mQl->skewMed = v/(DOUBLE) mpiVar.nPrcs;  
/*...................................................................*/

/*... skewMax */
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  MPI_Allreduce(&mQl->skewMax   ,&v      ,1          , MPI_DOUBLE
            ,MPI_MAX            , mpiVar.comm);
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mQl->skewMax = v;  
/*...................................................................*/
#endif

}
/*********************************************************************/

/********************************************************************* 
 * DGLABALCEL :globaliza os valores das celulas (DOUBLE)             *
 * no master(myId=0)                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m       -> memoria                                                *
 * pMesh   -> particionamento atualizado                             *
 * uG      -> nao definido                                           *
 * uL      -> variavel local                                         *
 * numelNov-> numero de elementos sem sobreposicoes                  *
 * nfd1    -> graus de liberdade linha  (tensor)                     *
 * nfd2    -> graus de liberdade coluna (tensor)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * uG      -> variavel global                                        *
 *-------------------------------------------------------------------* 
 * OBS: uG so esta alocado no master                                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void dGlobalCel(Memoria *m         ,PartMesh *pMesh
               ,DOUBLE *restrict uG,DOUBLE *restrict uL
               ,INT const numelNov 
               ,short const ndf1   ,short const ndf2)
            
{           
#ifdef _MPICH_
  INT i,lEl;
  INT elG = pMesh->elG ;
  INT *elLG  = pMesh->elLG;
  INT size= elG*ndf1*ndf2;
  DOUBLE *x1=NULL; 
#endif

/*...*/  
  tm.overHeadGCelMpi = getTimeC() - tm.overHeadGCelMpi;
/*...................................................................*/
   
/*...*/
  if(mpiVar.nPrcs < 2) return;
/*...................................................................*/

#ifdef _MPICH_

/*...*/
  if(!mpiVar.myId)
    zero(uG,size,DOUBLEC);
/*...................................................................*/

/*...*/
  HccaAlloc(DOUBLE,m,x1,size,"x1Gmpi",false);
  zero(x1,size,DOUBLEC);
/*...................................................................*/

/*...*/
  if(ndf1 == 1){
    for(i=0;i<numelNov;i++){
      lEl     = elLG[i];
      x1[lEl] = uL[i]; 
    }
  }
  else if(ndf1 == 2){
    for(i=0;i<numelNov;i++){
      lEl    = elLG[i];
      MAT2D(lEl,0,x1,2) = MAT2D(i,0,uL,2); 
      MAT2D(lEl,1,x1,2) = MAT2D(i,1,uL,2); 
    }
  }
  else{
    for(i=0;i<numelNov;i++){
      lEl    = elLG[i];
      MAT2D(lEl,0,x1,3) = MAT2D(i,0,uL,3); 
      MAT2D(lEl,1,x1,3) = MAT2D(i,1,uL,3); 
      MAT2D(lEl,2,x1,3) = MAT2D(i,2,uL,3); 
    }
  }  
/*...................................................................*/

/*...*/
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  MPI_Reduce(x1,uG,size,MPI_DOUBLE,MPI_SUM,0,mpiVar.comm);
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*...................................................................*/


/*...*/  
  HccaDealloc(m,x1,"x1Gmpi",false);
/*...................................................................*/

/*...*/  
  tm.overHeadGCelMpi = getTimeC() - tm.overHeadGCelMpi;
/*...................................................................*/
#endif
}  
/*********************************************************************/
  
/********************************************************************* 
 * DGLABALNODE:globaliza os valores nodais (DOUBLE) no master(myId=0)* 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m       -> memoria                                                *
 * pMesh   -> particionamento atualizado                             *
 * uG      -> nao definido                                           *
 * uL      -> variavel local                                         *
 * nfd1    -> graus de liberdade linha  (tensor)                     *
 * nfd2    -> graus de liberdade coluna (tensor)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * uG      -> variavel global                                        *
 *-------------------------------------------------------------------* 
 * OBS: uG so esta alocado no master                                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void dGlobalNode(Memoria *m         ,PartMesh *pMesh
                ,DOUBLE *restrict uG,DOUBLE *restrict uL
                ,short const ndf1   ,short const ndf2)
            
{           
#ifdef _MPICH_
  INT i,lNo;
  INT nno1   = pMesh->nno1;
  INT nnG    = pMesh->nnG;
  INT *noLG  = pMesh->noLG;
  INT size=nnG*ndf1*ndf2;
  DOUBLE *x1=NULL; 
#endif

/*...*/  
  tm.overHeadGNodMpi = getTimeC() - tm.overHeadGNodMpi;
/*...................................................................*/
   
/*...*/
  if(mpiVar.nPrcs < 2) return;
/*...................................................................*/

#ifdef _MPICH_

/*...*/
  if(!mpiVar.myId)
    zero(uG,size,DOUBLEC);
/*...................................................................*/

/*...*/
  HccaAlloc(DOUBLE,m,x1,size,"x1Gmpi",false);
  zero(x1,size,DOUBLEC);
/*...................................................................*/

/*...*/
  if(ndf1 == 1){
    for(i=0;i<nno1;i++){
      lNo    = noLG[i];
      x1[lNo] = uL[i]; 
    }
  }
  else if(ndf1 == 2){
    for(i=0;i<nno1;i++){
      lNo    = noLG[i];
      MAT2D(lNo,0,x1,2) = MAT2D(i,0,uL,2); 
      MAT2D(lNo,1,x1,2) = MAT2D(i,1,uL,2); 
    }
  }
  else{
    for(i=0;i<nno1;i++){
      lNo    = noLG[i];
      MAT2D(lNo,0,x1,3) = MAT2D(i,0,uL,3); 
      MAT2D(lNo,1,x1,3) = MAT2D(i,1,uL,3); 
      MAT2D(lNo,2,x1,3) = MAT2D(i,2,uL,3); 
    }
  }  
/*...................................................................*/

/*...*/
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  MPI_Reduce(x1,uG,size,MPI_DOUBLE,MPI_SUM,0,mpiVar.comm);
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*...................................................................*/


/*...*/  
  HccaDealloc(m,x1,"x1Gmpi",false);
/*...................................................................*/

/*...*/  
  tm.overHeadGNodMpi = getTimeC() - tm.overHeadGNodMpi;
/*...................................................................*/
#endif
}  
/*********************************************************************/

/********************************************************************* 
 * GETBUFFER  : obetem os valores do buffer de comunicacao           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x       -> nao definido                                           *
 * xb      -> buffer de recebimento                                  *
 * fMap    -> mapa de equacoes                                       *
 * nRcvs   -> numero de equacoes no buffer                           *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * x       -> vetor com os valores atualizado                        *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void getBuffer(DOUBLE *restrict x    ,DOUBLE *restrict xb
              ,INT *restrict fMap    ,INT const nRcvs){
  
  INT i,lNeq;

  for(i=0;i<nRcvs;i++){
    lNeq    = fMap[i];
    x[lNeq] = xb[i];
  }
    
} 
/*********************************************************************/

/********************************************************************* 
 * MAKEBUFFER  : gera buffer de comunicacao                          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x       -> vetor com os valores                                   *
 * xb      -> nao definido                                           *
 * fMap    -> mapa de equacoes                                       *
 * nSends  -> numero de equacoes no buffer                           *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * xb      -> buffer de envio                                        *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void makeBuffer(DOUBLE *restrict x    ,DOUBLE *restrict xb
               ,INT *restrict fMap    ,INT const nSends){
    
  INT i,lNeq;

  for(i=0;i<nSends;i++){
    lNeq   = fMap[i];
    xb[i]  = x[lNeq];
  }
  
} 
/*********************************************************************/

/********************************************************************* 
 * GETBUFFERCEL : obetem os valores do buffer de comunicacao para as * 
 * celulas                                                           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x       -> nao definido                                           *
 * xb      -> buffer de recebimento                                  *
 * fMap    -> mapa de equacoes                                       *
 * nRcvs   -> numero de equacoes no buffer                           *
 * nfd1    -> graus de liberdade linha  (tensor)                     *
 * nfd2    -> graus de liberdade coluna (tensor)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * x       -> vetor com os valores atualizado                        *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void getBufferCel(DOUBLE *restrict x    ,DOUBLE *restrict xb
              ,INT *restrict fMap       ,INT const nRcvs
              ,short const ndf1         ,short const ndf2){

  INT i,lCel;

/*... ndf = 1 */  
  if( ndf2 == 1){
/*... ndm = 3*/
    if( ndf1 == 3)
      for(i=0;i<nRcvs;i++){
        lCel    = fMap[i];
        MAT2D(lCel,0,x,3) = MAT2D(i,0,xb,3);
        MAT2D(lCel,1,x,3) = MAT2D(i,1,xb,3);
        MAT2D(lCel,2,x,3) = MAT2D(i,2,xb,3);
      }
/*... ndm = 2*/
    else if( ndf1 == 2)
      for(i=0;i<nRcvs;i++){
        lCel    = fMap[i];
        MAT2D(lCel,0,x,2) = MAT2D(i,0,xb,2);
        MAT2D(lCel,1,x,2) = MAT2D(i,1,xb,2);
      }
/*... ndm = 1*/
    else if(ndf1 == 1) 
      for(i=0;i<nRcvs;i++){
        lCel    = fMap[i];
        x[lCel] = xb[i];
      }
  }
    
} 
/*********************************************************************/

/********************************************************************* 
 * MAKEBUFFERCEL  : gera buffer de comunicacao para as celulas       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x       -> vetor com os valores                                   *
 * xb      -> nao definido                                           *
 * fMap    -> mapa das celulas                                       *
 * nSends  -> numero de celulas no buffer                            *
 * nfd1    -> graus de liberdade linha  (tensor)                     *
 * nfd2    -> graus de liberdade coluna (tensor)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * xb      -> buffer de envio                                        *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void makeBufferCel(DOUBLE *restrict x    ,DOUBLE *restrict xb
                  ,INT *restrict fMap    ,INT const nSends
                  ,short const ndf1          ,short const ndf2){
    
  INT i,lCel;

/*... ndf = 1 */  
  if( ndf2 == 1){
/*... ndm = 3*/
    if( ndf1 == 3)
      for(i=0;i<nSends;i++){
        lCel    = fMap[i];
        MAT2D(i,0,xb,3) = MAT2D(lCel,0,x,3);
        MAT2D(i,1,xb,3) = MAT2D(lCel,1,x,3);
        MAT2D(i,2,xb,3) = MAT2D(lCel,2,x,3);
      }
/*... ndm = 2*/
    else if( ndf1 == 2)
      for(i=0;i<nSends;i++){
        lCel    = fMap[i];
        MAT2D(i,0,xb,2) = MAT2D(lCel,0,x,2);
        MAT2D(i,1,xb,2) = MAT2D(lCel,1,x,2);
      }
/*... ndm = 1*/
    else if(ndf1 == 1) 
      for(i=0;i<nSends;i++){
        lCel    = fMap[i];
        xb[i]   = x[lCel];
      }
  }
/*..................................................................*/  
} 
/********************************************************************/

/********************************************************************* 
 * DGETBUFFERNOD: obetem os valores do buffer de comunicacao para os * 
 * nos (DOUBLE)                                                      * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x       -> nao definido                                           *
 * xb      -> buffer de recebimento                                  *
 * fMap    -> mapa de equacoes                                       *
 * nRcvs   -> numero de equacoes no buffer                           *
 * nfd1    -> graus de liberdade linha  (tensor)                     *
 * nfd2    -> graus de liberdade coluna (tensor)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * x       -> vetor com os valores atualizado                        *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void dGetBufferNod(DOUBLE *restrict x    ,DOUBLE *restrict xb
              ,INT *restrict fMap        ,INT const nRcvs
              ,short const ndf1          ,short const ndf2){

  INT i,lNod;

/*... ndf2 = 1 */  
  if( ndf2 == 1){
/*... ndf1 = 3*/
    if( ndf1 == 3)
      for(i=0;i<nRcvs;i++){
        lNod    = fMap[i];
        MAT2D(lNod,0,x,3) += MAT2D(i,0,xb,3);
        MAT2D(lNod,1,x,3) += MAT2D(i,1,xb,3);
        MAT2D(lNod,2,x,3) += MAT2D(i,2,xb,3);
      }
/*... ndf1 = 2*/
    else if( ndf1 == 2)
      for(i=0;i<nRcvs;i++){
        lNod    = fMap[i];
        MAT2D(lNod,0,x,2) += MAT2D(i,0,xb,2);
        MAT2D(lNod,1,x,2) += MAT2D(i,1,xb,2);
      }
/*... ndf1 = 1*/
    else if(ndf1 == 1) 
      for(i=0;i<nRcvs;i++){
        lNod    = fMap[i];
        x[lNod] += xb[i];
      }
  }
    
} 
/*********************************************************************/

/********************************************************************* 
 * DMAKEBUFFERNOD : gera buffer de comunicacao para os nos (DOUBLE)  * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x       -> vetor com os valores                                   *
 * xb      -> nao definido                                           *
 * fMap    -> mapa das celulas                                       *
 * nSends  -> numero de celulas no buffer                            *
 * nfd1    -> graus de liberdade linha  (tensor)                     *
 * nfd2    -> graus de liberdade coluna (tensor)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * xb      -> buffer de envio                                        *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void dMakeBufferNod(DOUBLE *restrict x    ,DOUBLE *restrict xb
                  ,INT *restrict fMap    ,INT const nSends
                  ,short const ndf1          ,short const ndf2){
    
  INT i,lNo;

/*... ndf2 = 1 */  
  if( ndf2 == 1){
/*... ndf1 = 3*/
    if( ndf1 == 3)
      for(i=0;i<nSends;i++){
        lNo     = fMap[i];
        MAT2D(i,0,xb,3) = MAT2D(lNo,0,x,3);
        MAT2D(i,1,xb,3) = MAT2D(lNo,1,x,3);
        MAT2D(i,2,xb,3) = MAT2D(lNo,2,x,3);
      }
/*... ndf2 = 2*/
    else if( ndf1 == 2)
      for(i=0;i<nSends;i++){
        lNo    = fMap[i];
        MAT2D(i,0,xb,2) = MAT2D(lNo,0,x,2);
        MAT2D(i,1,xb,2) = MAT2D(lNo,1,x,2);
      }
/*... ndf3 = 1*/
    else if(ndf1 == 1) 
      for(i=0;i<nSends;i++){
        lNo    = fMap[i];
        xb[i]  = x[lNo];
      }
  }
/*..................................................................*/  
} 
/********************************************************************/

/********************************************************************* 
 * IGETBUFFERNOD : obetem os valores do buffer de comunicacao para os* 
 * nos (INT)                                                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x       -> nao definido                                           *
 * xb      -> buffer de recebimento                                  *
 * fMap    -> mapa de equacoes                                       *
 * nRcvs   -> numero de equacoes no buffer                           *
 * nfd1    -> graus de liberdade linha  (tensor)                     *
 * nfd2    -> graus de liberdade coluna (tensor)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * x       -> vetor com os valores atualizado                        *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void iGetBufferNod(INT *restrict x     ,INT *restrict xb
                 ,INT *restrict fMap  ,INT const nRcvs
                 ,short const ndf1    ,short const ndf2){

  INT i,lNod;

/*... ndf2 = 1 */  
  if( ndf2 == 1){
/*... ndf1 = 3*/
    if( ndf1 == 3)
      for(i=0;i<nRcvs;i++){
        lNod    = fMap[i];
        MAT2D(lNod,0,x,3) += MAT2D(i,0,xb,3);
        MAT2D(lNod,1,x,3) += MAT2D(i,1,xb,3);
        MAT2D(lNod,2,x,3) += MAT2D(i,2,xb,3);
      }
/*... ndf1 = 2*/
    else if( ndf1 == 2)
      for(i=0;i<nRcvs;i++){
        lNod    = fMap[i];
        MAT2D(lNod,0,x,2) += MAT2D(i,0,xb,2);
        MAT2D(lNod,1,x,2) += MAT2D(i,1,xb,2);
      }
/*... ndf1 = 1*/
    else if(ndf1 == 1) 
      for(i=0;i<nRcvs;i++){
        lNod    = fMap[i];
        x[lNod] += xb[i];
      }
  }
    
} 
/*********************************************************************/

/********************************************************************* 
 * IMAKEBUFFERNOD  : gera buffer de comunicacao para os nos (INT)    * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x       -> vetor com os valores                                   *
 * xb      -> nao definido                                           *
 * fMap    -> mapa das celulas                                       *
 * nSends  -> numero de celulas no buffer                            *
 * nfd1    -> graus de liberdade linha  (tensor)                     *
 * nfd2    -> graus de liberdade coluna (tensor)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * xb      -> buffer de envio                                        *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void iMakeBufferNod(INT *restrict x       ,INT *restrict xb
                   ,INT *restrict fMap    ,INT const nSends
                   ,short const ndf1      ,short const ndf2){
    
  INT i,lNo;

/*... ndf = 1 */  
  if( ndf2 == 1){
/*... ndf = 3*/
    if( ndf1 == 3)
      for(i=0;i<nSends;i++){
        lNo     = fMap[i];
        MAT2D(i,0,xb,3) = MAT2D(lNo,0,x,3);
        MAT2D(i,1,xb,3) = MAT2D(lNo,1,x,3);
        MAT2D(i,2,xb,3) = MAT2D(lNo,2,x,3);
      }
/*... ndf = 2*/
    else if( ndf1 == 2)
      for(i=0;i<nSends;i++){
        lNo    = fMap[i];
        MAT2D(i,0,xb,2) = MAT2D(lNo,0,x,2);
        MAT2D(i,1,xb,2) = MAT2D(lNo,1,x,2);
      }
/*... ndf = 1*/
    else if(ndf1 == 1) 
      for(i=0;i<nSends;i++){
        lNo    = fMap[i];
        xb[i]  = x[lNo];
      }
  }
/*..................................................................*/  
} 
/********************************************************************/

/********************************************************************* 
 * COMUNICATENEQ : comunicao das equacoes de interface               * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * iNeq    -> mapa de equacoes de interface                          *
 * x       -> vetor a ser comunicado                                 *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * x       -> atualizado                                             *
 *-------------------------------------------------------------------* 
 * OBS: equacoes numeradas de primeiro por vizinho e depois por      *
 * numeracao global                                                  * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
  void comunicateNeq(Interface *iNeq,DOUBLE *restrict x){

#ifdef _MPICH_
  INT nRcvs  = iNeq->nRcvs;
  INT nSends = iNeq->nSends;
  INT k,kk;
  unsigned short partId,nPart,nVizParts = iNeq->nVizPart; 

/*...*/
  tm.overHeadNeqMpi    = getTimeC() - tm.overHeadNeqMpi;
/*...................................................................*/    
  
  if(mpiVar.nPrcs < 2 ) return;

/*... gerando o buffer de envio*/
   makeBuffer(x,&iNeq->xb[nRcvs],&iNeq->fMap[nRcvs],nSends);
/*...................................................................*/    

/*...*/  
  for(nPart=0;nPart < nVizParts;nPart++){
    partId = iNeq->vizPart[nPart];
/*... enviando para as particoes vizinhas*/
    k      = iNeq->iaSends[nPart];
    kk     = iNeq->iaSends[nPart+1] - iNeq->iaSends[nPart];
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
    MPI_Isend(&iNeq->xb[k]    ,kk
             ,MPI_DOUBLE      ,partId    
             ,mpiVar.myId     ,mpiVar.comm
             ,&mpiVar.sendRequest[nPart]);
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*... recebimentos das particoes vizinhas*/
    k      = iNeq->iaRcvs[nPart];
    kk     = iNeq->iaRcvs[nPart+1] - iNeq->iaRcvs[nPart];
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
    MPI_Irecv(&iNeq->xb[k]   ,kk
             ,MPI_DOUBLE     ,partId    
             ,partId         ,mpiVar.comm
             ,&mpiVar.recvRequest[nPart]);
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  }  
/*...................................................................*/    
    
/*... espera o recebimentos dos dados*/
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.recvRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*...................................................................*/    
  
/*...*/
  getBuffer(x,iNeq->xb,iNeq->fMap,nRcvs);
/*...................................................................*/    

/*... espera o envio de todos os dados*/
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.sendRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*...................................................................*/    

/*...*/
  tm.overHeadNeqMpi    = getTimeC() - tm.overHeadNeqMpi;
/*...................................................................*/    

#endif
}
/*********************************************************************/

/********************************************************************* 
 * COMUNICATECEL : comunicao valores das celulas nas interface       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * iEl     -> mapa das celulas nas interfaces                        *
 * x       -> vetor a ser comunicado                                 *
 * nfd1    -> graus de liberdade linha  (tensor)                     *
 * nfd2    -> graus de liberdade coluna (tensor)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * x       -> atualizado                                             *
 *-------------------------------------------------------------------* 
 * OBS: celulas numeradas primeiro por vizinho e depois por numeracao* 
 * global                                                            * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
  void comunicateCel(Interface *iCel ,DOUBLE *restrict x
                    ,short const ndf1,short const ndf2){

#ifdef _MPICH_
  INT nRcvs  = iCel->nRcvs;
  INT nSends = iCel->nSends;
  INT k,kk,nst=ndf1*ndf2;
  unsigned short partId,nPart,nVizParts = iCel->nVizPart; 

/*...*/
  tm.overHeadCelMpi    = getTimeC() - tm.overHeadCelMpi;
/*...................................................................*/    
  
  if(mpiVar.nPrcs < 2 ) return;

/*... gerando o buffer de envio*/
   makeBufferCel(x                     ,&iCel->xb[nRcvs*nst]
                ,&iCel->fMap[nRcvs]    ,nSends
                ,ndf1                  ,ndf2);
/*...................................................................*/    

/*...*/   
  for(nPart=0;nPart < nVizParts;nPart++){
    partId = iCel->vizPart[nPart];
/*... enviando para as particoes vizinhas*/
    k      = iCel->iaSends[nPart];
    kk     = iCel->iaSends[nPart+1] - iCel->iaSends[nPart];
    k      *= nst;
    kk     *= nst;
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
    MPI_Isend(&iCel->xb[k]    ,kk
             ,MPI_DOUBLE      ,partId    
             ,mpiVar.myId     ,mpiVar.comm
             ,&mpiVar.sendRequest[nPart]);
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*... recebimentos das particoes vizinhas*/
    k      = iCel->iaRcvs[nPart];
    kk     = iCel->iaRcvs[nPart+1] - iCel->iaRcvs[nPart];
    k      *= nst;
    kk     *= nst;
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
    MPI_Irecv(&iCel->xb[k]   ,kk
             ,MPI_DOUBLE     ,partId    
             ,partId         ,mpiVar.comm
             ,&mpiVar.recvRequest[nPart]);
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  }  
/*...................................................................*/    

/*... espera o recebimentos dos dados*/
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.recvRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*...................................................................*/    
  
/*...*/
  getBufferCel(x         ,iCel->xb
              ,iCel->fMap,nRcvs
              ,ndf1      ,ndf2);
/*...................................................................*/    

/*... espera o envio de todos os dados*/
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.sendRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*...................................................................*/    

/*...*/
  tm.overHeadCelMpi    = getTimeC() - tm.overHeadCelMpi;
/*...................................................................*/    

#endif
}
/*********************************************************************/

/********************************************************************* 
 * DCOMUNICATENOD : comunicao valores das nos nas interface (DOUBLE) * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * iNo     -> mapa das celulas nas interfaces                        *
 * x       -> vetor a ser comunicado                                 *
 * nfd1    -> graus de liberdade linha  (tensor)                     *
 * nfd2    -> graus de liberdade coluna (tensor)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * x       -> atualizado                                             *
 *-------------------------------------------------------------------* 
 * OBS: nos numerados primeiro por vizinho e depois por numeracao    * 
 * global                                                            * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
  void dComunicateNod(InterfaceNo *iNo ,DOUBLE *restrict x
                     ,short const ndf1 ,short const ndf2){

#ifdef _MPICH_
  INT nCom   = iNo->nCom;
  INT k,kk,nst=ndf1*ndf2;
  unsigned short partId,nPart,nVizParts = iNo->nVizPart; 

/*...*/
  tm.overHeadNodMpi    = getTimeC() - tm.overHeadNodMpi;
/*...................................................................*/    

  if(mpiVar.nPrcs < 2 ) return;

/*... gerando o buffer de envio*/
   dMakeBufferNod(x                    ,&iNo->xb[nCom*nst]
                ,iNo->fMap             ,nCom   
                ,ndf1                  ,ndf2);
/*...................................................................*/    

/*...*/   
  for(nPart=0;nPart < nVizParts;nPart++){
    k      = iNo->iaComNo[nPart];
    kk     = iNo->iaComNo[nPart+1] - iNo->iaComNo[nPart];
    partId = iNo->vizPart[nPart];
    k      *= nst;
    kk     *= nst;
/*... enviando para as particoes vizinhas*/
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
    MPI_Isend(&iNo->xb[nCom*nst+k]     ,kk
             ,MPI_DOUBLE               ,partId    
             ,mpiVar.myId              ,mpiVar.comm
             ,&mpiVar.sendRequest[nPart]);
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*... recebimentos das particoes vizinhas*/
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
    MPI_Irecv(&iNo->xb[k]   ,kk
             ,MPI_DOUBLE     ,partId    
             ,partId         ,mpiVar.comm
             ,&mpiVar.recvRequest[nPart]);
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  }  
/*...................................................................*/    

/*... espera o recebimentos dos dados*/
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.recvRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*...................................................................*/    
  
/*...*/
  dGetBufferNod(x         ,iNo->xb
              ,iNo->fMap ,nCom
              ,ndf1      ,ndf2);
/*...................................................................*/    

/*... espera o envio de todos os dados*/
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.sendRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*...................................................................*/    

/*...*/
  tm.overHeadNodMpi    = getTimeC() - tm.overHeadNodMpi;
/*...................................................................*/    

#endif
}
/*********************************************************************/

/********************************************************************* 
 * ICOMUNICATENOD : comunicao valores das nos nas interface (INT)    * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * iNo     -> mapa das celulas nas interfaces                        *
 * x       -> vetor a ser comunicado                                 *
 * nfd1    -> graus de liberdade linha  (tensor)                     *
 * nfd2    -> graus de liberdade coluna (tensor)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * x       -> atualizado                                             *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
  void iComunicateNod(InterfaceNo *iNo ,INT *restrict x
                     ,short const ndf1,short const ndf2){

#ifdef _MPICH_
  INT nCom   = iNo->nCom;
  INT k,kk,nst=ndf1*ndf2;
  unsigned short partId,nPart,nVizParts = iNo->nVizPart; 

/*...*/
  tm.overHeadNodMpi    = getTimeC() - tm.overHeadNodMpi;
/*...................................................................*/    

  if(mpiVar.nPrcs < 2 ) return;

/*... gerando o buffer de envio*/
  iMakeBufferNod(x                     ,&iNo->xi[nCom*nst]
                ,iNo->fMap             ,nCom   
                ,ndf1                  ,ndf2);
/*...................................................................*/    

/*...*/   
  for(nPart=0;nPart < nVizParts;nPart++){
    k      = iNo->iaComNo[nPart];
    kk     = iNo->iaComNo[nPart+1] - iNo->iaComNo[nPart];
    partId = iNo->vizPart[nPart];
    k      *= nst;
    kk     *= nst;
/*... enviando para as particoes vizinhas*/
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
    MPI_Isend(&iNo->xi[nCom*nst+k]     ,kk
             ,MPI_INT                  ,partId    
             ,mpiVar.myId              ,mpiVar.comm
             ,&mpiVar.sendRequest[nPart]);
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*... recebimentos das particoes vizinhas*/
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
    MPI_Irecv(&iNo->xi[k]    ,kk
             ,MPI_INT        ,partId    
             ,partId         ,mpiVar.comm
             ,&mpiVar.recvRequest[nPart]);
    tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  }  
/*...................................................................*/    

/*... espera o recebimentos dos dados*/
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.recvRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*...................................................................*/    
  
/*...*/
  iGetBufferNod(x         ,iNo->xi
              ,iNo->fMap ,nCom
              ,ndf1      ,ndf2);
/*...................................................................*/    

/*... espera o envio de todos os dados*/
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.sendRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadTotalMpi = getTimeC() - tm.overHeadTotalMpi;
/*...................................................................*/    

/*...*/
  tm.overHeadNodMpi    = getTimeC() - tm.overHeadNodMpi;
/*...................................................................*/    

#endif
}
/*********************************************************************/

/********************************************************************* 
 * COMUNICATEMESH : gera o mapa de inteface dos elementes e comunica * 
 * as particoes para os vizinhos                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m       -> vetor de memoria                                       *
 * mesh0   -> malha original                                         *
 * mesh    -> nao definido                                           *
 * pMesh   -> particionamento                                        *
 * loadsD1 -> cargas nas faces (D1)                                  *
 * loadsT1 -> cargas nas faces (T1)                                  *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * mesh    -> malha particionada e comunicada para o escravos        *
 * pMesh   -> particionamento atualizado                             *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void comunicateMesh(Memoria *m
                   ,Mesh *mesh0     ,Mesh *mesh
                   ,PartMesh *pMesh 
                   ,Loads *loadsD1  ,Loads *loadsT1) 
{

#ifdef _MPICH_
  bool *aux1 = NULL,*aux2 = NULL,*aux3=NULL;
  INT i,numelNov,numelOv,lNel,lNnode,nNodeNov,nNodeOv,nno1;
  INT *elLG=NULL,*elGL=NULL,*noLG=NULL,*noGL=NULL;
  INT *lEl=NULL,*fMap=NULL,*iaRcvs=NULL,*iaSends=NULL;
  INT *fMapNo=NULL,*iaComNo=NULL;
  INT size,nRcvs,nSends,nComNo,kk;
  INT maxGrade;
  DOUBLE *nD=NULL;
  short  *nS=NULL;
  short sSize=sizeof(short); 
  short iSize=sizeof(INT); 
  short dSize=sizeof(DOUBLE); 
  short *maxVizPart = NULL,*maxVizPartNo = NULL;
  unsigned short nPart,nVizPart,nVizPartNo;
  unsigned nPrcs=mpiVar.nPrcs,myId=mpiVar.myId; 
  short st=0;
  short ndm,maxViz,maxNo,numat;
  short ndfD[MAX_DIF_EQ],ndfT[MAX_TRANS_EQ],ndfF;
  short maxNdf;

/*...*/
  ndfF = 0;
  for(i=0;i<MAX_TRANS_EQ;i++)
    ndfT[i]  = 0;
  for(i=0;i<MAX_DIF_EQ;i++)
    ndfD[i]  = 0;
/*...................................................................*/
  
/*... alocando varaivel auxiliar*/  
  if(!myId){

/*... gera o incidencia do elementos*/
    mesh0->noIncid.nincid = (INT *) malloc(mesh0->nnode*iSize);
    zero(mesh0->noIncid.nincid,mesh->nnode,INTC);
    nodeGrade(mesh0->elm.node,mesh0->noIncid.nincid
             ,mesh0->elm.nen ,&maxGrade
             ,mesh0->nnode   ,mesh0->numel
             ,mesh0->maxNo);
/*...................................................................*/
 
/*...*/
    mesh0->noIncid.incid = (INT *) malloc(mesh0->nnode*maxGrade*iSize);
    zero(mesh0->noIncid.incid,mesh->nnode*maxGrade,INTC);
    elmIncid(mesh0->elm.node      ,mesh0->noIncid.incid
            ,mesh0->noIncid.nincid,mesh0->elm.nen
            ,mesh0->nnode         ,mesh0->numel
            ,maxGrade             ,mesh0->maxNo);
/*...................................................................*/  
 
    aux1 = (bool *) malloc(mesh0->nnode*sizeof(bool));
    ERRO_MALLOC_MPI(aux1,"fNode",__LINE__,__FILE__,__func__,st);

    aux2 = (bool *) malloc(mesh0->numel*sizeof(bool)); 
    ERRO_MALLOC_MPI(aux2,"fEl"  ,__LINE__,__FILE__,__func__,st);
    
    aux3 = (bool *) malloc(nPrcs*sizeof(bool)); 
    ERRO_MALLOC_MPI(aux2,"fEl"  ,__LINE__,__FILE__,__func__,st);
    
    noGL = (INT *) malloc(mesh0->nnode*iSize); 
    ERRO_MALLOC_MPI(noGL,"noGLTemp"  ,__LINE__,__FILE__,__func__,st);
    
    elGL = (INT *) malloc(mesh0->numel*iSize); 
    ERRO_MALLOC_MPI(elGL,"elGLTemp"  ,__LINE__,__FILE__,__func__,st);
    
    maxVizPart = (short *) malloc(nPrcs*sSize); 
    ERRO_MALLOC_MPI(maxVizPart,"mvp",__LINE__,__FILE__,__func__,st);
    
    maxVizPartNo = (short *) malloc(nPrcs*sSize); 
    ERRO_MALLOC_MPI(maxVizPartNo
                   ,"mvpNo",__LINE__,__FILE__,__func__,st);
      
    fMap = (INT *) malloc(2*mesh0->numel*iSize); 
    ERRO_MALLOC(fMap,"fMap"  ,__LINE__,__FILE__,__func__);
    
    fMapNo = (INT *) malloc(2*mesh0->nnode*iSize); 
    ERRO_MALLOC(fMap,"fMapNo"  ,__LINE__,__FILE__,__func__);
    
    iaRcvs = (INT *) malloc((nPrcs+1)*iSize); 
    ERRO_MALLOC(iaRcvs,"iaRcvs"  ,__LINE__,__FILE__,__func__);
    
    iaSends= (INT *) malloc((nPrcs+1)*iSize); 
    ERRO_MALLOC(iaSends,"iaSends"  ,__LINE__,__FILE__,__func__);
    
    iaComNo = (INT *) malloc((nPrcs+1)*iSize); 
    ERRO_MALLOC(iaRcvs,"iaComNo"  ,__LINE__,__FILE__,__func__);
    
  }
/*...................................................................*/

    
  MPI_Bcast(&st,1,MPI_SHORT,0,mpiVar.comm);
  if( st == -1){
    mpiStop();
    exit(EXIT_FAILURE); 
  }

/*... variaveis compartilhada*/
  if(!myId){ 
    pMesh->nnG  = mesh0->nnode;
    pMesh->elG  = mesh0->numel;
    mesh->ndm   = mesh0->ndm;
    mesh->numat = mesh0->numat;
    mesh->maxViz= mesh0->maxViz;
    mesh->maxNo = mesh0->maxNo;
    for(i=0;i<MAX_TRANS_EQ;i++)
      mesh->ndfT[i]  = mesh0->ndfT[i];
    for(i=0;i<MAX_DIF_EQ;i++)
      mesh->ndfD[i]  = mesh0->ndfD[i];
  }
  

  MPI_Bcast(&pMesh->nnG  ,1           ,MPI_INT  ,0,mpiVar.comm);
  MPI_Bcast(&pMesh->elG  ,1           ,MPI_INT  ,0,mpiVar.comm);
  MPI_Bcast(&mesh->ndm   ,1           ,MPI_SHORT,0,mpiVar.comm);
  MPI_Bcast(&mesh->numat ,1           ,MPI_SHORT,0,mpiVar.comm);
  MPI_Bcast(&mesh->maxViz,1           ,MPI_SHORT,0,mpiVar.comm);
  MPI_Bcast(&mesh->maxNo ,1           ,MPI_SHORT,0,mpiVar.comm);
  MPI_Bcast(&mesh->ndfT  ,MAX_TRANS_EQ,MPI_SHORT,0,mpiVar.comm);
  MPI_Bcast(&mesh->ndfD  ,MAX_DIF_EQ  ,MPI_SHORT,0,mpiVar.comm);
  MPI_Bcast(&mesh->ndfD  ,MAX_DIF_EQ  ,MPI_SHORT,0,mpiVar.comm);
/*... loadsD1*/
  for(i=0;i<MAXLOADD1;i++){
    MPI_Bcast(&loadsD1[i].type  ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&loadsD1[i].np    ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&loadsD1[i].par   ,MAXLOADPARAMETER ,MPI_DOUBLE
              ,0,mpiVar.comm);
  }
/*...................................................................*/

/*... loadsT1*/
  for(i=0;i<MAXLOADT1;i++){
    MPI_Bcast(&loadsT1[i].type  ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&loadsT1[i].np    ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&loadsT1[i].par   ,MAXLOADPARAMETER ,MPI_DOUBLE
              ,0,mpiVar.comm);
  }
/*...................................................................*/

/*...*/
  ndm    = mesh->ndm;
  maxViz = mesh->maxViz;
  maxNo  = mesh->maxNo;
  numat  = mesh->numat;
  for(i=0;i<MAX_TRANS_EQ;i++)
    ndfT[i]  = mesh->ndfT[i];
  for(i=0;i<MAX_DIF_EQ;i++)
    ndfD[i]  = mesh->ndfD[i];
/*...................................................................*/

/*... alocando materiais*/

/*... Prop*/ 
  HccaAlloc(DOUBLE,m,mesh->elm.material.prop,MAXPROP*numat     
           ,"propP" ,_AD_);
/*... type*/ 
  HccaAlloc(short   ,m,mesh->elm.material.type,numat     
           ,"typeP" ,_AD_);
  
/*... zerando os variavies*/
  zero(mesh->elm.material.prop,MAXPROP*numat,DOUBLEC);
  zero(mesh->elm.material.type,numat,"short");
/*...................................................................*/

  if(!myId){ 
    for(i=0;i<MAXPROP*numat;i++)
      mesh->elm.material.prop[i]  = mesh0->elm.material.prop[i];
    for(i=0;i<numat;i++)
      mesh->elm.material.type[i]  = mesh0->elm.material.type[i];
  }
  
  MPI_Bcast(mesh->elm.material.prop ,MAXPROP*numat
           ,MPI_DOUBLE,0,mpiVar.comm);
  MPI_Bcast(mesh->elm.material.type ,numat        
           ,MPI_SHORT,0,mpiVar.comm);

/*... dividindo a malha*/
  for(nPart=0;nPart<nPrcs;nPart++){
/*... so o processo master gera os mapas*/
    if(!myId){ 
      for(i=0;i<mesh0->nnode;i++)
        aux1[i] = false;
    
      for(i=0;i<mesh0->numel;i++)
        aux2[i] = false;
/*... obtendo o numero local de elementos e nos*/    
      getNumberLocalMesh(pMesh->ep      ,pMesh->np              
                        ,aux1           ,aux2
                        ,mesh0->elm.node,mesh0->elm.adj.nelcon
                        ,mesh0->elm.nen ,mesh0->elm.adj.nViz
                        ,mesh0->nnode   ,mesh0->numel    
                        ,maxNo          ,maxViz            
                        ,&numelNov      ,&numelOv
                        ,&nNodeNov      ,&nNodeOv
                        ,&nno1      
                        ,nPart);
/*...................................................................*/  

/*... numero de elementos totais*/
      lNel   = numelNov + numelOv;
      lNnode = nNodeNov + nNodeOv;
/*... map de elemento*/
      for(i=0;i<mesh0->numel;i++)
        aux2[i] = false;
      
      elLG = (INT *) malloc(lNel*iSize); 
      ERRO_MALLOC(elLG,"elLGtemp"  ,__LINE__,__FILE__,__func__);
      getMapElm(pMesh->ep             ,elLG
              ,elGL                   ,aux2
              ,mesh0->elm.node        ,mesh0->elm.adj.nelcon
              ,mesh0->elm.adj.nViz
              ,mesh0->numel           ,lNel
              ,maxViz                 ,nPart);
/*...................................................................*/

/*... map de nos*/
      for(i=0;i<mesh0->nnode;i++)
        aux1[i] = false;
      
      noLG = (INT *) malloc(lNnode*iSize); 
      ERRO_MALLOC(noLG,"noLGtemp"  ,__LINE__,__FILE__,__func__);
      getMapNode(pMesh->ep         ,pMesh->np  
              ,noLG                ,noGL    
              ,aux1
              ,mesh0->elm.node     ,mesh0->elm.adj.nelcon
              ,mesh0->elm.nen      ,mesh0->elm.adj.nViz
              ,mesh0->nnode        ,mesh0->numel       
              ,maxNo               ,maxViz            
              ,lNnode              ,nPart);
/*...................................................................*/

/*... mapa de vizinhos de elementos*/
      nVizPart = getMapViz(pMesh->ep  ,elLG
                          ,maxVizPart 
                          ,lNel       ,numelNov
                          ,nPart      ,nPrcs );
/*...................................................................*/

/*... mapa de vizinhos de no*/
      nVizPartNo =
      getMapVizNo(pMesh->ep            ,noLG
                 ,maxVizPartNo         ,aux3
                 ,mesh0->noIncid.nincid,mesh0->noIncid.incid
                 ,nNodeNov             ,maxGrade              
                 ,nPart                ,nPrcs);
/*...................................................................*/

/*... mapa da interface de elementos*/
      getMapInterfaceEl(pMesh->ep
                       ,mesh0->elm.adj.nelcon   ,mesh0->elm.adj.nViz
                       ,elLG                    ,aux2                    
                       ,&nRcvs                  ,&nSends
                       ,iaRcvs                  ,iaSends 
                       ,maxVizPart              ,fMap 
                       ,lNel                    ,numelNov
                       ,nPart                   ,nPrcs
                       ,nVizPart                ,maxViz);
/*...................................................................*/

/*... mapa da interface de nos*/
      getMapInterfaceNo(pMesh->ep
                     ,noLG                    ,aux1  
                     ,&nComNo                 ,iaComNo
                     ,mesh0->noIncid.nincid   ,mesh0->noIncid.incid
                     ,maxVizPartNo            ,fMapNo
                     ,nNodeNov                           
                     ,nPart                   ,nPrcs
                     ,nVizPartNo              ,maxGrade);
/*...................................................................*/

//    printf("rank %d nNodeNov %d nNodeOv %d\n",nPart,nNodeNov,nNodeOv);
//    for(i=0;i<nRcvs;i++)
//      printf("fMapEl %d %d %d\n",i+1,fMap[i]+1,elLG[fMap[i]]+1);
    
 //   printf("\n");
 //   for(i=nRcvs;i<nRcvs+nSends;i++)
 //     printf("fMapEl %d %d %d\n",i+1,fMap[i]+1,elLG[fMap[i]]+1);
      
 //   printf("\n");
 //   for(i=0;i<nComNo ;i++)
 //     printf("fMapNo %d %d %d\n",i+1,fMapNo[i]+1,noLG[fMapNo[i]]+1);
  
/*...................................................................*/

/*... partionamento local para o processo master(mater)*/
      if(!nPart){
        mesh->nnode     = lNnode;
        mesh->numel     = lNel;
        mesh->numelNov  = numelNov;
        mesh->nnodeOv   = nNodeOv;
        mesh->nnodeNov  = nNodeNov;
        pMesh->nno1     = nno1;
/*... particionamento*/     

/*... partViz (El)*/
        pMesh->iEl.nVizPart = nVizPart;
        HccaAlloc(short,m,pMesh->iEl.vizPart , nVizPart   
                 ,"nVizPart"                 ,false);
        zero(pMesh->iEl.vizPart,nVizPart,"short");


        for(i=0;i<nVizPart;i++)
          pMesh->iEl.vizPart[i] = maxVizPart[i];
        
/*... partViz (No)*/
        pMesh->iNo.nVizPart = nVizPartNo;
        HccaAlloc(short,m,pMesh->iNo.vizPart , nVizPartNo   
                 ,"nVizPartNo"              ,false);
        zero(pMesh->iNo.vizPart,nVizPartNo,"short");
        
        for(i=0;i<nVizPartNo;i++)
          pMesh->iNo.vizPart[i] = maxVizPartNo[i];
/*..................................................................*/
        
        pMesh->iEl.nRcvs  = nRcvs;      
        pMesh->iEl.nSends = nSends;      
        pMesh->iNo.nCom   = nComNo;      

/*... fMapEl(buffer de comunicacao)*/ 
        HccaAlloc(INT   ,m   ,pMesh->iEl.fMap , nRcvs + nSends   
                 ,"fMap",_AD_);
        zero(pMesh->iEl.fMap     , nRcvs+nSends,INTC);
        
        for(i=0;i<nRcvs+nSends;i++)
          pMesh->iEl.fMap[i] = fMap[i];

/*... fMapNo(buffer de comunicacao)*/ 
        HccaAlloc(INT   ,m   ,pMesh->iNo.fMap , nComNo   
                 ,"fMapNo",_AD_);
        zero(pMesh->iNo.fMap     , nComNo,INTC);
        
        for(i=0;i<nComNo;i++)
          pMesh->iNo.fMap[i] = fMapNo[i];
/*..................................................................*/
 
/*... send e rcvs*/       
        kk = nVizPart+1;

        HccaAlloc(INT   ,m   ,pMesh->iEl.iaSends , kk        
                 ,"iaSends",_AD_);
        zero(pMesh->iEl.iaSends   , kk,INTC);
        
        HccaAlloc(INT   ,m   ,pMesh->iEl.iaRcvs  , kk        
                 ,"iaRcvs",_AD_);
        zero(pMesh->iEl.iaRcvs   , kk,INTC);
        
        for(i=0;i<kk;i++){
          pMesh->iEl.iaSends[i] = iaSends[i];
          pMesh->iEl.iaRcvs[i]  = iaRcvs[i];
        }

/*... iaComNo*/
        kk = nVizPartNo+1;

        HccaAlloc(INT   ,m   ,pMesh->iNo.iaComNo , kk        
                 ,"iaComNo",_AD_);
        zero(pMesh->iNo.iaComNo   , kk,INTC);
        
        for(i=0;i<kk;i++){
          pMesh->iNo.iaComNo[i] = iaComNo[i];
        }
/*..................................................................*/

/*... elLG( numeracao local-global de elementos)*/
        HccaAlloc(INT   ,m   ,pMesh->elLG , lNel        
                 ,"elLG"   ,_AD_);
        zero(pMesh->elLG   , lNel,INTC);
        for(i=0;i<lNel;i++)
          pMesh->elLG[i] = elLG[i];
/*..................................................................*/

/*... noLG( numeracao local-global de nos)*/
        HccaAlloc(INT   ,m   ,pMesh->noLG , lNnode        
                 ,"noLG"   ,_AD_);
        zero(pMesh->noLG   , lNnode,INTC);
        for(i=0;i<lNnode;i++)
          pMesh->noLG[i] = noLG[i];
/*..................................................................*/

        
/*... elementos locais*/
        size = lNel*maxNo;
        lEl  = (INT *) malloc(size*iSize); 
        ERRO_MALLOC(lEl,"lEl"  ,__LINE__,__FILE__,__func__);
        
        HccaAlloc(INT,m,mesh->elm.node  ,size        
                 ,"elnodeP"  ,_AD_);
        zero(mesh->elm.node         ,size ,INTC);
        
        getLocalEl(mesh0->elm.node,lEl
                  ,elLG           ,noGL
                  ,mesh0->elm.nen         
                  ,lNel           ,maxNo       );
        

        for(i=0;i<size;i++)
          mesh->elm.node[i] = lEl[i];

        free(lEl);
/*...................................................................*/

/*... mat locais*/
        size = lNel;
        nS   = (short *) malloc(size*sSize); 
        ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

        HccaAlloc(short,m,mesh->elm.mat  ,size        
                 ,"elMatP"     ,_AD_);
        zero(mesh->elm.mat          ,size         ,"short"  );

        sGetLocalV(mesh0->elm.mat,nS
                  ,elLG
                  ,lNel          ,1); 
        
        for(i=0;i<size;i++)
          mesh->elm.mat[i] = nS[i];

        free(nS);  
/*...................................................................*/

/*... nen locais*/
        size = lNel;
        nS   = (short *) malloc(size*sSize); 
        ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

        HccaAlloc(short,m,mesh->elm.nen  ,size        
                 ,"elNenP"     ,_AD_);
        zero(mesh->elm.nen          ,size         ,"short"  );

        sGetLocalV(mesh0->elm.nen,nS
                  ,elLG
                  ,lNel          ,1); 
        
        for(i=0;i<size;i++)
          mesh->elm.nen[i] = nS[i];

        free(nS);  
/*...................................................................*/

/*... geomType locais*/
        size = lNel;
        nS   = (short *) malloc(size*sSize); 
        ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

        HccaAlloc(short,m,mesh->elm.geomType  ,size        
                 ,"elGtP"      ,_AD_);
        zero(mesh->elm.geomType     ,size          ,"short"  );
        
        sGetLocalV(mesh0->elm.geomType,nS
                  ,elLG
                  ,size              ,1); 
        
        for(i=0;i<size;i++)
          mesh->elm.geomType[i] = nS[i];

        free(nS);  
/*...................................................................*/

/*... transporte e fluido*/
        if(ndfT[0] > 0 || ndfF > 0) {     
/*... eVel*/
          size = lNel*ndm;
          nD   = (DOUBLE *) malloc(size*dSize); 
          ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

          HccaAlloc(DOUBLE,m,mesh->elm.vel      ,size        
                   ,"eVelP"         ,_AD_);
          zero(mesh->elm.vel       ,size,DOUBLEC);

          dGetLocalV(mesh0->elm.vel    ,nD
                    ,elLG
                    ,lNel              ,ndm); 
        
          for(i=0;i<size;i++){
            mesh->elm.vel[i] = nD[i];
          }

          free(nD);  
/*...................................................................*/
        }
/*...................................................................*/

/*... */
        if(ndfD[0] > 0){
/*... faceRd1 locais*/
          size = lNel*(maxViz+1)*ndfD[0];
          HccaAlloc(short,m,mesh->elm.faceRd1  ,size        
                   ,"faceRd1P"      ,_AD_);
          zero(mesh->elm.faceRd1   ,size,"short"  );
/*...................................................................*/
  
/*... faceLoadsD1 locais*/
          size = lNel*(maxViz+1)*ndfD[0];
          HccaAlloc(short,m,mesh->elm.faceLoadD1  ,size        
                   ,"faceLoadsD1P"      ,_AD_);
          zero(mesh->elm.faceLoadD1   ,size,"short"  );
/*...................................................................*/
  
/*... eU0d1 locais*/
          size = lNel*ndfD[0];
          HccaAlloc(DOUBLE,m,mesh->elm.u0D1  ,size        
                   ,"eU0d1P"             ,_AD_);
          zero(mesh->elm.u0D1            ,size,DOUBLEC);
/*...................................................................*/
  
/*... eUd1 locais*/
          size = lNel*ndfD[0];
          HccaAlloc(DOUBLE,m,mesh->elm.uD1  ,size        
                   ,"eUd1P"             ,_AD_);
          zero(mesh->elm.uD1            ,size,DOUBLEC);
/*...................................................................*/
  
/*... density locais*/
          size = lNel*2;
          HccaAlloc(DOUBLE,m,mesh->elm.densityUd1  ,size        
                   ,"densityP"            ,_AD_);
          zero(mesh->elm.densityUd1,size,DOUBLEC);
/*...................................................................*/

/*... */
          comunicate2(mesh0->elm.faceRd1   ,mesh->elm.faceRd1
                     ,mesh0->elm.faceLoadD1,mesh->elm.faceLoadD1
                     ,mesh0->elm.u0D1      ,mesh->elm.u0D1
                     ,mesh0->elm.uD1       ,mesh->elm.uD1
                     ,mesh0->elm.densityUd1,mesh->elm.densityUd1
                     ,lNel                 ,elLG
                     ,maxViz               ,ndfD[0]
                     ,0                    ,1);
/*...................................................................*/
        }
/*...................................................................*/

/*... */
        if(ndfT[0] > 0){ 
/*... faceRt1 locais*/
          size = lNel*(maxViz+1)*ndfT[0];
          HccaAlloc(short,m,mesh->elm.faceRt1  ,size        
                   ,"faceRt1P"      ,_AD_);
          zero(mesh->elm.faceRt1   ,size,"short"  );
/*...................................................................*/

/*... faceLoadsT1 locais*/
          size = lNel*(maxViz+1)*ndfT[0];
          HccaAlloc(short,m,mesh->elm.faceLoadT1  ,size        
                   ,"faceLoadsT1P"      ,_AD_);
          zero(mesh->elm.faceLoadT1   ,size,"short"  );
/*...................................................................*/

/*... eU0t1 locais*/
          size = lNel*ndfT[0];
          HccaAlloc(DOUBLE,m,mesh->elm.u0T1  ,size        
                   ,"eU0t1P"             ,_AD_);
          zero(mesh->elm.u0T1            ,size,DOUBLEC);
/*...................................................................*/

/*... eUt1 locais*/
          size = lNel*ndfT[0];
          HccaAlloc(DOUBLE,m,mesh->elm.uT1  ,size        
                   ,"eUt1P"             ,_AD_);
          zero(mesh->elm.uT1            ,size,DOUBLEC);
/*...................................................................*/

/*... density locais*/
          size = lNel*2;
          HccaAlloc(DOUBLE,m,mesh->elm.densityUt1  ,size        
                   ,"densityT1P"          ,_AD_);
          zero(mesh->elm.densityUt1,size,DOUBLEC);
/*...................................................................*/

/*... */
          comunicate2(mesh0->elm.faceRt1   ,mesh->elm.faceRt1
                     ,mesh0->elm.faceLoadT1,mesh->elm.faceLoadT1
                     ,mesh0->elm.u0T1      ,mesh->elm.u0T1
                     ,mesh0->elm.uT1       ,mesh->elm.uT1
                     ,mesh0->elm.densityUt1,mesh->elm.densityUt1
                     ,lNel                 ,elLG
                     ,maxViz               ,ndfT[0]
                     ,0                    ,1);
/*...................................................................*/
        }  
/*...................................................................*/

/*... nelcon locais*/
        size = lNel*maxViz;
        lEl  = (INT *) malloc(size*iSize); 
        ERRO_MALLOC(lEl,"lEl"   ,__LINE__,__FILE__,__func__);

        HccaAlloc(INT,m,mesh->elm.adj.nelcon  ,size        
                 ,"adjP"      ,_AD_);
        zero(mesh->elm.adj.nelcon,size,INTC);


        getLocalAdj(mesh0->elm.adj.nelcon ,lEl
                   ,elLG                  ,elGL
                   ,mesh0->elm.adj.nViz 
                   ,lNel                  ,maxViz); 
        
        for(i=0;i<size;i++)
          mesh->elm.adj.nelcon[i] = lEl[i];

        free(lEl);  
/*...................................................................*/

/*... nViz locais*/
        size = lNel;
        nS   = (short *) malloc(size*sSize); 
        ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

        HccaAlloc(short,m,mesh->elm.adj.nViz  ,size        
                 ,"nVizP"             ,_AD_);
        zero(mesh->elm.adj.nViz  ,size      ,"short");

        sGetLocalV(mesh0->elm.adj.nViz  ,nS
                  ,elLG
                  ,size                 ,1); 
        
        for(i=0;i<size;i++)
          mesh->elm.adj.nViz[i] = nS[i];

        free(nS);  
/*...................................................................*/

/*... no locais*/
        size = lNnode*ndm;
        nD   = (DOUBLE *) malloc(size*dSize); 
        ERRO_MALLOC(nD,"nD"  ,__LINE__,__FILE__,__func__);
        
        HccaAlloc(DOUBLE,m,mesh->node.x,size 
                 ,"xnodeP",_AD_);   
        zero(mesh->node.x,size,DOUBLEC);
        
        dGetLocalV(mesh0->node.x,nD
                  ,noLG
                  ,lNnode      ,ndm); 

        for(i=0;i<size;i++)
          mesh->node.x[i] = nD[i];
  
        free(nD);  
/*...................................................................*/
      }
/*...................................................................*/

/*... comunicando o partionamento local para o processos escravos 
      (master)*/
      else{
        
        MPI_Send(&lNnode, 1, MPI_INT,nPart,nPart
                ,mpiVar.comm);
        MPI_Send(&lNel  , 1, MPI_INT,nPart,nPart
                ,mpiVar.comm);
        MPI_Send(&numelOv, 1, MPI_INT,nPart,nPart
                ,mpiVar.comm);
        MPI_Send(&nNodeOv, 1   , MPI_INT,nPart,nPart
                ,mpiVar.comm);
        MPI_Send(&nno1   , 1   , MPI_INT,nPart,nPart
                ,mpiVar.comm);

/*... particionamento*/     

/*... partViz(El)*/
        MPI_Send(&nVizPart   , 1          , MPI_SHORT,nPart,nPart
                ,mpiVar.comm);
        MPI_Send(maxVizPart  , nVizPart   , MPI_SHORT,nPart,nPart    
                ,mpiVar.comm);
/*... partViz(No)*/
        MPI_Send(&nVizPartNo , 1          , MPI_SHORT,nPart,nPart
                ,mpiVar.comm);
        MPI_Send(maxVizPartNo, nVizPartNo , MPI_SHORT,nPart,nPart    
                ,mpiVar.comm);
/*...................................................................*/

/*...*/
        MPI_Send(&nRcvs      , 1          , MPI_INT  ,nPart,nPart
                ,mpiVar.comm);
        MPI_Send(&nSends     , 1          , MPI_INT  ,nPart,nPart
                ,mpiVar.comm);
        MPI_Send(&nComNo     , 1          , MPI_INT  ,nPart,nPart
                ,mpiVar.comm);
/*...................................................................*/
 
/*... fMapEl(buffer de comunicacao)*/ 
        MPI_Send(fMap        ,nRcvs+nSends, MPI_INT  ,nPart,nPart
                ,mpiVar.comm);
/*... fMapNo(buffer de comunicacao)*/ 
        MPI_Send(fMapNo      ,nComNo      , MPI_INT  ,nPart,nPart
                ,mpiVar.comm);
/*..................................................................*/
        
/*... send e rcvs*/       
        kk = nVizPart + 1;
        MPI_Send(iaSends     ,kk     , MPI_INT  ,nPart,nPart
                ,mpiVar.comm);

        MPI_Send(iaRcvs      ,kk     , MPI_INT  ,nPart,nPart
                ,mpiVar.comm);

/*... iaComNo*/       
        kk = nVizPartNo + 1;
        MPI_Send(iaComNo    ,kk     , MPI_INT  ,nPart,nPart
                ,mpiVar.comm);
/*..................................................................*/

/*... elLG( numeracao local-global de elementos)*/
        MPI_Send(elLG                ,lNel   , MPI_INT  ,nPart,nPart
                ,mpiVar.comm);
/*...................................................................*/

/*... noLG( numeracao local-global de elementos)*/
        MPI_Send(noLG                ,lNnode  , MPI_INT  ,nPart,nPart
                ,mpiVar.comm);
/*...................................................................*/

/*... elementos locais*/
        size = maxNo*lNel;
        lEl  = (INT *) malloc(size*iSize); 
        ERRO_MALLOC(lEl,"lEl"  ,__LINE__,__FILE__,__func__);
        
        getLocalEl(mesh0->elm.node,lEl
                  ,elLG           ,noGL 
                  ,mesh0->elm.nen         
                  ,lNel           ,maxNo);
        
        MPI_Send(lEl   ,    size, MPI_INT,nPart,1    
                ,mpiVar.comm);
        
        free(lEl);
/*...................................................................*/

/*... mat locais*/
        size = lNel;
        nS   = (short *) malloc(size*sSize); 
        ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

        sGetLocalV(mesh0->elm.mat,nS
                  ,elLG
                  ,size          ,1); 
        
        MPI_Send(nS    ,       size, MPI_SHORT,nPart,2    
                ,mpiVar.comm);
        
        free(nS);  
/*...................................................................*/

/*... nen locais*/
        size = lNel;
        nS   = (short *) malloc(size*sSize); 
        ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

        sGetLocalV(mesh0->elm.nen,nS
                  ,elLG
                  ,size          ,1); 
        
        MPI_Send(nS    ,      size, MPI_SHORT,nPart,3    
                ,mpiVar.comm);
        
        free(nS);
/*...................................................................*/

/*... geomType locais*/
        size = lNel;
        nS   = (short *) malloc(size*sSize); 
        ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

        sGetLocalV(mesh0->elm.geomType,nS
                  ,elLG
                  ,size               ,1); 
        
        MPI_Send(nS    ,      size, MPI_SHORT,nPart,4     
                ,mpiVar.comm);
        
        free(nS);
/*...................................................................*/

/*... transporte e fluido*/
        if(ndfT[0] > 0 || ndfF > 0) {     
/*... eVel*/
          size = lNel*ndm;
          nD   = (DOUBLE *) malloc(size*dSize); 
          ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

          dGetLocalV(mesh0->elm.vel    ,nD
                    ,elLG
                    ,lNel              ,ndm); 
          
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,5    
                  ,mpiVar.comm);
        
          free(nD);  
/*...................................................................*/
        }
/*...................................................................*/

/*... */
        if( ndfD[0] > 0){
          comunicate2(mesh0->elm.faceRd1   ,mesh->elm.faceRd1
                     ,mesh0->elm.faceLoadD1,mesh->elm.faceLoadD1
                     ,mesh0->elm.u0D1      ,mesh->elm.u0D1
                     ,mesh0->elm.uD1       ,mesh->elm.uD1
                     ,mesh0->elm.densityUd1,mesh->elm.densityUd1
                     ,lNel                 ,elLG
                     ,maxViz               ,ndfD[0]
                     ,nPart                ,2);
/*...................................................................*/
        }
/*...................................................................*/

/*... */
        if( ndfT[0] > 0){
          comunicate2(mesh0->elm.faceRt1   ,mesh->elm.faceRt1
                     ,mesh0->elm.faceLoadT1,mesh->elm.faceLoadT1
                     ,mesh0->elm.u0T1      ,mesh->elm.u0T1
                     ,mesh0->elm.uT1       ,mesh->elm.uT1
                     ,mesh0->elm.densityUt1,mesh->elm.densityUt1
                     ,lNel                 ,elLG
                     ,maxViz               ,ndfT[0]
                     ,nPart                ,2);
/*...................................................................*/
        }
/*...................................................................*/

/*... nelcon locais*/
        size = lNel*maxViz;
        lEl   = (INT *) malloc(size*iSize);

        ERRO_MALLOC(lEl,"lEl"   ,__LINE__,__FILE__,__func__);

        getLocalAdj(mesh0->elm.adj.nelcon ,lEl
                   ,elLG                  ,elGL
                   ,mesh0->elm.adj.nViz 
                   ,lNel                  ,maxViz); 
        
        MPI_Send(lEl   ,       size, MPI_INT,nPart,10    
                ,mpiVar.comm);
        
        free(lEl);  
/*...................................................................*/

/*... nViz locais*/
        size = lNel;
        nS   = (short *) malloc(size*sSize); 
        ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

        sGetLocalV(mesh0->elm.adj.nViz  ,nS
                  ,elLG
                  ,size                 ,1); 
        
        MPI_Send(nS    ,       size, MPI_SHORT,nPart,11    
                ,mpiVar.comm);
        
        free(nS);  
/*...................................................................*/

/*... no locais*/
        size = lNnode*ndm;
        nD   = (DOUBLE *) malloc(size*dSize); 
        ERRO_MALLOC(nD,"nD"  ,__LINE__,__FILE__,__func__);
        
        dGetLocalV(mesh0->node.x,nD
                  ,noLG
                  ,lNnode      ,ndm); 
        
        MPI_Send(nD    ,size , MPI_DOUBLE,nPart,12   
                ,mpiVar.comm);
        
        free(nD);
/*...................................................................*/
      } 
/*...................................................................*/
    }
/*...................................................................*/
    
/*... recebendo o partionamento local do processo master (escravos)*/
    if(nPart == myId && myId){
      
      MPI_Recv(&lNnode    , 1, MPI_INT,0,nPart,
               mpiVar.comm,MPI_STATUS_IGNORE); 
      MPI_Recv(&lNel      , 1, MPI_INT,0,nPart,
               mpiVar.comm,MPI_STATUS_IGNORE); 
      MPI_Recv(&numelOv   , 1, MPI_INT,0,nPart,
               mpiVar.comm,MPI_STATUS_IGNORE); 
      MPI_Recv(&nNodeOv   , 1, MPI_INT,0,nPart,
               mpiVar.comm,MPI_STATUS_IGNORE); 
      MPI_Recv(&nno1      , 1, MPI_INT,0,nPart,
               mpiVar.comm,MPI_STATUS_IGNORE); 
     
      mesh->nnode    = lNnode;
      mesh->numel    = lNel;
      mesh->numelNov = lNel-numelOv;
      mesh->nnodeOv  = nNodeOv ;
      mesh->nnodeNov = lNnode-nNodeOv;
      pMesh->nno1    = nno1;

/*... particionamento*/     

/*... partViz(El)*/
      MPI_Recv(&nVizPart  , 1, MPI_SHORT,0,nPart,
               mpiVar.comm,MPI_STATUS_IGNORE); 
      
      pMesh->iEl.nVizPart = nVizPart;

      HccaAlloc(short,m,pMesh->iEl.vizPart , nVizPart   
               ,"nVizPart"             ,_AD_);
      zero(pMesh->iEl.vizPart               ,nVizPart,"short");
      
      MPI_Recv(pMesh->iEl.vizPart,nVizPart, MPI_SHORT,0,nPart,
               mpiVar.comm,MPI_STATUS_IGNORE);
/*... partViz(No)*/
      MPI_Recv(&nVizPartNo, 1, MPI_SHORT,0,nPart,
               mpiVar.comm,MPI_STATUS_IGNORE); 
      
      pMesh->iNo.nVizPart = nVizPartNo;

      HccaAlloc(short,m,pMesh->iNo.vizPart   , nVizPartNo   
               ,"nVizPartNo"           ,_AD_);
      zero(pMesh->iNo.vizPart               ,nVizPartNo,"short");
      
      MPI_Recv(pMesh->iNo.vizPart,nVizPartNo, MPI_SHORT,0,nPart,
               mpiVar.comm,MPI_STATUS_IGNORE);
 
/*..................................................................*/

/*...*/     
      MPI_Recv(&nRcvs        ,       1, MPI_INT,0,nPart 
              ,mpiVar.comm,MPI_STATUS_IGNORE); 
      MPI_Recv(&nSends       ,       1, MPI_INT,0,nPart 
              ,mpiVar.comm,MPI_STATUS_IGNORE); 
      MPI_Recv(&nComNo       ,       1, MPI_INT,0,nPart 
              ,mpiVar.comm,MPI_STATUS_IGNORE); 

      pMesh->iEl.nRcvs  = nRcvs;
      pMesh->iEl.nSends = nSends;
      pMesh->iNo.nCom   = nComNo;
/*..................................................................*/
      
      
/*... fMapEl(buffer de comunicacao)*/ 
      size = nRcvs+nSends;
      HccaAlloc(INT  ,m,pMesh->iEl.fMap    , size   
               ,"fMap"                 ,_AD_);
      zero(pMesh->iEl.fMap,size,INTC);

      MPI_Recv(pMesh->iEl.fMap ,    size, MPI_INT,0,nPart 
              ,mpiVar.comm,MPI_STATUS_IGNORE); 
/*... fMapNo (buffer de comunicacao)*/ 
      size = nComNo;
      HccaAlloc(INT  ,m,pMesh->iNo.fMap    , size   
               ,"fMapNo"               ,_AD_);
      zero(pMesh->iNo.fMap,size,INTC);

      MPI_Recv(pMesh->iNo.fMap ,    size, MPI_INT,0,nPart 
              ,mpiVar.comm,MPI_STATUS_IGNORE); 
/*...................................................................*/
      
      
/*... send e rcvs*/       
      kk = nVizPart + 1;

      HccaAlloc(INT   ,m   ,pMesh->iEl.iaSends , kk        
                 ,"iaSends",_AD_);
      zero(pMesh->iEl.iaSends   , kk,INTC);
        
      HccaAlloc(INT   ,m   ,pMesh->iEl.iaRcvs , kk       
                 ,"iaRcvs",_AD_);
      zero(pMesh->iEl.iaRcvs   , kk,INTC);
        
      MPI_Recv(pMesh->iEl.iaSends,  kk, MPI_INT,0,nPart 
              ,mpiVar.comm,MPI_STATUS_IGNORE); 
      
      MPI_Recv(pMesh->iEl.iaRcvs,  kk, MPI_INT,0,nPart 
              ,mpiVar.comm,MPI_STATUS_IGNORE); 
/*... iaComNo    */       
      kk = nVizPartNo + 1;

      HccaAlloc(INT   ,m   ,pMesh->iNo.iaComNo , kk        
                 ,"iaSendsNo",_AD_);
      zero(pMesh->iNo.iaComNo   , kk,INTC);
        
      MPI_Recv(pMesh->iNo.iaComNo,  kk, MPI_INT,0,nPart 
              ,mpiVar.comm,MPI_STATUS_IGNORE); 
/*...................................................................*/

/*... elLG( numeracao local-global de elementos)*/
      HccaAlloc(INT   ,m   ,pMesh->elLG , lNel        
                 ,"elLG"   ,_AD_);
      zero(pMesh->elLG   , lNel,INTC);
      MPI_Recv(pMesh->elLG,  lNel, MPI_INT,0,nPart 
              ,mpiVar.comm,MPI_STATUS_IGNORE); 
/*..................................................................*/

/*... noLG( numeracao local-global de nos)*/
      HccaAlloc(INT   ,m   ,pMesh->noLG , lNnode       
                 ,"noLG"   ,_AD_);
      zero(pMesh->noLG   , lNnode,INTC);
      MPI_Recv(pMesh->noLG,  lNnode, MPI_INT,0,nPart 
              ,mpiVar.comm,MPI_STATUS_IGNORE); 
/*..................................................................*/

/*... elementos locais*/
      size = lNel*maxNo;
      
      HccaAlloc(INT,m,mesh->elm.node  ,size        
               ,"elnodeP"  ,_AD_);
      zero(mesh->elm.node         ,size ,INTC);
      
      MPI_Recv(mesh->elm.node,size, MPI_INT,0,1    
              ,mpiVar.comm,MPI_STATUS_IGNORE);
      
/*...................................................................*/

/*... mat locais*/
      size = lNel;

      HccaAlloc(short,m,mesh->elm.mat  ,size        
               ,"elMatP"     ,_AD_);
      zero(mesh->elm.mat          ,size         ,"short"  );
      
      MPI_Recv(mesh->elm.mat,lNel, MPI_SHORT,0,2    
              ,mpiVar.comm,MPI_STATUS_IGNORE);
      
/*...................................................................*/

/*... nen locais*/
      size = lNel;

      HccaAlloc(short,m,mesh->elm.nen  ,size        
               ,"elNenP"     ,_AD_);
      zero(mesh->elm.nen          ,size         ,"short"  );
      
      MPI_Recv(mesh->elm.nen,lNel, MPI_SHORT,0,3    
              ,mpiVar.comm,MPI_STATUS_IGNORE);
      
/*...................................................................*/

/*... geomType locais*/
      size = lNel;

      HccaAlloc(short,m,mesh->elm.geomType  ,size        
               ,"elGtP"     ,_AD_);
      zero(mesh->elm.geomType     ,size          ,"short"  );
      
      MPI_Recv(mesh->elm.geomType,lNel, MPI_SHORT,0,4    
              ,mpiVar.comm,MPI_STATUS_IGNORE);
      
/*...................................................................*/

/*... transporte e fluido*/
        if(ndfT[0] > 0 || ndfF > 0) {     
/*... eVel*/
          size = lNel*ndm;
          
          HccaAlloc(DOUBLE,m,mesh->elm.vel      ,size        
               ,"evelP"   ,_AD_);
          zero(mesh->elm.vel       ,size,DOUBLEC);

          
          MPI_Recv(mesh->elm.vel,       size, MPI_DOUBLE,0,5    
                  ,mpiVar.comm,MPI_STATUS_IGNORE);
        
          free(nD);  
/*...................................................................*/
        }
/*...................................................................*/

/*...*/
      if(ndfD[0] > 0) {
/*... faceRd1 locais*/
        size = lNel*(maxViz+1)*ndfD[0];
        HccaAlloc(short,m,mesh->elm.faceRd1  ,size        
               ,"faceRd1P"  ,_AD_);
        zero(mesh->elm.faceRd1   ,size,"short"  );
/*...................................................................*/

/*... faceLoadsD1 locais*/
        size = lNel*(maxViz+1)*ndfD[0];
        HccaAlloc(short,m,mesh->elm.faceLoadD1 ,size        
                 ,"faceLoadsD1P"  ,_AD_);
        zero(mesh->elm.faceLoadD1   ,size,"short"  );
/*...................................................................*/

/*... eU0d1 locais*/
        size = lNel*ndfD[0];
        HccaAlloc(DOUBLE,m,mesh->elm.u0D1        ,size        
                 ,"eU0d1P"         ,_AD_);
        zero(mesh->elm.u0D1          ,size,DOUBLEC);;
/*...................................................................*/

/*... eUd1 locais*/
        size = lNel*ndfD[0];
        HccaAlloc(DOUBLE,m,mesh->elm.uD1        ,size        
                 ,"eUd1P"         ,_AD_);
        zero(mesh->elm.uD1          ,size,DOUBLEC);;
/*...................................................................*/

/*... density locais*/
        size = lNel*2;
        HccaAlloc(DOUBLE,m,mesh->elm.densityUd1  ,size        
                 ,"densityP"            ,_AD_);
        zero(mesh->elm.densityUd1,size,DOUBLEC);
/*...................................................................*/

/*...*/
        comunicate2(mesh0->elm.faceRd1   ,mesh->elm.faceRd1
                   ,mesh0->elm.faceLoadD1,mesh->elm.faceLoadD1
                   ,mesh0->elm.u0D1      ,mesh->elm.u0D1
                   ,mesh0->elm.uD1       ,mesh->elm.uD1
                   ,mesh0->elm.densityUd1,mesh->elm.densityUd1
                   ,lNel                 ,elLG
                   ,maxViz               ,ndfD[0]
                   ,nPart                ,3);
/*...................................................................*/
      }
/*...................................................................*/

/*... faceRt1 locais*/
      if(ndfT[0] > 0) {
        size = lNel*(maxViz+1)*ndfT[0];
        HccaAlloc(short,m,mesh->elm.faceRt1  ,size        
               ,"faceRt1P"  ,_AD_);
        zero(mesh->elm.faceRt1   ,size,"short"  );
/*...................................................................*/

/*... faceLoadsT1 locais*/ 
        size = lNel*(maxViz+1)*ndfT[0];
        HccaAlloc(short,m,mesh->elm.faceLoadT1 ,size        
                 ,"faceLoadsT1P"  ,_AD_);
        zero(mesh->elm.faceLoadT1   ,size,"short"  );
/*...................................................................*/

/*... eU0t1 locais*/
        size = lNel*ndfT[0];
        HccaAlloc(DOUBLE,m,mesh->elm.u0T1        ,size        
                 ,"eU0t1P"         ,_AD_);
        zero(mesh->elm.u0T1          ,size,DOUBLEC);;
/*...................................................................*/

/*... eUt1 locais*/
        size = lNel*ndfT[0];
        HccaAlloc(DOUBLE,m,mesh->elm.uT1        ,size        
                 ,"eUt1P"         ,_AD_);
        zero(mesh->elm.uT1          ,size,DOUBLEC);;
/*...................................................................*/

/*... densityt1 locais*/
        size = lNel*2;
        HccaAlloc(DOUBLE,m,mesh->elm.densityUt1  ,size        
                 ,"densityT1P"            ,_AD_);
        zero(mesh->elm.densityUt1,size,DOUBLEC);
/*...................................................................*/

/*...*/
        comunicate2(mesh0->elm.faceRt1   ,mesh->elm.faceRt1
                   ,mesh0->elm.faceLoadT1,mesh->elm.faceLoadT1
                   ,mesh0->elm.u0T1      ,mesh->elm.u0T1
                   ,mesh0->elm.uT1       ,mesh->elm.uT1
                   ,mesh0->elm.densityUt1,mesh->elm.densityUt1
                   ,lNel                 ,elLG
                   ,maxViz               ,ndfT[0]
                   ,nPart                ,3);
/*...................................................................*/
      }
/*...................................................................*/

/*... nelvon locais*/
      size = lNel*maxViz;

      HccaAlloc(INT,m,mesh->elm.adj.nelcon  ,size        
                 ,"adjP"      ,_AD_);
      zero(mesh->elm.adj.nelcon,size,INTC);
      
      MPI_Recv( mesh->elm.adj.nelcon,size, MPI_INT,0,10   
              ,mpiVar.comm,MPI_STATUS_IGNORE);
      
/*...................................................................*/

/*... nViz locais*/           
      size = lNel;

      HccaAlloc(short,m,mesh->elm.adj.nViz ,size        
               ,"nVizP"         ,_AD_);
      zero(mesh->elm.adj.nViz  ,size      ,"short");
      
      MPI_Recv(mesh->elm.adj.nViz,size, MPI_SHORT,0,11   
              ,mpiVar.comm,MPI_STATUS_IGNORE);
      
/*...................................................................*/

/*... no locais*/
      size = lNnode*ndm;
      
      HccaAlloc(DOUBLE,m,mesh->node.x  ,size        
               ,"xnodeP",_AD_);   
      zero(mesh->node.x,size,DOUBLEC);
      
      MPI_Recv(mesh->node.x    ,size, MPI_DOUBLE,0,12   
              ,mpiVar.comm,MPI_STATUS_IGNORE);
      
/*...................................................................*/
    }
/*...................................................................*/
    free(elLG);
    free(noLG);
  }
/*...................................................................*/
  
  
/*... desalocando varaiveis auxiliares*/  
  if(!mpiVar.myId){
    free(mesh0->noIncid.incid); 
    free(mesh0->noIncid.nincid); 
    free(aux1);
    free(aux2);
    free(aux3);
    free(noGL);
    free(elGL);
    free(maxVizPart);
    free(fMap);
    free(fMapNo);
    free(iaRcvs); 
    free(iaSends);
    free(iaComNo); 
  }
/*...................................................................*/
  
  lNel     = mesh->numel;
  lNnode   = mesh->nnode;
  
/*... centroide */
  HccaAlloc(DOUBLE,m,mesh->elm.geom.cc ,lNel*ndm,"elCcP"    ,_AD_);
/*... vetor que une os centroides dos elementos */
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.ksi
         ,lNel*ndm*maxViz,"elksiP"  ,_AD_);
/*... modulo do vetor que une os centroides dos elementos */
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.mksi
        ,lNel*maxViz     ,"elmksiP",_AD_);
/*... vetor paralelo a face da celula */
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.eta
         ,lNel*ndm*maxViz,"eletaP"  ,_AD_);
/*... modulo do vetor paralelo a face da celula */
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.fArea
        ,lNel*maxViz     ,"elfAreaP",_AD_);
/*... volume da celula*/                           
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.volume
        ,lNel            ,"elVolP",_AD_);
/*... vetor normal a face da celula*/                           
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.normal
         ,lNel*maxViz*ndm       ,"elnormP",_AD_);
/*... ponto medio da face*/                           
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.xm
         ,lNel*maxViz*ndm       ,"elxmP",_AD_);
/*... vetor que une o centroide ao ponto medio*/                           
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.xmcc
         ,lNel*maxViz*ndm       ,"elxmccP",_AD_);
/*... vetor entre o ponto medio da face e ponto de intercao 
      entre a linha entre os centroides e a face*/         
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.vSkew  
         ,lNel*ndm*maxViz         ,"elvSkewP" ,_AD_);
/*... distancia entro o ponto medio da face e ponto de intercao 
      entre a linha entre os centroides e a face*/         
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.mvSkew  
         ,lNel*maxViz           ,"elmvSkewP" ,_AD_);
/*... distancia entro o ponto medio da face e ponto de intercao 
      entre a linha entre os centroides e a face*/         
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.dcca 
         ,lNel*maxViz           ,"eldccaP",_AD_);
/*... zerando os variavies*/
  zero(mesh->elm.geom.cc      ,lNel*ndm       ,DOUBLEC);
  zero(mesh->elm.geom.ksi     ,lNel*ndm*maxViz,DOUBLEC);
  zero(mesh->elm.geom.mksi    ,lNel*maxViz    ,DOUBLEC);
  zero(mesh->elm.geom.eta     ,lNel*ndm*maxViz,DOUBLEC);
  zero(mesh->elm.geom.fArea   ,lNel*maxViz    ,DOUBLEC);
  zero(mesh->elm.geom.volume  ,lNel           ,DOUBLEC);
  zero(mesh->elm.geom.normal  ,lNel*ndm*maxViz,DOUBLEC);
  zero(mesh->elm.geom.xm      ,lNel*ndm*maxViz,DOUBLEC);
  zero(mesh->elm.geom.xmcc    ,lNel*ndm*maxViz,DOUBLEC);
  zero(mesh->elm.geom.mvSkew  ,lNel*maxViz    ,DOUBLEC);
  zero(mesh->elm.geom.mvSkew  ,lNel*ndm*maxViz,DOUBLEC);
  zero(mesh->elm.geom.dcca    ,lNel*maxViz    ,DOUBLEC);
/*...................................................................*/

/*... transporte e fluido*/
  if(mesh->ndfT[0] > 0 || mesh->ndfF > 0) {     
/*... nVel*/
     HccaAlloc(DOUBLE,m,mesh->node.vel 
              ,lNnode*ndm    ,"nVelP"             ,_AD_);
     zero(mesh->node.vel     ,lNnode*ndm          ,DOUBLEC);
/*...................................................................*/
  }
/*...................................................................*/

/*... problema de difusao pura*/
  if(ndfD[0] > 0){
/*... uD1*/
    HccaAlloc(DOUBLE,m,mesh->node.uD1 
              ,lNnode*ndfD[0] ,"nUd1P"              ,_AD_);
/*... gradU1*/
    HccaAlloc(DOUBLE,m,mesh->node.gradUd1  
              ,lNnode*ndfD[0]*ndm ,"nGradUd1P"      ,_AD_);
    HccaAlloc(DOUBLE,m,mesh->elm.gradUd1 
              ,lNel*ndm*ndfD[0],"eTGradUd1P"     ,_AD_);
/*... rCell*/
    HccaAlloc(DOUBLE,m,mesh->elm.rCellUd1  
             ,lNel*ndm*ndfD[0],"rCellUd1P"      ,_AD_);
/*...................................................................*/
    zero(mesh->node.uD1      ,lNnode*ndfD[0]       ,DOUBLEC);
    zero(mesh->node.gradUd1  ,lNnode*ndm*ndfD[0]   ,DOUBLEC);
    zero(mesh->elm.gradUd1   ,lNel*ndm*ndfD[0]     ,DOUBLEC);
    zero(mesh->elm.rCellUd1  ,lNel*ndm*ndfD[0]     ,DOUBLEC);
  }
/*...................................................................*/

/*... problema de transporte*/
  if(ndfT[0] > 0){
/*... uD1*/
    HccaAlloc(DOUBLE,m,mesh->node.uT1 
              ,lNnode*ndfT[0] ,"nUt1P"              ,_AD_);
/*... gradU1*/
    HccaAlloc(DOUBLE,m,mesh->node.gradUt1  
              ,lNnode*ndfT[0]*ndm ,"nGradUt1P"      ,_AD_);
    HccaAlloc(DOUBLE,m,mesh->elm.gradUt1 
              ,lNel*ndm*ndfT[0],"eTGradUt1P"     ,_AD_);
/*... rCell*/
    HccaAlloc(DOUBLE,m,mesh->elm.rCellUt1  
             ,lNel*ndm*ndfT[0],"rCellUt1P"      ,_AD_);
/*...................................................................*/
    zero(mesh->node.uT1      ,lNnode*ndfT[0]       ,DOUBLEC);
    zero(mesh->node.gradUt1  ,lNnode*ndm*ndfT[0]   ,DOUBLEC);
    zero(mesh->elm.gradUt1   ,lNel*ndm*ndfT[0]     ,DOUBLEC);
    zero(mesh->elm.rCellUt1  ,lNel*ndm*ndfT[0]     ,DOUBLEC);
  }
/*...................................................................*/

/*...*/
  maxNdf = 0;
  maxNdf = max(ndm    ,maxNdf);
  maxNdf = max(ndfD[0],maxNdf);
  maxNdf = max(ndfD[1],maxNdf);
  maxNdf = max(ndfD[2],maxNdf);
  maxNdf = max(ndfT[0],maxNdf);
  maxNdf = max(ndfT[1],maxNdf);
  maxNdf = max(ndfT[2],maxNdf);
/*...................................................................*/

/*... buffer de comunicacao para variaveis de elementos*/
  size = (pMesh->iEl.nRcvs+pMesh->iEl.nSends)*ndm*maxNdf;
  HccaAlloc(DOUBLE,m,pMesh->iEl.xb 
           ,size ,"xBufferMpi"              ,false);
/*... buffer de comunicacao para variaveis de nos*/
  size = 2.0e0*pMesh->iNo.nCom*ndm*maxNdf;
  HccaAlloc(DOUBLE,m,pMesh->iNo.xb 
           ,size ,"xBufferMpiNo"            ,false);
  HccaAlloc(INT   ,m,pMesh->iNo.xi 
           ,size ,"xiBufferMpiNo"           ,false);
/*...................................................................*/
#endif
}
/*********************************************************************/ 

/*********************************************************************/ 
void comunicate2(short *m0faceR     ,short *faceR
                ,short *m0faceL     ,short *faceL
                ,DOUBLE *m0u0       ,DOUBLE *u0 
                ,DOUBLE *m0u        ,DOUBLE *u  
                ,DOUBLE *m0density  ,DOUBLE *density  
                ,INT const lNel     ,INT *elLG
                ,short const maxViz ,short const ndf
                ,short const nPart  ,short const iCod)
{
  INT i,size;
  DOUBLE *nD=NULL;
  short  *nS=NULL;
  short sSize=sizeof(short); 
//  short iSize=sizeof(INT); 
  short dSize=sizeof(DOUBLE); 

  switch(iCod){
    case 1: 
/*... faceR locais*/
      size = lNel*(maxViz+1)*ndf;
      nS   = (short *) malloc(size*sSize); 
      ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

      sGetLocalV(m0faceR           ,nS
                ,elLG
                ,lNel              ,(maxViz+1)*ndf); 
      
      for(i=0;i<size;i++) 
        faceR[i] = nS[i];
      
      free(nS);  
/*...................................................................*/

/*... faceL locais*/
      size = lNel*(maxViz+1)*ndf;
      nS   = (short *) malloc(size*sSize); 
      ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

      sGetLocalV(m0faceL,nS
                ,elLG
                ,lNel   ,(maxViz+1)*ndf); 
        
      for(i=0;i<size;i++)
        faceL[i] = nS[i];

      free(nS);
/*...................................................................*/

/*... eU0 locais*/
      size = lNel*ndf;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

      dGetLocalV(m0u0              ,nD
                ,elLG
                ,lNel              ,ndf); 
        
      for(i=0;i<size;i++)
        u0[i] = nD[i];

      free(nD);
/*...................................................................*/

/*... eU locais*/
      size = lNel*ndf;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

      dGetLocalV(m0u    ,nD
                ,elLG
                ,lNel   ,ndf); 
        
      for(i=0;i<size;i++)
        u[i] = nD[i];

      free(nD);
/*...................................................................*/

/*... density locais*/
      size = lNel*2;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

      dGetLocalV(m0density  ,nD
                ,elLG
                ,lNel       ,2); 
        
      for(i=0;i<size;i++)
        density[i] = nD[i];

      free(nD);
/*...................................................................*/
   break;

   case 2:

/*... faceR locais*/
     size = lNel*(maxViz+1)*ndf;
     nS   = (short *) malloc(size*sSize); 
     ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

     sGetLocalV(m0faceR,nS
               ,elLG
               ,lNel   ,(maxViz+1)*ndf); 
        
     MPI_Send(nS    ,       size, MPI_SHORT,nPart,5    
             ,mpiVar.comm);
        
     free(nS);  
/*...................................................................*/

/*... faceL locais*/
     size = lNel*(maxViz+1)*ndf;
     nS   = (short *) malloc(size*sSize); 
     ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

     sGetLocalV(m0faceL,nS
               ,elLG
               ,lNel   ,(maxViz+1)*ndf); 
        
     MPI_Send(nS    ,       size, MPI_SHORT,nPart,6    
             ,mpiVar.comm);
        
     free(nS);  
/*...................................................................*/

/*... eU0 locais*/
     size = lNel*ndf;
     nD   = (DOUBLE *) malloc(size*dSize); 
     ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

     dGetLocalV(m0u0     ,nD
               ,elLG
               ,lNel     ,ndf); 
        
     MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,7    
             ,mpiVar.comm);
        
     free(nD);  
/*...................................................................*/

/*... eU locais*/
     size = lNel*ndf;
     nD   = (DOUBLE *) malloc(size*dSize); 
     ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

     dGetLocalV(m0u      ,nD
               ,elLG
               ,lNel     ,ndf); 
        
     MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,8    
             ,mpiVar.comm);
        
     free(nD);  
/*...................................................................*/

/*... density locais*/
     size = lNel*2;
     nD   = (DOUBLE *) malloc(size*dSize); 
     ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

     dGetLocalV(m0density,nD
               ,elLG
               ,lNel     ,2); 
        
     MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,9    
               ,mpiVar.comm);
        
     free(nD);  
/*...................................................................*/
   break; 
/*...................................................................*/

/*...*/
   case 3:
   
/*... faceR locais*/
     size = lNel*(maxViz+1)*ndf;
     MPI_Recv(faceR,size, MPI_SHORT,0,5    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/
     
/*... faceL locais*/
     size = lNel*(maxViz+1)*ndf;
     MPI_Recv(faceL,size, MPI_SHORT,0,6    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... eU0 locais*/
     size = lNel*ndf;
     MPI_Recv(u0  ,size, MPI_DOUBLE,0,7    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/
     
/*... eU locais*/
     size = lNel*ndf;
     MPI_Recv(u,size, MPI_DOUBLE,0,8    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/
     
/*... density locais*/
     size = lNel*2;
     MPI_Recv(density,size, MPI_DOUBLE,0,9    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

   break; 
/*...................................................................*/


  }
/*...................................................................*/

}
/*********************************************************************/ 

