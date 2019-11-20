#include<ParallelMpi.h>
/********************************************************************* 
 * Data de criacao    : 16/08/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * allocFuild :                                                     * 
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
static void allocFuild(Memoria *m        , Mesh *mesh
                     , short const maxViz, short const lNel
                     , short const ndfFt)
{
  INT size;

/*... faceRvel locais*/
  size = lNel*(maxViz+1);
  HccaAlloc(short,m,mesh->elm.faceRvel  ,size        
       ,"faceRvelP"  ,_AD_);
  zero(mesh->elm.faceRvel   ,size,"short"  );
/*...................................................................*/

/*... faceRpres locais*/
  size = lNel*(maxViz+1);
  HccaAlloc(short,m,mesh->elm.faceRpres ,size        
         ,"faceRpresP"  ,_AD_);
  zero(mesh->elm.faceRpres  ,size,"short"  );
/*...................................................................*/

/*... pres0 locais*/
  size = lNel;
  HccaAlloc(DOUBLE,m,mesh->elm.pressure0     ,size        
           ,"ePres0P"           ,_AD_);
  zero(mesh->elm.pressure0,size,DOUBLEC);
/*...................................................................*/

/*... pres locais*/
  size = lNel;
  HccaAlloc(DOUBLE,m,mesh->elm.pressure     ,size        
           ,"ePresP",_AD_);
  zero(mesh->elm.pressure,size,DOUBLEC);
/*...................................................................*/

/*... energia*/
  if(ndfFt > 0)
  {
/*... faceRenergy locais*/
    size = lNel*(maxViz+1);
    HccaAlloc(short,m,mesh->elm.faceRenergy ,size        
             ,"faceRenP"  ,_AD_);
    zero(mesh->elm.faceRenergy,size,"short"  );
/*...................................................................*/

/*...*/
    size = lNel*(maxViz+1);
    HccaAlloc(short, m, mesh->elm.faceRrho    
                , size, "faceRrhoP", _AD_);
    zero(mesh->elm.faceRrho, size, "short");
/*...................................................................*/

/*... Energy0 locais*/
    size = lNel;
    HccaAlloc(DOUBLE,m,mesh->elm.energy0     ,size        
             ,"eEnergy0P"         ,_AD_);
    zero(mesh->elm.energy0         ,size,DOUBLEC);
/*...................................................................*/

/*... Energy locais*/
    size = lNel;
    HccaAlloc(DOUBLE,m,mesh->elm.energy     ,size        
             ,"eEnergyP"         ,_AD_);
    zero(mesh->elm.energy         ,size,DOUBLEC);
/*...................................................................*/

/*... Temp0 locais*/
    size = lNel;
    HccaAlloc(DOUBLE,m,mesh->elm.temp0       ,size        
             ,"eTemp0P"           ,_AD_);
    zero(mesh->elm.temp0           ,size,DOUBLEC);
/*...................................................................*/

/*... Temp locais*/
    size = lNel;
    HccaAlloc(DOUBLE,m,mesh->elm.temp       ,size        
             ,"eTempP"           ,_AD_);
    zero(mesh->elm.temp           ,size,DOUBLEC);
/*...................................................................*/

/*... densitytFluid locais*/
    size = lNel;
    HccaAlloc(DOUBLE,m,mesh->elm.densityFluid.t00  ,size        
           ,"densityFluidP00"            ,_AD_);
    HccaAlloc(DOUBLE,m,mesh->elm.densityFluid.t0  ,size        
           ,"densityFluidP0"            ,_AD_);
    HccaAlloc(DOUBLE,m,mesh->elm.densityFluid.t  ,size        
           ,"densityFluidP"            ,_AD_);
    zero(mesh->elm.densityFluid.t00,size,DOUBLEC);
    zero(mesh->elm.densityFluid.t0 ,size,DOUBLEC);
    zero(mesh->elm.densityFluid.t  ,size,DOUBLEC);   
/*...................................................................*/

/*... calor especifico locais*/
    size = lNel;
    HccaAlloc(DOUBLE,m,mesh->elm.specificHeat.t00 ,size        
           ,"sHeatP00",_AD_);
    HccaAlloc(DOUBLE,m,mesh->elm.specificHeat.t0  ,size         
           ,"sHeatP0",_AD_);
    HccaAlloc(DOUBLE,m,mesh->elm.specificHeat.t   ,size        
           ,"sHeatP",_AD_);
    zero(mesh->elm.specificHeat.t00,size,DOUBLEC);
    zero(mesh->elm.specificHeat.t0 ,size,DOUBLEC);
    zero(mesh->elm.specificHeat.t  ,size,DOUBLEC);
/*...................................................................*/

/*... calor especifico locais*/
    size = lNel;
    HccaAlloc(DOUBLE,m,mesh->elm.dViscosity  ,size        
           ,"dViscP",_AD_);
    zero(mesh->elm.dViscosity,size,DOUBLEC);
/*...................................................................*/

/*... calor especifico locais*/
    size = lNel;
    HccaAlloc(DOUBLE,m,mesh->elm.tConductivity  ,size        
           ,"tConductivityP",_AD_);
    zero(mesh->elm.tConductivity,size,DOUBLEC);
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 16/08/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * allocComb :                                                       * 
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
static void allocComb(Memoria *m          , Mesh *mesh
                     , short const maxViz , short const lNel
                     , short const ndfComb, short const nSp)
{
  INT size;

/*... faceResZcomb locais*/
  size = lNel*(maxViz+1);
  HccaAlloc(short,m,mesh->elm.faceResZcomb ,size        
           ,"faceRZP"   ,_AD_);
  zero(mesh->elm.faceResZcomb,size,"short"  );
/*...................................................................*/

/*... zComb0 locais*/
  size = lNel*ndfComb;
  HccaAlloc(DOUBLE,m,mesh->elm.zComb0     ,size        
           ,"zComb0P"           ,_AD_);
  zero(mesh->elm.zComb0         ,size,DOUBLEC);
/*...................................................................*/

/*... zComb  locais*/
  size = lNel*ndfComb;
  HccaAlloc(DOUBLE,m,mesh->elm.zComb      ,size        
           ,"zCombP"            ,_AD_);
  zero(mesh->elm.zComb          ,size,DOUBLEC);
/*...................................................................*/


/*... yFrac0 locais*/
  size = lNel*nSp;
  HccaAlloc(DOUBLE,m,mesh->elm.yFrac0     ,size        
           ,"zFrac0P"           ,_AD_);
  zero(mesh->elm.yFrac0         ,size,DOUBLEC);
/*...................................................................*/

/*... yFrac locais*/
  size = lNel*nSp;
  HccaAlloc(DOUBLE,m,mesh->elm.yFrac      ,size        
           ,"zFracP"            ,_AD_);
  zero(mesh->elm.yFrac          ,size,DOUBLEC);
/*...................................................................*/

/*... yFrac locais*/
  size = lNel*nSp;
  HccaAlloc(DOUBLE,m,mesh->elm.cDiffComb,size        
           ,"cDiffCombP"        ,_AD_);
  zero(mesh->elm.cDiffComb,size,DOUBLEC);
/*...................................................................*/

/*... mMolar*/
  HccaAlloc(DOUBLE, m, mesh->elm.mMolar.t00, lNel, "mMolar00P", _AD_);
  HccaAlloc(DOUBLE, m, mesh->elm.mMolar.t0 , lNel, "mMolar0P" , _AD_);
  HccaAlloc(DOUBLE, m, mesh->elm.mMolar.t  , lNel, "mMolarP"  , _AD_);
  zero(mesh->elm.mMolar.t00, lNel, DOUBLEC);
  zero(mesh->elm.mMolar.t0 , lNel, DOUBLEC);
  zero(mesh->elm.mMolar.t  , lNel, DOUBLEC);
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 15/11/2019                                   *
 * ------------------------------------------------------------------*
 * comunicate2 :                                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * iCod - 1 inicializada do proprio processo master                  *
 * iCod - 2 o master envia para os escravos                          *
 * iCod - 3 os escravos recebem                                      *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/ 
static void comunicateCD(short *m0faceR,short *faceR
                ,DOUBLE *m0u0          ,DOUBLE *u0 
                ,DOUBLE *m0u           ,DOUBLE *u  
                ,LevelTime m0density   ,LevelTime density  
                ,INT const lNel        ,INT *elLG
                ,short const maxViz    ,short const ndf
                ,short const nPart     ,short const iCod)
{
#ifdef _MPI_
  INT i,size;
  DOUBLE *nD=NULL;
  short  *nS=NULL;
  short sSize=sizeof(short); 
//  short iSize=sizeof(INT); 
  short dSize=sizeof(DOUBLE); 
  switch(iCod){
    case 1: 
    case 2:
/*... faceR locais*/
      size = lNel*(maxViz+1)*ndf;
      nS   = (short *) malloc(size*sSize); 
      ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

      sGetLocalV(m0faceR           ,nS
                ,elLG
                ,lNel              ,(maxViz+1)*ndf); 
      if(iCod == 1)
        for(i=0;i<size;i++) 
          faceR[i] = nS[i];
      else if(iCod == 2)
        MPI_Send(nS    ,       size, MPI_SHORT,nPart,5    
                ,mpiVar.comm); 
      
      free(nS);  
/*...................................................................*/

/*... eU0 locais*/
      size = lNel*ndf;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

      dGetLocalV(m0u0              ,nD
                ,elLG
                ,lNel              ,ndf); 
      if(iCod == 1)  
        for(i=0;i<size;i++)
          u0[i] = nD[i];
      else if(iCod == 2)
        MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,7    
                ,mpiVar.comm);
      free(nD);
/*...................................................................*/

/*... eU locais*/
      size = lNel*ndf;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

      dGetLocalV(m0u    ,nD
                ,elLG
                ,lNel   ,ndf); 
      if(iCod == 1)   
        for(i=0;i<size;i++)
          u[i] = nD[i];
      else if(iCod == 2)
        MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,8    
                ,mpiVar.comm);;
      free(nD);
/*...................................................................*/

/*... density locais*/
      size = lNel;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);
/*... t00*/
      dGetLocalV(m0density.t00,nD
                ,elLG
                ,lNel       ,1); 
       
      if(iCod == 1)    
        for(i=0;i<size;i++)
          density.t00[i] = nD[i];
      else if(iCod == 2)
        MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,9    
               ,mpiVar.comm);
/*... t0*/
      dGetLocalV(m0density.t0,nD
                ,elLG
                ,lNel       ,1); 
        
      if(iCod == 1)    
        for(i=0;i<size;i++)
          density.t0[i] = nD[i];
      else if(iCod == 2)
        MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,9    
               ,mpiVar.comm);

/*... t*/
      dGetLocalV(m0density.t,nD
                ,elLG
                ,lNel       ,1); 
        
      if(iCod == 1)    
        for(i=0;i<size;i++)
          density.t[i] = nD[i];
      else if(iCod == 2)
        MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,9    
               ,mpiVar.comm);


      free(nD);
/*...................................................................*/
   break;

/*...*/
   case 3:
   
/*... faceR locais*/
     size = lNel*(maxViz+1)*ndf;
     MPI_Recv(faceR,size, MPI_SHORT,0,5    
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
     size = lNel;
     MPI_Recv(density.t00,size, MPI_DOUBLE,0,9    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
     MPI_Recv(density.t0 ,size, MPI_DOUBLE,0,9    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
     MPI_Recv(density.t  ,size, MPI_DOUBLE,0,9    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

   break; 
/*...................................................................*/

  }
/*...................................................................*/
#endif
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 15/08/2019                                   *
 * Data de modificaco : 23/08/2019                                   *
 * ------------------------------------------------------------------*
 * comunicateComb :                                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * iCod - 1 inicializada do proprio processo master                  *
 * iCod - 2 o master envia para os escravos                          *
 * iCod - 3 os escravos recebem                                      *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/ 
static void comunicateComb(short *m0faceR     ,short *faceR
                ,DOUBLE *m0z0       ,DOUBLE *z0 
                ,DOUBLE *m0z        ,DOUBLE *z  
                ,DOUBLE *m0y0       ,DOUBLE *y0 
                ,DOUBLE *m0y        ,DOUBLE *y
                ,DOUBLE *m0Diff     ,DOUBLE *mDiff
                ,LevelTime m0mMolar ,LevelTime mmMolar
                ,INT const lNel     ,INT *elLG
                ,short const maxViz ,short const ndfComb
                ,short const nSp  
                ,short const nPart  ,short const iCod)
{
#ifdef _MPI_
  INT i,size;
  DOUBLE *nD=NULL;
  short  *nS=NULL;
  short sSize=sizeof(short); 
//  short iSize=sizeof(INT); 
  short dSize=sizeof(DOUBLE); 

  switch(iCod)
  {
    case 1: 
    case 2:

/*... faceRzComb locais*/
      size = lNel*(maxViz+1);
      nS   = (short *) malloc(size*sSize); 
      ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

      sGetLocalV(m0faceR           ,nS
                ,elLG
                ,lNel              ,maxViz+1); 
      if(iCod == 1)
        for(i=0;i<size;i++) 
          faceR[i] = nS[i];
      else if(iCod == 2)
        MPI_Send(nS    ,       size, MPI_SHORT,nPart,22   
                 ,mpiVar.comm);

      free(nS); 
/*...................................................................*/

/*... eZcomb0 locais*/
      size = lNel*ndfComb;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

      dGetLocalV(m0z0              ,nD
                ,elLG
                ,lNel              ,ndfComb); 
      if(iCod == 1)  
        for(i=0;i<size;i++)
          z0[i] = nD[i];
      else if(iCod == 2)
        MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,24   
                ,mpiVar.comm);

      free(nD);
/*...................................................................*/

/*... eZcomb locais*/
      size = lNel*ndfComb;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

      dGetLocalV(m0z    ,nD
                ,elLG
                ,lNel   ,ndfComb); 
        
      if(iCod == 1)  
        for(i=0;i<size;i++)
          z[i] = nD[i];
      else if(iCod == 2)
        MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,25   
                ,mpiVar.comm);

      free(nD);
/*...................................................................*/

/*... eY0 locais*/
      size = lNel*nSp;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

      dGetLocalV(m0y0              ,nD
                ,elLG
                ,lNel              ,nSp); 
        
      if(iCod == 1)  
        for(i=0;i<size;i++)
          y0[i] = nD[i];
      else if(iCod == 2)
        MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,26   
                ,mpiVar.comm);

      free(nD);
/*...................................................................*/

/*... eY locais*/
      size = lNel*nSp;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

      dGetLocalV(m0y    ,nD
                ,elLG
                ,lNel   ,nSp); 
        
      if(iCod == 1)  
        for(i=0;i<size;i++)
          y[i] = nD[i];
      else if(iCod == 2)
        MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,27    
                ,mpiVar.comm);

      free(nD);
/*...................................................................*/

/*... diff locais*/
      size = lNel*nSp;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);
  
      dGetLocalV(m0Diff ,nD
                ,elLG
                ,lNel   ,nSp); 
        
      if(iCod == 1)  
        for(i=0;i<size;i++)
          mDiff[i] = nD[i];
      else if(iCod == 2)
        MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,28    
                ,mpiVar.comm);

      free(nD);
/*...................................................................*/

/*... molar mass locais*/
        size = lNel;
        nD   = (DOUBLE *) malloc(size*dSize); 
        ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);
/*... t00*/
        dGetLocalV(m0mMolar.t00  ,nD
                  ,elLG
                  ,lNel         ,1); 
          
        if(iCod == 1)
          for(i=0;i<size;i++) 
            m0mMolar.t00[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,29   
                  ,mpiVar.comm);
/*... t0*/
        dGetLocalV(m0mMolar.t0  ,nD
                  ,elLG
                  ,lNel         ,1); 
          
        if(iCod == 1)
          for(i=0;i<size;i++) 
            m0mMolar.t0[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,29   
                  ,mpiVar.comm);

/*... t*/
        dGetLocalV(m0mMolar.t  ,nD
                  ,elLG
                  ,lNel         ,1); 
          
        if(iCod == 1)
          for(i=0;i<size;i++) 
            m0mMolar.t[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,29   
                  ,mpiVar.comm);
        
        free(nD);
/*...................................................................*/


   break;
/*...................................................................*/

/*...*/
   case 3:
 
/*... faceRzComb locais*/
     size = lNel*(maxViz+1);
     MPI_Recv(faceR,size, MPI_SHORT,0,22   
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/
     
/*... eZ0 locais*/
     size = lNel*ndfComb;
     MPI_Recv(z0  ,size, MPI_DOUBLE,0,24   
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/
     
/*... eZ locais*/
     size = lNel*ndfComb;
     MPI_Recv(z,size, MPI_DOUBLE,0,25   
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... eY0 locais*/
     size = lNel*nSp;
     MPI_Recv(y0  ,size, MPI_DOUBLE,0,26   
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/
     
/*... eY locais*/
     size = lNel*nSp;
     MPI_Recv(y,size, MPI_DOUBLE,0,27    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... diff locais*/
     size = lNel*nSp;
     MPI_Recv(mDiff,size, MPI_DOUBLE,0,28    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/
     
/*... massa molar*/
     size = lNel;
     MPI_Recv(mmMolar.t00,size, MPI_DOUBLE,0,29    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
     MPI_Recv(mmMolar.t0,size, MPI_DOUBLE,0,29    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
     MPI_Recv(mmMolar.t,size, MPI_DOUBLE,0,29    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

   break; 
/*...................................................................*/

  }
/*...................................................................*/
#endif
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 20/11/2019                                   *
 * ------------------------------------------------------------------*
 * comunicateFluid :                                                 * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * iCod - 1 inicializada do proprio processo master                  *
 * iCod - 2 o master envia para os escravos                          *
 * iCod - 3 os escravos recebem                                      *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/ 
static void comunicateFluid(short *m0faceRvel     ,short *faceRvel
                ,short *m0faceRpres    ,short *faceRpres
                ,short *m0faceRenergy  ,short *faceRenergy
                ,short *m0faceRrho     ,short *faceRrho     
                ,DOUBLE *m0pres0       ,DOUBLE *pres0
                ,DOUBLE *m0pres        ,DOUBLE *pres
                ,DOUBLE *m0energy0     ,DOUBLE *energy0
                ,DOUBLE *m0energy      ,DOUBLE *energy
                ,DOUBLE *m0temp0       ,DOUBLE *temp0
                ,DOUBLE *m0temp        ,DOUBLE *temp
                ,LevelTime m0density   ,LevelTime density
                ,LevelTime m0sHeat     ,LevelTime sHeat 
                ,DOUBLE *m0dVisc       ,DOUBLE *dVisc
                ,DOUBLE *m0tCond       ,DOUBLE *tCond 
                ,INT const lNel        ,INT *elLG
                ,short const maxViz    ,short const ndfF
                ,short const ndfFt
                ,short const nPart     ,short const iCod)
{
#ifdef _MPI_
  INT i,size;
  DOUBLE *nD=NULL;
  short  *nS=NULL;
  short sSize=sizeof(short); 
//  short iSize=sizeof(INT); 
  short dSize=sizeof(DOUBLE); 

/*...*/
  switch(iCod)
  {

/*...*/
    case 1: 
    case 2:
/*... faceRvel locais*/
      size = lNel*(maxViz+1);
      nS   = (short *) malloc(size*sSize); 
      ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

      sGetLocalV(m0faceRvel        ,nS
                ,elLG
                ,lNel              ,maxViz+1); 

      if(iCod == 1)
        for(i=0;i<size;i++) 
          faceRvel[i] = nS[i];
      else if(iCod == 2)
        MPI_Send(nS    ,       size, MPI_SHORT,nPart,6    
                ,mpiVar.comm);
      
      free(nS);  
/*...................................................................*/

/*... faceRpres locais*/
      size = lNel*(maxViz+1);
      nS   = (short *) malloc(size*sSize); 
      ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

      sGetLocalV(m0faceRpres       ,nS
                ,elLG
                ,lNel              ,maxViz+1); 
      
      if(iCod == 1)
        for(i=0;i<size;i++) 
          faceRpres[i] = nS[i];
      else if(iCod == 2)
        MPI_Send(nS    ,       size, MPI_SHORT,nPart,8    
                ,mpiVar.comm);
      
      free(nS);  
/*...................................................................*/

/*... pres0 locais*/
      size = lNel;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

      dGetLocalV(m0pres0           ,nD
                ,elLG
                ,lNel              ,1  ); 

      if(iCod == 1)
        for(i=0;i<size;i++) 
          pres0[i] = nD[i];
      else if(iCod == 2)
         MPI_Send(nD    ,       size, MPI_SHORT,nPart,10  
                 ,mpiVar.comm);

      free(nD);
/*...................................................................*/

/*... pres locais*/
      size = lNel;
      nD   = (DOUBLE *) malloc(size*dSize); 
      ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

      dGetLocalV(m0pres0    ,nD
                ,elLG
                ,lNel       ,1); 
        
      if(iCod == 1)
        for(i=0;i<size;i++) 
          pres[i] = nD[i];
      else if(iCod == 2)
         MPI_Send(nD    ,       size, MPI_SHORT,nPart,11   
                 ,mpiVar.comm);

      free(nD);
/*...................................................................*/

/*...*/
      if(ndfFt > 0)
      {
/*... faceRenergy locais*/
        size = lNel*(maxViz+1);
        nS   = (short *) malloc(size*sSize); 
        ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

        sGetLocalV(m0faceRenergy     ,nS
                  ,elLG
                  ,lNel              ,maxViz+1); 
      
        if(iCod == 1)
          for(i=0;i<size;i++) 
            faceRenergy[i] = nS[i];
        else if(iCod == 2)
          MPI_Send(nS    ,       size, MPI_SHORT,nPart,12   
                ,mpiVar.comm);
      
        free(nS);
/*...................................................................*/

/*... faceRrho locais*/
        size = lNel*(maxViz+1);
        nS   = (short *) malloc(size*sSize); 
        ERRO_MALLOC(nS,"nS"   ,__LINE__,__FILE__,__func__);

        sGetLocalV(m0faceRrho        ,nS
                  ,elLG
                  ,lNel              ,maxViz+1); 
      
        if(iCod == 1)
          for(i=0;i<size;i++) 
            faceRrho[i] = nS[i];
        else if(iCod == 2)
          MPI_Send(nS    ,       size, MPI_SHORT,nPart,13   
                ,mpiVar.comm);
      
        free(nS);
/*...................................................................*/

/*... energy0 locais*/
        size = lNel;
        nD   = (DOUBLE *) malloc(size*dSize); 
        ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

        dGetLocalV(m0energy0         ,nD
                  ,elLG
                  ,lNel              ,1  ); 
        
        if(iCod == 1)
          for(i=0;i<size;i++) 
            energy0[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,14   
                 ,mpiVar.comm);

        free(nD);
/*...................................................................*/

/*... energy locais*/
        size = lNel;
        nD   = (DOUBLE *) malloc(size*dSize); 
        ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

        dGetLocalV(m0energy   ,nD
                  ,elLG
                  ,lNel       ,1); 
        
        if(iCod == 1)
          for(i=0;i<size;i++) 
            energy[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,15   
                 ,mpiVar.comm);

        free(nD);
/*...................................................................*/

/*... temp0 locais*/
        size = lNel;
        nD   = (DOUBLE *) malloc(size*dSize); 
        ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

        dGetLocalV(m0temp0           ,nD
                  ,elLG
                  ,lNel              ,1  ); 
        
        if(iCod == 1)
          for(i=0;i<size;i++) 
            temp0[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,16   
                 ,mpiVar.comm);

        free(nD);
/*...................................................................*/

/*... temp locais*/
        size = lNel;
        nD   = (DOUBLE *) malloc(size*dSize); 
        ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);

        dGetLocalV(m0temp            ,nD
                  ,elLG
                  ,lNel              ,1  ); 
        
        if(iCod == 1)
          for(i=0;i<size;i++) 
            temp[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,17   
                 ,mpiVar.comm);

        free(nD);
/*...................................................................*/

/*... density locais*/
        size = lNel;
        nD   = (DOUBLE *) malloc(size*dSize); 
        ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);
/*... t00*/
        dGetLocalV(m0density.t00  ,nD
                  ,elLG
                  ,lNel         ,1); 
          
        if(iCod == 1)
          for(i=0;i<size;i++) 
            density.t00[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,18   
                  ,mpiVar.comm);
/*... t0*/
        dGetLocalV(m0density.t0  ,nD
                  ,elLG
                  ,lNel         ,1); 
          
        if(iCod == 1)
          for(i=0;i<size;i++) 
            density.t0[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,18   
                  ,mpiVar.comm);

/*... t*/
        dGetLocalV(m0density.t  ,nD
                  ,elLG
                  ,lNel         ,1); 
           
        if(iCod == 1)
          for(i=0;i<size;i++) 
            density.t[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,18   
                  ,mpiVar.comm);
        
        free(nD);
/*...................................................................*/

/*... calor especifico*/
        size = lNel;
        nD   = (DOUBLE *) malloc(size*dSize); 
        ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);
/*... t*/        
        dGetLocalV(m0sHeat.t  ,nD
                  ,elLG
                  ,lNel       ,1); 
          
        if(iCod == 1)
          for(i=0;i<size;i++) 
            sHeat.t[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,19   
                  ,mpiVar.comm);
/*... t0*/

        dGetLocalV(m0sHeat.t0 ,nD
                  ,elLG
                  ,lNel       ,1); 
          
        if(iCod == 1)
          for(i=0;i<size;i++) 
            sHeat.t0[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,19   
                  ,mpiVar.comm);
/*... t00*/
        dGetLocalV(m0sHeat.t00,nD
                  ,elLG
                  ,lNel       ,1); 
          
        if(iCod == 1)
          for(i=0;i<size;i++) 
            sHeat.t00[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,19   
                  ,mpiVar.comm);        

        free(nD);
/*...................................................................*/

/*... viscosidade dinamica locais*/
        size = lNel;
        nD   = (DOUBLE *) malloc(size*dSize); 
        ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);
        
        dGetLocalV(m0dVisc    ,nD
                  ,elLG
                  ,lNel       ,1); 
        
        if(iCod == 1)
          for(i=0;i<size;i++) 
            dVisc[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,20   
                  ,mpiVar.comm);  
        
        free(nD);
/*...................................................................*/

/*... condutividade termica locais*/
        size = lNel;
        nD   = (DOUBLE *) malloc(size*dSize); 
        ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);
        
        dGetLocalV(m0tCond            ,nD
                  ,elLG
                  ,lNel               ,1); 
        
        if(iCod == 1)
          for(i=0;i<size;i++) 
            tCond[i] = nD[i];
        else if(iCod == 2)
          MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,21   
                  ,mpiVar.comm); 
        
        free(nD);
/*...................................................................*/
      }
/*...................................................................*/

   break;

/*...*/
   case 3:
   
/*... faceRvel locais*/
     size = lNel*(maxViz+1);
     MPI_Recv(faceRvel,size, MPI_SHORT,0,6    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/
     
/*... faceRpres locais*/
     size = lNel*(maxViz+1);
     MPI_Recv(faceRpres,size, MPI_SHORT,0,8    
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/
     
/*... pres0 locais*/
     size = lNel;
     MPI_Recv(pres0      ,size, MPI_DOUBLE,0,10  
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... pres locais*/
     size = lNel;
     MPI_Recv(pres       ,size, MPI_DOUBLE,0,11   
             ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/
     
/*...*/
    if(ndfFt > 0)
    {
/*... faceRenergy locais*/
        size = lNel*(maxViz+1);
        MPI_Recv(faceRenergy,size, MPI_SHORT,0,12   
                ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... faceRpres locais*/
       size = lNel*(maxViz+1);
       MPI_Recv(faceRrho,size, MPI_SHORT,0,13   
               ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/


/*... energy0 locais*/
       size = lNel;
       MPI_Recv(energy0    ,size, MPI_DOUBLE,0,14   
               ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... energy locais*/
       size = lNel;
       MPI_Recv(energy     ,size, MPI_DOUBLE,0,15   
               ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... temp0 locais*/
       size = lNel;
       MPI_Recv(temp0      ,size, MPI_DOUBLE,0,16   
               ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... temp locais*/
       size = lNel;
       MPI_Recv(temp       ,size, MPI_DOUBLE,0,17   
               ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... massa especifica locais*/
       size = lNel;
       MPI_Recv(density.t,size, MPI_DOUBLE,0,18    
               ,mpiVar.comm,MPI_STATUS_IGNORE);
       MPI_Recv(density.t0,size, MPI_DOUBLE,0,18    
               ,mpiVar.comm,MPI_STATUS_IGNORE);
       MPI_Recv(density.t00,size, MPI_DOUBLE,0,18    
               ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... calor especifico locais*/
       size = lNel; 
       MPI_Recv(sHeat.t,size, MPI_DOUBLE,0,19    
               ,mpiVar.comm,MPI_STATUS_IGNORE);
       MPI_Recv(sHeat.t0,size, MPI_DOUBLE,0,19    
               ,mpiVar.comm,MPI_STATUS_IGNORE);
       MPI_Recv(sHeat.t00,size, MPI_DOUBLE,0,19    
               ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... viscosidade dinamica locais*/
       size = lNel;
       MPI_Recv(dVisc,size, MPI_DOUBLE,0,20    
               ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/

/*... condutividade termica locais*/
       size = lNel;
       MPI_Recv(tCond,size, MPI_DOUBLE,0,21    
               ,mpiVar.comm,MPI_STATUS_IGNORE);
/*...................................................................*/
     }
/*...................................................................*/
   break; 
/*...................................................................*/


  }
/*...................................................................*/
#endif
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 12/08/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * prMaster : inicializa da particao no processo master              * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 * não ha comunicao pois e o maste que gerou a particao              *
 *-------------------------------------------------------------------* 
 *********************************************************************/
static void prMaster(Memoria *m    , Mesh *mesh0
    , Mesh *mesh                   , PartMesh *pMesh
    , INT *fMap                    , INT *fMapNo
    , INT *iaSends                 , INT *iaRcvs
    , INT *iaComNo                 
    , INT *elLG                    , INT *elGL
    , INT *noLG                    , INT *noGL
    , short *maxVizPart            , short *maxVizPartNo
    , INT const lNel               , INT const lNnode
    , INT const nRcvs              , INT const nSends
    , INT const nComNo             , INT const maxNo  
    , short *ndfD                  , short *ndfT
    , short const ndfF             , short const ndfFt
    , short const ndfComb          , short const nSp  
    , short const ndm              , short const maxViz
    , unsigned short const nVizPart, unsigned short const nVizPartNo)
{

  short sSize=sizeof(short), 
        iSize=sizeof(INT), 
        dSize=sizeof(DOUBLE);
  short *nS=NULL;
  INT i,kk,size;
  INT *lEl=NULL;
  DOUBLE *nD=NULL;

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
  
  for(i=0;i<kk;i++)
  {
    pMesh->iEl.iaSends[i] = iaSends[i];
    pMesh->iEl.iaRcvs[i]  = iaRcvs[i];
  }

/*... iaComNo*/
  kk = nVizPartNo+1;
  
  HccaAlloc(INT   ,m   ,pMesh->iNo.iaComNo , kk        
           ,"iaComNo",_AD_);
  zero(pMesh->iNo.iaComNo   , kk,INTC);
  
  for(i=0;i<kk;i++)
    pMesh->iNo.iaComNo[i] = iaComNo[i];
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
  if(ndfT[0] > 0 || ndfF > 0 || ndfFt > 0)
  {   
/*... eVel0*/
    size = lNel*ndm;
    nD   = (DOUBLE *) malloc(size*dSize); 
    ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);
    
    HccaAlloc(DOUBLE,m,mesh->elm.vel0     ,size        
             ,"eVel0P"         ,_AD_);
    zero(mesh->elm.vel0      ,size,DOUBLEC);
    
    dGetLocalV(mesh0->elm.vel0   ,nD
              ,elLG
              ,lNel              ,ndm); 
    
    for(i=0;i<size;i++)
      mesh->elm.vel0[i] = nD[i];
    
    free(nD);  
/*...................................................................*/
  
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
    
    for(i=0;i<size;i++)
      mesh->elm.vel[i] = nD[i];
    
    free(nD);  
/*...................................................................*/
  }
/*...................................................................*/

/*... */
  if(ndfD[0] > 0)
  {
/*... faceRd1 locais*/
     size = lNel*(maxViz+1)*ndfD[0];
     HccaAlloc(short,m,mesh->elm.faceRd1  ,size        
              ,"faceRd1P"      ,_AD_);
     zero(mesh->elm.faceRd1   ,size,"short"  );
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
     size = lNel;
     HccaAlloc(DOUBLE,m,mesh->elm.densityUd1.t00  ,size        
              ,"densityP00"            ,_AD_); 
     HccaAlloc(DOUBLE,m,mesh->elm.densityUd1.t0  ,size        
              ,"densityP0"            ,_AD_); 
     HccaAlloc(DOUBLE,m,mesh->elm.densityUd1.t  ,size        
              ,"densityP"            ,_AD_); 
     zero(mesh->elm.densityUd1.t00,size,DOUBLEC);
     zero(mesh->elm.densityUd1.t0,size,DOUBLEC)
     zero(mesh->elm.densityUd1.t,size,DOUBLEC)
/*...................................................................*/

/*... */ 
     comunicateCD(mesh0->elm.faceRd1   ,mesh->elm.faceRd1
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
  if(ndfT[0] > 0)
  { 
/*... faceRt1 locais*/
    size = lNel*(maxViz+1)*ndfT[0];
    HccaAlloc(short,m,mesh->elm.faceRt1  ,size        
             ,"faceRt1P"      ,_AD_);
    zero(mesh->elm.faceRt1   ,size,"short"  );
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
    size = lNel;
    HccaAlloc(DOUBLE,m,mesh->elm.densityUt1.t00  ,size          
             ,"densityT1P00"          ,_AD_);
    HccaAlloc(DOUBLE,m,mesh->elm.densityUt1.t0  ,size          
             ,"densityT1P0"          ,_AD_);
    HccaAlloc(DOUBLE,m,mesh->elm.densityUt1.t  ,size          
             ,"densityT1P"          ,_AD_);
    zero(mesh->elm.densityUt1.t00,size,DOUBLEC);
    zero(mesh->elm.densityUt1.t0,size,DOUBLEC);
    zero(mesh->elm.densityUt1.t,size,DOUBLEC);
/*...................................................................*/

/*... */
    comunicateCD(mesh0->elm.faceRt1   ,mesh->elm.faceRt1
               ,mesh0->elm.u0T1      ,mesh->elm.u0T1
               ,mesh0->elm.uT1       ,mesh->elm.uT1
               ,mesh0->elm.densityUt1,mesh->elm.densityUt1
               ,lNel                 ,elLG
               ,maxViz               ,ndfT[0]
               ,0                    ,1);
  }
/*...................................................................*/

/*...*/
  if(ndfF > 0 || ndfFt > 0)
  {
    allocFuild(m     , mesh
              ,maxViz, lNel
              ,ndfFt);

/*...*/
    comunicateFluid(mesh0->elm.faceRvel  ,mesh->elm.faceRvel
              ,mesh0->elm.faceRpres      ,mesh->elm.faceRpres
              ,mesh0->elm.faceRenergy    ,mesh->elm.faceRenergy
              ,mesh0->elm.faceRrho       ,mesh->elm.faceRrho     
              ,mesh0->elm.pressure0      ,mesh->elm.pressure0
              ,mesh0->elm.pressure       ,mesh->elm.pressure
              ,mesh0->elm.energy0        ,mesh->elm.energy0
              ,mesh0->elm.energy         ,mesh->elm.energy
              ,mesh0->elm.temp0          ,mesh->elm.temp0  
              ,mesh0->elm.temp           ,mesh->elm.temp       
              ,mesh0->elm.densityFluid   ,mesh->elm.densityFluid
              ,mesh0->elm.specificHeat   ,mesh->elm.specificHeat 
              ,mesh0->elm.dViscosity     ,mesh->elm.dViscosity 
              ,mesh0->elm.tConductivity  ,mesh->elm.tConductivity 
              ,lNel                      ,elLG
              ,maxViz                    ,ndfF
              ,ndfFt
              ,0                         ,1);
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if(ndfComb>0)
  {
/*...*/
    allocComb(m      , mesh
            , maxViz , lNel
            , ndfComb, nSp);  
/*...................................................................*/

/*...*/
    comunicateComb(mesh0->elm.faceResZcomb ,mesh->elm.faceResZcomb 
                  ,mesh0->elm.zComb0       ,mesh->elm.zComb0
                  ,mesh0->elm.zComb        ,mesh->elm.zComb
                  ,mesh0->elm.yFrac0       ,mesh->elm.yFrac0
                  ,mesh0->elm.yFrac        ,mesh->elm.yFrac 
                  ,mesh0->elm.cDiffComb    ,mesh->elm.cDiffComb
                  ,mesh0->elm.mMolar       ,mesh->elm.mMolar
                  ,lNel                    ,elLG
                  ,maxViz                  ,ndfComb
                  ,nSp
                  ,0                       ,1);
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
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 12/08/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * sendPart : envia as particoes para os outros processos            * 
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
static void sendPart(Memoria *m    , Mesh *mesh0
    , Mesh *mesh                   , PartMesh *pMesh
    , INT *fMap                    , INT *fMapNo
    , INT *iaSends                 , INT *iaRcvs
    , INT *iaComNo                 
    , INT *elLG                    , INT *elGL
    , INT *noLG                    , INT *noGL
    , INT const numelOv            , INT const nNodeOv
    , INT const nno1 
    , short *maxVizPart            , short *maxVizPartNo
    , INT const lNel               , INT const lNnode
    , INT const nRcvs              , INT const nSends
    , INT const nComNo             , INT const maxNo  
    , short *ndfD                  , short *ndfT
    , short const ndfF             , short const ndfFt
    , short const ndfComb          , short const nSp
    , short const ndm              , short const maxViz
    , unsigned short const nVizPart, unsigned short const nVizPartNo
    , unsigned short const nPart )
{
#if _MPI_
  short sSize=sizeof(short), 
        iSize=sizeof(INT), 
        dSize=sizeof(DOUBLE);
  short *nS=NULL;
  INT i,kk,size;
  INT *lEl=NULL;
  DOUBLE *nD=NULL;
  
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
  if(ndfT[0] > 0 || ndfF > 0 || ndfFt > 0)
  {     
/*... eVel0*/
    size = lNel*ndm;
    nD   = (DOUBLE *) malloc(size*dSize); 
    ERRO_MALLOC(nD,"nD"   ,__LINE__,__FILE__,__func__);
    
    dGetLocalV(mesh0->elm.vel0   ,nD
              ,elLG
              ,lNel              ,ndm); 
    
    MPI_Send(nD    ,       size, MPI_DOUBLE,nPart,5    
            ,mpiVar.comm);
    
    free(nD);  
/*...................................................................*/

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
  if( ndfD[0] > 0)
    comunicateCD(mesh0->elm.faceRd1   ,mesh->elm.faceRd1
               ,mesh0->elm.u0D1      ,mesh->elm.u0D1
               ,mesh0->elm.uD1       ,mesh->elm.uD1
               ,mesh0->elm.densityUd1,mesh->elm.densityUd1
               ,lNel                 ,elLG
               ,maxViz               ,ndfD[0]
               ,nPart                ,2);
/*...................................................................*/

/*... */
  if( ndfT[0] > 0)
    comunicateCD(mesh0->elm.faceRt1   ,mesh->elm.faceRt1
               ,mesh0->elm.u0T1      ,mesh->elm.u0T1
               ,mesh0->elm.uT1       ,mesh->elm.uT1
               ,mesh0->elm.densityUt1,mesh->elm.densityUt1
               ,lNel                 ,elLG
               ,maxViz               ,ndfT[0]
               ,nPart                ,2);
/*...................................................................*/

/*...*/
  if(ndfF > 0 || ndfFt > 0)
    comunicateFluid(mesh0->elm.faceRvel  ,mesh->elm.faceRvel
              ,mesh0->elm.faceRpres      ,mesh->elm.faceRpres
              ,mesh0->elm.faceRenergy    ,mesh->elm.faceRenergy
              ,mesh0->elm.faceRrho       ,mesh->elm.faceRrho       
              ,mesh0->elm.pressure0      ,mesh->elm.pressure0
              ,mesh0->elm.pressure       ,mesh->elm.pressure
              ,mesh0->elm.energy0        ,mesh->elm.energy0
              ,mesh0->elm.energy         ,mesh->elm.energy
              ,mesh0->elm.temp0          ,mesh->elm.temp0  
              ,mesh0->elm.temp           ,mesh->elm.temp       
              ,mesh0->elm.densityFluid   ,mesh->elm.densityFluid
              ,mesh0->elm.specificHeat   ,mesh->elm.specificHeat 
              ,mesh0->elm.dViscosity     ,mesh->elm.dViscosity 
              ,mesh0->elm.tConductivity  ,mesh->elm.tConductivity 
              ,lNel                      ,elLG
              ,maxViz                    ,ndfF
              ,ndfFt
              ,nPart                     ,2);
/*...................................................................*/

/*...*/
  if(ndfComb>0)
    comunicateComb(mesh0->elm.faceResZcomb ,mesh->elm.faceResZcomb
                  ,mesh0->elm.zComb0       ,mesh->elm.zComb0
                  ,mesh0->elm.zComb        ,mesh->elm.zComb
                  ,mesh0->elm.yFrac0       ,mesh->elm.yFrac0
                  ,mesh0->elm.yFrac        ,mesh->elm.yFrac 
                  ,mesh0->elm.cDiffComb    ,mesh->elm.cDiffComb 
                  ,mesh0->elm.mMolar       ,mesh->elm.mMolar       
                  ,lNel                    ,elLG
                  ,maxViz                  ,ndfComb
                  ,nSp
                  ,nPart                   ,2);
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
#endif
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 14/08/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * recvPart : recebe as particoes para os outros processos           * 
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
static void recvPart(Memoria *m    , Mesh *mesh0
    , Mesh *mesh                   , PartMesh *pMesh
    , INT *fMap                    , INT *fMapNo
    , INT *iaSends                 , INT *iaRcvs
    , INT *iaComNo                 
    , INT *elLG                    , INT *elGL
    , INT *noLG                    , INT *noGL
    , INT *numelOv                 , INT *nNodeOv
    , INT *nno1 
    , short *maxVizPart            , short *maxVizPartNo
    , INT *lNel                    , INT *lNnode
    , INT *nRcvs                   , INT *nSends
    , INT *nComNo                  , INT const maxNo  
    , short *ndfD                  , short *ndfT
    , short const ndfF             , short const ndfFt
    , short const ndfComb          , short const nSp
    , short const ndm              , short const maxViz
    , short *nVizPart              , short *nVizPartNo
    , unsigned short const nPart )
{

#if _MPI_
  short sSize=sizeof(short), 
        iSize=sizeof(INT), 
        dSize=sizeof(DOUBLE);
  short *nS=NULL;
  INT i,kk,size;
  INT *lEl=NULL;
  DOUBLE *nD=NULL;

  MPI_Recv(lNnode    , 1, MPI_INT,0,nPart,
           mpiVar.comm,MPI_STATUS_IGNORE); 
  MPI_Recv(lNel      , 1, MPI_INT,0,nPart,
           mpiVar.comm,MPI_STATUS_IGNORE); 
  MPI_Recv(numelOv   , 1, MPI_INT,0,nPart,
           mpiVar.comm,MPI_STATUS_IGNORE); 
  MPI_Recv(nNodeOv   , 1, MPI_INT,0,nPart,
           mpiVar.comm,MPI_STATUS_IGNORE); 
  MPI_Recv(nno1      , 1, MPI_INT,0,nPart,
           mpiVar.comm,MPI_STATUS_IGNORE); 

  mesh->nnode    = *lNnode;
  mesh->numel    = *lNel;
  mesh->numelNov = *lNel-*numelOv;
  mesh->nnodeOv  = *nNodeOv ;
  mesh->nnodeNov = *lNnode-*nNodeOv;
  pMesh->nno1    = *nno1;
/*... particionamento*/     

/*... partViz(El)*/
  MPI_Recv(nVizPart  , 1, MPI_SHORT,0,nPart,
           mpiVar.comm,MPI_STATUS_IGNORE); 
  
  pMesh->iEl.nVizPart = *nVizPart;
  
  size = *nVizPart;
  HccaAlloc(short,m,pMesh->iEl.vizPart , size   
           ,"nVizPart"             ,_AD_);
  zero(pMesh->iEl.vizPart               ,size,"short");
  
  MPI_Recv(pMesh->iEl.vizPart,size, MPI_SHORT,0,nPart,
           mpiVar.comm,MPI_STATUS_IGNORE);
/*... partViz(No)*/
  MPI_Recv(nVizPartNo, 1, MPI_SHORT,0,nPart,
           mpiVar.comm,MPI_STATUS_IGNORE); 
  
  pMesh->iNo.nVizPart = *nVizPartNo;
  size = *nVizPartNo;
  HccaAlloc(short,m,pMesh->iNo.vizPart   , size   
           ,"nVizPartNo"           ,_AD_);
  zero(pMesh->iNo.vizPart               ,size,"short");
  
  MPI_Recv(pMesh->iNo.vizPart,size, MPI_SHORT,0,nPart,
           mpiVar.comm,MPI_STATUS_IGNORE); 
/*..................................................................*/

/*...*/     
  MPI_Recv(nRcvs        ,       1, MPI_INT,0,nPart 
          ,mpiVar.comm,MPI_STATUS_IGNORE); 
  MPI_Recv(nSends       ,       1, MPI_INT,0,nPart 
          ,mpiVar.comm,MPI_STATUS_IGNORE); 
  MPI_Recv(nComNo       ,       1, MPI_INT,0,nPart 
          ,mpiVar.comm,MPI_STATUS_IGNORE); 
  
  pMesh->iEl.nRcvs  = *nRcvs;
  pMesh->iEl.nSends = *nSends;
  pMesh->iNo.nCom   = *nComNo; 
/*..................................................................*/     
      
/*... fMapEl(buffer de comunicacao)*/ 
  size = *nRcvs + *nSends;
  HccaAlloc(INT  ,m,pMesh->iEl.fMap    , size   
           ,"fMap"                 ,_AD_);
  zero(pMesh->iEl.fMap,size,INTC);
  
  MPI_Recv(pMesh->iEl.fMap ,    size, MPI_INT,0,nPart 
          ,mpiVar.comm,MPI_STATUS_IGNORE); 
/*... fMapNo (buffer de comunicacao)*/ 
  size = *nComNo;
  HccaAlloc(INT  ,m,pMesh->iNo.fMap    , size   
           ,"fMapNo"               ,_AD_);
  zero(pMesh->iNo.fMap,size,INTC);
  
  MPI_Recv(pMesh->iNo.fMap ,    size, MPI_INT,0,nPart 
          ,mpiVar.comm,MPI_STATUS_IGNORE); 
/*...................................................................*/   
      
/*... send e rcvs*/       
  kk = *nVizPart + 1;
  
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
  kk = *nVizPartNo + 1;
  
  HccaAlloc(INT   ,m   ,pMesh->iNo.iaComNo , kk        
             ,"iaSendsNo",_AD_);
  zero(pMesh->iNo.iaComNo   , kk,INTC);
    
  MPI_Recv(pMesh->iNo.iaComNo,  kk, MPI_INT,0,nPart 
          ,mpiVar.comm,MPI_STATUS_IGNORE); 
/*...................................................................*/

/*... elLG( numeracao local-global de elementos)*/
  size = *lNel;
  HccaAlloc(INT   ,m   ,pMesh->elLG , size        
             ,"elLG"   ,_AD_);
  zero(pMesh->elLG   , size,INTC);
  MPI_Recv(pMesh->elLG,  size, MPI_INT,0,nPart 
          ,mpiVar.comm,MPI_STATUS_IGNORE); 
/*..................................................................*/

/*... noLG( numeracao local-global de nos)*/
  size = *lNnode;
  HccaAlloc(INT   ,m   ,pMesh->noLG , size       
             ,"noLG"   ,_AD_);
  zero(pMesh->noLG   , size,INTC);
  MPI_Recv(pMesh->noLG,  size, MPI_INT,0,nPart 
          ,mpiVar.comm,MPI_STATUS_IGNORE); 
/*..................................................................*/

/*... elementos locais*/
  size = (*lNel)*maxNo;
  
  HccaAlloc(INT,m,mesh->elm.node  ,size        
           ,"elnodeP"  ,_AD_);
  zero(mesh->elm.node         ,size ,INTC);
  
  MPI_Recv(mesh->elm.node,size, MPI_INT,0,1    
          ,mpiVar.comm,MPI_STATUS_IGNORE);      
/*...................................................................*/

/*... mat locais*/
  size = *lNel;
  
  HccaAlloc(short,m,mesh->elm.mat  ,size        
           ,"elMatP"     ,_AD_);
  zero(mesh->elm.mat          ,size         ,"short"  );
  
  MPI_Recv(mesh->elm.mat,size, MPI_SHORT,0,2    
          ,mpiVar.comm,MPI_STATUS_IGNORE);      
/*...................................................................*/

/*... nen locais*/
  size = *lNel;
  
  HccaAlloc(short,m,mesh->elm.nen  ,size        
           ,"elNenP"     ,_AD_);
  zero(mesh->elm.nen          ,size         ,"short"  );
  
  MPI_Recv(mesh->elm.nen,size, MPI_SHORT,0,3    
          ,mpiVar.comm,MPI_STATUS_IGNORE);      
/*...................................................................*/

/*... geomType locais*/
  size = *lNel;

  HccaAlloc(short,m,mesh->elm.geomType  ,size        
           ,"elGtP"     ,_AD_);
  zero(mesh->elm.geomType     ,size          ,"short"  );
  
  MPI_Recv(mesh->elm.geomType,size, MPI_SHORT,0,4    
          ,mpiVar.comm,MPI_STATUS_IGNORE);      
/*...................................................................*/

/*... transporte e fluido*/
  if(ndfT[0] > 0 || ndfF > 0 || ndfFt > 0)
  {     
/*... eVel0*/
    size = (*lNel)*ndm;
      
    HccaAlloc(DOUBLE,m,mesh->elm.vel0      ,size        
         ,"evel0P"   ,_AD_);
    zero(mesh->elm.vel0       ,size,DOUBLEC);
    
    MPI_Recv(mesh->elm.vel0,       size, MPI_DOUBLE,0,5    
            ,mpiVar.comm,MPI_STATUS_IGNORE);
    
    free(nD);  
/*...................................................................*/

/*... eVel*/
    size = (*lNel)*ndm;
      
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
  if(ndfD[0] > 0) 
  {
/*... faceRd1 locais*/
    size = (*lNel)*(maxViz+1)*ndfD[0];
    HccaAlloc(short,m,mesh->elm.faceRd1  ,size        
           ,"faceRd1P"  ,_AD_);
    zero(mesh->elm.faceRd1   ,size,"short"  );
/*...................................................................*/

/*... eU0d1 locais*/
     size = (*lNel)*ndfD[0];
     HccaAlloc(DOUBLE,m,mesh->elm.u0D1        ,size        
              ,"eU0d1P"         ,_AD_);
     zero(mesh->elm.u0D1          ,size,DOUBLEC);;
/*...................................................................*/

/*... eUd1 locais*/
     size = (*lNel)*ndfD[0];
     HccaAlloc(DOUBLE,m,mesh->elm.uD1        ,size        
              ,"eUd1P"         ,_AD_);
     zero(mesh->elm.uD1          ,size,DOUBLEC);;
/*...................................................................*/

/*... density locais*/
     size = (*lNel);
     HccaAlloc(DOUBLE,m,mesh->elm.densityUd1.t00  ,size        
              ,"densityP00"            ,_AD_);
     HccaAlloc(DOUBLE,m,mesh->elm.densityUd1.t0  ,size        
              ,"densityP0"            ,_AD_);
     HccaAlloc(DOUBLE,m,mesh->elm.densityUd1.t  ,size        
              ,"densityP"            ,_AD_);
     zero(mesh->elm.densityUd1.t00,size,DOUBLEC);
     zero(mesh->elm.densityUd1.t0 ,size,DOUBLEC);
     zero(mesh->elm.densityUd1.t  ,size,DOUBLEC);
/*...................................................................*/

/*...*/
     comunicateCD(mesh0->elm.faceRd1   ,mesh->elm.faceRd1
                ,mesh0->elm.u0D1      ,mesh->elm.u0D1
                ,mesh0->elm.uD1       ,mesh->elm.uD1
                ,mesh0->elm.densityUd1,mesh->elm.densityUd1
                ,*lNel                ,elLG
                ,maxViz               ,ndfD[0]
                ,0                    ,3);
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if(ndfT[0] > 0)
  {
/*... faceRt1 locais*/
    size = (*lNel)*(maxViz+1)*ndfT[0];
    HccaAlloc(short,m,mesh->elm.faceRt1  ,size        
           ,"faceRt1P"  ,_AD_);
    zero(mesh->elm.faceRt1   ,size,"short"  );
/*...................................................................*/

/*... eU0t1 locais*/
    size = (*lNel)*ndfT[0];
    HccaAlloc(DOUBLE,m,mesh->elm.u0T1        ,size        
             ,"eU0t1P"         ,_AD_);
    zero(mesh->elm.u0T1          ,size,DOUBLEC);;
/*...................................................................*/

/*... eUt1 locais*/
    size = (*lNel)*ndfT[0];
    HccaAlloc(DOUBLE,m,mesh->elm.uT1        ,size        
             ,"eUt1P"         ,_AD_);
    zero(mesh->elm.uT1          ,size,DOUBLEC);;
/*...................................................................*/

/*... densityt1 locais*/
     size = (*lNel);
     HccaAlloc(DOUBLE,m,mesh->elm.densityUt1.t00  ,size        
              ,"densityP00"            ,_AD_);
     HccaAlloc(DOUBLE,m,mesh->elm.densityUt1.t0  ,size        
              ,"densityP0"            ,_AD_);
     HccaAlloc(DOUBLE,m,mesh->elm.densityUt1.t  ,size        
              ,"densityP"            ,_AD_);
     zero(mesh->elm.densityUt1.t00,size,DOUBLEC);
     zero(mesh->elm.densityUt1.t0 ,size,DOUBLEC);
     zero(mesh->elm.densityUt1.t  ,size,DOUBLEC);
/*...................................................................*/

/*...*/
    comunicateCD(mesh0->elm.faceRt1   ,mesh->elm.faceRt1
               ,mesh0->elm.u0T1      ,mesh->elm.u0T1
               ,mesh0->elm.uT1       ,mesh->elm.uT1
               ,mesh0->elm.densityUt1,mesh->elm.densityUt1
               ,*lNel                ,elLG
               ,maxViz               ,ndfT[0]
               ,0                    ,3);
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if(ndfF > 0 || ndfFt > 0)
  {
/*...*/
    allocFuild(m     , mesh
              ,maxViz,*lNel
              ,ndfFt);
/*...................................................................*/

/*...*/
    comunicateFluid(mesh0->elm.faceRvel    ,mesh->elm.faceRvel
                ,mesh0->elm.faceRpres      ,mesh->elm.faceRpres
                ,mesh0->elm.faceRenergy    ,mesh->elm.faceRenergy
                ,mesh0->elm.faceRrho       ,mesh->elm.faceRrho     
                ,mesh0->elm.pressure0      ,mesh->elm.pressure0
                ,mesh0->elm.pressure       ,mesh->elm.pressure
                ,mesh0->elm.energy0        ,mesh->elm.energy0
                ,mesh0->elm.energy         ,mesh->elm.energy
                ,mesh0->elm.temp0          ,mesh->elm.temp0  
                ,mesh0->elm.temp           ,mesh->elm.temp       
                ,mesh0->elm.densityFluid   ,mesh->elm.densityFluid
                ,mesh0->elm.specificHeat   ,mesh->elm.specificHeat 
                ,mesh0->elm.dViscosity     ,mesh->elm.dViscosity 
                ,mesh0->elm.tConductivity  ,mesh->elm.tConductivity 
                ,*lNel                     ,elLG
                ,maxViz                    ,ndfF
                ,ndfFt
                ,0                         ,3);
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if(ndfComb>0)
  {
    allocComb(m      , mesh
            , maxViz , *lNel
            , ndfComb, nSp);

/*...*/
    comunicateComb(mesh0->elm.faceResZcomb ,mesh->elm.faceResZcomb
                  ,mesh0->elm.zComb0       ,mesh->elm.zComb0
                  ,mesh0->elm.zComb        ,mesh->elm.zComb
                  ,mesh0->elm.yFrac0       ,mesh->elm.yFrac0
                  ,mesh0->elm.yFrac        ,mesh->elm.yFrac 
                  ,mesh0->elm.cDiffComb    ,mesh->elm.cDiffComb 
                  ,mesh0->elm.mMolar       ,mesh->elm.mMolar     
                  ,*lNel                   ,elLG
                  ,maxViz                  ,ndfComb
                  ,nSp 
                  ,0                       ,3);
/*...................................................................*/
  }
/*...................................................................*/

/*... nelvon locais*/
  size = (*lNel)*maxViz;
  
  HccaAlloc(INT,m,mesh->elm.adj.nelcon  ,size        
             ,"adjP"      ,_AD_);
  zero(mesh->elm.adj.nelcon,size,INTC);
  
  MPI_Recv( mesh->elm.adj.nelcon,size, MPI_INT,0,10   
          ,mpiVar.comm,MPI_STATUS_IGNORE);    
/*...................................................................*/

/*... nViz locais*/           
  size = *lNel;
  
  HccaAlloc(short,m,mesh->elm.adj.nViz ,size        
           ,"nVizP"         ,_AD_);
  zero(mesh->elm.adj.nViz  ,size      ,"short");
  
  MPI_Recv(mesh->elm.adj.nViz,size, MPI_SHORT,0,11   
          ,mpiVar.comm,MPI_STATUS_IGNORE);      
/*...................................................................*/

/*... no locais*/
  size = (*lNnode)*ndm;
  
  HccaAlloc(DOUBLE,m,mesh->node.x  ,size        
           ,"xnodeP",_AD_);   
  zero(mesh->node.x,size,DOUBLEC);
  
  MPI_Recv(mesh->node.x    ,size, MPI_DOUBLE,0,12   
          ,mpiVar.comm,MPI_STATUS_IGNORE);      
/*...................................................................*/
#endif
}
/*********************************************************************/

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
 
  
#ifdef _MPI_
  mpiVar.comm = 0;
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
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 16/08/2019                                   * 
 *-------------------------------------------------------------------* 
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

#ifdef _MPI_
  mpiWait();
  MPI_Finalize();
#endif
}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 16/08/2019                                   * 
 *-------------------------------------------------------------------* 
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

#ifdef _MPI_
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
  MPI_Barrier(mpiVar.comm);
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
#endif

}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 27/08/2019                                   * 
 *-------------------------------------------------------------------*
 * GOBBALMESHQUALITY : obtem as estatistica globais da malha         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * mQl     -> estatisticas da malha particionada                     *
 * mQl0    -> estatisticas da malha global                           *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * mQl0    -> estatisticas da malha global atualizada                *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void globalMeshQuality(MeshQuality *mQl,MeshQuality *mQl0){

#ifdef _MPI_
  DOUBLE v;
/*... volume */
  MPI_Allreduce(&mQl->volume,&v ,1           ,MPI_DOUBLE
               ,MPI_SUM     ,mpiVar.comm);
  mQl0->volume = v;
/*...................................................................*/

/*... nonOrthMed */
  MPI_Allreduce(&mQl->nonOrthMed,&v      ,1          , MPI_DOUBLE
               ,MPI_SUM         , mpiVar.comm);
  mQl0->nonOrthMed = v/(DOUBLE) mpiVar.nPrcs;  
/*...................................................................*/

/*... nonOrthMAX */
  MPI_Allreduce(&mQl->nonOrthMax,&v      , 1          , MPI_DOUBLE
            ,MPI_MAX            , mpiVar.comm);
  mQl0->nonOrthMax = v;  
/*...................................................................*/

/*... skewMed*/
  MPI_Allreduce(&mQl->skewMed   ,&v      ,1          , MPI_DOUBLE
            ,MPI_SUM            , mpiVar.comm);
  mQl0->skewMed = v/(DOUBLE) mpiVar.nPrcs;  
/*...................................................................*/

/*... skewMax */
  MPI_Allreduce(&mQl->skewMax   ,&v      ,1          , MPI_DOUBLE
            ,MPI_MAX            , mpiVar.comm);
  mQl0->skewMax = v;  
/*...................................................................*/
#endif

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 16/08/2019                                   * 
 *-------------------------------------------------------------------* 
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
 * | u11    u12 ...    u1ndf1 |                                      *
 * |                          | = y(i,ndf1,ndf2)                     * 
 * | undf11 u12 ... undf1ndf2 |                                      *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void dGlobalCel(Memoria *m         ,PartMesh *pMesh
               ,DOUBLE *RESTRICT uG,DOUBLE *RESTRICT uL
               ,INT const numelNov 
               ,short const ndf1   ,short const ndf2)
            
{           
#ifdef _MPI_
  short j,k;
  INT i,lEl;
  INT elG = pMesh->elG ;
  INT *elLG  = pMesh->elLG;
  INT size= elG*ndf1*ndf2;
  DOUBLE *x1=NULL; 
#endif

/*...*/
  if(mpiVar.nPrcs < 2) return;
/*...................................................................*/

/*...*/  
  tm.overHeadGCelMpi = getTimeC() - tm.overHeadGCelMpi;
/*...................................................................*/
   
#ifdef _MPI_

/*...*/
  if(!mpiVar.myId)
    zero(uG,size,DOUBLEC);
/*...................................................................*/

/*...*/
  HccaAlloc(DOUBLE,m,x1,size,"x1Gmpi",false);
  zero(x1,size,DOUBLEC);
/*...................................................................*/

/*... ndf2 = 1*/
  if(ndf2 == 1)
  {
/*... ndf2 = 1*/
    if(ndf1 == 1)
    {
      for(i=0;i<numelNov;i++)
      {
        lEl     = elLG[i];
        x1[lEl] = uL[i]; 
      }
    }
/*...................................................................*/

/*...*/
    else
    {
      for(i=0;i<numelNov;i++)
      {
        lEl    = elLG[i];
        for(j=0;j<ndf1;j++)
          MAT2D(lEl,j,x1,ndf1) = MAT2D(i,j,uL,ndf1); 
      }
    }
/*...................................................................*/
  }  
/*...................................................................*/

/*...*/
  else
  {
    for(i=0;i<numelNov;i++)
    {
      lEl    = elLG[i];
      for(j=0;j<ndf1;j++)
        for(k=0;k<ndf2;k++)
        MAT3D(lEl,j,k,x1,ndf1,ndf2) = MAT3D(i,j,k,uL,ndf1,ndf2); 
    }
  }
/*...................................................................*/

/*...*/
  MPI_Reduce(x1,uG,size,MPI_DOUBLE,MPI_SUM,0,mpiVar.comm);
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
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 19/11/2019                                   * 
 *-------------------------------------------------------------------*
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
                ,DOUBLE *RESTRICT uG,DOUBLE *RESTRICT uL
                ,short const ndf1   ,short const ndf2)
            
{           
#ifdef _MPI_  
  INT i,j,k,lNo;
  INT nno    = pMesh->nno1;
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

#ifdef _MPI_

/*...*/
  if(!mpiVar.myId)
    zero(uG,size,DOUBLEC);
/*...................................................................*/

/*...*/
  HccaAlloc(DOUBLE,m,x1,size,"x1Gmpi",false);
  zero(x1,size,DOUBLEC);
/*...................................................................*/

/*... ndf2 = 1*/
  if(ndf2 == 1)
  {
/*...*/
    if(ndf1 == 1)
    {
      for(i=0;i<nno;i++){
        lNo    = noLG[i];
        x1[lNo] = uL[i]; 
      }
    }
/*...................................................................*/

/*...*/
    else
    {
      for(i=0;i<nno;i++)
      {
        lNo    = noLG[i];
        for(j=0;j<ndf1;j++)
          MAT2D(lNo,j,x1,ndf1) = MAT2D(i,j,uL,ndf1); 
      }
    } 
/*...................................................................*/
  } 
/*...................................................................*/
  else
  {
    for(i=0;i<nno;i++)
    {
      lNo    = noLG[i];
      for(j=0;j<ndf1;j++)
        for(k=0;k<ndf2;k++)
          MAT3D(lNo,j,k,x1,ndf1,ndf2) = MAT3D(i,j,k,uL,ndf1,ndf2); 
    }
  }
/*...................................................................*/

/*...*/
  MPI_Reduce(x1,uG,size,MPI_DOUBLE,MPI_SUM,0,mpiVar.comm);
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
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 16/08/2019                                   * 
 *-------------------------------------------------------------------*
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
void getBuffer(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
              ,INT *RESTRICT fMap    ,INT const nRcvs){
  
  INT i,lNeq;

  for(i=0;i<nRcvs;i++){
    lNeq    = fMap[i];
    x[lNeq] = xb[i];
  }
    
} 
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 16/08/2019                                   * 
 *-------------------------------------------------------------------* 
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
void makeBuffer(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
               ,INT *RESTRICT fMap    ,INT const nSends){
    
  INT i,lNeq;

  for(i=0;i<nSends;i++)
  {
    lNeq   = fMap[i];
    xb[i]  = x[lNeq];
  }
  
} 
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 16/08/2019                                   * 
 *-------------------------------------------------------------------*
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
void getBufferCel(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
              ,INT *RESTRICT fMap       ,INT const nRcvs
              ,short const ndf1         ,short const ndf2)
{
  short j,k;
  INT i,lCel;

/*... ndf = 1 */  
  if( ndf2 == 1)
  {
/*... ndm = 3*/
    if( ndf1 == 3)
    {
      for(i=0;i<nRcvs;i++)
      {
        lCel    = fMap[i];
        MAT2D(lCel,0,x,3) = MAT2D(i,0,xb,3);
        MAT2D(lCel,1,x,3) = MAT2D(i,1,xb,3);
        MAT2D(lCel,2,x,3) = MAT2D(i,2,xb,3);
      }
    }
/*..................................................................*/

/*... ndm = 2*/
    else if( ndf1 == 2)
    {
      for(i=0;i<nRcvs;i++)
      {
        lCel    = fMap[i];
        MAT2D(lCel,0,x,2) = MAT2D(i,0,xb,2);
        MAT2D(lCel,1,x,2) = MAT2D(i,1,xb,2);
      }
    }
/*..................................................................*/

/*... ndm = 1*/
    else if(ndf1 == 1)
    {
      for(i=0;i<nRcvs;i++)
      {
        lCel    = fMap[i];
        x[lCel] = xb[i];
      }
    }
/*..................................................................*/

/*... */
    else
    {
      for(i=0;i<nRcvs;i++)
      {
        lCel    = fMap[i];
        for(j=0;j<ndf1;j++)
          MAT2D(lCel,j,x,ndf1) = MAT2D(i,j,xb,ndf1);
      }
    }
/*.....................................................................*/

  }
/*..................................................................*/
   
/*...*/
  else
  {
    for(i=0;i<nRcvs;i++)
    {
      lCel    = fMap[i];
      for(j=0;j<ndf1;j++)
        for(k=0;k<ndf2;k++)
          MAT3D(lCel,j,k,x,ndf1,ndf2) = MAT3D(i,j,k,xb,ndf1,ndf2);
    }
  }
/*..................................................................*/
 
} 
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 22/08/2019                                   * 
 *-------------------------------------------------------------------* 
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
void makeBufferCel(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
                  ,INT *RESTRICT fMap    ,INT const nSends
                  ,short const ndf1      ,short const ndf2){
  short j,k;
  INT i,lCel;

/*... ndf = 1 */  
  if( ndf2 == 1)
  {
/*... ndm = 3*/
    if( ndf1 == 3)
    {
      for(i=0;i<nSends;i++)
      {
        lCel    = fMap[i];
        MAT2D(i,0,xb,3) = MAT2D(lCel,0,x,3);
        MAT2D(i,1,xb,3) = MAT2D(lCel,1,x,3);
        MAT2D(i,2,xb,3) = MAT2D(lCel,2,x,3);
      }
    }
/*.....................................................................*/

/*... ndm = 2*/
    else if( ndf1 == 2)
    {
      for(i=0;i<nSends;i++)
      {
        lCel    = fMap[i];
        MAT2D(i,0,xb,2) = MAT2D(lCel,0,x,2);
        MAT2D(i,1,xb,2) = MAT2D(lCel,1,x,2);
      }
    }
/*.....................................................................*/

/*... ndm = 1*/
    else if(ndf1 == 1)
    {
      for(i=0;i<nSends;i++)
      {
        lCel    = fMap[i];
        xb[i]   = x[lCel];
      }
    }
/*.....................................................................*/

/*... */
    else
    {
      for(i=0;i<nSends;i++)
      {
        lCel    = fMap[i];
        for(j=0;j<ndf1;j++)
          MAT2D(i,j,xb,ndf1) = MAT2D(lCel,j,x,ndf1);
      }
    }
/*.....................................................................*/

  }
/*.....................................................................*/ 

/*...*/
  else
  {
    for(i=0;i<nSends;i++)
    {
      lCel    = fMap[i];
      for(j=0;j<ndf1;j++)
        for(k=0;k<ndf2;k++)
          MAT3D(i,j,k,xb,ndf1,ndf2) = MAT3D(lCel,j,k,x,ndf1,ndf2);
    }
  }
/*....................................................................*/ 

} 
/********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 22/08/2019                                   * 
 *-------------------------------------------------------------------*
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
void dGetBufferNod(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
              ,INT *RESTRICT fMap        ,INT const nRcvs
              ,short const ndf1          ,short const ndf2){

  short j,k;
  INT i,lNod;

/*... ndf2 = 1 */  
  if( ndf2 == 1)
  {
/*... ndf1 = 3*/
    if( ndf1 == 3)
    {
      for(i=0;i<nRcvs;i++)
      {
        lNod    = fMap[i];
        MAT2D(lNod,0,x,3) += MAT2D(i,0,xb,3);
        MAT2D(lNod,1,x,3) += MAT2D(i,1,xb,3);
        MAT2D(lNod,2,x,3) += MAT2D(i,2,xb,3);
      }
    }
/*...................................................................*/

/*... ndf1 = 2*/
    else if( ndf1 == 2)
    {
      for(i=0;i<nRcvs;i++)
      {
        lNod    = fMap[i];
        MAT2D(lNod,0,x,2) += MAT2D(i,0,xb,2);
        MAT2D(lNod,1,x,2) += MAT2D(i,1,xb,2);
      }
    }
/*...................................................................*/
    
/*... ndf1 = 1*/
    else if(ndf1 == 1)
    {
      for(i=0;i<nRcvs;i++)
      {
        lNod    = fMap[i];
        x[lNod] += xb[i];
      }
    }
/*...................................................................*/
    
/*... */
    else
    {
      for(i=0;i<nRcvs;i++)
      {
        lNod    = fMap[i];
        for(j=0;j<ndf1;j++)
          MAT2D(lNod,j,x,ndf1) += MAT2D(i,j,xb,ndf1);
      }
    }
/*.....................................................................*/

  }
/*.....................................................................*/ 

/*...*/
  else
  {
    for(i=0;i<nRcvs;i++)
    {
      lNod    = fMap[i];
      for(j=0;j<ndf1;j++)
        for(k=0;k<ndf2;k++)
          MAT3D(lNod,j,k,x,ndf1,ndf2) += MAT3D(i,j,k,xb,ndf1,ndf2);
    }
  }
/*....................................................................*/ 


} 
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 22/08/2019                                   * 
 *-------------------------------------------------------------------* 
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
void dMakeBufferNod(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
                  ,INT *RESTRICT fMap    ,INT const nSends
                  ,short const ndf1          ,short const ndf2){
  short j,k;  
  INT i,lNo;

/*... ndf2 = 1 */  
  if( ndf2 == 1)
  {
/*... ndf1 = 3*/
    if( ndf1 == 3)
    {
      for(i=0;i<nSends;i++)
      {
        lNo     = fMap[i];
        MAT2D(i,0,xb,3) = MAT2D(lNo,0,x,3);
        MAT2D(i,1,xb,3) = MAT2D(lNo,1,x,3);
        MAT2D(i,2,xb,3) = MAT2D(lNo,2,x,3);
      }
    }
/*... ndf2 = 2*/
    else if( ndf1 == 2)
    {
      for(i=0;i<nSends;i++)
      {
        lNo    = fMap[i];
        MAT2D(i,0,xb,2) = MAT2D(lNo,0,x,2);
        MAT2D(i,1,xb,2) = MAT2D(lNo,1,x,2);
      }
    }
/*... ndf3 = 1*/
    else if(ndf1 == 1)
    { 
      for(i=0;i<nSends;i++)
      {
        lNo    = fMap[i];
        xb[i]  = x[lNo];
      } 
    }
/*..................................................................*/ 

/*... */
    else
    {
      for(i=0;i<nSends;i++)
      {
        lNo    = fMap[i];
        for(j=0;j<ndf1;j++)
          MAT2D(i,j,xb,ndf1) = MAT2D(lNo,j,x,ndf1);
      }
    }
/*.....................................................................*/

  }
/*.....................................................................*/ 

/*...*/
  else
  {
    for(i=0;i<nSends;i++)
    {
      lNo    = fMap[i];
      for(j=0;j<ndf1;j++)
        for(k=0;k<ndf2;k++)
          MAT3D(i,j,k,xb,ndf1,ndf2) = MAT3D(lNo,j,k,x,ndf1,ndf2);
    }
  }
/*....................................................................*/
 
} 
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 16/08/2019                                   * 
 *-------------------------------------------------------------------* 
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
void iGetBufferNod(INT *RESTRICT x     ,INT *RESTRICT xb
                 ,INT *RESTRICT fMap  ,INT const nRcvs
                 ,short const ndf1    ,short const ndf2){

  INT i,lNod;

/*... ndf2 = 1 */  
  if( ndf2 == 1){
/*... ndf1 = 3*/
    if( ndf1 == 3)
    {
      for(i=0;i<nRcvs;i++)
      {
        lNod    = fMap[i];
        MAT2D(lNod,0,x,3) += MAT2D(i,0,xb,3);
        MAT2D(lNod,1,x,3) += MAT2D(i,1,xb,3);
        MAT2D(lNod,2,x,3) += MAT2D(i,2,xb,3);
      }
    }
/*...................................................................*/

/*... ndf1 = 2*/
    else if( ndf1 == 2)
    {
      for(i=0;i<nRcvs;i++)
      {
        lNod    = fMap[i];
        MAT2D(lNod,0,x,2) += MAT2D(i,0,xb,2);
        MAT2D(lNod,1,x,2) += MAT2D(i,1,xb,2);
      }
    }
/*...................................................................*/

/*... ndf1 = 1*/
    else if(ndf1 == 1)
    {
      for(i=0;i<nRcvs;i++)
      {
        lNod    = fMap[i];
        x[lNod] += xb[i];
      }
    } 
/*...................................................................*/

/*...*/
    else
    {
        printf("Nao implementado!!");
        mpiStop();
        exit(EXIT_FAILURE);
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  else
  {
    printf("Nao implementado!!");
    mpiStop();
    exit(EXIT_FAILURE);
  }
/*...................................................................*/
} 
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 16/08/2019                                   * 
 *-------------------------------------------------------------------* 
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
void iMakeBufferNod(INT *RESTRICT x       ,INT *RESTRICT xb
                   ,INT *RESTRICT fMap    ,INT const nSends
                   ,short const ndf1      ,short const ndf2){
    
  INT i,lNo;

/*... ndf = 1 */  
  if( ndf2 == 1){
/*... ndf = 3*/
    if( ndf1 == 3)
    {
      for(i=0;i<nSends;i++)
      {
        lNo     = fMap[i];
        MAT2D(i,0,xb,3) = MAT2D(lNo,0,x,3);
        MAT2D(i,1,xb,3) = MAT2D(lNo,1,x,3);
        MAT2D(i,2,xb,3) = MAT2D(lNo,2,x,3);
      }
    }
/*... ndf = 2*/
    else if( ndf1 == 2)
    {
      for(i=0;i<nSends;i++)
      {
        lNo    = fMap[i];
        MAT2D(i,0,xb,2) = MAT2D(lNo,0,x,2);
        MAT2D(i,1,xb,2) = MAT2D(lNo,1,x,2);
      }
    }
/*... ndf = 1*/
    else if(ndf1 == 1) 
    {
      for(i=0;i<nSends;i++)
      {
        lNo    = fMap[i];
        xb[i]  = x[lNo];
      }
    }
    else
    {
        printf("Nao implementado!!");
        mpiStop();
        exit(EXIT_FAILURE);
    }
  }
  else
  {
    printf("Nao implementado!!");
    mpiStop();
    exit(EXIT_FAILURE);
  }
/*..................................................................*/  
} 
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
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
  void comunicateNeq(Interface *iNeq,DOUBLE *RESTRICT x){

#ifdef _MPI_
  INT nRcvs  = iNeq->nRcvs;
  INT nSends = iNeq->nSends;
  INT k,kk;
  unsigned short partId,nPart,nVizParts = iNeq->nVizPart; 

  if(mpiVar.nPrcs < 2 ) return;

/*...*/
  tm.overHeadNeqMpi = getTimeC() - tm.overHeadNeqMpi;
/*...................................................................*/    

/*... gerando o buffer de envio*/
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
  makeBuffer(x,&iNeq->xb[nRcvs],&iNeq->fMap[nRcvs],nSends);
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
/*...................................................................*/    

/*...*/  
  for(nPart=0;nPart < nVizParts;nPart++){
    partId = iNeq->vizPart[nPart];
/*... enviando para as particoes vizinhas*/
    k      = iNeq->iaSends[nPart];
    kk     = iNeq->iaSends[nPart+1] - iNeq->iaSends[nPart];
    tm.overHeadSendMpi = getTimeC() - tm.overHeadSendMpi;
    MPI_Isend(&iNeq->xb[k]    ,kk
             ,MPI_DOUBLE      ,partId    
             ,mpiVar.myId     ,mpiVar.comm
             ,&mpiVar.sendRequest[nPart]);
    tm.overHeadSendMpi = getTimeC() - tm.overHeadSendMpi;
/*... recebimentos das particoes vizinhas*/
    k      = iNeq->iaRcvs[nPart];
    kk     = iNeq->iaRcvs[nPart+1] - iNeq->iaRcvs[nPart];
    tm.overHeadRecvMpi = getTimeC() - tm.overHeadRecvMpi;
    MPI_Irecv(&iNeq->xb[k]   ,kk
             ,MPI_DOUBLE     ,partId    
             ,partId         ,mpiVar.comm
             ,&mpiVar.recvRequest[nPart]);
    tm.overHeadRecvMpi = getTimeC() - tm.overHeadRecvMpi;
  }  
/*...................................................................*/    
    
/*... espera o recebimentos dos dados*/
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.recvRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
/*...................................................................*/    
  
/*...*/
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
  getBuffer(x,iNeq->xb,iNeq->fMap,nRcvs);
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
/*...................................................................*/    

/*... espera o envio de todos os dados*/
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.sendRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
/*...................................................................*/    

/*...*/
  tm.overHeadNeqMpi    = getTimeC() - tm.overHeadNeqMpi;
/*...................................................................*/    

#endif
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*  
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
  void comunicateCel(Interface *iCel ,DOUBLE *RESTRICT x
                    ,short const ndf1,short const ndf2){

#ifdef _MPI_
  INT nRcvs  = iCel->nRcvs;
  INT nSends = iCel->nSends;
  INT k,kk,nst=ndf1*ndf2;
  unsigned short partId,nPart,nVizParts = iCel->nVizPart; 

/*...*/
  if(mpiVar.nPrcs < 2 ) return;
/*...................................................................*/

/*...*/
  tm.overHeadCelMpi    = getTimeC() - tm.overHeadCelMpi;
/*...................................................................*/    

/*... gerando o buffer de envio*/
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
  makeBufferCel(x                     ,&iCel->xb[nRcvs*nst]
                ,&iCel->fMap[nRcvs]    ,nSends
                ,ndf1                  ,ndf2);
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
/*...................................................................*/    

/*...*/   
  for(nPart=0;nPart < nVizParts;nPart++){
    partId = iCel->vizPart[nPart];
/*... enviando para as particoes vizinhas*/
    k      = iCel->iaSends[nPart];
    kk     = iCel->iaSends[nPart+1] - iCel->iaSends[nPart];
    k      *= nst;
    kk     *= nst;
    tm.overHeadSendMpi = getTimeC() - tm.overHeadSendMpi;
    MPI_Isend(&iCel->xb[k]    ,kk
             ,MPI_DOUBLE      ,partId    
             ,mpiVar.myId     ,mpiVar.comm
             ,&mpiVar.sendRequest[nPart]);
    tm.overHeadSendMpi = getTimeC() - tm.overHeadSendMpi;
/*... recebimentos das particoes vizinhas*/
    k      = iCel->iaRcvs[nPart];
    kk     = iCel->iaRcvs[nPart+1] - iCel->iaRcvs[nPart];
    k      *= nst;
    kk     *= nst;
    tm.overHeadRecvMpi = getTimeC() - tm.overHeadRecvMpi;
    MPI_Irecv(&iCel->xb[k]   ,kk
             ,MPI_DOUBLE     ,partId    
             ,partId         ,mpiVar.comm
             ,&mpiVar.recvRequest[nPart]);
    tm.overHeadRecvMpi = getTimeC() - tm.overHeadRecvMpi;
  }  
/*...................................................................*/    

/*... espera o recebimentos dos dados*/
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.recvRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
/*...................................................................*/    
  
/*...*/
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
  getBufferCel(x         ,iCel->xb
              ,iCel->fMap,nRcvs
              ,ndf1      ,ndf2);
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
/*...................................................................*/    

/*... espera o envio de todos os dados*/
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.sendRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
/*...................................................................*/    

/*...*/
  tm.overHeadCelMpi    = getTimeC() - tm.overHeadCelMpi;
/*...................................................................*/    

#endif
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 27/08/2019                                   *
 *-------------------------------------------------------------------*
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
  void dComunicateNod(InterfaceNo *iNo ,DOUBLE *RESTRICT x
                     ,short const ndf1 ,short const ndf2){

#ifdef _MPI_
  INT nCom   = iNo->nCom;
  INT k,kk,nst=ndf1*ndf2;
  unsigned short partId,nPart,nVizParts = iNo->nVizPart; 

/*...*/
  if(mpiVar.nPrcs < 2 ) return;
/*...................................................................*/

/*...*/
  tm.overHeadNodMpi    = getTimeC() - tm.overHeadNodMpi;
/*...................................................................*/    

/*... gerando o buffer de envio*/
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
  dMakeBufferNod(x                     ,&iNo->xb[nCom*nst]
                ,iNo->fMap             ,nCom   
                ,ndf1                  ,ndf2);
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
/*...................................................................*/    

/*...*/   
  for(nPart=0;nPart < nVizParts;nPart++){
    k      = iNo->iaComNo[nPart];
    kk     = iNo->iaComNo[nPart+1] - iNo->iaComNo[nPart];
    partId = iNo->vizPart[nPart];
    k      *= nst;
    kk     *= nst;
/*... enviando para as particoes vizinhas*/
    tm.overHeadSendMpi = getTimeC() - tm.overHeadSendMpi;
    MPI_Isend(&iNo->xb[nCom*nst+k]     ,kk
             ,MPI_DOUBLE               ,partId    
             ,mpiVar.myId              ,mpiVar.comm
             ,&mpiVar.sendRequest[nPart]);
    tm.overHeadSendMpi = getTimeC() - tm.overHeadSendMpi;
/*... recebimentos das particoes vizinhas*/
    tm.overHeadRecvMpi = getTimeC() - tm.overHeadRecvMpi;
    MPI_Irecv(&iNo->xb[k]   ,kk
             ,MPI_DOUBLE     ,partId    
             ,partId         ,mpiVar.comm
             ,&mpiVar.recvRequest[nPart]);
    tm.overHeadRecvMpi = getTimeC() - tm.overHeadRecvMpi;
  }  
/*...................................................................*/    

/*... espera o recebimentos dos dados*/
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.recvRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
/*...................................................................*/    
  
/*...*/
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
  dGetBufferNod(x        ,iNo->xb
              ,iNo->fMap ,nCom
              ,ndf1      ,ndf2);
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
/*...................................................................*/    

/*... espera o envio de todos os dados*/
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.sendRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
/*...................................................................*/    

/*...*/
  tm.overHeadNodMpi    = getTimeC() - tm.overHeadNodMpi;
/*...................................................................*/    

#endif
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 27/08/2019                                   *
 *-------------------------------------------------------------------*
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
  void iComunicateNod(InterfaceNo *iNo ,INT *RESTRICT x
                     ,short const ndf1,short const ndf2){

#ifdef _MPI_
  INT nCom   = iNo->nCom;
  INT k,kk,nst=ndf1*ndf2;
  unsigned short partId,nPart,nVizParts = iNo->nVizPart; 

/*...*/
  if(mpiVar.nPrcs < 2 ) return;
/*...................................................................*/

/*...*/
  tm.overHeadNodMpi    = getTimeC() - tm.overHeadNodMpi;
/*...................................................................*/

/*... gerando o buffer de envio*/
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
  iMakeBufferNod(x                     ,&iNo->xi[nCom*nst]
                ,iNo->fMap             ,nCom   
                ,ndf1                  ,ndf2);
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
/*...................................................................*/    

/*...*/   
  for(nPart=0;nPart < nVizParts;nPart++){
    k      = iNo->iaComNo[nPart];
    kk     = iNo->iaComNo[nPart+1] - iNo->iaComNo[nPart];
    partId = iNo->vizPart[nPart];
    k      *= nst;
    kk     *= nst;
/*... enviando para as particoes vizinhas*/
    tm.overHeadSendMpi = getTimeC() - tm.overHeadSendMpi;
    MPI_Isend(&iNo->xi[nCom*nst+k]     ,kk
             ,MPI_INT                  ,partId    
             ,mpiVar.myId              ,mpiVar.comm
             ,&mpiVar.sendRequest[nPart]);
    tm.overHeadSendMpi = getTimeC() - tm.overHeadSendMpi;
/*... recebimentos das particoes vizinhas*/
    tm.overHeadRecvMpi = getTimeC() - tm.overHeadRecvMpi;
    MPI_Irecv(&iNo->xi[k]    ,kk
             ,MPI_INT        ,partId    
             ,partId         ,mpiVar.comm
             ,&mpiVar.recvRequest[nPart]);
    tm.overHeadRecvMpi = getTimeC() - tm.overHeadRecvMpi;
  }  
/*...................................................................*/    

/*... espera o recebimentos dos dados*/
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.recvRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
/*...................................................................*/    
  
/*...*/
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
  iGetBufferNod(x         ,iNo->xi
              ,iNo->fMap ,nCom
              ,ndf1      ,ndf2);
  tm.overHeadBufferMpi = getTimeC() - tm.overHeadBufferMpi;
/*...................................................................*/    

/*... espera o envio de todos os dados*/
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
  mpiVar.ierr = MPI_Waitall(nVizParts,mpiVar.sendRequest,mpiVar.status);
  if (mpiVar.ierr != MPI_SUCCESS) {
    MPI_Error_string(mpiVar.ierr, mpiVar.errBuffer,&mpiVar.lString);
    printf("!%s!\n",mpiVar.errBuffer);
  }
  tm.overHeadWaitMpi = getTimeC() - tm.overHeadWaitMpi;
/*...................................................................*/    

/*...*/
  tm.overHeadNodMpi    = getTimeC() - tm.overHeadNodMpi;
/*...................................................................*/    

#endif
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 12/08/2019                                   *
 * ------------------------------------------------------------------*
 * COMUNICATEMESH : gera o mapa de inteface dos elementes e comunica * 
 * as particoes para os vizinhos                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * m       -> vetor de memoria                                       *
 * tModel  -> modelo de turubulencia                                 *
 * mesh0   -> malha original                                         *
 * mesh    -> nao definido                                           *
 * pMesh   -> particionamento                                        *
 * ldD1    -> cargas nas faces                                       *
 * ldT1    -> cargas nas faces                                       *
 * ldD1    -> cargas nas faces                                       *
 * ldVel   -> cargas nas faces                                       *
 * ldPres  -> cargas nas faces                                       *
 * ldPresC -> cargas nas faces                                       *
 * ldEnergy-> cargas nas faces                                       *
 * ldTemp  -> cargas nas faces                                       *
 * ldKturb -> cargas nas faces                                       *
 * ldZcomb -> cargas nas faces                                       *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * mesh    -> malha particionada e comunicada para o escravos        *
 * pMesh   -> particionamento atualizado                             *
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 * todos os processos tem tModel
 *-------------------------------------------------------------------* 
 *********************************************************************/
void comunicateMesh(Memoria *m        ,Combustion *cModel
                   ,Turbulence *tModel
                   ,Mesh *mesh0     ,Mesh *mesh
                   ,PartMesh *pMesh 
                   ,Loads *ldD1     ,Loads *ldT1
                   ,Loads *ldVel    ,Loads *ldPres
                   ,Loads *ldPresC  ,Loads *ldEnergy
                   ,Loads *ldTemp   ,Loads *ldKturb
                   ,Loads *ldZcomb) 
{

#ifdef _MPI_
  bool *aux1 = NULL,*aux2 = NULL,*aux3=NULL;
  bool fTurb          = tModel->fTurb,
       fWall          = tModel->fWall, 
       fTurbStruct    = tModel->fTurbStruct,
       fDynamic       = tModel->fDynamic;
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
  unsigned short nPrcs=mpiVar.nPrcs,myId=mpiVar.myId; 
  short st=0;
  short ndm,maxViz,maxNo,numat;
  short ndfD[MAX_DIF_EQ],ndfT[MAX_TRANS_EQ],ndfF,ndfFt,ntn;
  short maxNdf,ndfVel,ndfComb,nSp;

/*...*/
  ndfF    = 0;
  ndfFt   = 0;
  ndfComb = cModel->nComb;
  nSp     = cModel->nOfSpecies;
  for(i=0;i<MAX_TRANS_EQ;i++)
    ndfT[i]  = 0;
  for(i=0;i<MAX_DIF_EQ;i++)
    ndfD[i]  = 0;
/*...................................................................*/

/*... alocando varaivel auxiliar*/  
  if(!myId)
  {
/*... gera o incidencia do elementos*/
    mesh0->noIncid.nincid = (INT *) malloc(mesh0->nnode*iSize);
    zero(mesh0->noIncid.nincid,mesh0->nnode,INTC);
    nodeGrade(mesh0->elm.node,mesh0->noIncid.nincid
             ,mesh0->elm.nen ,&maxGrade
             ,mesh0->nnode   ,mesh0->numel
             ,mesh0->maxNo);
/*...................................................................*/
 
/*...*/
    mesh0->noIncid.incid = (INT *) malloc(mesh0->nnode*maxGrade*iSize);
    zero(mesh0->noIncid.incid,mesh0->nnode*maxGrade,INTC);
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

/*...*/ 
  MPI_Bcast(&st,1,MPI_SHORT,0,mpiVar.comm);
  if( st == -1)
  {
    mpiStop();
    exit(EXIT_MESH_COMM); 
  }
/*...................................................................*/

/*... variaveis compartilhadas*/
  if(!myId)
  { 
    pMesh->nnG  = mesh0->nnode;
    pMesh->elG  = mesh0->numel;
    mesh->ntn   = mesh0->ntn;
    mesh->ndm   = mesh0->ndm;
    mesh->numat = mesh0->numat;
    mesh->maxNo = mesh0->maxNo;
    mesh->maxViz= mesh0->maxViz;
    for(i=0;i<MAX_TRANS_EQ;i++)
      mesh->ndfT[i]  = mesh0->ndfT[i];
    for(i=0;i<MAX_DIF_EQ;i++)
      mesh->ndfD[i]  = mesh0->ndfD[i];
    mesh->ndfF  = mesh0->ndfF;
    mesh->ndfFt = mesh0->ndfFt;
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
  MPI_Bcast(&mesh->ndfF  ,1           ,MPI_SHORT,0,mpiVar.comm);
  MPI_Bcast(&mesh->ndfFt ,1           ,MPI_SHORT,0,mpiVar.comm);
  MPI_Bcast(&mesh->ntn   ,1           ,MPI_SHORT,0,mpiVar.comm);

/*... loadsD1*/
  for(i=0;i<MAXLOAD;i++)
  {
    MPI_Bcast(&ldD1[i].fUse  ,1 ,MPI_C_BOOL,0,mpiVar.comm);
    MPI_Bcast(&ldD1[i].type  ,1 ,MPI_SHORT ,0,mpiVar.comm);
    MPI_Bcast(&ldD1[i].np    ,1 ,MPI_SHORT ,0,mpiVar.comm);
    MPI_Bcast(&ldD1[i].par   ,MAXLOADPARAMETER ,MPI_DOUBLE
              ,0,mpiVar.comm);
  }
/*...................................................................*/

/*... loadsT1*/
  for(i=0;i<MAXLOAD;i++)
  {
    MPI_Bcast(&ldT1[i].fUse  ,1 ,MPI_C_BOOL,0,mpiVar.comm);
    MPI_Bcast(&ldT1[i].type  ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldT1[i].np    ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldT1[i].par   ,MAXLOADPARAMETER ,MPI_DOUBLE
             ,0,mpiVar.comm);
  }
/*...................................................................*/

/*... Fluid*/
  for(i=0;i<MAXLOAD;i++)
  {
/*... vel*/
    MPI_Bcast(&ldVel[i].fUse  ,1 ,MPI_C_BOOL,0,mpiVar.comm);
    MPI_Bcast(&ldVel[i].type  ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldVel[i].np    ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldVel[i].par   ,MAXLOADPARAMETER ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldVel[i].vel    ,3               ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldVel[i].density,1              ,MPI_DOUBLE
              ,0,mpiVar.comm);
/*... pres*/
    MPI_Bcast(&ldPres[i].fUse  ,1 ,MPI_C_BOOL,0,mpiVar.comm);
    MPI_Bcast(&ldPres[i].type  ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldPres[i].np    ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldPres[i].par   ,MAXLOADPARAMETER ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldPres[i].vel    ,3               ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldPres[i].density,1              ,MPI_DOUBLE
              ,0,mpiVar.comm);
/*... presC*/
    MPI_Bcast(&ldPresC[i].fUse  ,1 ,MPI_C_BOOL,0,mpiVar.comm);
    MPI_Bcast(&ldPresC[i].type  ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldPresC[i].np    ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldPresC[i].par   ,MAXLOADPARAMETER ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldPresC[i].vel    ,3               ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldPresC[i].density,1              ,MPI_DOUBLE
              ,0,mpiVar.comm);
/*... Energy*/
    MPI_Bcast(&ldEnergy[i].fUse  ,1 ,MPI_C_BOOL,0,mpiVar.comm);
    MPI_Bcast(&ldEnergy[i].type  ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldEnergy[i].np    ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldEnergy[i].par   ,MAXLOADPARAMETER ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldEnergy[i].vel    ,3               ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldEnergy[i].density,1              ,MPI_DOUBLE
              ,0,mpiVar.comm);
/*... Temp*/
    MPI_Bcast(&ldTemp[i].fUse  ,1 ,MPI_C_BOOL,0,mpiVar.comm);
    MPI_Bcast(&ldTemp[i].type  ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldTemp[i].np    ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldTemp[i].par   ,MAXLOADPARAMETER ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldTemp[i].vel    ,3               ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldTemp[i].density,1              ,MPI_DOUBLE
              ,0,mpiVar.comm);
/*... Kturb*/
    MPI_Bcast(&ldKturb[i].fUse  ,1 ,MPI_C_BOOL,0,mpiVar.comm);
    MPI_Bcast(&ldKturb[i].type  ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldKturb[i].np    ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldKturb[i].par   ,MAXLOADPARAMETER ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldKturb[i].vel    ,3               ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldKturb[i].density,1              ,MPI_DOUBLE
              ,0,mpiVar.comm);
/*... Zcomb*/
    MPI_Bcast(&ldZcomb[i].fUse  ,1 ,MPI_C_BOOL,0,mpiVar.comm);
    MPI_Bcast(&ldZcomb[i].type  ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldZcomb[i].np    ,1 ,MPI_SHORT,0,mpiVar.comm);
    MPI_Bcast(&ldZcomb[i].par   ,MAXLOADPARAMETER ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldZcomb[i].vel    ,3               ,MPI_DOUBLE
              ,0,mpiVar.comm);
    MPI_Bcast(&ldZcomb[i].density,1              ,MPI_DOUBLE
              ,0,mpiVar.comm);
  }
/*...................................................................*/

/*...*/
  ntn    = mesh->ntn;
  ndm    = mesh->ndm;
  ndfF   = mesh->ndfF;
  ndfFt  = mesh->ndfFt;
  numat  = mesh->numat;
  maxNo  = mesh->maxNo;
  maxViz = mesh->maxViz;
  for(i=0;i<MAX_TRANS_EQ;i++)
    ndfT[i]  = mesh->ndfT[i];
  for(i=0;i<MAX_DIF_EQ;i++)
    ndfD[i]  = mesh->ndfD[i];

/*...*/
  ndfVel = max(mesh->ndfF - 1,mesh->ndfFt - 2);
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
  for(nPart=0;nPart<nPrcs;nPart++)
  {
/*... so o processo master gera os mapas*/
    if(!myId)
    { 
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
      nVizPartNo = getMapVizNo(pMesh->ep            ,noLG
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

//    fprintf(stderr,"rank %d nNodeNov %d nNodeOv %d\n",nPart,nNodeNov,nNodeOv);
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
      if(!nPart)
      {
        mesh->nnode     = lNnode;
        mesh->numel     = lNel;
        mesh->numelNov  = numelNov;
        mesh->nnodeOv   = nNodeOv;
        mesh->nnodeNov  = nNodeNov;
        pMesh->nno1     = nno1;
/*... particionamento*/     

/*... partViz (El)*/
        prMaster(m                 , mesh0
                , mesh             , pMesh
                , fMap             , fMapNo
                , iaSends          , iaRcvs
                , iaComNo            
                , elLG             , elGL
                , noLG             , noGL
                , maxVizPart       , maxVizPartNo
                , lNel             , lNnode
                , nRcvs            , nSends
                , nComNo           , maxNo  
                , ndfD             , ndfT
                , ndfF             , ndfFt
                , ndfComb          , nSp
                , ndm              , maxViz
                , nVizPart         , nVizPartNo);  
      }
/*...................................................................*/

/*... comunicando o partionamento local para o processos escravos 
      (master)*/
      else
      {
        sendPart(m          , mesh0
                , mesh      , pMesh
                , fMap      , fMapNo
                , iaSends   , iaRcvs
                , iaComNo                   
                , elLG      , elGL
                , noLG      , noGL
                , numelOv   , nNodeOv
                , nno1 
                , maxVizPart, maxVizPartNo
                , lNel      , lNnode
                , nRcvs     , nSends
                , nComNo    , maxNo  
                , ndfD      , ndfT
                , ndfF      , ndfFt
                , ndfComb   , nSp
                , ndm       , maxViz
                , nVizPart  , nVizPartNo
                , nPart );      
      } 
/*...................................................................*/
    }
/*...................................................................*/
      
/*... recebendo o partionamento local do processo master (escravos)*/
    if(nPart == myId && myId) 
    {
      recvPart(m           , mesh0
             , mesh        , pMesh
             , fMap        , fMapNo
             , iaSends     , iaRcvs
             , iaComNo            
             , elLG        , elGL
             , noLG        , noGL
             , &numelOv    , &nNodeOv
             , &nno1 
             , maxVizPart  , maxVizPartNo
             , &lNel       , &lNnode
             , &nRcvs      , &nSends
             , &nComNo     , maxNo  
             , ndfD        , ndfT
             , ndfF        , ndfFt
             , ndfComb     , nSp
             , ndm         , maxViz
             , &nVizPart   , &nVizPartNo
             , nPart );        
    }
/*...................................................................*/

/*...*/
    if(!mpiVar.myId)
    {
      free(elLG);
      free(noLG);
    } 
/*...................................................................*/ 
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
  
/*... alocao de variaveis que não precisao de comunicacao*/
  lNel     = mesh->numel;
  numelNov = mesh->numelNov;
  lNnode   = mesh->nnode;
  
/*... centroide */
  HccaAlloc(DOUBLE,m,mesh->elm.geom.cc ,lNel*ndm,"elCcP"    ,_AD_);
/*... face por elemento*/
  HccaAlloc(INT, m, mesh->elm.cellFace, lNel*maxViz,"cellfaceP", _AD_);
/*... volume da celula*/                           
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.volume
        ,lNel            ,"elVolP",_AD_);
/*... vetor que une o centroide ao ponto medio*/                           
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.xmcc
         ,lNel*maxViz*ndm       ,"elxmccP",_AD_);
/*... distancia entro o ponto medio da face e ponto de intercao 
      entre a linha entre os centroides e a face*/         
  HccaAlloc(DOUBLE               ,m       ,mesh->elm.geom.dcca 
         ,lNel*maxViz           ,"eldccaP",_AD_);

/*... zerando os variavies*/
  zero(mesh->elm.geom.cc      ,lNel*ndm       ,DOUBLEC);
  zero(mesh->elm.cellFace     ,lNel*maxViz
      ,DOUBLEC);
  zero(mesh->elm.geom.volume  ,lNel           ,DOUBLEC);
  zero(mesh->elm.geom.xmcc    ,lNel*ndm*maxViz,DOUBLEC);
  zero(mesh->elm.geom.dcca    ,lNel*maxViz    ,DOUBLEC);
/*...................................................................*/

/*... transporte e fluido*/
  if(mesh->ndfT[0] > 0 || mesh->ndfF > 0 || mesh->ndfFt > 0) {     
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
             ,numelNov*ndm*ndfD[0],"rCellUd1P"      ,_AD_);
/*...................................................................*/
    zero(mesh->node.uD1      ,lNnode*ndfD[0]       ,DOUBLEC);
    zero(mesh->node.gradUd1  ,lNnode*ndm*ndfD[0]   ,DOUBLEC);
    zero(mesh->elm.gradUd1   ,lNel*ndm*ndfD[0]     ,DOUBLEC);
    zero(mesh->elm.rCellUd1  ,numelNov*ndm*ndfD[0] ,DOUBLEC);
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
             ,numelNov*ndm*ndfT[0],"rCellUt1P"      ,_AD_);
/*...................................................................*/
    zero(mesh->node.uT1      ,lNnode*ndfT[0]       ,DOUBLEC);
    zero(mesh->node.gradUt1  ,lNnode*ndm*ndfT[0]   ,DOUBLEC);
    zero(mesh->elm.gradUt1   ,lNel*ndm*ndfT[0]     ,DOUBLEC);
    zero(mesh->elm.rCellUt1  ,numelNov*ndm*ndfT[0] ,DOUBLEC);
  }
/*...................................................................*/

/*... problema de transporte*/
  if(fTurb)
  {
/*... eddyViscosity*/
    HccaAlloc(DOUBLE,m,mesh->elm.eddyViscosity
              ,lNel    ,"eddyViscP"              ,_AD_);
    zero(mesh->elm.eddyViscosity,lNel,DOUBLEC);
/*...................................................................*/

/*... stressR*/
    if(fTurbStruct)
    {
      HccaAlloc(DOUBLE   ,m,mesh->elm.stressR  
                ,lNel*ntn  ,"stressRP"      ,_AD_);
      zero(mesh->elm.eddyViscosity,lNel*ntn,DOUBLEC);
    }
/*...................................................................*/

/*... paramentros de parede*/
     if(fWall)
     {
        HccaAlloc(DOUBLE, m, mesh->elm.wallParameters
                 , lNel*NWALLPAR   , "wallParmP", _AD_);
        zero(mesh->elm.wallParameters, lNel*NWALLPAR, DOUBLEC);
     }
/*...................................................................*/

/*... coeficientes dinamicos locais*/
     if(fDynamic)
     {
       HccaAlloc(DOUBLE, m, mesh->elm.cd
              , lNel*2   , "cDynamicP"     , _AD_);
       zero(mesh->elm.cd, lNel*2, DOUBLEC);
     } 

  }
/*...................................................................*/

/*... problema de escoamento*/
  if(ndfF > 0 || ndfFt > 0)
  {
/*... pres*/
    HccaAlloc(DOUBLE,m,mesh->node.pressure
              ,lNnode         ,"npresP"               ,_AD_);
/*... gradPres*/
    HccaAlloc(DOUBLE,m,mesh->node.gradPres  
              ,lNnode*ndm ,"nGradPresP"      ,_AD_);
    HccaAlloc(DOUBLE,m,mesh->elm.gradPres  
              ,lNel*ndm   ,"eGradPresP"     ,_AD_);
/*... gradVel*/
    HccaAlloc(DOUBLE,m,mesh->node.gradVel  
              ,lNnode*ndm*ndfVel ,"nGradVelP"      ,_AD_);
    HccaAlloc(DOUBLE,m,mesh->elm.gradVel  
              ,lNel*ndm*ndfVel   ,"eGradVelP"     ,_AD_);

    zero(mesh->node.pressure,lNnode,DOUBLEC);
    zero(mesh->node.gradPres,lNnode*ndm,DOUBLEC);
    zero(mesh->elm.gradPres ,lNel*ndm  ,DOUBLEC);
    zero(mesh->node.gradVel,lNnode*ndm*ndfVel,DOUBLEC);
    zero(mesh->elm.gradVel ,lNel*ndm*ndfVel  ,DOUBLEC);

/*...*/
    if(ndfFt > 0)
    {
/*... temp*/
      HccaAlloc(DOUBLE,m,mesh->node.temp,lNnode,"nTempP",_AD_);
/*... gradTemp*/
      HccaAlloc(DOUBLE,m,mesh->node.gradTemp
              ,lNnode*ndm, "nGradTempP", _AD_);
      HccaAlloc(DOUBLE,m,mesh->elm.gradTemp  
              ,lNel*ndm   ,"eGradTempP"     ,_AD_);
      
      zero(mesh->node.temp     ,lNnode    ,DOUBLEC);
      zero(mesh->node.gradTemp ,lNnode*ndm,DOUBLEC);
      zero(mesh->elm.gradTemp  ,lNel*ndm  ,DOUBLEC);
/*...................................................................*/

/*... energy*/
      HccaAlloc(DOUBLE,m,mesh->node.energy   
              ,lNnode         ,"nEngP"               ,_AD_);
/*... gradEnergy*/
      HccaAlloc(DOUBLE,m,mesh->node.gradEnergy
              ,lNnode*ndm ,"nGradEngP"      ,_AD_);
      HccaAlloc(DOUBLE,m,mesh->elm.gradEnergy  
              ,lNel*ndm   ,"eGradEngP"     ,_AD_);
/*... rCellEnergy*/
      HccaAlloc(DOUBLE     ,m   ,mesh->elm.rCellEnergy
               ,numelNov   ,"nCellEngP"      ,_AD_);
/*...*/
      zero(mesh->node.energy     ,lNnode    ,DOUBLEC);
      zero(mesh->node.gradEnergy ,lNnode*ndm,DOUBLEC);
      zero(mesh->elm.gradEnergy  ,lNel*ndm  ,DOUBLEC);
      zero(mesh->elm.rCellEnergy ,numelNov  ,DOUBLEC);

/*... eGradRho*/
       HccaAlloc(DOUBLE ,m            ,mesh->elm.gradRhoFluid
                ,lNel*ndm,"eGradFluidP",_AD_);
       zero(mesh->elm.gradRhoFluid, lNel*ndm, DOUBLEC);

/*... nRhoFluid*/
       HccaAlloc(DOUBLE, m, mesh->node.rhoFluid,lNnode, "nRhoFluidP", _AD_);
       zero(mesh->node.rhoFluid, lNnode, DOUBLEC);

    }
/*...................................................................*/

/*... rCellVel*/
    size = lNel*ndm*ndfVel;
    HccaAlloc(DOUBLE,m,mesh->elm.rCellVel  
             ,size       ,"rCellVelP",_AD_);
    zero(mesh->elm.rCellVel  ,size,DOUBLEC);
/*...................................................................*/

/*... rCellPres*/
    HccaAlloc(DOUBLE,m,mesh->elm.rCellPres  
             ,numelNov,"rCellPresP"  ,_AD_);
    zero(mesh->elm.rCellPres  ,numelNov         ,DOUBLEC);
/*...................................................................*/

/*...*/
    if(tModel->fOneEq)
    {
      HccaAlloc(DOUBLE,m,mesh->elm.rCellKturb  
               ,numelNov          ,"rCellKturbP"  ,_AD_);
      zero(mesh->elm.rCellKturb,numelNov,DOUBLEC);
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  if(ndfComb > 0)
  {
/*... temp*/
    HccaAlloc(DOUBLE,m,mesh->node.zComb
             ,lNnode*ndfComb,"nZcombP",_AD_);
    zero(mesh->node.zComb     ,lNnode*ndfComb    ,DOUBLEC);
/*... gradTemp*/
    HccaAlloc(DOUBLE,m,mesh->node.gradZcomb
             ,lNnode*ndm*ndfComb, "nGradZcombP", _AD_);
    zero(mesh->node.gradZcomb ,lNnode*ndm*ndfComb,DOUBLEC);
    HccaAlloc(DOUBLE,m,mesh->elm.gradZcomb  
             ,lNel*ndm*ndfComb  ,"eGradZcombP"     ,_AD_);
    zero(mesh->elm.gradZcomb  ,lNel*ndm*ndfComb  ,DOUBLEC);
/*...................................................................*/

/*... wK*/
    HccaAlloc(DOUBLE,m, mesh->elm.wk
            ,lNel*nSp ,"wkP",_AD_);
    zero(mesh->elm.wk, lNel* nSp, DOUBLEC);
/*... HeatRe*/
    HccaAlloc(DOUBLE, m, mesh->elm.wT, lNel, "wT",_AD_);
    zero(mesh->elm.wT, lNel, DOUBLEC);

/*... enthalpyK*/
    HccaAlloc(DOUBLE, m, mesh->elm.enthalpyk
            , lNel*nSp, "enthalpykP", _AD_);
    zero(mesh->elm.enthalpyk, lNel*nSp, DOUBLEC);

/*... eGradY*/
    HccaAlloc(DOUBLE, m, mesh->elm.gradY
            , lNel*ndm*nSp, "eGradYP", _AD_);
    zero(mesh->elm.gradY, lNel*ndm*nSp, DOUBLEC);

/*... timeReactor*/
    HccaAlloc(DOUBLE, m, mesh->elm.tReactor
            , lNel*N_TERMS_REACTOR, "tReactorP", _AD_);
    zero(mesh->elm.tReactor, lNel*N_TERMS_REACTOR, DOUBLEC);

/*... rCellComb*/
    HccaAlloc(DOUBLE, m, mesh->elm.rCellComb
             , numelNov*ndfComb, "rCellCombP", _AD_);
    zero(mesh->elm.rCellComb, numelNov*ndfComb, DOUBLEC);
/*...................................................................*/

  }
/*...................................................................*/

/*...*/
  maxNdf = 0;
  maxNdf = max(ndm     ,maxNdf);
  maxNdf = max(ndfD[0] ,maxNdf);
  maxNdf = max(ndfD[1] ,maxNdf);
  maxNdf = max(ndfD[2] ,maxNdf);
  maxNdf = max(ndfT[0] ,maxNdf);
  maxNdf = max(ndfT[1] ,maxNdf);
  maxNdf = max(ndfT[2] ,maxNdf);
  maxNdf = max(ndfF    ,maxNdf);
  maxNdf = max(ndfFt   ,maxNdf);
  maxNdf = max(nSp     ,maxNdf);
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

void writeMeshPart(Mesh *mesh,Combustion *cModel)
{
  short j,nD;
  INT i;

/*... densityFluid*/
  fprintf(fileLogDebug,"Densisty\n");
  for(i=0;i<mesh->numel;i++)
  {
    fprintf(fileLogDebug,"%d ",i);
    fprintf(fileLogDebug," %e %e %e\n"
                               ,mesh->elm.densityFluid.t[i]
                               ,mesh->elm.densityFluid.t0[i]
                               ,mesh->elm.densityFluid.t00[i]);
    fprintf(fileLogDebug,"\n");
  }
/*......................................................................*/

/*... specific heat*/
  fprintf(fileLogDebug,"Specific heat\n");
  for(i=0;i<mesh->numel;i++)
  {
    fprintf(fileLogDebug,"%d ",i);
    fprintf(fileLogDebug," %e %e %e\n"
                               ,mesh->elm.specificHeat.t[i]
                               ,mesh->elm.specificHeat.t0[i]
                               ,mesh->elm.specificHeat.t00[i]);
    fprintf(fileLogDebug,"\n");
  }
/*......................................................................*/

/*... viscosidade dinamica*/
  fprintf(fileLogDebug,"Viscosidade dinamica\n");
  for(i=0;i<mesh->numel;i++)
    fprintf(fileLogDebug,"%d %e\n",i,mesh->elm.dViscosity[i]);
  fprintf(fileLogDebug,"\n");
/*......................................................................*/

/*... condutivade termica*/
  fprintf(fileLogDebug,"Condutividade termica\n");
  for(i=0;i<mesh->numel;i++)
    fprintf(fileLogDebug,"%d %e\n",i,mesh->elm.tConductivity[i]);
  fprintf(fileLogDebug,"\n");
/*......................................................................*/

/*... diffusion species*/
  fprintf(fileLogDebug,"Difusao massica\n");
  for(i=0;i<mesh->numel;i++)
  {
    fprintf(fileLogDebug,"%d ",i);
    for(j=0;j<nD;j++)
      fprintf(fileLogDebug," %e ",MAT2D(i,j,mesh->elm.cDiffComb,nD));
    fprintf(fileLogDebug,"\n");
  }
  fprintf(fileLogDebug,"\n");
/*......................................................................*/

/*... yFrac0*/
  nD = cModel->nOfSpecies;
  fprintf(fileLogDebug,"yFrac0\n");
  for(i=0;i<mesh->numel;i++)
  {
    fprintf(fileLogDebug,"%d ",i);
    for(j=0;j<nD;j++)
      fprintf(fileLogDebug," %e ",MAT2D(i,j,mesh->elm.yFrac0,nD));
    fprintf(fileLogDebug,"\n");
  }
  fprintf(fileLogDebug,"\n");
/*......................................................................*/

/*... yFrac*/
  nD = cModel->nOfSpecies;
  fprintf(fileLogDebug,"yFrac\n");
  for(i=0;i<mesh->numel;i++)
  {
    fprintf(fileLogDebug,"%d ",i);
    for(j=0;j<nD;j++)
      fprintf(fileLogDebug," %e ",MAT2D(i,j,mesh->elm.yFrac,nD));
    fprintf(fileLogDebug,"\n");
  }
  fprintf(fileLogDebug,"\n");
/*......................................................................*/

/*... z0*/
  nD = cModel->nComb;
  fprintf(fileLogDebug,"z0\n");
  for(i=0;i<mesh->numel;i++)
  {
    fprintf(fileLogDebug,"%d ",i);
    for(j=0;j<nD;j++)
      fprintf(fileLogDebug," %e ",MAT2D(i,j,mesh->elm.zComb0,nD));
    fprintf(fileLogDebug,"\n");
  }
  fprintf(fileLogDebug,"\n");
/*......................................................................*/

/*... z*/
  nD = cModel->nComb;
  fprintf(fileLogDebug,"z\n");
  for(i=0;i<mesh->numel;i++)
  {
    fprintf(fileLogDebug,"%d ",i);
    for(j=0;j<nD;j++)
      fprintf(fileLogDebug," %e ",MAT2D(i,j,mesh->elm.zComb,nD));
    fprintf(fileLogDebug,"\n");
  }
  fprintf(fileLogDebug,"\n");
/*......................................................................*/

/*... temp0*/
  fprintf(fileLogDebug,"temp0\n");
  for(i=0;i<mesh->numel;i++)
    fprintf(fileLogDebug,"%d %e\n",i,mesh->elm.temp0[i]);
  fprintf(fileLogDebug,"\n");
/*......................................................................*/

/*... temp*/
  fprintf(fileLogDebug,"temp\n");
  for(i=0;i<mesh->numel;i++)
    fprintf(fileLogDebug,"%d %e\n",i,mesh->elm.temp[i]);
  fprintf(fileLogDebug,"\n");
/*......................................................................*/

}
/*********************************************************************/