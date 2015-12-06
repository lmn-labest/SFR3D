#ifndef _PARALLELMPI_H_
  #define _PARALLELMPI_H_
/*...*/
  #include<mpi.h>
  #include<stdio.h>
  #include<stdlib.h>
/*...*/
  #include<Erro.h>
  #include<Define.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<PartMesh.h>
  #include<ElIncid.h>
  #include<HccaTime.h>
  

/*...*/
  void mpiStart(int *argc,char **argv);
  void mpiStop(void);
  void mpiWait(void);
  void comunicateMesh(Memoria *m
                   ,Mesh *mesh0     ,Mesh *mesh
                   ,PartMesh *pMesh
                   ,Loads *loadsD1  ,Loads *loadsT1);
/*... equacoes nas interfaces*/  
  void getBuffer(DOUBLE *restrict x    ,DOUBLE *restrict xb
              ,INT *restrict fMap    ,INT const nRcvs);
  
  void makeBuffer(DOUBLE *restrict x    ,DOUBLE *restrict xb
               ,INT *restrict fMap    ,INT const nSends);
  
  void comunicateNeq(Interface *iNeq,DOUBLE *restrict x);
/*...................................................................*/

/*... celulas nas interfaces*/  
  void makeBufferCel(DOUBLE *restrict x    ,DOUBLE *restrict xb
                  ,INT *restrict fMap    ,INT const nSends
                  ,short const ndf       ,short const ndm);
  
  void getBufferCel(DOUBLE *restrict x    ,DOUBLE *restrict xb
              ,INT *restrict fMap       ,INT const nRcvs
              ,short const ndf          ,short const ndm);
  
  void comunicateCel(Interface *iCel,DOUBLE *restrict x
                    ,short const ndf,short const ndm);
/*...................................................................*/


/*... nos nas interfaces*/  
  void dGetBufferNod(DOUBLE *restrict x    ,DOUBLE *restrict xb
                    ,INT *restrict fMap       ,INT const nRcvs
                    ,short const ndf          ,short const ndm);

  void dMakeBufferNod(DOUBLE *restrict x    ,DOUBLE *restrict xb
                  ,INT *restrict fMap    ,INT const nSends
                  ,short const ndf       ,short const ndm);


  void dComunicateNod(InterfaceNo *iNo ,DOUBLE *restrict x
                     ,short const ndf,short const ndm);


  void iGetBufferNod(INT *restrict x    ,INT *restrict xb
                    ,INT *restrict fMap ,INT const nRcvs
                    ,short const ndf    ,short const ndm);

  void iMakeBufferNod(INT *restrict x    ,INT *restrict xb
                     ,INT *restrict fMap ,INT const nSends
                     ,short const ndf    ,short const ndm);


  void iComunicateNod(InterfaceNo *iNo ,INT *restrict x
                     ,short const ndf,short const ndm);
/*...................................................................*/

/*... globaliza os valores nodas no processo master(0)*/
  void dGlobalNode(Memoria *m         ,PartMesh *pMesh
                  ,DOUBLE *restrict uG,DOUBLE *restrict uL
                  ,short const ndf1   ,short const ndf2);
/*...................................................................*/

/*... globaliza os valores das celulas no processo master(0)*/
  void dGlobalCel(Memoria *m         ,PartMesh *pMesh
               ,DOUBLE *restrict uG,DOUBLE *restrict uL
               ,INT const numelNov 
               ,short const ndf1   ,short const ndf2);
/*...................................................................*/

  void globalMeshQuality(MeshQuality *mQl,MeshQuality *mQl0);
/*...................................................................*/

/*...*/
  void comunicate2(short *m0faceR     ,short *faceR
                ,short *m0faceL     ,short *faceL
                ,DOUBLE *m0u0       ,DOUBLE *u0 
                ,DOUBLE *m0u        ,DOUBLE *u  
                ,DOUBLE *m0density  ,DOUBLE *density  
                ,INT const lNel     ,INT *elLG
                ,short const maxViz ,short const ndf
                ,short const npart  ,short const iCod);
/*...................................................................*/

/*... do Mpi*/
  typedef struct{
    unsigned short nPrcs;
    unsigned short myId;
    MPI_Comm comm;
    char errBuffer[MPI_MAX_ERROR_STRING];
    int lString;
    int ierr;
    MPI_Request sendRequest[MAX_MPI_PROCESS];
    MPI_Request recvRequest[MAX_MPI_PROCESS];
    MPI_Status  status[MAX_MPI_PROCESS];

  }Mpi;

  Mpi mpiVar;
/*...................................................................*/
 
#endif /*... _PARALLELMPI_H*/
