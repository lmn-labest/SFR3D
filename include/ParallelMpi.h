#ifndef _PARALLELMPI_H_
  #define _PARALLELMPI_H_
/*...*/
  #ifdef _MPI_
    #include<mpi.h>
  #endif
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
  void comunicateMesh(Memoria *m    ,Combustion *cModel
                   ,Turbulence *tModel
                   ,Mesh *mesh0     ,Mesh *mesh
                   ,PartMesh *pMesh
                   ,Loads *ldD1     ,Loads *ldT1
                   ,Loads *ldVel    ,Loads *ldPres
                   ,Loads *ldPresC  ,Loads *ldEnergy
                   ,Loads *ldTemp   ,Loads *ldKturb
                   ,Loads *ldZcomb);
/*... equacoes nas interfaces*/
  void getBuffer(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
              ,INT *RESTRICT fMap    ,INT const nRcvs);

  void makeBuffer(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
               ,INT *RESTRICT fMap    ,INT const nSends);

  void comunicateNeq(Interface *iNeq,DOUBLE *RESTRICT x);
/*...................................................................*/

/*... celulas nas interfaces*/
  void makeBufferCel(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
                  ,INT *RESTRICT fMap    ,INT const nSends
                  ,short const ndf       ,short const ndm);

  void getBufferCel(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
              ,INT *RESTRICT fMap       ,INT const nRcvs
              ,short const ndf          ,short const ndm);

  void comunicateCel(Interface *iCel,DOUBLE *RESTRICT x
                    ,short const ndf,short const ndm);
/*...................................................................*/


/*... nos nas interfaces*/
  void dGetBufferNod(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
                    ,INT *RESTRICT fMap       ,INT const nRcvs
                    ,short const ndf          ,short const ndm);

  void dMakeBufferNod(DOUBLE *RESTRICT x    ,DOUBLE *RESTRICT xb
                  ,INT *RESTRICT fMap    ,INT const nSends
                  ,short const ndf       ,short const ndm);


  void dComunicateNod(InterfaceNo *iNo ,DOUBLE *RESTRICT x
                     ,short const ndf,short const ndm);


  void iGetBufferNod(INT *RESTRICT x    ,INT *RESTRICT xb
                    ,INT *RESTRICT fMap ,INT const nRcvs
                    ,short const ndf    ,short const ndm);

  void iMakeBufferNod(INT *RESTRICT x    ,INT *RESTRICT xb
                     ,INT *RESTRICT fMap ,INT const nSends
                     ,short const ndf    ,short const ndm);


  void iComunicateNod(InterfaceNo *iNo ,INT *RESTRICT x
                     ,short const ndf,short const ndm);
/*...................................................................*/

/*... globaliza os valores nodas no processo master(0)*/
  void dGlobalNode(Memoria *m         ,PartMesh *pMesh
                  ,DOUBLE *RESTRICT uG,DOUBLE *RESTRICT uL
                  ,short const ndf1   ,short const ndf2);
/*...................................................................*/

/*... globaliza os valores das celulas no processo master(0)*/
  void dGlobalCel(Memoria *m         ,PartMesh *pMesh
               ,DOUBLE *RESTRICT uG,DOUBLE *RESTRICT uL
               ,INT const numelNov
               ,short const ndf1   ,short const ndf2);
/*...................................................................*/

/*...*/
  void globalMeshQuality(MeshQuality *mQl,MeshQuality *mQl0);
/*...................................................................*/

/*...*/
  void writeMeshPart(Mesh *mesh,Combustion *cModel);
/*...................................................................*/

/*... do Mpi*/
  typedef struct{
    unsigned short nPrcs;
    unsigned short myId;
	int lString;
	int ierr;
#ifdef _MPI_
    MPI_Comm comm;
	char errBuffer[MPI_MAX_ERROR_STRING];
    MPI_Request sendRequest[MAX_MPI_PROCESS];
    MPI_Request recvRequest[MAX_MPI_PROCESS];
    MPI_Status  status[MAX_MPI_PROCESS];
#endif

  }Mpi;

  extern Mpi mpiVar;
/*...................................................................*/

#endif /*_PARALLELMPI_H_*/
