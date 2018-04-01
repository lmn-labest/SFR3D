#ifndef _MESH_
  #define _MESH_
//  #include<Adjcency.h>
  #include<Reord.h>
  #include<HccaStdBool.h>
  #include<Define.h>

/*... Material*/
  typedef struct{
    short  *type;      /*tipo de calculo da celula*/ 
    double *prop;      /*valores numericos da propriedade*/    
  } Material;
/*...................................................................*/

/*...*/
  typedef struct{
    short *nViz;
    INT  *nelcon;
  }Adjacency;
/*...................................................................*/  

/*...*/
  typedef struct{
    bool fDensity;
    bool fSpecificHeat;
    bool fDynamicViscosity;
    bool fThermalCondutivty;
  }PropVar;
/*...................................................................*/
 
/*...*/
  typedef struct{
    bool fRes;
    INT maxIt;
    DOUBLE ck,ce,sk,tol;
  }EqK;
/*...................................................................*/  


/*...*/
  typedef struct{  
    bool fTurb;
    bool fWall;
    bool fDynamic;
    short wallType;
    short type;          /* 0 - LES*/
    short typeLes;       /* 0 - funcional
                            1 - estrutural
                            2 - misto*/ 
    short typeMixed[2];  /*[0] - estrutural - [1] funcional*/
    short typeDynamic;    /* 1 - um paremetro local 
                            2 - um parametro global padrao
                            3 - um parametro global modicado
                            4 - 2 paramento local*/

    DOUBLE cs,cf,c; /*constante    
                           1 -        
                           */   
    DOUBLE PrandltTwall;   /*Prandtl turbulento */
    DOUBLE PrandltTsgs;    /*Prandtl de sub-grid */

    EqK eK;
    
  }Turbulence; 
/*...................................................................*/

/*...*/
  typedef struct{
    bool fPresWork;
    bool fDissipation;
    bool fRes;
    bool fTemperature;
    bool fKelvin;
  }EnergyModel;
/*...................................................................*/  

/*...*/
  typedef struct{
    bool RhsDensity,LhsDensity;
  }MassEqModel;
/*...................................................................*/

/*...*/
  typedef struct{    
    bool fRes,fRhieChowInt,fViscosity,fDiv; 
    short iCodBuoyant;
  }MomentumModel;
/*...................................................................*/

/*...*/
  typedef struct{
    bool flag,fDynamic;
    short iCod,type;
    INT timeStep;
    DOUBLE total,dt[3],dtInicial,t;   
  }Temporal;
/*...................................................................*/  

/*...*/
  typedef struct{
    bool fRes;
    DOUBLE cEq[3];
  }kModel;
/*...................................................................*/  

/*...*/
/*...................................................................*/  

/*...*/
  typedef struct{
    short  iCod1;
    short  iCod2;
    DOUBLE par[NPADV];
  }Advection;
/*...................................................................*/  

/*...*/
  typedef struct{
    short  iCod;
  }Diffusion;
/*...................................................................*/  

/*...*/
  typedef struct{
    DOUBLE *cc;    /*centroide da celulas*/
    DOUBLE *ksi;   /*vetor que une os centroides da celulas*/
    DOUBLE *mksi;  /*modulo do vetor que une os centroides da celulas*/
    DOUBLE *eta;   /*vetor paralelo a aresta*/
    DOUBLE *fArea; /*area da face compartilhada*/
    DOUBLE *normal;/*vetor normal a face*/       
    DOUBLE *volume;/*volume da celula*/
    DOUBLE *xm;    /*ponto medio da face*/
    DOUBLE *xmcc;  /*vetor entre o centroide a ponto médio da aresta*/
    DOUBLE *dcca;  /*menor distancia entre o centroide a aresta*/
    DOUBLE *vSkew; /*vetor vSkew*/                               
    DOUBLE *mvSkew;/*modulo do vSkew*/                               
  }Geom;
/*...................................................................*/

/*...*/
  typedef struct{
    DOUBLE *vel;    /*Velocidade na face*/                             
  }Face;
/*...................................................................*/    

/*... interpol*/
  typedef struct{
    unsigned short np;                       /*numero de particoes*/  
    DOUBLE *x,*y;                            /*interpolacao*/
  }Interpol;
  Interpol iPol[MAXINTERPOL];
/*...................................................................*/  

/*... loads*/
  typedef struct{
    short nTypeVar;               /*0 - constante
                                    1 - funcao parabolica*/  
    short type;                     /*tipo*/
    short np;                       /*numero de particoes*/  
    DOUBLE par[MAXLOADPARAMETER];
    Interpol *intPol;               /*interpolacao*/
  }Loads;
  Loads  loadsD1[MAXLOADD1]         /*tipo de cargas (difusao pura)*/
        ,loadsT1[MAXLOADT1]         /*tipo de cargas (difusao-transporte)*/
        ,loadsVel[MAXLOADFLUID]     /*tipo de cargas (fluid-Vel)*/
        ,loadsPres[MAXLOADFLUID]    /*tipo de cargas (fluid-Pres)*/
        ,loadsPresC[MAXLOADFLUID]   /*tipo de cargas (fluid-Pres-correcao)*/
        ,loadsEnergy[MAXLOADFLUID]  /*tipo de cargas (fluid-energia)*/
        ,loadsTemp[MAXLOADFLUID]    /*tipo de cargas (fluid-temperatura)*/
        ,loadsKturb[MAXLOADFLUID];  /*tipo de cargas (fluid-turbulencia)*/
/*...................................................................*/

/*...*/
  typedef struct{
    int maxIt;
    DOUBLE tol;                               
  }NonLinear;
/*...................................................................*/  

/*...*/
  typedef struct{
    bool fInit;
    bool fMedia;
    bool fVel;
    int startSample,endSample;
    DOUBLE *mVel,*sVel;  
    DOUBLE t0;                             
  }Mean;
/*...................................................................*/  

/*... Elementos*/
  typedef struct{
    short  *mat;       /*materiais do elementos*/   
    short  *nen;       /*numero de no por elemento*/
    short  *geomType;  /*tipo geometrio do elemento*/
    short  *rNum;      /*renumeracao dos elementos*/ 
/*... */               
    short *faceRvel;       /*condicao  contorno na face (fluido)*/
    short *faceLoadVel;    /*tipo de carga contorno na face (fluido)*/
    short *faceRpres;      /*condicao  contorno na face (fluido)*/
    short *faceLoadPres;   /*tipo de carga contorno na face (fluido)*/
    short *faceRenergy;    /*condicao  contorno na face (fluido)*/
    short *faceLoadEnergy; /*tipo de carga contorno na face (fluido)*/
/*...*/
    short *faceReKturb;     /*condicao  contorno na face (fluido)*/
    short *faceLoadKturb;  /*tipo de carga contorno na face (fluido)*/
/*... */               
    short *faceRt1;    /*condicao  contorno na face (transporte)*/
    short *faceLoadT1; /*tipo de carga contorno na face (transporte)*/
/*...*/
    short *faceRd1;    /*condicao  contorno na face (difusa pura)*/
    short *faceLoadD1; /*tipo de carga contorno na face (difusa pura)*/
/*...................................................................*/
    INT    *node;       /*conectividades*/
    Geom   geom;  
/*...*/
    DOUBLE *energy;       /*energia*/
    DOUBLE *energy0;      /*energia*/
    DOUBLE *pressure;     /*pressao (n+1)*/
    DOUBLE *pressure0;    /*pressao (n)*/
    DOUBLE *temp;         /*temperatura (n+1)*/
    DOUBLE *temp0;         /*temperatura (n)*/
    DOUBLE *vel;          /*velocidade do fluido*/
    DOUBLE *vel0;         /*velocidade do fluido*/
    DOUBLE *densityFluid; /*massa especifica do fluido*/
    DOUBLE *specificHeat; /*calor especifico*/
    DOUBLE *dViscosity;   /*viscosidade dinamica*/
    DOUBLE *tConductivity;/*condutividae termica*/
/*...*/
    DOUBLE *gradVel;    /*gradiente do campo de velocidade*/    
    DOUBLE *gradPres;   /*gradiente do campo de pressao*/    
    DOUBLE *gradEnergy; /*gradiente do campo de energia*/
    DOUBLE *gradTemp;   /*gradiente da temperatura*/
/*...*/
    DOUBLE *rCellVel;   /*residuo da celula*/
    DOUBLE *rCellPres;  /*residuo da celula*/
    DOUBLE *rCellEnergy;/*residuo da celula*/
/*... turbulencia*/
    DOUBLE *kTurb;
    DOUBLE *kTurb0;
    DOUBLE *stressR;
    DOUBLE *eddyViscosity; 
    DOUBLE *wallParameters;
    DOUBLE *cd; 
    DOUBLE *gradKturb;    /*gradiente da difusao pura uT1*/
    DOUBLE *rCellKturb;   /*residuo da celula*/         
  
/*...*/
    DOUBLE *densityUd1; /*massa especifica do material uD1*/
    DOUBLE *uD1 ;       /*difusao pura uD1*/
    DOUBLE *u0D1;       /*difusao pura uD1*/
    DOUBLE *gradUd1;    /*gradiente da difusao pura uD1*/
    DOUBLE *rCellUd1;   /*residuo da celula*/
/*...*/
    DOUBLE *densityUt1; /*massa especifica do material uT1*/
    DOUBLE *uT1 ;       /*difusao pura uT1*/
    DOUBLE *u0T1;       /*difusao pura uT1*/
    DOUBLE *gradUt1;    /*gradiente da difusao pura uT1*/
    DOUBLE *rCellUt1;   /*residuo da celula*/
/*...*/
    DOUBLE *leastSquare; /*matriz de aproxima leastSquare*/
    DOUBLE *leastSquareR;/*fatoracao QR*/
/*...*/
    Material material;
    Adjacency adj;
  }Elmt;
/*...................................................................*/
  

/*... nos*/
  typedef struct{
    INT    *nno; 
    DOUBLE *x;         /*coordenadas*/
    DOUBLE *vel;       /*velocidades*/
    DOUBLE *energy;    /*energia*/
    DOUBLE *kTurb ;    /*energia cinetica turbulenta*/
    DOUBLE *pressure;  /*pressao*/
    DOUBLE *temp;      /*temperatura*/
    DOUBLE *uD1;       /*difusao pura uD1*/
    DOUBLE *uT1;       /*difusao pura uT1*/
    DOUBLE *gradTemp;  /*gradiente da temperatura*/
    DOUBLE *gradUd1 ;  /*gradiente da difusao pura uD1*/
    DOUBLE *gradUt1 ;  /*gradiente da difusao pura uT1*/
    DOUBLE *gradVel ;  /*gradiente da da velocidad ( Matriz Jacobiana) 
                            | du1dx1 du1dx2 du1dx3 |   
                            | du2dx1 du2dx2 du2dx3 |   
                            | du3dx1 du3dx2 du3dx3 |*/                          
    DOUBLE *gradPres;     /*gradiente da Pressao*/
    DOUBLE *gradEnergy;   /*gradiente da Energia*/
  }Node;
/*...................................................................*/

/*... meshQuality*/
  typedef struct{
    DOUBLE volume;
    DOUBLE nonOrthMed,nonOrthMax;
    DOUBLE skewMed,skewMax;
    DOUBLE aspectRaMax,aspectRaMin;
  }MeshQuality;
/*...................................................................*/

/*... Simple*/
  typedef struct{
    bool   sPressure;
    unsigned short faceInterpolVel;   /*tecnica de interpolacao
                                        das velocidades nas faces para
                                        evitar o problema checkboard
                                       */
    short  type;
 
    unsigned short kZeroVel;          /*iteracao com o qual o residuo
                                        e normalizado*/
    unsigned short kZeroPres;         /*iteracao com o qual o residuo
                                        e normalizado*/
    unsigned short kZeroEnergy;       /*iteracao com o qual o residuo
                                     e normalizado*/
    int    pSimple;
    int    maxIt;
    int nNonOrth;
    DOUBLE alphaPres,alphaVel,alphaEnergy
          ,alphaDensity;                  /*under-relaxation*/
    DOUBLE *ePresC,*nPresC  ,*eGradPresC;/*Pressao de correcao*/
    DOUBLE *ePresC1;                     /*Pressao de correcao 1*/
    DOUBLE *d;
    DOUBLE tolPres,tolVel[3],tolEnergy;
  }Simple;
/*...................................................................*/

/*... Prime*/
  typedef struct {
    bool   sPressure;
    unsigned short faceInterpolVel;   /*tecnica de interpolacao
                                      das velocidades nas faces para
                                      evitar o problema checkboard
                                     */
    short  type;
    unsigned short kZeroVel;          /*iteracao com o qual o residuo
                                      e normalizado*/
    unsigned short kZeroPres;         /*iteracao com o qual o residuo
                                      e normalizado*/
    int nNonOrth;
    int    pPrime;          /*iteracao que serao imprisas na tela*/ 
    int    maxIt;
    DOUBLE alphaPres, alphaVel;           /*under-relaxation*/
    DOUBLE *ePresC, *nPresC, *eGradPresC;/*Pressao de correcao*/
    DOUBLE *ePresC1;                     /*Pressao de correcao 1*/
    DOUBLE *d,*velUp,*aD,*bTemporal;
    DOUBLE tolPres, tolVel;
  }Prime;
/*...................................................................*/

 
/*...*/
  typedef struct{
    INT *nincid;
    INT *incid;
  }Pnode;
 
/*...................................................................*/
 
/*... Malha*/
  typedef struct{
    bool fOpen;   /*dominio aberto*/    
    short ntn;    /*dimenso do tensor simetrico*/
    short ndm;     /*dimensao*/
    short ndfF;    /*fluido*/    
    short ndfFt;   /*fluido termo ativado*/  
    short ndfT[MAX_TRANS_EQ];   /*graus de liberdade 
                                  para o problema de transporte*/
    short ndfD[MAX_DIF_EQ];   /*graus de liberdade 
                                  para o problema de difusa pura*/
    short numat;    /*numero maximo de materias no dominio*/
    short maxNo;    /*numero maximo de nos por elemento*/
    short maxViz;   /*numero maximo de vizinhos que um elemento possui*/
    INT nnode;      /*numero de nos*/
    INT numel;      /*numero de elementos*/
    INT numelNov;   /*numero de elementos sobrepostos*/
    INT nnodeOv;    /*numero de nos em elementos sobrepostos*/
    INT nnodeNov;   /*numero de nos em elementos nao sobrepostos*/
    DOUBLE mass[3]; /* mass inicial do sistema 
                       massa atual calculo incremental
                       massa atual calculo direto*/
    DOUBLE xRef[3]; /*... ponto de referencia*/
/*...*/    
    Elmt elm;     
    Node node;
    Pnode noIncid;
/*...*/
    MeshQuality mQuality;
  }Mesh;
/*...................................................................*/

/*...*/
  typedef struct{
    short rcGrad; /*tipo de rescontrucao de gradiente*/                     
/*...*/
    Temporal ddt;
/*...*/
    NonLinear nlTemp;
    NonLinear nlD1;
    NonLinear nlT1;
/*... equacao de transporte*/
    Advection advT1;
    Diffusion diffT1;
/*... equacao de velocidade*/
    Advection advVel;
    Diffusion diffVel;
/*... equacao de energia*/
    Advection advEnergy;
    Diffusion diffEnergy;
/*... equacao de energia*/
    Advection advKturb;
    Diffusion diffKturb;
/*... equacao de pressao*/
    Diffusion diffPres;
  }Scheme;
/*...................................................................*/

 typedef struct{
    INT  nRcvs;   /*numero de valores a serem recebidos*/
    INT  nSends;  /*numero de valores a serem enviados */
    INT *iaSends; /*ponteiros para o arranjo no buffer de envio*/
    INT *iaRcvs;  /*ponteiros para o arranjo buffer de recebimento*/              
    INT *fMap;    /*numeracao do buffer de comunicacao*/    
    DOUBLE *xb;   /*valores no buffer de comunicacao*/
    unsigned short nVizPart;
    short *vizPart;
    
  }Interface;
 
  typedef struct{
    unsigned short nVizPart;
    short *vizPart;
    INT  nCom;     /*numero de valores a serem comunicados*/
    INT  *iaComNo; /*ponteiros para o arranjo buffer de comunicacao*/              
    INT   *fMap;   /*numeracao do buffer de comunicacao*/    
    INT     *xi;   /*valores no buffer de comunicacao(INT)*/
    DOUBLE  *xb;   /*valores no buffer de comunicacao(DOUBLE)*/    
  }InterfaceNo;

/* variaveis do preparticionador*/
  typedef struct{
    bool fPartMesh;
    bool fPrintMesh;
    bool fPrintMeshPart;
    INT nno1;      /*numero de nos apenas da particao local*/
    INT nnG ;      /*numero de nos total da malha global*/
    INT elG ;      /*numero de elementos total da malha global*/
    INT *elLG;     /*mapa de local-> global de celulas*/
    INT *noLG;     /*mapa de local-> global de nos*/
    INT *np;       /*master apenas*/
    INT *ep;       /*master apenas*/
    Interface   iEl;
    InterfaceNo iNo;
  }PartMesh;
/*...................................................................*/
#endif/*_MESH_*/
