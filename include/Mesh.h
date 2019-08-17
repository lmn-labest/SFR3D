#ifndef _MESH_
  #define _MESH_

  #include<Reord.h>
  #include<HccaStdBool.h>
  #include<Define.h>
  #include<StructDefine.h>

/*...*/
  typedef struct{
    DOUBLE *cc;    /*centroide da celulas*/
    DOUBLE *volume;/*volume da celula*/
    DOUBLE *xmcc;  /*vetor entre o centroide a ponto médio da aresta*/
    DOUBLE *dcca;  /*menor distancia entre o centroide a aresta*/
  }Geom;
/*...................................................................*/

/*...*/
  typedef struct{
    INT *node;
    INT *owner;
    DOUBLE *ksi;   /*vetor que une os centroides da celulas*/
    DOUBLE *mksi;  /*modulo do vetor que une os centroides da celulas*/
    DOUBLE *eta;   /*vetor paralelo a aresta*/
    DOUBLE *area; /*area da face compartilhada*/
    DOUBLE *normal;/*vetor normal a face*/
    DOUBLE *xm;    /*ponto medio da face*/
    DOUBLE *vSkew; /*vetor vSkew*/
    DOUBLE *mvSkew;/*modulo do vSkew*/
  }Face;
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

/*...*/
    short *faceResZcomb;
    short *faceLoadZcomb;
/*...................................................................*/

    INT    *node;       /*conectividades*/
    INT    *cellFace;   /*faces que gera as celulas*/
    Geom   geom;  
/*...*/
    DOUBLE *enthalpyk;    /*Entapia por especie*/
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
/*... combustao*/
    DOUBLE *zComb0;
    DOUBLE *zComb;
    DOUBLE *yFrac0;
    DOUBLE *yFrac;
    DOUBLE *gradZcomb;
    DOUBLE *gradY;
    DOUBLE *rCellComb; 
    DOUBLE *cDiffComb;
    DOUBLE *wk;         /*taxa de consumo das especies kg/(m3 s)*/
    DOUBLE *rateHeatReComb;
    DOUBLE *tReactor;
/*...*/
    DOUBLE *densityUd1; /*massa especifica do material uD1*/
    DOUBLE *cDiffD1;    /*coefieiente de difusao D1*/
    DOUBLE *uD1 ;       /*difusao pura uD1*/
    DOUBLE *u0D1;       /*difusao pura uD1*/
    DOUBLE *gradUd1;    /*gradiente da difusao pura uD1*/
    DOUBLE *rCellUd1;   /*residuo da celula*/
/*...*/
    DOUBLE *densityUt1; /*massa especifica do material uT1*/
    DOUBLE *cDiffT1;    /*coefieiente de difusao T1*/
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
/*...*/
    Combustion comb;    /*Combustao*/
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
    
    DOUBLE *zComb;
    DOUBLE *gradZcomb;
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
    INT nFaces;     /*numero de faces totais*/
    INT numelNov;   /*numero de elementos sobrepostos*/
    INT nnodeOv;    /*numero de nos em elementos sobrepostos*/
    INT nnodeNov;   /*numero de nos em elementos nao sobrepostos*/
    DOUBLE mass[3]; /* mass inicial do sistema 
                       massa atual calculo incremental
                       massa atual calculo direto*/
    DOUBLE massInOut[2]; /* mass que adentra o dominio 
                            mass que sai do dominio*/ 
    DOUBLE xRef[3]; /*... ponto de referencia*/
    DOUBLE tempMax,tempMed;
/*...*/    
    Elmt elm;     
    Node node;
    Face face;
    Pnode noIncid;
/*...*/
    MeshQuality mQuality;
  }Mesh;
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
