#ifndef _MESH_
  #define _MESH_
  #include<Adjcency.h>
  #include<Reord.h>
  #include<HccaStdBool.h>
  #include<Define.h>

/*... Material*/
  typedef struct{
    double *prop;      /*valores numericos da propriedade*/
    short  *type;      /*tipo de calculo da celula*/ 
  } Material;
/*...................................................................*/

/*...*/
  typedef struct{
    INT  *nelcon;
    short *nViz;
  }Adjacency;
/*...................................................................*/  

/*...*/
  typedef struct{
    DOUBLE    dt;
    DOUBLE     t;
    short   type;
    bool    flag;
    DOUBLE total;
    INT timeStep;
  }Temporal;
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

/*... loads*/
  typedef struct{
    short type;                     /*tipo*/
    short np;                       /*numero de particoes*/  
    DOUBLE par[MAXLOADPARAMETER];
  }Loads;
  Loads  loadsD1[MAXLOAD1]; /*tipo de cargas (difusao pura)*/
/*...................................................................*/

/*...*/
  typedef struct{
    int maxIt;
    DOUBLE tol;                               
  }NonLinear;
/*...................................................................*/  

/*... Elementos*/
  typedef struct{
    INT    *node;       /*conectividades*/
    short  *mat;       /*materiais do elementos*/   
    short  *nen;       /*numero de no por elemento*/
    short  *geomType;  /*tipo geometrio do elemento*/
    short  *rNum;      /*renumeracao dos elementos*/ 
/*... */               
    short  *faceRt1;   /*condicao de contorno na face (transporte)*/
    DOUBLE *faceSt1;   /*valor da condicao de 
                        contorno na face (transporte)*/
/*...*/
    short *faceRd1;    /*condicao  contorno na face (difusa pura)*/
    short *faceLoadD1; /*tipo de carga contorno na face (difusa pura)*/
/*...................................................................*/
    Geom   geom;       
    DOUBLE *pressure;
/*...*/
    DOUBLE *temp;       /*temperatura*/
    DOUBLE *gradTemp;   /*gradiente da temperatura*/
    DOUBLE *rCellTemp;  /*residuo da celula*/
/*...*/
    DOUBLE *densityUd1; /*massa especifica do material uD1*/
    DOUBLE *uD1 ;       /*difusao pura uD1*/
    DOUBLE *u0D1;       /*difusao pura uD1*/
    DOUBLE *gradUd1;    /*gradiente da difusao pura uD1*/
    DOUBLE *rCellUd1;   /*residuo da celula*/
/*...*/
    DOUBLE *leastSquare; /*matriz de aproxima leastSquare*/
    DOUBLE *leastSquareR;/*fatoracao QR*/
    Material material;
    Adjacency adj;
  }Elmt;
/*...................................................................*/
  

/*... nos*/
  typedef struct{
    DOUBLE *x;      /*coordenadas*/
    DOUBLE *w;
    DOUBLE *pressure;
    DOUBLE *temp;      /*temperatura*/
    DOUBLE *uD1;       /*difusao pura uD1*/
    DOUBLE *gradTemp;  /*gradiente da temperatura*/
    DOUBLE *gradUd1 ;  /*gradiente da difusao pura uD1*/
    INT    *nno; 
  }Node;
/*...................................................................*/

/*... meshQuality*/
  typedef struct{
    DOUBLE volume;
    DOUBLE nonOrthMed;
    DOUBLE nonOrthMax;
    DOUBLE skewMed;
    DOUBLE skewMax;
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
    INT nnode;     /*numero de nos*/
    INT numel;     /*numero de elementos*/
    INT numelNov;  /*numero de elementos sobrepostos*/
    INT nnodeOv;   /*numero de nos em elementos sobrepostos*/
    INT nnodeNov;  /*numero de nos em elementos nao sobrepostos*/
    short ndm;     /*dimensao*/
    short ndfT[MAX_TRANS_EQ];   /*graus de liberdade 
                                  para o problema de transporte*/
    short ndfD[MAX_DIF_EQ];   /*graus de liberdade 
                                  para o problema de difusa pura*/
    short numat;   /*numero maximo de materias no dominio*/
    short maxNo;   /*numero maximo de nos por elemento*/
    short maxViz;  /*numero maximo de vizinhos que um elemento possui*/
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
    INT  nCom;     /*numero de valores a serem comunicados*/
    INT  *iaComNo; /*ponteiros para o arranjo buffer de comunicacao*/                        
    INT   *fMap;  /*numeracao do buffer de comunicacao*/    
    DOUBLE  *xb;    /*valores no buffer de comunicacao(DOUBLE)*/
    INT     *xi;    /*valores no buffer de comunicacao(INT)*/
    unsigned short nVizPart;
    short *vizPart;
    
  }InterfaceNo;

/* variaveis do preparticionador*/
  typedef struct{
    bool fPartMesh;
    bool fPrintMesh;
    bool fPrintMeshPart;
    Interface   iEl;
    InterfaceNo iNo;
    INT nno1;      /*numero de nos apenas da particao local*/
    INT nnG ;      /*numero de nos total da malha global*/
    INT elG ;      /*numero de elementos total da malha global*/
    INT *elLG;     /*mapa de local-> global de celulas*/
    INT *noLG;     /*mapa de local-> global de nos*/
    INT *np;       /*master apenas*/
    INT *ep;       /*master apenas*/
  }PartMesh;
/*...................................................................*/
#endif/*_MESH_*/
