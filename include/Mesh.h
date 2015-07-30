#ifndef _MESH_
  #define _MESH_
  #include<time.h>
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
  Loads  loadsD1[MAXLOAD1]; /*tipo de cargas (difusa pura)*/
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
    short  *faceRd1;  /*condicao  contorno na face (difusa pura)*/
    short *faceLoadD1; /*tipo de carga contorno na face (difusa pura)*/
/*...................................................................*/
    Geom   geom;       
    DOUBLE *pressure;
    DOUBLE *temp;       /*temperatura*/
    DOUBLE *gradTemp;   /*gradiente da temperatura*/
    DOUBLE *rCellTemp;  /*residuo da celula*/
    DOUBLE *uD1 ;       /*difusao pura uD1*/
    DOUBLE *gradUd1;    /*gradiente da difusao pura uD1*/
    DOUBLE *rCellUd1;   /*residuo da celula*/
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
    INT   *nno; 
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
  
/*... Malha*/
  typedef struct{
    short rcGrad; /*tipo de rescontrucao de gradiente*/                     
    INT nnode;/*numero de nos*/
    INT numel;/*numero de elementos*/
    short ndm;     /*dimensao*/
    short ndfT[MAX_TRANS_EQ];   /*graus de liberdade 
                                  para o problema de transporte*/
    short ndfD[MAX_DIF_EQ];   /*graus de liberdade 
                                  para o problema de difusa pura*/
    short numat;   /*numero maximo de materias no dominio*/
    short maxNo;   /*numero maximo de nos por elemento*/
    short maxViz;  /*numero maximo de vizinhos que um elemento possui*/
    Elmt elm;     
    Node node;
    NonLinear nlTemp;
    NonLinear nlD1;
    MeshQuality mQuality;
  }Mesh;
/*...................................................................*/

/* variaveis do preparticionador*/
  typedef struct{
    time_t tmetis;
    int *np;
    int *ep;
    int *noGL;
    int *noLG;
    int *elGL;
    int *elLG;
    int *ranks;
    int *fmap;
    int  nnof[2];
    int *rcvs;
    int *dspl;
    int  my_nnode;
    int  my_nno1;
    int  my_nno2;
    int  my_nno3;
    int  my_nno4;
    int  my_nno1a;
    int  my_numel;
    int  my_numel_nov;
    int  my_numel_ov;
    int  sizes;
    int  nprcs;
    bool ovlp;
  }PartMesh;
/*...................................................................*/
#endif/*_MESH_*/
