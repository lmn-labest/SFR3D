#ifndef _MESH_
  #define _MESH_
  #include<time.h>
  #include<Adjcency.h>
  #include<Reord.h>
  #include<Mystdbool.h>
  #include<Define.h>

/*... Material*/
  typedef struct Material{
    double *prop;      /*valores numericos da propriedade*/
    short  *type;      /*tipo da celula*/ 
  } Material;
/*...................................................................*/

/*...*/
  typedef struct Adjacency{
    INT  *nelcon;
    short *nViz;
  }Adjacency;
/*...................................................................*/  
/*... Elementos*/
  typedef struct Elmt{
    INT   *node;      /*conectividades*/
    short  *mat;       /*materiais do elementos*/   
    short  *nen;       /*numero de no por elemento*/
    short  *geomType;  /*tipo geometrio do elemento*/
    short  *faceRt1;   /*tipo de condicao de contorno na face*/
    short  *rNum;       /*renumeracao dos elementos*/                
    double *faceSt1;   /*valor da condicao de contorno na face*/
  }Elmt;
/*...................................................................*/
  
/*... nos*/
  typedef struct Node{
    double *x;      /*coordenadas*/
    double *w;
    double *pressure;
    double *temp;
    INT   *nno; 
  }Node;
/*...................................................................*/
  
/*... Malha*/
  typedef struct Mesh{
    INT nnode;/*numero de nos*/
    INT numel;/*numero de elementos*/
    short ndm;     /*dimensao*/
    short ndfT[MAX_TRANS_EQ];   /*graus de liberdade para o problema*/
    short numat;   /*numero maximo de materias no dominio*/
    short maxNo;   /*numero maximo de nos por elemento*/
    short maxViz;  /*numero maximo de vizinhos que um elemento possui*/
    Elmt elm;     
    Node node;
    Material material;
    Adjacency adj;
  }Mesh;
/*...................................................................*/

/* variaveis do preparticionador*/
  typedef struct PartMesh{
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
