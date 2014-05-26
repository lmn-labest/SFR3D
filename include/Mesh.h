#ifndef _MESH_
  #define _MESH_
  #include<time.h>
  #include<Mystdbool.h>
/*numero maximo de propriedade*/
  #define MAXPROP      15 
/*numero maximo de materiais*/
  #define MAXMAT      200
/*vtk elmentos*/
  #define VTK_TRIA  5
  #define VTK_QUAD  9
  #define VTK_TETR  10
  #define VTK_HEXA  12
  #define MAX_TRANS_EQ 3 
/*...................................................................*/

/*... Material*/
  typedef struct Material{
    double *prop;      /*valores numericos da propriedade*/
    short  *type;      /*tipo da celula*/ 
  } Material;
/*...................................................................*/
  
/*... Elementos*/
  typedef struct Elmt{
    int    *node;      /*conectividades*/
    short  *mat;       /*materiais do elementos*/   
    short  *nen;       /*numero de no por elemento*/
    short  *geomType;  /*tipo geometrio do elemento*/
    short  *faceRt1;   /*tipo de condicao de contorno na face*/
    double *faceSt1;   /*valor da condicao de contorno na face*/
    long   *nel; 
  }Elmt;
/*...................................................................*/
  
/*... nos*/
  typedef struct Node{
    double *x;      /*coordenadas*/
    double *w;
    double *pressure;
    double *temp;
    long   *nno; 
  }Node;
/*...................................................................*/
  
/*... Malha*/
  typedef struct Mesh{
    long int nnode;/*numero de nos*/
    long int numel;/*numero de elementos*/
    short ndm;   /*dimensao*/
    short ndfT[MAX_TRANS_EQ];   /*graus de liberdade para o problema*/
    short numat; /*numero maximo de materias no dominio*/
    short maxno; /*numero maximo de nos por elemento*/
    Elmt elm;     
    Node node;
    Material material;
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
