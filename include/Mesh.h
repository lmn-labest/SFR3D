#ifndef _MESH_
  #define _MESH_
  #include<time.h>
  #include<Mystdbool.h>
/*numero maximo de propriedade*/
  #define MAXPROP      15 
/*numero maximo de materiais*/
  #define MAXMAT      200
/*...................................................................*/

/*... Material*/
  typedef struct Material{
    double *prop;
    short  *type;/*tipo de elemento */
  } Material;
/*...................................................................*/
  
/*... Elementos*/
  typedef struct Elmt{
    int *node;      /*conectividades*/
    int *mat;       /*materiais do elementos*/   
    int *type;      /*tipo do geometrico do elemento*/
    int *nen;       /*numero de no por elemento*/
    int *eloads;    /*cargas nos elementos*/
  }Elmt;
/*...................................................................*/
  
/*... nos*/
  typedef struct Node{
    double *x;      /*coordenadas*/
  }Node;
/*...................................................................*/
  
/*... Malha*/
  typedef struct Mesh{
    long int nnode;/*numero de nos*/
    long int numel;/*numero de elementos*/
    short dim;   /*dimensao*/
    short ndf;   /*graus de liberdade para o problema acustico*/
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
