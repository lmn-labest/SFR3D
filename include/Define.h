#ifndef _DEFINE_H_
  #define _DEFINE_H_

/*...*/
  #include<stdlib.h>
  #include<stdio.h>
/*...................................................................*/

/*...*/
  #define PI 3.14159265358979323846
/*...................................................................*/

/*... trasiente*/
  #define EULER    0
  #define BACKWARD 1
/*...................................................................*/

/*...*/
  #define COEFDIF 0
  #define DENSITY 1
/*...................................................................*/

/*...*/
  #define MAX_MPI_PROCESS 2048
/*...................................................................*/

/*... */
  #define KELVINCONV 273.15e0
/*...................................................................*/

/*...*/
  #define DENSITY_LEVEL 2
/*...................................................................*/

/*... conv radianos para graus*/
  #define radToDeg(teta) 180.e0*teta/PI
/*...................................................................*/

/*... tipo de cargr no volume*/
  #define PCCELL       1 /*valor pescrito*/
  #define SCCELL       2 /*fonte*/
/*...................................................................*/

/*... tipos de CC (faces)*/
  #define DIRICHLETBC  1
  #define NEUMANNBC    2
  #define ROBINBC      3
  #define SINBC        4
  #define CONST        5
/*...................................................................*/

/*...*/
  #define MAXLOADPARAMETER 10
  #define MAXLOAD1         200
/*...................................................................*/

/*...*/
  #define DIFPROP      5  /*numero de propriedade de 
                            problemas difusao pura*/
  #define MAXPROP      5  /*numero maximo de propriedades*/
  #define MAXMAT      200 /*numero maximo de materias*/
  #define MAX_TRANS_EQ 3 /*numero maximo de equacoes de transporte*/ 
  #define MAX_DIF_EQ   3 /*numero maximo de equacoes de difusa*/ 
/*...................................................................*/

/*... Tipo geometrica da celula*/
  #define LINECELL  1
  #define TRIACELL  2
  #define QUADCELL  3  
  #define TETRCELL  4  
  #define HEXACELL  5  
/*...................................................................*/

/*...cellLoop*/
  #define  MAX_NUM_RES       28
  #define  MAX_NUM_NODE       8
  #define  MAX_NUM_PONT     168
  #define  MAX_NUM_FACE       6
  #define  MAX_NUM_NODE_FACE  4
  #define  MAX_SN            24 
  #define  MAX_NDM            3 
/*...................................................................*/

/*... reconstrucao de gradiente*/
  #define  RCGRADGAUSSC 1 
  #define  RCGRADGAUSSN 2 
  #define  RCLSQUARE    3 
  #define  RCLSQUAREQR  4    
/*...................................................................*/

/*...*/
  #define MAX_NDF 5
/*...................................................................*/

/*... solver*/
  #define PCG        1
  #define PBICGSTAB  2
/*...................................................................*/

/*... STORAGE*/
  #define BANDCSRMAX     0
  #define BANDCSRMED     1
  #define BANDCSRMIN     2
/*...................................................................*/

/*... STORAGE*/
  #define CSR     1
  #define CSRD    2
  #define CSRC    3
  #define ELLPACK 4
  #define CSRDCOO 5
/*...................................................................*/

/*... vtk elmentos*/
  #define VTK_LINE      3
  #define VTK_TRIA      5
  #define VTK_QUAD      9
  #define VTK_TETR     10
  #define VTK_HEXA     12
/*...................................................................*/

/*... definicao do tipo de inteiros usados*/
  #define INT         int 
  #define INTC       "int"
  #define DOUBLE    double
  #define DOUBLEC  "double"
/*...................................................................*/

/*...*/
  DOUBLE oneDivTree;
/*...................................................................*/

/*... macro para acesso matricial em vetores*/
  #define   MAT2D(i,j,vector,col)           vector[i*col+j]
  #define   MAT3D(i,j,k,vector,col1,col2)   vector[i*col1*col2+col2*j+k]
/*...................................................................*/

/*... definicao de funcoes*/
  #define min(a, b)  (((a) < (b)) ? (a) : (b))
  #define max(a, b)  (((a) > (b)) ? (a) : (b))
  #define vectorPlusOne(v,n,i)  for(i=0;i<n;i++) v[i]++ 
  #define vectorMinusOne(v,n,i) for(i=0;i<n;i++) v[i]--  
/*...................................................................*/

#endif/*_DEFINE_H_*/
