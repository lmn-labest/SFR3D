#ifndef _DEFINE_H_
  #define _DEFINE_H_
	
/*...*/
  #include<stdlib.h>
  #include<stdio.h>
/*...................................................................*/

/*...*/
  #define PI     3.14159265358979323846e0
  #define D5DIV6 8.33333333333333333333e-1
  #define D1DIV6 1.66666666666666666666e-1
  #define D1DIV3 3.33333333333333333333e-1
  #define D2DIV3 6.66666666666666666666e-1
/*...................................................................*/

/*... residuos*/
  #define RSCALED    1
  #define RSQRT      2 
  #define RSCALEDSUM 3

/*... trasiente*/
  #define EULER    0
  #define BACKWARD 1
/*...................................................................*/

/*...*/
  #define TDF 1          /*Discretizacao temporal em diferenca finita*/
  #define TVF 2          /*Discretizacao temporal em volume finito". 
/*...................................................................*/

/*...*/
  #define COEFDIF                    0
  #define DINAMICVISCOSITY           0
  #define DENSITY                    1
  #define COEFDIFFFLUID              2
  #define SPECIFICHEATCAPACITYFLUID  3
/*...................................................................*/

/*...*/
  #define SIMPLE  1
  #define SIMPLEC 2
  #define SIMPLER 3
  #define PISO    4
/*...................................................................*/

/*...*/
  #define MAX_MPI_PROCESS 2048
/*...................................................................*/

/*... */
  #define KELVINCONV 273.15e0
/*...................................................................*/

/*...*/
  #define DENSITY_LEVEL 3
  #define SHEAT_LEVEL   3
/*...................................................................*/

/*... conv radianos para graus*/
  #define radToDeg(teta) 180.e0*teta/PI
/*...................................................................*/

/*... tipo de cargr no volume*/
  #define PCCELL       1 /*valor pescrito*/
  #define SCCELL       2 /*fonte*/
/*...................................................................*/
  
/*... ADVECTION*/
  #define NFUNCLIMTFACE 6
  #define NFUNCNVD      7
/*...................................................................*/

/*... ADVECTION*/
  #define FOUP          1
/*...*/
  #define CD            2
/*...*/
  #define SOUP          3
/*...*/
  #define TVD           4
/*...*/
  #define NVD           5

/*... ADVECTION - TVD EDGEBASE LIMIT*/
  #define VANLEERFACE   1
  #define VANALBADAFACE 2
  #define MIDMODFACE    3
  #define OSHERFACE     4
  #define SUPERBEEFACE  5
  #define SPLITUPCFACE  6
/*...................................................................*/

/*... ADVECTION - NVD EDGEBASE LIMIT*/
  #define BCD_NVD       1
  #define MUSCL_NVD     2
  #define SMART_NVD     3
  #define MSMART_NVD    4
  #define SUPERBEE_NVD  5    
  #define MSUPERBEE_NVD 6
  #define STOIC_NVD     7    
/*...................................................................*/

/*... ADVECTION*/
  #define FBASE 1 /*limitacao por face*/
  #define VBASE 2 /*limitacao por volume*/
/*...................................................................*/

/*... DIFUSSION*/
  #define ORTHOGONAL  1
  #define MINIMAL     2
  #define ORTHOGONALC 3
  #define OVERRELAXED 4
/*...................................................................*/

/*... tipos de CC (faces)*/
  #define DIRICHLETBC      1
  #define NEUMANNBC        2  /*para caso de fluxo nao nulo*/
  #define ROBINBC          3
  #define SINBC            4
  #define CONST            5
  #define INLET            6
  #define OUTLET           7
  #define MOVEWALL         8
  #define INLETSTAICTPRES  9
  #define OPEN            10 
  #define STATICWALL      -1  /*parede impermevel para escoamento*/
/*...................................................................*/

/*...*/
  #define MAXLOADPARAMETER 10
  #define MAXLOADD1        200
  #define MAXLOADT1        200
  #define MAXLOADFLUID     200
/*...................................................................*/

/*...*/
  #define DIFPROP        5  /*numero de propriedade de 
                            problemas difusao pura*/
  #define MAXPROP        5  /*numero maximo de propriedades*/
  #define MAXMAT       200 /*numero maximo de materias*/
  #define MAX_TRANS_EQ   3 /*numero maximo de equacoes de transporte*/ 
  #define MAX_DIF_EQ     3 /*numero maximo de equacoes de difusa*/ 
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

/*... simple*/
  #define SZERO 1.e-16
/*...................................................................*/

/*...*/
  #define MAX_NDF 5
/*...................................................................*/

/*... solver*/
  #define PCG          1
  #define PBICGSTAB    2
  #define PBICGSTABL2  3
  #define GMRES        4
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
  #define CSRCCOO 6
/*...................................................................*/

/*... produto interno*/
  #define DOT      1
	#define DOTI2    2
	#define DOTI4    3
	#define DOTI6    4
	#define DOTI8    5
	#define DOTO2    6
	#define DOTO4    7
	#define DOTO6    8
	#define DOTO8    9
	#define DOTO2I2 10
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
//typedef int INT;
  #define INTC       "int"
  #define DOUBLE    double
//typedef double DOUBLE;
  #define DOUBLEC  "double"
  
  #ifdef _MSC_VER
/*  #define LONG_INT __int64*/
    typedef __int64 LONG_INT;
    #define RESTRICT __restrict
  #else
/*  #define LONG_INT long*/
    typedef __int64 LONG_INT;
    #define RESTRICT restrict
  #endif
/*...................................................................*/

/*... macro para acesso matricial em vetores*/
  #define   MAT2D(i,j,vector,col)           vector[(i)*(col)+(j)]
  #define   MAT3D(i,j,k,vector,col1,col2)   vector[(i)*(col1)*(col2)+(col2)*(j)+(k)]
/*...................................................................*/

/*... definicao de funcoes*/
  #define min(a, b)  (((a) < (b)) ? (a) : (b))
  #define max(a, b)  (((a) > (b)) ? (a) : (b))
  #define vectorPlusOne(v,n,i)  for(i=0;i<n;i++) v[i]++ 
  #define vectorMinusOne(v,n,i) for(i=0;i<n;i++) v[i]--  
/*...................................................................*/

/*...*/
  DOUBLE gravity[3];
/*...................................................................*/


#endif/*_DEFINE_H_*/
