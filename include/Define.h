#ifndef _DEFINE_H_
  #define _DEFINE_H_
	
/*...*/
  short debug_flag;
/*...*/
  #include<stdlib.h>
  #include<stdio.h>
/*...................................................................*/

/*...*/
  #define E_WALLMODEL   9.793e0
  #define VANDRIEST    26.e0
  #define VONKARMAN     0.4187e0
  #define AKADER        0.01e0
  #define BKADER        5.e0
/*...................................................................*/

/*...*/
  #define SMAGORINSKY 1
  #define WALEMODEL   2
  #define VREMAN      3
  #define DYNAMIC     4
  #define SIGMAMODEL  5
  #define MIXED       6
  #define BARDINA     7
  #define CLARK       8 
  #define BARDINAMOD  9
  #define ONEEQK     10
/*...................................................................*/
  
/*...*/
  #define LESFUNCMODEL       0
  #define LESSTRUMODEL       1
  #define LESMIXEDMODEL      2
  #define LESMIXEDTWOMODEL   3
  #define LESFUNCMODELONEEQK 4
/*...................................................................*/

/*...*/
  #define LES 0
/*...*/
  #define ESTMODEL 0
  #define FUNMODEL 1
/*...................................................................*/

/*... NEAR-WALL_MODEL*/
  #define STANDARDWALL 1
  #define ENHANCEDWALL 2
/*...................................................................*/

/*...*/
  #define NWALLPAR     4
/*...................................................................*/

/*...*/
  #define LDYNAMIC      1
  #define GDYNAMIC      2
  #define GDYNAMICMOD   3
  #define TWOPARDYNAMIC 4 
/*...................................................................*/

/*...*/
  #define SL_FUEL 0
  #define SL_AIR  1
  #define SL_PROD 2
/*...................................................................*/

/*...*/
  #define OUTPUT_FOR_FILE   1
  #define OUTPUT_FOR_SCREEN 2
/*...................................................................*/

/*... combustion*/
  #define MAXELEMENT   10
  #define MAXNAMELENSP 10
  #define MAXREAC       5
  #define MAXSPECIES   10
  #define HFORMATION    1
  #define HCOMBUSTION   2

/*...*/
  #define ARRHENIUS   1
  #define EDC         2 
  #define EDM         3 
/*...*/
  #define N_TERMS_REACTOR 5
/*...*/
  #define FDS_EDC                 1 
  #define FLUENT_EDC              2
  #define FLUENT_CONST_TMIX_EDC   3
  #define PANJWANI_EDC            4
  #define PANJWANI_CONST_TMIX_EDC 5
/*... Kg/kmol*/
  #define MW_O 15.9994e0
  #define MW_H  1.00794e0
  #define MW_N 14.00674e0 
  #define MW_C 12.01100e0
/*...................................................................*/

/*...*/
  #define POL            1
  #define NASAPOL7       2 
/*... viscosidade*/
  #define SUTHERLAND     2
  #define FDSVISCOSITY   3
  #define FDSTHERMALCOND 3
  #define WILKELAW       4
  #define FDSDIFF        2
  #define HIRSCHDIFF     3  
/*... densidade*/
  #define IDEALGAS   2
/*... densidade*/
  #define INCIDEALGAS   3
/*...................................................................*/

/*...*/
  #define TRAPEZIO 1
  #define SIMPSON  2 
/*...................................................................*/

/*...*/
  #define TREF       298.15e+00    /*Kelvin          */
  #define PREREF     1.01325e+05   /*Pa              */
  #define IDEALGASR  8.3144621e+03 /*J/(kmol.kelvin) */
  #define IDEALGASRC 1.9872063e+00  /*Cal/(mol.kelvin) */
//#define IDEALGASR  8.3144598e+03 /*J/(Kmol.kelvin) */
  #define MMOLARAR   2.896e+00     /*kg/Kmol         */
/*...................................................................*/

/*...*/
  #define TEMP_FOR_ENTHALPY(cp,t,tr) ( (cp)*(t-tr) )
  #define ENTHALPY_FOR_TEMP(cp,hs,tr) ( (hs/cp) + (tr) )
/*...................................................................*/

/*...*/
  #define PRESREF(dRef,R,T,Mg) ((dRef)*(R)*(T)/(Mg)) 
/*...................................................................*/

/*...*/
  #define PROP_UPDATE_NL_LOOP     0
  #define PROP_UPDATE_OLD_TIME    1
/*...................................................................*/

/*...*/
  #define MAXPLODEG  15  
/*...................................................................*/

/*...*/
  #define PI     3.14159265358979323846E0
  #define D5DIV6 8.33333333333333333333E-1
  #define D1DIV6 1.66666666666666666666E-1
  #define D1DIV3 3.33333333333333333333E-1
  #define D2DIV3 6.66666666666666666666E-1
/*...................................................................*/

/*... residuos*/
  #define RSCALED       1
  #define RSQRT         2 
  #define RSCALEDSUM    3
  #define RSCALEDSUMMAX 4
/*...................................................................*/

/*... trasiente*/
  #define EULER    0
  #define BACKWARD 1
/*...................................................................*/

/*...*/
  #define TDF 1          /*Discretizacao temporal em diferenca finita*/
  #define TVF 2          /*Discretizacao temporal em volume finito"*/ 
/*...................................................................*/

/*...*/
  #define COEFDIF                    0
  #define DYNAMICVISCOSITY           0
  #define DENSITY                    1
  #define THERMALCONDUCTIVITY        2
  #define SPECIFICHEATCAPACITYFLUID  3
  #define SPECIEDIFUSSIONFUEL        4
  #define SPECIEDIFUSSIONO2          5
  #define SPECIEDIFUSSIONCO2         6
  #define SPECIEDIFUSSIONH2O         7
  #define SPECIEDIFUSSIONN2          8
/*...................................................................*/

/*...*/
  #define TIME_N_MINUS_2 0
  #define TIME_N_MINUS_1 1
  #define TIME_N         2
/*...................................................................*/

/*...*/
  #define INTPOLFACELINEAR 0
  #define INTPOLFACEVOL    1
  #define INTPOLFACEMED    2
/*...................................................................*/

/*...*/
  #define SIMPLE  1
  #define SIMPLEC 2
  #define SIMPLER 3
  #define PISO    4
/*...................................................................*/

/*...*/
  #define CFL 1
/*...................................................................*/

/*...*/
  #define MAX_MPI_PROCESS 2048
/*...................................................................*/

/*...*/
  #define OWNER(cell,nel) ((cell) == (nel) ? 1:-1)
/*...................................................................*/

/*... */            
  #define KELVINCONV 273.15E+00
  #define CELSIUS_FOR_KELVIN(t) ((t)+(KELVINCONV)) 
  #define KELVIN_FOR_CELSIUS(t) ((t)-(KELVINCONV)) 
  #define TEMP(tc,t,fKelvin)  if((fKelvin)){ (tc) = (t);}\
                              else{(tc) = CELSIUS_FOR_KELVIN((t));}
/*...................................................................*/

/*...*/
  #define DVISCOSITY_LEVEL      1
  #define DENSITY_LEVEL         3
  #define SHEAT_LEVEL           3
  #define TCONDUCTIVITY_LEVEL   1
  #define COEFDIFF_LEVEL        1
  #define SPECIEDIFUSSION_LEVEL 1   
/*...................................................................*/

/*... conv radianos para graus*/
  #define radToDeg(teta) (180.e0*(teta)/(PI))
/*...................................................................*/

/*... tipo de cargr no volume*/
  #define PCCELL       1 /*valor pescrito*/
  #define SCCELL       2 /*fonte*/
/*...................................................................*/
  
/*... ADVECTION*/
  #define NFUNCLIMTFACE 5
  #define NFUNCNVD      9
/*...................................................................*/

/*... ADVECTION*/
  #define FOUP          1
  #define CD            2
  #define SOUP          3
  #define TVD           4
  #define NVD           5
  #define LUST          6

/*...*/
  #define NPADV         1

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
  #define MINMOD_NVD    8   
  #define MBCD_NVD      9 
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

/*... BUOYANT*/
  #define BUOYANT_HYDROSTATIC 1
  #define BUOYANT_PRGH        2
  #define BUOYANT_RHOREF      3      
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
  #define CONVECTIONHEAT  11
  #define STATICWALL      -1  /*parede impermevel para escoamento*/
  #define SLIP            12 /*parede perfeitamente lisa*/
  #define INLETTOTALPRES  13
/*...................................................................*/

/*...*/
  #define MAXLOADPARAMETER 20
  #define MAXLOADD1        200
  #define MAXLOADT1        200
  #define MAXLOADFLUID     200
  #define MAXINTERPOL      100
/*...................................................................*/

/*...*/
  #define LVARCONST          0
  #define LFUNCPARABOLIC     1
/*...................................................................*/

/*...*/
  #define DIFPROP        5  /*numero de propriedade de 
                            problemas difusao pura*/
  #define MAXPROP       20  /*numero maximo de propriedades*/
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
  #define PIRACELL  6
  #define DOTCELL   7
/*...................................................................*/

/*...cellLoop*/
  #define  MAX_NUM_RES       28
  #define  MAX_NUM_NODE       8
  #define  MAX_NUM_PONT     168
  #define  MAX_NUM_FACE       6
  #define  MAX_NUM_NODE_FACE  4
  #define  MAX_SN            24 
  #define  MAX_NDM            3 
  #define  MAX_COMB           7
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
  #define MAX_NDF 10
/*...................................................................*/

/*... solver*/
  #define PCG          1
  #define PBICGSTAB    2
  #define PBICGSTABL2  3
  #define GMRES        4
  #define MINRES       5
  #define PARDISO      6
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

/*... definicao do tipo de inteiros usados*/
  #define INT         int   
//typedef int INT;
  #define INTC       "int"
  #define DOUBLE    double
//typedef double DOUBLE;
  #define DOUBLEC  "double"
  #define LDOUBLE   long double
  #define LDOUBLEC  "long double"  

  #ifdef _MSC_VER
/*  #define LONG_INT __int64*/
    typedef __int64 LONG_INT;
    #define RESTRICT __restrict
  #else
    #define LONG_INT long  
/*  typedef __int64 LONG_INT;*/
    #define RESTRICT restrict
  #endif
/*...................................................................*/

/*... macro para acesso matricial em vetores*/
  #define   MAT2D(i,j,vector,col)           (vector[(i)*(col)+(j)])
  #define   MAT3D(i,j,k,vector,col1,col2)   (vector[(i)*(col1)*(col2)+(col2)*(j)+(k)])
/*...................................................................*/

/*... definicao de funcoes*/
  #define min(a, b)  (((a) < (b)) ? (a) : (b))
  #define max(a, b)  (((a) > (b)) ? (a) : (b))
  #define vectorPlusOne(v,n,i)  for(i=0;i<n;i++) (v[i]++) 
  #define vectorMinusOne(v,n,i) for(i=0;i<n;i++) (v[i]--)  
/*...................................................................*/

/*...*/
  #define LDMT           1
  #define LDLT           2
  #define GGT            3
  #define LUKIJ          4
  #define LUIKJ          5
  #define LUKJI          6
  #define LUJKI          7
  #define DOOLITTLEL     8
  #define DOOLITTLELC    9
  #define DOOLITTLECL   10
  #define LUKIJPP       11
  #define DOOLITTLECLPP 12
/*...................................................................*/

/*...*/
  DOUBLE gravity[3];
  FILE *fileLogExc,*fileLogDebug;
/*...................................................................*/


#endif/*_DEFINE_H_*/
