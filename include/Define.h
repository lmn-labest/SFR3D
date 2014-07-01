#ifndef _DEFINE_
  #define _DEFINE_
/*numero maximo de propriedade*/
  #define MAXPROP      15 /*numero maximo de propriedades*/
  #define MAXMAT      200 /*numero maximo de materias*/
  #define MAX_TRANS_EQ 3 /*numero maximo de equacoes de transporte*/ 
/*...................................................................*/

/*... vtk elmentos*/
  #define VTK_TRIA      5
  #define VTK_QUAD      9
  #define VTK_TETR     10
  #define VTK_HEXA     12
/*...................................................................*/

/*... definicao do tipo de inteiros usados*/
  #define INT   long
  #define INTC "long"
/*...................................................................*/

/*... macro para acesso matricial em vetores*/
  #define    VET(i,j,vector,col)        vector[i*col+j]
/*...................................................................*/

/*... definicao de funcoes*/
  #define min(a, b)  (((a) < (b)) ? (a) : (b))
  #define max(a, b)  (((a) > (b)) ? (a) : (b))
  #define vectorPlusOne(v,n,i)  for(i=0;i<n;i++) v[i] = v[i] + 1; 
  #define vectorMinusOne(v,n,i) for(i=0;i<n;i++) v[i] = v[i] - 1; 
/*...................................................................*/

/*... Saida de Erro*/                                                  
  #define ERRO_RCM fprintf(stderr,"\nrcm - fatal error!\n")
/*...................................................................*/
#endif/*_DEFINE_*/
