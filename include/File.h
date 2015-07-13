#ifndef _FILE_
  #define _FILE_
  #include<stdio.h>
  #include<stdlib.h>
  #include<HccaStdBool.h>
  #include<math.h>
  #include<string.h>
  #include<Define.h>
  #define MAX_STR_LEN_IN    129 /*tamanho do arquivo de entrada maximo*/
  #define MAX_STR_LEN_SUFIXO 80  /*tamanho do arquivo de entrada maximo*/
  #define SIZEMAX           200 /*tamanho do arquivo de saida maximo*/
  #define MAX_STR_NUMBER     65 /*INT integer*/
  #define MAX_EXT            20 /*numero maximo da extencao*/ 
  #define WORD_SIZE         120
  #define LINE_SIZE         120
  #define MAX_LINE          200

  #define FITPLOTD1           1
  #define FITPLOTD2           2
  #define FITPLOTT1           3
  #define FITPLOTT2           4
  #define FITPLOTTEMP         5
   
/*...................................................................*/  
  typedef struct{
    bool  bVtk;
    bool  fItPlotRes;
    bool  fItPlot;
    FILE *fileItPlot[6];
  }FileOpt;

  void  iota(INT,char*);
  FILE* openFile(char *,char *);
  void  fName(char*,INT,int,int,char**);
  void  readMacro(FILE*,char*,bool);
  void  clearLine(char *);
  int   rl(FILE *,char *);
  int   getNumProp(FILE *);
/*...*/  
  char macros[MAX_LINE][WORD_SIZE];/*todas as macros lidas no arquivo*/
  int  nmacro;
/*....................................................................*/  
#endif
