#ifndef _FILE_
  #define _FILE_
  #include<stdio.h>
  #include<stdlib.h>
  #include<Mystdbool.h>
  #include<math.h>
  #include<string.h>
  #define MAX_STR_LEN_IN 129 /*tamanho do arquivo de entrada maximo*/
  #define MAX_STR_LEN_SUFIXO 80  /*tamanho do arquivo de entrada maximo*/
  #define SIZEMAX 200 /*tamanho do arquivo de saida maximo*/
  #define MAX_STR_NUMBER 65 /*long integer*/
  #define MAX_EXT        20 /*numero maximo da extencao*/ 
  #define WORD_SIZE   120
  #define LINE_SIZE   120
  #define MAX_LINE  200
/*... uso nadepuração de erros*/
  #ifdef SRC_NAME
    #undef  SRC_NAME
  #endif  
  #define SRC_NAME "scr/File.c"
/*...................................................................*/  
  void  iota(long,char*);
  FILE* open_file(char *,char *);
  void  fname(char*,long,int,char**);
  void  readmacro(FILE*,char*,bool);
  void  clear_line(char *);
  int   rl(FILE *,char *);
  int   getnumprop(FILE *);
/*...*/  
  char macros[MAX_LINE][WORD_SIZE];/*todas as macros lidas no arquivo*/
  int nmacro;
/*....................................................................*/  
#endif
