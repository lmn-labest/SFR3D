#ifndef _ERRO_H_
  #define _ERRO_H_
  
  #define EXIT_HELP              1
  #define EXIT_MESH_COMM        -1
  #define EXIT_SOLVER           -2
  #define EXIT_SOLVER_CONFIG    -3
  #define EXIT_NRP_ET           -4
  #define EXIT_READ_LOOP        -5
  #define EXIT_PROG           -111

/*... Saida de Erro*/                                                  
  #define ERRO_RCM fprintf(stderr,"\nrcm - fatal error!\n")

  #define ERRO_OP(file,func,op)\
    {fprintf(stderr,"Opcao %d e invalida!!\n",op);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\n",file,func);\
    exit(EXIT_FAILURE);}
 
  #define ERRO_OP_NEW(file,func,line,str,op)\
    {fprintf(stderr,"%s\nOpcao %d e invalida!!\n",str,op);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
            ,file,func,line);\
    exit(EXIT_FAILURE);}

  #define ERRO_OP_WORD(file,fileSrc,func,line,str,op,cod)\
    {fprintf(stderr,"Erro !!\nWriting to file .. ");\
    fprintf(file,"%s\nOpcao %s e invalida!!\n",str,op);\
    fprintf(file,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
            ,fileSrc,func,line);\
    exit(cod);}

  #define ERRO_GERAL(file,fileSrc,func,line,str,cod)\
    {fprintf(stderr,"Erro !!\nWriting to file .. ");\
    fprintf(file,"Erro: %s!!\n",str);\
    fprintf(file,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
           ,fileSrc,func,line);\
    exit(cod);}

  #define ERRO_MALLOC(point,str,line,file,func)\
     if(point == NULL){\
     fprintf(stderr,"Erro na alocacao do vetor %s\n",str);\
     fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
            ,file,func,line);\
     exit(EXIT_FAILURE);}
  
  #define ERRO_MALLOC_MPI(point,str,line,file,func,a)\
     if(point == NULL){\
     fprintf(stderr,"Erro na alocacao do vetor %s\n",str);\
     fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
            ,file,func,line);\
     a = - 1;}

  #define ERRO_POL_READ(a,b,file,func,line)\
    if(a > b)\
    {\
      fprintf(stderr,"Numero de grau do polinomials execedido!!");\
      fprintf(stderr,"%d > %d",a,b);\
      fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
             ,file,func,line);\
      exit(EXIT_FAILURE);\
    }

  #define ERRO_NORM(line,file,func,lin,col)\
    {fprintf(stderr,"Erro no calulo da norma!!\n");\
    fprintf(stderr,"Numero de linha  : %d\n"\
                   "Numero de colunas: %d\n",lin,col);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
          ,file,func,line);\
    exit(EXIT_FAILURE);}

  #define ERRO_READ_LOOP(line,file,func,cod)\
    {fprintf(stderr,"Erro de leitura!!\n");\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
          ,file,func,line);\
    exit(cod);}

/*...................................................................*/
#endif /*_ERRO_H_*/
