#ifndef _ERRO_H_
  #define _ERRO_H_
/*... Saida de Erro*/                                                  
  #define ERRO_RCM fprintf(stderr,"\nrcm - fatal error!\n")

  #define ERRO_OP(file,func,op)\
    fprintf(stderr,"Opcao %d e invalida!!\n",op);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\n",file,func);\
    exit(EXIT_FAILURE);
 
  #define ERRO_OP_NEW(file,func,line,str,op)\
    fprintf(stderr,"%s\nOpcao %d e invalida!!\n",str,op);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
            ,file,func,line);\
    exit(EXIT_FAILURE);

  #define ERRO_OP_WORD(file,func,line,str,op)\
    fprintf(stderr,"%s\nOpcao %s e invalida!!\n",str,op);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
            ,file,func,line);\
    exit(EXIT_FAILURE);

  #define ERRO_GERAL(file,func,line,str)\
    {fprintf(stderr,"Erro: %s!!\n",str);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
           ,file,func,line);\
    exit(EXIT_FAILURE);}

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
      fprintf(stderr,"Numero de grau do polinmios execedido!!");\
      fprintf(stderr,"%d > %d",a,b);\
      fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
             ,file,func,line);\
      exit(EXIT_FAILURE);\
    }

  #define ERRO_NORM(line,file,func,lin,col)\
    fprintf(stderr,"Erro no calulo da norma!!\n");\
    fprintf(stderr,"Numero de linha  : %d\n"\
                   "Numero de colunas: %d\n",lin,col);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
          ,file,func,line);\
    exit(EXIT_FAILURE);

/*...................................................................*/
#endif /*_ERRO_H_*/
