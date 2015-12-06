#include<File.h>
/*********************************************************************
 * openFile : abre um arquivo                                        *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * name  -> nome do arquivo                                          * 
 * mod   -> mode de abertura                                         *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * FILE * ponterio par o arquivo de saida                            * 
 * ------------------------------------------------------------------*
 * *******************************************************************/
FILE* openFile(char *name,char *mod){
FILE *aux;
  if(( aux =fopen(name,mod))==NULL){
      fprintf(stderr,"Erro na abertura do arquivo %s.\n",name);
      exit(EXIT_FAILURE);
    }
  return aux;  
}
/********************************************************************/

/*********************************************************************
 * fname: add as extencoes dos arquivos de saida                     *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * out   -> nome do arquivo de saida                                 * 
 * num1  -> numeracao 1                                              *
 * num2  -> numeracao 2                                              *
 * cod   -> codigo                                                   *
 *          1 -> particianomento vtk novlp                           *
 *          2 -> particianomento vtk  ovlp                           *
 *          3 -> malha por particao vtk                              *
 *          4 -> arquivo de dados mefpar                             *      
 *          6 -> carrgamentos formato vtk                            *      
 *          7 -> arquivo de log de excucao                           *      
 *          8 -> resultados da temperatura formato vtk               *      
 *         10 -> log da iteracao nao linear                          *      
 *         11 -> log do solv                                         *      
 *         12 -> matrix de coeficientes no formato .mtx              *      
 *         13 -> matrix de coeficientes no formato .mtx  (0;1)       *      
 *         14 -> vetor de forcas no formato .mtx                     *      
 *         15 -> resultados da temperatura por iteracao formato vtk  *      
 *         16 -> resultados no formato .csv                          *      
 *         17 -> caregamento nas faces no formato vtk                *      
 *         18 -> caregamento nas faces no formato vtk                *      
 *         60 -> media dos tempos (MPI)                              *      
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * out -> aruivo de saida com a extencao                             * 
 * ------------------------------------------------------------------*
 * *******************************************************************/
void fName(char *name,INT num1,INT num2, int cod ,char **out ){
/*===*/
  char st[MAX_STR_NUMBER];
  char ext[MAX_EXT];
  int size1,size2;
/*===================================================================*/
/**/
/*===*/
  strcpy(ext,"");
  strcpy(*out,"");
  switch( cod ) {
/*.... vtk no-overlanping*/
    case 1:
      
      iota(num1,st);
      strcpy(ext,"_n_");
      strcat(ext,st);
      strcat(ext,".vtk");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
	      exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;     
/*...................................................................*/

/*.... vtk no-overlanping*/
    case 2:
      iota(num2,st);
      strcpy(ext,"_n_");
      strcat(ext,st);
      strcat(ext,"_part_");
      iota(num1,st);
      strcat(ext,st);
      strcat(ext,".vtk");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
	      exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;     
/*...................................................................*/

/*...*/
    case 3:
      iota(num1,st);
      strcpy(ext,"_n_");
      strcat(ext,st);
      strcat(ext,"_part.dat");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*...*/
    case 4:
      iota(num1,st);
      strcpy(ext,"_n_");
      strcat(ext,st);
      strcat(ext,"_part.vtk");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
	       exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*...*/
    case 6:
      strcat(ext,"_pgeo.vtk");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... arquivo com o log de excucao*/
    case 7:
      iota(num1,st);
      strcpy(ext,"_n_");
      strcat(ext,st);
      strcat(ext,"_log_");
      iota(num2,st);
      strcat(ext,st);
      strcat(ext,".txt");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... resultados da temperatura*/
    case 8:
      iota(num1,st);
      strcpy(ext,"_D1_step_");
      strcat(ext,st);
      strcat(ext,".vtk");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... arquivo com o log do solv*/
    case 10:
      iota(num1,st);
      strcpy(ext,"_n_");
      strcat(ext,st);
      strcat(ext,"_it_log.txt");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... arquivo com o log do solv*/
    case 11:
      iota(num1,st);
      strcpy(ext,"_n_");
      strcat(ext,st);
      strcat(ext,"_solv_log.txt");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... arquivo da matriz no formato coo*/
    case 12:
      strcat(ext,".mtx");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... arquivo da matriz no formato coo binario*/
    case 13:
      strcat(ext,"_bin.mtx");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... arquivo da matriz no formato coo*/
    case 14:
      strcat(ext,"_b.mtx");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... arquivo do vtk com iteraçoes intenas*/
    case 15:
      strcpy(ext,"_D1_step_");
      iota(num1,st);
      strcat(ext,st);
      strcat(ext,"_it_");
      iota(num2,st);
      strcat(ext,st);
      strcat(ext,".vtk");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... arquivo csv */                        
    case 16:
      iota(num1,st);
      strcat(ext,st);
      strcat(ext,".csv");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... arquivo vtk */                        
    case 17:
      strcat(ext,"_face.vtk");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... arquivo map */                        
    case 18:
      iota(num2,st);
      strcat(ext,"_map_");
      strcat(ext,st);
      strcat(ext,".dat");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... arquivo do vtk com iteraçoes intenas*/
    case 19:
      strcpy(ext,"_T1_step_");
      iota(num1,st);
      strcat(ext,st);
      strcat(ext,"_it_");
      iota(num2,st);
      strcat(ext,st);
      strcat(ext,".vtk");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... resultados Transporte 1*/
    case 20:
      iota(num1,st);
      strcpy(ext,"_T1_step_");
      strcat(ext,st);
      strcat(ext,".vtk");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*... medias do tempos MPI*/
    case 60:
      strcpy(ext,"_mTime_");
      iota(num1,st);
      strcat(ext,st);
      strcat(ext,".txt");
      size1 = strlen(name);
      size2 = strlen(ext);
      if( (size1+size2)  > SIZEMAX){
        fprintf(stderr,"Nome do arquivo muito extenso.\n"
	               "name : \"%s\"\n"
		       "Name maximo : %d\n"
		       "Funcao %s, arquivo fonte \"%s\"\n" 
		       ,name,SIZEMAX,__func__,__FILE__);
        exit(EXIT_FAILURE);      
      }
      strcpy(*out,name);
      strcat(*out,ext);
      break;
/*...................................................................*/

/*...*/
    default:
      ERRO_OP(__FILE__,__func__,cod);
    break;
/*...................................................................*/
  }
/*===================================================================*/
}
/*********************************************************************/

/*********************************************************************
 * ito:pegar o valor de um interio e trasnforma em uma cadeia de     *
 * caracter.                                                         *
 *********************************************************************
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*
 * istep -> posso no tempo                                           *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * st -> string com a volar de istep                                 *
 * ------------------------------------------------------------------*
 *********************************************************************/
void iota( INT t , char* st ){
/*...*/
  char c;
  int n , dec , j , k;
/*...................................................................*/
/**/
/*...*/
  n = ( t > 0 )? (int)( log10((double)t) + 1 ): 1;
  st[n] = '\0'; 
/*...................................................................*/
/**/
/*...*/
  dec = 10 ;   
  for( j = n - 1 ; j >= 0 ; --j ){
    k = ( t % dec ) / ( dec / 10 );
    c = (char)(k + '0');
    st[j] = c;
    dec = dec * 10;
  }
/*...................................................................*/
}
/*********************************************************************/ 


/*********************************************************************
 * READMACRO: le uma carcteres de um arquivo                         *
 *********************************************************************
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*
 * mc     -> string de caracter                                      *
 * allline-> leitura da linha toda                                   *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * mc -> string com lida                                             *
 * ------------------------------------------------------------------*
 *********************************************************************/
void readMacro(FILE* file,char *mc,bool allline)
{
 static char line[LINE_SIZE];
 static char word[WORD_SIZE];
 char c;
 int z;

#ifdef _DEBUG_READ 
 printf("DEGUB INFO: Antes funcao readmacro.\n"
        "linha  : \"%s\"\n",line);
#endif

 /*lendo um linha inteira*/

 if(allline){
  z = 0;
  while((c=getc(file))!='\n' && z<WORD_SIZE-1 && c!=EOF){
    line[z++] = c;
  }
  line[z]='\0';
  strcpy(mc,line);
 }
 
 /*lendo apenas uma paralvra*/
 else{
   fscanf(file,"%s",word);
   sscanf(word,"%s",mc);
 }   
 
/*... se encontrado '//' discarta a linha lida*/
 if(mc[0] =='/' && mc[1] == '/'){
   mc[0] = '\0'; 
   while((c=getc(file))!='\n'&& c!=EOF);
//   printf("linha discartada\n");
 } 

#ifdef _DEBUG_READ 
 printf("DEGUB INFO: Depois funcao readmacro.\n"
        "linha  : \"%s\"\n"
        "palavra: \"%s\"\n",line,mc);
#endif
}
/*********************************************************************/

int rl(FILE *f,char *st){
  
  long posicao;
  char s[LINE_SIZE];
  
  posicao=ftell(f);
  clearLine(s);
  readMacro(f,s,false);

#ifdef _DEBUG_READ 
   fprintf(stderr,"\nLinha lida\nlinha=!%s!",s);
   fprintf(stderr,"\nMACRO procurada   !%s!",st);
#endif
/*...*/  
  if(!strcmp(s,st)){
/*-------------------------------------------------------------------*/
#ifdef _DEBUG_READ 
      fprintf(stderr,"\nMacro !%s! achada\n",st); 
#endif
/*-------------------------------------------------------------------*/
    return 0;
  }
/*====== =============================================================*/
 /*...*/ 
  else{
/*-------------------------------------------------------------------*/
#ifdef _DEBUG_READ 
     fprintf(stderr,"\nRetonando linha");
#endif
/*-------------------------------------------------------------------*/
     fseek(f,posicao,0);
/*-------------------------------------------------------------------*/
#ifdef _DEBUG_READ 
     fprintf(stderr,"\nLinha retrocedida");
#endif
/*-------------------------------------------------------------------*/
  return 1;
  }
/*===================================================================*/
}
/*********************************************************************/
void clearLine(char *s){
   unsigned int i;
   i=0;
/*-------------------------------------------------------------------*/
#ifdef _DEBUG_READ 
           fprintf(stderr,"Funcao clear_line\n");
#endif
/*-------------------------------------------------------------------*/
   while(i<WORD_SIZE-1)
       s[i++]='\0';         
}

