#include<Vtk.h>

static void write_int(int,bool,FILE*);
static void write_double(double,bool,FILE*);
static void force_big_endian(unsigned char *,bool,int);
static void new_section(bool,FILE*); 

/**********************************************************************
 * head_vtk : escreve o cabecalho do arquivo vtk                      *
 * -------------------------------------------------------------------*
 * Parametro de entrada :                                             *
 * -------------------------------------------------------------------*
 *  f    -> ponteiro para o arquivo de saida                          *
 *  s*   -> informacao do cabecalho do arquivo vtk                    *
 *  cod  -> true binary format false aciss                            *
 * -------------------------------------------------------------------*
 * Parametro de saida :                                               *
 * -------------------------------------------------------------------*
 * -------------------------------------------------------------------*
 *********************************************************************/
void headVtk(char *s,bool cod,FILE *f)
{
/*=== Saida debug*/
  #ifdef _DEBUG_VTK_
     fprintf(stderr,"Escrevendo cabecalho do arquivo vtk...\n");
  #endif   
/*===================================================================*/
/**/
/*=== Escrevendo no arquivo*/
   fprintf(f,"# vtk DataFile Version 3.0\n");
   fprintf(f,"%s\n",s);
   if(cod)
     fprintf(f,"BINARY\n");
   else  
     fprintf(f,"ASCII\n");
   fprintf(f,"DATASET UNSTRUCTURED_GRID\n");
/*===================================================================*/
/**/
/*=== Saida debug*/
  #ifdef _DEBUG_VTK_
     fprintf(stderr,"Cabecalho escrito...\n");
  #endif   
/*===================================================================*/
}
/*********************************************************************/

/**********************************************************************
 * Data de criacao    : 00/00/0000                                    *
 * Data de modificaco : 15/02/2017                                    * 
 * -------------------------------------------------------------------*
 * timeVtk : informacao temporal do arquivo vtk                       *
 * -------------------------------------------------------------------*
 * Parametro de entrada :                                             *
 * -------------------------------------------------------------------*
 *  t    -> tempo                                                     *
 * iStep -> passo de tempo                                            *
 *  f    -> ponteiro para o arquivo de saida                          *
 * -------------------------------------------------------------------*
 * Parametro de saida :                                               *
 * -------------------------------------------------------------------*
 * -------------------------------------------------------------------*
 * OBS:                                                               *
 * -------------------------------------------------------------------*
 *********************************************************************/
void timeVtk(double t,int iStep,bool cod,FILE *f){

   new_section(cod,f);
   fprintf(f,"FIELD FieldData 3\n");

/*... passo de tempo*/
   fprintf(f,"CYCLE 1 1 int\n");
   write_int(iStep,cod,f);
   new_section(cod,f);

/*... tempo em segundos*/   
   fprintf(f,"TIME_S 1 1 double\n");
   write_double(t,cod,f);
   new_section(cod,f);
/*... tempo em horas*/ 
   fprintf(f,"TIME_H 1 1 double\n");
   write_double(t/3600.e0,cod,f);
   new_section(cod,f);
   
} 
/*********************************************************************/

/**********************************************************************
 * write_vtk_coor : escreve as coordenadas(arquivo vtk)               *
 * -------------------------------------------------------------------*
 * Parametro de entrada :                                             *
 * -------------------------------------------------------------------*
 *  f -> ponteiro para o arquivo de saida                             *
 *  x -> coordenadas x(x1,y1,z1,x2,y2,z2,...)                         *
 *  cod  -> true binary format false aciss                            *
 *  nnode - numero de nos                                             *
 *  ndm dimencao                                                      *
 * -------------------------------------------------------------------*
 * Parametro de saida :                                               *
 * -------------------------------------------------------------------*
 * -------------------------------------------------------------------*
 *********************************************************************/
void writeVtkCoor(double *x,INT nnode,int ndm,bool cod,FILE *f)
{
    INT i,j;
    double dum;
/*===*/
   new_section(cod,f); 
   fprintf(f,"POINTS %ld double\n",(long) nnode);
/*... 1D(x1,0,0,x2,0,0,...)*/  
   if(ndm == 1){
     j = 0;
     for(i=0;i<nnode;i++){
       write_double(x[j++],cod,f);
       dum = 0.0;
       write_double(dum,cod,f);
       write_double(dum,cod,f);
       new_section(cod,f); 
     }
   }
/*... 2D(x1,y1,0,x2,y2,0,...)*/  
   else if(ndm == 2){
     j = 0;
     for(i=0;i<nnode;i++){
       write_double(x[j++],cod,f);
       write_double(x[j++],cod,f);
       dum = 0.0;
       write_double(dum,cod,f);
       new_section(cod,f); 
     }
   }  
/*... 3D(x1,y1,z1,x2,y2,z2,...)*/  
   else{
     j = 0; 
     for(i=0;i<nnode;i++){
       write_double(x[j++],cod,f);
       write_double(x[j++],cod,f);
       write_double(x[j++],cod,f);
       new_section(cod,f);
     } 
   }
   new_section(cod,f); 
/*===================================================================*/
}
/**********************************************************************/

/**********************************************************************
 * writeVtkCell : escreve elementos no formato vtk                    *
 * -------------------------------------------------------------------*
 * Parametro de entrada :                                             *
 * -------------------------------------------------------------------*
 * el    -> conectividade                                             *
 * nen   -> numero de nos por elemento                                *
 * type  -> tipo geometrico de elemento                               *
 * numel -> tipo numerico elementos                                   *
 * maxno -> numero maximo de nos por elemento                         *
 * cod   -> true binary format false aciss                            *
 * f     -> ponteiro para o arquivo de saida                          *
 * -------------------------------------------------------------------*
 * Parametro de saida :                                               *
 * -------------------------------------------------------------------*
 *********************************************************************/
void writeVtkCell(int *el       ,short *nen ,short *type
                 ,INT numel     ,short maxno,bool cod
                 ,FILE *f)
{
  INT i,k;
  long aux=0;
  int dum;
  int j;

/*=== saida para a usuario*/
  #ifdef _DEBUG_VTK_
    fprintf(stderr,"Escrevendo elementos...\n");
  #endif  
/*===================================================================*/

  for(i=0;i<numel;i++)
    aux += nen[i] + 1;
/*...*/
  new_section(cod,f); 
  fprintf(f,"CELLS %ld %ld\n",(long) numel,aux);
  for(i=0;i<numel;i++){
    dum = nen[i];
    write_int(dum,cod,f);
    for(j=0;j<dum;j++){
      k = i*maxno + j;
      write_int(el[k],cod,f);
    }
  new_section(cod,f); 
  }	
/*...................................................................*/       

/*...*/
  new_section(cod,f); 
  fprintf(f,"CELL_TYPES %ld\n",(long) numel);
  for(i=0;i<numel;i++){
    if(type[i] == LINECELL) 
      dum = VTK_LINE;
    else if(type[i] == TRIACELL) 
      dum = VTK_TRIA;
    else if(type[i] == QUADCELL)
      dum = VTK_QUAD;
    else if(type[i] == TETRCELL)
      dum = VTK_TETR;
    else if(type[i] == HEXACELL)
      dum = VTK_HEXA;
    else if (type[i] == PIRACELL)
      dum = VTK_PIRA;
    else if (type[i] == DOTCELL)
      dum = VTK_DOT;
    else
      dum = 0;
    write_int(dum,cod,f);
    new_section(cod,f); 
  }	
  new_section(cod,f); 
/*...................................................................*/       

/*=== saida para o usuario*/
  #ifdef _DEBUG_VTK_
   fprintf(stderr,"Elementos escritos.\n");
  #endif 
/*===================================================================*/
}
/*********************************************************************/

/**********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 25/11/2017                                   * 
 *-------------------------------------------------------------------* 
 * writeVtkProp :  escreve propriedades                               *
 * -------------------------------------------------------------------*
 * Parametro de entrada :                                             *
 * -------------------------------------------------------------------*
 * iprop -> propriedades (int)                                        *
 * dprop -> propriedades (double)                                     *
 * n     -> numero de linhas                                          *
 * gdl   -> graus de liberdade                                        *
 * s     -> nome do campo                                             *
 * cod1  -> true BINARY vtk false ASCII vtk                           *
 * cod2  -> 1 int; 2 double                                           *
 * cod3  -> 1 SCALAR ; 2 VECTORS; 3 TENSOR                            *
 *  f    -> ponteiro para o arquivo de saida                          *
 * -------------------------------------------------------------------*
 * Parametro de saida :                                               *
 * -------------------------------------------------------------------*
 *********************************************************************/
void writeVtkProp(int *iprop,double *dprop,INT n     ,int gdl
                 ,char *s   ,bool cod1    ,short cod2,short cod3
                 ,FILE *f)
{
/*===*/
   INT i;
   int j;
/*====================================================================*/
/**/
/*===*/
   switch(cod2){
/*... int*/
     case INTEGER_VTK:   
       fprintf(f,"SCALARS %s int %d\n",s,gdl); 
       fprintf(f,"LOOKUP_TABLE default\n"); 
       for(i=0;i<n;i++){
         for(j=0;j<gdl;j++)
           write_int(iprop[i*gdl+j],cod1,f);
         new_section(cod1,f); 
       }  
       new_section(cod1,f);
       break; 
/*...................................................................*/

/*... double*/
     case DOUBLE_VTK:
/*...escalar*/
       if(cod3 == SCALARS_VTK){     
         fprintf(f,"SCALARS %s double %d\n",s,gdl); 
         fprintf(f,"LOOKUP_TABLE default\n"); 
         for(i=0;i<n;i++){
           for(j=0;j<gdl;j++)
             write_double(dprop[i*gdl+j],cod1,f);
           new_section(cod1,f); 
         }  
         new_section(cod1,f);
       }
/*...................................................................*/

/*... vetorial*/
       else if(cod3 == VECTORS_VTK){     
         fprintf(f,"VECTORS %s double\n",s); 
         for(i=0;i<n;i++){
/*...*/
           if(gdl == 2){ 
             for(j=0;j<gdl;j++)
               write_double(dprop[i*gdl+j],cod1,f);            
             write_double(0.0e0,cod1,f);
            }
/*...................................................................*/

/*...*/
           else
             for(j=0;j<gdl;j++)
               write_double(dprop[i*gdl+j],cod1,f);
/*...................................................................*/
           new_section(cod1,f); 
         }
/*...................................................................*/
       }
/*...................................................................*/

/*... tensor*/
       else if(cod3 == TENSORS_VTK){     
         fprintf(f,"TENSORS %s double\n",s); 
         for(i=0;i<n;i++){
           for(j=0;j<gdl*gdl;j++)
             write_double(dprop[i*gdl*gdl+j],cod1,f);
/*...................................................................*/
           new_section(cod1,f); 
         }   
         new_section(cod1,f);
       }
       break;
/*...................................................................*/
   } 
/*====================================================================*/
}
/*********************************************************************/
/**********************************************************************
 * write_vtk_node_prop: escreve propriedades nodais                   *
 * -------------------------------------------------------------------*
 * Parametro de entrada :                                             *
 * -------------------------------------------------------------------*
 * iprop-> propriedades (int     )                                    *
 * dprop-> propriedades (double  )                                    *
 * cod1-> 1-inteiro                                                   *
 *        2-double                                                    *
 * cod2-> 1-escalar                                                   *
 *        2-vetorial                                                  *
 * ndf -> graus de liberdade                                          *
 *   s -> campo a ser escrito                                         *
 * cod -> true BINARY vtk false ACISS vtk                             *
 *  f   -> ponteiro para o arquivo de saida                           *
 * -------------------------------------------------------------------*
 * Parametro de saida :                                               *
 * -------------------------------------------------------------------*
 * -------------------------------------------------------------------*
 *********************************************************************/
void writeVtkNodeProp(int *iprop,double *dprop,short cod1,short cod2
                        ,int nnode,short ndf,char *s,bool cod,FILE *f)
{
/*===*/
   INT i;
   int j;
   double ddum=0.0;
/*====================================================================*/
/**/
/*===*/

/*... escalar*/
   if( cod2 == INTEGER_VTK){
     if( cod1 == SCALARS_VTK ){
       fprintf(f,"SCALARS %s int %d\n",s,ndf); 
       fprintf(f,"LOOKUP_TABLE default\n"); 
       for(i=0;i<nnode;i++){
         for(j=0;j<ndf;j++)
           write_int(iprop[i*ndf+j],cod,f);
         new_section(cod,f);
       }  
       new_section(cod,f);
     }  
     else if( cod1 == DOUBLE_VTK ){
       fprintf(f,"SCALARS %s double %d\n",s,ndf); 
       fprintf(f,"LOOKUP_TABLE default\n"); 
       for(i=0;i<nnode;i++){
         for(j=0;j<ndf;j++)
           write_double(dprop[i*ndf+j],cod,f);
         new_section(cod,f);
       }  
       new_section(cod,f);
     }
   }
/*...................................................................*/

/*... vetorial*/
   else if(cod2 == VECTORS_VTK){
     if( cod1 == INTEGER_VTK ){
     }
     if( cod1 == DOUBLE_VTK ){
       fprintf(f,"VECTORS %s double \n",s); 
       for(i=0;i<nnode;i++){
         if(ndf == 2){
	         j=i*ndf;
           write_double(dprop[j],cod,f);
           write_double(dprop[j+1],cod,f);
           write_double(ddum,cod,f);
	      }
	      else if(ndf ==3){
          for(j=0;j<ndf;j++)
            write_double(dprop[i*ndf+j],cod,f);
          new_section(cod,f);
        }  
        new_section(cod,f);
       }
     }
   }
/*...................................................................*/

/*====================================================================*/
}
/*********************************************************************/

/* *******************************************************************
 *  Function: write_double
 *
 *  Purpose:
 *                                                                     
 *  Programmer:              
 *  Creation:                    
 * 
 *  Modifications:
 *
 * *******************************************************************/

static void write_double(double val,bool cod,FILE *f)
{
    if (cod)
    {
        force_big_endian((unsigned char *) &val,cod,8);
        fwrite(&val, sizeof(double), 1, f);
    }
    else
    {
/*      char str[128];
        sprintf(str, "%20.12e ", val);
        fprintf(f, str);*/
        fprintf(f," %.15e ", val);
    }
}

/* ********************************************************************
 *  Function: force_big_endian
 *
 *  Purpose:
 *      Determines if the machine is little-endian.  If so, then
 *      , for binary data, it will force the data to be big-endian.
 *
 *  Note:       This assumes that all inputs are 4 bytes INT.
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 * 
 * *******************************************************************/

static void force_big_endian(unsigned char *bytes,bool cod,int nbytes)
{
    static int doneTest = 0;
    static int shouldSwap = 0;
    if (!doneTest)
    {
        int tmp1 = 1;
        unsigned char *tmp2 = (unsigned char *) &tmp1;
        if (*tmp2 != 0)
            shouldSwap = 1;
        doneTest = 1;
    }

    if (shouldSwap & cod)
    {
        
	if(nbytes==8){
/*          fprintf(stderr,"8 bytes\n");*/
          unsigned char tmp = bytes[0];
          bytes[0] = bytes[7];
          bytes[7] = tmp;
          tmp = bytes[1];
          bytes[1] = bytes[6];
          bytes[6] = tmp;
          tmp = bytes[2];
          bytes[2] = bytes[5];
          bytes[5] = tmp;
          tmp = bytes[3];
          bytes[3] = bytes[4];
          bytes[4] = tmp;

	}
        else if(nbytes==4){	
/*	  fprintf(stderr,"4 bytes\n");*/
          unsigned char tmp = bytes[0];
          bytes[0] = bytes[3];
          bytes[3] = tmp;
          tmp = bytes[1];
          bytes[1] = bytes[2];
          bytes[2] = tmp;
       }
    }
}

/* *******************************************************************
 *  Function: write_int
 *
 *  Purpose:
 *      Writes an integer to the currently open file.  This routine
 *      takes care of ASCII vs binary issues.
 *
 *  Programmer: Hank Childs
 *  Creation:   September 3, 2004
 * 
 * *******************************************************************/

static void write_int(int val,bool cod,FILE *f)
{
    if (cod)
    {
        force_big_endian((unsigned char *) &val,cod,4);
        fwrite(&val, sizeof(int), 1, f);
    }
    else
    {
        char str[128];
        sprintf(str, "%d ", val);
        fprintf(f, str);

    }
}


static void new_section(bool cod,FILE *f){
 
  if(cod);
  else
    fprintf(f,"\n");

}


