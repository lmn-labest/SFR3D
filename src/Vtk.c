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
   fprintf(f,"# vtk DataFile Version 2.0\n");
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
void writeVtkCoor(double *x,long nnode,int ndm,bool cod,FILE *f)
{
    long i,j;
    double dum;
/*===*/
   new_section(cod,f); 
   fprintf(f,"POINTS %ld double\n",nnode);
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
     for(i=0;i<nnode*3;i++)
       write_double(x[i],cod,f);
     new_section(cod,f); 
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
void writeVtkCell(int *el        ,short int *nen , short int *type
                 ,long numel     ,short int maxno,bool cod
                 ,FILE *f)
{
  long i,k,aux=0;
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
  fprintf(f,"CELLS %ld %ld\n",numel,aux);
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
  fprintf(f,"CELL_TYPES %ld\n",numel);
  for(i=0;i<numel;i++){
    if(type[i] == 2) 
      dum = VTK_TRIA;
    else if(type[i] == 3)
      dum = VTK_QUAD;
    else if(type[i] == 4)
      dum = VTK_TETR;
    else if(type[i] == 5)
      dum = VTK_HEXA;
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
 * writeVtkCellProp :  escreve propriedades por celulas               *
 * -------------------------------------------------------------------*
 * Parametro de entrada :                                             *
 * -------------------------------------------------------------------*
 * iprop -> propriedades (int)                                        *
 * dprop -> propriedades (int)                                        *
 * numel -> numero de elementos                                       *
 * gdl   -> graus de liberdade                                        *
 * s     -> nome do campo                                             *
 * cod1  -> true BINARY vtk false ACISS vtk                           *
 * cod2  -> 1 int; 2 double                                           *
 *  f    -> ponteiro para o arquivo de saida                          *
 * -------------------------------------------------------------------*
 * Parametro de saida :                                               *
 * -------------------------------------------------------------------*
 *********************************************************************/
void writeVtkCellProp(int *iprop,double *dprop,long int numel,int gdl
                     ,char *s   ,bool cod1    ,short int cod2,FILE *f)
{
/*===*/
   long i;
   int j;
/*====================================================================*/
/**/
/*===*/
   switch(cod2){
/*... int*/
     case 1:   
       fprintf(f,"SCALARS %s int %d\n",s,gdl); 
       fprintf(f,"LOOKUP_TABLE default\n"); 
       for(i=0;i<numel;i++){
         for(j=0;j<gdl;j++)
           write_int(iprop[i*gdl+j],cod1,f);
         new_section(cod1,f); 
       }  
       new_section(cod1,f);
       break; 
/*...................................................................*/

/*... double*/
     case 2:   
       fprintf(f,"SCALARS %s double %d\n",s,gdl); 
       fprintf(f,"LOOKUP_TABLE default\n"); 
       for(i=0;i<numel;i++){
         for(j=0;j<gdl;j++)
           write_double(dprop[i*gdl+j],cod1,f);
         new_section(cod1,f); 
       }  
       new_section(cod1,f);
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
   long i;
   int j;
   double ddum=0.0;
/*====================================================================*/
/**/
/*===*/

/*... escalar*/
   if( cod2 == 1){
     if( cod1 == 1 ){
       fprintf(f,"SCALARS %s int %d\n",s,ndf); 
       fprintf(f,"LOOKUP_TABLE default\n"); 
       for(i=0;i<nnode;i++){
         for(j=0;j<ndf;j++)
           write_int(iprop[i*ndf+j],cod,f);
         new_section(cod,f);
       }  
       new_section(cod,f);
     }  
     else if( cod1 == 2 ){
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
   else if(cod2 == 2){
     if( cod1 == 1 ){
     }
     if( cod1 == 2 ){
       fprintf(f,"VECTORS %s double \n",s); 
       for(i=0;i<nnode;i++){
         if(ndf == 2){
	  j=i*ndf;
          write_double(dprop[j],cod,f);
          write_double(dprop[j+1],cod,f);
          write_double(ddum,cod,f);
	 }
	 else if(ndf ==3)
           for(j=0;j<ndf;j++)
             write_double(dprop[i*ndf+j],cod,f);
        new_section(cod,f);
       }  
       new_section(cod,f);
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
        char str[128];
        sprintf(str, "%20.12e ", val);
        fprintf(f, str);
    }
}

/* ********************************************************************
 *  Function: force_big_endian
 *
 *  Purpose:
 *      Determines if the machine is little-endian.  If so, then
 *      , for binary data, it will force the data to be big-endian.
 *
 *  Note:       This assumes that all inputs are 4 bytes long.
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
          bytes[3] = bytes[5];
          bytes[5] = tmp;

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


