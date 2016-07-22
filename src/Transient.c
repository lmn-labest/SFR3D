#include<Transient.h>  
/********************************************************************* 
 * SETTRANSIENTSCHEME : set a discretizacao temporal                 * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * word -> str com a discretizacao                                   * 
 * type -> nao definido                                              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * type -> tipo de discretizacao                                     * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void setTransientScheme(char *word,short *type){

  if(!strcmp(word,"EULER"))
   *type = EULER;
  else if(!strcmp(word,"BACKWARD"))
   *type =  BACKWARD;

}
/*********************************************************************/


