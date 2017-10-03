#include<Properties.h>

/*********************************************************************
 * Data de criacao    : 29/08/2017                                   *
 * Data de modificaco : 15/09/2017                                   *
 *-------------------------------------------------------------------*
 * IDEALGASDENSITY: kg/(m^3)                                         *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * t - temperatura em kelvin                                         *
 * presRef - pressao de referencia ou termomecanica                  *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE airDensity(DOUBLE const t,DOUBLE const presRef
                 ,bool const fKelvin) {

  DOUBLE tc,d;
  
  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  d = (MMOLARAR*presRef)/(IDEALGASR*tc);

  return d;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 29/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * AIRSPECIFIHEAT: kJ/(kg.K)                                         *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * t - temperatura em Kelvin                                         *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * http://www.engineeringtoolbox.com/air-properties-d_156.html   	   *
 *                                                                   *
 * regressao com polinomio de ordem 5 obtido pelo excel              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE airSpecifiHeat(DOUBLE const t,bool const fKelvin) {

  short i,n=sHeat.nPol;  
  DOUBLE a[10],y,d;
  DOUBLE tc;

  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  for (i = 0; i < n; i++)
    a[i] = sHeat.a[i];
  
/*... polinomio*/
  y = a[0];
  for (i = 1; i < n; i++)
    y += a[i]*pow(tc,i);
/*.....................................................................*/

/*...*/
  d = 1.e+0;
/*.....................................................................*/

 if (y < 0) {
    printf("Calor especifico negativo!!"
           "Calor especifico = %e\n"
           "Temperatura      = %lf\n!!",d*y,tc);
    exit(EXIT_FAILURE);
  }


  return d*y;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 29/08/2017                                   *
 * Data de modificaco : 16/09/2017                                   *
 *-------------------------------------------------------------------*
 * AIRDYNAMICVISCOSITY: kJ/(kg.K)                                    *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * t - temperatura em Kelvin                                         *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * http://www.engineeringtoolbox.com/air-properties-d_156.html       *
 *                                                                   *
 * regressao com polinomio de ordem 5 obtido pelo excel              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE airDynamicViscosity(DOUBLE const t,bool const fKelvin) {

  short i,n=dVisc.nPol;  
  DOUBLE a[6],x[5],y,d;
  DOUBLE tc;

  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  switch (dVisc.type) {
/*... polinomio*/
    case POL:
      for (i = 0; i < n; i++)
        a[i] = dVisc.a[i];
/*... t */  
        x[0] = tc;
/*... t^2 */
        x[1] = tc*tc;
/*... t^3 */
        x[2] = tc*x[1];
/*... t^4 */
        x[3] = tc*x[2];
/*... t^5 */
        x[4] = tc*x[3];
  
/*... polinomio*/
        y = a[0] + a[1]*x[0] + a[2]*x[1] 
          + a[3]*x[2] + a[4]*x[3] + a[5]*x[4];
        d = 1.e-05;
/*.....................................................................*/
        break;
/*.....................................................................*/

/*... polinomio*/
    case SUTHERLAND:
      a[0] = dVisc.a[0]; /*viscosidade de referencia*/
      a[1] = dVisc.a[1]; /*temperatura de referencia*/ 
      a[2] = dVisc.a[2]; /*constante de Sutherland*/

      x[0] = a[1]+ a[2];
      x[1] = tc  + a[2];

      y = a[0]*(x[0]/x[1])*pow(tc/a[1],1.5);

      d = 1.e0;
/*.....................................................................*/
      break;
/*.....................................................................*/

/*...*/
    default:  
      ERRO_OP(__FILE__,__func__,dVisc.type);
      break;
/*.....................................................................*/
  }
/*.....................................................................*/
 
  if (y < 0) {
    printf("Viscosidade dinamica negativa!!\n"
           "Viscosidade dinamica = %e\n"
           "Temperatura          = %lf\n!!",d*y,tc);
    exit(EXIT_FAILURE);
  }

  return d*y;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 28/08/2017                                   *
 * Data de modificaco : 19/09/2017                                   *
 *-------------------------------------------------------------------*
 * AIRTHERMALCONDUCTITY: [KW/m.K]                                    *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * t - temperatura em kelvin                                         *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * http://www.engineeringtoolbox.com/air-properties-d_156.html   	   *
 *                                                                   *
 * regressao com polinomio de ordem 5 obtido pelo excel              *
 *                                                                   *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE airThermalConductvity(DOUBLE const t,bool const fKelvin) {

  short i,n=thCond.nPol;  
  DOUBLE a[6],x[5],y,d;
  DOUBLE tc;

  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  
  switch (thCond.type) {
/*... polinomio*/
    case POL:
       for (i = 0; i < n; i++)
        a[i] = thCond.a[i];
/*... t */  
        x[0] = tc;
/*... t^2 */
        x[1] = tc*tc;
/*... t^3 */
        x[2] = tc*x[1];
/*... t^4 */
        x[3] = tc*x[2];
/*... t^5 */
        x[4] = tc*x[3];
  
/*... polinomio*/
        y = a[0] + a[1]*x[0] + a[2]*x[1] 
          + a[3]*x[2] + a[4]*x[3] + a[5]*x[4];
        d = 1.e-05;
/*.....................................................................*/
        break;
/*.....................................................................*/
 /*...*/
    default:  
      ERRO_OP(__FILE__,__func__,dVisc.type);
      break;
/*.....................................................................*/
  }
/*.....................................................................*/

/*...*/ 
  if (y < 0) {
    printf("Condutividade termica negativa!!\n"
           "Condutividade termica = %e\n"
           "Temperatura           = %lf\n!!",d*y,tc);
    exit(EXIT_FAILURE);
  }
/*.....................................................................*/

  return d*y;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 28/08/2017                                   *
 * Data de modificaco : 24/09/2017                                   *
 *-------------------------------------------------------------------*
 * TEMPFORSPECIFICENTHALPY: calcula a entalpia espeficia apartir da  *
 * temperatura                                                       *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * t        - temperatura (°C/K)                                     *
 * sHeatRef - calor especifico de referencia constante com temp      *
 * fSheat   - calor especifico com variacao com a Temperatura        *
 * fKelvin  - temperatura dada em kelvin                             *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * entalpia sensivel                                                 *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE tempForSpecificEnthalpy(DOUBLE const t   , DOUBLE const sHeatRef
                             , bool const fSheat, bool const fKelvin) {

  short i,n=sHeat.nPol;
  DOUBLE a[6],d,dt,tmp;
  DOUBLE tc,tRef= TREF ;

  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  if(fSheat){
    for (i = 0; i < n; i++)
      a[i] = sHeat.a[i];

    tmp = 0.0;
    for (i = 0; i < n; i++) {
      d    = (double) (i + 1);
      dt   = tc-tRef;
      tmp += a[i]*pow(dt,d)/d;
    }
  }

  else 
    tmp = TEMP_FOR_ENTHALPY(sHeatRef,tc,TREF);

  return tmp;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 28/08/2017                                   *
 * Data de modificaco : 24/09/2017                                   *
 *-------------------------------------------------------------------*
 * SPECIFICENTHALPYFORTEMP:  calcula a temperatura apartir da        *
 * entalpia especifica                                               * 
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * hs       - entalpia sensivel                                      *
 * sHeatRef - calor especifico de referencia constante com temp      *
 * fSheat   - calor especifico com variacao com a Temperatura        *
 * fKelvin  - temperatura dada em kelvin                             *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * temperatura (°C/Kelvin)                                          *
 *-------------------------------------------------------------------*
 * OBS:                                                              * 
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE specificEnthalpyForTemp(DOUBLE const hs  , DOUBLE const sHeatRef
                             , bool const fSheat, bool const fKelvin) 
{
  unsigned short i;
  bool flag = false;
  DOUBLE f,fl,t,conv,tol=1e-06;
 
/*...*/
  if(fSheat){
/*... chute inicial usando a massa espeficia constante*/
    t = ENTHALPY_FOR_TEMP(sHeatRef,hs,TREF);
/*...*/
    conv = (hs-tempForSpecificEnthalpy(t,sHeatRef,fSheat,true))*tol;
    conv = fabs(conv);
/*... Newton-Raphson*/
    for(i=0;i<60000;i++){
      f  = hs-tempForSpecificEnthalpy(t,sHeatRef,fSheat,true);
      if(fabs(f) < conv) {
        flag = true;
        break;
      }
    
      fl = airSpecifiHeat(t,true);
      t += f/fl;   
    }
/*...................................................................*/

    if(!flag){
      printf("%i %e %e %e\n",i,t,f,conv);
      ERRO_GERAL(__FILE__,__func__,__LINE__,
      "sEnthalpy->temperature:\n Newton-raphson did not converge !!");
    }
  }
/*...................................................................*/

/*...*/
  else
    t = ENTHALPY_FOR_TEMP(sHeatRef,hs,TREF);
/*...................................................................*/

  if(!fKelvin)
    t = KELVIN_FOR_CELSIUS(t);  

  return t;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 27/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * WATERDESNITY: kg/(m^3)                                            *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * t - temperatura em °C                                             *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE waterDensity(DOUBLE const t) {

  DOUBLE a,b,c,d;
  DOUBLE tmp1,tmp2,tmp3,tmp4;

  a = 288.9414e0;
  b = 508929.2e0;
  c = 68.12963e0;
  d = 3.9863e0;
  
  tmp1 = t + a;
  tmp2 = b*(t+c);
  tmp3 = (t-d)*(t-d);  

  tmp4 =  1.e0 - tmp1/tmp2*tmp3;
  
  if (tmp4 < 0) {
    printf("Massa especifica negativa\n!!");
    exit(EXIT_FAILURE);
  }

  return 1000.e0*tmp4;


}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 27/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * WATERSPECIFIHEAT: kJ/(kg.°C)                                      *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * t - temperatura em °C                                             *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * www.engineeringtoolbox.com/water-thermal-properties-d_162.html	   *
 *                                                                   *
 * regressao com polinomio de ordem 5 obtido pelo excel              *
 *                                                                   *
 * range  T=]0,100[                                                  * 
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE waterSpecifiHeat(DOUBLE const t) {

  DOUBLE a[5],x[4],y,d;

  a[0] = 4.21470622853361e+00;
  a[1] =-2.87563217630463e-03;
  a[2] = 7.57667401959410e-05;
  a[3] =-7.93659500485961e-07;
  a[4] = 3.28899958033620e-09;

/*... t */  
  x[0] = t;
/*... t*t */
  x[1] = t*t;
/*... t*t*t */
  x[2] = t*x[1];
/*... t*t*t*t */
  x[3] = t*x[2];
  
/*... polinomio*/
  y = a[0] + a[1]*x[0] + a[2]*x[1] + a[3]*x[2] + a[4]*x[3];
/*.....................................................................*/

  if (y < 0) {
    printf("Calor especifico negativo\n!!");
    exit(EXIT_FAILURE);
  }

/*...*/
  y = d = 1.e0;
/*.....................................................................*/

  return d*y;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 27/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * WATERDINAMICVICOSITY: Kg/(m.s)                                    *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * t - temperatura em °C                                             *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * www.engineeringtoolbox.com/water-thermal-properties-d_162.html	   *
 *                                                                   *
 * regressao com polinomio de ordem 5 obtido pelo excel              *
 *                                                                   *
 * range  T=]0,100[                                                  * 
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE waterDynamicViscosity(DOUBLE const t) {

  DOUBLE a[5],x[4],y,d;

  a[0] = 1.75373877709050e+00;
  a[1] =-5.20548110460391e-02;
  a[2] = 8.73350059860241e-04;
  a[3] =-7.60339598706936e-06;
  a[4] = 2.96188559332536E-16;

/*... t */  
  x[0] = t;
/*... t*t */
  x[1] = t*t;
/*... t*t*t */
  x[2] = t*x[1];
/*... t*t*t*t */
  x[3] = t*x[2];
  
/*... polinomio*/
  y = a[0] + a[1]*x[0] + a[2]*x[1] + a[3]*x[2] + a[4]*x[3];
/*.....................................................................*/
  
   if (y < 0) {
    printf("Visosidade dinamica negativa\n!!");
    exit(EXIT_FAILURE);
  }

/*...*/
  d = 1.e-3;
/*.....................................................................*/

  return d*y;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 28/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * WATERTHERMALCONDUCTITY: [W/m.K]                                   *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * t - temperatura em °C                                             *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 * www.engineeringtoolbox.com/water-thermal-properties-d_162.html	   *
 *                                                                   *
 * regressao com polinomio de ordem 5 obtido pelo excel              *
 *                                                                   *
 * range  T=]0,100[                                                  * 
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE waterThermalConductvity(DOUBLE const t) {

  DOUBLE a[5],x[4],y,d;

  a[0] = 0.5706e+00;
  a[1] = 1.756e-03;
  a[2] = 6.46e-06;

/*... t */  
  x[0] = t;
/*... t*t */
  x[1] = t*t;
  
/*... polinomio*/
  y = a[0] + a[1]*x[0] + a[2]*x[1];
/*.....................................................................*/
  
   if (y < 0) {
    printf("Condutividade termica negativa\n!!");
    exit(EXIT_FAILURE);
  }

/*...*/
  d = 1.e-3;
/*.....................................................................*/

  return d*y;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 27/08/2017                                   *
 * Data de modificaco : 18/09/2017                                   *
 *-------------------------------------------------------------------*
 * UPDATEDENSITY:                                                    *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void updateDensity(DOUBLE *RESTRICT temp , DOUBLE *RESTRICT density
                 , DOUBLE const alpha    , bool const iKelvin    
                 , INT const nEl         , char  const iCod)

{
  short nD = DENSITY_LEVEL;
  INT i;
  DOUBLE den,den0;
/*...*/
  switch (iCod){
    case PROP_UPDATE_SIMPLE_LOOP:
      for(i=0;i<nEl;i++){
        den0 =  MAT2D(i,2 ,density ,nD);         
        den = airDensity(temp[i], thDynamic.pTh[2], iKelvin);
/*...*/           
        MAT2D(i,2 ,density ,nD) =  alpha*den + (1.e0-alpha)*den0;
      }
/*..................................................................*/
    break;  

  case PROP_UPDATE_OLD_TIME:
    for(i=0;i<nEl;i++){
/*...t(n-2) = t(n-1)*/
      MAT2D(i,0 ,density ,nD) = MAT2D(i,1 ,density ,nD);
/*...t(n-1) = t(n)*/           
      MAT2D(i,1 ,density ,nD) = MAT2D(i,2 ,density ,nD);
    }
/*..................................................................*/
    break;
  }
/*..................................................................*/


}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 27/08/2017                                   *
 * Data de modificaco : 18/09/2017                                   *
 *-------------------------------------------------------------------*
 * UPDATESPECIFICHEAT:                                               *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void updateSpecificHeat(DOUBLE *RESTRICT temp, DOUBLE *RESTRICT sHeat
                       , bool const iKelvin 
                       ,INT const nEl        , char  const iCod)

{
  short nD = SHEAT_LEVEL;
  INT i;  
  
/*...*/
  switch (iCod){
    case PROP_UPDATE_SIMPLE_LOOP:
      for(i=0;i<nEl;i++)
/*...*/           
        MAT2D(i,2 ,sHeat ,nD) = airSpecifiHeat(temp[i],iKelvin);
/*..................................................................*/
    break;  

  case PROP_UPDATE_OLD_TIME:
    for(i=0;i<nEl;i++){
/*...*/
      MAT2D(i,0 ,sHeat ,nD) = MAT2D(i,1 ,sHeat ,nD);           
      MAT2D(i,1 ,sHeat ,nD) = MAT2D(i,2 ,sHeat ,nD);
    }
/*..................................................................*/
    break;
  }
/*..................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 27/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * UPDATESPECIFICHEAT:                                               *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void updateThermalCondutivty(DOUBLE *RESTRICT t,DOUBLE *RESTRICT thCond   
                            ,bool const iKelvin,INT const nEl)

{
  INT i;

  for(i=0;i<nEl;i++)         
    thCond[i] = airThermalConductvity(t[i],iKelvin);

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 01/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * UPDATESPECIFICHEAT:                                               *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void updateDynamicViscosity(DOUBLE *RESTRICT temp,DOUBLE *RESTRICT visc    
                          ,bool const iKelvin    ,INT const nEl)

{
  INT i;

  for(i=0;i<nEl;i++)         
    visc[i] = airDynamicViscosity(temp[i],iKelvin);

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 27/08/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * INITPROPTEMP: inicializao de propriedades com variacao temporal   *
 * dependentes da temperatura                                        *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * prop    -> nao definido                                           * 
 * t       -> temperatura                                            *
 * propMat -> propriedade de referencia por material                 * 
 * mat     -> material por celula                                    * 
 * np      -> numero niveis de tempos                                * 
 * nCell   -> numero de celulas                                      * 
 * iKelvin -> temperatura em kelvin (true|false)                     *
 * iProp   -> numero da propriedade                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * prop    -> propriedade iniciacializada                            * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
  *********************************************************************/
void initPropTemp(DOUBLE *RESTRICT prop   ,DOUBLE *RESTRICT t 
                 ,DOUBLE *RESTRICT propMat,short *RESTRICT mat
                 ,short const np          ,INT    const nCell 
                 ,bool const iKelvin      ,short const iProp)
{    
  INT i;
  unsigned short j,lMat;         
  for(i=0;i<nCell;i++){    

/*...*/
    lMat               = mat[i]-1;
/*...................................................................*/

/*...*/
    if( iProp == DENSITY )
      MAT2D(lMat, iProp, propMat, MAXPROP) 
        = airDensity(t[i], thDynamic.pTh[2], iKelvin);
/*...................................................................*/

/*...*/
    else if( iProp == SPECIFICHEATCAPACITYFLUID)  
      MAT2D(lMat, iProp, propMat, MAXPROP) 
      = airSpecifiHeat(t[i],iKelvin);
/*...................................................................*/

/*...*/
    else if( iProp == DYNAMICVISCOSITY)  
      MAT2D(lMat, iProp, propMat, MAXPROP) 
      = airDynamicViscosity(t[i],iKelvin);
/*...................................................................*/

/*...*/
    else if( iProp == THERMALCONDUCTIVITY)  
      MAT2D(lMat,iProp,propMat,MAXPROP) 
      = airThermalConductvity(t[i],iKelvin);
/*...................................................................*/

/*...*/
    for(j=0;j<np;j++)      
      MAT2D(i,j,prop,np) = MAT2D(lMat,iProp,propMat,MAXPROP); 
/*...................................................................*/
  
  }
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 04/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * INITSHEATPOL: inicializao a estrutura para o calculo dp calor     *
 * especifico em funcao da temperauta via polinomio                  *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * http://www.engineeringtoolbox.com/air-properties-d_156.html  	   *
 *                                                                   *
 * regressao com polinomio de ordem 5 obtido pelo excel              *
 *                                                                   *
 *********************************************************************/
void initSheatPol(void) {

  sHeat.nPol = 6;

  sHeat.a[0] = 1.048057862882000E+00;
  sHeat.a[1] =-4.283505419516730E-04;
  sHeat.a[2] = 1.196253372934490E-06;
  sHeat.a[3] =-9.784803154493030E-10;
  sHeat.a[4] = 3.497377099477960E-13;
  sHeat.a[5] =-4.663062674693290E-17;

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 16/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * INITVISCOSITYPOL: inicializao a estrutura para o calculo da       *
 * viscosidade dinamica via polinomio                                *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * http://www.engineeringtoolbox.com/air-properties-d_156.html  	   *
 *                                                                   *
 * regressao com polinomio de ordem 5 obtido pelo excel              *
 *                                                                   *
 *********************************************************************/
void initDviscosityPol(char *s) {

  if(!strcmp(s,"polinomio")){
    dVisc.type = POL;
    dVisc.nPol = 6;
    dVisc.a[0] = 7.90864535186124E-03;
    dVisc.a[1] = 7.83417856886777E-03;
    dVisc.a[2] =-7.04677460750563E-06;
    dVisc.a[3] = 4.96383517580279E-09;
    dVisc.a[4] =-1.90489010749328E-12;
    dVisc.a[5] = 2.96188559332536E-16;
  }

  else if(!strcmp(s,"Sutherland")){
    dVisc.type = SUTHERLAND;
    dVisc.a[0] = 1.789e-05; /*viscosidade de referencia*/
    dVisc.a[1] = 273.11e0;  /*temperatura de referencia*/
    dVisc.a[2] = 110.56e0;  /*constante de Sutherland*/
  }
  else {
    ERRO_GERAL(__FILE__,__func__,__LINE__,s);
  }
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 16/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * INITTHCONDPOL: inicializao a estrutura para o calculo da          *
 * condutividade termica via polinomio                               *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * http://www.engineeringtoolbox.com/air-properties-d_156.html  	   *
 *                                                                   *
 * regressao com polinomio de ordem 5 obtido pelo excel              *
 *                                                                   *
 *********************************************************************/
void initThCondPol(char *s) {

  if(!strcmp(s,"polinomio")){
    thCond.type = POL;
    thCond.nPol = 6;
    thCond.a[0] =-1.08258923946027E-01;
    thCond.a[1] = 1.07023532477279E-02;
    thCond.a[2] =-6.17983109198392E-06;
    thCond.a[3] = 3.19214156395719E-09;
    thCond.a[4] =-9.80535866543648E-13;
    thCond.a[5] = 1.29223908056440E-16;
  }
  else if(!strcmp(s,"Sutherland")){
    dVisc.type = SUTHERLAND;
    dVisc.a[0] = 1.789e-05; /*viscosidade de referencia*/
    dVisc.a[1] = 273.11e0;  /*temperatura de referencia*/
    dVisc.a[2] = 110.56e0;  /*constante de Sutherland*/
  }
  else {
    ERRO_GERAL(__FILE__,__func__,__LINE__,s);
  }
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 06/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * GETTEMPFORENERGY: obtem a temperatura aprtir da entalpia sensivel *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * temp   - nao definido                                             *
 * energy - entalpia sensivel                                        * 
 * prop   - propriedades por material                                *
 * mat    - material da celula                                       *
 * nCell  - numero da celulas                                        *
 * fTemp  - equaca da energia na forma de temperatura (true|false)   *
 * fsHeat - variacao do calor especifico em funcao da Temperatura    *
 *          (true|false)                                             *
 * fKelnvin - temperatura em kelvin (true|false)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*
 * temp   - temperatura (C ou kelvin)                                * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void getTempForEnergy(DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                     ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                     ,INT const nCell      ,bool const fTemp
                     ,bool const fSheat    ,bool const fKelvin){
  
  short lMat;
  INT i;  
  DOUBLE sHeat;

/*... resolucao da eq da energia na forma de temperatura*/ 
  if(fTemp)
    for (i = 0; i < nCell; i++)
      temp[i] = energy[i]; 
/*...................................................................*/ 

/*... resolucao da eq da energia na forma de entalpia sensivel*/  
  else{
/*...*/
    for (i = 0; i < nCell; i++) {
      lMat  = mat[i] - 1;
      sHeat = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
      temp[i] = specificEnthalpyForTemp(energy[i], sHeat
                                       , fSheat  , fKelvin);
    }
/*...................................................................*/ 
  }
/*...................................................................*/ 

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 06/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * GETENERGYFORTEMP: obtem a entalpia sensivel apartir da temp       *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * temp   - temp                                                     *
 * energy - nao definido                                             * 
 * prop   - propriedades por material                                *
 * mat    - material da celula                                       *
 * nCell  - numero da celulas                                        *
 * fTemp  - equaca da energia na forma de temperatura (true|false)   *
 * fsHeat - variacao do calor especifico em funcao da Temperatura    *
 *          (true|false)                                             *
 * fKelnvin - temperatura em kelvin (true|false)                     *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*
 * energy - entalpia sensivel                                        * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void getEnergyForTemp(DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                     ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                     ,INT const nCell     
                     ,bool const fSheat    ,bool const fKelvin) {
  short lMat;
  INT i;  
  DOUBLE sHeatRef;

/*...*/ 
  for (i = 0; i < nCell; i++) {
    lMat  = mat[i] - 1;
    sHeatRef = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
    energy[i] = tempForSpecificEnthalpy(temp[i],sHeatRef
                                       ,fSheat,fKelvin);
  }
/*...................................................................*/ 

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 15/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * SPECIFICMASSREF : calcula a massa especifica de referencia        *
 * atraves da media do valores nas celulas                           * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * temp   - temp                                                     *
 * volume - volume das celulas                                       *
 * prop   - propriedades por material                                *
 * mat    - material da celula                                       *
 * nCell  - numero da celulas                                        *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*
 * densidade de referencia - entalpia sensivel                       * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void specificMassRef(DOUBLE *RESTRICT density, DOUBLE *RESTRICT volume                  
                  , DOUBLE *RESTRICT prop    , short  *RESTRICT mat
                  , INT const nCell)
{
  short nD = DENSITY_LEVEL;
  INT i;  
  DOUBLE dm,vm;

  dm = vm = 0.e0;

  for (i = 0; i < nCell; i++) {
/*...*/   
   dm += MAT2D(i,2 ,density ,nD)*volume[i];
   vm += volume[i];
/*...................................................................*/ 
  }

  MAT2D(0,DENSITY,prop,MAXPROP) = dm/vm;  

  printf("densityRef :%e\n",MAT2D(0,DENSITY,prop,MAXPROP));

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 15/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * PRESREF : calcula da pressao de referencia atualizada             *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * temp0  - temperatura no passo n                                   *
 * temp   - temperatura no passo n + 1                               *
 * volume - volume das celulas                                       *
 * pTh    - pressao termodinamica de referencia                      *
 * nCell  - numero da celulas                                        *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*
 * pTh    - pressao termodinamica de referencia atualizada           *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void presRef(DOUBLE *RESTRICT temp0         , DOUBLE *RESTRICT temp  
                   , DOUBLE *RESTRICT volume, DOUBLE *pTh                               
                   , INT const nCell        , bool const fKelvin)
{
  INT i;  
  DOUBLE dm,vm,t,t0;

  dm = vm = 0.e0;

  for (i = 0; i < nCell; i++) {
/*...*/
    if(fKelvin){
      t  = temp[i];
      t0 = temp0[i];
    }  
    else{
      t  = CELSIUS_FOR_KELVIN(temp[i]); 
      t0 = CELSIUS_FOR_KELVIN(temp0[i]);  
    }
/*...................................................................*/ 

/*...*/   
    dm += t/t0*volume[i];
    vm += volume[i];
/*...................................................................*/ 
  }

  pTh[2] = pTh[1]* dm/vm;
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 15/09/2017                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * INIPRESREF : incializa a pressao ref atrazes da massa especifica  * 
 * de referencia e temperatura media do dominio                      *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * temp0  - temperatura no passo n                                   *
 * temp   - temperatura no passo n + 1                               *
 * volume - volume das celulas                                       *
 * pTh    - pressao termodinamica de referencia                      *
 * nCell  - numero da celulas                                        *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*
 * pTh    - pressao termodinamica de referencia atualizada           *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void initPresRef(DOUBLE *RESTRICT temp  
               , DOUBLE *RESTRICT volume, DOUBLE *pTh   
               , DOUBLE *RESTRICT prop  , short  *RESTRICT mat                     
               , INT const nCell        , bool const fKelvin)
{
  INT i;  
  DOUBLE dRef,tm,dm,vm;

  dm = vm = 0.e0;

  for (i = 0; i < nCell; i++) {

/*...*/   
    dm += temp[i]*volume[i];
    vm += volume[i];
/*...................................................................*/ 
  }
  
  dm = dm/vm;

/*...*/
  if(fKelvin)
    tm  = dm;
  else
    tm  = CELSIUS_FOR_KELVIN(dm); 
/*...................................................................*/ 

  dRef = MAT2D(0,DENSITY,prop,MAXPROP);  
  vm   = PRESREF(dRef, IDEALGASR, tm, MMOLARAR);

  pTh[0] = pTh[1] = pTh[2] = vm;

}
/*********************************************************************/