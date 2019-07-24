#include<Reaction.h>


void reactor( Combustion *cModel      , Prop *sHeatPol
            , DOUBLE *RESTRICT yFrac
            , DOUBLE *RESTRICT dy     , DOUBLE const density
            , DOUBLE const tMix       , DOUBLE const temp 
            , bool const fKelvin)
{
 
}

/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * arrheniusLaw: Lei de arrhenius                                    *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * T     -> temperatura em kelvin                                    *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*  
 *********************************************************************/
static DOUBLE arrheniusLaw(Arrhenius *arr,DOUBLE const T)
{
  DOUBLE A = arr->A, beta = arr->beta, Ta = arr->Ta;   

  return A*pow(T,beta)*exp(-Ta/T);    

}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 21/07/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * molaRateReaction: taxa de reacao molar                            *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * chem  -> concentracao molar das especies                          *
 * Q     -> taxe de reacao molar o reacao (kmol/m3s)                 *
 * w     -> nao definido                                             *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * w     -> taxa de reacao das especies em kg/m3 s                   *
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*  
 *********************************************************************/
void massRateReaction(Chemical *chem  ,DOUBLE *RESTRICT Q
                     ,DOUBLE *RESTRICT w)
{
  unsigned short i, j, nReac =chem->nReac, nSp = chem->nSp;
  DOUBLE de;
  
  for(j=0;j<nSp;j++)
    w[j] = 0.e0;

  for(i=0;i<nReac;i++)
  { 
    for(j=0;j<nSp;j++)
    {
/*... */
      de = chem->reac[i].stch[2][j];
/*... kg/m3s*/
      w[j] += de*Q[i]*chem->sp[j].mW;
    }
  }
}
/********************************************************************/


/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * massActionMass:                                                   *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * c     -> concentracao molar das especies                          *
 * w     -> nao definido                                             *
 * T     -> temperatura em kelvin                                    *
 * nSp   -> numero de especies                                       *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS: mol/(cm3 s)                                                  *
 * ------------------------------------------------------------------*  
 *********************************************************************/
DOUBLE massActionMass(Reaction *reac     ,Prop *sHeatPol 
                     ,DOUBLE  *RESTRICT c
                     ,DOUBLE const T     ,unsigned short const nSp)
{
  unsigned short i;
  DOUBLE kf,kr,qf,qr;  

/*... reacao direta*/
  kf = arrheniusLaw(&reac->ArrF, T);

  for(i=0,qf=kf;i<nSp;i++)
    qf *= pow(c[i],reac->exp[0][i]);

/*... reacao inversao*/
  qr = 0.e0;
  if(reac->reverse)
  {
/*  kr = reverseK(reac       ,sHeatPol 
                 ,T          ,R
                 ,kf         ,1.0     
                 ,nSp);*/
    kr = arrheniusLaw(&reac->ArrR, T);

    for(i=0,qr=kr;i<nSp;i++)
      qr *= pow(c[i],reac->exp[1][i]);
  }

  return qf-qr;
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 26/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * arrhenius :                                                       *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * y1       -> fracao massica                                        *
 * e1       -> expoente                                              *
 * mW1      -> massa molar                                           *
 * y2       -> fracao massica                                        *
 * e2       ->  expoente                                             *
 * mW2      -> massa molar                                           *
 * t        -> temperatura                                           *
 * alpha    -> coeficiente da temperatura                            *
 * desnity  -> densidade do fluido                                   *
 * tA       -> temperatura de ativacao                               *
 * coefA    -> coeficiente                                           *
 * fKelvin  ->                                                       *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * omega :kmol/m^3s                                                  *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 * Reacao quimica a1A1 + a2A2 => a3A3 + a4A4                         * 
 * Unidade de Coef                                                   *
 * (N/L^3)^(1-(e1+e2)*(1/theta^alpha)(1/T)                           *
 * Exemplo:                                                          *
 * N     = mol                                                       *
 * L     = cm                                                        *
 * theta = Kelvin                                                    *
 * T     = segundos                                                  *
 * A = ((mol/cm^3)^(1-(e1+e2))/((K^alpha)(1/s))                      *
 *********************************************************************/
DOUBLE arrhenius(DOUBLE const y1     ,DOUBLE const y2
                ,DOUBLE const e1     ,DOUBLE const e2
                ,DOUBLE const mW1    ,DOUBLE const mW2
                ,DOUBLE const t      ,DOUBLE const alpha
                ,DOUBLE const density,DOUBLE const tA    
                ,DOUBLE const coefA  ,bool const fKelvin)
{
  short iCod=3;
  DOUBLE tc;
  DOUBLE omega,al,c1,c2,k,prodC,d;


  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);
/*... A = ((kmol/m^3)^(1-(e1+e2))/((K^alpha)(1/s))  */
  if(iCod == 1)
  {
    c1 = density*y1/mW1;
    c2 = density*y2/mW2;
    d  = 1.e0;
  }
/*..................................................................*/

/*... mol/m3*/
/*... A = ((mol/m^3)^(1-(e1+e2))/((K^alpha)(1/s))  */
  else if(iCod == 2)
  {
/*... kmol/m3 - > mol/m3*/
    c1 = 1.e+03*density*y1/mW1;
    c2 = 1.e+03*density*y2/mW2;
/*... mol/m^3s -> kmol/m^3s*/
    d  = 1.e-03;
  }
/*..................................................................*/

/*... A = ((mol/cm^3)^(1-(e1+e2))/((K^alpha)(1/s))  */
  else if(iCod == 3)
  {
/*... kmol/m3 - > mol/cm3*/
    c1 = 1.e-03*density*y1/mW1;
    c2 = 1.e-03*density*y2/mW2;
/*... mol/cm^3s -> kmol/m^3s*/
    d  = 1.e+03;
  }
/*..................................................................*/

  c1 = max(c1,0.0);
  c2 = max(c2,0.0);

/*... (c1^a1)x(c1^a2)*/
  prodC = pow(c1,e1)*pow(c2,e2);
/*... exp(Ea/RT)*/
  k = coefA*pow(tc,alpha)*exp(-tA/tc);

  omega = d*prodC*k;

  return omega; 
}
/*********************************************************************/ 