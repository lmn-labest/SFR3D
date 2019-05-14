#include<Properties.h>

/*******************   MISTURA GASOSA  ******************************/
static DOUBLE pol(DOUBLE *RESTRICT a, DOUBLE const x,short const n)
{
  int i; 
  DOUBLE tmp = a[0];

  for(i=1;i<n;i++)
    tmp += a[i]*pow(x,(double) i);

  return tmp;
}

/********************************************************************/

/*********************************************************************
 * Data de criacao    : 29/08/2017                                   *
 * Data de modificaco : 04/05/2019                                   *
 *-------------------------------------------------------------------*
 * airDensity: kg/(m^3)                                              *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * t - temperatura em kelvin                                         *
 * presRef - pressao de referencia ou termomecanica                  *
 * p       - pressao ( pressao do modelo)                            *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE mixtureSpeciesDensity(PropPol *den        ,DOUBLE const malorMassMix
                            ,DOUBLE const t      ,DOUBLE const p
                            ,DOUBLE const presRef,bool const fKelvin)
{
  short i,n=den->nPol[0];
  DOUBLE a[MAXPLODEG],tc,y,d;

  for (i = 0; i < MAXPLODEG; i++)
    a[i] = 0.0e0;
  
  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  
  
  switch (den->type)
  {
/*...*/
    case IDEALGAS:
      y = (malorMassMix*p)/(IDEALGASR*tc);
      d = 1.e+00; 
      break;
/*.....................................................................*/

/*...*/
    case INCIDEALGAS:
      y = (malorMassMix*presRef)/(IDEALGASR*tc);
      d = 1.e+00;
      break;
/*.....................................................................*/

/*...*/
    default:  
      ERRO_OP(__FILE__,__func__,den->type);
      break;
/*.....................................................................*/
  }

  return y*d;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 25/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * updateMixDensity:                                                 *
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
void updateMixDensity(PropPol *pDen         , Combustion *cModel
                 , DOUBLE *RESTRICT temp    , DOUBLE *RESTRICT pressure
                 , DOUBLE *RESTRICT density , DOUBLE *RESTRICT yComb
                 , DOUBLE const alpha       , bool const iKelvin    
                 , INT const nEl            , char  const iCod)

{
  short nD = DENSITY_LEVEL,ns = cModel->nOfSpecies;
  INT i;
  DOUBLE den,den0,molarMassMix,*y=NULL;
/*...*/
  switch (iCod){
    case PROP_UPDATE_NL_LOOP:
      for(i=0;i<nEl;i++)
      {
        y = &MAT2D(i,0,yComb,ns);
        molarMassMix =  mixtureMolarMass(cModel,y); 
        den0 =  MAT2D(i,TIME_N ,density ,nD);         
        den =  mixtureSpeciesDensity(pDen            ,molarMassMix
                                    ,temp[i]         ,pressure[i]
                                    ,thDynamic.pTh[2],iKelvin);
/*...*/           
        MAT2D(i,TIME_N ,density ,nD) =  alpha*den + (1.e0-alpha)*den0;
      }
/*..................................................................*/
    break;  

    case PROP_UPDATE_OLD_TIME:
      for(i=0;i<nEl;i++)
      {
/*...t(n-2) = t(n-1)*/
        MAT2D(i,TIME_N_MINUS_2 ,density ,nD) = MAT2D(i,1 ,density ,nD);
/*...t(n-1) = t(n)*/           
        MAT2D(i,TIME_N_MINUS_1 ,density ,nD) = MAT2D(i,2 ,density ,nD);
      }
/*..................................................................*/
    break;
  }
/*..................................................................*/


}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 25/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * initMixtureSpeciesfiHeat:                                         *
 * ----------------------------------------------------------------- *
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void initMixtureSpeciesfiHeat(PropPol *prop, char *s, FILE *file)
{
  FILE *fileOut;
  
  char word[WORD_SIZE];
  char species[][WORD_SIZE] = { "fuel","n2"
                               ,"o2"  ,"co2" 
                               ,"h2o"};
  char nameAux[1000];
  short i,j,k;
  int nSpecies;
  double x[MAXPLODEG],g;


  if (!strcmp(s, "polinomio")) 
  {
    prop->type = POL;
  
    fscanf(file, "%s", nameAux);
    fileOut = openFile(nameAux, "r");

    fscanf(fileOut, "%d", &nSpecies);    

    for(j=0;j<nSpecies;j++)
    {

      readMacro(fileOut,word,false);
      convStringLower(word);
      
      fscanf(fileOut, "%lf", &g);

      for(k=0;k<MAXSPECIES;k++)
        if(!strcmp(word,species[k]))
          break;

      prop->nPol[j] = readFileLineSimple(x, fileOut);
      ERRO_POL_READ(prop->nPol[j], MAXPLODEG, __FILE__, __func__, __LINE__);

      for (i = 0; i < prop->nPol[j]; i++)
        MAT2D(k,i,prop->a,MAXPLODEG) = x[i]/g;
    }

    fclose(fileOut);

    fprintf(fileLogExc, "%-25s: %s\n", "Type", s);

  }
  
  printf("Write SpeciesfiHeat cp(T):\n");
  fileOut = openFile("species_cp.out", "w");
  for(i=0;i<nSpecies;i++)
  {
    fprintf(fileOut,"%s\n",species[i]);
    for (j = 0, g =200; j < 40; j++) 
    {
      fprintf(fileOut,"%lf %lf\n",g,pol(&MAT2D(i,0,prop->a,10)
                                   ,g,prop->nPol[i]));
      g += 200;
    }
  }
  fclose(fileOut);

}
/********************************************************************/

/********************************************************************* 
 * Data de criacao    : 20/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * etEnergyForTempMix: obtem a entalpia sensivel apartir da temp     *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * sHeatPol - estrutra par o polinoimio do calor especifico          *
 * yFrac    - fracao de massa da especies primitivas                 * 
 * temp   - temp                                                     *
 * energy - nao definido                                             * 
 * prop   - propriedades por material                                *
 * mat    - material da celula                                       *
 * nCell  - numero da celulas                                        *
 * nOfPrSp  - numero de especies primitivas                          *
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
void getEnergyForTempMix(PropPol *sHeatPol    ,DOUBLE *RESTRICT yFrac 
                        ,DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                        ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                        ,INT const nCell      ,short const nOfPrSp
                        ,bool const fSheat    ,bool const fKelvin
                        ,bool const fOmp      ,short const nThreads ) 
{
  short lMat;
  INT i;  
  DOUBLE sHeatRef,*y;

/*...*/
  if (fOmp) 
  {
/*...*/ 
#pragma omp parallel  for default(none) num_threads(nThreads)\
    private(i,lMat,sHeatRef,y) shared(prop,mat,energy,temp,sHeatPol,nOfPrSp,yFrac)
    for (i = 0; i < nCell; i++) 
    {
      lMat  = mat[i] - 1;
      sHeatRef = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
      y = &MAT2D(i,0,yFrac,nOfPrSp);
      energy[i] = tempForSpecificEnthalpyMix(sHeatPol ,yFrac
                                            ,temp[i]  ,sHeatRef
                                            ,nOfPrSp
                                            ,fSheat  ,fKelvin);  
    }
/*...................................................................*/ 
  }
/*...................................................................*/ 

  else
  {
/*...*/ 
    for (i = 0; i < nCell; i++) 
    {
      lMat  = mat[i] - 1;
      sHeatRef = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
      y = &MAT2D(i,0,yFrac,nOfPrSp);
      energy[i] = tempForSpecificEnthalpyMix(sHeatPol ,yFrac
                                            ,temp[i]  ,sHeatRef
                                            ,nOfPrSp
                                            ,fSheat  ,fKelvin); 
    }
/*...................................................................*/ 
  }
/*...................................................................*/ 
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 20/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * getTempForEnergyMix : obtem a temperatura apartir da entalpia     *  
 * sensivel                                                          *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sHeatPol - estrutra par o polinoimio do calor especifico          *
 * yFrac    - fracao de massa da especies primitivas                 * 
 * temp   - nao definido                                             *
 * energy - entalpia sensivel                                        * 
 * prop   - propriedades por material                                *
 * mat    - material da celula                                       *
 * nCell  - numero da celulas                                        *
 * nOfPrSp  - numero de especies primitivas                          *
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
void getTempForEnergyMix(PropPol *sHeatPol    ,DOUBLE *RESTRICT yFrac
                        ,DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                        ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                        ,INT const nCell      ,short const nOfPrSp 
                        ,bool const fTemp     ,bool const fSheat    
                        ,bool const fKelvin
                        ,bool const fOmp      ,short const nThreads )
{
  
  short lMat;
  INT i;  
  DOUBLE sHeatRef,*y;

/*... resolucao da eq da energia na forma de temperatura*/ 
  if(fTemp)
    for (i = 0; i < nCell; i++)
      temp[i] = energy[i]; 
/*...................................................................*/ 

/*... resolucao da eq da energia na forma de entalpia sensivel*/  
  else{
    if(fOmp){
#pragma omp parallel  for default(none) num_threads(nThreads)\
      private(i,lMat,sHeatRef,y) shared(prop,mat,energy,temp,sHeatPol,nOfPrSp,yFrac)
      for (i = 0; i < nCell; i++) {
        lMat  = mat[i] - 1;
        sHeatRef = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
        y = &MAT2D(i,0,yFrac,nOfPrSp);
        temp[i] = specificEnthalpyForTempOfMix(sHeatPol
                                         , energy[i], y    
                                         , sHeatRef , nOfPrSp
                                         , fSheat   , fKelvin
                                         , i);
      }
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    else{
      for (i = 0; i < nCell; i++) {
        lMat  = mat[i] - 1;
        sHeatRef = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
        y = &MAT2D(i,0,yFrac,nOfPrSp);
        temp[i] = specificEnthalpyForTempOfMix(sHeatPol
                                         , energy[i], y    
                                         , sHeatRef , nOfPrSp
                                         , fSheat   , fKelvin
                                         , i);
      }
/*...................................................................*/
    }
/*...................................................................*/ 
  }
/*...................................................................*/ 

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * updateMixSpecificHeat:                                            *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * sHeatPol - estrutra par o polinoimio do calor especifico          *
 * temp     - temperatura (°C/K)                                     *
 * yFrac    - fracao de massa da especies primitivas                 *
 * sHeat    - calor especifico por celula                            *
 * nOfPrSp  - numero de especies primitivas                          * 
 * fKelvin  - temperatura dada em kelvin                             *
 * nEl      - numero total de elmentos                               *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * sHeat    - calor especifico por celula atualizado                 *
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void updateMixSpecificHeat(PropPol *sHeatPol
                         , DOUBLE *RESTRICT temp  , DOUBLE *RESTRICT yFrac  
                         , DOUBLE *RESTRICT sHeat , short const nOfPrSp
                         , bool const iKelvin
                         , INT const nEl          , char  const iCod)

{
  short nD = SHEAT_LEVEL;
  INT i;
  DOUBLE *y;  
  
/*...*/
  switch (iCod)
  {
    case PROP_UPDATE_NL_LOOP:
      for(i=0;i<nEl;i++)
      {
        y = &MAT2D(i,0,yFrac,nOfPrSp);
/*...*/           
        MAT2D(i,TIME_N ,sHeat ,nD) = 
        mixtureSpecifiHeat(sHeatPol ,y,temp[i]  ,nOfPrSp,iKelvin);
      }
/*..................................................................*/
    break;  

  case PROP_UPDATE_OLD_TIME:
    for(i=0;i<nEl;i++){
/*...*/
      MAT2D(i, TIME_N_MINUS_2, sHeat, nD) = MAT2D(i,1 ,sHeat ,nD);           
      MAT2D(i, TIME_N_MINUS_1, sHeat, nD) = MAT2D(i,2 ,sHeat ,nD);
    }
/*..................................................................*/
    break;
  }
/*..................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 29/08/2017                                   *
 * Data de modificaco : 07/05/2019                                   *
 *-------------------------------------------------------------------*
 * mixtureSpecifiHeat: kJ/(kg.K)                                     *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * sHeatPol - estrutra par o polinoimio do calor especifico          *
 * yFrac    - fracao de massa da especies primitivas                 * 
 * t - temperatura em Kelvin                                         *
 * nOfPrSp  - numero de especies primitivas                          *
 * fKelvin  - temperatura dada em kelvin                             *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE mixtureSpecifiHeat(PropPol *sHeat    , DOUBLE *yFrac
                         , DOUBLE const t    , short const mOfPrSp
                         , bool const fKelvin) 
{

  short i,k,n;  
  DOUBLE a[MAXPLODEG],cpk,cp,d;
  DOUBLE tc;

  a[0] = 0.0e0;
  
  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  for(k=0,cp=0.e0;k<mOfPrSp;k++)
  {
    n = sHeat->nPol[k];
    for (i = 0; i < n; i++)
      a[i] = MAT2D(k,i,sHeat->a,MAXPLODEG);
  
/*... polinomio*/
    cpk = a[0];
    for (i = 1; i < n; i++)
      cpk += a[i]*pow(tc,i);
/*.....................................................................*/

    if (cpk < 0.e0 || yFrac[k] < 0.e0)
    {
      printf("Calor especifico negativo!!\n"
             "Y                = %e\n"
             "Species          = %d\n"   
             "Calor especifico = %e\n"
             "Temperatura      = %lf\n!!",yFrac[k], k,cpk,tc);
      exit(EXIT_FAILURE);
    }

/*...*/
    cp += yFrac[k]*cpk;
/*.....................................................................*/
  }
/*.....................................................................*/

/*...*/
  d = 1.e0;
/*.....................................................................*/

 if (cp < 0.e0) {
    printf("Calor especifico negativo!!"
           "Calor especifico = %e\n"
           "Temperatura      = %lf\n!!",d*cp,tc);
    exit(EXIT_FAILURE);
  }

  return d*cp;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 06/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * specieSpecifiHeat: calor especifico da especie k kJ/(kg.K)        *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * sHeatPol - estrutra par o polinoimio do calor especifico          *
 * kSpecie  - fracao de massa da especies primitivas                 * 
 * t - temperatura em Kelvin                                         *
 * fKelvin  - temperatura dada em kelvin                             *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE specieSpecifiHeat(PropPol *sHeat     , short const kSpecie
                        , DOUBLE const t    , bool const fKelvin) 
{

  short i,n;
  DOUBLE a[MAXPLODEG],cp,d,tc;

  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  n=sHeat->nPol[kSpecie];
  printf("%d %d\n",n,kSpecie);
  for (i = 0; i < n; i++)
    a[i] = MAT2D(kSpecie,i,sHeat->a,MAXPLODEG);
  
/*... polinomio*/
  cp = a[0];
  for (i = 1; i < n; i++)
    cp += a[i]*pow(tc,i);
/*.....................................................................*/

  if (cp < 0.e0 )
  {
    printf("Calor especifico negativo!!"
           "Species          = %hd "   
           "Calor especifico = %e\n"
           "Temperatura      = %lf\n!!",kSpecie, cp, tc);
    exit(EXIT_FAILURE);
  }

/*...*/
  d = 1.e0;
/*.....................................................................*/

 if (cp < 0.e0) {
    printf("Calor especifico negativo!!"
           "Calor especifico = %e\n"
           "Temperatura      = %lf\n!!",d*cp,tc);
    exit(EXIT_FAILURE);
  }

  return d*cp;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/08/2018                                   *
 * Data de modificaco : 07/05/2019                                   *
 *-------------------------------------------------------------------*
 * SPECIFICENTHALPYFORTEMPMIX:  calcula a temperatura apartir da     *
 * entalpia especifica                                               * 
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * sHeatPol - estrutra par o polinoimio do calor especifico          *
 * hs       - entalpia sensivel                                      *
 * yFrac    - fracao de massa da especies primitivas                 * 
 * sHeatRef - calor especifico de referencia constante com temp      *
 * nOfPrSp  - numero de especies primitivas                          *
 * fSheat   - calor especifico com variacao com a Temperatura        *
 * fKelvin  - temperatura dada em kelvin                             *
 * nel      - numero da celula  
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * temperatura (°C/Kelvin)                                           *
 *-------------------------------------------------------------------*
 * OBS:                                                              * 
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE specificEnthalpyForTempOfMix(PropPol *sHeatPol
                             , DOUBLE const hs        , DOUBLE *yFrac
                             , DOUBLE const sHeatRef  , short const nOfPrSp
                             , bool const fSheat      , bool const fKelvin
                             , INT const nel ) 
{
  unsigned short i;
  bool flag = false;
  DOUBLE f,fl,t,conv,tol=1e-04;
 
/*...*/
  if(fSheat)
  {
/*... chute inicial usando a massa espeficia constante*/

    t = 0.5e0*(TREF + ENTHALPY_FOR_TEMP(sHeatRef,hs,TREF))+tol;
    if(t > 4000.0) t = 4000.0;
 
/*...*/
    conv = tol*(hs-tempForSpecificEnthalpyMix(sHeatPol,yFrac
                                             ,t       ,sHeatRef
                                             ,nOfPrSp  
                                             ,fSheat  ,true));
    conv = fabs(conv);
/*... Newton-Raphson*/
    for(i=0;i<60000;i++)
    {
      f  = hs-tempForSpecificEnthalpyMix(sHeatPol,yFrac
                                        ,t       ,sHeatRef
                                        ,nOfPrSp 
                                        ,fSheat  ,true);
      if(fabs(f) < conv || fabs(f) == 0.e0) 
      {
        flag = true;
        break;
      }
    
      fl = mixtureSpecifiHeat(sHeatPol,yFrac
                             ,t       ,nOfPrSp
                             ,true);
      t += f/fl;   
    }
/*...................................................................*/

    if(!flag)
    {
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
 * Data de criacao    : 19/08/2018                                   *
 * Data de modificaco : 07/05/2019                                   *
 *-------------------------------------------------------------------*
 * TEMPFORSPECIFICENTHALPY: calcula a entalpia espeficia apartir da  *
 * temperatura                                                       *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * sHeat    - estrutra par o polinoimio do calor especifico          *
 * yFrac    - fracao de massa da especies primitivas                 *
 * t        - temperatura (°C/K)                                     *
 * sHeatRef - calor especifico de referencia constante com temp      *
 * nOfPrSp  - numero de especies primitivas                          *
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
DOUBLE tempForSpecificEnthalpyMix(PropPol *sHeat    , DOUBLE *yFrac
                                , DOUBLE const t    , DOUBLE const sHeatRef
                                , short const nOfPrSp
                                , bool const fSheat , bool const fKelvin) 
{

  short i,k,n;
  DOUBLE a[MAXPLODEG],d,dt,hk,hs;
  DOUBLE tc,tRef= TREF ;

  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  hs = 0.e0;
  if(fSheat)
  {
    for(k = 0; k < nOfPrSp; k++)
    {
      n = sHeat->nPol[k];
      for (i = 0; i < n; i++)
        a[i] = MAT2D(k,i,sHeat->a,MAXPLODEG);
  
      for (i = 0, hk =0.e0; i < n; i++) 
      {
        d    = (double) (i + 1);
        dt   =pow(tc,d) - pow(tRef,d);
        hk += a[i]*dt/d;
      }

      hs += yFrac[k]*hk; 
    }
  }

  else 
    hs = TEMP_FOR_ENTHALPY(sHeatRef,tc,TREF);

  return hs;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 07/05/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * TEMPFORSPECIFICENTHALPYSPECIES: calcula a entalpia espeficia da   *
 * especie k partir da temperatura                                   *                   *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * sHeat    - estrutra par o polinoimio do calor especifico          *
 * yFrac    - fracao de massa da especies primitivas                 *
 * t        - temperatura (°C/K)                                     *
 * nOfPrSp  - numero de especies primitivas                          *
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
DOUBLE tempForSpecificEnthalpySpecies(PropPol *sHeat, short const kSpecie
                               , DOUBLE const t     , DOUBLE const sHeatRef
                               , bool const fSheat  , bool const fKelvin) 
{

  short i,n;
  DOUBLE a[MAXPLODEG],d,dt,hk;
  DOUBLE tc,tRef= TREF;

  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  if(fSheat)
  {
    n=sHeat->nPol[kSpecie];
    for (i = 0; i < n; i++)
      a[i] = MAT2D(kSpecie,i,sHeat->a,MAXPLODEG);
  
    for (i = 0, hk =0.e0; i < n; i++) 
    {
      d    = (double) (i + 1);
      dt   =pow(tc,d) - pow(tRef,d);
      hk += a[i]*dt/d;
    }
  }

  else
    hk       = TEMP_FOR_ENTHALPY(sHeatRef,tc,TREF);
  
  return hk;

}
/**********************************************************************/

/********************************************************************* 
 * Data de criacao    : 20/08/2018                                   *
 * Data de modificaco : 25/08/2018                                   *
 *-------------------------------------------------------------------*
 * initPropTempMix: inicializao de propriedades com variacao temporal*
 * dependentes da temperatura                                        *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * propFluid -> estrutura de dados das propriedades varaiveis        *
 * prop      -> nao definido                                         * 
 * t         -> temperatura                                          *
 * pressure  -> pressao                                              *
 * yFrac    - fracao de massa da especies primitivas                 * 
 * propMat   -> propriedade de referencia por material               * 
 * mat       -> material por celula                                  * 
 * nOfPrSp  - numero de especies primitivas                          *
 * np        -> numero niveis de tempos                              * 
 * nCell     -> numero de celulas                                    * 
 * iKelvin   -> temperatura em kelvin (true|false)                   *
 * iProp     -> numero da propriedade                                * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * prop    -> propriedade iniciacializada                            * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
  *********************************************************************/
void initPropTempMix(PropVarFluid *propFluid, Combustion *cModel
                 ,DOUBLE *RESTRICT prop 
                 ,DOUBLE *RESTRICT t        ,DOUBLE *RESTRICT pressure
                 ,DOUBLE *RESTRICT yFrac    ,DOUBLE *RESTRICT propMat
                 ,short *RESTRICT mat       ,short const nOfPrSp 
                 ,short const np            ,INT    const nCell 
                 ,bool const iKelvin        ,short const iProp)
{    
  INT i;
  unsigned short j,lMat;
  DOUBLE *y,molarMassMix;
         
  for(i=0;i<nCell;i++){    

/*...*/
    lMat               = mat[i]-1;
/*...................................................................*/

/*...*/
    if( iProp == DENSITY )
    {
      y = &MAT2D(i,0,yFrac, nOfPrSp);
      molarMassMix =  mixtureMolarMass(cModel,y); 
      MAT2D(lMat, iProp, propMat, MAXPROP) 
        = mixtureSpeciesDensity(&propFluid->den,molarMassMix
                      ,t[i]                    , pressure[i]
                      ,thDynamic.pTh[2]        , iKelvin);
    }
/*...................................................................*/

/*...*/
    else if( iProp == SPECIFICHEATCAPACITYFLUID)
    {
      y = &MAT2D(i,0,yFrac, nOfPrSp);
      MAT2D(lMat, iProp, propMat, MAXPROP) 
      = mixtureSpecifiHeat(&propFluid->sHeat, y
                           ,t[i]            , nOfPrSp
                           ,iKelvin);
    }
/*...................................................................*/

/*...*/
    else if( iProp == DYNAMICVISCOSITY)  
      MAT2D(lMat, iProp, propMat, MAXPROP) 
      = airDynamicViscosity(&propFluid->dVisc,t[i],iKelvin);
/*...................................................................*/

/*...*/
    else if( iProp == THERMALCONDUCTIVITY)  
      MAT2D(lMat,iProp,propMat,MAXPROP) 
      = airThermalConductvity(&propFluid->thCond,t[i],iKelvin);
/*...................................................................*/

/*...*/
    for(j=0;j<np;j++)      
      MAT2D(i,j,prop,np) = MAT2D(lMat,iProp,propMat,MAXPROP); 
/*...................................................................*/
  
  }
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 29/08/2017                                   *
 * Data de modificaco : 12/07/2018                                   *
 *-------------------------------------------------------------------*
 * airDensity: kg/(m^3)                                              *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * t - temperatura em kelvin                                         *
 * presRef - pressao de referencia ou termomecanica                  *
 * p       - pressao ( pressao do modelo)                            *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
DOUBLE airDensity(PropPol *den
                 ,DOUBLE const t      ,DOUBLE const p
                 ,DOUBLE const presRef,bool const fKelvin) {
  short i,n=den->nPol[0];
  DOUBLE a[MAXPLODEG],tc,y,d;

  for (i = 0; i < MAXPLODEG; i++)
    a[i] = 0.0e0;
  
  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  switch (den->type) {
/*... polinomio*/
    case POL:
      for (i = 0; i < n; i++)
        a[i] = den->a[i];

/*... polinomio*/
      y = a[0];
      for (i = 1; i < n; i++)
        y += a[i]*pow(tc,i);
/*.....................................................................*/

/*.....................................................................*/
      d = 1.e0; 
      break;
/*.....................................................................*/
 
/*...*/
    case IDEALGAS:
      y = (MMOLARAR*p)/(IDEALGASR*tc);
      d = 1.e0; 
      break;
/*.....................................................................*/

/*...*/
    case INCIDEALGAS:
      y = (MMOLARAR*presRef)/(IDEALGASR*tc);
      d = 1.e0; 
      break;
/*.....................................................................*/


/*...*/
    default:  
      ERRO_OP(__FILE__,__func__,den->type);
      break;
/*.....................................................................*/
  }

  return y*d;

}
/**********************************************************************/

/*********************************************************************
* Data de criacao    : 12/05/2018                                   *
* Data de modificaco : 07/05/2019                                   *
*-------------------------------------------------------------------*
* diffProp:                                                         *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* u - valor da propriedade                                          *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* retorna a valor da massa especifica                               *
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
DOUBLE diffProp(PropPol *pol  , DOUBLE u) 
{
  short i, n = pol->nPol[0];
  DOUBLE a[MAXPLODEG], y;

  for (i = 0; i < MAXPLODEG; i++)
    a[i] = 0.0e0;

  switch (pol->type)
  {
/*... polinomio*/
    case POL:
      for (i = 0; i < n; i++)
        a[i] = pol->a[i];

/*... polinomio*/
      y = a[0];
      for (i = 1; i < n; i++)
        y += a[i] * pow(u, i);
/*.....................................................................*/

/*.....................................................................*/
      break;
/*.....................................................................*/


/*...*/
    default:
      ERRO_OP(__FILE__, __func__, pol->type);
      break;
/*.....................................................................*/
  }

  return y;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 29/08/2017                                   *
 * Data de modificaco : 07/05/2019                                   *
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
DOUBLE airSpecifiHeat(PropPol *sHeat, DOUBLE const t
                     ,bool const fKelvin) {

  short i,n=sHeat->nPol[0];  
  DOUBLE a[MAXPLODEG],y,d;
  DOUBLE tc;

  a[0] = 0.0e0;
  
  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  for (i = 0; i < n; i++)
    a[i] = sHeat->a[i];
  
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
 * Data de modificaco : 07/05/2019                                   *
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
DOUBLE airDynamicViscosity(PropPol *dVisc,DOUBLE const t
                          ,bool const fKelvin) {

  short i,n=dVisc->nPol[0];
  DOUBLE a[MAXPLODEG],x[MAXPLODEG-1],y,d;
  DOUBLE tc;

  a[0] = 0.0e0;
  
  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  switch (dVisc->type) {
/*... polinomio*/
    case POL:
      for (i = 0; i < n; i++)
        a[i] = dVisc->a[i];
/*.....................................................................*/
  
/*... polinomio*/
      y = a[0];
      for (i = 1; i < n; i++)
        y += a[i]*pow(tc,i);
      d = 1.e-05;
/*.....................................................................*/
      break;
/*.....................................................................*/

/*... polinomio*/
    case SUTHERLAND:
      a[0] = dVisc->a[0]; /*viscosidade de referencia*/
      a[1] = dVisc->a[1]; /*temperatura de referencia*/ 
      a[2] = dVisc->a[2]; /*constante de Sutherland*/

      x[0] = a[1]+ a[2];
      x[1] = tc  + a[2];

      y = a[0]*(x[0]/x[1])*pow(tc/a[1],1.5);

      d = 1.e0;
/*.....................................................................*/
      break;
/*.....................................................................*/

/*...*/
    default:  
      ERRO_OP(__FILE__,__func__,dVisc->type);
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
 * Data de modificaco : 07/05/2019                                   *
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
DOUBLE airThermalConductvity(PropPol *thCond, DOUBLE const t 
                            ,bool const fKelvin) {

  short i,n=thCond->nPol[0];  
  DOUBLE a[MAXPLODEG],y,d;
  DOUBLE tc;

  a[0] = 0.0e0;

  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  
  switch (thCond->type) {
/*... polinomio*/
    case POL:
      for (i = 0; i < n; i++)
        a[i] = thCond->a[i];
/*.....................................................................*/

/*... polinomio*/
      y = a[0];
      for (i = 1; i < n; i++)
        y += a[i]*pow(tc,i);
      d = 1.e-05;
/*.....................................................................*/
      break;
/*.....................................................................*/
 /*...*/
    default:  
      ERRO_OP(__FILE__,__func__,thCond->type);
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
 * Data de modificaco : 07/05/2019                                   *
 *-------------------------------------------------------------------*
 * TEMPFORSPECIFICENTHALPY: calcula a entalpia espeficia apartir da  *
 * temperatura                                                       *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * sHeat    - estrutra par o polinoimio do calor especifico          *
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
DOUBLE tempForSpecificEnthalpy(PropPol *sHeat
                             , DOUBLE const t   , DOUBLE const sHeatRef
                             , bool const fSheat, bool const fKelvin) {

  short i,n=sHeat->nPol[0];
  DOUBLE a[6],d,dt,tmp;
  DOUBLE tc,tRef= TREF ;

  if(fKelvin)
    tc = t;  
  else
    tc = CELSIUS_FOR_KELVIN(t);  

  if(fSheat){
    for (i = 0; i < n; i++)
      a[i] = sHeat->a[i];

    tmp = 0.0;
    for (i = 0; i < n; i++) {
      d    = (double) (i + 1);
      dt   =pow(tc,d) - pow(tRef,d);
      tmp += a[i]*dt/d;
    }
  }

  else 
    tmp = TEMP_FOR_ENTHALPY(sHeatRef,tc,TREF);

  return tmp;

}
/**********************************************************************/

/*********************************************************************
 * Data de criacao    : 28/08/2017                                   *
 * Data de modificaco : 15/07/2018                                   *
 *-------------------------------------------------------------------*
 * SPECIFICENTHALPYFORTEMP:  calcula a temperatura apartir da        *
 * entalpia especifica                                               * 
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * sHeatPol - estrutra par o polinoimio do calor especifico          *
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
DOUBLE specificEnthalpyForTemp(PropPol *sHeatPol
                             , DOUBLE const hs  , DOUBLE const sHeatRef
                             , bool const fSheat, bool const fKelvin) 
{
  unsigned short i;
  bool flag = false;
  DOUBLE f,fl,t,conv,tol=1e-04;
 
/*...*/
  if(fSheat){
/*... chute inicial usando a massa espeficia constante*/
    t = ENTHALPY_FOR_TEMP(sHeatRef,hs,TREF);
/*...*/
    conv = (hs-tempForSpecificEnthalpy(sHeatPol,t,sHeatRef,fSheat,true))*tol;
    conv = fabs(conv);
/*... Newton-Raphson*/
    for(i=0;i<60000;i++){
      f  = hs-tempForSpecificEnthalpy(sHeatPol,t,sHeatRef,fSheat,true);
      if(fabs(f) < conv) {
        flag = true;
        break;
      }
    
      fl = airSpecifiHeat(sHeatPol,t,true);
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

  DOUBLE a[MAXPLODEG],x[MAXPLODEG-1],y,d;

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

  DOUBLE a[MAXPLODEG],x[MAXPLODEG-1],y,d;

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

  DOUBLE a[MAXPLODEG],x[MAXPLODEG-1],y,d;

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
 * Data de modificaco : 12/07/2018                                   *
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
void updateDensity(PropPol *pDen
                 , DOUBLE *RESTRICT temp    , DOUBLE *RESTRICT pressure
                 , DOUBLE *RESTRICT density                 
                 , DOUBLE const alpha       , bool const iKelvin    
                 , INT const nEl            , char  const iCod)

{
  short nD = DENSITY_LEVEL;
  INT i;
  DOUBLE den,den0;
/*...*/
  switch (iCod){
    case PROP_UPDATE_NL_LOOP:
      for(i=0;i<nEl;i++){
        den0 =  MAT2D(i,2 ,density ,nD);         
        den = airDensity(pDen,temp[i],pressure[i],thDynamic.pTh[2]
                        ,iKelvin);
/*...*/           
        MAT2D(i,TIME_N ,density ,nD) =  alpha*den + (1.e0-alpha)*den0;
      }
/*..................................................................*/
    break;  

  case PROP_UPDATE_OLD_TIME:
    for(i=0;i<nEl;i++){
/*...t(n-2) = t(n-1)*/
      MAT2D(i,TIME_N_MINUS_2 ,density ,nD) = MAT2D(i,1 ,density ,nD);
/*...t(n-1) = t(n)*/           
      MAT2D(i,TIME_N_MINUS_1 ,density ,nD) = MAT2D(i,2 ,density ,nD);
    }
/*..................................................................*/
    break;
  }
/*..................................................................*/

}
/*********************************************************************/ 

/*********************************************************************
 * Data de criacao    : 27/08/2017                                   *
 * Data de modificaco : 15/07/2018                                   *
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
void updateSpecificHeat( PropPol *sHeatPol
                       , DOUBLE *RESTRICT temp, DOUBLE *RESTRICT sHeat
                       , bool const iKelvin 
                       ,INT const nEl        , char  const iCod)

{
  short nD = SHEAT_LEVEL;
  INT i;  
  
/*...*/
  switch (iCod)
  {
    case PROP_UPDATE_NL_LOOP:
      for(i=0;i<nEl;i++)
/*...*/           
        MAT2D(i,TIME_N ,sHeat ,nD) = airSpecifiHeat(sHeatPol,temp[i],iKelvin);
/*..................................................................*/
    break;  

  case PROP_UPDATE_OLD_TIME:
    for(i=0;i<nEl;i++){
/*...*/
      MAT2D(i, TIME_N_MINUS_2, sHeat, nD) = MAT2D(i,1 ,sHeat ,nD);           
      MAT2D(i, TIME_N_MINUS_1, sHeat, nD) = MAT2D(i,2 ,sHeat ,nD);
    }
/*..................................................................*/
    break;
  }
/*..................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 27/08/2017                                   *
 * Data de modificaco : 12/07/2018                                   *
 *-------------------------------------------------------------------*
 * UPDATESPECIFICHEAT:                                               *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * thC     -> intepolacao a condutividade termica                    *
 * t       -> temperatura                                            *
 * thCond  -> vetor de condutividade termica por celula              *
 * iKelvin -> kelvin ou celsus                                       *
 * nEl     -> numero total de celulas                                *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void updateThermalconductivity(PropPol *thC
                              ,DOUBLE *RESTRICT t,DOUBLE *RESTRICT thCond   
                              ,bool const iKelvin,INT const nEl)
{
  INT i;

  for(i=0;i<nEl;i++)         
    thCond[i] = airThermalConductvity(thC,t[i],iKelvin);

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 01/09/2017                                   *
 * Data de modificaco : 15/07/2018                                   *
 *-------------------------------------------------------------------*
 * updateDynamicViscosity:                                           *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * dVisc   -> intepolacao da viscosidade dinamica                    *
 * t       -> temperatura                                            *
 * visc    -> vetor de viscosidades dianamica por celula             *
 * iKelvin -> kelvin ou celsus                                       *
 * nEl     -> numero total de celulas                                *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void updateDynamicViscosity(PropPol *dVisc
                          ,DOUBLE *RESTRICT temp ,DOUBLE *RESTRICT visc    
                          ,bool const iKelvin    ,INT const nEl)
{
  INT i;

  for(i=0;i<nEl;i++)         
    visc[i] = airDynamicViscosity(dVisc,temp[i],iKelvin);

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 13/05/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * updateProp:                                                       *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * pol     - polinomico de interpolacao com u                        *
 * u       - variavel                                                *
 * coef    - coeficiente a ser atualizado                            *
 * nEl     - numero de elementos                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void updateProp(PropPol *pol         , DOUBLE *RESTRICT u
              , DOUBLE *RESTRICT coef, INT nEl)

{
  INT i;

  for (i = 0; i<nEl; i++)
    coef[i] = diffProp(pol, u[i]);

}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 12/05/2018                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* UPDATESPECIFICHEAT:                                               *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* pol     - polinomico de interpolacao com u                        *
* u       - variavel                                                *
* density - massa especifica                                        *
* nEl     - numero de elementos                                     *
* iCod    - codigo                                                  *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
void updateDensityCD(PropPol *pol            , DOUBLE *RESTRICT u
                   , DOUBLE *RESTRICT density, INT nEl    
                   , char  iCod)

{
  short nD = SHEAT_LEVEL;
  INT i;

/*...*/
  switch (iCod) {
  case PROP_UPDATE_NL_LOOP:
    for (i = 0; i<nEl; i++)
/*...*/
      MAT2D(i, TIME_N, density, nD) = diffProp(pol,u[i]);
/*..................................................................*/
    break;

  case PROP_UPDATE_OLD_TIME:
    for (i = 0; i<nEl; i++) {
/*...*/
      MAT2D(i,TIME_N_MINUS_2,density, nD) = MAT2D(i,1,density, nD);
      MAT2D(i,TIME_N_MINUS_1,density, nD) = MAT2D(i,2,density, nD);
    }
/*..................................................................*/
    break;
  }
/*..................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 27/08/2017                                   *
 * Data de modificaco : 12/07/2018                                   *
 *-------------------------------------------------------------------*
 * INITPROPTEMP: inicializao de propriedades com variacao temporal   *
 * dependentes da temperatura                                        *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * propFluid -> estrutura de dados das propriedades varaiveis        *
 * prop      -> nao definido                                         * 
 * t         -> temperatura                                          *
 * pressure  -> pressao                                              *
 * propMat   -> propriedade de referencia por material               * 
 * mat       -> material por celula                                  * 
 * np        -> numero niveis de tempos                              * 
 * nCell     -> numero de celulas                                    * 
 * iKelvin   -> temperatura em kelvin (true|false)                   *
 * iProp     -> numero da propriedade                                * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * prop    -> propriedade iniciacializada                            * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
  *********************************************************************/
void initPropTemp(PropVarFluid *propFluid
                 ,DOUBLE *RESTRICT prop   ,DOUBLE *RESTRICT t 
                 ,DOUBLE *RESTRICT pressure,DOUBLE *RESTRICT propMat
                 ,short *RESTRICT mat
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
        = airDensity(&propFluid->den,t[i], pressure[i]
                    ,thDynamic.pTh[2], iKelvin);
/*...................................................................*/

/*...*/
    else if( iProp == SPECIFICHEATCAPACITYFLUID)  
      MAT2D(lMat, iProp, propMat, MAXPROP) 
      = airSpecifiHeat(&propFluid->sHeat,t[i],iKelvin);
/*...................................................................*/

/*...*/
    else if( iProp == DYNAMICVISCOSITY)  
      MAT2D(lMat, iProp, propMat, MAXPROP) 
      = airDynamicViscosity(&propFluid->dVisc,t[i],iKelvin);
/*...................................................................*/

/*...*/
    else if( iProp == THERMALCONDUCTIVITY)  
      MAT2D(lMat,iProp,propMat,MAXPROP) 
      = airThermalConductvity(&propFluid->thCond,t[i],iKelvin);
/*...................................................................*/

/*...*/
    for(j=0;j<np;j++)      
      MAT2D(i,j,prop,np) = MAT2D(lMat,iProp,propMat,MAXPROP); 
/*...................................................................*/
  
  }
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 12/05/2018                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* initPropCD: inicializao de propriedades com variacao de           *
* temperatura                                                       *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* pol     -> polinomio de baixa ordem                               *
* prop    -> nao definido                                           *
* u       -> temperatura                                            *
* propMat -> propriedade de referencia por material                 *
* mat     -> material por celula                                    *
* np      -> numero niveis de tempos                                *
* nCell   -> numero de celulas                                      *
* iProp   -> numero da propriedade                                  *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* prop    -> propriedade iniciacializada                            *
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
void initPropCD(PropPol *pol            , DOUBLE *RESTRICT prop   
              , DOUBLE *RESTRICT u      , DOUBLE *RESTRICT propMat
              , short *RESTRICT mat
              , short np                , INT    nCell
              , short iProp)
{
  INT i;
  unsigned short j, lMat;
/*...*/
  for (i = 0; i<nCell; i++) 
  {
/*...*/
    lMat = mat[i] - 1;
/*...................................................................*/

/*...*/
    MAT2D(lMat, iProp, propMat, MAXPROP) = diffProp(pol, u[i]);
/*...................................................................*/

/*...*/
    for (j = 0; j<np; j++)
      MAT2D(i, j, prop, np) = MAT2D(lMat, iProp, propMat, MAXPROP);
/*...................................................................*/
  }
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 04/09/2017                                   *
 * Data de modificaco : 07/05/2019                                   *
 *-------------------------------------------------------------------*
 * INITSHEATPOL: inicializao a estrutura para o calculo do calor     *
 * especifico em funcao da temperatura via polinomio                 *
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
void initSheatPol(PropPol *prop, char *s, FILE *file) {

  FILE *fileOut;
  char nameAux[1000];
  short i;
  double x[MAXPLODEG];

  if (!strcmp(s, "polinomio")) {
    prop->type = POL;
    fscanf(file, "%s", nameAux);
    fileOut = openFile(nameAux, "r");

    prop->nPol[0] = readFileLineSimple(x, fileOut);
    ERRO_POL_READ(prop->nPol[0], MAXPLODEG, __FILE__, __func__, __LINE__);

    for (i = 0; i < prop->nPol[0]; i++)
      prop->a[i] = x[i];

    fclose(fileOut);

    fprintf(fileLogExc, "%-25s: %s\n", "Type", s);

  }

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 16/09/2017                                   *
 * Data de modificaco : 07/05/2019                                   *
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
void initDviscosityPol(PropPol *prop, char *s, FILE *file) {

  FILE *fileOut;
  char nameAux[1000];
  short i;
  double x[MAXPLODEG];

  if (!strcmp(s, "polinomio")) {
    prop->type = POL;
    fscanf(file, "%s", nameAux);
    fileOut = openFile(nameAux, "r");

    prop->nPol[0] = readFileLineSimple(x, fileOut);
    ERRO_POL_READ(prop->nPol[0], MAXPLODEG, __FILE__, __func__, __LINE__);

    for (i = 0; i < prop->nPol[0]; i++)
      prop->a[i] = x[i];

    fclose(fileOut);

    fprintf(fileLogExc, "%-25s: %s\n", "Type", s);

  }

  else if(!strcmp(s,"sutherland")){
    prop->type = SUTHERLAND;
    prop->a[0] = 1.789e-05; /*viscosidade de referencia*/
    prop->a[1] = 273.11e0;  /*temperatura de referencia*/
    prop->a[2] = 110.56e0;  /*constante de Sutherland*/
  }
  else {
    ERRO_GERAL(__FILE__,__func__,__LINE__,s);
  }
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 05/11/2017                                   *
 * Data de modificaco : 07/05/2019                                   *
 *-------------------------------------------------------------------*
 * INITDENSITY: inicializao a estrutura para o calculo da            *
 * densidade via polinomio                                           *
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
void initDensityPol(PropPol *prop, char *s, FILE *file) {

  FILE *fileOut;
  char nameAux[1000];
  short i;
  double x[MAXPLODEG];

  if(!strcmp(s,"polinomio")){
    prop->type = POL;
    fscanf(file, "%s", nameAux);
    fileOut = openFile(nameAux, "r");

    prop->nPol[0] = readFileLineSimple(x, fileOut);
    ERRO_POL_READ(prop->nPol[0], MAXPLODEG, __FILE__, __func__, __LINE__);

    for (i = 0; i < prop->nPol[0]; i++)
      prop->a[i] = x[i];

    fclose(fileOut);

    fprintf(fileLogExc, "%-25s: %s\n", "Type", s);

  }
/*... ideal gas (p)*/
  else if(!strcmp(s,"idealgas")){
    prop->type = IDEALGAS;
    fprintf(fileLogExc,"%-25s: %s\n","Type",s);
  }
/*... ideal gas incompressivel (PRef)*/
  else if(!strcmp(s,"incidealgas")){
    prop->type = INCIDEALGAS;
    fprintf(fileLogExc,"%-25s: %s\n","Type",s);
  }
  else {
    ERRO_GERAL(__FILE__,__func__,__LINE__,s);
  }
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 05/11/2017                                   *
* Data de modificaco : 07/05/2019                                   *
*-------------------------------------------------------------------*
* initCdPol: inicializao a estrutura para o calculo da              *
* propriedade via polinomio                                         *
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
void initCdPol(PropPol *prop,char *s,FILE *file)
{

  FILE *fileOut;  
  char nameAux[1000];
  short i;
  double x[MAXPLODEG];

  if (!strcmp(s, "polinomio")) 
  {
    prop->type = POL;
    fscanf(file, "%s", nameAux);
    fileOut = openFile(nameAux, "r");

    prop->nPol[0] = readFileLineSimple(x,fileOut);
    ERRO_POL_READ(prop->nPol[0], MAXPLODEG, __FILE__,__func__, __LINE__);

    for (i = 0; i < prop->nPol[0]; i++)
      prop->a[i] = x[i];

    fclose(fileOut);

    fprintf(fileLogExc, "%-25s: %s\n", "Type", s);
  }
  else 
  {
    ERRO_GERAL(__FILE__, __func__, __LINE__, s);
  }
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 16/09/2017                                   *
 * Data de modificaco : 07/05/2019                                   *
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
void initThCondPol(PropPol *prop, char *s, FILE *file) {

  FILE *fileOut;
  char nameAux[1000];
  short i;
  double x[MAXPLODEG];

  if(!strcmp(s,"polinomio"))
  {

    prop->type = POL;
    fscanf(file, "%s", nameAux);
    fileOut = openFile(nameAux, "r");

    prop->nPol[0] = readFileLineSimple(x, fileOut);
    ERRO_POL_READ(prop->nPol[0], MAXPLODEG, __FILE__, __func__, __LINE__);

    for (i = 0; i < prop->nPol[0]; i++)
      prop->a[i] = x[i];

    fclose(fileOut);

    fprintf(fileLogExc, "%-25s: %s\n", "Type", s);
  }
  else if(!strcmp(s,"sutherland"))
  {
    prop->type = SUTHERLAND;
    prop->a[0] = 1.789e-05; /*viscosidade de referencia*/
    prop->a[1] = 273.11e0;  /*temperatura de referencia*/
    prop->a[2] = 110.56e0;  /*constante de Sutherland*/
  }
  else {
    ERRO_GERAL(__FILE__,__func__,__LINE__,s);
  }
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 06/09/2017                                   *
 * Data de modificaco : 15/07/2018                                   *
 *-------------------------------------------------------------------*
 * GETTEMPFORENERGY: obtem a temperatura aprtir da entalpia sensivel *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * sHeatPol - estrutra par o polinoimio do calor especifico          *
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
void getTempForEnergy(PropPol *sHeatPol
                     ,DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                     ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                     ,INT const nCell      ,bool const fTemp
                     ,bool const fSheat    ,bool const fKelvin
                     ,bool const fOmp      ,short const nThreads ){
  
  short lMat;
  INT i;  
  DOUBLE sHeatRef;

/*... resolucao da eq da energia na forma de temperatura*/ 
  if(fTemp)
    for (i = 0; i < nCell; i++)
      temp[i] = energy[i]; 
/*...................................................................*/ 

/*... resolucao da eq da energia na forma de entalpia sensivel*/  
  else{
    if(fOmp){
#pragma omp parallel  for default(none) num_threads(nThreads)\
      private(i,lMat,sHeatRef) shared(prop,mat,energy,temp,sHeatPol)
      for (i = 0; i < nCell; i++) {
        lMat  = mat[i] - 1;
        sHeatRef = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
        temp[i] = specificEnthalpyForTemp(sHeatPol
                                         , energy[i], sHeatRef
                                         , fSheat  , fKelvin);
      }
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    else{
      for (i = 0; i < nCell; i++) {
        lMat  = mat[i] - 1;
        sHeatRef = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
        temp[i] = specificEnthalpyForTemp(sHeatPol
                                         , energy[i], sHeatRef
                                         , fSheat   , fKelvin);
      }
/*...................................................................*/
    }
/*...................................................................*/ 
  }
/*...................................................................*/ 

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 20/08/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * etEnergyForTemp   : obtem a entalpia sensivel apartir da temp     *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * sHeatPol - estrutra par o polinoimio do calor especifico          *
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
void getEnergyForTemp(PropPol *sHeatPol
                     ,DOUBLE *RESTRICT temp,DOUBLE *RESTRICT energy
                     ,DOUBLE *RESTRICT prop,short  *RESTRICT mat 
                     ,INT const nCell     
                     ,bool const fSheat    ,bool const fKelvin
                     ,bool const fOmp      ,short const nThreads ) 
{
  short lMat;
  INT i;  
  DOUBLE sHeatRef;

/*...*/
  if (fOmp) 
  {
/*...*/ 
#pragma omp parallel  for default(none) num_threads(nThreads)\
    private(i,lMat,sHeatRef) shared(prop,mat,energy,temp,sHeatPol)
    for (i = 0; i < nCell; i++) 
    {
      lMat  = mat[i] - 1;
      sHeatRef = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
      energy[i] = tempForSpecificEnthalpy(sHeatPol
                                         ,temp[i]  ,sHeatRef
                                         ,fSheat  ,fKelvin);
    }
/*...................................................................*/ 
  }
/*...................................................................*/ 

  else
  {
/*...*/ 
    for (i = 0; i < nCell; i++) 
    {
      lMat  = mat[i] - 1;
      sHeatRef = MAT2D(lMat, SPECIFICHEATCAPACITYFLUID, prop, MAXPROP);
      energy[i] = tempForSpecificEnthalpy(sHeatPol
                                         ,temp[i] ,sHeatRef
                                         ,fSheat  ,fKelvin);
    }
/*...................................................................*/ 
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
 * Data de modificaco : 26/08/2018                                   *
 *-------------------------------------------------------------------*
 * INIPRESREF : incializa a pressao ref atrazes da massa especifica  * 
 * de referencia e temperatura media do dominio                      *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * temp   - temperatura no passo n + 1                               *
 * volume - volume das celulas                                       *
 * pTh    - pressao termodinamica de referencia                      *
 * prop                                                              *
 * mat                                                               *
 * molarMass                                                         *
 * nCell  - numero da celulas                                        *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------*
 * pTh    - pressao termodinamica de referencia atualizada           *
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 *-------------------------------------------------------------------*
 *********************************************************************/
void initPresRef(DOUBLE *RESTRICT temp  , DOUBLE *RESTRICT volume
               , DOUBLE *pTh            , DOUBLE *RESTRICT prop  
               , short  *RESTRICT mat   , DOUBLE const molarMass                  
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
  vm   = PRESREF(dRef, IDEALGASR, tm, molarMass);

  pTh[0] = pTh[1] = pTh[2] = vm;

}
/*********************************************************************/

/**********************************************************************
 * Data de criacao    : 12/05/2018                                    *
 * Data de modificaco : 00/00/0000                                    *
 *--------------------------------------------------------------------*
 * setReGrad:                                                         *
 *--------------------------------------------------------------------*
 * Parametros de entrada:                                             *
 *--------------------------------------------------------------------*
 * x       -> nao definido                                            *
 * file    -> arquivo de arquivo                                      *
 *--------------------------------------------------------------------*
 * Parametros de saida:                                               *
 *--------------------------------------------------------------------*
 * x       -> valores lidos                                           *
 * nTerms  -> retornas o numero de valores lidos                      *
 *--------------------------------------------------------------------*
 * OBS:                                                               *
 *--------------------------------------------------------------------*
 **********************************************************************/
int readFileLineSimple(DOUBLE *a, FILE *file)
{

  int i, nTerms;

  fscanf(file, "%d", &nTerms);

  for (i = 0; i<nTerms; i++)
    fscanf(file, "%lf", a + i);

  return nTerms;

}
/**********************************************************************/

/**********************************************************************
 * Data de criacao    : 02/06/2018                                    *
 * Data de modificaco : 00/00/0000                                    *
 *--------------------------------------------------------------------*
 * initPropCD : inicializacao                                         *
 *--------------------------------------------------------------------*
 * Parametros de entrada:                                             *
 *--------------------------------------------------------------------*
 * propVar -> nao definido                                            *
 * n       -> numero de termors                                       *
 *--------------------------------------------------------------------*
 * Parametros de saida:                                               *
 *--------------------------------------------------------------------*
 * propVar -> inicializado                                            *
 *--------------------------------------------------------------------*
 * OBS:                                                               *
 *--------------------------------------------------------------------*
 **********************************************************************/
void initPropStructCD(PropVarCD *propVar, short const n)
{
  short i;

  for(i=0;i<n;i++)
  {
    propVar[i].fDensity        = false;
    propVar[i].fCeofDiff       = false;
    propVar[i].ceofDiff.type   = -1;
    propVar[i].ceofDiff.nPol[0] = 0;
    propVar[i].den.nPol[0]     = 0;
    propVar[i].den.type        = -1;
  }

}
/**********************************************************************/

