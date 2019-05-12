#ifndef _STRUCTDEFINEH_
  #define _STRUCTDEFINEH_

  #include<HccaStdBool.h>
  #include<Define.h>

/*...*/
  typedef struct{
    DOUBLE energyAtivation; /* KJ/Kmol*/
    DOUBLE alpha;
    DOUBLE a;
  }ArrheniusLaw;
/*...................................................................*/  

/*...*/
  typedef struct{
    bool fDensityRef;
    bool fPresTh;
    DOUBLE pTh[3];
  }ThermoDynamic;
   ThermoDynamic thDynamic;
/*...................................................................*/  

/*...*/
   typedef struct {
     unsigned char type;
     short nPol[MAXSPECIES];
     DOUBLE a[MAXSPECIES*MAXPLODEG];
   }PropPol;
/*...................................................................*/

/*...*/
   typedef struct {
     bool fDensity;
     bool fSpecificHeat;
     bool fDynamicViscosity;
     bool fThermalconductivity;
     PropPol den, thCond, dVisc, sHeat;
     DOUBLE molarMass; 
   }PropVarFluid;
/*...................................................................*/

/*...*/
   typedef struct {
     bool fDensity;
     bool fCeofDiff;
     PropPol den, ceofDiff;
   }PropVarCD;
/*...................................................................*/

/*... interpol*/
  typedef struct{
    unsigned short np;                       /*numero de particoes*/  
    DOUBLE *x,*y;                            /*interpolacao*/
  }Interpol;
  Interpol iPol[MAXINTERPOL];
/*...................................................................*/  

/*... interpol*/
  typedef struct{
    short c,h,o;
  }Fuel;
/*...................................................................*/


/*...*/
  typedef struct {
    bool fRes;
    bool fCombustion;
    bool fLump;
    short nOfSpecies;     /* numero total de especies*/
    short nOfSpeciesLump; /* numero de especies agrupadas*/
    short nComb;          /* numero especies transportadas*/
    short typeHeatRealese;
    short reactionKinetic;
    Fuel fuel;
    DOUBLE sMassAir ,sMolar,tMix;
    DOUBLE sMassO2  ,sMassN2;
    DOUBLE sMassCO2p,sMassH2Op,sMassN2p;
    DOUBLE stoichO2,stoichN2;
    DOUBLE stoichCO2p,stoichH2Op,stoichN2p; 
    DOUBLE lumpedMatrix[21];
/*... massa molar*/
    DOUBLE mW_Fuel,mW_N2,mW_O2,mW_CO2,mW_CO,mW_H2O,mW_C,mW_Air;
/*... entalpia de formacao*/    
    DOUBLE entalphyOfForm[MAXSPECIES]; /*0 - Fuel
                                         1 - O2   
                                         2 - N2
                                         3 - CO2     
                                         4 - H2O*/      

    DOUBLE entalphyOfCombustion;        /* Entalpia de combustao calculada 
                                          pelas especies primitivas*/
    DOUBLE entalphyOfCombustionGrouped;  /* Entalpia de combustao calculada 
                                          pelas especies agrupadas*/
    DOUBLE totalHeat;     /* Calor total liberado pela reacao de combustao*/
    DOUBLE totalMassFuel;  /* massa total de combustivel consumido*/
/*... composicao do ar*/
    DOUBLE O2InAir,N2InAir;

    DOUBLE CO2InProd,H2OInProd,N2InProd;

    ArrheniusLaw arrhenius; 

  } Combustion;
/*...................................................................*/

/*... Material*/
  typedef struct{
    short  *type;      /*tipo de calculo da celula*/ 
    double *prop;      /*valores numericos da propriedade*/    
  } Material;
/*...................................................................*/

/*...*/
  typedef struct{
    short *nViz;
    INT  *nelcon;
  }Adjacency;
/*...................................................................*/  

/*...*/
  typedef struct{
    bool fRes;
    INT maxIt;
    DOUBLE ck,ce,sk,tol;
  }EqK;
/*...................................................................*/  


/*...*/
  typedef struct{  
    bool fTurb;
    bool fWall;
    bool fDynamic;
    short wallType;
    short type;          /* 0 - LES*/
    short typeLes;       /* 0 - funcional
                            1 - estrutural
                            2 - misto*/ 
    short typeMixed[2];  /*[0] - estrutural - [1] funcional*/
    short typeDynamic;    /* 1 - um paremetro local 
                            2 - um parametro global padrao
                            3 - um parametro global modicado
                            4 - 2 paramento local*/

    DOUBLE cs,cf,c; /*constante    
                           1 -        
                           */   
    DOUBLE PrandltTwall;   /*Prandtl turbulento */
    DOUBLE PrandltTsgs;    /*Prandtl de sub-grid */

    EqK eK;
    
  }Turbulence; 
/*...................................................................*/

/*...*/
  typedef struct{
    bool fPresWork;
    bool fDissipation;
    bool fRes;
    bool fTemperature;
    bool fKelvin;
  }EnergyModel;
/*...................................................................*/  

/*...*/
  typedef struct{
    bool RhsDensity,LhsDensity;
  }MassEqModel;
/*...................................................................*/

/*...*/
  typedef struct{    
    bool fRes,fRhieChowInt,fViscosity,fDiv; 
    short iCodBuoyant;
  }MomentumModel;
/*...................................................................*/

/*...*/
  typedef struct{
    bool flag,fDynamic;
    short iCod,type;
    INT timeStep;
    DOUBLE total,dt[3],dtInicial,t;   
  }Temporal;
/*...................................................................*/  

/*...*/
  typedef struct{
    bool fRes;
    DOUBLE cEq[3];
  }kModel;
/*...................................................................*/  

/*...*/
  typedef struct{
    bool fRes;
  }DiffModel;
/*...................................................................*/

/*...*/
  typedef struct {
    bool fRes;
  }TransModel;
/*...................................................................*/

/*...*/
  typedef struct{
    short  iCod1;
    short  iCod2;
    DOUBLE par[NPADV];
  }Advection;
/*...................................................................*/  

/*...*/
  typedef struct{
    short  iCod;
  }Diffusion;
/*...................................................................*/ 

/*...*/
  typedef struct{
    int maxIt;
    int pPlot;
    DOUBLE tol;                               
  }NonLinear;
/*...................................................................*/  

/*...*/
  typedef struct{
    bool fInit;
    bool fMedia;
    bool fVel;
    int startSample,endSample;
    DOUBLE *mVel,*sVel;  
    DOUBLE t0;                             
  }Mean;
/*...................................................................*/  

/*... loads*/
  typedef struct{
    short nTypeVar;               /*0 - constante
                                    1 - funcao parabolica*/  
    short type;                     /*tipo*/
    short np;                       /*numero de particoes*/  
    DOUBLE par[MAXLOADPARAMETER];
    DOUBLE vel[3],density;
    Interpol *intPol;               /*interpolacao*/
  }Loads;
  Loads  loadsD1[MAXLOADD1]         /*tipo de cargas (difusao pura)*/
        ,loadsT1[MAXLOADT1]         /*tipo de cargas (difusao-transporte)*/
        ,loadsVel[MAXLOADFLUID]     /*tipo de cargas (fluid-Vel)*/
        ,loadsPres[MAXLOADFLUID]    /*tipo de cargas (fluid-Pres)*/
        ,loadsPresC[MAXLOADFLUID]   /*tipo de cargas (fluid-Pres-correcao)*/
        ,loadsEnergy[MAXLOADFLUID]  /*tipo de cargas (fluid-energia)*/
        ,loadsTemp[MAXLOADFLUID]    /*tipo de cargas (fluid-temperatura)*/
        ,loadsKturb[MAXLOADFLUID]   /*tipo de cargas (fluid-turbulencia)*/
        ,loadsZcomb[MAXLOADFLUID];   /*tipo de cargas (combustivel      )*/
/*...................................................................*/

/*...*/
  typedef struct{
    short rcGrad; /*tipo de rescontrucao de gradiente*/                     
/*...*/
    Temporal ddt;
/*...*/
    NonLinear nlTemp;
    NonLinear nlD1;
    NonLinear nlT1;
/*... equacao de difusao*/
    Diffusion diffD1;
/*... equacao de transporte*/
    Advection advT1;
    Diffusion diffT1;
/*... equacao de velocidade*/
    Advection advVel;
    Diffusion diffVel;
/*... equacao de energia*/
    Advection advEnergy;
    Diffusion diffEnergy;
/*... equacao de energia*/
    Advection advKturb;
    Diffusion diffKturb;
/*... equacao de pressao*/
    Diffusion diffPres;
/*... equacao de fracao massica*/
    Advection advComb;
    Diffusion diffComb;
  }Scheme;
/*...................................................................*/

/*... Simple*/
  typedef struct{
    bool   sPressure;
    unsigned short faceInterpolVel;   /*tecnica de interpolacao
                                        das velocidades nas faces para
                                        evitar o problema checkboard
                                       */
    short  type;
 
    unsigned short kZeroVel;          /*iteracao com o qual o residuo
                                        e normalizado*/
    unsigned short kZeroPres;         /*iteracao com o qual o residuo
                                        e normalizado*/
    unsigned short kZeroEnergy;       /*iteracao com o qual o residuo
                                     e normalizado*/
    unsigned short kZeroComb;         /*iteracao com o qual o residuo
                                       e normalizado*/
    int    pSimple;
    int    maxIt;
    int nNonOrth;
    DOUBLE alphaPres,alphaVel,alphaEnergy
          ,alphaDensity,alphaComb;       /*under-relaxation*/   
    DOUBLE *ePresC,*nPresC  ,*eGradPresC;/*Pressao de correcao*/
    DOUBLE *ePresC1;                     /*Pressao de correcao 1*/
    DOUBLE *d;
    DOUBLE tolPres,tolVel[3],tolEnergy,tolComb;
  }Simple;
/*...................................................................*/

/*... Prime*/
  typedef struct {
    bool   sPressure;
    unsigned short faceInterpolVel;   /*tecnica de interpolacao
                                      das velocidades nas faces para
                                      evitar o problema checkboard
                                     */
    short  type;
    unsigned short kZeroVel;          /*iteracao com o qual o residuo
                                      e normalizado*/
    unsigned short kZeroPres;         /*iteracao com o qual o residuo
                                      e normalizado*/
    int nNonOrth;
    int    pPrime;          /*iteracao que serao imprisas na tela*/ 
    int    maxIt;
    DOUBLE alphaPres, alphaVel;           /*under-relaxation*/
    DOUBLE *ePresC, *nPresC, *eGradPresC;/*Pressao de correcao*/
    DOUBLE *ePresC1;                     /*Pressao de correcao 1*/
    DOUBLE *d,*velUp,*aD,*bTemporal;
    DOUBLE tolPres, tolVel;
  }Prime;
/*...................................................................*/
#endif /*_MESH_*/