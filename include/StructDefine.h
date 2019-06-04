#ifndef _STRUCTDEFINEH_
  #define _STRUCTDEFINEH_

  #include<HccaStdBool.h>
  #include<Define.h>

  DOUBLE gStep;

/*...*/
  typedef struct{
    DOUBLE energyAtivation; /* KJ/Kmol*/
    DOUBLE alpha;
    DOUBLE a,e1,e2;
  }ArrheniusLaw;
/*...................................................................*/  

/*...*/
  typedef struct{
    short type;
    DOUBLE cGamma,cTau,tMix; /* KJ/Kmol*/    
  }Edc;
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
     short nPol;
     DOUBLE range[2];
     DOUBLE a[MAXPLODEG];
   }Pol;
/*...................................................................*/

/*...*/
   typedef struct {
     short type;
     DOUBLE range[3][2];
     DOUBLE a1[2][9],a2[2][9];  /*(KJ/KmolK) (KJ/KGK)*/
   }PolNasa;
/*...................................................................*/

/*...*/
   typedef struct {
     unsigned char type;
     Pol     pol[MAXSPECIES];  
     PolNasa nasa[MAXSPECIES];
     DOUBLE surtherland[3];
   }Prop;
/*...................................................................*/

/*...*/
   typedef struct {
     bool fDensity;
     bool fSpecificHeat;
     bool fDynamicViscosity;
     bool fThermalConductivity;
     bool fDiffusion;
     Prop den, thCond, dVisc, sHeat, diff;
     DOUBLE molarMass;
     DOUBLE sHeatRef,dViscosityRef,ThermalConductivityRef,densityRef; 
   }PropVarFluid;
/*...................................................................*/

/*...*/
   typedef struct {
     bool fDensity;
     bool fCeofDiff;
     Prop den, ceofDiff;
   }PropVarCD;
/*...................................................................*/

/*... interpol*/
  typedef struct{
    unsigned short np;                       /*numero de particoes*/  
    DOUBLE *x,*y;                            /*interpolacao*/
  }Interpol;
  Interpol iPol[MAXINTERPOL];
/*...................................................................*/  

/*...*/
  typedef struct{
    short c,h,o;
    char name[20];
    DOUBLE hf,lj[2];
  }Fuel;
/*...................................................................*/

/*...*/
  typedef struct {
    bool fRes;
    bool fCombustion;
    bool fLump;
    bool fCorrectVel;
    short nOfSpecies;     /* numero total de especies*/
    short nOfSpeciesLump; /* numero de especies agrupadas*/
    short nComb;          /* numero especies transportadas*/
    short typeHeatRealese;
    short reactionKinetic;
    short nReac; 
    short sp_fuel[2],sp_CO2,sp_H2O,sp_O2,sp_N2;
    short speciesPart[2][2][MAXSPECIES],nSpeciesPart[2][2]; /*especie que participam da reacao i*/
    DOUBLE sMass[2][2][MAXSPECIES],sMassAir;    
    DOUBLE stoich[2][2][MAXSPECIES],stoichAir;/*(i,j,k) 
                                              i - reacao quimica
                                              j - 0 reagente - 1 produtos
                                              k - especies de 0 a N - 1
                                              */
           
    DOUBLE lumpedMatrix[MAXSPECIES*3];
/*... massa molar*/
    DOUBLE mW[MAXSPECIES],mW_Air;
/*... Leornad-Jone parametros*/
    DOUBLE leornadJones[MAXSPECIES][2]; /*col 1 - sigma, col 2 -e/k
/*... entalpia de formacao*/    
    DOUBLE entalphyOfForm[MAXSPECIES]; /*0   - Fuel
                                         1   - O2   
                                         3   - CO2     
                                         4   - H2O
                                         5   - CO
                                         6   - ... 
                                         N-1 - N2*/      

    DOUBLE entalphyOfCombustion;        /* Entalpia de combustao calculada 
                                          pelas especies primitivas*/
    DOUBLE entalphyOfCombustionGrouped;  /* Entalpia de combustao calculada 
                                          pelas especies agrupadas*/
    DOUBLE totalHeat;     /* Calor total liberado pela reacao de combustao*/
    DOUBLE totalMassFuel;  /* massa total de combustivel consumido*/
/*... composicao do ar*/
    DOUBLE O2InAir,N2InAir;

    DOUBLE CO2InProd,H2OInProd,N2InProd;

    ArrheniusLaw arrhenius[2]; 
    Edc          edc;
    Fuel         fuel[2];
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
    DOUBLE SchmidtTsgs;    /*Chmidt  de sub-grid */

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