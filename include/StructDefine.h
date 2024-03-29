#ifndef _STRUCTDEFINEH_
  #define _STRUCTDEFINEH_

  #include<HccaStdBool.h>
  #include<Define.h>
  #include<File.h>

  extern INT gStep;

/*...*/
  typedef struct
  {
    DOUBLE *t00,*t0,*t;
  } LevelTime;
/*...................................................................*/

/*...*/
  typedef struct
  {
    bool flWord;
    short  jLoop,kLoop;
    char loopWord[100][MAX_LINE];
  } Macros;
/*...................................................................*/

/*...*/
  typedef struct
  {
    bool fLimiter;
    short  type,func;
    DOUBLE beta;
  } RcGrad;
/*...................................................................*/

/*...*/
  typedef struct
  {
    short  type;
    INT maxIt;
    DOUBLE tol;
  } Edo;
/*...................................................................*/

/*...*/
  typedef struct
  {
    DOUBLE A,E,Ta,beta;
  } Arrhenius;
/*...................................................................*/

/*...*/
  typedef struct
  {
    bool reverse;
    short nPartSp[2];
    short partSp[2][MAXSPECIES];
    DOUBLE stch[3][MAXSPECIES];
    DOUBLE exp[2][MAXSPECIES];
    DOUBLE sO2;
    Arrhenius ArrF,ArrR;
  } Reaction;
/*...................................................................*/

/*...*/
  typedef struct
  {
    char  name[MAXNAMELENSP];
    unsigned short nO,nN,nC,nH;
    DOUBLE mW;
    DOUBLE leornadJones[2]; /*col 1 - sigma, col 2 -e/k*/
/*... entalpia de formacao*/
    DOUBLE entalphyOfForm;       /*kj/kg*/
    DOUBLE entalphyOfFormMolar;  /*kj/kmol*/

  } Specie;
/*...................................................................*/

/*...*/
  typedef struct
  {
    unsigned short nReac,nSp,nEp;
    unsigned short eO,eN,eC,eH;
    unsigned short sCO2,sCO,sH2O,sO2,sN2,sCH4,sC3H8;
    short fuel[MAXSPECIES],nFuel;
    short ox[MAXSPECIES],nOx;
    short prod[MAXSPECIES],nProd;
    DOUBLE mE[MAXELEMENT];
    Specie sp[MAXSPECIES];
    Reaction reac[MAXREAC];
  } Chemical;
/*...................................................................*/

/*...*/
  typedef struct{
    short type;
    DOUBLE cGamma,cTau,tMix; /* KJ/Kmol*/
    Edo edo;
  }Edc;
/*...................................................................*/

/*...*/
  typedef struct{
    bool tMixConst,fProd;
    DOUBLE coef[3]; /* 0 - A
                       1 - B
                       2 - tMix*/
  }Edm;
/*...................................................................*/

/*...*/
  typedef struct{
    bool fDensityRef;
    bool fPresTh;
    DOUBLE pTh[3];
  }ThermoDynamic;
  extern ThermoDynamic thDynamic;
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
     DOUBLE unit;
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
  extern Interpol iPol[MAXINTERPOL];
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

    DOUBLE lumpedMatrix[MAXSPECIES*3];

    DOUBLE entalphyOfCombustion;        /* Entalpia de combustao calculada
                                          pelas especies primitivas*/
    DOUBLE entalphyOfCombustionGrouped;  /* Entalpia de combustao calculada
                                          pelas especies agrupadas*/
    DOUBLE totalHeat;     /* Calor total liberado pela reacao de combustao*/
/*... composicao do ar*/
    DOUBLE O2InAir,N2InAir;

    DOUBLE CO2InProd,H2OInProd,N2InProd;

    Edc      edc;
    Edm      edm;
    Chemical chem;
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
    bool fOneEq;
    bool fTurbStruct;

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
    bool fDiffEnergy;
  }EnergyModel;
/*...................................................................*/

/*...*/
  typedef struct{
    bool RhsDensity,LhsDensity;
  }MassEqModel;
/*...................................................................*/

/*...*/
  typedef struct{
    bool fRes;
    bool fRhieChowInt;
    bool fViscosity;
    bool fDiv;
    bool fSoPressure;
    short iCodBuoyant;
  }MomentumModel;
/*...................................................................*/

/*...*/
  typedef struct{
    bool flag,fDynamic;
    short iCod,type,typeReal;
    INT timeStep;
    DOUBLE total,dt[3],dtInicial,t,t0,dtMax,dtMin,cfl,chem;
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
    bool fUse;
    short nTypeVar;               /*0 - constante
                                    1 - funcao*/
    short iCod[2];               /* conjunto de codigos*/
    short type;                     /*tipo*/
    short np;                       /*numero de particoes*/
    DOUBLE par[MAXLOADPARAMETER];
    DOUBLE vel[3],density;
    Interpol *intPol;               /*interpolacao*/
  }Loads;

  extern Loads
        loadsD1[MAXLOAD],      /*tipo de cargas (difusao pura)*/
        loadsT1[MAXLOAD],      /*tipo de cargas (difusao-transporte)*/
        loadsVel[MAXLOAD],     /*tipo de cargas (fluid-Vel)*/
        loadsPres[MAXLOAD],    /*tipo de cargas (fluid-Pres)*/
        loadsPresC[MAXLOAD],   /*tipo de cargas (fluid-Pres-correcao)*/
        loadsEnergy[MAXLOAD],  /*tipo de cargas (fluid-energia)*/
        loadsTemp[MAXLOAD],    /*tipo de cargas (fluid-temperatura)*/
        loadsKturb[MAXLOAD],   /*tipo de cargas (fluid-turbulencia)*/
        loadsZcomb[MAXLOAD],    /*tipo de cargas (combustivel      )*/
        loadsRhoFluid[MAXLOAD];
/*...................................................................*/

/*...*/
  typedef struct{
    RcGrad rcGrad; /*tipo de rescontrucao de gradiente*/
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
    bool fRel;
    unsigned short type;
    unsigned short k;         /*iteracao com o qual o residuo*/
    DOUBLE tol[3];
    char name[50];
  }Residual;
/*...................................................................*/

/*... Simple*/
  typedef struct{
    bool   sPressure;
    unsigned short faceInterpolVel;   /*tecnica de interpolacao
                                        das velocidades nas faces para
                                        evitar o problema checkboard
                                       */
    short  type;

    int    pSimple;
    int    maxIt;
    int nNonOrth;
    DOUBLE alphaPres,alphaVel,alphaEnergy
          ,alphaDensity,alphaComb;       /*under-relaxation*/
    DOUBLE *ePresC,*nPresC  ,*eGradPresC;/*Pressao de correcao*/
    DOUBLE *ePresC1;                     /*Pressao de correcao 1*/
    DOUBLE *d;
    Residual vel,z,mass,energy,kTurb;
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

/*... interpolacao da variaveis temporais*/
  typedef struct
  {
    DOUBLE *vel0,*vel, *velI  ,*velG, *nVelI,*nVelG;
    DOUBLE *gradVelI,*gradVelG,*nGradVelI, *nGradVelG;
    DOUBLE *p0  ,*p  , *pI    ,*pG, *nPresI,*nPresG;
    DOUBLE *gradPresI,*gradPresG,*nGradPresI,*nGradPresG;
    DOUBLE *temp0,*temp,*tempI,*tempG,*nTempI,*nTempG;
    DOUBLE *gradTempI,*gradTempG,*nGradTempI,*nGradTempG;
    DOUBLE *y0 ,*y,*yI,*yG,*nYI,*nYG;
    DOUBLE *gradYI,*gradYG,*nGradYI,*nGradYG;
    DOUBLE *wT0  ,*wT  ,*wTI  ,*wTG,*nWTI,*nWTG;
    DOUBLE *rho0,*rho,*rhoI,*rhoG,*nRhoI,*nRhoG;
    DOUBLE *sHeat0,*sHeat,*sHeatI,*sHeatG,*nsHeatI,*nsHeatG;
    DOUBLE *dVisc0,*dVisc,*dViscI,*dViscG,*ndViscI,*ndViscG;
    DOUBLE *tCond0,*tCond,*tCondI,*tCondG,*ntCondI,*ntCondG;
    DOUBLE *cDiff0,*cDiff,*cDiffI,*cDiffG,*ncDiffI,*ncDiffG;
    DOUBLE *mMolar0,*mMolar,*mMolarI,*mMolarG,*nmMolarI,*nmMolarG;
    DOUBLE *eddyVisc0,*eddyVisc,*eddyViscI,*eddyViscG,*neddyViscI,*neddyViscG;
  }TimeInterpol;
/*...................................................................*/

#endif /*_MESH_*/
