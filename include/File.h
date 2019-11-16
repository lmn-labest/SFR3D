#ifndef _FILE_H_
  #define _FILE_H_
/*...*/
  #include<math.h>
  #include<stdio.h>
  #include<stdlib.h>
  #include<ctype.h>
  #include<string.h>
/*...*/
  #include<Erro.h>
  #include<HccaStdBool.h>  
  #include<Define.h>
/*...*/
  #define MAX_STR_LEN_IN    129 /*tamanho do arquivo de entrada maximo*/
  #define MAX_STR_LEN_SUFIXO 80  /*tamanho do arquivo de entrada maximo*/
  #define SIZEMAX           200 /*tamanho do arquivo de saida maximo*/
  #define MAX_STR_NUMBER     65 /*INT integer*/
  #define MAX_EXT            50 /*numero maximo da extencao*/ 
  #define WORD_SIZE         120
  #define LINE_SIZE         120
  #define MAX_LINE          200

  #define FITPLOTD1           1
  #define FITPLOTD2           2
  #define FITPLOTT1           3
  #define FITPLOTT2           4
  #define FITPLOTTEMP         5
  #define FITPLOTSIMPLE       6
  
/*...................................................................*/  
  typedef struct{
    bool fPolt;
    bool bVtk;           /*escrita do arquivo binario*/
    bool fItPlotRes; 
    bool fItPlot;
    bool fCell;
    bool fNode;        
    bool gradPres;       
    bool gradVel;        
    bool gradEnergy; 
    bool gradTemp; 
    bool graduD1; 
    bool graduT1;
    bool gradZcomb;
    bool uD1;
    bool uT1;
    bool vel;            
    bool pres;
    bool presTotal;           
    bool energy;   
    bool temp;      
    bool eddyViscosity;  
    bool densityFluid;   
    bool specificHeat;   
    bool dViscosity;     
    bool tConductivity; 
    bool densityD1;
    bool coefDiffD1;
    bool densityT1;
    bool coefDiffT1;
    bool coefDiffSp;
    bool vorticity;
    bool wallParameters;
    bool stress; 
    bool kinetic;
    bool stressR;  
    bool cDynamic; 
    bool Qcriterion; 
    bool kTurb; 
    bool zComb;
    bool wk;         
    bool yFrac;
    bool rateHeatComb; 
    bool enthalpyk;
    bool gradY;
    bool tReactor;
    bool gradRho;
    bool mMolar;

    bool bconditions;    /*insere as condicoes de contorno nos valores nodais*/
    bool cc;             /*centro geomentrico da celula*/       
    bool pKelvin;        /*plot em kelvin*/  


    short timeFile;
    short stepPlot[2],nextStepPlot[2];
    bool fStepPlot;

    DOUBLE t,ta,tNext;
    bool fTimePlot;

    FILE *fileItPlot[8];
    FILE *fileParameters;
    
  }FileOpt;

  void  iota(INT t, char* st);
  FILE* openFile(const char* const name, const char* const mod);
  FILE* openFileBuffer(const char* const name, const char* const mod, bool buffer);
  void fName(const char* const name, INT num1, INT num2, int cod, char *out);
  void readMacro(FILE* file, char *mc, bool allline);
  void readMacroV2(FILE* file, char *mc, bool allline, bool lower);
  void  clearLine(char *s);
  int   rl(FILE *f, char *st);
  void convStringLower(char *s);
//int   getNumProp(FILE *f);
/*...*/  
  char macros[MAX_LINE][WORD_SIZE];/*todas as macros lidas no arquivo*/
  int  nmacro;
/*....................................................................*/  
#endif /*_FILE_H*/
