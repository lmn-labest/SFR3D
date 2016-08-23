#ifndef _CELLLOOP_H_
  #define _CELLLOOP_H_
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<Assbly.h>
  #include<Define.h>
  #include<HccaBlas.h>
  #include<HccaStdBool.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<ParallelMpi.h>
  #include<Loads.h>
/*...................................................................*/

/*...*/
  void convTempForKelvin(DOUBLE *restrict uC,INT const n
                        ,bool const fKelvin);
/*...................................................................*/

/* ... calculo das propriedade geometicas*/
  void pGeomForm(DOUBLE *restrict x    ,INT    *restrict el
              ,INT    *restrict nelcon ,short  *restrict nen    
              ,short  *restrict nFace  ,short  *restrict geomType
              ,DOUBLE *restrict gCc    ,DOUBLE *restrict gKsi  
              ,DOUBLE *restrict gmKsi  ,DOUBLE *restrict gEta 
              ,DOUBLE *restrict gmEta  ,DOUBLE *restrict gNormal
              ,DOUBLE *restrict gVolume,DOUBLE *restrict gXm    
              ,DOUBLE *restrict gXmcc  
              ,DOUBLE *restrict gvSkew ,DOUBLE *restrict gmSkew
              ,DOUBLE *restrict gDcca               
              ,short maxNo             ,short maxViz
              ,short ndm               ,INT numel);
/*...................................................................*/

/*.. reconstrucao de gradiente*/
  void rcGradU(Memoria *m             ,Loads *loads
         ,INT    *restrict el         ,INT    *restrict nelcon
         ,DOUBLE *restrict cc         ,DOUBLE *restrict x     
         ,short  *restrict nen        ,short  *restrict nFace
         ,short  *restrict geomType   ,DOUBLE *restrict prop 
         ,short  *restrict mat 
         ,DOUBLE *restrict leastSquare,DOUBLE *restrict leastSquareR   
         ,DOUBLE *restrict gKsi       ,DOUBLE *restrict gmKsi 
         ,DOUBLE *restrict gEta       ,DOUBLE *restrict gmEta 
         ,DOUBLE *restrict gNormal    ,DOUBLE *restrict gVolume
         ,DOUBLE *restrict gvSkew     
         ,DOUBLE *restrict gXm        ,DOUBLE *restrict gXmcc
         ,DOUBLE *restrict gDcca 
         ,short  *restrict faceR      ,short *restrict faceL  
         ,DOUBLE *restrict u          ,DOUBLE *restrict gradU               
         ,DOUBLE *restrict nU         ,short const lib
         ,short const maxNo           ,short const maxViz
         ,short const ndf             ,short const ndm
         ,InterfaceNo *iNo            ,Interface *iCel 
         ,INT const numelNov          ,INT const numel 
         ,INT const nNodeNov          ,INT const nNode);
/*...................................................................*/

/*... */
  void interCellNode(Memoria *m           ,Loads *loads
                   ,DOUBLE *restrict noU  ,DOUBLE *restrict elU
                   ,INT *restrict el      ,short  *restrict geomType 
                   ,DOUBLE *restrict cc   ,DOUBLE *restrict x
                   ,DOUBLE *restrict xm   
                   ,short *restrict nen   ,short *restrict nFace
                   ,short  *restrict faceR,short *restrict faceL  
                   ,InterfaceNo *iNo        
                   ,INT const numelNov    ,INT const numel       
                   ,INT const nNodeNov    ,INT const nNode 
                   ,short const maxNo     ,short const maxViz    
                   ,short const ndf1      ,short const ndf2  
                   ,short const ndm
                   ,bool const fBc        ,short const type);
/*...................................................................*/

/* ... montagem do sistemas de equacoes (difusao)*/
 void systFormDif(Loads *loads
             ,INT    *restrict el  ,INT    *restrict nelcon 
             ,short  *restrict nen     ,short  *restrict nFace
             ,short  *restrict geomType,DOUBLE *restrict prop 
             ,short  *restrict calType ,short  *restrict mat  
             ,DOUBLE *restrict gKsi   ,DOUBLE *restrict gmKsi
             ,DOUBLE *restrict gEta   ,DOUBLE *restrict gfArea
             ,DOUBLE *restrict gNormal,DOUBLE *restrict gVolume
             ,DOUBLE *restrict gXm    ,DOUBLE *restrict gXmcc
             ,DOUBLE *restrict gvSkew ,DOUBLE *restrict gmvSkew
             ,DOUBLE *restrict gDcca  ,DOUBLE *restrict density 
             ,INT    *restrict ia     ,INT    *restrict ja
             ,DOUBLE *restrict al     ,DOUBLE *restrict ad
             ,DOUBLE *restrict b      ,INT    *restrict id
             ,short  *restrict faceR  ,short  *restrict faceL
             ,DOUBLE *restrict u0     ,DOUBLE *restrict gradU0
             ,DOUBLE *restrict rCell  ,Temporal const ddt
             ,INT const nEq           ,INT const nEqNov
             ,INT const nAd           ,INT const nAdR
             ,short const maxNo       ,short const maxViz
             ,short const ndm         ,INT const numel
             ,short const ndf         ,short const storage
             ,bool const forces       ,bool const matrix 
             ,bool const calRcell     ,bool const unsym); 
/*...................................................................*/
 
/* ... montagem do sistemas de equacoes (Transporte)*/
  void systFormTrans(Loads *loads      
             ,Advection advT           ,Diffusion diffT
             ,INT    *restrict el      ,INT    *restrict nelcon 
             ,short  *restrict nen     ,short  *restrict nFace
             ,short  *restrict geomType,DOUBLE *restrict prop 
             ,short  *restrict calType ,short  *restrict mat  
             ,DOUBLE *restrict gCc
             ,DOUBLE *restrict gKsi   ,DOUBLE *restrict gmKsi
             ,DOUBLE *restrict gEta   ,DOUBLE *restrict gfArea
             ,DOUBLE *restrict gNormal,DOUBLE *restrict gVolume
             ,DOUBLE *restrict gXm    ,DOUBLE *restrict gXmcc   
             ,DOUBLE *restrict gvSkew ,DOUBLE *restrict gmvSkew   
             ,DOUBLE *restrict gDcca  ,DOUBLE *restrict density
             ,INT    *restrict ia     ,INT    *restrict ja
             ,DOUBLE *restrict al     ,DOUBLE *restrict ad 
             ,DOUBLE *restrict b      ,INT    *restrict id 
             ,short  *restrict faceR  ,short  *restrict faceL
             ,DOUBLE *restrict u0     ,DOUBLE *restrict gradU0 
             ,DOUBLE *restrict vel 
             ,DOUBLE *restrict rCell  ,Temporal const ddt
             ,INT const nEq           ,INT const nEqNov
             ,INT const nAd           ,INT const nAdR          
             ,short const maxNo       ,short const maxViz
             ,short const ndm         ,INT const numel
             ,short const ndf         ,short const storage
             ,bool const forces       ,bool const matrix 
             ,bool const calRcell     ,bool const unsym); 
/*...................................................................*/

/* ... montagem do sistemas de equacoes (Simple - VEl)*/
  void systFormSimpleVel(Loads *loadsVel,Loads *loadsPres 
             ,Advection advVel          ,Diffusion diffVel      
             ,short const typeSimple 
             ,INT    *restrict el       ,INT    *restrict nelcon 
             ,short  *restrict nen      ,short  *restrict nFace
             ,short  *restrict geomType ,DOUBLE *restrict prop 
             ,short  *restrict calType  ,short  *restrict mat  
             ,DOUBLE *restrict gCc
             ,DOUBLE *restrict gKsi     ,DOUBLE *restrict gmKsi
             ,DOUBLE *restrict gEta     ,DOUBLE *restrict gfArea
             ,DOUBLE *restrict gNormal  ,DOUBLE *restrict gVolume
             ,DOUBLE *restrict gXm      ,DOUBLE *restrict gXmcc
             ,DOUBLE *restrict gvSkew   ,DOUBLE *restrict gmvSkew 
             ,DOUBLE *restrict gDcca    ,DOUBLE *restrict density
             ,INT    *restrict ia       ,INT    *restrict ja
             ,DOUBLE *restrict al       ,DOUBLE *restrict ad
             ,DOUBLE *restrict b        ,INT    *restrict id 
             ,short  *restrict faceVelR ,short  *restrict faceVelL 
             ,short  *restrict facePresR,short  *restrict facePresL 
             ,DOUBLE *restrict pres     ,DOUBLE *restrict gradPres
             ,DOUBLE *restrict vel      ,DOUBLE *restrict gradVel
             ,DOUBLE *restrict dField   ,DOUBLE const underU
             ,DOUBLE *restrict rCell    ,Temporal const ddt
             ,INT const nEq             ,INT const nEqNov
             ,INT const nAd             ,INT const nAdR 
             ,short const maxNo         ,short const maxViz
             ,short const ndm           ,INT const numel
             ,short const ndf           ,short const storage
             ,bool const forces         ,bool const matrix 
             ,bool const calRcell       ,bool const unsym
             ,const bool sPressure);  
/*...................................................................*/

/* ... montagem do sistemas de equacoes (Simple - Pres)*/
  void systFormSimplePres(Loads *loadsVel,Loads *loadsPres 
						 ,Diffusion diffPres
             ,INT    *restrict el      ,INT    *restrict nelcon 
             ,short  *restrict nen     ,short  *restrict nFace
             ,short  *restrict geomType,DOUBLE *restrict prop 
             ,short  *restrict calType ,short  *restrict mat  
             ,DOUBLE *restrict gKsi   ,DOUBLE *restrict gmKsi
             ,DOUBLE *restrict gEta   ,DOUBLE *restrict gfArea
             ,DOUBLE *restrict gNormal,DOUBLE *restrict gVolume
             ,DOUBLE *restrict gXm    ,DOUBLE *restrict gXmcc   
             ,DOUBLE *restrict gvSkew ,DOUBLE *restrict gmvSkew   
             ,DOUBLE *restrict gDcca  ,DOUBLE *restrict density
             ,INT    *restrict ia     ,INT    *restrict ja
             ,DOUBLE *restrict al     ,DOUBLE *restrict ad
             ,DOUBLE *restrict b      ,INT    *restrict id
             ,short  *restrict faceVelR ,short  *restrict faceVelL
             ,short  *restrict facePresR,short  *restrict facePresL
             ,DOUBLE *restrict pres   ,DOUBLE *restrict gradPres
             ,DOUBLE *restrict vel    ,DOUBLE *restrict dField 
             ,DOUBLE *restrict rCell  ,Temporal const ddt
             ,INT const nEq           ,INT const nEqNov
             ,INT const nAd           ,INT const nAdR
             ,short const maxNo       ,short const maxViz
             ,short const ndm         ,INT const numel
             ,short const ndf         ,short const storage
             ,bool const forces       ,bool const matrix 
             ,bool const calRcell     ,bool const unsym); 
/*...................................................................*/

/*... correcao nao ortogonal para correcao de pressao*/
  void simpleNonOrthPres(Diffusion diffPres
               ,INT *restrict el           ,INT *restrict nelcon
							 ,short  *restrict nen       ,short  *restrict nFace
               ,short  *restrict geomType  ,DOUBLE *restrict prop 
               ,short  *restrict calType   ,short  *restrict mat
               ,DOUBLE *restrict cc
               ,DOUBLE *restrict gKsi      ,DOUBLE *restrict gmKsi 
               ,DOUBLE *restrict gEta      ,DOUBLE *restrict gfArea 
               ,DOUBLE *restrict gNormal   ,DOUBLE *restrict gVolume
               ,DOUBLE *restrict gXm       ,DOUBLE *restrict gXmcc 
               ,DOUBLE *restrict gvSkew    ,DOUBLE *restrict gmvSkew 
               ,DOUBLE *restrict gDcca     ,DOUBLE *restrict density
               ,DOUBLE *restrict b         ,INT    *restrict id
               ,short  *restrict facePresR ,DOUBLE *restrict pres 
               ,DOUBLE *restrict gradPres  ,DOUBLE *restrict dField 
               ,short const maxNo          ,short const maxViz
               ,short const ndm            ,INT const numel);
/*...................................................................*/

/*... chamada da biblioteca de elementos (difusao)*/
  void cellLibDif(Loads *loads
                 ,short *restrict lGeomType,DOUBLE *restrict lprop
                 ,INT   *restrict lViz     ,INT    *restrict lId
                 ,DOUBLE *restrict ksi     ,DOUBLE *restrict mksi    
                 ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea 
                 ,DOUBLE *restrict normal  ,DOUBLE *restrict volume 
                 ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
                 ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity 
                 ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
                 ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
                 ,DOUBLE *restrict lRcell  ,Temporal const ddt
                 ,short  *restrict lFaceR  ,short  *restrict lFaceL
                 ,DOUBLE *restrict u0      ,DOUBLE *restrict lGradU0 
                 ,short const nEn          ,short const nFace
                 ,short const ndm          ,short const lib 
                 ,INT const nel);
/*...................................................................*/

/*... chamada da biblioteca de elementos (transporte)*/
  void cellLibTrans(Loads *loads           
                 ,Advection advT           ,Diffusion diffT
                 ,short *restrict lGeomType,DOUBLE *restrict lprop
                 ,INT   *restrict lViz     ,INT    *restrict lId
                 ,DOUBLE *restrict ksi     ,DOUBLE *restrict mksi
                 ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea 
                 ,DOUBLE *restrict normal  ,DOUBLE *restrict volume 
                 ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
                 ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity 
                 ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
                 ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
                 ,DOUBLE *restrict lRcell  ,Temporal const ddt
                 ,short  *restrict lFaceR  ,short  *restrict lFaceL
                 ,DOUBLE *restrict lu0     ,DOUBLE *restrict lGradU0 
                 ,DOUBLE *restrict lVel    ,DOUBLE *restrict cc
                 ,short const nEn          ,short const nFace
                 ,short const ndm          ,short const lib 
                 ,INT const nel);
/*...................................................................*/


/*.......................... SIMPLE .................................*/
/*... chamada da biblioteca de elementos (escoamento-vel)*/
  void cellLibSimpleVel(Loads *loadsVel     ,Loads *loadsPres 
          ,Advection advVel           ,Diffusion diffVel
          ,short const typeSimple 
          ,short *restrict lGeomType  ,DOUBLE *restrict lprop
          ,INT   *restrict lViz       ,INT *restrict lId  
          ,DOUBLE *restrict ksi       ,DOUBLE *restrict mKsi
          ,DOUBLE *restrict eta       ,DOUBLE *restrict fArea
          ,DOUBLE *restrict normal    ,DOUBLE *restrict volume
          ,DOUBLE *restrict xm        ,DOUBLE *restrict xmcc
          ,DOUBLE *restrict dcca      ,DOUBLE *restrict lDensity
          ,DOUBLE *restrict vSkew     ,DOUBLE *restrict mvSkew
          ,DOUBLE *restrict lA        ,DOUBLE *restrict lB
          ,DOUBLE *restrict lRcell    ,Temporal const ddt
          ,short  *restrict lFaceVelR ,short  *restrict lFaceVelL
          ,short  *restrict lFacePresR,short  *restrict lFacePresL
          ,DOUBLE *restrict pres      ,DOUBLE *restrict gradPres 
          ,DOUBLE *restrict vel       ,DOUBLE *restrict gradVel 
          ,DOUBLE *restrict dField    ,DOUBLE *restrict cc 
          ,DOUBLE const underU        ,const bool sPressure
          ,short const nEn            ,short  const nFace     
          ,short const ndm            ,short const lib    
          ,INT const nel);

/*... chamada da biblioteca de elementos (escoamento-pres)*/
  void cellLibSimplePres(Loads *loadsVel     ,Loads *loadsPres
	           	 ,Diffusion diffPres
               ,short *restrict lGeomType,DOUBLE *restrict lprop
               ,INT   *restrict lViz     ,INT *restrict lId  
               ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
               ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea
               ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
               ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
               ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
               ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
               ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
               ,DOUBLE *restrict lRcell  ,Temporal const ddt
               ,short  *restrict lFaceVelR ,short  *restrict lFaceVelL
               ,short  *restrict lFacePresR,short  *restrict lFacePresL
               ,DOUBLE *restrict pres    ,DOUBLE *restrict gradPres 
               ,DOUBLE *restrict vel     ,DOUBLE *restrict dField 
               ,short const nEn          ,short  const nFace     
               ,short const ndm          ,short const lib    
               ,INT const nel);
/*...................................................................*/

/*... chamada da biblioteca de correcao nao ortogonal para correcao 
      de pressao*/
  void cellLibSimpleNonOrthPres(Diffusion diffPres
               ,short *restrict lGeomType
               ,DOUBLE *restrict lprop   ,INT   *restrict lViz
               ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
               ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea
               ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
               ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
               ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
               ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
               ,DOUBLE *restrict lB      
               ,DOUBLE *restrict pres    ,DOUBLE *restrict gradPres 
               ,DOUBLE *restrict dField  ,DOUBLE *restrict cc
               ,short const nEn          ,short  const nFace     
               ,short const ndm          ,INT const nel);
/*...................................................................*/

/*...................................................................*/

/*... carga por elmento e condicoes pescritas por celula*/
  void cellPload(Loads *loads           ,DOUBLE *restrict gCc
                ,short  *restrict faceR ,short *restrict faceS
                ,DOUBLE *restrict volume,INT *restrict id 
                ,DOUBLE *restrict u     ,DOUBLE *restrict f
                ,INT const numel        ,short const ndf
                ,short const ndm        ,short const maxViz);

/*... carga por elmento e condicoes pescritas por celula no 
      metodo simple*/
  void pLoadSimple(DOUBLE *restrict sP  ,DOUBLE *restrict p
          ,DOUBLE *restrict tA          ,DOUBLE *restrict velC
          ,DOUBLE *restrict n
          ,DOUBLE *restrict gradVel     ,DOUBLE *restrict xmcc
          ,DOUBLE const viscosityC      ,DOUBLE const densityC
          ,DOUBLE const fArea           ,DOUBLE const dcca
          ,Loads ld                     ,short const ndm
          ,bool const fCal1             ,bool const fCal2);

  void pLoadSimplePres(DOUBLE *restrict sP  ,DOUBLE *restrict p
          ,DOUBLE *restrict tA
          ,DOUBLE const viscosityC,DOUBLE const densityC
          ,DOUBLE const wfn
          ,DOUBLE const fArea     ,DOUBLE const dcca 
          ,Loads ld               ,bool const fCal);
/*...................................................................*/

/*... carga por elmento e condicoes pescritas por celula no metodo 
      simple*/
  void cellPloadSimple(Loads *loads           ,DOUBLE *restrict cc 
                    ,short  *restrict faceR ,short *restrict faceL
                    ,DOUBLE *restrict volume
                    ,INT *restrict idVel    ,INT *restrict idPres
                    ,DOUBLE *restrict vel   ,DOUBLE *restrict pres
                    ,DOUBLE *restrict fVel  ,DOUBLE *restrict fPres
                    ,INT const numel        ,short const ndf
                    ,short const ndm        ,short const maxViz);
/*...................................................................*/

/*... */
  void updateCellValue(DOUBLE *restrict u,DOUBLE *restrict x
                      ,INT *restrict id  ,Interface *iNeq
                      ,INT const numel   ,short const ndf
                      ,bool const fAdd   ,short const fCom);
  
  void updateCellValueSimple(DOUBLE *restrict u,DOUBLE *restrict x
                      ,INT *restrict id  ,Interface *iNeq
                      ,INT const numel   ,INT const nEq 
                      ,short const ndf
                      ,bool const fAdd   ,short const fCom);
/*...................................................................*/

/*... funcoes geometricas*/
  void cellGeom2D(DOUBLE *restrict lx       ,short *restrict lnFace
                 ,short  *restrict lGeomType,DOUBLE *restrict xc
                 ,DOUBLE *restrict ksi      ,DOUBLE *restrict mksi
                 ,DOUBLE *restrict eta      ,DOUBLE *restrict mEta
                 ,DOUBLE *restrict normal   ,DOUBLE *restrict volume
                 ,DOUBLE *restrict xm       ,DOUBLE *restrict xmcc
                 ,DOUBLE *restrict dcca
                 ,DOUBLE *restrict vSkew    ,DOUBLE *restrict mvSkew
                 ,short  *restrict sn       ,short const maxNo   
                 ,short const maxViz        ,short const ndm
                 ,INT const nel);

   void cellGeom3D(DOUBLE *restrict lx      ,short *restrict lGeomType
                 ,short *restrict lnFac     ,short *restrict lnEn
                 ,DOUBLE *restrict xc
                 ,DOUBLE *restrict ksi      ,DOUBLE *restrict mksi
                 ,DOUBLE *restrict eta      ,DOUBLE *restrict fArea
                 ,DOUBLE *restrict normal   ,DOUBLE *restrict volume
                 ,DOUBLE *restrict xm       ,DOUBLE *restrict xmcc
                 ,DOUBLE *restrict dcca     
                 ,DOUBLE *restrict vSkew    ,DOUBLE *restrict mvSkew
                 ,short  *restrict sn
                 ,short const maxNo         ,short const maxViz
                 ,short const ndm           ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (difusao)*/
  void cellDif2D(Loads *loads
                ,short *restrict lGeomType,DOUBLE *restrict lprop
                ,INT   *restrict lViz     ,INT *restrict lId
                ,DOUBLE *restrict ksi     ,DOUBLE *restrict mksi
                ,DOUBLE *restrict eta     ,DOUBLE *restrict mEta
                ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
                ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
                ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
                ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
                ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
                ,DOUBLE *restrict lRcell  ,Temporal const ddt
                ,short  *restrict lFaceR  ,short *restrict lFaceL
                ,DOUBLE *restrict u0      ,DOUBLE *restrict lGradU0
                ,short const nen          ,short const nFace
                ,short const ndm          ,INT const nel);
  
  void cellDif3D(Loads *loads
                ,short *restrict lGeomType,DOUBLE *restrict prop
                ,INT *restrict lViz       ,INT *restrict lId  
                ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
                ,DOUBLE *restrict eta     ,DOUBLE *restrict fArea
                ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
                ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
                ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity 
                ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
                ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
                ,DOUBLE *restrict lRcell  ,Temporal const ddt   
                ,short  *restrict lFaceR  ,short  *restrict lFaceL
                ,DOUBLE *restrict u0      ,DOUBLE *restrict gradU0
                ,const short nEn          ,short const nFace    
                ,const short ndm          ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (transporte)*/
  void cellTrans2D(Loads *loads           
                ,Advection advT           ,Diffusion diffT
                ,short *restrict lGeomType,DOUBLE *restrict lprop
                ,INT   *restrict lViz     ,INT *restrict lId
                ,DOUBLE *restrict ksi     ,DOUBLE *restrict mksi
                ,DOUBLE *restrict eta     ,DOUBLE *restrict mEta
                ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
                ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
                ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity   
                ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
                ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
                ,DOUBLE *restrict lRcell  ,Temporal const ddt
                ,short  *restrict lFaceR  ,short *restrict lFaceL
                ,DOUBLE *restrict u0      ,DOUBLE *restrict lGradU0
                ,DOUBLE *restrict vel     ,DOUBLE *restrict cc
                ,short const nen          ,short const nFace
                ,short const ndm          ,INT const nel);
  
  void cellTrans3D(Loads *loads           
                ,Advection advT           , Diffusion diffT
                ,short *restrict lGeomType,DOUBLE *restrict lprop
                ,INT   *restrict lViz     ,INT *restrict lId
                ,DOUBLE *restrict ksi     ,DOUBLE *restrict mksi
                ,DOUBLE *restrict eta     ,DOUBLE *restrict mEta
                ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
                ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
                ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
                ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
                ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
                ,DOUBLE *restrict lRcell  ,Temporal const ddt
                ,short  *restrict lFaceR  ,short *restrict lFaceL
                ,DOUBLE *restrict u0      ,DOUBLE *restrict lGradU0
                ,DOUBLE *restrict vel     ,DOUBLE *restrict cc
                ,short const nen          ,short const nFace    
                ,short const ndm          ,INT const nel);
/*...................................................................*/

/*.......................... SIMPLE .................................*/

/*... biblioteca de celulas (simple - vel)*/
  void cellSimpleVel2D(Loads *loadsVel    ,Loads *loadsPres 
              ,Advection  advVel          ,Diffusion diffVel
              ,short const typeSimple 
              ,short *restrict lGeomType  ,DOUBLE *restrict lprop
              ,INT   *restrict lViz       ,INT *restrict lId
              ,DOUBLE *restrict ksi       ,DOUBLE *restrict mksi
              ,DOUBLE *restrict eta       ,DOUBLE *restrict mEta
              ,DOUBLE *restrict normal    ,DOUBLE *restrict volume
              ,DOUBLE *restrict xm        ,DOUBLE *restrict xmcc
              ,DOUBLE *restrict dcca      ,DOUBLE *restrict lDensity
              ,DOUBLE *restrict vSkew     ,DOUBLE *restrict mvSkew
              ,DOUBLE *restrict lA        ,DOUBLE *restrict lB
              ,DOUBLE *restrict lRcell    ,Temporal const ddt
              ,short  *restrict lFaceVelR ,short  *restrict lFaceVelL
              ,short  *restrict lFacePresR,short  *restrict lFacePresL
              ,DOUBLE *restrict pres      ,DOUBLE *restrict gradPres 
              ,DOUBLE *restrict vel       ,DOUBLE *restrict gradVel
              ,DOUBLE *restrict dField    ,DOUBLE *restrict cc
              ,DOUBLE const underU        ,bool const sPressure
              ,short const nen            ,short const nFace
              ,short const ndm            ,INT const nel);

/*... biblioteca de celulas (simple - vel)*/
  void cellSimpleVel3D(Loads *loadsVel    ,Loads *loadsPres
              ,Advection  advVel          ,Diffusion  diffVel
							,short const typeSimple 
              ,short *restrict lGeomType  ,DOUBLE *restrict lprop
              ,INT   *restrict lViz       ,INT *restrict lId
              ,DOUBLE *restrict ksi       ,DOUBLE *restrict mksi
              ,DOUBLE *restrict eta       ,DOUBLE *restrict mEta
              ,DOUBLE *restrict normal    ,DOUBLE *restrict volume
              ,DOUBLE *restrict xm        ,DOUBLE *restrict xmcc
              ,DOUBLE *restrict dcca      ,DOUBLE *restrict lDensity
              ,DOUBLE *restrict vSkew     ,DOUBLE *restrict mvSkew
              ,DOUBLE *restrict lA        ,DOUBLE *restrict lB
              ,DOUBLE *restrict lRcell    ,Temporal const ddt
              ,short  *restrict lFaceVelR ,short  *restrict lFaceVelL
              ,short  *restrict lFacePresR,short  *restrict lFacePresL
              ,DOUBLE *restrict pres      ,DOUBLE *restrict gradPres 
              ,DOUBLE *restrict vel       ,DOUBLE *restrict gradVel 
              ,DOUBLE *restrict dField    ,DOUBLE *restrict cc 
              ,DOUBLE const underU        ,bool const sPressure
              ,short const nen            ,short const nFace
              ,short const ndm            ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (simple - pres)*/
  void cellSimplePres2D(Loads *loadsVel    ,Loads *loadsPres   
							,Diffusion diffPres        
              ,short *restrict lGeomType,DOUBLE *restrict prop
              ,INT *restrict lViz       ,INT *restrict lId  
              ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
              ,DOUBLE *restrict eta     ,DOUBLE *restrict mEta
              ,DOUBLE *restrict normal  ,DOUBLE *restrict area   
              ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
              ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
              ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
              ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
              ,DOUBLE *restrict lRcell  
              ,short  *restrict lFaceVelR ,short  *restrict lFaceVelL
              ,short  *restrict lFacePresR,short  *restrict lFacePresL
              ,DOUBLE *restrict pres    ,DOUBLE *restrict gradPres 
              ,DOUBLE *restrict vel     ,DOUBLE *restrict dField  
              ,const short nEn          ,short const nFace    
              ,const short ndm          ,INT const nel);
  
/*... biblioteca de celulas (simple - pres)*/
  void cellSimplePres3D(Loads *loadsVel    ,Loads *loadsPres
							,Diffusion diffPres
		          ,short *restrict lGeomType,DOUBLE *restrict prop
              ,INT *restrict lViz       ,INT *restrict lId  
              ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
              ,DOUBLE *restrict eta     ,DOUBLE *restrict mEta
              ,DOUBLE *restrict normal  ,DOUBLE *restrict area   
              ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
              ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
              ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
              ,DOUBLE *restrict lA      ,DOUBLE *restrict lB
              ,DOUBLE *restrict lRcell  
              ,short  *restrict lFaceVelR ,short  *restrict lFaceVelL
              ,short  *restrict lFacePresR,short  *restrict lFacePresL
              ,DOUBLE *restrict pres    ,DOUBLE *restrict gradPres 
              ,DOUBLE *restrict vel     ,DOUBLE *restrict dField  
              ,const short nEn          ,short const nFace    
              ,const short ndm          ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (simple - pres-correcao non-ortogonal)*/
  void cellSimpleNonOrthPres2D(Diffusion diffPres
              ,short *restrict lGeomType
              ,DOUBLE *restrict prop    ,INT *restrict lViz
              ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
              ,DOUBLE *restrict eta     ,DOUBLE *restrict mEta
              ,DOUBLE *restrict normal  ,DOUBLE *restrict area   
              ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
              ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
              ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
              ,DOUBLE *restrict lB      
              ,DOUBLE *restrict pres    ,DOUBLE *restrict gradPres 
              ,DOUBLE *restrict dField  ,DOUBLE *restrict cc
              ,const short nEn          ,short const nFace
              ,const short ndm          ,INT const nel);
/*...................................................................*/

/*... biblioteca de celulas (simple - pres-correcao non-ortogonal)*/
  void cellSimpleNonOrthPres3D(Diffusion diffPres
              ,short *restrict lGeomType
              ,DOUBLE *restrict prop    ,INT *restrict lViz
              ,DOUBLE *restrict ksi     ,DOUBLE *restrict mKsi
              ,DOUBLE *restrict eta     ,DOUBLE *restrict mEta
              ,DOUBLE *restrict normal  ,DOUBLE *restrict area   
              ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc
              ,DOUBLE *restrict dcca    ,DOUBLE *restrict lDensity
              ,DOUBLE *restrict vSkew   ,DOUBLE *restrict mvSkew
              ,DOUBLE *restrict lB      
              ,DOUBLE *restrict pres    ,DOUBLE *restrict gradPres 
              ,DOUBLE *restrict dField  ,DOUBLE *restrict cc
              ,const short nEn          ,short const nFace
              ,const short ndm          ,INT const nel);
/*...................................................................*/

/*...................................................................*/

/*...*/
  void cellTransient(DOUBLE *restrict volume  ,INT *restrict id 
                   ,DOUBLE *restrict u0      ,DOUBLE *restrict u
                   ,DOUBLE *restrict density ,DOUBLE *restrict f
                   ,DOUBLE const dt
                   ,INT const numel          ,short const ndf
                   ,short const type         ,bool const fAdd);
/*...................................................................*/

/*...*/
  void cellTransientSimple(DOUBLE *restrict volume  ,INT *restrict id 
                        ,DOUBLE *restrict u0      ,DOUBLE *restrict u
                        ,DOUBLE *restrict density ,DOUBLE *restrict f
                        ,DOUBLE const dt          ,INT const nEq
                        ,INT const numel          ,short const ndf
                        ,short const type         ,bool const fAdd);
/*...................................................................*/

/*... biblioteca de reconstrucao de gradiente*/
  void cellLibRcGrad(Loads *loads
           ,INT   *restrict lViz,DOUBLE *restrict lProp
           ,DOUBLE *restrict leastSquare ,DOUBLE *restrict leastSquareR
           ,DOUBLE *restrict ksi         ,DOUBLE *restrict mKsi
           ,DOUBLE *restrict eta         ,DOUBLE *restrict fArea
           ,DOUBLE *restrict normal      ,DOUBLE *restrict volume
           ,DOUBLE *restrict vSkew       
           ,DOUBLE *restrict xm          ,DOUBLE *restrict xmcc 
           ,DOUBLE *restrict lDcca  
           ,short  *restrict lFaceR      ,short *restrict lFaceL
           ,DOUBLE *restrict u           ,DOUBLE *restrict gradU 
           ,DOUBLE *restrict nU          ,short const ty 
           ,short const nFace            ,short const ndm  
           ,short const lib              ,short const ndf  
           ,short *restrict  isNodIN    ,INT const nel);

  void rcLeastSquare(DOUBLE *restrict gKsi    ,DOUBLE *restrict gmKsi
                  ,DOUBLE *restrict lSquare   ,DOUBLE *restrict lSquareR
                  ,short *restrict nFace       
                  ,INT const numel          ,short const maxViz
                  ,short const type         ,short const ndm);
/*...................................................................*/

/*... least square*/
  void leastSquareMatrix(DOUBLE *restrict lKsi    ,DOUBLE *restrict lmKsi
                    ,DOUBLE *restrict lLsquare,DOUBLE *restrict lSquareR
                    ,short const type         
                    ,short const lnFace       ,short const ndm);


  void leastSquare(Loads *loads
                 ,DOUBLE *restrict lLsquare,INT *restrict lViz 
                 ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc 
                 ,DOUBLE *restrict lProp   ,DOUBLE *restrict lDcca
                 ,DOUBLE *restrict u       ,DOUBLE *restrict gradU
                 ,short  *restrict lFaceR  ,short *restrict lFaceL
                 ,short const nFace        ,short const ndf
                 ,short const ndm          ,INT const nel);

  void leastSquareQR(Loads *loads
                 ,DOUBLE *restrict lLs     ,DOUBLE *restrict lLsR
                 ,DOUBLE *restrict lProp   ,DOUBLE *restrict lDcca
                 ,INT *restrict lViz       
                 ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc  
                 ,DOUBLE *restrict u       ,DOUBLE *restrict gradU
                 ,short  *restrict lFaceR  ,short *restrict lFaceL
                 ,short const nFace        ,short const ndf
                 ,short const ndm);

/*...................................................................*/

  
/*... green-gauss*/
  void greenGaussCell(Loads *loads
               ,INT *restrict lViz ,DOUBLE *restrict mKsi
               ,DOUBLE *restrict lProp   ,DOUBLE *restrict lDcca
               ,DOUBLE *restrict eta     ,DOUBLE *restrict mfArea
               ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
               ,DOUBLE *restrict mvSkew
               ,DOUBLE *restrict xm      ,DOUBLE *restrict xmcc    
               ,short  *restrict lFaceR  ,short  *restrict lFaceS
               ,DOUBLE *restrict u       ,DOUBLE *restrict gradU 
               ,short const nFace        ,short const ndm   
               ,short const ndf          ,INT const nel);
  void greenGaussNode(INT *restrict lViz   ,DOUBLE *restrict mEta
               ,DOUBLE *restrict normal  ,DOUBLE *restrict volume
               ,DOUBLE *restrict u       ,DOUBLE *restrict gradU 
               ,short *restrict isNod       
               ,short const nFace        ,short const ndm   
               ,short const ndf          ,short const ty);
/*...................................................................*/

/*....*/
  void meshQuality(MeshQuality *mq
                ,short  *restrict nFace   ,DOUBLE *restrict volume
                ,DOUBLE *restrict gKsi    ,DOUBLE *restrict gNormal
                ,DOUBLE *restrict gmvSkew
                ,short const maxViz      ,short const ndm
                ,INT const numel); 
/*...................................................................*/

/*...*/
  void wallFluid(short *restrict faceR ,INT *restrict nelcon
              ,short *restrict nFace     
              ,INT const nEl         ,short const maxViz);
/*...................................................................*/

/*... funcoes limitadoras de fluxo*/
  DOUBLE limitFaceBase(DOUBLE const rr,short const iCod);
  
  DOUBLE faceBaseTvd(short const nAresta    ,short const idCell
                 ,DOUBLE *restrict u0
                 ,DOUBLE *restrict gradUv,DOUBLE *restrict gradUp
                 ,DOUBLE *restrict lKsi  ,DOUBLE const lModKsi 
                 ,DOUBLE const cv
                 ,short const iCod       ,short const ndm);

  DOUBLE faceBaseTvdV1(DOUBLE const uC        ,DOUBLE const uV
                 ,DOUBLE *restrict gradUc,DOUBLE *restrict gradUv
                 ,DOUBLE *restrict lKsi  ,DOUBLE const lModKsi 
                 ,DOUBLE const cv
                 ,short const iCod       ,short const ndm);
/*...................................................................*/

/*... NVD-HR-CDC*/
  DOUBLE faceBaseNvd(DOUBLE const uP, DOUBLE const uV
                     , DOUBLE *restrict gradUp, DOUBLE *restrict gradUv
                     , DOUBLE *restrict lKsi, DOUBLE const lModKsi
                     , DOUBLE const m
                     , short const iCod, short const ndm);
  DOUBLE nvd(DOUBLE const phiTil, short const iCod);
/*...................................................................*/

/*... parametro fisicos do escoamento*/  
  void parameterCell(DOUBLE *restrict vel  ,DOUBLE *restrict prop
                ,DOUBLE *restrict density,DOUBLE *restrict volume
                ,short  *restrict mat 
                ,DOUBLE *cfl             ,DOUBLE *reynolds
                ,bool *fParameter        ,DOUBLE const dt
                ,INT const nEl           ,short const ndm);
/*...................................................................*/

/*... funcoes de correcao atrasada*/
  void difusionScheme(DOUBLE *restrict n ,DOUBLE *restrict ksi 
                    ,DOUBLE const fArea ,DOUBLE const lModKsi  
                    ,DOUBLE *restrict e ,DOUBLE *restrict t   
                    ,short const ndm    ,short const iCod);

  void difusionSchemeAnisotropic(DOUBLE *restrict s,DOUBLE *restrict ksi
                                ,DOUBLE *restrict e,DOUBLE *restrict t
                                ,short const ndm   ,short const iCod);
  void advectiveScheme(DOUBLE *restrict velC       ,DOUBLE *restrict velV
                     ,DOUBLE *restrict gradVelC   ,DOUBLE *restrict gradVelV
                     ,DOUBLE *restrict gradVelComp,DOUBLE *restrict vSkew
                     ,DOUBLE *restrict xmcc       ,DOUBLE const wfn
                     ,DOUBLE *restrict ksi        ,DOUBLE const modKsi
                     ,DOUBLE const m              ,DOUBLE *restrict cvc
                     ,short const ndm
                     ,short const iCod1, short const iCod2);
   
  void advectiveSchemeScalar(DOUBLE const uC     ,DOUBLE const uV
                     ,DOUBLE *restrict gradUc   ,DOUBLE *restrict gradUv
                     ,DOUBLE *restrict gradUcomp,DOUBLE *restrict vSkew
                     ,DOUBLE *restrict rC       ,DOUBLE *restrict rV
                     ,DOUBLE *restrict ksi      ,DOUBLE const modKsi
                     ,DOUBLE const m            ,DOUBLE *cvc
                     ,short const ndm
                     ,short const iCod1, short const iCod2);

  DOUBLE deferredCd(DOUBLE const velC,DOUBLE const velV
                   ,DOUBLE const wfn);

  DOUBLE upwindLinearV1(DOUBLE const uC  ,DOUBLE const uV
                 ,DOUBLE *restrict gradUc,DOUBLE *restrict gradUv
                 ,DOUBLE *restrict rC    ,DOUBLE *restrict rV
                 ,DOUBLE const m         ,short const ndm);
/*...................................................................*/

/*... funcoes de apoio*/ 
  void gradFaceNull(DOUBLE *restrict gradVelFace
                   ,DOUBLE *restrict gradVelCell
                   ,DOUBLE *restrict xmcc       ,short const ndm);

  short sn(short *s,short ty, INT nel);
  DOUBLE areaQuadCell(DOUBLE *restrict eta,short ndm);
  DOUBLE areaTriaCell(DOUBLE *restrict eta,short ndm);
  DOUBLE areaCell(DOUBLE *eta,short ty,short ndm,INT nel);
  DOUBLE volume3DGreenGauss(DOUBLE *restrict xm,DOUBLE *restrict normal
                         ,DOUBLE *restrict fArea
                         ,short const nFace);
  void vectorKm2d(DOUBLE *restrict x      ,DOUBLE *restrict xc
                 ,DOUBLE *restrict xm
                  ,DOUBLE *restrict vSkew ,DOUBLE *restrict mvSkew
                 ,short  *restrict sn     ,short const nFace
                 ,short const maxViz      ,short const maxNo       
                 ,short const ndm         ,INT const nel);

  void vectorKm3d(DOUBLE *restrict xc     ,DOUBLE *restrict xm
                 ,DOUBLE *restrict ksi    ,DOUBLE *restrict normal
                 ,DOUBLE *restrict vSkew ,DOUBLE *restrict mvSkew
                 ,short const nFace      ,short const ndm
                 ,short const maxViz     ,INT nel);
  
  void pLoad(DOUBLE *restrict sP  ,DOUBLE *restrict p
            ,DOUBLE *restrict tA
            ,DOUBLE const coefDifC,DOUBLE const densityC
            ,DOUBLE const wfn     ,DOUBLE *restrict xm                   
            ,DOUBLE const fArea   ,DOUBLE const dcca
            ,Loads ld             ,bool const fCal);
  void setNvd(char *word,short *iCod);
  void setTvd(char *word,short *iCod);
  void  setDiffusionScheme(char *word,short *iCod);
  void  setAdvectionScheme(char *word, Advection *adv,FILE *fileIn);
  DOUBLE sizeCar(DOUBLE const volume,short const ndm);
/*...................................................................*/

#endif/*_CELLLOOP_H_*/
