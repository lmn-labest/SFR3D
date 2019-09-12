#include<CellLib.h>
/*********************************************************************
 * Data de criacao    : 10/09/2017                                   *
 * Data de modificaco : 12/12/2017                                   *
 *-------------------------------------------------------------------*
 * CELLLIBSIMPLETURBULENCE: chamada de bibliotecas de celulas para   *
 * o calculo da viscosidae turbelenta                                *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * tModel    -> modelo de turbulencia                                *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * fArea     -> area das faces                                       *
 * normal    -> vetores normais as faces das celulas                 *
 * volume    -> volume celula central                                *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os *
 *            centrois compartilhado nessa face da celula central    *
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * vel       -> campo de velocidade conhecido                        *
 * lDensity  -> massa especifica com variacao temporal               *
 * lViscosity-> viscosidade dinamica com variacao temporal           *
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 * dynamic   -> M:L e M:M do les dinamico                            *
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * lib       -> numero da biblioteca                                 *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lViscosity-> viscosidae turbulenta                                *
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellLibTurbulence(Loads *lVel    , Turbulence tModel 
      , short *RESTRICT lGeomType     , DOUBLE *RESTRICT lprop
      , INT   *RESTRICT lViz          , DOUBLE *RESTRICT ksi
      , DOUBLE *RESTRICT mKsi         , DOUBLE *RESTRICT eta
      , DOUBLE *RESTRICT fArea        , DOUBLE *RESTRICT normal    
      , DOUBLE *RESTRICT volume       , DOUBLE *RESTRICT xm        
      , DOUBLE *RESTRICT xmcc         , DOUBLE *RESTRICT dcca    
      , DOUBLE *RESTRICT cc           , DOUBLE *RESTRICT vSkew   
      , DOUBLE *RESTRICT mvSkew       , short  *RESTRICT faceVelR 
      , short *RESTRICT faceVelL      , DOUBLE *RESTRICT vel    
      , DOUBLE *RESTRICT gradVel      , DOUBLE *RESTRICT lDensity  
      , DOUBLE const dViscosity       , DOUBLE *viscosity          
      , DOUBLE *RESTRICT wallPar      , DOUBLE const cDyn
      , short const nEn               , short  const nFace
      , short const ndm               , short const lib
      , INT const nel)
{

/*... */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2){
      cellLes(lVel       , tModel       
            , lGeomType  , lprop
            , lViz       
            , ksi        , mKsi
            , eta        , fArea  
            , normal     , volume
            , xm         , xmcc
            , dcca       , cc
            , vSkew      , mvSkew 
            , faceVelR   , faceVelL      
            , vel        , gradVel
            , lDensity   , dViscosity
            , viscosity  
            , nEn        , nFace 
            , ndm        , nel);  
    }
/*..................................................................*/

/*... 3D*/
    else if(ndm == 3){
      eddyViscosity3D(lVel          , tModel       
              , lGeomType     , lprop 
              , lViz          , fArea  
              , normal        , volume
              , dcca          , faceVelR  
              , faceVelL      , vel  
              , gradVel       , lDensity 
              , dViscosity    , viscosity 
              , wallPar       , cDyn
              , nEn           , nFace 
              , ndm           , nel); 
    }  
/*..................................................................*/
  }

}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 22/08/2017                                   *
* Data de modificaco : 05/06/2019                                   *
*-------------------------------------------------------------------*
* CELLLIBENERGY: chamada de bibliotecas de celulas para             *
* problema de escoamento de fluidos (Energy)                        *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* loads     -> definicoes de cargas                                 *
* loads     -> definicoes de cargas para velocidade                 *
* adv       -> tecnica da discretizacao do termo advecao            *
* diff      -> tecnica da discretizacao do termo difusivo           *
* tModel    -> modelo de turbulencia                                *
* eModel    -> modelo da equacao de energia                         *
* cModel    -> modelo de combustao                                  *
* vProp     -> propedades variaveis                                 *
* lGeomType -> tipo geometrico da celula central e seus vizinhos    *
* lViz      -> viznhos da celula central                            *
* lId       -> equa da celula                                       *
* Ksi       -> vetores que unem centroide da celula central aos     *
*            vizinhos destas                                        *
* mKsi      -> modulo do vetor ksi                                  *
* eta       -> vetores paralelos as faces das celulas               *
* fArea     -> area das faces                                       *
* normal    -> vetores normais as faces das celulas                 *
* volume    -> volume celula central                                *
* xm        -> pontos medios das faces da celula central            *
* xmcc      -> vetores que unem o centroide aos pontos medios das   *
*            faces da celula central                                *
* vSkew     -> vetor entre o ponto medio a intersecao que une os    *
*            centrois compartilhado nessa face da celula central    *
* mvSkew    -> distacia entre o ponto medio a intersecao que une os *
*            centrois compartilhado nessa face da celula central    *
* dcca      -> menor distacia do centroide central a faces desta    *
*              celula                                               *
* cc        -> centroides da celula centra e seus vizinhos          *
* lA        -> nao definido                                         *
* lB        -> nao definido                                         *
* lRcell    -> nao definido                                         *
* ddt       -> discretizacao temporal                               *
* faceVelR  -> restricoes por elemento de velocidades               *
* faceVelL  -> carga por elemento de velocidades                    *
* facePresR -> restricoes por elemento de pressao                   *
* facePresL -> carga por elemento de pressao                        *
* vel       -> campo de velocidade conhecido                        *
* dField    -> matriz D do metodo simple                            *
* lDensity  -> massa especifica com variacao temporal               *
* lSheat    -> calor especifico com variacao temporal               *
* lDviscosity-> viscosidade dinamica com variacao temporal          *
* ltConductvity -> condutividade termica com variacao temporal      *
* lEnthalpyk -> Entalpia das especies                               *
* lGradY    -> gradiente das especies                               *
* yFac      -> fracoes massicas                                     *
* rateHeat  -> taxa de liberacao de energia por reacao quimica      *
* dField    -> matriz D do metodo simple                            * 
* wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
* underU    -> fator underrelaxtion sinple                          *
* nEn       -> numero de nos da celula central                      *
* nFace     -> numero de faces da celula central                    *
* ndm       -> numero de dimensoes                                  *
* lib       -> numero da biblioteca                                 *
* nel       -> numero da celula                                     *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* lA        -> coeficiente da linha i                               *
* lB        -> vetor de forca da linha i                            *
*-------------------------------------------------------------------*
*********************************************************************/
void cellLibEnergy(Loads *lEnergy  , Loads *lVel                      
     , Advection  *adv             , Diffusion *diff      
     , Turbulence *tModel          , EnergyModel *model  
     , Combustion *cModel          , PropVarFluid *vProp         
     , short *RESTRICT lGeomType   
     , INT   *RESTRICT lViz        , INT *RESTRICT lId
     , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
     , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT fArea
     , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume
     , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
     , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
     , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
     , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
     , DOUBLE *RESTRICT lRcell     , Temporal *ddt
     , short  *RESTRICT lFaceR     , short  *RESTRICT lFaceL
     , short  *RESTRICT lFaceVelR  , short  *RESTRICT lFaceVelL     
     , DOUBLE *RESTRICT u          , DOUBLE *RESTRICT gradU
     , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
     , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres  
     , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lSheat
     , DOUBLE *RESTRICT lDviscosity, DOUBLE *RESTRICT lTconductvity
     , DOUBLE *RESTRICT lEnthalpyk , DOUBLE *RESTRICT lGradY 
     , DOUBLE *RESTRICT lDiffY     , DOUBLE *RESTRICT lYfrac 
     , DOUBLE const lRateHeat  
     , DOUBLE *RESTRICT dField     , DOUBLE *RESTRICT wallPar
     , DOUBLE const underU           
     , short const nEn             , short  const nFace
     , short const ndm             , short const lib
     , INT const nel)
{

/*... */
  if (lib == 1) {
/*... 2D*/
    if (ndm == 2)  
      cellEnergy2D(lEnergy    , lVel
                 , adv        , diff
                 , tModel     , model 
                 , vProp 
                 , lGeomType  
                 , lViz       , lId
                 , ksi        , mKsi
                 , eta        , fArea
                 , normal     , volume
                 , xm         , xmcc
                 , dcca       , cc
                 , vSkew      , mvSkew 
                 , lA         , lB
                 , lRcell     , *ddt
                 , lFaceR     , lFaceL
                 , lFaceVelR  , lFaceVelL
                 , u          , gradU    
                 , vel        , gradVel
                 , pres       , gradPres  
                 , lDensity   , lSheat
                 , lDviscosity, lTconductvity
                 , dField     , wallPar
                 , underU       
                 , nEn        , nFace
                 , ndm        , nel);    
/*..................................................................*/

/*... 3D*/
      else if (ndm == 3) {
        cellEnergy3D(lEnergy  , lVel
                 , adv        , diff
                 , tModel     , model 
                 , cModel     , vProp 
                 , lGeomType  
                 , lViz       , lId
                 , ksi        , mKsi
                 , eta        , fArea
                 , normal     , volume
                 , xm         , xmcc
                 , dcca       , cc
                 , vSkew      , mvSkew 
                 , lA         , lB
                 , lRcell     , ddt
                 , lFaceR     , lFaceL
                 , lFaceVelR  , lFaceVelL
                 , u          , gradU    
                 , vel        , gradVel
                 , pres       , gradPres  
                 , lDensity   , lSheat
                 , lDviscosity, lTconductvity
                 , lEnthalpyk , lGradY           
                 , lDiffY     , lYfrac
                 , lRateHeat
                 , dField     , wallPar
                 , underU       
                 , nEn        , nFace
                 , ndm        , nel);    
/*..................................................................*/
    }  
/*..................................................................*/
  }

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 18/01/2018                                   *
 * Data de modificaco : 22/01/2018                                   *
 *-------------------------------------------------------------------*
 * cellLibOneEqK: chamada de bibliotecas de celulas para             *
 * problema de escoamento de fluidos (Energia cinetica turbulenta)   * 
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loads   -> definicoes de cargas K                                 *
 * ldsVel  -> definicoes de cargas para velocidades                  *
 * model     -> modelo da equacao de energia                         *
 * adv       -> tecnica da discretizacao do termo advecao            *
 * diff      -> tecnica da discretizacao do termo difusivo           *
 * vProp     -> propedades variaveis (true|false)                    *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lId       -> equa da celula                                       *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * fArea     -> area das faces                                       *
 * normal    -> vetores normais as faces das celulas                 *
 * volume    -> volume celula central                                *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os *
 *            centrois compartilhado nessa face da celula central    *
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceReK -> restricoes por elemento                                *
 * faceLdK -> carga por elemento                                     *
 * faceReVel-> restricoes por elemento                               *
 * faceLdVel-> carga por elemento                                    *
 * vel       -> campo de velocidade conhecido                        *
 * dField    -> matriz D do metodo simple                            *
 * lDensity  -> massa especifica com variacao temporal               *
 * lDviscosity-> viscosidade dinamica com variacao temporal          *
 * dField    -> matriz D do metodo simple                            * 
 * wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  * 
 * cDyn    -> coeficiente dinanmicos ( Ck, Ce)                       * 
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * lib       -> numero da biblioteca                                 *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellLibOneEqK(Loads *ldsK     , Loads *ldsVel
     , Turbulence *tModel                       
     , Advection  *adv             , Diffusion *diff              
     , short *RESTRICT lGeomType   , DOUBLE *RESTRICT lprop
     , INT   *RESTRICT lViz        , INT *RESTRICT lId
     , DOUBLE *RESTRICT ksi        , DOUBLE *RESTRICT mKsi
     , DOUBLE *RESTRICT eta        , DOUBLE *RESTRICT fArea
     , DOUBLE *RESTRICT normal     , DOUBLE *RESTRICT volume
     , DOUBLE *RESTRICT xm         , DOUBLE *RESTRICT xmcc
     , DOUBLE *RESTRICT dcca       , DOUBLE *RESTRICT cc
     , DOUBLE *RESTRICT vSkew      , DOUBLE *RESTRICT mvSkew
     , DOUBLE *RESTRICT lA         , DOUBLE *RESTRICT lB
     , DOUBLE *RESTRICT lRcell     , Temporal *ddt
     , short  *RESTRICT lFaceReK   , short  *RESTRICT lFaceLdK 
     , short  *RESTRICT lFaceReVel , short  *RESTRICT lFaceLdVel   
     , DOUBLE *RESTRICT u          , DOUBLE *RESTRICT gradU
     , DOUBLE *RESTRICT vel        , DOUBLE *RESTRICT gradVel
     , DOUBLE *RESTRICT pres       , DOUBLE *RESTRICT gradPres  
     , DOUBLE *RESTRICT lDensity   , DOUBLE *RESTRICT lDviscosity
     , DOUBLE *RESTRICT dField     , DOUBLE *RESTRICT wallPar    
     , DOUBLE *RESTRICT cDyn
     , short const nEn             , short  const nFace
     , short const ndm             , INT const nel)
{

/*... 2D*/
  if (ndm == 2);  
 /*..................................................................*/

/*... 3D*/
  else if (ndm == 3) 
    cellKinectTurb3D(ldsK        , ldsVel 
                    , tModel
                    , adv        , diff                      
                    , lGeomType  , lprop
                    , lViz       , lId
                    , ksi        , mKsi
                    , eta        , fArea
                    , normal     , volume
                    , xm         , xmcc
                    , dcca       , cc
                    , vSkew      , mvSkew 
                    , lA         , lB
                    , lRcell     , ddt
                    , lFaceReK   , lFaceLdK
                    , lFaceReVel , lFaceLdVel
                    , u          , gradU    
                    , vel        , gradVel
                    , pres       , gradPres  
                    , lDensity   , lDviscosity
                    , dField     , wallPar
                    , cDyn
                    , nEn        , nFace
                    , ndm        , nel);    
/*..................................................................*/
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 30/06/2016                                   *
 * Data de modificaco : 19/07/2018                                   *
 *-------------------------------------------------------------------*
 * CELLLIBSIMPLEVEl: chamada de bibliotecas de celulas para          *
 * problema de escoamento de fluidos (VEL)                           *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  *
 * loadsPres -> definicoes de cargas de pressao                      *
 * advVel    -> tecnica da discretizacao do termo advecao            *
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * typeSimple-> tipo do metodo simple                                *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lId       -> equa da celula                                       *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * fArea     -> area das faces                                       *
 * normal    -> vetores normais as faces das celulas                 *
 * volume    -> volume celula central                                *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os *
 *            centrois compartilhado nessa face da celula central    *
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               *
 * lDensity  -> massa especifica com variacao temporal               *
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceVelR  -> restricoes por elemento de velocidades               *
 * faceVelL  -> carga por elemento de velocidades                    *
 * facePresR -> restricoes por elemento de pressao                   *
 * facePresL -> carga por elemento de pressao                        *
 * pres      -> campo de pressao conhecido                           *
 * gradPes   -> gradiente reconstruido da pressao                    *
 * vel       -> campo de velocidade conhecido                        *
 * dField    -> matriz D do metodo simple                            *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * underU    -> fator underrelaxtion sinple                          *
 * sPressure -> reconstrucao de segunda ordem para pressoes nas      *
 *              faces                                                *
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * lib       -> numero da biblioteca                                 *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellLibSimpleVel(Loads *lVel        ,Loads *lPres     
             ,Advection  *advVel         ,Diffusion *diffVel    
             ,short const typeSimple 
             ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT lprop
             ,INT   *RESTRICT lViz       ,INT *RESTRICT lId  
             ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mKsi
             ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT fArea
             ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT volume
             ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
             ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity
             ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
             ,DOUBLE *RESTRICT lA        ,DOUBLE *RESTRICT lB
             ,DOUBLE *RESTRICT lRcell    ,Temporal *ddt
             ,short  *RESTRICT lFaceVelR ,short  *RESTRICT lFaceVelL
             ,short  *RESTRICT lFacePresR,short  *RESTRICT lFacePresL
             ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres 
             ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT gradVel
             ,DOUBLE *RESTRICT dField    ,DOUBLE *RESTRICT cc
             ,DOUBLE const underU        ,const bool sPressure
             ,short const nEn            ,short  const nFace
             ,short const ndm            ,short const lib
             ,INT const nel)
{

/*... */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2){
      cellSimpleVel2D(lVel    ,lPres       
                 ,advVel      ,diffVel
                 ,typeSimple
                 ,lGeomType   ,lprop
                 ,lViz        ,lId
                 ,ksi         ,mKsi
                 ,eta         ,fArea
                 ,normal      ,volume
                 ,xm          ,xmcc
                 ,dcca        ,lDensity 
                 ,vSkew       ,mvSkew
                 ,lA          ,lB
                 ,lRcell      ,*ddt 
                 ,lFaceVelR   ,lFaceVelL
                 ,lFacePresR  ,lFacePresL
                 ,pres        ,gradPres 
                 ,vel         ,gradVel
                 ,dField      ,cc
                 ,underU      ,sPressure
                 ,nEn         ,nFace 
                 ,ndm         ,nel);  
    }
/*..................................................................*/

/*... 3D*/
    else if(ndm == 3){
      cellSimpleVel3D(lVel    ,lPres       
                 ,advVel      ,diffVel     
                 ,typeSimple
                 ,lGeomType   ,lprop
                 ,lViz        ,lId
                 ,ksi         ,mKsi
                 ,eta         ,fArea
                 ,normal      ,volume
                 ,xm          ,xmcc
                 ,dcca        ,lDensity 
                 ,vSkew       ,mvSkew
                 ,lA          ,lB
                 ,lRcell      ,ddt 
                 ,lFaceVelR   ,lFaceVelL
                 ,lFacePresR  ,lFacePresL
                 ,pres        ,gradPres 
                 ,vel         ,gradVel
                 ,dField      ,cc
                 ,underU      ,sPressure
                 ,nEn         ,nFace 
                 ,ndm         ,nel);  
    } 
/*..................................................................*/
  }
/*..................................................................*/
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 05/08/2018                                   *
* Data de modificaco : 27/08/2019                                   *
*-------------------------------------------------------------------*
* cellLibCombustion: chamada de bibliotecas de celulas para         *
* problema de escoamento de fluidos (Energy)                        *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* loads     -> definicoes de cargas                                 *
* model     -> modelo da equacao de energia                         *
* adv       -> tecnica da discretizacao do termo advecao            *
* diff      -> tecnica da discretizacao do termo difusivo           *
* tModel    -> modelo de turbulencia                                *
* cModel    -> modelo de combustao                                  *
* vProp     -> propedades variaveis (true|false)                    *
* lGeomType -> tipo geometrico da celula central e seus vizinhos    *
* lprop     -> propriedade fisicas das celulas                      *
* lViz      -> viznhos da celula central                            *
* lId       -> equa da celula                                       *
* Ksi       -> vetores que unem centroide da celula central aos     *
*            vizinhos destas                                        *
* mKsi      -> modulo do vetor ksi                                  *
* eta       -> vetores paralelos as faces das celulas               *
* fArea     -> area das faces                                       *
* normal    -> vetores normais as faces das celulas                 *
* volume    -> volume celula central                                *
* xm        -> pontos medios das faces da celula central            *
* xmcc      -> vetores que unem o centroide aos pontos medios das   *
*            faces da celula central                                *
* vSkew     -> vetor entre o ponto medio a intersecao que une os    *
*            centrois compartilhado nessa face da celula central    *
* mvSkew    -> distacia entre o ponto medio a intersecao que une os *
*            centrois compartilhado nessa face da celula central    *
* dcca      -> menor distacia do centroide central a faces desta    *
*              celula                                               *
* cc        -> centroides da celula centra e seus vizinhos          *
* lA        -> nao definido                                         *
* lB        -> nao definido                                         *
* lRcell    -> nao definido                                         *
* ddt       -> discretizacao temporal                               *
* faceVelR  -> restricoes por elemento de velocidades               *
* faceVelL  -> carga por elemento de velocidades                    *
* facePresR -> restricoes por elemento de pressao                   *
* facePresL -> carga por elemento de pressao                        *
* vel       -> campo de velocidade conhecido                        *
* dField    -> matriz D do metodo simple                            *
* lDensity  -> massa especifica com variacao temporal               *
* lDiff     -> coeficiente de difusao das especies                  *
* lEddyVisc -> Viscosidade turbulenta                               *  
* dField    -> matriz D do metodo simple                            *
* wallPar   -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)  *
* underU    -> fator underrelaxtion sinple                          *
* nEn       -> numero de nos da celula central                      *
* nFace     -> numero de faces da celula central                    *
* ndm       -> numero de dimensoes                                  *
* lib       -> numero da biblioteca                                 *
* nel       -> numero da celula                                     *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* lA        -> coeficiente da linha i                               *
* lB        -> vetor de forca da linha i                            *
*-------------------------------------------------------------------*
*********************************************************************/
void cellLibCombustion(Loads *lComb        , Loads *lVel
               , Advection  *adv           , Diffusion *diff
               , Turbulence *tModel        , Combustion *cModel
               , PropVarFluid *vProp
               , short *RESTRICT lGeomType 
               , INT   *RESTRICT lViz      , INT *RESTRICT lId
               , DOUBLE *RESTRICT ksi      , DOUBLE *RESTRICT mKsi
               , DOUBLE *RESTRICT eta      , DOUBLE *RESTRICT fArea
               , DOUBLE *RESTRICT normal   , DOUBLE *RESTRICT volume
               , DOUBLE *RESTRICT xm       , DOUBLE *RESTRICT xmcc
               , DOUBLE *RESTRICT dcca     , DOUBLE *RESTRICT cc
               , DOUBLE *RESTRICT vSkew    , DOUBLE *RESTRICT mvSkew
               , DOUBLE *RESTRICT lA       , DOUBLE *RESTRICT lB
               , DOUBLE *RESTRICT lRcell   , Temporal *ddt
               , short  *RESTRICT lFaceR   , short  *RESTRICT lFaceL
               , short  *RESTRICT lFaceVelR, short  *RESTRICT lFaceVelL
               , DOUBLE *RESTRICT u        , DOUBLE *RESTRICT gradU
               , DOUBLE *RESTRICT Q        , DOUBLE *RESTRICT vel      
               , DOUBLE *RESTRICT pres     , DOUBLE *RESTRICT gradPres
               , DOUBLE *RESTRICT lDensity , DOUBLE *RESTRICT lDiff 
               , DOUBLE *RESTRICT lEddyVisc
               , DOUBLE *RESTRICT dField   , DOUBLE *RESTRICT wallPar
               , DOUBLE const underU
               , short const nEn           , short  const nFace
               , short const ndm           , short const lib
               , INT const nel)
{

/*... */
  if (lib == 1) {
/*... 2D*/
    if (ndm == 2);
/*    cellEnergy2D(lEnergy, lVel
        , adv, diff
        , tModel, model
        , vProp
        , lGeomType, lprop
        , lViz, lId
        , ksi, mKsi
        , eta, fArea
        , normal, volume
        , xm, xmcc
        , dcca, cc
        , vSkew, mvSkew
        , lA, lB
        , lRcell, ddt
        , lFaceR, lFaceL
        , lFaceVelR, lFaceVelL
        , u, gradU
        , vel, gradVel
        , pres, gradPres
        , lDensity, lSheat
        , lDviscosity, lTconductvity
        , dField, wallPar
        , underU
        , nEn, nFace
        , ndm, nel);*/
/*..................................................................*/

/*... 3D*/
    else if (ndm == 3) 
    {
      cellCombustion3D(lComb    , lVel
                     , adv       , diff
                     , tModel    , cModel
                     , vProp
                     , lGeomType 
                     , lViz      , lId
                     , ksi       , mKsi
                     , eta       , fArea
                     , normal    , volume
                     , xm        , xmcc
                     , dcca      , cc
                     , vSkew     , mvSkew
                     , lA        , lB
                     , lRcell    , ddt
                     , lFaceR    , lFaceL
                     , lFaceVelR , lFaceVelL
                     , u         , gradU
                     , Q         , vel       
                     , pres      , gradPres
                     , lDensity  , lDiff
                     , lEddyVisc
                     , dField    , wallPar
                     , underU
                     , nEn       , nFace
                     , ndm       , nel);
/*..................................................................*/
    }
/*..................................................................*/
  }

}
/*********************************************************************/


/*********************************************************************
 * Data de criacao    : 27/08/2017                                   *
 * Data de modificaco : 21/05/2019                                   *
 *-------------------------------------------------------------------*
 * CELLLIBSIMPLEVEl: chamada de bibliotecas de celulas para          *
 * problema de escoamento de fluidos (VEL -low mach)                 *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  *
 * loadsPres -> definicoes de cargas de pressao                      *
 * advVel    -> tecnica da discretizacao do termo advecao            *
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * ModelMomentum -> termos/modelos da equacao de momento linear      *
 * typeSimple-> tipo do metodo simple                                *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lViz      -> viznhos da celula central                            *
 * lId       -> equa da celula                                       *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * fArea     -> area das faces                                       *
 * normal    -> vetores normais as faces das celulas                 *
 * volume    -> volume celula central                                *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os *
 *            centrois compartilhado nessa face da celula central    *
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceVelR  -> restricoes por elemento de velocidades               *
 * faceVelL  -> carga por elemento de velocidades                    *
 * facePresR -> restricoes por elemento de pressao                   *
 * facePresL -> carga por elemento de pressao                        *
 * pres      -> campo de pressao conhecido                           *
 * gradPes   -> gradiente reconstruido da pressao                    *
 * vel       -> campo de velocidade conhecido                        *
 * lDensity  -> massa especifica com variacao temporal               *
 * lDviscosity-> viscosidade dinamica com variacao temporal          *
 * dField    -> matriz D do metodo simple                            *
 * stressR   -> tensor residual                                      *
 * densityMed-> densidade media do meio                              *
 * underU    -> fator underrelaxtion sinple                          *
 * sPressure -> reconstrucao de segunda ordem para pressoes nas      *
 *              faces                                                *
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * lib       -> numero da biblioteca                                 *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellLibSimpleVelLm(Loads *lVel     , Loads *lPres
           , Advection  *advVel         , Diffusion *diffVel    
           , Turbulence *tModel         , MomentumModel *ModelMomentum
           , short const typeSimple 
           , short *RESTRICT lGeomType  
           , INT   *RESTRICT lViz       , INT *RESTRICT lId   
           , DOUBLE *RESTRICT ksi       , DOUBLE *RESTRICT mKsi 
           , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT fArea 
           , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume 
           , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc 
           , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT cc 
           , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew 
           , DOUBLE *RESTRICT lA        , DOUBLE *RESTRICT lB 
           , DOUBLE *RESTRICT lRcell    , Temporal *ddt 
           , short  *RESTRICT lFaceVelR , short  *RESTRICT lFaceVelL 
           , short  *RESTRICT lFacePresR, short  *RESTRICT lFacePresL 
           , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres  
           , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT gradVel 
           , DOUBLE *RESTRICT lDensity  , DOUBLE *RESTRICT lViscosity 
           , DOUBLE *RESTRICT dField    , DOUBLE *RESTRICT stressR   
           , DOUBLE *RESTRICT wallPar   , DOUBLE const densityMed
           , DOUBLE const underU        , const bool sPressure 
           , short const nEn            , short  const nFace 
           , short const ndm            , short const lib 
           , INT const nel)
{

/*... */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2){
      cellSimpleVel2DLm(lVel      , lPres    
                     , *advVel    , *diffVel 
                     , *tModel    , *ModelMomentum
                     , typeSimple , lGeomType  
                     , lViz       , lId 
                     , ksi        , mKsi 
                     , eta        , fArea 
                     , normal     , volume 
                     , xm         , xmcc 
                     , dcca       , cc 
                     , vSkew      , mvSkew 
                     , lA         , lB 
                     , lRcell     , *ddt  
                     , lFaceVelR  , lFaceVelL 
                     , lFacePresR , lFacePresL 
                     , pres       , gradPres  
                     , vel        , gradVel 
                     , lDensity   , lViscosity 
                     , dField     , stressR
                     , wallPar    , densityMed
                     , underU     , sPressure 
                     , nEn        , nFace  
                     , ndm        , nel);  
    }
/*..................................................................*/

/*... 3D*/
    else if(ndm == 3)
    { 
      cellSimpleVel3DLm(lVel      , lPres    
                     , advVel     , diffVel 
                     , tModel     , ModelMomentum
                     , typeSimple , lGeomType  
                     , lViz       , lId 
                     , ksi        , mKsi 
                     , eta        , fArea 
                     , normal     , volume 
                     , xm         , xmcc 
                     , dcca       , cc 
                     , vSkew      , mvSkew 
                     , lA         , lB 
                     , lRcell     , ddt  
                     , lFaceVelR  , lFaceVelL 
                     , lFacePresR , lFacePresL 
                     , pres       , gradPres  
                     , vel        , gradVel 
                     , lDensity   , lViscosity 
                     , dField     , stressR
                     , wallPar    , densityMed
                     , underU     , sPressure 
                     , nEn        , nFace  
                     , ndm        , nel); 
    }  
/*..................................................................*/
  }

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 01/07/2016                                   *
 * Data de modificaco : 09/07/2016                                   *
 *-------------------------------------------------------------------*
 * CELLLIBSIMPLEPRES : chamada de bibliotecas de celulas para        *
 * problema de escoamento de fluidos (PRES)                          *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  *
 * loadsPres -> definicoes de cargas de pressao                      *
 * diffPres  -> tecnica da discretizacao do termo difusivo           *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lId       -> equa da celula                                       *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * fArea     -> area das faces                                       *
 * normal    -> vetores normais as faces das celulas                 *
 * volume    -> volume celula central                                *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os *
 *            centrois compartilhado nessa face da celula central    *
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               *
 * lDensity  -> massa especifica com variacao temporal               *
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceVelR  -> restricoes por elemento de velocidades               *
 * faceVelL  -> carga por elemento de velocidades                    *
 * facePresR -> restricoes por elemento de pressao                   *
 * facePresL -> carga por elemento de pressao                        *
 * pres      -> campo de pressao conhecido                           *
 * gradPes   -> gradiente reconstruido da pressao                    *
 * vel       -> campo de velocidade conhecido                        *
 * dField    -> matriz D do metodo simple                            *
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * lib       -> numero da biblioteca                                 *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellLibSimplePres(Loads *lVel       ,Loads *lPres
	             ,Diffusion *diffPres
               ,short *RESTRICT lGeomType,DOUBLE *RESTRICT lprop
               ,INT   *RESTRICT lViz     ,INT *RESTRICT lId  
               ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
               ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea
               ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
               ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
               ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
               ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
               ,DOUBLE *RESTRICT lRcell  ,Temporal *ddt
               ,short  *RESTRICT lFaceVelR ,short  *RESTRICT lFaceVelL
               ,short  *RESTRICT lFacePresR,short  *RESTRICT lFacePresL
               ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres 
               ,DOUBLE *RESTRICT vel     ,DOUBLE *RESTRICT dField 
               ,short const nEn          ,short  const nFace     
               ,short const ndm          ,short const lib    
               ,INT const nel)
{

/*... */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2){
      cellSimplePres2D(lVel,lPres
								 ,diffPres	
                 ,lGeomType,lprop
                 ,lViz     ,lId
                 ,ksi      ,mKsi
                 ,eta      ,fArea
                 ,normal   ,volume
                 ,xm       ,xmcc
                 ,dcca     ,lDensity 
                 ,vSkew    ,mvSkew
                 ,lA       ,lB
                 ,lRcell    
                 ,lFaceVelR   ,lFaceVelL
                 ,lFacePresR  ,lFacePresL
                 ,pres     ,gradPres 
                 ,vel      ,dField
                 ,nEn      ,nFace  
                 ,ndm      ,nel);  
    }
/*..................................................................*/

/*... 3D*/
    else if(ndm == 3){
      cellSimplePres3D(lVel,lPres
								 ,diffPres
                 ,lGeomType,lprop
                 ,lViz     ,lId
                 ,ksi      ,mKsi
                 ,eta      ,fArea
                 ,normal   ,volume
                 ,xm       ,xmcc
                 ,dcca     ,lDensity 
                 ,vSkew    ,mvSkew
                 ,lA       ,lB
                 ,lRcell    
                 ,lFaceVelR   ,lFaceVelL
                 ,lFacePresR  ,lFacePresL
                 ,pres     ,gradPres
                 ,vel      ,dField
                 ,nEn      ,nFace
                 ,ndm      ,nel);
    } 
/*..................................................................*/
  }

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 17/09/2017                                   *
 * Data de modificaco : 17/07/2018                                   *
 *-------------------------------------------------------------------*
 * CELLLIBSIMPLEPRESLM: chamada de bibliotecas de celulas para       *
 * problema de escoamento de fluidos a baixo Mach (PRES)             *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  *
 * loadsPres -> definicoes de cargas de pressao                      *
 * diffPres  -> tecnica da discretizacao do termo difusivo           *
 * eMass     -> termos/modelos da equacao de mass                    *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lId       -> equa da celula                                       *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * fArea     -> area das faces                                       *
 * normal    -> vetores normais as faces das celulas                 *
 * volume    -> volume celula central                                *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os *
 *            centrois compartilhado nessa face da celula central    *
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               *
 * lDensity  -> massa especifica com variacao temporal               *
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceVelR  -> restricoes por elemento de velocidades               *
 * faceVelL  -> carga por elemento de velocidades                    *
 * facePresR -> restricoes por elemento de pressao                   *
 * facePresL -> carga por elemento de pressao                        *
 * pres      -> campo de pressao conhecido                           *
 * gradPes   -> gradiente reconstruido da pressao                    *
 * vel       -> campo de velocidade conhecido                        *
 * dField    -> matriz D do metodo simple                            *
 * temp    -> temperatura                                            *   
 * wallPar -> parametros de parede  ( yPlus, uPlus, uFri,sTressW)    *
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * lib       -> numero da biblioteca                                 *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellLibSimplePresLm(Loads *lVel        , Loads *lPres
	             , Diffusion *diffPres        , MassEqModel *eMass  
               , short *RESTRICT lGeomType  
               , INT   *RESTRICT lViz       , INT *RESTRICT lId  
               , DOUBLE *RESTRICT ksi       , DOUBLE *RESTRICT mKsi
               , DOUBLE *RESTRICT eta       , DOUBLE *RESTRICT fArea
               , DOUBLE *RESTRICT normal    , DOUBLE *RESTRICT volume
               , DOUBLE *RESTRICT xm        , DOUBLE *RESTRICT xmcc
               , DOUBLE *RESTRICT dcca      , DOUBLE *RESTRICT lDensity
               , DOUBLE *RESTRICT vSkew     , DOUBLE *RESTRICT mvSkew
               , DOUBLE *RESTRICT lA        , DOUBLE *RESTRICT lB
               , DOUBLE *RESTRICT lRcell    , Temporal *ddt
               , short  *RESTRICT lFaceVelR , short  *RESTRICT lFaceVelL
               , short  *RESTRICT lFacePresR, short  *RESTRICT lFacePresL
               , DOUBLE *RESTRICT pres      , DOUBLE *RESTRICT gradPres 
               , DOUBLE *RESTRICT vel       , DOUBLE *RESTRICT dField 
               , DOUBLE *RESTRICT temp      , DOUBLE *RESTRICT wallPar
               , short const nEn            , short  const nFace     
               , short const ndm            , short const lib    
               , INT const nel)
{

/*... */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2){
      cellSimplePres2DLm(lVel, lPres
								 , diffPres	 , eMass 
                 , lGeomType 
                 , lViz      , lId
                 , ksi       , mKsi
                 , eta       , fArea
                 , normal    , volume
                 , xm        , xmcc
                 , dcca      , lDensity 
                 , vSkew     , mvSkew
                 , lA        , lB
                 , lRcell    , *ddt 
                 , lFaceVelR , lFaceVelL
                 , lFacePresR, lFacePresL
                 , pres      , gradPres 
                 , vel       , dField
                 , temp
                 , nEn       , nFace  
                 , ndm       , nel);  
    }
/*..................................................................*/

/*... 3D*/
    else if(ndm == 3){
      cellSimplePres3DLm(lVel, lPres
								 , diffPres	 , eMass 
                 , lGeomType 
                 , lViz      , lId
                 , ksi       , mKsi
                 , eta       , fArea
                 , normal    , volume
                 , xm        , xmcc
                 , dcca      , lDensity 
                 , vSkew     , mvSkew
                 , lA        , lB
                 , lRcell    , ddt 
                 , lFaceVelR , lFaceVelL
                 , lFacePresR, lFacePresL
                 , pres      , gradPres 
                 , vel       , dField
                 , temp      , wallPar
                 , nEn       , nFace  
                 , ndm       , nel);  
    } 
/*..................................................................*/
  }

}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 05/09/2016                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * CELLLIBVELEXP:chamada de bibliotecas de celulas para              *
 * problema de escoamento de fluidos (Vel - Explicito)               *
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * loadsVel  -> definicoes de cargas de velocidades                  *
 * loadsPres -> definicoes de cargas de pressao                      *
 * advVel    -> tecnica da discretizacao do termo advecao            *
 * diffVel   -> tecnica da discretizacao do termo difusivo           *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lId       -> equa da celula                                       *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * fArea     -> area das faces                                       *
 * normal    -> vetores normais as faces das celulas                 *
 * volume    -> volume celula central                                *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os *
 *            centrois compartilhado nessa face da celula central    *
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               *
 * lDensity  -> massa especifica com variacao temporal               *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceVelR  -> restricoes por elemento de velocidades               *
 * faceVelL  -> carga por elemento de velocidades                    *
 * facePresR -> restricoes por elemento de pressao                   *
 * facePresL -> carga por elemento de pressao                        *
 * pres      -> campo de pressao conhecido                           *
 * gradPes   -> gradiente reconstruido da pressao                    *
 * vel       -> campo de velocidade conhecido                        *
 * gradVel   -> gradiente rescontruido das velocidades               *
 * dField    -> matriz D do metodo simple                            *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * lBt       -> discretizacao temporal                               *
 * underU    -> parametro de sobre relaxamento                       *
 * nEn       -> numero de nos da celula central                      *
 * sPressure -> reconstrucao de segunda ordem para pressoes nas      *
 *              faces                                                *
 * fResidual -> calculo do residuo                                   *
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lB        -> velocidade (calculo explicito) ou residuo            *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellLibVelExp(Loads *loadsVel, Loads *loadsPres
         ,Advection advVel           ,Diffusion diffVel
         ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT lprop
         ,INT   *RESTRICT lViz      
         ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mKsi
         ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT fArea
         ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT volume
         ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
         ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity
         ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
         ,DOUBLE *RESTRICT lB        ,Temporal const ddt
         ,short  *RESTRICT lFaceVelR ,short  *RESTRICT lFaceVelL
         ,short  *RESTRICT lFacePresR,short  *RESTRICT lFacePresL
         ,DOUBLE *RESTRICT pres      ,DOUBLE *RESTRICT gradPres
         ,DOUBLE *RESTRICT vel       ,DOUBLE *RESTRICT gradVel
         ,DOUBLE *RESTRICT dField    ,DOUBLE *RESTRICT cc
         ,DOUBLE *RESTRICT lBt       ,DOUBLE const underU
         ,const bool sPressure       ,const bool fResidual
         ,const short nEn            ,short const nFace
         ,const short ndm            ,INT const nel)
         
{

/*... 2D*/
  if (ndm == 2) 
    cellVelExp2D(loadsVel ,loadsPres
                ,advVel    ,diffVel
                ,lGeomType ,lprop
                ,lViz
                ,ksi       ,mKsi
                ,eta       ,fArea
                ,normal    ,volume
                ,xm        ,xmcc
                ,dcca      ,lDensity
                ,vSkew     ,mvSkew
                ,lB        ,ddt
                ,lFaceVelR ,lFaceVelL
                ,lFacePresR,lFacePresL
                ,pres      ,gradPres
                ,vel       ,gradVel
                ,dField    ,cc
                ,lBt       ,underU   
                ,sPressure ,fResidual   
                ,nEn       ,nFace
                ,ndm       ,nel);
/*..................................................................*/

/*... 3D*/
    else if (ndm == 3) 
      cellVelExp3D(loadsVel ,loadsPres
                ,advVel    ,diffVel
                ,lGeomType ,lprop
                ,lViz
                ,ksi       ,mKsi
                ,eta       ,fArea
                ,normal    ,volume
                ,xm        ,xmcc
                ,dcca      ,lDensity
                ,vSkew     ,mvSkew
                ,lB        ,ddt
                ,lFaceVelR ,lFaceVelL
                ,lFacePresR,lFacePresL
                ,pres      ,gradPres
                ,vel       ,gradVel
                ,dField    ,cc
                ,lBt       ,underU   
                ,sPressure ,fResidual   
                ,nEn       ,nFace
                ,ndm       ,nel);
/*..................................................................*/
 
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 02/08/2016                                   *
 * Data de modificaco : 09/08/2016                                   *
 *-------------------------------------------------------------------*
 * CELLLIBSIMPLENONORTHPRES : chamada de bibliotecas de celulas para *
 * correcao nao ortogonal do problema de escoamento de fluidos (PRES)*
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 * diffPres  -> tecnica da discretizacao do termo difusivo           *
 * loadsPres -> definicoes de cargas de pressao                      *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    *
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lId       -> equa da celula                                       *
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        *
 * mKsi      -> modulo do vetor ksi                                  *
 * eta       -> vetores paralelos as faces das celulas               *
 * fArea     -> area das faces                                       *
 * normal    -> vetores normais as faces das celulas                 *
 * volume    -> volume celula central                                *
 * xm        -> pontos medios das faces da celula central            *
 * xmcc      -> vetores que unem o centroide aos pontos medios das   *
 *            faces da celula central                                *
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    *
 *            centrois compartilhado nessa face da celula central    *
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os *
 *            centrois compartilhado nessa face da celula central    *
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               *
 * lDensity  -> massa especifica com variacao temporal               *
 * lB        -> nao definido                                         *
 * gradPes   -> gradiente reconstruido da pressao                    *
 * dField    -> matriz D do metodo simple                            *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * nEn       -> numero de nos da celula central                      *
 * nFace     -> numero de faces da celula central                    *
 * ndm       -> numero de dimensoes                                  *
 * nel       -> numero da celula                                     *
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------*
 *********************************************************************/
void cellLibSimpleNonOrthPres(Diffusion diffPres
               ,short *RESTRICT lGeomType
               ,DOUBLE *RESTRICT lprop   ,INT   *RESTRICT lViz
               ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
               ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea
               ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
               ,DOUBLE *RESTRICT dcca    ,DOUBLE *RESTRICT lDensity
               ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
               ,DOUBLE *RESTRICT lB      
               ,DOUBLE *RESTRICT pres    ,DOUBLE *RESTRICT gradPres
               ,DOUBLE *RESTRICT dField  ,DOUBLE *RESTRICT cc
               ,short const nEn          ,short  const nFace     
               ,short const ndm          ,INT const nel)
{

/*... 2D*/
  if(ndm == 2){
    cellSimpleNonOrthPres2D(diffPres
                 ,lGeomType
                 ,lprop    ,lViz     
                 ,ksi      ,mKsi
                 ,eta      ,fArea
                 ,normal   ,volume
                 ,xm       ,xmcc
                 ,dcca     ,lDensity 
                 ,vSkew    ,mvSkew
                 ,lB       
                 ,pres     ,gradPres 
                 ,dField   ,cc
                 ,nEn      ,nFace  
                 ,ndm      ,nel);  
  }
/*..................................................................*/   

/*... 3D*/
  else if(ndm == 3){
    cellSimpleNonOrthPres3D(diffPres
                 ,lGeomType
                 ,lprop    ,lViz     
                 ,ksi      ,mKsi
                 ,eta      ,fArea
                 ,normal   ,volume
                 ,xm       ,xmcc
                 ,dcca     ,lDensity 
                 ,vSkew    ,mvSkew
                 ,lB       
                 ,pres     ,gradPres 
                 ,dField   ,cc
                 ,nEn      ,nFace  
                 ,ndm      ,nel);  
  } 
/*..................................................................*/   

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 22/08/2016                                   * 
 *-------------------------------------------------------------------* 
 * CELLLIBTRANS : chamada de bibliotecas de celulas para o problema  *
 * de transporte.                                                    * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads     -> definicoes de cargas                                 * 
 * advT      -> tecnica da discretizacao do termo advecao            *
 * diffT     -> tecnica da discretizacao do termo difusivo           *
 * tModel    -> configuracoes do modelo de transporte                *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            * 
 * lId       -> equa da celula                                       * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area das faces                                       * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> volume celula central                                * 
 * xm        -> pontos medios das faces da celula central            * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               * 
 * lDensity  -> massa especifica com variacao temporal               * 
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * faceR     -> restricoes por elemento                              * 
 * faceL     -> carga por elemento                                   * 
 * u0        -> solucao conhecida                                    * 
 * gradU0    -> gradiente rescontruido da solucao conhecida          * 
 * vel       -> campo de velocidade conhecido                        *
 * cc        -> centroides da celula centra e seus vizinhos          *
 * nEn       -> numero de nos da celula central                      * 
 * nFace     -> numero de faces da celula central                    * 
 * ndm       -> numero de dimensoes                                  * 
 * lib       -> numero da biblioteca                                 * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellLibTrans(Loads *loads           
               ,Advection *advT          ,Diffusion *diffT
               ,TransModel *tModel
               ,short *RESTRICT lGeomType,DOUBLE *RESTRICT lprop
               ,INT   *RESTRICT lViz     ,INT *RESTRICT lId  
               ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
               ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea
               ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
               ,DOUBLE *RESTRICT dcca    
               ,DOUBLE *RESTRICT lDensity,DOUBLE *RESTRICT lCoefDiff
               ,DOUBLE *RESTRICT vSkew   ,DOUBLE *RESTRICT mvSkew
               ,DOUBLE *RESTRICT lA      ,DOUBLE *RESTRICT lB
               ,DOUBLE *RESTRICT lRcell  ,Temporal *ddt
               ,short  *RESTRICT lFaceR  ,short  *RESTRICT lFaceL       
               ,DOUBLE *RESTRICT u0      ,DOUBLE *RESTRICT gradU0
               ,DOUBLE *RESTRICT vel     ,DOUBLE *RESTRICT cc
               ,short const nEn          ,short  const nFace     
               ,short const ndm          ,short const lib    
               ,INT const nel)
{

/*... */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2){
/*    cellTrans2D(loads    
                 ,advT     ,diffT
                 ,lGeomType,lprop
                 ,lViz     ,lId
                 ,ksi      ,mKsi
                 ,eta      ,fArea
                 ,normal   ,volume
                 ,xm       ,xmcc
                 ,dcca     ,lDensity 
                 ,vSkew    ,mvSkew
                 ,lA       ,lB
                 ,lRcell   ,ddt 
                 ,lFaceR   ,lFaceL
                 ,u0       ,gradU0      
                 ,vel      ,cc             
                 ,nEn      ,nFace 
                 ,ndm      ,nel);*/
    }
/*..................................................................*/   

/*... 3D*/
    else if(ndm == 3){
      cellTrans3D(loads    
                 ,advT     ,diffT
                 ,tModel
                 ,lGeomType,lprop
                 ,lViz     ,lId
                 ,ksi      ,mKsi
                 ,eta      ,fArea
                 ,normal   ,volume
                 ,xm       ,xmcc
                 ,dcca     
                 ,lDensity ,lCoefDiff
                 ,vSkew    ,mvSkew
                 ,lA       ,lB
                 ,lRcell   ,ddt 
                 ,lFaceR   ,lFaceL
                 ,u0       ,gradU0      
                 ,vel      ,cc             
                 ,nEn      ,nFace 
                 ,ndm      ,nel);  
    } 
/*..................................................................*/   
  }

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 13/05/2018                                   * 
 *-------------------------------------------------------------------* 
 * CELLLIBDIF : chamada de bibliotecas de celulas para o problema    *
 * de difusao.                                                       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads     -> definicoes de cargas                                 * 
 * diff      -> tecnica da discretizacao do termo difusivo           *
 * dModel   -> configuracoes do modelo difusivo                      *
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            * 
 * lId       -> equa da celula                                       * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area das faces                                       * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> volume celula central                                * 
 * xm        -> pontos medios das faces da celula central            * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               * 
 * lDensity  -> massa especifica com variacao temporal               * 
 * lCeofDiffD-> coeficiente de difusao com variacao temporal         *
 * lA        -> nao definido                                         *
 * lB        -> nao definido                                         *
 * lRcell    -> nao definido                                         *
 * ddt       -> discretizacao temporal                               *
 * u0        -> solucao conhecida                                    * 
 * gradU0    -> gradiente rescontruido da solucao conhecida          * 
 * faceR     -> restricoes por elemento                              * 
 * faceLd1   -> carga por elemento                                   * 
 * nEn       -> numero de nos da celula central                      * 
 * nFace     -> numero de faces da celula central                    * 
 * ndm       -> numero de dimensoes                                  * 
 * lib       -> numero da biblioteca                                 * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lA        -> coeficiente da linha i                               *
 * lB        -> vetor de forca da linha i                            *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellLibDif(Loads *loads               ,Diffusion *diff
               ,DiffModel *dModel          
               ,short *RESTRICT lGeomType  ,DOUBLE *RESTRICT lprop
               ,INT   *RESTRICT lViz       ,INT *RESTRICT lId  
               ,DOUBLE *RESTRICT ksi       ,DOUBLE *RESTRICT mKsi
               ,DOUBLE *RESTRICT eta       ,DOUBLE *RESTRICT fArea
               ,DOUBLE *RESTRICT normal    ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT xm        ,DOUBLE *RESTRICT xmcc
               ,DOUBLE *RESTRICT dcca      ,DOUBLE *RESTRICT lDensity
               ,DOUBLE *RESTRICT lCeofDiffD
               ,DOUBLE *RESTRICT vSkew     ,DOUBLE *RESTRICT mvSkew
               ,DOUBLE *RESTRICT lA        ,DOUBLE *RESTRICT lB
               ,DOUBLE *RESTRICT lRcell    ,Temporal *ddt
               ,short  *RESTRICT lFaceR    ,short  *RESTRICT lFaceL       
               ,DOUBLE *RESTRICT u0        ,DOUBLE *RESTRICT gradU0
               ,short const nEn            ,short  const nFace     
               ,short const ndm            ,short const lib    
               ,INT const nel)
{

/*... difusao pura */
  if(lib == 1){
/*... 2D*/
    if(ndm == 2){
/*    cellDif2D(loads
               ,lGeomType,lprop
               ,lViz     ,lId
               ,ksi      ,mKsi
               ,eta      ,fArea
               ,normal   ,volume
               ,xm       ,xmcc
               ,dcca     ,lDensity 
               ,vSkew    ,mvSkew
               ,lA       ,lB
               ,lRcell   ,ddt 
               ,lFaceR   ,lFaceL
               ,u0       ,gradU0      
               ,nEn      ,nFace 
               ,ndm      ,nel);*/
    }
/*..................................................................*/   

/*... 3D*/
    else if(ndm == 3){
      cellDif3D(loads     ,diff
               ,dModel    
               ,lGeomType ,lprop
               ,lViz      ,lId
               ,ksi       ,mKsi
               ,eta       ,fArea
               ,normal    ,volume
               ,xm        ,xmcc
               ,dcca      ,lDensity  
               ,lCeofDiffD
               ,vSkew     ,mvSkew
               ,lA        ,lB
               ,lRcell    ,ddt
               ,lFaceR    ,lFaceL 
               ,u0        ,gradU0      
               ,nEn       ,nFace 
               ,ndm       ,nel);
    } 
/*..................................................................*/   
  }

}
/*********************************************************************/
 
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * CELLGEOM2D : calculo geometrico de propriedade de celula 2D       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lx        -> coordenadas dos nos da celula central e seus viznhos * 
 * lnFace    -> numero da faces da celula central e seus vizinhos    * 
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * xc        -> indefenido                                           * 
 * ksi       -> indefinido                                           * 
 * mksi      -> indefinido                                           * 
 * eta       -> indefinido                                           * 
 * meta      -> indefinido                                           * 
 * normal    -> indefinido                                           * 
 * area      -> indefinido                                           * 
 * xm        -> indefinido                                           * 
 * xmcc      -> indefinido                                           * 
 * dcca      -> indefinido                                           * 
 * vSkew     -> indefinido                                           * 
 * mvSkew    -> indefinido                                           * 
 * sn        -> numeracao dos nos por aresta                         * 
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * numel     -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * xc     -> centroide das celulas                                   * 
 * Ksi    -> vetores que unem centroide da celula central aos        *
 *            vizinhos destas                                        * 
 * mKsi   -> modulo do vetor ksi                                     * 
 * eta    -> vetores paralelos as faces das celulas                  * 
 * fArea  -> modulo do vetor eta                                     * 
 * normal -> vetores normais as faces das celulas                    * 
 * area   -> area da celula central                                  * 
 * xm     -> pontos medios das faces da celula cenral                * 
 * xmcc   -> vetores que unem o centroide aos pontos medios das      * 
 *            faces da celula central                                * 
 * vSkew  -> vetor entre o ponto medio a intersecao que une os       * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew -> distacia entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca   -> menor distacia do centroide central a faces desta celula* 
 *-------------------------------------------------------------------* 
 * OBS: xc | x1 x2 x3 | -> celula vizinha da aresta 1                * 
 *         | x1 x2 x3 | -> celula vizinha da aresta 2                * 
 *         | x1 x2 x3 | -> celula vizinha da aresta 3                * 
 *         |   ...    | ->   ...                                     * 
 *         | x1 x2 x3 | -> celula central                            * 
 *                                                                   * 
 *     lx(vizinho,numero do no do vizinho, dimensao)                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellGeom2D(DOUBLE *RESTRICT lx       ,short *RESTRICT lnFace
               ,short  *RESTRICT lGeomType,DOUBLE *RESTRICT xc
               ,DOUBLE *RESTRICT ksi      ,DOUBLE *RESTRICT mksi
               ,DOUBLE *RESTRICT eta      ,DOUBLE *RESTRICT meta
               ,DOUBLE *RESTRICT normal   ,DOUBLE *RESTRICT area
               ,DOUBLE *RESTRICT xm       ,DOUBLE *RESTRICT xmcc
               ,DOUBLE *RESTRICT dcca     
               ,DOUBLE *RESTRICT vSkew    ,DOUBLE *RESTRICT mvSkew
               ,short  *RESTRICT sn
               ,short const maxNo         ,short const maxViz
               ,short const ndm           ,INT const nel)
{
  short i,j,k,cCell = maxViz;
  short no1,no2;
  double x1,x2;
/*...*/
  for(i=0;i<lnFace[cCell];i++){
    mksi[i]  = 0.0e0;
    meta[i]  = 0.0e0;
    dcca[i]  = 0.0e0;
    MAT2D(i,0,ksi,ndm)= 0.0e0; 
    MAT2D(i,1,ksi,ndm)= 0.0e0; 
    MAT2D(i,0,eta,ndm)= 0.0e0; 
    MAT2D(i,1,eta,ndm)= 0.0e0; 
  }
/*...................................................................*/

/*... calculo do centro geometrico das celulas*/      
  for(i=0;i<=maxViz;i++){
    for(j=0;j<ndm;j++){
      MAT2D(i,j,xc,ndm) = 0.0e0;
      if(lnFace[i]){
        for(k=0;k<lnFace[i];k++){
          MAT2D(i,j,xc,ndm) += MAT3D(i,k,j,lx,maxNo,ndm);
        }
        MAT2D(i,j,xc,ndm)*=(1.0e0/((DOUBLE)lnFace[i]));
      }
    }
  }
/*...................................................................*/

/*... vetor que une o centroide da celula ao vizinho(ksi)*/
  for(i=0;i<lnFace[cCell];i++){
    if(lnFace[i]){
      for(j=0;j<ndm;j++){
        MAT2D(i,j,ksi,ndm) = 
        MAT2D(i,j,xc,ndm) - MAT2D(cCell,j,xc,ndm);
        mksi[i] += MAT2D(i,j,ksi,ndm)*MAT2D(i,j,ksi,ndm); 
      }
      mksi[i] = sqrt(mksi[i]);
    }   
  }
/*...................................................................*/

/*... vetor paralelo as arestas (eta)*/
  for(i=0;i<lnFace[cCell];i++){
    no1 = MAT2D(i,0,sn,2); 
    no2 = MAT2D(i,1,sn,2);
    for(j=0;j<ndm;j++){
      x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
      x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
      MAT2D(i,j,eta,ndm) = x2 - x1;
    }
    *area = areaCell(eta,lGeomType[cCell],ndm,nel);
  }
/*...................................................................*/

/*... modulo do vetor das arestas (eta)*/
  for(i=0;i<lnFace[cCell];i++){
    for(j=0;j<ndm;j++){
      meta[i] += MAT2D(i,j,eta,ndm)*MAT2D(i,j,eta,ndm); 
    }
    meta[i] = sqrt(meta[i]);
  }
/*...................................................................*/

/*... calculo dos vetores unitarios*/
  for(i=0;i<lnFace[cCell];i++){
    for(j=0;j<ndm;j++){
      MAT2D(i,j,eta,ndm) = MAT2D(i,j,eta,ndm) / meta[i];
      if(lnFace[i]) 
        MAT2D(i,j,ksi,ndm) = MAT2D(i,j,ksi,ndm) / mksi[i];
    } 
  }
/*...................................................................*/

/*... calculo do vetor normal unitario*/
  for(i=0;i<lnFace[cCell];i++){
    MAT2D(i,0,normal,ndm) =  MAT2D(i,1,eta,ndm);
    MAT2D(i,1,normal,ndm) = -MAT2D(i,0,eta,ndm);
  }
/*...................................................................*/  

/*... ponto medio da aresta(xm) e vetor que une este ponto ao 
      centroide(mxcc) */
  for(i=0;i<lnFace[cCell];i++){
    no1 = MAT2D(i,0,sn,2); 
    no2 = MAT2D(i,1,sn,2);
    for(j=0;j<ndm;j++){
      x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
      x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
/*...*/
      MAT2D(i,j,xm  ,ndm) = 0.5e0*(x1+x2);
/*...................................................................*/  

/*...*/
      MAT2D(i,j,xmcc,ndm) = MAT2D(i,j,xm,ndm) - MAT2D(cCell,j,xc,ndm);   
/*...................................................................*/  

/*...*/
      dcca[i] += MAT2D(i,j,xmcc,ndm)*MAT2D(i,j,normal,ndm); 
/*...................................................................*/  
    }
  }
/*...................................................................*/  

/*... calculo do centro da celula fantasma*/      
/*  for(i=0;i<=maxViz;i++){
    for(j=0;j<ndm;j++){
      if(!lnFace[i]){
        MAT2D(i,j,xc,ndm) = MAT2D(cCell,j,xc,ndm) 
                          + 2.0e0*dcca[i]*MAT2D(i,j,normal,ndm);
      }
    }
  }*/
/*...................................................................*/

/*... vetor que une o centroide da celula ao vizinho fantasma(ksi)*/
  for(i=0;i<lnFace[cCell];i++){
    if(!lnFace[i]){
      for(j=0;j<ndm;j++){
        MAT2D(i,j,ksi,ndm) = MAT2D(i,j,xmcc,ndm);
        mksi[i] += MAT2D(i,j,ksi,ndm)*MAT2D(i,j,ksi,ndm); 
      }
      mksi[i] = sqrt(mksi[i]);
      for(j=0;j<ndm;j++)
        MAT2D(i,j,ksi,ndm) /= mksi[i];
    }   
  }
/*...................................................................*/

/*...*/
  vectorKm2d(lx     , xc
            ,xm     
            ,vSkew  ,mvSkew
            ,sn     ,lnFace[cCell]
            ,maxViz ,maxNo       
            ,ndm    ,nel);
/*...................................................................*/  
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 21/08/2019                                   * 
 *-------------------------------------------------------------------* 
 * CELLGEOM3D : calculo geometrico de propriedade de celula 3D       * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lx        -> coordenadas dos nos da celula central e seus viznhos * 
 * lGeomType -> tipo geometrico da celula central e seus vizinhos    * 
 * lnEn      -> numero da de nos da celula central e seus vizinhos   * 
 * lnFace    -> numero da faces da celula central e seus vizinhos    * 
 * xc        -> indefenido                                           * 
 * ksi       -> indefinido                                           * 
 * mksi      -> indefinido                                           * 
 * eta       -> indefinido                                           * 
 * fArea     -> indefinido                                           * 
 * normal    -> indefinido                                           * 
 * area      -> indefinido                                           * 
 * xm        -> indefinido                                           * 
 * xmcc      -> indefinido                                           * 
 * dcca      -> indefinido                                           * 
 * mvSkew    -> indefinido                                           * 
 * sn        -> numeracao dos nos por aresta                         * 
 * maxNo     -> numero de nos por celula maximo da malha             * 
 * maxViz    -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * numel     -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * xc     -> centroide das celulas                                   * 
 * Ksi    -> vetores que unem centroide da celula central aos        *
 *            vizinhos destas                                        * 
 * mKsi   -> modulo do vetor ksi                                     * 
 * eta    -> vetores paralelos as faces das celulas                  * 
 * fArea  -> area da face                                            * 
 * normal -> vetores normais as faces das celulas                    * 
 * volume -> volume da celula central                                * 
 * xm     -> pontos medios das faces da celula cenral                * 
 * xmcc   -> vetores que unem o centroide aos pontos medios das      * 
 *            faces da celula central                                * 
 * vSkew  -> vetor entre o ponto medio a intersecao que une os       * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew -> distacia entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * dcca   -> menor distacia do centroide central a faces desta celula* 
 *-------------------------------------------------------------------* 
 * OBS: xc | x1 x2 x3 | -> celula vizinha da aresta 1                * 
 *         | x1 x2 x3 | -> celula vizinha da aresta 2                * 
 *         | x1 x2 x3 | -> celula vizinha da aresta 3                * 
 *         |   ...    | ->   ...                                     * 
 *         | x1 x2 x3 | -> celula central                            * 
 *                                                                   * 
 *     lx(vizinho,numero do no do vizinho, dimensao)                 * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellGeom3D(DOUBLE *RESTRICT lx       ,short  *RESTRICT lGeomType
               ,short *RESTRICT lnFace    ,short *RESTRICT lnEn      
               ,DOUBLE *RESTRICT xc        
               ,DOUBLE *RESTRICT ksi      ,DOUBLE *RESTRICT mksi
               ,DOUBLE *RESTRICT eta      ,DOUBLE *RESTRICT fArea
               ,DOUBLE *RESTRICT normal   ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT xm       ,DOUBLE *RESTRICT xmcc
               ,DOUBLE *RESTRICT dcca     
               ,DOUBLE *RESTRICT vSkew    ,DOUBLE *RESTRICT mvSkew
               ,short  *RESTRICT sn                         
               ,short const maxNo         ,short const maxViz
               ,short const ndm           ,INT const nel)
{
  short i,j,k,cCell = maxViz;
  short no1,no2,no3,no4,nenFaces=0;
  DOUBLE x1,x2,x3,x4,v1[3],v2[3],n[3],dot;

  for(i=0;i<lnFace[cCell];i++){
    mksi[i]  = 0.0e0;
    fArea[i] = 0.0e0;
    dcca[i]  = 0.0e0;
    MAT2D(i,0,ksi,ndm)= 0.0e0; 
    MAT2D(i,1,ksi,ndm)= 0.0e0; 
    MAT2D(i,2,ksi,ndm)= 0.0e0; 
    MAT2D(i,0,eta,ndm)= 0.0e0; 
    MAT2D(i,1,eta,ndm)= 0.0e0; 
    MAT2D(i,2,eta,ndm)= 0.0e0; 
  }

/*...*/
  if( lGeomType[cCell] == HEXACELL || lGeomType[cCell] == PIRACELL)
    nenFaces = 4;  
  else if( lGeomType[cCell] == TETRCELL )  
    nenFaces = 3;   
/*...................................................................*/

/*... calculo do centro geometrico das celulas*/
  for(i=0;i<= maxViz;i++){
    for(j=0;j<ndm;j++){
      MAT2D(i,j,xc,ndm) = 0.0e0;
      if(lnFace[i]){
        for(k=0;k<lnEn[i];k++){
          MAT2D(i,j,xc,ndm) += MAT3D(i,k,j,lx,maxNo,ndm);
        }
        MAT2D(i,j,xc,ndm)*=(1.0e0/((DOUBLE)lnEn[i]));
      }
    }
  }
/*...................................................................*/

/*... vetor que une o centroide da celula ao vizinho(ksi)*/
  for(i=0;i<lnFace[cCell];i++){
    if(lnFace[i]){
      for(j=0;j<ndm;j++){
        MAT2D(i,j,ksi,ndm) = 
        MAT2D(i,j,xc,ndm) - MAT2D(cCell,j,xc,ndm);
        mksi[i] += MAT2D(i,j,ksi,ndm)*MAT2D(i,j,ksi,ndm); 
      }
      mksi[i] = sqrt(mksi[i]);
      for(j=0;j<ndm;j++)
        MAT2D(i,j,ksi,ndm) /=mksi[i];
    }   
  }
/*...................................................................*/

/*... calculo do vetor normal e area da face*/
  if( lGeomType[cCell] == HEXACELL ){
    for(i=0;i<lnFace[cCell];i++){
      no1 = MAT2D(i,0,sn,nenFaces);
      no2 = MAT2D(i,1,sn,nenFaces);
      no3 = MAT2D(i,2,sn,nenFaces);
      no4 = MAT2D(i,3,sn,nenFaces);
      for(j=0;j<ndm;j++){
        x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
        x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
        x3 = MAT3D(cCell,no3,j,lx,maxNo,ndm);
        x4 = MAT3D(cCell,no4,j,lx,maxNo,ndm);
        v1[j] = MAT2D(0,j,eta,ndm) = x2 - x1;
        MAT2D(1,j,eta,ndm) = x3 - x2;
        MAT2D(2,j,eta,ndm) = x4 - x3;
        v2[j] = MAT2D(3,j,eta,ndm) = x1 - x4;
      }
/*... area da face*/
      fArea[i] = areaCell(eta,QUADCELL,ndm,nel);
/*... normal a face*/
      prodVet(v2,v1,n);
/*... normal do vetor normal*/
      dot = 0.e0;
      for(j=0;j<ndm;j++){
        dot += n[j]*n[j];                             
      }
      dot = sqrt(dot);
/*... vetor unitario normal*/
      for(j=0;j<ndm;j++)
        MAT2D(i,j,normal,ndm) =  n[j]/dot;
    }
  }
/*...................................................................*/

/*... calculo do vetor normal e area da face*/
  else if( lGeomType[cCell] == TETRCELL ){
    for(i=0;i<lnFace[cCell];i++){
      no1 = MAT2D(i,0,sn,nenFaces);
      no2 = MAT2D(i,1,sn,nenFaces);
      no3 = MAT2D(i,2,sn,nenFaces);
      for(j=0;j<ndm;j++){
        x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
        x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
        x3 = MAT3D(cCell,no3,j,lx,maxNo,ndm);
        v1[j] = MAT2D(0,j,eta,ndm) = x2 - x1;
        v2[j] = MAT2D(1,j,eta,ndm) = x3 - x1;
      }
/*... area da face*/
      fArea[i] = areaCell(eta,TRIACELL,ndm,nel);
/*... normal a face*/
      prodVet(v1,v2,n);
/*... normal do vetor normal*/
      dot = 0.e0;
      for(j=0;j<ndm;j++){
        dot += n[j]*n[j];                             
      }
      dot = sqrt(dot);
/*... vetor unitario normal*/
      for(j=0;j<ndm;j++)
        MAT2D(i,j,normal,ndm) =  n[j]/dot;
    }
  }
/*...................................................................*/

/*... calculo do vetor normal e area da face*/
  else if (lGeomType[cCell] == PIRACELL)
  {
    for (i = 0; i<lnFace[cCell]; i++) {
/*... face triangular*/
      if (i)
      {
        no1 = MAT2D(i, 0, sn, nenFaces);
        no2 = MAT2D(i, 1, sn, nenFaces);
        no3 = MAT2D(i, 2, sn, nenFaces);
        for (j = 0; j<ndm; j++) {
          x1 = MAT3D(cCell, no1, j, lx, maxNo, ndm);
          x2 = MAT3D(cCell, no2, j, lx, maxNo, ndm);
          x3 = MAT3D(cCell, no3, j, lx, maxNo, ndm);
          v1[j] = MAT2D(0, j, eta, ndm) = x2 - x1;
          v2[j] = MAT2D(1, j, eta, ndm) = x3 - x1;
        }
/*... area da face*/
        fArea[i] = areaCell(eta, TRIACELL, ndm, nel);
/*... normal a face*/
        prodVet(v1, v2, n);
/*... normal do vetor normal*/
        dot = 0.e0;
        for (j = 0; j<ndm; j++) {
          dot += n[j] * n[j];
        }
        dot = sqrt(dot);
/*... vetor unitario normal*/
        for (j = 0; j<ndm; j++)
          MAT2D(i, j, normal, ndm) = n[j] / dot;
      }
/*... face quandrigular*/
      else
      {
        no1 = MAT2D(i, 0, sn, nenFaces);
        no2 = MAT2D(i, 1, sn, nenFaces);
        no3 = MAT2D(i, 2, sn, nenFaces);
        no4 = MAT2D(i, 3, sn, nenFaces);
        for (j = 0; j<ndm; j++) 
        {
          x1 = MAT3D(cCell, no1, j, lx, maxNo, ndm);
          x2 = MAT3D(cCell, no2, j, lx, maxNo, ndm);
          x3 = MAT3D(cCell, no3, j, lx, maxNo, ndm);
          x4 = MAT3D(cCell, no4, j, lx, maxNo, ndm);
          v1[j] = MAT2D(0, j, eta, ndm) = x2 - x1;
          MAT2D(1, j, eta, ndm) = x3 - x2;
          MAT2D(2, j, eta, ndm) = x4 - x3;
          v2[j] = MAT2D(3, j, eta, ndm) = x1 - x4;
        }
/*... area da face*/
        fArea[i] = areaCell(eta, QUADCELL, ndm, nel);
/*... normal a face*/
        prodVet(v2, v1, n);
/*... normal do vetor normal*/
        dot = 0.e0;
        for (j = 0; j<ndm; j++) {
          dot += n[j] * n[j];
        }
        dot = sqrt(dot);
/*... vetor unitario normal*/
        for (j = 0; j<ndm; j++)
          MAT2D(i, j, normal, ndm) = n[j] / dot;
      }
    }
  }
/*...................................................................*/

/*... ponto medio da aresta(xm) e vetor que une este ponto ao 
      centroide(mxcc) */
  if( lGeomType[cCell] == HEXACELL ){
    for(i=0;i<lnFace[cCell];i++){
      no1 = MAT2D(i,0,sn,nenFaces);
      no2 = MAT2D(i,1,sn,nenFaces);
      no3 = MAT2D(i,2,sn,nenFaces);
      no4 = MAT2D(i,3,sn,nenFaces);
      for(j=0;j<ndm;j++){
        x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
        x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
        x3 = MAT3D(cCell,no3,j,lx,maxNo,ndm);
        x4 = MAT3D(cCell,no4,j,lx,maxNo,ndm);
/*...*/
        MAT2D(i,j,xm  ,ndm) = 0.25e0*(x1+x2+x3+x4);
/*...................................................................*/  

/*...*/
        MAT2D(i,j,xmcc,ndm) = MAT2D(i,j,xm,ndm) - MAT2D(cCell,j,xc,ndm);   
/*...................................................................*/  

/*...*/
        dcca[i] += MAT2D(i,j,xmcc,ndm)*MAT2D(i,j,normal,ndm); 
/*...................................................................*/  
      }
    }
/*...................................................................*/  
  }
/*...................................................................*/  

/*...*/  
  else if( lGeomType[cCell] == TETRCELL ){
    for(i=0;i<lnFace[cCell];i++){
      no1 = MAT2D(i,0,sn,nenFaces);
      no2 = MAT2D(i,1,sn,nenFaces);
      no3 = MAT2D(i,2,sn,nenFaces);
      for(j=0;j<ndm;j++){
        x1 = MAT3D(cCell,no1,j,lx,maxNo,ndm);
        x2 = MAT3D(cCell,no2,j,lx,maxNo,ndm);
        x3 = MAT3D(cCell,no3,j,lx,maxNo,ndm);
/*...*/
        MAT2D(i,j,xm  ,ndm) =  D1DIV3*(x1+x2+x3);
/*...................................................................*/  

/*...*/
        MAT2D(i,j,xmcc,ndm) = MAT2D(i,j,xm,ndm) - MAT2D(cCell,j,xc,ndm);   
/*...................................................................*/  

/*...*/
        dcca[i] += MAT2D(i,j,xmcc,ndm)*MAT2D(i,j,normal,ndm); 
/*...................................................................*/  
      }
    }
/*...................................................................*/  
  }
/*...................................................................*/  

/*...*/
  else if (lGeomType[cCell] == PIRACELL) {
    for (i = 0; i<lnFace[cCell]; i++) {
/*... faces triangulares*/
      if(i)
      {
        no1 = MAT2D(i, 0, sn, nenFaces);
        no2 = MAT2D(i, 1, sn, nenFaces);
        no3 = MAT2D(i, 2, sn, nenFaces);
        for (j = 0; j<ndm; j++) {
          x1 = MAT3D(cCell, no1, j, lx, maxNo, ndm);
          x2 = MAT3D(cCell, no2, j, lx, maxNo, ndm);
          x3 = MAT3D(cCell, no3, j, lx, maxNo, ndm);
/*...*/
          MAT2D(i, j, xm, ndm) = D1DIV3 * (x1 + x2 + x3);
/*...................................................................*/

/*...*/
          MAT2D(i, j, xmcc, ndm) = MAT2D(i, j, xm, ndm) 
                                 - MAT2D(cCell, j, xc, ndm);
/*...................................................................*/

/*...*/
          dcca[i] += MAT2D(i, j, xmcc, ndm)*MAT2D(i, j, normal, ndm);
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/

/*... face quadrada*/
      else
      {
        no1 = MAT2D(i, 0, sn, nenFaces);
        no2 = MAT2D(i, 1, sn, nenFaces);
        no3 = MAT2D(i, 2, sn, nenFaces);
        no4 = MAT2D(i, 3, sn, nenFaces);
        for (j = 0; j<ndm; j++)
        {
          x1 = MAT3D(cCell, no1, j, lx, maxNo, ndm);
          x2 = MAT3D(cCell, no2, j, lx, maxNo, ndm);
          x3 = MAT3D(cCell, no3, j, lx, maxNo, ndm);
          x4 = MAT3D(cCell, no4, j, lx, maxNo, ndm);
/*...*/
          MAT2D(i, j, xm, ndm) = 0.25e0*(x1 + x2 + x3 + x4);
/*...................................................................*/

/*...*/
          MAT2D(i, j, xmcc, ndm) = MAT2D(i, j, xm, ndm) 
                                 - MAT2D(cCell, j, xc, ndm);
/*...................................................................*/

/*...*/
          dcca[i] += MAT2D(i, j, xmcc, ndm)*MAT2D(i, j, normal, ndm);
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
    
/*... calculo do centro da celula fantasma*/      
/*  for(i=0;i<=maxViz;i++){
    for(j=0;j<ndm;j++){
      if(!lnFace[i]){
        MAT2D(i,j,xc,ndm) = MAT2D(cCell,j,xc,ndm) 
                          + 2.0e0*dcca[i]*MAT2D(i,j,normal,ndm);
      }
    }
  }*/
/*...................................................................*/

/*... vetor que une o centroide da celula ao vizinho fantasma(ksi)*/
  for(i=0;i<lnFace[cCell];i++){
    if(!lnFace[i]){
      mksi[i] = 0.e0;
      for(j=0;j<ndm;j++){
        MAT2D(i,j,ksi,ndm) = MAT2D(i,j,xmcc,ndm);
        mksi[i] += MAT2D(i,j,ksi,ndm)*MAT2D(i,j,ksi,ndm); 
      }
      mksi[i] = sqrt(mksi[i]);
      for(j=0;j<ndm;j++)
        MAT2D(i,j,ksi,ndm) /= mksi[i];
    }   
  }
/*...................................................................*/

/*... volume da celula*/
  if( lGeomType[cCell] == HEXACELL )
    *volume = volumeHexa(&MAT3D(cCell,0,0,lx,maxNo,ndm));
  else if( lGeomType[cCell] == TETRCELL )
    *volume = volumeTetra(&MAT3D(cCell,0,0,lx,maxNo,ndm));
  else if (lGeomType[cCell] = PIRACELL)
    *volume = volume3DGreenGauss(xm,normal,fArea,lnFace[cCell]);
/*...................................................................*/

  if (*volume <= 0.e0)
  {
    fprintf(stderr,"Volume negativo da celula %d!",nel);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"
           ,__FILE__,__func__,__LINE__);
    exit(EXIT_FAILURE);
  }

/*...*/
   vectorKm3d(xc            ,xm
             ,ksi           ,normal
             ,vSkew         ,mvSkew     
             ,lnFace[cCell] ,ndm
             ,maxViz        ,nel);
/*...................................................................*/  
}
/*********************************************************************/ 
    
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * CELLRCGRAD : chamada de bibliotecas de celulas para a reconstrucao*
 * de gradiente.                                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads     -> definicoes de cargas                                 * 
 * lprop     -> propriedade fisicas das celulas                      *
 * lViz      -> viznhos da celula central                            *
 * lSquare -> matriz para a reconstrucao least Square                * 
 * lSquareR-> fatoracao R (RCLSQUAREQR)                              * 
 * Ksi       -> vetores que unem centroide da celula central aos     *
 *            vizinhos destas                                        * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area da face                                         * 
 * normal    -> vetores normais as faces das celulas                 * 
 * area      -> area da celula central                               * 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face                      * 
 * xm        -> pontos medios das faces das celulas                  * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * dcca      -> menor distacia do centroide central a faces desta    *
 *              celula                                               * 
 * u         -> solucao conhecida na celula                          * 
 * gradU     -> gradiente rescontruido da solucao conhecida          * 
 * nU        -> solucao conhecida no no                              * 
 * faceR     -> restricoes por elmento                               * 
 * faceL     -> carga por elemento                                   * 
 * nFace     -> numero de faces                                      * 
 * ndm       -> numero de dimensoes                                  * 
 * lib       -> tipo de reconstrucao de gradiente                    * 
 * ndf       -> grauss de liberdade                                  * 
 * isNod     -> numeracao dos nos por aresta                         * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * gradU     -> gradiente calculodo                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void cellLibRcGrad(Loads *loads
                 ,INT   *RESTRICT lViz    ,DOUBLE *RESTRICT lProp    
                 ,DOUBLE *RESTRICT lLsquare,DOUBLE *RESTRICT lLsquareR
                 ,DOUBLE *RESTRICT ksi     ,DOUBLE *RESTRICT mKsi
                 ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea
                 ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
                 ,DOUBLE *RESTRICT vSkew   
                 ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
                 ,DOUBLE *RESTRICT lDcca 
                 ,short  *RESTRICT lFaceR  ,short *RESTRICT lFaceL
                 ,DOUBLE *RESTRICT u       ,DOUBLE *RESTRICT gradU 
                 ,DOUBLE *RESTRICT lnU     ,short const ty                
                 ,short const nFace        ,short const ndm      
                 ,short const lib          ,short const ndf
                 ,short *RESTRICT  isNod   ,INT const nel){
  long aux;
    
  switch(lib){
/*... green-Gauss linear baseado na celula*/  
    case RCGRADGAUSSC:
      greenGaussCell(loads
                    ,lViz    ,mKsi
                    ,lProp   ,lDcca
                    ,eta     ,fArea
                    ,normal  ,volume
                    ,vSkew   
                    ,xm      ,xmcc
                    ,lFaceR  ,lFaceL
                    ,u       ,gradU 
                    ,nFace   ,ndm   
                    ,ndf     ,nel);
    break;
/*...................................................................*/ 
    
/*... green-Gauss linear baseado no no*/  
    case RCGRADGAUSSN:
       greenGaussNode(lViz    ,fArea
                     ,normal  ,volume
                     ,lnU     ,gradU 
                     ,isNod   
                     ,nFace   ,ndm   
                     ,ndf     ,ty);
    break;
/*...................................................................*/ 
   
/*... minimos quadrados*/  
    case RCLSQUARE:
      leastSquare(loads
                 ,lLsquare,lViz
                 ,xm      ,xmcc
                 ,lProp   ,lDcca
                 ,u       ,gradU
                 ,lFaceR  ,lFaceL
                 ,nFace   ,ndf
                 ,ndm     ,nel);
    break;
/*...................................................................*/ 

/*... minimos quadrados QR*/  
    case RCLSQUAREQR:
      leastSquareQR(loads
                   ,lLsquare,lLsquareR
                   ,lProp   ,lDcca
                   ,lViz    
                   ,xm      ,xmcc
                   ,u       ,gradU
                   ,lFaceR  ,lFaceL
                   ,nFace   ,ndf
                   ,ndm);
    break;
/*...................................................................*/ 

/*...*/
    default: 
      aux = (long) nel;
      printf("Erro: tipo de reconstrucao nao exitentes. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
               ,aux,__FILE__,__func__);
      exit(EXIT_FAILURE);
    break; 
/*...................................................................*/ 

  }  

}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 11/05/2019                                   * 
 *-------------------------------------------------------------------* 
 * GRREENGAUSSCELL: reconstrucao de gradiente green-gauss linear por * 
 * celula                                                            *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * loads     -> definicoes de cargas                                 * 
 * lViz      -> viznhos da celula central                            * 
 * lProp     -> propriedades dos material                            *
 * lDcca     -> menor distancia do centroide a faces desta celula    * 
 * mKsi      -> modulo do vetor ksi                                  * 
 * eta       -> vetores paralelos as faces das celulas               * 
 * fArea     -> area da face                                         * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> area da celula central                               * 
 * vSkew      -> vetor entre o ponto medio a intersecao que une os   * 
 *            centrois compartilhado nessa face                      * 
 * xm        -> pontos medios das faces das celulas (x1,x2,x3,t)     * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * lFaceR    -> restricoes por elmento                               * 
 * lFaceS    -> carga por elemento                                   * 
 * u         -> solucao conhecida                                    * 
 * gradU     -> gradiente rescontruido da solucao conhecida          * 
 * nFace     -> numero vizinhos por celula maximo da malha           * 
 * ndf       -> grauss de liberdade                                  * 
 * ndm       -> numero de dimensoes                                  * 
 * nel       -> numero da celula                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * gradU     -> gradiente calculado                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void greenGaussCell(Loads *loads
               ,INT *RESTRICT lViz   ,DOUBLE *RESTRICT mKsi
               ,DOUBLE *RESTRICT lProp   ,DOUBLE *RESTRICT lDcca 
               ,DOUBLE *RESTRICT eta     ,DOUBLE *RESTRICT fArea
               ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT vSkew   
               ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc 
               ,short  *RESTRICT lFaceR  ,short *RESTRICT lFaceL
               ,DOUBLE *RESTRICT u       ,DOUBLE *RESTRICT gradU 
               ,short const nFace        ,short const ndm   
               ,short const ndf          ,INT const nel)
{
  DOUBLE v,dPviz,coefDif,tmp;
  DOUBLE uf[MAX_NUM_FACE*MAX_NDF],uC[MAX_NDF],par[MAXLOADPARAMETER];
  DOUBLE lModKsi,alpha,alphaMenosUm,invVol,xx[4];
  INT vizNel;
  short idCell = nFace,nCarg,type;
  short i,j,k;

	invVol = 1.e0/volume[idCell];
/*...*/

/*... */
  if(ndf == 1){
    uC[0] = u[idCell];
    for(i=0;i<nFace;i++){
      vizNel = lViz[i];
/*... dominio*/
      if(vizNel > -1){
        lModKsi    = mKsi[i];
/*...*/
        dPviz = 0.e0;
        for(j=0;j<ndm;j++){
          v          = MAT2D(i,j,vSkew,ndm) + MAT2D(i,j,xmcc,ndm);
          dPviz     += v*v;
        }
/*...................................................................*/

/*...*/
        dPviz        = sqrt(dPviz);
        alpha        = dPviz/lModKsi;
        alphaMenosUm = 1.e0 - alpha;
/*...................................................................*/

/*...*/
        uf[i] = alphaMenosUm*uC[0] + alpha*u[i];
/*...................................................................*/
      }
/*...................................................................*/

/*... contorno*/
      else{
        xx[0] = MAT2D(i,0,xm,ndm);
        xx[1] = MAT2D(i,1,xm,ndm);
        xx[2] = 0.e0;
        xx[3] = 0.e0;
        if(ndm == 3) xx[2] = MAT2D(i,2,xm,ndm);
/*... temperatura prescrita na face(extrapolacao linear)*/
        if(lFaceR[i]){
          nCarg= lFaceL[i]-1;
          type = loads[nCarg].type;
/*... valor prescrito*/
          if (type == DIRICHLETBC) {
            getLoads(par,&loads[nCarg],xx);
            uf[i] =par[0];
          }
/*...................................................................*/

/*...*/
          else if (type == INLET) {
            getLoads(par,&loads[nCarg],xx);
            uf[i] =par[0];
          }
/*...................................................................*/

/*... fluxo prescrito*/
          else if (type == NEUMANNBC ){
            coefDif = lProp[COEFDIF];
            if (coefDif != 0.e0) {
              getLoads(par,&loads[nCarg],xx);
              uf[i] = uC[0] - (par[0]/coefDif)*lDcca[i];
            }              
          }
/*...................................................................*/

/*... condicao de robin*/
          else if(type == ROBINBC || type == CONVECTIONHEAT ){
            uf[i] =uC[0];
            for(j=0;j<ndm;j++)
              uf[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
          }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
          else if (type == OUTLET){
            uf[i] = uC[0];
          }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
          else if(type == SINBC){
            loadSenProd(&uf[i],loads[nCarg].par,xx);
          }
/*...................................................................*/
        } 
/*...................................................................*/

/*...*/
        else{
          uf[i] =uC[0];
          for(j=0;j<ndm;j++)
            uf[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/
    for(k=0;k<ndm;k++){
      tmp = 0.e0;
      for(i=0;i<nFace;i++){
        tmp += uf[i]*fArea[i]*MAT2D(i,k,normal,ndm);; 
      gradU[k] = tmp*invVol; 
      }
    }
/*...................................................................*/
  } 
/*...................................................................*/

/*... graus de liberdade maior que 1*/
  else{
/*...*/
    for(k=0;k<ndf;k++)
      uC[k] = MAT2D(idCell,k,u,ndf);
/*...................................................................*/
    for(i=0;i<nFace;i++){
      vizNel = lViz[i];
/*... dominio*/
      if(vizNel > -1){
        lModKsi    = mKsi[i];
/*...*/
        dPviz = 0.e0;
        for(j=0;j<ndm;j++){
          v          = MAT2D(i,j,vSkew,ndm) + MAT2D(i,j,xmcc,ndm);
          dPviz     += v*v;
        }
/*...................................................................*/

/*...*/
        dPviz        = sqrt(dPviz);
        alpha        = dPviz/lModKsi;
        alphaMenosUm = 1.e0 - alpha;
/*...................................................................*/

/*...*/
        for(k=0;k<ndf;k++)
          MAT2D(i,k,uf,ndf) = alphaMenosUm*uC[k] 
                            + alpha*MAT2D(i,k,u,ndf);
/*...................................................................*/
      }
/*...................................................................*/

/*... contorno*/
      else{
/*... temperatura prescrita na face(extrapolacao linear)*/
        if(lFaceR[i] > 0){
          nCarg= lFaceL[i]-1;
          type = loads[nCarg].type;
/*... valor prescrito*/
          if( type == DIRICHLETBC || type == MOVEWALL){
            getLoads(par,&loads[nCarg],xx);
            for(k=0;k<ndf;k++)
              MAT2D(i,k,uf,ndf) = par[k];            
          }
/*...................................................................*/

/*...*/          
          else if( type == INLET ){
            getLoads(par,&loads[nCarg],xx);
            for(k=0;k<ndf;k++)
              MAT2D(i,k,uf,ndf) = par[k];      
          }
/*...................................................................*/

/*... fluxo prescrito*/
          else if (type == NEUMANNBC ){
            coefDif = lProp[COEFDIF];
            if (coefDif != 0.e0) {
              getLoads(par,&loads[nCarg],xx);
              for(k=0;k<ndf;k++)
                MAT2D(i,k,uf,ndf) = uC[k] 
                             - (par[k]/coefDif)*lDcca[i];
            }
              
          }
/*...................................................................*/

/*... condicao de robin*/
          else if(type == ROBINBC || type == CONVECTIONHEAT)
          {
            ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__
            ,"Condicao de robin na implementada",EXIT_FAILURE);
          }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
          else if (type == OUTLET || type == INLETSTAICTPRES 
                || type == OPEN){
            for(k=0;k<ndf;k++)
              MAT2D(i,k,uf,ndf) = uC[k];
          }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
          else if(type == SINBC){
            ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__
            ,"Condicao SINBC robin nao implementada",EXIT_PROG);
          }
/*...................................................................*/
        } 
/*...................................................................*/

/*... parede static impermeavel*/
        else if(lFaceR[i]==STATICWALL){
          for(k=0;k<ndf;k++)
            MAT2D(i,k,uf,ndf) = 0.e0;
        }
/*...................................................................*/

/*... fluxo*/
        else{
          for(k=0;k<ndf;k++){
            MAT2D(i,k,uf,ndf) = uC[k];
            for(j=0;j<ndm;j++)
              MAT2D(i,k,uf,ndf) += 
              MAT2D(k,j,gradU,ndm)*MAT2D(i,j,xmcc,ndm);
          }
        }
/*...................................................................*/
      }
/*...................................................................*/
    }
/*...................................................................*/

/*...*/    
    for(k=0;k<ndf;k++){
      for(j=0;j<ndm;j++){
        tmp = 0.e0;
        for(i=0;i<nFace;i++)
          tmp += MAT2D(i,k,uf,ndf)*fArea[i]*MAT2D(i,j,normal,ndm); 
        MAT2D(k,j,gradU,ndm) = tmp*invVol; 
      }
    }
/*...................................................................*/
  }
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 30/04/2018                                   * 
 *-------------------------------------------------------------------* 
 * GRREENGAUSSNODE: reconstrucao de gradiente green-gauss linear por * 
 * celula                                                            *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * lViz      -> viznhos da celula central                            * 
 * fArea     -> area do face                                         * 
 * normal    -> vetores normais as faces das celulas                 * 
 * volume    -> area da celula central                               * 
 * lFaceR    -> restricoes por elmento                               * 
 * lFaceS    -> carga por elemento                                   * 
 * u         -> solucao conhecida nodal                              * 
 * gradU     -> gradiente rescontruido da solucao conhecida por      *
 *              celula                                               * 
 * isNod     -> numeracao dos nos por aresta                         * 
 * nFace     -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * ndf       -> grauss de liberdade                                  * 
 * ty        -> tipo geometrico da celula                            * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * gradU     -> gradiente calculado                                  *
 *-------------------------------------------------------------------* 
 *********************************************************************/
void greenGaussNode(INT *RESTRICT lViz   ,DOUBLE *RESTRICT fArea
               ,DOUBLE *RESTRICT normal  ,DOUBLE *RESTRICT volume
               ,DOUBLE *RESTRICT u       ,DOUBLE *RESTRICT gradU 
               ,short *RESTRICT isNod       
               ,short const nFace        ,short const ndm   
               ,short const ndf          ,short const ty)
{
  DOUBLE lNormal[MAX_NDM],aux1;
  DOUBLE uf[MAX_NDF];
  DOUBLE invVol;
  INT no[4];
  short idCell = nFace;
  short i,j,k;
  
  invVol = 1.e0/volume[idCell];
  
  for(i=0;i<ndf*ndm;i++)
    gradU[i]  = 0.e0;
  
  for(i=0;i<ndm;i++)
    lNormal[i]  = 0.e0;
  
  for(i=0;i<MAX_NDF;i++)
    uf[i]  = 0.e0;
  
/*...*/
  for(i=0;i<nFace;i++)
  {    
    for(j=0;j<ndm;j++)
      lNormal[j] = MAT2D(i,j,normal,ndm);
/*... triangulos e quadrilateros*/    
    if(ty == TRIACELL || ty == QUADCELL)
    {
      no[0] = MAT2D(i,0,isNod,2);
      no[1] = MAT2D(i,1,isNod,2);
      for(j=0;j<ndf;j++)
        uf[j] = 0.5e0*(MAT2D(no[0],j,u,ndf) + MAT2D(no[1],j,u,ndf));
    }
/*...................................................................*/

/*... tetraedro*/    
    else if(ty == TETRCELL)
    {
      no[0] = MAT2D(i,0,isNod,3);
      no[1] = MAT2D(i,1,isNod,3);
      no[2] = MAT2D(i,2,isNod,3);
      for(j=0;j<ndf;j++)
        uf[j] = D1DIV3*(MAT2D(no[0],j,u,ndf) 
                       +MAT2D(no[1],j,u,ndf)
                       +MAT2D(no[2],j,u,ndf));
    }
/*...................................................................*/

/*... hexaedro*/    
    else if(ty == HEXACELL)
    {
      no[0] = MAT2D(i,0,isNod,4);
      no[1] = MAT2D(i,1,isNod,4);
      no[2] = MAT2D(i,2,isNod,4);
      no[3] = MAT2D(i,3,isNod,4);
      for(j=0;j<ndf;j++)
        uf[j] = 0.25e0*(MAT2D(no[0],j,u,ndf) 
                       +MAT2D(no[1],j,u,ndf)
                       +MAT2D(no[2],j,u,ndf)
                       +MAT2D(no[3],j,u,ndf));
    }
/*...................................................................*/

/*... piramide*/
    else if (ty == PIRACELL)
    {
/*... face triangular*/
      if (i)
      {
        no[0] = MAT2D(i, 0, isNod, 4);
        no[1] = MAT2D(i, 1, isNod, 4);
        no[2] = MAT2D(i, 2, isNod, 4);
        for (j = 0; j<ndf; j++)
          uf[j] = D1DIV3 * (MAT2D(no[0], j, u, ndf)
            + MAT2D(no[1], j, u, ndf)
            + MAT2D(no[2], j, u, ndf));
      }
/*...................................................................*/

/*... face quadrangular*/
      else
      {
        no[0] = MAT2D(i, 0, isNod, 4);
        no[1] = MAT2D(i, 1, isNod, 4);
        no[2] = MAT2D(i, 2, isNod, 4);
        no[3] = MAT2D(i, 3, isNod, 4);
        for (j = 0; j<ndf; j++)
          uf[j] = 0.25e0*(MAT2D(no[0], j, u, ndf)
                        + MAT2D(no[1], j, u, ndf)
                        + MAT2D(no[2], j, u, ndf)
                        + MAT2D(no[3], j, u, ndf));
      }
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    for(j=0;j<ndf;j++)
    {
      aux1 = uf[j]*fArea[i]; 
      for(k=0;k<ndm;k++)
        MAT2D(j,k,gradU,ndm) += aux1*lNormal[k]; 
    }
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
  for(i=0;i<ndf;i++)
    for(j=0;j<ndm;j++)
      MAT2D(i,j,gradU,ndm) *= invVol; 
/*...................................................................*/

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 11/09/2019                                   * 
 *-------------------------------------------------------------------* 
 * LEASTSQUARE : calcula o gradiente por minimos quadrados           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * loads     -> definicoes de cargas                                 * 
 * lLsquare  -> matriz para a reconstrucao least Square              * 
 * xm        -> pontos medios das faces das celulas                  * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * lProp     -> propriedades dos material                            *
 * lDcca     -> menor distancia do centroide a faces desta celula    * 
 * lViz      -> viznhos da celula central                            * 
 * u         -> solucao conhecida                                    * 
 * gradU     -> gradiente rescontruido da solucao conhecida          * 
 * lFaceR    -> restricoes por elmento                               * 
 * lFaceS    -> carga por elemento                                   * 
 * nFace     -> numero vizinhos por celula maximo da malha           * 
 * ndf       -> grauss de liberdade                                  * 
 * ndm       -> numero de dimensoes                                  * 
 * nel       -> numero do elemento                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lSquare   -> matriz para a reconstrucao least Square              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void  leastSquare(Loads *loads
                 ,DOUBLE *RESTRICT lLsquare,INT *RESTRICT lViz
                 ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc 
                 ,DOUBLE *RESTRICT lProp   ,DOUBLE *RESTRICT lDcca 
                 ,DOUBLE *RESTRICT u       ,DOUBLE *RESTRICT gradU
                 ,short  *RESTRICT lFaceR  ,short *RESTRICT lFaceL
                 ,short const nFace        ,short const ndf
                 ,short const ndm          ,INT const nel){

  DOUBLE du[MAX_NUM_FACE*MAX_NDF],uC[MAX_NDF];
  DOUBLE uT[MAX_NDF],xx[4],par[MAXLOADPARAMETER];
  DOUBLE tmp,coefDif;
  INT vizNel;
  short idCell = nFace,nCarg,type;
  short i,j,k,l,it=1,itNumber=8;
  
  for(j=0;j<MAX_NDF;j++)
    uC[j] = 0.e0;
 
  for(l=0;l<it;l++){    
/*... um grau de liberdade*/  
    if(ndf == 1){
      uC[0] = u[idCell];  
      for(i=0;i<nFace;i++){
        du[i]  = 0.e0;
        vizNel = lViz[i];
/*... dominio*/
        if(vizNel > -1)
          du[i] = u[i] - uC[0];
/*... contorno*/
        else{
          xx[0] = MAT2D(i,0,xm,ndm);
          xx[1] = MAT2D(i,1,xm,ndm);
          xx[2] = 0.e0;
          if(ndm == 3) xx[2] = MAT2D(i,2,xm,ndm);
          xx[3] = 0.e0;          
/*... temperatura prescrita na face(extrapolacao linear)*/
          if(lFaceR[i])
          {
            nCarg = lFaceL[i]-1;
            type  = loads[nCarg].type;
/*... valor prescrito*/
            if (type == DIRICHLETBC || type == MOVEWALL) {
              getLoads(par,&loads[nCarg],xx);
              du[i] = par[0] - uC[0];
            }              
/*...................................................................*/

/*...*/
            else if (type == INLET) {
              getLoads(par,&loads[nCarg],xx);
              du[i] = par[0] - uC[0];
            }                          
/*...................................................................*/

/*... fluxo prescrito*/
            else if (type == NEUMANNBC){
              coefDif = lProp[COEFDIF];
              if (coefDif != 0.e0) {
                getLoads(par,&loads[nCarg],xx);
                du[i] = (par[0]/coefDif)*lDcca[i];
              }
                
            }
/*...................................................................*/

/*... condicao de robin*/
            else if( type == ROBINBC || type == CONVECTIONHEAT )
            {
              if(l==0) it = itNumber;
              du[i] =0.e0;
              for(j=0;j<ndm;j++)
                du[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
            }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
            else if (type == OUTLET || type == INLETSTAICTPRES
                  || type == OPEN){
              du[i] = 0.e0;
            }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
            else if( type == SINBC ){
              loadSenProd(uT,loads[nCarg].par,xx);
              du[i] = uT[0] - uC[0];
            }
/*...................................................................*/
          } 
/*...................................................................*/

/*... fluxo nulo*/
          else
          {
            if(l==0) it = itNumber;
            du[i] =0.e0;
            for(j=0;j<ndm;j++)
              du[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/

/*... gradU = G du*/
      for(i=0;i<ndm;i++){
        tmp = 0.e0;
        for(j=0;j<nFace;j++)
          tmp += MAT2D(i,j,lLsquare,nFace)*du[j]; 
        gradU[i] = tmp;
      }
/*...................................................................*/
    }
/*...................................................................*/

/*... graus de liberdade maior que 1 */
    else
    {
      for(k=0;k<ndf;k++)
        uC[k] = MAT2D(idCell,k,u,ndf);
/*... loop nas faces*/
      for(i=0;i<nFace;i++)
        {
        for(k=0;k<ndf;k++)
          MAT2D(i,k,du,ndf) = 0.e0;    
        vizNel = lViz[i];
/*... dominio*/
        if(vizNel > -1)
          for(k=0;k<ndf;k++)
            MAT2D(i,k,du,ndf) =  MAT2D(i,k,u,ndf) - uC[k];
/*... contorno*/
        else
        {
          xx[0] = MAT2D(i,0,xm,ndm);
          xx[1] = MAT2D(i,1,xm,ndm);
          xx[2] = 0.e0;
          if(ndm == 3) xx[2] = MAT2D(i,2,xm,ndm);
          xx[3] = 0.e0;
/*... potencial prescrito na face(extrapolacao linear)*/
          if(lFaceR[i] > 0)
          {
            nCarg= lFaceL[i]-1;
            type = loads[nCarg].type;
/*... valor prescrito*/
            if( type == DIRICHLETBC || type == MOVEWALL)
            {
              getLoads(par,&loads[nCarg],xx);
              for(k=0;k<ndf;k++)
                MAT2D(i,k,du,ndf) = par[k] - uC[k];
            }
/*...................................................................*/

/*...*/
            else if( type == INLET)
            {
              getLoads(par,&loads[nCarg],xx);
              for(k=0;k<ndf;k++)
                MAT2D(i,k,du,ndf) = par[k] - uC[k];
            }
/*...................................................................*/

/*... condicao de robin*/
            else if(type == ROBINBC || type == CONVECTIONHEAT){
              ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__
              ,"Condicao de robin na implementada",EXIT_PROG);
            }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
            else if(type == SINBC){
              ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__
              ,"Condicao SINBC robin nao implementada",EXIT_PROG);
            }
/*...................................................................*/
         
/*... derivada nula (condicao localmente parabolica saida)*/
            else if (type == OUTLET || type == INLETSTAICTPRES
                  || type == OPEN)
              for(k=0;k<ndf;k++)
                MAT2D(i,k,du,ndf) = 0.e0;                          
/*...................................................................*/
          } 
/*...................................................................*/

/*... parede static impermeavel*/
          else if(lFaceR[i]==STATICWALL)
          {
            for(k=0;k<ndf;k++)
              MAT2D(i,k,du,ndf) = - uC[k];
           }
/*...................................................................*/

/*... fluxo nulo*/
          else
          {
            if(l==0) it = itNumber;
            for(k=0;k<ndf;k++)
            {
              MAT2D(i,k,du,ndf) = 0.e0;
              for(j=0;j<ndm;j++)
                MAT2D(i,k,du,ndf) += 
                MAT2D(k,j,gradU,ndm)*MAT2D(i,j,xmcc,ndm);
            }
          }
/*...................................................................*/
        } 
/*...................................................................*/
      }
/*...................................................................*/

/*... gradU = G du*/
      for(k=0;k<ndf;k++)
        for(i=0;i<ndm;i++){
          tmp = 0.e0;
          for(j=0;j<nFace;j++)
            tmp += MAT2D(i,j,lLsquare,nFace)*MAT2D(j,k,du,ndf); 
          MAT2D(k,i,gradU,ndm) = tmp;
        }
/*...................................................................*/
    }
/*...................................................................*/
  }
/*...................................................................*/
} 
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 11/05/2019                                   * 
 *-------------------------------------------------------------------* 
 * LEASTSQUAREQR: calcula o gradiente por minimos quadrados           *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * loads     -> definicoes de cargas                                 * 
 * lLsQt     -> matriz para a reconstrucao least Square              * 
 * lLsR      -> fatoracao R (RCLSQUAREQR)                            * 
 * lProp     -> propriedades dos material                            *
 * lDcca     -> menor distancia do centroide a faces desta celula    * 
 * xm        -> pontos medios das faces das celulas (x1,x2,x3)       * 
 * xmcc      -> vetores que unem o centroide aos pontos medios das   * 
 *            faces da celula central                                * 
 * lViz      -> viznhos da celula central                            * 
 * u         -> solucao conhecida                                    * 
 * gradU     -> gradiente rescontruido da solucao conhecida          * 
 * lFaceR    -> restricoes por elmento                               * 
 * lFaceS    -> carga por elemento                                   * 
 * nFace     -> numero vizinhos por celula maximo da malha           * 
 * ndm       -> numero de dimensoes                                  * 
 * ndf       -> grauss de liberdade                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lSquare   -> matriz para a reconstrucao least Square              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void  leastSquareQR(Loads *loads
                   ,DOUBLE *RESTRICT lLsQt   ,DOUBLE *RESTRICT lLsR
                   ,DOUBLE *RESTRICT lProp   ,DOUBLE *RESTRICT lDcca 
                   ,INT *RESTRICT lViz       
                   ,DOUBLE *RESTRICT xm      ,DOUBLE *RESTRICT xmcc
                   ,DOUBLE *RESTRICT u       ,DOUBLE *RESTRICT gradU
                   ,short  *RESTRICT lFaceR  ,short *RESTRICT lFaceL
                   ,short const nFace        ,short const ndf
                   ,short const ndm){

  DOUBLE du[MAX_NUM_FACE*MAX_NDF],uC[MAX_NDF],coefDif;
  DOUBLE uT[MAX_NDF],xx[4],par[MAXLOADPARAMETER];
  DOUBLE b[MAX_NUM_FACE*MAX_NDF],b0,b1,b2;
  INT vizNel;
  short idCell = nFace,type;
  short i,j,k,l,nCarg;

  for(j=0;j<MAX_NDF;j++)
    uC[j] = 0.e0;
 
  for(l=0;l<1;l++){
/*... um grau de liberdade*/  
    if(ndf == 1){
      uC[0] = u[idCell];  
      for(i=0;i<nFace;i++){
        du[i]  = 0.e0;
        vizNel = lViz[i];
/*... dominio*/
        if(vizNel > -1)
          du[i] = u[i] - uC[0];
/*... contorno*/
        else{
          xx[0] = MAT2D(i,0,xm,ndm);
          xx[1] = MAT2D(i,1,xm,ndm);
          xx[2] = 0.e0;
          if(ndm == 3) xx[2] = MAT2D(i,2,xm,ndm);
          xx[3] = 0.e0;
/*... temperatura prescrita na face(extrapolacao linear)*/
          if(lFaceR[i]){
            nCarg=lFaceL[i]-1;
            type = loads[nCarg].type;
/*... valor prescrito*/
            if (type == DIRICHLETBC || type == MOVEWALL) {
              getLoads(par,&loads[nCarg],xx);
              du[i] = par[0] - uC[0];
            }              
/*...................................................................*/

/*...*/
            else if (type == INLET) {
              getLoads(par,&loads[nCarg],xx);
              du[i] = par[0] - uC[0];
            }                          
/*...................................................................*/

/*... fluxo prescrito*/
            else if (type == NEUMANNBC){
              coefDif = lProp[COEFDIF];
              if (coefDif != 0.e0) {
                getLoads(par,&loads[nCarg],xx);
                du[i] = (par[0]/coefDif)*lDcca[i];
              }
                
            }
/*...................................................................*/

/*... condicao de robin*/
            else if( type == ROBINBC || type == CONVECTIONHEAT ){
              du[i] =0.e0;
              for(j=0;j<ndm;j++)
                du[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
            }
/*...................................................................*/

/*... derivada nula (condicao localmente parabolica saida)*/
            else if (type == OUTLET || type == OPEN ){
              du[i] = 0.e0;
            }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
            else if( type == SINBC ){
              loadSenProd(uT,loads[nCarg].par,xx);
              du[i] = uT[0] - uC[0];
            }
/*...................................................................*/
          } 
/*...................................................................*/

/*... fluxo nulo*/
          else{
            du[i] =0.e0;
            for(j=0;j<ndm;j++)
              du[i] += gradU[j]*MAT2D(i,j,xmcc,ndm);
          }
/*...................................................................*/
        }
/*...................................................................*/
      }
/*...................................................................*/

/*... b = Qtdu*/
      for(i=0;i<ndm;i++){
        b[i]      = 0.e0;
        for(j=0;j<nFace;j++){
          b[i] += MAT2D(i,j,lLsQt,nFace)*du[j]; 
        } 
      }
/*...................................................................*/
    
/*... R grad = Qtdu*/
      if(ndm == 2){
/*... retrosubstituicao*/
        gradU[1] = b[1]/lLsR[2];
        gradU[0] = (b[0] - lLsR[1]*gradU[1])/lLsR[0];
      }
      else if(ndm == 3){
/*... retrosubstituicao*/
        gradU[2] = b[2]/lLsR[5];
        gradU[1] = (b[1] - lLsR[4]*gradU[2])/lLsR[3];
        gradU[0] = (b[0] 
                 - lLsR[1]*gradU[1] 
                 - lLsR[2]*gradU[2])/lLsR[0];
      }
/*...................................................................*/
    }
/*...................................................................*/
  
/*... graus de liberdade maior que 1 */
    else{
      for(k=0;k<ndf;k++)
        uC[k] = MAT2D(idCell,k,u,ndf);
/*... loop nas faces*/
      for(i=0;i<nFace;i++){
        vizNel = lViz[i];
/*... dominio*/
        if(vizNel > -1)
          for(k=0;k<ndf;k++)
            MAT2D(i,k,du,ndf) =  MAT2D(i,k,u,ndf) - uC[k];
/*... contorno*/
        else{
          xx[0] = MAT2D(i,0,xm,ndm);
          xx[1] = MAT2D(i,1,xm,ndm);
          xx[2] = 0.e0;
          if(ndm == 3) xx[2] = MAT2D(i,2,xm,ndm);
          xx[3] = 0.e0;
/*... potencial prescrito na face(extrapolacao linear)*/
          if(lFaceR[i] > 0){
            nCarg= lFaceL[i]-1;
            type = loads[nCarg].type;
/*... valor prescrito*/
            if( type == DIRICHLETBC || type == MOVEWALL){
              getLoads(par,&loads[nCarg],xx);
              for(k=0;k<ndf;k++)
                MAT2D(i,k,du,ndf) = par[k] - uC[k];
            }
/*...................................................................*/

/*...*/
            else if( type == INLET){
              getLoads(par,&loads[nCarg],xx);
              for(k=0;k<ndf;k++)
                MAT2D(i,k,du,ndf) = par[k] - uC[k];
            }
/*...................................................................*/

/*... condicao de robin*/
            else if(type == ROBINBC || type == CONVECTIONHEAT){
              ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__
              ,"Condicao de robin na implementada",EXIT_PROG);
            }
/*...................................................................*/

/*... potencial senoidal prescrito 
      ((a1*sin(w1*x1+c1))*(a2*sin(w2*x2+c2)*(a3*sin(w3*x3+c3))*/
            else if(type == SINBC){
              ERRO_GERAL(fileLogDebug,__FILE__,__func__,__LINE__
              ,"Condicao SINBC robin nao implementada",EXIT_PROG);
            }
/*...................................................................*/
         
/*... derivada nula (condicao localmente parabolica saida)*/
            else if (type == OUTLET || type == INLETSTAICTPRES
                    || type == OPEN)
              for(k=0;k<ndf;k++)
                MAT2D(i,k,du,ndf) = 0.e0;                          
/*...................................................................*/
          } 
/*...................................................................*/

/*... parede static impermeavel*/
          else if(lFaceR[i]==STATICWALL){
            for(k=0;k<ndf;k++)
              MAT2D(i,k,du,ndf) = - uC[k];
           }
/*...................................................................*/

/*... fluxo nulo*/
          else{
            for(k=0;k<ndf;k++){
              MAT2D(i,k,du,ndf) = 0.e0;
              for(j=0;j<ndm;j++)
                MAT2D(i,k,du,ndf) += 
                MAT2D(k,j,gradU,ndm)*MAT2D(i,j,xmcc,ndm);
            }
          }
/*...................................................................*/
        } 
/*...................................................................*/
      }
/*...................................................................*/

/*... b = Qtdu*/
      for(k=0;k<ndf;k++)
        for(i=0;i<ndm;i++){
          MAT2D(k,i,b,ndf) = 0.e0;
          for(j=0;j<nFace;j++)
            MAT2D(k,i,b,ndf) += MAT2D(i,j,lLsQt,nFace)*MAT2D(j,k,du,ndf); 
        }
/*...................................................................*/
    
/*... R grad = Qtdu*/
      for(k=0;k<ndf;k++){
        if(ndm == 2){
          b0 =  MAT2D(k,0,b,ndf);
          b1 =  MAT2D(k,1,b,ndf);
/*... retrosubstituicao*/
          MAT2D(k,1,gradU,ndf) = b1/lLsR[2];
          MAT2D(k,0,gradU,ndf) = 
                           (b0  -lLsR[1]*MAT2D(k,1,gradU,ndf))/lLsR[0];
        }
        else if(ndm == 3){
          b0 =  MAT2D(k,0,b,ndf);
          b1 =  MAT2D(k,1,b,ndf);
          b2 =  MAT2D(k,2,b,ndf);
/*... retrosubstituicao*/
          MAT2D(k,2,gradU,ndf) = b2/lLsR[5];
          MAT2D(k,1,gradU,ndf)= (b1 
                   - lLsR[4]* MAT2D(k,2,gradU,ndf))/lLsR[3];
          MAT2D(k,0,gradU,ndf) = (b0 
                   - lLsR[1]* MAT2D(k,1,gradU,ndf) 
                   - lLsR[2]* MAT2D(k,2,gradU,ndf))/lLsR[0];
        }
/*...................................................................*/
      }
/*...................................................................*/
    } 
/*...................................................................*/
  }
/*...................................................................*/
} 
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * LEASTSQUAREMATRIX : calcula a matriz de minimos quadrados para a  * 
 * rescontrucao de gradiente                                         *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * lKsi      -> vetores que unem centroide da celula central aos     *
 * mKsi      -> modulo do vetor ksi                                  * 
 * lSquare   -> nao definido                                         * 
 * lSquareR  -> nao definido                                         * 
 * type      -> tecnica de resolucao                                 * 
 * lnFace    -> numero vizinhos da celula                            * 
 * ndm       -> numero de dimensoes                                  * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * lSquare   -> matriz para a reconstrucao least Square              * 
 * lSquareR  -> fatoracao R (RCLSQUAREQR)                            * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 void leastSquareMatrix(DOUBLE *RESTRICT lKsi    ,DOUBLE *RESTRICT lmKsi
                  ,DOUBLE *RESTRICT lLsquare,DOUBLE *RESTRICT lLsquareR
                  ,short const type        
                  ,short const lnFace       ,short const ndm){

  DOUBLE dx[MAX_NUM_FACE*MAX_NDM],w[MAX_NUM_FACE*MAX_NDM];
  DOUBLE aTwA[2*MAX_NDM],inv[2*MAX_NDM];
  DOUBLE r[MAX_NDM*MAX_NDM];
  DOUBLE aux1,detA;
  short i,j,k,nf;
 
  if(ndm == 2){
    aTwA[0] = 0.e0;
    aTwA[1] = 0.e0;
    aTwA[2] = 0.e0;
  }
  else{
    aTwA[0] = 0.e0;
    aTwA[1] = 0.e0;
    aTwA[2] = 0.e0;
    aTwA[3] = 0.e0;
    aTwA[4] = 0.e0;
    aTwA[5] = 0.e0;
  }
  
  for(i=0;i<2*MAX_NDM;i++){
    inv[i] = 0.e0;
    r[i]   = 0.e0;
  }
   

/*...G=inv(aTwA)wAt */
  if(type == RCLSQUARE){
    for(nf=0;nf<lnFace;nf++){
      w[nf]  = 1.e0/lmKsi[nf];
      w[nf] *= w[nf];
/*... dimensao 2*/
      if(ndm == 2){
        MAT2D(nf,0,dx,ndm) = MAT2D(nf,0,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,1,dx,ndm) = MAT2D(nf,1,lKsi,ndm)*lmKsi[nf];
/*atTwa(0,0)*/
        aTwA[0] +=  w[nf]*MAT2D(nf,1,dx,ndm)*MAT2D(nf,1,dx,ndm);
/*atTwa(0,1)*/
        aTwA[1] +=  w[nf]*MAT2D(nf,0,dx,ndm)*MAT2D(nf,1,dx,ndm);
/*atTwa(1,1)*/
        aTwA[2] +=  w[nf]*MAT2D(nf,0,dx,ndm)*MAT2D(nf,0,dx,ndm);
      }
/*...................................................................*/

/*... dimensao 3*/
      else if(ndm == 3){
        MAT2D(nf,0,dx,ndm) = MAT2D(nf,0,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,1,dx,ndm) = MAT2D(nf,1,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,2,dx,ndm) = MAT2D(nf,2,lKsi,ndm)*lmKsi[nf];
/*atTwa(0,0)*/
        aTwA[0] +=  w[nf]*MAT2D(nf,0,dx,ndm)*MAT2D(nf,0,dx,ndm);
/*atTwa(0,1)*/
        aTwA[1] +=  w[nf]*MAT2D(nf,0,dx,ndm)*MAT2D(nf,1,dx,ndm);
/*atTwa(0,2)*/
        aTwA[2] +=  w[nf]*MAT2D(nf,0,dx,ndm)*MAT2D(nf,2,dx,ndm);
/*atTwa(1,1)*/
        aTwA[3] +=  w[nf]*MAT2D(nf,1,dx,ndm)*MAT2D(nf,1,dx,ndm);
/*atTwa(1,2)*/
        aTwA[4] +=  w[nf]*MAT2D(nf,1,dx,ndm)*MAT2D(nf,2,dx,ndm);
/*atTwa(2,2)*/
        aTwA[5] +=  w[nf]*MAT2D(nf,2,dx,ndm)*MAT2D(nf,2,dx,ndm);
      }
/*...................................................................*/
    }

/*... dimensao 2*/
    if(ndm == 2){
      detA   = 1.e0/(aTwA[2]*aTwA[0] - aTwA[1]*aTwA[1]);
      inv[0] =  aTwA[2]*detA;
      inv[1] = -aTwA[1]*detA;
      inv[2] =  aTwA[0]*detA;
    }
/*...................................................................*/
 
/*... dimensao 3*/
    else if(ndm == 3){
      detA = aTwA[0]*aTwA[3]*aTwA[5] + 2.e0*aTwA[1]*aTwA[4]*aTwA[2] 
             -( aTwA[2]*aTwA[2]*aTwA[3] 
              + aTwA[4]*aTwA[4]*aTwA[0]    
              + aTwA[1]*aTwA[1]*aTwA[5]); 
      detA = 1.e0/detA;
/*... inversao com cofatores*/
      inv[0] = (aTwA[3]*aTwA[5] - aTwA[4]*aTwA[4])*detA;
      inv[1] =-(aTwA[1]*aTwA[5] - aTwA[2]*aTwA[4])*detA;
      inv[2] = (aTwA[1]*aTwA[4] - aTwA[2]*aTwA[3])*detA;

      inv[3] = (aTwA[0]*aTwA[5] - aTwA[2]*aTwA[2])*detA;
      inv[4] =-(aTwA[0]*aTwA[4] - aTwA[1]*aTwA[2])*detA;

      inv[5] = (aTwA[0]*aTwA[3] - aTwA[1]*aTwA[1])*detA;
    }

    for(nf=0;nf<lnFace;nf++){
/*... dimensao 2*/
      if(ndm == 2){
        aux1 = w[nf];
        MAT2D(0,nf,lLsquare,lnFace) = aux1*(MAT2D(nf,0,dx,ndm)*inv[0] 
                                          + MAT2D(nf,1,dx,ndm)*inv[1]);
        MAT2D(1,nf,lLsquare,lnFace) = aux1*(MAT2D(nf,0,dx,ndm)*inv[1] 
                                          + MAT2D(nf,1,dx,ndm)*inv[2]);
      }
/*...................................................................*/

/*... dimensao 3*/

      else if(ndm == 3){
        aux1 = w[nf];
        MAT2D(0,nf,lLsquare,lnFace) = aux1*(MAT2D(nf,0,dx,ndm)*inv[0] 
                                        + MAT2D(nf,1,dx,ndm)*inv[1]
                                        + MAT2D(nf,2,dx,ndm)*inv[2]);
        MAT2D(1,nf,lLsquare,lnFace) = aux1*(MAT2D(nf,0,dx,ndm)*inv[1] 
                                        + MAT2D(nf,1,dx,ndm)*inv[3]
                                        + MAT2D(nf,2,dx,ndm)*inv[4]);
        MAT2D(2,nf,lLsquare,lnFace) = aux1*(MAT2D(nf,0,dx,ndm)*inv[2] 
                                        + MAT2D(nf,1,dx,ndm)*inv[4]
                                        + MAT2D(nf,2,dx,ndm)*inv[5]);
      }
/*...................................................................*/
    }
  }
/*...................................................................*/

/*... fatoracao QR-MGS*/
  else if(type == RCLSQUAREQR){
    for(nf=0;nf<lnFace;nf++){
      w[nf]  = 1.e0/lmKsi[nf];
      w[nf] *= w[nf];
/*... dimensao 2*/
      if(ndm == 2){
        MAT2D(nf,0,dx,ndm) = MAT2D(nf,0,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,1,dx,ndm) = MAT2D(nf,1,lKsi,ndm)*lmKsi[nf];
      }
/*...................................................................*/

/*... dimensao 3*/
      else if(ndm == 3){
        MAT2D(nf,0,dx,ndm) = w[nf]*MAT2D(nf,0,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,1,dx,ndm) = w[nf]*MAT2D(nf,1,lKsi,ndm)*lmKsi[nf];
        MAT2D(nf,2,dx,ndm) = w[nf]*MAT2D(nf,2,lKsi,ndm)*lmKsi[nf];
      }
/*...................................................................*/
    }

/*... MGS*/
    for(k=0;k<ndm;k++){
      aux1 = 0.e0;
      for(i=0;i<lnFace;i++)
        aux1 += MAT2D(i,k,dx,ndm)* MAT2D(i,k,dx,ndm);
      
      MAT2D(k,k,r,ndm) = sqrt(aux1);
      
      for(i=0;i<lnFace;i++)
        MAT2D(i,k,dx,ndm) /= MAT2D(k,k,r,ndm); 
      
      for(j=k+1;j<ndm;j++){
        aux1 = 0.e0;
        for(i=0;i<lnFace;i++)
          aux1 += MAT2D(i,k,dx,ndm)*MAT2D(i,j,dx,ndm);
        
        MAT2D(k,j,r,ndm) = aux1; 
        
        for(i=0;i<lnFace;i++)
          MAT2D(i,j,dx,ndm) -=  MAT2D(i,k,dx,ndm)*MAT2D(k,j,r,ndm); 
      }
    }
/*...................................................................*/

/*... dimensao 2*/
    if(ndm == 2){
        lLsquareR[0] = MAT2D(0,0,r,ndm);
        lLsquareR[1] = MAT2D(0,1,r,ndm);
        lLsquareR[2] = MAT2D(1,1,r,ndm);
    }
/*...................................................................*/

/*... R dimensao 3*/
    else if(ndm == 3){
        lLsquareR[0] = MAT2D(0,0,r,ndm);
        lLsquareR[1] = MAT2D(0,1,r,ndm);
        lLsquareR[2] = MAT2D(0,2,r,ndm);
        lLsquareR[3] = MAT2D(1,1,r,ndm);
        lLsquareR[4] = MAT2D(1,2,r,ndm);
        lLsquareR[5] = MAT2D(2,2,r,ndm);
    }
/*...................................................................*/

    for(nf=0;nf<lnFace;nf++){
/*... dimensao 2*/
      if(ndm == 2){
        MAT2D(0,nf,lLsquare,lnFace) = MAT2D(nf,0,dx,ndm);
        MAT2D(1,nf,lLsquare,lnFace) = MAT2D(nf,1,dx,ndm);
      }
/*...................................................................*/

/*... Qt - dimensao 3*/
      else if(ndm == 3){
        MAT2D(0,nf,lLsquare,lnFace) = w[nf]*MAT2D(nf,0,dx,ndm);
        MAT2D(1,nf,lLsquare,lnFace) = w[nf]*MAT2D(nf,1,dx,ndm);
        MAT2D(2,nf,lLsquare,lnFace) = w[nf]*MAT2D(nf,2,dx,ndm);
      }
/*...................................................................*/
    }
/*...................................................................*/
  } 
/*...................................................................*/
   
  }
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 28/04/2018                                   * 
 *-------------------------------------------------------------------* 
 * SN: numera dos nos das faces das celulas                          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * s   -> indefinido                                                 * 
 * ty  -> numero do tipo geometrico da celula                        * 
 * nel -> numero da celula                                           * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * s   -> numeracao dos nos por aresta                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
short sn(short *s,short ty,INT nel)
{
    
  long aux;
  short i;
  short isNodTria[6] = {0,1
                       ,1,2
                       ,2,0};
  short isNodQuad[8] = {0,1
                       ,1,2
                       ,2,3
                       ,3,0};
  short isNodTetr[12] = {1,2,3
                        ,0,3,2
                        ,0,1,3
                        ,0,2,1};

/*... face 1 - 1 4 3 2
      face 2 - 5 6 7 8
      face 3 - 1 2 6 5
      face 4 - 2 3 7 6
      face 5 - 2 3 8 7
      face 6 - 4 1 4 7   */
  short isNodHex[24] = {0,3,2,1
                       ,4,5,6,7
                       ,0,1,5,4
                       ,1,2,6,5
                       ,2,3,7,6
                       ,3,0,4,7};

  short isNodPira[20] = { 0,3,2, 1
                         ,0,1,4,-1
                         ,1,2,4,-1
                         ,2,3,4,-1
                         ,3,0,4,-1};

  switch(ty){

/* ... triangulo*/  
    case TRIACELL:
      for(i=0;i< 6;i++)
        s[i] = isNodTria[i];
      return 2;
    break;
/*...................................................................*/

/* ... quadrilateros*/  
    case QUADCELL:
      for(i=0;i< 8;i++)
        s[i] = isNodQuad[i];
      return 2;
    break;
/*...................................................................*/

/* ... tetraedros*/  
    case TETRCELL:
      for(i=0;i< 12;i++)
        s[i] = isNodTetr[i];
      return 3;
    break;
/*...................................................................*/

/* ... hexaedros*/  
    case HEXACELL:
      for(i=0;i<24;i++)
        s[i] = isNodHex[i];
      return 4;
    break;
/*...................................................................*/

/* ... hexaedros*/
    case PIRACELL:
      for (i = 0; i<20; i++)
        s[i] = isNodPira[i];
      return 4;
    break;
/*...................................................................*/

/*...*/    
    default: 
      aux = (long) nel;
      printf("Erro: tipo de celula nao exitentes. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
               ,aux,__FILE__,__func__);
      exit(EXIT_FAILURE);
    break; 
/*...................................................................*/
  }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * VOLUMECELL : volume(area) da celula                               * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * eta -> vetores paralelos as faces da celula                       * 
 * ty  -> tipo geometrico da celula                                  * 
 * ndm -> numero de dimensao                                         * 
 * nel -> numero da celula                                           * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * retorna o valores da volume (area) da celula                      * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE areaCell(DOUBLE *eta,short ty,short ndm,INT nel)
{
  
  long aux;

  switch(ty){

/* ... triangulo*/  
    case TRIACELL: 
      return areaTriaCell(eta,ndm);
    break;
/*...................................................................*/

/* ... quadrilateros*/  
    case QUADCELL:
      return areaQuadCell(eta,ndm);
    break;
/*...................................................................*/

/*...*/    
    default: 
      aux = (long) nel;
			printf("Erro: tipo de celula nao exitentes. Nel = %ld.\n"
				"Arquivo fonte:  \"%s\".\n"
				"Nome da funcao: \"%s\".\n"
				, aux, __FILE__, __func__);
      exit(EXIT_FAILURE);
    break; 
/*...................................................................*/
  }
}
/*********************************************************************/ 


/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * AREATRIACELL: calculo da area da celula triangular                * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * eta -> vetores paralelos as arestas da celula                     * 
 * ndm -> numero de dimensao                                         * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * retorna o valores da area da celula                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE areaTriaCell(DOUBLE *RESTRICT eta, short ndm)
{
  double v1[3],v2[3],v3[3],a,dot;

/*... aresta 1*/
  v1[0] = MAT2D(0,0,eta,ndm);
  v1[1] = MAT2D(0,1,eta,ndm);
  if( ndm == 2)
    v1[2] = 0.0e0;               
  else
    v1[2] = MAT2D(0,2,eta,ndm);
  

/*... aresta 2*/
  v2[0] = MAT2D(1,0,eta,ndm);
  v2[1] = MAT2D(1,1,eta,ndm);
  if( ndm == 2)
    v2[2] = 0.0e0;               
  else
    v2[2] = MAT2D(1,2,eta,ndm);

/*...*/
  prodVet(v1,v2,v3);
  dot = v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2];
  a = 0.5e0*sqrt(dot); 
  return a;
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * AREAQUADCELL: calculo da area da celula quadrangular              * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * eta -> vetores paralelos as arestas da celula                     * 
 * ndm -> numero de dimensao                                         * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * retorna o valores da area da celula                               * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE areaQuadCell(DOUBLE *RESTRICT eta,short ndm)
{

  double v1[3],v2[3],c[3];
  double a,dot;

/*...*/
  v1[0] = MAT2D(0,0,eta,ndm);
  v1[1] = MAT2D(0,1,eta,ndm);
  if(ndm == 2) 
    v1[2] = 0.0e0;
  else
    v1[2] = MAT2D(0,2,eta,ndm);
    
  v2[0] = MAT2D(1,0,eta,ndm);
  v2[1] = MAT2D(1,1,eta,ndm);
  if(ndm == 2) 
    v2[2] = 0.0e0;
  else
    v2[2] = MAT2D(1,2,eta,ndm);
  
  prodVet(v1,v2,c);
  dot = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
  a = 0.5e0*sqrt(dot);
/*...................................................................*/

/*...*/
  v1[0] = MAT2D(2,0,eta,ndm);
  v1[1] = MAT2D(2,1,eta,ndm);
  if(ndm == 2) 
    v1[2] = 0.0e0;
  else
    v1[2] = MAT2D(2,2,eta,ndm);
  

  v2[0] = MAT2D(3,0,eta,ndm);
  v2[1] = MAT2D(3,1,eta,ndm);
  if(ndm == 2) 
    v2[2] = 0.0e0;
  else
    v2[2] = MAT2D(3,2,eta,ndm);
  
  prodVet(v1,v2,c);
  dot = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
  a += 0.5e0*sqrt(dot);
/*...................................................................*/


/*...*/
  v1[0] = MAT2D(0,0,eta,ndm);
  v1[1] = MAT2D(0,1,eta,ndm);
  if(ndm == 2) 
    v1[2] = 0.0e0;
  else
    v1[2] = MAT2D(0,2,eta,ndm);
  
  v2[0] = MAT2D(3,0,eta,ndm);
  v2[1] = MAT2D(3,1,eta,ndm);
  if(ndm == 2) 
    v2[2] = 0.0e0;
  else
    v2[2] = MAT2D(3,2,eta,ndm);
  
  prodVet(v1,v2,c);
  dot = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
  a += 0.5e0*sqrt(dot);
/*...................................................................*/


/*...*/
  v1[0] = MAT2D(1,0,eta,ndm);
  v1[1] = MAT2D(1,1,eta,ndm);
  if(ndm == 2) 
    v1[2] = 0.0e0;
  else
    v1[2] = MAT2D(1,2,eta,ndm);
  
  v2[0] = MAT2D(2,0,eta,ndm);
  v2[1] = MAT2D(2,1,eta,ndm);
  if(ndm == 2) 
    v2[2] = 0.0e0;
  else
    v2[2] = MAT2D(2,2,eta,ndm);
  
  prodVet(v1,v2,c);
  dot = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
  a += 0.5e0*sqrt(dot);
/*...................................................................*/
  return 0.5e0*a;
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * VECTORKM2D : calculo da distancia entre a intersecao entre a reta *
 * que une o centroide da celulas e a face compartilhada             *  
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x      -> coordenadas dos nos da celula central e seus viznhos    * 
 * xc     -> centroide das celulas                                   * 
 * xm     -> pontos medios das faces da celula cenral                * 
 * vSkew  -> indefinido                                              * 
 * mvSkew -> indefinido                                              * 
 * sn     -> numeracao dos nos por aresta                            * 
 * nFace  -> numero de face da celula                                * 
 * sn     -> numeracao dos nos por aresta                            * 
 * maxNo  -> numero de nos por celula maximo da malha                * 
 * maxViz -> numero vizinhos por celula maximo da malha              * 
 * ndm    -> numero de dimensoes                                     * 
 * nel    -> numero da celula                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * vSkew  -> vetor entre o ponto medio a intersecao que une os       * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew -> distacia entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void vectorKm2d(DOUBLE *RESTRICT x     ,DOUBLE *RESTRICT xc
               ,DOUBLE *RESTRICT xm
               ,DOUBLE *RESTRICT vSkew ,DOUBLE *RESTRICT mvSkew
               ,short  *RESTRICT sn    ,short const nFace
               ,short const maxViz     ,short const maxNo       
               ,short const ndm        ,INT const nel)
{
  short i,cCell=maxViz;
  double a1,b1,c1,a2,b2,c2;
  double a[2][2],f[2],xi[2],lvSkew[2];
  short no1,no2;
  double det = 0.0;
  long aux;
/*... reta 1: a1x + b1y + c1 = 0
      reta 2: a2x + b2y + c2 = 0

      | a1 b1 | |x| = |-c1| 
      | a2 b2 | |y| = |-c2| 
    
  a1 = i e b1 = cCell (celula central)
  a2 e b2 pontos da aresta da celula central 
.......................................................................*/
  for(i=0;i<nFace;i++){
/*...ya - yb*/
    a1 = MAT2D(i,1,xc,ndm)     - MAT2D(cCell,1,xc,ndm);  
/*...xb - xa*/
    b1 = MAT2D(cCell,0,xc,ndm) - MAT2D(i,0,xc,ndm);  
/*...xayb - xbya*/
    c1 = MAT2D(i    ,0,xc,ndm)*MAT2D(cCell,1,xc,ndm) 
       - MAT2D(cCell,0,xc,ndm)*MAT2D(i    ,1,xc,ndm);
/*.....................................................................*/

    no1 = MAT2D(i,0,sn,ndm);
    no2 = MAT2D(i,1,sn,ndm);

/*...ya - yb*/
    a2 = MAT3D(cCell,no1,1,x,maxNo,ndm) - MAT3D(cCell,no2,1,x,maxNo,ndm);  
/*...xb - xa*/
    b2 = MAT3D(cCell,no2,0,x,maxNo,ndm) - MAT3D(cCell,no1,0,x,maxNo,ndm);  
/*...xayb - xbya*/
    c2 = MAT3D(cCell,no1,0,x,maxNo,ndm)*MAT3D(cCell,no2,1,x,maxNo,ndm)  
       - MAT3D(cCell,no2,0,x,maxNo,ndm)*MAT3D(cCell,no1,1,x,maxNo,ndm);  
/*.....................................................................*/

/*... verifica se ha intersecao entre duas retas*/
    det = a1*b2 - a2*b1;
    if(fabs(det) < sqrt(a1*a1+a2*a2)*1.e-16){
      aux = (long) nel;
      printf("Erro: As retas nao sao concorrentes. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
               "Determinante  : \"%e\".\n"
               ,aux,__FILE__,__func__,fabs(det));
      exit(EXIT_FAILURE);
    }
/*... matriz inversa A*/
    a[0][0] = (1.0/det)*b2;
    a[0][1] =-(1.0/det)*b1;
    a[1][0] =-(1.0/det)*a2;
    a[1][1] = (1.0/det)*a1;
/*... f*/
    f[0] = -c1;
    f[1] = -c2;
/*... x = A-1f*/
    xi[0] = a[0][0]*f[0] + a[0][1]*f[1];
    xi[1] = a[1][0]*f[0] + a[1][1]*f[1];
/*...vertor vSkew*/
    lvSkew[0] = xi[0] - MAT2D(i,0,xm,ndm);
    lvSkew[1] = xi[1] - MAT2D(i,1,xm,ndm);
/*...................................................................*/
    mvSkew[i] = sqrt(lvSkew[0]*lvSkew[0] + lvSkew[1]*lvSkew[1]);
/*...*/
    MAT2D(i,0,vSkew,ndm) = lvSkew[0];
    MAT2D(i,1,vSkew,ndm) = lvSkew[1];
  }
/*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * VECTORKM3D : calculo da distancia entre a intersecao entre a reta *
 * que une o centroide da celulas e a face compartilhada             *  
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x      -> coordenadas dos nos da celula central e seus viznhos    * 
 * xc     -> centroide das celulas                                   * 
 * xm     -> pontos medios das faces da celula cenral                * 
 * normal -> vetores normais as faces das celulas                    * 
 * ksi    -> vetores que unem centroide da celula central aos        *
 *           vizinhos destas                                         * 
 * vSkew  -> indefinido                                              * 
 * mvSkew -> indefinido                                              * 
 * nFace  -> numero de faces da celula                               * 
 * ndm    -> numero de dimensoes                                     * 
 * maxViz -> numero vizinhos por celula maximo da malha              * 
 * nel    -> numero da celula                                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * vSkew     -> vetor entre o ponto medio a intersecao que une os    * 
 *            centrois compartilhado nessa face da celula central    * 
 * mvSkew    -> distacia entre o ponto medio a intersecao que une os * 
 *            centrois compartilhado nessa face da celula central    * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void vectorKm3d(DOUBLE *RESTRICT xc   ,DOUBLE *RESTRICT xm
               ,DOUBLE *RESTRICT ksi  ,DOUBLE *RESTRICT normal
               ,DOUBLE *RESTRICT vSkew,DOUBLE *RESTRICT mvSkew
               ,short const nFace     ,short const ndm
               ,short const maxViz    ,INT nel)
{
  short i,cCell=maxViz;
  DOUBLE n[3],lKsi[3],xP[3],xR[3],xi[3],lvSkew[3],dot1,dot2,t;
  long aux;
/*... eq do plano :  nx*(x-xm) + ny*(y-ym) nz*(z-zm) = 0
      eq parametrica (x,y,z) = (xc,yc,zc) + t*(ksix,kisy,kisz)

  t = (nx(xm-xc) + ny(ym - yc) + nz(zm - zc))/(nx*ksix + ny*ksiy + nz*ksiz)

  t = n * dx / n*ksi  
  
.......................................................................*/
  xR[0] = MAT2D(cCell,0,xc,ndm);
  xR[1] = MAT2D(cCell,1,xc,ndm);
  xR[2] = MAT2D(cCell,2,xc,ndm);
  for(i=0;i<nFace;i++){
/*... vetor normal que define o plano e ponto*/
    n[0]  = MAT2D(i,0,normal,ndm);
    n[1]  = MAT2D(i,1,normal,ndm);
    n[2]  = MAT2D(i,2,normal,ndm);
    xP[0] = MAT2D(i,0,xm,ndm);
    xP[1] = MAT2D(i,1,xm,ndm);
    xP[2] = MAT2D(i,2,xm,ndm);
    
/*... vetor paralelo a reta*/
    lKsi[0] = MAT2D(i,0,ksi,ndm);
    lKsi[1] = MAT2D(i,1,ksi,ndm);
    lKsi[2] = MAT2D(i,2,ksi,ndm);
/*... */
    dot2 = n[0]*(xP[0]-xR[0]) +  n[1]*(xP[1]-xR[1]) + n[2]*(xP[2]-xR[2]);

/*... verifica se ha intersecao entre duas retas*/
    dot1 = n[0]*lKsi[0] + n[1]*lKsi[1] + n[2]*lKsi[2];
    if(fabs(dot1) < fabs(dot2)*1.e-16){
      aux = (long) nel;
      printf("Erro: A reta e o plano nao cruzam. Nel = %ld.\n"
               "Arquivo fonte:  \"%s\".\n"
               "Nome da funcao: \"%s\".\n"
               "Determinante  : \"%e\".\n"
               ,aux,__FILE__,__func__,fabs(dot1));
      exit(EXIT_FAILURE);
    }

/*... valor de t*/
    t = dot2/dot1;
/*...*/
    xi[0] = xR[0] + t*lKsi[0];
    xi[1] = xR[1] + t*lKsi[1];
    xi[2] = xR[2] + t*lKsi[2];

/*...vertor vSkew*/
    lvSkew[0] = xi[0] - MAT2D(i,0,xm,ndm);
    lvSkew[1] = xi[1] - MAT2D(i,1,xm,ndm);
    lvSkew[2] = xi[2] - MAT2D(i,2,xm,ndm);
/*...................................................................*/
    mvSkew[i] = 
     sqrt(lvSkew[0]*lvSkew[0]
        + lvSkew[1]*lvSkew[1] 
        + lvSkew[2]*lvSkew[2]);
/*...*/
    MAT2D(i,0,vSkew,ndm) = lvSkew[0];
    MAT2D(i,1,vSkew,ndm) = lvSkew[1];
    MAT2D(i,2,vSkew,ndm) = lvSkew[2];
  }
/*...................................................................*/
}
/*********************************************************************/
  
/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * VOLUME3DGREENGAUSS : volume calculado pelo teoreima do divergente *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * xm     -> pontos medios das faces da celula central               * 
 * normal -> vetor normal as faces                                   * 
 * fArea  -> area das faces                                          * 
 * nFace  -> numero de face da celula                                * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE volume3DGreenGauss(DOUBLE *RESTRICT xm,DOUBLE *RESTRICT normal
                         ,DOUBLE *RESTRICT fArea
                         ,short const nFace)
{
  
  DOUBLE dot = 0.e0;
  short i;
 
  for(i=0;i<nFace;i++){
    dot += fArea[i]
         * (MAT2D(i,0,xm,3)*MAT2D(i,0,normal,3) 
         +  MAT2D(i,1,xm,3)*MAT2D(i,1,normal,3) 
         +  MAT2D(i,2,xm,3)*MAT2D(i,2,normal,3)); 
  }   
  
  return dot*D1DIV3;

}
/*********************************************************************/
   
/********************************************************************* 
 * Data de criacao    : 02/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * volumeHexa: volume do hexaedro                                    *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x -> coordenadas dos vertices                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *                         | v3x v3y v3z |                           *
 * vol = v3 * ( v1 x v2) = | v1x v2y v2z |                           *
 *                         | v2x v2y v2z |                           *     
 *********************************************************************/
DOUBLE volumeHexa(DOUBLE *RESTRICT x)
{
  
  DOUBLE d[3][3],det;
 
/*... v3 - no5 - no1*/
  d[0][0] = MAT2D(4,0,x,3) - MAT2D(0,0,x,3); 
  d[0][1] = MAT2D(4,1,x,3) - MAT2D(0,1,x,3); 
  d[0][2] = MAT2D(4,2,x,3) - MAT2D(0,2,x,3);

/*... v1 - no2 - no1*/
  d[1][0] = MAT2D(1,0,x,3) - MAT2D(0,0,x,3); 
  d[1][1] = MAT2D(1,1,x,3) - MAT2D(0,1,x,3); 
  d[1][2] = MAT2D(1,2,x,3) - MAT2D(0,2,x,3); 

/*... v2 - no4 - no1*/
  d[2][0] = MAT2D(3,0,x,3) - MAT2D(0,0,x,3); 
  d[2][1] = MAT2D(3,1,x,3) - MAT2D(0,1,x,3); 
  d[2][2] = MAT2D(3,2,x,3) - MAT2D(0,2,x,3); 


  det  = d[0][0]*d[1][1]*d[2][2] + d[0][1]*d[1][2]*d[2][0] +
         d[0][2]*d[1][0]*d[2][1] - d[2][0]*d[1][1]*d[0][1] -
         d[0][1]*d[1][0]*d[2][2] - d[0][0]*d[2][1]*d[1][2];
  
  return det;

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 02/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * volumeTetra: volume do tetraedro                                  *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * x -> coordenadas dos vertices                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *             | x1 - x4   y1 - y4   z1 - z4 |                        *
 * vol = (1/6) | x2 - x4   y2 - y4   z2 - z4 |                        *
 *             | x3 - x4   z3 - y4   z3 - z4 |                        *        
 *********************************************************************/
DOUBLE volumeTetra(DOUBLE *RESTRICT x)
{
  
  DOUBLE d[3][3],det;
 
/*... v3 - no1 - no4*/
  d[0][0] = MAT2D(0,0,x,3) - MAT2D(3,0,x,3); 
  d[0][1] = MAT2D(0,1,x,3) - MAT2D(3,1,x,3); 
  d[0][2] = MAT2D(0,2,x,3) - MAT2D(3,2,x,3);

/*... v1 - no2 - no4*/
  d[1][0] = MAT2D(1,0,x,3) - MAT2D(3,0,x,3); 
  d[1][1] = MAT2D(1,1,x,3) - MAT2D(3,1,x,3); 
  d[1][2] = MAT2D(1,2,x,3) - MAT2D(3,2,x,3); 

/*... v2 - no3 - no1*/
  d[2][0] = MAT2D(2,0,x,3) - MAT2D(3,0,x,3); 
  d[2][1] = MAT2D(2,1,x,3) - MAT2D(3,1,x,3); 
  d[2][2] = MAT2D(2,2,x,3) - MAT2D(3,2,x,3); 


  det  = d[0][0]*d[1][1]*d[2][2] + d[0][1]*d[1][2]*d[2][0] +
         d[0][2]*d[1][0]*d[2][1] - d[2][0]*d[1][1]*d[0][1] -
         d[0][1]*d[1][0]*d[2][2] - d[0][0]*d[2][1]*d[1][2];
  printf("Volume tetra nao testado ainda!!!\n");
  exit(0);
  return det/6.e0;

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * LIMITFACEBASE : funcao limitadora de fluxo baseado na face        *
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * r       -> parametro utilizado                                    * 
 * iCod    -> tipo de funcao limitadora                              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE limitFaceBase(DOUBLE const r,short const iCod)
{
  #define NFUNCLIMTEB 5
  DOUBLE a,b,c;
  short i;
  char word[][WORD_SIZE]=
       {"VanLeer","VanAlbada"
       ,"MidMod ","Osher"   ,"SuperBee"};

  switch(iCod){

/*... Van Leer - TVD*/
    case VANLEERFACE: 
      a   = fabs(r);
      return (r + a ) / (1.0e0 + a);
    break;
/*...................................................................*/

/*... Van Albada  - TVD*/
    case VANALBADAFACE:
      if( r > 0.e0){ 
        a = r*r;
        return (r + a ) / (1.0e0 + a);
      }
      else 
        return 0.e0;
    break;
/*...................................................................*/

/*... Mid-Mod - TVD*/
    case MIDMODFACE:
      return max(0.e0,min(r,1.e0));
    break;
/*...................................................................*/

/*... OSHER - TVD*/
    case OSHERFACE:
      return  max(0.e0,min(r,2.e0));
    break;
/*...................................................................*/

/*... SUPERBEE - TVD*/
    case SUPERBEEFACE:
      a   = min(2.e0*r,1.e0);
      b   = min(r,2.e0);
      c   = max(a,b);
      return max(0.e0,c);
    break;
/*...................................................................*/

/*...*/
    default:
      printf("Erro: tipo de funcao limitadora de fluxo invalida.\n"
             "Arquivo fonte:  \"%s\".\n"
             "Nome da funcao: \"%s\".\n"
             "Linha         : \"%d\".\n"
             ,__FILE__,__func__,__LINE__);
      printf("Funcoes disponiveis:\n");
      for(i=0;i<NFUNCLIMTEB;i++)
        printf("%s\n",word[i]);
      exit(EXIT_FAILURE);
/*...................................................................*/
  }
  return 0.e0;
}
/*********************************************************************/ 

/*********************************************************************
* Data de criacao    : 22/08/2016                                   *
* Data de modificaco : 19/11/2017                                   *
*-------------------------------------------------------------------*
* NVD :                                                             *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* phiTil  -> parametro utilizado                                    *
* iCod    -> tipo de funcao limitadora                              *
* p       -> nao definido                                           *  
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* p[0] -> alfa                                                      *
* p[1] -> beta                                                      *
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
* phiFtil = alfa*phiCtil + beta                                     *
*********************************************************************/
void nvd(DOUBLE const phiTil,DOUBLE *p, short const iCod)
{

  short i; 

  char scheme_nvd[][20] ={"BCD"        ,"MUSCL"   ,"Smart"
                         ,"ModSmart"   ,"SuperBee","ModSuperBee"
                         ,"Stoic"      ,"MinMod"  ,"ModBCD"};

  p[1] = p[0] = 0.0e0;

  switch (iCod) {

/*... Bounded CD - NVD*/
    case BCD_NVD:
/*...*/
      if (phiTil >= 0 && phiTil <= 1.e0) {
        p[0] = 0.5e0;
        p[1] = 0.5e0;
      }
      else {
        p[0] = 1.0e0;
        p[1] = 0.0e0;
      }
      
    break;
/*..................................................................*/

/*... MUSCL - NVD*/
  case MUSCL_NVD:
/*...*/
//  if (phiTil >= 0 && phiTil <= 0.25e0)
//    a = 2.e0*phiTil;
//  else if (phiTil > 0.25e0 && phiTil <= 0.75e0)
//    a = phiTil + 2.5e-1;
//  else if (phiTil > 0.75e0 && phiTil <= 1.e0)
//    a = 1.e0;
//  else
//    a = phiTil;
//  return a;
  break;
/*...................................................................*/

/*... SMART - NVD*/
  case SMART_NVD:
/*...*/
      if (phiTil >= 0 && phiTil <= 5.0/6.0) {
        p[0] = 3.0/4.0;
        p[1] = 3.0/8.0;
      }
      else if (phiTil > 5.0/6.0 && phiTil <= 1.0) {
        p[0] = 0.0e0;
        p[1] = 1.0e0;
      }
      else {
        p[0] = 1.0e0;
        p[1] = 0.0e0;
      }      
    break;
/*..................................................................*/

/*... SMART - NVD*/
  case MSMART_NVD:
/*...*/
      if (phiTil >= 0 && phiTil <= 1.0/6.0) {
        p[0] = 3.0;
        p[1] = 0.0;
      }
      else if (phiTil > 1.0/6.0 && phiTil <= 0.7) {
        p[0] = 3.0e0/4.e0;
        p[1] = 3.0e0/8.0e0;
      }
      else if (phiTil > 0.7 && phiTil <= 1.0) {
        p[0] = D1DIV3;
        p[1] = D2DIV3;
      }
      else {
        p[0] = 1.0e0;
        p[1] = 0.0e0;
      }      
    break;
/*..................................................................*/

/*... SUPERBBE - NVD*/
    case SUPERBEE_NVD:
/*...*/
      if (phiTil >= 0 && phiTil <= 5.e-1) {
        p[0] = 0.5e0;
        p[1] = 0.5e0;
      }
      else if (phiTil > 5e-1 && phiTil <= D2DIV3) {
        p[0] = 1.5e0;
        p[1] = 0.0e0;
      }
      else if (phiTil > D2DIV3 && phiTil <= 1.e0) {
        p[0] = 0.0e0;
        p[1] = 1.0e0;
      }
      else {
        p[0] = 1.0e0;
        p[1] = 0.0e0;
      }      
    break;
/*..................................................................*/

/*... MSUPERBBE - NVD*/
    case MSUPERBEE_NVD:
/*...*/
      if (phiTil >= 0 && phiTil <= D1DIV3) {
        p[0] = 2.0e0;
        p[1] = 0.0e0;
      }
      else if (phiTil > D1DIV3 && phiTil <= 0.5e0) {
        p[0] = 0.5e0;
        p[1] = 0.5e0;
      }
      else if (phiTil > 0.5e0 && phiTil <= D2DIV3) {
        p[0] = 1.5e0;
        p[1] = 0.0e0;
      }
      else if (phiTil >  D2DIV3 && phiTil <= 1) {
        p[0] = 0.0e0;
        p[1] = 1.0e0;
      }
      else {
        p[0] = 1.0e0;
        p[1] = 0.0e0;
      }  
    break;
/*..................................................................*/

/*...STOIC - NVD*/
    case STOIC_NVD:
/*...*/
//  if (phiTil > 0.e0 && phiTil <= 5.e-1)
//    a = 5.e-1*(phiTil + 1.e0);
//  else if (phiTil > 5e-1 && phiTil <= D5DIV6)
//    a = 7.5e-1*phiTil + 3.75e-1;
//  else if (phiTil > D5DIV6 && phiTil <= 1.e0)
//    a = 1.e0;
//  else
//    a = phiTil;
//  return a;
    break;
/*..................................................................*/

/*...MINMOD - NVD*/
    case MINMOD_NVD:
/*...*/
      if (phiTil >= 0 && phiTil <= 5.e-1) {
        p[0] = 1.5e0;
        p[1] = 0.0e0;
      }
      else if (phiTil > 5e-1 && phiTil <= 1.0) {
        p[0] = 0.5e0;
        p[1] = 0.5e0;
      }
      else {
        p[0] = 1.0e0;
        p[1] = 0.0e0;
      }      
    break;
/*..................................................................*/

/*...MBCD - NVD*/
    case MBCD_NVD:
/*...*/
      if (phiTil >= 0.e0 && phiTil <= 1.0/4.0) {
        p[0] = 2.5;
        p[1] = 0.0e0;
      }
      else if (phiTil > 1.0/4.0 && phiTil <= 1.e0) {
        p[0] = 0.5e0;
        p[1] = 0.5e0;
      }
      else {
        p[0] = 1.0e0;
        p[1] = 0.0e0;
      }
    break;
/*..................................................................*/


    default:
      printf("Erro: tipo de funcao NVD invalida.\n"
             "Arquivo fonte:  \"%s\".\n"
             "Nome da funcao: \"%s\".\n"
             "Linha         : \"%d\".\n"
            , __FILE__, __func__, __LINE__);
      printf("Funcoes disponiveis:\n");
      for (i = 0; i<NFUNCNVD; i++)
        printf("%s\n", scheme_nvd[i]);
      exit(EXIT_FAILURE);
/*...................................................................*/
  }

}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * FACEBASETVD : metodo TVD baseado na face                          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE faceBaseTvd(short const nAresta    ,short const idCell
                 ,DOUBLE *RESTRICT u0
                 ,DOUBLE *RESTRICT gradUv,DOUBLE *RESTRICT gradUp
                 ,DOUBLE *RESTRICT lKsi  ,DOUBLE const lModKsi 
                 ,DOUBLE const cv
                 ,short const iCod       ,short const ndm)
{                    
  
  DOUBLE du,r,cvc=0.0e0,gf;
  DOUBLE eps=1.e-14;
 
  if( cv < 0.0e0) {
    du    = u0[idCell] - u0[nAresta] + eps;
    gf    = gradUv[0]*lKsi[0] + gradUv[1]*lKsi[1];
    if(ndm == 3) gf += gradUv[2]*lKsi[2];
    gf  = -2.e0*gf*lModKsi;
    r   = gf/du - 1.0e0;
    cvc = 0.5e0*limitFaceBase(r ,iCod)*du;
  }
  else{
    du    = u0[nAresta] - u0[idCell] + eps;
    gf    = gradUp[0]*lKsi[0] + gradUp[1]*lKsi[1];
    if(ndm == 3) gf += gradUp[2]*lKsi[2];
    gf  = 2.e0*gf*lModKsi;
    r   = gf/du - 1.0e0;
    cvc  = 0.5e0*limitFaceBase(r,iCod)*du;
  }
 
  return cvc;
} 
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 04/07/2016                                   *
 * Data de modificaco : 19/11/2017                                   * 
 *-------------------------------------------------------------------* 
 * FACEBASETVDV1 : metodo TVD baseado na face                        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * uC     -> valor central                                           * 
 * uV     -> valor do vizinho                                        * 
 * gradUc -> gradiente central                                       * 
 * gradUv -> gradiente do vizinho                                    * 
 * lKsi   -> vetor unitario entre os centroides                      * 
 * lModKsi-> distancia entre os centroides                           * 
 * cv     -> fluxo de massa                                          * 
 * iCod   -> tecnica TVD                                             * 
 * ndm    -> numero de dimensoes                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 * Tface = (TfaceUp)Imp + limte*(TfaceDown - TfaceUp)Exp             *   
 *********************************************************************/
DOUBLE faceBaseTvdV1(DOUBLE const uC     ,DOUBLE const uV
                 ,DOUBLE *RESTRICT gradUc,DOUBLE *RESTRICT gradUv
                 ,DOUBLE *RESTRICT lKsi  ,DOUBLE const lModKsi 
                 ,DOUBLE const cv
                 ,short const iCod       ,short const ndm)
{                    
  
  DOUBLE du,r,cvc=0.0e0,gf;
  DOUBLE eps=1.e-14;
 
  if( cv < 0.0e0) {
    du  = uC - uV + eps;
    gf  = gradUv[0]*lKsi[0] + gradUv[1]*lKsi[1];
    if(ndm == 3) gf += gradUv[2]*lKsi[2];
    gf  =-2.e0*gf*lModKsi;
    r   = gf/du - 1.0e0;
    cvc = 0.5e0*limitFaceBase(r ,iCod)*du;
  }
  else{
    du    = uV - uC + eps;
    gf    = gradUc[0]*lKsi[0] + gradUc[1]*lKsi[1];
    if(ndm == 3) gf += gradUc[2]*lKsi[2];
    gf    = 2.e0*gf*lModKsi;
    r     = gf/du - 1.0e0;
    cvc   = 0.5e0*limitFaceBase(r,iCod)*du;
  }

  return cvc;
} 
/*********************************************************************/ 

/*********************************************************************
* Data de criacao    : 22/08/2016                                   *
* Data de modificaco : 19/11/2017                                   *
*-------------------------------------------------------------------*
* FACEBANVD :                                                       *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* uP     -> valor central                                           *
* uV     -> valor do vizinho                                        *
* gradUc -> gradiente central                                       *
* gradUv -> gradiente do vizinho                                    *
* lKsi   -> vetor unitario entre os centroides                      *
* lModKsi-> distancia entre os centroides                           *
* m      -> fluxo de massa                                          *
* iCod   -> tecnica TVD                                             *
* ndm    -> numero de dimensoes                                     *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
DOUBLE faceBaseNvd(DOUBLE const uP        ,DOUBLE const uV
                  ,DOUBLE *RESTRICT gradUp,DOUBLE *RESTRICT gradUv
                  ,DOUBLE *RESTRICT lKsi  ,DOUBLE const lModKsi
                  ,DOUBLE const m
                  ,short const iCod       ,short const ndm)
{

  DOUBLE phiCtil, du, cvc = 0.0e0, gf,par[2];
  DOUBLE eps = 1.e-32;

  
/*... phiC = uV; phiD = uP*/
  if (m < 0.0e0) {
    gf      = gradUv[0]*lKsi[0] + gradUv[1]*lKsi[1];
    if (ndm == 3) gf += gradUv[2]*lKsi[2];
    gf = -2.e0*gf*lModKsi + eps;
    du      = uV - uP;
    phiCtil = du/gf + 1.e0;
    nvd(phiCtil,par,iCod);
    cvc     = (par[0] - 1.0)*uV + par[1]*uP
            + (1.0-par[0]-par[1])*(uP - gf);
//  printf("%e %e %e %e\n",par[0],par[1],gf,phiCtil);
//  cvc     = uP + (phiFtil - 1.e0)*gf - uV;
  }
/*... phiC = uP; phiD = uV*/
  else {
    gf      = gradUp[0]*lKsi[0] + gradUp[1]*lKsi[1];
    if (ndm == 3) gf += gradUp[2]*lKsi[2];
    gf = 2.e0*gf*lModKsi + eps;
    du      = uP - uV;
    phiCtil = du/gf + 1.e0;
    nvd(phiCtil,par, iCod);
    cvc     = (par[0] - 1.0)*uP + par[1]*uV
            + (1.0-par[0]-par[1])*(uV - gf); 
//  printf("%e %e %e %e\n",par[0],par[1],gf,phiCtil);
//  cvc     = uV + (phiFtil - 1.e0)*gf - uP;
  }

  return cvc;
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 22/08/2016                                   *
* Data de modificaco : 19/11/2017                                   *
*-------------------------------------------------------------------*
* SETNVD :                                                          *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* word ->                                                           *
* iCod -> nao definido                                              *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* iCod -> codigo da tecnica do termo advectivo                      *
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
void setNvd(char *word, short *iCod)
{
  short i;
  char scheme_nvd[][20] ={"BCD"        ,"MUSCL"   ,"Smart"
                         ,"ModSmart"   ,"SuperBee","ModSuperBee"
                         ,"Stoic"      ,"MinMod"  ,"ModBCD"};

/*...*/
  if (!strcmp(word,scheme_nvd[0])) {
    *iCod = BCD_NVD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : BCD\n");
  }
/*...................................................................*/

/*...*/
  else if (!strcmp(word, scheme_nvd[1])) {
    *iCod = MUSCL_NVD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : Muscl\n");
  }
/*...................................................................*/

/*...*/
  else if (!strcmp(word, scheme_nvd[2])) {
    *iCod = SMART_NVD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : Smart\n");
  }
/*...................................................................*/

/*...*/
  else if (!strcmp(word, scheme_nvd[3])) {
    *iCod = MSMART_NVD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : ModSmart\n");
  }
/*...................................................................*/

/*...*/
  else if (!strcmp(word, scheme_nvd[4])) {
    *iCod = SUPERBEE_NVD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : SuperBee\n");
  }
/*...................................................................*/

/*...*/
  else if (!strcmp(word, scheme_nvd[5])) {
    *iCod = MSUPERBEE_NVD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : ModSuperBee\n");
  }
/*...................................................................*/

/*...*/
  else if (!strcmp(word, scheme_nvd[6])) {
    *iCod = STOIC_NVD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : Stoic\n");
  }
/*...................................................................*/

/*...*/
  else if (!strcmp(word, scheme_nvd[7])) {
    *iCod = MINMOD_NVD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : MinMod\n");
  }
/*...................................................................*/

/*...*/
  else if (!strcmp(word, scheme_nvd[8])) {
    *iCod = MBCD_NVD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : MBCD\n");
  }
/*...................................................................*/

/*...*/
  else {
    printf("Erro: tipo de funcao limitadora fluxo invalida.\n"
           "Arquivo fonte:  \"%s\".\n"
           "Nome da funcao: \"%s\".\n"
           "Linha         : \"%d\".\n"
           , __FILE__, __func__, __LINE__);
    printf("Funcoes disponiveis:\n");
    for (i = 0; i< NFUNCNVD; i++)
      printf("%s\n", scheme_nvd[i]);
    exit(EXIT_FAILURE);
  }
  /*...................................................................*/
}
/*********************************************************************/

/********************************************************************* 
 * Data de criacao    : 16/07/2016                                   *
 * Data de modificaco : 20/11/2017                                   * 
 *-------------------------------------------------------------------* 
 * UPWINDLINEARV1: upwind linear com resconstrucao linear            * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * uC     -> valor central                                           * 
 * uV     -> valor do vizinho                                        * 
 * gradUc -> gradiente central                                       * 
 * gradUv -> gradiente do vizinho                                    * 
 * rC     -> distancia ate o ponto central da face (Celula central)  * 
 * rV     -> distancia ate o ponto central da face (Celula vizinha)  * 
 * wfn    -> velociade normal a face                                 * 
 * ndm    -> numero de dimensoes                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              *
 * Tface = (TfaceUp)Imp + (TfaceUp + grad(T)*r - TFaceUp)Exp         *   
 *       = (TfaceUp)Imp + (grad(T)*r)Exp                             *   
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE upwindLinearV1(DOUBLE const uC     ,DOUBLE const uV
                 ,DOUBLE *RESTRICT gradUc ,DOUBLE *RESTRICT gradUv
                 ,DOUBLE *RESTRICT rC     ,DOUBLE *RESTRICT rV
                 ,DOUBLE const wfn        ,short const ndm)
{                    
  
  DOUBLE cvc=0.e0;
 
  if (wfn < 0.0e0) {
    cvc = gradUv[0]*rV[0] + gradUv[1]*rV[1];
    if(ndm == 3) cvc += gradUv[2]*rV[2];
  }      
  else {
    cvc = gradUc[0]*rC[0] + gradUc[1]*rC[1];
    if(ndm == 3) cvc += gradUc[2]*rC[2];
  }
 
  return cvc;
} 
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 16/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * DEFERREDCD : Diferencao central com correcao atrasada             * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * uC     -> valor central                                           * 
 * uV     -> valor do vizinho                                        * 
 * wfn    -> velocidade normal                                       * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------*
 * valor na face = (Upwind)implicito + (Central - Upwind)explicito   * 
 *********************************************************************/
DOUBLE deferredCd(DOUBLE const uC,DOUBLE const uV,DOUBLE const wfn)
{                    
  
  DOUBLE cvc=0.0e0;
 
  if(wfn < 0.0e0) {
    cvc = 0.5e0*(uC+uV) - uV;
  }
  else{
    cvc = 0.5e0*(uC+uV) - uC;
  }
 
  return cvc;
} 
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 20/11/2017                                   *
 * Data de modificaco : 05/02/2018                                   * 
 *-------------------------------------------------------------------* 
 * deferredLust : tecnica do openFOAM                                * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * uC     -> valor central                                           * 
 * uV     -> valor do vizinho                                        * 
 * gradUc -> gradiente central                                       * 
 * gradUv -> gradiente do vizinho                                    * 
 * rC     -> distancia ate o ponto central da face (Celula central)  * 
 * rV     -> distancia ate o ponto central da face (Celula vizinha)  * 
 * alphaMenosUm -> interpolacao linear na face (1.e0 - alfa)         *
 * alpha        -> interpolacao linear na face (alfa)                *
 * beta   -> para metro do lust ( beta 0 -1 )                        *
 *         0 - CD                                                    *
 *         1 - UPWIND                                                *
 * wfn    -> velociade normal a face                                 * 
 * ndm    -> numero de dimensoes                                     * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *                                                                   * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------*
 * valor na face = (Upwind)implicito + (Lustl - Upwind)explicito     *
 * Lust = beta*UpwindLinear + (1 - beta)*central                     *
 * UpwindLinear = Upwind + gradU*r 
 *********************************************************************/
DOUBLE deferredLust(DOUBLE const uC          ,DOUBLE const uV
                   ,DOUBLE *RESTRICT gradUc  ,DOUBLE *RESTRICT gradUv
                   ,DOUBLE *RESTRICT rC      ,DOUBLE *RESTRICT rV
                   ,DOUBLE const alphaMenosUm,DOUBLE const alpha
                   ,DOUBLE const beta 
                   ,DOUBLE const wfn         ,short const ndm)
{                    
  
  DOUBLE gf,cd,up,ul,uf1,uf2;

/*... upwind linear *SouUp*/   
  if (wfn < 0.0e0) {
    gf = gradUv[0]*rV[0] + gradUv[1]*rV[1];
    if(ndm == 3) gf += gradUv[2]*rV[2];
    up = uV;
  }      
  else {
    gf = gradUc[0]*rC[0] + gradUc[1]*rC[1];
    if(ndm == 3) gf += gradUc[2]*rC[2];
    up = uC;
  }  
  ul = up + gf;
/*....................................................................*/

/*... centrado*/
//cd = alphaMenosUm*uC+alpha*uV;
  uf1 = uC + gradUc[0]*rC[0] + gradUc[1]*rC[1];
  uf2 = uV + gradUv[0]*rV[0] + gradUv[1]*rV[1];
  if (ndm == 3) {
    uf1 += gradUc[2]*rC[2];
    uf2 += gradUv[2]*rV[2];
  }
  cd = alphaMenosUm*uf1+alpha*uf2;
/*....................................................................*/
 
  return  beta*ul + (1.e0 - beta)*cd - up ;

} 
/*********************************************************************/ 


/********************************************************************* 
 * Data de criacao    : 00/00/2015                                   *
 * Data de modificaco : 21/08/2016                                   * 
 *-------------------------------------------------------------------* 
 * SETTVD :                                                          * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * word ->                                                           * 
 * iCod -> nao definido                                              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * iCod -> codigo da tecnica do termo advectivo                      * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void setTvd(char *word,short *iCod)
{
  short i;
  char fBase[][WORD_SIZE]=
       {"VanLeer" ,"VanAlbada"
       ,"MidMod ","Osher"   ,"SuperBee"};
/*...*/
  if(!strcmp(word,"VanLeer")){
    *iCod = VANLEERFACE; 
    if(!mpiVar.myId ) fprintf(fileLogExc,"iCod  : VanLeer\n");
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,"VanAlbada")){
    *iCod =  VANALBADAFACE; 
    if(!mpiVar.myId ) fprintf(fileLogExc,"iCod  : VanAlbade\n");
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,"MidMod")){
    *iCod =  MIDMODFACE; 
    if(!mpiVar.myId ) fprintf(fileLogExc,"iCod  : MidMod\n");
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,"Osher")){
    *iCod =  OSHERFACE; 
    if(!mpiVar.myId ) fprintf(fileLogExc,"iCod  : Osher\n");
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,"SuperBee")){
    *iCod =  SUPERBEEFACE; 
    if(!mpiVar.myId )fprintf(fileLogExc,"iCod  : SuperBee\n");
  }
/*...................................................................*/

/*...*/
  else{
    printf("Erro: tipo de funcao limitadora fluxo invalida.\n"
           "Arquivo fonte:  \"%s\".\n"
           "Nome da funcao: \"%s\".\n"
           "Linha         : \"%d\".\n"
           ,__FILE__,__func__,__LINE__);
    printf("Funcoes disponiveis:\n");
    for(i=0;i<NFUNCLIMTFACE;i++)
      printf("%s\n",fBase[i]);
    exit(EXIT_FAILURE);
  }
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 08/08/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * SETDIFFUSIONSCHME;                                                * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * word ->                                                           * 
 * iCod -> nao definido                                              * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * iCod -> codigo da tecnica do termo difusivo                       * 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void  setDiffusionScheme(char *word,short *iCod)
{
  short i;
  char fDif[][WORD_SIZE]=
       {"Orthogonal"  ,"Minimal"
       ,"OrthogonalC" ,"OverRelaxed"};
/*...*/  
  if(!strcmp(word,fDif[0])){
    *iCod = ORTHOGONAL; 
    if(!mpiVar.myId ) fprintf(fileLogExc,"iCod  : %s\n",fDif[0]);
  }
/*...................................................................*/
  
/*...*/
  else if(!strcmp(word,fDif[1])){
    *iCod = MINIMAL; 
    if(!mpiVar.myId ) fprintf(fileLogExc,"iCod  : %s\n",fDif[1]);
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,fDif[2])){
    *iCod = ORTHOGONALC; 
    if(!mpiVar.myId ) fprintf(fileLogExc,"iCod  : %s\n",fDif[2]);
  }
/*...................................................................*/

/*...*/
  else if(!strcmp(word,fDif[3])){
    *iCod = OVERRELAXED; 
    if(!mpiVar.myId ) fprintf(fileLogExc,"iCod  : %s\n",fDif[3]);
  }
/*...................................................................*/

/*...*/
  else{
    printf("Erro: tipo de correcao nao ortogonal invalida.\n"
           "Arquivo fonte:  \"%s\".\n"
           "Nome da funcao: \"%s\".\n"
           "Linha         : \"%d\".\n"
           ,__FILE__,__func__,__LINE__);
    printf("Funcoes disponiveis:\n");
    for(i=0;i<4;i++)
      printf("%s\n",fDif[i]);
    exit(EXIT_FAILURE);
  }
/*...................................................................*/
}
/*********************************************************************/ 

/*********************************************************************
* Data de criacao    : 08/08/2016                                    *
* Data de modificaco : 20/11/2017                                    *
*------------------------------------------------------------------- *
* SETADVECTIONSCHEME :                                               *
*------------------------------------------------------------------- *
* Parametros de entrada:                                             *
*------------------------------------------------------------------- *
* word   ->                                                          *
* iCod   -> nao definido                                             *
* FileIn -> arquivo                                                  *
*------------------------------------------------------------------- *
* Parametros de saida:                                               *
*------------------------------------------------------------------- *
* adv -> codigo da tecnica do termo advectivo                        *
*------------------------------------------------------------------- *
* OBS:                                                               *
*------------------------------------------------------------------- *
*********************************************************************/
void  setAdvectionScheme(char *word, Advection *adv,FILE *fileIn)
{
  short i,nTerm=6;
  char fAdv[][WORD_SIZE] = { "FoUp","CD"   ,"SoUp"  /*0,1,2*/
                            ,"TVD" ,"NVD"  ,"LUST"};/*3,4,5*/
  /*...*/
  if (!strcmp(word, fAdv[0])) {
    adv->iCod1 = FOUP;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : %s\n", fAdv[0]);
  }
  /*...................................................................*/

  /*...*/
  else if (!strcmp(word, fAdv[1])) {
    adv->iCod1 = CD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : %s\n", fAdv[1]);
  }
  /*...................................................................*/

/*...*/
  else if (!strcmp(word, fAdv[2])) {
    adv->iCod1 = SOUP;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : %s\n", fAdv[2]);
  }
/*...................................................................*/

/*...*/
  else if (!strcmp(word, fAdv[3])) {
    adv->iCod1 = TVD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : %s\n", fAdv[3]);
    readMacro(fileIn, word, false);
    setTvd(word, &adv->iCod2);
  }
/*...................................................................*/

/*...*/
  else if (!strcmp(word, fAdv[4])) {
    adv->iCod1 = NVD;
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : %s\n", fAdv[4]);
    readMacro(fileIn, word, false);
    setNvd(word, &adv->iCod2);
  }
/*...................................................................*/

/*...*/
  else if (!strcmp(word, fAdv[5])) {
    adv->iCod1 = LUST;
    readMacro(fileIn, word, false);
    adv->par[0] = atof(word);
    if (!mpiVar.myId) fprintf(fileLogExc,"iCod  : %s(%lf)\n"
                                        , fAdv[5],adv->par[0]);
  }
/*...................................................................*/

/*...*/
  else {
    printf("Erro: tipo de tecnica advectiva invalida.\n"
           "Arquivo fonte:  \"%s\".\n"
           "Nome da funcao: \"%s\".\n"
           "Linha         : \"%d\".\n"
           , __FILE__, __func__, __LINE__);
    printf("Funcoes disponiveis:\n");
    for (i = 0; i<nTerm; i++)
      printf("%s\n", fAdv[i]);
    exit(EXIT_FAILURE);
  }
  /*...................................................................*/
}
/*********************************************************************/
      
/********************************************************************* 
 * Data de criacao    : 18/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * SIZECAR: tamanho caracteristico da celula                         * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------* 
 * volume -> area ou voluma                                          * 
 * ndm    -> dimensao                                                * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
DOUBLE sizeCar(DOUBLE const volume,short const ndm)
{

  if(ndm == 2)      return sqrt(volume);
  else if(ndm == 3) return pow(volume,D1DIV3);

  return 0.0;
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 22/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * GRADFACENULL :  garante o gradiente nulo da face usando uma extra * 
 * polacao no centroide                                              * 
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
void gradFaceNull(DOUBLE *RESTRICT gradVelFace
                 ,DOUBLE *RESTRICT gradVelCell
                 ,DOUBLE *RESTRICT xmcc       ,short const ndm)
{

  DOUBLE prod[3],tensor[3][3],eFace[3],mod;

/*...................................................................*/
  if(ndm == 2){
/*... vetor unitario na direcao do centroide ate o ponto medio da face
      eFace */
    mod      = sqrt(xmcc[0]*xmcc[0] + xmcc[1]*xmcc[1]);
    eFace[0] = xmcc[0]/mod;
    eFace[1] = xmcc[1]/mod;
/*...................................................................*/

/*... (Grad(VelP),eFace)*/
    prod[0] = MAT2D(0,0,gradVelCell,2)*eFace[0]
            + MAT2D(0,1,gradVelCell,2)*eFace[1];

    prod[1] = MAT2D(1,0,gradVelCell,2)*eFace[0]
            + MAT2D(1,1,gradVelCell,2)*eFace[1];
/*...................................................................*/

/*... produto tensorial entre dois vetores*/
    tensor[0][0] = prod[0]*eFace[0];
    tensor[0][1] = prod[0]*eFace[1];
    
    tensor[1][0] = prod[1]*eFace[0];
    tensor[1][1] = prod[1]*eFace[1];

/*...................................................................*/

/*... Grad(VelFace) = Grad(VelP) - (Grad(VelP),eFace)eFace*/
     MAT2D(0,0,gradVelFace,2) =  MAT2D(0,0,gradVelCell,2)
                              -  tensor[0][0];        

     MAT2D(0,1,gradVelFace,2) =  MAT2D(0,1,gradVelCell,2)
                              -  tensor[0][1];        

     MAT2D(1,0,gradVelFace,2) =  MAT2D(1,0,gradVelCell,2)
                              -  tensor[1][0];        

     MAT2D(1,1,gradVelFace,2) =  MAT2D(1,1,gradVelCell,2)
                              -  tensor[1][1];        
/*...................................................................*/

  }
/*...................................................................*/

/*...*/  
  else if(ndm == 3){
/*... vetor unitario na direcao do centroide ate o ponto medio da face
      eFace */
    mod      = sqrt(xmcc[0]*xmcc[0] 
                  + xmcc[1]*xmcc[1] 
                  + xmcc[2]*xmcc[2]);
    eFace[0] = xmcc[0]/mod;
    eFace[1] = xmcc[1]/mod;
    eFace[2] = xmcc[2]/mod;
/*...................................................................*/

/*... (Grad(VelP),eFace)*/
    prod[0] = MAT2D(0,0,gradVelCell,3)*eFace[0]
            + MAT2D(0,1,gradVelCell,3)*eFace[1] 
            + MAT2D(0,2,gradVelCell,3)*eFace[2];

    prod[1] = MAT2D(1,0,gradVelCell,3)*eFace[0]
            + MAT2D(1,1,gradVelCell,3)*eFace[1] 
            + MAT2D(1,2,gradVelCell,3)*eFace[2];
    
    prod[2] = MAT2D(2,0,gradVelCell,3)*eFace[0]
            + MAT2D(2,1,gradVelCell,3)*eFace[1] 
            + MAT2D(2,2,gradVelCell,3)*eFace[2];
/*...................................................................*/

/*... produto tensorial entre dois vetores*/
    tensor[0][0] = prod[0]*eFace[0];
    tensor[0][1] = prod[0]*eFace[1];
    tensor[0][2] = prod[0]*eFace[2];
    
    tensor[1][0] = prod[1]*eFace[0];
    tensor[1][1] = prod[1]*eFace[1];
    tensor[1][2] = prod[1]*eFace[2];
    
    tensor[2][0] = prod[2]*eFace[0];
    tensor[2][1] = prod[2]*eFace[1];
    tensor[2][2] = prod[2]*eFace[2];

/*...................................................................*/

/*... Grad(VelFace) = Grad(VelP) - (Grad(VelP),eFace)eFace*/
     MAT2D(0,0,gradVelFace,3) =  MAT2D(0,0,gradVelCell,3)
                              -  tensor[0][0];        
     MAT2D(0,1,gradVelFace,3) =  MAT2D(0,1,gradVelCell,3)
                              -  tensor[0][1];        
     MAT2D(0,2,gradVelFace,3) =  MAT2D(0,2,gradVelCell,3)
                              -  tensor[0][2];        

     MAT2D(1,0,gradVelFace,3) =  MAT2D(1,0,gradVelCell,3)
                              -  tensor[1][0];        
     MAT2D(1,1,gradVelFace,3) =  MAT2D(1,1,gradVelCell,3)
                              -  tensor[1][1];        
     MAT2D(1,2,gradVelFace,3) =  MAT2D(1,2,gradVelCell,3)
                              -  tensor[1][2];        
     
     MAT2D(2,0,gradVelFace,3) =  MAT2D(2,0,gradVelCell,3)
                              -  tensor[2][0];        
     MAT2D(2,1,gradVelFace,3) =  MAT2D(2,1,gradVelCell,3)
                              -  tensor[2][1];        
     MAT2D(2,2,gradVelFace,3) =  MAT2D(2,2,gradVelCell,3)
                              -  tensor[2][2];        
/*...................................................................*/

  }
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 08/07/2016                                   *
 * Data de modificaco : 00/00/0000                                   * 
 *-------------------------------------------------------------------* 
 * difusionScheme :  correcao nao ortogonal do termo difusivo        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * n       -> vetor normal unitarios                                 *
 * ksi     -> vetor que une os centroide central e o vizinho         *   
 * lArea   -> area da face                                           *   
 * modKsi  -> distancia entre o centroide central e o vizinho        *   
 * e       -> nao definido                                           *   
 * t       -> nao definido                                           *   
 * ndm     -> dimensao                                               *   
 * iCodDif -> codigo da tecnica nao ortogonal                        * 
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * e       -> parcela implicita                                      *   
 * t       -> parcela explicita                                      *   
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 * grad(phi)*S = grad(phi)*E + grad(phi)*T                           * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 void difusionScheme(DOUBLE *RESTRICT n ,DOUBLE *RESTRICT ksi 
                    ,DOUBLE const fArea ,DOUBLE const lModKsi  
                    ,DOUBLE *RESTRICT e ,DOUBLE *RESTRICT t   
                    ,short const ndm    ,short const iCod){

  DOUBLE nk,tmp;
  short i;
  char word[][WORD_SIZE]=
       {"Orthogonal"  ,"Minimal"
       ,"OrthogonalC" ,"OverRelaxed"};

  switch(iCod){

/*... ORTHOGONAL*/
    case ORTHOGONAL:
      if(ndm == 2){
/*... vetor E*/
        e[0] = fArea*ksi[0];        
        e[1] = fArea*ksi[1];
/*...................................................................*/

/*... vetor T*/
				t[0] = 0.e0;        
				t[1] = 0.e0;
/*...................................................................*/
      }        
/*...................................................................*/

/*...*/
      else if(ndm == 3){
/*... vetor E*/
        e[0] = fArea*ksi[0];        
        e[1] = fArea*ksi[1];
        e[2] = fArea*ksi[2];
/*...................................................................*/

/*... vetor T*/
				t[0] = 0.e0;        
				t[1] = 0.e0;        
				t[2] = 0.e0;
/*...................................................................*/
      }        
    break;
/*...................................................................*/

/*... Minimal*/
    case MINIMAL: 
      if(ndm == 2){
/*...*/
        nk = n[0]*ksi[0] + n[1]*ksi[1];
/*...................................................................*/

/*... vetor E*/
        tmp  = fArea*nk;
        e[0] = tmp*ksi[0];        
        e[1] = tmp*ksi[1];
/*...................................................................*/

/*... vetor T*/
        t[0] = fArea*n[0] - e[0];        
        t[1] = fArea*n[1] - e[1];        
/*...................................................................*/
      }        
/*...................................................................*/

/*...*/
      else if(ndm == 3){
/*...*/
        nk = n[0]*ksi[0] + n[1]*ksi[1] +  n[2]*ksi[2];
/*...................................................................*/

/*... vetor E*/
        tmp  = fArea*nk;
        e[0] = tmp*ksi[0];        
        e[1] = tmp*ksi[1];
        e[2] = tmp*ksi[2];
/*...................................................................*/

/*... vetor T*/
        t[0] = fArea*n[0] - e[0];        
        t[1] = fArea*n[1] - e[1];        
        t[2] = fArea*n[2] - e[2];        
/*...................................................................*/
      }        
    break;
/*...................................................................*/

/*... Orthogonal*/
    case ORTHOGONALC:
      if(ndm == 2){
/*... vetor E*/
        e[0] = fArea*ksi[0];        
        e[1] = fArea*ksi[1];
/*...................................................................*/

/*... vetor T*/
        t[0] = fArea*n[0] - e[0];        
        t[1] = fArea*n[1] - e[1];        
/*...................................................................*/
      }        
/*...................................................................*/

/*...*/
      else if(ndm == 3){
/*... vetor E*/
        e[0] = fArea*ksi[0];        
        e[1] = fArea*ksi[1];
        e[2] = fArea*ksi[2];
/*...................................................................*/

/*... vetor T*/
        t[0] = fArea*n[0] - e[0];        
        t[1] = fArea*n[1] - e[1];        
        t[2] = fArea*n[2] - e[2];        
/*...................................................................*/
      }        
    break;
/*...................................................................*/

/*... OverRelaxed*/
    case OVERRELAXED:
      if(ndm == 2){
/*...*/
        nk = n[0]*ksi[0] + n[1]*ksi[1];
        nk = 1.e0/nk;
/*...................................................................*/

/*... vetor E*/
        tmp  = fArea*nk;
        e[0] = tmp*ksi[0];        
        e[1] = tmp*ksi[1];        
/*...................................................................*/

/*... vetor T*/
        t[0] = fArea*n[0] - e[0];        
        t[1] = fArea*n[1] - e[1];        
/*...................................................................*/
      }
/*...................................................................*/

/*...*/
      else if(ndm == 3){
/*...*/
        nk = n[0]*ksi[0] + n[1]*ksi[1]  + n[2]*ksi[2];
        nk = 1.e0/nk;
/*...................................................................*/

/*... vetor E*/
        tmp  = fArea*nk;
        e[0] = tmp*ksi[0];        
        e[1] = tmp*ksi[1];        
        e[2] = tmp*ksi[2];        
/*...................................................................*/

/*... vetor T*/
        t[0] = fArea*n[0] - e[0];        
        t[1] = fArea*n[1] - e[1];        
        t[2] = fArea*n[2] - e[2];        
/*...................................................................*/
      }
    break;
/*...................................................................*/

/*...*/
    default:
      printf("Erro: tipo de correcao nao ortogonal invalida.\n"
             "Arquivo fonte:  \"%s\".\n"
             "Nome da funcao: \"%s\".\n"
             "Linha         : \"%d\".\n"
             ,__FILE__,__func__,__LINE__);
      printf("Funcoes disponiveis:\n");
      for(i=0;i<4;i++)
        printf("%s\n",word[i]);
      exit(EXIT_FAILURE);
/*...................................................................*/
  }
}
/*********************************************************************/ 

/********************************************************************* 
 * Data de criacao    : 08/07/2016                                   *
 * Data de modificaco : 30/01/2018                                   * 
 *-------------------------------------------------------------------* 
 * difusionScheme :  correcao nao ortogonal do termo difusivo        * 
 *-------------------------------------------------------------------* 
 * Parametros de entrada:                                            * 
 *-------------------------------------------------------------------*
 * s       -> vetor de area modificado                               *
 * ksi     -> vetor que une os centroide central e o vizinho         *
 * e       -> nao definido                                           *
 * t       -> nao definido                                           *
 * ndm     -> dimensao                                               *
 * iCodDif -> codigo da tecnica nao ortogonal                        *
 *-------------------------------------------------------------------* 
 * Parametros de saida:                                              * 
 *-------------------------------------------------------------------* 
 * e       -> parcela implicita                                      *   
 * t       -> parcela explicita                                      *   
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 * grad(phi)*S = grad(phi)*E + grad(phi)*T                           * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
 void difusionSchemeNew(DOUBLE *RESTRICT s , DOUBLE *RESTRICT ksi 
                      , DOUBLE *RESTRICT e , DOUBLE *RESTRICT t   
                      , short const ndm    , short const iCod){

  DOUBLE nk,tmp, ss;
  short i;
  char word[][WORD_SIZE]=
       {"Orthogonal"  ,"Minimal"
       ,"OrthogonalC" ,"OverRelaxed"};

  switch(iCod){

/*... ORTHOGONAL*/
    case ORTHOGONAL:
      if(ndm == 2){
/*...*/
        ss = sqrt(s[0]*s[0]+s[1]*s[1]);
/*...................................................................*/

/*... vetor E*/
        e[0] = ss*ksi[0];        
        e[1] = ss*ksi[1];
/*...................................................................*/

/*... vetor T*/
				t[0] = 0.e0;        
				t[1] = 0.e0;
/*...................................................................*/
      }        
/*...................................................................*/

/*...*/
      else if(ndm == 3){
/*...*/
        ss = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
/*...................................................................*/

/*... vetor E*/
        e[0] = ss*ksi[0];        
        e[1] = ss*ksi[1];
        e[2] = ss*ksi[2];
/*...................................................................*/

/*... vetor T*/
				t[0] = 0.e0;        
				t[1] = 0.e0;        
				t[2] = 0.e0;
/*...................................................................*/
      }        
    break;
/*...................................................................*/

/*... Minimal*/
    case MINIMAL: 
      if(ndm == 2){
 /*...*/
        ss = sqrt(s[0]*s[0]+s[1]*s[1]);
/*...................................................................*/

/*...*/
        nk = s[0]*ksi[0] + s[1]*ksi[1];
/*...................................................................*/

/*... vetor E*/
        tmp  = nk*ss;
        e[0] = tmp*ksi[0];        
        e[1] = tmp*ksi[1];
/*...................................................................*/

/*... vetor T*/
        t[0] = s[0] - e[0];        
        t[1] = s[1] - e[1];        
/*...................................................................*/
      }        
/*...................................................................*/

/*...*/
      else if(ndm == 3){
/*...*/
        ss = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
/*...................................................................*/

/*...*/
        nk = s[0]*ksi[0] + s[1]*ksi[1] + s[2]*ksi[2];
/*...................................................................*/

/*... vetor E*/
        tmp  = nk*ss;
        e[0] = tmp*ksi[0];        
        e[1] = tmp*ksi[1];
        e[2] = tmp*ksi[2];
/*...................................................................*/

/*... vetor T*/
        t[0] = s[0] - e[0];        
        t[1] = s[1] - e[1];        
        t[2] = s[2] - e[2];        
/*...................................................................*/
      }        
    break;
/*...................................................................*/

/*... Orthogonal*/
    case ORTHOGONALC:
      if(ndm == 2){
 /*...*/
        ss = sqrt(s[0]*s[0]+s[1]*s[1]);
/*...................................................................*/

/*... vetor E*/
        e[0] = ss*ksi[0];        
        e[1] = ss*ksi[1];
/*...................................................................*/

/*... vetor T*/
        t[0] = s[0] - e[0];        
        t[1] = s[1] - e[1];        
/*...................................................................*/
      }        
/*...................................................................*/

/*...*/
      else if(ndm == 3){
/*...*/
        ss = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
/*...................................................................*/

/*... vetor E*/
        e[0] = ss*ksi[0];        
        e[1] = ss*ksi[1];
        e[2] = ss*ksi[2];
/*...................................................................*/

/*... vetor T*/
        t[0] = s[0] - e[0];        
        t[1] = s[1] - e[1];        
        t[2] = s[2] - e[2];        
/*...................................................................*/
      }     
    break;
/*...................................................................*/

/*... OverRelaxed*/
    case OVERRELAXED:
      if(ndm == 2){
/*...*/
        ss = s[0]*s[0] + s[1]*s[1];
/*...................................................................*/

/*...*/
        nk = s[0]*ksi[0] + s[1]*ksi[1];
        nk = 1.e0/nk;
/*...................................................................*/

/*... vetor E*/
        tmp  = ss*nk;
        e[0] = tmp*ksi[0];        
        e[1] = tmp*ksi[1];        
/*...................................................................*/

/*... vetor T*/
        t[0] = s[0] - e[0];        
        t[1] = s[1] - e[1];        
/*...................................................................*/
      }
/*...................................................................*/

/*...*/
      else if(ndm == 3){
/*...*/
        ss = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
/*...................................................................*/

/*...*/
        nk = s[0]*ksi[0] + s[1]*ksi[1] + s[2]*ksi[2];
/*...................................................................*/

/*... vetor E*/
        tmp  = ss/nk;
        e[0] = tmp*ksi[0];        
        e[1] = tmp*ksi[1];        
        e[2] = tmp*ksi[2];        
/*...................................................................*/

/*... vetor T*/
        t[0] = s[0] - e[0];        
        t[1] = s[1] - e[1];        
        t[2] = s[2] - e[2];        
/*...................................................................*/
      }
    break;
/*...................................................................*/

/*...*/
    default:
      printf("Erro: tipo de correcao nao ortogonal invalida.\n"
             "Arquivo fonte:  \"%s\".\n"
             "Nome da funcao: \"%s\".\n"
             "Linha         : \"%d\".\n"
             ,__FILE__,__func__,__LINE__);
      printf("Funcoes disponiveis:\n");
      for(i=0;i<4;i++)
        printf("%s\n",word[i]);
      exit(EXIT_FAILURE);
/*...................................................................*/
  }
}
/*********************************************************************/ 

/*********************************************************************
* Data de criacao    : 15/08/2016                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* difusionSchemeAnisotropic : correcao nao ortogonal do termo       * 
* difusivo anisotropico                                             *                                                      *  
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* s       -> vetor de area modificado                               *
* ksi     -> vetor que une os centroide central e o vizinho         *
* e       -> nao definido                                           *
* t       -> nao definido                                           *
* ndm     -> dimensao                                               *
* iCodDif -> codigo da tecnica nao ortogonal                        *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* e       -> parcela implicita                                      *
* t       -> parcela explicita                                      *
*-------------------------------------------------------------------*
* OBS:                                                              *
* grad(phi)*Sl = grad(phi)*E + grad(phi)*T                          *
*             | k11 k21 k31 | | S1 |   | k11S1 + k21S2 + k31S3 |    *
* Sl = kt*S = | k12 k22 k32 | | S2 | = | k12S1 + k22S2 + k32S3 |    *
*             | k13 k23 k33 | | S3 |   | k13S1 + k23S2 + k33S3 |    *
*-------------------------------------------------------------------*
*********************************************************************/
void difusionSchemeAnisotropic(DOUBLE *RESTRICT s,DOUBLE *RESTRICT ksi
                              ,DOUBLE *RESTRICT e,DOUBLE *RESTRICT t
                              ,short const ndm   ,short const iCod) {

  DOUBLE nk, tmp,ss;
  short i;
  char word[][WORD_SIZE] =
  { "Orthogonal"  ,"Minimal"
    ,"OrthogonalC" ,"OverRelaxed" };

  switch (iCod) {

/*... Minimal*/
  case ORTHOGONAL:
    if (ndm == 2) {
/*...*/
      ss = sqrt(s[0]*s[0]+s[1]*s[1]);
/*...................................................................*/

/*... vetor E*/
      e[0] = ss*ksi[0];
      e[1] = ss*ksi[1];
/*...................................................................*/

/*... vetor T*/
      t[0] = 0.e0;
      t[1] = 0.e0;
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    else if (ndm == 3) {
/*...*/
      ss = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
/*...................................................................*/

/*... vetor E*/
      e[0] = ss*ksi[0];
      e[1] = ss*ksi[1];
      e[2] = ss*ksi[2];
/*...................................................................*/

/*... vetor T*/
      t[0] = 0.e0;
      t[1] = 0.e0;
      t[2] = 0.e0;
/*...................................................................*/
    }
/*...................................................................*/
  break;
/*...................................................................*/

/*... Minimal*/
  case MINIMAL:
  break;
/*...................................................................*/

/*... Orthogonal*/
  case ORTHOGONALC:
  break;
/*...................................................................*/

/*... OverRelaxed*/
  case OVERRELAXED:
    if (ndm == 2) {
/*...*/
      ss = s[0]*s[0] + s[1]*s[1];
/*...................................................................*/

/*...*/
      nk = s[0]*ksi[0] + s[1]*ksi[1];
/*...................................................................*/

/*... vetor E*/
      tmp  = ss/nk;
      e[0] = tmp*ksi[0];
      e[1] = tmp*ksi[1];
/*...................................................................*/

/*... vetor T*/
      t[0] = s[0] - e[0];
      t[1] = s[1] - e[1];
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    else if(ndm == 3) {

/*...*/
      ss = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
/*...................................................................*/

/*...*/
      nk = s[0]*ksi[0] + s[1]*ksi[1] + s[2]*ksi[2]; 
/*...................................................................*/
 
/*... vetor E*/
      tmp = ss / nk;
      e[0] = tmp*ksi[0];
      e[1] = tmp*ksi[1];
      e[2] = tmp*ksi[2];
/*...................................................................*/

/*... vetor T*/
      t[0] = s[0] - e[0];
      t[1] = s[1] - e[1];
      t[2] = s[2] - e[2];
/*...................................................................*/
    }
/*...................................................................*/
  break;
/*...................................................................*/

/*...*/
  default:
  printf("Erro: tipo de correcao nao ortogonal invalida.\n"
         "Arquivo fonte:  \"%s\".\n"
         "Nome da funcao: \"%s\".\n"
         "Linha         : \"%d\".\n"
         , __FILE__, __func__, __LINE__);
  printf("Funcoes disponiveis:\n");
  for (i = 0; i<4; i++)
    printf("%s\n", word[i]);
  exit(EXIT_FAILURE);
/*...................................................................*/
  }
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 04/05/2019                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* advecticeShcemeNdim: discretizacao do termo advectivo             *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* uC        -> valor da celula central                              *
* UF        -> valor da celula vizinha                              *
* gradUc    -> gradiente da celula central                          *
* gradUf    -> gradiente da celula vizinha                          *
* gradUcomp -> gradiente compacto na face                           *
* vSkew       -> vetor entre o ponto medio a intersecao que une os  *
*                centrois compartilhado nessa face da celula central*
* rC          -> distancia do centroide central ate o ponto medio da*
*                face                                               * 
* rV          -> distancia do centroide vizinho ate o ponto medio da*
*                face                                               *
* Ksi         -> vetores que unem centroide da celula central aos   *
*                vizinhos destas                                    *
* mKsi        -> modulo do vetor ksi                                *
* m           -> velocidade normal                                  *
* cvc         -> nao definido                                       *
* alphaMenosUm -> interpolocao linear                               *      
* alfa      -> interpolocao linear                                  *   
* parameters  -> paramentros das tecnicas de advecao                *
* ndm         -> dimensao                                           *
* ndf         -> graus de liberdade                                 *
* iCod1       -> tipo de tecnica                                    *
* iCod2       -> funcoes NVD ou TVD                                 *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* cvc         -> inteporlacao                                       *
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
void advectiveSchemeNdim(DOUBLE *RESTRICT uC ,DOUBLE *RESTRICT uV
                ,DOUBLE *RESTRICT gradUc     ,DOUBLE *RESTRICT gradUv
                ,DOUBLE *RESTRICT gradUComp  ,DOUBLE *RESTRICT vSkew
                ,DOUBLE *RESTRICT rC         ,DOUBLE *RESTRICT rV
                ,DOUBLE *RESTRICT ksi        ,DOUBLE const modKsi
                ,DOUBLE const wfn            ,DOUBLE *RESTRICT cvc 
                ,DOUBLE const alphaMenosUm   ,DOUBLE const alpha
                ,DOUBLE *RESTRICT parameters ,short const ndm
                ,short const ndf               
                ,short const iCod1           ,short const iCod2) {

  short i,j;
  DOUBLE tmp;
  char word[][WORD_SIZE] =
  { "FoUp","Cd"
   ,"SoUp","Tvd"
   ,"Nvd" ,"LUST"};
  switch(iCod1){
/*... Upwind de primeira ordem*/
    case FOUP:
      for(i=0;i<ndf;i++)
        cvc[i] = 0.e0;
    break;
/*...................................................................*/

/*... metodo centrado  atraso( up(implicito) + (ucd - up)explicito) */
  case CD:
/*...*/
    for(i=0;i<ndf;i++)
    {
      cvc[i] = deferredCd(uC[i], uV[i], wfn);

/*... interpolacao undirecional*/
/*    cvc[i] -= MAT2D(i,0,gradVelComp,3)*vSkew[0]
              + MAT2D(i,1,gradVelComp,3)*vSkew[1] 
              + MAT2D(i,2,gradVelComp,3)*vSkew[2];
/*...................................................................*/
    }
/*...................................................................*/
 
  break;
/*...................................................................*/

/*... metodo upwind linear =  up(implicito) + gradU*r*/
  case SOUP:
/*...*/
    for(i=0;i<ndf;i++)
      cvc[i] = upwindLinearV1(uC[i]         , uV[i]
                             ,&gradUc[i*ndm],&gradUv[i*ndm]
                             ,rC            ,rV 
                             ,wfn           ,ndm);   

/*...................................................................*/
  break;
/*...................................................................*/

/*... upwind + termo anti-difusivos*/
  case TVD:
/*...*/
    for(i=0;i<ndf;i++)
    {
      cvc[i] = faceBaseTvdV1(uC[i]         ,uV[i]
                           ,&gradUc[i*ndm] ,&gradUv[i*ndm]
                            ,ksi           ,modKsi
                            ,wfn 
                            ,iCod2         ,ndm);
/*... interpolacao undirecional*/
/*    for(j=0,tmp = 0.e0;j<ndm;j++)
        tmp += MAT2D(i,j,gradUComp,ndm)*vSkew[j];
      cvc[i] -= tmp;*/
    } 
/*...................................................................*/
  break;
/*...................................................................*/

/*... NVD*/
  case NVD:
/*...*/
    for(i=0;i<ndf;i++)
    {
      cvc[i] = faceBaseNvd(uC[i]          , uV[i]
                          ,&gradUc[i*ndm] ,&gradUv[i*ndm]
                          ,ksi, modKsi
                          ,wfn
                          ,iCod2, ndm);
/*... interpolacao undirecional*/
      tmp = 0.e0;
      for(j=0;j<ndm;j++)
        tmp += MAT2D(i,j,gradUComp,ndm)*vSkew[j];
      cvc[i] -= tmp;
    }
/*...................................................................*/
  break;
/*...................................................................*/

/*... NVD*/
  case LUST:
/*...*/
    for(i=0;i<ndf;i++)
      cvc[i] = deferredLust(uC[i]            ,uV[i]
                            ,&gradUc[i*ndm]  ,&gradUv[i*ndm]
                            ,rC              ,rV  
                            ,alphaMenosUm    ,alpha   
                            ,parameters[0]   
                            ,wfn             ,ndm);  
/*...................................................................*/

  break;
/*...................................................................*/

/*...*/
  default:
    printf("Erro: tipo de tecnica de adveccao nao existente.\n"
         "Arquivo fonte:  \"%s\".\n"
         "Nome da funcao: \"%s\".\n"
         "Linha         : \"%d\".\n"
         , __FILE__, __func__, __LINE__);
    printf("Funcoes disponiveis:\n");
    for (i = 0; i<4; i++)
      printf("%s\n", word[i]);
    exit(EXIT_FAILURE);
/*...................................................................*/
  }
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 08/07/2016                                   *
* Data de modificaco : 05/02/2018                                   *
*-------------------------------------------------------------------*
* advecticeShceme: discretizacao do termo advectivo                 *
*-------------------------------------------------------------------*
* Parametros de entrada:                                            *
*-------------------------------------------------------------------*
* velC        -> velocidade da celula central                       *
* velF        -> velocidade da celula vizinha                       *
* gradVelC    -> gradiente da velocidade da celula central          * 
* gradVelV    -> gradiente da velocidade da celula vizinha          *
* gradVelComp -> gradiente de velocidade compacto na face           *
* vSkew       -> vetor entre o ponto medio a intersecao que une os  *
*                centrois compartilhado nessa face da celula central*
* rC          -> distancia do centroide central ate o ponto medio da*
*                face                                               * 
* rV          -> distancia do centroide vizinho ate o ponto medio da*
*                face                                               *
* Ksi         -> vetores que unem centroide da celula central aos   *
*                vizinhos destas                                    *
* mKsi        -> modulo do vetor ksi                                *
* m           -> velocidade normal                                  *
* cvc         -> nao definido                                       *
* alphaMenosUm -> interpolocao linear                               *      
* alfa      -> interpolocao linear                                  *   
* parameters  -> paramentros das tecnicas de advecao                *
* ndm         -> dimensao                                           *
* iCod1       -> tipo de tecnica                                    *
* iCod2       -> funcoes NVD ou TVD                                 *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* cvc         -> inteporlacao                                       *
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
*********************************************************************/
void advectiveScheme(DOUBLE *RESTRICT velC   ,DOUBLE *RESTRICT velV
                ,DOUBLE *RESTRICT gradVelC   ,DOUBLE *RESTRICT gradVelV
                ,DOUBLE *RESTRICT gradVelComp,DOUBLE *RESTRICT vSkew
                ,DOUBLE *RESTRICT rC         ,DOUBLE *RESTRICT rV
                ,DOUBLE *RESTRICT ksi        ,DOUBLE const modKsi
                ,DOUBLE const wfn            ,DOUBLE *RESTRICT cvc 
                ,DOUBLE const alphaMenosUm   ,DOUBLE const alpha
                ,DOUBLE *RESTRICT parameters ,short const ndm             
                ,short const iCod1           ,short const iCod2) {

  short i;
  char word[][WORD_SIZE] =
  { "FoUp","Cd"
   ,"SoUp","Tvd"
   ,"Nvd" ,"LUST"};
  switch(iCod1){
/*... Upwind de primeira ordem*/
    case FOUP:
      cvc[0] = 0.e0;
      cvc[1] = 0.e0;
      if (ndm == 3) cvc[2] = 0.e0;
    break;
/*...................................................................*/

/*... metodo centrado  atraso( up(implicito) + (ucd - up)explicito) */
  case CD:
/*...*/
    if (ndm == 2) {
      cvc[0] = deferredCd(velC[0], velV[0], wfn);
      cvc[1] = deferredCd(velC[1], velV[1], wfn);
/*... interpolacao undirecional*/
      cvc[0] -= MAT2D(0,0,gradVelComp,2)*vSkew[0]
              + MAT2D(0,1,gradVelComp,2)*vSkew[1];
   
      cvc[1] -= MAT2D(1,0,gradVelComp,2)*vSkew[0]
              + MAT2D(1,1,gradVelComp,2)*vSkew[1];
/*...................................................................*/
  }
/*...................................................................*/

/*...*/
    else if (ndm == 3) {
      cvc[0] = deferredCd(velC[0], velV[0], wfn);
      cvc[1] = deferredCd(velC[1], velV[1], wfn);
      cvc[2] = deferredCd(velC[2], velV[2], wfn);
/*... interpolacao undirecional*/
/*    cvc[0] -= MAT2D(0,0,gradVelComp,3)*vSkew[0]
              + MAT2D(0,1,gradVelComp,3)*vSkew[1] 
              + MAT2D(0,2,gradVelComp,3)*vSkew[2];

      cvc[1] -= MAT2D(1,0,gradVelComp,3)*vSkew[0]
              + MAT2D(1,1,gradVelComp,3)*vSkew[1]
              + MAT2D(1,2,gradVelComp,3)*vSkew[2];
      
      cvc[2] -= MAT2D(2,0,gradVelComp,3)*vSkew[0]
              + MAT2D(2,1,gradVelComp,3)*vSkew[1]
              + MAT2D(2,2,gradVelComp,3)*vSkew[2];*/
/*...................................................................*/
    }
/*...................................................................*/
 
  break;
/*...................................................................*/

/*... metodo upwind linear =  up(implicito) + gradU*r*/
  case SOUP:
/*...*/
    if (ndm == 2) {
      cvc[0] = upwindLinearV1(velC[0]    ,velV[0]
                             ,gradVelC   ,gradVelV
                             ,rC         ,rV 
                             ,wfn        ,ndm);
/*...*/
      cvc[1] = upwindLinearV1(velC[1]      ,velV[1]
                             ,&gradVelC[2],&gradVelV[2]
                             ,rC          ,rV 
                             ,wfn         ,ndm);
    }
/*...................................................................*/

/*...*/
    else if (ndm == 3) {
      cvc[0] = upwindLinearV1(velC[0]    ,velV[0]
                             ,gradVelC   ,gradVelV
                             ,rC         ,rV 
                             ,wfn        ,ndm);
/*...*/
      cvc[1] = upwindLinearV1(velC[1]     ,velV[1]
                             ,&gradVelC[3],&gradVelV[3]
                             ,rC          ,rV 
                             ,wfn         ,ndm);
/*...*/
      cvc[2] = upwindLinearV1(velC[2]     ,velV[2]
                             ,&gradVelC[6],&gradVelV[6]
                             ,rC          ,rV  
                             ,wfn         ,ndm);
  }
/*...................................................................*/
  break;
/*...................................................................*/

/*... upwind + termo anti-difusivos*/
  case TVD:
/*...*/
    if (ndm == 2) {
      cvc[0] = faceBaseTvdV1(velC[0]     ,velV[0]
                            ,&gradVelC[0],&gradVelV[0]
                            ,ksi         ,modKsi
                            ,wfn 
                            ,iCod2       ,ndm);
 
      cvc[1] = faceBaseTvdV1(velC[1]     ,velV[1]
                            ,&gradVelC[2],&gradVelV[2]
                            ,ksi         ,modKsi
                            ,wfn 
                            ,iCod2       ,ndm);

/*... interpolacao undirecional*/
      cvc[0] -= MAT2D(0,0,gradVelComp,2)*vSkew[0]
              + MAT2D(0,1,gradVelComp,2)*vSkew[1];
   
      cvc[1] -= MAT2D(1,0,gradVelComp,2)*vSkew[0]
              + MAT2D(1,1,gradVelComp,2)*vSkew[1];
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    else if (ndm == 3) {
      cvc[0] = faceBaseTvdV1(velC[0]      ,velV[0]
                            ,&gradVelC[0] ,&gradVelV[0]
                            ,ksi          ,modKsi
                            ,wfn 
                            ,iCod2        ,ndm);
      cvc[1] = faceBaseTvdV1(velC[1]     ,velV[1]
                            ,&gradVelC[3],&gradVelV[3]
                            ,ksi         ,modKsi
                            ,wfn  
                            ,iCod2      ,ndm);
      cvc[2] = faceBaseTvdV1(velC[2]     ,velV[2]
                            ,&gradVelC[6],&gradVelV[6]
                            ,ksi         ,modKsi
                            ,wfn 
                            ,iCod2       ,ndm);
/*... interpolacao undirecional*/
      cvc[0] -= MAT2D(0,0,gradVelComp,3)*vSkew[0]
              + MAT2D(0,1,gradVelComp,3)*vSkew[1]
              + MAT2D(0,2,gradVelComp,3)*vSkew[2];

      cvc[1] -= MAT2D(1,0,gradVelComp,3)*vSkew[0]
              + MAT2D(1,1,gradVelComp,3)*vSkew[1]
              + MAT2D(1,2,gradVelComp,3)*vSkew[2];

      cvc[2] -= MAT2D(2,0,gradVelComp,3)*vSkew[0]
              + MAT2D(2,1,gradVelComp,3)*vSkew[1]
              + MAT2D(2,2,gradVelComp,3)*vSkew[2];
/*...................................................................*/
    }
/*...................................................................*/
  break;
/*...................................................................*/

/*... NVD*/
  case NVD:
/*...*/
    if (ndm == 2) {
      cvc[0] = faceBaseNvd(velC[0], velV[0]
                          ,&gradVelC[0], &gradVelV[0]
                          ,ksi, modKsi
                          ,wfn
                          ,iCod2, ndm);
      cvc[1] = faceBaseNvd(velC[1], velV[1]
                          , &gradVelC[2], &gradVelV[2]
                          , ksi, modKsi
                          , wfn
                          , iCod2, ndm);

/*... interpolacao undirecional*/
      cvc[0] -= MAT2D(0, 0, gradVelComp, 2)*vSkew[0]
              + MAT2D(0, 1, gradVelComp, 2)*vSkew[1];

      cvc[1] -= MAT2D(1, 0, gradVelComp, 2)*vSkew[0]
              + MAT2D(1, 1, gradVelComp, 2)*vSkew[1];
/*...................................................................*/
    }
/*...................................................................*/

/*...*/
    else if (ndm == 3) {
      cvc[0] = faceBaseNvd(velC[0]     ,velV[0]
                          ,&gradVelC[0],&gradVelV[0]
                          ,ksi         ,modKsi
                          ,wfn
                          ,iCod2       ,ndm);
      cvc[1] = faceBaseNvd(velC[1]     ,velV[1]
                          ,&gradVelC[3],&gradVelV[3]
                          ,ksi         ,modKsi
                          ,wfn
                          ,iCod2       ,ndm);
      cvc[2] = faceBaseNvd(velC[2]     ,velV[2]
                          ,&gradVelC[6],&gradVelV[6]
                          ,ksi         ,modKsi
                          ,wfn
                          ,iCod2       ,ndm);
/*... interpolacao undirecional*/
/*    cvc[0] -= MAT2D(0, 0, gradVelComp, 3)*vSkew[0]
              + MAT2D(0, 1, gradVelComp, 3)*vSkew[1]
              + MAT2D(0, 2, gradVelComp, 3)*vSkew[2];

      cvc[1] -= MAT2D(1, 0, gradVelComp, 3)*vSkew[0]
              + MAT2D(1, 1, gradVelComp, 3)*vSkew[1]
              + MAT2D(1, 2, gradVelComp, 3)*vSkew[2];

      cvc[2] -= MAT2D(2, 0, gradVelComp, 3)*vSkew[0]
              + MAT2D(2, 1, gradVelComp, 3)*vSkew[1]
              + MAT2D(2, 2, gradVelComp, 3)*vSkew[2];*/
/*...................................................................*/
    }
/*...................................................................*/
  break;
/*...................................................................*/

/*... NVD*/
  case LUST:
/*...*/
    cvc[0] = deferredLust(velC[0]     ,velV[0]
                          ,&gradVelC[0],&gradVelV[0]
                          ,rC          ,rV  
                          ,alphaMenosUm,alpha   
                          ,parameters[0]   
                          ,wfn         ,ndm);
    
    cvc[1] = deferredLust(velC[1]     ,velV[1]
                          ,&gradVelC[3],&gradVelV[3]
                          ,rC          ,rV  
                          ,alphaMenosUm,alpha  
                          ,parameters[0]     
                          ,wfn         ,ndm);
    if (ndm == 3)
      cvc[2] = deferredLust(velC[2]     ,velV[2]
                            ,&gradVelC[6],&gradVelV[6]
                            ,rC          ,rV  
                            ,alphaMenosUm,alpha     
                            ,parameters[0]   
                            ,wfn         ,ndm);
/*...................................................................*/

  break;
/*...................................................................*/

/*...*/
  default:
    printf("Erro: tipo de tecnica de adveccao nao existente.\n"
         "Arquivo fonte:  \"%s\".\n"
         "Nome da funcao: \"%s\".\n"
         "Linha         : \"%d\".\n"
         , __FILE__, __func__, __LINE__);
    printf("Funcoes disponiveis:\n");
    for (i = 0; i<4; i++)
      printf("%s\n", word[i]);
    exit(EXIT_FAILURE);
/*...................................................................*/
  }
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 22/08/2016                                   *
* Data de modificaco : 07/02/2018                                   *
*-------------------------------------------------------------------*
* advecticeShcemeScalar: discretizacao do termo advectivo           *
*-------------------------------------------------------------------*
* uC        -> valor da celula central                              *
* UF        -> valor da celula vizinha                              *
* gradUc    -> gradiente da celula central                          *
* gradUf    -> gradiente da celula vizinha                          *
* gradUcomp -> gradiente compacto na face                           *
* vSkew     -> vetor entre o ponto medio a intersecao que une os    *
*              centrois compartilhado nessa face da celula central  *
* rC        -> distancia do centroide central ate o ponto medio da  *
*              face                                                 *
* rV        -> distancia do centroide vizinho ate o ponto medio da  *
*              face                                                 *
* Ksi       -> vetores que unem centroide da celula central aos     *
*              vizinhos destas                                      *
* mKsi      -> modulo do vetor ksi                                  *
* m         -> velocidade normal a face                             *
* cvc       -> nao definido                                         *
* alphaMenosUm -> interpolocao linear                               *      
* ndm       -> dimensao                                             *
* alphaMenosUm -> interpolocao linear                               *      
* alfa      -> interpolocao linear                                  *                                    *
* iCod1     -> tipo de tecnica                                      *
* iCod2     -> funcoes NVD ou TVD                                   *
*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* cvc         -> inteporlacao                                       *
*-------------------------------------------------------------------*
* OBS:                                                              *
*-------------------------------------------------------------------*
* Uface = (alpha-1)*Up + alpha*Uv                                   * 
*********************************************************************/
void advectiveSchemeScalar(DOUBLE const uC, DOUBLE const uV
              ,DOUBLE *RESTRICT gradUc    ,DOUBLE *RESTRICT gradUv
              ,DOUBLE *RESTRICT gradUcomp ,DOUBLE *RESTRICT vSkew
              ,DOUBLE *RESTRICT rC        ,DOUBLE *RESTRICT rV
              ,DOUBLE *RESTRICT ksi       ,DOUBLE const modKsi
              ,DOUBLE const wfn           ,DOUBLE *cvc
              ,DOUBLE const alphaMenosUm  ,DOUBLE const alpha          
              ,DOUBLE *RESTRICT parameters ,short const ndm
              ,short const iCod1          ,short const iCod2) {

  short i;
  char word[][WORD_SIZE] =
  { "FoUp","Cd"
   ,"SoUp","Tvd" 
   ,"Nvd" ,"LUST"};

  *cvc = 0.e0;
  switch (iCod1) {
/*... Upwind de primeira ordem*/
    case FOUP:
      *cvc = 0.e0;
  break;
/*...................................................................*/

/*... metodo centrado  atraso( up(implicito) + (ucd - up)explicito) */
  case CD:
/*...*/
    *cvc = deferredCd(uC,uV,wfn);
/*...................................................................*/

/*... interpolacao undirecional*/
    if(ndm == 2)
      *cvc -= gradUcomp[0]*vSkew[0] + gradUcomp[1]*vSkew[1];
    else if(ndm == 3)
      *cvc -= gradUcomp[0]*vSkew[0]   
            + gradUcomp[1]*vSkew[1] 
            + gradUcomp[2]*vSkew[2];
/*...................................................................*/
  break;
/*...................................................................*/

/*... metodo upwind linear =  up(implicito) + gradU*r*/
  case SOUP:
    *cvc = upwindLinearV1(uC       ,uV
                        ,gradUc    ,gradUv
                        ,rC        ,rV       
                        ,wfn       ,ndm);
  break;
/*...................................................................*/

/*... upwind + termo anti-difusivos*/
  case TVD:
/*...*/
    *cvc = faceBaseTvdV1(uC     ,uV
                       ,gradUc ,gradUv
                       ,ksi    ,modKsi
                       ,wfn
                       ,iCod2  ,ndm);
/*...................................................................*/

/*... interpolacao undirecional*/
    if (ndm == 2)
      *cvc -= gradUcomp[0]*vSkew[0] + gradUcomp[1]*vSkew[1];
    else if (ndm == 3)
      *cvc -= gradUcomp[0]*vSkew[0]
            + gradUcomp[1]*vSkew[1]
            + gradUcomp[2]*vSkew[2];
/*...................................................................*/
  break;
/*...................................................................*/

/*... NVD*/
  case NVD:
/*...*/
    *cvc = faceBaseNvd(uC    ,uV
                      ,gradUc,gradUv
                      ,ksi   ,modKsi
                      ,wfn
                      ,iCod2 ,ndm);
/*...................................................................*/

/*... interpolacao undirecional*/
/*  if (ndm == 2)
      *cvc -= gradUcomp[0]*vSkew[0] + gradUcomp[1]*vSkew[1];
    else if (ndm == 3)
      *cvc -= gradUcomp[0]*vSkew[0]
            + gradUcomp[1]*vSkew[1]
            + gradUcomp[2]*vSkew[2];*/
/*...................................................................*/
  break;
/*...................................................................*/

/*... LUST*/
  case LUST:
/*...*/
    *cvc =  deferredLust(uC           ,uV
                        ,gradUc       ,gradUv
                        ,rC           ,rV  
                        ,alphaMenosUm ,alpha      
                        ,parameters[0]
                        ,wfn          ,ndm);
/*...................................................................*/

  break;
/*...................................................................*/

  break;
/*...................................................................*/


/*...*/
  default:
    printf("Erro: tipo de tecnica de adveccao nao existente.\n"
         "Arquivo fonte:  \"%s\".\n"
         "Nome da funcao: \"%s\".\n"
         "Linha         : \"%d\".\n"
         , __FILE__, __func__, __LINE__);
    printf("Funcoes disponiveis:\n");
    for (i = 0; i<4; i++)
      printf("%s\n", word[i]);
    exit(EXIT_FAILURE);
/*...................................................................*/
  }
}
/*********************************************************************/

/*********************************************************************
* Data de criacao    : 25/08/2017                                   *
* Data de modificaco : 00/00/0000                                   *
*-------------------------------------------------------------------*
* INTERPOLFACEVEL : interpolacao das velocidae na faces             *
*-------------------------------------------------------------------*
* velC      -> velocidade da celula central                         *
* velV      -> velocidade celula vizinha                            *
* presC     -> pressao da celula central                            *
* presV     -> pressao da celula vizinha                            *
* gradPesC  -> gradiente reconstruido da pressao da celula central  *
* gradPesV  -> gradiente reconstruido da pressao da celula vizinha  *
* vKsi      -> vetores que unem centroide da celula central aos     *
*              vizinhos destas                                      *
* mKsi      -> modulo do vetor ksi                                  *
* dField    -> matriz D do metodo simple                            *
* alphaMenosUm -> interpolacao linear                               *
* alpha        -> interpolacao linear                               *

*-------------------------------------------------------------------*
* Parametros de saida:                                              *
*-------------------------------------------------------------------*
* wfn -> velocidade normal a face                                   *
*-------------------------------------------------------------------*
* OBS:                                                              *
* 2016 - Darwish-Moukalled  
*-------------------------------------------------------------------*
*********************************************************************/
DOUBLE interpolFaceVel(DOUBLE *RESTRICT velC,DOUBLE *RESTRICT velV
              ,DOUBLE const presC            ,DOUBLE const presV
              ,DOUBLE *RESTRICT gradPresC    ,DOUBLE *RESTRICT gradPresV
              ,DOUBLE *RESTRICT lNormal      ,DOUBLE *RESTRICT vKsi
              ,DOUBLE const mKsi             ,DOUBLE *RESTRICT dFieldF
              ,DOUBLE const alphaMenosUm     ,DOUBLE const alpha
              ,short const ndm)
{

  DOUBLE wfn,nk,wf[3],gf[3],gfp[3];

  wfn = gfp[0] = gfp[1] = gfp[2] = 0.e0;
  if( ndm == 2 ){
/*...*/
    wf[0] = alphaMenosUm*velC[0] + alpha*velV[0];
    wf[1] = alphaMenosUm*velC[1] + alpha*velV[1]; 
/*...................................................................*/

/*...................................................................*/
    wfn = wf[0] * lNormal[0] + wf[1] * lNormal[1];
/*...................................................................*/

/*... 2016 - Darwish-Moukalled*/
/*... interpolacao linear dos gradientes das pressoes*/
    gf[0] = alphaMenosUm*gradPresC[0] + alpha*gradPresV[0];
    gf[1] = alphaMenosUm*gradPresC[1] + alpha*gradPresV[1];
/*...................................................................*/

/*...*/
    gfp[0] = (presV - presC) / mKsi;
/*...................................................................*/

/*... produtos interno*/
    nk = dFieldF[0] * vKsi[0] * lNormal[0]
       + dFieldF[1] * vKsi[1] * lNormal[1];
    gfp[1] = gf[0] * vKsi[0] + gf[1] * vKsi[1];
/*...................................................................*/

/*...*/
    gfp[1] = (gfp[1] - gfp[0])*nk;
/*...................................................................*/
  }
/*...................................................................*/

  else if( ndm == 3 ){
/*...*/
    wf[0] = alphaMenosUm*velC[0] + alpha*velV[0];
    wf[1] = alphaMenosUm*velC[1] + alpha*velV[1]; 
    wf[2] = alphaMenosUm*velC[2] + alpha*velV[2]; 
/*...................................................................*/

/*...................................................................*/
    wfn = wf[0] * lNormal[0] 
        + wf[1] * lNormal[1]
        + wf[2] * lNormal[2];
/*...................................................................*/

/*... 2016 - Darwish-Moukalled*/
/*... interpolacao linear dos gradientes das pressoes*/
    gf[0] = alphaMenosUm*gradPresC[0] + alpha*gradPresV[0];
    gf[1] = alphaMenosUm*gradPresC[1] + alpha*gradPresV[1];
    gf[2] = alphaMenosUm*gradPresC[2] + alpha*gradPresV[2];
/*...................................................................*/

/*...*/
    gfp[0] = (presV - presC) / mKsi;
/*...................................................................*/

/*... produtos interno*/
    nk = dFieldF[0] * vKsi[0] * lNormal[0]
       + dFieldF[1] * vKsi[1] * lNormal[1]
       + dFieldF[2] * vKsi[2] * lNormal[2];
    gfp[1] = gf[0] * vKsi[0] + gf[1] * vKsi[1] + gf[2] * vKsi[2];
/*...................................................................*/

/*...*/
    gfp[1] = (gfp[1] - gfp[0])*nk;
/*...................................................................*/
  }
/*...................................................................*/

/*...*/ 
  wfn += gfp[1];
/*...................................................................*/

/*...*/
  return wfn;
/*...................................................................*/
}
/*********************************************************************/ 

/********************************************************************** 
 * Data de criacao    : 11/11/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 *------------------------------------------------------------------- *
 * vorticity: calcula a vorticidade                                   *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * w       -> nao definido                                            * 
 * gradVel -> gradienta das velocidades                               * 
 * ndm     -> dimensao                                                * 
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * w       -> vorticidade                                             *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *------------------------------------------------------------------- *
 *                                                                    *
 * gradVel(ndf,ndm) -> gradVel(3,3)                                   *
 *                                                                    *
 *             | du1dx1 du1dx2 du1dx3 |                               * 
 * grad(*,*) = | du2dx1 du2dx2 du2dx3 |                               *
 *             | du3dx1 du3dx2 du3dx3 |                               *
 *                                                                    *
 **********************************************************************/
void vorticity(DOUBLE *RESTRICT w,DOUBLE *RESTRICT gradVel
              ,const short ndm) {

  if (ndm == 3) {
/*... dv3/dx2 - dv2/dx3*/
    w[0] = MAT2D(2,1,gradVel,ndm) - MAT2D(1,2,gradVel,ndm);
/*... dv1/dx3 - dv3/dx1*/ 
    w[1] = MAT2D(0,2,gradVel,ndm) - MAT2D(2,0,gradVel,ndm);
/*... dv2/dx1 - dv1/dx2*/
    w[2] = MAT2D(1,0,gradVel,ndm) - MAT2D(0,1,gradVel,ndm);
  }

}
/*********************************************************************/ 

/********************************************************************** 
 * Data de criacao    : 25/11/2017                                    *
 * Data de modificaco : 09/12/2017                                    *
 *------------------------------------------------------------------- *
 * stress: tensor de forcas viscosas                                  *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * s       -> nao definido                                            * 
 * gradVel -> gradienta das velocidades                               * 
 * nu      -> viscosidade molecular                                   * 
 * lambda  -> coeficiente de formacao linear                          * 
 * ndm     -> dimensao                                                * 
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * s       -> tensor                                                  *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *------------------------------------------------------------------- *
 *                                                                    *
 * gradVel(ndf,ndm) -> gradVel(3,3)                                   *
 *                                                                    *
 *             | du1dx1 du1dx2 du1dx3 |                               * 
 * grad(*,*) = | du2dx1 du2dx2 du2dx3 |                               *
 *             | du3dx1 du3dx2 du3dx3 |                               *
 *                                                                    *
 **********************************************************************/
void stress(DOUBLE *RESTRICT s,DOUBLE *RESTRICT gradVel
           ,DOUBLE const nu   ,DOUBLE const lambda
           ,short const ndm) {

  DOUBLE tmp;

  if (ndm == 3) {
    tmp = MAT2D(0,0,gradVel,ndm)
        + MAT2D(1,1,gradVel,ndm) 
        + MAT2D(2,2,gradVel,ndm);
/*... s11 - s0*/
    s[0] = 2.e0*nu*MAT2D(0,0,gradVel,ndm) + lambda*tmp;
/*... s22 - s1*/
    s[1] = 2.e0*nu*MAT2D(1,1,gradVel,ndm) + lambda*tmp;
/*... s33 - s2*/
    s[2] = 2.e0*nu*MAT2D(2,2,gradVel,ndm) + lambda*tmp;
/*... s12 - s3*/
    s[3] = nu*(MAT2D(1,0,gradVel,ndm) + MAT2D(0,1,gradVel,ndm));
/*... s23 - s4*/
    s[4] = nu*(MAT2D(2,1,gradVel,ndm) + MAT2D(1,2,gradVel,ndm));
/*... s13 - s5*/
    s[5] = nu*(MAT2D(2,0,gradVel,ndm) + MAT2D(0,2,gradVel,ndm));

  }

}
/*********************************************************************/

/********************************************************************** 
 * Data de criacao    : 09/12/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 *------------------------------------------------------------------- *
 * stressEddyViscosity : tensor de forcas viscosas turbulenta         *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * s       -> nao definido                                            * 
 * gradVel -> gradienta das velocidades                               * 
 * nu      -> viscosidade turbulenta                                  * 
 * ndm     -> dimensao                                                * 
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * s       -> tensor                                                  *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *------------------------------------------------------------------- *
 *                                                                    *
 * gradVel(ndf,ndm) -> gradVel(3,3)                                   *
 *                                                                    *
 *             | du1dx1 du1dx2 du1dx3 |                               * 
 * grad(*,*) = | du2dx1 du2dx2 du2dx3 |                               *
 *             | du3dx1 du3dx2 du3dx3 |                               *
 *                                                                    *
 **********************************************************************/
void stressEddyViscosity(DOUBLE *RESTRICT s,DOUBLE *RESTRICT gradVel
                       , DOUBLE const nut  , short const ndm) {

  DOUBLE tmp;

  if (ndm == 3) {
     tmp =D1DIV3*(MAT2D(0,0,gradVel,ndm)
                + MAT2D(1,1,gradVel,ndm) 
                + MAT2D(2,2,gradVel,ndm));
/*... s11 - s0*/
    s[0] = -2.e0*nut*(MAT2D(0,0,gradVel,ndm) - tmp);
/*... s22 - s1*/
    s[1] = -2.e0*nut*(MAT2D(1,1,gradVel,ndm) - tmp);
/*... s33 - s2*/
    s[2] = -2.e0*nut*(MAT2D(2,2,gradVel,ndm) - tmp);
/*... s12 - s3*/
    s[3] = -nut*(MAT2D(1,0,gradVel,ndm) + MAT2D(0,1,gradVel,ndm));
/*... s23 - s4*/
    s[4] = -nut*(MAT2D(2,1,gradVel,ndm) + MAT2D(1,2,gradVel,ndm));
/*... s13 - s5*/
    s[5] = -nut*(MAT2D(2,0,gradVel,ndm) + MAT2D(0,2,gradVel,ndm));

  }

}
/*********************************************************************/  

/**********************************************************************
 * Data de criacao    : 13/12/2017                                    *
 * Data de modificaco : 00/00/0000                                    *
 *------------------------------------------------------------------- * 
 * qCriterion : calculo do parametro Q                                *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * gradVel -> gradienbte de velocidades                               *
 * ndm     -> dimensao                                                * 
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 *                                                                    *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *       | 0   w1  w2 |            | s1  s3  s6 |                     *
 *  W =  |-w1   0  w3 |        S = | s3  s2  s4 |                     *
 *       |-w2 -w3   0 |            | s6  s4  s3 |                     * 
 *------------------------------------------------------------------- *
 **********************************************************************/
DOUBLE qCriterion(DOUBLE *RESTRICT gradVel, short const ndm) {

  DOUBLE q=0.0e0,ss,ww,s[6],w[3];


  if (ndm == 3) {

    s[0] = MAT2D(0, 0, gradVel, 3);         /*s11*/
    s[1] = MAT2D(1, 1, gradVel, 3);         /*s22*/
    s[2] = MAT2D(2, 2, gradVel, 3);         /*s33*/ 
    s[3] = 0.5e0*(MAT2D(0, 1, gradVel, 3) + MAT2D(1, 0, gradVel, 3)); /*s12*/
    s[4] = 0.5e0*(MAT2D(1, 2, gradVel, 3) + MAT2D(2, 1, gradVel, 3)); /*s23*/
    s[5] = 0.5e0*(MAT2D(0, 2, gradVel, 3) + MAT2D(2, 0, gradVel, 3)); /*s13*/

    ss = doubleDotSym(s);


    w[0] = 0.5e0*(MAT2D(0, 1, gradVel, 3) - MAT2D(1, 0, gradVel, 3)); /*w12*/
    w[1] = 0.5e0*(MAT2D(0, 2, gradVel, 3) - MAT2D(2, 0, gradVel, 3)); /*w13*/
    w[2] = 0.5e0*(MAT2D(1, 2, gradVel, 3) - MAT2D(1, 2, gradVel, 3)); /*w23*/

    ww = 2.e0*(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

    q = ww - ss;
  }

  return 0.25*q; 
}
/*********************************************************************/  

/**********************************************************************
 * Data de criacao    : 28/02/2018                                    *
 * Data de modificaco : 00/00/0000                                    *
 *------------------------------------------------------------------- * 
 * interpolFace : interpolacao linear de valores na face              *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * vSkew  -> vetor entre o ponto medio a intersecao que une os        *
 *            centrois compartilhado nessa face da celula central     *
 * xmcc      -> vetores que unem o centroide aos pontos medios das    *
 *            faces da celula central                                 * 
 * volP    -> volume da celula central                                *
 * volV    -> volume da celula vizinha                                *
 * mKsi    -> modulo do vetor ksi ( vetor que une os dois centroiodes)*
 * ndm     -> dimensao                                                * 
 * iCod    -> 0 LINEAR                                                *
 *            1 por volume                                            *
 *            2 media                                                 *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 *                                                                    *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *------------------------------------------------------------------- *
 **********************************************************************/
DOUBLE interpolFace(DOUBLE *RESTRICT vSkew, DOUBLE *RESTRICT xmcc
                  , DOUBLE volP            , DOUBLE volV
                  , DOUBLE lModKsi         , short ndm           
                  , short iCod    ) 
{

  DOUBLE v[3],alpha=0.e0,dPviz;


  switch(iCod){

/*... interpolacao usando distancias lineares*/
    case INTPOLFACELINEAR: 
      if (ndm == 3) {
        v[0]         = vSkew[0] + xmcc[0];
        v[1]         = vSkew[1] + xmcc[1];
        v[2]         = vSkew[2] + xmcc[2];
        dPviz        = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        alpha        = dPviz/lModKsi;
      }
      break;
/*....................................................................*/    

/*... interpolacao usando os volumes das celulas*/
    case INTPOLFACEVOL: 
      if (ndm == 3) {
        alpha        = volV/(volP+volV);
      }
      break;
/*....................................................................*/ 

/*... interpolacao usando distancias lineares*/
    case INTPOLFACEMED: 
      alpha        = 0.5e0;
      break;
/*....................................................................*/ 

  }

  return alpha; 
}
/*********************************************************************/   

/**********************************************************************
 * Data de criacao    : 12/05/2019                                    *
 * Data de modificaco : 00/00/0000                                    *
 *------------------------------------------------------------------- * 
 * velCorrectCombustion : corecao da velocidade de difusao no modelo  * 
 * de transporte de especies (combustao)                              *  
 * ------------------------------------------------------------------ *
 * parametros de entrada:                                             * 
 * ------------------------------------------------------------------ *
 * diff   -> coef de difusao das especies                             *
 *            centrois compartilhado nessa face da celula central     *
 * gradZ  -> gradiente das especies                                   *
 * volC    -> nao definido                                            *
 * ndm     -> dimensao                                                * 
 * ns      -> numero de especies                                      *
 * ------------------------------------------------------------------ *
 * parametros de saida  :                                             * 
 * ------------------------------------------------------------------ *
 * velC   -> velocidade de correcao                                   *
 * ------------------------------------------------------------------ *
 * OBS:                                                               *
 *           | du1/dx1 du1/dx2 du1/dx3 |                              *
 * gradZ  =  | du2/dx1 du2/dx2 du2/dx3 |                              *
 *           | du3/dx1 du3/dx2 du3/dx3 |                              *
 *           | du4/dx1 du4/dx2 du4/dx3 |                              *
 *                     ...                                            *
 * Livro : Theoretical and numerical Combustion                       *
 *------------------------------------------------------------------- *
 **********************************************************************/
void velCorrectCombustion(DOUBLE *RESTRICT diff, DOUBLE *RESTRICT gradZ
                  , DOUBLE *RESTRICT velC      , short const ndm           
                  , short const ns    ) 
{
  short i;
  DOUBLE vel[3];


  if(ndm == 3)
  {
    vel[0] = vel[1] = vel[2] = 0.e0;
  
    for(i=0;i<ns;i++)
    {
      vel[0] += diff[i]*MAT2D(i,0,gradZ,3);
      vel[1] += diff[i]*MAT2D(i,1,gradZ,3);
      vel[2] += diff[i]*MAT2D(i,2,gradZ,3);

    }

    velC[0] += vel[0];
    velC[1] += vel[1];
    velC[2] += vel[2];
  }
}
/*********************************************************************/   