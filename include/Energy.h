#ifndef _ENERGY_H_
  #define  _ENERGY_H_
/*...*/
  #include<math.h>
/*...................................................................*/

/*...*/
  #include<Define.h>
  #include<Memoria.h>
  #include<Mesh.h>
  #include<Sisteq.h>
  #include<Solv.h>
  #include<Properties.h>
  #include<CellLoop.h>
/*...................................................................*/

/*...*/
bool energyEquation(Memoria *m               , PropVarFluid *prop 
                   , Loads *loadsVel         , Loads *loadsEnergy  
                   , EnergyModel *eModel     , Turbulence *tModel  
                   , ThermoDynamic *thDynamic, Mesh *mesh          
                   , SistEq *sistEqEnergy    , Solv *solvEnergy
                   , Simple *sp  
                   , Scheme *sc              , PartMesh *pMesh);
/*...................................................................*/

#endif /*_ENERGY_H_*/