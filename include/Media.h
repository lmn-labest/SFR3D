#ifndef _MEDIA_H_
  #define _MEDIA_H_
/*...*/ 
  #include<stdlib.h>
  #include<stdio.h>
 
/*...*/
  #include<Mesh.h>
/*...................................................................*/

/*...*/
  void calMean(Mean *media    , Mesh *mesh
             , DOUBLE const dt, DOUBLE const t
             , INT const timeStep );
  void integralTemp(Mean *media, Mesh *mesh, DOUBLE const dt); 
  void  calMediaTempFinal(Mean *media, DOUBLE const t
                        , INT const n, short const ndf) ;
  void  prime2(DOUBLE *RESTRICT mPhi , DOUBLE *RESTRICT phi
             , DOUBLE *RESTRICT p2Phi
             , INT const n           , short const ndf); 
/*...................................................................*/

#endif/*_MEDIA_H_*/
