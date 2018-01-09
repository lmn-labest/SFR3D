#ifndef _RANDOM_H_
  #define _RANDOM_H_
/*...*/
  #include<stdlib.h>
  #include<Time.h>
/*...................................................................*/

/*...*/ 
  #include<Define.h>
  #include<HccaStdBool.h>
/*...................................................................*/

/*...*/
  void initRand(bool const flag);
  DOUBLE doubleRand(DOUBLE min, DOUBLE max );
/*...................................................................*/

#endif/*_RANDOM_H_*/