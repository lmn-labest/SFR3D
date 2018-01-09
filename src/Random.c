#include<Random.h>

void initRand(bool const flag){

/* initialize random seed: */
  if(flag)
    srand (time(NULL));
  else
    srand (1);

}


DOUBLE doubleRand(DOUBLE min, DOUBLE max )
{
    DOUBLE scale = rand() / (DOUBLE) RAND_MAX; /* [0, 1.0] */
    return min + scale * ( max - min );      /* [min, max] */
}