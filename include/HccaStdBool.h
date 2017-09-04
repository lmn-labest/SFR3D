#ifndef _MYSTDBOOL_H_
  #ifdef _MSC_VER  
    #if (_MSC_VER < 1800)
      #define bool short
      #define true  1
      #define false 0 
    #else
      #include<stdbool.h>
    #endif
  #else
    #include<stdbool.h>
  #endif
#endif  /*_MYSTDBOOL_H_*/
