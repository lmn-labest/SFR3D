#ifndef _GRAPH_H  
  #define _GRAPH_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
  #include<Define.h>
  #include<HccaStdBool.h>
  #include<HccaSort.h>
  #include<Mesh.h>
/*...*/
  void sortGraphCsr(INT *restrict ia,INT *restrict ja,INT const n);
/*...*/  
  void convGraph(INT *restrict xAdj    ,INT *restrict adjncy
                ,INT const *adj        ,short const *nViz
                ,short const maxViz    ,INT const numel 
                ,bool const xAdjFlag   ,bool const adjFlag );

  void convGraphPart(INT *restrict xAdj,INT *restrict adjncy 
                  ,INT const *adj      ,short const *nViz     
                  ,short const maxViz  ,INT const numelNov
                  ,bool const xAdjFlag ,bool const adjFlag);
/*...................................................................*/

#endif/*_CSR_H*/
