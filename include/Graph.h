#ifndef _GRAPH_H  
  #define _GRAPH_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
  #include<Define.h>
  #include<HccaStdBool.h>
  #include<HccaSort.h>
//  #include<Mesh.h>
/*...*/
  void sortGraphCsr(INT *RESTRICT ia,INT *RESTRICT ja,INT const n);
/*...*/  
  void convGraph(INT *RESTRICT xAdj    ,INT *RESTRICT adjncy
                ,INT const *adj        ,short const *nViz
                ,short const maxViz    ,INT const numel 
                ,bool const xAdjFlag   ,bool const adjFlag );

  void convGraphPart(INT *RESTRICT xAdj,INT *RESTRICT adjncy 
                  ,INT const *adj      ,short const *nViz     
                  ,short const maxViz  ,INT const numelNov
                  ,bool const xAdjFlag ,bool const adjFlag);
/*...................................................................*/

#endif/*_CSR_H*/
