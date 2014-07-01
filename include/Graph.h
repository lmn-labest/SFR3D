#ifndef _GRAPH_H  
  #define _GRAPH_H
   #ifdef _AD_
    #undef  _AD_
  #endif  
  #define _AD_ false
  #include<Mystdbool.h>
  #include<Mesh.h>
/*...*/
  void sortGraphCsr(long *ia,long *ja,long n);
/*...*/  
  void convGraph(long *xAdj    ,long *adjncy,long const *adj,short const *nViz
                ,short maxViz  ,long numel  ,bool xAdjFlag  ,bool adjFlag );
/*...*/
  void bubblesort(long *ja,long n);

#endif/*_CSR_H*/
