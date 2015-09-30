#ifndef _RCM_H  
  #define _RCM_H
/*...*/
  #include<stdio.h>
  #include<stdlib.h>
/*...*/
  #include<Define.h>
  #include<Erro.h>
/*...*/
  void genrcm(INT node_num, INT *adj_row, INT *adj
             ,INT *perm);

  void root_find( INT *root , INT *adj_row  , INT *adj
                , char *mask , INT *level_num, INT *level_row
                , INT *level, INT node_num );

  void level_set( INT root  ,INT *adj_row, INT *adj
                , char *mask ,INT *level_num, INT *level_row
                , INT *level,INT node_num );
  
  void degree ( INT root, INT *adj_row, INT *adj, char *mask,
                INT *deg, INT *iccsze, INT *ls, INT node_num );
  
  void ivec_reverse ( INT n, INT *a );

  void rcm ( INT root, INT *adj_row, INT *adj, char *mask,
            INT *perm, INT *iccsze , INT node_num );
 


#endif/*_RCM_H*/
