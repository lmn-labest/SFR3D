#ifndef _RCM_H  
  #define _RCM_H
  #define ERRO_RCM fprintf(stderr,"\nrcm - fatal error!\n")
  #define vectorPlusOne(v,n,i)  for(i=0;i<n;i++) v[i] = v[i] + 1; 
  #define vectorMinusOne(v,n,i) for(i=0;i<n;i++) v[i] = v[i] - 1; 
  #include<stdio.h>
  #include<stdlib.h>
/*...*/
  void genrcm(long node_num, long *adj_row, long *adj
             ,long *perm);

  void root_find( long *root , long *adj_row  , long *adj
                , char *mask , long *level_num, long *level_row
                , long *level, long node_num );

  void level_set( long root  ,long *adj_row, long *adj
                , char *mask ,long *level_num, long *level_row
                , long *level,long node_num );
  
  void degree ( long root, long *adj_row, long *adj, char *mask,
                long *deg, long *iccsze, long *ls, long node_num );
  
  void ivec_reverse ( long n, long *a );

  void rcm ( long root, long *adj_row, long *adj, char *mask,
            long *perm, long *iccsze , long node_num );
 


#endif/*_RCM_H*/
