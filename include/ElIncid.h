#ifndef _ELINCIDH_
  #define _ELINCIDH_
/*...*/
  #include<Define.h>
  #include<Mesh.h>
  

  void mkElIncid(Memoria *m,Pnode *pn,Mesh *mesh);

  void nodeGrade(INT *restrict el       ,INT *restrict nincid     
                ,short *restrict nen     ,INT *maxGrade
                ,INT const nNode        ,INT const numel
                ,short const maxNo);

  void elmIncid(INT *restrict el    ,INT *restrict incid
               ,INT *restrict nincid,short *restrict nen  
               ,INT const nnode     ,INT const numel
               ,INT const maxGrade  ,short const maxNo);


#endif/*_ELINCIDH_*/
