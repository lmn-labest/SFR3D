#ifndef _ELINCID_H_
  #define _ELINCID_H_
/*...*/
  #include<Define.h>
  #include<Mesh.h>
  

  void mkElIncid(Memoria *m,Pnode *pn,Mesh *mesh);

  void nodeGrade(INT *RESTRICT el       ,INT *RESTRICT nincid     
                ,short *RESTRICT nen     ,INT *maxGrade
                ,INT const nNode        ,INT const numel
                ,short const maxNo);

  void elmIncid(INT *RESTRICT el    ,INT *RESTRICT incid
               ,INT *RESTRICT nincid,short *RESTRICT nen  
               ,INT const nnode     ,INT const numel
               ,INT const maxGrade  ,short const maxNo);


#endif/*_ELINCID_H_*/
