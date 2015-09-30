#ifndef _PARTMESH_
  #define _PARTMESH_
/*...*/
  #include<Mesh.h>
  #include<Memoria.h>
  #include<HccaTime.h>
  #include<HccaStdBool.h>
  #include<Define.h>
  #include<File.h>
  #ifdef _METIS_
    #include<metis/metis.h>
  #endif

  void partMesh(Memoria *m      
             ,DOUBLE *restrict x ,INT *restrict el
             ,short *restrict nen
             ,INT const nNode    ,INT const nNumel
             ,PartMesh *pMesh  
             ,short const ndm    ,short const maxViz
             ,short const nDiv);

/*... impressao do mapa*/
  void printMap(PartMesh pMesh
               ,INT const nNode ,INT const nEl
               ,short const myId
               ,char *nameOut   ,FILE *f);
/*...................................................................*/

/*... divisao da malha na coordenadas XY*/
  void divCoorXY(DOUBLE *restrict x ,INT *restrict el
                ,short *restrict nen      
                ,INT const nNode   ,INT const nEl   
                ,INT *restrict np  ,INT *restrict eq   
                ,short const ndm   ,short const maxNo 
                ,short const nDiv  ,bool const fC);
/*...................................................................*/

/*... convertendo a malha para o formato CSR*/
  INT meshToCsrIa(short *restrict nen,INT *restrict eptr
                 ,INT const numel    );
  void meshToCsrJa(INT *restrict el  ,short *restrict nen
                ,INT *restrict eptr,INT *restrict eind
                ,INT const numel   ,short const maxNo);
/*...................................................................*/

/*... map*/
  unsigned short getMapViz(INT *restrict ep        ,INT *restrict  elLG
                ,short *restrict mVizPart 
                ,INT const  lNel          ,INT const numelNov
                ,short const rank         ,short const nPrcs);
  
  unsigned short getMapVizNo(INT *restrict ep,INT *restrict  noLG
                  ,short *restrict vizPart  ,bool *restrict contViz
                  ,INT *restrict nincid     ,INT *restrict incid 
                  ,INT const  nNodeNov      ,INT const maxGrade              
                  ,short const rank         ,short const nPrcs);



  void getMapInterfaceEl(INT *restrict ep
                        ,INT *restrict nelcon   ,short *restrict nFace
                        ,INT *restrict elLG     ,bool *restrict fEl
                        ,INT *nRcvs             ,INT *nSends 
                        ,INT *restrict iaRcvs   ,INT *restrict iaSends
                        ,short *restrict vizPart,INT *restrict fMap 
                        ,INT const lNel         ,INT const numelNov
                        ,short const rank       ,short const nPrcs
                        ,short const nVizPart   ,short const maxViz);
  
  void  getMapInterfaceNo(INT *restrict ep
                         ,INT *restrict noLG      ,bool *restrict fNod
                         ,INT *nComNo             ,INT *restrict iaComNo
                         ,INT *restrict nincid    ,INT *restrict incid         
                         ,short *restrict vizPart ,INT *restrict fMap 
                         ,INT const lNnode         
                         ,short const rank        ,short const nPrcs
                         ,short const nVizPart    ,short const maxGrade);


  void getNumberLocalMesh(INT *restrict ep   ,INT *restrict np                   
                       ,bool *restrict fNode ,bool *restrict fEl
                       ,INT *restrict el     ,INT *restrict nelcon 
                       ,short *restrict nen  ,short *restrict nFace
                       ,INT const nNode      ,INT const nEl       
                       ,short const maxNo    ,short const maxViz
                       ,INT *numelNov        ,INT *numelOv       
                       ,INT *nNodeNov        ,INT *nNodeOv
                       ,INT *nno1                                   
                       ,short const rank );

  void getMapElm(INT *restrict ep                      
                ,INT *restrict elLG    ,INT *restrict elGL 
                ,bool *restrict fEp
                ,INT *restrict el      ,INT *restrict nelcon 
                ,short *restrict nFace
                ,INT const nEl         ,INT const lNel  
                ,short const maxViz    ,short const rank );

  void getMapNode(INT *restrict ep        ,INT *restrict np                      
                 ,INT *restrict noLG      ,INT *restrict noGL
                 ,bool *restrict fNode
                 ,INT *restrict el        ,INT *restrict nelcon 
                 ,short *restrict nen     ,short *restrict nFace
                 ,INT const nNode         ,INT const nEl           
                 ,short const maxNo       ,short const maxViz  
                 ,INT const lNode         ,short const rank );

  void getLocalEl(INT *restrict el   ,INT *restrict lEl
                 ,INT *restrict elLG ,INT *restrict noGL
                 ,short *restrict nen   
                 ,INT const nLel     ,short const maxNo);
    
  void dGetLocalV(DOUBLE *restrict vG ,DOUBLE *restrict vL
                 ,INT *restrict mapLG
                 ,INT const nLin      ,INT const nCol);

  void sGetLocalV(short *restrict vG ,short *restrict vL
                 ,INT *restrict mapLG
                 ,INT const nLin     ,INT const nCol);

  void getLocalAdj(INT *restrict adj  ,INT *restrict lAdj
                 ,INT *restrict elLG ,INT *restrict elGL
                 ,short *restrict nen   
                 ,INT const nLel     ,short const maxViz);
/*...................................................................*/


#endif/*_PARTMESH_*/
