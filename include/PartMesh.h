#ifndef _PARTMESH_H_
  #define _PARTMESH_H_
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

  void fPartMesh(Memoria *m       ,DOUBLE *RESTRICT x 
             ,INT *RESTRICT el   ,short *RESTRICT nen
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
  void divCoorXY(DOUBLE *RESTRICT x ,INT *RESTRICT el
                ,short *RESTRICT nen      
                ,INT const nNode   ,INT const nEl   
                ,INT *RESTRICT np  ,INT *RESTRICT eq   
                ,short const ndm   ,short const maxNo 
                ,short const nDiv  ,bool const fC);
/*...................................................................*/

/*... divisao da malha na coordenadas XYZ*/
  void divCoorXYZ(DOUBLE *RESTRICT coor   
              ,INT *RESTRICT el      ,short  *RESTRICT nen                          
              ,INT const nNode       ,INT const nEl   
              ,INT *RESTRICT np      ,INT *RESTRICT ep   
              ,short const ndm       ,short const maxNo
              ,short const nDiv      ,bool const fC);
/*...................................................................*/

/*... convertendo a malha para o formato CSR*/
  INT meshToCsrIa(short *RESTRICT nen,INT *RESTRICT eptr
                 ,INT const numel    );
  void meshToCsrJa(INT *RESTRICT el  ,short *RESTRICT nen
                ,INT *RESTRICT eptr,INT *RESTRICT eind
                ,INT const numel   ,short const maxNo);
/*...................................................................*/

/*... map*/
  unsigned short getMapViz(INT *RESTRICT ep        ,INT *RESTRICT  elLG
                ,short *RESTRICT mVizPart 
                ,INT const  lNel          ,INT const numelNov
                ,short const rank         ,short const nPrcs);
  
  unsigned short getMapVizNo(INT *RESTRICT ep,INT *RESTRICT  noLG
                  ,short *RESTRICT vizPart  ,bool *RESTRICT contViz
                  ,INT *RESTRICT nincid     ,INT *RESTRICT incid 
                  ,INT const  nNodeNov      ,INT const maxGrade              
                  ,short const rank         ,short const nPrcs);



  void getMapInterfaceEl(INT *RESTRICT ep
                        ,INT *RESTRICT nelcon   ,short *RESTRICT nFace
                        ,INT *RESTRICT elLG     ,bool *RESTRICT fEl
                        ,INT *nRcvs             ,INT *nSends 
                        ,INT *RESTRICT iaRcvs   ,INT *RESTRICT iaSends
                        ,short *RESTRICT vizPart,INT *RESTRICT fMap 
                        ,INT const lNel         ,INT const numelNov
                        ,short const rank       ,short const nPrcs
                        ,short const nVizPart   ,short const maxViz);
  
  void  getMapInterfaceNo(INT *RESTRICT ep
                         ,INT *RESTRICT noLG      ,bool *RESTRICT fNod
                         ,INT *nComNo             ,INT *RESTRICT iaComNo
                         ,INT *RESTRICT nincid    ,INT *RESTRICT incid         
                         ,short *RESTRICT vizPart ,INT *RESTRICT fMap 
                         ,INT const lNnode         
                         ,short const rank        ,short const nPrcs
                         ,short const nVizPart    ,short const maxGrade);


  void getNumberLocalMesh(INT *RESTRICT ep   ,INT *RESTRICT np                   
                       ,bool *RESTRICT fNode ,bool *RESTRICT fEl
                       ,INT *RESTRICT el     ,INT *RESTRICT nelcon 
                       ,short *RESTRICT nen  ,short *RESTRICT nFace
                       ,INT const nNode      ,INT const nEl       
                       ,short const maxNo    ,short const maxViz
                       ,INT *numelNov        ,INT *numelOv       
                       ,INT *nNodeNov        ,INT *nNodeOv
                       ,INT *nno1                                   
                       ,short const rank );

  void getMapElm(INT *RESTRICT ep                      
                ,INT *RESTRICT elLG    ,INT *RESTRICT elGL 
                ,bool *RESTRICT fEp
                ,INT *RESTRICT el      ,INT *RESTRICT nelcon 
                ,short *RESTRICT nFace
                ,INT const nEl         ,INT const lNel  
                ,short const maxViz    ,short const rank );

  void getMapNode(INT *RESTRICT ep        ,INT *RESTRICT np                      
                 ,INT *RESTRICT noLG      ,INT *RESTRICT noGL
                 ,bool *RESTRICT fNode
                 ,INT *RESTRICT el        ,INT *RESTRICT nelcon 
                 ,short *RESTRICT nen     ,short *RESTRICT nFace
                 ,INT const nNode         ,INT const nEl           
                 ,short const maxNo       ,short const maxViz  
                 ,INT const lNode         ,short const rank );

  void getLocalEl(INT *RESTRICT el   ,INT *RESTRICT lEl
                 ,INT *RESTRICT elLG ,INT *RESTRICT noGL
                 ,short *RESTRICT nen   
                 ,INT const nLel     ,short const maxNo);
    
  void dGetLocalV(DOUBLE *RESTRICT vG ,DOUBLE *RESTRICT vL
                 ,INT *RESTRICT mapLG
                 ,INT const nLin      ,INT const nCol);

  void sGetLocalV(short *RESTRICT vG ,short *RESTRICT vL
                 ,INT *RESTRICT mapLG
                 ,INT const nLin     ,INT const nCol);

  void getLocalAdj(INT *RESTRICT adj  ,INT *RESTRICT lAdj
                 ,INT *RESTRICT elLG ,INT *RESTRICT elGL
                 ,short *RESTRICT nen   
                 ,INT const nLel     ,short const maxViz);
/*...................................................................*/


#endif/*_PARTMESH_H_*/
