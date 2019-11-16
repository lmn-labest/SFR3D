#include<SaveLoad.h>

/*********************************************************************
 * Data de criacao    : 07/02/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * wSave :                                                           * 
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void wSave(PropVarFluid *prop           ,Turbulence *tModel
         ,ThermoDynamic *thDynamic
         ,Mesh *mesh              ,Temporal *ddt
         ,Mean *media
         ,char *preName           ,char *nameOut) {

  short k,j,ndm=mesh->ndm,nD;
  INT i;
  FILE *file;

/*...*/
  fName(preName,ddt->timeStep,0,24,nameOut);
  file = openFile(nameOut,"w");
/*...................................................................*/

/*...*/
  fprintf(file,"/dt/\n");
  fprintf(file,"%d %.15e %.15e %.15e %.15e\n",ddt->timeStep
                                             ,ddt->t
                                             ,ddt->dt[0] 
                                             ,ddt->dt[1]
                                             ,ddt->dt[2]);
/*.....................................................................*/

/*... velocidade*/
  fprintf(file,"/Vel/\n");
  for (i = 0; i < mesh->numel; i++) {
    for (j = 0; j < mesh->ndm; j++)
      fprintf(file,"%.15e ",MAT2D(i,j,mesh->elm.vel0,ndm));  
    for (j = 0; j < mesh->ndm; j++)
      fprintf(file,"%.15e ",MAT2D(i,j,mesh->elm.vel,ndm));  
    fprintf(file,"\n");  
  }
/*.....................................................................*/

/*... gradVel*/
  fprintf(file,"/gradVel/\n");
  for (i = 0; i < mesh->numel; i++) {
    for (j = 0; j < ndm; j++)
      for (k = 0; k < ndm; k++)
        fprintf(file,"%.15e ",MAT3D(i,j,k,mesh->elm.gradVel,ndm,ndm));  
    fprintf(file,"\n");  
  }
/*.....................................................................*/

/*... pressao*/
  fprintf(file,"/pres/\n");
  for (i = 0; i < mesh->numel; i++)
    fprintf(file,"%.15e %.15e\n",mesh->elm.pressure0[i]
                               ,mesh->elm.pressure[i]);
/*.....................................................................*/

/*... gradPres*/
  fprintf(file,"/gradPres/\n");
  for (i = 0; i < mesh->numel; i++) {
    for (j = 0; j < ndm; j++)
        fprintf(file,"%.15e ",MAT2D(i,j,mesh->elm.gradVel,ndm));  
    fprintf(file,"\n");  
  }
/*.....................................................................*/

/*... temperatura0*/
  fprintf(file,"/temp/\n");
  for (i = 0; i < mesh->numel; i++){
    fprintf(file,"%.15e %.15e\n",mesh->elm.temp0[i],mesh->elm.temp[i]);
  }
/*.....................................................................*/

/*... gradTemp*/
  fprintf(file,"/gradTemp/\n");
  for (i = 0; i < mesh->numel; i++) {
    for (j = 0; j < ndm; j++)
        fprintf(file,"%.15e ",MAT2D(i,j,mesh->elm.gradTemp,ndm));  
    fprintf(file,"\n");  
  }
/*.....................................................................*/

/*... energy0*/
  fprintf(file,"/energy/\n");
  for (i = 0; i < mesh->numel; i++)
    fprintf(file,"%.15e %.15e\n",mesh->elm.energy0[i]
                                ,mesh->elm.energy[i]);
/*.....................................................................*/

/*... gradEnergy*/
  fprintf(file,"/gradEnergy/\n");
  for (i = 0; i < mesh->numel; i++) {
    for (j = 0; j < ndm; j++)
        fprintf(file,"%.15e ",MAT2D(i,j,mesh->elm.gradEnergy,ndm));  
    fprintf(file,"\n");  
  }
/*.....................................................................*/

/*... desnsity*/
  fprintf(file,"/densityFluid/\n");
  for (i = 0; i < mesh->numel; i++){
    fprintf(file,"%.15e %.15e %.15e\n"
            ,mesh->elm.densityFluid.t[i]
            ,mesh->elm.densityFluid.t0[i]
            ,mesh->elm.densityFluid.t00[i]);
     fprintf(file,"\n");
  }
/*.....................................................................*/

/*... specificHeat*/
  fprintf(file,"/specificHeat/\n");
  for (i = 0; i < mesh->numel; i++){
    fprintf(file,"%.15e %.15e %.15e\n"
            ,mesh->elm.specificHeat.t[i]
            ,mesh->elm.specificHeat.t0[i]
            ,mesh->elm.specificHeat.t00[i]);
  }
/*.....................................................................*/

/*... condutividade*/
  fprintf(file,"/tConductivity/\n");
  for (i = 0; i < mesh->numel; i++){
    fprintf(file,"%.15e\n",mesh->elm.tConductivity[i]);
  }
/*.....................................................................*/

/*... viscosidade*/
  fprintf(file,"/dViscosity/\n");
  for (i = 0; i < mesh->numel; i++){
    fprintf(file,"%.15e\n",mesh->elm.dViscosity[i]);
  }
/*.....................................................................*/

/*... eddyViscosity*/
  if(tModel->fTurb){
    fprintf(file,"/eddyViscosity/\n");
    for (i = 0; i < mesh->numel; i++){
      fprintf(file,"%.15e\n",mesh->elm.eddyViscosity[i]);
    }
  }
/*.....................................................................*/

/*... cDynamic*/
  if(tModel->fDynamic){
    fprintf(file,"/cDynamic/\n");
    for (i = 0; i < mesh->numel; i++){
      fprintf(file,"%.15e %.15e\n",MAT2D(i,0,mesh->elm.cd,2)
                                  ,MAT2D(i,1,mesh->elm.cd,2));
    }
  }
/*.....................................................................*/

/*... stressR*/
  if(tModel->typeLes == LESSTRUMODEL){
    nD = mesh->ntn;
    fprintf(file,"/stressR/\n");
    for (i = 0; i < mesh->numel; i++){
      for (j = 0; j < nD; j++)
        fprintf(file," %.15e ",MAT2D(i,j,mesh->elm.stressR,nD));
      fprintf(file,"\n");
    }
  }
/*.....................................................................*/

/*... condutividade*/
  if(media->fMedia){
    fprintf(file,"/media/ %e %d %d %d\n",media->t0
                                        ,media->startSample
                                        ,media->endSample
                                        ,mesh->numel);
    if(media->fVel){
      fprintf(file,"/mVel/\n");
      for (i = 0; i < mesh->numel; i++){
        for (j = 0; j < mesh->ndm; j++)
          fprintf(file,"%.15e ",MAT2D(i,j,media->mVel,ndm));  
        for (j = 0; j < mesh->ndm; j++)
          fprintf(file,"%.15e ",MAT2D(i,j,media->sVel,ndm));  
        fprintf(file,"\n");
      }
    }
  }
/*..................................................................*/
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 07/02/2018                                   *
 * Data de modificaco : 00/00/0000                                   *
 *-------------------------------------------------------------------*
 * load :                                                            * 
 *-------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------*
 * Parametros de saida:                                              *
 *-------------------------------------------------------------------*
 *-------------------------------------------------------------------* 
 * OBS:                                                              * 
 *-------------------------------------------------------------------* 
 *********************************************************************/
void load(PropVarFluid *prop           ,Turbulence *tModel
         ,ThermoDynamic *thDynamic
         ,Mesh *mesh              ,Temporal *ddt
         ,Mean *media
         ,FILE *fileIn) {
  
  char word[WORD_SIZE], str[1024];
  short k,j,ndm=mesh->ndm,nD;
  INT i;
  FILE *file;

  readMacro(fileIn, word, false);

/*...*/  
  file = openFile(word,"r");
/*...................................................................*/

/*...*/
  fscanf(file,"%s",str);
  fscanf(file,"%d %lf %lf %lf %lf",&ddt->timeStep
                                  ,&ddt->t 
                                  ,&ddt->dt[0] 
                                  ,&ddt->dt[1]
                                  ,&ddt->dt[2]);
/*.....................................................................*/

/*... velocidade*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++) {
    for (j = 0; j < mesh->ndm; j++)
      fscanf(file,"%lf",&MAT2D(i,j,mesh->elm.vel0,ndm));  
    for (j = 0; j < mesh->ndm; j++)
      fscanf(file,"%lf",&MAT2D(i,j,mesh->elm.vel,ndm));  
  }
/*.....................................................................*/

/*... gradVel*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++) {
    for (j = 0; j < ndm; j++)
      for (k = 0; k < ndm; k++)
        fscanf(file,"%lf",&MAT3D(i,j,k,mesh->elm.gradVel,ndm,ndm));  
  }
/*.....................................................................*/

/*... pressao*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++)
    fscanf(file,"%lf %lf",&mesh->elm.pressure0[i]
                         ,&mesh->elm.pressure[i]);
/*.....................................................................*/

/*... gradPres*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++) {
    for (j = 0; j < ndm; j++)
        fscanf(file,"%lf",&MAT2D(i,j,mesh->elm.gradPres,ndm));  
  }
/*.....................................................................*/

/*... temperatura0*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++){
    fscanf(file,"%lf %lf",&mesh->elm.temp0[i]
                         ,&mesh->elm.temp[i]);
  }
/*.....................................................................*/

/*... gradTemp*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++) {
    for (j = 0; j < ndm; j++)
      fscanf(file,"%lf",&MAT2D(i,j,mesh->elm.gradTemp,ndm));  
  }
/*.....................................................................*/

/*... energy0*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++)
    fscanf(file,"%lf %lf",&mesh->elm.energy0[i]
                         ,&mesh->elm.energy[i]);
/*.....................................................................*/

/*... gradEnergy*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++) {
    for (j = 0; j < ndm; j++)
      fscanf(file,"%lf",&MAT2D(i,j,mesh->elm.gradEnergy,ndm));  
  }
/*.....................................................................*/

/*... desnsity*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++){
    fscanf(file,"%lf %lf %lf"
               ,&mesh->elm.densityFluid.t[i]
               ,&mesh->elm.densityFluid.t0[i]
               ,&mesh->elm.densityFluid.t00[i]);
  }
/*.....................................................................*/

/*... specificHeat*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++){
    fscanf(file,"%lf %lf %lf"
               ,&mesh->elm.specificHeat.t[i]
               ,&mesh->elm.specificHeat.t0[i]
               ,&mesh->elm.specificHeat.t00[i]);
  }
/*.....................................................................*/

/*... condutividade*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++){
    fscanf(file,"%lf",&mesh->elm.tConductivity[i]);
  }
/*.....................................................................*/

/*... condutividade*/
  fscanf(file,"%s",str);
  for (i = 0; i < mesh->numel; i++){
    fscanf(file,"%lf",&mesh->elm.dViscosity[i]);
  }
/*.....................................................................*/

/*... eddyViscosity*/
  if(tModel->fTurb){
    fscanf(file,"%s",str);
    for (i = 0; i < mesh->numel; i++){
      fscanf(file,"%lf\n",&mesh->elm.eddyViscosity[i]);
    }
  }
/*.....................................................................*/

/*... cDynamic*/
  if(tModel->fDynamic){
    fscanf(file,"%s",str);
    for (i = 0; i < mesh->numel; i++){
      fscanf(file,"%lf %lf",&MAT2D(i,0,mesh->elm.cd,2)
                            ,&MAT2D(i,1,mesh->elm.cd,2));
    }
  }
/*.....................................................................*/

/*... stressR*/
  if(tModel->typeLes == LESSTRUMODEL){
    nD = mesh->ntn;
    fscanf(file,"%s",str);
    for (i = 0; i < mesh->numel; i++){
      for (j = 0; j < nD; j++)
        fscanf(file,"%lf",&MAT2D(i,j,mesh->elm.stressR,nD));
    }
  }
/*.....................................................................*/

/*... condutividade*/
  if(media->fMedia){
    fscanf(file,"%s %lf %d %d\n",str
                                    ,&media->t0
                                    ,&media->startSample
                                    ,&media->endSample);
    if(media->fVel){
      fscanf(file,"%s",str);
      for (i = 0; i < mesh->numel; i++){
        for (j = 0; j < mesh->ndm; j++)
          fscanf(file,"%lf",&MAT2D(i,j,media->mVel,ndm));  
        for (j = 0; j < mesh->ndm; j++)
          fscanf(file,"%lf",&MAT2D(i,j,media->sVel,ndm)); 
       }
    }
  }
/*.....................................................................*/


}
