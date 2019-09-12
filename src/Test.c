#include<Define.h>
#include<Mesh.h>
#include<Math.h>

/*############### Teste do gradiente ####################*/
double func1(double x, double y, double z)
{
  double pi = PI,v1,v2,v3;
  v1 = cos(pi*x);
  v2 = sin(pi*y);
  v3 = z*z*z + 1;
  return v1+v2+v3;
}

void gradFunc1(double x, double y, double z,double *gradZ)
{
  double pi = PI;

  gradZ[0] = -pi * sin(pi*x);
  gradZ[1] =  pi * cos(pi*y);
  gradZ[2] = 3*z*z;

}  


void gradErro(Mesh *mesh)
{
  int i;
  double grad[3],x,y,z,err[3],errR=0.0, sumGrad=0.0;
  for (i = 0; i<mesh->numel; i++)
  {
    x = MAT2D(i, 0, mesh->elm.geom.cc, 3);
    y = MAT2D(i, 1, mesh->elm.geom.cc, 3);
    z = MAT2D(i, 2, mesh->elm.geom.cc, 3);
    
    gradFunc1(x,y,z,grad);

    err[0] = grad[0] - MAT2D(i, 0, mesh->elm.gradTemp, 3);
    err[1] = grad[1] - MAT2D(i, 1, mesh->elm.gradTemp, 3);
    err[2] = grad[2] - MAT2D(i, 2, mesh->elm.gradTemp, 3);
    errR   += err[0]*err[0]+ err[1]* err[1]+ err[2]*err[2];
    sumGrad = grad[0]* grad[0] + grad[1]*grad[1] + grad[2]*grad[2];
  }
  
  printf("Err A %e\n", sqrt(errR));
  printf("Err R %e\n", sqrt(errR)/ sqrt(sumGrad));

}

void grad(Mesh *mesh)
{
  int i;
  DOUBLE x, y, z;
  fprintf(fileLogDebug, "initialTemp\n");
  for (i = 0; i<mesh->numel; i++)
  {
    x = MAT2D(i, 0, mesh->elm.geom.cc, 3);
    y = MAT2D(i, 1, mesh->elm.geom.cc, 3);
    z = MAT2D(i, 2, mesh->elm.geom.cc, 3);
    fprintf(fileLogDebug, "%9d %12.6lf\n", i + 1, func1(x, y, z));
  }
  fprintf(fileLogDebug, "endInitialTemp\nreturn\n");
}
/*#######################################################*/