#!/bin/sh
INPUTRCM='retangulo.dat retangulo_rcm.dat'
INPUTVTK='retangulo.dat'

#regime estacionario difusao
INPUTDIF1='ex3_60_0_gglc.dat    ex3_60_3_gglc.dat
           ex3_60_0_ggln.dat    ex3_60_3_ggln.dat
           ex3_60_0_lSquare.dat ex3_60_3_lSquare.dat'

INPUTDIF2='ex2_1_gglc.dat       ex2_2_gglc.dat
           ex2_1_ggln.dat       ex2_2_ggln.dat
           ex2_1_lSquare.dat    ex2_2_lSquare.dat' 

INPUTDIF3='ex3_60_0_gglc.dat    ex3_60_3_gglc.dat
           ex3_60_0_ggln.dat    ex3_60_3_ggln.dat
           ex3_60_0_lSquare.dat ex3_60_3_lSquare.dat'


INPUTDIF4='ex4_1_ls.dat    ex4_2_ls.dat 
           ex4_1_ggln.dat  ex4_2_ggln.dat  
           ex4_1_gglc.dat  ex4_2_gglc.dat' 

INPUTDIF5='ex5_1_ls.dat    ex5_2_ls.dat 
           ex5_1_ggln.dat  ex5_2_ggln.dat  
           ex5_1_gglc.dat  ex5_2_gglc.dat' 

INPUTDIF3DEX1='ex1_3d_hexa_0_lsqr.dat ex1_3d_hexa_1_lsqr.dat  
               ex1_3d_hexa_0_ls.dat   ex1_3d_hexa_1_ls.dat    
               ex1_3d_hexa_0_ggln.dat ex1_3d_hexa_1_ggln.dat  
               ex1_3d_hexa_0_gglc.dat ex1_3d_hexa_1_gglc.dat  
               ex1_3d_tetra_0_ls.dat ex1_3d_tetra_1_ls.dat  
               ex1_3d_tetra_0_lsqr.dat ex1_3d_tetra_1_lsqr.dat   
               ex1_3d_tetra_0_ggln.dat ex1_3d_tetra_1_ggln.dat   
               ex1_3d_tetra_0_gglc.dat ex1_3d_tetra_1_gglc.dat'  

INPUTDIF3DEX2='ex2_3D_0_ls.dat    ex2_3D_1_ls.dat  
               ex2_3D_0_ggln.dat  ex2_3D_1_ggln.dat  
               ex2_3D_0_gglc.dat  ex2_3D_1_gglc.dat 
               ex2_3D_0_lsqr.dat  ex2_3D_1_lsqr.dat'

INPUTDIF3DEX2NONORT='ex2_3d_tetra_0_ggln.dat
                     ex2_3d_tetra_0_gglc.dat 
                     ex2_3d_tetra_0_ls.dat   
                     ex2_3d_tetra_0_lsqr.dat' 

INPUTDIF3DEX3='ex3_3d_hexa_0_ls.dat    ex3_3d_hexa_1_ls.dat
               ex3_3d_hexa_0_ggln.dat  ex3_3d_hexa_1_ggln.dat 
               ex3_3d_hexa_0_gglc.dat  ex3_3d_hexa_1_gglc.dat
               ex3_3d_hexa_0_lsqr.dat    ex3_3d_hexa_1_lsqr.dat
               ex3_3d_tetra_0_ggln.dat ex3_3d_tetra_1_ggln.dat 
               ex3_3d_tetra_0_gglc.dat ex3_3d_tetra_1_gglc.dat 
               ex3_3d_tetra_0_ls.dat   ex3_3d_tetra_1_ls.dat   
               ex3_3d_tetra_0_lsqr.dat   ex3_3d_tetra_1_lsqr.dat'  

INPUTDIF3DEX4='ex4_3d_tetra_0_gglc.dat   ex4_3d_tetra_1_gglc.dat
               ex4_3d_tetra_0_ggln.dat   ex4_3d_tetra_1_ggln.dat
               ex4_3d_tetra_0_lsqr.dat   ex4_3d_tetra_1_lsqr.dat
               ex4_3d_tetra_0_ls.dat     ex4_3d_tetra_1_ls.dat
               ex4_3d_hexa_0_ls.dat      ex4_3d_hexa_1_ls.dat
               ex4_3d_hexa_0_lsqr.dat    ex4_3d_hexa_1_lsqr.dat
               ex4_3d_hexa_0_ggln.dat    ex4_3d_hexa_1_ggln.dat
               ex4_3d_hexa_0_gglc.dat    ex4_3d_hexa_1_gglc.dat'     

# transiente difusao
INPUTDIF2TIME='ex2_1_gglc.dat ex2_1_gglc_e.dat'

INPUTDIF4TIME='ex4_1_gglc.dat ex4_1_gglc_e.dat' 

INPUTDIF5TIME='ex5_1_gglc.dat  ex5_1_gglc_e.dat' 

INPUTDIF3DEX2TIME='ex2_3D_1_gglc.dat ex2_3D_1_gglc_e.dat'

# transporte estacionario
INPUTTRANS2D='ex2_1_transporte_gglc.dat ex2_1_transporte_ggln.dat ex2_1_transporte_ls.dat' 

INPUTTRANS3D='ex3_3D_1_transporte.dat ex3_3D_1_transporte_tetra.dat' 


NAMEBIN=mfvCell_gnu_mpi_O3

#difusao
ex1=true 
ex2=true 
ex3=false
ex4=true
ex5=true
ex1_3D=true 
ex2_3D=true
ex2_3D_NONORT=true 
ex3_3D=true 
ex4_3D=true  

# dinamico
ex2t=true
ex4t=true
ex5t=true
ex2_3Dt=true

#transporte
ex2Transporte2D=true
ex3Transporte3D=true

#gerando o executavel
echo `make clean` > /dev/null 
echo `make` > /dev/null       
mv bin/$NAMEBIN test/binTest
cp test/input/*.dat test/binTest
#rodando os teste COO
#test/binTest/run_coo.sh "$INPUTRCM" "$NAMEBIN" 

#rodando os teste VTK
#cp test/input/vtk/*.dat test/binTest
#test/binTest/run_vtk.sh "$INPUTVTK" "$NAMEBIN" 
#if [ $? == 1 ];then
#  exit 1
#fi

#rodando os teste de transporte
if  $ex2Transporte2D ; then
  cp test/input/transport/orthogonal/2D/gaussGreenCell/*.dat test/binTest
  cp test/input/transport/orthogonal/2D/gaussGreenNode/*.dat test/binTest
  cp test/input/transport/orthogonal/2D/LeastSquare/*.dat test/binTest
#
  test/binTest/run_trans_ex2_2D.sh "$INPUTTRANS2D" "$NAMEBIN" 
  if [ $? == 1 ]; then
    exit 1
  fi
fi

if  $ex3Transporte3D ; then
  cp test/input/transport/orthogonal/3D/gaussGreenCell/*.dat test/binTest
  cp test/input/transport/non_orthogonal/3D/gaussGreenCell/*.dat test/binTest
#
  test/binTest/run_trans_ex3_3D.sh "$INPUTTRANS3D" "$NAMEBIN" 
  if [ $? == 1 ]; then
    exit 1
  fi
fi

#rodando os teste de difusao
if  $ex2 ; then
  cp test/input/dif/orthogonal/2D/gaussGreenCell/*.dat test/binTest
  cp test/input/dif/orthogonal/2D/gaussGreenNode/*.dat test/binTest
  cp test/input/dif/orthogonal/2D/LeastSquare/*.dat test/binTest
#
  test/binTest/run_dif_exato.sh "$INPUTDIF2" "$NAMEBIN" 
  if [ $? == 1 ]; then
    exit 1
  fi
fi

if  $ex3 ; then
  cp test/input/dif/non_orthogonal/2D/60/gaussGreenCell/*.dat test/binTest
  cp test/input/dif/non_orthogonal/2D/60/gaussGreenNode/*.dat test/binTest
  cp test/input/dif/non_orthogonal/2D/60/LeastSquare/*.dat test/binTest
  test/binTest/run_dif.sh "$INPUTDIF3" "$NAMEBIN" 

  if [ $? == 1 ];then
    exit 1
  fi
fi

if $ex4 ; then
  cp test/input/dif/orthogonal/2D/gaussGreenCell/*.dat test/binTest
  cp test/input/dif/orthogonal/2D/gaussGreenNode/*.dat test/binTest
  cp test/input/dif/orthogonal/2D/LeastSquare/*.dat test/binTest
  test/binTest/run_dif_exato_Ex4.sh "$INPUTDIF4" "$NAMEBIN" 
  if [ $? == 1 ];then
    exit 1
  fi
fi

#
if $ex5 ; then
  cp test/input/dif/orthogonal/2D/gaussGreenCell/*.dat test/binTest
  cp test/input/dif/orthogonal/2D/gaussGreenNode/*.dat test/binTest
  cp test/input/dif/orthogonal/2D/LeastSquare/*.dat test/binTest
  test/binTest/run_dif_exato_Ex5.sh "$INPUTDIF5" "$NAMEBIN" 
  if [ $? == 1 ];then
    exit 1
  fi
fi

# difusao 3D
if  $ex2_3D ; then
  cp test/input/dif/orthogonal/3D/GreenGaussCell/*.dat test/binTest
  cp test/input/dif/orthogonal/3D/GreenGaussNode/*.dat test/binTest
  cp test/input/dif/orthogonal/3D/LeastSquare/*.dat test/binTest
  cp test/input/dif/orthogonal/3D/LeastSquareQR/*.dat test/binTest
  test/binTest/run_dif_exato_3D_Ex2.sh "$INPUTDIF3DEX2" "$NAMEBIN" 
  if [ $? == 1 ];then
    exit 1
  fi
fi

# difusao 3D
if  $ex3_3D ; then 
  for f in test/input/dif/non_orthogonal/3D/*.tar.gz ; do tar -xvzf $f -C test/binTest ; done
  cp test/input/dif/non_orthogonal/3D/GreenGaussCell/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/GreenGaussNode/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/LeastSquare/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/LeastSquareQR/*.dat test/binTest
  test/binTest/run_dif_exato_3D_Ex3.sh "$INPUTDIF3DEX3" "$NAMEBIN" 
  if [ $? == 1 ];then
    exit 1
  fi
fi

# difusao 3D
if  $ex2_3D_NONORT ; then 
  for f in test/input/dif/non_orthogonal/3D/*.tar.gz ; do tar -xvzf $f -C test/binTest ; done
  cp test/input/dif/non_orthogonal/3D/GreenGaussCell/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/GreenGaussNode/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/LeastSquare/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/LeastSquareQR/*.dat test/binTest
  test/binTest/run_dif_exato_3D_Ex2.sh "$INPUTDIF3DEX2NONORT" "$NAMEBIN" 
  if [ $? == 1 ];then
    exit 1
  fi
fi

# difusao 3D
if  $ex1_3D ; then 
  for f in test/input/dif/orthogonal/3D/*.tar.gz ; do tar -xvzf $f -C test/binTest ; done
  for f in test/input/dif/non_orthogonal/3D/*.tar.gz ; do tar -xvzf $f -C test/binTest ; done
  cp test/input/dif/orthogonal/3D/GreenGaussCell/*.dat test/binTest
  cp test/input/dif/orthogonal/3D/GreenGaussNode/*.dat test/binTest
  cp test/input/dif/orthogonal/3D/LeastSquare/*.dat test/binTest
  cp test/input/dif/orthogonal/3D/LeastSquareQR/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/GreenGaussCell/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/GreenGaussNode/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/LeastSquare/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/LeastSquareQR/*.dat test/binTest
  test/binTest/run_dif_exato_3D_Ex1.sh "$INPUTDIF3DEX1" "$NAMEBIN" 
  if [ $? == 1 ];then
    exit 1
  fi
fi

# difusao 3D
if  $ex1_3D ; then 
  for f in test/input/dif/orthogonal/3D/*.tar.gz ; do tar -xvzf $f -C test/binTest ; done
  for f in test/input/dif/non_orthogonal/3D/*.tar.gz ; do tar -xvzf $f -C test/binTest ; done
  cp test/input/dif/orthogonal/3D/GreenGaussCell/*.dat test/binTest
  cp test/input/dif/orthogonal/3D/GreenGaussNode/*.dat test/binTest
  cp test/input/dif/orthogonal/3D/LeastSquare/*.dat test/binTest
  cp test/input/dif/orthogonal/3D/LeastSquareQR/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/GreenGaussCell/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/GreenGaussNode/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/LeastSquare/*.dat test/binTest
  cp test/input/dif/non_orthogonal/3D/LeastSquareQR/*.dat test/binTest
  test/binTest/run_dif_exato_3D_Ex4.sh "$INPUTDIF3DEX4" "$NAMEBIN" 
  if [ $? == 1 ];then
    exit 1
  fi
fi

if  $ex2t ; then
  cp test/input/dif/orthogonal/2D_temp/gaussGreenCell/*.dat test/binTest
# cp test/input/dif/orthogonal/2D/gaussGreenNode/*.dat test/binTest
# cp test/input/dif/orthogonal/2D/LeastSquare/*.dat test/binTest
#
  test/binTest/run_dif_exato_temp.sh "$INPUTDIF2TIME" "$NAMEBIN" 
  if [ $? == 1 ]; then
    exit 1
  fi
fi

if  $ex4t ; then
  cp test/input/dif/orthogonal/2D_temp/gaussGreenCell/*.dat test/binTest
# cp test/input/dif/orthogonal/2D/gaussGreenNode/*.dat test/binTest
# cp test/input/dif/orthogonal/2D/LeastSquare/*.dat test/binTest
#
  test/binTest/run_dif_exato_Ex4_temp.sh "$INPUTDIF4TIME" "$NAMEBIN" 
  if [ $? == 1 ]; then
    exit 1
  fi
fi

if  $ex5t ; then
  cp test/input/dif/orthogonal/2D_temp/gaussGreenCell/*.dat test/binTest
# cp test/input/dif/orthogonal/2D/gaussGreenNode/*.dat test/binTest
# cp test/input/dif/orthogonal/2D/LeastSquare/*.dat test/binTest
#
  test/binTest/run_dif_exato_Ex5_temp.sh "$INPUTDIF5TIME" "$NAMEBIN" 
  if [ $? == 1 ]; then
    exit 1
  fi
fi

# difusao 3D
if  $ex2_3Dt ; then
  cp test/input/dif/orthogonal/3D_temp/GreenGaussCell/*.dat test/binTest
  test/binTest/run_dif_exato_3D_Ex2_temp.sh "$INPUTDIF3DEX2TIME" "$NAMEBIN" 
  if [ $? == 1 ];then
    exit 1
  fi
fi


DIR="test/binTest"
rm  $DIR/*.dat $DIR/*.vtk  $DIR/*.txt   $DIR/*.csv $DIR/$NAMEBIN 
#rm  $DIR/*.mtx 

exit 0
