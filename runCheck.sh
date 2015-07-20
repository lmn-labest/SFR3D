#!/bin/sh
INPUTRCM='retangulo.dat retangulo_rcm.dat'
INPUTVTK='retangulo.dat'
INPUTDIF1='ex3_60_0_gglc.dat    ex3_60_3_gglc.dat
          ex3_60_0_ggln.dat    ex3_60_3_ggln.dat
          ex3_60_0_lSquare.dat ex3_60_3_lSquare.dat'

INPUTDIF2='ex2_1_gglc.dat       ex2_2_gglc.dat
           ex2_1_ggln.dat       ex2_2_ggln.dat
           ex2_1_lSquare.dat    ex2_2_lSquare.dat' 


INPUTDIF4='ex4_1_ls.dat    ex4_2_ls.dat 
           ex4_1_ggln.dat  ex4_2_ggln.dat  
           ex4_1_gglc.dat  ex4_2_gglc.dat' 

INPUTDIF5='ex5_1_ls.dat    ex5_2_ls.dat 
           ex5_1_ggln.dat  ex5_2_ggln.dat  
           ex5_1_gglc.dat  ex5_2_gglc.dat' 

INPUTDIF3DEX2='ex2_3D_0_ls.dat    ex2_3D_1_ls.dat  
               ex2_3D_0_ggln.dat  ex2_3D_1_ggln.dat  
               ex2_3D_0_gglc.dat  ex2_3D_1_gglc.dat' 

INPUTDIF3DEX3='  ex3_3d_hexa_0_ls.dat   ex3_3d_hexa_1_ls.dat ex3_3d_hexa_2_ls.dat  
               ex3_3d_hexa_0_ggln.dat ex3_3d_hexa_1_ggln.dat ex3_3d_hexa_2_ggln.dat 
               ex3_3d_hexa_0_gglc.dat ex3_3d_hexa_1_gglc.dat ex3_3d_hexa_2_gglc.dat'  


NAMEBIN=mfvCell_gnu_O3

#gerando o executavel
echo `make clean` > /dev/null 
echo `make` > /dev/null       
mv bin/$NAMEBIN test/binTest
cp test/input/*.dat test/binTest
#rodando os teste COO
#test/binTest/run_coo.sh "$INPUTRCM" "$NAMEBIN" 

#rodando os teste VTK
cp test/input/vtk/*.dat test/binTest
test/binTest/run_vtk.sh "$INPUTVTK" "$NAMEBIN" 
if [ $? == 1 ];then
  exit 1
fi

#rodando os teste difusao
#cp test/input/dif/non_orthogonal/2D/60/gaussGreenCell/*.dat test/binTest
#cp test/input/dif/non_orthogonal/2D/60/gaussGreenNode/*.dat test/binTest
#cp test/input/dif/non_orthogonal/2D/60/LeastSquare/*.dat test/binTest
#test/binTest/run_dif.sh "$INPUTDIF1" "$NAMEBIN" 

#if [ $? == 1 ];then
#  exit 1
#fi

cp test/input/dif/orthogonal/2D/gaussGreenCell/*.dat test/binTest
cp test/input/dif/orthogonal/2D/gaussGreenNode/*.dat test/binTest
cp test/input/dif/orthogonal/2D/LeastSquare/*.dat test/binTest
#
test/binTest/run_dif_exato.sh "$INPUTDIF2" "$NAMEBIN" 
if [ $? == 1 ]; then
  exit 1
fi


test/binTest/run_dif_exato_Ex4.sh "$INPUTDIF4" "$NAMEBIN" 
if [ $? == 1 ];then
  exit 1
fi

#
test/binTest/run_dif_exato_Ex5.sh "$INPUTDIF5" "$NAMEBIN" 
if [ $? == 1 ];then
  exit 1
fi

# difusao 3D
cp test/input/dif/orthogonal/3D/GreenGaussCell/*.dat test/binTest
cp test/input/dif/orthogonal/3D/GreenGaussNode/*.dat test/binTest
cp test/input/dif/orthogonal/3D/LeastSquare/*.dat test/binTest
test/binTest/run_dif_exato_3D_Ex2.sh "$INPUTDIF3DEX2" "$NAMEBIN" 
if [ $? == 1 ];then
  exit 1
fi

# difusao 3D
tar -xvzf test/input/dif/non_orthogonal/3D/*.tar.gz -C test/binTest
cp test/input/dif/non_orthogonal/3D/GreenGaussCell/*.dat test/binTest
cp test/input/dif/non_orthogonal/3D/GreenGaussNode/*.dat test/binTest
cp test/input/dif/non_orthogonal/3D/LeastSquare/*.dat test/binTest
test/binTest/run_dif_exato_3D_Ex3.sh "$INPUTDIF3DEX3" "$NAMEBIN" 
if [ $? == 1 ];then
  exit 1
fi

DIR="test/binTest"
rm  $DIR/*.dat $DIR/*.vtk  $DIR/*.txt   $DIR/*.csv $DIR/$NAMEBIN $DIR/*.mtx 

exit 0
