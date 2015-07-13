#!/bin/sh
INPUTRCM='retangulo.dat retangulo_rcm.dat'
INPUTVTK='retangulo.dat'
INPUTDIF1='ex3_60_0_gglc.dat    ex3_60_3_gglc.dat
          ex3_60_0_ggln.dat    ex3_60_3_ggln.dat
          ex3_60_0_lSquare.dat ex3_60_3_lSquare.dat'

INPUTDIF2='ex2_1_gglc.dat       ex2_2_gglc.dat
          ex2_1_ggln.dat       ex2_2_ggln.dat
          ex2_1_lSquare.dat    ex2_2_lSquare.dat' 


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

#rodando os teste difusao
cp test/input/dif/non_orthogonal/2D/60/gaussGreenCell/*.dat test/binTest
cp test/input/dif/non_orthogonal/2D/60/gaussGreenNode/*.dat test/binTest
cp test/input/dif/non_orthogonal/2D/60/LeastSquare/*.dat test/binTest
test/binTest/run_dif.sh "$INPUTDIF1" "$NAMEBIN" 

cp test/input/dif/orthogonal/2D/gaussGreenCell/*.dat test/binTest
cp test/input/dif/orthogonal/2D/gaussGreenNode/*.dat test/binTest
cp test/input/dif/orthogonal/2D/LeastSquare/*.dat test/binTest
test/binTest/run_dif_exato.sh "$INPUTDIF2" "$NAMEBIN" 

DIR="test/binTest"
rm  $DIR/*.dat $DIR/*.vtk  $DIR/*.txt   $DIR/*.csv $DIR/$NAMEBIN $DIR/*.mtx

exit 1
