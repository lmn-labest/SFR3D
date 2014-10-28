#!/bin/sh
INPUTRCM='retangulo.dat retangulo_rcm.dat'
INPUTVTK='retangulo.dat'
NAMEBIN=mfvCell_gnu_O3

#gerando o executavel
echo `make clean` > /dev/null 
echo `make` > /dev/null       
mv bin/$NAMEBIN test/binTest
cp test/input/*.dat test/binTest
#rodando os teste COO
test/binTest/run_coo.sh "$INPUTRCM" "$NAMEBIN" 

#rodando os teste VTK
test/binTest/run_vtk.sh "$INPUTVTK" "$NAMEBIN" 

DIR="test/binTest"
#rm $DIR/*.mtx $DIR/*.dat $DIR/*.vtk $DIR/$NAMEBIN

exit 1
