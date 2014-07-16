#!/bin/sh
INPUT='retangulo.dat retangulo_rcm.dat'
NAMEBIN=mfvCell_gnu_O0

#gerando o executavel
echo `make clean` > /dev/null 
echo `make` > /dev/null       
mv bin/$NAMEBIN test/binTest
cp test/input/*.dat test/binTest
#rodando os teste
test/binTest/run.sh "$INPUT" "$NAMEBIN"

exit 1
