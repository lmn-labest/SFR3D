#!/bin/sh
INPUT=$1
NAMEBIN=$2

usage(){
  if [ "$1" -lt 2 ]; then
    echo "Usage $0: [exc name] [input] [clean true|false]"
    exit 1
  fi
}

usage "$#"
DIR="test/binTest"

#executando os testes
for name in $INPUT; do
  echo --------------------------------------------------------------
  OUT=`echo $name | sed -e "s/.dat//"`
  echo  teste: $NAMEBIN $name $OUT
  echo `$DIR/$NAMEBIN $DIR/$name $DIR/$OUT` >/dev/null
  echo --------------------------------------------------------------
  OUT2=`echo $name | sed -e "s/.dat/_D1_cell_1.csv/"`
  $DIR/comp_3D_ex1.py $DIR/$OUT2 
  if [ $? == 1 ]; then
    echo --------------------------------------------------------------
    echo teste falhou !!! : $OUT2 
    echo --------------------------------------------------------------
    exit 1
  fi 
done
#----------------------------------------------------------------------

