echo off
cls

set arg1=%1
set arg2=%2

echo Modificando os nomes dos arquivos

set nameFile=%1
set n=%2

:: Nome dos arquivos
FOR /L %%I IN (0,1,%n%) DO (

	mv %nameFile%_%%I_0.vtu %nameFile%.%%I.vtu
)