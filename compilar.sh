#!/bin/bash
# Compila REDONDITO con gfortran en el login de sakura: la fftp de
# SPECTER que se va convirtiendo en sphp (src/) y sus pruebas. Los
# módulos y los objetos quedan en build/ y los ejecutables en bin/.
set -e
cd $(dirname $0)
module purge; module load gnu15 openmpi5 fftw/3.3.11
FFTW=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11
SHTNS=$HOME/SHTNS
FC="mpifort -cpp -O3 -DGDOUBLE_PRECISION -DIKIND_=8 -DCSIZE_=32"
FC="$FC -DNSTRIP_=1 -I../src -I$FFTW/include -I$SHTNS/include"
LIBS="$SHTNS/lib/libshtns_omp.a -L$FFTW/lib -lfftw3_omp -lfftw3"
LIBS="$LIBS -lgomp -lm"
mkdir -p build bin && cd build
# La biblioteca, una sola vez
$FC -c ../src/gtimer_mod.f90 ../src/fprecision_mod.F90 \
  ../src/sphp_mod.F90 ../src/fcgram_mod.f90 ../src/sphp_sca.F90 \
  ../src/sphp_vec.F90
OBJ="gtimer_mod.o fprecision_mod.o sphp_mod.o fcgram_mod.o sphp_sca.o"
OBJ="$OBJ sphp_vec.o"
# Una prueba por ejecutable
for p in prueba_sphp prueba_vec prueba_z prueba_3d; do
  $FC ../src/$p.f90 $OBJ $LIBS -o ../bin/$p
done
