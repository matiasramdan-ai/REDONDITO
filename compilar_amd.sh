#!/bin/bash
# Compila REDONDITO para las GPU AMD MI210 (gfx90a) con amdflang y
# OpenMP offload, como build/build_amd.sh de SPECTER: sphp-gpu
# (src/sphp-gpu/), el manejo del dispositivo y de la memoria
# (src/memmgt/), el timer (src/utils/), las rutinas pseudoespectrales
# (src/pseudo/) y la prueba de la GPU contra la CPU (src/). Los
# módulos y los objetos quedan en build/amd/ y el ejecutable en
# bin/amd/; corre solo en un nodo con GPU.
set -e
cd $(dirname $0)
module purge; module load rocm/7.2
module load llvm/22.0
module load openmpi5-rocm/5.0.7
module load fftw-gnu/3.3.11
FFTW=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11
SHTNS=$HOME/SHTNS
SP=../../src/sphp-gpu
FC="mpifort -cpp -O3 -funroll-loops -DGDOUBLE_PRECISION -DIKIND_=8"
FC="$FC -DCSIZE_=32 -DNSTRIP_=1 -DGHOST_GPU -DGPU_AMD"
FC="$FC -fopenmp --offload-arch=gfx90a -fopenmp-version=52"
FC="$FC -Wno-incompatible-implicit-interfaces"
FC="$FC -I$SP -I$FFTW/include -I$SHTNS/include"
LIBS="$SHTNS/lib/libshtns_cuda_omp.a -L$FFTW/lib -lfftw3_omp -lfftw3"
LIBS="$LIBS -L$ROCM_PATH/lib -lhipfft -lamdhip64 -lhiprtc -lstdc++ -lm"
mkdir -p build/amd bin/amd && cd build/amd
# La biblioteca, una sola vez
$FC -c ../../src/utils/gtimer_mod.f90 ../../src/memmgt/gdevice_mod.f90 \
  $SP/fprecision_mod.F90 $SP/gfft_bindings.F90 $SP/sphp_mod.F90 \
  $SP/fcgram_mod.f90 $SP/sphp_sca.F90 $SP/sphp_vec.F90 \
  ../../src/memmgt/gmem_mod.f90 ../../src/pseudo/pseudospec_mod.f90 \
  ../../src/memmgt/workspace.f90
OBJ="gtimer_mod.o gdevice_mod.o fprecision_mod.o gfft_bindings.o"
OBJ="$OBJ sphp_mod.o fcgram_mod.o sphp_sca.o sphp_vec.o gmem_mod.o"
OBJ="$OBJ pseudospec_mod.o workspace.o"
$FC ../../src/prueba_gpu.f90 $OBJ $LIBS -o ../../bin/amd/prueba_gpu
