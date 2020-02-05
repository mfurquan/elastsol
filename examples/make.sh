#!/bin/bash
MKLROOT="/opt/intel/mkl/"
F95ROOT=/opt/intel/mkl/
gfortran -o build/beam vtk_legacy.f90 csr3s.f90 hexa_lagr.f90 elastsol.f90 beam.f90 -fdefault-integer-8 -I${F95ROOT}/include/intel64/ilp64 -m64 -I${MKLROOT}/include ${F95ROOT}/lib/intel64/libmkl_blas95_ilp64.a ${F95ROOT}/lib/intel64/libmkl_lapack95_ilp64.a -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl \
-Wall -g -fbacktrace -fimplicit-none -fcheck=all -fcheck=bounds
