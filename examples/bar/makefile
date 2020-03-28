SRCDIR = ../../src
MODDIR = ../../include
OBJDIR = ../../build
MKLROOT = /opt/intel/mkl

EXE = $(OBJDIR)/bar

FC = gfortran
FFLAGS = -J$(MODDIR) -fdefault-integer-8 \
	 -I${MKLROOT}/include/intel64/ilp64 -m64 -I${MKLROOT}/include \
	 -Ofast
#	 -ffree-line-length-0 -finit-real=nan -fcheck=all -fcheck=bounds \
#	 -g -O0 -Wall -Wconversion -fimplicit-none -fbacktrace
#	 -Ofast

LD = gfortran

LFLAGS = ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a \
	 ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a \
	 -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 \
	 -lmkl_sequential -lmkl_core -lpthread -lm -ldl \
	 -Ofast
#	 -g -fbounds-check -Wall -fbacktrace
#	 -Ofast

$(EXE) : $(OBJDIR)/csr3s.o $(OBJDIR)/hexa_lagr.o $(OBJDIR)/elastsol.o \
	$(OBJDIR)/vtk_legacy.o $(OBJDIR)/bar.o
	$(LD) -o $(EXE) $(OBJDIR)/csr3s.o $(OBJDIR)/hexa_lagr.o \
		$(OBJDIR)/elastsol.o $(OBJDIR)/vtk_legacy.o \
		$(OBJDIR)/bar.o $(LFLAGS)

$(OBJDIR)/csr3s.o $(MODDIR)/csr3s.mod : $(SRCDIR)/csr3s.f90
	$(FC) -c -o $(OBJDIR)/csr3s.o $< $(FFLAGS)

$(OBJDIR)/hexa_lagr.o $(MODDIR)/hexa_lagr.mod : $(SRCDIR)/hexa_lagr.f90
	$(FC) -c -o $(OBJDIR)/hexa_lagr.o $< $(FFLAGS)

$(OBJDIR)/elastsol.o $(MODDIR)/elastsol.mod : $(SRCDIR)/elastsol.f90 \
	$(MODDIR)/csr3s.mod $(MODDIR)/hexa_lagr.mod
	$(FC) -c -o $(OBJDIR)/elastsol.o $< $(FFLAGS)

$(OBJDIR)/vtk_legacy.o $(MODDIR)/vtk_legacy.mod : $(SRCDIR)/vtk_legacy.f90
	$(FC) -c -o $(OBJDIR)/vtk_legacy.o $< $(FFLAGS)

$(OBJDIR)/bar.o : bar.f90 $(MODDIR)/elastsol.mod $(MODDIR)/vtk_legacy.mod
	$(FC) -c -o $(OBJDIR)/bar.o $< $(FFLAGS)

.PHONY:clean cleanall
clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod
cleanall:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(EXE)
