SRCDIR = ../src
MODDIR = ../include
OBJDIR = ../build
MKLROOT = /opt/intel/mkl

EXE = $(OBJDIR)/test1

FC = gfortran
FFLAGS = -g -O0 -Wall -Wconversion -fimplicit-none -fbacktrace \
	 -ffree-line-length-0 -finit-real=nan -fcheck=all -fcheck=bounds \
	 -J$(MODDIR) -fdefault-integer-8 \
	 -I${MKLROOT}/include/intel64/ilp64 -m64 -I${MKLROOT}/include

LD = gfortran
LFLAGS = -g -fbounds-check -Wall -fbacktrace \
	 ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a \
	 ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a \
	 -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 \
	 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

$(EXE) : $(OBJDIR)/csr3s.o $(OBJDIR)/hexa_lagr.o $(OBJDIR)/elastsol.o \
	$(OBJDIR)/test1.o
	$(LD) -o $(EXE) $(OBJDIR)/csr3s.o $(OBJDIR)/hexa_lagr.o $(OBJDIR)/elastsol.o \
		$(OBJDIR)/test1.o $(LFLAGS)

$(OBJDIR)/csr3s.o $(MODDIR)/csr3s.mod : $(SRCDIR)/csr3s.f90
	$(FC) -c -o $(OBJDIR)/csr3s.o $< $(FFLAGS)

$(OBJDIR)/hexa_lagr.o $(MODDIR)/hexa_lagr.mod : $(SRCDIR)/hexa_lagr.f90
	$(FC) -c -o $(OBJDIR)/hexa_lagr.o $< $(FFLAGS)

$(OBJDIR)/elastsol.o $(MODDIR)/elastsol.mod : $(SRCDIR)/elastsol.f90 \
	$(MODDIR)/csr3s.mod $(MODDIR)/hexa_lagr.mod
	$(FC) -c -o $(OBJDIR)/elastsol.o $< $(FFLAGS)

$(OBJDIR)/test1.o : test1.f90 $(MODDIR)/elastsol.mod
	$(FC) -c -o $(OBJDIR)/test1.o $< $(FFLAGS)

.PHONY:clean cleanall
clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod
cleanall:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(EXE)
