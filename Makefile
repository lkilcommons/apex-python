AR=ar

####################
# IBM XLF compiler #
####################
#FCOMPILER=ibm
#FC=xlf
#OPT=-O3
#FFLAGS=$(OPT)

##################
# Intel Compiler #
##################
#FCOMPILER=intelem
#FC=ifort
#OPT=-O3
#FFLAGS=$(OPT)
#Debugging flags (leave these commented out unless you know what you're doing):
##FFLAGS+= -debug full -traceback -check bounds -check format -check output_conversion -check pointers -check uninit

################
# GNU compiler #
################
FCOMPILER=gnu95
FC=gfortran
OPT=-g
FFLAGS=$(OPT)
#Debugging flags (leave these commented out unless you know what you're doing):
##FFLAGS+=-Wall -Wextra -pedantic -fimplicit-none -fbounds-check

###########################
# Portland Group compiler #
###########################
#FCOMPILER=pg
#FC=pgf90
#OPT=-O3
#FFLAGS=$(OPT)
#Debugging flags (leave these commented out unless you know what you're doing):
##FFLAGS+= -C -Mchkfpstk -Mchkstk -Mpgicoff -traceback -Ktrap=fp -Minform=inform
##FFLAGS+= -fpic -Minform=inform -Mnosecond_underscore -fast
## -Ktrap=align,denorm,divz,fp,inexact,inv,ovf,unf

########################################################################
# You shouldn't need to touch anything below this line.
########################################################################
.SUFFIXES : .o .f90 .f
%.o: %.f90
	$(FC) $(FFLAGS) -c  -o $@ $<

default:
	@ echo "   ==============================================="
	@ echo "   =    APEX Coordinate System Transformation    ="
	@ echo "   =    http://www.hao.ucar.edu/apex             ="
	@ echo "   ==============================================="
	@ echo "   "
	@ echo "    To build, type"
	@ echo "   "
	@ echo "       make <target>"
	@ echo "   "
	@ echo "   where <target> is \"libapex.a\", \"apex.so\", \"test\", or \"test_f\""
	@ echo ""
	@ echo "   libapex.a:  Static APEX library. Installed to lib/"
	@ echo "   apex.so:  Python module to Apex routines. Installed to lib/"
	@ echo "   test:  Execute Python unit tests which validate the Apex installation"
	@ echo "   test_f: Fortran test program.  Installed to bin/"
	@ echo "   "
	@ echo "   To build everything, type"
	@ echo "   "
	@ echo "       make all"

BINDIR=bin
LIBDIR=lib
SRCDIR=apexpy
ARFLAGS = cvru

all: libapex.a apex.so test

libapex.a: $(addprefix $(SRCDIR)/, apex.o apxntrp.o apxntrp_legacy.o cossza.o ggrid.o magfld.o magloctm.o subsol.o)
	$(AR) $(ARFLAGS) $(LIBDIR)/$@ $?
	ranlib $(LIBDIR)/$@

apex.so:
	@ echo "Warning: You might have to clean up apex.pyf by hand.  Check the stdout carefully!"
	cd $(SRCDIR) && f2py *.f *.f90 -m apex -h apex.pyf
	#f2py --overwrite-signature apex.f apxntrp.f90 cossza.f ggrid.f magfld.f magloctm.f subsol.f -m apex -h apex.pyf 2>&1 | tee apex.pyf.out
	#cd $(SRCDIR) && f2py --fcompiler=$(FCOMPILER) -c apex.pyf *.f *.f90
	python setup.py config_fc --fcompiler=$(FCOMPILER) --opt=$(OPT) install --install-lib $(LIBDIR)

test: apex.so
	PYTHONPATH=lib/:$$PYTHONPATH nosetests -w src

test_f: libapex.a
	$(FC) $(SRCDIR)/test_apxntrp.f -o $(BINDIR)/test_apxntrp -L $(LIBDIR) -lapex
	$(FC) $(SRCDIR)/test_foster.f -o $(BINDIR)/test_foster -L $(LIBDIR) -lapex

clean:
	rm -rf build
	rm -f $(SRCDIR)/*.o \
		$(LIBDIR)/libapex.a $(LIBDIR)/apex.so \
		$(BINDIR)/test_apxntrp $(BINDIR)/test_foster
