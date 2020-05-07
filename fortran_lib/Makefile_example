.SUFFIXES: .f90

include ./make.sys

PROG = atomic
WORKDIR = EXEdir
SRCDIR = src
VPATH = $(SRCDIR):$(WORKDIR)


modn = atomic_constants.o atomic_control.o atomic_context.o
lev1 = atomic_clebsch.o atomic_angular.o atomic_util.o atomic_pfact.o gw_make_new.o
lev2 = atomic_stream.o atomic_basis.o atomic_hmat.o atomic_symm.o atomic_fmat.o atomic_tran_bd.o
dump = atomic_dumper.o atomic_print.o
main = atomic_main.o

objects = $(modn) $(lev1) $(lev2) $(dump) $(main)

default: headers $(PROG)
headers:
	@if [ ! -d $(WORKDIR)  ];then \
	mkdir $(WORKDIR);\
	echo "create workdir" $(WORKDIR); \
	fi

all: $(PROG)

$(PROG): $(objects)
	(cd $(WORKDIR);$(LINKER) $(objects) -o $(PROG) $(LFLAGS) $(LIBS);cd ..) #-Wl,-stack_size,0x100000000

.SUFFIXES: .o .f90
.f90.o:
	(cd $(WORKDIR);$(F90) $(FFLAGS) -c  ../$< -o $@;cd ..)

clobber:
	-rm -rf $(WORKDIR)

.PHONY : clean
clean: clobber

