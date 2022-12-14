#-------------------------------------------------------------------------------
# HNU_BASE Makefile
#-------------------------------------------------------------------------------
include ../../Makefile.option
INSTALL = $(realpath $(CURDIR)/../..)
AR_TARGET = libsparsebase.a

#-------------------------------------------------------------------------------

default: base

base:
	make install INSTALL=$(INSTALL)

install: $(AR_TARGET)
	@mkdir -p $(INSTALL)/lib
	@mkdir -p $(INSTALL)/include
	rm *.o
	mv -f *.a $(INSTALL)/lib

uninstall: 
	rm $(INSTALL)/lib/$(AR_TARGET)

#-------------------------------------------------------------------------------

# CC = gcc -std=c99 #-g
# CF = -O2 -fexceptions -fPIC -fopenmp
# I = -I../../include -I../../include/tpsm
# C = $(CC) $(CF) $(I)

INC = ../../include/SparseBase_config.h ../../include/xerbla.h ../../include/colamd.h \
	  ../../include/amd.h ../../include/amd_internal.h ../../include/tpsm/tpsm.h \
	  ../../include/tpsm/tpsm_base.h

OBJ = SparseBase_config.o xerbla.o amd.o colamd.o base.o distribution.o \
	  buffer.o synchronization.o tcb.o barrier.o threads.o auxiliary.o tpsm.o
	  #threadbase.o threadpool_manager.o

$(OBJ): $(INC)

$(AR_TARGET): $(OBJ)
	ar rv $@ $^
	- ranlib $@

#-------------------------------------------------------------------------------

SparseBase_config.o: SparseBase_config.c
	$(C) -c $< -o $@

xerbla.o: xerbla.c
	$(C) -c $< -o $@

amd.o: amd.c
	$(C) $(LONG_FLAGS) -c $< -o $@

colamd.o: colamd.c
	$(C) $(LONG_FLAGS) -c $< -o $@

#  ------------------------------
#           线程池
#  ------------------------------
INCdir = ../../include/tpsm

base.o : ./tpsm_base.c $(INCdir)/tpsm_base.h
	$(C) -c $< -o $@

distribution.o : ./tpsm_distribution.c $(INCdir)/tpsm_distribution.h \
						 $(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h \
						 base.o
	$(C) -c $< -o $@

buffer.o : ./tpsm_buffer.c $(INCdir)/tpsm_buffer.h \
		   $(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
		   base.o
	$(C) -c $< -o $@

synchronization.o : ./tpsm_synchronization.c $(INCdir)/tpsm_synchronization.h \
					$(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
					base.o
	$(C) -c $< -o $@

tcb.o : ./tpsm_tcb.c $(INCdir)/tpsm_tcb.h \
		$(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
		base.o
	$(C) -c $< -o $@

barrier.o : ./tpsm_barrier.c $(INCdir)/tpsm_barrier.h \
		$(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
		base.o
	$(C) -c $< -o $@

threads.o : ./tpsm_threads.c $(INCdir)/tpsm_threads.h \
			$(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
			base.o tcb.o
	$(C) -c $< -o $@

auxiliary.o:./tpsm_auxiliary.c $(INCdir)/tpsm_auxiliary.h
	$(C) -c $< -o $@

tpsm.o : ./tpsm.c $(INCdir)/tpsm.h \
		 $(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
		 base.o buffer.o synchronization.o tcb.o threads.o distribution.o auxiliary.o
	$(C) -c $< -o $@