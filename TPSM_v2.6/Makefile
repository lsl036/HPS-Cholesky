# ########################################################################
# VERSION: 2.6-Release
# AUTHOR:  蔡沁耘@湖南大学
# EMAIL:   hnutsai@hnu.edu.cn
# DATE:    2020年12月16日
# FILE:    tpsm_threads.c
# BRIEF:   TPSM并行框架Release版本(version 2.6),为TPSM_v2.5版本的精炼版本
# ########################################################################

CC = gcc -std=c99

TPSMflag = #-DSETSTACKSIZE #-DPREVENTION

CFLAG = -g -fPIC -O2 $(TPSMflag)

INCdir = ./include
SRCdir = ./src
LIBdir = ./lib
INCLUDE = -I$(INCdir)


libtpsm.a : main.o buffer.o synchronization.o tcb.o barrier.o threads.o distribution.o base.o auxiliary.o
	ar -rc $(LIBdir)/$@ $^
	rm $^
	
base.o : $(SRCdir)/tpsm_base.c $(INCdir)/tpsm_base.h
	$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@

distribution.o : $(SRCdir)/tpsm_distribution.c $(INCdir)/tpsm_distribution.h \
						 $(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h \
						 base.o
	$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@

buffer.o : $(SRCdir)/tpsm_buffer.c $(INCdir)/tpsm_buffer.h \
		   $(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
		   base.o
	$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@

synchronization.o : $(SRCdir)/tpsm_synchronization.c $(INCdir)/tpsm_synchronization.h \
					$(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
					base.o
	$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@

tcb.o : $(SRCdir)/tpsm_tcb.c $(INCdir)/tpsm_tcb.h \
		$(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
		base.o
	$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@

barrier.o : $(SRCdir)/tpsm_barrier.c $(INCdir)/tpsm_barrier.h \
		$(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
		base.o
	#$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@ -DDEBUG
	$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@

threads.o : $(SRCdir)/tpsm_threads.c $(INCdir)/tpsm_threads.h \
			$(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
			base.o tcb.o
	#$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@ -DDEBUG
	$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@ 

auxiliary.o:$(SRCdir)/tpsm_auxiliary.c $(INCdir)/tpsm_auxiliary.h
	$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@

main.o : $(SRCdir)/tpsm_main.c $(INCdir)/tpsm_main.h \
		 $(INCdir)/tpsm_sysinfo.h $(INCdir)/tpsm_default.h $(INCdir)/tpsm_structure.h \
		 base.o buffer.o synchronization.o tcb.o threads.o distribution.o auxiliary.o
	#$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@ -DDEBUG
	$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@ 




# ----------------------------------------------
# TEST
# ----------------------------------------------
test.o : ./test.c
	$(CC) $(CFLAG) $(INCLUDE) -c $< -o $@

xdtst: test.o  $(LIBdir)/libtpsm.a
	$(CC) $(CFLAG) -o $@  $^ -lpthread -lnuma




.PHONY : clean
clean :
	rm *.o $(LIBdir)/libtpsm.a