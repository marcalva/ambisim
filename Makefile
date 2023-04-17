
CC = gcc
# CFLAGS = -g -O1 -Wall -Wextra -Wfloat-equal -Wno-unused-function -fsanitize=address -fno-omit-frame-pointer
CFLAGS = -g -O2 -Wall -Wextra -Wfloat-equal -Wno-unused-function -Wpointer-arith -Wshadow

all : ambisim

HTSDIR = htslib
include $(HTSDIR)/htslib.mk
include $(HTSDIR)/htslib_static.mk
HTSLIB = $(HTSDIR)/libhts.a

CPPFLAGS = -I. -I$(HTSDIR)

OBJS = main.o rvd.o array_util.o bins.o gtf_anno.o str_util.o overlap.o \
	   variants.o counts.o region.o sam_read.o gex_prob.o atac_prob.o bc_sim.o

LDFLAGS = -L$(HTSDIR)
LIBS = -lm -lhts -lpthread $(HTSLIB_static_LIBS)

ambisim : $(OBJS) $(HTSLIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $< 

cleano :
	rm -f *o
	rm -f ambisim

clean:
	rm -f *o
	rm -f ambimux
	cd $(HTSDIR) && make clean && cd ../

