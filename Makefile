
# CC = $(shell which gcc)
# CFLAGS = -g -O1 -Wall -Wextra -Wfloat-equal -Wno-unused-function -fsanitize=address -fno-omit-frame-pointer
CFLAGS += -g -O2 -Wall -Wextra -Wfloat-equal -Wno-unused-function -Wpointer-arith -Wshadow

HTSDIR = htslib

CPPFLAGS += -I. -I$(HTSDIR)

OBJS = main.o rvd.o array_util.o bins.o gtf_anno.o str_util.o overlap.o \
	   variants.o counts.o region.o sam_read.o gex_prob.o atac_prob.o sam_prob.o bc_sim.o

ambsim_make : ambisim

# $(HTSLIB):
# 	+cd $(HTSDIR) && $(MAKE) lib-static

-include $(HTSDIR)/htslib.mk
-include $(HTSDIR)/htslib_static.mk


LDFLAGS += -L$(HTSDIR)
LIBS += -lm -l:libhts.a -lpthread $(HTSLIB_static_LIBS)

hts :
	echo "building htslib"
	rm -rf $(HTSDIR)
	git clone --branch 1.17 https://github.com/samtools/htslib.git $(HTSDIR)
	cd $(HTSDIR) && git submodule update --init --recursive
	cd $(HTSDIR) && autoreconf -i && ./configure
	cd $(HTSDIR) && make lib-static

check_lib :
	if [ ! -d "$(HTSDIR)" ]; then \
		echo "htslib subdirectory not found, run 'make hts' first"; \
		exit 1; \
	fi

check_static :
	if [ ! -s "$(HTSDIR)/htslib_static.mk" ]; then \
		echo "htslib_static.mk not found, run 'make hts' first"; \
		exit 1; \
	fi

ambisim : check_lib check_static $(HTSDIR)/libhts.a $(OBJS)
	echo "building ambisim"
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $< 

cleano :
	rm -f *o
	rm -f ambisim

clean:
	rm -f *o
	rm -f ambimux
	cd $(HTSDIR) && make clean && cd ../

