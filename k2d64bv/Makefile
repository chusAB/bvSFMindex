# Makefile written by Jesus Alastruey Benede (jalastru@unizar.es)
# October 2017

# Source file
VERSION=$(shell basename $(CURDIR))
# VERSION=k2d64bv

# Select the compiler,
#    c=0 corresponds icc
#    c=1 corresponds gcc
#    c=2 corresponds gcc-6
#    c=3 corresponds gcc-7
# You can override the default value from the command line!
# The command
#    "make c=2"
# will the use value c=2, regardless of the default.
CC=gcc
c=1

ifeq ($(c),0)
    CC=icc
else ifeq ($(c),1)
    CC=gcc
else ifeq ($(c),2)
    CC=gcc-6
else ifeq ($(c),3)
    CC=gcc-7
endif

# Specify number of overlapped sequence searches (s) per thread
# The command
#    "make s=4"
# will the use value s=4
s=4
OVERLAP_FLAG = -DNSEQS=$(s)

# Select the target architecture,
#    a=0 corresponds to generic architecture (gen)
#    a=1 corresponds to native architecture (native)
#    a=2 corresponds to Knights Landing (knl)
#    a=3 corresponds to Ivy Bridge architecture (ivb)
#    a=4 corresponds to Haswell architecture (hsw)
#    a=5 corresponds to Broadwell architecture (bdw)
#    a=6 corresponds to Skylake client architecture (skl)
#    a=7 corresponds to Skylake server architecture (skx)
# You can override the default value from the command line!
# The command
#    "make a=1"
# will the use value a=1, regardless of the default.
a=0

ifeq ($(a),0)
    ARCH=gen
else ifeq ($(a),1)
    ARCH=nat
else ifeq ($(a),2)
    ARCH=knl
else ifeq ($(a),3)
    ARCH=ivb
else ifeq ($(a),4)
    ARCH=hsw
else ifeq ($(a),5)
    ARCH=bdw
else ifeq ($(a),6)
    ARCH=skl
else ifeq ($(a),7)
    ARCH=skx
else
    ARCH=gen
endif

# Compiler dependent flags
# machine flags
ifeq ($(c),0)
    # icc
    ifeq ($(a),1)
        ARCH_FLAG= -xHost
    else ifeq ($(a),2)
        ARCH_FLAG = -xMIC-AVX512 -DKNL
    else ifeq ($(a),3)
        ARCH_FLAG = -xAVX
    else ifeq ($(a),4)
        ARCH_FLAG = -xCORE-AVX2
    else ifeq ($(a),5)
        ARCH_FLAG = -xCORE-AVX2
    else ifeq ($(a),6)
        ARCH_FLAG = -march=skylake
        # ARCH_FLAG = -xCORE-AVX2
    else ifeq ($(a),7)
        ARCH_FLAG = -xCORE-AVX512
    else
        ARCH_FLAG =
    endif
else
	# gcc
    ifeq ($(a),1)
        ARCH_FLAG = -march=native
    else ifeq ($(a),2)
        ARCH_FLAG = -march=knl  -DKNL
    else ifeq ($(a),3)
        ARCH_FLAG = -march=ivybridge
    else ifeq ($(a),4)
        ARCH_FLAG = -march=haswell
    else ifeq ($(a),5)
        ARCH_FLAG = -march=broadwell
    else ifeq ($(a),6)
        ARCH_FLAG = -march=skylake
    else ifeq ($(a),7)
        ARCH_FLAG = -march=skylake-avx512
    else
        ARCH_FLAG =
    endif
endif
# knl:            -mavx512f -mavx512cd -mavx512er -mavx512pf
# skylake-avx512: -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -mavx512ifma -mavx512vbmi

# More compiler dependent flags
# optimization and report flags
# -std=gnu11 is the default mode, it is needed by the scan-build clang static analyzer
ifeq ($(c),0)
    # icc
    CFLAGS = -O3  -W -Wall -Winline -ip -qopenmp $(ARCH_FLAG) -std=gnu11  -qno-opt-prefetch 
    # CFLAGS = -O0 -g  ...
    REPORT_FLAGS = -qopt-report=4 -qopt-report-phase ipo -qopt-report-file=stdout
else
	# gcc
	CFLAGS = -O3 -W -Wall -Wextra -Wshadow -Winline -fopenmp $(ARCH_FLAG) -std=gnu11 -fno-prefetch-loop-arrays
    # CFLAGS = -O0 -g  ...
	REPORT_FLAGS = -fopt-info-optimized
endif

# huge page support
CFLAGS := $(CFLAGS) -DHUGEPAGES

# single (L2) or dual prefetch (L2+L1)
p=dp
ifeq ($(p),dp)
    CFLAGS := $(CFLAGS) -DDUAL_PREFETCH
endif

# libnuma support
n=0
ifeq ($(n),1)
    CFLAGS := $(CFLAGS) -DLIBNUMA
endif

# perf counters
w=0
ifeq ($(w),1)
    CFLAGS := $(CFLAGS) -DPERF
endif

# debug threads
d=0
ifeq ($(w),1)
    CFLAGS := $(CFLAGS) -DDEBUG_THREADS
endif

# iaca support
i=0
ifeq ($(i),1)
    CFLAGS := $(CFLAGS) -DIACA
endif

# required libraries
CLIBS = -lm

# target dependent libraries
# Determine if the target architecture is Knights Landing or not
ifeq ($(a),1)
    # Knights Landing (Xeon Phi)
    CLIBS := $(CLIBS) -lmemkind  -L../lib
    # if libdivsufsort library is not system-wide installed:
    # CLIBS := $(CLIBS) -lmemkind -Wl,-rpath=/home/user/[path]/libdivsufsort/build/lib -L/home/user/[path]/libdivsufsort/build/lib/
# else
	# Skylake, Broadwell, Ivy Bridge, Native ...
    # CLIBS := $(CLIBS) -Wl,-rpath=/home/user/[path]/libdivsufsort/build/lib -L/home/user/[path]/libdivsufsort/build/lib/
endif

# optional library
# libnuma support
ifeq ($(n),1)
    CLIBS := $(CLIBS) -lnuma
endif


#
# Directories
#
# VPATH = ..
vpath %.c ..
OBJDIR := obj
BINDIR := ../bin
REPDIR := cc_report

#
# all sources
#
SRCS = $(wildcard *.c)

#
# common, build and count sources and objects
#
COMMON_SRCS = aux.c bit_mng.c file_mng.c $(VERSION).c
COMMON_OBJS = $(patsubst %.c, $(OBJDIR)/%.o, $(COMMON_SRCS))

BUILD_SRCS = BWT.c 
BUILD_OBJS = $(patsubst %.c, $(OBJDIR)/%.o, $(BUILD_SRCS))

COUNT_SRCS = mem.c perf.c
COUNT_OBJS = $(patsubst %.c, $(OBJDIR)/%.o, $(COUNT_SRCS))

# replace .c by .o and prepend $OBJDIR

#
# all target
#
all: $(BINDIR)/$(VERSION)_build.$(ARCH).$(CC)  $(BINDIR)/$(VERSION)_fcount.$(ARCH).$(CC).$(s)seq.$(p)

#
# alias to avoid errors when executing, for instance, $ make k2d64bv_fcount
#
$(VERSION)_fcount: $(BINDIR)/$(VERSION)_fcount.$(ARCH).$(CC).$(s)seq.$(p)
	@echo "HOLA" > /dev/null
$(VERSION)_fcount.$(ARCH).$(CC).$(s)seq.$(p): $(BINDIR)/$(VERSION)_fcount.$(ARCH).$(CC).$(s)seq.$(p)
	@echo "HOLA" > /dev/null

$(VERSION)_build: $(BINDIR)/$(VERSION)_build.$(ARCH).$(CC)
	@echo "HOLA" > /dev/null
$(VERSION)_build.$(ARCH).$(CC): $(BINDIR)/$(VERSION)_build.$(ARCH).$(CC)
	@echo "HOLA" > /dev/null

#
# dependencies to force the creation of the $(OBJDIR) and $(BINDIR) directories
# @: suppress the echoing of the command
#
$(OBJDIR):
	@mkdir -p $@
$(BINDIR):
	@mkdir -p $@
$(REPDIR):
	@mkdir -p $@	

#
# there are three different kinds of compilation linesth
# https://www.gnu.org/software/make/manual/html_node/Static-Usage.html
#
# $<: name of the first prerequisite
# $@: file name of the target of the rule
# $(*F): the file-within-directory part of the stem. If the value of ‘$@’ is dir/foo.o then ‘$(*F)’ is foo 
# http://www.gnu.org/software/make/manual/make.html#Automatic-Variables
#
SRCS1 = aux.c mem.c file_mng.c BWT.c bit_mng.c perf.c $(VERSION).c $(VERSION)_build.c
$(patsubst %.c, $(OBJDIR)/%.o, $(SRCS1)): $(OBJDIR)/%.o: %.c | $(OBJDIR) $(REPDIR)
	$(CC)  $(CFLAGS)  -c $<  -o $@  | tee $(REPDIR)/$(VERSION).$(*F).$(ARCH).$(CC).txt 2>&1
#	@$(CC)  $(CFLAGS)  -c $<  $(CLIBS)  -o $@  > $(REPDIR)/$(VERSION).$(*F).$(ARCH).$(CC).txt 2>&1

SRCS2 = $(VERSION)_fcount.c 
$(patsubst %.c, $(OBJDIR)/%.o, $(SRCS2)): $(OBJDIR)/%.o: %.c | $(OBJDIR) $(REPDIR)
	$(CC)  $(CFLAGS) $(OVERLAP_FLAG) $(REPORT_FLAGS) -g -c $<  -o $@  | tee $(REPDIR)/$(VERSION)_fcount.$(ARCH).$(CC).$(s)seq.$(p).report.txt 2>&1
#	@$(CC)  $(CFLAGS) $(OVERLAP_FLAG) $(REPORT_FLAGS) -g -c $<  $(CLIBS) -o $@  > $(REPDIR)/$(VERSION)_fcount.$(ARCH).$(CC).$(s)seq.$(p).report.txt 2>&1
#	$(CC)  $(CFLAGS) $(OVERLAP_FLAG) $(REPORT_FLAGS) -g -c $<  $(CLIBS) -o $@  -Wa,-adghln=$(VERSION)_fcount.$(ARCH).$(CC).s   > $(VERSION)_fcount.$(ARCH).$(CC).report.txt 2>&1

#
# add this line to generate assembly code
#	-Wa,-adghln=$(VERSION)_fcount.$(ARCH).$(CC).s
#

# create list of auto dependencies
AUTODEPS:= $(patsubst %.c, $(OBJDIR)/%.d, $(SRCS))
# https://www.gnu.org/software/make/manual/html_node/Include.html#Include
-include $(AUTODEPS)

#
# binaries
# $^: names of all the prerequisites
# http://www.gnu.org/software/make/manual/make.html#Automatic-Variables

$(BINDIR)/$(VERSION)_build.$(ARCH).$(CC): $(COMMON_OBJS) $(BUILD_OBJS) $(OBJDIR)/$(VERSION)_build.o | $(BINDIR)
	$(CC)  $(CFLAGS)  $^  $(CLIBS) -ldivsufsort64 -o $@  && strip $@

$(BINDIR)/$(VERSION)_fcount.$(ARCH).$(CC).$(s)seq.$(p): $(COMMON_OBJS) $(COUNT_OBJS) $(OBJDIR)/$(VERSION)_fcount.o | $(BINDIR)
	$(CC)  $(CFLAGS)  $^  $(CLIBS)  -o $@  && strip $@

#
# other targets
#
clean:
	@rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d $(BINDIR)/$(VERSION)_build.$(ARCH).$(CC) $(VERSION)_fcount.$(ARCH).$(CC).$(s)seq.$(p)

flags:
	@$(CC) $(CFLAGS) -E -v - </dev/null 2>&1 | grep cc1

gccversion:
	@$(CC) -v 2>&1 | tail -1

.PHONY: clean all flags gccversion
