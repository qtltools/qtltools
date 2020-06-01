##########################################
# SET THESE 6 PATHS CORRECTLY TO COMPILE #
##########################################
BOOST_INC=
BOOST_LIB=
RMATH_INC=
RMATH_LIB=
HTSLD_INC=
HTSLD_LIB=
#########################################################
# EXAMPLES:                                             #
# BOOST_INC=/usr/include/                               #
# BOOST_LIB=/usr/lib/x86_64-linux-gnu/                  #
# RMATH_INC=$(HOME)/Tools/R-3.6.1/src/include           #
# RMATH_LIB=$(HOME)/Tools/R-3.6.1/src/nmath/standalone  #
# HTSLD_INC=$(HOME)/Tools/htslib-1.9                    #
# HTSLD_LIB=$(HOME)/Tools/htslib-1.9                    #
#########################################################
define n


endef

#COMPILER MODE C++11
CXX=g++ -std=c++0x

#COMPILER FLAGS
CXXFLAG_REL=-O3
CXXFLAG_DBG=-g
CXXFLAG_WRN=-Wall -Wextra -Wno-sign-compare -Wno-unused-local-typedefs -Wno-deprecated -Wno-unused-parameter

#BASE LIBRARIES
LIB_FLAGS=-lz -lgsl -lbz2 -llzma -lgslcblas -lm -lpthread

#FILE LISTS
BFILE=bin/QTLtools
HFILE=$(shell find src -name *.h)
TFILE=$(shell find lib -name *.h)
CFILE=$(shell find src -name *.cpp | LC_ALL=C sort)
OFILE=$(shell for file in `find src -name *.cpp | LC_ALL=C sort`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp | LC_ALL=C sort`; do echo $$(dirname $$file); done)
ifneq (, $(shell which git))
 GITVS=$(shell git describe --tags --long --abbrev=10 2>/dev/null)
 ifneq (, $(GITVS))
  $(info Compiling version $(GITVS))
  CXXFLAG_REL+= -DQTLTOOLS_VERSION=\"$(GITVS)\"
  CXXFLAG_DBG+= -DQTLTOOLS_VERSION=\"$(GITVS)\"
 endif
endif

#STATICLY LINKED LIBS
LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
#INCLUDE DIRS
IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)

#ONLY FOR MAC STATIC LINKING, ARCHIVES ASSUMED TO BE INSTALLED WITH BREW.IF YOU HAVE THESE IN OTHER LOCATIONS MODIFY THE NEXT 5 LINES 
MZ=/usr/local/opt/zlib/lib/libz.a
MCBLAS=/usr/local/lib/libgslcblas.a
MGSL=/usr/local/lib/libgsl.a
MBZ2=/usr/local/opt/bzip2/lib/libbz2.a
MLZMA=/usr/local/lib/liblzma.a

#MAC SPECIFIC STUFF
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
 CXXFLAG_REL+= -fvisibility=hidden -fvisibility-inlines-hidden
 CXXFLAG_DBG+= -fvisibility=hidden -fvisibility-inlines-hidden
 ifeq ($(MAKECMDGOALS),static)
   ifeq ("$(wildcard  $(MZ))","")
     $(error Cannot find $(MZ)! $nTry:	brew install zlib$nOr edit the makefile with the correct location of libz.a)
   endif
   ifeq ("$(wildcard  $(MCBLAS))","")
     $(error Cannot find $(MCBLAS)! $nTry:	brew install gsl$nOr edit the makefile with the correct location of libgslcblas.a)
   endif
   ifeq ("$(wildcard  $(MGSL))","")
     $(error Cannot find $(MGSL)! $nTry:	brew install gsl$nOr edit the makefile with the correct location of libgsl.a)
   endif
   ifeq ("$(wildcard  $(MBZ2))","")
     $(error Cannot find $(MBZ2)! $nTry:	brew install bzip2$nOr edit the makefile with the correct location of libbz2.a)
   endif
   ifeq ("$(wildcard  $(MLZMA))","")
     $(error Cannot find $(MLZMA)! $nTry:	brew install xz$nOr edit the makefile with the correct location of liblzma.a)
   endif        
 endif
endif

#DEFAULT VERSION (SET UP THE VARIABLES IN THE BEGINING OF THE MAKEFILE)
all: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
all: LDFLAG=$(CXXFLAG_REL)
all: $(BFILE)

#DEBUG VERSION (SET UP THE VARIABLES IN THE BEGINING OF THE MAKEFILE)
debug: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
debug: LDFLAG=$(CXXFLAG_DBG)
debug: $(BFILE)

#STATIC VERSION (SET UP THE VARIABLES IN THE BEGINING OF THE MAKEFILE)
static: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
static: LDFLAG=$(CXXFLAG_REL)
ifeq ($(UNAME),Darwin)
#ASSUMES YOU INSTALLED REQUIRED LIBRARIES WITH BREW. SEE ABOVE WHERE THESE VARIABLES ARE SET
static: LIB_FILES+= $(MZ) $(MCBLAS) $(MGSL) $(MBZ2) $(MLZMA)
static: LIB_FLAGS=-lm -lpthread
else
static: LIB_FLAGS=-Wl,-Bstatic -lz -lgsl -lbz2 -llzma -lgslcblas -Wl,-Bdynamic -lm -lpthread
endif
static: $(BFILE)

#STATIC DEBUG VERSION (SET UP THE VARIABLES IN THE BEGINING OF THE MAKEFILE)
static-dbg: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
static-dbg: LDFLAG=$(CXXFLAG_DBG)
ifeq ($(UNAME),Darwin)
#ASSUMES YOU INSTALLED REQUIRED LIBRARIES WITH BREW, SEE ABOVE WHERE THESE VARIABLES ARE SET
static-dbg: LIB_FILES+= $(MZ) $(MCBLAS) $(MGSL) $(MBZ2) $(MLZMA)
static-dbg: LIB_FLAGS=-lm -lpthread
else
static-dbg: LIB_FLAGS=-Wl,-Bstatic -lz -lgsl -lbz2 -llzma -lgslcblas -Wl,-Bdynamic -lm -lpthread
endif
static-dbg: $(BFILE)

personal: BOOST_INC=/usr/include/
personal: BOOST_LIB=/usr/lib/x86_64-linux-gnu/
personal: RMATH_INC=$(HOME)/Tools/R-3.6.1/src/include
personal: RMATH_LIB=$(HOME)/Tools/R-3.6.1/src/nmath/standalone
personal: HTSLD_INC=$(HOME)/Tools/htslib-1.9
personal: HTSLD_LIB=$(HOME)/Tools/htslib-1.9
personal: all

personal: BOOST_INC=/usr/include/
personal: BOOST_LIB=/usr/lib/x86_64-linux-gnu/
personal: RMATH_INC=$(HOME)/Tools/R-3.6.1/src/include
personal: RMATH_LIB=$(HOME)/Tools/R-3.6.1/src/nmath/standalone
personal: HTSLD_INC=$(HOME)/Tools/htslib-1.9
personal: HTSLD_LIB=$(HOME)/Tools/htslib-1.9
personal: debug


baobab: BOOST_INC=/srv/beegfs/scratch/groups/funpopgen/Tools/boost_1_71_0/
baobab: BOOST_LIB=/srv/beegfs/scratch/groups/funpopgen/Tools/boost_1_71_0/stage/lib/
baobab: RMATH_INC=/srv/beegfs/scratch/groups/funpopgen/Tools/R-3.6.1/src/include/
baobab: RMATH_LIB=/srv/beegfs/scratch/groups/funpopgen/Tools/R-3.6.1/src/nmath/standalone/
baobab: HTSLD_INC=/srv/beegfs/scratch/groups/funpopgen/Tools/htslib-1.9/
baobab: HTSLD_LIB=/srv/beegfs/scratch/groups/funpopgen/Tools/htslib-1.9/
baobab: static

baobab-dbg: BOOST_INC=/srv/beegfs/scratch/groups/funpopgen/Tools/boost_1_71_0/
baobab-dbg: BOOST_LIB=/srv/beegfs/scratch/groups/funpopgen/Tools/boost_1_71_0/stage/lib/
baobab-dbg: RMATH_INC=/srv/beegfs/scratch/groups/funpopgen/Tools/R-3.6.1/src/include/
baobab-dbg: RMATH_LIB=/srv/beegfs/scratch/groups/funpopgen/Tools/R-3.6.1/src/nmath/standalone/
baobab-dbg: HTSLD_INC=/srv/beegfs/scratch/groups/funpopgen/Tools/htslib-1.9/
baobab-dbg: HTSLD_LIB=/srv/beegfs/scratch/groups/funpopgen/Tools/htslib-1.9/
baobab-dbg: static-dbg

#COMPILATION RULES
$(BFILE): $(OFILE)
	$(CXX) $^ $(LIB_FILES) -o $@ $(LIB_FLAGS) $(LDFLAG)

obj/QTLtools.o: src/QTLtools.cpp $(HFILE) $(TFILE) $(CFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

obj/data.o: src/common/data.cpp src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

obj/cis_%.o: cis_%.cpp cis_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/correct_%.o: correct_%.cpp correct_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/fenrich_%.o: fenrich_%.cpp fenrich_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/match_%.o: match_%.cpp match_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/trans_%.o: trans_%.cpp trans_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/rtc_%.o: rtc_%.cpp rtc_data.h  src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

obj/pca_%.o: pca_%.cpp pca_data.h pca_pca.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/genrich_%.o: genrich_%.cpp genrich_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/union_%.o: union_%.cpp union_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/extract_%.o: extract_%.cpp extract_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

obj/quan_%.o: quan_%.cpp quan_data.h quan_xxhash.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/ase_%.o: ase_%.cpp ase_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/bamstat_%.o: bamstat_%.cpp bamstat_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/fdensity_%.o: fdensity_%.cpp fdensity_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/rep_%.o: rep_%.cpp rep_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/gwas_%.o: gwas_%.cpp gwas_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
clean: 
	rm -f obj/*.o $(BFILE)

clean-cis:
	rm -f obj/cis_*.o $(BFILE)

clean-correct:
	rm -f obj/correct_*.o $(BFILE)

clean-fenrich:
	rm -f obj/fenrich_*.o $(BFILE)

clean-genrich:
	rm -f obj/genrich_*.o $(BFILE)

clean-match:
	rm -f obj/match_*.o $(BFILE)

clean-trans:
	rm -f obj/trans_*.o $(BFILE)

clean-rtc:
	rm -f obj/rtc_*.o $(BFILE)

clean-pca:
	rm -f obj/pca_*.o $(BFILE)
	
clean-extract:
	rm -f obj/extract_*.o $(BFILE)
	
clean-ase:
	rm -f obj/ase_*.o $(BFILE)

clean-union:
	rm -f obj/union_*.o $(BFILE)

clean-quan:
	rm -f obj/quan_*.o $(BFILE)

clean-bamstat:
	rm -f obj/bamstat_*.o $(BFILE)
		
clean-fdensity:
	rm -f obj/fdensity_*.o $(BFILE)
		
clean-rep:
	rm -f obj/rep_*.o $(BFILE)
		
clean-gwas:
	rm -f obj/gwas_*.o $(BFILE)
		
		