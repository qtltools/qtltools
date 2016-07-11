#COMPILER MODE C++11
CXX=g++ -std=c++11

#COMPILER FLAGS
CXXFLAG_REL=-O2
CXXFLAG_DBG=-g
CXXFLAG_WRN=-Wall -Wextra -Wno-sign-compare -Wno-unused-local-typedefs -Wno-deprecated -Wno-unused-parameter

#LINKER FLAGS
LDFLAG_REL=-O2
LDFLAG_DBG=-g

#BASE LIBRARIES
LIB_FLAGS=-lm -lz -lgsl -lblas -lbz2 -lpthread

#FILE LISTS
BFILE=bin/QTLtools
HFILE=$(shell find src -name *.h)
TFILE=$(shell find lib -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#DEFAULT VERSION (I.E. UNIGE DESKTOP RELEASE VERSION)
all: default

#FOR USER DEFINED COMPILATION [YOU NEED TO UPDATE THE SIX PATHS TO COMPILE QTLTOOLS]
default: RMATH_INC=
default: RMATH_LIB=
default: HTSLD_INC=
default: HTSLD_LIB=
default: BOOST_INC=
default: BOOST_LIB=
default: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
default: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
default: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
default: LDFLAG=$(LDFLAG_REL)
default: $(BFILE)

#UNIGE DESKTOP RELEASE VERSION
desktop: RMATH_INC=$(HOME)/Tools/R-3.2.2/src/include
desktop: RMATH_LIB=$(HOME)/Tools/R-3.2.2/src/nmath/standalone
desktop: HTSLD_INC=$(HOME)/Tools/htslib-1.3
desktop: HTSLD_LIB=$(HOME)/Tools/htslib-1.3
desktop: BOOST_INC=/usr/include
desktop: BOOST_LIB=/usr/lib/x86_64-linux-gnu
desktop: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
desktop: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
desktop: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
desktop: LDFLAG=$(LDFLAG_REL)
desktop: $(BFILE)

#UNIGE DESKTOP DEBUG VERSION
desktop-dbg: RMATH_INC=$(HOME)/Tools/R-3.2.2/src/include
desktop-dbg: RMATH_LIB=$(HOME)/Tools/R-3.2.2/src/nmath/standalone
desktop-dbg: HTSLD_INC=$(HOME)/Tools/htslib-1.3
desktop-dbg: HTSLD_LIB=$(HOME)/Tools/htslib-1.3
desktop-dbg: BOOST_INC=/usr/include
desktop-dbg: BOOST_LIB=/usr/lib/x86_64-linux-gnu
desktop-dbg: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
desktop-dbg: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
desktop-dbg: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
desktop-dbg: LDFLAG=$(LDFLAG_DBG)
desktop-dbg: $(BFILE)

#UNIGE SERVER RELEASE VERSION (=BINARY RELEASED TO PUBLIC)
server: RMATH_INC=/home/popgen/delaneau/SOFT/R-3.2.1/src/include
server: RMATH_LIB=/home/popgen/delaneau/SOFT/R-3.2.1/src/nmath/standalone
server: HTSLD_INC=/home/popgen/delaneau/SOFT/htslib-1.2.1
server: HTSLD_LIB=/home/popgen/delaneau/SOFT/htslib-1.2.1
server: BOOST_INC=/home/popgen/delaneau/SOFT/boost/boost_1_59_0/build/include
server: BOOST_LIB=/home/popgen/delaneau/SOFT/boost/boost_1_59_0/build/lib
server: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
server: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
server: LIB_FLAGS=$(LIB_FLAGS) -lgslcblas
server: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
server: LDFLAG=$(LDFLAG_REL)
server: $(BFILE)

#DELL LAPTOP RELEASE VERSION
laptop: RMATH_INC=$(HOME)/Libraries/R-3.2.2/src/include
laptop: RMATH_LIB=$(HOME)/Libraries/R-3.2.2/src/nmath/standalone
laptop: HTSLD_INC=$(HOME)/Libraries/htslib-1.2.1
laptop: HTSLD_LIB=$(HOME)/Libraries/htslib-1.2.1
laptop: BOOST_INC=/usr/include
laptop: BOOST_LIB=/usr/lib/x86_64-linux-gnu
laptop: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
laptop: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
laptop: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
laptop: LDFLAG=$(LDFLAG_REL)
laptop: $(BFILE)

#DELL LAPTOP DEBUG VERSION
laptop-dbg: RMATH_INC=$(HOME)/Libraries/R-3.2.2/src/include
laptop-dbg: RMATH_LIB=$(HOME)/Libraries/R-3.2.2/src/nmath/standalone
laptop-dbg: HTSLD_INC=$(HOME)/Libraries/htslib-1.2.1
laptop-dbg: HTSLD_LIB=$(HOME)/Libraries/htslib-1.2.1
laptop-dbg: BOOST_INC=/usr/include
laptop-dbg: BOOST_LIB=/usr/lib/x86_64-linux-gnu
laptop-dbg: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
laptop-dbg: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
laptop-dbg: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
laptop-dbg: LDFLAG=$(LDFLAG_DBG)
laptop-dbg: $(BFILE)

#VITAL-IT RELEASE VERSION
cluster: RMATH_INC=/software/R/3.1.1/include
cluster: RMATH_LIB=/software/R/3.1.1/lib64
cluster: HTSLD_INC=/software/UHTS/Analysis/samtools/1.2/include
cluster: HTSLD_LIB=/software/UHTS/Analysis/samtools/1.2/lib64
cluster: BOOST_INC=/software/include
cluster: BOOST_LIB=/software/lib64
cluster: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
cluster: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
cluster: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
cluster: LDFLAG=$(LDFLAG_REL)
cluster: $(BFILE)

#VITAL-IT DEBUG VERSION
cluster-dbg: RMATH_INC=/software/R/3.1.1/include
cluster-dbg: RMATH_LIB=/software/R/3.1.1/lib64
cluster-dbg: HTSLD_INC=/software/UHTS/Analysis/samtools/1.2/include
cluster-dbg: HTSLD_LIB=/software/UHTS/Analysis/samtools/1.2/lib64
cluster-dbg: BOOST_INC=/software/include
cluster-dbg: BOOST_LIB=/software/lib64
cluster-dbg: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
cluster-dbg: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
cluster-dbg: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
cluster-dbg: LDFLAG=$(LDFLAG_DBG)
cluster-dbg: $(BFILE)

#UBUNTU RELEASE VERSION
ubuntu: RMATH_INC=$(HOME)/R-3.2.2/src/include
ubuntu: RMATH_LIB=$(HOME)/R-3.2.2/src/nmath/standalone
ubuntu: HTSLD_INC=/usr/local/include/
ubuntu: HTSLD_LIB=/usr/local/lib
ubuntu: BOOST_INC=/usr/include
ubuntu: BOOST_LIB=/usr/lib/x86_64-linux-gnu
ubuntu: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
ubuntu: IFLAG=-Ilib/OTools -Ilib/ -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
ubuntu: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
ubuntu: LDFLAG=$(LDFLAG_REL)
ubuntu: $(BFILE)

#UBUNTU DEBUG VERSION
ubuntu-dbg: RMATH_INC=$(HOME)/R-3.2.2/src/include
ubuntu-dbg: RMATH_LIB=$(HOME)/R-3.2.2/src/nmath/standalone
ubuntu-dbg: HTSLD_INC=/usr/local/include/
ubuntu-dbg: HTSLD_LIB=/usr/local/lib
ubuntu-dbg: BOOST_INC=/usr/include
ubuntu-dbg: BOOST_LIB=/usr/lib/x86_64-linux-gnu
ubuntu-dbg: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
ubuntu-dbg: IFLAG=-Ilib/OTools -Ilib/ -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
ubuntu-dbg: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
ubuntu-dbg: LDFLAG=$(LDFLAG_DBG)
ubuntu-dbg: $(BFILE)

#MAC RELEASE VERSION
mac: RMATH_INC=$(HOME)/Libraries/R-3.2.2/src/include
mac: RMATH_LIB=$(HOME)/Libraries/R-3.2.2/src/nmath/standalone
mac: HTSLD_INC=$(HOME)/Libraries/htslib-1.2.1
mac: HTSLD_LIB=$(HOME)/Libraries/htslib-1.2.1
mac: BOOST_INC=/opt/local/include
mac: BOOST_LIB=/opt/local/lib
mac: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
mac: IFLAG=-Ilib/OTools -Ilib/ -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
mac: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams-mt.a $(BOOST_LIB)/libboost_program_options-mt.a
mac: LDFLAG=$(LDFLAG_REL) -L /opt/local/lib
mac: $(BFILE)

#MAC DEBUG VERSION
mac-dbg: RMATH_INC=$(HOME)/Libraries/R-3.2.2/src/include
mac-dbg: RMATH_LIB=$(HOME)/Libraries/R-3.2.2/src/nmath/standalone
mac-dbg: HTSLD_INC=$(HOME)/Libraries/htslib-1.2.1
mac-dbg: HTSLD_LIB=$(HOME)/Libraries/htslib-1.2.1
mac-dbg: BOOST_INC=/opt/local/include
mac-dbg: BOOST_LIB=/opt/local/lib
mac-dbg: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
mac-dbg: IFLAG=-Ilib/OTools -Ilib/ -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
mac-dbg: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams-mt.a $(BOOST_LIB)/libboost_program_options-mt.a
mac-dbg: LDFLAG=$(LDFLAG_DBG) -L /opt/local/lib
mac-dbg: $(BFILE)

#COMPILATION RULES
$(BFILE): $(OFILE)
	$(CXX) $^ $(LIB_FILES) -o $@ $(LIB_FLAGS) $(LDFLAG)

obj/QTLtools.o: src/QTLtools.cpp $(HFILE) $(TFILE)
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

obj/quan_%.o: quan_%.cpp quan_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/ase_%.o: ase_%.cpp ase_data.h src/common/data.h src/common/filter.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/bamstat_%.o: bamstat_%.cpp bamstat_data.h src/common/data.h src/common/filter.h $(TFILE)
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
		
