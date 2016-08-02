# Please edit Makefile for your system
SHELL = /bin/sh
CC    = g++ 
F77   = gfortran
#CXX   = mpicxx
#CXX   = /Users/saksena/gcc-4.8.1/bin/g++
CXX   = g++
#F90   = mpif90
F90   = gfortran
MAKE  = make --no-print-directory
#LIBS  = -lmpi
#MPIRUN= mpirun
#NP_ARG= -np
LD_FLAGS=

all: ZTM generateDataset

#CFLAGS = -O0 -g -fopenmp -I/Users/saksena/boost_1_53_0/include -L/Users/saksena/boost_1_53_0/lib -lboost_regex
#CXXFLAGS = -O0 -g -fopenmp -I/Users/saksena/boost_1_53_0/include -L/Users/saksena/boost_1_53_0/lib -lboost_regex

# CFLAGS = -O2 -fopenmp -I/Users/saksena/boost_1_53_0/include -L/Users/saksena/boost_1_53_0/lib -lboost_regex
# CXXFLAGS = -O2 -fopenmp -I/Users/saksena/boost_1_53_0/include -L/Users/saksena/boost_1_53_0/lib -lboost_regex

# CFLAGS = -O2 -fopenmp -I/home/insong/boost_1_53_0/include -L/home/insong/boost_1_53_0/lib -lboost_regex
# CXXFLAGS = -O2 -fopenmp -I/home/insong/boost_1_53_0/include -L/home/insong/boost_1_53_0/lib -lboost_regex

## use this in della
CFLAGS = -O2 -fopenmp -I/home/insong/boost_1_53_0/include -L/home/insong/boost_1_53_0/lib -lboost_regex
CXXFLAGS = -O2 -fopenmp -I/home/insong/boost_1_53_0/include -L/home/insong/boost_1_53_0/lib -lboost_regex

# FLAGS = -O2 -fopenmp -lboost_regex
# CXXFLAGS = -O2 -fopenmp -lboost_regex

.SUFFIXES: .o .c .f .cc .cpp .f90
.c:
	${CC} ${CFLAGS} -o $* $< ${LIBS}
.c.o:
	${CC} ${CFLAGS} -c $<
.f:
	${F77} ${FFLAGS} -o $* $< ${LIBS}
.f.o:
	${F77} ${FFLAGS} -c $<
.f90:
	${F90} ${F90FLAGS} -o $* $< ${LIBS}
.cc:
	${CXX} ${CXXFLAGS} -o $* $< ${LIBS}
.cpp:
	${CXX} ${CXXFLAGS} -o $* $< ${LIBS}

OBJS = Analytics.o ChkpBase.o ZeroTradeModel.o ZeroTradeModelIO.o main.o
OBJSF = Analytics.o ChkpBase.o ZeroTradeModel.o ZeroTradeModelIO.o mainF.o
OBJS2 = generateDataset.o

DEFAULT: ZTM ZTMF generateDataset
clean:
	rm -f *.o
	rm -f *~
	rm -f core
	rm -f ZeroTradeModel
	rm -f ZeroTradeModelIO
	rm -f ZTM

ZTM: $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(LD_FLAGS) -fopenmp

ZTMF: $(OBJSF)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJSF) $(LD_FLAGS) -fopenmp

generateDataset: $(OBJS2)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS2) $(LD_FLAGS) -fopenmp

#sources:
#	@echo ${SOURCES}
