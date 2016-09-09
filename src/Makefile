# Please edit Makefile for your system
SHELL = /bin/sh
CC    = g++ 
CXX   = g++
LIBS  = -m64 -shared -L/usr/lib64/R/lib -lR
LD_FLAGS= 

all: ZTM dynCluster.so

FLAGS = -O2 -fopenmp -lboost_regex
CXXFLAGS = -DNDEBUG -fopenmp -lboost_regex -fpic -O2 -g -pipe -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector --param=ssp-buffer-size=4 -m64 -mtune=generic
RCPPFLAGS = -I"/home/hubertj/R/x86_64-redhat-linux-gnu-library/3.3/Rcpp/include" -I/usr/include/R -L/usr/lib64/R/lib -lR

.SUFFIXES: .o .c .cc .cpp 
.c:
	${CC} ${CFLAGS} -o $* $< ${LIBS}
.cc:
	${CXX} ${CXXFLAGS} -o $* $< ${LIBS}
.cpp:
	${CXX} ${CXXFLAGS} -o $* $< ${LIBS}

OBJS = Analytics.o ChkpBase.o ZeroTradeModel.o ZeroTradeModelIO.o
OBJ1 = main.o
OBJ2 = mainRcpp.o

DEFAULT: ZTM dynCluster.so
clean:
	rm -f *.o 
	rm -f *.so
	rm -f *~
	rm -f core
	rm -f ZeroTradeModel
	rm -f ZeroTradeModelIO
	rm -f ZTM

ZTM: $(OBJS) $(OBJ1)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(OBJ1) $(LD_FLAGS)

mainRcpp.o: 
	$(CXX) $(RCPPFLAGS) $(CXXFLAGS) -c mainRcpp.cpp -o mainRcpp.o ${LIBS}

dynCluster.so: $(OBJS) $(OBJ2)
	$(CXX) $(FLAGS) $(OBJS) $(OBJ2) -o dynCluster.so ${LIBS}
	cp -f dynCluster.so ../libs/.
