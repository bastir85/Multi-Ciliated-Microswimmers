RM=rm -f
DEFAULTFLAGS=-O3 -std=c++0x -std=c++11 -Wall -ffast-math -I libcilia/\
		 -DCILIA_PHASE_POWER=1600.0 -DSIM_MEASUREINTERVAL=50 -g

MPCDFLAGS=-DOMP_PARALLEL_MDOBJ -DPARALLEL  -fopenmp
MDFLAGS=  -DSIM_TAU=0.05 -DSIM_SMALL_TAU=0.0001 -DBD_SIMULATION -DBD_INTEGRATOR #-DOSEEN

CPPFLAGS=$(DEFAULTFLAGS) $(MPCDFLAGS)
##CPPFLAGS=$(DEFAULTFLAGS) $(MDFLAGS)

#-DPARALLEL -fopenmp#-DMPCD_RANDOM
#-DXY_WALL -DBOUNDS -DMPCD_CURVATURE=10

LDFLAGS=-fopenmp -g #-openmp-link dynamic
LDLIBS=-lm -lconfig++

SRCS=main.cpp simulation.cpp cilia.cpp libcilia/fluid.cpp libcilia/mpcd.cpp\
	 sphere.cpp libcilia/mdobject.cpp libcilia/globalforces.cpp clami.cpp ConstantBeatCilia.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: cilia

cilia: $(OBJS)
	$(CXX) $(LDFLAGS) -o cilia $(OBJS) $(LDLIBS)

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)
	$(RM) cilia

dist-clean: clean
	$(RM) *~ .dependtool

include .depend


