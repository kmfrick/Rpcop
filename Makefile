SHELL = /bin/sh

SRCS = espai_rba_f.cpp fileops.cpp ll_flt.cpp ll_p.cpp ll_pnt.cpp ll_q.cpp ma.cpp main.cpp mb.cpp pila.cpp
OBJS=$(subst .cc,.o,$(SRCS))

CXXFLAG = -Wall -g
CXX = g++
INCLUDE =
LIBS = -lm 

pcop: $(OBJS)
	${CXX} $(LDFLAGS) -o pcop $(OBJS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<
