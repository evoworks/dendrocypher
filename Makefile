#CXX = /usr/bin/clang 
#CXX = /usr/bin/g++
#CXX = /opt/local/bin/g++-mp-4.6
CXX = /opt/local/bin/g++-apple-4.2

#INC = -I./ -I/usr/local/include -I/home/jrm/git/proteus/optare/opt/include
INC = -I./ 

#LIB = -L/usr/local/lib -L/home/jrm/git/proteus/optare/opt/lib
LIB = -L/usr/local/lib

#CXXFLAGS = -DHAVE_STD -DHAVE_NAMESPACES $(INC) $(LIB) -Wall -Wl,-rpath,/usr/local/lib/gcc46
#CXXFLAGS = $(INC) $(LIB) -Wall -Wl,-rpath,/opt/local/lib/gcc46
#CXXFLAGS = $(INC) $(LIB) -Wall -Wl,-rpath,/usr/local/lib/gcc46
CXXFLAGS = $(INC) $(LIB) -Wall -Wl,-rpath,/usr/local/lib/apple-gcc42



FFLAGS = 
ifeq ($(MAKECMDGOALS), debug)
	CXXFLAGS += -O0 -DDEBUG -g -ggdb
	FFLAGS += -O0 -DDEBUG -ggdb
else 
	#CXXFLAGS += -O4 -funroll-loops -fomit-frame-pointer -finline-functions
	#CXXFLAGS += -Ofast -funroll-loops -fomit-frame-pointer -finline-functions
	CXXFLAGS += -O3 -g
	FFLAGS += -O3
endif


SRCS = \
dendrology/forestry.cpp \
dendrology/node.cpp \
dendrology/tree.cpp \
matrices/eigen.cpp \
matrices/int_matrix.cpp \
matrices/matrix.cpp \
matrices/str_matrix.cpp \
main.cpp \
tools.cpp \

OBJS = $(SRCS:.cpp=.o)
#OBJS += $(FSRCS:.f=.o)
DEPS = $(SRCS:.cpp=.d)
#DEPS += $(FSRCS:.f=.d)
PROG = DendroCypher


$(PROG):	$(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) -lm

%.d: %.cpp %.f
	@set -e; rm -f $@; \
	$(CXX) -MM $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

all: $(PROG)
debug: $(PROG)

-include $(DEPS)

.PHONY: clean
clean:
	rm -f $(DEPS) $(OBJS) $(PROG) $(PROG).core
