#---------------------------------------------------------------------
# Detect system type using uname command.
#
# On Linux,    uname -s returns 'Linux'
# On Linux 64, uname -m additionally returns x86_64
# On Mac OS X, uname -s returns 'Darwin'
#
# See discussion at http://osdir.com/ml/gnu.make.bugs/2002-09/msg00003.html
# Conclusion: instead of using OSTYPE (env var defined in bash), manually
# set OSTYPE := $(shell uname -msr)
# In order to avoid confusion, I'm going to just use OSNAME := $(shell uname -s).
#---------------------------------------------------------------------

OSNAME := $(shell uname -s)
MACHINE := $(shell uname -m)

ifeq ($(OSNAME),Linux)
  ifeq ($(MACHINE),x86_64)
    PLATFORM := LINUX64
  else
    PLATFORM := LINUX
  endif
else 
  ifeq ($(OSNAME),Darwin)
    PLATFORM := MAC_OS_X
  else
    PLATFORM := OTHER
  endif
endif

#---------------------------------------------------------------------
# Choose a compiler & its options
#---------------------------------------------------------------------
g++ = g++
CXX = g++

OPTS = -std=c++11 -pipe -O3 -MMD -Wno-write-strings  -Wno-unused-result -Wno-deprecated

#-------------------------------------------------------------------
# Base parameters for INCLUDE and LIB.
#
# The next sections after this one add items to INCLUDE and LIB using
# the += operator. 
#-------------------------------------------------------------------

INCLUDE = -I.
LIB = $(LIBS)


#--------------------------------------------------------------------
# Add Xlib and OpenGL
#--------------------------------------------------------------------

ifneq ($(GL_ON),0)

	ifeq ($(PLATFORM),MAC_OS_X)
		OPTS += -DMACOS
   		GL_LIB = -framework OpenGL -framework GLUT -framework Foundation -lGLEW.osx
		LIB += $(GL_LIB) -L./GL 
	else
		GL_LIB = -lglut  -lGL -lGLU -lGLEW.linux64
		LIB += $(GL_LIB) -L./GL 
		INCLUDE += -I/usr/include
	endif
endif


#--------------------------------------------------------------------

DIRS= . mathtool GL modelgraph

SRCS=$(wildcard $(addsuffix /*.cpp,$(DIRS))) 
CSRCS=$(wildcard $(addsuffix /*.c,$(DIRS)))
OBJS=${SRCS:.cpp=.o} $(CSRCS:.c=.o)

#--------------------------------------------------------------------
# Set CFLAGS
#--------------------------------------------------------------------
INCLUDE    +=   $(addprefix -I,$(DIRS))
CFLAGS     = $(OPTS) $(INCLUDE) 
CXXFLAGS   = $(CFLAGS)

TARGET = raytracer 

all: $(TARGET)

#--------------------------------------------------------------------
$(TARGET): $(OBJS) 
	${CXX} ${CXXFLAGS} -o $@ $(OBJS) $(LIB)

clean:
	rm -f *.o $(OBJS) Dependencies $(TARGET)

#--------------------------------------------------------------------
# Build Rules
#--------------------------------------------------------------------

.SUFFIXES: .cpp .c

.cpp.o: 
	${CXX} ${CXXFLAGS} -c $< -o $@
	cat $*.d >> Dependencies
	-rm -f $*.d

.c.o: 
	${CXX} ${CXXFLAGS} -c $< -o $@ 
	cat $*.d >> Dependencies
	-rm -f $*.d

Dependencies:
	touch Dependencies

-include Dependencies

