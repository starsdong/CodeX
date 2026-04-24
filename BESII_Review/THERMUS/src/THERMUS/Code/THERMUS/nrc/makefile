#########################
ARCH          = $(shell root-config --arch)

CXX           =
ObjSuf        = o
SrcSuf        = c
ExeSuf        =
ifeq ($(ARCH),macosx)
DllSuf        = dylib
else
DllSuf        = so
endif
OutPutOpt     = -o 

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)


CXX           = g++
CXXFLAGS      = -O  -fPIC
LD            = g++
LDFLAGS       = -O  
ifeq ($(ARCH),macosx)
SOFLAGS       = -dynamiclib -single_module -undefined dynamic_lookup
else
SOFLAGS       = -shared
endif


CXXFLAGS     += $(ROOTCFLAGS)
CXXFLAGS     += -I$(THERMUS)/main -I$(THERMUS)/include -I$(THERMUS)/nrc
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------


FNCSO        	= broydn.o \
		  zbrent.o \
		  fdjac.o \
		  fmin.o \
		  rsolv.o \
		  qrdcmp.o \
		  qrupdt.o \
		  lnsrch.o \
		  nrutil.o \
		  rotate.o \

FNCSS        	= broydn.c \
		  zbrent.c \
		  fdjac.c \
		  fmin.c \
		  rsolv.c \
		  qrdcmp.c \
		  qrupdt.c \
		  lnsrch.c \
		  nrutil.c \
		  rotate.c \

FNCSSO		= $(THERMUS)/lib/libNRCFunctions.$(DllSuf)

OBJS          	= $(FNCSO) 
PROGRAMS      	= $(FNCSSO)
 
#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:   Functions

ifeq ($(ARCH),macosx)
all:            $(PROGRAMS)
		@ln -sf $(THERMUS)/lib/libNRCFunctions.$(DllSuf) \
		       $(THERMUS)/lib/libNRCFunctions.so
else
all:            $(PROGRAMS)
endif


clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.root *.ps *.so .def so_locations
		@rm -rf cxx_repository

###

Functions:	$(FNCSSO)	

$(FNCSSO):	$(FNCSO)
	$(LD) $(SOFLAGS) $(LIBS) $(LDFLAGS) $^ $(OutPutOpt) $@


$(FNCSO): 	$(FNCSS)
	$(CXX) $(CXXFLAGS) -c $(FNCSS) 

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
