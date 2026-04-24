#########################
ARCH          = $(shell root-config --arch)

CXX           =
ObjSuf        = o
SrcSuf        = cxx
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


SEP_CLASSESH	= $(THERMUS)/main/TTMParticle.h \
		  $(THERMUS)/main/TTMDecay.h \
		  $(THERMUS)/main/TTMDecayChannel.h \
		  $(THERMUS)/main/TTMParticleSet.h \
		  $(THERMUS)/main/TTMParameter.h \
		  $(THERMUS)/main/TTMParameterSet.h \
		  $(THERMUS)/main/TTMParameterSetBSQ.h \
		  $(THERMUS)/main/TTMParameterSetBQ.h \
		  $(THERMUS)/main/TTMParameterSetCanBSQ.h \
		  $(THERMUS)/main/TThermalParticle.h \
		  $(THERMUS)/main/TThermalParticleBSQ.h \
		  $(THERMUS)/main/TThermalParticleBQ.h \
		  $(THERMUS)/main/TThermalParticleCanBSQ.h \
		  $(THERMUS)/main/TThermalModel.h \
		  $(THERMUS)/main/TTMYield.h \
		  $(THERMUS)/main/TTMDensObj.h \
		  $(THERMUS)/main/TThermalModelBSQ.h \
		  $(THERMUS)/main/TThermalModelBQ.h \
		  $(THERMUS)/main/TThermalModelCanBSQ.h \
		  $(THERMUS)/main/TThermalFit.h \
		  $(THERMUS)/main/TThermalFitBSQ.h \
		  $(THERMUS)/main/TThermalFitBQ.h \
		  $(THERMUS)/main/TThermalFitCanBSQ.h \
		  $(THERMUS)/main/TTMIDObj.h \

FNCSO        	= BQConstrainQ.o \
		  BSQConstrainS.o \
		  BSQConstrainQ.o \
		  BSQConstrainSQ.o \
		  BSQConstrainSQEN.o \
		  BSQConstrainSEN.o \
		  BSQConstrainSQST3.o \
		  BSQConstrainSST3.o \
		  BSQConstrainSQPercolation.o \
		  BSQConstrainSPercolation.o \
		  BQConstrainEN.o \
		  BQConstrainQEN.o \
		  BQConstrainST3.o \
   		  BQConstrainBDens.o \
		  BQConstrainQST3.o \
	  	  BQConstrainQBDens.o \
	  	  BQConstrainQNetBDens.o \
		  BQConstrainQPercolation.o \
		  BSQConstrainSBDens.o \
		  BSQConstrainSQBDens.o \
		  BSQConstrainSQNetBDens.o \
		  BQfuncEN.o \
                  BQfuncQEN.o \
                  BQfuncST3.o \
                  BQfuncBDens.o \
                  BQfuncQST3.o \
                  BQfuncQBDens.o \
                  BQfuncQNetBDens.o \
		  BQfuncQ.o \
		  BQfuncQPercolation.o \
		  BSQfuncSQPercolation.o \
		  BSQfuncSPercolation.o \
	          BSQfuncS.o \
	          BSQfuncQ.o \
		  BSQfuncSQ.o \
		  BSQfuncSQEN.o \
		  BSQfuncSEN.o \
		  BSQfuncSQST3.o \
		  BSQfuncSST3.o \
		  BSQfuncSBDens.o \
		  BSQfuncSQBDens.o \
		  BSQfuncSQNetBDens.o \
		  fit_function.o \
		  Minuit_function.o \
		  Functions.o \
		  BSQCanonicalFncs.o \
		  IntegrationRoutines.o \
		  IntegrationRoutines2D.o \
		  Int_2_String.o \
		  FindExclVolPressure.o \
		  ExclVolPressureFunc.o \

FNCSS        	= BQConstrainQ.c \
		  BSQConstrainS.c \
		  BSQConstrainQ.c \
		  BSQConstrainSQ.c \
		  BSQConstrainSQEN.c \
		  BSQConstrainSEN.c \
		  BSQConstrainSQST3.c \
		  BSQConstrainSST3.c \
		  BSQConstrainSQPercolation.c \
		  BSQConstrainSPercolation.c \
		  BQConstrainEN.c \
                  BQConstrainQEN.c \
                  BQConstrainST3.c \
                  BQConstrainBDens.c \
                  BQConstrainQST3.c \
                  BQConstrainQBDens.c \
                  BQConstrainQNetBDens.c \
		  BQConstrainQPercolation.c \
		  BSQConstrainSBDens.c \
		  BSQConstrainSQBDens.c \
		  BSQConstrainSQNetBDens.c \
		  BQfuncEN.c \
		  BQfuncQEN.c \
		  BQfuncST3.c \
		  BQfuncBDens.c \
		  BQfuncQST3.c \
		  BQfuncQBDens.c \
		  BQfuncQNetBDens.c \
		  BQfuncQ.c \
		  BQfuncQPercolation.c \
	          BSQfuncS.c \
	          BSQfuncQ.c \
		  BSQfuncSQ.c \
		  BSQfuncSQEN.c \
		  BSQfuncSEN.c \
		  BSQfuncSQST3.c \
		  BSQfuncSST3.c \
		  BSQfuncSBDens.c \
		  BSQfuncSQBDens.c \
		  BSQfuncSQNetBDens.c \
		  BSQfuncSQPercolation.c \
		  BSQfuncSPercolation.c \
		  fit_function.c \
		  Minuit_function.c \
		  Functions.c \
		  BSQCanonicalFncs.c \
		  IntegrationRoutines.c \
		  IntegrationRoutines2D.c \
		  Int_2_String.c \
		  FindExclVolPressure.c \
		  ExclVolPressureFunc.c \

FNCSSO		= $(THERMUS)/lib/libFunctions.$(DllSuf)

OBJS          	= $(FNCSO) 
PROGRAMS      	= $(FNCSSO)
 
#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:   Functions

ifeq ($(ARCH),macosx)
all:            $(PROGRAMS)
		@ln -sf $(THERMUS)/lib/libFunctions.$(DllSuf) \
		        $(THERMUS)/lib/libFunctions.so
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


$(FNCSO): 	$(FNCSS) $(SEP_CLASSESH) 
	$(CXX) $(CXXFLAGS) -c $(FNCSS) 

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
