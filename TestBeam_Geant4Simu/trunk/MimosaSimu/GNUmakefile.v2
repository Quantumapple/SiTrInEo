# $Id: GNUmakefile 68058 2013-03-13 14:47:43Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := MimosaSimu
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

EXTRALIBS := $(DIGI_DIR)/lib/digit_b2.o

#EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --glibs)
#CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags) 


include $(G4INSTALL)/config/binmake.gmk

# Root (exlude libNew and libpthread from library list)
ROOTINC       = -DTIARA_USEROOT -I $(ROOTSYS)/include/ -I$(DIGI_DIR)/src/

ROOTLIBS      = -L$(ROOTSYS)/lib  -lCore -lCint -lMathCore -lRIO -lNet -lHist -lGraf -lGraf3d  -lTree -lRint  -lMatrix -lPhysics -lThread -pthread -lMinuit -lHtml  -lGeom -lEG -lGui -lGed -lTreePlayer  -lFoam  -lPostscript  -lGpad  -lm -ldl -rdynamic

# Extra flags for G4 
CPPFLAGS += -m64 -fPIC $(ROOTINC)
LDLIBS   += -m64 $(ROOTLIBS)

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

