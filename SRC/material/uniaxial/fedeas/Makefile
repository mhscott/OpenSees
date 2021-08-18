include ../../../../Makefile.def

OBJS       = FedeasHardeningMaterial.o Hard1.o \
	FedeasBond1Material.o \
	FedeasBond2Material.o \
	FedeasConcr1Material.o \
	FedeasConcr2Material.o Concr2.o \
	FedeasConcr3Material.o \
	FedeasHyster1Material.o Hyster1.o \
	FedeasHyster2Material.o Hyster2.o \
	FedeasSteel1Material.o Steel1.o \
	FedeasSteel2Material.o Steel2f.o \
	PD1.o PlasticDamageMaterial.o \
	TclFedeasMaterialCommand.o

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
