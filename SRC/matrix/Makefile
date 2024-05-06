#  Written: fmk 
#  Created: 10/97
# updated Boris Jeremic@ucdavis.edu Aug2001

include ../../Makefile.def

OBJS       = ID.o Vector.o Matrix.o

################### TARGETS ########################
all: $(OBJS) 

test: $(OBJS) crossProduct.o
	$(LINKER) crossProduct.o $(OBJS) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(FE_LIBRARY) \
	-o matrix_tst

pi: $(OBJS) pi.o
	$(LINKER) pi.o $(OBJS) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(FE_LIBRARY)  \
	-o pi_tst

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core test

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core test
	@$(RM) $(RMFLAGS) test.\$\$\$ 
	@$(RM) $(RMFLAGS) *.ti 
	@$(RM) $(RMFLAGS) tca.map 
	@$(RM) $(RMFLAGS) .inslog*

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
