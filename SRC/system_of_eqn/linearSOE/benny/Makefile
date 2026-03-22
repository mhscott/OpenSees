include ../../../../Makefile.def

OBJS       = BennySparseLinSOE.o BennySparseLinSolver.o \
	benny_sparse.o benny_symbolic.o benny_numeric.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS)

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
