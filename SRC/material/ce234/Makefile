include ../../../Makefile.def

OBJS       = CE234Model01.o CE234Model02.o \
	CE234Model03.o CE234PlaneStress.o \
	TclCE234MaterialCommand.o

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
