# Update the following lines to your system requirements
AURYNPATH=$(HOME)/auryn

# Path to Auryn include files
AURYNINC=$(AURYNPATH)/src

# Path to Auryn library
AURYNLIB=$(AURYNPATH)/build/home


# The following should not require updating in most cases 
CC = mpic++
CFLAGS=-ansi -pipe -O3 -march=native -ffast-math -pedantic -I/usr/include -I$(AURYNINC)
LDFLAGS=-lboost_program_options -lboost_mpi -lboost_serialization -L$(AURYNLIB) -lauryn

SIMFILES=sim_rc_p10c sim_rc_p10d sim_rc_p11

simfiles: $(SIMFILES)

sim_%: sim_%.o 
	$(CC) $(CFLAGS) $< -L. -lauryn $(LDFLAGS) -o $(subst .o,,$<)

%.o : %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *~ *.o core a.out *.log *.a $(SIMFILES)

