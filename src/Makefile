# Update the following lines to your system requirements
AURYNPATH=$(HOME)/auryn

# Path to Auryn include files
AURYNINC=$(AURYNPATH)/src

# Path to Auryn library
AURYNLIB=$(AURYNPATH)/build/release/src



# The following should not require updating in most cases 
CC = mpic++
CFLAGS=-ansi -pipe -O3 -march=native -ffast-math -pedantic -I/usr/include -I$(AURYNINC)
LDFLAGS=-lboost_program_options -lboost_mpi -lboost_serialization -lboost_system -lboost_filesystem -L$(AURYNLIB) -lauryn

SIMFILES=sim_rc_p10c sim_rc_p10d sim_rc_p11 sim_rc_p12
OBJFILES=P10Connection.o P11Connection.o P12Connection.o GlobalPFConnection.o

simfiles: $(SIMFILES)

sim_%: sim_%.o $(OBJFILES) 
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $(subst .o,,$<)

%.o : %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *~ *.o core a.out *.log *.a $(SIMFILES)

