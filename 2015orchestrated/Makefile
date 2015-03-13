AURYNPATH=$(HOME)/auryn
AURYNINC=$(AURYNPATH)/src
AURYNLIB=$(AURYNPATH)/build/home

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

