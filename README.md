# On-line learning with orchestrated plasticity rules

This repository contains example  simulations which reproduce the key results
from our paper as described in Figure 3:

```
Zenke, F., Agnes, E.J., and Gerstner, W. (2015). Diverse synaptic plasticity
mechanisms orchestrated to form and retrieve memories in spiking neural
networks. Nat Commun 6.
http://www.nature.com/ncomms/2015/150421/ncomms7922/full/ncomms7922.html
```


To run the code you need to install and compile Auryn. The present version of
the coded uses slightly different paramters which lead to shorter training
times. Note, that due to ongoing development of the simulator and specifically
changes to the random seeding of SparseConnections since v0.4.1 the network you
are simulating is not identical at the individual connection level to the
simulation shown in the paper. To replicate the original simulation bit-by-bit
please refer to the commit history of this repository and the simulator.


This development version was last tested with Auryn v0.8.1-dev (commit
94892708c3).  


## Included Auryn classes

```
P10Connection 
       Implements the excitatory connection object for orchestrated
       plasticity without homeostatic sliding threshold.

P11Connection 
       Implements the excitatory connection object for orchestrated
       plasticity with homeostatic sliding threshold.

P12Connection 
       Implements variation of P11Connection with two limits of the
       sliding threshold.
```



## Running an example

In the following I will assume that you have git installed and up and running
on your system. Moreover, you have all dependencies to compile Auryn installed.
First, download and compile Auryn (see https://github.com/fzenke/auryn for
instructions). Should you have difficulties compiling the simulator please
refer to the installation and troubleshooting section in the manual
(https://fzenke.net/auryn).

Now go to the installation directory of the simulation code (when you are
reading this, chances are you are already in this directory) and run `make`
there (you might have to update the auryn path in the Makefile if you are using
a different install directory). This will build the necessary Auryn libraries
that implement plasticity and the simulation libraries. 

Finally, update the output path in 'globalvars.sh' to point to a portion of
disk with sufficient space. If you want to run distributed simulations, make
sure this path is accessible from all nodes of your cluster.

Invoking `$ ./run_orchestrated_stdp.sh` will run the three scripts
1run_init.sh, 2run_learn.sh and 3run_cued.sh. To regenerate the file containing
the responsive cells, run `bootstrap.sh`. 

The simulation file `sim_rc_p10c` is the downstream version of the simulation
behind Figure 3 in our paper.  The simulation `sim_rc_p11` is a version with a
homeostatic sliding threshold. The current default  is to run the homeostatic
one, but this can be changed easily in `globalvars.sh`. 


## Example analysis

You will find some example analysis files (including a Jupyter notebook) in the
`ana` subdirectory.
