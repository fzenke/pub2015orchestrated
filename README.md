# Auryn classes and simulations to reproduce on-line learning with orchestrated plasticity rules.

These simulations are to reproduce the findings from Zenke, F., Agnes, E.J., and Gerstner, W. (2015). 
Diverse synaptic plasticity mechanisms orchestrated to form and retrieve memories in spiking neural networks. Nat Commun 6. 
http://www.nature.com/ncomms/2015/150421/ncomms7922/full/ncomms7922.html

A similiar version to this document is available at
https://www.fzenke.net/auryn/doku.php?id=examples:orchestrated_plasticity

## Auryn classes
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

To run the orchestrated plasticity simulation as described in Figure 3 of you
need to install and compile Auryn v0.5 (the develop branch of this repository 
might be tracking more recent versions of Auryn). 
In the following we will assume that you have git installed and up and running
on your system. Moreover, you have all dependencies to compile Auryn installed.

To download and compile Auryn do the following:
```
$ cd ~
$ git clone https://github.com/fzenke/auryn.git
$ cd auryn/build/home
$ make
```

Should you have difficulties compiling the simulator please refer to the
installation and troubleshooting section in the manual (www.fzenke.net/auryn).

Now go to the installation directory of the `src/` simulation code (when you are
reading this, chances are you are already in this directory). And run make
there. This should build the necessary Auryn libraries that implement plasticity
and the simulation libraries. For instance the binary file sim_rc_p10c is the
one behind Figure 3 in Zenke, F., Agnes, E.J., Gerstner, W., 2015. Diverse
synaptic plasticity mechanisms orchestrated to form and retrieve memories in
spiking neural networks. Nature Communications.

Finally update the output path in globalvars.sh to point to a portion of disk
with sufficient space. If you want to run distributed simulations make sure this
path is accessible from all nodes of your cluster.

Invoking
`$ ./run_orchestrated_stdp.sh`
will run the the three scripts 1run_init.sh, 2run_learn.sh and 3run_cued.sh in
sequence which corresponds to the two-fold learning and recall protocol shown in
Figure 3.
