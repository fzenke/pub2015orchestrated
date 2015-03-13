/* 
* Copyright 2014 Friedemann Zenke
*
* This file is part of Auryn, a simulation package for plastic
* spiking neural networks.
* 
* Auryn is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* Auryn is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with Auryn.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "auryn.h"
#include "P10Connection.h"

#define N_EXEC_WEIGHTS 600

using namespace std;

namespace po = boost::program_options;
namespace mpi = boost::mpi;

int main(int ac, char* av[]) 
{

	string dir = "./";
	string file_prefix = "testbench";
	string wmat = "";

	char strbuf [255];
	string msg;

	NeuronID size = 800;
	const NeuronID n_neurons = 100;

	NeuronID seed = 1;
	double kappa = 5;
	double poststim = 0.0;

	double wmax = 10.0;

	bool record_rates = true;


	double eta = 1e-3;
	double beta = 0.05;
	double weightc = 0.5;
	double sparseness = 0.1;
	AurynWeight we = 0.3;
	AurynWeight wctl = 0.1;

	
	double simtime = 100.;

	int errcode = 0;

    try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("kappa", po::value<double>(), "exc input rate")
            ("poststim", po::value<double>(), "postsynaptic stim")
            ("we", po::value<double>(), "input weight (exc)")
            ("simtime", po::value<double>(), "simulation time")
            ("eta", po::value<double>(), "the learning rate")
            ("beta", po::value<double>(), "the reset rate")
            ("weightc", po::value<double>(), "weightc value")
            ("size", po::value<int>(), "simulation size")
            ("rates", po::value<bool>(), "record single unit rates")
            ("wmat", po::value<string>(), "wmat to load")
            ("prefix", po::value<string>(), "file prefix")
            ("dir", po::value<string>(), "file dir")
        ;

        po::variables_map vm;        
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }


        if (vm.count("kappa")) {
            cout << "kappa set to " 
                 << vm["kappa"].as<double>() << ".\n";
			kappa = vm["kappa"].as<double>();
        } 

        if (vm.count("poststim")) {
            cout << "poststim set to " 
                 << vm["poststim"].as<double>() << ".\n";
			poststim = vm["poststim"].as<double>();
        } 

        if (vm.count("simtime")) {
            cout << "simtime set to " 
                 << vm["simtime"].as<double>() << ".\n";
			simtime = vm["simtime"].as<double>();
        } 

        if (vm.count("we")) {
            cout << "we set to " 
                 << vm["we"].as<double>() << ".\n";
			we = vm["we"].as<double>();
        } 

        if (vm.count("eta")) {
            cout << "eta set to " 
                 << vm["eta"].as<double>() << ".\n";
			eta = vm["eta"].as<double>();
        } 

        if (vm.count("beta")) {
            cout << "beta set to " 
                 << vm["beta"].as<double>() << ".\n";
			beta = vm["beta"].as<double>();
        } 

        if (vm.count("weightc")) {
            cout << "weightc set to " 
                 << vm["weightc"].as<double>() << ".\n";
			weightc = vm["weightc"].as<double>();
        } 

        if (vm.count("size")) {
            cout << "size set to " 
                 << vm["size"].as<int>() << ".\n";
			size = vm["size"].as<int>();
        } 

        if (vm.count("rates")) {
            cout << "rates set to " 
                 << vm["rates"].as<bool>() << ".\n";
			record_rates = vm["rates"].as<bool>();
        } 

        if (vm.count("wmat")) {
            cout << "wmat set to " 
                 << vm["wmat"].as<string>() << ".\n";
			wmat = vm["wmat"].as<string>();
        } 

        if (vm.count("prefix")) {
            cout << "prefix set to " 
                 << vm["prefix"].as<string>() << ".\n";
			file_prefix = vm["prefix"].as<string>();
        } 

        if (vm.count("dir")) {
            cout << "dir set to " 
                 << vm["dir"].as<string>() << ".\n";
			dir = vm["dir"].as<string>();
        } 

        if (vm.count("seed")) {
            cout << "seed set to " 
                 << vm["seed"].as<int>() << ".\n";
			seed = vm["seed"].as<int>();
        } 
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

	// BEGIN Global stuff
	mpi::environment env(ac, av);
	mpi::communicator world;
	communicator = &world;

	sprintf(strbuf, "%s/%s.%d.log", dir.c_str(), file_prefix.c_str(), world.rank());
	string logfile = strbuf;
	logger = new Logger(logfile,world.rank(),PROGRESS,EVERYTHING);

	sys = new System(&world);
	// END Global stuff
	

	sys->quiet = true;


	PoissonGroup * poisson_e = new PoissonGroup(size, kappa);

	PoissonGroup * poisson_ctl = new PoissonGroup(size, 1.0);

	PoissonGroup * poisson_post = new PoissonGroup(n_neurons, poststim);

	IFGroup * neurons = new IFGroup(n_neurons);

	if ( poststim > 0.0 ) {
		IdentityConnection * con_post = new IdentityConnection(poisson_post,neurons,1.0);
		con_post->set_transmitter(MEM);
	} 

	P10Connection * con_e = new P10Connection(
			poisson_e,
			neurons,
			we,
			sparseness,
			eta,
			1.0, // not used
			wmax
			);

	con_e->set_min_weight(0.01);
	con_e->set_max_weight(1000);
	con_e->set_tau_d(0.2);
	con_e->set_tau_f(0.6);
	con_e->set_ujump(0.2);
	con_e->set_urest(0.2);
	con_e->set_beta(beta);
	con_e->set_weight_a(0.0);
	con_e->set_weight_c(weightc);

	P10Connection * con_ctl = new P10Connection(
			poisson_ctl,
			neurons,
			0.0,
			sparseness,
			eta,
			1.0, // not used
			wmax
			);

	con_ctl->set_min_weight(0.01);
	con_ctl->set_max_weight(1000);
	con_ctl->set_tau_d(0.2);
	con_ctl->set_tau_f(0.6);
	con_ctl->set_ujump(0.2);
	con_ctl->set_urest(0.2);
	con_ctl->set_beta(beta);
	con_ctl->set_weight_a(0.0);
	con_ctl->set_weight_c(weightc);

	sprintf(strbuf, "%s/%s.%d.ras", dir.c_str(), file_prefix.c_str(), world.rank() );
	SpikeMonitor * smon_e = new SpikeMonitor( neurons, strbuf, size);

	// sprintf(strbuf, "%s/%s.%d.p.ras", dir.c_str(), file_prefix.c_str(), world.rank() );
	// SpikeMonitor * smon_p = new SpikeMonitor( poisson_e, strbuf, size);

	// sprintf(strbuf, "%s/%s.%d.mem", dir.c_str(), file_prefix.c_str(), world.rank() );
	// VoltageMonitor * vmon = new VoltageMonitor( neurons, 0, strbuf, 10 );

	// sprintf(strbuf, "%s/%s.%d.ampa", dir.c_str(), file_prefix.c_str(), world.rank() );
	// AmpaMonitor * ampamon = new AmpaMonitor( neurons, 0, strbuf, 10 );

	// sprintf(strbuf, "%s/%s.%d.gaba", dir.c_str(), file_prefix.c_str(), world.rank() );
	// GabaMonitor * gabamon = new GabaMonitor( neurons, 0, strbuf, 10 );

	sprintf(strbuf, "%s/%s.%d.prate", dir.c_str(), file_prefix.c_str(), world.rank() );
	PopulationRateMonitor * pmon = new PopulationRateMonitor( neurons, strbuf, 10 );

	sprintf(strbuf, "%s/%s.%d.rates", dir.c_str(), file_prefix.c_str(), world.rank() );
	RateMonitor * rmon = new RateMonitor( neurons, strbuf, 2 );

	sprintf(strbuf, "%s/%s.%d.we", dir.c_str(), file_prefix.c_str(), world.rank() );
	WeightMonitor * wmone = new WeightMonitor( con_e, strbuf );
	for ( int i = 0 ; i < 100 ; ++i )
		wmone->add_to_list(i,0);

	sprintf(strbuf, "%s/%s.%d.wectl", dir.c_str(), file_prefix.c_str(), world.rank() );
	WeightMonitor * wmonctl = new WeightMonitor( con_ctl, strbuf );
	for ( int i = 0 ; i < 100 ; ++i )
		wmonctl->add_to_list(i,0);

	RateChecker * chk = new RateChecker( neurons , -1 , 20. , 1);

	if ( !wmat.empty() ) {
		logger->msg("Loading weight matrix ...",PROGRESS,true);
		con_e->load_from_complete_file(wmat);
		con_e->set_all(we);

		con_ctl->load_from_complete_file(wmat);
		con_ctl->set_all(we);
		con_ctl->set_all(wctl);
	}

	con_e->stdp_active = false;
	con_ctl->stdp_active = false;

	if (!sys->run(60,false)) 
			errcode = 1;

	con_e->stdp_active = true;
	con_ctl->stdp_active = true;

	if (!sys->run(simtime,false)) 
			errcode = 1;

	// sprintf(strbuf, "%s/%s.%d.wmat", dir.c_str(), file_prefix.c_str(), world.rank() );
	// con_e->write_to_file(strbuf);
	sprintf(strbuf, "%s/%s", dir.c_str(), file_prefix.c_str() );
	sys->save_network_state(strbuf);

	logger->msg("Freeing ...",PROGRESS,true);
	delete sys;

	if (errcode)
		env.abort(errcode);
	return errcode;
}
