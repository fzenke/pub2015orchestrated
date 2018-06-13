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


#define N_EXEC_WEIGHTS 800

using namespace auryn;

namespace po = boost::program_options;
namespace mpi = boost::mpi;

int main(int ac, char* av[]) 
{

	string dir = "./";
	string file_prefix = "p10f";
	string wmat = "";

	char strbuf [255];
	string msg;

	NeuronID size = 1000;
	const NeuronID n_neurons = 1;

	NeuronID seed = 1;
	double kappa = 5;

	double on = 5.0;
	double off = 20.0;

	double wmax = 10.0;


	double eta = 1e-3;
	double beta = 0.05;
	double sparseness = 0.1;
	AurynWeight we = 0.3;
	
	double simtime = 100.;

	int errcode = 0;

    try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("kappa", po::value<double>(), "exc input rate")
            ("we", po::value<double>(), "input weight (exc)")
            ("simtime", po::value<double>(), "simulation time")
            ("eta", po::value<double>(), "the learning rate")
            ("beta", po::value<double>(), "the reset rate")
            ("on", po::value<double>(), "mean on time")
            ("off", po::value<double>(), "mean off time")
            ("size", po::value<int>(), "simulation size")
            ("wmat", po::value<string>(), "wmat to load")
            ("prefix", po::value<string>(), "file prefix")
            ("dir", po::value<string>(), "file dir")
        ;

        po::variables_map vm;        
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
			std::cout << desc << "\n";
            return 1;
        }


        if (vm.count("kappa")) {
			std::cout << "kappa set to " 
                 << vm["kappa"].as<double>() << ".\n";
			kappa = vm["kappa"].as<double>();
        } 

        if (vm.count("simtime")) {
			std::cout << "simtime set to " 
                 << vm["simtime"].as<double>() << ".\n";
			simtime = vm["simtime"].as<double>();
        } 

        if (vm.count("we")) {
			std::cout << "we set to " 
                 << vm["we"].as<double>() << ".\n";
			we = vm["we"].as<double>();
        } 

        if (vm.count("eta")) {
			std::cout << "eta set to " 
                 << vm["eta"].as<double>() << ".\n";
			eta = vm["eta"].as<double>();
        } 

        if (vm.count("beta")) {
			std::cout << "beta set to " 
                 << vm["beta"].as<double>() << ".\n";
			beta = vm["beta"].as<double>();
        } 

        if (vm.count("on")) {
			std::cout << "on set to " 
                 << vm["on"].as<double>() << ".\n";
			on = vm["on"].as<double>();
        } 

        if (vm.count("off")) {
			std::cout << "off set to " 
                 << vm["off"].as<double>() << ".\n";
			off = vm["off"].as<double>();
        } 

        if (vm.count("size")) {
			std::cout << "size set to " 
                 << vm["size"].as<int>() << ".\n";
			size = vm["size"].as<int>();
        } 

        if (vm.count("wmat")) {
			std::cout << "wmat set to " 
                 << vm["wmat"].as<string>() << ".\n";
			wmat = vm["wmat"].as<string>();
        } 

        if (vm.count("prefix")) {
			std::cout << "prefix set to " 
                 << vm["prefix"].as<string>() << ".\n";
			file_prefix = vm["prefix"].as<string>();
        } 

        if (vm.count("dir")) {
			std::cout << "dir set to " 
                 << vm["dir"].as<string>() << ".\n";
			dir = vm["dir"].as<string>();
        } 

        if (vm.count("seed")) {
			std::cout << "seed set to " 
                 << vm["seed"].as<int>() << ".\n";
			seed = vm["seed"].as<int>();
        } 
    }
    catch(std::exception& e) {
		std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
		std::cerr << "Exception of unknown type!\n";
    }

	auryn_init(ac, av, dir);
	

	IFGroup * neurons = new IFGroup(n_neurons);

	sprintf(strbuf, "%s/%s.%d.stimtimes", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	string stimtimefile = strbuf;

	double width = 0.05;
	double ampl  = 20;
	MovingBumpGroup * stimgroup = new MovingBumpGroup(size,on,width,ampl,stimtimefile);
	stimgroup->set_floor(0.05);
	// stimgroup->set_mean_on_period(on);
	// stimgroup->set_mean_off_period(off);
	// stimgroup->scale = 20;
	// stimgroup->set_next_action_time(100); // let network settle for some time
	// stimgroup->background_rate = kappa;
	// stimgroup->background_during_stimulus = true;
	// stimgroup->load_patterns("/home/zenke/stim/1k/gauss10xsigma50.pat");


	P10Connection * con_e = new P10Connection(
			stimgroup,
			neurons,
			0.0,
			1.0,
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
	con_e->set_weight_c(0.3);
	con_e->set_all(we);

	sprintf(strbuf, "%s/%s.%d.prate", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	PopulationRateMonitor * pmon = new PopulationRateMonitor( neurons, strbuf, 1 );

	// sprintf(strbuf, "%s/%s.%d.rates", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	// RateMonitor * rmon = new RateMonitor( neurons, strbuf, 2 );
	
	sprintf(strbuf, "%s/%s.%d.ras", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	SpikeMonitor * rmon = new SpikeMonitor( stimgroup, strbuf );

	sprintf(strbuf, "%s/%s.%d.syn", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	WeightMonitor * wmone = new WeightMonitor( con_e, strbuf );
	wmone->add_equally_spaced(100);

	RateChecker * chk = new RateChecker( neurons , -1 , 20. , 1);

	if ( !wmat.empty() ) {
		logger->msg("Loading weight matrix ...",PROGRESS,true);
		con_e->load_from_file(wmat);
		con_e->set_all(we);
	}

	if (!sys->run(simtime,false)) 
			errcode = 1;

	// sprintf(strbuf, "%s/%s.%d.wmat", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	// con_e->write_to_file(strbuf);
	// sprintf(strbuf, "%s/%s", dir.c_str(), file_prefix.c_str() );
	// sys->save_network_state(strbuf);

	if (errcode)
		auryn_abort(errcode);

	logger->msg("Freeing ...",PROGRESS,true);
	auryn_free();
	return errcode;
}
