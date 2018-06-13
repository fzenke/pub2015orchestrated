/* 
* Copyright 2015 Friedemann Zenke
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
#include "GlobalPFConnection.h"
#include "P11Connection.h"

#define N_EXEC_WEIGHTS 800

using namespace auryn;

namespace po = boost::program_options;
namespace mpi = boost::mpi;

int main(int ac, char* av[]) 
{

	string dir = "/lcncluster/zenke/reset/";
	string file_prefix = "rc";
	string infilename = "";

	char strbuf [255];
	string msg;

	bool save = false;
	bool chain = false;
	bool prime = false;
	bool consolidation = true;
	bool isp_active = true;

	bool inh_input = false;
	bool noisy_initial_weights = false;
	bool consolidate_initial_weights = false;
	bool quiet = false;

	NeuronID size = 4096;
	NeuronID seed = 1;
	double alpha = 3;
	double kappa = 10;
	double tauf = 200e-3;
	double ujump = 0.2;
	double taud = 200e-3;
	double eta = 1e-3;

	double tauh = 600.0;

	double beta = 5.0e-2;
	double delta = 2.0e-2;
	double weight_a = 0.1;
	double weight_c = 0.5;
	double adapt = 0.0;

	double pot_strength = 0.1; 
	
	double ontime = 1.0;
	double offtime = 5.0;

	double scale = 35;

	double wmax = 5.0;
	double wmin = 0.0;
	double wmaxi = 5.0;
	double wmini = 0.0;

	double bgrate = 10.0;

	int preferred = -1;

	string stimfile = ""; // stimulus patterns file
	string prefile = ""; // preload patters file
	string recfile = ""; // file with receptive fields
	string monfile = ""; // patternsto monitor file
	string premonfile = ""; // patternsto monitor file


	AurynWeight wee = 0.1;
	AurynWeight wei = 0.2;
	AurynWeight wie = 0.2;
	AurynWeight wii = 0.2;

	AurynWeight chi = 1.0;
	AurynWeight xi = 0.5;

	double sparseness = 0.1;
	double sparseness_ext = 0.05;
	double wext = 0.2;
	
	double simtime = 3600.;

	int errcode = 0;

    try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("load", po::value<string>(), "input weight matrix")
            ("prefix", po::value<string>(), "set file prefix")
            ("save", "save network state at end of sim")
            ("chain", "chain mode for pattern loader")
            ("prime", "prime network with a burn-in phase")
            ("noconsolidation", "switches off consolidation")
            ("noisp", "switches off isp")
            ("noisyweights", "switches noisy initial weights on")
            ("consolidateweights", "initialize weights as consolidated")
            ("inhinput", "switches external input to inh on")
            ("quiet", "quiet mode")
            ("alpha", po::value<double>(), "exc input rate")
            ("kappa", po::value<double>(), "hom parameter")
            ("taud", po::value<double>(), "time constant of synaptic depression")
            ("tauf", po::value<double>(), "time constant of synaptic facilitation")
            ("tauh", po::value<double>(), "time constant of homeostasis")
            ("ujump", po::value<double>(), "u jump STP constant")
            ("chi", po::value<double>(), "chi factor - pattern preload strength")
            ("xi", po::value<double>(), "xi factor - stimulation strength")
            ("wext", po::value<double>(), "recurrent weight (ext)")
            ("wee", po::value<double>(), "recurrent weight (wee)")
            ("wei", po::value<double>(), "recurrent weight (wei)")
            ("wii", po::value<double>(), "recurrent weight (wii)")
            ("wie", po::value<double>(), "recurrent weight (wie)")
            ("extsparse", po::value<double>(), "external sparseness")
            ("intsparse", po::value<double>(), "external sparseness")
            ("simtime", po::value<double>(), "simulation time")
            ("ontime", po::value<double>(), "simulation time")
            ("offtime", po::value<double>(), "simulation time")
            ("dir", po::value<string>(), "output dir")
            ("eta", po::value<double>(), "the learning rate")
            ("beta", po::value<double>(), "decay parameter")
            ("potstrength", po::value<double>(), "potential strength parameter")
            ("delta", po::value<double>(), "growth parameter")
            ("weight_a", po::value<double>(), "weight_a")
            ("weight_c", po::value<double>(), "weight_c")
            ("size", po::value<int>(), "simulation size")
            ("seed", po::value<int>(), "random seed ")
            ("stimfile", po::value<string>(), "stimulus file")
            ("prefile", po::value<string>(), "preload file")
            ("recfile", po::value<string>(), "receptive field file")
            ("scale", po::value<double>(), "stimulus strength")
            ("adapt", po::value<double>(), "adaptation jump size for long time constant")
            ("bgrate", po::value<double>(), "background rate of input")
            ("preferred", po::value<int>(), "num of preferred stim")
            ("monfile", po::value<string>(), "monitor file")
            ("premonfile", po::value<string>(), "presynaptic monitor file")
        ;

        po::variables_map vm;        
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
			std::cout << desc << "\n";
            return 1;
        }

        if (vm.count("load")) {
			std::cout << "load from matrix " 
                 << vm["load"].as<string>() << ".\n";
			infilename = vm["load"].as<string>();
        } 

        if (vm.count("prefix")) {
			std::cout << "prefix from matrix " 
                 << vm["prefix"].as<string>() << ".\n";
			file_prefix = vm["prefix"].as<string>();
        } 

        if (vm.count("save")) {
			save = true;
        } 

        if (vm.count("chain")) {
			chain = true;
        } 

        if (vm.count("prime")) {
			prime = true;
        } 

        if (vm.count("noconsolidation")) {
			consolidation = false;
        } 

        if (vm.count("noisp")) {
			isp_active = false;
        } 

        if (vm.count("noisyweights")) {
			noisy_initial_weights = true;
        } 

        if (vm.count("consolidateweights")) {
			consolidate_initial_weights = true;
        } 

        if (vm.count("inhinput")) {
			inh_input = true;
        } 

        if (vm.count("quiet")) {
			quiet = true;
        } 
        if (vm.count("alpha")) {
			std::cout << "alpha set to " 
                 << vm["alpha"].as<double>() << ".\n";
			alpha = vm["alpha"].as<double>();
        } 

        if (vm.count("kappa")) {
			std::cout << "kappa set to " 
                 << vm["kappa"].as<double>() << ".\n";
			kappa = vm["kappa"].as<double>();
        } 

        if (vm.count("taud")) {
			std::cout << "taud set to " 
                 << vm["taud"].as<double>() << ".\n";
			taud = vm["taud"].as<double>();
        } 

        if (vm.count("tauf")) {
			std::cout << "tauf set to " 
                 << vm["tauf"].as<double>() << ".\n";
			tauf = vm["tauf"].as<double>();
        } 

        if (vm.count("tauh")) {
			std::cout << "tauh set to " 
                 << vm["tauh"].as<double>() << ".\n";
			tauh = vm["tauh"].as<double>();
        } 

        if (vm.count("ujump")) {
			std::cout << "ujump set to " 
                 << vm["ujump"].as<double>() << ".\n";
			ujump = vm["ujump"].as<double>();
        } 

        if (vm.count("simtime")) {
			std::cout << "simtime set to " 
                 << vm["simtime"].as<double>() << ".\n";
			simtime = vm["simtime"].as<double>();
        } 

        if (vm.count("ontime")) {
			std::cout << "ontime set to " 
                 << vm["ontime"].as<double>() << ".\n";
			ontime = vm["ontime"].as<double>();
        } 

        if (vm.count("offtime")) {
			std::cout << "offtime set to " 
                 << vm["offtime"].as<double>() << ".\n";
			offtime = vm["offtime"].as<double>();
        } 

        if (vm.count("dir")) {
			std::cout << "dir set to " 
                 << vm["dir"].as<string>() << ".\n";
			dir = vm["dir"].as<string>();
        } 

        if (vm.count("chi")) {
			std::cout << "chi set to " 
                 << vm["chi"].as<double>() << ".\n";
			chi = vm["chi"].as<double>();
        } 

        if (vm.count("xi")) {
			std::cout << "xi set to " 
                 << vm["xi"].as<double>() << ".\n";
			xi = vm["xi"].as<double>();
        } 

        if (vm.count("wext")) {
			std::cout << "wext set to " 
                 << vm["wext"].as<double>() << ".\n";
			wext = vm["wext"].as<double>();
        } 

        if (vm.count("wee")) {
			std::cout << "wee set to " 
                 << vm["wee"].as<double>() << ".\n";
			wee = vm["wee"].as<double>();
        } 

        if (vm.count("wei")) {
			std::cout << "wei set to " 
                 << vm["wei"].as<double>() << ".\n";
			wei = vm["wei"].as<double>();
        } 

        if (vm.count("wii")) {
			std::cout << "wii set to " 
                 << vm["wii"].as<double>() << ".\n";
			wii = vm["wii"].as<double>();
        } 

        if (vm.count("wie")) {
			std::cout << "wie set to " 
                 << vm["wie"].as<double>() << ".\n";
			wie = vm["wie"].as<double>();
        } 

        if (vm.count("extsparse")) {
			std::cout << "extsparse set to " 
                 << vm["extsparse"].as<double>() << ".\n";
			sparseness_ext = vm["extsparse"].as<double>();
        } 

        if (vm.count("intsparse")) {
			std::cout << "intsparse set to " 
                 << vm["intsparse"].as<double>() << ".\n";
			sparseness = vm["intsparse"].as<double>();
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

        if (vm.count("potstrength")) {
			std::cout << "potstrength set to " 
                 << vm["potstrength"].as<double>() << ".\n";
			pot_strength = vm["potstrength"].as<double>();
        } 

        if (vm.count("delta")) {
			std::cout << "delta set to " 
                 << vm["delta"].as<double>() << ".\n";
			delta = vm["delta"].as<double>();
        } 

        if (vm.count("weight_a")) {
			std::cout << "weight_a set to " 
                 << vm["weight_a"].as<double>() << ".\n";
			weight_a = vm["weight_a"].as<double>();
        } 

        if (vm.count("weight_c")) {
			std::cout << "weight_c set to " 
                 << vm["weight_c"].as<double>() << ".\n";
			weight_c = vm["weight_c"].as<double>();
        } 

        if (vm.count("size")) {
			std::cout << "size set to " 
                 << vm["size"].as<int>() << ".\n";
			size = vm["size"].as<int>();
        } 

        if (vm.count("stimfile")) {
			std::cout << "stimfile set to " 
                 << vm["stimfile"].as<string>() << ", ";
			stimfile = vm["stimfile"].as<string>();
			monfile = stimfile;
			premonfile = stimfile;
        } 

        if (vm.count("prefile")) {
			std::cout << "prefile set to " 
                 << vm["prefile"].as<string>() << ", ";
			prefile = vm["prefile"].as<string>();
        } 

        if (vm.count("recfile")) {
			std::cout << "recfile set to " 
                 << vm["recfile"].as<string>() << ", ";
			recfile = vm["recfile"].as<string>();
        } 

        if (vm.count("scale")) {
			std::cout << "scale set to " 
                 << vm["scale"].as<double>() << ".\n";
			scale = vm["scale"].as<double>();
        } 

        if (vm.count("adapt")) {
			std::cout << "adapt set to " 
                 << vm["adapt"].as<double>() << ".\n";
			adapt = vm["adapt"].as<double>();
        } 

        if (vm.count("bgrate")) {
			std::cout << "bgrate set to " 
                 << vm["bgrate"].as<double>() << ".\n";
			bgrate = vm["bgrate"].as<double>();
        } 

        if (vm.count("preferred")) {
			std::cout << "preferred set to " 
                 << vm["preferred"].as<int>() << ".\n";
			preferred = vm["preferred"].as<int>();
        } 

        if (vm.count("monfile")) {
			std::cout << "monfile set to " 
                 << vm["monfile"].as<string>() << ", ";
			monfile = vm["monfile"].as<string>();
        } 

        if (vm.count("premonfile")) {
			premonfile = vm["premonfile"].as<string>();
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


	auryn_init( ac, av, dir, "sim_rc_p11", file_prefix );
	sys->set_master_seed(42);
	logger->set_logfile_loglevel(VERBOSE);
	
	//log params
	logger->parameter("alpha",alpha);
	logger->parameter("beta",beta);
	logger->parameter("delta",delta);
	logger->parameter("eta",eta);
	logger->parameter("wee",wee);
	logger->parameter("wext",wext);
	logger->parameter("chi",chi);
	logger->parameter("xi",xi);

	logger->parameter("stimfile",stimfile);
	logger->parameter("monfile",monfile);
	logger->parameter("offtime",offtime);
	logger->parameter("ontime",ontime);

	logger->parameter("taud",taud);
	logger->parameter("tauf",tauf);
	logger->parameter("ujump",ujump);

	logger->parameter("tauh",tauh);

	AIF2Group * neurons_e = new AIF2Group(size);
	// AIFSclGroup * neurons_e = new AIFSclGroup(size);
	// neurons_e->randomize_state_vector_gauss("scaling_weight",1.0,0.2,2351);
	// neurons_e->set_eta(eta); // exc scaling
	// neurons_e->set_kappa(alpha); // exc scaling


	neurons_e->dg_adapt1  = 0.1;
	neurons_e->dg_adapt2  = adapt;
	
	// AIFGroup * neurons_e = new AIFGroup(size);
	// neurons_e->dg_adapt1  = 0.1;
	// neurons_e->random_adapt(1,1);

	// IFGroup * neurons_e = new IFGroup(size);

	// IF2Group * neurons_e = new IF2Group(size);
	neurons_e->set_tau_ampa(5e-3);
	neurons_e->set_tau_gaba(10e-3);
	neurons_e->set_tau_nmda(100e-3);
	neurons_e->set_ampa_nmda_ratio(0.2);

	IFGroup * neurons_i2 = new IFGroup(size/4);
	neurons_i2->set_tau_ampa(5e-3);
	neurons_i2->set_tau_gaba(10e-3);
	neurons_i2->set_tau_nmda(100e-3);
	neurons_i2->set_ampa_nmda_ratio(0.3);


	StimulusGroup * stimgroup;
	sprintf(strbuf, "%s/%s.%d.stimtimes", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	string stimtimefile = strbuf;
	// stimgroup is initialized here. If stimfile is empty no patterns are loaded
	// and it acts simply as PoissonGroup
	stimgroup = new StimulusGroup(size,stimtimefile);
	stimgroup->set_mean_on_period(ontime);
	stimgroup->set_mean_off_period(offtime);
	stimgroup->binary_patterns = true;
	stimgroup->scale = scale;
	stimgroup->background_rate = bgrate;
	stimgroup->background_during_stimulus = true;
	// stimgroup->randomintensities = true;
	if (seed!=1) stimgroup->seed(seed);





	double raw_delta = delta*eta/1e-3;

	P11Connection * con_ee;
	con_ee = new P11Connection(neurons_e,neurons_e,
			wee,sparseness,
			eta,
			kappa,
			wmax
			);
	con_ee->set_transmitter(AMPA);
	con_ee->set_name("EE");
	con_ee->set_weight_a(weight_a); 
	con_ee->set_weight_c(weight_c); 
	con_ee->consolidation_active = consolidation;
	double wtmax = 1.0/4*(weight_c-weight_a);
	double normalization_factor = (wtmax-weight_a)*(wtmax-(weight_a+weight_c)/2)*(wtmax-weight_c); 
	con_ee->pot_strength = pot_strength/normalization_factor;
	logger->parameter("normalized pot_strength",con_ee->pot_strength);
	if ( noisy_initial_weights ) 
		con_ee->random_data(wee,wee);
	if ( consolidate_initial_weights )
		con_ee->consolidate();
	// con_ee->random_data(wee,wee);
	// con_ee->consolidate();
	// STP parameters 
	con_ee->set_tau_d(taud);
	con_ee->set_tau_f(tauf);
	con_ee->set_ujump(ujump);
	// con_ee->set_urest(ujump);
	con_ee->set_beta(beta);
	con_ee->delta = raw_delta*eta;
	con_ee->set_min_weight(wmin);
	con_ee->set_tau_hom(tauh);
	// con_ee->sparse_set_data(0.01,weight_c);
	// con_ee->w_min = 0.01;
	//
	





	// SparseConnection * con_ei2 = new SparseConnection(neurons_e,neurons_i2,wei,sparseness,GLUT);
	// SparseConnection * con_ei2 = new SparseConnection(neurons_e,neurons_i2,wei,sparseness,GLUT);
	STPConnection * con_ei2 = new STPConnection(neurons_e,neurons_i2,3*wei,sparseness,GLUT);
	con_ei2->set_tau_d(taud);
	con_ei2->set_tau_f(0.6);
	con_ei2->set_ujump(0.2);
	// con_ei2->set_urest(0.2);

	// P11Connection * con_ei2 = NULL;
	// con_ei2 = new P11Connection( neurons_e, neurons_i2,
	// 	3*wei,sparseness,
	// 	eta,
	// 	kappa, // supposedly deprecated
	// 	wmax,
	// 	GLUT
	// 	);

	// con_ei2->set_weight_a(weight_a);
	// con_ei2->set_weight_c(weight_c);
	// con_ei2->set_tau_d(taud);
	// con_ei2->set_tau_f(tauf);
	// con_ei2->set_ujump(ujump);
	// con_ei2->set_urest(ujump);
	// con_ei2->set_beta(beta);
	// con_ei2->A2_minus = eta;
	// con_ei2->delta = raw_delta*eta;
	// con_ei2->set_min_weight(wmin);
	// con_ei2->set_name("E->I");

	// P11Connection * con_ei2 = new P11Connection(neurons_e,neurons_i2,
	// 		wei,sparseness,
	// 		eta,
	// 		alpha,
	// 		wmax
	// 		);
	// con_ei2->set_name("EI2");
	// con_ei2->set_tau_d(taud);
	// con_ei2->set_tau_f(tauf);
	// con_ei2->set_ujump(0.2);
	// con_ei2->set_urest(0.2);
	// con_ei2->set_beta(beta);
	// con_ei2->weight_a = weight_a;


	double geta = -eta*1e-4;
	SparseConnection * con_i2i2 = new SparseConnection(neurons_i2,neurons_i2,wii,sparseness,GABA);
	// RateModulatedConnection * con_i2i2;
	// con_i2i2 = new RateModulatedConnection(neurons_i2,neurons_i2,
	// 		wii,sparseness,
	// 		GABA
	// 		);
	// con_i2i2->eta = geta;
	// con_i2i2->rate_target   = 3*alpha;
	// con_i2i2->rate_estimate = 3*alpha;
	// con_i2i2->random_data(wii,wii);
	con_i2i2->set_name("I2->I2");

	GlobalPFConnection * con_i2e;
	con_i2e = new GlobalPFConnection(neurons_i2,neurons_e,
			wie,sparseness,
			10.0,
			eta/50,
			alpha, 
			wmaxi,
			GABA
			);
	con_i2e->set_name("I2E");
	// con_i2e->set_eta(geta);
	// con_i2e->rate_target = alpha;
	// con_i2e->rate_estimate = alpha;
	// con_i2e->random_data(wie,wie);

	if ( inh_input ) {
		STPConnection * con_stim_i 
			= new STPConnection( stimgroup, 
					neurons_i2, 
					wext, 
					sparseness_ext, 
					GLUT);
		con_stim_i->set_tau_d(taud);
		con_stim_i->set_tau_f(0.6);
		con_stim_i->set_ujump(0.2);
		// con_stim_i->set_urest(0.2);
	}

	// External input
	P11Connection * con_stim_e = NULL;
	con_stim_e = new P11Connection( stimgroup, neurons_e,
		wext,sparseness_ext,
		eta,
		kappa, // supposedly deprecated
		wmax,
		GLUT
		);


	con_stim_e->set_weight_a(weight_a);
	con_stim_e->set_weight_c(weight_c);
	con_stim_e->set_tau_d(taud);
	con_stim_e->set_tau_f(tauf);
	con_stim_e->set_ujump(ujump);
	// con_stim_e->set_urest(ujump);
	con_stim_e->set_beta(beta);
	con_stim_e->delta = raw_delta*eta;
	con_stim_e->set_min_weight(wmin);
	// con_stim_e->random_data(wext,wext);
	if ( noisy_initial_weights ) 
		con_stim_e->random_data(wext,wext);
	con_stim_e->set_name("Stim->E");
	con_stim_e->consolidation_active = consolidation;
	con_stim_e->pot_strength = pot_strength/normalization_factor;
	con_stim_e->set_tau_hom(tauh);
	if ( consolidate_initial_weights )
		con_stim_e->consolidate();
	// con_stim_e->sparse_set_data(0.01,weight_c);
	// con_stim_e->consolidate();

	
	// External input
	// P11Connection * con_stim_i = NULL;
	// con_stim_i = new P11Connection( stimgroup, neurons_i2,
	// 	wext,sparseness_ext,
	// 	eta,
	// 	kappa, // supposedly deprecated
	// 	wmax,
	// 	GLUT
	// 	);

	// con_stim_i->set_weight_a(weight_a);
	// con_stim_i->set_weight_c(weight_c);
	// con_stim_i->random_data(wext,wext);
	// // con_stim_i->sparse_set_data(0.1,1.0);
	// con_stim_i->set_tau_d(taud);
	// con_stim_i->set_tau_f(tauf);
	// con_stim_i->set_ujump(ujump);
	// con_stim_i->set_urest(ujump);
	// con_stim_i->set_beta(beta);
	// con_stim_i->delta = raw_delta*eta;
	// con_stim_i->set_min_weight(wmin);
	// // con_stim_e->random_data(wext,wext);
	// con_stim_i->set_name("Stim->E");
	// // con_stim_e->random_data(wext,wext);


	if (!stimfile.empty()) {
		logger->msg("Setting up stimulus ...",PROGRESS,true);
		stimgroup->load_patterns(stimfile.c_str());
		stimgroup->set_next_action_time(10); // let network settle for some time

		sprintf(strbuf, "%s/%s.%d.s.spk", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
		BinarySpikeMonitor * smon_s = new BinarySpikeMonitor( stimgroup, string(strbuf), size );

		// gives the first 3 patterns half of the probability
		if ( preferred > 0 ) { 
			std::vector<double> dist = stimgroup->get_distribution();
			int r = preferred;
			for ( int i = 0 ; i < dist.size() ; ++i ) {
				if ( i == r ) 
					dist[i] = 1.0;
				else
					dist[i] = 1.0/dist.size();
			}
			stimgroup->set_distribution(dist);
		}

		sprintf(strbuf, "%s/%s.%d.wse", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
		new WeightStatsMonitor( con_stim_e, string(strbuf) );
	}


	// load if necessary
	if (!infilename.empty()) {
		logger->msg("Loading from file ...",PROGRESS,true);
		sys->load_network_state(infilename.c_str());
	}


		
	if ( !prefile.empty() && chi > 0.0 ) {
		con_ee->patterns_ignore_gamma = true;
		con_ee->load_patterns(prefile,chi);
		// con_ee->consolidate();
	}

	if ( !prefile.empty() && xi > 0.0 ) {
		con_stim_e->patterns_ignore_gamma = true;
		con_stim_e->load_patterns(prefile,xi,false);
		if ( consolidate_initial_weights )
			con_stim_e->consolidate();
	}

	sprintf(strbuf, "%s/%s.%d.wse", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	new WeightStatsMonitor( con_stim_e, string(strbuf) );

	if ( !recfile.empty() && xi > 0.0 ) {
		con_stim_e->load_fragile_matrix(recfile); // TODO
		con_stim_e->scale_all(xi);
	}

	sprintf(strbuf, "%s/%s.%d.sse", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	WeightMonitor * wmon_s = new WeightMonitor( con_stim_e, string(strbuf), 1.0 ); 
	wmon_s->add_equally_spaced(50);
	if ( !monfile.empty() ) {
		sprintf(strbuf, "%s/%s.%d.pact", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
		PatternMonitor * patmon = new PatternMonitor( neurons_e, string(strbuf) , monfile.c_str(), 100);

		if ( !stimfile.empty() ) // 
			wmon_s->load_pattern_connections(stimfile,monfile,20,20,ASSEMBLIES_ONLY); // true for assemblies only
		else 
			wmon_s->load_pattern_connections(monfile,20,20,ASSEMBLIES_ONLY); // true for assemblies only
	}


	// sprintf(strbuf, "%s/%s.%d.sei", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	// WeightMonitor * wmon_ei = new WeightMonitor( con_ei2, 0, 100, strbuf, 1.0, DATARANGE); 
	
	sprintf(strbuf, "%s/%s.%d.see", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	WeightMonitor * wmon = new WeightMonitor( con_ee, string(strbuf), 1.0); 
	wmon->add_equally_spaced(50);

	if ( !monfile.empty() ) 
		wmon->load_pattern_connections(monfile,10,10,ASSEMBLIES_ONLY); // true for assemblies only

	sprintf(strbuf, "%s/%s.%d.hom", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	StateMonitor * stmon_hom = new StateMonitor( con_ee->hom, 12, string(strbuf), 1 ); 

	sprintf(strbuf, "%s/%s.%d.mem", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	VoltageMonitor * stmon_mem = new VoltageMonitor( neurons_e, 3, string(strbuf) ); 
	stmon_mem->record_for(10); // stops recording after 10s

	sprintf(strbuf, "%s/%s.%d.imem", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	VoltageMonitor * stmon_imem = new VoltageMonitor( neurons_i2, 3, string(strbuf) ); 
	stmon_imem->record_for(10); // stops recording after 10s

	// sprintf(strbuf, "%s/%s.%d.si1e", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	// WeightMonitor * wmon_i1e = new WeightMonitor( con_i1e, string(strbuf) ); 
	// wmon_i1e->add_equally_spaced(50);

	sprintf(strbuf, "%s/%s.%d.si2e", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	WeightMonitor * wmon_i2e = new WeightMonitor( con_i2e, string(strbuf) ); 
	wmon_i2e->add_equally_spaced(50);


	// sprintf(strbuf, "%s/%s.%d.ipe", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	// WeightMonitor * wmon_ise = new WeightMonitor( con_se_inh, 0, 100, strbuf, 1.0, DATARANGE); 

	// sprintf(strbuf, "%s/%s.%d.sii", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	// WeightMonitor * wmon_ii = new WeightMonitor( con_ii, 0, 100, strbuf, 1.0, DATARANGE); 
	
	sprintf(strbuf, "%s/%s.%d.wee", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	new WeightStatsMonitor( con_ee, string(strbuf) );

	sprintf(strbuf, "%s/%s.%d.wi2e", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	new WeightStatsMonitor( con_i2e, string(strbuf) );

	if ( !monfile.empty() ) {
		sprintf(strbuf, "%s/%s.%d.wprec", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
		WeightPatternMonitor * wpmon = new WeightPatternMonitor( con_ee, string(strbuf), 60 );
		wpmon->load_patterns(monfile);
	}

	if ( !premonfile.empty() && !monfile.empty() ) {
		sprintf(strbuf, "%s/%s.%d.wpin", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
		WeightPatternMonitor * wpmon = new WeightPatternMonitor( con_stim_e, string(strbuf), 60 );
		wpmon->load_pre_patterns(premonfile);
		wpmon->load_post_patterns(monfile);
	}

	// sprintf(strbuf, "%s/%s.%d.wi1e", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	// new WeightStatsMonitor( con_i1e, string(strbuf) );


	// sprintf(strbuf, "%s/%s.%d.wpe", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	// new WeightStatsMonitor( con_se, string(strbuf) );

	sprintf(strbuf, "%s/%s.%d.e.spk", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	BinarySpikeMonitor * smon_e = new BinarySpikeMonitor( neurons_e, string(strbuf), size );


	sprintf(strbuf, "%s/%s.%d.i2.spk", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	BinarySpikeMonitor * smon_i2 = new BinarySpikeMonitor( neurons_i2, string(strbuf), size );

	sprintf(strbuf, "%s/%s.%d.e.prate", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	PopulationRateMonitor * pmon_e = new PopulationRateMonitor( neurons_e, string(strbuf), 0.1 );


	sprintf(strbuf, "%s/%s.%d.i2.prate", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	PopulationRateMonitor * pmon_i2 = new PopulationRateMonitor( neurons_i2, string(strbuf), 0.1 );

	RateChecker * chk = new RateChecker( neurons_e , -1 , 20. , 0.1);


	// prime
	if ( prime ) {
		// tries to decrease training time by initially rapidly decorrelating patterns
		// Note, that this was not done in the original publication
		logger->msg("High intensity priming ...",PROGRESS,true);
		stimgroup->set_mean_on_period(1.0);
		stimgroup->set_mean_off_period(0.2);

		if (!sys->run(200.0,false))
			errcode = 1;
	}

	stimgroup->set_mean_on_period(ontime);
	stimgroup->set_mean_off_period(offtime);

	if ( eta > 0 ) {
		con_ee->stdp_active = true;
		con_stim_e->stdp_active = true;
		con_i2e->stdp_active = true;
	} else {
		con_ee->stdp_active = false;
		con_stim_e->stdp_active = false;
		con_i2e->stdp_active = false;
	}

	con_i2e->stdp_active = isp_active;


	logger->msg("Main simtime ...",PROGRESS,true);
	if (!sys->run(simtime,false)) 
			errcode = 1;

	logger->msg("Saving ...",PROGRESS,true);
	if ( save ) {
		sys->set_output_dir(dir);
		sys->save_network_state(file_prefix);
	}

	logger->msg("Writing connectivity matrices ...",PROGRESS,true);
	sprintf(strbuf, "%s/%s.%d.ee.wmat", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	con_ee->write_to_file(strbuf);

	sprintf(strbuf, "%s/%s.%d.ext.wmat", dir.c_str(), file_prefix.c_str(), sys->mpi_rank() );
	con_stim_e->write_to_file(strbuf);

	if (errcode) auryn_abort(errcode);

	logger->msg("Freeing ...",PROGRESS,true);
	auryn_free();

	return errcode;

}

