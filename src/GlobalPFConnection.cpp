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

#include "GlobalPFConnection.h"

using namespace auryn;

void GlobalPFConnection::init(AurynFloat tau_hom, AurynFloat eta, AurynFloat kappa, AurynFloat maxweight)
{
	if ( dst->get_post_size() == 0 ) return;

	tau_post = tau_hom;
	post_factor_mul = exp(-auryn_timestep/tau_post);

	target_rate = kappa;
	expected_spikes = kappa*dst->get_size();
	post_factor = expected_spikes;

	learning_rate = eta/expected_spikes;

	zeta = 0.0;
	xi = 1.0;

	double tau_stdp = 20e-3;
	tr_pre = src->get_pre_trace(tau_stdp);
	tr_post = dst->get_post_trace(tau_stdp);

	set_min_weight(0.0);
	set_max_weight(maxweight);

	stdp_active = true;

}

void GlobalPFConnection::init_shortcuts() 
{
	if ( dst->get_post_size() == 0 ) return; // if there are no target neurons on this rank

	fwd_ind = w->get_row_begin(0); 
	fwd_data = w->get_data_begin();

	bkw_ind = bkw->get_row_begin(0); 
	bkw_data = bkw->get_data_begin();
}

void GlobalPFConnection::finalize() {
	DuplexConnection::finalize();
	init_shortcuts();
}

void GlobalPFConnection::free()
{
}

GlobalPFConnection::GlobalPFConnection(SpikingGroup * source, NeuronGroup * destination, TransmitterType transmitter) : DuplexConnection(source, destination, transmitter)
{
}

GlobalPFConnection::GlobalPFConnection(SpikingGroup * source, NeuronGroup * destination, 
		const char * filename, 
		AurynFloat tau_hom, 
		AurynFloat eta, 
		AurynFloat kappa, AurynFloat maxweight , 
		TransmitterType transmitter) 
: DuplexConnection(source, 
		destination, 
		filename, 
		transmitter)
{
	init(tau_hom, eta, kappa, maxweight);
	init_shortcuts();
}

GlobalPFConnection::GlobalPFConnection(SpikingGroup * source, NeuronGroup * destination, 
		AurynWeight weight, AurynFloat sparseness, 
		AurynFloat tau_hom, 
		AurynFloat eta, 
		AurynFloat kappa, AurynFloat maxweight , 
		TransmitterType transmitter,
		string name) 
: DuplexConnection(source, 
		destination, 
		weight, 
		sparseness, 
		transmitter, 
		name)
{
	init(tau_hom, eta, kappa, maxweight);
	if ( name.empty() ) 
		set_name("GlobalPFConnection");
	init_shortcuts();
}

GlobalPFConnection::~GlobalPFConnection()
{
	if ( dst->get_post_size() > 0 ) 
		free();
}

inline AurynWeight GlobalPFConnection::dw_pre(NeuronID post)
{
	double dw = (tr_post->get(post)+xi);
	return dw;
}

inline AurynWeight GlobalPFConnection::dw_post(NeuronID pre)
{
	double dw = tr_pre->get(pre);
	return dw;
}



void GlobalPFConnection::propagate_forward()
{
	AurynDouble global_modulation = learning_rate*(post_factor-expected_spikes);
	for (SpikeContainer::const_iterator spike = src->get_spikes()->begin() ; // spike = pre_spike
			spike != src->get_spikes()->end() ; ++spike ) {
		for (NeuronID * c = w->get_row_begin(*spike) ; c != w->get_row_end(*spike) ; ++c ) { // c = post index
			AurynWeight * value = &(fwd_data[c-fwd_ind]); 
			transmit( *c , *value );
			if ( stdp_active ) {
				NeuronID translated_spike = dst->global2rank(*c);
			    *value += global_modulation*dw_pre(translated_spike) - zeta;
			    // *value += learning_rate*(dw_het);
			    if ( *value < get_min_weight() ) 
					*value = get_min_weight();
			}
		}
	}
}

inline void GlobalPFConnection::propagate_backward()
{
	if ( !stdp_active ) return;
	AurynDouble global_modulation = learning_rate*(post_factor-expected_spikes);
	NeuronID * ind = bkw->get_row_begin(0); // first element of index array
	AurynWeight ** data = bkw->get_data_begin();
	SpikeContainer::const_iterator spikes_end = dst->get_spikes_immediate()->end();
	for (SpikeContainer::const_iterator spike = dst->get_spikes_immediate()->begin() ; // spike = post_spike
			spike != spikes_end ; ++spike ) {
		for (NeuronID * c = bkw->get_row_begin(*spike) ; c != bkw->get_row_end(*spike) ; ++c ) {
			AurynWeight * value = data[c-ind]; 
			*value += global_modulation*dw_post(*c);
			if ( *value > get_max_weight() ) 
				*value = get_max_weight();
		}
	}
}


void GlobalPFConnection::propagate()
{
	propagate_forward();
	propagate_backward();
}

void GlobalPFConnection::evolve()
{
	post_factor *= post_factor_mul;
	post_factor += dst->get_spikes()->size()/tau_post;
}

bool GlobalPFConnection::write_to_file(string filename)
{

	std::stringstream oss;
	oss << filename << ".cstate";

	std::ofstream outfile;
	outfile.open(oss.str().c_str(),std::ios::out);
	if (!outfile) {
		std::cerr << "Can't open output file " << filename << std::endl;
	  throw AurynOpenFileException();
	}

	boost::archive::text_oarchive oa(outfile);
	oa << post_factor ;

	outfile.close();

	return SparseConnection::write_to_file(filename);
}

bool GlobalPFConnection::load_from_file(string filename)
{

	std::stringstream oss;
	oss << filename << ".cstate";
	std::ifstream infile (oss.str().c_str());

	if (!infile) {
		std::stringstream oes;
		oes << "Can't open input file " << filename;
		logger->msg(oes.str(),ERROR);
		throw AurynOpenFileException();
	}

	boost::archive::text_iarchive ia(infile);
	ia >> post_factor;
	logger->parameter("loaded: post_factor",post_factor);

	infile.close();

	return SparseConnection::load_from_file(filename);
}

