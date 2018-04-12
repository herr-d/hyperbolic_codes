#ifndef BASIC_INCLUDES_HPP
#define BASIC_INCLUDES_HPP

#include <cassert>
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <set>
#include <vector>
#include <map>
#include <utility>
#include <unordered_map>

using size_type = uint32_t;
using ParityCheck = std::set<size_type>;
using Graph = std::map<size_type,ParityCheck>;
using StabilizerContainer = std::vector<ParityCheck>;


/*
 * This struct stores all the information for the simulation.
 * An instance of this object called "code" needs to be defined globally (I know that's not pretty)
 * such that autotune has access to the information.
 */
struct Code_info{
	// basic simulation info
	bool boot;
	char *dir;
	uint16_t r; // first schaeffli symbol
	uint16_t s; // second schaeffli symbol
	uint16_t l; // number of layers
	uint16_t t_check; // frequency of t_checks
	uint16_t t_delete;
	size_type big_t_max; // maximum numbers of error correction rounds
	double probability; // probability for errors
	size_type seed0;
	size_type seed1;
	size_type max_num_X; // maximum logical X errors
	size_type max_num_Z; // maximum logical Z errors

	// generate code structure

	StabilizerContainer X_stabilizer, Z_stabilizer; // list of stabilizers
	StabilizerContainer X_boundary, Z_boundary; // stabilizer boundary
	size_type num_qubits; // total number of qubits
	StabilizerContainer X_operator, Z_operator; //operator chains for logical qubits (not necessarily optimal)
	ParityCheck XBoundaryQubits, ZBoundaryQubits; // Data qubits on the X boundary and Z boundary


	/*
	 * Probably not needed...
	 */
	size_type GetBoundaryQubitFromStabilizer(ParityCheck& boundary_qubits, ParityCheck& stabilizer){
		for(auto i : stabilizer){
			if(boundary_qubits.find(i) != boundary_qubits.end())
				return i;
		}
	std::cerr << "Could not find Boundary Qubit" << std::endl;
	throw 0;
	}
};


#endif