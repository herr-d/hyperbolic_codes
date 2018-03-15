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



struct Code_info{
	StabilizerContainer X_stabilizer, Z_stabilizer; // list of stabilizers
	StabilizerContainer X_boundary, Z_boundary; // stabilizer boundary
	size_type num_qubits; // total number of qubits
	StabilizerContainer X_operator, Z_operator; //operator chains for logical qubits (not necessarily optimal)
	ParityCheck XBoundaryQubits, ZBoundaryQubits; // Data qubits on the X boundary and Z boundary

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