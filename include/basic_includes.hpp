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
	StabilizerContainer X_stabilizer, Z_stabilizer;
	StabilizerContainer X_boundary, Z_boundary;
	size_type num_qubits;
	StabilizerContainer X_operator, Z_operator;
};


#endif