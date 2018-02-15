#ifndef BASIC_INCLUDES_HPP
#define BASIC_INCLUDES_HPP

#include <cassert>
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <set>
#include <vector>
#include <unordered_map>

using size_type = uint32_t;
using ParityCheck = std::set<size_type>;
using Graph = std::unordered_map<size_type,ParityCheck>;
using StabilizerContainer = std::vector<ParityCheck>;


#endif