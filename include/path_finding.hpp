#ifndef PATH_FINDING
#define PATH_FINDING

#include <basic_includes.hpp>
#include <unordered_set>

void Dijkstra(StabilizerContainer stabilizer, ParityCheck& result_container, size_type num_qubits);
void DataQubits(StabilizerContainer& stabilizers, size_type begin, ParityCheck& end, ParityCheck& result_container);

#endif