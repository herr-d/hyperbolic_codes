#include <path_finding.hpp>



void Dijkstra(StabilizerContainer stabilizer, ParityCheck& result_container, size_type num_qubits){
	std::vector<size_type> dist(num_qubits+1, num_qubits);
	std::vector<size_type> prev(num_qubits+1, num_qubits);
	assert(result_container.size() == 2);

	size_type begin, end;


	//very ugly !!!!!
	bool second = false;
	for(auto element : result_container){
		if (second){
			end = element;
		}
		else{
			begin = element;
		}
		second = true;
	}

	std::unordered_set<size_type> open_set;

	for(int i=0; i < num_qubits; ++i){
		open_set.insert(i);
	}

	dist.at(begin) = 0;

	while(open_set.size() >0){
		size_type min_element = num_qubits;
		size_type cur_dist = num_qubits;
		for(auto element : open_set){
			if (dist.at(element) < cur_dist){
				cur_dist = dist.at(element);
				min_element = element;
			}
		}

		open_set.erase(min_element);

		if(min_element == end)
			break;

		for(auto stabil : stabilizer){
			if(stabil.find(min_element) != stabil.end()){
				for(size_type neighbor : stabil){
					if (dist.at(neighbor) > cur_dist+1){
						dist.at(neighbor) = cur_dist + 1;
						prev.at(neighbor) = min_element;
					}
				}
			}
		}

	}

	//generate path by tracking backwards
	assert(dist.at(end) < num_qubits); // checks if the end was found in the graph

	
	size_type cur_element = end;
	ParityCheck path_qubits;
	path_qubits.insert(cur_element);

	while(cur_element != begin){
		cur_element = prev[cur_element];
		path_qubits.insert(cur_element);
	}
	std::swap(path_qubits, result_container);
}


// Path finding between two stabilizers. Data qubits along the path are returned.
void DataQubits(StabilizerContainer& stabilizers, size_type begin, ParityCheck& end, ParityCheck& result_container){
	std::vector<size_type> dist(stabilizers.size()+1, stabilizers.size());
	std::vector<size_type> prev(stabilizers.size()+1, stabilizers.size());

	std::unordered_set<size_type> open_set;

	for(int i=0; i < stabilizers.size(); ++i){
		open_set.insert(i);
	}

	dist.at(begin) = 0;

	//TODO: TOO SLOW!!!!
	while(open_set.size() > 0){
		size_type min_element = stabilizers.size();
		size_type cur_dist = stabilizers.size();
		for(auto element : open_set){
			if (dist.at(element) < cur_dist){
				cur_dist = dist.at(element);
				min_element = element;
			}
		}

		open_set.erase(min_element);

		if(end.find(min_element) != end.end())
			break;

		//how to get neighboring stabilizers efficiently?
		for(auto i : open_set){
			for(auto s : stabilizers.at(min_element)){
				if(stabilizers.at(i).find(s) != stabilizers.at(i).end()){
					dist.at(s) = cur_dist + 1;
					prev.at(s) = min_element;
				}
			}

		}
	}

	size_type cur_element;
	for(auto i : end){
		cur_element = i;
		if (dist.at(i) < stabilizers.size()+1)
			break;
	}

	end.clear();
	end.insert(cur_element);

	ParityCheck path_qubits;

	while(cur_element != begin){
		for(auto i : stabilizers.at(cur_element)){
			if(stabilizers.at(prev[cur_element]).find(i) != stabilizers.at(prev[cur_element]).end())
				path_qubits.insert(i);
		}
		cur_element = prev[cur_element];
	}
	std::swap(path_qubits, result_container);
}