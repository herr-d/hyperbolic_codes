#include <basic_includes.hpp>
#include <generate_code.hpp>
#include <path_finding.hpp>




size_type CSSCodes::get_number_qubits(){
	return num_qubits_;
}
size_type CSSCodes::get_number_stabilizers(){
	return Z_stabilizers_.size() + graph_.size();
}

size_type CSSCodes::get_number_X_stabilizers(){
	//return X_stabilizers_.size();
	return graph_.size();
}
size_type CSSCodes::get_number_Z_stabilizers(){
	return Z_stabilizers_.size();
}

size_type get_logical_start(ParityCheck& boundary_data_qubits, StabilizerContainer& stabilizer, ParityCheck& boundary){
	//find the start inside the boundary
	for(auto stabil_pos : boundary){
		for(auto qubit : stabilizer.at(stabil_pos)){
			if (boundary_data_qubits.find(qubit) != boundary_data_qubits.end())
			{
				return qubit;
			}
		}
	}
	throw 1;
}


void CSSCodes::generate_state(Code_info & code){
	//create qubits and a temporary lookup table
	std::map<std::pair<size_type,size_type>, size_type> tmp_lookup;
	ParityCheck XDataQubits;
	ParityCheck ZDataQubits;

	size_type qubit_id = 0;
	for(auto node : graph_){
		size_type center_node = node.first;
		for(size_type neighbor : node.second){
			if(center_node < neighbor){
				//add new qubit to lookup table
				tmp_lookup[std::make_pair(center_node, neighbor)] = qubit_id;
				++qubit_id;
			}

		}
	}
	code.num_qubits = tmp_lookup.size();

	//generate X checks with physical qubits
	for(auto X_check : graph_){
		if (X_check.second.size() > 1){
			ParityCheck tmp;
			for(auto neighbor : X_check.second){
				tmp.insert(tmp_lookup[std::make_pair(std::min(X_check.first,neighbor),std::max(X_check.first,neighbor))]);
			}
			code.X_stabilizer.push_back(tmp);
		}
		else if(X_check.second.size() == 1){
			XDataQubits.insert(tmp_lookup[std::make_pair(std::min(X_check.first,*X_check.second.begin()),std::max(X_check.first,*X_check.second.begin()))]);
		}
	}


	//generate Z-checks with physical qubits
	for(auto Z_check : Z_stabilizers_){
		ParityCheck tmp;
		for(size_type node : Z_check){
			for(auto neighbor : graph_.at(node)){
				if(neighbor > node)
					continue;
				if (Z_check.find(neighbor) != Z_check.end()){
					size_type qubit = tmp_lookup[std::make_pair(std::min(node,neighbor),std::max(node,neighbor))];
					tmp.insert(qubit);
					if(ZDataQubits.find(qubit) != ZDataQubits.end()){
						ZDataQubits.erase(qubit);
					}
					else{
						ZDataQubits.insert(qubit);
					}
				}
			}
		}
		code.Z_stabilizer.push_back(tmp);
	}

	std::swap(code.X_boundary, X_bound_result);
	std::swap(code.Z_boundary, Z_bound_result);
	std::swap(code.XBoundaryQubits, XDataQubits);
	std::swap(code.ZBoundaryQubits, ZDataQubits);

	//generate logical operator chains

	//find endpoints for X operator chain

	for(int i = 0; i < code.X_boundary.size()/2; ++i){
		size_type start_boundary = (2*i)%code.X_boundary.size();
		size_type end_boundary = (2*i + 1)%code.X_boundary.size();
		ParityCheck tmp;
		tmp.insert(get_logical_start(code.XBoundaryQubits, code.X_stabilizer, code.X_boundary[start_boundary]));
		tmp.insert(get_logical_start(code.XBoundaryQubits, code.X_stabilizer, code.X_boundary[end_boundary]));
		Dijkstra(code.X_stabilizer, tmp, code.num_qubits);
		code.X_operator.push_back(tmp);
	}

	//find endpoints for Z operator chain
	for(int i = 0; i < code.Z_boundary.size()/2; ++i){
		size_type start_boundary = (2*i)%code.Z_boundary.size();
		size_type end_boundary = (2*i + 1)%code.Z_boundary.size();
		ParityCheck tmp;
		tmp.insert(get_logical_start(code.ZBoundaryQubits, code.Z_stabilizer, code.Z_boundary[start_boundary] ));
		tmp.insert(get_logical_start(code.ZBoundaryQubits, code.Z_stabilizer, code.Z_boundary[end_boundary] ));
		Dijkstra(code.Z_stabilizer, tmp, code.num_qubits);
		code.Z_operator.push_back(tmp);
	}

	//resolve dual
	if(dual_){
		std::swap(code.X_stabilizer,code.Z_stabilizer);
		std::swap(code.X_boundary, code.Z_boundary);
		std::swap(code.X_operator, code.Z_operator);
		std::swap(code.XBoundaryQubits, code.ZBoundaryQubits);
	}

	return;
}



void CSSCodes::printout_graph(){
	std::cout << "----------- GRAPH -----------" << std::endl;
	for(auto element :graph_){
		std::cout << " "<<  element.first << " : ";
		for (auto k : element.second){
			std::cout << k <<" ";
		}
		std::cout << std::endl;
	}
}

size_type CSSCodes::add_vertex(){
	size_type id = graph_.size();
	graph_[id] = ParityCheck();
	return id;
}

void CSSCodes::merge_vertices(const size_type& qubit1, const size_type& qubit2){
	//merge all neighbors to smaller node
	size_type min = std::min(qubit1,qubit2);
	size_type max = std::max(qubit1,qubit2);

	for(auto n : graph_.at(max)){
		graph_.at(min).insert(n);
		graph_.at(n).insert(min);
		//remove the edge
		remove_edge(max, n);
	}
}

void CSSCodes::add_edge(const size_type& qubit1, const size_type& qubit2){
	graph_[qubit1].insert(qubit2);
	graph_[qubit2].insert(qubit1);
	++num_qubits_;
}

void CSSCodes::remove_edge(const size_type& qubit1, const size_type& qubit2){
	graph_[qubit1].erase(qubit2);
	graph_[qubit2].erase(qubit1);
}

void CSSCodes::add_face(ParityCheck qubits){
	Z_stabilizers_.push_back(qubits);
}

//
// Interesting part starts HERE
// Generation of the lattice
//

void Hyperbolic::build_lattice(){
	generate_center_face();
	for(int i = 0; i < size_; ++i){
		new_layer();
	}
	return;
}

void Hyperbolic::generate_center_face(){
	//this will create a single center face from which the patch grows
	assert(graph_.size() == 0);
	ParityCheck face_check;
	size_type id = add_vertex();
	face_check.insert(id);
	for(int i=0; i < r_-1; ++i){
		id = add_vertex();
		add_edge(id,id-1);
		face_check.insert(id);
	}
	add_edge(0,id);
	add_face(face_check);
}

void Hyperbolic::generate_center_vertex(){
	//This will create the patches around a center vertex.
	//This function will result in the dual lattice of generate_center_faces.
	assert(graph_.size() == 0);
	size_type id = add_vertex();
	Graph frontier;
	frontier[id] = ParityCheck();

	add_vertex_neighbors(frontier,id);
	gen_faces_from_vertex(frontier,id);
	// add last face
	ParityCheck face_check;
	face_check.insert(id);
	size_type last = *(std::prev(frontier[id].end(),1));
	face_check.insert(last);
	for (int i = 3; i < r_; ++i){
		id = add_vertex();
		add_edge(last,id);
		face_check.insert(id);
		last = id;
	}
	add_edge(id, *frontier[0].begin());
	face_check.insert(*frontier[0].begin());
	add_face(face_check);
}

void Hyperbolic::add_vertex_neighbors(Graph& frontier, size_type vertex){
	for (int i = graph_[vertex].size(); i < s_; ++i){

		size_type id = add_vertex();
		frontier[vertex].insert(id);
		add_edge(vertex, id);
	}
	return;
}

void Hyperbolic::get_frontier(Graph & frontier){
	for(auto element : graph_){
		//only elments on the boundary
		if (element.second.size() >= s_)
			continue;
		if(frontier.find(element.first)==frontier.end()){
			auto tmp = element.first;
			frontier[tmp] = ParityCheck();
		}
	}
	return;
}

void Hyperbolic::generate_frontier_loop(std::vector<size_type> & original_loop){
	Graph frontier;
	get_frontier(frontier);

	std::vector<size_type> loop(frontier.size());
	loop.at(0) = frontier.begin()->first;

	for(int i=1; i < frontier.size(); ++i){
		//find the neighbors
		std::vector<size_type> tmp;
		for(auto neigh : graph_[loop.at(i-1)])
			if(frontier.find(neigh) != frontier.end())
				tmp.push_back(neigh);

		assert(tmp.size() == 2);


		if(tmp.at(0) != loop.at(0) && tmp.at(0) != loop.at(abs(i-2))){
				loop.at(i) = tmp.at(0);
			}
		else{
			assert(tmp.at(0) == loop.at(0) || tmp.at(1) !=  loop.at(abs(i-2)));
			loop.at(i) = tmp.at(1);
		}
	}
	std::swap(loop, original_loop);
	return;
}



void Hyperbolic::new_layer(){
	//adds another layer around the code that is existing already
	assert(graph_.size()!=0);
	Graph frontier;
	get_frontier(frontier);
	
	//for this lets first create a path along the frontier
	//TODO ugly... fix this
	std::vector<size_type> loop;
	generate_frontier_loop(loop);

	//add new vertices around all frontier vertices
	for(auto element : frontier){
		add_vertex_neighbors(frontier, element.first);
	}
	//now close all unfinished faces
	for(auto f : frontier){
		gen_faces_from_vertex(frontier, f.first);
	}

	for(int i = 0; i < loop.size(); ++i){
		gen_faces_from_edge(frontier, loop.at(i), loop.at((i+1)%loop.size()));
	}
}


void Hyperbolic::gen_faces_from_vertex(Graph& frontier, size_type vertex){
	//assert(frontier[vertex].size() > 1);
	ParityCheck::iterator it = frontier[vertex].begin();
	size_type last = *it;
	it = std::next(it, 1);
	while(it != frontier[vertex].end()){
		//now add as many vertices to get an r_ sided face
		ParityCheck face_check;

		//first add the vertex on the frontier
		face_check.insert(vertex);

		//now add the first already generated vertex
		face_check.insert(last);
		//now add any additional vertices in between
		for (int i = 3; i < r_; ++i){
			size_type id = add_vertex();
			face_check.insert(id);
			add_edge(last, id);
			last = id;
		}
		face_check.insert(*it);
		add_edge(last, *it);
		add_face(face_check);
		last = *it;
		it = std::next(it, 1);
	}
}


void Hyperbolic::gen_faces_from_edge(Graph& frontier, size_type vertex1, size_type vertex2){
	//generate remaining edge for face between vertex1 and vertex2 and their children

	ParityCheck face_check({vertex1,vertex2});
	//chose the one 
	auto last_it = std::prev(frontier[vertex1].end(),1);
	if (graph_[*last_it].size() > graph_[*frontier[vertex1].begin()].size()){
		last_it =  frontier[vertex1].begin();
	}
	auto end_it = frontier[vertex2].begin();
	if (graph_[*end_it].size() > graph_[*std::prev(frontier[vertex2].end(),1)].size()){
		end_it =  std::prev(frontier[vertex2].end(),1);
	}

	//now add as many vertices to get an r_ sided face
	auto last = *last_it;
	face_check.insert(last);
	for (int i = 4; i < r_; ++i){
		size_type id = add_vertex();
		add_edge(last, id);
		face_check.insert(id);
		last = id;
	}
	add_edge(last, *end_it);
	face_check.insert(*end_it);
	add_face(face_check);
	return;
}


void Hyperbolic::generate_rough_edges(short number){
	StabilizerContainer X_boundary;
	StabilizerContainer Z_boundary;

	size_type shift = 0; // rotate the rough/smooth edges around the boundary

	number = 2* number;
	// add loop
	std::vector<size_type> loop;
	generate_frontier_loop(loop);

	printout_graph();	



	if(loop.size()%number != 0){
		std::cerr << "Cannot partition into equal sizes for rough edges. Last partition will be smaller." << std::endl;
	}

	// iterate through loop and divide into pieces (interested in qubits = edges not vertices)
	for(int i = 0; i < number; ++i){
		if (i%2){
			X_boundary.push_back(ParityCheck());
			for(int j=0; j < loop.size()/number; ++j){
				//two vertices
				size_type v1 = ((i* loop.size()/number + j) + shift)%loop.size();
				size_type v2 = (v1 + 1)%loop.size();
				// remove qubits not needed
				--num_qubits_;
				remove_edge(loop.at(v2),loop.at(v1));
			}
			/*// remove nodes with no neighbors. Nodes with one neighbor will be ignored when exporting
			// but are still needed to determine the qubit ids
			for(int j=0; j < loop.size()/number; ++j){
				size_type v1 = (i* loop.size()/number + j + shift)%loop.size();
				if(graph_.at(loop.at(v1)).size() == 0 ){
					graph_.erase(loop.at(v1));
				}
			}*/
			// add neighbors of vertices with single neighbors to the boundary stabilizers
			for(int j=0; j <= loop.size()/number; ++j){
				size_type v1 = ((i* loop.size()/number + j) + shift) % loop.size();
				if(graph_[loop.at(v1)].size() == 1 ){
					X_boundary.back().insert(*graph_[loop.at(v1)].begin());
				}
			}
		}
		else{
			Z_boundary.push_back(ParityCheck());
			for(int j=0; j <= loop.size()/number; ++j){
				size_type v1 = (i* loop.size()/number + j + shift) % loop.size();
				size_type v2 = (v1 + 1)%loop.size();
				// insert boundary Z stabilizer
				for(uint i = 0; i<Z_stabilizers_.size(); ++i){
					if(Z_stabilizers_.at(i).find(loop.at(v1)) != Z_stabilizers_.at(i).end() && Z_stabilizers_.at(i).find(loop.at(v2)) != Z_stabilizers_.at(i).end()){
						Z_boundary.back().insert(i);
					}
				}
			}
		}
	}

	
	// Some of the X-stabilizers get deleted and thus the boundary ids need to change
	// Z stabilizers do not get deleted and thus do not need to be modified
	// boundaries
	std::set<size_type> deleted_vector;
	
	for(auto i : graph_){
		if(i.second.size() == 0){
			deleted_vector.insert(i.first);
		}
	}


	//TODO

	for(const auto& bound : Z_boundary){
		std::vector<size_type> tmp;
		for(const auto& stabil : bound){
			for(auto node : stabil){
				td::cout << "bound "<< eighbor << std::endl;
				for(auto stabil :Z_stabilizer)
				if(Z_stabilizers.
					graph_.at(neighbor).size() == 1){
					tmp.push_back(neighbor);
				}
			}
		}
		std::cout << std::endl;
		if (tmp.size()>1){
			assert(tmp.size() == 2);
			merge_vertices(tmp.at(0),tmp.at(1));
			deleted_vector.insert(tmp.at(1));
		}
	}


	exit(0);
	//only rarely executed... but still expensive
	for(auto & stabilizer : Z_stabilizers_){
		for(auto deleted : deleted_vector)
			stabilizer.erase(deleted);
	}



	for(auto i : graph_){
		if(i.second.size() == 2 || i.second.size()==1){
			deleted_vector.insert(i.first);
		}
	}

	// since not all vertices are used in the output, the ids need to be adjusted
	// probably not needed but better safe than sorry


	for(auto & edge : X_boundary){
		ParityCheck new_edge;
		for(auto id : edge){
			size_type count = 0;
			for(auto deleted : deleted_vector){
				if(id >= deleted){
					++count;
				}
			}
			new_edge.insert(id-count);
		}
		std::swap(new_edge, edge);
	}

	// remove border artifacts, where one stabilizer has more than one final qubit
	// these lead to errors in the code



	std::swap(X_bound_result, X_boundary);
	std::swap(Z_bound_result, Z_boundary);
	return;
}
