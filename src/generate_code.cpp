#include <basic_includes.hpp>
#include <generate_code.hpp>


size_type CSSCodes::get_number_qubits(){
	return num_qubits_;
}
size_type CSSCodes::get_number_stabilizers(){
	return Z_stabilizers_.size() + X_stabilizers_.size() + Y_stabilizers_.size();
}

size_type CSSCodes::get_number_X_stabilizers(){
	//return X_stabilizers_.size();
	return graph_.size();
}
size_type CSSCodes::get_number_Z_stabilizers(){
	return Z_stabilizers_.size();
}
size_type CSSCodes::get_number_Y_stabilizers(){
	return Y_stabilizers_.size();
}



const StabilizerContainer& CSSCodes::get_X_stabilizer_list(){
	return X_stabilizers_;
}
const StabilizerContainer& CSSCodes::get_Z_stabilizer_list(){
	return Z_stabilizers_;
}
const StabilizerContainer& CSSCodes::get_Y_stabilizer_list(){
	return Y_stabilizers_;
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


void CSSCodes::generate_stabilizers(){
	/*//X_stabilizer
	StabilizerContainer tmp;
	for (auto element : graph_){
		tmp.push_back(element.second);
	}
	std::swap(tmp,X_stabilizers_);*/
}

void CSSCodes::add_vertex(){
	graph_[num_vertices_] = ParityCheck();
	++num_vertices_;
}

void CSSCodes::add_edge(const size_type& qubit1, const size_type& qubit2){
	graph_[qubit1].insert(qubit2);
	graph_[qubit2].insert(qubit1);
	++num_qubits_;
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
	assert(num_vertices_ == 0);
	ParityCheck face_check;
	add_vertex();
	face_check.insert(num_vertices_-1);
	for(int i=0; i < r_-1; ++i){
		add_vertex();
		add_edge(num_vertices_-1,num_vertices_-2);
		face_check.insert(num_vertices_-1);
	}
	add_edge(0,num_vertices_-1);
	add_face(face_check);
}

void Hyperbolic::generate_center_vertex(){
	//This will create the patches around a center vertex.
	//This function will result in the dual lattice of generate_center_faces.
	assert(num_vertices_ == 0);
	add_vertex();
	Graph frontier;
	frontier[0] = ParityCheck();

	add_vertex_neighbors(frontier,0);
	gen_faces_from_vertex(frontier,0);
	// add last face
	ParityCheck face_check;
	face_check.insert(0);
	size_type last = *(std::prev(frontier[0].end(),1));
	face_check.insert(last);
	for (int i = 3; i < r_; ++i){
		add_vertex();
		add_edge(last, num_vertices_-1);
		last = num_vertices_-1;
		face_check.insert(last);
	}
	add_edge(last, *frontier[0].begin());
	face_check.insert(*frontier[0].begin());
	add_face(face_check);
}

void Hyperbolic::add_vertex_neighbors(Graph& frontier, size_type vertex){
	for (int i = graph_[vertex].size(); i < s_; ++i){

		add_vertex();
		frontier[vertex].insert(num_vertices_-1);
		add_edge(vertex, num_vertices_-1);
	}
	return;
}

void Hyperbolic::new_layer(){
	//adds another layer around the code that is existing already
	assert(num_vertices_!=0);
	Graph frontier;
	for(auto element : graph_){
		//only elments on the boundary
		if (element.second.size() >= s_)
			continue;
		if(frontier.find(element.first)==frontier.end()){
			auto tmp = element.first;
			frontier[tmp] = ParityCheck();
		}
	}
	//add new vertices around all frontier vertices
	for(auto element : frontier){
		add_vertex_neighbors(frontier, element.first);
	}
	//now close all unfinished faces
	for(auto f : frontier){
		gen_faces_from_vertex(frontier, f.first);
	}

	//for this lets first create a path along the frontier
	std::vector<size_type> loop;
	loop.push_back(frontier.begin()->first);
	for(int i=loop.size(); i < frontier.size(); ++i){
		//find the neighbors
		std::vector<size_type> tmp;
		for(auto neigh : graph_[loop.back()])
			if(frontier.find(neigh) != frontier.end())
				tmp.push_back(neigh);

		assert(tmp.size() == 2);


		if(tmp.at(0) != *loop.begin() && tmp.at(0) != loop.at(abs((loop.size()-2)%loop.size()))){
				loop.push_back(tmp.at(0));
			}
		else{
			assert(tmp.at(0) == *loop.begin() || tmp.at(1) !=  loop.at(abs((loop.size()-2)%loop.size())));
			loop.push_back(tmp.at(1));
		}
	}
	for(int i = 0; i < loop.size(); ++i){
		gen_faces_from_edge(frontier, loop.at(i), loop.at((i+1)%loop.size()));
	}
}

void Hyperbolic::generate_rough_edges(){

}

void Hyperbolic::gen_faces_from_vertex(Graph& frontier, size_type vertex){
	//assert(frontier[vertex].size() > 1);
	ParityCheck::iterator it = frontier[vertex].begin();
	size_type last = *it;
	it = std::next(it, 1);
	while(it != frontier[vertex].end()){
		//now add as many vertices to get an r_ sided face
		ParityCheck face_check;
		face_check.insert(vertex);
		for (int i = 3; i < r_; ++i){
			add_vertex();
			face_check.insert(num_vertices_-1);
			add_edge(last, num_vertices_-1);
			last = num_vertices_-1;
		}
		add_face(face_check);
		add_edge(last, *it);
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
		add_vertex();
		add_edge(last, num_vertices_-1);
		face_check.insert(num_vertices_-1);
		last = num_vertices_-1;
	}
	add_edge(last, *end_it);
	face_check.insert(*end_it);
	add_face(face_check);
	return;
}