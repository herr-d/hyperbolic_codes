#ifndef GENERATE_CODE_HPP
#define GENERATE_CODE_HPP


#include <basic_includes.hpp>


class CSSCodes
{
protected:
	StabilizerContainer Z_stabilizers_;
	StabilizerContainer X_stabilizers_;
	StabilizerContainer Y_stabilizers_;
	size_type num_qubits_;
	size_type num_vertices_;
	Graph graph_;
	
public:
	CSSCodes() : num_vertices_(0), num_qubits_(0), graph_(), Z_stabilizers_(),
				X_stabilizers_(), Y_stabilizers_() {};
	~CSSCodes(){};

	void generate_patch();
	void generate_stabilizers();

	size_type get_number_qubits();
	size_type get_number_stabilizers();
	size_type get_number_X_stabilizers();
	size_type get_number_Y_stabilizers();
	size_type get_number_Z_stabilizers();

	const StabilizerContainer& get_X_stabilizer_list();
	const StabilizerContainer& get_Y_stabilizer_list();
	const StabilizerContainer& get_Z_stabilizer_list();

	void printout_graph();
	void add_vertex();
	void add_edge(const size_type& qubit1, const size_type& qubit2);
	void add_face(ParityCheck qubits);
};


class Hyperbolic : public CSSCodes{
private:
	size_type r_; //number of face-neighbors
	size_type s_; //number of vertex-neighbors
	size_type size_; //how many layers will be created
	bool dual_;

	void add_vertex_neighbors(Graph& frontier, size_type vertex);
	void gen_faces_from_edge(Graph &frontier, size_type vertex1, size_type vertex2);
	void gen_faces_from_vertex(Graph& frontier, size_type vertex);

public:
	Hyperbolic(size_type r, size_type s, size_type size) : CSSCodes() {
		// the construction relies on vertex neighbors > 3
		// if ==3 then create the dual lattice and take the dual of the result
		assert(s>3);
		assert(r>3);
		r_ = r;
		s_ = s;
		size_=size;
	}

	~Hyperbolic(){};
	void build_lattice();
	void generate_rough_edges();
	void new_layer();
	void generate_center_vertex();
	void generate_center_face();
	void add_rough_edges();
};


#endif