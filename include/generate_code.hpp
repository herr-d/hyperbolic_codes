#ifndef GENERATE_CODE_HPP
#define GENERATE_CODE_HPP

#include <basic_includes.hpp>
#include <autotune_wrapper.hpp>

/* Parent class for Hyperbolic codes
 * Any CSS code can be set up manually using
 * functions defined. 
 * MORE INTERESTING: Derived class
*/
class CSSCodes
{
protected:
	StabilizerContainer Z_stabilizers_; // container for the face parity checks
	StabilizerContainer X_bound_result; // container for the boundaries
	StabilizerContainer Z_bound_result; // container for the boundaries
	size_type num_qubits_; // total number of qubits
	Graph graph_; // stores the graph - container for node parity checks
	bool dual_; // in the final output swap X and Z-stabilizers
	
public:
	CSSCodes() : num_qubits_(0), dual_(false), Z_stabilizers_(), graph_(), X_bound_result(), Z_bound_result(){};
	~CSSCodes(){};


	/* returns the number of qubits the current layout has
	 * [in]: no input
	 * [out]: number of qubits
	*/
	size_type get_number_qubits();
	
	/* returns the total number of stabilizers
	 * [in]: no input
	 * [out]: count of all stabilizers
	*/
	size_type get_number_stabilizers();
	
	/* returns the number of X_stabilizers
	 * [in]: no input
	 * [out]: count of all X_stabilizers
	*/
	size_type get_number_X_stabilizers();
	
	/* returns the number of Z_stabilizers
	 * [in]: no input
	 * [out]: count of all Z_stabilizers
	*/
	size_type get_number_Z_stabilizers();

	/* Prints out the graph. This only makes sense for small codes/ for manual checking
	 * [in]: no input
	 * [out]: no output
	*/
	void printout_graph();

	/* adds a vertex to the graph and assings it the proper id
	 * [in]: no input
	 * [out]: Id of the vertex
	*/
	size_type add_vertex();

	/* adds an edges between the given two vertices
	 * [in]: id of vertex2
	 *		 id of vertex1
	 * [out]: no output
	*/
	void add_edge(const size_type& qubit1, const size_type& qubit2);

	/* adds a face to the graph
	 * [in]: a set of vertices that constitute a face
	 * [out]: no output
	*/
	void add_face(ParityCheck qubits);

	/* generates the physical qubits, and stabilizers and writes that information in the code reference
	 * [in]: container of the type code
	 * [out]: no output
	*/
	void generate_state(Code_info & code);
};


class Hyperbolic : public CSSCodes{
private:
	size_type r_; //number of face-neighbors
	size_type s_; //number of vertex-neighbors
	size_type size_; //how many layers will be created

	/* 
	 * [in]: 
	 * [out]:
	*/
	void add_vertex_neighbors(Graph& frontier, size_type vertex);
	
	/*
	 * [in]: 
	 * [out]:
	*/
	void gen_faces_from_edge(Graph &frontier, size_type vertex1, size_type vertex2);
	
	/* 
	 * [in]: 
	 * [out]:
	*/
	void gen_faces_from_vertex(Graph& frontier, size_type vertex);
	
	/* Infers which vertices are not fully connected yet and adds them to the frontiers reference
	 * [in]: reference to a graph object
	 * [out]: no return value
	*/
	void get_frontier(Graph & frontier);
	
	/* creates a closed loop from the frontier of the code
	 * [in]: reference to vector that will contain the information
	 * [out]: no return
	*/
	void generate_frontier_loop(std::vector<size_type> & original_loop);

public:
	/* Here all important variables are set.
	 * [in]: r: number of neighbors for each face
	 		 s: number of neighbors for each vertex
	 		 size: how large a lattice is generated (number of iterations where the lattice grows by one layer)
	 * [out]: no output
	*/
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


	/* This function completely builds the lattice given {r,s} and size.
	 * NOTE: to encode qubits, rough edges still need to be generated
	 * [in]: no input
	 * [out]: no output
	*/
	void build_lattice();

	/* generates rough edges
	 * [in]: the number of rough edges that should be created (ex. surface code: 2)
	 * [out]: no output
	*/
	void generate_rough_edges(short number);

	/*  grows a new layer around the existing lattice
	 * [in]: no input
	 * [out]: no output
	*/
	void new_layer();

	/* generates a single vertex surrounded by its adjacient faces: starting point for the creation of arbitrary sized patches
	 * [in]: no input
	 * [out]: no output
	*/
	void generate_center_vertex();

	/* generates a single face and its nodes: starting point for arbitrary sized patches
	 * [in]: no input
	 * [out]: no output
	*/
	void generate_center_face();
};

#endif
