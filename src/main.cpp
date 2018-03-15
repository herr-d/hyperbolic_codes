#include <generate_code.hpp>
#include <autotune_wrapper.hpp>
#include <debug.hpp>

Code_info code;


void write_usage(char *arg0){
	std::cout << "Usage: " << arg0 << "num_face_neighbors " << "num_vertex_neighbors "  << "layers" << std::endl;
	std::cout << " num_face_neighbors: number of neighboring faces for the lattice\n";
	std::cout << " num_vertex_neighbors: number of neighboring vertices for the lattice\n";
	std::cout << " layers: determines the size of the lattice\n\n";
	std::cout << "Example: Surface code with distance 5/6 (currently: X error distance + 1 = Z error distance)\n";
	std::cout << arg0 << " 4 4 2" << std::endl;
}


void generate_hyperbolic_code(size_type r, size_type s, size_type size){
	Hyperbolic hyperbolic(r, s, size);
	hyperbolic.build_lattice();

	#ifdef DEBUG
		std::cout << "Created code with " << hyperbolic.get_number_qubits() << " data qubits" << std::endl;
	#endif

	hyperbolic.generate_rough_edges(2);

	#ifdef DEBUG
		std::cout << "Generated rough edges "<< std::endl;
	#endif

	hyperbolic.generate_state(code);
	#ifdef DEBUG
		std::cout << "Finalized the creation of the lattice -> variable \"code\" is now valid" << std::endl;
		//printout_code();
	#endif
}




int main(int argc, char **argv){


	double error_probability = 0.0001;

	generate_hyperbolic_code(4,4,5);

	#ifdef DEBUG
	#endif
	
	// initialize error simulation
	Wrapper simulation(error_probability);
	simulation.generate_recipe();
	//simulation.run_simulation(4);
	simulation.calculate_t_check(error_probability);

	// run individual experiments


	return 0;
}