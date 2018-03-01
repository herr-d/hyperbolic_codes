#include <generate_code.hpp>
#include <autotune_wrapper.hpp>


Code_info code;


void generate_hyperbolic_code(size_type r, size_type s, size_type size){
		Hyperbolic hyperbolic(r, s, size);
		hyperbolic.build_lattice();
		
		#ifdef DEBUG
		std::cout << "Created code with " << hyperbolic.get_number_qubits() << " data qubits" << std::endl;
		#endif

		//hyperbolic.generate_rough_edges(2);

		hyperbolic.generate_state(code);
}



int main(int argc, char **argv){

	double error_probability = 0.0001;

	generate_hyperbolic_code(5,4,5);

	#ifdef DEBUG
	std::cout << "X stabilizers: " << code.X_stabilizer.size() << std::endl;
	std::cout << "Z stabilizers: " << code.Z_stabilizer.size() << std::endl;
	#endif
	
	// initialize error simulation
	Wrapper simulation(code, error_probability);

	// run individual experiments


	// initialize RNG & error class
	return 0;
}