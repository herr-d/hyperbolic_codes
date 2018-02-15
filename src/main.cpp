#include <generate_code.hpp>
#include <simulate_errors.hpp>


int main(int argc, char **argv){
	Hyperbolic small(5, 4, 10);
	small.build_lattice();
	small.generate_stabilizers();
	std::cout << small.get_number_Z_stabilizers() << std::endl;
	return 0;
}