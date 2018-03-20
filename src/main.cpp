#include <generate_code.hpp>
#include <autotune_wrapper.hpp>
#include <debug.hpp>
#include <getopt.h>
#include <iomanip>
#include <fstream>
#include <cstring>

Code_info code;


void write_usage(char *arg0){


	std::cout << "\nUsage:\n" << "\t" << arg0 << "\t[--OPTIONS]\n\n" << std::endl;
	std::cout << "Options:\n";
	std::cout << "\t--layers     \tnumber of layers and thus the size of the code\n";
	std::cout << "\t--k          \tnumber of neighbors each surface has\n";
	std::cout << "\t--r          \tnumber of neighbors each vertex has\n";
	std::cout << "\t--t_check    \tThe number of timesteps between each round of perfect stabilizer measurements\n";
	std::cout << "\t--big_t_max  \tThe maximum value of big_t\n";
	std::cout << "\t--probability\tThe probability of a random error\n";
	std::cout << "\t--boot       \tFlag if boot up phase should be performed to optimize t_check\n";
	std::cout << "\t--seed0      \tThe first random seed\n";
	std::cout << "\t--seed1      \tThe second random seed\n";
	std::cout << "\t--help       \tPrints this help message\n";

	std::cout << "\nStandard parameters:\n\t" <<  arg0 << " --layers " << code.l << " -k " << code.s
	<< " -r " << code.r << " --t_check " << code.t_check << " --big_t_max " << code.big_t_max
	<< " --probability " << code.probability << " --seed0 " << code.seed0 << " --seed1 " << code.seed1 << "\n" << std::endl;
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


void process_arguments(int argc, char **argv){

	//first set default arguments in the code
	code.boot = false;
	code.r = 4;
	code.s = 4;
	code.l = 4;
	code.t_check = 10;
	code.t_delete = 100;
	code.big_t_max = 20000;
	code.probability = 0.001;
	code.seed0 = 42;
	code.seed1 = 42;
	code.max_num_X = 5000;
	code.max_num_Z = 5000;
	code.dir = (char*)malloc(100*sizeof(char));
	strcpy(code.dir, "../../autotune/ex/ems/");

	bool boot = code.boot;
	uint16_t r = code.r;
	uint16_t s = code.s;
	uint16_t l = code.l;
	uint16_t t_check = code.t_check;
	size_type big_t_max = code.big_t_max;
	double probability = code.probability;
	size_type seed0 = code.seed0;
	size_type seed1 = code.seed1;
	size_type t_delete = code.t_delete;

	int c;

	while(1){
		struct option long_options[] = 
		{
			{"layers", required_argument, 0, 'l'},
			{"k", required_argument, 0, 'k'},
			{"r", required_argument, 0, 'r'},
			{"t_check", required_argument, 0, 't'},
			{"big_t_max", required_argument, 0, 'T'},
			{"probability", required_argument, 0, 'p'},
			{"seed0", required_argument, 0, 's'},
			{"seed1", required_argument, 0, 'S'},
			{"help", no_argument, 0, 'h'},
			{"boot", no_argument, 0, 'b'},
			{"t_delete", required_argument, 0, 'x'},
			{"ems", required_argument, 0, 'e'},
			{"d", required_argument, 0, 'd'}
		};
		int option_index = 0;
		c = getopt_long(argc, argv, "a:s:r:t:T:p:s:S:", long_options, &option_index);
		if(c == -1)
			break;

		switch(c)
		{
			case 'l':{
				l = atoi(optarg);
				break;
			}
			case 'k':{
				s = atoi(optarg);
				break;
			}
			case 'r':{
				r = atoi(optarg);
				break;
			}
			case 'b':{
				boot = true;
				break;
			}
			case 't':{
				t_check = atoi(optarg);
				break;
			}
			case 'T':{
				big_t_max = atoi(optarg);
				break;
			}
			case 'p':{
				probability = atof(optarg);
				break;
			}
			case 's':{
				seed0 = atoi(optarg);
				break;
			}
			case 'S':{
				seed1 = atoi(optarg);
				break;
			}
			case 'x':{
				t_delete = atoi(optarg);
				break;
			}
			case 'e':{
				strcpy(code.dir, optarg);
				break;
			}
			case 'd':{
				//do nothing: only here to be compatible with the same execution scripts
				break;
			}
			case 'h':{
				write_usage(argv[0]);
				exit(0);			}
			default:{
				write_usage(argv[0]);
				exit(1);
			}
		}

	}

	code.boot = boot;
	code.r = r;
	code.s = s;
	code.l = l;
	code.t_check = t_check;
	code.big_t_max = big_t_max;
	code.probability = probability;
	code.seed0 = seed0;
	code.seed1 = seed1;
	code.t_delete = t_delete;


	return;
}


int main(int argc, char **argv){

	process_arguments(argc, argv);

	generate_hyperbolic_code(code.r,code.s,code.l);

	// initialize error simulation
	Wrapper simulation(code.probability);
	simulation.generate_recipe();

	if(code.boot)
		simulation.calculate_t_check();

	simulation.run_simulation(code.big_t_max);
	return 0;
}