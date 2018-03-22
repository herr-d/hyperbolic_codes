#include <autotune_wrapper.hpp>
#include <debug.hpp>
#include <getopt.h>
#include <iomanip>
#include <fstream>
#include <cstring>

Code_info code;

size_type get_data_qubit_id(int x, int y){
	//calculate id
	int id = y/2 * (2*code.r-1);
	id += y%2 * code.r;
	id += x/2;
	return id;
}



void generate_surface_code(){
	code.num_qubits = 2*code.r*code.r - 2*code.r + 1;

	ParityCheck Xbound_one;
	ParityCheck Xbound_two;
	ParityCheck Zbound_one;
	ParityCheck Zbound_two;


	ParityCheck Xdata_one;
	ParityCheck Xdata_two;
	ParityCheck Zdata_one;
	ParityCheck Zdata_two;

	//add stabilizers
	for(int x = 0; x < 2*code.r -1; ++x){
		for(int y = 0; y < 2*code.r -1; ++y){
			if((x%2 == 0 && y%2 == 0) || x%2 == 1 && y%2 ==1){ //data qubit
				if(x == 0){
					code.XBoundaryQubits.insert(get_data_qubit_id(x,y));
				}
				if(x == 2*code.r-2){
					code.XBoundaryQubits.insert(get_data_qubit_id(x,y));
				}
				if(y == 0){
					code.ZBoundaryQubits.insert(get_data_qubit_id(x,y));
				}
				if(y == 2*code.r-2){
					code.ZBoundaryQubits.insert(get_data_qubit_id(x,y));
				}
			}
			else if(x%2 == 1 && y%2 == 0){ // X stabilizer
				ParityCheck tmp;
	
				for(int dx=-1; dx <=1; dx=dx+2)
					if(x+dx >= 0 && x+dx < 2*code.r - 1)
						tmp.insert(get_data_qubit_id(x+dx,y));

					for(int dy=-1; dy <=1; dy=dy+2)
						if(y+dy >= 0 && y+dy < 2*code.r - 1)
							tmp.insert(get_data_qubit_id(x,y+dy));

				code.X_stabilizer.push_back(tmp);
				if(y == 0){
					Xbound_one.insert(code.X_stabilizer.size()-1);
				}
				if(y == 2*code.r-2){
					Xbound_two.insert(code.X_stabilizer.size()-1);
				}
			}

			else{ // Z stabilizer
				ParityCheck tmp;
				for(int dx=-1; dx <=1; dx=dx+2)
					if(x+dx >= 0 && x+dx < 2*code.r - 1)
							tmp.insert(get_data_qubit_id(x+dx,y));

				for(int dy=-1; dy <=1; dy=dy+2)
					if(y+dy >= 0 && y+dy < 2*code.r - 1)
							tmp.insert(get_data_qubit_id(x,y+dy));

				code.Z_stabilizer.push_back(tmp);
				if(x == 0){
					Zbound_one.insert(code.Z_stabilizer.size()-1);
				}
				if(x == 2*code.r-2){
					Zbound_two.insert(code.Z_stabilizer.size()-1);
				}
			}
		}
	}

	code.X_boundary.push_back(Xbound_one);
	code.X_boundary.push_back(Xbound_two);
	code.Z_boundary.push_back(Zbound_one);
	code.Z_boundary.push_back(Zbound_two);

	//add logical operators
	ParityCheck X_logical;
	ParityCheck Z_logical;
	for(int i = 0; i < code.r; ++i){
		X_logical.insert(get_data_qubit_id(2*i,code.r-1));
		Z_logical.insert(get_data_qubit_id(code.r-1,2*i));
	}
	code.X_operator.push_back(X_logical);
	code.Z_operator.push_back(Z_logical);

}





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


void process_arguments(int argc, char **argv){

	//first set default arguments in the code
	code.boot = false;
	code.r = 3;
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
				r = atoi(optarg);
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

	generate_surface_code();
	// initialize error simulation
	Wrapper simulation(code.probability);
	simulation.generate_recipe();

	if(code.boot)
		simulation.calculate_t_check();

	simulation.run_simulation(code.big_t_max);
	return 0;
}