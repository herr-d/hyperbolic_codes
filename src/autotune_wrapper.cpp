#include <autotune_wrapper.hpp>
#include <string.h>
#include <path_finding.hpp>
#include <climits>


#define DE_HT_FACTOR 100
#define STICK_HT_FACTOR 200
#define QUBIT_HT_SIZE 1000
// This rearranges example 1 such that it can be used with hyperbolic codes


// DEBUG FUNCTIONS

//
// This function only works for codes generated by generate_surface_code()
//
void autotune_coords_to_real(int i, int j, int t){
	int dist = code.r;
	int x = 0;
	int y = 0;
	if(j<0){
		//boundaries
		if(i < 0){
			//temporal boundaries

		}
		else{
			//physical boundaries
			x = dist;
			if( i == 0)
				y =  -1;
			else
				y = 2*dist;
			if(j == -2)
				std::swap(x,y);
		}
	}
	else{
		if(i < code.num_qubits){
			//data qubits
			y = 2*(i / (2*dist-1));
			int rem = i % (2*dist - 1);
			y += rem/dist;
			x = (rem / dist) ? 1 + 2*(rem % dist): 2*(rem % dist);
		}
		else if(i < code.num_qubits + code.X_stabilizer.size()){
			std::cout << " X ";
			i -= code.num_qubits;
			x = (2 * (i/dist)) + 1;
			y = 2* (i % dist);
		}
		else{
			// Z stabilizer
			std::cout << " Z ";
			i -= code.num_qubits + code.X_stabilizer.size();
			x = 2*(i / (dist-1));
			y = 1 + 2*(i % (dist - 1));
		}
	}

	std::cout << " x: " << x << " y: " << y;
	return;
}



int infer_test_boundary(size_type pos, bool is_dual){
	//infer from the stabilizer id which boundary is needed
	if(is_dual){
		for(int i=0; i < code.Z_boundary.size(); ++i){
			if (code.Z_boundary.at(i).find(pos) != code.Z_boundary.at(i).end()){
				return i;
			}
		}
	} else{
		for(int i=0; i < code.X_boundary.size(); ++i){
			if (code.X_boundary.at(i).find(pos) != code.X_boundary.at(i).end()){
				return i;
			}
		}
	}
	return -1;
}


// HERE BEGINS THE REAL CODE


SET* Wrapper::create_boundary_set(QC *qc, int id, int type, int i, int j, int t) {
	SET *set;
	set = qc_create_set_adv(qc, type, i, j, t, NULL);
	boundaries[id] = set->ball;

	return set;
}

BALL *get_boundary(int i, int j, long int big_t, int type, void *boundaries) {
	BALL **bdys;


	bdys = (BALL **)boundaries;

	// Boundaries do not have a big_t set, they are initialised to LONG_MAX.
	if (big_t != LONG_MAX) {
		std::cout << "hmm" << std::endl;
		std::cout << "long max: " << LONG_MAX << std::endl;
		std::cout<< "get_boundaries: i: " << i <<" j: " << j  << " big_t: " << big_t << " type: " << type << std::endl;
		return NULL;
	}

	if(j == -1) {
		if(i == -1){
			assert(type == PRIMAL_BOUNDARY);
			return bdys[0];
		}else if(i == -2){
			assert(type == PRIMAL_BOUNDARY);
			return bdys[1];
		}
	} else if(j == -2){
		if(i == -1){
			assert(type == DUAL_BOUNDARY);
			return bdys[2];
		}else if(i == -2){
			assert(type == DUAL_BOUNDARY);
			return bdys[3];
		}
	} else if(j == -3){
		if(i == -1){
			assert(type == PRIMAL_BOUNDARY);
			return bdys[4];
		}else if(i == -2){
			assert(type == PRIMAL_BOUNDARY);
			return bdys[5];
		}
	} else if(j == -4){
		if(i == -1){
			assert(type == DUAL_BOUNDARY);
			return bdys[6];
		}else if(i == -2){
			assert(type == DUAL_BOUNDARY);
			return bdys[7];
		}
	}
	// There is no boundary at this coordinate
	return NULL;
}


SET* Wrapper::infer_boundary(size_type pos, bool is_dual){
	//infer from the stabilizer id which boundary is needed
	if(is_dual){
		for(int i=0; i < code.Z_boundary.size(); ++i){
			if (code.Z_boundary.at(i).find(pos) != code.Z_boundary.at(i).end()){
				return dual_boundary[i];
			}
		}
	} else{
		for(int i=0; i < code.X_boundary.size(); ++i){
			if (code.X_boundary.at(i).find(pos) != code.X_boundary.at(i).end()){
				return primal_boundary[i];
			}
		}
	}
	return NULL;
}

void Wrapper::create_initial_syndrome_and_set_arrays(){
	for(int i = 0; i < code.X_stabilizer.size(); ++i){
		primal_syndrome[i] = qc_create_syndrome();
		qc_insert_syndrome(dp_qc->qc, primal_syndrome[i]);

		//debug
		//std::cout << "create set with boundary "<< infer_test_boundary(i,false) << " and type: " << 1 << " at id "<< i << " and position " << i << std::endl;

		primal_set[i][1] = qc_create_set_adv(dp_qc->qc, PRIMAL, code.num_qubits + i, 0, 1, infer_boundary(i,false));
		primal_set[i][0] = qc_create_set_adv(dp_qc->qc, PRIMAL, code.num_qubits + i, 0, 2, infer_boundary(i,false));

		qc_insert_set(dp_qc->qc, primal_set[i][1]);
		qc_associate_syndrome(primal_set[i][1], primal_syndrome[i]);
		qc_insert_set(dp_qc->qc, primal_set[i][0]);
	}

	for(int i = 0; i < code.Z_stabilizer.size(); ++i){
		dual_syndrome[i] = qc_create_syndrome();
		qc_insert_syndrome(dp_qc->qc, dual_syndrome[i]);

		//debug
		//std::cout << "create set with boundary "<< infer_test_boundary(i,true) << " and type:" << -1 << " at id "<< i <<" and position " << i << std::endl;

		dual_set[i][1] = qc_create_set_adv(dp_qc->qc, DUAL, code.num_qubits + code.X_stabilizer.size() + i, 0, 1, infer_boundary(i,true));
		dual_set[i][0] = qc_create_set_adv(dp_qc->qc, DUAL, code.num_qubits + code.X_stabilizer.size() + i, 0, 2, infer_boundary(i,true));

		qc_insert_set(dp_qc->qc, dual_set[i][1]);
		qc_associate_syndrome(dual_set[i][1], dual_syndrome[i]);
		qc_insert_set(dp_qc->qc, dual_set[i][0]);
	}
}


void Wrapper::init_X_stabilizers(){
	// init all X stabilizer qubits
	for(int i = code.num_qubits; i < code.num_qubits + code.X_stabilizer.size(); ++i){
		dp_init_Z(dp_qc,qubit_array[i]); //first initialize in Z
		dp_H(dp_qc,qubit_array[i]); // apply a Hadamard 
	}
}


void Wrapper::init_Z_stabilizers(){
	//init all Z stabilizer qubits
	for(int i = code.num_qubits + code.X_stabilizer.size(); i < code.num_qubits + code.X_stabilizer.size() + code.Z_stabilizer.size(); ++i){
		dp_dead_H(dp_qc,qubit_array[i]); // apply a dead Hadamard 
		dp_init_Z(dp_qc,qubit_array[i]); //first initialize in Z
	}
}

void Wrapper::data_error_during_init(){
	for (int i = 0; i < code.num_qubits; ++i){
		dp_iden_init_Z(dp_qc, qubit_array[i]);
		dp_dead_H(dp_qc, qubit_array[i]);
	}	
}


void Wrapper::apply_cnots(){
	StabilizerContainer X_queue(code.X_stabilizer);
	StabilizerContainer Z_queue(code.Z_stabilizer);
	bool finished = false;
	int timestep;
	while(finished==false){
		timestep = qubit_array[0]->t;
		finished = true;

		for(int pos = 0; pos < code.X_stabilizer.size(); ++pos){
			for(auto qubit : X_queue.at(pos)){
				if(qubit_array[qubit]->t == timestep){
					dp_cnot(dp_qc, qubit_array[code.num_qubits + pos], qubit_array[qubit]);
					X_queue.at(pos).erase(qubit);
					finished = false;
					break;
				}
			}
		}
		for(int pos = 0; pos < code.Z_stabilizer.size(); ++pos){
			for(auto qubit : Z_queue.at(pos)){
				if(qubit_array[qubit]->t == timestep){
					dp_cnot(dp_qc, qubit_array[qubit], qubit_array[pos + code.num_qubits + code.X_stabilizer.size()]);
					Z_queue.at(pos).erase(qubit);
					finished = false;
					break;
				}
			}
		}

		//bring everything to the same timestep
		if (finished)
			break;
		// Bring all qubits to the same time by applying identities
		for(int i=0; i < code.num_qubits + code.Z_stabilizer.size() + code.X_stabilizer.size(); ++i){
			if(qubit_array[i]->t <= timestep){
				dp_iden_cnot(dp_qc,qubit_array[i]);
			}
		}

	} // end while
}


void Wrapper::measure_X_stabilizers(size_type big_t){
	size_type lay1 = (big_t+1)%2;
	size_type lay2 = (big_t)%2;
	// measurement
	for(int pos = 0; pos < code.X_stabilizer.size(); ++pos){
		dp_meas_Z(dp_qc, qubit_array[pos + code.num_qubits], primal_set[pos][lay1], primal_set[pos][lay2]);

		qc_unassociate_syndrome(primal_syndrome[pos]);
		qc_associate_syndrome(primal_set[pos][lay2], primal_syndrome[pos]);
		//create new sets and check for boundaries
		primal_set[pos][lay1] = qc_create_set_adv(dp_qc->qc, PRIMAL, code.num_qubits + pos, 0, 2, infer_boundary(pos,false));

		qc_insert_set(dp_qc->qc, primal_set[pos][lay1]);
	}

	for(int i = 0; i < code.num_qubits; ++i){
		dp_iden_meas_Z(dp_qc,qubit_array[i]);
	}
}



void Wrapper::measure_Z_stabilizers(size_type big_t){
	size_type lay1 = (big_t+1)%2;
	size_type lay2 = (big_t)%2;

	for(int pos = 0; pos < code.Z_stabilizer.size(); ++pos){
		dp_meas_Z(dp_qc, qubit_array[code.num_qubits + code.X_stabilizer.size() + pos], dual_set[pos][lay1], dual_set[pos][lay2]);
		qc_unassociate_syndrome(dual_syndrome[pos]);
		qc_associate_syndrome(dual_set[pos][lay2], dual_syndrome[pos]);
		
		//create new sets
		dual_set[pos][lay1] = qc_create_set_adv(dp_qc->qc, DUAL, code.num_qubits + code.X_stabilizer.size() + pos, 0, 2, infer_boundary(pos,true));
		qc_insert_set(dp_qc->qc, dual_set[pos][lay1]);
	}

	// Apply Hadamard for X stabilizers
	for(int i = code.num_qubits; i < code.num_qubits + code.X_stabilizer.size(); ++i){
		dp_H(dp_qc,qubit_array[i]);
	}
	for(int i = 0; i < code.num_qubits; ++i){
		dp_iden_H(dp_qc,qubit_array[i]);
	}
	// apply dead syndrome due to earlier measuremnt
	for(int i = code.num_qubits + code.X_stabilizer.size(); i < code.num_qubits + code.X_stabilizer.size() + code.Z_stabilizer.size(); ++i){
		dp_dead_H(dp_qc,qubit_array[i]);
	}


}

void Wrapper::measure_stabilizers(size_type big_t){
	// loop through syndromes
	// unfortunately lattice is frustrated in general, so assume worst case
	// where all syndrome block each other
	init_X_stabilizers();
	init_Z_stabilizers();

	data_error_during_init();
	
	apply_cnots();

	measure_Z_stabilizers(big_t);
	measure_X_stabilizers(big_t);

	assert(dp_qc->qc->big_t == big_t+1);
	return;
}


Wrapper::Wrapper(){
	frame = NULL;
	qubit_array	= NULL;
	primal_boundary = NULL;
	dual_boundary = NULL;
	primal_syndrome = NULL;
	dual_syndrome = NULL;
	temporal_primal_boundary = NULL;
	temporal_dual_boundary = NULL;
	boundaries = NULL;
	primal_set = NULL;
	dual_set = NULL;
	dp_qc = NULL;
	copy = false;
	distance = 10;
	num_checks = 0;
	probab = 0;
}

Wrapper::Wrapper(const Wrapper & a): Wrapper(){
	copy = true;

	distance = a.distance;
	num_checks = a.num_checks;
	num_X_changes = a.num_X_changes;
	num_Z_changes = a.num_Z_changes;
	last_X_check = a.last_X_check;
	last_Z_check = a.last_Z_check;


	const size_type num_ancillae = code.Z_stabilizer.size() + code.X_stabilizer.size();

	//allocate all data
	frame = (int *)my_malloc((code.num_qubits + num_ancillae) * sizeof(int));
	qubit_array = (QUBIT**) my_calloc(code.num_qubits + num_ancillae, sizeof(QUBIT *));

	primal_syndrome = (SYNDROME**)my_calloc(code.X_stabilizer.size(), sizeof(SYNDROME *));
	dual_syndrome = (SYNDROME**)my_calloc(code.Z_stabilizer.size(), sizeof(SYNDROME *));

	primal_set = (SET ***)my_2d_calloc(code.X_stabilizer.size(), 2, sizeof(SET *));
	dual_set = (SET ***)my_2d_calloc(code.Z_stabilizer.size(), 2, sizeof(SET *));

	//copy the frame
	for(int i = 0; i < code.num_qubits + num_ancillae; ++i){
		frame[i] = a.frame[i];
	}

	//copy the depolarizing quantum computer
	dp_qc = dp_copy_dp_qc(a.dp_qc);

	//copy the qubits
	for(int i = 0; i < code.num_qubits + num_ancillae; ++i){
		qubit_array[i] = a.qubit_array[i]->copy;
	}

	//copy the PRIMAL syndrome and set arrays
	for(int i = 0; i < code.X_stabilizer.size(); ++i){
		primal_syndrome[i] = a.primal_syndrome[i]->copy;
		primal_set[i][0] = a.primal_set[i][0]->copy;
		primal_set[i][1] = a.primal_set[i][1]->copy;
	}
	//copy the DUAL syndrome and set arrays
	for(int i = 0; i < code.Z_stabilizer.size(); ++i){	
		dual_syndrome[i] = a.dual_syndrome[i]->copy;
		dual_set[i][0] = a.dual_set[i][0]->copy;
		dual_set[i][1] = a.dual_set[i][1]->copy;
	}

	boundaries = a.boundaries;
	primal_boundary = a.primal_boundary;
	dual_boundary = a.dual_boundary;
	temporal_primal_boundary = a.temporal_primal_boundary;
	temporal_dual_boundary = a.temporal_dual_boundary;
}


Wrapper::Wrapper(double error_probability) : Wrapper(){
	probab = error_probability;
	copy = false;

	//create recipe-- whatever that does
	recipe = qc_create_recipe_adv(RECIPE_INFINITE, code.num_qubits + code.Z_stabilizer.size() + code.X_stabilizer.size(), 1, false);

	//allocate data containers
	allocate();
	//init for the first time
	init(NULL);
}

void Wrapper::allocate(){
	const size_type num_ancillae = code.Z_stabilizer.size() + code.X_stabilizer.size();
	delete_data();
	//allocation
	frame = (int *)my_malloc((code.num_qubits + num_ancillae) * sizeof(int));
	qubit_array = (QUBIT**) my_calloc(code.num_qubits + num_ancillae, sizeof(QUBIT *));
	primal_boundary = (SET**) my_calloc(code.X_boundary.size(), sizeof(SET *));
	dual_boundary = (SET**) my_calloc(code.Z_boundary.size(), sizeof(SET *));

	primal_syndrome = (SYNDROME**)my_calloc(code.X_stabilizer.size(), sizeof(SYNDROME *));
	dual_syndrome = (SYNDROME**)my_calloc(code.Z_stabilizer.size(), sizeof(SYNDROME *));

	primal_set = (SET ***)my_2d_calloc(code.X_stabilizer.size(), 2, sizeof(SET *));
	dual_set = (SET ***)my_2d_calloc(code.Z_stabilizer.size(), 2, sizeof(SET *));

	temporal_primal_boundary = (SET**) my_calloc(2,sizeof(SET *));
	temporal_dual_boundary = (SET**) my_calloc(2,sizeof(SET *));

	boundaries = (BALL **)my_malloc((code.X_boundary.size() + code.Z_boundary.size() + 4) * sizeof(BALL *));
}

void Wrapper::init(RECIPE_ADV *rec){
	num_X_changes = 0;
	num_Z_changes = 0;
	last_X_check = 0;
	last_Z_check = 0;
	const size_type num_ancillae = code.Z_stabilizer.size() + code.X_stabilizer.size();
	//initialization
	dp_qc = dp_create_dp_qc_adv(code.seed0, code.seed1, DE_HT_FACTOR*(code.num_qubits+num_ancillae)*(code.num_qubits+num_ancillae),
		STICK_HT_FACTOR*distance*distance, probab, code.t_delete, rec, code.dir);
	//create qubits
	for (int i = 0; i < code.num_qubits+num_ancillae; ++i) {
		qubit_array[i] = qc_create_and_insert_qubit(dp_qc->qc, i, 0, 0, QUBIT_HT_SIZE);
	}

	// first primal
	for(int i = 0; i < code.X_boundary.size(); ++i){
		primal_boundary[i] = create_boundary_set(dp_qc->qc, i, PRIMAL_BOUNDARY, -1 - i, -1, 0);
	}
	// second dual
	for(int i = 0; i < code.Z_boundary.size(); ++i){
		dual_boundary[i] = create_boundary_set(dp_qc->qc, i + code.X_boundary.size(), DUAL_BOUNDARY, -1 - i, -2, 0);
	}

	//temporal
	temporal_primal_boundary[0] = create_boundary_set(dp_qc->qc, code.Z_boundary.size() + code.X_boundary.size(), PRIMAL_BOUNDARY, -1, -3, 0);
	temporal_primal_boundary[1] = create_boundary_set(dp_qc->qc, code.Z_boundary.size() + code.X_boundary.size() + 1, PRIMAL_BOUNDARY, -2, -3, 0);

	temporal_dual_boundary[0] = create_boundary_set(dp_qc->qc, code.Z_boundary.size() + code.X_boundary.size()+2, DUAL_BOUNDARY, -1, -4, 0);
	temporal_dual_boundary[1] = create_boundary_set(dp_qc->qc, code.Z_boundary.size() + code.X_boundary.size()+3, DUAL_BOUNDARY, -2, -4, 0);

	dp_qc->qc->get_boundary = get_boundary;
	dp_qc->qc->boundaries = (void*)boundaries;

	create_initial_syndrome_and_set_arrays();
}

void Wrapper::delete_data(){
	if(dp_qc){
		if(copy){
			dp_free_dp_qc_copy(dp_qc);
		}
		else{
			dp_free_dp_qc(dp_qc);
		}
	}


	if(frame)
		free(frame);
	if(qubit_array);
		free(qubit_array);

	if(primal_syndrome){
		for(int i = 0; i < code.X_stabilizer.size(); ++i){
			qc_free_syndrome(primal_syndrome[i]);
		}
		free(primal_syndrome);
	}

	if(dual_syndrome){
		for(int i=0; i< code.Z_stabilizer.size(); ++i){
			qc_free_syndrome(dual_syndrome[i]);
		}
		free(dual_syndrome);
	}
	
	
	if(primal_set)
		my_2d_free(code.X_stabilizer.size(), (void**)primal_set);
	if(dual_set)
		my_2d_free(code.Z_stabilizer.size(), (void**)dual_set);


	//no need to free boundaries for copy: these are not copied
	if(!copy){
		if(primal_boundary){
			for(int i=0; i < code.X_boundary.size(); ++i)
				qc_free_bdy(primal_boundary[i]);
		}

		if(dual_boundary){
			for(int i=0; i < code.Z_boundary.size(); ++i)
				qc_free_bdy(dual_boundary[i]);
		}


		if(temporal_primal_boundary){
			qc_free_bdy(temporal_primal_boundary[0]);
			qc_free_bdy(temporal_primal_boundary[1]);
		}
		if(temporal_dual_boundary){
			qc_free_bdy(temporal_dual_boundary[0]);
			qc_free_bdy(temporal_dual_boundary[1]);		
		}
	
		if(boundaries)
			free(boundaries);
	}


	dp_qc = NULL;
	frame = NULL;
	qubit_array	= NULL;
	primal_syndrome = NULL;
	dual_syndrome = NULL;
	temporal_primal_boundary = NULL;
	temporal_dual_boundary = NULL;
	primal_set = NULL;
	dual_set = NULL;

	if(!copy){
		primal_boundary = NULL;
		dual_boundary = NULL;
		boundaries = NULL;
	}
}

Wrapper::~Wrapper(){
	delete_data();
	if(!copy)
		qc_free_recipe_adv(recipe);
}




void Wrapper::generate_recipe(){
	size_type big_t = 0;
	do{
		measure_stabilizers(big_t++);
	}while (qc_boot_up_adv(dp_qc->qc, recipe, 2, 12) != DONE);
	return;
}


void Wrapper::run_simulation(size_type big_t_max){
	allocate();
	init(recipe);
	
	//simulation parameters

	//output basic data
	std::cout << "p = " << dp_qc->p << "\n"
			<< "d = " << std::min(code.Z_operator.at(0).size(), code.X_operator.at(0).size()) << "\n"
			<< "big_t_max = " << big_t_max << "\n"
			<< "new t_check = " << code.t_check << "\n"
			<< "t_delete = " << code.t_delete << "\n"
			<< "max_num_X = " << code.max_num_X << "\n"
			<< "max_num_Z = " << code.max_num_Z << "\n"
			<< "verbose = " << false << "\n"
			<< "s0 = " << code.seed0 << "\n"
			<< "s1 = " << code.seed1 << "\n"
			<< "-s0 " << code.seed0 << " -s1 " << code.seed1 <<std::endl;

	size_type big_t = 0;

	// Disable error tracking
	dp_qc->qc->track = false;
	
	// Let QC know about the boundaries
	dp_qc->qc->get_boundary = get_boundary;
	dp_qc->qc->boundaries = (void *)boundaries;

	
	while (big_t <= big_t_max && (num_X_changes < code.max_num_X || num_Z_changes < code.max_num_Z)) {
		measure_stabilizers(big_t);

		qc_convert_nests(dp_qc->qc, false);

		qc_mwpm(dp_qc->qc, false);


		correct_mts();
		// not quite sure why following needs -2 rather than -1
		qc_trim_nests(dp_qc->qc, dp_qc->qc->unfinalized_big_t - 20);
		m_time_delete(dp_qc->qc->m_pr);
		m_time_delete(dp_qc->qc->m_du);
		
		//TODO uncomment
		if (big_t > 0 && big_t%code.t_check == 0) {
			test_correct();
		}
		
		assert(qc->big_t == big_t+1);
		big_t++;
	}
}


void Wrapper::calculate_t_check(){
	allocate();
	init(recipe);
	dp_qc->qc->track = false;
	dp_qc->qc->get_boundary = get_boundary;
	
	size_type big_t = 0;
	size_type t_check = code.t_check;
	size_type t_check_scale = 10;
	size_type boot_num_X = 15;
	size_type boot_num_Z = 15;


	// basic data
	std::cout << "p = " << dp_qc->p << "\n"
		<< "d = " << std::min(code.Z_operator.at(0).size(), code.X_operator.at(0).size()) << "\n"
		<< "big_t_max = " << code.big_t_max << "\n"
		<< "t_check = " << code.t_check << "\n"
		<< "t_delete = " << code.t_delete << "\n"
		<< "max_num_X = " << code.max_num_X << "\n"
		<< "max_num_Z = " << code.max_num_Z << "\n"
		<< "boot = " << true << "\n"
		<< "boot_num_X = " << boot_num_X << "\n"
		<< "boot_num_Z = " << boot_num_Z << "\n"
		<< "t_check_scale = " << t_check_scale << "\n"
		<< "verbose = " << false << "\n"
		<< "s0 = " << code.seed0 << "\n"
		<< "s1 = " << code.seed1 << "\n"
		<< "-s0 " << code.seed0 << " -s1 " << code.seed1 <<std::endl;



	while (true) {
		measure_stabilizers(big_t);
		qc_convert_nests(dp_qc->qc, false);
		qc_mwpm(dp_qc->qc, false);
		correct_mts();
		qc_trim_nests(dp_qc->qc, dp_qc->qc->unfinalized_big_t - 20);
		m_time_delete(dp_qc->qc->m_pr);
		m_time_delete(dp_qc->qc->m_du);


		if (big_t > 0 && (big_t % t_check) == 0) {
			test_correct();	

			// If we have gotten to 50 time checks and have yet to find a change, 
			// double t_check or we'll be here all year.
			if (num_X_changes == 0 && num_Z_changes == 0 && num_checks >= 50) {
				t_check <<= 1;
				std::cout << "Doubling tcheck: " << t_check << std::endl;
				num_checks = 0;
			}

			if (num_X_changes >= boot_num_X || num_Z_changes >= boot_num_Z) {
				int changes = (num_X_changes >= num_Z_changes) ? num_X_changes : num_Z_changes; 
				std::cout << "Calculating tcheck: " << t_check << " "  << num_checks << " " << changes << " " << t_check_scale <<std::endl;

				code.t_check = t_check * num_checks / changes / t_check_scale;

				if (code.t_check < 1) { 
					code.t_check = 1;
				}
				
				break;
			}
		}

		// Ensure that the qc has been properly advanced in big_t
		assert(dp_qc->qc->big_t == big_t+1);
		
		big_t++;
	}
	// Free the sc_dp_qc now that the t_check has been calculated
	delete_data();
	qc_reset_recipe_adv(recipe);
}

void Wrapper::reset_recipe(){

}


void Wrapper::test_correct() {
	int count;
	long big_t = 0;
	Wrapper second(*this);

	// Enable perfect gates and perform a perfect round of stabilizer
	// measurements, followed by matching, then correction.
	second.dp_qc->qc->perfect_gates = true;

	second.dp_qc->qc->m_pr->undo_flag = true;
	second.dp_qc->qc->m_du->undo_flag = true;

	second.measure_stabilizers(second.dp_qc->qc->big_t);

	qc_finalize_nests(second.dp_qc->qc, second.dp_qc->qc->big_t);
	qc_convert_nests(second.dp_qc->qc, true);
	
	
	qc_mwpm(second.dp_qc->qc, true);

	second.correct_mts();

	// Track how many test corrects have been performed for this simulation
	num_checks++;

	// Count the number of Z errors that have occured along the vertical
	// boundary. An Z-error is a disagreement between the state of the qubit in
	// the qc and the state of the correction in the frame. 
	//
	// If there is a Z error on the qubit, but no Z correction, or the reverse
	// with a Z correction but no Z error on the qubit, then an error along
	// the boundary has occured.

	// calculate the parity along a path that defines a logical qubit
	for(auto logical_qubit : code.Z_operator){
		count = 0;
		for(auto id : logical_qubit){
			if ((dp_contains_Z(second.qubit_array[id]->e) && ((second.frame[id]&Z) != Z)) ||
				(!dp_contains_Z(second.qubit_array[id]->e) && ((second.frame[id]&Z) == Z))) {
				++count;
			}
		}
		//check if the parity is the same
		if (count % 2 != last_Z_check){
			last_Z_check = count%2;
			num_Z_changes++;
			print_stats();
		}
	}


	for(auto logical_qubit : code.X_operator){
		count = 0;
		for(auto id : logical_qubit){
			if ((dp_contains_X(second.qubit_array[id]->e) && ((second.frame[id]&X) != X)) ||
				(!dp_contains_X(second.qubit_array[id]->e) && ((second.frame[id]&X) == X))) {
				++count;
			}
		}
		//check if the parity is the same
		if (count % 2 != last_X_check){
			last_X_check = count%2;
			num_X_changes++;
			print_stats();
		}
	}


	qc_undo_mwpm(second.dp_qc->qc);

	dp_qc->qc->m_pr->undo_flag = false;
	dp_qc->qc->m_du->undo_flag = false;
}


void Wrapper::print_stats() {
	std::cout << ">>> " << 0 << " | " << dp_qc->qc->big_t << " | Last: "
				<< last_X_check << " " << last_Z_check << " | Checks: " << num_checks << " " << num_checks
				<< " | Changes: " << num_X_changes << " " << num_Z_changes << std::endl;
}



void Wrapper::correct_mts() {
	process_aug_edges(dp_qc->qc->m_pr);
	process_aug_edges(dp_qc->qc->m_du);
}

void Wrapper::process_aug_edge(AUG_EDGE *ae) {
	// Performs correction given two endpoints of a matched error chain
	// Get the coordinates at the start and end of the edge (come from boundary coordinates and set coordinates)
	// ae->va: source vertex of the edge
	// ae->vb: destination vertex of the edge
	int i1 = ae->va->i;
	int j1 = ae->va->j;
	int	i2 = ae->vb->i;
	int	j2 = ae->vb->j;
	int error = 0;
	bool boundary = false;
	ParityCheck end; //list of end positions (only multiple elments when boundary is involved)
	// first resolve boundary to proper coordinates
	if(j1 < 0){
		//only one boundary possible -> always have the boundary as destinaltion
		std::swap(i1,i2);
		std::swap(j1,j2);
	}
	
	// now resolve the boundary
	if(j2 < 0){
		// boundary
		if(j2 <= -3){ // temporal boundary
			return;
		}
		boundary = true;
		if(j2 == -1){ // primal boundary: set to closest qubit position
			for(auto element : code.X_boundary.at(-i2-1)){
				end.insert(element);
			}
			error = Z; // primal has X-stabilizers -> detection of Z errors
		}
		if(j2 == -2){ // dual boundary: set to closest qubit position
			for(auto element : code.Z_boundary.at(-i2-1)){
				end.insert(element);
			}
			error = X; // dual has Z-stabilizers -> detection of X errors
		}
	}
	else{
		// not on the boundary
		if (i2 >= code.num_qubits && i2 < code.num_qubits + code.X_stabilizer.size()){
			end.insert(i2-code.num_qubits);
			error = X;
		} else if(i2 >= code.num_qubits + code.X_stabilizer.size()){ // dual
			end.insert(i2 - code.num_qubits - code.X_stabilizer.size());
			error = Z;
		}
		else{
			std::cout << "i2: " << i2 <<  std::endl;
			std::cout << "process_edge: no proper ids" <<std::endl;
			throw(1);
		}
	}

	// starting postion
	if (i1 >= code.num_qubits && i1 < code.num_qubits + code.X_stabilizer.size()){
		i1 = i1-code.num_qubits;
		assert(error == X);
	} else if(i1 >= code.num_qubits + code.X_stabilizer.size()){ // dual
		i1 = i1 - code.num_qubits - code.X_stabilizer.size();
		assert(error == Z);
	}
	else{
		std::cout << "i1: " << i1 << " i2: " << i2 << std::endl;
		std::cout << "process_edge: no proper ids" <<std::endl;
		throw(1);
	}


	// perform pathfinding and obtain data qubits with an error
	ParityCheck path;
	if(error == X){
		DataQubits(code.X_stabilizer, i1, end, path);
	} else{ //dual syndromes
		DataQubits(code.Z_stabilizer, i1, end, path);
	}
	
	//add last element on boundary
	if(boundary){
		if(j1 == -1){
			path.insert(code.GetBoundaryQubitFromStabilizer(code.XBoundaryQubits,code.X_stabilizer.at(*end.begin())));
		}
		if(j1 == -2){
			path.insert(code.GetBoundaryQubitFromStabilizer(code.ZBoundaryQubits,code.Z_stabilizer.at(*end.begin())));
		}
	}

	// apply an error all elements of the frame along the path
	for(auto id : path){
		frame[id] ^= error;
	}
	return;
}


void Wrapper::process_aug_edges(MATCHING *m) {
	AUG_EDGE *ae;
	ae = m_get_aug_edge(m);
	while (ae != NULL) {
		process_aug_edge(ae);
		m_delete_aug_edge(ae);
		ae = m_get_aug_edge(m);
	}
}
