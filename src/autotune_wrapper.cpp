#include <autotune_wrapper.hpp>
#include <string.h>
#include <path_finding.hpp>


#define LONG_MAX 100
#define DE_HT_FACTOR 100
#define STICK_HT_FACTOR 200
#define QUBIT_HT_SIZE 1000
// This rearranges example 1 such that it can be used with hyperbolic codes

SET* Wrapper::create_boundary_set(QC *qc, int id, int type, int i, int j, int t) {
	SET *set;

	set = qc_create_set_adv(qc, type, i, j, t, NULL);
	boundaries[id] = set->ball;

	return set;
}

BALL *get_boundary(int i, int j, long int big_t, int type, void *boundaries) {
	BALL **bdys;
	std::cout<< "get_boundaries... is this function called?" << std::endl;
	// look up boundaries
	// TODO


	bdys = (BALL **)boundaries;

	// Boundaries do not have a big_t set, they are initialised to LONG_MAX.
	if (big_t != LONG_MAX) {
		return NULL;
	}



	if(j == 1){
		if(i ==0 || i == 1){
			if(type == PRIMAL_BOUNDARY){
				return bdys[code.X_boundary.size() + code.Z_boundary.size() + i];
			}
			return bdys[code.X_boundary.size() + code.Z_boundary.size() + 2 + i];
		}
	} else{
		if(type == PRIMAL_BOUNDARY){
			if (i < code.Z_boundary.size())
				return bdys[i];
		}
		else{
			if (i < code.X_boundary.size())
				return bdys[code.Z_boundary.size() + i];
		}
	}

	// There is no boundary at this coordinate
	return NULL;
}

void Wrapper::correct_mts() {
	process_aug_edges(dp_qc->qc->m_pr);
	process_aug_edges(dp_qc->qc->m_du);
}


SET* Wrapper::infer_boundary(size_type pos, bool is_dual){
	//infer from the stabilizer id which boundary is needed
	if(is_dual){
		for(int i=0; i < code.Z_boundary.size(); ++i){
			if (code.Z_boundary.at(i).find(pos) != code.Z_boundary.at(i).end()){
				return primal_boundary[i];
			}
		}
	} else{
		for(int i=0; i < code.X_boundary.size(); ++i){
			if (code.X_boundary.at(i).find(pos) != code.X_boundary.at(i).end()){
				return dual_boundary[i];
			}
		}
	}
	return NULL;
}


void Wrapper::create_initial_syndrome_and_set_arrays(){

	for(int i = 0; i < code.X_stabilizer.size(); ++i){
		primal_syndrome[i] = qc_create_syndrome();
		qc_insert_syndrome(dp_qc->qc, primal_syndrome[i]);

		primal_set[i][1] = qc_create_set_adv(dp_qc->qc, PRIMAL, i, 0, 1, infer_boundary(i,true));
		primal_set[i][0] = qc_create_set_adv(dp_qc->qc, PRIMAL, i, 0, 2, infer_boundary(i,true));


		qc_insert_set(dp_qc->qc, primal_set[i][1]);
		qc_associate_syndrome(primal_set[i][1], primal_syndrome[i]);
		qc_insert_set(dp_qc->qc, primal_set[i][0]);
	}

	for(int i = 0; i < code.Z_stabilizer.size(); ++i){
		dual_syndrome[i] = qc_create_syndrome();
		qc_insert_syndrome(dp_qc->qc, dual_syndrome[i]);

		dual_set[i][1] = qc_create_set_adv(dp_qc->qc, DUAL, i, 1, 1, infer_boundary(i,false));
		dual_set[i][0] = qc_create_set_adv(dp_qc->qc, DUAL, i, 1, 2, infer_boundary(i,false));

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

		// for worst case each ancilla qubit needs to wait 2 steps
		// each qubit has 2 X and 2 Z-stabilizers, due to the hadamard operation they are shifted by 1
		// thus 1 step definitely works, for a clash one of the choices also works thus the worst case has a delay of 2
		for(int j=0; j < num_padding; ++j)
			dp_iden_cnot(dp_qc,qubit_array[i]);
	}
}


void Wrapper::init_Z_stabilizers(){
	//init all Z stabilizer qubits
	for(int i = code.num_qubits + code.X_stabilizer.size(); i < code.num_qubits + code.X_stabilizer.size() + code.Z_stabilizer.size(); ++i){
		dp_init_Z(dp_qc,qubit_array[i]); //first initialize in Z
		dp_dead_H(dp_qc,qubit_array[i]); // apply a Hadamard 

		// for worst case each ancilla qubit needs to wait 2 steps
		// each qubit has 2 X and 2 Z-stabilizers, due to the hadamard operation they are shifted by 1
		// thus 1 step definitely works, for a clash one of the choices also works thus the worst case has a delay of 2
		for(int j=0; j < num_padding; ++j)
			dp_iden_cnot(dp_qc,qubit_array[i]);
	}
}

void Wrapper::data_error_during_init(){
	// all data qubits need to wait 1 time step due to initialization
	// only 1 means a shifted initialization for X and Z stabilizers, such that the data qubits are best distributed
	for (int i = 0; i < code.num_qubits; ++i){
		dp_iden_init_Z(dp_qc, qubit_array[i]);
		dp_dead_H(dp_qc, qubit_array[i]);

		for(int j=0; j < num_padding; ++j){
			dp_dead_cnot(dp_qc,qubit_array[i]);
		}
	}	
}



void Wrapper::perform_X_stabilizers(){
	// now apply CNOTS for X measurements
	for (int pos = 0; pos < code.X_stabilizer.size(); ++pos){
		for(auto qubit : code.X_stabilizer.at(pos)){
			dp_cnot(dp_qc, qubit_array[qubit], qubit_array[code.num_qubits + pos]);
			
			//unfortunately autotune requires the same timestep... so we need to add a dead_identity operation to all other qubits
			for(int i=0; i < code.num_qubits + code.Z_stabilizer.size() + code.X_stabilizer.size(); ++i){
				if(i != qubit && i != pos+code.num_qubits)
					dp_dead_cnot(dp_qc,qubit_array[i]);
			}
		}
	}

	// preparation for the measurment of X_stabilizers
	for(int i = code.num_qubits; i < code.num_qubits + code.X_stabilizer.size(); ++i){
		dp_H(dp_qc,qubit_array[i]);
	}
	// apply dead identity to all other qubits
	for(int i = 0; i < code.num_qubits; ++i){
		dp_dead_H(dp_qc,qubit_array[i]);
	}
	for(int i = code.num_qubits + code.X_stabilizer.size(); i < code.num_qubits + code.X_stabilizer.size() + code.Z_stabilizer.size(); ++i){
		dp_dead_H(dp_qc,qubit_array[i]);
	}
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
		primal_set[pos][lay1] = qc_create_set_adv(dp_qc->qc, PRIMAL, 0, 0, 2, infer_boundary(pos,false));

		qc_insert_set(dp_qc->qc, primal_set[pos][lay1]);
	}
}



void Wrapper::perform_Z_stabilizers(){
	// apply CNOTS for Z measurements
	for (int pos = 0; pos < code.Z_stabilizer.size(); ++pos){
		for(auto qubit : code.Z_stabilizer.at(pos)){
			dp_cnot(dp_qc, qubit_array[pos + code.num_qubits + code.X_stabilizer.size()], qubit_array[qubit]);
			
			//unfortunately autotune requires the same timestep... so we need to add a dead_identity operation to all other qubits
			for(int i=0; i < code.num_qubits + code.Z_stabilizer.size() + code.X_stabilizer.size(); ++i){
				if(i != qubit && i != pos + code.num_qubits + code.X_stabilizer.size())
					dp_dead_cnot(dp_qc,qubit_array[i]);
			}
		}
	}
}



void Wrapper::measure_Z_stabilizers(size_type big_t){
	size_type lay1 = (big_t+1)%2;
	size_type lay2 = (big_t)%2;

	// measurement
	// TODO: change to here

	for(int pos = 0; pos < code.Z_stabilizer.size(); ++pos){
		dp_meas_Z(dp_qc, qubit_array[code.num_qubits + code.X_stabilizer.size() + pos], dual_set[pos][lay1], dual_set[pos][lay2]);
		qc_unassociate_syndrome(dual_syndrome[pos]);
		qc_associate_syndrome(dual_set[pos][lay2], dual_syndrome[pos]);
		
		//create new sets
		dual_set[pos][lay1] = qc_create_set_adv(dp_qc->qc, DUAL, 0, 0, 2, infer_boundary(pos,true));

		qc_insert_set(dp_qc->qc, dual_set[pos][lay1]);
	}
}

void Wrapper::measure_stabilizers(size_type big_t){
	// loop through syndromes
	// unfortunately lattice is fustrated in general, so assume worst case
	// where all syndrome block each other

	init_X_stabilizers();
	init_Z_stabilizers();
	data_error_during_init(); // add dead operators for all qubits to be at the same time step

	perform_X_stabilizers();
	
	perform_Z_stabilizers();

	measure_X_stabilizers(big_t);
	measure_Z_stabilizers(big_t);

	//one dead idendity for measurements (the errors have been applied in the data_error_during_init() step)
	for(int i = 0; i < code.num_qubits; ++i){
		dp_dead_meas_Z(dp_qc,qubit_array[i]);
	}

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

	seed0 = 42;
	seed1 = 43;
	distance = 10;
	t_delete = 100; //delete after 100 steps
	num_checks = 0;
}

Wrapper::Wrapper(const Wrapper & a): Wrapper(){
	//copy directory for the error model
	dir = (char*)malloc(100*sizeof(char));
	strcpy(dir, a.dir);

	const size_type num_ancillae = code.Z_stabilizer.size() + code.X_stabilizer.size();

	//allocate all data
	allocate();

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
		primal_syndrome[i] = a.primal_syndrome[i];
		primal_set[i][0] = a.primal_set[i][0];
		primal_set[i][1] = a.primal_set[i][1];
	}
	//copy the DUAL syndrome and set arrays
	for(int i = 0; i < code.Z_stabilizer.size(); ++i){	
		dual_syndrome[i] = a.dual_syndrome[i];
		dual_set[i][0] = a.primal_set[i][0];
		dual_set[i][1] = a.primal_set[i][1];
	}
	//copy the boundaries
	for(int i = 0; i < code.X_boundary.size(); ++i){
		primal_boundary[i] = a.primal_boundary[i];
	}
	// second dual
	for(int i = 0; i < code.Z_boundary.size(); ++i){
		dual_boundary[i] = dual_boundary[i];
	}

	temporal_primal_boundary[0] = a.temporal_primal_boundary[0];
	temporal_primal_boundary[1] = a.temporal_primal_boundary[1];

	temporal_dual_boundary[0] = a.temporal_dual_boundary[0];
	temporal_dual_boundary[1] = a.temporal_dual_boundary[1];
}


Wrapper::Wrapper(double error_probability) : Wrapper(){
	dir = (char*)malloc(100*sizeof(char));
	strcpy(dir, "../../autotune/ex/ems/");


	//create recipe-- whatever that does
	recipe = qc_create_recipe_adv(RECIPE_INFINITE, 100, 100, false);

	//allocate data containers
	allocate();
	//init for the first time
	init(NULL, error_probability);
}

void Wrapper::allocate(){
	const size_type num_ancillae = code.Z_stabilizer.size() + code.X_stabilizer.size();
	delete_data();
	//allocation
	frame = (int *)my_malloc((code.num_qubits + num_ancillae) * sizeof(int));
	qubit_array = (QUBIT**) my_calloc(code.num_qubits + num_ancillae, sizeof(QUBIT *));
	primal_boundary = (SET**) my_calloc(2*code.X_boundary.size(), sizeof(SET *));
	dual_boundary = (SET**) my_calloc(2*code.Z_boundary.size(), sizeof(SET *));

	primal_syndrome = (SYNDROME**)my_calloc(code.X_stabilizer.size(), sizeof(SYNDROME *));
	dual_syndrome = (SYNDROME**)my_calloc(code.Z_stabilizer.size(), sizeof(SYNDROME *));

	primal_set = (SET ***)my_2d_calloc(code.X_stabilizer.size(), 2, sizeof(SET *));
	dual_set = (SET ***)my_2d_calloc(code.Z_stabilizer.size(), 2, sizeof(SET *));

	temporal_primal_boundary = (SET**) my_calloc(2,sizeof(SET *));
	temporal_dual_boundary = (SET**) my_calloc(2,sizeof(SET *));
}

void Wrapper::init(RECIPE_ADV *rec, double error_probability){
	num_X_changes = 0;
	num_Z_changes = 0;
	last_X_check = 0;
	last_Z_check = 0;
	const size_type num_ancillae = code.Z_stabilizer.size() + code.X_stabilizer.size();
	//initialization
	dp_qc = dp_create_dp_qc_adv(42, 42,DE_HT_FACTOR*(code.num_qubits+num_ancillae)*(code.num_qubits+num_ancillae),
		STICK_HT_FACTOR*distance*distance, error_probability, t_delete, rec, dir);
	//create qubits
	for (int j = 0; j < code.num_qubits+num_ancillae; ++j) {
		qubit_array[j] = qc_create_and_insert_qubit(dp_qc->qc, j, 0, 0, QUBIT_HT_SIZE);
	}

	boundaries = (BALL **)my_malloc((code.X_boundary.size() + code.Z_boundary.size() + 4) * sizeof(BALL *));

	// first primal
	for(int i = 0; i < code.X_boundary.size(); ++i){
		primal_boundary[i] = create_boundary_set(dp_qc->qc, i, PRIMAL_BOUNDARY, -i, -1, 0);
	}
	// second dual
	for(int i = 0; i < code.Z_boundary.size(); ++i){
		dual_boundary[i] = create_boundary_set(dp_qc->qc, i + code.X_boundary.size(), DUAL_BOUNDARY, - i - code.X_boundary.size(), -2, 0);
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
	if(frame)
		free(frame);
	if(qubit_array);
		free(qubit_array);
	if(primal_boundary)
		free(primal_boundary);
	if(dual_boundary)
		free(dual_boundary);
	if(primal_syndrome)
		free(primal_syndrome);
	if(dual_syndrome)
		free(dual_syndrome);
	if(temporal_primal_boundary)
		free(temporal_primal_boundary);
	if(temporal_dual_boundary)
		free(temporal_dual_boundary);
	if(boundaries)
		free(boundaries);
	if(primal_set)
		free(primal_set);
	if(dual_set)
		free(dual_set);

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
}

Wrapper::~Wrapper(){
	free(dir);
	delete_data();
	qc_free_recipe_adv(recipe);
	dp_free_dp_qc(dp_qc);
}



//simulation related

void Wrapper::generate_recipe(){
	size_type big_t = 0;
	do{
		measure_stabilizers(big_t++);
	}while (qc_boot_up_adv(dp_qc->qc, recipe, 2, 12) != DONE);

	return;
}


void Wrapper::run_simulation(size_type num_runs){

}


void Wrapper::calculate_t_check(double error_probability){
	allocate();
	init(recipe,error_probability);
	dp_qc->qc->track = false;
	dp_qc->qc->get_boundary = get_boundary;
	
	size_type big_t = 0;
	size_type t_check = 1;
	size_type t_check_scale = 10;
	size_type boot_num_X = 50;
	size_type boot_num_Z = 50;

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
				std::cout << "Calculating tcheck: " << t_check << " "  << num_checks << " " << changes << " " << t_check_scale;

				t_check = t_check * num_checks / changes / t_check_scale;

				if (t_check < 1) { 
					t_check = 1;
				}
				
				break;
			}
		}

		// If we are verbose, then print out the current SC_DP_QC statistics
		#ifdef DEBUG
			print_stats();
		#endif

		// Ensure that the qc has been properly advanced in big_t
		assert(dp_qc->qc->big_t == big_t+1);
		
		big_t++;
	}
	// Free the sc_dp_qc now that the t_check has been calculated
	delete_data();
}

void Wrapper::reset_recipe(){

}


void Wrapper::test_correct() {
	int count;
	long big_t = 0;
	//printf("Test correct\n");
	Wrapper second(*this);


	// Enable perfect gates and perform a perfect round of stabilizer
	// measurements, followed by matching, then correction.
	second.dp_qc->qc->perfect_gates = true;

	second.dp_qc->qc->m_pr->undo_flag = true;
	second.dp_qc->qc->m_du->undo_flag = true;

	second.measure_stabilizers(big_t);

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
	//printf("End test correct\n");
}


void Wrapper::print_stats() {
	std::cout << "TODO" << std::endl;
}


void Wrapper::process_aug_edge(AUG_EDGE *ae) {
	// Get the coordinates at the start and end of the edge (come from boundary coordinates and set coordinates)
	// ae->va: source vertex of the edge
	// ae->vb: destination vertex of the edge
	int i1 = ae->va->i;
	int j1 = ae->va->j;
	int	i2 = ae->vb->i;
	int	j2 = ae->vb->j;
	int error = 0;
	assert(!(j1<0 && j2<0));
	// first resolve boundary to proper coordinates
	if(j1 < 0){
		// boundary
		if(j1 <= -3){ // temporal boundary
			return;
		}
		if(j1 == -1){ // primal boundary: set to closest qubit position
			i1 =
			j1 =
			error = Z; // primal has X-stabilizers -> detection of Z errors
		}
		if(j1 == -2){ // dual boundary: set to closest qubit position
			i1 = 
			j1 = 
			error = X; // dual has Z-stabilizers -> detection of X errors
		}
	}
	// same for the destination
	if(j2 < 0){
		// boundary
		if(j2 <= -3){ // temporal boundary
			return;
		}
		if(j2 == -1){ // primal boundary: set to closest qubit position
			i2 =
			j2 =
			error = Z; // primal has X-stabilizers -> detection of Z errors
		}
		if(j2 == -2){ // dual boundary: set to closest qubit position
			i2 = 
			j2 = 
			error = X; // dual has Z-stabilizers -> detection of X errors
		}
	}


	// perform pathfinding
	ParityCheck path;
	// from the path deduce the data qubits

	// apply an error all elements of the frame along the path
	for(auto id : path){
		frame[id] ^= error;
	}

	return;
}

/**
 * \brief Process augmented edges from a \ref matching
 * 
 * \param[in] m The \ref matching containing the augmented edges
 * \param[in] d The distance of the matching
 * \param[in] frame The Pauli frame to process the edges on
 */
void Wrapper::process_aug_edges(MATCHING *m) {
	AUG_EDGE *ae;

	ae = m_get_aug_edge(m);
	while (ae != NULL) {
		process_aug_edge(ae);
		m_delete_aug_edge(ae);
		ae = m_get_aug_edge(m);
	}
}
