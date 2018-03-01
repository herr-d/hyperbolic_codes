#include <autotune_wrapper.hpp>
#include <string.h>

// This rearranges example 1 such that it can be used with hyperbolic codes



SET* Wrapper::create_boundary_set(QC *qc, int id, int type, int i, int j, int t) {
	SET *set;

	set = qc_create_set_adv(qc, type, i, j, t, NULL);
	boundaries[id] = set->ball;

	return set;
}


SET* Wrapper::on_boundary(size_type pos, bool X_not_Z){
	return NULL;
}


void Wrapper::create_initial_syndrome_and_set_arrays(){
	for(int i = 0; i < X_stabilizers.size(); ++i){
		primal_syndrome[i] = qc_create_syndrome();
		qc_insert_syndrome(dp_qc->qc, primal_syndrome[i]);
		if(true) // TODO check if on boundary
		{
			primal_set[i][1] = qc_create_set_adv(dp_qc->qc, PRIMAL, 0, 0, 1, spacial_boundary[0]);
			primal_set[i][0] = qc_create_set_adv(dp_qc->qc, PRIMAL, 0, 0, 2, spacial_boundary[0]);
		}
		else{
			primal_set[i][1] = qc_create_set_adv(dp_qc->qc, PRIMAL, 0, 0, 1, NULL);
			primal_set[i][0] = qc_create_set_adv(dp_qc->qc, PRIMAL, 0, 0, 2, NULL);
		}

		qc_insert_set(dp_qc->qc, primal_set[i][1]);
		qc_associate_syndrome(primal_set[i][1], primal_syndrome[i]);
		qc_insert_set(dp_qc->qc, primal_set[i][0]);
	}
	for(int i = 0; i < Z_stabilizers.size(); ++i){
		dual_syndrome[i] = qc_create_syndrome();
		qc_insert_syndrome(dp_qc->qc, dual_syndrome[i]);
		if(on_boundary(i,false)!=NULL) // TODO check if on boundary
		{
			dual_set[i][1] = qc_create_set_adv(dp_qc->qc, DUAL, 0, 0, 1, on_boundary(i,false));
			dual_set[i][0] = qc_create_set_adv(dp_qc->qc, DUAL, 0, 0, 2, on_boundary(i,false));
		}
		else{
			dual_set[i][1] = qc_create_set_adv(dp_qc->qc, DUAL, 0, 0, 1, NULL);
			dual_set[i][0] = qc_create_set_adv(dp_qc->qc, DUAL, 0, 0, 2, NULL);
		}
		qc_insert_set(dp_qc->qc, dual_set[i][1]);
		qc_associate_syndrome(dual_set[i][1], dual_syndrome[i]);
		qc_insert_set(dp_qc->qc, dual_set[i][0]);
	}
}


void Wrapper::perform_X_stabilizer( size_type pos){
	// init
	size_type ancilla_qubit = num_data_qubits + pos;
	dp_init_Z(dp_qc,qubit_array[ancilla_qubit]);
	dp_H(dp_qc,qubit_array[ancilla_qubit]);

	// for worst case each ancilla qubit needs to wait 2 steps
	// each qubit has 2 X and 2 Z-stabilizers, due to the hadamard operation they are shifted by 1
	// thus 1 step definitely works, for a clash one of the choices also works thus the worst case has a delay of 2
	for(int i=0; i < 2; ++i)
		dp_iden_H(dp_qc,qubit_array[ancilla_qubit]);
	
	// CNOTS
	for(auto qubit : X_stabilizers.at(pos)){
		dp_cnot(dp_qc, qubit_array[qubit], qubit_array[ancilla_qubit]);
	}

	// measurement
	dp_H(dp_qc,qubit_array[ancilla_qubit]);
}


void Wrapper::perform_Z_stabilizer(size_type pos){
	size_type lay1, lay2;

	// init
	size_type ancilla_qubit = num_data_qubits + pos; 
	dp_init_Z(dp_qc,qubit_array[ancilla_qubit]);

	// for worst case each ancilla qubit needs to wait 2 steps
	// each qubit has 2 X and 2 Z-stabilizers, due to the hadamard operation they are shifted by 1
	// thus 1 step definitely works' For a clash, one of the choices also works thus the worst case has a delay of 2
	for(int i=0; i < 2; ++i)
		dp_iden_init_Z(dp_qc, qubit_array[ancilla_qubit]);
	
	// CNOTS
	for(auto qubit : X_stabilizers.at(pos)){
		dp_cnot(dp_qc, qubit_array[ancilla_qubit], qubit_array[qubit]);
	}

	// measurement
	// TODO: change to here
	dp_meas_Z(dp_qc, qubit_array[ancilla_qubit], primal_set[pos][lay1], primal_set[pos][lay2]);
	qc_unassociate_syndrome(primal_syndrome[pos]);
	qc_associate_syndrome(primal_set[pos][lay2], primal_syndrome[pos]);
	//create new sets and check for boundaries
	if (on_boundary(pos,true)!=NULL) {
		primal_set[pos][lay1] = qc_create_set_adv(dp_qc->qc, PRIMAL, 0, 0, 2, on_boundary(pos,true));
	}
	else {
		primal_set[pos][lay1] = qc_create_set_adv(dp_qc->qc, PRIMAL, 0, 0, 2, NULL);
	}
	qc_insert_set(dp_qc->qc, primal_set[pos][lay1]);
}

void Wrapper::add_idle_data_qubit_errors(){
	// due to initialization and Hadamards on syndrome qubits
	// identity operations need to be applied
	for(int i = 0; i < num_data_qubits; ++i){
		dp_iden_init_Z(dp_qc, qubit_array[i]);
		dp_iden_H(dp_qc, qubit_array[i]);
		dp_iden_H(dp_qc, qubit_array[i]);
		dp_iden_meas_Z(dp_qc, qubit_array[i]);

	}
}



void Wrapper::syndrome_checks(){
	// loop through syndromes
	// unfortunately lattice is fustrated in general, so assume worst case
	// where all syndrome block each other

	add_idle_data_qubit_errors();

	for(int i =0; i < X_stabilizers.size(); ++i){
		perform_X_stabilizer(i);
	}
	for(int i =0; i < Z_stabilizers.size(); ++i){
		perform_Z_stabilizer(i);
	}

	return;
}

	

Wrapper::Wrapper(const Code_info & code, double error_probability)
	 				: X_stabilizers(code.X_stabilizer), Z_stabilizers(code.Z_stabilizer), //set const references to the generated graph
	 				X_boundary(code.X_boundary), Z_boundary(code.Z_boundary), // references to the boundary stabilizers
	 				num_data_qubits(code.num_qubits){

	strcpy(dir, "../../autotune/ex/ems/");

	//allocate memory
	frame = (int *)my_malloc((num_data_qubits + num_ancillae) * sizeof(int));
	qubit_array = (QUBIT**) my_calloc(num_data_qubits+num_ancillae, sizeof(QUBIT *));
	spacial_boundary = (SET**) my_calloc(2*code.X_boundary.size(), sizeof(SET *));
	temporal_boundary = (SET**) my_calloc(2*code.Z_boundary.size(), sizeof(SET *));

	primal_syndrome = (SYNDROME**)my_calloc(code.Z_stabilizer.size(), sizeof(SYNDROME *));
	dual_syndrome = (SYNDROME**)my_calloc(code.X_stabilizer.size(), sizeof(SYNDROME *));

	primal_set = (SET ***)my_2d_calloc(code.Z_stabilizer.size(), 2, sizeof(SET *));
	dual_set = (SET ***)my_2d_calloc(code.X_stabilizer.size(), 2, sizeof(SET *));


	//create recipe-- whatever that does
	recipe = qc_create_recipe_adv(RECIPE_INFINITE, 100, 100, FALSE);

	//create dp_qc
	dp_qc = dp_create_dp_qc_adv(seed0, seed1,100*(num_data_qubits+num_ancillae), 200*distance*distance, error_probability, t_delete, recipe, dir);

	//create qubits
	for (int j = 0; j < num_data_qubits+num_ancillae; j++) {
		qubit_array[j] = qc_create_and_insert_qubit(dp_qc->qc, 0, 0, 0, 1000);
	}

	// first primal
	for(int i = 0; i < code.X_boundary.size()*2; ++i){
		spacial_boundary[i] = create_boundary_set(dp_qc->qc, i, PRIMAL_BOUNDARY, -1, 0, 0);
	}
	// second dual
	for(int i = 0; i < code.Z_boundary.size()*2; ++i){
		temporal_boundary[i] = create_boundary_set(dp_qc->qc, i, DUAL_BOUNDARY, -1, 0, 0);
	}
}


Wrapper::~Wrapper(){
	free(frame);
	free(qubit_array);
	free(spacial_boundary);
	free(temporal_boundary);
	free(primal_syndrome);
	free(dual_syndrome);

	free(primal_set);
	free(dual_set);
	qc_free_recipe_adv(recipe);
	dp_free_dp_qc(dp_qc);
}

void Wrapper::run_simulation(size_type num_runs){
	return;
}
