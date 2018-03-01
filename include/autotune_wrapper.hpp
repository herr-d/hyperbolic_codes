#ifndef AUTOTUNE_WRAPPER
#define AUTOTUNE_WRAPPER

#include <basic_includes.hpp>
#include <depolar/depolar.h>
#include <memory/memory.h>
#include <my_time/my_time.h>



/*
 * Wrapper class for the autotune library.
 * 
 *
*/
class Wrapper{
private:
	//these will be replaced
	int seed0 = 42;
	int seed1 = 43;
	int distance = 10;
	int t_delete = 100; //delete after 100 steps
	char dir[];


	///////////////////////////////////////////////////////////////////////////////////////////
	// Varibles needed for Hyperbolic codes
	///////////////////////////////////////////////////////////////////////////////////////////

	const StabilizerContainer& X_stabilizers, Z_stabilizers; // Stabilizers of the code
	const StabilizerContainer& X_boundary, Z_boundary; // List of stabilizers on the boundary
	DP_QC *dp_qc; //depolarizing quantum computer
	size_type num_data_qubits; // number of data qubits of the code
	size_type num_ancillae; // how many syndrome qubits are needed


	///////////////////////////////////////////////////////////////////////////////////////////
	// autotune specific variables
	///////////////////////////////////////////////////////////////////////////////////////////
	RECIPE_ADV *recipe; // recipe... whatever that is

	int* frame; // frame: stores the changes from the syndrome measurements
	
	QUBIT** qubit_array; // qubits needed for the actual simulation

	SYNDROME** primal_syndrome; // primal syndomes used in the actual simulation
	SYNDROME** dual_syndrome; // dual syndromes used in the actual simulation
	SET*** primal_set; // sets are associated with each primal syndrome on each time step
	SET*** dual_set; // sets are associated with each dual syndrome on each time step

	BALL** boundaries; // ball of boundaries defined above used for the simulation
	SET** spacial_boundary; // i don't know
	SET** temporal_boundary; // i don't know




	///////////////////////////////////////////////////////////////////////////////////////////
	// Functions adapted from ex1 of autotune
	///////////////////////////////////////////////////////////////////////////////////////////

	/*
	 * Creates the boundary sets for the code and stores them in spacial_ or temporal_boundary
	 * [in]:
	 *
	 * [out]: pointer to created SET
	*/
	SET *create_boundary_set(QC *qc, int id, int type, int i, int j, int t);

	/*
	 * Generates the sets and syndromes needed for the simulation
	 * [in]: no argument, but relies on defined X_stabilzers, etc.
	 * [out]: no return, but primal_set, primal_syndrome, etc. are initialized
	*/	
	void create_initial_syndrome_and_set_arrays();

	/*
	 * Performes all syndrome checks for the entire code once
	 * [in]: no argument, but everything needs to be initialized
	 *
	 * [out]: adds results to frame such that autotune can process it
	*/
	void syndrome_checks();

	/*
	 * performs the Z_stabilizer defined in Z_stabilizers.at(pos)
	 * [in]:
	 * 		pos: indicates which stabilizer of the array should be performed
	 * [out]: no return, but result is added to frame
	*/
	void perform_Z_stabilizer(size_type pos);

	/*
	 * performs the X_stabilizer defined in X_stabilizers.at(pos)
	 * [in]:
	 * 		pos: indicates which stabilizer of the array should be performed
	 * [out]: no return, but result is added to frame
	*/
	void perform_X_stabilizer(size_type pos);

	/*
	 * During initialization of syndrome qubits, data qubits are idle and subject to errors
	 * [in]:
	 *
	 * [out]: no return, errors are tracked by autotune
	*/
	void add_idle_data_qubit_errors();

	/*
	 * Given a stabilizer, is it on the boundary?
	 * [in]:
	 * 		pos: indicates which stabilizer of the array
	 * 		X_not_Z: flag to switch between X and Z stabilizers
	 * [out]: pointer to boundary SET, if not on boundary return NULL
	*/
	SET *on_boundary(size_type pos, bool X_not_Z);

public:
	/*
	 * Initialization of everything
	 * [in]:
	 * 		X_stabil: const reference to all X stabilizers
	 * 		Z_stabil: const reference to all Z stabilizers
 	 * 		X_bound: const reference to all X boundaries
	 * 		Z_bound: const reference to all X boundaries
	 * 		qubits: number of data qubits in the code
	 * 		error_probability: sets a constant probability for errors to occur
	 				(the weight of different errors can be adjusted in the EMS directory)
	 * [out]: no return, but all memory is allocated and all variables are initialized
	*/
	Wrapper(const Code_info & code, double error_probability);

	/*
	 * Runs the simulation, to get an estimate on the logical error rate
	 * [in]:
	 * 		num_runs: number of times the simulation runs
	 * [out]: no return, but average on logical error rate is written as standard output
	*/
	void run_simulation(size_type num_runs);


	/*
	 * Savely frees all memory associated with the simulation
	 * [in]: no input
	 * [out]: no return
	*/
	~Wrapper();
};


#endif
//AUTOTUNE_WRAPPER