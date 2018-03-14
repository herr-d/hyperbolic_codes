#ifndef AUTOTUNE_WRAPPER
#define AUTOTUNE_WRAPPER

#include <basic_includes.hpp>
#include <depolar/depolar.h>
#include <memory/memory.h>
#include <my_time/my_time.h>


//stupid function pointer... I cant define it in the class and thus need access to variables... thus: extern Code_info code
BALL *get_boundary(int i, int j, long int big_t, int type, void *boundaries);
extern Code_info code;

/*
 * Wrapper class for the autotune library.
 * 
 *
*/
class Wrapper{
private:
	char *dir;

	const int num_padding = 2; //total number of steps that need to be added due to scheduling conflicts: worst case estimate;
	int seed0;
	int seed1;
	int distance;
	int t_delete; //delete after 100 steps
	int num_checks;
	int num_X_changes;
	int num_Z_changes;
	int last_X_check;
	int last_Z_check;

	///////////////////////////////////////////////////////////////////////////////////////////
	// Varibles needed for Hyperbolic codes
	///////////////////////////////////////////////////////////////////////////////////////////
	DP_QC *dp_qc; //depolarizing quantum computer


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
	SET** primal_boundary; // i don't know
	SET** dual_boundary; // i don't know
	SET** temporal_primal_boundary;
	SET** temporal_dual_boundary;



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

	void correct_mts();

	void process_aug_edge(AUG_EDGE *ae);
	void process_aug_edges(MATCHING *m);


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
	void measure_stabilizers(size_type big_t);


	void init_X_stabilizers();
	void init_Z_stabilizers();
	void data_error_during_init();
	void delete_data();


	/*
	 * performs the Z_stabilizer defined in Z_stabilizers.at(pos)
	 * [in]:
	 * 		pos: indicates which stabilizer of the array should be performed
	 * [out]: no return, but result is added to frame
	*/
	void perform_Z_stabilizers();

	/*
	 * performs the X_stabilizer defined in X_stabilizers.at(pos)
	 * [in]:
	 * 		pos: indicates which stabilizer of the array should be performed
	 * [out]: no return, but result is added to frame
	*/
	void perform_X_stabilizers();

	void measure_X_stabilizers(size_type big_t);
	void measure_Z_stabilizers(size_type big_t);
	void allocate();
	void init(RECIPE_ADV *rec, double error_probability);

	/*
	 * Given a stabilizer, is it on the boundary?
	 * [in]:
	 * 		pos: indicates which stabilizer of the array
	 * 		X_not_Z: flag to switch between X and Z stabilizers
	 * [out]: pointer to boundary SET, if not on boundary return NULL
	*/
	SET *infer_boundary(size_type pos, bool X_not_Z);

	void print_stats();

	void test_correct();

public:
	/*
	 * Initialization of everything
	 * [in]:
	 * 		
	 * [out]: no return, but all memory is allocated and all variables are initialized
	*/
	Wrapper(double error_probability);
	Wrapper();
	Wrapper(const Wrapper & a);


	/*
	 *
	 *
	 *
	*/
	void generate_recipe();


	/*
	 *
	 *
	 *
	*/
	void calculate_t_check(double error_probability);

	/*
	 *
	 *
	 *
	*/
	void reset_recipe();


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