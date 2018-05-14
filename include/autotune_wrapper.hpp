#ifndef AUTOTUNE_WRAPPER
#define AUTOTUNE_WRAPPER

#ifdef GTEST
#define private public
#define protected public
#endif



#include <basic_includes.hpp>
#include <depolar/depolar.h>
#include <memory/memory.h>
#include <my_time/my_time.h>


//stupid function pointer... I cant define it in the class and thus need access to variables... thus: extern Code_info code
BALL *get_boundary(int i, int j, long int big_t, int type, void *boundaries);
extern Code_info code;

/*
 * Wrapper class for the autotune library.
 * This wrapper is a modification of ex1 from the Autotune library.
 * It should always work, as long as the Code_info container code is properly initialized.
 * Any changes throughout this wrapper class should be unnecessary.
 *
 * The user should only need to execute the following public functions to be able to run the threshold calculation:
 * Constructor: Wrapper simulation(code.probability);
 * generate_recipe();
 * calculate_t_check(); //optional
 * run_simulation(code.big_t_max);
 *
*/
class Wrapper{
private:
	bool copy; // is it a copy of anothec element?

	const int num_padding = 2; //total number of steps that need to be added due to scheduling conflicts: worst case estimate;
	int distance;
	int num_checks;
	std::vector<size_type> num_X_changes;
	std::vector<size_type> num_Z_changes;
	std::vector<size_type> last_X_check;
	std::vector<size_type> last_Z_check;
	double probab;

	///////////////////////////////////////////////////////////////////////////////////////////
	// autotune specific variables
	///////////////////////////////////////////////////////////////////////////////////////////
	DP_QC *dp_qc; //depolarizing quantum computer
	RECIPE_ADV *recipe; // recipe... whatever that is

	int* frame; // frame: stores the changes from the syndrome measurements
	
	QUBIT** qubit_array; // qubits needed for the actual simulation

	SYNDROME** primal_syndrome;
	SYNDROME** dual_syndrome;
	SET*** primal_set; // sets are associated with each primal syndrome + has two time steps
	SET*** dual_set; // sets are associated with each dual syndrome + has two time steps

	BALL** boundaries; // ball of boundaries defined above used for the simulation
	SET** primal_boundary; // primal sets associated with the boundaries
	SET** dual_boundary; // dual sets associated with the boundaries
	SET** temporal_primal_boundary; // sets that indicate a primal boundary in time
	SET** temporal_dual_boundary; // sets that indicate a dual boundary in time direction



	///////////////////////////////////////////////////////////////////////////////////////////
	// Functions adapted from ex1 of autotune
	///////////////////////////////////////////////////////////////////////////////////////////

	/*
	 * Creates the boundary sets for the code and stores them in spacial_ or temporal_boundary
	 * [in]:
	 *		id:		unique id of the boundary
	 *		type:	type of the boundary PRIMAL_BOUNDARY, DUAL_BOUNDARY
	 *		i,j,t:	coordinates of the boundary
	 * [out]: pointer to created SET
	*/
	SET *create_boundary_set(QC *qc, int id, int type, int i, int j, int t);


	/* Applies corrections due to the syndrome measurements by processing
	 *	augmented edges from primal and dual matchings 
	 * [in]:
	 * [out]:
	*/
	void correct_mts();

	/*
	 * Creates the boundary sets for the code and stores them in spacial_ or temporal_boundary
	 * [in]:
	 *		ae: The augmented edge to process
	 * [out]: pointer to created SET
	*/
	void process_aug_edge(AUG_EDGE *ae);

	/*
	 * Creates the boundary sets for the code and stores them in spacial_ or temporal_boundary
	 * [in]:
	 *
	 * [out]: pointer to created SET
	*/
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


	/*
	 * Initializes the ancilla qubits for the X stabilizers
	 * [in]: no argument, information is provided in the code object
	 * [out]:
	*/
	void init_X_stabilizers();

	/*
	 * Initializes the ancilla qubits for the Z stabilizers
	 * [in]: no argument, information is provided in the code object
	 * [out]:
	*/
	void init_Z_stabilizers();

	/*
	 * Adds errors to idle data qubits during initialization of syndrome qubits
	 * [in]: no argument, information is provided in the code object
	 * [out]:
	*/
	void data_error_during_init();

	/*
	 * This function cleans up all objects created except for the recipe
	 * [in]:
	 * [out]:
	*/
	void delete_data();


	/*
	 * performs all Z_stabilizers defined in the code object
	 * [in]:
	 * [out]:
	*/
	void perform_Z_stabilizers();

	/*
	 * performs all X stabilizers defined in the code object
	 * [in]:
	 * [out]:
	*/
	void perform_X_stabilizers();


	void apply_cnots();
	
	/*
	 * Adds the measurement of all X stabilizers of this time step.
	 * [in]:
	 * 		pos:
	 * [out]: no return, but result are stored internally
	*/
	void measure_X_stabilizers(size_type big_t);

	/*
	 * Adds the measurement of all Z stabilizers of this time step.
	 * [in]:
	 * 		pos:
	 * [out]: no return, but result are stored internally
	*/
	void measure_Z_stabilizers(size_type big_t);


	/*
	 * Memory for the simulation is allocated.
	 * [in]:
	 * [out]: no return, but result are stored internally
	*/
	void allocate();

	/*
	 * All allocated objects are initialized. The simulation can now be run.
	 * [in]:
	 * 		pos:
	 * [out]: no return, but result are stored internally
	*/
	void init(RECIPE_ADV *rec);

	/*
	 * Given a stabilizer, is it on the boundary?
	 * [in]:
	 * 		pos: indicates which stabilizer of the array
	 * 		X_not_Z: flag to switch between X and Z stabilizers
	 * [out]: pointer to boundary SET, if not on boundary return NULL
	*/
	SET *infer_boundary(size_type pos, bool X_not_Z);

	/*
	 * Debbuging tool, that outputs the statistics of the current run.
	 * [in]:
	 * [out]: pointer to boundary SET, if not on boundary return NULL
	*/
	void print_stats();

	/*
	 * Debbuging tool, that outputs the statistics of the current run.
	 * [in]:
	 * [out]: pointer to boundary SET, if not on boundary return NULL
	*/
	void print_frame();
	void print_qubit_array();


	/*
	 * Test the attempted correction. If it is incorrect a logical error happens.
	 * [in]:
	 * [out]:
	*/
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
	 * A recipe is generated by performing several rounds of error correction.
	 * Internally autotune detects the proper cycle length and sets up internal variables
	 * [in]:
	 * [out]:
	*/
	void generate_recipe();


	/*
	 * How often should a t check be performed/Detection of error events.
	 * This adjusts the provided t_check value to improve the execution time of the threshold calculation
	 * [in]:
	 * [out]:
	*/
	void calculate_t_check();

	/*
	 * Resets the recipe. This should not be needed, except if two completely different simulations are run with the same instance.
	 * Afterwards a recipe needs to be generated again.
	 * [in]:
	 * [out]:
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