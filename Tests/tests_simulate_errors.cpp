#include <gtest/gtest.h>

#define private public


#include <generate_code.hpp>
#include <path_finding.hpp>
#include <autotune_wrapper.hpp>
#include <path_finding.hpp>
#include <gen_surface_code.hpp>


Code_info code;


TEST(autotune_wrapper, infer_boundary){

	Wrapper simulation(code.probability);

	SET* b_p0 = simulation.primal_boundary[0];
	SET* b_p1 = simulation.primal_boundary[1];
	SET* b_d0 = simulation.dual_boundary[0];
	SET* b_d1 = simulation.dual_boundary[1];

	//Data qubits
	std::vector<SET*> results_primal = {b_p0, NULL, b_p1, b_p0, NULL, b_p1, b_p0, NULL, b_p1, b_p0, NULL, b_p1};
	std::vector<SET*> results_dual = {b_d0, b_d0, b_d0, b_d0, NULL, NULL, NULL, NULL, b_d1, b_d1, b_d1, b_d1};
	for(int i=0; i < 12; ++i){
		EXPECT_EQ(results_primal.at(i), simulation.infer_boundary(i,false));
		EXPECT_EQ(results_dual.at(i), simulation.infer_boundary(i,true));
	}
}


TEST(autotune_wrapper, get_boundary){

	Wrapper simulation(code.probability);
	
	//Data qubits
	for(int i = 0; i <= code.num_qubits; ++i){
		EXPECT_EQ(NULL, get_boundary(i, 0, LONG_MAX, PRIMAL_BOUNDARY, simulation.boundaries));
		EXPECT_EQ(NULL, get_boundary(i, 0, LONG_MAX, DUAL_BOUNDARY, simulation.boundaries));
	}

	//boundaries
	EXPECT_EQ(simulation.boundaries[0], get_boundary(0, -1, LONG_MAX, PRIMAL_BOUNDARY, simulation.boundaries));
	EXPECT_EQ(simulation.boundaries[1], get_boundary(1, -1, LONG_MAX, PRIMAL_BOUNDARY, simulation.boundaries));
	EXPECT_EQ(simulation.boundaries[2], get_boundary(0, -2, LONG_MAX, DUAL_BOUNDARY, simulation.boundaries));
	EXPECT_EQ(simulation.boundaries[3], get_boundary(1, -2, LONG_MAX, DUAL_BOUNDARY, simulation.boundaries));
	EXPECT_EQ(simulation.boundaries[4], get_boundary(-1, -3, LONG_MAX, PRIMAL_BOUNDARY, simulation.boundaries));
	EXPECT_EQ(simulation.boundaries[5], get_boundary(-2, -3, LONG_MAX, PRIMAL_BOUNDARY, simulation.boundaries));
	EXPECT_EQ(simulation.boundaries[6], get_boundary(-1, -4, LONG_MAX, DUAL_BOUNDARY, simulation.boundaries));
	EXPECT_EQ(simulation.boundaries[7], get_boundary(-2, -4, LONG_MAX, DUAL_BOUNDARY, simulation.boundaries));
}


int main(int argc, char* argv[])
{
	code.dir = (char*)malloc(100*sizeof(char));
	code.boot = false;
	code.r = 4;  // distance
	code.s = 3;
	code.l = 3;
	code.t_check = 1;
	code.t_delete = 100;
	code.big_t_max= 2000;
	code.probability = 0.0001;
	code.seed0 = 42;
	code.seed1 = 42;
	code.max_num_X = 2000;
	code.max_num_Z = 2000;
	generate_surface_code();


	strcpy(code.dir, EMS_DIR);


	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
