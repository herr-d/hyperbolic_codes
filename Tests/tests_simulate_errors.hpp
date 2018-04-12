#include <gtest/gtest.h>

#include <autotune_wrapper.hpp>
#include <path_finding.h>
#include <generate_code.hpp>

Code_info code;

TEST(Pathfinding, DataQubits){
	code.boot = false;
	strcpy(code.dir, "../../autotune/ex/ems/");
	code.r = 3;  // distance
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

	// Test pathfinding algorithm

	// NO NEED SO FAR... THE ERROR CORRECTION ROUTINE IS NOT YET CALLED	
}

TEST(autotune_wrapper, infer_boundaries){
	code.boot = false;
	strcpy(code.dir, "../../autotune/ex/ems/");
	code.r = 3;  // distance
	code.s = 3;
	code.l = 1;
	code.t_check = 1;
	code.t_delete = 100;
	code.big_t_max= 2000;
	code.probability = 0.0001;
	code.seed0 = 42;
	code.seed1 = 42;
	code.max_num_X = 2000;
	code.max_num_Z = 2000;

	generate_surface_code();

	Wrapper simulation(code.probability);

	simulation.generate_recipe();

	// Test pathfinding algorithm
	EXPECT_EQ(simulation.primal_boundary[0],simulation.infer_boundary(0,false));
	EXPECT_EQ(NULL,simulation.infer_boundary(1,false));
	EXPECT_EQ(simulation.primal_boundary[1],simulation.infer_boundary(2,false));
	EXPECT_EQ(simulation.primal_boundary[0],simulation.infer_boundary(3,false));
	EXPECT_EQ(NULL,simulation.infer_boundary(4,false));
	EXPECT_EQ(simulation.primal_boundary[1],simulation.infer_boundary(5,false));

	EXPECT_EQ(simulation.dual_boundary[0],simulation.infer_boundary(0,true));
	EXPECT_EQ(simulation.dual_boundary[0],simulation.infer_boundary(1,true));
	EXPECT_EQ(NULL,simulation.infer_boundary(2,true));
	EXPECT_EQ(NULL,simulation.infer_boundary(3,true));
	EXPECT_EQ(simulation.dual_boundary[1],simulation.infer_boundary(4,true));
	EXPECT_EQ(simulation.dual_boundary[2],simulation.infer_boundary(5,true));
}


TEST(autotune_wrapper, measure_stabilizers){
	code.boot = false;
	strcpy(code.dir, "../../autotune/ex/ems/");
	code.r = 3;  // distance
	code.s = 3;
	code.l = 1;
	code.t_check = 1;
	code.t_delete = 100;
	code.big_t_max= 2000;
	code.probability = 0.0001;
	code.seed0 = 42;
	code.seed1 = 42;
	code.max_num_X = 2000;
	code.max_num_Z = 2000;

	generate_surface_code();

	Wrapper simulation(code.probability);

	simulation.generate_recipe();

	//Test the logic

	

}
