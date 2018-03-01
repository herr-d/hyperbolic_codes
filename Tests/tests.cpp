
#include <gtest/gtest.h>
#include <generate_code.hpp>

//test the creation of the hyperbolic lattice
//#include "tests_generate_code.cpp"




//Tests

// tests for the creation of the code
#include "tests_generate_code.hpp"

// tests for the error simulation
#include "tests_simulate_errors.hpp"


int main(int argc, char* argv[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
