#include <gtest/gtest.h>
#include <generate_code.hpp>


TEST(Hyperbolic, Creation_Surface){
	//surface code
	Hyperbolic small(4, 4, 1);
	small.build_lattice();
	small.printout_graph();
	EXPECT_EQ(9,small.get_number_Z_stabilizers());
	EXPECT_EQ(16,small.get_number_X_stabilizers());
	EXPECT_EQ(24,small.get_number_qubits());
	//make rough edges

/*	small.add_rough_edges(2);
	EXPECT_EQ(9,small.get_number_Z_stabilizers());
	EXPECT_EQ(8,small.get_number_X_stabilizers());
	EXPECT_EQ(18,small.get_number_qubits());
	*/
}


TEST(Hyperbolic, Creation_Hyperbolic){
	//hyperbolic code
	Hyperbolic small(5, 4, 1);
	small.build_lattice();
	EXPECT_EQ(11,small.get_number_Z_stabilizers());
	EXPECT_EQ(30,small.get_number_X_stabilizers());
	EXPECT_EQ(40,small.get_number_qubits());
	//make rough edges
	//needs to be thought over again
	/*
	small.add_rough_edges(2);
	EXPECT_EQ(0,small.get_number_Z_stabilizers());
	EXPECT_EQ(0,small.get_number_X_stabilizers());
	EXPECT_EQ(0,small.get_number_qubits());*/
}

TEST(Hyperbolic, new_layer){
	Hyperbolic test(4,4,1);
	test.generate_center_face();
	test.new_layer();
	EXPECT_EQ(9,test.get_number_Z_stabilizers());
	EXPECT_EQ(16,test.get_number_X_stabilizers());
	EXPECT_EQ(24,test.get_number_qubits());
}


TEST(Hyperbolic, generate_center_face){
	Hyperbolic test(8,8,2);
	test.generate_center_face();
	//EXPECT_EQ(1,test.get_number_Z_stabilizers());
	EXPECT_EQ(8,test.get_number_X_stabilizers());
	EXPECT_EQ(8,test.get_number_qubits());
}

TEST(Hyperbolic, generate_center_vertex){
	Hyperbolic test(4,4,2);
	test.generate_center_vertex();
	EXPECT_EQ(4,test.get_number_Z_stabilizers());
	EXPECT_EQ(9,test.get_number_X_stabilizers());
	EXPECT_EQ(12,test.get_number_qubits());
}


TEST(Hyperbolic, generate_rough_edges){

}


int main(int argc, char* argv[]){
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
