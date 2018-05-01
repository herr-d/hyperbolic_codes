
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

add_executable(test_generate tests_generate_code.cpp)
add_executable(test_simulate tests_simulate_errors.cpp)

target_link_libraries(test_generate ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} Generate_Code Autotune_Wrapper Path_Finding autotune ${BLOSSOM_PM} ${BLOSSOM_MISC} ${BLOSSOM_MIN_COST} ${BLOSSOM_GEOM})
target_link_libraries(test_simulate ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} Generate_Code Autotune_Wrapper Path_Finding autotune ${BLOSSOM_PM} ${BLOSSOM_MISC} ${BLOSSOM_MIN_COST} ${BLOSSOM_GEOM})


gtest_add_tests(TARGET	test_generate
				TEST_SUFFIX	.noArgs)
gtest_add_tests(TARGET	test_simulate
				TEST_SUFFIX	.noArgs)