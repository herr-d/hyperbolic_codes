
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

add_executable(test_generate tests_generate_code.cpp)

target_link_libraries(test_generate ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} Generate_Code Autotune_Wrapper Path_Finding)

gtest_add_tests(TARGET      test_generate
                TEST_SUFFIX .noArgs
                TEST_LIST   noArgsTests
)

set_tests_properties(${noArgsTests}   PROPERTIES TIMEOUT 10)
