
add_library(Generate_Code generate_code.cpp)

add_library(Autotune_Wrapper autotune_wrapper.cpp)

add_library(Path_Finding path_finding.cpp)


add_executable(surface_code_test_sim surface_code_test_sim.cpp)
target_link_libraries(surface_code_test_sim Autotune_Wrapper Path_Finding autotune)


add_executable(main main.cpp)
target_link_libraries(main  Generate_Code Path_Finding Autotune_Wrapper autotune)

