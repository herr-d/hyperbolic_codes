
add_library(Generate_Code generate_code.cpp)

add_library(Autotune_Wrapper autotune_wrapper.cpp)

add_library(Path_Finding path_finding.cpp)


add_executable(old_ex old_ex.cpp)
target_link_libraries(old_ex autotune ${BLOSSOM_PM} ${BLOSSOM_MISC} ${BLOSSOM_MIN_COST} ${BLOSSOM_GEOM})


add_executable(surface_code_test_sim surface_code_test_sim.cpp)
target_link_libraries(surface_code_test_sim Autotune_Wrapper Path_Finding autotune ${BLOSSOM_PM} ${BLOSSOM_MISC} ${BLOSSOM_MIN_COST} ${BLOSSOM_GEOM})


add_executable(main main.cpp)
target_link_libraries(main  Generate_Code Path_Finding Autotune_Wrapper autotune ${BLOSSOM_PM} ${BLOSSOM_MISC} ${BLOSSOM_MIN_COST} ${BLOSSOM_GEOM})

