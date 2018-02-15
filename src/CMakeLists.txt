add_library(Error_Class error_class.cpp)

add_library(Generate_Code generate_code.cpp)

add_library(Simulate_Errors simulate_errors.cpp)

add_executable(main main.cpp)

target_link_libraries(main Generate_Code Error_Class GraphSim LocCliff)
