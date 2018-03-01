
add_library(Generate_Code generate_code.cpp)

add_library(Autotune_Wrapper autotune_wrapper.cpp)

add_executable(main main.cpp)

target_link_libraries(main Generate_Code Autotune_Wrapper autotune)
