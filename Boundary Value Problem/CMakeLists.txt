cmake_minimum_required(VERSION 3.17)
project(bv_problem)

enable_language(Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

file(GLOB TASK1 "task1.f90" "modules/*.f90")
file(GLOB TASK2 "task2.f90" "modules/*.f90")
file(GLOB TASK3 "task3.f90" "modules/*.f90")
file(GLOB TASK4 "task4.f90" "modules/*.f90")

add_executable(bv_problem_task1 ${TASK1})
add_executable(bv_problem_task2 ${TASK2})
add_executable(bv_problem_task3 ${TASK3})
add_executable(bv_problem_task4 ${TASK4})













