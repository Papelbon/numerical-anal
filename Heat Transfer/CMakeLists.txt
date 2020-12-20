cmake_minimum_required(VERSION 3.17)
project(heat_transfer)

enable_language(Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

file(GLOB TASK1 "task1.f90" "modules/*.f90")
file(GLOB TASK2 "task2.f90" "modules/*.f90")
file(GLOB TASK3 "task3.f90" "modules/*.f90")
file(GLOB TASK4 "task4.f90" "modules/*.f90")
file(GLOB TASK5 "task5.f90" "modules/*.f90")

add_executable(heat_transfer_task1 ${TASK1})
add_executable(heat_transfer_task2 ${TASK2})
add_executable(heat_transfer_task3 ${TASK3})
add_executable(heat_transfer_task4 ${TASK4})
add_executable(heat_transfer_task5 ${TASK5})













