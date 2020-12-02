cmake_minimum_required(VERSION 3.17)

project(lab11)

enable_language(Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

add_executable(
        lab11.x
        modules/ogpf.f90
        modules/base.f90
        modules/qr.f90
        modules/lu.f90
        modules/integrate.f90
        modules/collocation.f90
        modules/least_squares.f90
        modules/galerkin.f90
        main.f90
)
















