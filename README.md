## Numerical Analysis Problems and Solutions

### Boundary Value Problem
Solve boundary value problem of the form

![bv_problem](https://user-images.githubusercontent.com/62307154/102728720-7848c380-433e-11eb-8f7c-a20760905ace.png)
1. Collocation method
```fortran
subroutine collocation_method(x_min, x_max, n, solution, solution_ext, f, a, p, q, basic, &
            x_points, is_print, is_draw, info)
```
2. Integral Least Squares method
```fortran
subroutine int_least_squares_method(x_min, x_max, n, solution, solution_ext, f, a, p, q, basic, &
            x_points, is_print, is_draw, info)
```
3. Discrete Least Squares method 
```fortran
subroutine disc_least_squares_method(x_min, x_max, n, solution, solution_ext, f, a, p, q, basic, &
            x_points, is_print, is_draw, info)
```
4. Galerkin method
```fortran
subroutine galerkin_method(x_min, x_max, n, solution, solution_ext, f, a, p, q, basic, &
            x_points, is_print, is_draw, info)
```
5. Finite Difference method
```fortran
subroutine finite_difference_method(AA, BB, x_min, x_max, n, f, a, p, q, alpha1, beta1, &
            alpha2, beta2, y_coef, y_coef_0, y_coef_n, x_grid, sol, is_print, is_draw, info)
```
#### Examples
![bv_task1](https://user-images.githubusercontent.com/62307154/102730439-d9c06080-4345-11eb-971c-030807730a08.png)

- [example1](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task1.f90)

![bv_task2](https://user-images.githubusercontent.com/62307154/102729016-d1652700-433f-11eb-9ee4-e2ea0b139f18.png)

- [example2](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task2.f90)

![bv_task3](https://user-images.githubusercontent.com/62307154/102729230-eb533980-4340-11eb-97f7-25084775cbdc.png)

- [example3](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task3.f90)

![bv_task4](https://user-images.githubusercontent.com/62307154/102729860-83eab900-4343-11eb-948a-ab2a26330f8d.png)

- [example4](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task4.f90)

### Heat Transfer
Estimate the solution of the time-dependent (time-independent) heat equation over a one dimensional region.

Time-independent form:

![time-independent](https://user-images.githubusercontent.com/62307154/102730324-69194400-4345-11eb-9b79-c7f54b7daf08.png)

Time-dependent form:

![time-dependent](https://user-images.githubusercontent.com/62307154/102730899-44be6700-4347-11eb-894f-a5c407dacb8e.png)
1. Finite Difference methods (1D Steady State Heat Equation)
```fortran
subroutine fd_centered_heat_transfer_steady(n, a, b, Ua, Ub, k, f, x_grid, sol, &
            is_print, is_draw, info)
```
```fortran
subroutine fd_balance_heat_transfer_steady(n, a, b, Ua, Ub, k, f, x_points, sol, &
            is_print, is_draw, info, psources_)
```
2. Finite Difference method (Time Dependent 1D Heat Equation using Explicit Time Stepping)
```fortran
subroutine fd_heat_transfer_explicit(n, m, a, b, t, g1, g2, phi, k, f, x_list, t_list, &
            matrix, is_check_cfl, bc, is_print, is_draw, info)
```
3. Finite Difference method (Time Dependent 1D Heat Equation using Implicit Time Stepping)
```fortran
subroutine fd_heat_transfer_implicit(n, m, a, b, t, g1, g2, phi, k, f, x_list, t_list, &
            matrix, bc, is_print, is_draw, info)
```
#### Examples
1. Modeling of steady-state thermal conductivity processes.

- [example1](https://github.com/Papelbon/numerical-anal/blob/main/Heat%20Transfer/task1.f90)
- [example2](https://github.com/Papelbon/numerical-anal/blob/main/Heat%20Transfer/task2.f90)

2. Modeling of unsteady-state thermal conductivity processes.

![task3_3d_2](https://user-images.githubusercontent.com/62307154/102731937-4e959980-434a-11eb-8c02-4cdda2ddb8ce.gif)

![task3_3d_3](https://user-images.githubusercontent.com/62307154/102731951-5c4b1f00-434a-11eb-8b00-1bc313605605.gif)

- [example3](https://github.com/Papelbon/numerical-anal/blob/main/Heat%20Transfer/task3.f90)
- [example4](https://github.com/Papelbon/numerical-anal/blob/main/Heat%20Transfer/task4.f90)
- [example5](https://github.com/Papelbon/numerical-anal/blob/main/Heat%20Transfer/task5.f90)
