## Numerical Analysis Problems and Solutions

### Boundary Value Problem
Solve boundary value problem of the form

![bv_problem](https://user-images.githubusercontent.com/62307154/101572047-fc1eb980-39e8-11eb-866a-cdca0f3c25f2.png)
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
- [example1](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task1.f90)

- [example2](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task2.f90)

- [example3](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task3.f90)

- [example4](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task4.f90)
