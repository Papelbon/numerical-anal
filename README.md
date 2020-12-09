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
![task1](https://user-images.githubusercontent.com/62307154/101655478-b8f93080-3a52-11eb-9785-2bd8ec22d53a.png)

- [example1](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task1.f90)

![task2](https://user-images.githubusercontent.com/62307154/101656386-b0552a00-3a53-11eb-88e2-002531193f12.png)

- [example2](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task2.f90)

![task3](https://user-images.githubusercontent.com/62307154/101656825-33768000-3a54-11eb-8aee-81c412ef1dd7.png)

- [example3](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task3.f90)

![task4](https://user-images.githubusercontent.com/62307154/101659674-62dabc00-3a57-11eb-9fad-0abfd7747e7a.png)

- [example4](https://github.com/Papelbon/numerical-anal/blob/main/Boundary%20Value%20Problem/task4.f90)
