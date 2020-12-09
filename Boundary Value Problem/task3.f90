! ---------------------------------------------------------------------------------------------------------------------
! Задача 3.
!
! Краевая задача вида
! y'' + (-3)* y' + (8 * x) * y = 8,
! x_min = 1.8, x_max = 3.8.
!
! Граничные условия имеют вид
! y(x_min)  + (-0.5) * y'(x_min) = 2,
! y(x_max) = 5.
!
! Для аппроксимации производных в граничных условиях использовать разностные отношения вида
! y'(0) = (-y(2) + 4 * y(1) - 3 * y(0)) / 2 / h
! y'(n) = (3 * y(n) - 4 * y(n-1) + y(n-2)) / 2 / h
! ---------------------------------------------------------------------------------------------------------------------
! Задача 3(test).
!
! Краевая задача вида
! y'' + (x^2)* y' + (-x) * y = (6 / x^4) - 3 / x,
! x_min = 1, x_max = 2.
!
! Граничные условия имеют вид
! y(x_min) = 1,
! 3 * y(x_max) + y'(x_max) = 0.5.
!
! Точное решение: y = 1 / x^2
! ---------------------------------------------------------------------------------------------------------------------

program task3

    use ogpf
    use base

    use finite_difference

    ! Ввод начальных данных -------------------------------------------------------------------------------------------

    implicit none

    type(gpf) :: gp

    integer, parameter :: dp = kind(0.d0)
    real(dp), dimension(:), allocatable :: fd_sol, x_grid, fd_sol2, x_grid2, fd_sol3, x_grid3, fd_sol4, x_grid4, y_exact
    integer :: i, j, n, n2, n3, n4

    ! Метод конечных разностей ----------------------------------------------------------------------------------------
    ! Задача 3 --------------------------------------------------------------------------------------------------------

    n = 100; n2 = 200; n3 = 5; n4 = 10
    allocate(fd_sol(1:n+1), x_grid(1:n+1), fd_sol2(1:n2+1), x_grid2(1:n2+1), fd_sol3(1:n3+1), x_grid3(1:n3+1), &
        fd_sol4(1:n4+1), x_grid4(1:n4+1))

    call finite_difference_method(2._dp, 5._dp, 1.8_dp, 3.8_dp, n, f3, a1, p3, q3, 1._dp, -0.5_dp, 1._dp, 0._dp, &
            y_coef_dr1_bp1, y_coef_dr1_bc0, y_coef_dr1_bcn, x_grid, fd_sol, .false., .false., &
            'Task 3. Finite Difference method')

    call finite_difference_method(2._dp, 5._dp, 1.8_dp, 3.8_dp, n2, f3, a1, p3, q3, 1._dp, -0.5_dp, 1._dp, 0._dp, &
            y_coef_dr1_bp1, y_coef_dr1_bc0, y_coef_dr1_bcn, x_grid2, fd_sol2, .false., .false., &
            'Task 3. Finite Difference method')

    call finite_difference_method(2._dp, 5._dp, 1.8_dp, 3.8_dp, n3, f3, a1, p3, q3, 1._dp, -0.5_dp, 1._dp, 0._dp, &
            y_coef_dr1_bp1, y_coef_dr1_bc0, y_coef_dr1_bcn, x_grid3, fd_sol3, .true., .false., &
            'Task 3. Finite Difference method')

    call finite_difference_method(2._dp, 5._dp, 1.8_dp, 3.8_dp, n4, f3, a1, p3, q3, 1._dp, -0.5_dp, 1._dp, 0._dp, &
            y_coef_dr1_bp1, y_coef_dr1_bc0, y_coef_dr1_bcn, x_grid4, fd_sol4, .false., .false., &
            'Task 3. Finite Difference method')

    ! Сравнительный график --------------------------------------------------------------------------------------------

    call gp%title("Task 3. Finite Difference method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('X grid')
    call gp%ylabel('Y solution')
    call gp%plot(x1=x_grid, y1=fd_sol + 0.1, ls1='with lines title "n=100, +0.1" lc rgb "#000"', &
            x2=x_grid, y2=fd_sol - 0.1, ls2='with lines title "n=100, -0.1" lc rgb "#000"', &
            x3=x_grid2, y3=fd_sol2, ls3='with lines title "n=200" lc rgb "#FFC300"')

    call gp%title("Task 3. Finite Difference method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('X grid')
    call gp%ylabel('Y solution')
    call gp%plot(x1=x_grid3, y1=fd_sol3, ls1='with lines title "n=5"', &
            x2=x_grid4, y2=fd_sol4, ls2='with lines title "n=10"', &
            x3=x_grid, y3=fd_sol, ls3='with lines title "n=100"')

    deallocate(fd_sol, x_grid, fd_sol2, x_grid2, fd_sol3, x_grid3, fd_sol4, x_grid4)

    ! Метод конечных разностей ----------------------------------------------------------------------------------------
    ! Задача 3(test) --------------------------------------------------------------------------------------------------

    n = 5
    allocate(fd_sol(1:n+1), x_grid(1:n+1), y_exact(1:n+1))

    call finite_difference_method(1._dp, 0.5_dp, 1._dp, 2._dp, n, f4, a1, p4, q4, 1._dp, 0._dp, 3._dp, 1._dp, &
            y_coef_dr1_bp1, y_coef_dr1_bc0, y_coef_dr1_bcn, x_grid, fd_sol, .false., .false., &
            'Task 3(test). Finite Difference method')

    ! Сравнительный график --------------------------------------------------------------------------------------------

    y_exact = 1 / (x_grid) ** 2

    call gp%title("Task 3(test). Finite Difference method")
    call gp%options('set key bottom left; set grid')
    call gp%xlabel('X grid')
    call gp%ylabel('Y solution')
    call gp%plot(x1=x_grid, y1=y_exact, ls1='with lines title "exact solution"', &
            x2=x_grid, y2=fd_sol, ls2='with lines title "numerical solution"')

    deallocate(fd_sol, x_grid, y_exact)

end program

