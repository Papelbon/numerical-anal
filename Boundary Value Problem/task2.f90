! ---------------------------------------------------------------------------------------------------------------------
! Задача 2.
!
! Краевая задача вида
! y'' + (-0.5 + sin(x)) * y' + (8 / (1 + 0.25 * x^2)) * y = 5 * (1 - x^2),
! x_min = 0, x_max = 2.
!
! Граничные условия однородны и имеют вид
! y(x_min) = 0,
! y(x_max) = 0.
! ---------------------------------------------------------------------------------------------------------------------

program task2

    use ogpf
    use base

    use collocation
    use least_squares
    use galerkin
    use finite_difference

    ! Ввод начальных данных -------------------------------------------------------------------------------------------

    implicit none

    type(gpf) :: gp

    integer, parameter :: dp = kind(0.d0), X_POINTS = 100
    real(dp), dimension(:), allocatable :: col_sol, ils_sol, dls_sol, gal_sol, &
            fd_sol, x_grid, fd_sol2, x_grid2, fd_sol3, x_grid3
    real(dp), dimension(1:X_POINTS) :: col_sol_ext, ils_sol_ext, dls_sol_ext, gal_sol_ext, x_list
    integer :: i, j, n, n2, n3

    ! Метод коллокаций ------------------------------------------------------------------------------------------------

    n = 7
    allocate(col_sol(1:n))

    call collocation_method(x_min=0.000001_dp, x_max=2._dp, n=n, solution=col_sol, solution_ext=col_sol_ext,   &
            f=f1, a=a1, p=p1, q=q1, basic=basic_v3, x_points=X_POINTS, is_print=.false., is_draw=.false., &
            info='Task 2. Collocation method')

    ! Интегральный метод наименьших квадратов -------------------------------------------------------------------------

    n = 7
    allocate(ils_sol(1:n))

    call int_least_squares_method(0.000001_dp, 2._dp, n, ils_sol, ils_sol_ext, f1, a1, p1, q1, basic_v3, &
            X_POINTS, .false., .false., 'Task 2. Integral Least Squares method')

    ! Дискретный метод наименьших квадратов ---------------------------------------------------------------------------

    n = 7
    allocate(dls_sol(1:n))

    call disc_least_squares_method(0.000001_dp, 2._dp, n, dls_sol, dls_sol_ext, f1, a1, p1, q1, basic_v3, &
            X_POINTS, .false., .false., 'Task 2. Discrete Least Squares method')

    ! Метод Галеркина -------------------------------------------------------------------------------------------------

    n = 7
    allocate(gal_sol(1:n))

    call galerkin_method(0.000001_dp, 2._dp, n, gal_sol, gal_sol_ext, f1, a1, p1, q1, basic_v3, &
            X_POINTS, .false., .false., 'Task 2. Galerkin method')

    ! Метод конечных разностей ----------------------------------------------------------------------------------------

    n = 30; n2 = 60; n3 = 5
    allocate(fd_sol(1:n+1), x_grid(1:n+1), fd_sol2(1:n2+1), x_grid2(1:n2+1), fd_sol3(1:n3+1), x_grid3(1:n3+1))

    call finite_difference_method(0._dp, 0._dp, 0._dp, 2._dp, n, f1, a1, p1, q1, 1._dp, 0._dp, 1._dp, 0._dp, &
            y_coef_dr1_bp1, y_coef_dr1_bc0, y_coef_dr1_bcn, x_grid, fd_sol, .false., .false., &
            'Task 2. Finite Difference method')

    call finite_difference_method(0._dp, 0._dp, 0._dp, 2._dp, n2, f1, a1, p1, q1, 1._dp, 0._dp, 1._dp, 0._dp, &
            y_coef_dr1_bp1, y_coef_dr1_bc0, y_coef_dr1_bcn, x_grid2, fd_sol2, .false., .false., &
            'Task 2. Finite Difference method')

    call finite_difference_method(0._dp, 0._dp, 0._dp, 2._dp, n3, f1, a1, p1, q1, 1._dp, 0._dp, 1._dp, 0._dp, &
            y_coef_dr1_bp1, y_coef_dr1_bc0, y_coef_dr1_bcn, x_grid3, fd_sol3, .true., .false., &
            'Task 2. Finite Difference method')

    ! Сравнительные графики -------------------------------------------------------------------------------------------

    call gp%title("Task 2. Finite Difference method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('X grid')
    call gp%ylabel('Y solution')
    call gp%plot(x1=x_grid, y1=fd_sol + 1e-2, ls1='with lines title "n=30, +1e-2" lc rgb "#000"', &
            x2=x_grid, y2=fd_sol - 1e-2, ls2='with lines title "n=30, -1e-2" lc rgb "#000"', &
            x3=x_grid2, y3=fd_sol2, ls3='with lines title "n=60" lc rgb "#FFC300"',          &
            x4=x_grid3, y4=fd_sol3, ls4='with lines title "n=5" lc rgb "#00AD09"')

    call noised_linspace(0._dp, 2._dp, x_list, noise=0._dp)

    call gp%title("Task 2. All methods")
    call gp%options('set key top left; set grid')
    call gp%xlabel('X')
    call gp%ylabel('Y')
    call gp%plot(x1=x_list, y1=col_sol_ext, ls1='with lines title "Collocation method"',            &
            x2=x_list, y2=dls_sol_ext, ls2='with lines title "Discrete Least Squares method"', &
            x3=x_list, y3=gal_sol_ext, ls3='with lines title "Galerkin method"',               &
            x4=x_grid2, y4=fd_sol2, ls4='with lines title "Finite Difference method"')


    deallocate(fd_sol, x_grid, fd_sol2, x_grid2, fd_sol3, x_grid3, col_sol, ils_sol, dls_sol, gal_sol)

end program

