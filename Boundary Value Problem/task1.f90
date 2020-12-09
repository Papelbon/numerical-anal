! ---------------------------------------------------------------------------------------------------------------------
! Задача 1.
!
! Краевая задача вида
! sin(12) * y'' + (1 + cos(12) * x^2) * y = -1,
! x_min = -1, x_max = 1.
!
! Граничные условия однородны и имеют вид
! y(x_min) = 0,
! y(x_max) = 0.
! ---------------------------------------------------------------------------------------------------------------------

program task1

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

    n = 3
    allocate(col_sol(1:n))

    call collocation_method(x_min=-1._dp, x_max=1._dp, n=n, solution=col_sol, solution_ext=col_sol_ext,   &
            f=f2, a=a2, p=p2, q=q2, basic=basic_v2, x_points=X_POINTS, is_print=.true., is_draw=.false., &
            info='Task 1. Collocation method')

    ! Интегральный метод наименьших квадратов -------------------------------------------------------------------------

    n = 3
    allocate(ils_sol(1:n))

    call int_least_squares_method(-1._dp, 1._dp, n, ils_sol, ils_sol_ext, f2, a2, p2, q2, basic_v2, &
            X_POINTS, .false., .false., 'Task 1. Integral Least Squares method')

    ! Дискретный метод наименьших квадратов ---------------------------------------------------------------------------

    n = 3
    allocate(dls_sol(1:n))

    call disc_least_squares_method(-1._dp, 1._dp, n, dls_sol, dls_sol_ext, f2, a2, p2, q2, basic_v2, &
            X_POINTS, .false., .false., 'Task 1. Discrete Least Squares method')

    ! Метод Галеркина -------------------------------------------------------------------------------------------------

    n = 3
    allocate(gal_sol(1:n))

    call galerkin_method(-1._dp, 1._dp, n, gal_sol, gal_sol_ext, f2, a2, p2, q2, basic_v2, &
            X_POINTS, .false., .false., 'Task 1. Galerkin method')

    ! Метод конечных разностей ----------------------------------------------------------------------------------------

    n = 30; n2 = 60; n3 = 5
    allocate(fd_sol(1:n+1), x_grid(1:n+1), fd_sol2(1:n2+1), x_grid2(1:n2+1), fd_sol3(1:n3+1), x_grid3(1:n3+1))

    call finite_difference_method(0._dp, 0._dp, -1._dp, 1._dp, n, f2, a2, p2, q2, 1._dp, 0._dp, 1._dp, 0._dp, &
            y_coef_dr1_bp1, y_coef_dr1_bc0, y_coef_dr1_bcn, x_grid, fd_sol, .false., .false., &
            'Task 1. Finite Difference method')

    call finite_difference_method(0._dp, 0._dp, -1._dp, 1._dp, n2, f2, a2, p2, q2, 1._dp, 0._dp, 1._dp, 0._dp, &
            y_coef_dr1_bp1, y_coef_dr1_bc0, y_coef_dr1_bcn, x_grid2, fd_sol2, .false., .false., &
            'Task 1. Finite Difference method')

    call finite_difference_method(0._dp, 0._dp, -1._dp, 1._dp, n3, f2, a2, p2, q2, 1._dp, 0._dp, 1._dp, 0._dp, &
            y_coef_dr1_bp1, y_coef_dr1_bc0, y_coef_dr1_bcn, x_grid3, fd_sol3, .true., .true., &
            'Task 1. Finite Difference method')

    ! Сравнительные графики -------------------------------------------------------------------------------------------

    call gp%title("Task 1. Finite Difference method")
    call gp%options('set key bottom left; set grid')
    call gp%xlabel('X grid')
    call gp%ylabel('Y solution')
    call gp%plot(x1=x_grid, y1=fd_sol + 1e-3, ls1='with lines title "n=30, +1e-3" lc rgb "#000"', &
                 x2=x_grid, y2=fd_sol - 1e-3, ls2='with lines title "n=30, -1e-3" lc rgb "#000"', &
                 x3=x_grid2, y3=fd_sol2, ls3='with lines title "n=60" lc rgb "#FFC300"',          &
                 x4=x_grid3, y4=fd_sol3, ls4='with lines title "n=5" lc rgb "#00AD09"')

    call noised_linspace(-1._dp, 1._dp, x_list, noise=0._dp)

    call gp%title("Task 1. All methods")
    call gp%options('set key top left; set grid')
    call gp%xlabel('X')
    call gp%ylabel('Y')
    call gp%plot(x1=x_list, y1=col_sol_ext, ls1='with lines title "Collocation method"',            &
                 x2=x_list, y2=ils_sol_ext, ls2='with lines title "Integral Least Squares method"', &
                 x3=x_list, y3=gal_sol_ext, ls3='with lines title "Galerkin method"',               &
                 x4=x_grid2, y4=fd_sol2, ls4='with lines title "Finite Difference method"')


    deallocate(fd_sol, x_grid, fd_sol2, x_grid2, fd_sol3, x_grid3, col_sol, ils_sol, dls_sol, gal_sol)

end program

