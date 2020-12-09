! ---------------------------------------------------------------------------------------------------------------------
! Рассматривается краевая задача вида
! -(k(x) * u')' + q(x) * u = f(x), x_min = a, x_max = b,
! -k(a) * u'(a) + 0.5 * u(a) = A,
! k(b) * u'(b) + 0.5 * u(b) = B.
! ---------------------------------------------------------------------------------------------------------------------
! Задача 4.
!
! a = 0, b = 2, c = 1.125, A = 0, B = 0,
! if a < x < c then k(x) = 1.8 and q(x) = 6.5,
! if c < x < b then k(x) = 0.6 and q(x) = 7.8,
! f(x) = 8 * x * (2.5 - x)
!
! Для аппроксимации производных в граничных условиях использовать метод баланса
! ---------------------------------------------------------------------------------------------------------------------

program task4

    use ogpf
    use base

    use finite_difference

    ! Ввод начальных данных -------------------------------------------------------------------------------------------

    implicit none

    type(gpf) :: gp

    integer, parameter :: dp = kind(0.d0)
    real(dp), dimension(:), allocatable :: fd_sol_1, fd_sol_2, x_grid_1, x_grid_2, fd_sol2_1, fd_sol2_2, &
            x_grid2_1, x_grid2_2, fd_sol3, x_grid3
    integer :: i, j, n_1, n_2, n2_1, n2_2, n3

    ! Метод конечных разностей ----------------------------------------------------------------------------------------

    n_1 = 500; n2_1 = 1000; n_2 = 70; n2_2 = 140; n3 = 5
    allocate(fd_sol_1(1:n_1+1), x_grid_1(1:n_1+1), fd_sol2_1(1:n2_1+1), x_grid2_1(1:n2_1+1), &
            fd_sol_2(1:n_2+1), x_grid_2(1:n_2+1), fd_sol2_2(1:n2_2+1), x_grid2_2(1:n2_2+1), &
            fd_sol3(1:n3+1), x_grid3(1:n3+1))

    call finite_difference_method(0._dp,0._dp,0._dp,2._dp,n_1,f6,a3,p2,q6,0.5_dp,-k(0._dp),0.5_dp,k(2._dp), &
            y_coef_dr1_bp1,y_coef_dr2_bc0,y_coef_dr2_bcn,x_grid_1,fd_sol_1,.false.,.false., &
            'Task 4. Finite Difference method')

    call finite_difference_method(0._dp,0._dp,0._dp,2._dp,n2_1,f6,a3,p2,q6,0.5_dp,-k(0._dp),0.5_dp,k(2._dp), &
            y_coef_dr1_bp1,y_coef_dr2_bc0,y_coef_dr2_bcn,x_grid2_1,fd_sol2_1,.false.,.false., &
            'Task 4. Finite Difference method')

    call finite_difference_method(0._dp,0._dp,0._dp,2._dp,n_2,f6,a3,p2,q6,0.5_dp,-k(0._dp),0.5_dp,k(2._dp), &
            y_coef_dr1_bp1,y_coef_dr1_bc0,y_coef_dr1_bcn,x_grid_2,fd_sol_2,.false.,.false., &
            'Task 4. Finite Difference method')

    call finite_difference_method(0._dp,0._dp,0._dp,2._dp,n2_2,f6,a3,p2,q6,0.5_dp,-k(0._dp),0.5_dp,k(2._dp), &
            y_coef_dr1_bp1,y_coef_dr1_bc0,y_coef_dr1_bcn,x_grid2_2,fd_sol2_2,.false.,.false., &
            'Task 4. Finite Difference method')

    call finite_difference_method(0._dp,0._dp,0._dp,2._dp,n3,f6,a3,p2,q6,0.5_dp,-k(0._dp),0.5_dp,k(2._dp), &
            y_coef_dr1_bp1,y_coef_dr2_bc0,y_coef_dr2_bcn,x_grid3,fd_sol3,.true.,.true., &
            'Task 4. Finite Difference method')

    ! Сравнительный график --------------------------------------------------------------------------------------------

    call gp%title("Task 4. Finite Difference method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('X grid')
    call gp%ylabel('Y solution')
    call gp%plot(x1=x_grid_1, y1=fd_sol_1 + 1e-3, ls1='with lines title "n=500, +1e-3" lc rgb "#000"', &
            x2=x_grid_1, y2=fd_sol_1 - 1e-3, ls2='with lines title "n=500, -1e-3" lc rgb "#000"', &
            x3=x_grid2_1, y3=fd_sol2_1, ls3='with lines title "balance (n=1000)" lc rgb "#FFC300"')

    call gp%title("Task 4. Finite Difference method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('X grid')
    call gp%ylabel('Y solution')
    call gp%plot(x1=x_grid_2, y1=fd_sol_2 + 1e-3, ls1='with lines title "n=70, +1e-3" lc rgb "#000"', &
            x2=x_grid_2, y2=fd_sol_2 - 1e-3, ls2='with lines title "n=70, -1e-3" lc rgb "#000"', &
            x3=x_grid2_2, y3=fd_sol2_2, ls3='with lines title "sec-order diff (n=140)" lc rgb "#FFC300"')


    deallocate(fd_sol_1, fd_sol_2, x_grid_1, x_grid_2, fd_sol2_1, fd_sol2_2, &
            x_grid2_1, x_grid2_2, fd_sol3, x_grid3)


end program

