program task5

    use ogpf
    use base
    use heat_transfer

    ! Ввод начальных данных -------------------------------------------------------------------------------------------

    implicit none

    type(gpf) :: gp
    integer, parameter :: dp = kind(0.d0)
    real(dp), parameter :: a = -1._dp, b = 1._dp, t = 0.4_dp, t1 = 0.1_dp, t2 = 0.2_dp

    real(dp), allocatable :: x_list(:), t_list(:), u(:,:), x_list2(:), t_list2(:), u2(:,:), err_t1(:), err_t2(:), &
            errors_t1(:), errors_t2(:), dt_array(:), dx_array(:)

    integer :: i, n, n2, m, m2
    real(dp) :: dt, dx

    ! Finite Difference Explicit Heat Transfer method -----------------------------------------------------------------
    ! Явная схема. Левая разность (h - фиксированное) -----------------------------------------------------------------

    print *, 'Task 5. Finite Difference Explicit Heat Transfer method. Left difference (h is fixed):'
    print *, '  N             τ                   s(t=t1)                   s(t=t2)                 ' // &
            'max(t=t1)                 max(t=t2)'
    print *, '---  ------------  ------------------------  ------------------------  ' // &
            '------------------------  ------------------------'

    n = 10; m = 20
    allocate(errors_t1(5), errors_t2(5), dt_array(5))
    do i = 1, 5
        m = m * 2
        m2 = m * 2

        allocate(u(m+1,n+1), x_list(n+1), t_list(m+1), u2(m2+1,n+1), x_list2(n+1), t_list2(m2+1), err_t1(n+1), err_t2(n+1))

        call fd_heat_transfer_explicit(n, m, a, b, t, g1, g2, phi, k, f, x_list, t_list, u, bc=1, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Left difference'))

        call fd_heat_transfer_explicit(n, m2, a, b, t, g1, g2, phi, k, f, x_list2, t_list2, u2, bc=1, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Left difference'))

        dt = t / m; dt_array(i) = dt
        err_t1 = u(int(t1/dt),:) - u2(2*int(t1/dt),:); errors_t1(i) = norm2(err_t1)
        err_t2 = u(int(t2/dt),:) - u2(2*int(t2/dt),:); errors_t2(i) = norm2(err_t2)

        print '(i4,a2,f12.7,a2,f24.16,a2,f24.16,a2,f24.16,a2,f24.16)',n,'  ',dt,'  ',std(err_t1),'  ',std(err_t2), &
                '  ',maxval(abs(err_t1)),'  ',maxval(abs(err_t2))
        deallocate(u, x_list, t_list, u2, x_list2, t_list2, err_t1, err_t2)
    end do

    call gp%title("Task 5. Finite Difference Explicit Heat Transfer method. Left difference (h is fixed)")
    call gp%options('set key top left; set grid')
    call gp%xlabel('τ')
    call gp%ylabel('error')
    call gp%plot(x1=dt_array, y1=errors_t1, ls1='with lines title "t=t1"', &
            x2=dt_array, y2=errors_t2, ls2='with lines title "t=t2"')

    deallocate(errors_t1, errors_t2, dt_array)

    ! Finite Difference Explicit Heat Transfer method -----------------------------------------------------------------
    ! Явная схема. Левая разность (τ - фиксированное) -----------------------------------------------------------------

    print *, ''
    print *, 'Task 5. Finite Difference Explicit Heat Transfer method. Left difference (τ is fixed):'
    print *, '  N             τ                   s(t=t1)                   s(t=t2)                 ' // &
            'max(t=t1)                 max(t=t2)'
    print *, '---  ------------  ------------------------  ------------------------  ' // &
            '------------------------  ------------------------'

    n = 3; m = 4000
    allocate(errors_t1(5), errors_t2(5), dx_array(5))
    do i = 1, 5
        n = n * 2
        n2 = n * 2

        allocate(u(m+1,n+1), x_list(n+1), t_list(m+1), u2(m+1,n2+1), x_list2(n2+1), t_list2(m+1), err_t1(n+1), err_t2(n+1))

        call fd_heat_transfer_explicit(n, m, a, b, t, g1, g2, phi, k, f, x_list, t_list, u, bc=1, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Left difference'))

        call fd_heat_transfer_explicit(n2, m, a, b, t, g1, g2, phi, k, f, x_list2, t_list2, u2, bc=1, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Left difference'))

        dx = (b - a) / n; dt = t / m; dx_array(i) = dx
        err_t1 = u(int(t1/dt),:) - u2(int(t1/dt),1:n2+1:2); errors_t1(i) = norm2(err_t1)
        err_t2 = u(int(t2/dt),:) - u2(int(t2/dt),1:n2+1:2); errors_t2(i) = norm2(err_t2)

        print '(i4,a2,f12.7,a2,f24.16,a2,f24.16,a2,f24.16,a2,f24.16)',n,'  ',dt,'  ',std(err_t1),'  ',std(err_t2), &
                '  ',maxval(abs(err_t1)),'  ',maxval(abs(err_t2))
        deallocate(u, x_list, t_list, u2, x_list2, t_list2, err_t1, err_t2)
    end do

    call gp%title("Task 5. Finite Difference Explicit Heat Transfer method. Left difference (τ is fixed)")
    call gp%options('set key top left; set grid')
    call gp%xlabel('h')
    call gp%ylabel('error')
    call gp%plot(x1=dx_array, y1=errors_t1, ls1='with lines title "t=t1"', &
            x2=dx_array, y2=errors_t2, ls2='with lines title "t=t2"')

    deallocate(errors_t1, errors_t2, dx_array)

    ! Finite Difference Explicit Heat Transfer method -----------------------------------------------------------------
    ! Явная схема. Левая разность (τ = h^2 / 6) -----------------------------------------------------------------------

    print *, ''
    print *, 'Task 5. Finite Difference Explicit Heat Transfer method. Left difference (τ = h^2 / 6):'
    print *, '  N             τ                   s(t=t1)                   s(t=t2)                 ' // &
            'max(t=t1)                 max(t=t2)'
    print *, '---  ------------  ------------------------  ------------------------  ' // &
            '------------------------  ------------------------'

    n = 3
    allocate(errors_t1(5), errors_t2(5), dx_array(5))
    do i = 1, 5
        n = n * 2; n2 = n * 2
        dx = (b - a) / n; dt = dx * dx / 6._dp; m = t / dt

        allocate(u(m+1,n+1), x_list(n+1), t_list(m+1), u2(m+1,n2+1), x_list2(n2+1), t_list2(m+1), err_t1(n+1), err_t2(n+1))

        call fd_heat_transfer_explicit(n, m, a, b, t, g1, g2, phi, k, f, x_list, t_list, u, bc=1, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Left difference'))

        call fd_heat_transfer_explicit(n2, m, a, b, t, g1, g2, phi, k, f, x_list2, t_list2, u2, bc=1, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Left difference'))

        dx_array(i) = dx
        err_t1 = u(int(t1/dt),:) - u2(int(t1/dt),1:n2+1:2); errors_t1(i) = norm2(err_t1)
        err_t2 = u(int(t2/dt),:) - u2(int(t2/dt),1:n2+1:2); errors_t2(i) = norm2(err_t2)

        print '(i4,a2,f12.10,a2,f24.16,a2,f24.16,a2,f24.16,a2,f24.16)',n,'  ',dt,'  ',std(err_t1),'  ',std(err_t2), &
                '  ',maxval(abs(err_t1)),'  ',maxval(abs(err_t2))
        deallocate(u, x_list, t_list, u2, x_list2, t_list2, err_t1, err_t2)
    end do

    call gp%title("Task 5. Finite Difference Explicit Heat Transfer method. Left difference (τ = h^2 / 6)")
    call gp%options('set key top left; set grid')
    call gp%xlabel('h')
    call gp%ylabel('error')
    call gp%plot(x1=dx_array, y1=errors_t1, ls1='with lines title "t=t1"', &
            x2=dx_array, y2=errors_t2, ls2='with lines title "t=t2"')

    deallocate(errors_t1, errors_t2, dx_array)

    ! Finite Difference Explicit Heat Transfer method -----------------------------------------------------------------
    ! Явная схема. Центральная разность (h - фиксированное) -----------------------------------------------------------

    print *, ''
    print *, 'Task 5. Finite Difference Explicit Heat Transfer method. Central difference (h is fixed):'
    print *, '  N             τ                   s(t=t1)                   s(t=t2)                 ' // &
            'max(t=t1)                 max(t=t2)'
    print *, '---  ------------  ------------------------  ------------------------  ' // &
            '------------------------  ------------------------'

    n = 10; m = 20
    allocate(errors_t1(5), errors_t2(5), dt_array(5))
    do i = 1, 5
        m = m * 2
        m2 = m * 2

        allocate(u(m+1,n+1), x_list(n+1), t_list(m+1), u2(m2+1,n+1), x_list2(n+1), t_list2(m2+1), err_t1(n+1), err_t2(n+1))

        call fd_heat_transfer_explicit(n, m, a, b, t, g1, g2, phi, k, f, x_list, t_list, u, bc=2, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Central difference'))

        call fd_heat_transfer_explicit(n, m2, a, b, t, g1, g2, phi, k, f, x_list2, t_list2, u2, bc=2, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Central difference'))

        dt = t / m; dt_array(i) = dt
        err_t1 = u(int(t1/dt),:) - u2(2*int(t1/dt),:); errors_t1(i) = norm2(err_t1)
        err_t2 = u(int(t2/dt),:) - u2(2*int(t2/dt),:); errors_t2(i) = norm2(err_t2)

        print '(i4,a2,f12.7,a2,f24.16,a2,f24.16,a2,f24.16,a2,f24.16)',n,'  ',dt,'  ',std(err_t1),'  ',std(err_t2), &
                '  ',maxval(abs(err_t1)),'  ',maxval(abs(err_t2))
        deallocate(u, x_list, t_list, u2, x_list2, t_list2, err_t1, err_t2)
    end do

    call gp%title("Task 5. Finite Difference Explicit Heat Transfer method. Central difference (h is fixed)")
    call gp%options('set key top left; set grid')
    call gp%xlabel('τ')
    call gp%ylabel('error')
    call gp%plot(x1=dt_array, y1=errors_t1, ls1='with lines title "t=t1"', &
            x2=dt_array, y2=errors_t2, ls2='with lines title "t=t2"')

    deallocate(errors_t1, errors_t2, dt_array)

    ! Finite Difference Explicit Heat Transfer method -----------------------------------------------------------------
    ! Явная схема. Центральная разность (τ - фиксированное) -----------------------------------------------------------

    print *, ''
    print *, 'Task 5. Finite Difference Explicit Heat Transfer method. Central difference (τ is fixed):'
    print *, '  N             τ                   s(t=t1)                   s(t=t2)                 ' // &
            'max(t=t1)                 max(t=t2)'
    print *, '---  ------------  ------------------------  ------------------------  ' // &
            '------------------------  ------------------------'

    n = 3; m = 4000
    allocate(errors_t1(5), errors_t2(5), dx_array(5))
    do i = 1, 5
        n = n * 2
        n2 = n * 2

        allocate(u(m+1,n+1), x_list(n+1), t_list(m+1), u2(m+1,n2+1), x_list2(n2+1), t_list2(m+1), err_t1(n+1), err_t2(n+1))

        call fd_heat_transfer_explicit(n, m, a, b, t, g1, g2, phi, k, f, x_list, t_list, u, bc=2, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Central difference'))

        call fd_heat_transfer_explicit(n2, m, a, b, t, g1, g2, phi, k, f, x_list2, t_list2, u2, bc=2, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Central difference'))

        dx = (b - a) / n; dt = t / m; dx_array(i) = dx
        err_t1 = u(int(t1/dt),:) - u2(int(t1/dt),1:n2+1:2); errors_t1(i) = norm2(err_t1)
        err_t2 = u(int(t2/dt),:) - u2(int(t2/dt),1:n2+1:2); errors_t2(i) = norm2(err_t2)

        print '(i4,a2,f12.7,a2,f24.16,a2,f24.16,a2,f24.16,a2,f24.16)',n,'  ',dt,'  ',std(err_t1),'  ',std(err_t2), &
                '  ',maxval(abs(err_t1)),'  ',maxval(abs(err_t2))
        deallocate(u, x_list, t_list, u2, x_list2, t_list2, err_t1, err_t2)
    end do

    call gp%title("Task 5. Finite Difference Explicit Heat Transfer method. Central difference (τ is fixed)")
    call gp%options('set key top left; set grid')
    call gp%xlabel('h')
    call gp%ylabel('error')
    call gp%plot(x1=dx_array, y1=errors_t1, ls1='with lines title "t=t1"', &
            x2=dx_array, y2=errors_t2, ls2='with lines title "t=t2"')

    deallocate(errors_t1, errors_t2, dx_array)

    ! Finite Difference Explicit Heat Transfer method -----------------------------------------------------------------
    ! Явная схема. Центральная разность (τ = h^2 / 6) -----------------------------------------------------------------

    print *, ''
    print *, 'Task 5. Finite Difference Explicit Heat Transfer method. Central difference (τ = h^2 / 6):'
    print *, '  N             τ                   s(t=t1)                   s(t=t2)                 ' // &
            'max(t=t1)                 max(t=t2)'
    print *, '---  ------------  ------------------------  ------------------------  ' // &
            '------------------------  ------------------------'

    n = 3
    allocate(errors_t1(5), errors_t2(5), dx_array(5))
    do i = 1, 5
        n = n * 2; n2 = n * 2
        dx = (b - a) / n; dt = dx * dx / 6._dp; m = t / dt

        allocate(u(m+1,n+1), x_list(n+1), t_list(m+1), u2(m+1,n2+1), x_list2(n2+1), t_list2(m+1), err_t1(n+1), err_t2(n+1))

        call fd_heat_transfer_explicit(n, m, a, b, t, g1, g2, phi, k, f, x_list, t_list, u, bc=2, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Central difference'))

        call fd_heat_transfer_explicit(n2, m, a, b, t, g1, g2, phi, k, f, x_list2, t_list2, u2, bc=2, &
                is_check_cfl=.true., is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Explicit Heat Transfer method. Central difference'))

        dx_array(i) = dx
        err_t1 = u(int(t1/dt),:) - u2(int(t1/dt),1:n2+1:2); errors_t1(i) = norm2(err_t1)
        err_t2 = u(int(t2/dt),:) - u2(int(t2/dt),1:n2+1:2); errors_t2(i) = norm2(err_t2)

        print '(i4,a2,f12.10,a2,f24.16,a2,f24.16,a2,f24.16,a2,f24.16)',n,'  ',dt,'  ',std(err_t1),'  ',std(err_t2), &
                '  ',maxval(abs(err_t1)),'  ',maxval(abs(err_t2))
        deallocate(u, x_list, t_list, u2, x_list2, t_list2, err_t1, err_t2)
    end do

    call gp%title("Task 5. Finite Difference Explicit Heat Transfer method. Central difference (τ = h^2 / 6)")
    call gp%options('set key top left; set grid')
    call gp%xlabel('h')
    call gp%ylabel('error')
    call gp%plot(x1=dx_array, y1=errors_t1, ls1='with lines title "t=t1"', &
            x2=dx_array, y2=errors_t2, ls2='with lines title "t=t2"')

    deallocate(errors_t1, errors_t2, dx_array)

    ! Finite Difference Implicit Heat Transfer method -----------------------------------------------------------------
    ! Неявная схема. Левая разность (h - фиксированное) ---------------------------------------------------------------

    print *, ''
    print *, 'Task 5. Finite Difference Implicit Heat Transfer method. Left difference (h is fixed):'
    print *, '  N             τ                   s(t=t1)                   s(t=t2)                 ' // &
            'max(t=t1)                 max(t=t2)'
    print *, '---  ------------  ------------------------  ------------------------  ' // &
            '------------------------  ------------------------'

    n = 10; m = 20
    allocate(errors_t1(5), errors_t2(5), dt_array(5))
    do i = 1, 5
        m = m * 2
        m2 = m * 2

        allocate(u(m+1,n+1), x_list(n+1), t_list(m+1), u2(m2+1,n+1), x_list2(n+1), t_list2(m2+1), err_t1(n+1), err_t2(n+1))

        call fd_heat_transfer_implicit(n, m, a, b, t, g1, g2, phi, k(1._dp), f, x_list, t_list, u, bc=1, &
                is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Implicit Heat Transfer method. Left difference'))

        call fd_heat_transfer_implicit(n, m2, a, b, t, g1, g2, phi, k(1._dp), f, x_list2, t_list2, u2, bc=1, &
                is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Implicit Heat Transfer method. Left difference'))

        dt = t / m; dt_array(i) = dt
        err_t1 = u(int(t1/dt),:) - u2(2*int(t1/dt),:); errors_t1(i) = norm2(err_t1)
        err_t2 = u(int(t2/dt),:) - u2(2*int(t2/dt),:); errors_t2(i) = norm2(err_t2)

        print '(i4,a2,f12.7,a2,f24.16,a2,f24.16,a2,f24.16,a2,f24.16)',n,'  ',dt,'  ',std(err_t1),'  ',std(err_t2), &
                '  ',maxval(abs(err_t1)),'  ',maxval(abs(err_t2))
        deallocate(u, x_list, t_list, u2, x_list2, t_list2, err_t1, err_t2)
    end do

    call gp%title("Task 5. Finite Difference Implicit Heat Transfer method. Left difference (h is fixed)")
    call gp%options('set key top left; set grid')
    call gp%xlabel('τ')
    call gp%ylabel('error')
    call gp%plot(x1=dt_array, y1=errors_t1, ls1='with lines title "t=t1"', &
            x2=dt_array, y2=errors_t2, ls2='with lines title "t=t2"')

    deallocate(errors_t1, errors_t2, dt_array)

    ! Finite Difference Implicit Heat Transfer method -----------------------------------------------------------------
    ! Неявная схема. Левая разность (τ - фиксированное) ---------------------------------------------------------------

    print *, ''
    print *, 'Task 5. Finite Difference Implicit Heat Transfer method. Left difference (τ is fixed):'
    print *, '  N             τ                   s(t=t1)                   s(t=t2)                 ' // &
            'max(t=t1)                 max(t=t2)'
    print *, '---  ------------  ------------------------  ------------------------  ' // &
            '------------------------  ------------------------'

    n = 3; m = 4000
    allocate(errors_t1(5), errors_t2(5), dx_array(5))
    do i = 1, 5
        n = n * 2
        n2 = n * 2

        allocate(u(m+1,n+1), x_list(n+1), t_list(m+1), u2(m+1,n2+1), x_list2(n2+1), t_list2(m+1), err_t1(n+1), err_t2(n+1))

        call fd_heat_transfer_implicit(n, m, a, b, t, g1, g2, phi, k(1._dp), f, x_list, t_list, u, bc=1, &
                is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Implicit Heat Transfer method. Left difference'))

        call fd_heat_transfer_implicit(n2, m, a, b, t, g1, g2, phi, k(1._dp), f, x_list2, t_list2, u2, bc=1, &
                is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Implicit Heat Transfer method. Left difference'))

        dx = (b - a) / n; dt = t / m; dx_array(i) = dx
        err_t1 = u(int(t1/dt),:) - u2(int(t1/dt),1:n2+1:2); errors_t1(i) = norm2(err_t1)
        err_t2 = u(int(t2/dt),:) - u2(int(t2/dt),1:n2+1:2); errors_t2(i) = norm2(err_t2)

        print '(i4,a2,f12.7,a2,f24.16,a2,f24.16,a2,f24.16,a2,f24.16)',n,'  ',dt,'  ',std(err_t1),'  ',std(err_t2), &
                '  ',maxval(abs(err_t1)),'  ',maxval(abs(err_t2))
        deallocate(u, x_list, t_list, u2, x_list2, t_list2, err_t1, err_t2)
    end do

    call gp%title("Task 5. Finite Difference Implicit Heat Transfer method. Left difference (τ is fixed)")
    call gp%options('set key top left; set grid')
    call gp%xlabel('h')
    call gp%ylabel('error')
    call gp%plot(x1=dx_array, y1=errors_t1, ls1='with lines title "t=t1"', &
            x2=dx_array, y2=errors_t2, ls2='with lines title "t=t2"')

    deallocate(errors_t1, errors_t2, dx_array)

    ! Finite Difference Implicit Heat Transfer method -----------------------------------------------------------------
    ! Неявная схема. Центральная разность (h - фиксированное) ---------------------------------------------------------

    print *, ''
    print *, 'Task 5. Finite Difference Implicit Heat Transfer method. Central difference (h is fixed):'
    print *, '  N             τ                   s(t=t1)                   s(t=t2)                 ' // &
            'max(t=t1)                 max(t=t2)'
    print *, '---  ------------  ------------------------  ------------------------  ' // &
            '------------------------  ------------------------'

    n = 10; m = 20
    allocate(errors_t1(5), errors_t2(5), dt_array(5))
    do i = 1, 5
        m = m * 2
        m2 = m * 2

        allocate(u(m+1,n+1), x_list(n+1), t_list(m+1), u2(m2+1,n+1), x_list2(n+1), t_list2(m2+1), err_t1(n+1), err_t2(n+1))

        call fd_heat_transfer_implicit(n, m, a, b, t, g1, g2, phi, k(1._dp), f, x_list, t_list, u, bc=2, &
                is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Implicit Heat Transfer method. Central difference'))

        call fd_heat_transfer_implicit(n, m2, a, b, t, g1, g2, phi, k(1._dp), f, x_list2, t_list2, u2, bc=2, &
                is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Implicit Heat Transfer method. Central difference'))

        dt = t / m; dt_array(i) = dt
        err_t1 = u(int(t1/dt),:) - u2(2*int(t1/dt),:); errors_t1(i) = norm2(err_t1)
        err_t2 = u(int(t2/dt),:) - u2(2*int(t2/dt),:); errors_t2(i) = norm2(err_t2)

        print '(i4,a2,f12.7,a2,f24.16,a2,f24.16,a2,f24.16,a2,f24.16)',n,'  ',dt,'  ',std(err_t1),'  ',std(err_t2), &
                '  ',maxval(abs(err_t1)),'  ',maxval(abs(err_t2))
        deallocate(u, x_list, t_list, u2, x_list2, t_list2, err_t1, err_t2)
    end do

    call gp%title("Task 5. Finite Difference Implicit Heat Transfer method. Central difference (h is fixed)")
    call gp%options('set key top left; set grid')
    call gp%xlabel('τ')
    call gp%ylabel('error')
    call gp%plot(x1=dt_array, y1=errors_t1, ls1='with lines title "t=t1"', &
            x2=dt_array, y2=errors_t2, ls2='with lines title "t=t2"')

    deallocate(errors_t1, errors_t2, dt_array)

    ! Finite Difference Implicit Heat Transfer method -----------------------------------------------------------------
    ! Неявная схема. Центральная разность (τ - фиксированное) ---------------------------------------------------------

    print *, ''
    print *, 'Task 5. Finite Difference Implicit Heat Transfer method. Central difference (τ is fixed):'
    print *, '  N             τ                   s(t=t1)                   s(t=t2)                 ' // &
            'max(t=t1)                 max(t=t2)'
    print *, '---  ------------  ------------------------  ------------------------  ' // &
            '------------------------  ------------------------'

    n = 3; m = 4000
    allocate(errors_t1(5), errors_t2(5), dx_array(5))
    do i = 1, 5
        n = n * 2
        n2 = n * 2

        allocate(u(m+1,n+1), x_list(n+1), t_list(m+1), u2(m+1,n2+1), x_list2(n2+1), t_list2(m+1), err_t1(n+1), err_t2(n+1))

        call fd_heat_transfer_implicit(n, m, a, b, t, g1, g2, phi, k(1._dp), f, x_list, t_list, u, bc=2, &
                is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Implicit Heat Transfer method. Central difference'))

        call fd_heat_transfer_implicit(n2, m, a, b, t, g1, g2, phi, k(1._dp), f, x_list2, t_list2, u2, bc=2, &
                is_print=.false., is_draw=.false., &
                info=('Task 5. Finite Difference Implicit Heat Transfer method. Central difference'))

        dx = (b - a) / n; dt = t / m; dx_array(i) = dx
        err_t1 = u(int(t1/dt),:) - u2(int(t1/dt),1:n2+1:2); errors_t1(i) = norm2(err_t1)
        err_t2 = u(int(t2/dt),:) - u2(int(t2/dt),1:n2+1:2); errors_t2(i) = norm2(err_t2)

        print '(i4,a2,f12.7,a2,f24.16,a2,f24.16,a2,f24.16,a2,f24.16)',n,'  ',dt,'  ',std(err_t1),'  ',std(err_t2), &
                '  ',maxval(abs(err_t1)),'  ',maxval(abs(err_t2))
        deallocate(u, x_list, t_list, u2, x_list2, t_list2, err_t1, err_t2)
    end do

    call gp%title("Task 5. Finite Difference Implicit Heat Transfer method. Central difference (τ is fixed)")
    call gp%options('set key top left; set grid')
    call gp%xlabel('h')
    call gp%ylabel('error')
    call gp%plot(x1=dx_array, y1=errors_t1, ls1='with lines title "t=t1"', &
            x2=dx_array, y2=errors_t2, ls2='with lines title "t=t2"')

    deallocate(errors_t1, errors_t2, dx_array)

contains

    real(dp) function f(x,t)

        implicit none
        real(dp), intent(in) :: x, t

        f = x

    end function f

    real(dp) function k(x)

        implicit none
        real(dp), intent(in) :: x

        k = 0.5_dp

    end function k

    real(dp) function g1(t)

        implicit none
        real(dp), intent(in) :: t

        g1 = 1._dp

    end function g1

    real(dp) function g2(t)

        implicit none
        real(dp), intent(in) :: t

        g2 = 1._dp

    end function g2

    real(dp) function phi(x)

        implicit none
        real(dp), intent(in) :: x

        phi = x * x

    end function phi

end program