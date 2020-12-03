! Решение краевых задач. Методы коллокаций, наименьших квадратов и Галеркина.

! -Основное задание (по умолчанию k = 12):
!  Методами коллокаций, галеркина, интегральным и дискретным методами наименьших квадратов
!  получить численное решение краевой задачи:
!  a * y'' + (1 + b * x^2) * y = -1,  -1 <= x <= 1, a = sin(k), b = cos(k), y(-1) = 0, y(1) = 0.
!  Базисная система вида:
!  φ0 = 0, φi(x) = x^i * (1 - x^2), i = 1, 2, ...

! -Тестовый пример 2.1(2.2)
!  Схож с основным заданием, но a = 1, b = 1 и
!  Базисная система вида:
!  φ0 = 0, φ1(x) = 1 - x^2, φi(x) = x^i * (1 - x^2), i = 2, 3, ...

! -Тестовый пример 2.3
!  Краевая задача вида:
!  y'' + y = x,  0 <= x <= 1, y(0) = 0, y(1) = 0.
!  Базисная система вида:
!  φ0 = 0, φi(x) = x^i * (1 - x), i = 1, 2, ...


program main

    use ogpf
    use base
    use collocation
    use galerkin
    use least_squares

    ! Ввод начальных данных -------------------------------------------------------------------------------------------
    implicit none

    type(gpf) :: gp
    integer, parameter :: dp = kind(0.d0), N_TEST = 3, X_POINTS = 100

    character *16 :: info
    real(dp), parameter :: X_MIN = -1._dp, X_MAX = 1._dp, K = 12._dp
    real(dp) :: a, b

    real(dp), dimension(:), allocatable :: sol1, sol2, sol3, sol4

    real(dp), dimension(1:N_TEST) :: test_sol1, test_sol2, test_sol3, test_sol4
    real(dp), dimension(1:X_POINTS) :: sol1_ext, tsol1_ext, sol2_ext, tsol2_ext, &
            sol3_ext, tsol3_ext, sol4_ext, tsol4_ext, x_list

    integer :: n, i, j

    print '(a9)', "Input n: "
    read *, n

    allocate(sol1(1:n), sol2(1:n), sol3(1:n), sol4(1:n))

    call noised_linspace(X_MIN, X_MAX, x_list, noise=0._dp)

    ! Метод коллокаций ------------------------------------------------------------------------------------------------
    print "(/, a21)", "1. Collocation method"

    info = "Test sample 2.1."
    a = 1._dp; b = 1._dp
    print "(/, a16)", info
    call collocation_method(a, b, X_MIN, X_MAX, N_TEST, f1, test_sol1, tsol1_ext, x_points=X_POINTS, &
        is_test=.true., is_print=.true., is_draw=.true., info=info)

    info = "General sample."
    a = SIN(K); b = COS(K)
    print "(/, a15)", info
    call collocation_method(a, b, X_MIN, X_MAX, n, f1, sol1, sol1_ext, x_points=X_POINTS, &
            is_test=.true., is_print=.true., is_draw=.true., info=info)


    ! Интегральный метод наименьших квадратов -------------------------------------------------------------------------
    print "(/, a34)", "2.1. Integral least squares method"

    info = "Test sample 2.2."
    a = 1._dp; b = 1._dp
    print "(/, a16)", info
    call int_least_squares_method(a, b, X_MIN, X_MAX, N_TEST, f1, test_sol2, tsol2_ext, x_points=X_POINTS, &
            is_test=.true., is_print=.true., is_draw=.true., info=info)

    info = "General sample."
    a = SIN(K); b = COS(K)
    print "(/, a15)", info
    call int_least_squares_method(a, b, X_MIN, X_MAX, n, f1, sol2, sol2_ext, x_points=X_POINTS, &
            is_test=.true., is_print=.true., is_draw=.true., info=info)

    ! Дискретный метод наименьших квадратов ---------------------------------------------------------------------------
    print "(/, a34)", "2.2. Discrete least squares method"

    info = "Test sample 2.2."
    a = 1._dp; b = 1._dp
    print "(/, a16)", info
    call disc_least_squares_method(a, b, X_MIN, X_MAX, N_TEST, f1, test_sol3, tsol3_ext, x_points=X_POINTS, &
            is_test=.true., is_print=.true., is_draw=.true., info=info)

    info = "General sample."
    a = SIN(K); b = COS(K)
    print "(/, a15)", info
    call disc_least_squares_method(a, b, X_MIN, X_MAX, n, f1, sol3, sol3_ext, x_points=X_POINTS, &
            is_test=.true., is_print=.true., is_draw=.true., info=info)

    ! Метод Галеркина -------------------------------------------------------------------------------------------------
    print "(/, a18)", "3. Galerkin method"

! для правильного решения тестового задания в сабрутине замените test_basic_func на test_galerkin_basic_func

!    info = "Test sample 2.3."
!    a = 1._dp; b = 0
!    print "(/, a16)", info
!    call galerkin_method(a, b, 0.0000001_dp, X_MAX, N_TEST, f2, test_sol4, tsol4_ext, X_POINTS, &
!            is_test=.true., is_print=.true., is_draw=.true., info=info)

    info = "General sample."
    a = SIN(K); b = COS(K)
    print "(/, a15)", info
    call galerkin_method(a, b, X_MIN, X_MAX, n, f1, sol4, sol4_ext, X_POINTS, .true., .true., .true., info)

    ! Общий график ----------------------------------------------------------------------------------------------------
    call gp%title('General sample')
    call gp%options('set key top left; set grid')
    call gp%xlabel('-1 <= X <= 1')
    call gp%ylabel('Yn(x) = φ0(x) + a1 * φ1(x) + ... + an * φn(x)')
    call gp%plot(x1=x_list, y1=sol1_ext, ls1='with lines title "Collocation method"', &
            x2=x_list, y2=sol2_ext, ls2='with lines title "Integral least squares method"', &
            x3=x_list, y3=sol3_ext, ls3='with lines title "Discrete least squares method"', &
            x4=x_list, y4=sol4_ext, ls4='with lines title "Galerkin method"')
            
    deallocate(sol1, sol2, sol3, sol4)

end program

