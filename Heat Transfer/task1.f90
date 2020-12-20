! ---------------------------------------------------------------------------------------------------------------------
! Задача 1.
!
! Решить стационарное (не зависящее от времени) уравнение теплопроводности в одномерной области:
! -d/dx (k(x) du/dx) = f(x),
! u(a) = Ua, u(b) = Ub.
!
! Исходные наборы данных:
! k(x) = x, a = 1, b = 2, Ua = 3, Ub = 3.
!
! -------------------------------------------------------------------------------------------------------------
! параметры |   набор 1   |   набор 2   |   набор 3   |   набор 4   |   набор 5   |   набор 6   |   набор 7   |
! -------------------------------------------------------------------------------------------------------------
!    k(x)   |     k(x)    |   2 * k(x)  | 0.1 * k(x)  |  1 / k(x)   |    k(x)     |    k(x)     |    k(x)     |
! -------------------------------------------------------------------------------------------------------------
!     Ua    |     Ua      |     Ua      |     Ua      |      Ua     |    -Ua      |    Ua       |    -Ua      |
! -------------------------------------------------------------------------------------------------------------
!     Ub    |     Ub      |     Ub      |     Ub      |      Ub     |     Ub      |    -Ua      |    -Ua      |
! -------------------------------------------------------------------------------------------------------------
!
! ---------------------------------------------------------------------------------------------------------------------

program task1

    use ogpf
    use base
    use heat_transfer


    ! Ввод начальных данных -------------------------------------------------------------------------------------------

    implicit none

    type(gpf) :: gp
    integer, parameter :: dp = kind(0.d0), nk = 7
    real(dp), parameter :: a = 1._dp, b = 2._dp, Ua = 3._dp, Ub = 3._dp

    interface
        real(kind=8) function func(x)
            implicit none
            real(kind=8), intent (in) :: x

        end function func
    end interface

    type pp
        procedure(func), pointer, nopass :: f => null()
    end type pp

    type(pp) :: k_variants(nk)
    real(dp), dimension(1:nk) :: ua_array = (/Ua, Ua, Ua, Ua, -Ua, Ua, -Ua/), &
            ub_array = (/Ub, Ub, Ub, Ub, Ub, -Ub, -Ub/)
    real(dp), allocatable :: x_list(:), u_variants(:,:)
    integer :: i, n

    k_variants(1)%f => k
    k_variants(2)%f => k2
    k_variants(3)%f => k3
    k_variants(4)%f => k4
    k_variants(5)%f => k
    k_variants(6)%f => k
    k_variants(7)%f => k

    n = 150
    allocate(u_variants(1:nk,1:n+1), x_list(1:n+1))

    ! Finite Difference Centered Heat Transfer Steady method ----------------------------------------------------------
    ! Наборы 1-3 ------------------------------------------------------------------------------------------------------
    ! Наборы 1 и 4 ----------------------------------------------------------------------------------------------------
    ! Наборы 5-7 ------------------------------------------------------------------------------------------------------

    do i = 1, nk
        call fd_centered_heat_transfer_steady(n, a, b, ua_array(i), ub_array(i), k_variants(i)%f, f, x_list, &
                u_variants(i,1:n+1), .false., .false., 'Task 1. Finite Difference Centered Heat Transfer Steady method')
    end do

    ! 1-3
    call gp%title("Task 1. Finite Difference Centered Heat Transfer Steady method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('a <= X <= b')
    call gp%ylabel('Temperature')
    call gp%plot(x1=x_list, y1=u_variants(1,1:n+1), ls1='with lines title "set 1"', &
            x2=x_list, y2=u_variants(2,1:n+1), ls2='with lines title "set 2"',      &
            x3=x_list, y3=u_variants(3,1:n+1), ls3='with lines title "set 3"')

    ! 1 и 4
    call gp%title("Task 1. Finite Difference Centered Heat Transfer Steady method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('a <= X <= b')
    call gp%ylabel('Temperature')
    call gp%plot(x1=x_list, y1=u_variants(1,1:n+1), ls1='with lines title "set 1"', &
            x2=x_list, y2=u_variants(4,1:n+1), ls2='with lines title "set 4"')

    ! 5-7
    call gp%title("Task 1. Finite Difference Centered Heat Transfer Steady method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('a <= X <= b')
    call gp%ylabel('Temperature')
    call gp%plot(x1=x_list, y1=u_variants(5,1:n+1), ls1='with lines title "set 5"', &
            x2=x_list, y2=u_variants(6,1:n+1), ls2='with lines title "set 6"',      &
            x3=x_list, y3=u_variants(7,1:n+1), ls3='with lines title "set 7"')

    ! Finite Difference Balance Heat Transfer Steady method -----------------------------------------------------------
    ! Наборы 1-3 ------------------------------------------------------------------------------------------------------
    ! Наборы 1 и 4 ----------------------------------------------------------------------------------------------------
    ! Наборы 5-7 ------------------------------------------------------------------------------------------------------

    do i = 1, nk
        call fd_balance_heat_transfer_steady(n, a, b, ua_array(i), ub_array(i), k_variants(i)%f, f, x_list, &
                u_variants(i,1:n+1), .false., .false., 'Task 1. Finite Difference Balance Heat Transfer Steady method')
    end do

    ! 1-3
    call gp%title("Task 1. Finite Difference Balance Heat Transfer Steady method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('a <= X <= b')
    call gp%ylabel('Temperature')
    call gp%plot(x1=x_list, y1=u_variants(1,1:n+1), ls1='with lines title "set 1"', &
            x2=x_list, y2=u_variants(2,1:n+1), ls2='with lines title "set 2"', &
            x3=x_list, y3=u_variants(3,1:n+1), ls3='with lines title "set 3"')

    ! 1 и 4
    call gp%title("Task 1. Finite Difference Balance Heat Transfer Steady method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('a <= X <= b')
    call gp%ylabel('Temperature')
    call gp%plot(x1=x_list, y1=u_variants(1,1:n+1), ls1='with lines title "set 1"', &
            x2=x_list, y2=u_variants(4,1:n+1), ls2='with lines title "set 4"')

    ! 5-7
    call gp%title("Task 1. Finite Difference Balance Heat Transfer Steady method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('a <= X <= b')
    call gp%ylabel('Temperature')
    call gp%plot(x1=x_list, y1=u_variants(5,1:n+1), ls1='with lines title "set 5"', &
            x2=x_list, y2=u_variants(6,1:n+1), ls2='with lines title "set 6"', &
            x3=x_list, y3=u_variants(7,1:n+1), ls3='with lines title "set 7"')

    deallocate(u_variants, x_list)

contains
    real(dp) function f(x)

        implicit none
        real(dp), intent(in) :: x

        f = 1 / x

    end function f

    real(dp) function k(x)

        implicit none
        real(dp), intent(in) :: x

        k = x

    end function k

    real(dp) function k2(x)

        implicit none
        real(dp), intent(in) :: x

        k2 = 2 * x

    end function k2

    real(dp) function k3(x)

        implicit none
        real(dp), intent(in) :: x

        k3 = 0.1 * x

    end function k3

    real(dp) function k4(x)

        implicit none
        real(dp), intent(in) :: x

        k4 = 1 / x

    end function k4


end program


