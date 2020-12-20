! ---------------------------------------------------------------------------------------------------------------------
! Задача 2.
!
! Решить стационарное (не зависящее от времени) уравнение теплопроводности в одномерной области:
! -d/dx (k(x) du/dx) = f(x),
! u(a) = Ua, u(b) = Ub.
!
! Исходные наборы данных:
! a = 1, b = 2, Ua = 3, Ub = 3.
!
! a) стержень состоит из двух материалов с различными коэффициентами теплопроводности:
!
!           /  k1, a <= x <= (b + a) / 2,
! k(x) = --|                                    при a) k1 << k2; b) k2 << k1.
!           \  k2, (b + a) / 2 < x <= b
!
! b) стержень состоит из трех материалов с различными свойствами:
!
!           /  k1, a <= x <= a + (b - a) / 3,
! k(x) = --|   k2, a + (b - a) / 3 < x <= a + 2 * (b - a) / 3,
!           \  k3, a + 2 * (b - a) / 3 < x <= b
!
!  при a) k1 < k2 < k3;
!      b) k1 > k2 > k3;
!      c) k1 = k, k2 = 2k, k3 = k;
!      d) k1 = 20k, k2 = k, k3 = 20k;
!
! c) f(x) = c * δ(x - x0)     --    точечный источник тепла,
!
! c    - константа (мощность источника),
! δ(x) - дельта функция,
! x0   - точка из [a, b], где располагается точечный источник.
!
! рассмотреть a) x0 - середина отрезка [a,b];
!             b) два одинаковых по мощности источника поставлены в разные точки
!                отрезка, симметричные относительно середины отрезка;
!             c) два различных по мощности источника поставлены симметрично;
!             d) ...
! ---------------------------------------------------------------------------------------------------------------------

program task2

    use ogpf
    use base
    use heat_transfer


    ! Ввод начальных данных -------------------------------------------------------------------------------------------

    implicit none

    type(gpf) :: gp
    integer, parameter :: dp = kind(0.d0), nk_a = 2, nk_b = 4, nf = 4
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

    type(pp) :: k_variants_a(nk_a), k_variants_b(nk_b)
    real(dp), dimension(1:nf, 1:4) :: f_args
    real(dp), allocatable :: x_list(:), u_variants(:,:)

    integer :: i, n

    ! Task 2(a) -------------------------------------------------------------------------------------------------------

    k_variants_a(1)%f => ka1
    k_variants_a(2)%f => ka2

    ! Task 2(b) -------------------------------------------------------------------------------------------------------

    k_variants_b(1)%f => kb1
    k_variants_b(2)%f => kb2
    k_variants_b(3)%f => kb3
    k_variants_b(4)%f => kb4

    ! Task 2(c) -------------------------------------------------------------------------------------------------------

    f_args = 0._dp
    ! x0_1, c1 - источник в середине отрезка
    f_args(1,1) = a + (b - a) / 2; f_args(1,2) = 2._dp
    ! x0_1, c1, x0_2, c2 - два разных источника - одинаковой мощности - симметрично
    f_args(2,1) = a + (b - a) / 4; f_args(2,2) = 4._dp; f_args(2,3) = a + 3 * (b - a) / 4; f_args(2,4) = 4._dp
    ! x0_1, c1, x0_2, c2 - два разных источника - разной мощности - симметрично
    f_args(3,1) = a + (b - a) / 4; f_args(3,2) = 2._dp; f_args(3,3) = a + 3 * (b - a) / 4; f_args(3,4) = 5.5_dp
    ! x0_1, c1, x0_2, c2 - два разных источника - одинаковой мощности - несимметрично
    f_args(4,1) = a + (b - a) / 2; f_args(4,2) = 7._dp; f_args(4,3) = b - (b - a) / 3;     f_args(4,4) = 7._dp


    n = 150

    ! Finite Difference Balance Heat Transfer Steady method -----------------------------------------------------------
    ! Task 2(a) - стержень состоит из двух материалов с различными коэффициентами теплопроводности --------------------

    allocate(u_variants(1:nk_a,1:n+1), x_list(1:n+1))

    do i = 1, nk_a
        call fd_balance_heat_transfer_steady(n, a, b, Ua, Ub, k_variants_a(i)%f, f, x_list, u_variants(i,1:n+1), &
                .false., .false., 'Task 2(a). Finite Difference Balance Heat Transfer Steady method')
    end do

    call gp%title("Task 2(a). Finite Difference Balance Heat Transfer Steady method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('a <= X <= b')
    call gp%ylabel('Temperature')
    call gp%plot(x1=x_list, y1=u_variants(1,1:n+1), ls1='with lines title "set a"', &
            x2=x_list, y2=u_variants(2,1:n+1), ls2='with lines title "set b"')

    deallocate(u_variants, x_list)

    ! Finite Difference Balance Heat Transfer Steady method -----------------------------------------------------------
    ! Task 2(b) - стержень состоит из трех материалов с различными свойствами -----------------------------------------

    allocate(u_variants(1:nk_b,1:n+1), x_list(1:n+1))

    do i = 1, nk_b
        call fd_balance_heat_transfer_steady(n, a, b, Ua, Ub, k_variants_b(i)%f, f, x_list, u_variants(i,1:n+1), &
                .false., .false., 'Task 2(b). Finite Difference Balance Heat Transfer Steady method')
    end do

    call gp%title("Task 2(b). Finite Difference Balance Heat Transfer Steady method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('a <= X <= b')
    call gp%ylabel('Temperature')
    call gp%plot(x1=x_list, y1=u_variants(1,1:n+1), ls1='with lines title "set a"', &
            x2=x_list, y2=u_variants(2,1:n+1), ls2='with lines title "set b"', &
            x3=x_list, y3=u_variants(3,1:n+1), ls3='with lines title "set c"', &
            x4=x_list, y4=u_variants(4,1:n+1), ls4='with lines title "set d"')

    deallocate(u_variants, x_list)

    ! Finite Difference Balance Heat Transfer Steady method -----------------------------------------------------------
    ! Task 2(c) - точечный источник тепла -----------------------------------------------------------------------------

    allocate(u_variants(1:nf,1:n+1), x_list(1:n+1))

    do i = 1, nf
        call fd_balance_heat_transfer_steady(n, a, b, Ua, Ub, k, f, x_list, u_variants(i,1:n+1), .false., .false., &
            'Task 2(c). Finite Difference Balance Heat Transfer Steady method', f_args(i,1:4))
    end do

    call gp%title("Task 2(c). Finite Difference Balance Heat Transfer Steady method")
    call gp%options('set key top left; set grid')
    call gp%xlabel('a <= X <= b')
    call gp%ylabel('Temperature')
    call gp%plot(x1=x_list, y1=u_variants(1,1:n+1), ls1='with lines title "set a"', &
            x2=x_list, y2=u_variants(2,1:n+1), ls2='with lines title "set b"', &
            x3=x_list, y3=u_variants(3,1:n+1), ls3='with lines title "set c"', &
            x4=x_list, y4=u_variants(4,1:n+1), ls4='with lines title "set d"')

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

    real(dp) function ka1(x)

        implicit none
        real(dp), intent(in) :: x
        real(dp), parameter :: k1 = 2._dp, k2 = 100._dp

        if (a <= x .and. x <= (b + a) / 2) then
            ka1 = k1
        else
            ka1 = k2
        end if

    end function ka1

    real(dp) function ka2(x)

        implicit none
        real(dp), intent(in) :: x
        real(dp), parameter :: k1 = 100._dp, k2 = 2._dp

        if (a <= x .and. x <= (b + a) / 2) then
            ka2 = k1
        else
            ka2 = k2
        end if

    end function ka2

    real(dp) function kb1(x)

        implicit none
        real(dp), intent(in) :: x
        real(dp), parameter :: k1 = 2._dp, k2 = 15._dp, k3 = 30._dp

        if (a <= x .and. x <= a + (b - a) / 3) then
            kb1 = k1
        elseif (a + (b - a) / 3 < x .and. x <= a + 2 * (b - a) / 3) then
            kb1 = k2
        else
            kb1 = k3
        end if

    end function kb1

    real(dp) function kb2(x)

        implicit none
        real(dp), intent(in) :: x
        real(dp), parameter :: k1 = 30._dp, k2 = 15._dp, k3 = 2._dp

        if (a <= x .and. x <= a + (b - a) / 3) then
            kb2 = k1
        elseif (a + (b - a) / 3 < x .and. x <= a + 2 * (b - a) / 3) then
            kb2 = k2
        else
            kb2 = k3
        end if

    end function kb2

    real(dp) function kb3(x)

        implicit none
        real(dp), intent(in) :: x
        real(dp), parameter :: k1 = 5._dp, k2 = 10._dp, k3 = 5._dp

        if (a <= x .and. x <= a + (b - a) / 3) then
            kb3 = k1
        elseif (a + (b - a) / 3 < x .and. x <= a + 2 * (b - a) / 3) then
            kb3 = k2
        else
            kb3 = k3
        end if

    end function kb3

    real(dp) function kb4(x)

        implicit none
        real(dp), intent(in) :: x
        real(dp), parameter :: k1 = 100._dp, k2 = 5._dp, k3 = 100._dp

        if (a <= x .and. x <= a + (b - a) / 3) then
            kb4 = k1
        elseif (a + (b - a) / 3 < x .and. x <= a + 2 * (b - a) / 3) then
            kb4 = k2
        else
            kb4 = k3
        end if

    end function kb4

end program


