program task4

    use ogpf
    use base
    use heat_transfer

    ! Ввод начальных данных -------------------------------------------------------------------------------------------

    implicit none

    type(gpf) :: gp
    integer, parameter :: dp = kind(0.d0)
    real(dp), parameter :: a = -1._dp, b = 1._dp, t = 0.4_dp

    real(dp), allocatable :: x_list(:), t_list(:), u(:,:)

    integer :: i, n, m


    n = 10; m = 100

    ! Finite Difference Explicit Heat Transfer method -----------------------------------------------------------------

    allocate(u(m+1,n+1), x_list(n+1), t_list(m+1))


    call fd_heat_transfer_explicit(n, m, a, b, t, g1, g2, phi, k, f, x_list, t_list, u, bc=0, &
            is_check_cfl=.true., is_print=.true., is_draw=.true., &
            info=('Task 4. Finite Difference Explicit Heat Transfer method'))

    deallocate(u, x_list, t_list)

contains

    real(dp) function f(x,t)

        implicit none
        real(dp), intent(in) :: x, t

        f = 0._dp

    end function f

    real(dp) function k(x)

        implicit none
        real(dp), intent(in) :: x

        k = 2._dp

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


