program task3

    use ogpf
    use base
    use heat_transfer

    ! Ввод начальных данных -------------------------------------------------------------------------------------------

    implicit none

    type(gpf) :: gp
    character(len=21) :: phi_str
    integer, parameter :: dp = kind(0.d0), n_phi = 3
    real(dp), parameter :: a = 1._dp, b = 2._dp, t = 0.5_dp

    interface
        real(kind=8) function func(x)
            implicit none
            real(kind=8), intent (in) :: x

        end function func
    end interface

    type pp
        procedure(func), pointer, nopass :: f => null()
    end type pp

    type(pp) :: phi_variants(n_phi)
    real(dp), allocatable :: x_list(:), t_list(:), u_variants(:,:,:), u_steady(:), x_steady(:)

    integer :: i, n, m

    phi_variants(1)%f => phi1
    phi_variants(2)%f => phi2
    phi_variants(3)%f => phi3

    n = 10; m = 200

    ! Finite Difference Centered Heat Transfer Steady method ----------------------------------------------------------
    ! Пример (не зависящий от времени) с аналогичными исходными данными -----------------------------------------------

    allocate(u_steady(1:n+1), x_steady(1:n+1))

    call fd_centered_heat_transfer_steady(n, a, b, 3._dp, 3._dp, k, f1, x_steady, &
            u_steady, .false., .true., 'Task 3. Finite Difference Centered Heat Transfer Steady method')

    deallocate(u_steady, x_steady)

    ! Finite Difference Explicit Heat Transfer method -----------------------------------------------------------------

    allocate(u_variants(n_phi,m+1,n+1), x_list(n+1), t_list(m+1))

    do i = 1, n_phi
        if (i == 1) phi_str = " φ = 2 * x"
        if (i == 2) phi_str = " φ = sin(x)"
        if (i == 3) phi_str = " φ = 10 * (x ** 1/3)"

        call fd_heat_transfer_explicit(n=n, m=m, a=a, b=b, t=t, g1=g1, g2=g2, phi=phi_variants(i)%f, k=k, f=f, &
                x_list=x_list, t_list=t_list, matrix=u_variants(i,:,:), bc=0, is_check_cfl=.false., &
                is_print=.false., is_draw=.true., &
                info=('Task 3. Finite Difference Explicit Heat Transfer method:' //  phi_str))
    end do

    deallocate(u_variants, x_list, t_list)

contains

    real(dp) function f(x,t)

        implicit none
        real(dp), intent(in) :: x, t

        f = 1._dp / x
        f = f * (1._dp - exp(-t))

    end function f

    real(dp) function f1(x)

        implicit none
        real(dp), intent(in) :: x

        f1 = 1._dp / x

    end function f1

    real(dp) function k(x)

        implicit none
        real(dp), intent(in) :: x

        k = x

    end function k

    real(dp) function g1(t)

        implicit none
        real(dp), intent(in) :: t

        g1 = 3._dp

    end function g1

    real(dp) function g2(t)

        implicit none
        real(dp), intent(in) :: t

        g2 = 3._dp

    end function g2

    real(dp) function phi1(x)

        implicit none
        real(dp), intent(in) :: x

        phi1 = 2 * x

    end function phi1

    real(dp) function phi2(x)

        implicit none
        real(dp), intent(in) :: x

        phi2 = sin(x)

    end function phi2

    real(dp) function phi3(x)

        implicit none
        real(dp), intent(in) :: x

        phi3 = 10 * (x ** 1/3)

    end function phi3

end program


