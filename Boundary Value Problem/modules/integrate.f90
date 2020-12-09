module integrate
    implicit none

    private

    public trapezium, simpson

    integer, parameter:: dp = kind(0.d0)
contains
    ! -----------------------------------------------------------------------------------------------------------------
    ! Численное интегрирование методом трапеций
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function trapezium(func, x0, x1, i, j, f, a, p, q, basic, div_no) result(sum)

        implicit none

        real(dp), external :: func, f, p, q, basic, a
        real(dp), intent(in) :: x0, x1
        integer, intent(in) :: div_no, i, j

        real(dp) :: x, dx
        integer :: l

        dx = (x1 - x0) / div_no
        sum = func(x0, i, j, f, a, p, q, basic) + func(x1, i, j, f, a, p, q, basic)
        x = x0

        do l = 1, div_no - 1
            x = x + dx
            sum = sum + 2.0 * func(x, i, j, f, a, p, q, basic)
        end do
        sum = dx * sum / 2.0

    end function trapezium
    ! -----------------------------------------------------------------------------------------------------------------
    ! Численное интегрирование методом Симпсона
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function simpson(func, x0, x1, i, j, f, a, p, q, basic, div_no) result(sum)

        implicit none

        real(dp), external :: func, f, p, q, basic, a
        real(dp), intent(in) :: x0, x1
        integer, intent(in) :: div_no, i, j

        real(dp) :: x, dx
        integer :: l

        dx = (x1 - x0) / (2.0 * div_no)
        sum = func(x0, i, j, f, a, p, q, basic)  + func(x1, i, j, f, a, p, q, basic)
        x = x0

        do l = 1, 2 * div_no - 1
            x = x + dx

            if(mod(l, 2) /= 0) then
                sum = sum + 4.0 * func(x, i, j, f, a, p, q, basic)
            else
                sum = sum + 2.0 * func(x, i, j, f, a, p, q, basic)
            end if

        end do

        sum = dx * sum / 3.0

    end function simpson
end module integrate