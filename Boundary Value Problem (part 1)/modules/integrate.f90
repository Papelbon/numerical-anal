module integrate
    implicit none

    private

    public trapezium, simpson

    integer, parameter:: dp = kind(0.d0)
contains
    ! -----------------------------------------------------------------------------------------------------------------
    ! Численное интегрирование методом трапеций
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function trapezium(f, x0, x1, i, j, a, b, fm_val, div_no, is_test) result(sum)

        implicit none

        real(dp), external :: f
        real(dp), external :: fm_val

        real(dp), intent(in) :: x0, x1, a, b
        integer, intent(in) :: div_no, i, j
        logical, intent(in) :: is_test

        real(dp) :: x, dx
        integer :: l

        dx = (x1 - x0) / div_no
        sum = f(x0, i, j, a, b, fm_val, is_test) + f(x1, i, j, a, b, fm_val, is_test)
        x = x0

        do l = 1, div_no - 1
            x = x + dx
            sum = sum + 2.0 * f(x, i, j, a, b, fm_val, is_test)
        end do
        sum = dx * sum / 2.0

    end function trapezium
    ! -----------------------------------------------------------------------------------------------------------------
    ! Численное интегрирование методом Симпсона
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function simpson(f, x0, x1, i, j, a, b, fm_val, div_no, is_test) result(sum)

        implicit none

        real(dp), external :: f
        real(dp), external :: fm_val

        real(dp), intent(in) :: x0, x1, a, b
        integer, intent(in) :: div_no, i, j
        logical, intent(in) :: is_test

        real(dp) :: x, dx

        integer :: l

        dx = (x1 - x0) / (2.0 * div_no)
        sum = f(x0, i, j, a, b, fm_val, is_test) + f(x1, i, j, a, b, fm_val, is_test)
        x = x0

        do l = 1, 2 * div_no - 1
            x = x + dx

            if(mod(l, 2) /= 0) then
                sum = sum + 4.0 * f(x, i, j, a, b, fm_val, is_test)
            else
                sum = sum + 2.0 * f(x, i, j, a, b, fm_val, is_test)
            end if

        end do

        sum = dx * sum / 3.0

    end function simpson
end module integrate