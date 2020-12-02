! Модуль содержит основные функции и сабрутины

module base

    implicit none

    private

    public basic_func, test_basic_func, test_galerkin_basic_func, noised_linspace, l2_norm, residual2, f1, f2

    integer, parameter:: dp = kind(0.d0)

contains
    ! -----------------------------------------------------------------------------------------------------------------
    ! Базисная функция для основного задания
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function basic_func(x, i, deriv_lvl) result(basic)

        implicit none

        real(dp), intent(in) :: x
        integer, intent(in) :: i, deriv_lvl

        if (i == 0) then
            basic = 0._dp
        else
            select case(deriv_lvl)
            case(0)
                basic = x ** i * (1 - x * x)
            case(1)
                basic = x ** (i - 1) * (i - x * (i + 1))
            case(2)
                basic = x ** (i - 2) * ((i * i - i) - x * x * (i + 1) * (i + 2))
            end select
        end if

    end function basic_func
    ! -----------------------------------------------------------------------------------------------------------------
    ! Базисная функция для тестовых заданий
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function test_basic_func(x, i, deriv_lvl) result(basic)

        implicit none

        real(dp), intent(in) :: x
        integer, intent(in) :: i, deriv_lvl

        if (i == 0) then
            basic = 0._dp
        else if (i == 1) then
            select case(deriv_lvl)
            case(0)
                basic = 1 - x * x
            case(1)
                basic = -2 * x
            case(2)
                basic = -2
            end select
        else
            select case(deriv_lvl)
            case(0)
                basic = x ** i * (1 - x * x)
            case(1)
                basic = x ** (i - 1) * (i - x * (i + 1))
            case(2)
                basic = x ** (i - 2) * ((i * i - i) - x * x * (i + 1) * (i + 2))
            end select
        end if

    end function test_basic_func
    ! -----------------------------------------------------------------------------------------------------------------
    ! Базисная функция для тестового задания для метода Галеркина
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function test_galerkin_basic_func(x, i, deriv_lvl) result(basic)

        implicit none

        real(dp), intent(in) :: x
        integer, intent(in) :: i, deriv_lvl

        if (i == 0) then
            basic = 0._dp
        else
            select case(deriv_lvl)
            case(0)
                basic = x ** i * (1 - x)
            case(1)
                basic = i * (1 - x) * (x ** (i - 1)) - (x ** i)
            case(2)
                basic = (i - 1) * i * (1 - x) * (x ** (i - 2)) - (2 * i * x ** (i - 1))
            end select
        end if

    end function test_galerkin_basic_func
    ! -----------------------------------------------------------------------------------------------------------------
    ! Равномерное разбиение с учетом крайних шумов
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine noised_linspace(from, to, array, noise)

        implicit none

        real(dp), intent(in) :: from, to, noise
        real(dp), intent(out) :: array(:)

        real(dp) :: range
        integer :: n, i

        n = size(array)
        range = (to-noise) - (from+noise)

        if (n == 0) return

        if (n == 1) then
            array(1) = from+noise
            return
        end if

        do i = 1, n
            array(i) = (from+noise) + range * (i - 1) / (n - 1)

            ! во избежание деления на 0
            if (array(i) == 0) then
                array(i) = 1e-8
            end if

        end do

        return

    end
    ! -----------------------------------------------------------------------------------------------------------------
    ! L2 норма (погрешность, для вычисления интеграла используется метод трапеций)
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function l2_norm(f, x0, x1, n, a, b, a_sol, fm_val, div_no, basic) result(l2)

        implicit none

        real(dp), external :: f, basic, fm_val

        real(dp), intent(in) :: x0, x1, a, b, a_sol(1:n)
        integer, intent(in) :: div_no, n

        real(dp) :: x, dx
        integer :: l

        dx = (x1 - x0) / div_no
        l2 = f(x0, a, b, n, a_sol, fm_val, basic) + f(x1, a, b, n, a_sol, fm_val, basic)
        x = x0

        do l = 1, div_no - 1
            x = x + dx
            l2 = l2 + 2.0 * f(x, a, b, n, a_sol, fm_val, basic)
        end do
        l2 = sqrt(dx * l2 / 2.0)

    end function l2_norm
    ! -----------------------------------------------------------------------------------------------------------------
    ! Квадрат невязки
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function residual2(x, a, b, n, aa, f, basic) result(sum)

        implicit none

        real(dp), external :: basic, f

        real(dp), intent(in) :: x, a, b, aa(1:n)
        integer, intent(in) :: n

        integer :: i, j

        sum = 0._dp
        do i = 1, n
            sum = sum + aa(i) * (a * basic(x, i - 1, 2) + (1 + b * x * x) * basic(x, i - 1, 0))
        end do

        sum = (sum - f(x)) ** 2

    end function residual2
    ! -----------------------------------------------------------------------------------------------------------------
    ! f1(x) - правая часть уравнения
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function f1(x)

        implicit none
        real(dp), intent(in) :: x

        f1 = -1

    end function f1
    ! -----------------------------------------------------------------------------------------------------------------
    ! f2(x) - правая часть уравнения
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function f2(x)

        implicit none
        real(dp), intent(in) :: x

        f2 = x

    end function f2

end module base
