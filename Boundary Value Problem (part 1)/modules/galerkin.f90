
module galerkin

    use base
    use ogpf
    use lu
    use qr
    use integrate

    implicit none

    private

    type(gpf) :: gp

    public galerkin_method

    integer, parameter:: dp = kind(0.d0)
contains
    ! -----------------------------------------------------------------------------------------------------------------
    ! Метод Галеркина
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine galerkin_method(a, b, x_min, x_max, n, fm_val, solution, solution_ext, &
            x_points, is_test, is_print, is_draw, info)

        implicit none

        character *16, intent(in) :: info
        integer, intent(in) :: n, x_points
        real(dp), intent(in) :: a, b, x_min, x_max
        logical, intent(in) :: is_test, is_print, is_draw

        real(dp), intent(out) :: solution(:), solution_ext(:)

        real(dp), dimension(1:x_points) :: x_list
        real(dp), dimension(1:n) :: free_members_vector
        real(dp), dimension(1:n, 1:n) :: coef_matrix

        real(dp), external :: fm_val

        integer :: i, j, error
        real(dp) :: l2

        solution = 0._dp; x_list = 0._dp; solution_ext = 0._dp; free_members_vector = 0._dp; coef_matrix = 0._dp

        ! иниициализация вектора свободных членов
        do i = 1, n
            free_members_vector(i) = trapezium(fm_galerkin_func, x_min, x_max, i, j, a, b, fm_val, 100, is_test)
        end do

        ! инициализация матрицы коэффициентов
        do i = 1, n
            do j = 1, n
                coef_matrix(i, j) = trapezium(int_galerkin_func, x_min, x_max, i, j, a, b, fm_val, 100, is_test)
            end do
        end do

        ! решение системы линейных уравнений
        call lu_solve(coef_matrix, free_members_vector, solution, error)
        if (error > 0) then
            call qr_solve(n, n, coef_matrix, free_members_vector, solution)
        end if

        ! численное решение
        call noised_linspace(x_min, x_max, x_list, noise=0._dp)
        do i = 1, x_points ! заменить на test_galerkin_basic_func
            do j = 1, size(solution)
                if (is_test) then
                    solution_ext(i) = solution_ext(i) + solution(j) * test_basic_func(x_list(i), j - 1, 0)
                else
                    solution_ext(i) = solution_ext(i) + solution(j) * basic_func(x_list(i), j - 1, 0)
                end if
            end do
        end do

        ! вычисление погрешности
        if (is_test) then ! заменить на test_galerkin_basic_func
            l2 = l2_norm(residual2, x_min, x_max, n, a, b, solution, fm_val, 100, test_basic_func)
        else
            l2 = l2_norm(residual2, x_min, x_max, n, a, b, solution, fm_val, 100, basic_func)
        end if

        if (is_print) then

            print "(/, a19)", "Coefficient matrix:"
            do i = 1, n
                print *, (coef_matrix(i, j), j = 1, n)
            end do

            print "(/, a20)", "Free members vector:"
            print *, (free_members_vector(i), i = 1, n)

            if (error > 0) then
                print "(/, a42)", "LU_SOLVE: Matrix is special. Use QR_SOLVE."
            end if

            print "(/, a23)", "Linear system solution:"
            print *, (solution(i), i = 1, n)

            print "(/, a21)", "Collocation solution:"
            print *, (solution_ext(i), i = 1, x_points)

            print "(/, a8)", "L2 norm:"
            print *, l2

        end if

        if (is_draw) then
            call gp%title(info // " Galerkin method")
            call gp%xlabel('-1 <= X <= 1')
            call gp%ylabel('Yn(x) = φ0(x) + a1 * φ1(x) + ... + an * φn(x)')
            call gp%plot(x_list, solution_ext, 'with lines')
        end if

    end subroutine galerkin_method
    ! -----------------------------------------------------------------------------------------------------------------
    ! Вспомогательные функции для метода Галеркина
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function int_galerkin_func(x, i, j, a, b, fm_val, is_test) result(res)

        implicit none

        real(dp), external :: fm_val
        real(dp), intent(in) :: x, a, b
        integer, intent(in) :: i, j
        logical, intent(in) :: is_test

        if (is_test) then  ! заменить на test_galerkin_basic_func
            res = (a * test_basic_func(x, j - 1, 2) + (1 + b * x * x) * &
                    test_basic_func(x, j - 1, 0))  * test_basic_func(x, i - 1, 0)
        else
            res = (a * basic_func(x, j - 1, 2) + (1 + b * x * x) * basic_func(x, j - 1, 0)) * basic_func(x, i - 1, 0)
        end if

    end function int_galerkin_func

    real(dp) function fm_galerkin_func(x, i, j, a, b, fm_val, is_test) result(res)

        implicit none

        real(dp), external :: fm_val
        real(dp), intent(in) :: x, a, b
        integer, intent(in) :: i, j
        logical, intent(in) :: is_test

        if (is_test) then ! заменить на test_gakerkin_basic_func
            res = fm_val(x) * (1 * test_basic_func(x, i - 1, 0))
        else
            res = fm_val(x) * (1 * basic_func(x, i - 1, 0))
        end if

    end function fm_galerkin_func
end module galerkin