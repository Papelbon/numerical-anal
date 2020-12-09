
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
    subroutine galerkin_method(x_min, x_max, n, solution, solution_ext, f, a, p, q, basic,&
            x_points, is_print, is_draw, info)

        implicit none

        character(*), intent(in) :: info
        integer, intent(in) :: n, x_points
        real(dp), intent(in) :: x_min, x_max
        logical, intent(in) :: is_print, is_draw

        real(dp), intent(out) :: solution(:), solution_ext(:)

        real(dp), dimension(1:x_points) :: x_list
        real(dp), dimension(1:n) :: free_members_vector
        real(dp), dimension(1:n, 1:n) :: coef_matrix

        real(dp), external :: f, a, p, q, basic

        integer :: i, j, error
        real(dp) :: l2

        solution = 0._dp; x_list = 0._dp; solution_ext = 0._dp; free_members_vector = 0._dp; coef_matrix = 0._dp

        ! иниициализация вектора свободных членов
        do i = 1, n
            free_members_vector(i) = trapezium(fm_galerkin_func, x_min, x_max, i, j, f, a, p, q, basic, 100)
        end do

        ! инициализация матрицы коэффициентов
        do i = 1, n
            do j = 1, n
                coef_matrix(i, j) = trapezium(int_galerkin_func, x_min, x_max, i, j, f, a, p, q, basic, 100)
            end do
        end do

        ! решение системы линейных уравнений
        call lu_solve(coef_matrix, free_members_vector, solution, error)
        if (error > 0) then
            call qr_solve(n, n, coef_matrix, free_members_vector, solution)
        end if

        ! численное решение
        call noised_linspace(x_min, x_max, x_list, noise=0._dp)
        do i = 1, x_points
            solution_ext(i) = solution_ext(i) + basic(x_list(i), 0, 0)
            do j = 2, size(solution)
                solution_ext(i) = solution_ext(i) + solution(j) * basic(x_list(i), j-1, 0)
            end do
        end do

        ! вычисление погрешности
        l2 = l2_norm(residual2, x_min, x_max, n, solution, f, a, p, q, basic, 100)

        if (is_print) then

            print *, ''
            print *, info

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
            call gp%title(info)
            call gp%xlabel('X')
            call gp%ylabel('Yn(x) = φ0(x) + a1 * φ1(x) + ... + an * φn(x)')
            call gp%plot(x_list, solution_ext, 'with lines')
        end if

    end subroutine galerkin_method
    ! -----------------------------------------------------------------------------------------------------------------
    ! Вспомогательные функции для метода Галеркина
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function int_galerkin_func(x, i, j, f, a, p, q, basic) result(res)

        implicit none

        real(dp), external :: f, a, p, q, basic
        real(dp), intent(in) :: x
        integer, intent(in) :: i, j

        res = basic(x, i-1, 0) * (a(x) * basic(x, j-1, 2) + p(x) * basic(x, j-1, 1) + q(x) * basic(x, j-1, 0))

    end function int_galerkin_func

    real(dp) function fm_galerkin_func(x, i, j, f, a, p, q, basic) result(res)

        implicit none

        real(dp), external :: f, a, p, q, basic
        real(dp), intent(in) :: x
        integer, intent(in) :: i, j

        res = f(x) * basic(x, i-1, 0)

    end function fm_galerkin_func
end module galerkin