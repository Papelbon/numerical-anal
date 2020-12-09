
module collocation

    use base
    use lu
    use qr
    use ogpf

    implicit none

    private

    type(gpf) :: gp

    public collocation_method

    integer, parameter:: dp = kind(0.d0)

contains
    ! -----------------------------------------------------------------------------------------------------------------
    ! Метод коллокаций
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine collocation_method(x_min, x_max, n, solution, solution_ext, f, a, p, q, basic, &
            x_points, is_print, is_draw, info)

        implicit none

        character(*), intent(in) :: info
        integer, intent(in) :: n, x_points
        real(dp), intent(in) :: x_min, x_max
        logical, intent(in) :: is_print, is_draw

        real(dp), intent(out) :: solution(:), solution_ext(:)

        real(dp), dimension(1:x_points) :: x_list
        real(dp), dimension(1:n) :: collocation_points, free_members_vector
        real(dp), dimension(1:n, 1:n) :: coef_matrix

        integer :: i, j, error
        real(dp) :: l2

        real(dp), external :: f, a, p, q, basic

        solution = 0._dp; x_list = 0._dp; solution_ext = 0._dp
        free_members_vector = 0._dp; coef_matrix = 0._dp; collocation_points = 0._dp

        ! точки коллокации
        call noised_linspace(from=x_min, to=x_max, array=collocation_points, noise=0.01_dp)

        ! иниициализация вектора свободных членов
        do i = 1, n
            free_members_vector(i) = f(collocation_points(i))
        end do

        ! инициализация матрицы коэффициентов
        call init_linear_system(a, p, q, collocation_points, coef_matrix, basic)

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

            print "(/, a19)", "Collocation points:"
            print *, (collocation_points(i), i = 1, n)

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
            call gp%plot(x_list, solution_ext, 'with lines ')
        end if

    end subroutine
    ! -----------------------------------------------------------------------------------------------------------------
    ! Инициализация системы линейных уравнений
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine init_linear_system(a, p, q, col_points, matrix, basic)

        implicit none

        real(dp), intent(in) :: col_points(:)
        real(dp), external :: a, p, q, basic
        real(dp), intent(out) :: matrix(:,:)

        integer :: i, j, n
        real(dp) :: x

        n = size(col_points)
        do i = 1, n
            x = col_points(i)
            do j = 1, n
                matrix(i, j) = a(x) * basic(x, j-1, 2) + p(x) * basic(x, j-1, 1) + q(x) * basic(x, j-1, 0)
            end do
        end do

        return

    end subroutine init_linear_system

end module collocation