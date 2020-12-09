module finite_difference

    use ogpf
    use base
    use tridiagonal_matrix
    use qr

    implicit none

    private

    type(gpf) :: gp

    public finite_difference_method

    integer, parameter:: dp = kind(0.d0)

contains
    ! -----------------------------------------------------------------------------------------------------------------
    ! Метод разностных аппроксимаций
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine finite_difference_method(AA, BB, x_min, x_max, n, f, a, p, q, alpha1, beta1, alpha2, beta2, &
            y_coef, y_coef_0, y_coef_n, x_grid, sol, is_print, is_draw, info)

        implicit none

        real(dp), external :: f, a, p, q, y_coef, y_coef_0, y_coef_n

        character(*), intent(in) :: info
        integer, intent(in) :: n
        real(dp), intent(in) :: x_min, x_max, AA, BB, alpha1, beta1, alpha2, beta2
        logical, intent(in) :: is_print, is_draw

        real(dp), intent(out) :: sol(:), x_grid(1:n+1)

        real(dp), dimension(1:n+1) :: fm_vector, main_diag, low_diag, up_diag
        real(dp), dimension(1:n+1, 1:n+1) :: coef_matrix

        logical :: is_tridiagonal
        integer :: i, j, error
        real(dp) :: h

        sol = 0._dp; x_grid = 0._dp; coef_matrix = 0._dp; fm_vector = 0._dp
        low_diag = 0._dp; main_diag = 0._dp; up_diag = 0._dp; is_tridiagonal = .true.

        ! шаг
        h = (x_max - x_min) / n

        ! сетка x
        call noised_linspace(from=x_min, to=x_max, array=x_grid, noise=0._dp)

        ! иниициализация вектора свободных членов
        fm_vector(1) = AA
        fm_vector(n+1) = BB
        do i = 2, n
            fm_vector(i) = f(x_grid(i))
        end do

        ! инициализация матрицы коэффициентов
        call init_linear_system(x_grid, h, a, p, q, alpha1, beta1, alpha2, beta2, y_coef, y_coef_0, y_coef_n, &
                coef_matrix, is_tridiagonal)

        if (is_tridiagonal) then
            ! инициализация трех диагналей
            call init_tridiagonal(x_grid, h, a, p, q, alpha1, beta1, alpha2, beta2, y_coef, y_coef_0, y_coef_n, &
                    main_diag, up_diag, low_diag)
        end if

        ! решение системы линейных уравнений
        if (is_tridiagonal) then
            call tridiagonal_solve(main_diag, up_diag, low_diag, fm_vector, sol, error)
            if (error > 0) then
                sol = 0._dp
                call qr_solve(n+1, n+1, coef_matrix, fm_vector, sol)
            end if
        else
            call qr_solve(n+1, n+1, coef_matrix, fm_vector, sol)
        end if


        if (is_print) then

            print *, ''
            print *, info

            print "(/, a2)", "h:"
            print *, h

            print "(/, a7)", "X grid:"
            print *, (x_grid(i), i = 1, n+1)

            print "(/, a19)", "Coefficient matrix:"
            do i = 1, n+1
                print *, (coef_matrix(i, j), j = 1, n+1)
            end do

            if (is_tridiagonal) then
                print "(/, a9)", "Diagonals"
                print "(a5)", "Upper"
                print *, (up_diag(i), i = 1, n)
                print "(a4)", "Main"
                print *, (main_diag(i), i = 1, n+1)
                print "(a5)", "Lower"
                print *, (low_diag(i), i = 2, n+1)
            end if

            print "(/, a20)", "Free members vector:"
            print *, (fm_vector(i), i = 1, n+1)

            if (error > 0 .and. is_tridiagonal) then
                print "(/, a54)", "TRIDIAGONAL_SOLVE: Matrix is degenerate. Use QR_SOLVE."
            end if

            print "(/, a9)", "Solution:"
            print *, (sol(i), i = 1, n+1)

        end if

        if (is_draw) then
            call gp%title(info)
            call gp%xlabel('X grid')
            call gp%ylabel('Y solution')
            call gp%plot(x_grid, sol, 'with lines')
        end if

    end subroutine
    ! -----------------------------------------------------------------------------------------------------------------
    ! Инициализация системы линейных уравнений (для решения системы обычными методами)
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine init_linear_system(x_grid, h, a, p, q, alpha1, beta1, alpha2, beta2, y_coef, y_coef_0, y_coef_n, &
            matrix, is_tridiagonal)

        implicit none

        real(dp), external :: a, p, q, y_coef, y_coef_0, y_coef_n
        real(dp), intent(in) :: h, alpha1, beta1, alpha2, beta2, x_grid(:)
        logical, intent(out) :: is_tridiagonal
        real(dp), intent(out) :: matrix(:,:)

        integer :: i, j, t, n, lvl
        real(dp) :: x, temp

        t = 1; n = size(x_grid) - 1; matrix = 0._dp

        matrix(1,1) = y_coef_0(h, alpha1, beta1, 0)
        matrix(1,2) = y_coef_0(h, alpha1, beta1, 1)
        temp = y_coef_0(h, alpha1, beta1, 2)
        if (temp /= 0) then
            is_tridiagonal = .false.
        end if
        matrix(1,3) = temp

        temp = y_coef_n(h, alpha2, beta2, 0)
        if (temp /= 0) then
            is_tridiagonal = .false.
        end if
        matrix(n+1,n-1) = temp
        matrix(n+1,n)   = y_coef_n(h, alpha2, beta2, 1)
        matrix(n+1,n+1) = y_coef_n(h, alpha2, beta2, 2)

        do i = 2, n
            x = x_grid(i); lvl = 0
            do j = t, t + 2
                matrix(i, j) = y_coef(x, h, a, p, q, lvl)
                lvl = lvl + 1
            end do
            t = t + 1
        end do

        return

    end subroutine init_linear_system
    ! -----------------------------------------------------------------------------------------------------------------
    ! Инициализация трех векторов (главной, верхней и нижней диагоналей)
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine init_tridiagonal(x_grid, h, a, p, q, alpha1, beta1, alpha2, beta2, y_coef, y_coef_0, y_coef_n, &
            main_diag, up_diag, low_diag)

        implicit none

        real(dp), external :: a, p, q, y_coef, y_coef_0, y_coef_n
        real(dp), intent(in) :: h, alpha1, beta1, alpha2, beta2, x_grid(:)
        real(dp), intent(out) :: main_diag(:), up_diag(:), low_diag(:)

        integer :: i, j, n, lvl
        real(dp) :: x

        low_diag = 0._dp; main_diag = 0._dp; up_diag = 0._dp
        n = size(x_grid) - 1

        ! главная диагональ
        lvl = 1
        main_diag(1)   = y_coef_0(h, alpha1, beta1, 0)
        main_diag(n+1) = y_coef_n(h, alpha2, beta2, 2)
        do i = 2, n
            x = x_grid(i)
            main_diag(i) = y_coef(x, h, a, p, q, lvl)
        end do

        ! нижняя диагональ
        lvl = 0
        low_diag(n+1) = y_coef_n(h, alpha2, beta2, 1)
        do i = 2, n
            x = x_grid(i)
            low_diag(i) = y_coef(x, h, a, p, q, lvl)
        end do

        ! верхняя диагональ
        lvl = 2
        up_diag(1) = y_coef_0(h, alpha1, beta1, 1)
        do i = 2, n
            x = x_grid(i)
            up_diag(i) = y_coef(x, h, a, p, q, lvl)
        end do

        return

    end subroutine init_tridiagonal

end module finite_difference