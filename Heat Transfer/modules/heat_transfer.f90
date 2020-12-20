module heat_transfer

    use ogpf
    use base

    use tridiagonal_matrix
    use qr
    use integrate

    implicit none

    private

    type(gpf) :: gp

    public fd_centered_heat_transfer_steady, fd_balance_heat_transfer_steady, &
            fd_heat_transfer_explicit, fd_heat_transfer_implicit

    integer, parameter:: dp = kind(0.d0)

contains
    ! -----------------------------------------------------------------------------------------------------------------
    ! Метод конечных разностей (приближение центрированной разности)
    ! для численного решения стационарного (не зависящего от времени) уравнения теплопроводности в одномерной области:
    ! -d/dx (k(x) du/dx) = f(x),
    ! u(a) = Ua, u(b) = Ub.
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine fd_centered_heat_transfer_steady(n, a, b, Ua, Ub, k, f, x_grid, sol, is_print, is_draw, info)

        implicit none

        real(dp), external :: k, f

        character(*), intent(in) :: info
        integer, intent(in) :: n
        real(dp), intent(in) :: a, b, Ua, Ub
        logical, intent(in) :: is_print, is_draw

        real(dp), intent(out) :: x_grid(1:n+1), sol(:)

        real(dp), dimension(1:n+1) :: fm_vector, main_diag, low_diag, up_diag

        integer :: i, error
        real(dp) :: h, xm, xp

        sol = 0._dp; x_grid = 0._dp; fm_vector = 0._dp; low_diag = 0._dp; main_diag = 0._dp; up_diag = 0._dp

        ! шаг
        h = (b - a) / n

        ! сетка x
        call noised_linspace(from=a, to=b, array=x_grid, noise=0._dp)

        ! инициализация вектора свободных членов и трех диагоналей
        main_diag(1) = 1._dp
        fm_vector(1) = Ua
        do i = 2, n

            xm = (x_grid(i-1) + x_grid(i)) / 2._dp
            xp = (x_grid(i) + x_grid(i+1)) / 2._dp

            low_diag(i)  = -k(xm) / h / h
            main_diag(i) = (k(xm) + k(xp)) / h / h
            up_diag(i)   = -k(xp) / h / h

            fm_vector(i) = f(x_grid(i))

        end do
        fm_vector(n+1) = Ub
        main_diag(n+1) = 1._dp

        call tridiagonal_solve(main_diag, up_diag, low_diag, fm_vector, sol, error)
        if (error == 1) then
            print "(/, a42)", "TRIDIAGONAL_SOLVE: Matrix dimension N < 3."; return
        elseif (error == 65) then
            print "(/, a40)", "TRIDIAGONAL_SOLVE: Matrix is degenerate."; return
        end if


        if (is_print) then

            print *, ''
            print *, info

            print "(/, a2)", "h:"
            print *, h

            print "(/, a7)", "X grid:"
            print *, (x_grid(i), i = 1, n+1)

            print "(/, a9)", "Diagonals"
            print "(a5)", "Upper"
            print *, (up_diag(i), i = 1, n)
            print "(a4)", "Main"
            print *, (main_diag(i), i = 1, n+1)
            print "(a5)", "Lower"
            print *, (low_diag(i), i = 2, n+1)

            print "(/, a20)", "Free members vector:"
            print *, (fm_vector(i), i = 1, n+1)

            print "(/, a9)", "Solution:"
            print *, (sol(i), i = 1, n+1)

        end if

        if (is_draw) then
            call gp%title(info)
            call gp%xlabel('a <= X <= b')
            call gp%ylabel('Temperature')
            call gp%plot(x_grid, sol, 'with lines')
        end if

    end subroutine
    ! -----------------------------------------------------------------------------------------------------------------
    ! Метод конечных разностей (метод баланса)
    ! для численного решения стационарного (не зависящего от времени) уравнения теплопроводности в одномерной области:
    ! -d/dx (k(x) du/dx) = f(x),
    ! u(a) = Ua, u(b) = Ub.
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine fd_balance_heat_transfer_steady(n, a, b, Ua, Ub, k, f, x_points, sol, is_print, is_draw, info, psources_)

        implicit none

        real(dp), external :: k, f

        character(*), intent(in) :: info
        integer, intent(in) :: n
        real(dp), intent(in) :: a, b, Ua, Ub
        logical, intent(in) :: is_print, is_draw

        real(dp), intent(out) :: x_points(1:n+1), sol(:)

        ! опциональный параметр: psources_ = (x01, c1, x02, c2, ...)
        ! x0 - положение источника, c - мощность источника
        real(dp), optional, intent(in) :: psources_(:)

        real(dp), dimension(1:n+1) :: fm_vector, main_diag, low_diag, up_diag

        integer :: i, j, error
        real(dp) :: h, xm, xp

        sol = 0._dp; x_points = 0._dp; fm_vector = 0._dp; low_diag = 0._dp; main_diag = 0._dp; up_diag = 0._dp

        ! шаг
        h = (b - a) / n

        ! сетка x
        call noised_linspace(from=a, to=b, array=x_points, noise=0._dp)

        ! иниициализация вектора свободных членов и трех диагоналей
        main_diag(1) = 1._dp
        fm_vector(1) = Ua
        do i = 2, n

            low_diag(i)   = c(x_points(i-1), x_points(i))
            main_diag(i)  = -c(x_points(i-1), x_points(i)) - c(x_points(i), x_points(i+1))
            up_diag(i)    = c(x_points(i), x_points(i+1))

            if (present(psources_)) then
                do j = 1, size(psources_)-1, 2
                    fm_vector(i) = fm_vector(i) - get_psource_val(x_points(i), psources_(j), psources_(j+1), h)
                end do
            else
                fm_vector(i)  = -phi(x_points(i) - h/2, x_points(i) + h/2)
            end if

        end do
        fm_vector(n+1) = Ub
        main_diag(n+1) = 1._dp

        call tridiagonal_solve(main_diag, up_diag, low_diag, fm_vector, sol, error)
        if (error == 1) then
            print "(/, a42)", "TRIDIAGONAL_SOLVE: Matrix dimension N < 3."; return
        elseif (error == 65) then
            print "(/, a40)", "TRIDIAGONAL_SOLVE: Matrix is degenerate."; return
        end if


        if (is_print) then

            print *, ''
            print *, info

            print "(/, a2)", "h:"
            print *, h

            print "(/, a9)", "X points:"
            print *, (x_points(i), i = 1, n+1)

            print "(/, a9)", "Diagonals"
            print "(a5)", "Upper"
            print *, (up_diag(i), i = 1, n)
            print "(a4)", "Main"
            print *, (main_diag(i), i = 1, n+1)
            print "(a5)", "Lower"
            print *, (low_diag(i), i = 2, n+1)

            print "(/, a20)", "Free members vector:"
            print *, (fm_vector(i), i = 1, n+1)

            print "(/, a9)", "Solution:"
            print *, (sol(i), i = 1, n+1)

        end if

        if (is_draw) then
            call gp%title(info)
            call gp%xlabel('a <= X <= b')
            call gp%ylabel('Temperature')
            call gp%plot(x_points, sol, 'with lines')
        end if

    contains
        real(dp) function q(x)
            
            implicit none
            real(dp), intent(in) :: x

            q = 1 / k(x)

        end function q

        real(dp) function c(a, b)

            implicit none
            real(dp), intent(in) :: a, b

            c = 1 / trapezium(q, a, b, 200)

        end function c

        real(dp) function phi(a, b)

            implicit none
            real(dp), intent(in) :: a, b

            phi = trapezium(f, a, b, 200)

        end function phi

    end subroutine

    real(dp) function get_psource_val(x, x0, c, h) result(fm)

        implicit none
        real(dp), intent(in) :: x, x0, c, h

        if (abs(x - x0) - h/2 < 1e-5) then
            fm = c / 2
        elseif (abs(x - x0) < h/2) then
            fm = c
        else
            fm = 0._dp
        end if

    end function get_psource_val
    ! -----------------------------------------------------------------------------------------------------------------
    ! Метод конечных разностей (явная схема)
    ! для численного решения зависящего от времени уравнения теплопроводности в одномерной области:
    ! du/dt - d/dx (k(x) du/dx) = f(x,t), a < x < b, 0 < t <= T.
    !
    ! bc=0: u(a,t) = g1(t), u(b,t) = g2(t);
    ! bc=1: u(a,t) = g1(t), du/dt(b,t) = g2(t)  -  левая разность;
    ! bc=2: u(a,t) = g1(t), du/dt(b,t) = g2(t)  -  центральная разность.
    !
    ! u(x,0) = φ(x), a <= x <= b.
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine fd_heat_transfer_explicit(n, m, a, b, t, g1, g2, phi, k, f, x_list, t_list, matrix, is_check_cfl, bc, &
            is_print, is_draw, info)

        implicit none

        real(dp), external :: k, f, g1, g2, phi

        character(*), intent(in) :: info
        integer, intent(in) :: n, m, bc
        real(dp), intent(in) :: a, b, t
        logical, intent(in) :: is_print, is_draw, is_check_cfl

        real(dp), intent(out) :: x_list(:), t_list(:), matrix(:,:)

        integer :: i, j
        real(dp) :: dx, dt, cfl

        x_list = 0._dp; t_list = 0._dp; matrix = 0._dp

        ! шаг x
        dx = (b - a) / n

        ! сетка x
        call noised_linspace(from=a, to=b, array=x_list, noise=0._dp)

        ! шаг t
        dt = (t - 0) / m

        ! сетка t
        call noised_linspace(from=0._dp, to=t, array=t_list, noise=0._dp)

        if (is_check_cfl) then
            cfl = k(1) * dt / dx / dx
            if (cfl >= 0.5_dp) then
                print *, ''
                print *, info
                print *, 'fd_heat_transfer_explicit - fatal error!'
                print *, 'CFL condition failed'
                print *, '0.5 <= k * dt / dx / dx = CFL'
                return
            end if
        end if

        ! численное решение
        do i = 1, m+1
            matrix(i,1) = g1(t_list(i))
        end do
        do i = 1, n+1
            matrix(1,i) = phi(x_list(i))
        end do

        do i = 1, m
            do j = 2, n
                matrix(i+1,j) = dt * k(x_list(j) - dx / 2) / dx / dx * matrix(i,j-1) + &
                    (1 - dt * (k(x_list(j) - dx / 2) + k(x_list(j) + dx / 2)) / dx / dx) * matrix(i,j) + &
                    dt * k(x_list(j) + dx / 2) / dx / dx * matrix(i,j+1) + &
                    dt * f(x_list(j), t_list(i))
            end do
            if (bc == 0) then
                matrix(i+1,n+1) = g2(t_list(i))
            elseif (bc == 1) then
                matrix(i+1,n+1) = matrix(i+1,n) + dx * g2(t_list(i))
            elseif (bc == 2) then
                ! matrix(i+1,n+1) = matrix(i+1,n-1) + 2 * dx * g2(t_list(i))
                matrix(i+1,n+1) = dt * k(x_list(j) - dx / 2) / dx / dx * matrix(i,n) + &
                        (1 - dt * (k(x_list(j) - dx / 2) + k(x_list(j) + dx / 2)) / dx / dx) * matrix(i,n+1) + &
                        dt * k(x_list(j) + dx / 2) / dx / dx * (2 * dx * g2(t_list(i)) + matrix(i+1,n)) + &
                        dt * f(x_list(j), t_list(i))
            end if

        end do

        if (is_print) then

            print *, ''
            print *, info

            print "(/, a3)", "dx:"
            print *, dx

            print "(/, a3)", "dt:"
            print *, dt

            if (is_check_cfl) then
                print "(/, a4)", "cfl:"
                print *, cfl
            end if

            print "(/, a7)", "X list:"
            print *, (x_list(i), i = 1, n+1)

            print "(/, a7)", "T list:"
            print *, (t_list(i), i = 1, m+1)

            print "(/, a9)", "Solution:"
            do i = 1, m+1
                print *, (matrix(i,j), j = 1, n+1)
            end do

        end if

        if (is_draw) then
            call draw_animation(x_list, t_list, matrix, info)
        end if

    end subroutine
    ! -----------------------------------------------------------------------------------------------------------------
    ! Метод конечных разностей (неявная схема)
    ! для численного решения зависящего от времени уравнения теплопроводности в одномерной области:
    ! du/dt - d/dx (k(x) du/dx) = f(x,t), a < x < b, 0 < t <= T.
    !
    ! bc=1: u(a,t) = g1(t), du/dt(b,t) = g2(t)  -  левая разность;
    ! bc=2: u(a,t) = g1(t), du/dt(b,t) = g2(t)  -  центральная разность.
    ! else: u(a,t) = g1(t), u(b,t) = g2(t);
    !
    ! u(x,0) = φ(x), a <= x <= b.
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine fd_heat_transfer_implicit(n, m, a, b, t, g1, g2, phi, k, f, x_list, t_list, matrix, bc, &
            is_print, is_draw, info)

        implicit none

        real(dp), external :: f, g1, g2, phi

        character(*), intent(in) :: info
        integer, intent(in) :: n, m, bc
        real(dp), intent(in) :: a, b, t, k
        logical, intent(in) :: is_print, is_draw

        real(dp), intent(out) :: x_list(:), t_list(:), matrix(:,:)
        real(dp), allocatable :: temp_matrix(:,:), fm_vector(:), temp_result(:)

        integer :: i, j, nn
        real(dp) :: dx, dt

        x_list = 0._dp; t_list = 0._dp; matrix = 0._dp

        ! шаг x
        dx = (b - a) / n

        ! сетка x
        call noised_linspace(from=a, to=b, array=x_list, noise=0._dp)

        ! шаг t
        dt = (t - 0) / m

        ! сетка t
        call noised_linspace(from=0._dp, to=t, array=t_list, noise=0._dp)

        ! численное решение
        do i = 1, n+1
            matrix(1,i) = phi(x_list(i))
        end do
        do i = 1, m+1
            matrix(i,1) = g1(t_list(i))
        end do

        nn = n + 1
        if (bc == 2) nn = n + 2

        allocate(temp_matrix(nn,nn), fm_vector(nn), temp_result(nn))
        temp_matrix = 0._dp; fm_vector = 0._dp; temp_result = 0._dp

        temp_matrix(1,1) = 1._dp
        do i = 2, nn - 1
            temp_matrix(i,i-1) = -k * dt / dx / dx
            temp_matrix(i,i)   = 1 + 2 * k * dt / dx / dx
            temp_matrix(i,i+1) = -k * dt / dx / dx
        end do
        temp_matrix(nn,nn) = 1._dp

        if (bc == 1) temp_matrix(nn,nn-1) = -1._dp
        if (bc == 2) temp_matrix(nn,nn-2) = -1._dp

        do i = 2, m+1
            fm_vector(1) = g1(t_list(i))
            do j = 2, nn-1
                fm_vector(j) = dt * f(x_list(j), t_list(i) + dt) + matrix(i-1,j)
            end do
            fm_vector(nn) = g2(t_list(i))
            if (bc == 1) fm_vector(nn) = dx * fm_vector(nn)
            if (bc == 2) fm_vector(nn) = 2 * dx * fm_vector(nn)

            call qr_solve(nn, nn, temp_matrix, fm_vector, temp_result)
            if (bc == 2) then
                matrix(i,:) = temp_result(1:nn-1)
            else
                matrix(i,:) = temp_result(1:nn)
            end if
            temp_result = 0._dp
        end do
        deallocate(temp_matrix, fm_vector, temp_result)

        if (is_print) then

            print *, ''
            print *, info

            print "(/, a3)", "dx:"
            print *, dx

            print "(/, a3)", "dt:"
            print *, dt

            print "(/, a7)", "X list:"
            print *, (x_list(i), i = 1, n+1)

            print "(/, a7)", "T list:"
            print *, (t_list(i), i = 1, m+1)

            print "(/, a9)", "Solution:"
            do i = 1, m+1
                print *, (matrix(i,j), j = 1, n+1)
            end do

        end if

        if (is_draw) then
            call draw_animation(x_list, t_list, matrix, info)
        end if

    end subroutine

    subroutine draw_animation(x_list, t_list, matrix, info)

        implicit none

        character(*), intent(in) :: info
        real(dp), allocatable :: xx(:,:), tt(:,:)
        real(dp), intent(in) :: t_list(:), x_list(:), matrix(:,:)
        integer :: i

        ! 3D анимация
        call meshgrid(xx, tt, x_list, t_list)
        call gp%title(info)

        call gp%options('unset colorbox')
        call gp%options('set ticslevel 0')
        call gp%xlabel('--TIME--')
        call gp%ylabel('--X--')
        call gp%zlabel('--U(X,T)--')

        call gp%animation_start(1)
        do i = 1, size(t_list), 10
            call gp%surf(tt(1:i,:), xx(1:i,:), matrix(1:i,:), palette='jet')
        end do
        call gp%animation_show()

        ! 2D анимация
        call gp%xlabel('X')
        call gp%ylabel('U=TEMP')

        call gp%animation_start(1)
        do i = 1, size(t_list), 10
            call gp%plot(x_list, matrix(i,:), 'with lines lt 8')
        end do
        call gp%animation_show()

    end subroutine draw_animation

end module heat_transfer