module tridiagonal_matrix

    implicit none

    private

    public tridiagonal_solve

    integer, parameter:: dp = kind(0.d0)

contains
    ! -----------------------------------------------------------------------------------------------------------------
    ! Метод прогонки для решения систем линейных уравнений Mx=D, где M - тридиагональная матрица
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine tridiagonal_solve(A, B, C, D, X, error)

        implicit none

        real(dp), intent(in) :: A(:), B(:), C(:), D(:)
        real(dp), intent(out) :: X(:)
        integer, intent(out) :: error
        real(dp), allocatable:: P(:), Q(:), R(:), S(:), T(:)
        real(dp) :: W
        integer ::  i, ii, n

        n = size(A)
        if (n < 3) then
            error = 1; return
        end if

        allocate(P(n), Q(n), R(n), S(n), T(n))

        P(1)=0._dp; Q(1)=0._dp; R(1)=1._dp
        do i = 1, n - 1
            ii = i + 1
            W = A(i) + Q(i) * C(i)
            if (1._dp + W == 1._dp) then
                deallocate(P, Q, R, S, T)
                error = 65; return
            end if
            P(ii) = (D(i) - P(i) * C(i)) / W
            Q(ii) = -B(i) / W
            R(ii) = -R(i) * C(i) / W
        end do
        S(n) = 1._dp; T(n) = 0._dp
        do i = n - 1, 1, -1
            ii = i + 1
            S(i) = Q(ii) * S(ii) + R(ii)
            T(i) = Q(ii) * T(ii) + P(ii)
        end do
        W = A(n) + B(n) * S(1) + c(n) * S(n-1)
        if (1._dp + W == 1._dp) then
            deallocate(P, Q, R, S, T)
            error = 65; return
        end if
        X(n) = (D(n) - B(n) * T(1) - C(n) * T(n - 1)) / W
        do i = 1, n - 1
            X(i) = S(i) * X(n) + T(i)
        end do
        deallocate(P, Q, R, S, T)
        error=0; return
    end subroutine tridiagonal_solve

end module tridiagonal_matrix