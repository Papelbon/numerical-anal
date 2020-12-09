
module lu

    implicit none

    private

    public lu_solve

    integer, parameter:: dp = kind(0.d0)

contains
    ! -----------------------------------------------------------------------------------------------------------------
    ! Решение системы линейных уравнений: метод разложения LU
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine lu_solve(coef_matrix, free_members, solution, error)

        implicit none

        real(dp), intent(in) :: free_members(:), coef_matrix(:,:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: error

        integer :: i, j, m, n
        integer, dimension(1:size(free_members)) :: permutation_array
        real(dp), dimension(1:size(free_members), 1:size(free_members)) :: matrix_lu
        real(dp) :: sum, det

        call matrix_decomposition(coef_matrix, matrix_lu, permutation_array, det, error)
        if (error > 0) then
            return
        end if

        n = size(matrix_lu, 1)
        if (loc(free_members) /= loc(solution)) solution = free_members
        do i = 1, n
            m = permutation_array(i)
            sum = dble(solution(m))
            solution(m) = solution(i)
            do j = 1, i - 1
                sum = sum - matrix_lu(i, j) * dble(solution(j))
            end do
            solution(i) = sngl(sum)
        end do
        do i = n, 1, -1
            sum = dble(solution(i))
            do j = i + 1, n
                sum = sum - matrix_lu(i, j) * dble(solution(j))
            end do
            solution(i) = sngl(sum / matrix_lu(i, i))
        end do
        return
    end subroutine lu_solve
    ! -----------------------------------------------------------------------------------------------------------------
    ! Разложение матрицы LU
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine matrix_decomposition(matrix, matrix_lu, permutation_array, det, error)

        implicit none

        real(dp), intent(in):: matrix(:,:)
        real(dp), intent(out) :: matrix_lu(:,:), det
        integer, intent(out):: permutation_array(:), error
        real(dp) :: aa, ww, rn, sum
        integer i, j, k, p, max, n

        n = size(matrix, 1)
        rn = float(n)
        det = 1.0
        error = 0
        if (loc(matrix) /= loc(matrix_lu)) matrix_lu = matrix
        do p = 1, n-1
            if (p > 1) then
                do i = p, N
                    sum = dble(matrix_lu(i, p))
                    do k = 1, p - 1
                        sum = sum - matrix_lu(i, k) * dble(matrix_lu(k, p))
                    end do
                    matrix_lu(i, p) = sngl(sum)
                end do
            end if
            aa = abs(matrix_lu(p, p))
            max = p
            do i = p + 1, n
                if (aa < abs(matrix_lu(i, p))) then
                    aa = abs(matrix_lu(i, p))
                    max = i
                end if
            end do
            if (rn + aa == rn) then
                det = 0.0
                error = 65
                return
            end if
            if (max /= p) then
                do k = 1, N
                    aa = matrix_lu(max, k)
                    matrix_lu(max, k) = matrix_lu(p, k)
                    matrix_lu(p, k) = aa
                end do
                det = -det
            end if
            det = det * matrix_lu(p, p)
            permutation_array(p) = max
            do i = p + 1, n
                matrix_lu(i, p) = matrix_lu(i, p) / matrix_lu(p, p)
            end do
            if (p > 1) then
                do j = p + 1, n
                    sum = dble(matrix_lu(p, j))
                    do k = 1, p - 1
                        sum = sum - matrix_lu(p, k) * dble(matrix_lu(k, j))
                    end do
                    matrix_lu(p, j) = sngl(sum)
                end do
            end if
        end do
        sum = dble(matrix_lu(n, n))
        do k = 1, n - 1
            sum = sum - matrix_lu(n, k) * dble(matrix_lu(k, n))
        end do
        matrix_lu(n, n) = sngl(sum)
        det = det * sngl(sum)
        permutation_array(n) = n
        return
    end subroutine matrix_decomposition
end module lu