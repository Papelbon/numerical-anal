module gauss

    implicit none

    private

    public gauss_solve

    integer, parameter:: dp = kind(0.d0)

contains
    ! -----------------------------------------------------------------------------------------------------------------
    ! Метод Гаусса решения системы линейных уравнений A*x=b
    ! -----------------------------------------------------------------------------------------------------------------
    subroutine gauss_solve(matrix, free_members, x, n, is_singular)

        implicit none

        real(dp) :: a(1:n, 1:n), b(1:n), s(1:n)

        real(dp), intent(in) :: matrix(:, :), free_members(:)
        integer, intent(in) :: n

        real(dp), intent(out) :: x(:)
        logical, intent(out) :: is_singular

        real(dp) :: c, pivot, store
        integer :: i, j, k, l

        is_singular = .false.
        a = matrix; b = free_members

        ! step 1: begin forward elimination
        do k=1, n-1
            ! step 2: "scaling"
            ! s(i) will have the largest element from row i
            do i=k,n ! loop over rows
                s(i) = 0.0
                do j=k,n ! loop over elements of row i
                    s(i) = max(s(i),abs(a(i,j)))
                end do
            end do
            ! step 3: "pivoting 1"
            ! find a row with the largest pivoting element
            pivot = abs(a(k,k)/s(k))
            l=k
            do j=k+1,n
                if(abs(a(j,k)/s(j)) > pivot) then
                    pivot = abs(a(j,k)/s(j))
                    l=j
                end if
            end do
            ! Check if the system has a singular matrix
            if(pivot == 0.0) then
                is_singular = .true.
                return
            end if
            ! step 4: "pivoting 2" interchange rows k and l (if needed)
            if (l /= k) then
                do j=k,n
                    store = a(k,j)
                    a(k,j) = a(l,j)
                    a(l,j) = store
                end do
                store = b(k)
                b(k) = b(l)
                b(l) = store
            end if
            ! step 5: the elimination (after scaling and pivoting)
            do i=k+1,n
                c=a(i,k)/a(k,k)
                a(i,k) = 0.0
                b(i)=b(i)- c*b(k)
                do j=k+1,n
                    a(i,j) = a(i,j)-c*a(k,j)
                end do
            end do
        end do
        ! step 6: back substiturion
        x(n) = b(n)/a(n,n)
        do i=n-1,1,-1
            c=0.0
            do j=i+1,n
                c= c + a(i,j)*x(j)
            end do
            x(i) = (b(i)- c)/a(i,i)
        end do
    end subroutine gauss_solve

end module gauss