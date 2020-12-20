! Модуль содержит основные функции и сабрутины

module base

    implicit none

    private

    public noised_linspace, std

    integer, parameter:: dp = kind(0.d0)

contains
    ! -----------------------------------------------------------------------------------------------------------------
    ! Стандтартное отклонение
    ! -----------------------------------------------------------------------------------------------------------------
    real(dp) function std(data)

        implicit none

        real(dp), intent(in) :: data(:)
        real(dp) :: mean, variance
        integer :: i, n

        n = size(data)

        mean = 0._dp
        do i = 1, n
            mean = mean + data(i)
        end do
        mean = mean / n

        variance = 0._dp
        do i = 1, n
            variance = variance + (data(i) - mean) ** 2
        end do
        variance = variance / (n-1)

        std = sqrt(variance)
    end
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

end module base
