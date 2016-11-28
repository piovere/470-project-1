subroutine mat(width, sigma_r, D, matrix)
    use iso_fortran_env
    implicit none

    integer :: N
    real(real64), intent(in) :: sigma_r, D, width
    real(real64), intent(inout) :: matrix(:,:)
    real(real64), allocatable :: ipiv(:)
    integer :: i
    integer :: info
    real(real64) :: dx

    external DGETRF

    N = size(matrix, dim=1)
    allocate (ipiv(N))

    dx = width / real(N,real64)

    ! Top row of matrix
    matrix(1,1) = (0.5 * sigma_r) + (D / dx**2)
    matrix(1,2) = -D / dx**2

    ! Bottom row of matrix
    matrix(N,N) = sigma_r + 2.0 * D / dx**2
    matrix(N,N-1) = -D / dx**2

    ! All the other rows
    do i=2, N-1
        matrix(i,i) = sigma_r + 2.0 * D / dx**2
        matrix(i,i+1) = -D / dx**2
        matrix(i,i-1) = -D / dx**2
    enddo
    
    call DGETRF(N, N, matrix, N, ipiv, info)
    if (info /= 0) stop 'Matrix is numerically singular!'
end subroutine mat