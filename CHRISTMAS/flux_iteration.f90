subroutine iterate(matrix, S, flux, k, flux_error, k_error)
    use iso_fortran_env
    implicit none

    integer :: N ! size of matrix
    integer :: info ! error code from DGETRS
    real(real64), intent(inout) :: matrix(:,:)
    real(real64), intent(inout) :: k
    real(real64), intent(inout) :: flux(:), S(:)
    real(real64), intent(out) :: flux_error, k_error
    real(real64), allocatable :: ipiv(:)
    real(real64) :: m, m_old, k_old
    real(real64), allocatable :: flux_old(:)
    
    external DGETRS

    N = size(matrix, dim=1)
    allocate (ipiv(N), flux_old(N))

    flux_old = flux
    k_old = k

    flux = S * flux

    ! solve for the new flux
    call DGETRS('N', N, 1, matrix(:,:), N, ipiv(:), flux(:), N, info)
    if (info /= 0) stop 'Solution of the linear system failed'

    ! calculate convergence
    m = norm2(flux)
    m_old = norm2(flux_old)

    flux = flux / m

    flux_error = maxval(abs(flux - flux_old) / flux)
    k_error = abs(k - k_old) / k
end subroutine iterate
