function matrix(nodes, sigma_a, sigma_inscatter, sigma_outscatter)

implicit none
integer                     :: nodes
real(real64), allocatable   :: A(:,:)
real(real64), allocatable   :: ipiv(:,:)
real(real64)                :: sigma_a, sigma_inscatter, sigma_outscatter
real(real64)                :: D
real(real64)                :: w
real(real64)                :: dx, dx2 

end function matrix

! subroutine matrix
! ! Neutron Flux desperation Group 3
! ! DGETRS - solve a system of linear equations A * X = B or 
! ! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html

! use iso_fortran_env
! implicit none

! real(real64), allocatable :: A(:,:)                 ! matrix to LU decompose
! real(real64), allocatable :: b(:)                   ! b vector
! real(real64), allocatable :: ipiv(:)                ! pivoting array for LAPACK

! !--------------------------Matrix Declaration-----------------------------!

! A(1,1) = (0.5*sigma_a)+(D/(dx**2))
! A(1,2) = -D/(dx**2.0)
! A(n,n) = sigma_a+2.0*D/(dx**2)
! A(n,n-1) = -D/(dx**2)

! do i=2 , n-1

!     A(i,i) = sigma_a+2.0*D/(dx**2)
!     A(i,i+1) = -D/(dx**2)
!     A(i,i-1) = -D/(dx**2)

! enddo

! ! Pretty matrix printing
! ! From http://jblevins.org/log/array-write
! ! do i=1,n
! !     write(*,"(100g15.5)") ( A(i,j), j=1,n )
! ! enddo

! ! LU factorization of a general M-by-N matrix A
! ! See: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
! call DGETRF(n, n, A, n, ipiv, info)
! if (info /= 0) stop 'Matrix is numerically singular!'

! end subroutine matrix